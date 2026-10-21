// Tests of the exponential Runge-Kutta update of the first-harmonic blocks (dg_pkpm_harm_exp),
// driven through the three stages exactly as the SSP-RK3 driver of the app does: the nodal
// transforms are exact inverses, the update is exact for a constant field and forcing at any
// Omega dt, it is SSP-RK3 itself without a field, and it is third order with a field that changes
// magnitude and direction across the step.
#include <math.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_dg_pkpm_harm_exp_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// One phase cell of 8 blocks (the first-harmonic blocks are 2-7) and its configuration cell.
struct hx_cell {
  int cdim, np, nc;
  struct gkyl_basis cbasis;
  struct gkyl_range crange, prange;
  struct gkyl_dg_pkpm_harm_exp *up;
  struct gkyl_array *em, *coef, *S, *A;
};

static struct hx_cell
hx_cell_new(int cdim)
{
  struct hx_cell c = {.cdim = cdim};
  gkyl_cart_modal_serendip(&c.cbasis, cdim, 1);
  c.nc = c.cbasis.num_basis;
  c.np = 3 * c.nc;
  int lo[GKYL_MAX_DIM], hi[GKYL_MAX_DIM];
  for (int d = 0; d < cdim + 1; ++d) {
    lo[d] = hi[d] = 1;
  }
  gkyl_range_init(&c.crange, cdim, lo, hi);
  gkyl_range_init(&c.prange, cdim + 1, lo, hi);
  c.up = gkyl_dg_pkpm_harm_exp_new(&(struct gkyl_dg_pkpm_harm_exp_inp) {
    .conf_basis = &c.cbasis, .qbym = 1.0, .use_gpu = false});
  c.em = gkyl_array_new(GKYL_DOUBLE, 8 * c.nc, 1);
  c.coef = gkyl_array_new(GKYL_DOUBLE, gkyl_dg_pkpm_harm_exp_num_coef(c.up), 1);
  c.S = gkyl_array_new(GKYL_DOUBLE, gkyl_dg_pkpm_harm_exp_num_stage(c.up), 1);
  c.A = gkyl_array_new(GKYL_DOUBLE, gkyl_dg_pkpm_harm_exp_num_stage(c.up), 1);
  return c;
}

static void
hx_cell_release(struct hx_cell *c)
{
  gkyl_dg_pkpm_harm_exp_release(c->up);
  gkyl_array_release(c->em);
  gkyl_array_release(c->coef);
  gkyl_array_release(c->S);
  gkyl_array_release(c->A);
}

// A spatially uniform field B (only the constant mode of the configuration basis).
static void
hx_set_field(struct hx_cell *c, const double *B)
{
  double *em = gkyl_array_fetch(c->em, 0);
  for (int k = 0; k < 8 * c->nc; ++k) {
    em[k] = 0.0;
  }
  for (int a = 0; a < 3; ++a) {
    em[(3 + a) * c->nc] = B[a] * pow(sqrt(2.0), c->cdim);
  }
}

// The test problem, mode by mode of the first-harmonic blocks: dU/dt = U x omega(t) + S_k(t). A
// uniform omega makes the Galerkin product with U the cross product of the coefficient vectors.
struct hx_problem {
  double om0; // magnitude scale of omega
  bool vary; // omega and S change in time
};

static void
hx_omega(const struct hx_problem *pb, double t, double *om)
{
  if (!pb->vary) {
    double n = sqrt(0.29);
    om[0] = pb->om0 * 0.3 / n;
    om[1] = pb->om0 * -0.2 / n;
    om[2] = pb->om0 * 0.4 / n;
    return;
  }
  double ang = 0.3 * sin(0.7 * t), mag = pb->om0 * (1.0 + 0.2 * sin(0.7 * t + 1.0));
  om[0] = mag * sin(ang);
  om[1] = 0.0;
  om[2] = mag * cos(ang);
}

static double
hx_source(const struct hx_problem *pb, double t, int k, int a)
{
  double base = 0.3 + 0.1 * k - 0.2 * a;
  return pb->vary ? base * (1.0 + 0.3 * sin(0.5 * t + k + a)) : base;
}

// rhs = U x omega(t) + S(t) in the first-harmonic blocks; F and G get zero.
static void
hx_rhs(const struct hx_cell *c, const struct hx_problem *pb, double t, const double *f, double *rhs)
{
  double om[3];
  hx_omega(pb, t, om);
  for (int k = 0; k < 8 * c->np; ++k) {
    rhs[k] = 0.0;
  }
  for (int v = 0; v < 2; ++v) {
    for (int k = 0; k < c->np; ++k) {
      double u[3];
      for (int a = 0; a < 3; ++a) {
        u[a] = f[(2 + 3 * v + a) * c->np + k];
      }
      double uxo[3] = {u[1] * om[2] - u[2] * om[1], u[2] * om[0] - u[0] * om[2], u[0] * om[1] - u[1] * om[0]};
      for (int a = 0; a < 3; ++a) {
        rhs[(2 + 3 * v + a) * c->np + k] = uxo[a] + hx_source(pb, t, k + 3 * v, a);
      }
    }
  }
}

// One step of the app's SSP-RK3 driver (Shu-Osher form, forward Euler stages) with the
// exponential update written into the first-harmonic blocks of each stage output.
static void
hx_step(struct hx_cell *c, const struct hx_problem *pb, double t, double h, bool use_exp, double *f)
{
  int n = 8 * c->np;
  double f1[8 * 24], fnew[8 * 24], rhs[8 * 24];
  struct gkyl_array *fa = gkyl_array_new(GKYL_DOUBLE, n, 1), *f1a = gkyl_array_new(GKYL_DOUBLE, n, 1);
  struct gkyl_array *ra = gkyl_array_new(GKYL_DOUBLE, n, 1), *oa = gkyl_array_new(GKYL_DOUBLE, n, 1);
  double om[3];
  hx_omega(pb, t, om);
  hx_set_field(c, om); // qbym = 1
  gkyl_dg_pkpm_harm_exp_coef(c->up, &c->crange, h, c->em, c->coef);
  for (int k = 0; k < n; ++k) {
    ((double *)gkyl_array_fetch(fa, 0))[k] = f[k];
  }
  // stage 1
  hx_rhs(c, pb, t, f, rhs);
  for (int k = 0; k < n; ++k) {
    f1[k] = f[k] + h * rhs[k];
  }
  if (use_exp) {
    for (int k = 0; k < n; ++k) {
      ((double *)gkyl_array_fetch(ra, 0))[k] = rhs[k];
      ((double *)gkyl_array_fetch(oa, 0))[k] = f1[k];
    }
    gkyl_dg_pkpm_harm_exp_accum(c->up, 1, h, &c->crange, &c->prange, c->coef, fa, ra, c->S, c->A);
    gkyl_dg_pkpm_harm_exp_state(c->up, 1, h, &c->crange, &c->prange, c->coef, fa, c->S, c->A, oa);
    for (int k = 0; k < n; ++k) {
      f1[k] = ((double *)gkyl_array_fetch(oa, 0))[k];
    }
  }
  // stage 2
  hx_rhs(c, pb, t + h, f1, rhs);
  for (int k = 0; k < n; ++k) {
    fnew[k] = f1[k] + h * rhs[k];
  }
  if (use_exp) {
    for (int k = 0; k < n; ++k) {
      ((double *)gkyl_array_fetch(f1a, 0))[k] = f1[k];
      ((double *)gkyl_array_fetch(ra, 0))[k] = rhs[k];
      ((double *)gkyl_array_fetch(oa, 0))[k] = fnew[k];
    }
    gkyl_dg_pkpm_harm_exp_accum(c->up, 2, h, &c->crange, &c->prange, c->coef, f1a, ra, c->S, c->A);
    gkyl_dg_pkpm_harm_exp_state(c->up, 2, h, &c->crange, &c->prange, c->coef, fa, c->S, c->A, oa);
    for (int k = 0; k < n; ++k) {
      fnew[k] = ((double *)gkyl_array_fetch(oa, 0))[k];
    }
  }
  for (int k = 0; k < n; ++k) {
    f1[k] = 0.75 * f[k] + 0.25 * fnew[k];
  }
  // stage 3
  hx_rhs(c, pb, t + 0.5 * h, f1, rhs);
  for (int k = 0; k < n; ++k) {
    fnew[k] = f1[k] + h * rhs[k];
  }
  if (use_exp) {
    for (int k = 0; k < n; ++k) {
      ((double *)gkyl_array_fetch(f1a, 0))[k] = f1[k];
      ((double *)gkyl_array_fetch(ra, 0))[k] = rhs[k];
      ((double *)gkyl_array_fetch(oa, 0))[k] = fnew[k];
    }
    gkyl_dg_pkpm_harm_exp_accum(c->up, 3, h, &c->crange, &c->prange, c->coef, f1a, ra, c->S, c->A);
    gkyl_dg_pkpm_harm_exp_state(c->up, 3, h, &c->crange, &c->prange, c->coef, fa, c->S, c->A, oa);
    for (int k = 0; k < n; ++k) {
      fnew[k] = ((double *)gkyl_array_fetch(oa, 0))[k];
    }
  }
  for (int k = 0; k < n; ++k) {
    f[k] = f[k] / 3.0 + 2.0 / 3.0 * fnew[k];
  }
  gkyl_array_release(fa);
  gkyl_array_release(f1a);
  gkyl_array_release(ra);
  gkyl_array_release(oa);
}

static void
hx_init(const struct hx_cell *c, double *f)
{
  for (int k = 0; k < 8 * c->np; ++k) {
    f[k] = 0.1 * sin(1.3 * k + 0.4);
  }
}

static void
hx_run(struct hx_cell *c, const struct hx_problem *pb, double T, int nstep, bool use_exp, double *f)
{
  hx_init(c, f);
  double h = T / nstep;
  for (int s = 0; s < nstep; ++s) {
    hx_step(c, pb, s * h, h, use_exp, f);
  }
}

// Classical RK4 reference of the same ODE.
static void
hx_ref(struct hx_cell *c, const struct hx_problem *pb, double T, int nstep, double *f)
{
  int n = 8 * c->np;
  double k1[8 * 24], k2[8 * 24], k3[8 * 24], k4[8 * 24], tmp[8 * 24];
  hx_init(c, f);
  double h = T / nstep;
  for (int s = 0; s < nstep; ++s) {
    double t = s * h;
    hx_rhs(c, pb, t, f, k1);
    for (int k = 0; k < n; ++k) tmp[k] = f[k] + 0.5 * h * k1[k];
    hx_rhs(c, pb, t + 0.5 * h, tmp, k2);
    for (int k = 0; k < n; ++k) tmp[k] = f[k] + 0.5 * h * k2[k];
    hx_rhs(c, pb, t + 0.5 * h, tmp, k3);
    for (int k = 0; k < n; ++k) tmp[k] = f[k] + h * k3[k];
    hx_rhs(c, pb, t + h, tmp, k4);
    for (int k = 0; k < n; ++k) f[k] += h / 6.0 * (k1[k] + 2.0 * k2[k] + 2.0 * k3[k] + k4[k]);
  }
}

// Exact solution with constant omega and S: U_par + S_par t, d + e^(tL)(U - d) with d = (S x om)/|om|^2.
static void
hx_exact_const(struct hx_cell *c, const struct hx_problem *pb, double T, double *f)
{
  double om[3];
  hx_omega(pb, 0.0, om);
  double mag = sqrt(om[0] * om[0] + om[1] * om[1] + om[2] * om[2]);
  double oh[3] = {om[0] / mag, om[1] / mag, om[2] / mag};
  hx_init(c, f);
  for (int v = 0; v < 2; ++v) {
    for (int k = 0; k < c->np; ++k) {
      double u[3], s[3];
      for (int a = 0; a < 3; ++a) {
        u[a] = f[(2 + 3 * v + a) * c->np + k];
        s[a] = hx_source(pb, 0.0, k + 3 * v, a);
      }
      double ub = u[0] * oh[0] + u[1] * oh[1] + u[2] * oh[2], sb = s[0] * oh[0] + s[1] * oh[1] + s[2] * oh[2];
      double sp[3], d[3], w[3];
      for (int a = 0; a < 3; ++a) sp[a] = s[a] - sb * oh[a];
      d[0] = (sp[1] * om[2] - sp[2] * om[1]) / (mag * mag);
      d[1] = (sp[2] * om[0] - sp[0] * om[2]) / (mag * mag);
      d[2] = (sp[0] * om[1] - sp[1] * om[0]) / (mag * mag);
      for (int a = 0; a < 3; ++a) w[a] = u[a] - ub * oh[a] - d[a];
      double wxo[3] = {w[1] * oh[2] - w[2] * oh[1], w[2] * oh[0] - w[0] * oh[2], w[0] * oh[1] - w[1] * oh[0]};
      double cs = cos(mag * T), sn = sin(mag * T);
      for (int a = 0; a < 3; ++a) {
        f[(2 + 3 * v + a) * c->np + k] = (ub + sb * T) * oh[a] + d[a] + cs * w[a] + sn * wxo[a];
      }
    }
  }
}

static double
hx_rel_err(const struct hx_cell *c, const double *f, const double *g)
{
  double num = 0.0, den = 0.0;
  for (int k = 2 * c->np; k < 8 * c->np; ++k) {
    num += (f[k] - g[k]) * (f[k] - g[k]);
    den += g[k] * g[k];
  }
  return sqrt(num / den);
}

static void
test_transforms(int cdim)
{
  struct hx_cell c = hx_cell_new(cdim);
  double f[8 * 24], g[8 * 24], fn[24 * 6];
  for (int k = 0; k < 8 * c.np; ++k) {
    f[k] = sin(0.7 * k + 0.1);
    g[k] = -1.0;
  }
  pkpm_hx_to_nodes(cdim, f, fn);
  pkpm_hx_to_modal(cdim, fn, g);
  for (int k = 0; k < 2 * c.np; ++k) {
    TEST_CHECK(g[k] == -1.0); // F and G untouched
  }
  double err = 0.0;
  for (int k = 2 * c.np; k < 8 * c.np; ++k) {
    err = fmax(err, fabs(g[k] - f[k]));
  }
  TEST_CHECK(err < 1e-14);
  TEST_MSG("cdim %d: round trip error %g", cdim, err);
  hx_cell_release(&c);
}

static void test_transforms_1x(void) { test_transforms(1); }
static void test_transforms_2x(void) { test_transforms(2); }
static void test_transforms_3x(void) { test_transforms(3); }

// Constant field and forcing: exact at any Omega dt (the drift, the rotation and the parallel growth).
static void
test_constant(int cdim)
{
  struct hx_cell c = hx_cell_new(cdim);
  double f[8 * 24], g[8 * 24];
  double om_h[] = {0.1, 2.0, 10.0, 100.0, 1000.0};
  for (int i = 0; i < 5; ++i) {
    struct hx_problem pb = {.om0 = om_h[i] / 0.1, .vary = false};
    hx_run(&c, &pb, 2.0, 20, true, f);
    hx_exact_const(&c, &pb, 2.0, g);
    double err = hx_rel_err(&c, f, g);
    TEST_CHECK(err < 1e-12);
    TEST_MSG("cdim %d, Omega dt = %g: relative error %g", cdim, om_h[i], err);
  }
  hx_cell_release(&c);
}

static void test_constant_1x(void) { test_constant(1); }
static void test_constant_2x(void) { test_constant(2); }
static void test_constant_3x(void) { test_constant(3); }

// Without a field the update is SSP-RK3 (to round-off).
static void
test_no_field(void)
{
  struct hx_cell c = hx_cell_new(2);
  double f[8 * 24], g[8 * 24];
  struct hx_problem pb = {.om0 = 0.0, .vary = true};
  hx_run(&c, &pb, 2.0, 20, true, f);
  hx_run(&c, &pb, 2.0, 20, false, g);
  double err = hx_rel_err(&c, f, g);
  TEST_CHECK(err < 1e-14);
  TEST_MSG("relative difference to SSP-RK3 %g", err);
  hx_cell_release(&c);
}

// Zero first harmonic stays exactly zero.
static void
test_zero(void)
{
  struct hx_cell c = hx_cell_new(2);
  struct hx_problem pb = {.om0 = 50.0, .vary = true};
  double f[8 * 24], rhs[8 * 24];
  for (int k = 0; k < 8 * c.np; ++k) {
    f[k] = k < 2 * c.np ? 1.0 + 0.1 * k : 0.0;
  }
  // a zero right-hand side in the first-harmonic blocks: no source
  struct gkyl_array *fa = gkyl_array_new(GKYL_DOUBLE, 8 * c.np, 1);
  struct gkyl_array *ra = gkyl_array_new(GKYL_DOUBLE, 8 * c.np, 1);
  struct gkyl_array *oa = gkyl_array_new(GKYL_DOUBLE, 8 * c.np, 1);
  double B[3] = {1.0, 2.0, 3.0};
  hx_set_field(&c, B);
  gkyl_dg_pkpm_harm_exp_coef(c.up, &c.crange, 0.3, c.em, c.coef);
  for (int k = 0; k < 8 * c.np; ++k) {
    ((double *)gkyl_array_fetch(fa, 0))[k] = f[k];
    ((double *)gkyl_array_fetch(ra, 0))[k] = k < 2 * c.np ? 0.5 : 0.0;
    ((double *)gkyl_array_fetch(oa, 0))[k] = 7.0;
  }
  bool ok = true;
  for (int s = 1; s <= 3; ++s) {
    gkyl_dg_pkpm_harm_exp_accum(c.up, s, 0.3, &c.crange, &c.prange, c.coef, fa, ra, c.S, c.A);
    gkyl_dg_pkpm_harm_exp_state(c.up, s, 0.3, &c.crange, &c.prange, c.coef, fa, c.S, c.A, oa);
    const double *o = gkyl_array_cfetch(oa, 0);
    for (int k = 0; k < 8 * c.np; ++k) {
      ok = ok && (k < 2 * c.np ? o[k] == 7.0 : o[k] == 0.0);
    }
  }
  TEST_CHECK(ok);
  gkyl_array_release(fa);
  gkyl_array_release(ra);
  gkyl_array_release(oa);
  hx_cell_release(&c);
}

// Field changing magnitude (20%) and direction (0.3 rad): third order with the frozen
// exponential and the stage remainder.
static void
test_order(void)
{
  struct hx_cell c = hx_cell_new(1);
  double f[8 * 24], g[8 * 24];
  struct hx_problem pb = {.om0 = 4.0, .vary = true};
  hx_ref(&c, &pb, 3.0, 60000, g);
  double err[4];
  int nsteps[4] = {15, 30, 60, 120};
  for (int i = 0; i < 4; ++i) {
    hx_run(&c, &pb, 3.0, nsteps[i], true, f);
    err[i] = hx_rel_err(&c, f, g);
  }
  for (int i = 1; i < 4; ++i) {
    double order = log2(err[i - 1] / err[i]);
    TEST_CHECK(order > 2.8);
    TEST_MSG("Omega dt = %g -> %g: errors %g -> %g, order %g", 4.0 * 3.0 / nsteps[i - 1],
      4.0 * 3.0 / nsteps[i], err[i - 1], err[i], order);
  }
  hx_cell_release(&c);
}

TEST_LIST = {
  {"test_transforms_1x", test_transforms_1x},
  {"test_transforms_2x", test_transforms_2x},
  {"test_transforms_3x", test_transforms_3x},
  {"test_constant_1x", test_constant_1x},
  {"test_constant_2x", test_constant_2x},
  {"test_constant_3x", test_constant_3x},
  {"test_no_field", test_no_field},
  {"test_zero", test_zero},
  {"test_order", test_order},
  {NULL, NULL},
};
