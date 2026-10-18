// The collision operator of the thermal-velocity normalized PKPM solver as a
// self-contained update: the conservative LBO on [F, G] (drag and diffusion of F with the
// boundary corrections, the G source 2 nu vth^2 F - 2 nu G folded into the drag kernel)
// with primitive moments from the raw normalized moments (mass 1) must preserve the
// three constraints M0(F) = 1, M1(F) = 0 and M2(F) + 2 M0(G) = 3 of any normalized state
// to round-off, i.e. the moments [M0, M1, M2 + 2 M0_G] of its right-hand side vanish
// cell by cell. The state here is anisotropic and non-Maxwellian with a varying
// collisionality so that every term of the operator is active.
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_lbo_pkpm_diff.h>
#include <gkyl_dg_lbo_pkpm_drag.h>
#include <gkyl_dg_updater_lbo_pkpm.h>
#include <gkyl_dg_updater_moment_pkpm.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

struct lbo_ctx {
  int cdim;
  double nu_amp; // relative amplitude of the collisionality variation
};

// A two-stream, anisotropic normalized state: F a weighted pair of shifted Gaussians,
// G = theta(x) F with theta != 1 (perpendicular temperature differs from the parallel one).
static void
eval_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct lbo_ctx *c = ctx;
  double x = xn[0], w = xn[c->cdim];
  double a = 0.7 + 0.2 * sin(M_PI * x), s1 = 0.8, s2 = 1.3;
  double F = a * exp(-0.5 * (w - 0.6) * (w - 0.6) / (s1 * s1)) / (s1 * sqrt(2.0 * M_PI))
    + (1.0 - a) * exp(-0.5 * (w + 1.1) * (w + 1.1) / (s2 * s2)) / (s2 * sqrt(2.0 * M_PI));
  double theta = 0.6 + 0.3 * cos(M_PI * x);
  fout[0] = F;
  fout[1] = theta * F;
}

static void
eval_nu(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct lbo_ctx *c = ctx;
  fout[0] = 1.0 + c->nu_amp * sin(M_PI * xn[0]);
}

static void
test_lbo_constraints(int cdim, double nu_amp)
{
  struct lbo_ctx ctx = {.cdim = cdim, .nu_amp = nu_amp};
  // The velocity domain holds the tails of the state to below round-off: with the wider
  // Gaussian reaching 1e-6 at the boundary (wmax = 8) the boundary treatment of the
  // operator leaves a 1e-10 relative residual that decreases with the cell size.
  double wmax = 12.0;
  double lower[4] = {-1.0, -1.0, -1.0, 0.0}, upper[4] = {1.0, 1.0, 1.0, 0.0};
  int cells[4] = {4, 3, 2, 0}, ghost[4] = {1, 1, 1, 0};
  lower[cdim] = -wmax;
  upper[cdim] = wmax;
  cells[cdim] = 48;
  ghost[cdim] = 0;
  struct gkyl_rect_grid cgrid, pgrid, vgrid;
  gkyl_rect_grid_init(&cgrid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&pgrid, cdim + 1, lower, upper, cells);
  gkyl_rect_grid_init(&vgrid, 1, &lower[cdim], &upper[cdim], &cells[cdim]);
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext, vlocal, vlocal_ext;
  gkyl_create_grid_ranges(&cgrid, ghost, &clocal_ext, &clocal);
  gkyl_create_grid_ranges(&pgrid, ghost, &plocal_ext, &plocal);
  gkyl_create_grid_ranges(&vgrid, &ghost[cdim], &vlocal_ext, &vlocal);
  int nc = cbasis.num_basis, np = pbasis.num_basis;
  int vdim = 1;

  struct gkyl_array *f = mkarr(2 * np, plocal_ext.volume), *rhs = mkarr(2 * np, plocal_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, plocal_ext.volume);
  struct gkyl_array *nu = mkarr(nc, clocal_ext.volume);
  struct gkyl_array *nmoms = mkarr(6 * nc, clocal_ext.volume), *rhs_moms = mkarr(6 * nc, clocal_ext.volume);
  struct gkyl_array *bcorr = mkarr((vdim + 1) * nc, clocal_ext.volume);
  struct gkyl_array *prim = mkarr((vdim + 1) * nc, clocal_ext.volume);
  struct gkyl_array *nu_prim = mkarr((vdim + 1) * nc, clocal_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&pgrid, &pbasis, 8, 2, eval_f, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &plocal_ext, f);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 1, eval_nu, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, nu);
  gkyl_proj_on_basis_release(proj);

  // Raw normalized moments [M0, M2, M0_G, M1, M3, M1_G] (mass 1), the boundary
  // corrections, the primitive moments [u_par, vth^2] and nu times them.
  struct gkyl_dg_updater_moment *mcalc =
    gkyl_dg_updater_moment_pkpm_new(&pgrid, &cbasis, &pbasis, &clocal_ext, &vlocal, 1.0, false, false);
  gkyl_array_clear(nmoms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &plocal, &clocal, f, nmoms);
  double v_bounds[2] = {-wmax, wmax};
  struct gkyl_mom_calc_bcorr *bcorr_calc =
    gkyl_mom_calc_bcorr_lbo_pkpm_new(&pgrid, &cbasis, &pbasis, v_bounds, 1.0, false);
  gkyl_array_clear(bcorr, 0.0);
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &plocal, &clocal, f, bcorr);
  struct gkyl_prim_lbo_calc *pcalc = gkyl_prim_lbo_pkpm_calc_new(&pgrid, &cbasis, &pbasis, &clocal, false);
  gkyl_prim_lbo_calc_advance(pcalc, &clocal, nmoms, bcorr, nu, prim);
  for (int d = 0; d < vdim + 1; ++d) {
    gkyl_dg_mul_op(&cbasis, d, nu_prim, d, prim, 0, nu);
  }

  struct gkyl_dg_lbo_pkpm_drag_auxfields drag_inp = {.nuSum = nu, .nuPrimMomsSum = nu_prim};
  struct gkyl_dg_lbo_pkpm_diff_auxfields diff_inp = {.nuSum = nu, .nuPrimMomsSum = nu_prim};
  struct gkyl_dg_updater_collisions *lbo =
    gkyl_dg_updater_lbo_pkpm_new(&pgrid, &cbasis, &pbasis, &clocal, &drag_inp, &diff_inp, false);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_updater_lbo_pkpm_advance(lbo, &plocal, f, cflrate, rhs);

  gkyl_array_clear(rhs_moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &plocal, &clocal, rhs, rhs_moms);

  // The operator is active: the second moment of F and the zeroth of G change at a rate
  // comparable to nu; their weighted sum and the lower moments do not.
  double scale = 0.0, err = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &clocal);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&clocal_ext, iter.idx);
    const double *m = gkyl_array_cfetch(rhs_moms, loc);
    for (int k = 0; k < nc; ++k) {
      scale = fmax(scale, fmax(fabs(m[nc + k]), fabs(m[2 * nc + k])));
      err = fmax(err, fabs(m[k]));
      err = fmax(err, fabs(m[3 * nc + k]));
      err = fmax(err, fabs(m[nc + k] + 2.0 * m[2 * nc + k]));
    }
  }
  TEST_CHECK(scale > 0.1);
  TEST_MSG("cdim=%d: the collision operator is inactive (scale %g)", cdim, scale);
  TEST_CHECK(err < 1e-11 * scale);
  TEST_MSG("cdim=%d nu_amp=%g: constraint moments of the LBO right-hand side %g (scale %g)", cdim, nu_amp, err, scale);

  gkyl_dg_updater_lbo_pkpm_release(lbo);
  gkyl_prim_lbo_calc_release(pcalc);
  gkyl_mom_calc_bcorr_release(bcorr_calc);
  gkyl_dg_updater_moment_pkpm_release(mcalc);
  gkyl_array_release(nu_prim);
  gkyl_array_release(prim);
  gkyl_array_release(bcorr);
  gkyl_array_release(rhs_moms);
  gkyl_array_release(nmoms);
  gkyl_array_release(nu);
  gkyl_array_release(cflrate);
  gkyl_array_release(rhs);
  gkyl_array_release(f);
}

void test_lbo_constraints_1x() { test_lbo_constraints(1, 0.0); test_lbo_constraints(1, 0.5); }
void test_lbo_constraints_2x() { test_lbo_constraints(2, 0.0); test_lbo_constraints(2, 0.5); }

TEST_LIST = {
  {"test_lbo_constraints_1x", test_lbo_constraints_1x},
  {"test_lbo_constraints_2x", test_lbo_constraints_2x},
  {NULL, NULL},
};
