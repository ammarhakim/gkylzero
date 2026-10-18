// Tests of the momentum solve of the refactored PKPM solver: the exact modal volume
// kernel against quadrature, the face flux (dg_pkpm_fluid_flux_surf) against an
// independent nodal assembly sharing only the ten-moment Riemann solver, and the
// complete 1x momentum right-hand side against the legacy solver for b along x with
// uniform density, where the two formulations coincide.
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_dg_pkpm_fluid_flux_surf.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_wv_ten_moment_priv.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

struct field_ctx {
  int cdim;
  bool uniform_rho; // for the legacy comparison
  bool b_along_x;
  bool zero_u; // for the legacy comparison
  bool inadmissible; // make p_par negative in part of the domain to exercise the HLL fallback
};

static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct field_ctx *c = ctx;
  int cdim = c->cdim;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  for (int i = 0; i < 8; ++i) {
    fout[i] = 0.0;
  }
  if (c->b_along_x) {
    fout[3] = 2.0;
    return;
  }
  fout[3] = 1.0 + 0.4 * sin(x) * cos(y);
  fout[4] = 0.6 * cos(x) * sin(z) + 0.5 * sin(y);
  fout[5] = -0.5 + 0.3 * sin(y) * cos(z) + 0.2 * sin(x);
}

// Physical moments [rho, p_par, p_perp, M1, q_par, q_perp]
static void
eval_moms(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct field_ctx *c = ctx;
  int cdim = c->cdim;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = c->uniform_rho ? 2.0 : 2.0 + 0.5 * sin(x) * cos(y) + 0.2 * sin(z);
  fout[1] = 1.0 + 0.3 * cos(x) * sin(y);
  fout[2] = 0.7 + 0.25 * cos(y) * sin(x) + 0.1 * cos(z);
  fout[3] = 0.0;
  fout[4] = 0.3 * sin(x) * cos(z) + 0.1 * cos(y);
  fout[5] = -0.2 * cos(x) + 0.15 * sin(y) * sin(z);
  if (c->inadmissible && x > 1.0) {
    fout[1] = -0.2; // negative parallel pressure: P_nn <= 0 where b is along the normal
    fout[2] = 0.05;
  }
}

// Normalized moments [1, M2, theta, 0, M3, M1_G] of the same state with p = (p_par + 2 p_perp)/3
// and xi^2 = p/rho: M2 = p_par/p, theta = p_perp/p, and the third moments scaled by 1/(p xi).
static void
eval_nmoms(double t, const double *xn, double *restrict fout, void *ctx)
{
  double m[6];
  eval_moms(t, xn, m, ctx);
  double p = (m[1] + 2.0 * m[2]) / 3.0, xi = sqrt(fabs(p) / m[0]);
  fout[0] = 1.0;
  fout[1] = m[1] / p;
  fout[2] = m[2] / p;
  fout[3] = 0.0;
  fout[4] = m[4] / (p * xi);
  fout[5] = m[5] / (p * xi);
}

// Thermal speed and scalar pressure [xi, xi^2, p] of the same state.
static void
eval_xi(double t, const double *xn, double *restrict fout, void *ctx)
{
  double m[6];
  eval_moms(t, xn, m, ctx);
  double p = (m[1] + 2.0 * m[2]) / 3.0;
  fout[1] = p / m[0];
  fout[0] = sqrt(fabs(fout[1]));
  fout[2] = p;
}

// Momentum M = rho u for a smooth u that stays positive along x (no stagnation switch).
static void
eval_euler(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct field_ctx *c = ctx;
  int cdim = c->cdim;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double rho = c->uniform_rho ? 2.0 : 2.0 + 0.5 * sin(x) * cos(y) + 0.2 * sin(z);
  fout[0] = rho * (1.5 + 0.3 * sin(x) + 0.2 * cos(y));
  fout[1] = rho * (0.4 * cos(y) * sin(z) + 0.2 * sin(x));
  fout[2] = rho * (-0.3 * sin(z) + 0.25 * cos(x) * cos(y));
  if (c->zero_u) {
    fout[0] = fout[1] = fout[2] = 0.0;
  }
}

// Fluid state [rho, rho ux, rho uy, rho uz, E] of the energy-constrained update.
static void
eval_fluid(double t, const double *xn, double *restrict fout, void *ctx)
{
  double m[6], M[3];
  eval_euler(t, xn, M, ctx);
  eval_moms(t, xn, m, ctx);
  fout[0] = m[0];
  fout[1] = M[0];
  fout[2] = M[1];
  fout[3] = M[2];
  fout[4] = 0.5 * (M[0] * M[0] + M[1] * M[1] + M[2] * M[2]) / m[0] + 0.5 * (m[1] + 2.0 * m[2]);
}

// Kinetic coupling moments [p_par, p_perp, q_par, q_perp] of the same state.
static void
eval_cmoms(double t, const double *xn, double *restrict fout, void *ctx)
{
  double m[6];
  eval_moms(t, xn, m, ctx);
  fout[0] = m[1];
  fout[1] = m[2];
  fout[2] = m[4];
  fout[3] = m[5];
}


static const double gl3_x[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double gl3_w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
static const double gl4_x[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
  0.8611363115940526};
static const double gl4_w[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
  0.3478548451374538};

struct test_setup {
  int cdim;
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
  struct gkyl_array *moms, *cmoms, *euler, *nmoms, *fluid, *bvar, *u, *xi;
  struct gkyl_wv_eqn *ten;
  struct gkyl_wave_geom *geom;
  struct gkyl_dg_calc_pkpm_vars *vars;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
};

static void
setup(struct test_setup *s, struct field_ctx *ctx)
{
  int cdim = ctx->cdim;
  s->cdim = cdim;
  double lower[3] = {-M_PI, -M_PI, -M_PI}, upper[3] = {M_PI, M_PI, M_PI};
  int cells[3] = {8, 4, 3};
  int ghost[3] = {1, 1, 1};
  gkyl_rect_grid_init(&s->grid, cdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->basis, cdim, 1);
  gkyl_create_grid_ranges(&s->grid, ghost, &s->local_ext, &s->local);
  int nc = s->basis.num_basis;

  s->moms = mkarr(6 * nc, s->local_ext.volume);
  s->cmoms = mkarr(4 * nc, s->local_ext.volume);
  s->euler = mkarr(3 * nc, s->local_ext.volume);
  s->nmoms = mkarr(6 * nc, s->local_ext.volume);
  s->fluid = mkarr(5 * nc, s->local_ext.volume);
  s->xi = mkarr(3 * nc, s->local_ext.volume);

  s->bvar = mkarr(3 * nc, s->local_ext.volume);
  s->u = mkarr(3 * nc, s->local_ext.volume);
  struct gkyl_array *em = mkarr(8 * nc, s->local_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 6, eval_moms, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->moms);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 4, eval_cmoms, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->cmoms);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 3, eval_euler, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->euler);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 6, eval_nmoms, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->nmoms);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 5, eval_fluid, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->fluid);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 3, eval_xi, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->xi);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 8, eval_em, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, em);
  gkyl_proj_on_basis_release(proj);

  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->grid, &s->basis, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->local_ext, em, s->bvar);
  gkyl_array_release(em);

  s->ten = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  s->geom = gkyl_wave_geom_new(&s->grid, &s->local_ext, 0, 0, false);
  s->vars = gkyl_dg_calc_pkpm_vars_new(&s->grid, &s->basis, &s->local_ext, false);
  struct gkyl_array *cell_avg = gkyl_array_new(GKYL_INT, 1, s->local_ext.volume);
  gkyl_dg_calc_pkpm_vars_u(s->vars, s->moms, s->euler, cell_avg, s->u);
  gkyl_array_release(cell_avg);
}

static void
teardown(struct test_setup *s)
{
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_dg_calc_pkpm_vars_release(s->vars);
  gkyl_wave_geom_release(s->geom);
  gkyl_wv_eqn_release(s->ten);
  gkyl_array_release(s->u);
  gkyl_array_release(s->bvar);
  gkyl_array_release(s->euler);
  gkyl_array_release(s->moms);
  gkyl_array_release(s->cmoms);
  gkyl_array_release(s->nmoms);
  gkyl_array_release(s->fluid);
  gkyl_array_release(s->xi);
}

static struct gkyl_dg_eqn *
make_eqn(struct test_setup *s, const struct gkyl_array *flux)
{
  struct gkyl_dg_pkpm_fluid_inp inp = {
    .conf_basis = &s->basis,
    .conf_range = &s->local,
    .auxfields = {.vlasov_pkpm_moms = s->cmoms, .bvar = s->bvar, .u = s->u, .fluid_flux_surf = flux},
    .use_gpu = false,
  };
  return gkyl_dg_pkpm_fluid_inew(&inp);
}

// Pressure tensor P = p_perp I + (p_par - p_perp) b b, heat flux q = b (q_par + 2 q_perp)/2,
// density and energy at reference point z of cell cidx from the p=1 factors.
static void
pressure_at(
  const struct test_setup *s, long cidx, const double *z, double P[3][3], double q[3], double *rho,
  double *E
)
{
  int nc = s->basis.num_basis;
  const double *m = gkyl_array_cfetch(s->moms, cidx), *b = gkyl_array_cfetch(s->bvar, cidx);
  const double *fl = gkyl_array_cfetch(s->fluid, cidx);
  double bq[3];
  for (int i = 0; i < 3; ++i) {
    bq[i] = s->basis.eval_expand(z, &b[i * nc]);
  }
  *rho = s->basis.eval_expand(z, &fl[0]);
  *E = s->basis.eval_expand(z, &fl[4 * nc]);
  double ppar = s->basis.eval_expand(z, &m[nc]), pperp = s->basis.eval_expand(z, &m[2 * nc]);
  double h = s->basis.eval_expand(z, &m[4 * nc]) + 2.0 * s->basis.eval_expand(z, &m[5 * nc]);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      P[i][j] = (ppar - pperp) * bq[i] * bq[j] + (i == j ? pperp : 0.0);
    }
    q[i] = 0.5 * bq[i] * h;
  }
}

static void
test_vol(int cdim)
{
  struct field_ctx ctx = {.cdim = cdim, .uniform_rho = false, .b_along_x = false};
  struct test_setup s;
  setup(&s, &ctx);
  int nc = s.basis.num_basis;
  struct gkyl_dg_eqn *eqn = make_eqn(&s, NULL);

  double max_err = 0.0, max_ref = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *M = gkyl_array_cfetch(s.fluid, cidx), *u = gkyl_array_cfetch(s.u, cidx);
    double xc[3], out[64] = {0.0}, ref[64] = {0.0};
    gkyl_rect_grid_cell_center(&s.grid, iter.idx, xc);
    eqn->vol_term(eqn, xc, s.grid.dx, iter.idx, M, out);

    int nnodes = 1;
    for (int d = 0; d < cdim; ++d) {
      nnodes *= 4;
    }
    for (int n = 0; n < nnodes; ++n) {
      double z[3], w = 1.0;
      int m = n;
      for (int d = cdim - 1; d >= 0; --d) {
        z[d] = gl4_x[m % 4];
        w *= gl4_w[m % 4];
        m /= 4;
      }
      double P[3][3], q[3], rho, E;
      pressure_at(&s, cidx, z, P, q, &rho, &E);
      double Mq[3], uq[3];
      for (int i = 0; i < 3; ++i) {
        Mq[i] = s.basis.eval_expand(z, &M[(1 + i) * nc]);
        uq[i] = s.basis.eval_expand(z, &u[i * nc]);
      }
      for (int k = 0; k < nc; ++k) {
        double unit[8] = {0.0};
        unit[k] = 1.0;
        for (int d = 0; d < cdim; ++d) {
          double dphi = 2.0 / s.grid.dx[d] * s.basis.eval_grad_expand(d, z, unit);
          double Pu = 0.0;
          for (int i = 0; i < 3; ++i) {
            ref[(1 + i) * nc + k] += w * dphi * (Mq[i] * uq[d] + P[i][d]);
            Pu += P[i][d] * uq[i];
          }
          ref[0 * nc + k] += w * dphi * Mq[d];
          ref[4 * nc + k] += w * dphi * (uq[d] * E + Pu + q[d]);
        }
      }
    }
    for (int k = 0; k < 5 * nc; ++k) {
      double err = fabs(out[k] - ref[k]);
      max_err = err > max_err ? err : max_err;
      max_ref = fabs(ref[k]) > max_ref ? fabs(ref[k]) : max_ref;
    }
  }
  TEST_CHECK(max_ref > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_ref);
  TEST_MSG("cdim=%d: max volume error %g (max reference %g)", cdim, max_err, max_ref);

  gkyl_dg_eqn_release(eqn);
  teardown(&s);
}

// Independent assembly of the lower-face flux of cell idx (idx_l the left neighbor):
// HLL on [rho, M, E] with the speed bounds u_n -/+ sqrt((3 + 2 b_n^2) P_nn/rho).
static void
reference_face(
  struct test_setup *s, int dir, const int *idx_l, const int *idx_r, double *flux_ref,
  double *lambda_max
)
{
  int cdim = s->cdim, nc = s->basis.num_basis;
  long cl = gkyl_range_idx(&s->local_ext, idx_l), cr = gkyl_range_idx(&s->local_ext, idx_r);
  const double *M_l = gkyl_array_cfetch(s->fluid, cl), *M_r = gkyl_array_cfetch(s->fluid, cr);
  const double *u_l = gkyl_array_cfetch(s->u, cl), *u_r = gkyl_array_cfetch(s->u, cr);
  const struct gkyl_wave_cell_geom *g = gkyl_wave_geom_get(s->geom, idx_r);
  struct gkyl_basis sbasis;
  int ns = 1;
  if (cdim > 1) {
    gkyl_cart_modal_serendip(&sbasis, cdim - 1, 1);
    ns = sbasis.num_basis;
  }
  for (int k = 0; k < 5 * ns; ++k) {
    flux_ref[k] = 0.0;
  }
  *lambda_max = 0.0;
  int nnodes = 1;
  for (int d = 0; d < cdim - 1; ++d) {
    nnodes *= 3;
  }
  for (int n = 0; n < nnodes; ++n) {
    double zs[3], zl[3], zr[3], w = 1.0;
    int m = n;
    for (int d = cdim - 2; d >= 0; --d) {
      zs[d] = gl3_x[m % 3];
      w *= gl3_w[m % 3];
      m /= 3;
    }
    int t = 0;
    for (int d = 0; d < cdim; ++d) {
      if (d == dir) {
        zl[d] = 1.0;
        zr[d] = -1.0;
      } else {
        zl[d] = zr[d] = zs[t++];
      }
    }
    double Pl[3][3], Pr[3][3], ql[3], qr[3], rho_l, rho_r, E_l, E_r, Ml[3], Mr[3], ul[3], ur[3];
    pressure_at(s, cl, zl, Pl, ql, &rho_l, &E_l);
    pressure_at(s, cr, zr, Pr, qr, &rho_r, &E_r);
    const double *b_l = gkyl_array_cfetch(s->bvar, cl), *b_r = gkyl_array_cfetch(s->bvar, cr);
    double bn_l = s->basis.eval_expand(zl, &b_l[dir * nc]), bn_r = s->basis.eval_expand(zr, &b_r[dir * nc]);
    for (int i = 0; i < 3; ++i) {
      Ml[i] = s->basis.eval_expand(zl, &M_l[(1 + i) * nc]);
      Mr[i] = s->basis.eval_expand(zr, &M_r[(1 + i) * nc]);
      ul[i] = s->basis.eval_expand(zl, &u_l[i * nc]);
      ur[i] = s->basis.eval_expand(zr, &u_r[i * nc]);
    }
    bool ok_l = rho_l > 0.0 && Pl[dir][dir] > 0.0, ok_r = rho_r > 0.0 && Pr[dir][dir] > 0.0;
    double c_l = ok_l ? sqrt((3.0 + 2.0 * bn_l * bn_l) * Pl[dir][dir] / rho_l) : 0.0;
    double c_r = ok_r ? sqrt((3.0 + 2.0 * bn_r * bn_r) * Pr[dir][dir] / rho_r) : 0.0;
    *lambda_max = fmax(*lambda_max, fmax(fabs(ul[dir]) + c_l, fabs(ur[dir]) + c_r));
    double Fl[5], Fr[5], F[5], Ql[5], Qr[5], Pul = 0.0, Pur = 0.0;
    for (int i = 0; i < 3; ++i) {
      Fl[1 + i] = Ml[i] * ul[dir] + Pl[i][dir];
      Fr[1 + i] = Mr[i] * ur[dir] + Pr[i][dir];
      Ql[1 + i] = Ml[i];
      Qr[1 + i] = Mr[i];
      Pul += Pl[i][dir] * ul[i];
      Pur += Pr[i][dir] * ur[i];
    }
    Fl[0] = Ml[dir];
    Fr[0] = Mr[dir];
    Ql[0] = rho_l;
    Qr[0] = rho_r;
    Fl[4] = ul[dir] * E_l + Pul + ql[dir];
    Fr[4] = ur[dir] * E_r + Pur + qr[dir];
    Ql[4] = E_l;
    Qr[4] = E_r;
    {
      double S_L = fmin(ul[dir] - c_l, ur[dir] - c_r), S_R = fmax(ul[dir] + c_l, ur[dir] + c_r);
      for (int i = 0; i < 5; ++i) {
        F[i] = S_L >= 0.0 ? Fl[i] : (S_R <= 0.0 ? Fr[i] : (S_R * Fl[i] - S_L * Fr[i] + S_L * S_R * (Qr[i] - Ql[i])) / (S_R - S_L));
      }
    }
    for (int k = 0; k < ns; ++k) {
      double chi = 1.0;
      if (cdim > 1) {
        double unit[4] = {0.0};
        unit[k] = 1.0;
        chi = sbasis.eval_expand(zs, unit);
      }
      for (int i = 0; i < 5; ++i) {
        flux_ref[i * ns + k] += w * chi * F[i];
      }
    }
  }
}

static void
test_flux_ctx(int cdim, bool inadmissible)
{
  struct field_ctx ctx = {.cdim = cdim, .uniform_rho = false, .b_along_x = false, .inadmissible = inadmissible};
  struct test_setup s;
  setup(&s, &ctx);
  int nc = s.basis.num_basis;

  struct gkyl_dg_pkpm_fluid_flux_surf_inp inp = {
    .conf_grid = &s.grid,
    .conf_basis = &s.basis,
    .wv_eqn = s.ten,
    .geom = s.geom,
    .use_gpu = false,
  };
  struct gkyl_dg_pkpm_fluid_flux_surf *up = gkyl_dg_pkpm_fluid_flux_surf_inew(&inp);
  int ns = gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(up);
  TEST_CHECK(ns == nc / 2);
  struct gkyl_array *flux = mkarr(5 * cdim * ns, s.local_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, s.local_ext.volume);
  gkyl_array_clear(flux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_pkpm_fluid_flux_surf_advance(up, &s.local, &s.local_ext, s.cmoms, s.bvar, s.u, s.fluid, cflrate, flux);

  double flux_ref[20];
  double max_err = 0.0, max_flux = 0.0, max_cfl_err = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *fl = gkyl_array_cfetch(flux, cidx);
    double cfl_ref = 0.0;
    for (int dir = 0; dir < cdim; ++dir) {
      int idx_l[3], idx_r[3];
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      idx_l[dir] -= 1;
      double lam;
      reference_face(&s, dir, idx_l, iter.idx, flux_ref, &lam);
      cfl_ref += 1.5 * 2.0 / s.grid.dx[dir] * lam;
      for (int k = 0; k < 5 * ns; ++k) {
        double err = fabs(fl[dir * 5 * ns + k] - flux_ref[k]);
        max_err = err > max_err ? err : max_err;
        max_flux = fabs(flux_ref[k]) > max_flux ? fabs(flux_ref[k]) : max_flux;
      }
      if (iter.idx[dir] == s.local.upper[dir]) {
        gkyl_copy_int_arr(cdim, iter.idx, idx_r);
        idx_r[dir] += 1;
        const double *fg = gkyl_array_cfetch(flux, gkyl_range_idx(&s.local_ext, idx_r));
        reference_face(&s, dir, iter.idx, idx_r, flux_ref, &lam);
        cfl_ref += 1.5 * 2.0 / s.grid.dx[dir] * lam;
        for (int k = 0; k < 5 * ns; ++k) {
          double err = fabs(fg[dir * 5 * ns + k] - flux_ref[k]);
          max_err = err > max_err ? err : max_err;
        }
      }
    }
    const double *cfl = gkyl_array_cfetch(cflrate, cidx);
    double cerr = fabs(cfl[0] - cfl_ref);
    max_cfl_err = cerr > max_cfl_err ? cerr : max_cfl_err;
  }
  TEST_CHECK(max_flux > 0.0);
  TEST_CHECK(max_err < 1e-11 * max_flux);
  TEST_MSG("cdim=%d: max flux error %g (max flux %g)", cdim, max_err, max_flux);
  TEST_CHECK(max_cfl_err < 1e-10);
  TEST_MSG("cdim=%d: max cflrate error %g", cdim, max_cfl_err);

  gkyl_array_release(cflrate);
  gkyl_array_release(flux);
  gkyl_dg_pkpm_fluid_flux_surf_release(up);
  teardown(&s);
}

// Absorbing-boundary fluxes (dg_pkpm_fluid_flux_surf_absorb): for a uniform state with
// Maxwellian [F, G] the moments of the kinetic upwind flux through the wall have the
// closed form I_k = int_{a > 0} w^k a M(w) dw (a = xi b_n w + u_n), integrated here by
// Simpson's rule; u_n is chosen so that the switch a = 0 lies on a velocity cell edge,
// where the 3-point Gauss-Legendre integration of the projected Maxwellian is exact.
struct absorb_ctx {
  int cdim;
  double rho, xi, u[3], b[3];
};

static void
eval_absorb_fluid(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  double uu = c->u[0] * c->u[0] + c->u[1] * c->u[1] + c->u[2] * c->u[2];
  fout[0] = c->rho;
  for (int i = 0; i < 3; ++i) {
    fout[1 + i] = c->rho * c->u[i];
  }
  fout[4] = 0.5 * c->rho * uu + 1.5 * c->rho * c->xi * c->xi;
}

static void
eval_absorb_u(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  for (int i = 0; i < 3; ++i) {
    fout[i] = c->u[i];
  }
}

static void
eval_absorb_b(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  for (int i = 0; i < 3; ++i) {
    fout[i] = c->b[i];
  }
}

static void
eval_absorb_xi(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  fout[0] = c->xi;
  fout[1] = c->xi * c->xi;
  fout[2] = c->rho * c->xi * c->xi;
}

static void
eval_absorb_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct absorb_ctx *c = ctx;
  double w = xn[c->cdim];
  fout[0] = fout[1] = exp(-0.5 * w * w) / sqrt(2.0 * M_PI);
}

// I_k = int w^k phi(w) M(w) dw over w in [wa, wb] by Simpson's rule, phi = xi b_n w + u_n.
static double
absorb_moment(int k, double wa, double wb, double xibn, double un)
{
  int n = 20000;
  double h = (wb - wa) / n, sum = 0.0;
  for (int i = 0; i <= n; ++i) {
    double w = wa + i * h, wt = (i == 0 || i == n) ? 1.0 : (i % 2 ? 4.0 : 2.0);
    sum += wt * pow(w, k) * (xibn * w + un) * exp(-0.5 * w * w) / sqrt(2.0 * M_PI);
  }
  return sum * h / 3.0;
}

static void
test_absorb(int cdim)
{
  double wmax = 8.0;
  int nv = 96; // dv = 1/6
  struct absorb_ctx ctx = {.cdim = cdim, .rho = 2.0, .xi = 1.5};
  double theta = cdim > 1 ? 0.4 : 0.0;
  ctx.b[0] = cos(theta);
  ctx.b[1] = sin(theta);
  ctx.b[2] = 0.0;
  // a = 0 at w0 = 1/3 (a velocity cell edge): u_n = -xi b_n w0
  ctx.u[0] = -ctx.xi * ctx.b[0] / 3.0;
  ctx.u[1] = 0.1;
  ctx.u[2] = -0.2;

  double lower[4] = {-1.0, -1.0, -1.0, -wmax}, upper[4] = {1.0, 1.0, 1.0, wmax};
  int cells[4] = {4, 3, 2, nv}, ghost[4] = {1, 1, 1, 0};
  lower[cdim] = -wmax;
  upper[cdim] = wmax;
  cells[cdim] = nv;
  ghost[cdim] = 0;
  struct gkyl_rect_grid cgrid, pgrid;
  gkyl_rect_grid_init(&cgrid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&pgrid, cdim + 1, lower, upper, cells);
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext;
  gkyl_create_grid_ranges(&cgrid, ghost, &clocal_ext, &clocal);
  gkyl_create_grid_ranges(&pgrid, ghost, &plocal_ext, &plocal);
  int nc = cbasis.num_basis, np = pbasis.num_basis;

  struct gkyl_array *fluid = mkarr(5 * nc, clocal_ext.volume), *u = mkarr(3 * nc, clocal_ext.volume);
  struct gkyl_array *bvar = mkarr(3 * nc, clocal_ext.volume), *xi = mkarr(3 * nc, clocal_ext.volume);
  struct gkyl_array *f = mkarr(2 * np, plocal_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 5, eval_absorb_fluid, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, fluid);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 3, eval_absorb_u, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, u);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 3, eval_absorb_b, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, bvar);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&cgrid, &cbasis, 8, 3, eval_absorb_xi, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &clocal_ext, xi);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&pgrid, &pbasis, 8, 2, eval_absorb_f, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &plocal_ext, f);
  gkyl_proj_on_basis_release(proj);

  struct gkyl_wv_eqn *ten = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&cgrid, &clocal_ext, 0, 0, false);
  struct gkyl_dg_pkpm_fluid_flux_surf_inp inp = {
    .conf_grid = &cgrid, .conf_basis = &cbasis, .wv_eqn = ten, .geom = geom, .phase_grid = &pgrid, .use_gpu = false};
  struct gkyl_dg_pkpm_fluid_flux_surf *up = gkyl_dg_pkpm_fluid_flux_surf_inew(&inp);
  int ns = gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(up);
  struct gkyl_array *flux = mkarr(5 * cdim * ns, clocal_ext.volume);
  gkyl_array_clear(flux, 0.0);

  int dir = 0;
  double xibn = ctx.xi * ctx.b[dir], un = ctx.u[dir], w0 = -un / xibn;
  double bu = 0.0, uu = 0.0;
  for (int i = 0; i < 3; ++i) {
    bu += ctx.b[i] * ctx.u[i];
    uu += ctx.u[i] * ctx.u[i];
  }
  // modal coefficient of a constant on the surface basis: sqrt(2)^(cdim-1)
  double c0 = pow(sqrt(2.0), cdim - 1);
  for (int e = 0; e < 2; ++e) {
    enum gkyl_edge_loc edge = e == 0 ? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
    struct gkyl_range skin, ghost_r;
    gkyl_skin_ghost_ranges(&skin, &ghost_r, dir, edge, &clocal_ext, ghost);
    gkyl_dg_pkpm_fluid_flux_surf_absorb(
      up, dir, edge, &skin, &clocal_ext, &plocal, &plocal_ext, bvar, u, xi, fluid, f, flux
    );
    // outflow: a < 0 (w < w0) through the lower wall, a > 0 (w > w0) through the upper
    double I[3], J0;
    for (int k = 0; k < 3; ++k) {
      I[k] = e == 0 ? absorb_moment(k, -wmax, w0, xibn, un) : absorb_moment(k, w0, wmax, xibn, un);
    }
    J0 = I[0];
    double ref[5];
    ref[0] = ctx.rho * I[0];
    for (int i = 0; i < 3; ++i) {
      double delta = i == dir ? 1.0 : 0.0;
      ref[1 + i] = ctx.rho * (ctx.u[i] * I[0] + ctx.xi * ctx.b[i] * I[1])
        + ctx.rho * ctx.xi * ctx.xi * (delta - ctx.b[i] * ctx.b[dir]) * J0;
    }
    ref[4] = 0.5 * ctx.rho * (uu * I[0] + 2.0 * ctx.xi * bu * I[1] + ctx.xi * ctx.xi * I[2])
      + ctx.rho * ctx.xi * ctx.xi * J0;
    // the flux of the lower wall is stored in the skin cell, that of the upper wall in
    // the ghost cell beyond it
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &skin);
    while (gkyl_range_iter_next(&iter)) {
      int idx[3];
      gkyl_copy_int_arr(cdim, iter.idx, idx);
      if (e == 1) {
        idx[dir] += 1;
      }
      const double *fl = gkyl_array_cfetch(flux, gkyl_range_idx(&clocal_ext, idx));
      for (int c = 0; c < 5; ++c) {
        double got = fl[dir * 5 * ns + c * ns] / c0;
        TEST_CHECK(gkyl_compare(got, ref[c], 1e-4 * fabs(ref[4])));
        TEST_MSG("cdim=%d edge=%d comp=%d: got %g expected %g", cdim, e, c, got, ref[c]);
        for (int k = 1; k < ns; ++k) {
          TEST_CHECK(fabs(fl[dir * 5 * ns + c * ns + k]) < 1e-12 * fabs(ref[4]));
        }
      }
    }
    TEST_CHECK(e == 0 ? ref[0] < 0.0 : ref[0] > 0.0);
  }

  gkyl_array_release(flux);
  gkyl_dg_pkpm_fluid_flux_surf_release(up);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(ten);
  gkyl_array_release(f);
  gkyl_array_release(xi);
  gkyl_array_release(bvar);
  gkyl_array_release(u);
  gkyl_array_release(fluid);
}
void test_absorb_1x() { test_absorb(1); }
void test_absorb_2x() { test_absorb(2); }
void test_vol_1x() { test_vol(1); }
void test_vol_2x() { test_vol(2); }
void test_vol_3x() { test_vol(3); }
void test_flux_1x() { test_flux_ctx(1, false); }
void test_flux_2x() { test_flux_ctx(2, false); }
void test_flux_3x() { test_flux_ctx(3, false); }
// Part of the domain has P_nn <= 0: the HLL fallback must be taken there.
void test_flux_hll_1x() { test_flux_ctx(1, true); }
void test_flux_hll_2x() { test_flux_ctx(2, true); }

TEST_LIST = {
  {"test_vol_1x", test_vol_1x},
  {"test_vol_2x", test_vol_2x},
  {"test_vol_3x", test_vol_3x},
  {"test_flux_1x", test_flux_1x},
  {"test_flux_2x", test_flux_2x},
  {"test_flux_3x", test_flux_3x},
  {"test_flux_hll_1x", test_flux_hll_1x},
  {"test_flux_hll_2x", test_flux_hll_2x},
  {"test_absorb_1x", test_absorb_1x},
  {"test_absorb_2x", test_absorb_2x},
  {NULL, NULL},
};
