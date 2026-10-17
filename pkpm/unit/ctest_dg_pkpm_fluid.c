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
#include <gkyl_dg_euler_pkpm.h>
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

// [rho, p_par, p_perp, M1]
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
  if (c->inadmissible && x > 1.0) {
    fout[1] = -0.2; // negative parallel pressure: P_nn <= 0 where b is along the normal
    fout[2] = 0.05;
  }
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
  struct gkyl_array *moms, *euler, *bvar, *u;
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

  s->moms = mkarr(4 * nc, s->local_ext.volume);
  s->euler = mkarr(3 * nc, s->local_ext.volume);

  s->bvar = mkarr(3 * nc, s->local_ext.volume);
  s->u = mkarr(3 * nc, s->local_ext.volume);
  struct gkyl_array *em = mkarr(8 * nc, s->local_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 4, eval_moms, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->moms);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 3, eval_euler, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, s->euler);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->grid, &s->basis, 8, 8, eval_em, ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local_ext, em);
  gkyl_proj_on_basis_release(proj);

  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->grid, &s->basis, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->local_ext, em, s->bvar);
  gkyl_array_release(em);

  s->ten = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  s->geom = gkyl_wave_geom_new(&s->grid, &s->local_ext, 0, 0, false);
  s->vars = gkyl_dg_calc_pkpm_vars_new(&s->grid, &s->basis, &s->local_ext, s->ten, s->geom, 0.0, false);
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
}

static struct gkyl_dg_eqn *
make_eqn(struct test_setup *s, const struct gkyl_array *flux)
{
  struct gkyl_dg_pkpm_fluid_inp inp = {
    .conf_basis = &s->basis,
    .conf_range = &s->local,
    .auxfields = {.moms = s->moms, .bvar = s->bvar, .u = s->u, .fluid_flux_surf = flux},
    .use_gpu = false,
  };
  return gkyl_dg_pkpm_fluid_inew(&inp);
}

static void
pressure_at(const struct test_setup *s, long cidx, const double *z, double P[3][3], double *rho)
{
  int nc = s->basis.num_basis;
  const double *m = gkyl_array_cfetch(s->moms, cidx), *b = gkyl_array_cfetch(s->bvar, cidx);
  double bq[3];
  for (int i = 0; i < 3; ++i) {
    bq[i] = s->basis.eval_expand(z, &b[i * nc]);
  }
  double ppar = s->basis.eval_expand(z, &m[nc]), pperp = s->basis.eval_expand(z, &m[2 * nc]);
  *rho = s->basis.eval_expand(z, &m[0]);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      P[i][j] = (ppar - pperp) * bq[i] * bq[j] + (i == j ? pperp : 0.0);
    }
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
    const double *M = gkyl_array_cfetch(s.euler, cidx), *u = gkyl_array_cfetch(s.u, cidx);
    double xc[3], out[24] = {0.0}, ref[24] = {0.0};
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
      double P[3][3], rho;
      pressure_at(&s, cidx, z, P, &rho);
      double Mq[3], uq[3];
      for (int i = 0; i < 3; ++i) {
        Mq[i] = s.basis.eval_expand(z, &M[i * nc]);
        uq[i] = s.basis.eval_expand(z, &u[i * nc]);
      }
      for (int k = 0; k < nc; ++k) {
        double unit[8] = {0.0};
        unit[k] = 1.0;
        for (int d = 0; d < cdim; ++d) {
          double dphi = 2.0 / s.grid.dx[d] * s.basis.eval_grad_expand(d, z, unit);
          for (int i = 0; i < 3; ++i) {
            ref[i * nc + k] += w * dphi * (Mq[i] * uq[d] + P[i][d]);
          }
        }
      }
    }
    for (int k = 0; k < 3 * nc; ++k) {
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

// Independent assembly of the lower-face flux of cell idx (idx_l the left neighbor),
// sharing only the ten-moment Riemann solver with the kernel.
static void
reference_face(
  struct test_setup *s, int dir, const int *idx_l, const int *idx_r, double *flux_ref,
  double *mass_ref, double *lambda_max
)
{
  int cdim = s->cdim, nc = s->basis.num_basis;
  long cl = gkyl_range_idx(&s->local_ext, idx_l), cr = gkyl_range_idx(&s->local_ext, idx_r);
  const double *M_l = gkyl_array_cfetch(s->euler, cl), *M_r = gkyl_array_cfetch(s->euler, cr);
  const double *u_l = gkyl_array_cfetch(s->u, cl), *u_r = gkyl_array_cfetch(s->u, cr);
  const struct gkyl_wave_cell_geom *g = gkyl_wave_geom_get(s->geom, idx_r);
  struct gkyl_basis sbasis;
  int ns = 1;
  if (cdim > 1) {
    gkyl_cart_modal_serendip(&sbasis, cdim - 1, 1);
    ns = sbasis.num_basis;
  }
  for (int k = 0; k < 3 * ns; ++k) {
    flux_ref[k] = 0.0;
  }
  for (int k = 0; k < ns; ++k) {
    mass_ref[k] = 0.0;
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
    double Pl[3][3], Pr[3][3], rho_l, rho_r, Ml[3], Mr[3], ul[3], ur[3];
    pressure_at(s, cl, zl, Pl, &rho_l);
    pressure_at(s, cr, zr, Pr, &rho_r);
    for (int i = 0; i < 3; ++i) {
      Ml[i] = s->basis.eval_expand(zl, &M_l[i * nc]);
      Mr[i] = s->basis.eval_expand(zr, &M_r[i * nc]);
      ul[i] = s->basis.eval_expand(zl, &u_l[i * nc]);
      ur[i] = s->basis.eval_expand(zr, &u_r[i * nc]);
    }
    bool ok_l = rho_l > 0.0 && Pl[dir][dir] > 0.0, ok_r = rho_r > 0.0 && Pr[dir][dir] > 0.0;
    double c_l = ok_l ? sqrt(3.0 * Pl[dir][dir] / rho_l) : 0.0;
    double c_r = ok_r ? sqrt(3.0 * Pr[dir][dir] / rho_r) : 0.0;
    *lambda_max = fmax(*lambda_max, fmax(fabs(ul[dir]) + c_l, fabs(ur[dir]) + c_r));
    double Fl[3], Fr[3], F[3], F0;
    for (int i = 0; i < 3; ++i) {
      Fl[i] = Ml[i] * ul[dir] + Pl[i][dir];
      Fr[i] = Mr[i] * ur[dir] + Pr[i][dir];
    }
    if (!(ok_l && ok_r)) {
      double S_L = fmin(ul[dir] - c_l, ur[dir] - c_r), S_R = fmax(ul[dir] + c_l, ur[dir] + c_r);
      for (int i = 0; i < 3; ++i) {
        F[i] = S_L >= 0.0 ? Fl[i] : (S_R <= 0.0 ? Fr[i] : (S_R * Fl[i] - S_L * Fr[i] + S_L * S_R * (Mr[i] - Ml[i])) / (S_R - S_L));
      }
      F0 = S_L >= 0.0 ? Ml[dir] : (S_R <= 0.0 ? Mr[dir] : (S_R * Ml[dir] - S_L * Mr[dir] + S_L * S_R * (rho_r - rho_l)) / (S_R - S_L));
    } else {
      double ql[10], qr[10], qll[10], qrl[10], delta[10], waves[50], speeds[5], amdq[10], apdq[10], am[10], ap[10];
      ql[0] = rho_l;
      qr[0] = rho_r;
      for (int i = 0; i < 3; ++i) {
        ql[1 + i] = Ml[i];
        qr[1 + i] = Mr[i];
      }
      int ij[6][2] = {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
      for (int c = 0; c < 6; ++c) {
        int i = ij[c][0], j = ij[c][1];
        ql[4 + c] = Pl[i][j] + Ml[i] * Ml[j] / rho_l;
        qr[4 + c] = Pr[i][j] + Mr[i] * Mr[j] / rho_r;
      }
      rot_to_local(s->ten, g->tau1[dir], g->tau2[dir], g->norm[dir], ql, qll);
      rot_to_local(s->ten, g->tau1[dir], g->tau2[dir], g->norm[dir], qr, qrl);
      for (int i = 0; i < 10; ++i) {
        delta[i] = qrl[i] - qll[i];
      }
      wave_roe(s->ten, delta, qll, qrl, waves, speeds);
      for (int i = 0; i < 5; ++i) {
        speeds[i] *= g->lenr[dir];
      }
      qfluct_roe(s->ten, qll, qrl, waves, speeds, am, ap);
      rot_to_global(s->ten, g->tau1[dir], g->tau2[dir], g->norm[dir], am, amdq);
      rot_to_global(s->ten, g->tau1[dir], g->tau2[dir], g->norm[dir], ap, apdq);
      for (int i = 0; i < 3; ++i) {
        F[i] = 0.5 * (Fl[i] + Fr[i]) - 0.5 * (apdq[1 + i] - amdq[1 + i]);
      }
      F0 = 0.5 * (Ml[dir] + Mr[dir]) - 0.5 * (apdq[0] - amdq[0]);
    }
    for (int k = 0; k < ns; ++k) {
      double chi = 1.0;
      if (cdim > 1) {
        double unit[4] = {0.0};
        unit[k] = 1.0;
        chi = sbasis.eval_expand(zs, unit);
      }
      for (int i = 0; i < 3; ++i) {
        flux_ref[i * ns + k] += w * chi * F[i];
      }
      mass_ref[k] += w * chi * F0;
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
  struct gkyl_array *flux = mkarr(3 * cdim * ns, s.local_ext.volume);
  struct gkyl_array *mflux = mkarr(cdim * ns, s.local_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, s.local_ext.volume);
  gkyl_array_clear(flux, 0.0);
  gkyl_array_clear(mflux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_pkpm_fluid_flux_surf_advance(up, &s.local, &s.local_ext, s.moms, s.bvar, s.u, s.euler, cflrate, flux, mflux);

  double flux_ref[12], mass_ref[4];
  double max_err = 0.0, max_flux = 0.0, max_cfl_err = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *fl = gkyl_array_cfetch(flux, cidx);
    const double *ml = gkyl_array_cfetch(mflux, cidx);
    double cfl_ref = 0.0;
    for (int dir = 0; dir < cdim; ++dir) {
      int idx_l[3], idx_r[3];
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      idx_l[dir] -= 1;
      double lam;
      reference_face(&s, dir, idx_l, iter.idx, flux_ref, mass_ref, &lam);
      cfl_ref += 1.5 * 2.0 / s.grid.dx[dir] * lam;
      for (int k = 0; k < 3 * ns; ++k) {
        double err = fabs(fl[dir * 3 * ns + k] - flux_ref[k]);
        max_err = err > max_err ? err : max_err;
        max_flux = fabs(flux_ref[k]) > max_flux ? fabs(flux_ref[k]) : max_flux;
      }
      for (int k = 0; k < ns; ++k) {
        double err = fabs(ml[dir * ns + k] - mass_ref[k]);
        max_err = err > max_err ? err : max_err;
      }
      if (iter.idx[dir] == s.local.upper[dir]) {
        gkyl_copy_int_arr(cdim, iter.idx, idx_r);
        idx_r[dir] += 1;
        const double *fg = gkyl_array_cfetch(flux, gkyl_range_idx(&s.local_ext, idx_r));
        const double *mg = gkyl_array_cfetch(mflux, gkyl_range_idx(&s.local_ext, idx_r));
        reference_face(&s, dir, iter.idx, idx_r, flux_ref, mass_ref, &lam);
        cfl_ref += 1.5 * 2.0 / s.grid.dx[dir] * lam;
        for (int k = 0; k < 3 * ns; ++k) {
          double err = fabs(fg[dir * 3 * ns + k] - flux_ref[k]);
          max_err = err > max_err ? err : max_err;
        }
        for (int k = 0; k < ns; ++k) {
          double err = fabs(mg[dir * ns + k] - mass_ref[k]);
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
  gkyl_array_release(mflux);
  gkyl_dg_pkpm_fluid_flux_surf_release(up);
  teardown(&s);
}

// 1x, b along x, uniform density and zero flow (pressure-only dynamics): the legacy
// pressure projection is exact, the advective fluxes vanish and the kinetic mass flux
// correction is multiplied by u = 0, so the new and legacy right-hand sides coincide.
static void
test_legacy_1x()
{
  struct field_ctx ctx = {.cdim = 1, .uniform_rho = true, .b_along_x = true, .zero_u = true};
  struct test_setup s;
  setup(&s, &ctx);
  int cdim = 1, nc = s.basis.num_basis;

  // New right-hand side.
  struct gkyl_dg_pkpm_fluid_flux_surf_inp inp = {
    .conf_grid = &s.grid, .conf_basis = &s.basis, .wv_eqn = s.ten, .geom = s.geom};
  struct gkyl_dg_pkpm_fluid_flux_surf *up = gkyl_dg_pkpm_fluid_flux_surf_inew(&inp);
  int ns = gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(up);
  struct gkyl_array *flux = mkarr(3 * cdim * ns, s.local_ext.volume);
  struct gkyl_array *mflux = mkarr(cdim * ns, s.local_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, s.local_ext.volume);
  gkyl_array_clear(flux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_pkpm_fluid_flux_surf_advance(up, &s.local, &s.local_ext, s.moms, s.bvar, s.u, s.euler, cflrate, flux, mflux);
  gkyl_array_release(mflux);
  struct gkyl_dg_eqn *eqn = make_eqn(&s, flux);
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *slvr = gkyl_hyper_dg_new(&s.grid, &s.basis, eqn, cdim, up_dirs, zero_flux, 1, false);
  struct gkyl_array *rhs = mkarr(3 * nc, s.local_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_hyper_dg_advance(slvr, &s.local, s.euler, cflrate, rhs);

  // Legacy right-hand side from the legacy variable pipeline.
  struct gkyl_array *bvar_leg = mkarr(9 * nc, s.local_ext.volume);
  gkyl_array_clear(bvar_leg, 0.0);
  gkyl_array_shiftc(bvar_leg, sqrt(2.0), 0 * nc); // bx = 1
  gkyl_array_shiftc(bvar_leg, sqrt(2.0), 3 * nc); // bxbx = 1
  struct gkyl_array *p_ij = mkarr(6 * nc, s.local_ext.volume);
  struct gkyl_array *div_ppar = mkarr(nc, s.local_ext.volume);
  struct gkyl_array *prim = mkarr(9 * nc, s.local_ext.volume);
  int nsurf = nc / 2;
  struct gkyl_array *prim_surf = mkarr(2 * cdim * 4 * nsurf, s.local_ext.volume);
  struct gkyl_array *lax = mkarr(cdim * nsurf, s.local_ext.volume);
  struct gkyl_array *pen = mkarr(3 * cdim * nsurf, s.local_ext.volume);
  struct gkyl_array *cell_avg = gkyl_array_new(GKYL_INT, 1, s.local_ext.volume);
  gkyl_array_clear(div_ppar, 0.0);
  gkyl_array_clear(p_ij, 0.0);
  gkyl_array_clear(prim, 0.0);
  gkyl_array_clear(prim_surf, 0.0);
  gkyl_array_clear(lax, 0.0);
  gkyl_array_clear(pen, 0.0);
  gkyl_dg_calc_pkpm_vars_pressure(s.vars, &s.local_ext, bvar_leg, s.moms, p_ij);
  gkyl_dg_calc_pkpm_vars_advance(s.vars, s.moms, s.euler, p_ij, div_ppar, cell_avg, prim, prim_surf);
  gkyl_dg_calc_pkpm_vars_penalization(s.vars, &s.local, &s.local_ext, s.moms, p_ij, prim, s.euler, lax, pen);

  struct gkyl_dg_eqn *eqn_leg = gkyl_dg_euler_pkpm_new(&s.basis, &s.local, false);
  struct gkyl_dg_euler_pkpm_auxfields aux = {
    .vlasov_pkpm_moms = s.moms,
    .pkpm_prim = prim,
    .pkpm_prim_surf = prim_surf,
    .pkpm_p_ij = p_ij,
    .pkpm_lax = lax,
    .pkpm_penalization = pen,
  };
  gkyl_euler_pkpm_set_auxfields(eqn_leg, aux);
  gkyl_hyper_dg *slvr_leg = gkyl_hyper_dg_new(&s.grid, &s.basis, eqn_leg, cdim, up_dirs, zero_flux, 1, false);
  struct gkyl_array *rhs_leg = mkarr(3 * nc, s.local_ext.volume);
  struct gkyl_array *cfl_leg = mkarr(1, s.local_ext.volume);
  gkyl_array_clear(rhs_leg, 0.0);
  gkyl_array_clear(cfl_leg, 0.0);
  gkyl_hyper_dg_advance(slvr_leg, &s.local, s.euler, cfl_leg, rhs_leg);

  // The new u equals the legacy weak-division u.
  double max_u_err = 0.0, max_err = 0.0, max_rhs = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.local);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&s.local_ext, iter.idx);
    const double *u = gkyl_array_cfetch(s.u, cidx), *pr = gkyl_array_cfetch(prim, cidx);
    const double *r = gkyl_array_cfetch(rhs, cidx), *rl = gkyl_array_cfetch(rhs_leg, cidx);
    for (int k = 0; k < 3 * nc; ++k) {
      max_u_err = fmax(max_u_err, fabs(u[k] - pr[k]));
      max_err = fmax(max_err, fabs(r[k] - rl[k]));
      max_rhs = fmax(max_rhs, fabs(rl[k]));
    }
  }
  TEST_CHECK(max_u_err < 1e-13);
  TEST_CHECK(max_rhs > 0.0);
  TEST_CHECK(max_err < 1e-11 * max_rhs);
  TEST_MSG("legacy comparison: max rhs error %g (max rhs %g, max u error %g)", max_err, max_rhs, max_u_err);

  gkyl_array_release(cfl_leg);
  gkyl_array_release(rhs_leg);
  gkyl_hyper_dg_release(slvr_leg);
  gkyl_dg_eqn_release(eqn_leg);
  gkyl_array_release(cell_avg);
  gkyl_array_release(pen);
  gkyl_array_release(lax);
  gkyl_array_release(prim_surf);
  gkyl_array_release(prim);
  gkyl_array_release(div_ppar);
  gkyl_array_release(p_ij);
  gkyl_array_release(bvar_leg);
  gkyl_array_release(rhs);
  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_array_release(cflrate);
  gkyl_array_release(flux);
  gkyl_dg_pkpm_fluid_flux_surf_release(up);
  teardown(&s);
}

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
  {"test_legacy_1x", test_legacy_1x},
  {NULL, NULL},
};
