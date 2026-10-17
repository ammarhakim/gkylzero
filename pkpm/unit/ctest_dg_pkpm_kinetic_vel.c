// Tests of the velocity-space pass (stage B) of the refactored PKPM kinetic update:
// the v_par flux (dg_pkpm_vel_flux_surf) against an independent nodal evaluation,
// the stage-B volume kernel against quadrature, and the first-moment invariant of the
// complete two-pass right-hand side: M1(rhs_A + rhs_B) = 0 per cell and mode.
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_dg_pkpm_kinetic_vel.h>
#include <gkyl_dg_pkpm_vel_flux_surf.h>
#include <gkyl_dg_updater_moment_pkpm.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_wv_ten_moment.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  for (int i = 0; i < 8; ++i) {
    fout[i] = 0.0;
  }
  fout[3] = 1.0 + 0.4 * sin(x) * cos(y);
  fout[4] = 0.6 * cos(x) * sin(z) + 0.5 * sin(y);
  fout[5] = -0.5 + 0.3 * sin(y) * cos(z) + 0.2 * sin(x);
}

static void
eval_u(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 0.8 * sin(x) + 0.3 * cos(y);
  fout[1] = -0.7 * cos(y) * sin(z) + 0.2 * sin(x);
  fout[2] = 0.5 * sin(z) - 0.4 * cos(x) * cos(y);
}

static void
eval_conf_scalars(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  fout[0] = 0.7 + 0.2 * cos(x) * sin(y) + 0.1 * z; // theta
  fout[1] = 0.3 + 0.1 * sin(x); // nu_sum
  fout[2] = 0.2 + 0.05 * cos(y) + 0.03 * sin(z); // nu_vthsq_sum
  fout[3] = 0.4 * sin(x) - 0.3 * cos(y) * sin(z); // a stand-in a_P for the flux test
}

static void
eval_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double v = xn[cdim];
  double n = 1.0 + 0.3 * sin(x) * cos(y) + 0.1 * cos(z);
  double vt = 1.0 + 0.2 * cos(x);
  double F0 = n / sqrt(2.0 * M_PI * vt * vt) * exp(-v * v / (2.0 * vt * vt));
  fout[0] = F0 * (1.0 + 0.1 * v); // a deliberately nonzero first moment
  fout[1] = 0.5 * vt * vt * F0 * (1.0 - 0.1 * v);
}

static const double gl3_x[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double gl3_w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
static const double gl4_x[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
  0.8611363115940526};
static const double gl4_w[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
  0.3478548451374538};

struct test_setup {
  int cdim, pdim;
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext, vlocal, vlocal_ext;
  struct gkyl_array *bvar, *grad_b, *u, *grad_u, *theta, *nu_sum, *nu_vthsq_sum, *a_P, *f;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
};

// Weak gradient of u on the interior cells (u is projected on the extended range).
static void
compute_grad_u(struct test_setup *s)
{
  int cdim = s->cdim, nc = s->cbasis.num_basis;
  s->grad_u = mkarr(cdim * 3 * nc, s->clocal_ext.volume);
  gkyl_array_clear(s->grad_u, 0.0);
  struct gkyl_wv_eqn *ten = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&s->cgrid, &s->clocal_ext, 0, 0, false);
  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_new(&s->cgrid, &s->cbasis, &s->clocal_ext, ten, geom, 0.0, false);
  gkyl_dg_calc_pkpm_vars_grad_u(vars, &s->clocal, s->u, s->grad_u);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(ten);
}

static void
setup(struct test_setup *s, int cdim)
{
  int pdim = cdim + 1;
  s->cdim = cdim;
  s->pdim = pdim;
  double lower[4] = {-M_PI, -M_PI, -M_PI, -4.0}, upper[4] = {M_PI, M_PI, M_PI, 4.0};
  int cells[4] = {6, 4, 3, 8};
  int ccells[3] = {6, 4, 3};
  int cghost[3] = {1, 1, 1}, pghost[4] = {1, 1, 1, 0}, vghost[1] = {0};
  if (cdim < 3) {
    lower[cdim] = -4.0;
    upper[cdim] = 4.0;
    cells[cdim] = 8;
    pghost[cdim] = 0;
  }
  gkyl_rect_grid_init(&s->cgrid, cdim, lower, upper, ccells);
  gkyl_rect_grid_init(&s->pgrid, pdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&s->pbasis, cdim, 1);
  gkyl_create_grid_ranges(&s->cgrid, cghost, &s->clocal_ext, &s->clocal);
  gkyl_create_grid_ranges(&s->pgrid, pghost, &s->plocal_ext, &s->plocal);
  struct gkyl_rect_grid vgrid;
  gkyl_rect_grid_init(&vgrid, 1, &lower[cdim], &upper[cdim], &cells[cdim]);
  gkyl_create_grid_ranges(&vgrid, vghost, &s->vlocal_ext, &s->vlocal);

  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  struct gkyl_array *em = mkarr(8 * nc, s->clocal_ext.volume);
  struct gkyl_array *scal = mkarr(4 * nc, s->clocal_ext.volume);
  s->bvar = mkarr(3 * nc, s->clocal_ext.volume);
  s->u = mkarr(3 * nc, s->clocal_ext.volume);
  s->theta = mkarr(nc, s->clocal_ext.volume);
  s->nu_sum = mkarr(nc, s->clocal_ext.volume);
  s->nu_vthsq_sum = mkarr(nc, s->clocal_ext.volume);
  s->a_P = mkarr(nc, s->clocal_ext.volume);
  s->f = mkarr(2 * np, s->plocal_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, em);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_u, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->u);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 4, eval_conf_scalars, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, scal);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_set_offset(s->theta, 1.0, scal, 0 * nc);
  gkyl_array_set_offset(s->nu_sum, 1.0, scal, 1 * nc);
  gkyl_array_set_offset(s->nu_vthsq_sum, 1.0, scal, 2 * nc);
  gkyl_array_set_offset(s->a_P, 1.0, scal, 3 * nc);
  proj = gkyl_proj_on_basis_new(&s->pgrid, &s->pbasis, 8, 2, eval_f, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->plocal_ext, s->f);
  gkyl_proj_on_basis_release(proj);

  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->cgrid, &s->cbasis, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->clocal_ext, em, s->bvar);
  s->grad_b = mkarr(cdim * 3 * nc, s->clocal_ext.volume);
  gkyl_array_clear(s->grad_b, 0.0);
  gkyl_dg_calc_pkpm_bvar_grad(s->calc_b, &s->clocal, s->bvar, s->grad_b);
  gkyl_array_release(em);
  gkyl_array_release(scal);
  compute_grad_u(s);
}

static void
teardown(struct test_setup *s)
{
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_array_release(s->grad_b);
  gkyl_array_release(s->grad_u);
  gkyl_array_release(s->f);
  gkyl_array_release(s->a_P);
  gkyl_array_release(s->nu_vthsq_sum);
  gkyl_array_release(s->nu_sum);
  gkyl_array_release(s->theta);
  gkyl_array_release(s->u);
  gkyl_array_release(s->bvar);
}

// Configuration-space coefficients at reference point z of cell cidx.
static void
conf_coeffs(
  struct test_setup *s, long cidx, const double *z, const struct gkyl_array *a_P, double *K,
  double *beta, double *theta_q, double *aP_q
)
{
  int cdim = s->cdim, nc = s->cbasis.num_basis;
  const double *b = gkyl_array_cfetch(s->bvar, cidx), *gu = gkyl_array_cfetch(s->grad_u, cidx);
  const double *gb = gkyl_array_cfetch(s->grad_b, cidx);
  double bq[3], gradu[3][3];
  for (int i = 0; i < 3; ++i) {
    bq[i] = s->cbasis.eval_expand(z, &b[i * nc]);
    for (int d = 0; d < 3; ++d) {
      gradu[d][i] = d < cdim ? s->cbasis.eval_expand(z, &gu[(d * 3 + i) * nc]) : 0.0;
    }
  }
  *K = 0.0;
  *beta = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      *K += bq[i] * bq[j] * gradu[i][j];
    }
  }
  for (int d = 0; d < cdim; ++d) {
    *beta += s->cbasis.eval_expand(z, &gb[(d * 3 + d) * nc]);
  }
  *theta_q = s->cbasis.eval_expand(z, gkyl_array_cfetch(s->theta, cidx));
  *aP_q = s->cbasis.eval_expand(z, gkyl_array_cfetch(a_P, cidx));
}

// Independent reference of the stored v_par flux on the lower face of cell idx
// (idx_l is the cell below). Returns the CFL speed estimate at the nodes.
static double
reference_vel_face(
  struct test_setup *s, const struct gkyl_array *a_P, const int *idx_l, const int *idx,
  double *flux_ref
)
{
  int cdim = s->cdim, pdim = s->pdim;
  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  long cidx = gkyl_range_idx(&s->clocal_ext, idx);
  long pl = gkyl_range_idx(&s->plocal_ext, idx_l), pr = gkyl_range_idx(&s->plocal_ext, idx);
  const double *f_l = gkyl_array_cfetch(s->f, pl), *f_r = gkyl_array_cfetch(s->f, pr);
  double xc[4];
  gkyl_rect_grid_cell_center(&s->pgrid, idx, xc);
  double v_face = xc[cdim] - 0.5 * s->pgrid.dx[cdim];

  for (int k = 0; k < 2 * nc; ++k) {
    flux_ref[k] = 0.0;
  }
  double cfl_max = 0.0;
  int nnodes = 1;
  for (int d = 0; d < cdim; ++d) {
    nnodes *= 3;
  }
  for (int n = 0; n < nnodes; ++n) {
    double zs[4], zl[4], zr[4], w = 1.0;
    int m = n;
    for (int d = cdim - 1; d >= 0; --d) {
      zs[d] = gl3_x[m % 3];
      w *= gl3_w[m % 3];
      m /= 3;
    }
    for (int d = 0; d < cdim; ++d) {
      zl[d] = zr[d] = zs[d];
    }
    zl[cdim] = 1.0;
    zr[cdim] = -1.0;
    double K, beta, thq, aPq;
    conf_coeffs(s, cidx, zs, a_P, &K, &beta, &thq, &aPq);
    double alpha = aPq - v_face * K;
    double Fl = s->pbasis.eval_expand(zl, &f_l[0]), Fr = s->pbasis.eval_expand(zr, &f_r[0]);
    double Gl = s->pbasis.eval_expand(zl, &f_l[np]), Gr = s->pbasis.eval_expand(zr, &f_r[np]);
    double Hl = 4.0 * thq * Gl - 2.0 * thq * thq * Fl, Hr = 4.0 * thq * Gr - 2.0 * thq * thq * Fr;
    double ap = fmax(alpha, 0.0), am = fmin(alpha, 0.0), bp = fmax(beta, 0.0), bm = fmin(beta, 0.0);
    double ghat[2] = {ap * Fl + am * Fr + bp * Gl + bm * Gr, ap * Gl + am * Gr + bp * Hl + bm * Hr};
    cfl_max = fmax(cfl_max, fabs(alpha) + (2.0 + sqrt(2.0)) * fabs(thq * beta));
    for (int c = 0; c < 2; ++c) {
      for (int k = 0; k < nc; ++k) {
        double unit[8] = {0.0};
        unit[k] = 1.0;
        double znorm[3] = {0.4, 0.4, 0.4};
        double chi_k = s->cbasis.eval_expand(zs, unit) / s->cbasis.eval_expand(znorm, unit);
        flux_ref[c * nc + k] += w * chi_k * ghat[c];
      }
    }
  }
  return cfl_max;
}

static void
test_vel_flux(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis;

  struct gkyl_dg_pkpm_vel_flux_surf_inp inp = {
    .phase_grid = &s.pgrid,
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .use_gpu = false,
  };
  struct gkyl_dg_pkpm_vel_flux_surf *up = gkyl_dg_pkpm_vel_flux_surf_inew(&inp);
  int nk = gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(up);
  TEST_CHECK(nk == nc);
  struct gkyl_array *flux = mkarr(2 * nk, s.plocal_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, s.plocal_ext.volume);
  gkyl_array_clear(flux, 1.0e30); // the edge faces must be overwritten with zeros
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_pkpm_vel_flux_surf_advance(
    up, &s.clocal, &s.plocal, s.bvar, s.grad_b, s.u, s.grad_u, s.theta, s.a_P, s.f, cflrate, flux
  );

  double flux_ref[16];
  double max_err = 0.0, max_flux = 0.0, max_cfl_err = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    const double *fl = gkyl_array_cfetch(flux, pidx);
    const double *cfl = gkyl_array_cfetch(cflrate, pidx);
    if (iter.idx[cdim] == s.plocal.lower[cdim]) {
      for (int k = 0; k < 2 * nk; ++k) {
        TEST_CHECK(fl[k] == 0.0);
      }
      TEST_CHECK(cfl[0] == 0.0);
      continue;
    }
    int idx_l[4];
    gkyl_copy_int_arr(s.pdim, iter.idx, idx_l);
    idx_l[cdim] -= 1;
    double cfl_max = reference_vel_face(&s, s.a_P, idx_l, iter.idx, flux_ref);
    for (int k = 0; k < 2 * nk; ++k) {
      double err = fabs(fl[k] - flux_ref[k]);
      max_err = err > max_err ? err : max_err;
      max_flux = fabs(flux_ref[k]) > max_flux ? fabs(flux_ref[k]) : max_flux;
    }
    double cerr = fabs(cfl[0] - 2.5 * 2.0 / s.pgrid.dx[cdim] * cfl_max);
    max_cfl_err = cerr > max_cfl_err ? cerr : max_cfl_err;
  }
  TEST_CHECK(max_flux > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_flux);
  TEST_MSG("cdim=%d: max flux error %g (max flux %g)", cdim, max_err, max_flux);
  TEST_CHECK(max_cfl_err < 1e-11);
  TEST_MSG("cdim=%d: max cflrate error %g", cdim, max_cfl_err);

  gkyl_array_release(cflrate);
  gkyl_array_release(flux);
  gkyl_dg_pkpm_vel_flux_surf_release(up);
  teardown(&s);
}

static void
test_vel_vol(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int np = s.pbasis.num_basis;

  struct gkyl_dg_pkpm_kinetic_vel_inp inp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields = {.a_P = s.a_P, .vel_flux_surf = NULL},
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn = gkyl_dg_pkpm_kinetic_vel_inew(&inp);

  double max_err = 0.0, max_ref = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    long cidx = gkyl_range_idx(&s.clocal_ext, iter.idx);
    const double *f = gkyl_array_cfetch(s.f, pidx);
    const double *aP = gkyl_array_cfetch(s.a_P, cidx);
    double xc[4], out[128] = {0.0}, ref[128] = {0.0};
    gkyl_rect_grid_cell_center(&s.pgrid, iter.idx, xc);
    eqn->vol_term(eqn, xc, s.pgrid.dx, iter.idx, f, out);

    int nnodes = 1;
    for (int d = 0; d < s.pdim; ++d) {
      nnodes *= 4;
    }
    double dv1 = 2.0 / s.pgrid.dx[cdim];
    for (int n = 0; n < nnodes; ++n) {
      double z[4], w = 1.0;
      int m = n;
      for (int d = s.pdim - 1; d >= 0; --d) {
        z[d] = gl4_x[m % 4];
        w *= gl4_w[m % 4];
        m /= 4;
      }
      double aPq = s.cbasis.eval_expand(z, aP);
      double Fq = s.pbasis.eval_expand(z, &f[0]), Gq = s.pbasis.eval_expand(z, &f[np]);
      for (int k = 0; k < np; ++k) {
        double unit[64] = {0.0};
        unit[k] = 1.0;
        double dphi_v = dv1 * s.pbasis.eval_grad_expand(cdim, z, unit);
        ref[k] += w * dphi_v * aPq * Fq;
        ref[np + k] += w * dphi_v * aPq * Gq;
      }
    }
    for (int k = 0; k < 2 * np; ++k) {
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

// The complete two-pass right-hand side: configuration fluxes, stage A, M1 closure,
// v_par fluxes, stage B. The first v_par moment of the total must vanish in every
// cell and configuration mode, to round-off relative to the cancelling stage-A moment.
static void
test_m1_invariant(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis;

  // Stage A.
  struct gkyl_dg_pkpm_conf_flux_surf_inp cinp = {
    .phase_grid = &s.pgrid, .conf_basis = &s.cbasis, .phase_basis = &s.pbasis, .use_gpu = false};
  struct gkyl_dg_pkpm_conf_flux_surf *calc_cflux = gkyl_dg_pkpm_conf_flux_surf_inew(&cinp);
  int nkc = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(calc_cflux);
  struct gkyl_array *cflux = mkarr(2 * cdim * nkc, s.plocal_ext.volume);
  struct gkyl_array *cmflux = mkarr(cdim * gkyl_dg_pkpm_conf_flux_surf_num_face_basis(calc_cflux), s.plocal_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, s.plocal_ext.volume);
  struct gkyl_array *rhs = mkarr(2 * np, s.plocal_ext.volume);
  gkyl_array_clear(cflux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);
  gkyl_dg_pkpm_conf_flux_surf_advance(
    calc_cflux, &s.clocal, &s.plocal, &s.plocal_ext, s.bvar, s.u, s.theta, s.f, cflrate, cflux, cmflux
  );
  struct gkyl_dg_pkpm_kinetic_conf_inp ainp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields =
      {
        .bvar = s.bvar,
        .grad_b = s.grad_b,
        .u = s.u,
        .grad_u = s.grad_u,
        .theta = s.theta,
        .nu_sum = s.nu_sum,
        .nu_vthsq_sum = s.nu_vthsq_sum,
        .conf_flux_surf = cflux,
      },
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn_a = gkyl_dg_pkpm_kinetic_conf_inew(&ainp);
  int up_dirs_a[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux_a[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *slvr_a =
    gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn_a, cdim, up_dirs_a, zero_flux_a, 1, false);
  gkyl_hyper_dg_advance(slvr_a, &s.plocal, s.f, cflrate, rhs);

  // Closure.
  double mass = 2.0;
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_pkpm_new(
    &s.pgrid, &s.cbasis, &s.pbasis, &s.clocal, &s.vlocal, mass, false, false
  );
  struct gkyl_array *moms = mkarr(4 * nc, s.clocal_ext.volume);
  struct gkyl_array *rhs_moms = mkarr(4 * nc, s.clocal_ext.volume);
  gkyl_array_clear(moms, 0.0);
  gkyl_array_clear(rhs_moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, s.f, moms);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, rhs, rhs_moms);

  struct gkyl_wv_eqn *ten = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_new(&s.cgrid, &s.clocal_ext, 0, 0, false);
  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_new(&s.cgrid, &s.cbasis, &s.clocal_ext, ten, geom, 0.0, false);
  struct gkyl_array *a_P = mkarr(nc, s.clocal_ext.volume);
  gkyl_array_clear(a_P, 0.0);
  int num_null = gkyl_dg_calc_pkpm_vars_p_force(vars, &s.clocal, moms, rhs_moms, a_P);
  TEST_CHECK(num_null == 0);

  double max_m1_a = 0.0;
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &s.clocal);
  while (gkyl_range_iter_next(&citer)) {
    const double *rm = gkyl_array_cfetch(rhs_moms, gkyl_range_idx(&s.clocal_ext, citer.idx));
    for (int k = 0; k < nc; ++k) {
      max_m1_a = fabs(rm[3 * nc + k]) > max_m1_a ? fabs(rm[3 * nc + k]) : max_m1_a;
    }
  }

  // Stage B.
  struct gkyl_dg_pkpm_vel_flux_surf_inp vinp = {
    .phase_grid = &s.pgrid, .conf_basis = &s.cbasis, .phase_basis = &s.pbasis, .use_gpu = false};
  struct gkyl_dg_pkpm_vel_flux_surf *calc_vflux = gkyl_dg_pkpm_vel_flux_surf_inew(&vinp);
  struct gkyl_array *vflux = mkarr(2 * nc, s.plocal_ext.volume);
  gkyl_array_clear(vflux, 0.0);
  gkyl_dg_pkpm_vel_flux_surf_advance(
    calc_vflux, &s.clocal, &s.plocal, s.bvar, s.grad_b, s.u, s.grad_u, s.theta, a_P, s.f, cflrate,
    vflux
  );
  struct gkyl_dg_pkpm_kinetic_vel_inp binp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields = {.a_P = a_P, .vel_flux_surf = vflux},
    .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn_b = gkyl_dg_pkpm_kinetic_vel_inew(&binp);
  int up_dirs_b[GKYL_MAX_DIM] = {cdim}, zero_flux_b[2 * GKYL_MAX_DIM] = {0};
  zero_flux_b[cdim] = zero_flux_b[cdim + s.pdim] = 1;
  gkyl_hyper_dg *slvr_b =
    gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn_b, 1, up_dirs_b, zero_flux_b, 1, false);
  gkyl_hyper_dg_advance(slvr_b, &s.plocal, s.f, cflrate, rhs);

  // First moment of the total right-hand side.
  gkyl_array_clear(rhs_moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, rhs, rhs_moms);
  double max_m1 = 0.0;
  gkyl_range_iter_init(&citer, &s.clocal);
  while (gkyl_range_iter_next(&citer)) {
    const double *rm = gkyl_array_cfetch(rhs_moms, gkyl_range_idx(&s.clocal_ext, citer.idx));
    for (int k = 0; k < nc; ++k) {
      max_m1 = fabs(rm[3 * nc + k]) > max_m1 ? fabs(rm[3 * nc + k]) : max_m1;
    }
  }
  TEST_CHECK(max_m1_a > 0.0);
  TEST_CHECK(max_m1 < 1e-12 * max_m1_a);
  TEST_MSG("cdim=%d: max |M1(rhs)| %g after closure (stage-A moment %g)", cdim, max_m1, max_m1_a);

  gkyl_hyper_dg_release(slvr_b);
  gkyl_dg_eqn_release(eqn_b);
  gkyl_array_release(vflux);
  gkyl_dg_pkpm_vel_flux_surf_release(calc_vflux);
  gkyl_array_release(a_P);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(ten);
  gkyl_array_release(rhs_moms);
  gkyl_array_release(moms);
  gkyl_dg_updater_moment_pkpm_release(mcalc);
  gkyl_hyper_dg_release(slvr_a);
  gkyl_dg_eqn_release(eqn_a);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(cflux);
  gkyl_array_release(cmflux);
  gkyl_dg_pkpm_conf_flux_surf_release(calc_cflux);
  teardown(&s);
}

void test_vel_flux_1x() { test_vel_flux(1); }
void test_vel_flux_2x() { test_vel_flux(2); }
void test_vel_flux_3x() { test_vel_flux(3); }
void test_vel_vol_1x() { test_vel_vol(1); }
void test_vel_vol_2x() { test_vel_vol(2); }
void test_vel_vol_3x() { test_vel_vol(3); }
void test_m1_invariant_1x() { test_m1_invariant(1); }
void test_m1_invariant_2x() { test_m1_invariant(2); }
void test_m1_invariant_3x() { test_m1_invariant(3); }

TEST_LIST = {
  {"test_vel_flux_1x", test_vel_flux_1x},
  {"test_vel_flux_2x", test_vel_flux_2x},
  {"test_vel_flux_3x", test_vel_flux_3x},
  {"test_vel_vol_1x", test_vel_vol_1x},
  {"test_vel_vol_2x", test_vel_vol_2x},
  {"test_vel_vol_3x", test_vel_vol_3x},
  {"test_m1_invariant_1x", test_m1_invariant_1x},
  {"test_m1_invariant_2x", test_m1_invariant_2x},
  {"test_m1_invariant_3x", test_m1_invariant_3x},
  {NULL, NULL},
};
