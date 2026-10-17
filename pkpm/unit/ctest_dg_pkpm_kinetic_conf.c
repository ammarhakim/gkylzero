// Tests of the configuration-space pass (stage A) of the refactored PKPM kinetic
// update: the exact modal volume kernel against an independent Gauss-Legendre
// quadrature of the same terms, and the a_P force closure (weak residual of
// rho a_P + M1(rhs_A) = 0 in every mode).
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
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
  fout[0] = F0 * (1.0 + 0.1 * v);
  fout[1] = 0.5 * vt * vt * F0 * (1.0 - 0.1 * v);
}

static const double gl4_x[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
  0.8611363115940526};
static const double gl4_w[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
  0.3478548451374538};

struct test_setup {
  int cdim, pdim;
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext, vlocal, vlocal_ext;
  struct gkyl_array *bvar, *grad_b, *u, *grad_u, *theta, *nu_sum, *nu_vthsq_sum, *f;
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
  struct gkyl_array *scal = mkarr(3 * nc, s->clocal_ext.volume);
  s->bvar = mkarr(3 * nc, s->clocal_ext.volume);
  s->u = mkarr(3 * nc, s->clocal_ext.volume);
  s->theta = mkarr(nc, s->clocal_ext.volume);
  s->nu_sum = mkarr(nc, s->clocal_ext.volume);
  s->nu_vthsq_sum = mkarr(nc, s->clocal_ext.volume);
  s->f = mkarr(2 * np, s->plocal_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, em);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_u, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->u);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_conf_scalars, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, scal);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_set_offset(s->theta, 1.0, scal, 0 * nc);
  gkyl_array_set_offset(s->nu_sum, 1.0, scal, 1 * nc);
  gkyl_array_set_offset(s->nu_vthsq_sum, 1.0, scal, 2 * nc);
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
  gkyl_array_release(s->nu_vthsq_sum);
  gkyl_array_release(s->nu_sum);
  gkyl_array_release(s->theta);
  gkyl_array_release(s->u);
  gkyl_array_release(s->bvar);
}

static struct gkyl_dg_eqn *
make_eqn(struct test_setup *s, const struct gkyl_array *flux)
{
  struct gkyl_dg_pkpm_kinetic_conf_inp inp = {
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .conf_range = &s->clocal,
    .phase_range = &s->plocal,
    .auxfields =
      {
        .bvar = s->bvar,
        .grad_b = s->grad_b,
        .u = s->u,
        .grad_u = s->grad_u,
        .theta = s->theta,
        .nu_sum = s->nu_sum,
        .nu_vthsq_sum = s->nu_vthsq_sum,
        .conf_flux_surf = flux,
      },
    .use_gpu = false,
  };
  return gkyl_dg_pkpm_kinetic_conf_inew(&inp);
}

// Independent evaluation of the stage-A volume terms of one cell by 4-point
// Gauss-Legendre quadrature in every phase-space direction (exact for every
// polynomial product of the kernel).
static void
reference_vol(struct test_setup *s, const int *idx, double *out_ref)
{
  int cdim = s->cdim, pdim = s->pdim;
  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  long cidx = gkyl_range_idx(&s->clocal_ext, idx), pidx = gkyl_range_idx(&s->plocal_ext, idx);
  const double *b = gkyl_array_cfetch(s->bvar, cidx), *u = gkyl_array_cfetch(s->u, cidx);
  const double *gu = gkyl_array_cfetch(s->grad_u, cidx), *gb = gkyl_array_cfetch(s->grad_b, cidx);
  const double *th = gkyl_array_cfetch(s->theta, cidx), *nu = gkyl_array_cfetch(s->nu_sum, cidx);
  const double *nuvt = gkyl_array_cfetch(s->nu_vthsq_sum, cidx);
  const double *f = gkyl_array_cfetch(s->f, pidx);
  double xc[4];
  gkyl_rect_grid_cell_center(&s->pgrid, idx, xc);
  double wv = xc[cdim], dv = s->pgrid.dx[cdim];
  double dx1[4];
  for (int d = 0; d < pdim; ++d) {
    dx1[d] = 2.0 / s->pgrid.dx[d];
  }
  for (int k = 0; k < 2 * np; ++k) {
    out_ref[k] = 0.0;
  }
  int nnodes = 1;
  for (int d = 0; d < pdim; ++d) {
    nnodes *= 4;
  }
  for (int n = 0; n < nnodes; ++n) {
    double z[4], w = 1.0;
    int m = n;
    for (int d = pdim - 1; d >= 0; --d) {
      z[d] = gl4_x[m % 4];
      w *= gl4_w[m % 4];
      m /= 4;
    }
    double v = wv + 0.5 * dv * z[cdim];
    double bq[3], uq[3], gradu[3][3];
    for (int i = 0; i < 3; ++i) {
      bq[i] = s->cbasis.eval_expand(z, &b[i * nc]);
      uq[i] = s->cbasis.eval_expand(z, &u[i * nc]);
      for (int d = 0; d < 3; ++d) {
        gradu[d][i] = d < cdim ? s->cbasis.eval_expand(z, &gu[(d * 3 + i) * nc]) : 0.0;
      }
    }
    double K = 0.0, divu = 0.0, divb = 0.0;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        K += bq[i] * bq[j] * gradu[i][j];
      }
    }
    for (int d = 0; d < cdim; ++d) {
      divu += gradu[d][d];
      divb += s->cbasis.eval_expand(z, &gb[(d * 3 + d) * nc]);
    }
    double thq = s->cbasis.eval_expand(z, th), nuq = s->cbasis.eval_expand(z, nu);
    double nuvtq = s->cbasis.eval_expand(z, nuvt);
    double Fq = s->pbasis.eval_expand(z, &f[0]), Gq = s->pbasis.eval_expand(z, &f[np]);
    double Hq = 4.0 * thq * Gq - 2.0 * thq * thq * Fq;
    double src_G = (K - divu - 2.0 * nuq - v * divb) * Gq + 2.0 * nuvtq * Fq;

    for (int k = 0; k < np; ++k) {
      double unit[64] = {0.0};
      unit[k] = 1.0;
      double phi = s->pbasis.eval_expand(z, unit);
      double dphi_v = dx1[cdim] * s->pbasis.eval_grad_expand(cdim, z, unit);
      double stream_F = 0.0, stream_G = 0.0;
      for (int d = 0; d < cdim; ++d) {
        double dphi = dx1[d] * s->pbasis.eval_grad_expand(d, z, unit);
        stream_F += dphi * (v * bq[d] + uq[d]) * Fq;
        stream_G += dphi * (v * bq[d] + uq[d]) * Gq;
      }
      out_ref[k] += w * (stream_F + dphi_v * (-v * K * Fq + divb * Gq));
      out_ref[np + k] += w * (stream_G + dphi_v * (-v * K * Gq + divb * Hq) + phi * src_G);
    }
  }
}

static void
test_vol(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int np = s.pbasis.num_basis;

  struct gkyl_dg_eqn *eqn = make_eqn(&s, NULL);
  struct gkyl_array *rhs = mkarr(2 * np, s.plocal_ext.volume);
  gkyl_array_clear(rhs, 0.0);

  double out_ref[128];
  double max_err = 0.0, max_ref = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    double xc[4];
    gkyl_rect_grid_cell_center(&s.pgrid, iter.idx, xc);
    double *out = gkyl_array_fetch(rhs, pidx);
    double cfl = eqn->vol_term(eqn, xc, s.pgrid.dx, iter.idx, gkyl_array_cfetch(s.f, pidx), out);
    TEST_CHECK(cfl >= 0.0);
    reference_vol(&s, iter.idx, out_ref);
    for (int k = 0; k < 2 * np; ++k) {
      double err = fabs(out[k] - out_ref[k]);
      max_err = err > max_err ? err : max_err;
      max_ref = fabs(out_ref[k]) > max_ref ? fabs(out_ref[k]) : max_ref;
    }
  }
  TEST_CHECK(max_ref > 0.0);
  TEST_CHECK(max_err < 1e-11 * max_ref);
  TEST_MSG("cdim=%d: max volume error %g (max reference %g)", cdim, max_err, max_ref);

  gkyl_array_release(rhs);
  gkyl_dg_eqn_release(eqn);
  teardown(&s);
}

// Full stage-A right-hand side (fluxes + volume) through hyper_dg, then the a_P
// closure: the weak residual <psi_k, rho a_P + M1(rhs_A)> vanishes in every mode.
static void
test_p_force(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis;

  struct gkyl_dg_pkpm_conf_flux_surf_inp finp = {
    .phase_grid = &s.pgrid,
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .use_gpu = false,
  };
  struct gkyl_dg_pkpm_conf_flux_surf *calc_flux = gkyl_dg_pkpm_conf_flux_surf_inew(&finp);
  int nk = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(calc_flux);
  struct gkyl_array *flux = mkarr(2 * cdim * nk, s.plocal_ext.volume);
  struct gkyl_array *mflux = mkarr(cdim * gkyl_dg_pkpm_conf_flux_surf_num_face_basis(calc_flux), s.plocal_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, s.plocal_ext.volume);
  gkyl_array_clear(flux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_pkpm_conf_flux_surf_advance(
    calc_flux, &s.clocal, &s.plocal, &s.plocal_ext, s.bvar, s.u, s.theta, s.f, cflrate, flux, mflux
  );

  struct gkyl_dg_eqn *eqn = make_eqn(&s, flux);
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *slvr = gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn, cdim, up_dirs, zero_flux, 1, false);
  struct gkyl_array *rhs = mkarr(2 * np, s.plocal_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_hyper_dg_advance(slvr, &s.plocal, s.f, cflrate, rhs);

  // Moments of the stage state and of the partial right-hand side.
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

  // Weak residual: rho a_P (weak product) + M1(rhs) in every mode.
  struct gkyl_array *rho_aP = mkarr(nc, s.clocal_ext.volume);
  gkyl_array_clear(rho_aP, 0.0);
  gkyl_dg_mul_op_range(&s.cbasis, 0, rho_aP, 0, moms, 0, a_P, &s.clocal);
  double max_res = 0.0, max_m1 = 0.0, max_aP = 0.0;
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &s.clocal);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, citer.idx);
    const double *ra = gkyl_array_cfetch(rho_aP, cidx);
    const double *rm = gkyl_array_cfetch(rhs_moms, cidx);
    const double *ap = gkyl_array_cfetch(a_P, cidx);
    for (int k = 0; k < nc; ++k) {
      double res = fabs(ra[k] + rm[3 * nc + k]);
      max_res = res > max_res ? res : max_res;
      max_m1 = fabs(rm[3 * nc + k]) > max_m1 ? fabs(rm[3 * nc + k]) : max_m1;
      max_aP = fabs(ap[k]) > max_aP ? fabs(ap[k]) : max_aP;
    }
  }
  TEST_CHECK(max_m1 > 0.0);
  TEST_CHECK(max_aP > 0.0);
  TEST_CHECK(max_res < 1e-12 * max_m1);
  TEST_MSG("cdim=%d: weak residual %g (max |M1(rhs)| %g, max |a_P| %g)", cdim, max_res, max_m1, max_aP);

  gkyl_array_release(rho_aP);
  gkyl_array_release(a_P);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_wave_geom_release(geom);
  gkyl_wv_eqn_release(ten);
  gkyl_array_release(rhs_moms);
  gkyl_array_release(moms);
  gkyl_dg_updater_moment_pkpm_release(mcalc);
  gkyl_array_release(rhs);
  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_array_release(cflrate);
  gkyl_array_release(flux);
  gkyl_array_release(mflux);
  gkyl_dg_pkpm_conf_flux_surf_release(calc_flux);
  teardown(&s);
}

// The weak gradient kernels against an independent quadrature of the weak form
//   <psi_k, d_i u_j> = dx1_i [ <psi_k(xi_i=+1), ubar_r> - <psi_k(xi_i=-1), ubar_l> - <d_xi_i psi_k, u_c> ]
// with ubar the average of the two traces at the face.
static void
test_grad_u(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis;
  double max_err = 0.0, max_val = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.clocal);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, iter.idx);
    const double *gu = gkyl_array_cfetch(s.grad_u, cidx), *uc = gkyl_array_cfetch(s.u, cidx);
    for (int dir = 0; dir < cdim; ++dir) {
      int idx_l[3], idx_r[3];
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      gkyl_copy_int_arr(cdim, iter.idx, idx_r);
      idx_l[dir] -= 1;
      idx_r[dir] += 1;
      const double *ul = gkyl_array_cfetch(s.u, gkyl_range_idx(&s.clocal_ext, idx_l));
      const double *ur = gkyl_array_cfetch(s.u, gkyl_range_idx(&s.clocal_ext, idx_r));
      double dx1 = 2.0 / s.cgrid.dx[dir];
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < nc; ++k) {
          double unit[8] = {0.0};
          unit[k] = 1.0;
          double ref = 0.0;
          // Volume term: 4-point Gauss-Legendre in every direction.
          int nv = 1;
          for (int d = 0; d < cdim; ++d) {
            nv *= 4;
          }
          for (int n = 0; n < nv; ++n) {
            double z[3], w = 1.0;
            int m = n;
            for (int d = cdim - 1; d >= 0; --d) {
              z[d] = gl4_x[m % 4];
              w *= gl4_w[m % 4];
              m /= 4;
            }
            ref -= w * dx1 * s.cbasis.eval_grad_expand(dir, z, unit) * s.cbasis.eval_expand(z, &uc[j * nc]);
          }
          // Face terms: 4-point rule in the transverse directions.
          int nf = 1;
          for (int d = 0; d < cdim - 1; ++d) {
            nf *= 4;
          }
          for (int n = 0; n < nf; ++n) {
            double zp[3], zm[3], w = 1.0;
            int m = n, t = 0;
            double zt[3];
            for (int d = cdim - 2; d >= 0; --d) {
              zt[d] = gl4_x[m % 4];
              w *= gl4_w[m % 4];
              m /= 4;
            }
            for (int d = 0; d < cdim; ++d) {
              if (d == dir) {
                zp[d] = 1.0;
                zm[d] = -1.0;
              } else {
                zp[d] = zm[d] = zt[t++];
              }
            }
            // Upper face of this cell: average of this cell at +1 and the right cell at -1.
            double ubar_r = 0.5 * (s.cbasis.eval_expand(zp, &uc[j * nc]) + s.cbasis.eval_expand(zm, &ur[j * nc]));
            double ubar_l = 0.5 * (s.cbasis.eval_expand(zp, &ul[j * nc]) + s.cbasis.eval_expand(zm, &uc[j * nc]));
            ref += w * dx1 * (s.cbasis.eval_expand(zp, unit) * ubar_r - s.cbasis.eval_expand(zm, unit) * ubar_l);
          }
          double val = gu[(dir * 3 + j) * nc + k];
          max_err = fmax(max_err, fabs(val - ref));
          max_val = fmax(max_val, fabs(ref));
        }
      }
    }
  }
  TEST_CHECK(max_val > 0.0);
  TEST_CHECK(max_err < 1e-11 * max_val);
  TEST_MSG("cdim=%d: weak gradient max error %g (max value %g)", cdim, max_err, max_val);
  teardown(&s);
}

void test_grad_u_1x() { test_grad_u(1); }
void test_grad_u_2x() { test_grad_u(2); }
void test_grad_u_3x() { test_grad_u(3); }

void test_vol_1x() { test_vol(1); }
void test_vol_2x() { test_vol(2); }
void test_vol_3x() { test_vol(3); }
void test_p_force_1x() { test_p_force(1); }
void test_p_force_2x() { test_p_force(2); }
void test_p_force_3x() { test_p_force(3); }

TEST_LIST = {
  {"test_grad_u_1x", test_grad_u_1x},
  {"test_grad_u_2x", test_grad_u_2x},
  {"test_grad_u_3x", test_grad_u_3x},
  {"test_vol_1x", test_vol_1x},
  {"test_vol_2x", test_vol_2x},
  {"test_vol_3x", test_vol_3x},
  {"test_p_force_1x", test_p_force_1x},
  {"test_p_force_2x", test_p_force_2x},
  {"test_p_force_3x", test_p_force_3x},
  {NULL, NULL},
};
