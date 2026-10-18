// Tests of the configuration-space pass (stage A) of the refactored PKPM kinetic
// update: the exact modal volume kernel against an independent Gauss-Legendre
// quadrature of the same terms, and the a_P force closure (weak residual of
// M0(f) c_0 + M0(rhs_A) = 0 and M0(f) a_P + M1(rhs_A) + c_0 M1(f) = 0 in every mode).
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

// Identity configuration-space position map (the updaters always take one).
static struct gkyl_vlasov_position_map *
identity_pos_map(
  const struct gkyl_rect_grid *cgrid, const struct gkyl_range *clocal, const struct gkyl_range *clocal_ext,
  const struct gkyl_basis *cbasis
)
{
  struct gkyl_vlasov_position_map_inp inp[GKYL_MAX_CDIM] = {0};
  return gkyl_vlasov_position_map_new(cgrid, clocal, clocal_ext, cbasis, inp, false);
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

// Phase-space ranges with configuration ghosts on both sides and ONE velocity ghost cell
// above w_max (the empty state beyond the characteristic velocity edge, whose flux slot
// holds the upper-edge flux); plocal_ext_conf is extended in configuration space only.
static void
phase_ranges(
  const struct gkyl_rect_grid *pgrid, int cdim, struct gkyl_range *plocal_ext, struct gkyl_range *plocal,
  struct gkyl_range *plocal_ext_conf
)
{
  int pdim = pgrid->ndim, lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int elo[GKYL_MAX_DIM] = {0}, eup[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < pdim; ++d) {
    lower[d] = 1;
    upper[d] = pgrid->cells[d];
  }
  for (int d = 0; d < cdim; ++d) {
    elo[d] = eup[d] = 1;
  }
  eup[cdim] = 1;
  struct gkyl_range local;
  gkyl_range_init(&local, pdim, lower, upper);
  gkyl_range_extend(plocal_ext, &local, elo, eup);
  gkyl_sub_range_init(plocal, plocal_ext, lower, upper);
  int lo_c[GKYL_MAX_DIM], up_c[GKYL_MAX_DIM];
  for (int d = 0; d < pdim; ++d) {
    lo_c[d] = d < cdim ? plocal_ext->lower[d] : lower[d];
    up_c[d] = d < cdim ? plocal_ext->upper[d] : upper[d];
  }
  gkyl_sub_range_init(plocal_ext_conf, plocal_ext, lo_c, up_c);
}

struct test_setup {
  struct gkyl_vlasov_position_map *pos_map; // identity configuration-space map
  int cdim, pdim;
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext, plocal_ext_conf, vlocal, vlocal_ext;
  struct gkyl_array *bvar, *grad_b, *u, *grad_u, *grad_ln_n, *theta, *f;
  struct gkyl_array *xi, *grad_xi;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
};

// Thermal speed [xi, xi^2, p] (positive; p is not used by the kinetic kernels).
static void
eval_xi(double t, const double *xn, double *restrict fout, void *ctx)
{
  int cdim = *(int *)ctx;
  double x = xn[0], y = cdim > 1 ? xn[1] : 0.3, z = cdim > 2 ? xn[2] : -0.2;
  double xi = 1.0 + 0.2 * cos(x) * sin(y) + 0.05 * z;
  fout[0] = xi;
  fout[1] = xi * xi;
  fout[2] = 0.9 * xi * xi;
}

// Weak gradient of u on the interior cells (u is projected on the extended range).
static void
compute_grad_u(struct test_setup *s)
{
  int cdim = s->cdim, nc = s->cbasis.num_basis;
  s->grad_u = mkarr(cdim * 3 * nc, s->clocal_ext.volume);
  s->grad_xi = mkarr(cdim * 3 * nc, s->clocal_ext.volume);
  s->grad_ln_n = mkarr(cdim * nc, s->clocal_ext.volume);
  gkyl_array_clear(s->grad_ln_n, 0.0);
  gkyl_array_clear(s->grad_u, 0.0);
  gkyl_array_clear(s->grad_xi, 0.0);
  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_new(&s->cgrid, &s->cbasis, &s->clocal_ext, s->pos_map, false);
  gkyl_dg_calc_pkpm_vars_grad_u(vars, &s->clocal, s->u, s->grad_u);
  gkyl_dg_calc_pkpm_vars_grad_u(vars, &s->clocal, s->xi, s->grad_xi);
  gkyl_dg_calc_pkpm_vars_release(vars);
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
  int cghost[3] = {1, 1, 1}, vghost[1] = {0};
  if (cdim < 3) {
    lower[cdim] = -4.0;
    upper[cdim] = 4.0;
    cells[cdim] = 8;
  }
  gkyl_rect_grid_init(&s->cgrid, cdim, lower, upper, ccells);
  gkyl_rect_grid_init(&s->pgrid, pdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&s->pbasis, cdim, 1);
  gkyl_create_grid_ranges(&s->cgrid, cghost, &s->clocal_ext, &s->clocal);
  s->pos_map = identity_pos_map(&s->cgrid, &s->clocal, &s->clocal_ext, &s->cbasis);
  phase_ranges(&s->pgrid, cdim, &s->plocal_ext, &s->plocal, &s->plocal_ext_conf);
  struct gkyl_rect_grid vgrid;
  gkyl_rect_grid_init(&vgrid, 1, &lower[cdim], &upper[cdim], &cells[cdim]);
  gkyl_create_grid_ranges(&vgrid, vghost, &s->vlocal_ext, &s->vlocal);

  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  struct gkyl_array *em = mkarr(8 * nc, s->clocal_ext.volume);
  struct gkyl_array *scal = mkarr(nc, s->clocal_ext.volume);
  s->bvar = mkarr(3 * nc, s->clocal_ext.volume);
  s->u = mkarr(3 * nc, s->clocal_ext.volume);
  s->xi = mkarr(3 * nc, s->clocal_ext.volume);
  s->theta = mkarr(nc, s->clocal_ext.volume);
  s->f = mkarr(2 * np, s->plocal_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, em);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_u, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->u);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 1, eval_conf_scalars, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, scal);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_set_offset(s->theta, 1.0, scal, 0 * nc);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_xi, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->xi);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->pgrid, &s->pbasis, 8, 2, eval_f, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->plocal_ext_conf, s->f);
  gkyl_proj_on_basis_release(proj);

  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->cgrid, &s->cbasis, s->pos_map, false);
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
  gkyl_vlasov_position_map_release(s->pos_map);
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_array_release(s->grad_b);
  gkyl_array_release(s->grad_u);
  gkyl_array_release(s->f);
  gkyl_array_release(s->theta);
  gkyl_array_release(s->xi);
  gkyl_array_release(s->grad_xi);
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
        .xi = s->xi,
        .grad_xi = s->grad_xi,
        .grad_ln_n = s->grad_ln_n,
        .theta = s->theta,
        .conf_flux_surf = flux,
      },
    .pos_map = s->pos_map, .use_gpu = false,
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
  const double *gln = gkyl_array_cfetch(s->grad_ln_n, cidx), *gxi = gkyl_array_cfetch(s->grad_xi, cidx);
  const double *xi = gkyl_array_cfetch(s->xi, cidx);
  const double *th = gkyl_array_cfetch(s->theta, cidx);
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
    double K = 0.0, divb = 0.0, bgn = 0.0, bgx = 0.0;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        K += bq[i] * bq[j] * gradu[i][j];
      }
    }
    for (int d = 0; d < cdim; ++d) {
      divb += s->cbasis.eval_expand(z, &gb[(d * 3 + d) * nc]);
      bgn += bq[d] * s->cbasis.eval_expand(z, &gln[d * nc]);
      bgx += bq[d] * s->cbasis.eval_expand(z, &gxi[(d * 3 + 0) * nc]);
    }
    double xiq = s->cbasis.eval_expand(z, xi);
    double thq = s->cbasis.eval_expand(z, th);
    double Fq = s->pbasis.eval_expand(z, &f[0]), Gq = s->pbasis.eval_expand(z, &f[np]);
    double Hq = 4.0 * thq * Gq - 2.0 * thq * thq * Fq;
    double src_F = -xiq * v * bgn * Fq;
    double src_G = (K - xiq * v * divb - 2.0 * v * bgx - xiq * v * bgn) * Gq;

    for (int k = 0; k < np; ++k) {
      double unit[64] = {0.0};
      unit[k] = 1.0;
      double phi = s->pbasis.eval_expand(z, unit);
      double dphi_v = dx1[cdim] * s->pbasis.eval_grad_expand(cdim, z, unit);
      double stream_F = 0.0, stream_G = 0.0;
      for (int d = 0; d < cdim; ++d) {
        double dphi = dx1[d] * s->pbasis.eval_grad_expand(d, z, unit);
        stream_F += dphi * (xiq * v * bq[d] + uq[d]) * Fq;
        stream_G += dphi * (xiq * v * bq[d] + uq[d]) * Gq;
      }
      out_ref[k] += w * (stream_F + dphi_v * (-(v * K + v * v * bgx) * Fq + xiq * divb * Gq) + phi * src_F);
      out_ref[np + k] += w * (stream_G + dphi_v * (-(v * K + v * v * bgx) * Gq + xiq * divb * Hq) + phi * src_G);
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

// Full stage-A right-hand side (fluxes + volume) through hyper_dg, then the closures of
// the factored update: the weak residuals <psi_k, M0(f) c_0 + M0(rhs_A)> and
// <psi_k, M0(f) a_P + M1(rhs_A) + c_0 M1(f)> vanish in every mode (f need not be normalized
// for this: the closures are weak divisions by its zeroth moment).
static void
test_closure(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis;
  // No mirror force: its velocity-edge flux (xi div b) H is switched by the sign of
  // xi div b alone, which would spoil the edge-free closure checked here.
  gkyl_array_clear(s.grad_b, 0.0);

  struct gkyl_dg_pkpm_conf_flux_surf_inp finp = {
    .phase_grid = &s.pgrid,
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .pos_map = s.pos_map, .use_gpu = false,
  };
  struct gkyl_dg_pkpm_conf_flux_surf *calc_flux = gkyl_dg_pkpm_conf_flux_surf_inew(&finp);
  int nk = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(calc_flux);
  struct gkyl_array *flux = mkarr(2 * cdim * nk, s.plocal_ext.volume);
  struct gkyl_array *cflrate = mkarr(1, s.plocal_ext.volume);
  gkyl_array_clear(flux, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_dg_pkpm_conf_flux_surf_advance(
    calc_flux, &s.clocal, &s.plocal, &s.plocal_ext, s.bvar, s.u, s.xi, s.f, cflrate, flux
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
  struct gkyl_array *moms = mkarr(6 * nc, s.clocal_ext.volume);
  struct gkyl_array *rhs_moms = mkarr(6 * nc, s.clocal_ext.volume);
  gkyl_array_clear(moms, 0.0);
  gkyl_array_clear(rhs_moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, s.f, moms);
  gkyl_dg_updater_moment_pkpm_advance(mcalc, &s.plocal, &s.clocal, rhs, rhs_moms);

  struct gkyl_dg_calc_pkpm_vars *vars =
    gkyl_dg_calc_pkpm_vars_new(&s.cgrid, &s.cbasis, &s.clocal_ext, s.pos_map, false);
  struct gkyl_array *a_P = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c_h = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c_0 = mkarr(nc, s.clocal_ext.volume);
  // Previous closures alpha_0 = 0, c = -100 make both velocity edges inflow at every node
  // (the known forces are far smaller than 100 w_max), so the edge fluxes vanish and the
  // closure is the edge-free solve whose weak residuals are checked below.
  gkyl_array_clear(a_P, 0.0);
  gkyl_array_clear(c_h, 0.0);
  gkyl_array_shiftc(c_h, -100.0 * pow(sqrt(2.0), cdim), 0);
  gkyl_array_clear(c_0, 0.0);
  struct gkyl_array *vflux = mkarr(2 * nc, s.plocal_ext.volume);
  gkyl_array_clear(vflux, 0.0);
  int num_null = gkyl_dg_calc_pkpm_vars_closure(
    vars, &s.clocal, &s.plocal, &s.pgrid, moms, rhs_moms, s.bvar, s.grad_b, s.grad_u, s.xi, s.grad_xi, s.theta,
    s.f, c_0, a_P, c_h, vflux, cflrate
  );
  // both edge slots (lowest cell, velocity ghost cell) must hold no flux
  double edge_max = 0.0;
  struct gkyl_range_iter eit;
  gkyl_range_iter_init(&eit, &s.clocal);
  while (gkyl_range_iter_next(&eit)) {
    int idx[GKYL_MAX_DIM];
    gkyl_copy_int_arr(cdim, eit.idx, idx);
    for (int e = 0; e < 2; ++e) {
      idx[cdim] = e == 0 ? s.plocal.lower[cdim] : s.plocal.upper[cdim] + 1;
      const double *fe = gkyl_array_cfetch(vflux, gkyl_range_idx(&s.plocal_ext, idx));
      for (int k = 0; k < 2 * nc; ++k) edge_max = fmax(edge_max, fabs(fe[k]));
    }
  }
  TEST_CHECK(edge_max == 0.0);
  gkyl_array_release(vflux);
  TEST_CHECK(num_null == 0);

  // Weak residuals in every mode: M0 c_0 + r0, M0 alpha_0 + M1 c + r1 + c_0 M1 and
  // 2 M1 alpha_0 + 2 (M2 + 2 M0G) c + r2 + 2 r0G + c_0 M2.
  struct gkyl_array *n_c0 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *n_aP = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c0_m1 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *ch_m1 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *m1_aP = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *s2 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *ch_s2 = mkarr(nc, s.clocal_ext.volume);
  struct gkyl_array *c0_m2 = mkarr(nc, s.clocal_ext.volume);
  gkyl_array_clear(n_c0, 0.0);
  gkyl_array_clear(n_aP, 0.0);
  gkyl_array_clear(c0_m1, 0.0);
  gkyl_array_clear(ch_m1, 0.0);
  gkyl_array_clear(m1_aP, 0.0);
  gkyl_array_clear(ch_s2, 0.0);
  gkyl_array_clear(c0_m2, 0.0);
  gkyl_array_set_offset(s2, 1.0, moms, nc);
  gkyl_array_accumulate_offset(s2, 2.0, moms, 2 * nc);
  gkyl_dg_mul_op_range(&s.cbasis, 0, n_c0, 0, moms, 0, c_0, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, n_aP, 0, moms, 0, a_P, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, c0_m1, 0, c_0, 3, moms, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, ch_m1, 0, c_h, 3, moms, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, m1_aP, 0, a_P, 3, moms, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, ch_s2, 0, c_h, 0, s2, &s.clocal);
  gkyl_dg_mul_op_range(&s.cbasis, 0, c0_m2, 0, c_0, 1, moms, &s.clocal);
  double max_res = 0.0, max_m1 = 0.0, max_aP = 0.0, max_res0 = 0.0, max_m0 = 0.0;
  double max_res2 = 0.0, max_m2 = 0.0, max_ch = 0.0;
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &s.clocal);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(&s.clocal_ext, citer.idx);
    const double *na = gkyl_array_cfetch(n_aP, cidx);
    const double *nc0 = gkyl_array_cfetch(n_c0, cidx);
    const double *cm1 = gkyl_array_cfetch(c0_m1, cidx), *chm1 = gkyl_array_cfetch(ch_m1, cidx);
    const double *m1a = gkyl_array_cfetch(m1_aP, cidx), *chs = gkyl_array_cfetch(ch_s2, cidx);
    const double *cm2 = gkyl_array_cfetch(c0_m2, cidx);
    const double *rm = gkyl_array_cfetch(rhs_moms, cidx);
    const double *ap = gkyl_array_cfetch(a_P, cidx), *ch = gkyl_array_cfetch(c_h, cidx);
    for (int k = 0; k < nc; ++k) {
      double res0 = fabs(nc0[k] + rm[k]);
      max_res0 = res0 > max_res0 ? res0 : max_res0;
      max_m0 = fabs(rm[k]) > max_m0 ? fabs(rm[k]) : max_m0;
      double res = fabs(na[k] + chm1[k] + rm[3 * nc + k] + cm1[k]);
      max_res = res > max_res ? res : max_res;
      max_m1 = fabs(rm[3 * nc + k]) > max_m1 ? fabs(rm[3 * nc + k]) : max_m1;
      max_aP = fabs(ap[k]) > max_aP ? fabs(ap[k]) : max_aP;
      double r2 = rm[nc + k] + 2.0 * rm[2 * nc + k];
      double res2 = fabs(2.0 * m1a[k] + 2.0 * chs[k] + r2 + cm2[k]);
      max_res2 = res2 > max_res2 ? res2 : max_res2;
      max_m2 = fabs(r2) > max_m2 ? fabs(r2) : max_m2;
      max_ch = fabs(ch[k]) > max_ch ? fabs(ch[k]) : max_ch;
    }
  }
  TEST_CHECK(max_m0 > 0.0);
  TEST_CHECK(max_res0 < 1e-12 * max_m0);
  TEST_MSG("cdim=%d: max normalization residual %g (max M0(rhs) %g)", cdim, max_res0, max_m0);
  gkyl_array_release(n_c0);
  gkyl_array_release(n_aP);
  gkyl_array_release(c0_m1);
  gkyl_array_release(ch_m1);
  gkyl_array_release(m1_aP);
  gkyl_array_release(s2);
  gkyl_array_release(ch_s2);
  gkyl_array_release(c0_m2);
  gkyl_array_release(c_0);
  TEST_CHECK(max_m1 > 0.0);
  TEST_CHECK(max_aP > 0.0);
  TEST_CHECK(max_res < 1e-12 * max_m1);
  TEST_MSG("cdim=%d: first-moment residual %g (max |M1(rhs)| %g, max |alpha_0| %g)", cdim, max_res, max_m1, max_aP);
  TEST_CHECK(max_m2 > 0.0);
  TEST_CHECK(max_ch > 0.0);
  TEST_CHECK(max_res2 < 1e-12 * max_m2);
  TEST_MSG("cdim=%d: second-moment residual %g (max |M2 + 2 M0G (rhs)| %g, max |c| %g)", cdim, max_res2, max_m2, max_ch);

  gkyl_array_release(a_P);
  gkyl_array_release(c_h);
  gkyl_dg_calc_pkpm_vars_release(vars);
  gkyl_array_release(rhs_moms);
  gkyl_array_release(moms);
  gkyl_dg_updater_moment_pkpm_release(mcalc);
  gkyl_array_release(rhs);
  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_array_release(cflrate);
  gkyl_array_release(flux);
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
void test_closure_1x() { test_closure(1); }
void test_closure_2x() { test_closure(2); }
void test_closure_3x() { test_closure(3); }

TEST_LIST = {
  {"test_grad_u_1x", test_grad_u_1x},
  {"test_grad_u_2x", test_grad_u_2x},
  {"test_grad_u_3x", test_grad_u_3x},
  {"test_vol_1x", test_vol_1x},
  {"test_vol_2x", test_vol_2x},
  {"test_vol_3x", test_vol_3x},
  {"test_closure_1x", test_closure_1x},
  {"test_closure_2x", test_closure_2x},
  {"test_closure_3x", test_closure_3x},
  {NULL, NULL},
};
