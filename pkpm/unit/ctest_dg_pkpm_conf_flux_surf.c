// Tests of the PKPM configuration-space HLL surface flux (dg_pkpm_conf_flux_surf)
// against an independent nodal evaluation, and of the stored-flux scatter of the
// configuration pass (dg_pkpm_kinetic_conf) through hyper_dg: conservation on a
// periodic domain.
#include <math.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_kinetic_conf.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_hyper_dg.h>
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

// Periodic smooth fields on [-pi, pi]^cdim so that projection on the extended
// range fills the ghost cells with their periodic values.
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
  // Comparable to the thermal speed so that a_l and a_r change sign across nodes.
  fout[0] = 0.8 * sin(x) + 0.3 * cos(y);
  fout[1] = -0.7 * cos(y) * sin(z) + 0.2 * sin(x);
  fout[2] = 0.5 * sin(z) - 0.4 * cos(x) * cos(y);
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
  fout[0] = F0 * (1.0 + 0.1 * v); // F_0 (a nonzero first moment is fine here)
  fout[1] = 0.5 * vt * vt * F0 * (1.0 - 0.1 * v);
}

static const double gl3_x[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double gl3_w[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

// Reference HLL flux, long form, for the test oracle.
static double
hll_ref(double a_l, double a_r, double f_l, double f_r)
{
  double s_l = fmin(fmin(a_l, a_r), 0.0), s_r = fmax(fmax(a_l, a_r), 0.0);
  if (s_l == 0.0 && s_r == 0.0) {
    return 0.0;
  }
  if (s_l >= 0.0) {
    return a_l * f_l;
  }
  if (s_r <= 0.0) {
    return a_r * f_r;
  }
  return (s_r * a_l * f_l - s_l * a_r * f_r + s_l * s_r * (f_r - f_l)) / (s_r - s_l);
}

struct test_setup {
  struct gkyl_vlasov_position_map *pos_map; // identity configuration-space map
  int cdim, pdim;
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext;
  struct gkyl_array *bvar, *u, *xi, *f, *cflrate, *flux;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
  struct gkyl_dg_pkpm_conf_flux_surf *calc_flux;
};

static void
setup(struct test_setup *s, int cdim)
{
  int pdim = cdim + 1;
  s->cdim = cdim;
  s->pdim = pdim;
  double lower[4] = {-M_PI, -M_PI, -M_PI, -4.0}, upper[4] = {M_PI, M_PI, M_PI, 4.0};
  int cells[4] = {6, 4, 3, 8};
  int ccells[3] = {6, 4, 3};
  int cghost[3] = {1, 1, 1}, pghost[4] = {1, 1, 1, 0};
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
  s->pos_map = identity_pos_map(&s->cgrid, &s->clocal, &s->clocal_ext, &s->cbasis);
  gkyl_create_grid_ranges(&s->pgrid, pghost, &s->plocal_ext, &s->plocal);

  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis;
  struct gkyl_array *em = mkarr(8 * nc, s->clocal_ext.volume);
  s->bvar = mkarr(3 * nc, s->clocal_ext.volume);
  s->u = mkarr(3 * nc, s->clocal_ext.volume);
  // thermal speed [xi, xi^2, p] = [1, 1, 0]: the streaming speed is v b + u as in the reference
  s->xi = mkarr(3 * nc, s->clocal_ext.volume);
  gkyl_array_clear(s->xi, 0.0);
  gkyl_array_shiftc(s->xi, pow(sqrt(2.0), cdim), 0 * nc);
  gkyl_array_shiftc(s->xi, pow(sqrt(2.0), cdim), 1 * nc);
  s->f = mkarr(2 * np, s->plocal_ext.volume);
  s->cflrate = mkarr(1, s->plocal_ext.volume);

  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 8, eval_em, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, em);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 3, eval_u, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, s->u);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s->pgrid, &s->pbasis, 8, 2, eval_f, &cdim);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->plocal_ext, s->f);
  gkyl_proj_on_basis_release(proj);

  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->cgrid, &s->cbasis, s->pos_map, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->clocal_ext, em, s->bvar);
  gkyl_array_release(em);

  struct gkyl_dg_pkpm_conf_flux_surf_inp inp = {
    .phase_grid = &s->pgrid,
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .pos_map = s->pos_map, .use_gpu = false,
  };
  s->calc_flux = gkyl_dg_pkpm_conf_flux_surf_inew(&inp);
  int nk = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(s->calc_flux);
  TEST_CHECK(nk == np / 2);
  s->flux = mkarr(2 * cdim * nk, s->plocal_ext.volume);
  gkyl_array_clear(s->flux, 0.0);
  gkyl_array_clear(s->cflrate, 0.0);
  gkyl_dg_pkpm_conf_flux_surf_advance(
    s->calc_flux, &s->clocal, &s->plocal, &s->plocal_ext, s->bvar, s->u, s->xi, s->f, s->cflrate, s->flux
  );
}

static void
teardown(struct test_setup *s)
{
  gkyl_vlasov_position_map_release(s->pos_map);
  gkyl_dg_pkpm_conf_flux_surf_release(s->calc_flux);
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_array_release(s->flux);
  gkyl_array_release(s->cflrate);
  gkyl_array_release(s->f);
  gkyl_array_release(s->u);
  gkyl_array_release(s->xi);
  gkyl_array_release(s->bvar);
}

// Independent reference for the stored flux on the face between cells idx_l and
// idx_r in direction dir (the lower face of idx_r). Returns the maximum speed.
static double
reference_face(
  struct test_setup *s, int dir, const int *idx_l, const int *idx_r, const struct gkyl_basis *sbasis,
  double *flux_ref
)
{
  int cdim = s->cdim, pdim = s->pdim;
  int nc = s->cbasis.num_basis, np = s->pbasis.num_basis, nk = sbasis->num_basis;
  long cl = gkyl_range_idx(&s->clocal_ext, idx_l), cr = gkyl_range_idx(&s->clocal_ext, idx_r);
  long pl = gkyl_range_idx(&s->plocal_ext, idx_l), pr = gkyl_range_idx(&s->plocal_ext, idx_r);
  const double *b_l = gkyl_array_cfetch(s->bvar, cl), *b_r = gkyl_array_cfetch(s->bvar, cr);
  const double *u_l = gkyl_array_cfetch(s->u, cl), *u_r = gkyl_array_cfetch(s->u, cr);
  const double *f_l = gkyl_array_cfetch(s->f, pl), *f_r = gkyl_array_cfetch(s->f, pr);
  double xc[4];
  gkyl_rect_grid_cell_center(&s->pgrid, idx_r, xc);
  double wv = xc[cdim], dv = s->pgrid.dx[cdim];

  for (int k = 0; k < 2 * nk; ++k) {
    flux_ref[k] = 0.0;
  }
  double amax = 0.0;
  // Transverse configuration directions, then v_par: 3 nodes each.
  int ntr = cdim - 1;
  int nnodes = 1;
  for (int d = 0; d < ntr + 1; ++d) {
    nnodes *= 3;
  }
  for (int n = 0; n < nnodes; ++n) {
    int sub[4];
    int m = n;
    for (int d = ntr; d >= 0; --d) { // last surface variable fastest
      sub[d] = m % 3;
      m /= 3;
    }
    double w = 1.0, zs[4], zl[4], zr[4];
    int t = 0;
    for (int d = 0; d < cdim; ++d) {
      if (d == dir) {
        zl[d] = 1.0;
        zr[d] = -1.0;
      } else {
        zl[d] = zr[d] = zs[t] = gl3_x[sub[t]];
        w *= gl3_w[sub[t]];
        t += 1;
      }
    }
    zl[cdim] = zr[cdim] = zs[t] = gl3_x[sub[t]];
    w *= gl3_w[sub[t]];
    double v = wv + 0.5 * dv * zs[t];

    double a_l = v * s->cbasis.eval_expand(zl, &b_l[dir * nc]) + s->cbasis.eval_expand(zl, &u_l[dir * nc]);
    double a_r = v * s->cbasis.eval_expand(zr, &b_r[dir * nc]) + s->cbasis.eval_expand(zr, &u_r[dir * nc]);
    amax = fmax(amax, fmax(fabs(a_l), fabs(a_r)));
    for (int c = 0; c < 2; ++c) {
      double fl = s->pbasis.eval_expand(zl, &f_l[c * np]);
      double fr = s->pbasis.eval_expand(zr, &f_r[c * np]);
      double ghat = hll_ref(a_l, a_r, fl, fr);
      for (int k = 0; k < nk; ++k) {
        // Surface basis function k at the surface node, divided by its value at the
        // emitter's shape-normalization point (2/5 in every variable): the stored
        // flux is the shape coefficient, and the scatter maps fold the constant back in.
        double unit[64] = {0.0};
        unit[k] = 1.0;
        double znorm[4] = {0.4, 0.4, 0.4, 0.4};
        double chi_k = sbasis->eval_expand(zs, unit) / sbasis->eval_expand(znorm, unit);
        flux_ref[c * nk + k] += w * chi_k * ghat;
      }
    }
  }
  return amax;
}

static void
test_conf_flux(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nk = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(s.calc_flux);

  struct gkyl_basis sbasis;
  if (cdim == 1) {
    gkyl_cart_modal_serendip(&sbasis, 1, 2);
  } else {
    gkyl_cart_modal_hybrid(&sbasis, cdim - 1, 1);
  }
  TEST_CHECK(sbasis.num_basis == nk);

  double flux_ref[64];
  double max_flux = 0.0, max_err = 0.0, max_cfl_err = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    const double *flux = gkyl_array_cfetch(s.flux, pidx);
    double cfl_ref = 0.0;
    for (int dir = 0; dir < cdim; ++dir) {
      int idx_l[4], idx_r[4];
      gkyl_copy_int_arr(s.pdim, iter.idx, idx_l);
      idx_l[dir] -= 1;
      double amax = reference_face(&s, dir, idx_l, iter.idx, &sbasis, flux_ref);
      cfl_ref += 1.5 * 2.0 / s.pgrid.dx[dir] * amax;
      for (int k = 0; k < 2 * nk; ++k) {
        double err = fabs(flux[dir * 2 * nk + k] - flux_ref[k]);
        max_err = err > max_err ? err : max_err;
        max_flux = fabs(flux_ref[k]) > max_flux ? fabs(flux_ref[k]) : max_flux;
      }
      // Upper ghost-owned face.
      if (iter.idx[dir] == s.plocal.upper[dir]) {
        gkyl_copy_int_arr(s.pdim, iter.idx, idx_r);
        idx_r[dir] += 1;
        long pidx_r = gkyl_range_idx(&s.plocal_ext, idx_r);
        const double *flux_g = gkyl_array_cfetch(s.flux, pidx_r);
        amax = reference_face(&s, dir, iter.idx, idx_r, &sbasis, flux_ref);
        cfl_ref += 1.5 * 2.0 / s.pgrid.dx[dir] * amax;
        for (int k = 0; k < 2 * nk; ++k) {
          double err = fabs(flux_g[dir * 2 * nk + k] - flux_ref[k]);
          max_err = err > max_err ? err : max_err;
        }
      }
    }
    const double *cfl = gkyl_array_cfetch(s.cflrate, pidx);
    double cerr = fabs(cfl[0] - cfl_ref);
    max_cfl_err = cerr > max_cfl_err ? cerr : max_cfl_err;
  }
  TEST_CHECK(max_flux > 0.0);
  TEST_CHECK(max_err < 1e-12 * max_flux);
  TEST_MSG("cdim=%d: max flux error %g (max flux %g)", cdim, max_err, max_flux);
  TEST_CHECK(max_cfl_err < 1e-11);
  TEST_MSG("cdim=%d: max cflrate error %g", cdim, max_cfl_err);

  teardown(&s);
}

// hyper_dg over the configuration directions with the scatter surface term and no
// volume term: the increments telescope, so the cell-averaged RHS sums to zero on
// the periodic domain for both components, while individual cells are nonzero.
static void
test_scatter(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int nc = s.cbasis.num_basis, np = s.pbasis.num_basis;

  // The volume term needs theta and the collision coefficients; zeros here (the
  // conservation check below holds for any values).
  struct gkyl_array *zeros = mkarr(nc, s.clocal_ext.volume);
  gkyl_array_clear(zeros, 0.0);
  struct gkyl_array *zeros_g = mkarr(3 * s.cdim * nc, s.clocal_ext.volume);
  gkyl_array_clear(zeros_g, 0.0);
  struct gkyl_dg_pkpm_kinetic_conf_inp inp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields =
      {
        .bvar = s.bvar,
        .grad_b = zeros_g,
        .u = s.u,
        .grad_u = zeros_g,
        .xi = s.xi,
        .grad_xi = zeros_g,
        .grad_ln_n = zeros_g,
        .theta = zeros,
        .conf_flux_surf = s.flux,
      },
    .pos_map = s.pos_map, .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn = gkyl_dg_pkpm_kinetic_conf_inew(&inp);
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *slvr = gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn, cdim, up_dirs, zero_flux, 1, false);

  struct gkyl_array *rhs = mkarr(2 * np, s.plocal_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(s.cflrate, 0.0);
  gkyl_hyper_dg_advance(slvr, &s.plocal, s.f, s.cflrate, rhs);

  // F_0 has no volume sources (its volume terms are divergences), so its cell averages
  // change through the scattered face fluxes only, which telescope on the periodic
  // grid. G carries the (K - div u) G source and is not conserved.
  double sum = 0.0, mag = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    const double *r = gkyl_array_cfetch(rhs, gkyl_range_idx(&s.plocal_ext, iter.idx));
    sum += r[0];
    mag += fabs(r[0]);
  }
  TEST_CHECK(mag > 0.0);
  TEST_CHECK(fabs(sum) < 1e-12 * mag);
  TEST_MSG("cdim=%d: F_0 sum %g, |sum| %g", cdim, sum, mag);

  gkyl_array_release(rhs);
  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_array_release(zeros_g);
  gkyl_array_release(zeros);
  teardown(&s);
}

// The full stage-A right-hand side (hyper_dg: volume kernel + scattered face fluxes)
// against the volume kernel plus an independent DG surface term built from the stored
// fluxes of every configuration direction:
//   dx1_d [ <psi_m(xi_d = -1), Fhat_l> - <psi_m(xi_d = +1), Fhat_r> ]
// with Fhat = sum_k c_k chi_k over the face variables (the other configuration
// variables and v_par), c_k = g_k chi_k(2/5, ...) (shape to orthonormal coefficient).
static void
test_scatter_ref(int cdim)
{
  struct test_setup s;
  setup(&s, cdim);
  int pdim = s.pdim, nc = s.cbasis.num_basis, np = s.pbasis.num_basis, nk = np / 2;
  struct gkyl_basis sbasis;
  if (cdim == 1) {
    gkyl_cart_modal_serendip(&sbasis, 1, 2);
  } else {
    gkyl_cart_modal_hybrid(&sbasis, cdim - 1, 1);
  }
  TEST_CHECK(sbasis.num_basis == nk);

  struct gkyl_array *zeros = mkarr(nc, s.clocal_ext.volume);
  gkyl_array_clear(zeros, 0.0);
  struct gkyl_array *zeros_g = mkarr(3 * cdim * nc, s.clocal_ext.volume);
  gkyl_array_clear(zeros_g, 0.0);
  struct gkyl_dg_pkpm_kinetic_conf_inp inp = {
    .conf_basis = &s.cbasis,
    .phase_basis = &s.pbasis,
    .conf_range = &s.clocal,
    .phase_range = &s.plocal,
    .auxfields =
      {.bvar = s.bvar, .grad_b = zeros_g, .u = s.u, .grad_u = zeros_g, .xi = s.xi, .grad_xi = zeros_g, .grad_ln_n = zeros_g, .theta = zeros,
        .conf_flux_surf = s.flux},
    .pos_map = s.pos_map, .use_gpu = false,
  };
  struct gkyl_dg_eqn *eqn = gkyl_dg_pkpm_kinetic_conf_inew(&inp);
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2}, zero_flux[2 * GKYL_MAX_DIM] = {0};
  gkyl_hyper_dg *slvr = gkyl_hyper_dg_new(&s.pgrid, &s.pbasis, eqn, cdim, up_dirs, zero_flux, 1, false);
  struct gkyl_array *rhs = mkarr(2 * np, s.plocal_ext.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(s.cflrate, 0.0);
  gkyl_hyper_dg_advance(slvr, &s.plocal, s.f, s.cflrate, rhs);

  // Volume kernels of the three dimensionalities.
  double max_err = 0.0, max_rhs = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.plocal);
  while (gkyl_range_iter_next(&iter)) {
    long pidx = gkyl_range_idx(&s.plocal_ext, iter.idx);
    long cidx = gkyl_range_idx(&s.clocal_ext, iter.idx);
    double xc[4];
    gkyl_rect_grid_cell_center(&s.pgrid, iter.idx, xc);
    double ref[64] = {0.0};
    const double *bv = gkyl_array_cfetch(s.bvar, cidx), *uv = gkyl_array_cfetch(s.u, cidx);
    const double *zg = gkyl_array_cfetch(zeros_g, cidx), *z = gkyl_array_cfetch(zeros, cidx);
    const double *fv = gkyl_array_cfetch(s.f, pidx), *xv = gkyl_array_cfetch(s.xi, cidx);
    const double *jac = gkyl_array_cfetch(s.pos_map->jacob_pos, cidx); // identity map: ones
    if (cdim == 1) {
      pkpm_kinetic_conf_vol_1x1v_ser_p1(xc, s.pgrid.dx, jac, bv, zg, uv, zg, xv, zg, zg, z, fv, ref);
    } else if (cdim == 2) {
      pkpm_kinetic_conf_vol_2x1v_ser_p1(xc, s.pgrid.dx, jac, bv, zg, uv, zg, xv, zg, zg, z, fv, ref);
    } else {
      pkpm_kinetic_conf_vol_3x1v_ser_p1(xc, s.pgrid.dx, jac, bv, zg, uv, zg, xv, zg, zg, z, fv, ref);
    }
    for (int dir = 0; dir < cdim; ++dir) {
      int idx_r[4];
      gkyl_copy_int_arr(pdim, iter.idx, idx_r);
      idx_r[dir] += 1;
      long pidx_r = gkyl_range_idx(&s.plocal_ext, idx_r);
      const double *g_l = gkyl_array_cfetch(s.flux, pidx), *g_r = gkyl_array_cfetch(s.flux, pidx_r);
      double dx1 = 2.0 / s.pgrid.dx[dir];
      int nnodes = 1;
      for (int d = 0; d < cdim; ++d) {
        nnodes *= 3;
      }
      for (int c = 0; c < 2; ++c) {
        for (int m = 0; m < np; ++m) {
          double unit_m[64] = {0.0};
          unit_m[m] = 1.0;
          double surf = 0.0;
          for (int n = 0; n < nnodes; ++n) {
            // Surface nodes: the other configuration variables then v_par, 3 each.
            int sub[4], q = n;
            for (int d = cdim - 1; d >= 0; --d) {
              sub[d] = q % 3;
              q /= 3;
            }
            double w = 1.0, zs[4], zl[4], zr[4];
            int t = 0;
            for (int d = 0; d < cdim; ++d) {
              if (d == dir) {
                zl[d] = 1.0;
                zr[d] = -1.0;
              } else {
                zl[d] = zr[d] = zs[t] = gl3_x[sub[t]];
                w *= gl3_w[sub[t]];
                t += 1;
              }
            }
            zl[cdim] = zr[cdim] = zs[t] = gl3_x[sub[t]];
            w *= gl3_w[sub[t]];
            double fl = 0.0, fr = 0.0;
            for (int k = 0; k < nk; ++k) {
              double unit_k[64] = {0.0};
              unit_k[k] = 1.0;
              double znorm[4] = {0.4, 0.4, 0.4, 0.4};
              double chi_norm = sbasis.eval_expand(znorm, unit_k), chi = sbasis.eval_expand(zs, unit_k);
              fl += g_l[dir * 2 * nk + c * nk + k] * chi_norm * chi;
              fr += g_r[dir * 2 * nk + c * nk + k] * chi_norm * chi;
            }
            // psi_m at the lower face of this cell (lower flux enters with +) and at its
            // upper face (upper flux enters with -).
            double zlo[4], zhi[4];
            for (int d = 0; d < pdim; ++d) {
              zlo[d] = zr[d];
              zhi[d] = zl[d];
            }
            surf += w * (s.pbasis.eval_expand(zlo, unit_m) * fl - s.pbasis.eval_expand(zhi, unit_m) * fr);
          }
          ref[c * np + m] += dx1 * surf;
        }
      }
    }
    const double *r = gkyl_array_cfetch(rhs, pidx);
    for (int k = 0; k < 2 * np; ++k) {
      max_err = fmax(max_err, fabs(r[k] - ref[k]));
      max_rhs = fmax(max_rhs, fabs(ref[k]));
    }
  }
  TEST_CHECK(max_rhs > 0.0);
  TEST_CHECK(max_err < 1e-11 * max_rhs);
  TEST_MSG("cdim=%d: stage-A rhs vs reference: max error %g (max rhs %g)", cdim, max_err, max_rhs);

  gkyl_array_release(rhs);
  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
  gkyl_array_release(zeros_g);
  gkyl_array_release(zeros);
  teardown(&s);
}


void test_scatter_ref_1x() { test_scatter_ref(1); }
void test_scatter_ref_2x() { test_scatter_ref(2); }
void test_scatter_ref_3x() { test_scatter_ref(3); }

void test_conf_flux_1x() { test_conf_flux(1); }
void test_conf_flux_2x() { test_conf_flux(2); }
void test_conf_flux_3x() { test_conf_flux(3); }
void test_scatter_1x() { test_scatter(1); }
void test_scatter_2x() { test_scatter(2); }
void test_scatter_3x() { test_scatter(3); }

TEST_LIST = {
  {"test_conf_flux_1x", test_conf_flux_1x},
  {"test_conf_flux_2x", test_conf_flux_2x},
  {"test_conf_flux_3x", test_conf_flux_3x},
  {"test_scatter_ref_1x", test_scatter_ref_1x},
  {"test_scatter_ref_2x", test_scatter_ref_2x},
  {"test_scatter_ref_3x", test_scatter_ref_3x},
  {"test_scatter_1x", test_scatter_1x},
  {"test_scatter_2x", test_scatter_2x},
  {"test_scatter_3x", test_scatter_3x},
  {NULL, NULL},
};
