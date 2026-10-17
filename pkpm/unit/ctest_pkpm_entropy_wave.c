// Semi-discrete stability check of the coupled refactored PKPM update (kinetic two-pass
// update plus the momentum solve) on a 1x1v periodic entropy wave: a density pulse
// advected by a uniform flow with uniform pressure. The exact solution keeps u and p
// uniform; the test integrates the coupled zero-layer operators with SSP-RK3 and checks
// that the flow stays bounded. (The legacy split-form momentum flux with a Lax mass
// flux was unstable here at a rate ~ u/dx; the nodal Roe flux is not.) The kinetic and
// fluid halves can be frozen separately to localize an instability (diagnostics,
// printed with -v).
#include <math.h>
#include <stdio.h>

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_pkpm_bvar.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_fluid.h>
#include <gkyl_dg_pkpm_fluid_flux_surf.h>
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

struct wave_ctx {
  double u0, alpha, vt;
  double b_theta, b_perp; // static field: B = (sin theta, -cos theta, -b_perp cos(pi x))
  double uy_amp; // sheared perpendicular drift u_y = uy_amp cos(pi x)
  double ux_amp; // compressive flow u_x = u0 + ux_amp cos(pi x)
};

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

static void
eval_em(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct wave_ctx *c = ctx;
  for (int i = 0; i < 8; ++i) fout[i] = 0.0;
  if (c->b_perp == 0.0 && c->b_theta == 0.0) {
    fout[3] = 1.0; // B = x-hat
    return;
  }
  fout[3] = sin(c->b_theta);
  fout[4] = -cos(c->b_theta);
  fout[5] = -c->b_perp * cos(M_PI * xn[0]);
}

static void
eval_f(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct wave_ctx *c = ctx;
  double x = xn[0], v = xn[1];
  double n = 1.0 + c->alpha * sin(M_PI * x);
  double T = c->vt * c->vt / n; // uniform pressure p = n T
  double F0 = n / sqrt(2.0 * M_PI * T) * exp(-v * v / (2.0 * T));
  fout[0] = F0;
  fout[1] = T * F0;
}

static void
eval_euler(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct wave_ctx *c = ctx;
  double n = 1.0 + c->alpha * sin(M_PI * xn[0]);
  fout[0] = n * (c->u0 + c->ux_amp * cos(M_PI * xn[0]));
  fout[1] = n * c->uy_amp * cos(M_PI * xn[0]);
  fout[2] = 0.0;
}

struct run_opts {
  double u0, alpha, t_end;
  int Nx, Nv;
  bool freeze_fluid, freeze_kinetic;
  double b_theta, b_perp; // oblique, x-dependent static field (0, 0: b = x-hat)
  double uy_amp; // sheared perpendicular drift amplitude
  double ux_amp; // compressive flow amplitude
};

struct sim {
  struct gkyl_rect_grid cgrid, pgrid;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_range clocal, clocal_ext, plocal, plocal_ext, vlocal, vlocal_ext;
  struct gkyl_range cskin[2], cghost[2], pskin[2], pghost[2];
  int nc, np, nkc, nkf, nmf;
  struct gkyl_array *bvar, *grad_b, *u, *grad_u, *theta, *nu_sum, *nu_vthsq_sum, *a_P, *moms, *rhs_moms;
  struct gkyl_array *cflux, *vflux, *fflux, *cflrate_f, *cflrate_m, *cell_avg;
  struct gkyl_array *roe_mf, *kin_mf, *kin_mf_face;
  struct gkyl_dg_calc_pkpm_bvar *calc_b;
  struct gkyl_dg_calc_pkpm_vars *vars;
  struct gkyl_wv_eqn *ten;
  struct gkyl_wave_geom *geom;
  struct gkyl_dg_updater_moment *mcalc;
  struct gkyl_dg_pkpm_conf_flux_surf *calc_cflux;
  struct gkyl_dg_pkpm_vel_flux_surf *calc_vflux;
  struct gkyl_dg_pkpm_fluid_flux_surf *calc_fflux;
  struct gkyl_dg_eqn *eqn_a, *eqn_b, *eqn_m;
  gkyl_hyper_dg *slvr_a, *slvr_b, *slvr_m;
};

static void
periodic(struct gkyl_array *arr, const struct gkyl_range *skin, const struct gkyl_range *ghost)
{
  // lower ghost <- upper skin, upper ghost <- lower skin
  gkyl_array_copy_range_to_range(arr, arr, &ghost[0], &skin[1]);
  gkyl_array_copy_range_to_range(arr, arr, &ghost[1], &skin[0]);
}

static void
sim_init(struct sim *s, const struct run_opts *o)
{
  struct wave_ctx bctx = {.b_theta = o->b_theta, .b_perp = o->b_perp};
  double lower[2] = {-1.0, -6.0}, upper[2] = {1.0, 6.0};
  int cells[2] = {o->Nx, o->Nv};
  int cghost[2] = {1, 0}, vghost[1] = {0};
  gkyl_rect_grid_init(&s->cgrid, 1, lower, upper, cells);
  gkyl_rect_grid_init(&s->pgrid, 2, lower, upper, cells);
  gkyl_cart_modal_serendip(&s->cbasis, 1, 1);
  gkyl_cart_modal_hybrid(&s->pbasis, 1, 1);
  gkyl_create_grid_ranges(&s->cgrid, cghost, &s->clocal_ext, &s->clocal);
  gkyl_create_grid_ranges(&s->pgrid, cghost, &s->plocal_ext, &s->plocal);
  struct gkyl_rect_grid vgrid;
  gkyl_rect_grid_init(&vgrid, 1, &lower[1], &upper[1], &cells[1]);
  gkyl_create_grid_ranges(&vgrid, vghost, &s->vlocal_ext, &s->vlocal);
  gkyl_skin_ghost_ranges(&s->cskin[0], &s->cghost[0], 0, GKYL_LOWER_EDGE, &s->clocal_ext, cghost);
  gkyl_skin_ghost_ranges(&s->cskin[1], &s->cghost[1], 0, GKYL_UPPER_EDGE, &s->clocal_ext, cghost);
  gkyl_skin_ghost_ranges(&s->pskin[0], &s->pghost[0], 0, GKYL_LOWER_EDGE, &s->plocal_ext, cghost);
  gkyl_skin_ghost_ranges(&s->pskin[1], &s->pghost[1], 0, GKYL_UPPER_EDGE, &s->plocal_ext, cghost);

  int nc = s->nc = s->cbasis.num_basis, np = s->np = s->pbasis.num_basis;
  long cvol = s->clocal_ext.volume, pvol = s->plocal_ext.volume;
  s->bvar = mkarr(3 * nc, cvol);
  s->u = mkarr(3 * nc, cvol);
  s->grad_u = mkarr(3 * nc, cvol);
  gkyl_array_clear(s->grad_u, 0.0);
  s->theta = mkarr(nc, cvol);
  s->nu_sum = mkarr(nc, cvol);
  s->nu_vthsq_sum = mkarr(nc, cvol);
  s->a_P = mkarr(nc, cvol);
  s->moms = mkarr(4 * nc, cvol);
  s->rhs_moms = mkarr(4 * nc, cvol);
  s->cflrate_f = mkarr(1, pvol);
  s->cflrate_m = mkarr(1, cvol);
  s->cell_avg = gkyl_array_new(GKYL_INT, 1, cvol);
  gkyl_array_clear(s->nu_sum, 0.0);
  gkyl_array_clear(s->nu_vthsq_sum, 0.0);
  gkyl_array_clear(s->a_P, 0.0);

  struct gkyl_array *em = mkarr(8 * nc, cvol);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s->cgrid, &s->cbasis, 8, 8, eval_em, &bctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->clocal_ext, em);
  gkyl_proj_on_basis_release(proj);
  s->calc_b = gkyl_dg_calc_pkpm_bvar_new(&s->cgrid, &s->cbasis, false);
  gkyl_dg_calc_pkpm_bvar_advance(s->calc_b, &s->clocal_ext, em, s->bvar);
  s->grad_b = mkarr(3 * nc, cvol);
  gkyl_array_clear(s->grad_b, 0.0);
  gkyl_dg_calc_pkpm_bvar_grad(s->calc_b, &s->clocal, s->bvar, s->grad_b);
  gkyl_array_release(em);

  s->ten = gkyl_wv_ten_moment_new(0.0, false, false, 0, 0, false);
  s->geom = gkyl_wave_geom_new(&s->cgrid, &s->clocal_ext, 0, 0, false);
  s->vars =
    gkyl_dg_calc_pkpm_vars_new(&s->cgrid, &s->cbasis, &s->clocal_ext, s->ten, s->geom, 0.0, false);
  s->mcalc = gkyl_dg_updater_moment_pkpm_new(
    &s->pgrid, &s->cbasis, &s->pbasis, &s->clocal_ext, &s->vlocal, 1.0, false, false
  );

  struct gkyl_dg_pkpm_conf_flux_surf_inp cinp = {
    .phase_grid = &s->pgrid, .conf_basis = &s->cbasis, .phase_basis = &s->pbasis, .use_gpu = false};
  s->calc_cflux = gkyl_dg_pkpm_conf_flux_surf_inew(&cinp);
  s->nkc = gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(s->calc_cflux);
  s->cflux = mkarr(2 * s->nkc, pvol);
  s->nmf = gkyl_dg_pkpm_conf_flux_surf_num_face_basis(s->calc_cflux);
  s->kin_mf = mkarr(s->nmf, pvol);
  s->kin_mf_face = mkarr(s->nmf, cvol);
  s->roe_mf = mkarr(s->nmf, cvol);
  gkyl_array_clear(s->roe_mf, 0.0);
  struct gkyl_dg_pkpm_kinetic_conf_inp ainp = {
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
        .conf_flux_surf = s->cflux,
      },
    .use_gpu = false,
  };
  s->eqn_a = gkyl_dg_pkpm_kinetic_conf_inew(&ainp);
  int up_dirs_a[GKYL_MAX_DIM] = {0}, zero_flux_a[2 * GKYL_MAX_DIM] = {0};
  s->slvr_a = gkyl_hyper_dg_new(&s->pgrid, &s->pbasis, s->eqn_a, 1, up_dirs_a, zero_flux_a, 1, false);

  struct gkyl_dg_pkpm_vel_flux_surf_inp vinp = {
    .phase_grid = &s->pgrid, .conf_basis = &s->cbasis, .phase_basis = &s->pbasis, .use_gpu = false};
  s->calc_vflux = gkyl_dg_pkpm_vel_flux_surf_inew(&vinp);
  s->vflux = mkarr(2 * nc, pvol);
  struct gkyl_dg_pkpm_kinetic_vel_inp binp = {
    .conf_basis = &s->cbasis,
    .phase_basis = &s->pbasis,
    .conf_range = &s->clocal,
    .phase_range = &s->plocal,
    .auxfields = {.a_P = s->a_P, .vel_flux_surf = s->vflux},
    .use_gpu = false,
  };
  s->eqn_b = gkyl_dg_pkpm_kinetic_vel_inew(&binp);
  int up_dirs_b[GKYL_MAX_DIM] = {1}, zero_flux_b[2 * GKYL_MAX_DIM] = {0};
  zero_flux_b[1] = zero_flux_b[1 + 2] = 1;
  s->slvr_b = gkyl_hyper_dg_new(&s->pgrid, &s->pbasis, s->eqn_b, 1, up_dirs_b, zero_flux_b, 1, false);

  struct gkyl_dg_pkpm_fluid_flux_surf_inp finp = {
    .conf_grid = &s->cgrid,
    .conf_basis = &s->cbasis,
    .wv_eqn = s->ten,
    .geom = s->geom,
    .use_gpu = false,
  };
  s->calc_fflux = gkyl_dg_pkpm_fluid_flux_surf_inew(&finp);
  s->nkf = gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(s->calc_fflux);
  s->fflux = mkarr(3 * s->nkf, cvol);
  struct gkyl_dg_pkpm_fluid_inp minp = {
    .conf_basis = &s->cbasis,
    .conf_range = &s->clocal,
    .auxfields = {.moms = s->moms, .bvar = s->bvar, .u = s->u, .fluid_flux_surf = s->fflux},
    .use_gpu = false,
  };
  s->eqn_m = gkyl_dg_pkpm_fluid_inew(&minp);
  int up_dirs_m[GKYL_MAX_DIM] = {0}, zero_flux_m[2 * GKYL_MAX_DIM] = {0};
  s->slvr_m = gkyl_hyper_dg_new(&s->cgrid, &s->cbasis, s->eqn_m, 1, up_dirs_m, zero_flux_m, 1, false);
}

static void
sim_release(struct sim *s)
{
  gkyl_hyper_dg_release(s->slvr_m);
  gkyl_dg_eqn_release(s->eqn_m);
  gkyl_array_release(s->fflux);
  gkyl_dg_pkpm_fluid_flux_surf_release(s->calc_fflux);
  gkyl_hyper_dg_release(s->slvr_b);
  gkyl_dg_eqn_release(s->eqn_b);
  gkyl_array_release(s->vflux);
  gkyl_dg_pkpm_vel_flux_surf_release(s->calc_vflux);
  gkyl_hyper_dg_release(s->slvr_a);
  gkyl_dg_eqn_release(s->eqn_a);
  gkyl_array_release(s->cflux);
  gkyl_array_release(s->kin_mf);
  gkyl_array_release(s->kin_mf_face);
  gkyl_array_release(s->roe_mf);
  gkyl_dg_pkpm_conf_flux_surf_release(s->calc_cflux);
  gkyl_dg_updater_moment_pkpm_release(s->mcalc);
  gkyl_dg_calc_pkpm_vars_release(s->vars);
  gkyl_wave_geom_release(s->geom);
  gkyl_wv_eqn_release(s->ten);
  gkyl_dg_calc_pkpm_bvar_release(s->calc_b);
  gkyl_array_release(s->cell_avg);
  gkyl_array_release(s->cflrate_m);
  gkyl_array_release(s->cflrate_f);
  gkyl_array_release(s->rhs_moms);
  gkyl_array_release(s->moms);
  gkyl_array_release(s->a_P);
  gkyl_array_release(s->nu_vthsq_sum);
  gkyl_array_release(s->nu_sum);
  gkyl_array_release(s->theta);
  gkyl_array_release(s->grad_u);
  gkyl_array_release(s->u);
  gkyl_array_release(s->grad_b);
  gkyl_array_release(s->bvar);
}

// Primitive variables of the stage (ghost cells included), then the right-hand sides.
// Returns the maximum CFL rate.
static double
sim_rhs(
  struct sim *s, const struct run_opts *o, const struct gkyl_array *f, const struct gkyl_array *M,
  struct gkyl_array *rhs_f, struct gkyl_array *rhs_M
)
{
  gkyl_array_clear(s->moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(s->mcalc, &s->plocal_ext, &s->clocal_ext, f, s->moms);
  if (o->freeze_fluid) {
    // Frozen fluid: the flow is exactly the uniform u0 (the kinetic update alone).
    gkyl_array_clear(s->u, 0.0);
    gkyl_array_shiftc(s->u, o->u0 * sqrt(2.0), 0);
  } else {
    gkyl_dg_calc_pkpm_vars_u(s->vars, s->moms, M, s->cell_avg, s->u);
  }
  gkyl_dg_calc_pkpm_vars_theta(s->vars, &s->clocal_ext, s->moms, s->theta);
  gkyl_dg_calc_pkpm_vars_grad_u(s->vars, &s->clocal, s->u, s->grad_u);

  gkyl_array_clear(s->cflrate_f, 0.0);
  gkyl_array_clear(s->cflrate_m, 0.0);
  gkyl_array_clear(rhs_f, 0.0);
  gkyl_array_clear(rhs_M, 0.0);

  // Momentum face fluxes first (their mass flux constrains the kinetic flux).
  gkyl_array_clear(s->fflux, 0.0);
  gkyl_dg_pkpm_fluid_flux_surf_advance(
    s->calc_fflux, &s->clocal, &s->clocal_ext, s->moms, s->bvar, s->u, M, s->cflrate_m, s->fflux,
    s->roe_mf
  );
  if (!o->freeze_kinetic) {
    gkyl_array_clear(s->cflux, 0.0);
    gkyl_dg_pkpm_conf_flux_surf_advance(
      s->calc_cflux, &s->clocal, &s->plocal, &s->plocal_ext, s->bvar, s->u, s->theta, f, s->cflrate_f, s->cflux,
      s->kin_mf
    );
    gkyl_dg_pkpm_conf_flux_surf_constrain(
      s->calc_cflux, &s->clocal, &s->plocal, &s->plocal_ext, 1.0, s->moms, s->roe_mf, s->kin_mf,
      s->kin_mf_face, f, s->cflux
    );
    gkyl_hyper_dg_advance(s->slvr_a, &s->plocal, f, s->cflrate_f, rhs_f);
    gkyl_array_clear(s->rhs_moms, 0.0);
    gkyl_dg_updater_moment_pkpm_advance(s->mcalc, &s->plocal, &s->clocal, rhs_f, s->rhs_moms);
    gkyl_dg_calc_pkpm_vars_p_force(s->vars, &s->clocal, s->moms, s->rhs_moms, s->a_P);
    gkyl_array_clear(s->vflux, 0.0);
    gkyl_dg_pkpm_vel_flux_surf_advance(
      s->calc_vflux, &s->clocal, &s->plocal, s->bvar, s->grad_b, s->u, s->grad_u, s->theta, s->a_P,
      f, s->cflrate_f, s->vflux
    );
    gkyl_hyper_dg_advance(s->slvr_b, &s->plocal, f, s->cflrate_f, rhs_f);
  }
  if (!o->freeze_fluid) {
    gkyl_hyper_dg_advance(s->slvr_m, &s->clocal, M, s->cflrate_m, rhs_M);
  }
  double wf = 0.0, wm = 0.0;
  gkyl_array_reduce_range(&wf, s->cflrate_f, GKYL_MAX, &s->plocal);
  gkyl_array_reduce_range(&wm, s->cflrate_m, GKYL_MAX, &s->clocal);
  return fmax(wf, wm);
}

struct diag {
  double mass, mom, ke, ppar, pperp, du_max, dp_max;
};

struct run_result {
  double dev_max; // max |u - u0| (|p - 1| with the fluid frozen) over the run
  double denergy_max; // max |E(t) - E(0)| of the integrated energy rho u^2/2 + p_par/2 + p_perp
  double dinternal_max; // max drift of the internal energy p_par/2 + p_perp alone
};

static struct diag
sim_diag(struct sim *s, const struct run_opts *o, const struct gkyl_array *f, const struct gkyl_array *M)
{
  int nc = s->nc;
  double dx = s->cgrid.dx[0];
  struct diag d = {0};
  gkyl_array_clear(s->moms, 0.0);
  gkyl_dg_updater_moment_pkpm_advance(s->mcalc, &s->plocal, &s->clocal, f, s->moms);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->clocal);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(&s->clocal_ext, iter.idx);
    const double *m = gkyl_array_cfetch(s->moms, cidx);
    const double *Mx = gkyl_array_cfetch(M, cidx);
    double rho = m[0] / sqrt(2.0), ppar = m[nc] / sqrt(2.0), pperp = m[2 * nc] / sqrt(2.0);
    double mx = Mx[0] / sqrt(2.0), my = Mx[nc] / sqrt(2.0), mz = Mx[2 * nc] / sqrt(2.0);
    d.mass += rho * dx;
    d.mom += mx * dx;
    d.ke += (mx * mx + my * my + mz * mz) / rho * dx;
    d.ppar += ppar * dx;
    d.pperp += pperp * dx;
    double du = fabs(mx / rho - o->u0), dp = fabs(ppar - 1.0);
    d.du_max = du > d.du_max ? du : d.du_max;
    d.dp_max = dp > d.dp_max ? dp : d.dp_max;
  }
  return d;
}

// Runs the coupled system and returns the maximum deviation of the cell-average flow
// from u0 over the run (of the parallel pressure from 1 when the fluid is frozen) and
// the maximum drift of the integrated energy; the deviation is 1e300 if the run blew up.
static struct run_result
run(const struct run_opts *o, bool verbose)
{
  struct sim s;
  sim_init(&s, o);
  int np = s.np, nc = s.nc;
  long pvol = s.plocal_ext.volume, cvol = s.clocal_ext.volume;
  struct wave_ctx ctx = {
    .u0 = o->u0, .alpha = o->alpha, .vt = 1.0, .uy_amp = o->uy_amp, .ux_amp = o->ux_amp};

  struct gkyl_array *f = mkarr(2 * np, pvol), *f1 = mkarr(2 * np, pvol), *fnew = mkarr(2 * np, pvol);
  struct gkyl_array *M = mkarr(3 * nc, cvol), *M1 = mkarr(3 * nc, cvol), *Mnew = mkarr(3 * nc, cvol);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&s.pgrid, &s.pbasis, 8, 2, eval_f, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.plocal_ext, f);
  gkyl_proj_on_basis_release(proj);
  proj = gkyl_proj_on_basis_new(&s.cgrid, &s.cbasis, 8, 3, eval_euler, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &s.clocal_ext, M);
  gkyl_proj_on_basis_release(proj);
  periodic(f, s.pskin, s.pghost);
  periodic(M, s.cskin, s.cghost);

  double t = 0.0, du_run = 0.0, de_run = 0.0, di_run = 0.0, t_print = 0.0;
  int step = 0;
  struct diag d0 = sim_diag(&s, o, f, M);
  if (verbose) {
    printf("  Nx=%d Nv=%d u0=%g alpha=%g freeze_fluid=%d freeze_kinetic=%d\n", o->Nx, o->Nv, o->u0,
      o->alpha, o->freeze_fluid, o->freeze_kinetic);
    printf("  %8s %6s %10s %10s %10s %10s %10s %10s\n", "t", "step", "dmass", "dKE", "dppar", "dpperp",
      "max|u-u0|", "max|p-1|");
  }
  while (t < o->t_end) {
    // SSP-RK3.
    double omega = sim_rhs(&s, o, f, M, fnew, Mnew);
    double dt = fmin(1.0 / omega, o->t_end - t);
    gkyl_array_accumulate(gkyl_array_scale(fnew, dt), 1.0, f);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, dt), 1.0, M);
    periodic(fnew, s.pskin, s.pghost);
    periodic(Mnew, s.cskin, s.cghost);
    gkyl_array_copy(f1, fnew);
    gkyl_array_copy(M1, Mnew);

    sim_rhs(&s, o, f1, M1, fnew, Mnew);
    gkyl_array_accumulate(gkyl_array_scale(fnew, dt), 1.0, f1);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, dt), 1.0, M1);
    gkyl_array_accumulate(gkyl_array_scale(fnew, 0.25), 0.75, f);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, 0.25), 0.75, M);
    periodic(fnew, s.pskin, s.pghost);
    periodic(Mnew, s.cskin, s.cghost);
    gkyl_array_copy(f1, fnew);
    gkyl_array_copy(M1, Mnew);

    sim_rhs(&s, o, f1, M1, fnew, Mnew);
    gkyl_array_accumulate(gkyl_array_scale(fnew, dt), 1.0, f1);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, dt), 1.0, M1);
    gkyl_array_accumulate(gkyl_array_scale(fnew, 2.0 / 3.0), 1.0 / 3.0, f);
    gkyl_array_accumulate(gkyl_array_scale(Mnew, 2.0 / 3.0), 1.0 / 3.0, M);
    periodic(fnew, s.pskin, s.pghost);
    periodic(Mnew, s.cskin, s.cghost);
    gkyl_array_copy(f, fnew);
    gkyl_array_copy(M, Mnew);

    t += dt;
    step += 1;
    struct diag d = sim_diag(&s, o, f, M);
    // With the fluid frozen the flow is fixed, so the pressure deviation is the measure.
    double dev = o->freeze_fluid ? d.dp_max : d.du_max;
    du_run = dev > du_run ? dev : du_run;
    double de = fabs((0.5 * d.ke + 0.5 * d.ppar + d.pperp) - (0.5 * d0.ke + 0.5 * d0.ppar + d0.pperp));
    de_run = de > de_run ? de : de_run;
    double di = fabs((0.5 * d.ppar + d.pperp) - (0.5 * d0.ppar + d0.pperp));
    di_run = di > di_run ? di : di_run;
    if (verbose && (t >= t_print || t >= o->t_end)) {
      printf("  %8.4f %6d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e\n", t, step, d.mass - d0.mass,
        d.ke - d0.ke, d.ppar - d0.ppar, d.pperp - d0.pperp, d.du_max, d.dp_max);
      t_print += o->t_end / 10.0;
    }
    if (!(du_run < 10.0)) {
      du_run = 1e300;
      break;
    }
  }
  if (verbose) {
    // Cell averages and slopes of rho and u_x at the end: a slope much larger than the
    // neighbor difference / 2 sqrt(3/2) is a sub-cell sawtooth.
    gkyl_array_clear(s.moms, 0.0);
    gkyl_dg_updater_moment_pkpm_advance(s.mcalc, &s.plocal, &s.clocal, f, s.moms);
    gkyl_dg_calc_pkpm_vars_u(s.vars, s.moms, M, s.cell_avg, s.u);
    printf("  final rho avg / slope, ux avg / slope per cell:\n");
    struct gkyl_range_iter it;
    gkyl_range_iter_init(&it, &s.clocal);
    while (gkyl_range_iter_next(&it)) {
      long cidx = gkyl_range_idx(&s.clocal_ext, it.idx);
      const double *m = gkyl_array_cfetch(s.moms, cidx), *u = gkyl_array_cfetch(s.u, cidx);
      printf("   %3d  %9.5f %9.5f   %9.5f %9.5f\n", it.idx[0], m[0] / sqrt(2.0), m[1], u[0] / sqrt(2.0), u[1]);
    }
  }

  gkyl_array_release(Mnew);
  gkyl_array_release(M1);
  gkyl_array_release(M);
  gkyl_array_release(fnew);
  gkyl_array_release(f1);
  gkyl_array_release(f);
  sim_release(&s);
  return (struct run_result){.dev_max = du_run, .denergy_max = de_run, .dinternal_max = di_run};
}

static bool
verbose(void)
{
  return test_verbose_level_ >= 3;
}

// Coupled runs: the flow stays bounded (the collisionless temperature perturbation phase
// mixes into a sound wave of modest amplitude) and the energy drift is small.
static void
check_coupled(const struct run_opts *o)
{
  struct run_result r = run(o, verbose());
  TEST_CHECK(r.dev_max < 0.2);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
  TEST_CHECK(r.denergy_max < 1e-3);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

void
test_coupled_u1()
{
  struct run_opts o = {.u0 = 1.0, .alpha = 0.2, .t_end = 0.6, .Nx = 16, .Nv = 16};
  check_coupled(&o);
}

void
test_coupled_u0()
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.2, .t_end = 0.6, .Nx = 16, .Nv = 16};
  check_coupled(&o);
}

// Oblique, x-dependent static field (div b != 0): the mirror force exchanges parallel
// and perpendicular pressure in a plasma initially at rest; the total energy must be
// conserved up to the (small) dissipation of the fluxes.
void
test_mirror_oblique()
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.0, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
  TEST_CHECK(r.denergy_max < 1e-2);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

// Oblique varying field with a flowing density pulse: streaming along b, mirror
// forces and the K = bb:grad(u) work all act at once.
void
test_mirror_oblique_flow()
{
  struct run_opts o = {.u0 = 0.5, .alpha = 0.2, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
  TEST_CHECK(r.denergy_max < 1e-2);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

// Same field, kinetic update alone (u = u0 exactly): without compression work the
// internal energy p_par/2 + p_perp must be conserved by the mirror-force pairing.
void
test_mirror_oblique_frozen_fluid()
{
  struct run_opts o = {.u0 = 0.5, .alpha = 0.2, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5, .freeze_fluid = true};
  struct run_result r = run(&o, verbose());
  // The momentum is not evolved here, so only the internal energy is meaningful.
  TEST_CHECK(r.dinternal_max < 1e-5);
  TEST_MSG("max internal energy drift over the run: %g", r.dinternal_max);
}

// Sheared perpendicular drift through an oblique uniform field: K = b_x b_y du_y/dx
// exchanges perpendicular kinetic energy with the pressures through P_xy; the total
// energy must be conserved.
static void
check_shear(int Nx)
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.0, .t_end = 2.0, .Nx = Nx, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.0, .uy_amp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-2);
  TEST_MSG("Nx=%d: max energy drift over the run: %g", Nx, r.denergy_max);
}
void test_shear_nx16() { check_shear(16); }

// Same shear with a nearly perpendicular field: the streaming speed along x is
// small for every particle, so the slope modes get little upwind damping.
// Compressive flow varying along x across a nearly perpendicular field: the upwind
// switch w* = -u_x/b_x of the face flux sits inside different velocity cells at the
// two faces of a cell, so its quadrature error acts on the slope modes.
static void
check_compress(double theta, double amp, bool b_along_x)
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.0, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = b_along_x ? 0.0 : theta, .b_perp = 0.0, .ux_amp = amp};
  if (b_along_x) {
    o.b_theta = 0.0;
  }
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-2);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}
void test_compress_quasi_perp() { check_compress(0.0871, 0.1, false); }
void test_compress_parallel() { check_compress(0.0, 0.1, true); }

void
test_shear_quasi_perp()
{
  struct run_opts o = {.u0 = 0.0, .alpha = 0.0, .t_end = 2.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.0871, .b_perp = 0.0, .uy_amp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-2);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

// Nearly perpendicular field (b_x = 0.087, as in the kinetic Alfven wave input) with a
// slow flow comparable to the parallel streaming speed along x: a slope-mode check.
void
test_quasi_perp()
{
  struct run_opts o = {.u0 = 0.05, .alpha = 0.2, .t_end = 20.0, .Nx = 16, .Nv = 16,
    .b_theta = 0.0871, .b_perp = 0.0};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
  TEST_CHECK(r.denergy_max < 1e-2);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}
void test_shear_nx32() { check_shear(32); }

void
test_mirror_oblique_flow_nx32()
{
  struct run_opts o = {.u0 = 0.5, .alpha = 0.2, .t_end = 2.0, .Nx = 32, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-2);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

void
test_mirror_oblique_flow_nx64()
{
  struct run_opts o = {.u0 = 0.5, .alpha = 0.2, .t_end = 2.0, .Nx = 64, .Nv = 16,
    .b_theta = 0.3, .b_perp = 0.5};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.denergy_max < 1e-2);
  TEST_MSG("max energy drift over the run: %g", r.denergy_max);
}

// Finer grid: the inconsistent coupling grew at a rate ~ u/dx, so this run is the
// sharper instability check.
void
test_coupled_u1_nx64()
{
  struct run_opts o = {.u0 = 1.0, .alpha = 0.2, .t_end = 0.4, .Nx = 64, .Nv = 16};
  check_coupled(&o);
}

void
test_frozen_fluid_u1()
{
  struct run_opts o = {
    .u0 = 1.0, .alpha = 0.2, .t_end = 0.6, .Nx = 16, .Nv = 16, .freeze_fluid = true};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |p - 1| over the run: %g", r.dev_max);
}

void
test_frozen_kinetic_u1()
{
  struct run_opts o = {
    .u0 = 1.0, .alpha = 0.2, .t_end = 0.6, .Nx = 16, .Nv = 16, .freeze_kinetic = true};
  struct run_result r = run(&o, verbose());
  TEST_CHECK(r.dev_max < 1.0);
  TEST_MSG("max |u - u0| over the run: %g", r.dev_max);
}

TEST_LIST = {
  {"coupled_u1", test_coupled_u1},
  {"coupled_u0", test_coupled_u0},
  {"coupled_u1_nx64", test_coupled_u1_nx64},
  {"mirror_oblique", test_mirror_oblique},
  {"mirror_oblique_flow", test_mirror_oblique_flow},
  {"mirror_oblique_frozen_fluid", test_mirror_oblique_frozen_fluid},
  {"shear_nx16", test_shear_nx16},
  {"quasi_perp", test_quasi_perp},
  {"shear_quasi_perp", test_shear_quasi_perp},
  {"compress_quasi_perp", test_compress_quasi_perp},
  {"compress_parallel", test_compress_parallel},
  {"shear_nx32", test_shear_nx32},
  {"mirror_oblique_flow_nx32", test_mirror_oblique_flow_nx32},
  {"mirror_oblique_flow_nx64", test_mirror_oblique_flow_nx64},
  {"frozen_fluid_u1", test_frozen_fluid_u1},
  {"frozen_kinetic_u1", test_frozen_kinetic_u1},
  {NULL, NULL},
};
