#include <assert.h>
#include <gkyl_pkpm_priv.h>

// Evaluate a p=1 configuration-space expansion at the computational coordinate xn.
static double
eval_conf_field(const struct pkpm_bgk_eval_ctx *c, const struct gkyl_array *arr, int comp, const double *xn)
{
  int idx[GKYL_MAX_CDIM];
  double xc[GKYL_MAX_CDIM], eta[GKYL_MAX_CDIM];
  gkyl_rect_grid_coord_idx(c->conf_grid, xn, idx);
  for (int d = 0; d < c->cdim; ++d) {
    idx[d] = GKYL_MAX2(GKYL_MIN2(idx[d], c->local_ext->upper[d]), c->local_ext->lower[d]);
  }
  gkyl_rect_grid_cell_center(c->conf_grid, idx, xc);
  for (int d = 0; d < c->cdim; ++d) {
    eta[d] = 2.0 * (xn[d] - xc[d]) / c->conf_grid->dx[d];
  }
  const double *coeff = gkyl_array_cfetch(arr, gkyl_range_idx(c->local_ext, idx));
  return c->conf_basis->eval_expand(eta, &coeff[comp * c->conf_basis->num_basis]);
}

// The BGK target in the normalized variables: a Maxwellian of variance r = T_0/(m xi^2)
// in w (r = 1 when relaxing to the local temperature) and its isotropic partner G = r F.
static void
eval_bgk_target(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct pkpm_bgk_eval_ctx *c = ctx;
  double w = xn[c->cdim];
  double r = 1.0;
  if (c->fixed_temp_relax) {
    double xisq = eval_conf_field(c, c->xi, 1, xn), vtsq0 = eval_conf_field(c, c->fixed_vtsq, 0, xn);
    r = (xisq > 0.0 && vtsq0 > 0.0) ? vtsq0 / xisq : 1.0;
  }
  double F = exp(-0.5 * w * w / r) / sqrt(2.0 * M_PI * r);
  fout[0] = F;
  fout[1] = r * F;
}

void
pkpm_species_bgk_init(gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk)
{
  int nc = app->confBasis.num_basis, np = app->basis.num_basis;
  // Collision frequency (projected once, like the LBO's).
  bgk->nu = mkarr(app->use_gpu, nc, app->local_ext.volume);
  struct gkyl_array *nu_host = mkarr(false, nc, app->local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->confBasis,
    .num_quad = app->poly_order + 1,
    .num_ret_vals = 1,
    .eval = s->info.collisions.self_nu,
    .ctx = s->info.collisions.ctx,
    .c2p_func = pkpm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local_ext, nu_host);
  gkyl_proj_on_basis_release(proj);
  double nu_frac = s->info.collisions.nuFrac > 0.0 ? s->info.collisions.nuFrac : 1.0;
  gkyl_array_scale(nu_host, nu_frac);
  gkyl_array_copy(bgk->nu, nu_host);
  gkyl_array_release(nu_host);

  bgk->fixed_temp_relax = s->info.collisions.fixed_temp_relax;
  bgk->fixed_vtsq = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(bgk->fixed_vtsq, 0.0);
  bgk->fixed_vtsq_host = app->use_gpu ? mkarr(false, nc, app->local_ext.volume) : gkyl_array_acquire(bgk->fixed_vtsq);
  bgk->xi_host = app->use_gpu ? mkarr(false, 3 * nc, app->local_ext.volume) : 0;
  bgk->f_target = mkarr(app->use_gpu, 2 * np, s->local_ext.volume);
  bgk->f_target_host = app->use_gpu ? mkarr(false, 2 * np, s->local_ext.volume) : gkyl_array_acquire(bgk->f_target);
  gkyl_array_clear(bgk->f_target, 0.0);
  bgk->E_target = mkarr(app->use_gpu, nc, app->local_ext.volume);
  bgk->scratch = mkarr(app->use_gpu, nc, app->local_ext.volume);
  bgk->uu = mkarr(app->use_gpu, nc, app->local_ext.volume);
  bgk->eval_ctx = (struct pkpm_bgk_eval_ctx){
    .cdim = app->cdim,
    .conf_grid = &app->grid,
    .conf_basis = &app->confBasis,
    .local_ext = &app->local_ext,
    .fixed_temp_relax = bgk->fixed_temp_relax,
    .fixed_vtsq = bgk->fixed_vtsq_host,
  };
  bgk->proj_target = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &s->grid,
    .basis = &app->basis,
    .num_quad = 4,
    .num_ret_vals = 2,
    .eval = eval_bgk_target,
    .ctx = &bgk->eval_ctx,
  });
}

// The fixed temperature is the initial one: xi^2 = T_0/m of the initial fluid state (call
// after the initial condition, with the PKPM variables of that state computed).
void
pkpm_species_bgk_set_fixed_temp(gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk)
{
  int nc = app->confBasis.num_basis;
  gkyl_array_set_offset(bgk->fixed_vtsq, 1.0, s->pkpm_xi, nc);
  gkyl_array_copy(bgk->fixed_vtsq_host, bgk->fixed_vtsq);
}

// nu (f_target - f) into the configuration-pass right-hand side; the target is gyrotropic
// ([F, G] only), so the first harmonic [M0_perp, N_perp] relaxes to zero (before the closures:
// with a fixed temperature the target's energy differs from the state's and the closures
// keep the constraints exact while the fluid energy carries the change).
void
pkpm_species_bgk_kinetic_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  struct timespec tm = gkyl_wall_clock();
  if (app->use_gpu) {
    gkyl_array_copy(bgk->xi_host, s->pkpm_xi);
    bgk->eval_ctx.xi = bgk->xi_host;
  } else {
    bgk->eval_ctx.xi = s->pkpm_xi;
  }
  gkyl_proj_on_basis_advance(bgk->proj_target, app->tcurr, &s->local, bgk->f_target_host);
  if (app->use_gpu) {
    gkyl_array_copy(bgk->f_target, bgk->f_target_host);
  }
  gkyl_vlasov_position_map_rescale_jacobpos(app->pos_map, &app->basis, &s->local, bgk->f_target, bgk->f_target);
  pkpm_species_accumulate_conf_phase(app, s, rhs, 1.0, bgk->nu, bgk->f_target);
  pkpm_species_accumulate_conf_phase(app, s, rhs, -1.0, bgk->nu, fin);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
}

// With a fixed temperature the fluid energy relaxes to rho |u|^2/2 + 3 rho xi_0^2/2 at the
// rate nu (density and momentum are conserved by the BGK operator).
void
pkpm_species_bgk_fluid_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk, struct gkyl_array *rhs
)
{
  if (!bgk->fixed_temp_relax) {
    return;
  }
  int nc = app->confBasis.num_basis;
  gkyl_array_clear(bgk->uu, 0.0);
  for (int i = 0; i < 3; ++i) {
    gkyl_dg_mul_op_range(&app->confBasis, 0, bgk->scratch, i, s->pkpm_u, i, s->pkpm_u, &app->local);
    gkyl_array_accumulate(bgk->uu, 1.0, bgk->scratch);
  }
  gkyl_dg_mul_op_range(&app->confBasis, 0, bgk->E_target, 0, s->pkpm_rho, 0, bgk->uu, &app->local);
  gkyl_dg_mul_op_range(&app->confBasis, 0, bgk->scratch, 0, s->pkpm_rho, 0, bgk->fixed_vtsq, &app->local);
  gkyl_array_scale(bgk->E_target, 0.5);
  gkyl_array_accumulate(bgk->E_target, 1.5, bgk->scratch);
  // physical E of the stage: block 4 of fluid_no_J
  gkyl_array_accumulate_offset(bgk->E_target, -1.0, s->fluid_no_J, 4 * nc);
  // the evolved fluid carries the position-map Jacobian
  gkyl_vlasov_position_map_rescale_jacobpos_conf(app->pos_map, &app->local, bgk->E_target, bgk->E_target);
  gkyl_dg_mul_op_range(&app->confBasis, 0, bgk->scratch, 0, bgk->nu, 0, bgk->E_target, &app->local);
  gkyl_array_accumulate_offset_range(rhs, 1.0, bgk->scratch, 4 * nc, &app->local);
}

void
pkpm_species_bgk_release(const gkyl_pkpm_app *app, const struct pkpm_bgk_collisions *bgk)
{
  gkyl_proj_on_basis_release(bgk->proj_target);
  gkyl_array_release(bgk->nu);
  gkyl_array_release(bgk->fixed_vtsq);
  gkyl_array_release(bgk->fixed_vtsq_host);
  gkyl_array_release(bgk->f_target);
  gkyl_array_release(bgk->f_target_host);
  gkyl_array_release(bgk->E_target);
  gkyl_array_release(bgk->scratch);
  gkyl_array_release(bgk->uu);
  if (app->use_gpu) {
    gkyl_array_release(bgk->xi_host);
  }
}
