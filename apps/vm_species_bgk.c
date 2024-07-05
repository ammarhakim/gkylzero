#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_bgk_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_bgk_collisions *bgk)
{
  int cdim = app->cdim, vdim = app->vdim;
  // allocate nu and initialize it
  bgk->nu_sum = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  bgk->self_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  struct gkyl_array *self_nu = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_new(&app->grid, &app->confBasis,
    app->poly_order+1, 1, s->info.collisions.self_nu, s->info.collisions.ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local, self_nu);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_copy(bgk->self_nu, self_nu);
  gkyl_array_copy(bgk->nu_sum, self_nu);
  gkyl_array_release(self_nu);

  bgk->spitzer_calc = 0;
  bgk->normNu = false;
  if (s->info.collisions.normNu) {
    bgk->normNu = true;
    double nuFrac = s->info.collisions.nuFrac ? s->info.collisions.nuFrac : 1.0;
    double eps0 = 1.0;
    double hbar = 1.0;
    bgk->spitzer_calc = gkyl_spitzer_coll_freq_new(&app->confBasis, app->poly_order+1,
      nuFrac, eps0, hbar, app->use_gpu);
    // Create arrays for scaling collisionality by normalization factor
    // norm_nu is computed from Spitzer calc and is the normalization factor for the local
    // density and thermal velocity, norm_nu_sr = n/(vth_s^2 + vth_r^2)^(3/2)
    // nu_init is the inital collisionality profile, which must be stored so that at every time
    // time step the collisionality profile is properly scaled and the effects are not cumulative
    bgk->norm_nu = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    bgk->nu_init = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    gkyl_array_copy(bgk->nu_init, bgk->self_nu);
  }

  // Host-side copy for I/O
  bgk->nu_sum_host = bgk->nu_sum;
  if (app->use_gpu) {
    bgk->nu_sum_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  }

  bgk->model_id = s->model_id;

  int max_iter = s->info.collisions.max_iter > 0 ? s->info.collisions.max_iter : 100;
  double iter_eps = s->info.collisions.iter_eps > 0 ? s->info.collisions.iter_eps  : 1e-12;

  // Allocate everything needed to make f_lte
  vm_species_lte_init(app, s, &bgk->lte, s->info.collisions.correct_all_moms);
  
  bgk->nu_f_lte = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  // BGK updater (also computes stable timestep)
  bgk->up_bgk = gkyl_bgk_collisions_new(&app->confBasis, &app->basis, app->use_gpu);
}

// computes moments
void
vm_species_bgk_moms(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  vm_species_moment_calc(&bgk->lte.moms, species->local, app->local, fin);
  
  app->stat.species_coll_mom_tm += gkyl_time_diff_now_sec(wst);    
}

// updates the collision terms in the rhs
void
vm_species_bgk_rhs(gkyl_vlasov_app *app, struct vm_species *species,
  struct vm_bgk_collisions *bgk, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(bgk->nu_f_lte, 0.0);

  // Project the LTE distribution function
  vm_species_lte(app, species, &bgk->lte, fin);

  gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, bgk->lte.f_lte, 
    bgk->self_nu, bgk->lte.f_lte, &app->local, &species->local);
  gkyl_array_accumulate(bgk->nu_f_lte, 1.0, bgk->lte.f_lte);

  gkyl_bgk_collisions_advance(bgk->up_bgk, &app->local, &species->local, 
    bgk->nu_sum, bgk->nu_f_lte, fin, bgk->implicit_step, bgk->dt_implicit, rhs, species->cflrate);

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void 
vm_species_bgk_release(const struct gkyl_vlasov_app *app, const struct vm_bgk_collisions *bgk)
{
  gkyl_array_release(bgk->self_nu);
  gkyl_array_release(bgk->nu_sum);

  gkyl_array_release(bgk->nu_f_lte);

  if (app->use_gpu) {
    gkyl_array_release(bgk->nu_sum_host);
  }

  if (bgk->normNu) {
    gkyl_array_release(bgk->norm_nu);
    gkyl_array_release(bgk->nu_init);
    gkyl_spitzer_coll_freq_release(bgk->spitzer_calc);
  }

  vm_species_lte_release(app, &bgk->lte);

  gkyl_bgk_collisions_release(bgk->up_bgk);
}
