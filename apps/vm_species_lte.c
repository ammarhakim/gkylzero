#include <assert.h>
#include <gkyl_vlasov_priv.h>

void 
vm_species_lte_init(struct gkyl_vlasov_app *app, struct vm_species *s, struct vm_lte *lte, bool correct_all_moms)
{
  int cdim = app->cdim, vdim = app->vdim;

  // allocate moments needed for lte update
  vm_species_moment_init(app, s, &lte->moms, "LTEMoments");

  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .phase_basis_on_dev = app->basis_on_dev.basis, // pointer to (potentially) device-side basis
    .conf_basis_on_dev = app->basis_on_dev.confBasis, // pointer to (potentially) device-side conf basis
    .conf_range =  &app->local,
    .conf_range_ext = &app->local_ext,
    .vel_range = &s->local_vel,
    .p_over_gamma = s->p_over_gamma,
    .gamma = s->gamma,
    .gamma_inv = s->gamma_inv,
    .h_ij_inv = s->h_ij_inv,
    .det_h = s->det_h,
    .model_id = s->model_id,
    .mass = s->info.mass,
    .use_gpu = app->use_gpu,
  };
  lte->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew( &inp_proj );

  lte->correct_all_moms = false;
  int max_iter = s->info.max_iter > 0 ? s->info.max_iter : 100;
  double iter_eps = s->info.iter_eps > 0 ? s->info.iter_eps  : 1e-12;
  
  if (correct_all_moms) {
    lte->correct_all_moms = true;

    struct gkyl_vlasov_lte_correct_inp inp_corr = {
      .phase_grid = &s->grid,
      .conf_basis = &app->confBasis,
      .phase_basis = &app->basis,
      .phase_basis_on_dev = app->basis_on_dev.basis, // pointer to (potentially) device-side basis
      .conf_basis_on_dev = app->basis_on_dev.confBasis, // pointer to (potentially) device-side conf basis
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &s->local_vel,
      .p_over_gamma = s->p_over_gamma,
      .gamma = s->gamma,
      .gamma_inv = s->gamma_inv,
      .h_ij_inv = s->h_ij_inv,
      .det_h = s->det_h,
      .model_id = s->model_id,
      .mass = s->info.mass,
      .use_gpu = app->use_gpu,
      .max_iter = max_iter,
      .eps = iter_eps,
    };
    lte->self_niter = 0;
    lte->corr_lte = gkyl_vlasov_lte_correct_inew( &inp_corr );

    lte->corr_stat = gkyl_dynvec_new(GKYL_DOUBLE,7);
    lte->is_first_corr_status_write_call = true;
  }

  lte->f_lte = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
}

// computes moments
void
vm_species_lte_moms(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lte *lte, const struct gkyl_array *fin)
{
  struct timespec wst = gkyl_wall_clock();

  vm_species_moment_calc(&lte->moms, species->local, app->local, fin);
  
  app->stat.species_lte_tm += gkyl_time_diff_now_sec(wst);    
}

// updates f_lte to the current fin
void
vm_species_lte(gkyl_vlasov_app *app, const struct vm_species *species,
  struct vm_lte *lte, const struct gkyl_array *fin)
{

  // Always update the moments
  vm_species_lte_moms(app, species, lte, fin);

  struct timespec wst = gkyl_wall_clock();
  gkyl_array_clear(lte->f_lte, 0.0);

  // Project the LTE distribution function to obtain f_lte.
  // e.g., Maxwellian for non-relativistic and Maxwell-Juttner for relativistic.
  // Projection routine also corrects the density of the projected distribution function.
  gkyl_vlasov_lte_proj_on_basis_advance(lte->proj_lte, &species->local, &app->local, 
    lte->moms.marr, lte->f_lte);

  // Correct all the moments of the projected LTE distribution function.
  if (lte->correct_all_moms) {
    struct gkyl_vlasov_lte_correct_status status_corr = gkyl_vlasov_lte_correct_all_moments(lte->corr_lte, lte->f_lte, lte->moms.marr,
      &species->local, &app->local);
    double corr_vec[7];
    corr_vec[0] = status_corr.num_iter;
    corr_vec[1] = status_corr.iter_converged;
    corr_vec[2] = status_corr.error[0];
    corr_vec[3] = status_corr.error[1];
    corr_vec[4] = status_corr.error[2];
    corr_vec[5] = status_corr.error[3];
    corr_vec[6] = status_corr.error[4];
    gkyl_dynvec_append(lte->corr_stat,app->tcurr,corr_vec);

    lte->self_niter += status_corr.num_iter;
  } 

  app->stat.species_lte_tm += gkyl_time_diff_now_sec(wst);
}

void 
vm_species_lte_release(const struct gkyl_vlasov_app *app, const struct vm_lte *lte)
{
  gkyl_array_release(lte->f_lte);

  vm_species_moment_release(app, &lte->moms);

  gkyl_vlasov_lte_proj_on_basis_release(lte->proj_lte);
  if (lte->correct_all_moms) {
    gkyl_vlasov_lte_correct_release(lte->corr_lte);
    gkyl_dynvec_release(lte->corr_stat);
  }
}
