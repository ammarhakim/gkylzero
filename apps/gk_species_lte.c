#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_lte_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_lte *lte, 
  struct correct_all_moms_inp corr_inp)
{
  int cdim = app->cdim, vdim = app->vdim;

  // allocate moments needed for Maxwellian (LTE=local thermodynamic equilibrium) update
  gk_species_moment_init(app, s, &lte->moms, "MaxwellianMoments");

  struct gkyl_gk_maxwellian_proj_on_basis_inp inp_proj = {
    .phase_grid = &s->grid,
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .conf_range =  &app->local,
    .conf_range_ext = &app->local_ext,
    .vel_range = &s->local_vel,
    .gk_geom = app->gk_geom,
    .vel_map = s->vel_map,
    .mass = s->info.mass,
    .use_gpu = app->use_gpu,
  };
  lte->proj_max = gkyl_gk_maxwellian_proj_on_basis_inew( &inp_proj );

  lte->correct_all_moms = corr_inp.correct_all_moms;
  int max_iter = corr_inp.max_iter > 0 ? s->info.max_iter : 50;
  double iter_eps = corr_inp.iter_eps > 0 ? s->info.iter_eps  : 1e-10;
  bool use_last_converged = corr_inp.use_last_converged;
  
  if (lte->correct_all_moms) {
    struct gkyl_gk_maxwellian_correct_inp inp_corr = {
      .phase_grid = &s->grid,
      .conf_basis = &app->confBasis,
      .phase_basis = &app->basis,
      .conf_range =  &app->local,
      .conf_range_ext = &app->local_ext,
      .vel_range = &s->local_vel, 
      .gk_geom = app->gk_geom,
      .vel_map = s->vel_map,
      .mass = s->info.mass,
      .max_iter = max_iter,
      .eps = iter_eps,
      .use_last_converged = use_last_converged, 
      .use_gpu = app->use_gpu,
    };
    lte->niter = 0; // total number of iterations from correcting moments
    lte->num_corr = 0; // total number of times the correction updater is called
    lte->corr_max = gkyl_gk_maxwellian_correct_inew( &inp_corr );

    lte->corr_stat = gkyl_dynvec_new(GKYL_DOUBLE, 5);
    lte->is_first_corr_status_write_call = true;
  }

  lte->f_lte = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
}

// Compute f_lte from input Maxwellian (LTE=local thermodynamic equilibrium) moments
void
gk_species_lte_from_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lte *lte, const struct gkyl_array *moms_lte)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_array_clear(lte->f_lte, 0.0);

  // Project the Maxwellian distribution function to obtain f_lte.
  // Projection routine also corrects the density of the projected distribution function.
  gkyl_gk_maxwellian_proj_on_basis_advance(lte->proj_max, &species->local, &app->local, 
    moms_lte, false, lte->f_lte);

  // Correct all the moments of the projected Maxwellian distribution function.
  if (lte->correct_all_moms) {
    struct gkyl_gk_maxwellian_correct_status status_corr;
    status_corr = gkyl_gk_maxwellian_correct_all_moments(lte->corr_max, lte->f_lte, moms_lte,
      &species->local, &app->local);
    double corr_vec[5] = { 0.0 };
    corr_vec[0] = status_corr.num_iter;
    corr_vec[1] = status_corr.iter_converged;
    corr_vec[2] = status_corr.error[0];
    corr_vec[3] = status_corr.error[1];
    corr_vec[4] = status_corr.error[2];
    double corr_vec_global[5] = { 0.0 };
    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 5, corr_vec, corr_vec_global);    
    gkyl_dynvec_append(lte->corr_stat, app->tcurr, corr_vec_global);

    lte->niter += status_corr.num_iter;
    lte->num_corr += 1;
  } 

  app->stat.species_lte_tm += gkyl_time_diff_now_sec(wst);   
}

// Compute equivalent f_lte from fin
void
gk_species_lte(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_lte *lte, const struct gkyl_array *fin)
{
// compute needed Maxwellian moments (J*n, u_par, T/m)   
  gk_species_moment_calc(&lte->moms, species->local, app->local, fin);
  
  // divide out the Jacobian from the density
  gkyl_dg_div_op_range(lte->moms.mem_geo, app->confBasis, 
    0, lte->moms.marr, 0, lte->moms.marr, 0, 
    app->gk_geom->jacobgeo, &app->local);  

  gk_species_lte_from_moms(app, species, lte, lte->moms.marr);
}

void 
gk_species_lte_release(const struct gkyl_gyrokinetic_app *app, const struct gk_lte *lte)
{
  gkyl_array_release(lte->f_lte);

  gk_species_moment_release(app, &lte->moms);

  gkyl_gk_maxwellian_proj_on_basis_release(lte->proj_max);
  if (lte->correct_all_moms) {
    gkyl_gk_maxwellian_correct_release(lte->corr_max);
    gkyl_dynvec_release(lte->corr_stat);
  }
}