#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_const.h>

#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_calc_gk_rad_vars_priv.h>

void 
gk_species_radiation_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_rad_drag *rad)
{
  rad->radiation_id = s->info.radiation.radiation_id;
  int cdim = app->cdim, vdim = app->vdim;
  int pdim = cdim+vdim;
  // make appropriate reduced bases and surface bases for radiation variables
  struct gkyl_basis rad_basis, surf_rad_vpar_basis, surf_rad_mu_basis, surf_vpar_basis, surf_mu_basis;
  if (app->poly_order > 1) {
    // radiation variables have no dependence on y since B = B(x,z)
    if (cdim == 3) {
      gkyl_cart_modal_serendip(&rad_basis, pdim-1, app->poly_order);
      gkyl_cart_modal_serendip(&surf_rad_vpar_basis, pdim-2, app->poly_order);
      gkyl_cart_modal_serendip(&surf_rad_mu_basis, pdim-2, app->poly_order);
    }
    else {
      gkyl_cart_modal_serendip(&rad_basis, pdim, app->poly_order);
      gkyl_cart_modal_serendip(&surf_rad_vpar_basis, pdim-1, app->poly_order);
      gkyl_cart_modal_serendip(&surf_rad_mu_basis, pdim-1, app->poly_order);
    }
    gkyl_cart_modal_serendip(&surf_vpar_basis, pdim-1, app->poly_order);
    gkyl_cart_modal_serendip(&surf_mu_basis, pdim-1, app->poly_order);
  }
  else {
    // radiation variables have no dependence on y since B = B(x,z)
    if (cdim == 3) {
      gkyl_cart_modal_gkhybrid(&rad_basis, cdim-1, vdim);
      // constant vparallel surface, only depends on (x,z,mu)
      gkyl_cart_modal_serendip(&surf_rad_vpar_basis, pdim-2, app->poly_order);
      // constant mu surface, only depends on (x,z,vpar), vpar is poly_order = 2
      gkyl_cart_modal_gkhybrid(&surf_rad_mu_basis, cdim-1, 1);
    }
    else {
      gkyl_cart_modal_gkhybrid(&rad_basis, cdim, vdim);
      gkyl_cart_modal_serendip(&surf_rad_vpar_basis, pdim-1, app->poly_order);
      gkyl_cart_modal_gkhybrid(&surf_rad_mu_basis, cdim, 1);
    }
    gkyl_cart_modal_serendip(&surf_vpar_basis, pdim-1, app->poly_order);
    gkyl_cart_modal_gkhybrid(&surf_mu_basis, cdim, 1);   
  }

  // Updater to compute drag coefficients.
  rad->calc_gk_rad_vars = gkyl_dg_calc_gk_rad_vars_new(&s->grid, &app->confBasis,
	&app->basis, s->info.charge, s->info.mass, app->gk_geom, s->vel_map, app->use_gpu);

  // Fitting parameters
  double rad_fit_a[GKYL_MAX_RAD_DENSITIES], rad_fit_alpha[GKYL_MAX_RAD_DENSITIES], rad_fit_beta[GKYL_MAX_RAD_DENSITIES],
    rad_fit_gamma[GKYL_MAX_RAD_DENSITIES], rad_fit_v0[GKYL_MAX_RAD_DENSITIES], rad_fit_ne[GKYL_MAX_RAD_DENSITIES];
  struct all_radiation_states *rad_data = gkyl_radiation_read_rad_fit_params();

  rad->num_cross_collisions = s->info.radiation.num_cross_collisions;
  int num_dens_per_coll[rad->num_cross_collisions];
  for (int i=0; i<rad->num_cross_collisions; ++i) {
    num_dens_per_coll[i] = s->info.radiation.num_of_densities[i] ? s->info.radiation.num_of_densities[i] : 1;
    int status = gkyl_radiation_read_get_num_densities(*rad_data, s->info.radiation.z[i],
      s->info.radiation.charge_state[i], s->info.radiation.min_ne, s->info.radiation.max_ne, &num_dens_per_coll[i]);
  }

  // Allocate drag coefificents.
  rad->vnu_surf = gkyl_dg_calc_gk_rad_vars_drag_new(rad->num_cross_collisions, num_dens_per_coll,
    surf_rad_vpar_basis.num_basis, s->local_ext.volume, app->use_gpu);
  rad->vnu = gkyl_dg_calc_gk_rad_vars_drag_new(rad->num_cross_collisions, num_dens_per_coll,
    rad_basis.num_basis, s->local_ext.volume, app->use_gpu);
  rad->vsqnu_surf = gkyl_dg_calc_gk_rad_vars_drag_new(rad->num_cross_collisions, num_dens_per_coll,
    surf_rad_mu_basis.num_basis, s->local_ext.volume, app->use_gpu);
  rad->vsqnu = gkyl_dg_calc_gk_rad_vars_drag_new(rad->num_cross_collisions, num_dens_per_coll,
    rad_basis.num_basis, s->local_ext.volume, app->use_gpu);

  int max_num_densities = num_dens_per_coll[0];
  for (int i=1; i<rad->num_cross_collisions; i++) {
    if (num_dens_per_coll[i] > max_num_densities)
      max_num_densities = num_dens_per_coll[i];
  }
  
  // Make array for cutoff below which radiation is set to 0. Keep the radiation from driving Te negative.
  rad->vtsq_min_per_species = gkyl_malloc(rad->num_cross_collisions*sizeof(struct gkyl_array*));
  
  // Initialize drag coefficients.
  for (int i=0; i<rad->num_cross_collisions; ++i) {
    int num_densities = s->info.radiation.num_of_densities[i] ? s->info.radiation.num_of_densities[i] : 1;
    int status = gkyl_radiation_read_get_fit_params(*rad_data, s->info.radiation.z[i], s->info.radiation.charge_state[i],
      rad_fit_a, rad_fit_alpha, rad_fit_beta, rad_fit_gamma, rad_fit_v0, &num_densities, rad_fit_ne,
      s->info.radiation.reference_ne, s->info.radiation.min_ne, s->info.radiation.max_ne);
    assert(num_densities == num_dens_per_coll[i]); // Consistency check.
    rad->vtsq_min_per_species[i] = mkarr(app->use_gpu, 1, num_densities);
    rad->rad_fit_ne[i] = mkarr(app->use_gpu, 1, num_densities);
    struct gkyl_array *rad_fit_ne_host = mkarr(false, 1, num_densities);
    struct gkyl_array *vtsq_min_host = mkarr(false, 1, num_densities);
    memcpy(rad_fit_ne_host->data, rad_fit_ne, num_densities*sizeof(double));
    gkyl_array_copy(rad->rad_fit_ne[i], rad_fit_ne_host);
    gkyl_array_release(rad_fit_ne_host);

    // Fetch the species we are colliding with and the fitting parameters for that species
    rad->collide_with_idx[i] = gk_find_species_idx(app, s->info.radiation.collide_with[i]);
    if (rad->collide_with_idx[i] == -1) {
      rad->collide_with_idx[i] = gk_find_neut_species_idx(app, s->info.radiation.collide_with[i]);
      rad->collide_with_neut[i] = gk_find_neut_species(app, s->info.radiation.collide_with[i]);
      rad->is_neut_species[i] = true;
      gk_neut_species_moment_init(app, rad->collide_with_neut[i], &rad->moms[i], "M0");
    }
    else {
      rad->collide_with[i] = gk_find_species(app, s->info.radiation.collide_with[i]);
      rad->is_neut_species[i] = false;
      // allocate density calculation needed for radiation update
      gk_species_moment_init(app, rad->collide_with[i], &rad->moms[i], "M0");
    }

    if (status == 1) {
      char msg[100];
      sprintf(msg, "No radiation fits exist for z=%d, charge state=%d\n",s->info.radiation.z[i], s->info.radiation.charge_state[i]);
      gkyl_gyrokinetic_app_cout(app, stderr, msg);
      exit(EXIT_FAILURE);
    }

    for (int n=0; n<num_densities; n++) {
      // allocate drag coefficients in vparallel and mu for each collision, both surface and volume expansions
      // nu = nu(v) for both vparallel and mu updates, 
      // where |v| = sqrt(v_par^2 + 2 mu B/m)
      // Note that through the spatial variation of B = B(x,z), 
      // both these drag coefficients depend on phase space, but a reduced (x,z,vpar,mu) phase space
      gkyl_dg_calc_gk_rad_vars_nu_advance(rad->calc_gk_rad_vars, &app->local, &s->local,
        rad_fit_a[n], rad_fit_alpha[n], rad_fit_beta[n], rad_fit_gamma[n], rad_fit_v0[n],
        rad->vnu_surf[i].data[n].arr, rad->vnu[i].data[n].arr,
        rad->vsqnu_surf[i].data[n].arr, rad->vsqnu[i].data[n].arr);

      double Te_min_eV;
      if (s->info.radiation.te_min_model == GKYL_CONST_TE && s->info.radiation.Te_min) {
	// Turn off radiation below a constant temperature
	Te_min_eV = s->info.radiation.Te_min / GKYL_ELEMENTARY_CHARGE;
      }
      else if (s->info.radiation.te_min_model == GKYL_VARY_TE_AGGRESSIVE) {
	// Turn off radiation below 10^-4*max(Lz)
	Te_min_eV = 0.1372 * pow(rad_fit_v0[n], 1.867);
      }
      else {
	// (s->info.radiation.te_min_model == GKYL_VARY_TE_CONSERVATIVE) i.e. Turn off radiation below 3.16*10^-3*max(Lz)
	Te_min_eV = 0.2815 * pow(rad_fit_v0[n], 1.768);
      }
      double *vtsq_min_host_d = (double*) gkyl_array_fetch(vtsq_min_host, n);
      vtsq_min_host_d[0] = Te_min_eV * fabs(s->info.charge)/s->info.mass * pow(sqrt(2.0), cdim);
    }
    gkyl_array_copy(rad->vtsq_min_per_species[i], vtsq_min_host);
    gkyl_array_release(vtsq_min_host);
    // Allocate emissivity.
    rad->emissivity[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    rad->emissivity_host[i] = rad->emissivity[i];
    if (app->use_gpu) 
      rad->emissivity_host[i] = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  }

  gkyl_radiation_read_release_fit_params(rad_data);

  // Total vparallel and mu radiation drag including density scaling
  rad->nvnu_surf = mkarr(app->use_gpu, surf_vpar_basis.num_basis, s->local_ext.volume);
  rad->nvnu = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  rad->nvsqnu_surf = mkarr(app->use_gpu, surf_mu_basis.num_basis, s->local_ext.volume);
  rad->nvsqnu = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);

  // allocate moments needed for temperature update
  gk_species_moment_init(app, s, &rad->prim_moms, "MaxwellianMoments");
  rad->vtsq = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  rad->m0 = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  
  rad->nvnu_surf_host = rad->nvnu_surf;
  rad->nvnu_host = rad->nvnu;
  rad->nvsqnu_surf_host = rad->nvsqnu_surf;
  rad->nvsqnu_host = rad->nvsqnu;
  if (app->use_gpu) {
    rad->nvnu_surf_host = mkarr(false, surf_vpar_basis.num_basis, s->local_ext.volume);
    rad->nvnu_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
    rad->nvsqnu_surf_host = mkarr(false, surf_mu_basis.num_basis, s->local_ext.volume);
    rad->nvsqnu_host = mkarr(false, app->basis.num_basis, s->local_ext.volume);
  }

  // Allocate data and updaters for integrated moments.
  gk_species_moment_init(app, s, &rad->integ_moms, "Integrated");
  if (app->use_gpu) {
    rad->red_integ_diag = gkyl_cu_malloc(sizeof(double[vdim+2]));
    rad->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[vdim+2]));
  } 
  else {
    rad->red_integ_diag = gkyl_malloc(sizeof(double[vdim+2]));
    rad->red_integ_diag_global = gkyl_malloc(sizeof(double[vdim+2]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  rad->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, vdim+2);
  rad->is_first_integ_write_call = true;
  // Allocate rhs arry to be used for calculation of integrated moments
  rad->integrated_moms_rhs = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);

  // Arrays for emissivity
  rad->emissivity_rhs = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  rad->emissivity_denominator = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
  // allocate M2 for emissivity
  gk_species_moment_init(app, s, &rad->m2, "M2");

  // Radiation updater
  struct gkyl_dg_rad_gyrokinetic_auxfields drag_inp = { .nvnu_surf = rad->nvnu_surf, .nvnu = rad->nvnu,
    .nvsqnu_surf = rad->nvsqnu_surf, .nvsqnu = rad->nvsqnu};
  rad->drag_slvr = gkyl_dg_updater_rad_gyrokinetic_new(&s->grid, 
    &app->confBasis, &app->basis, &s->local, &app->local, s->vel_map, &drag_inp, app->use_gpu);
}

// computes density for computation of total radiation drag and primitive moments
void
gk_species_radiation_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_rad_drag *rad, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{

  // compute needed Maxwellian moments (n, u_par, T/m) (Jacobian factors already eliminated)
  gk_species_moment_calc(&rad->prim_moms, species->local, app->local, species->f);

  gkyl_array_set_offset(rad->m0, 1.0, rad->prim_moms.marr, 0*app->confBasis.num_basis);
  gkyl_array_set_offset(rad->vtsq, 1.0, rad->prim_moms.marr, 2*app->confBasis.num_basis);

  gkyl_array_clear(rad->nvnu_surf, 0.0);
  gkyl_array_clear(rad->nvnu, 0.0);
  gkyl_array_clear(rad->nvsqnu_surf, 0.0);
  gkyl_array_clear(rad->nvsqnu, 0.0);

  for (int i=0; i<rad->num_cross_collisions; ++i) {
    // compute needed moments
    if (rad->is_neut_species[i])
      gk_neut_species_moment_calc(&rad->moms[i], rad->collide_with_neut[i]->local,
        app->local, fin_neut[rad->collide_with_idx[i]]);
    else
      gk_species_moment_calc(&rad->moms[i], rad->collide_with[i]->local, app->local, fin[rad->collide_with_idx[i]]);
    // divide out Jacobian from ion density before computation of final drag coefficient
    gkyl_dg_div_op_range(rad->moms[i].mem_geo, app->confBasis, 0, rad->moms[i].marr, 0,
      rad->moms[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
    gkyl_dg_calc_gk_rad_vars_nI_nu_advance(
      rad->calc_gk_rad_vars, &app->local, &species->local, 
      &rad->vnu_surf[i], &rad->vnu[i], &rad->vsqnu_surf[i], &rad->vsqnu[i], 
      rad->rad_fit_ne[i], rad->m0, rad->moms[i].marr, 
      rad->nvnu_surf, rad->nvnu, rad->nvsqnu_surf, rad->nvsqnu,
      rad->vtsq_min_per_species[i], rad->vtsq);					   
  }
}

void
gk_species_radiation_integrated_moms(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_rad_drag *rad, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  gkyl_array_clear(rad->integrated_moms_rhs, 0.0);
  gkyl_dg_updater_rad_gyrokinetic_advance(rad->drag_slvr, &species->local,
    species->f, species->cflrate, rad->integrated_moms_rhs);
  gk_species_moment_calc(&rad->integ_moms, species->local, app->local, rad->integrated_moms_rhs);
}

// computes emissivity
void
gk_species_radiation_emissivity(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_rad_drag *rad, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{

  // compute needed Maxwellian moments (n, u_par, T/m) (Jacobian factors already eliminated)
  gk_species_moment_calc(&rad->prim_moms, species->local, app->local, species->f);

  gkyl_array_set_offset(rad->m0, 1.0, rad->prim_moms.marr, 0*app->confBasis.num_basis);
  gkyl_array_set_offset(rad->vtsq, 1.0, rad->prim_moms.marr, 2*app->confBasis.num_basis);

  //Calculate m2
  for (int i=0; i<rad->num_cross_collisions; ++i) {
    gkyl_array_clear(rad->nvnu_surf, 0.0);
    gkyl_array_clear(rad->nvnu, 0.0);
    gkyl_array_clear(rad->nvsqnu_surf, 0.0);
    gkyl_array_clear(rad->nvsqnu, 0.0);
    gkyl_array_clear(rad->emissivity_rhs, 0.0);
    gkyl_array_clear(rad->emissivity_denominator, 0.0);
    if (rad->is_neut_species[i])
      gk_neut_species_moment_calc(&rad->moms[i], rad->collide_with_neut[i]->local, app->local, fin_neut[rad->collide_with_idx[i]]);
    else
      gk_species_moment_calc(&rad->moms[i], rad->collide_with[i]->local, app->local, fin[rad->collide_with_idx[i]]);

    // divide out Jacobian from ion density before computation of final drag coefficient
    gkyl_dg_div_op_range(rad->moms[i].mem_geo, app->confBasis, 0, rad->moms[i].marr, 0,
      rad->moms[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
    gkyl_dg_calc_gk_rad_vars_nI_nu_advance(
      rad->calc_gk_rad_vars, 
      &app->local, &species->local, 
      &rad->vnu_surf[i], &rad->vnu[i],
      &rad->vsqnu_surf[i], &rad->vsqnu[i],
      rad->rad_fit_ne[i], rad->m0, rad->moms[i].marr, 
      rad->nvnu_surf, rad->nvnu, 
      rad->nvsqnu_surf, rad->nvsqnu,
      rad->vtsq_min_per_species[i], rad->vtsq);
    
    gkyl_dg_updater_rad_gyrokinetic_advance(rad->drag_slvr, &species->local,
      species->f, species->cflrate, rad->emissivity_rhs);
    gk_species_moment_calc(&rad->m2, species->local, app->local, rad->emissivity_rhs);

    gkyl_dg_mul_op(app->confBasis, 0, rad->emissivity_denominator, 0, rad->m0, 0, rad->moms[i].marr);
    gkyl_dg_mul_op(app->confBasis, 0, rad->emissivity_denominator, 0, rad->emissivity_denominator, 0, app->gk_geom->jacobgeo);
    gkyl_dg_div_op_range(rad->m2.mem_geo ,app->confBasis, 0, rad->emissivity[i], 0, rad->m2.marr, 0, rad->emissivity_denominator, &app->local);

    rad->emissivity[i] = gkyl_array_scale(rad->emissivity[i], -species->info.mass/2.0);
  }  
}

// updates the collision terms in the rhs
void
gk_species_radiation_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_rad_drag *rad, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  // accumulate update due to collisions onto rhs
  gkyl_dg_updater_rad_gyrokinetic_advance(rad->drag_slvr, &species->local,
    fin, species->cflrate, rhs);
  
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

void 
gk_species_radiation_release(const struct gkyl_gyrokinetic_app *app, const struct gk_rad_drag *rad)
{
  for (int i=0; i<rad->num_cross_collisions; ++i) {
    gkyl_array_release(rad->emissivity[i]);
    if (app->use_gpu)
      gkyl_array_release(rad->emissivity_host[i]);
    gkyl_array_release(rad->rad_fit_ne[i]);
    gkyl_array_release(rad->vtsq_min_per_species[i]);
    gk_species_moment_release(app, &rad->moms[i]);
  }
  gkyl_dg_calc_gk_rad_vars_release(rad->calc_gk_rad_vars);
  gk_species_moment_release(app, &rad->integ_moms); 
  if (app->use_gpu) {
    gkyl_cu_free(rad->red_integ_diag);
    gkyl_cu_free(rad->red_integ_diag_global);
  }
  else {
    gkyl_free(rad->red_integ_diag);
    gkyl_free(rad->red_integ_diag_global);
  }
  gkyl_dg_calc_gk_rad_vars_drag_release(rad->vnu, rad->num_cross_collisions, app->use_gpu);
  gkyl_dg_calc_gk_rad_vars_drag_release(rad->vnu_surf, rad->num_cross_collisions, app->use_gpu);
  gkyl_dg_calc_gk_rad_vars_drag_release(rad->vsqnu, rad->num_cross_collisions, app->use_gpu);
  gkyl_dg_calc_gk_rad_vars_drag_release(rad->vsqnu_surf, rad->num_cross_collisions, app->use_gpu);
  gkyl_dynvec_release(rad->integ_diag);
  gkyl_array_release(rad->integrated_moms_rhs);
  gk_species_moment_release(app, &rad->m2);
  gk_species_moment_release(app, &rad->prim_moms);
  gkyl_array_release(rad->emissivity_denominator);
  gkyl_array_release(rad->emissivity_rhs);
  gkyl_array_release(rad->nvnu_surf);
  gkyl_array_release(rad->nvnu);
  gkyl_array_release(rad->nvsqnu_surf);
  gkyl_array_release(rad->nvsqnu);
  gkyl_array_release(rad->vtsq);
  gkyl_array_release(rad->m0);
  gkyl_free(rad->vtsq_min_per_species);
  if (app->use_gpu) {
    gkyl_array_release(rad->nvnu_surf_host);
    gkyl_array_release(rad->nvnu_host);
    gkyl_array_release(rad->nvsqnu_surf_host);
    gkyl_array_release(rad->nvsqnu_host);    
  }
  gkyl_dg_updater_rad_gyrokinetic_release(rad->drag_slvr);
 }
