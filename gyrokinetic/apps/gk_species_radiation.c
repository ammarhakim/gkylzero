#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_const.h>

#include <gkyl_dg_calc_gk_rad_vars.h>

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
  rad->calc_gk_rad_vars = gkyl_dg_calc_gk_rad_vars_new(&s->grid, &app->basis,
	&s->basis, s->info.charge, s->info.mass, app->gk_geom, s->vel_map, app->use_gpu);

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
      gk_neut_species_moment_init(app, rad->collide_with_neut[i], &rad->moms[i], GKYL_F_MOMENT_M0, false);
    }
    else {
      rad->collide_with[i] = gk_find_species(app, s->info.radiation.collide_with[i]);
      rad->is_neut_species[i] = false;
      // allocate density calculation needed for radiation update
      gk_species_moment_init(app, rad->collide_with[i], &rad->moms[i], GKYL_F_MOMENT_M0, false);
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
    rad->emissivity[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    rad->emissivity_host[i] = app->use_gpu? mkarr(false, rad->emissivity[i]->ncomp, rad->emissivity[i]->size)
	                                  : gkyl_array_acquire(rad->emissivity[i]);
  }

  gkyl_radiation_read_release_fit_params(rad_data);

  // Total vparallel and mu radiation drag including density scaling
  rad->nvnu_surf = mkarr(app->use_gpu, surf_vpar_basis.num_basis, s->local_ext.volume);
  rad->nvnu = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
  rad->nvsqnu_surf = mkarr(app->use_gpu, surf_mu_basis.num_basis, s->local_ext.volume);
  rad->nvsqnu = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);

  // Allocate moments needed for temperature update.
  gk_species_moment_init(app, s, &rad->prim_moms, GKYL_F_MOMENT_MAXWELLIAN, false);

  rad->vtsq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  rad->m0 = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  
  if (app->use_gpu) {
    rad->nvnu_surf_host   = mkarr(false, rad->nvnu_surf  ->ncomp, rad->nvnu_surf  ->size); 
    rad->nvnu_host        = mkarr(false, rad->nvnu       ->ncomp, rad->nvnu       ->size); 
    rad->nvsqnu_surf_host = mkarr(false, rad->nvsqnu_surf->ncomp, rad->nvsqnu_surf->size); 
    rad->nvsqnu_host      = mkarr(false, rad->nvsqnu     ->ncomp, rad->nvsqnu     ->size); 
  }
  else {
    rad->nvnu_surf_host   = gkyl_array_acquire(rad->nvnu_surf);
    rad->nvnu_host        = gkyl_array_acquire(rad->nvnu);
    rad->nvsqnu_surf_host = gkyl_array_acquire(rad->nvsqnu_surf);
    rad->nvsqnu_host      = gkyl_array_acquire(rad->nvsqnu);
  }

  // Allocate data and updaters for integrated moments.
  gk_species_moment_init(app, s, &rad->integ_moms, GKYL_F_MOMENT_M0M1M2PARM2PERP, true);
  int num_mom = rad->integ_moms.num_mom;
  if (app->use_gpu) {
    rad->red_integ_diag = gkyl_cu_malloc(sizeof(double[num_mom]));
    rad->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[num_mom]));
  } 
  else {
    rad->red_integ_diag = gkyl_malloc(sizeof(double[num_mom]));
    rad->red_integ_diag_global = gkyl_malloc(sizeof(double[num_mom]));
  }
  // allocate dynamic-vector to store all-reduced integrated moments 
  rad->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, num_mom);
  rad->is_first_integ_write_call = true;

  // Arrays for emissivity
  rad->emissivity_rhs = mkarr(app->use_gpu, s->basis.num_basis, s->local_ext.volume);
  rad->emissivity_denominator = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  // allocate M2 for emissivity
  gk_species_moment_init(app, s, &rad->m2, GKYL_F_MOMENT_M2, false);

  // Radiation updater
  struct gkyl_dg_rad_gyrokinetic_auxfields drag_inp = { .nvnu_surf = rad->nvnu_surf, .nvnu = rad->nvnu,
    .nvsqnu_surf = rad->nvsqnu_surf, .nvsqnu = rad->nvsqnu};
  rad->drag_slvr = gkyl_dg_updater_rad_gyrokinetic_new(&s->grid, 
    &app->basis, &s->basis, &s->local, &app->local, s->vel_map, &drag_inp, app->use_gpu);
}

// computes density for computation of total radiation drag and primitive moments
void
gk_species_radiation_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_rad_drag *rad, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  struct timespec wst = gkyl_wall_clock(); 

  // compute needed Maxwellian moments (n, u_par, T/m) (Jacobian factors already eliminated)
  gk_species_moment_calc(&rad->prim_moms, species->local, app->local, species->f);

  gkyl_array_set_offset(rad->m0, 1.0, rad->prim_moms.marr, 0*app->basis.num_basis);
  gkyl_array_set_offset(rad->vtsq, 1.0, rad->prim_moms.marr, 2*app->basis.num_basis);

  gkyl_array_clear(rad->nvnu_surf, 0.0);
  gkyl_array_clear(rad->nvnu, 0.0);
  gkyl_array_clear(rad->nvsqnu_surf, 0.0);
  gkyl_array_clear(rad->nvsqnu, 0.0);

  for (int i=0; i<rad->num_cross_collisions; ++i) {
    // compute needed moments
    if (rad->is_neut_species[i]) {
      gk_neut_species_moment_calc(&rad->moms[i], rad->collide_with_neut[i]->local,
        app->local, fin_neut[rad->collide_with_idx[i]]);
    }
    else {
      gk_species_moment_calc(&rad->moms[i], rad->collide_with[i]->local, app->local, fin[rad->collide_with_idx[i]]);
    }
    // divide out Jacobian from ion density before computation of final drag coefficient
    gkyl_dg_div_op_range(rad->moms[i].mem_geo, app->basis, 0, rad->moms[i].marr, 0,
      rad->moms[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
    gkyl_dg_calc_gk_rad_vars_nI_nu_advance(
      rad->calc_gk_rad_vars, &app->local, &species->local, 
      &rad->vnu_surf[i], &rad->vnu[i], &rad->vsqnu_surf[i], &rad->vsqnu[i], 
      rad->rad_fit_ne[i], rad->m0, rad->moms[i].marr, 
      rad->nvnu_surf, rad->nvnu, rad->nvsqnu_surf, rad->nvsqnu,
      rad->vtsq_min_per_species[i], rad->vtsq);					   
  }

  app->stat.species_rad_mom_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_species_radiation_integrated_moms(gkyl_gyrokinetic_app *app, struct gk_species *species,
	struct gk_rad_drag *rad, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  gkyl_array_clear(rad->emissivity_rhs, 0.0);
  gkyl_dg_updater_rad_gyrokinetic_advance(rad->drag_slvr, &species->local,
    species->f, species->cflrate, rad->emissivity_rhs);
  gk_species_moment_calc(&rad->integ_moms, species->local, app->local, rad->emissivity_rhs);
}

// computes emissivity
void
gk_species_radiation_emissivity(gkyl_gyrokinetic_app *app, struct gk_species *species,
  struct gk_rad_drag *rad, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  struct timespec wst = gkyl_wall_clock(); 

  // Compute needed Maxwellian moments (n, u_par, T/m) (Jacobian factors already eliminated).
  gk_species_moment_calc(&rad->prim_moms, species->local, app->local, species->f);

  gkyl_array_set_offset(rad->m0, 1.0, rad->prim_moms.marr, 0*app->basis.num_basis);
  gkyl_array_set_offset(rad->vtsq, 1.0, rad->prim_moms.marr, 2*app->basis.num_basis);

  // Calculate m2.
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
    gkyl_dg_div_op_range(rad->moms[i].mem_geo, app->basis, 0, rad->moms[i].marr, 0,
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

    gkyl_dg_mul_op(app->basis, 0, rad->emissivity_denominator, 0, rad->m0, 0, rad->moms[i].marr);
    gkyl_dg_mul_op(app->basis, 0, rad->emissivity_denominator, 0, rad->emissivity_denominator, 0, app->gk_geom->jacobgeo);
    gkyl_dg_div_op_range(rad->m2.mem_geo ,app->basis, 0, rad->emissivity[i], 0, rad->m2.marr, 0, rad->emissivity_denominator, &app->local);

    rad->emissivity[i] = gkyl_array_scale(rad->emissivity[i], -species->info.mass/2.0);
  }  
  app->stat.species_rad_mom_tm += gkyl_time_diff_now_sec(wst);
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
  
  app->stat.species_rad_tm += gkyl_time_diff_now_sec(wst);
}

// write functions
void
gk_species_radiation_write_drag(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (gks->rad.radiation_id == GKYL_GK_RADIATION) {
    struct timespec wst = gkyl_wall_clock();
    // Compute radiation drag coefficients
    const struct gkyl_array *fin_neut[app->num_neut_species];
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i) {
      fin[i] = app->species[i].f;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      fin_neut[i] = app->neut_species[i].f;
    }

    gk_species_radiation_moms(app, gks, &gks->rad, fin, fin_neut);
    app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);

    struct timespec wtm = gkyl_wall_clock();
    struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->basis.id      
      }
    );

    // Construct the file handles for vparallel and mu drag
    const char *fmt_nvnu_surf = "%s-%s_radiation_nvnu_surf_%d.gkyl";
    int sz_nvnu_surf = gkyl_calc_strlen(fmt_nvnu_surf, app->name, gks->info.name, frame);
    char fileNm_nvnu_surf[sz_nvnu_surf+1]; // ensures no buffer overflow
    snprintf(fileNm_nvnu_surf, sizeof fileNm_nvnu_surf, fmt_nvnu_surf, app->name, gks->info.name, frame);

    const char *fmt_nvnu = "%s-%s_radiation_nvnu_%d.gkyl";
    int sz_nvnu = gkyl_calc_strlen(fmt_nvnu, app->name, gks->info.name, frame);
    char fileNm_nvnu[sz_nvnu+1]; // ensures no buffer overflow
    snprintf(fileNm_nvnu, sizeof fileNm_nvnu, fmt_nvnu, app->name, gks->info.name, frame);

    const char *fmt_nvsqnu_surf = "%s-%s_radiation_nvsqnu_surf_%d.gkyl";
    int sz_nvsqnu_surf = gkyl_calc_strlen(fmt_nvsqnu_surf, app->name, gks->info.name, frame);
    char fileNm_nvsqnu_surf[sz_nvsqnu_surf+1]; // ensures no buffer overflow
    snprintf(fileNm_nvsqnu_surf, sizeof fileNm_nvsqnu_surf, fmt_nvsqnu_surf, app->name, gks->info.name, frame);

    const char *fmt_nvsqnu = "%s-%s_radiation_nvsqnu_%d.gkyl";
    int sz_nvsqnu = gkyl_calc_strlen(fmt_nvsqnu, app->name, gks->info.name, frame);
    char fileNm_nvsqnu[sz_nvsqnu+1]; // ensures no buffer overflow
    snprintf(fileNm_nvsqnu, sizeof fileNm_nvsqnu, fmt_nvsqnu, app->name, gks->info.name, frame);

    // copy data from device to host before writing it out
    if (app->use_gpu) {
      gkyl_array_copy(gks->rad.nvnu_surf_host, gks->rad.nvnu_surf);
      gkyl_array_copy(gks->rad.nvnu_host, gks->rad.nvnu);
      gkyl_array_copy(gks->rad.nvsqnu_surf_host, gks->rad.nvsqnu_surf);
      gkyl_array_copy(gks->rad.nvsqnu_host, gks->rad.nvsqnu);
    }

    gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->rad.nvnu_surf_host, fileNm_nvnu_surf);
    gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->rad.nvnu_host, fileNm_nvnu);
    gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->rad.nvsqnu_surf_host, fileNm_nvsqnu_surf);
    gkyl_comm_array_write(gks->comm, &gks->grid, &gks->local, mt, gks->rad.nvsqnu_host, fileNm_nvsqnu);
    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.n_diag_io += 4;

    gk_array_meta_release(mt);   
  }
}

void
gk_species_radiation_write_emissivity(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm, int frame)
{
  if (gks->rad.radiation_id == GKYL_GK_RADIATION) {
    struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->basis.id
      }
    );
    
    struct timespec wst = gkyl_wall_clock();
    const struct gkyl_array *fin_neut[app->num_neut_species];
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i) {
      fin[i] = app->species[i].f;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      fin_neut[i] = app->neut_species[i].f;
    }

    gk_species_radiation_emissivity(app, gks, &gks->rad, fin, fin_neut);
    double diag_calc_tm = gkyl_time_diff_now_sec(wst);

    int num_charged_species = 0;
    int num_neut_species = 0;
    for (int i=0; i<gks->rad.num_cross_collisions; i++) {
      struct timespec wtm = gkyl_wall_clock();
      // copy data from device to host before writing it out
      if (app->use_gpu)
        gkyl_array_copy(gks->rad.emissivity_host[i], gks->rad.emissivity[i]);

      // Construct the file handles for vparallel and mu drag
      const char *fmt_emissivity = "%s-%s_radiation_emissivity_%s_%d.gkyl";  
      if (gks->rad.is_neut_species[i]) {
        int sz_emissivity = gkyl_calc_strlen(fmt_emissivity, app->name, gks->info.name,
          app->neut_species[gks->rad.collide_with_idx[i]].info.name, frame);
        char fileNm_emissivity[sz_emissivity+1]; // ensures no buffer overflow
        snprintf(fileNm_emissivity, sizeof fileNm_emissivity, fmt_emissivity, app->name,
          gks->info.name, app->neut_species[gks->rad.collide_with_idx[i]].info.name, frame);

        gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->rad.emissivity_host[i], fileNm_emissivity);
        app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
        num_neut_species++;
      } 
      else {
        int sz_emissivity = gkyl_calc_strlen(fmt_emissivity, app->name, gks->info.name,
          app->species[gks->rad.collide_with_idx[i]].info.name, frame);
        char fileNm_emissivity[sz_emissivity+1]; // ensures no buffer overflow
        snprintf(fileNm_emissivity, sizeof fileNm_emissivity, fmt_emissivity, app->name,
          gks->info.name, app->species[gks->rad.collide_with_idx[i]].info.name, frame);

        gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gks->rad.emissivity_host[i], fileNm_emissivity);
        app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
        num_charged_species++;
      }  
      app->stat.n_diag_io += 1;
    }

    app->stat.species_diag_calc_tm += (diag_calc_tm * num_charged_species) / gks->rad.num_cross_collisions;
    app->stat.neut_species_diag_calc_tm += (diag_calc_tm * num_neut_species) / gks->rad.num_cross_collisions;

    gk_array_meta_release(mt); 
  }
}

void
gk_species_radiation_calc_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks, double tm)
{
  struct timespec wst = gkyl_wall_clock(); 
  if (gks->rad.radiation_id == GKYL_GK_RADIATION) {
    int vdim = app->vdim;
    int num_mom = gks->rad.integ_moms.num_mom;
    double avals_global[num_mom];

    // Compute radiation drag coefficients
    const struct gkyl_array *fin_neut[app->num_neut_species];
    const struct gkyl_array *fin[app->num_species];
    for (int i=0; i<app->num_species; ++i)  {
      fin[i] = app->species[i].f;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      fin_neut[i] = app->neut_species[i].f;
    }
    gk_species_radiation_moms(app, gks, &gks->rad, fin, fin_neut);
    gk_species_radiation_integrated_moms(app, gks, &gks->rad, fin, fin_neut);
  
    // reduce to compute sum over whole domain, append to diagnostics
    gkyl_array_reduce_range(gks->rad.red_integ_diag, gks->rad.integ_moms.marr, GKYL_SUM, &app->local);
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_mom, 
      gks->rad.red_integ_diag, gks->rad.red_integ_diag_global);
    if (app->use_gpu) {
      gkyl_cu_memcpy(avals_global, gks->rad.red_integ_diag_global, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      memcpy(avals_global, gks->rad.red_integ_diag_global, sizeof(double[num_mom]));
    }
    gkyl_dynvec_append(gks->rad.integ_diag, tm, avals_global);
  }
  app->stat.species_diag_calc_tm += gkyl_time_diff_now_sec(wst);  
}

void
gk_species_radiation_write_integrated_mom(gkyl_gyrokinetic_app* app, struct gk_species *gks)
{
  if (gks->rad.radiation_id == GKYL_GK_RADIATION) {
    struct timespec wst = gkyl_wall_clock();

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s_radiation_%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name, "integrated_moms");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name, "integrated_moms");

      if (gks->rad.is_first_integ_write_call) {
        gkyl_dynvec_write(gks->rad.integ_diag, fileNm);
        gks->rad.is_first_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(gks->rad.integ_diag, fileNm);
      }
    }
    gkyl_dynvec_clear(gks->rad.integ_diag);

    app->stat.species_diag_io_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_diag_io += 1;
  }
}
  
void 
gk_species_radiation_release(const struct gkyl_gyrokinetic_app *app, const struct gk_rad_drag *rad)
{
  for (int i=0; i<rad->num_cross_collisions; ++i) {
    gkyl_array_release(rad->emissivity[i]);
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
  gkyl_array_release(rad->nvnu_surf_host);
  gkyl_array_release(rad->nvnu_host);
  gkyl_array_release(rad->nvsqnu_surf_host);
  gkyl_array_release(rad->nvsqnu_host);    
  gkyl_dg_updater_rad_gyrokinetic_release(rad->drag_slvr);
}
