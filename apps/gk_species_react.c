#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_species_react_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, 
  struct gkyl_gyrokinetic_react inp, struct gk_react *react, bool all_gk)
{
  react->num_react = inp.num_react; 
  react->all_gk = all_gk;
  // initialize information about reactions from input struct
  for (int i=0; i<react->num_react; ++i) 
    react->react_type[i] = inp.react_type[i];
}

static double
gk_species_react_get_vtsq_min(struct gkyl_gyrokinetic_app *app, struct gk_species *s)
{
  double bmag_mid = app->bmag_ref;

  int vdim = app->vdim;
  double dv_min[vdim];
  gkyl_velocity_map_reduce_dv_range(s->vel_map, GKYL_MIN, dv_min, s->vel_map->local_vel);

  double tpar_min = (s->info.mass/6.0)*pow(dv_min[0],2);
  double tperp_min = vdim>1 ? (bmag_mid/3.0)*dv_min[1] : tpar_min;
  return (tpar_min + 2.0*tperp_min)/(3.0*s->info.mass);
}

static double
gk_neut_species_react_get_vtsq_min(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  double bmag_mid = app->bmag_ref;

  int vdim = app->vdim+1; // neutral species are 3v otherwise
  double dv_min[vdim];
  gkyl_velocity_map_reduce_dv_range(s->vel_map, GKYL_MIN, dv_min, s->vel_map->local_vel);

  double t_min = 0.0;
  for (int i=0; i<vdim; i++)
    t_min += (s->info.mass/6.0)*pow(dv_min[0],2);
  return t_min/(3.0*s->info.mass);
}

void 
gk_species_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_react *react)
{
  // distribution function which holds update for each reaction
  // form depend on react->type_self, e.g., for ionization and react->type_self == GKYL_SELF_ELC
  // react->f_react = n_elc*coeff_react*(2*fmax(n_elc, upar_donor, vtiz^2) - f_elc)
  react->f_react = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  react->f_react_other = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  react->proj_max = gkyl_proj_maxwellian_on_basis_new(&s->grid,
    &app->confBasis, &app->basis, app->basis.poly_order+1, s->vel_map, app->use_gpu);

  int vdim = app->vdim;

  for (int i=0; i<react->num_react; ++i) {
    react->react_id[i] = react->react_type[i].react_id;
    react->type_self[i] = react->react_type[i].type_self;

    // Fetch pointers to species objects
    react->species_elc[i] = gk_find_species(app, react->react_type[i].elc_nm);
    react->species_ion[i] = gk_find_species(app, react->react_type[i].ion_nm);

    // Fetch index of species for indexing arrays
    react->elc_idx[i] = gk_find_species_idx(app, react->react_type[i].elc_nm);
    react->ion_idx[i] = gk_find_species_idx(app, react->react_type[i].ion_nm);

    // Compute a minimum representable temperature based on the smallest dv in the grid.
    react->ion_vtsq_min[i] = gk_species_react_get_vtsq_min(app, gk_find_species(app, react->react_type[i].ion_nm));

    gk_species_moment_init(app, &app->species[react->elc_idx[i]], &react->moms_elc[i], "ThreeMoments", false);
    gk_species_moment_init(app, &app->species[react->ion_idx[i]], &react->moms_ion[i], "ThreeMoments", false);

    // If all the reacting species are gyrokinetic species, need to use 
    // gk methods to fetch pointers and indices, otherwise use gk_neut methods
    // to get the necessary neutral species information
    if (react->all_gk && gk_find_species(app, react->react_type[i].donor_nm)) {
      react->donor_idx[i] = gk_find_species_idx(app, react->react_type[i].donor_nm);
      gk_species_moment_init(app, &app->species[react->donor_idx[i]], &react->moms_donor[i], "ThreeMoments", false);
    }
    else if (gk_find_neut_species(app, react->react_type[i].donor_nm)) {
      react->donor_idx[i] = gk_find_neut_species_idx(app, react->react_type[i].donor_nm);
      gk_neut_species_moment_init(app, &app->neut_species[react->donor_idx[i]], &react->moms_donor[i], "FiveMoments");   
    }
    else if (gk_find_neut_species(app, react->react_type[i].partner_nm)) {

      react->partner_idx[i] = gk_find_neut_species_idx(app, react->react_type[i].partner_nm);

      // Compute a minimum representable temperature based on the smallest dv in the grid.
      react->neut_vtsq_min[i] = gk_neut_species_react_get_vtsq_min(app, &app->neut_species[react->partner_idx[i]]);

      gk_neut_species_moment_init(app, &app->neut_species[react->partner_idx[i]], &react->moms_partner[i], "FiveMoments");

      react->prim_vars_cxi[i] = mkarr(app->use_gpu, (2+app->vdim)*app->confBasis.num_basis, app->local_ext.volume);
      react->prim_vars_cxn[i] = mkarr(app->use_gpu, (2+app->vdim)*app->confBasis.num_basis, app->local_ext.volume);
    }

    react->coeff_react[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    if(app->use_gpu)
      react->coeff_react_host[i] = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    else
      react->coeff_react_host[i] = react->coeff_react[i];

    react->vt_sq_iz1[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->vt_sq_iz2[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_elc[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_ion[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_donor[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_partner[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_mod[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->prim_vars[i] = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
    react->prim_vars_donor[i] = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
    react->prim_vars_proj_inp[i] = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    react->prim_vars_se_proj_inp[i] = mkarr(app->use_gpu, 3*app->confBasis.num_basis, app->local_ext.volume);
    
    if (react->react_id[i] == GKYL_REACT_IZ) {
      struct gkyl_dg_iz_inp iz_inp = {
        .grid = &s->grid,
        .cbasis = &app->confBasis,
        .pbasis = &app->basis,
        .conf_rng = &app->local,
        .conf_rng_ext = &app->local_ext,
        .phase_rng = &s->local,
        .mass_ion = react->react_type[i].ion_mass,
        .type_ion = react->react_type[i].ion_id,
        .charge_state = react->react_type[i].charge_state,
        .type_self = react->type_self[i],
        .all_gk = react->all_gk,
      };
      react->iz[i] = gkyl_dg_iz_new(&iz_inp, app->use_gpu);
    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      struct gkyl_dg_recomb_inp recomb_inp = {
        .grid = &s->grid,
        .cbasis = &app->confBasis,
        .pbasis = &app->basis,
        .conf_rng = &app->local,
        .conf_rng_ext = &app->local_ext,
        .phase_rng = &s->local,
        .mass_self = s->info.mass,
        .type_ion = react->react_type[i].ion_id,
        .charge_state = react->react_type[i].charge_state,
        .type_self = react->type_self[i],
        .all_gk = react->all_gk,
      };
      react->recomb[i] = gkyl_dg_recomb_new(&recomb_inp, app->use_gpu);
    }
    else if (react->react_id[i] == GKYL_REACT_CX) {
      struct gkyl_dg_cx_inp cx_inp = {
        .grid = &s->grid,
	.cbasis = &app->confBasis,
	.pbasis_gk = &app->basis,
	.pbasis_vl = &app->neut_basis,
	.conf_rng = &app->local,
	.conf_rng_ext = &app->local_ext,
	.phase_rng = &s->local,
	.mass_ion = react->react_type[i].ion_mass,
	.mass_neut = react->react_type[i].partner_mass,
	.type_ion = react->react_type[i].ion_id,
      };
      react->cx[i] = gkyl_dg_cx_new(&cx_inp, app->use_gpu);
    }
  }
}

// computes reaction coefficients
void
gk_species_react_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_react *react, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  for (int i=0; i<react->num_react; ++i) {
    if (react->react_id[i] == GKYL_REACT_IZ) {
      // compute needed moments
      gk_species_moment_calc(&react->moms_elc[i], app->species[react->elc_idx[i]].local,
        app->local, fin[react->elc_idx[i]]);
      for (int j=0; j<react->moms_elc[i].num_mom; ++j) {
        gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, j, react->moms_elc[i].marr, j,
          react->moms_elc[i].marr, 0, app->gk_geom->jacobgeo, &app->local);   
      }  
      gkyl_array_set_range(react->m0_elc[i], 1.0, react->moms_elc[i].marr, &app->local);

      if (react->all_gk) {
        gk_species_moment_calc(&react->moms_donor[i], app->species[react->donor_idx[i]].local,
          app->local, fin[react->donor_idx[i]]);
      }
      else {
        gk_neut_species_moment_calc(&react->moms_donor[i], app->neut_species[react->donor_idx[i]].local,
          app->local, fin_neut[react->donor_idx[i]]);
      }
      for (int j=0; j<react->moms_donor[i].num_mom; ++j) {
        gkyl_dg_div_op_range(react->moms_donor[i].mem_geo, app->confBasis, j, react->moms_donor[i].marr, j,
          react->moms_donor[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      }
      gkyl_array_set_range(react->m0_donor[i], 1.0, react->moms_donor[i].marr, &app->local);

      // compute ionization reaction rate
      gkyl_dg_iz_coll(react->iz[i], react->moms_elc[i].marr, react->moms_donor[i].marr,
    	app->gk_geom->b_i, react->prim_vars[i], react->prim_vars_donor[i],
        react->vt_sq_iz1[i], react->vt_sq_iz2[i], react->coeff_react[i], 0);
    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      // compute needed moments
      gk_species_moment_calc(&react->moms_elc[i], app->species[react->elc_idx[i]].local,
        app->local, fin[react->elc_idx[i]]);
      for (int j=0; j<react->moms_elc[i].num_mom; ++j) {
        gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, j, react->moms_elc[i].marr, j,
          react->moms_elc[i].marr, 0, app->gk_geom->jacobgeo, &app->local);   
      }  
      gkyl_array_set_range(react->m0_elc[i], 1.0, react->moms_elc[i].marr, &app->local);

      gk_species_moment_calc(&react->moms_ion[i], app->species[react->ion_idx[i]].local,
        app->local, fin[react->ion_idx[i]]);
      for (int j=0; j<react->moms_ion[i].num_mom; ++j) {
        gkyl_dg_div_op_range(react->moms_ion[i].mem_geo, app->confBasis, j, react->moms_ion[i].marr, j,
          react->moms_ion[i].marr, 0, app->gk_geom->jacobgeo, &app->local);   
      }  
      gkyl_array_set_range(react->m0_ion[i], 1.0, react->moms_ion[i].marr, &app->local);
      
      // compute recombination reaction rate
      gkyl_dg_recomb_coll(react->recomb[i], react->moms_elc[i].marr, react->moms_ion[i].marr,
        app->gk_geom->b_i, react->prim_vars[i], react->coeff_react[i], 0);
    }
    else if (react->react_id[i] == GKYL_REACT_CX) {
      // calc moms_ion, moms_neut
      gk_species_moment_calc(&react->moms_ion[i], app->species[react->ion_idx[i]].local,
        app->local, fin[react->ion_idx[i]]);
      for (int j=0; j<react->moms_ion[i].num_mom; ++j) {
        gkyl_dg_div_op_range(react->moms_ion[i].mem_geo, app->confBasis, j, react->moms_ion[i].marr, j,
          react->moms_ion[i].marr, 0, app->gk_geom->jacobgeo, &app->local);   
      }
      gkyl_array_set_range(react->m0_ion[i], 1.0, react->moms_ion[i].marr, &app->local);

      gk_neut_species_moment_calc(&react->moms_partner[i], app->neut_species[react->partner_idx[i]].local,
        app->local, fin_neut[react->partner_idx[i]]);
      
      for (int j=0; j<react->moms_partner[i].num_mom; ++j) {
        gkyl_dg_div_op_range(react->moms_partner[i].mem_geo, app->confBasis, j, react->moms_partner[i].marr, j,
          react->moms_partner[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      }
      gkyl_array_set_range(react->m0_partner[i], 1.0, react->moms_partner[i].marr, &app->local);

      // prim_vars_neut_gk is returned to prim_vars[i] here.
      gkyl_dg_cx_coll(react->cx[i], react->ion_vtsq_min[i], react->neut_vtsq_min[i],
        react->moms_ion[i].marr, react->moms_partner[i].marr, app->gk_geom->bcart, react->prim_vars_cxi[i],
        react->prim_vars_cxn[i], react->prim_vars[i], react->coeff_react[i], 0);
    }
  }
}

// updates the reaction terms in the rhs
void
gk_species_react_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<react->num_react; ++i) {
    gkyl_array_clear(react->f_react, 0.0);

    if (react->react_id[i] == GKYL_REACT_IZ) {
      if (react->type_self[i] == GKYL_SELF_ELC) {
	// electron update is n_elc*coeff_react*(fmax1(n_elc, upar_elc, vtiz1^2)
	//   + fmax2(n_elc, upar_donor, vtiz2^2) - f_elc)

	// primary elc
        gkyl_array_set_offset(react->prim_vars_proj_inp[i], 1.0, react->m0_elc[i], 0*app->confBasis.num_basis);
        gkyl_array_set_offset(react->prim_vars_proj_inp[i], 1.0, react->prim_vars[i], 1*app->confBasis.num_basis);        
        gkyl_array_set_offset(react->prim_vars_proj_inp[i], 1.0, react->vt_sq_iz1[i], 2*app->confBasis.num_basis);

        gkyl_proj_gkmaxwellian_on_basis_prim_mom(react->proj_max, &s->local, &app->local,
          react->prim_vars_proj_inp[i],
          app->gk_geom->bmag, app->gk_geom->jacobtot, s->info.mass, react->f_react);
	gkyl_array_scale_by_cell(react->f_react, s->vel_map->jacobvel);
				 
        // scale to correct m0
        gk_species_moment_calc(&s->m0, s->local_ext, app->local_ext, react->f_react);
        gkyl_dg_div_op_range(s->m0.mem_geo, app->confBasis, 0, react->m0_mod[i], 0,
          react->m0_elc[i], 0, s->m0.marr, &app->local);
        gkyl_dg_mul_op_range(app->confBasis, 0, react->m0_mod[i], 0, react->m0_mod[i], 0, app->gk_geom->jacobgeo, &app->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react, 
          react->m0_mod[i], react->f_react, &app->local_ext, &s->local_ext);
       
	// secondary elc (se)
	gkyl_array_set_offset(react->prim_vars_se_proj_inp[i], 1.0, react->m0_elc[i], 0*app->confBasis.num_basis);
	gkyl_array_set_offset(react->prim_vars_se_proj_inp[i], 1.0, react->prim_vars_donor[i], 1*app->confBasis.num_basis);
        gkyl_array_set_offset(react->prim_vars_se_proj_inp[i], 1.0, react->vt_sq_iz2[i], 2*app->confBasis.num_basis);

      	gkyl_array_set_offset(react->prim_vars_donor[i], 1.0, react->vt_sq_iz2[i], 1*app->confBasis.num_basis);
        gkyl_proj_gkmaxwellian_on_basis_prim_mom(react->proj_max, &s->local, &app->local,
          react->prim_vars_se_proj_inp[i],
          app->gk_geom->bmag, app->gk_geom->jacobtot, s->info.mass, react->f_react_other);
	gkyl_array_scale_by_cell(react->f_react_other, s->vel_map->jacobvel);

        // scale to correct m0
        gk_species_moment_calc(&s->m0, s->local_ext, app->local_ext, react->f_react_other);
        gkyl_dg_div_op_range(s->m0.mem_geo, app->confBasis, 0, react->m0_mod[i], 0,
          react->m0_elc[i], 0, s->m0.marr, &app->local);

	gkyl_dg_mul_op_range(app->confBasis, 0, react->m0_mod[i], 0, react->m0_mod[i], 0, app->gk_geom->jacobgeo, &app->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react_other,
          react->m0_mod[i], react->f_react_other, &app->local_ext, &s->local_ext);	

        gkyl_array_accumulate(react->f_react, 1.0, react->f_react_other);
        gkyl_array_accumulate(react->f_react, -1.0, fin);
	
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_donor[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, 1.0, react->f_react);
      }
      else if (react->type_self[i] == GKYL_SELF_ION) {
        gkyl_array_set_offset(react->prim_vars_proj_inp[i], 1.0, react->m0_donor[i], 0*app->confBasis.num_basis);
        gkyl_array_set_offset(react->prim_vars_proj_inp[i], 1.0, react->prim_vars_donor[i], 1*app->confBasis.num_basis);        
        gkyl_proj_gkmaxwellian_on_basis_prim_mom(react->proj_max, &s->local, &app->local,
	  react->prim_vars_proj_inp[i], 
          app->gk_geom->bmag, app->gk_geom->jacobtot, s->info.mass, react->f_react);
	gkyl_array_scale_by_cell(react->f_react, s->vel_map->jacobvel);

        // scale to correct m0
        gk_species_moment_calc(&s->m0, s->local_ext, app->local_ext, react->f_react);
        gkyl_dg_div_op_range(s->m0.mem_geo, app->confBasis, 0, react->m0_mod[i], 0,
          react->m0_donor[i], 0, s->m0.marr, &app->local);
        gkyl_dg_mul_op_range(app->confBasis, 0, react->m0_mod[i], 0, react->m0_mod[i], 0, app->gk_geom->jacobgeo, &app->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react, 
          react->m0_mod[i], react->f_react, &app->local_ext, &s->local_ext);

        // ion update is n_elc*coeff_react*fmax(n_donor, upar_donor, vt_donor^2)
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_elc[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, 1.0, react->f_react);
      }
      else {
        // donor update is -n_elc*coeff_react*f_donor
        gkyl_array_set(react->f_react, 1.0, fin);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_elc[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, -1.0, react->f_react);
      }
    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      if (react->type_self[i] == GKYL_SELF_ELC) {
        // update is -n_ion*coeff_react*f_elc
        gkyl_array_set(react->f_react, 1.0, fin);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_ion[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, -1.0, react->f_react);
      }
      else if (react->type_self[i] == GKYL_SELF_ION) {
        // update is -n_elc*coeff_react*f_ion
        gkyl_array_set(react->f_react, 1.0, fin);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_elc[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, -1.0, react->f_react);
      }
      else {
        gkyl_proj_gkmaxwellian_on_basis_lab_mom(react->proj_max, &s->local, &app->local,
          react->moms_ion[i].marr, app->gk_geom->bmag, app->gk_geom->jacobtot, s->info.mass, react->f_react);
	gkyl_array_scale_by_cell(react->f_react, s->vel_map->jacobvel);
        // scale to correct m0
        gk_species_moment_calc(&s->m0, s->local_ext, app->local_ext, react->f_react); 
        gkyl_dg_div_op_range(s->m0.mem_geo, app->confBasis, 0, react->m0_mod[i], 0,
          react->m0_ion[i], 0, s->m0.marr, &app->local);
        gkyl_dg_mul_op_range(app->confBasis, 0, react->m0_mod[i], 0, react->m0_mod[i], 0, app->gk_geom->jacobgeo, &app->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react, 
          react->m0_mod[i], react->f_react, &app->local_ext, &s->local_ext);

        // receiver update is n_elc*coeff_react*fmax(n_ion, upar_ion, vt_ion^2)
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
          react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
          react->m0_elc[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, 1.0, react->f_react);
      }
    }
    else if (react->react_id[i] == GKYL_REACT_CX) {
      // calculate RHS for ions
      // here prim_vars[i] is prim_vars_neut_gk
      gkyl_array_set_offset(react->prim_vars_proj_inp[i], 1.0, react->m0_partner[i], 0*app->confBasis.num_basis);
        gkyl_array_set_offset(react->prim_vars_proj_inp[i], 1.0, react->prim_vars[i], 1*app->confBasis.num_basis); 
      gkyl_proj_gkmaxwellian_on_basis_prim_mom(react->proj_max, &s->local, &app->local,
        react->prim_vars_proj_inp[i],
	app->gk_geom->bmag, app->gk_geom->jacobtot, s->info.mass, react->f_react);
      gkyl_array_scale_by_cell(react->f_react, s->vel_map->jacobvel);

      // scale to correct m0 and multiply f
      gk_species_moment_calc(&s->m0, s->local_ext, app->local_ext, react->f_react);
      gkyl_dg_div_op_range(s->m0.mem_geo, app->confBasis, 0, react->m0_mod[i], 0,
        react->m0_partner[i], 0, s->m0.marr, &app->local);
      gkyl_dg_mul_op_range(app->confBasis, 0, react->m0_mod[i], 0, react->m0_mod[i], 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react, 
        react->m0_mod[i], react->f_react, &app->local_ext, &s->local_ext);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
        react->m0_ion[i], react->f_react, &app->local, &s->local);
      
      // ion update is coeff_react*(n_ion*f_neut - n_neut*f_ion)
      gkyl_array_set(react->f_react_other, 1.0, fin);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react_other,
        react->m0_partner[i], react->f_react_other, &app->local, &s->local);
      gkyl_array_accumulate(react->f_react, -1.0, react->f_react_other);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
        react->coeff_react[i], react->f_react, &app->local, &s->local);
      gkyl_array_accumulate(rhs, 1.0, react->f_react);      
      
    }
  }

  app->stat.species_coll_tm += gkyl_time_diff_now_sec(wst);
}

// write functions
void
gk_species_react_write(gkyl_gyrokinetic_app* app, struct gk_species *gks, struct gk_react *gkr,
  int ridx, double tm, int frame)
{
  if (gkr->type_self[ridx] == GKYL_SELF_ION) {
    struct timespec wst = gkyl_wall_clock();
    struct timespec wtm; 
    struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->confBasis.id
      }
    );

    // Compute reaction rate
    const struct gkyl_array *fin[app->num_species];
    const struct gkyl_array *fin_neut[app->num_neut_species];
    for (int i=0; i<app->num_species; ++i) 
      fin[i] = app->species[i].f;
    for (int i=0; i<app->num_neut_species; ++i)
      fin_neut[i] = app->neut_species[i].f;
    gk_species_react_cross_moms(app, gks, gkr, fin, fin_neut);
    
    if (app->use_gpu) {
      gkyl_array_copy(gkr->coeff_react_host[ridx], gkr->coeff_react[ridx]);
    }
    
    if (gkr->react_id[ridx] == GKYL_REACT_IZ) {
      const char *fmt = "%s-%s_%s_%s_iz_react_%d.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name,
        gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].donor_nm, frame);
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
        gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].donor_nm, frame);
  
      wtm = gkyl_wall_clock();
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gkr->coeff_react_host[ridx], fileNm);
    }
   if (gkr->react_id[ridx] == GKYL_REACT_RECOMB) {
     const char *fmt = "%s-%s_%s_%s_recomb_react_%d.gkyl";
     int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name,
       gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].recvr_nm, frame);
     char fileNm[sz+1]; // ensures no buffer overflow
     snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
       gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].recvr_nm, frame);
  
     wtm = gkyl_wall_clock();
     gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gkr->coeff_react_host[ridx], fileNm);
   }
   if (gkr->react_id[ridx] == GKYL_REACT_CX) {
     const char *fmt = "%s-%s_%s_cx_react_%d.gkyl";
     int sz = gkyl_calc_strlen(fmt, app->name, gks->info.name,
       gkr->react_type[ridx].partner_nm, frame);
     char fileNm[sz+1]; // ensures no buffer overflow
     snprintf(fileNm, sizeof fileNm, fmt, app->name, gks->info.name,
       gkr->react_type[ridx].partner_nm, frame);
     
     wtm = gkyl_wall_clock();
     gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, gkr->coeff_react_host[ridx], fileNm);
   }
   app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
   app->stat.nio += 1;
   
   gk_array_meta_release(mt); 
   app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
   app->stat.ndiag += 1;
  }
}

void 
gk_species_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react)
{
  gkyl_array_release(react->f_react);
  gkyl_array_release(react->f_react_other);
  gkyl_proj_maxwellian_on_basis_release(react->proj_max);
  for (int i=0; i<react->num_react; ++i) {
    gk_species_moment_release(app, &react->moms_elc[i]);
    gk_species_moment_release(app, &react->moms_ion[i]);
    if (react->all_gk && gk_find_species(app, react->react_type[i].donor_nm))
      gk_species_moment_release(app, &react->moms_donor[i]);
    else if (gk_find_neut_species(app, react->react_type[i].donor_nm))
      gk_neut_species_moment_release(app, &react->moms_donor[i]);
    else if (gk_find_neut_species(app, react->react_type[i].partner_nm))
      gk_neut_species_moment_release(app, &react->moms_partner[i]);

    gkyl_array_release(react->coeff_react[i]);
    gkyl_array_release(react->vt_sq_iz1[i]);
    gkyl_array_release(react->vt_sq_iz2[i]);
    gkyl_array_release(react->m0_elc[i]);
    gkyl_array_release(react->m0_ion[i]);
    gkyl_array_release(react->m0_donor[i]);
    gkyl_array_release(react->m0_partner[i]);
    gkyl_array_release(react->m0_mod[i]);
    gkyl_array_release(react->prim_vars[i]);
    gkyl_array_release(react->prim_vars_donor[i]); 
    gkyl_array_release(react->prim_vars_cxi[i]);
    gkyl_array_release(react->prim_vars_cxn[i]);
    gkyl_array_release(react->prim_vars_proj_inp[i]); 
    gkyl_array_release(react->prim_vars_se_proj_inp[i]); 

    if (app->use_gpu)
      gkyl_array_release(react->coeff_react_host[i]);

    if (react->react_id[i] == GKYL_REACT_IZ) 
      gkyl_dg_iz_release(react->iz[i]);
    else if (react->react_id[i] == GKYL_REACT_RECOMB)  
      gkyl_dg_recomb_release(react->recomb[i]);
    else if (react->react_id[i] == GKYL_REACT_CX) { 
      gkyl_dg_cx_release(react->cx[i]);
      gkyl_array_release(react->prim_vars_cxi[i]);
      gkyl_array_release(react->prim_vars_cxn[i]);
    }
  }
}
