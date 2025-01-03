#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_neut_species_react_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, 
  struct gkyl_gyrokinetic_react inp, struct gk_react *react)
{
  react->num_react = inp.num_react; 
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
gk_neut_species_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s, struct gk_react *react)
{
  // distribution function which holds update for each reaction
  // form depend on react->type_self, e.g., for recombination and react->type_self == GKYL_SELF_RECVR
  // react->f_react = n_elc*coeff_react*fmax(n_ion, upar_ion b_i, vt_ion^2)
  react->f_react = mkarr(app->use_gpu, app->neut_basis.num_basis, s->local_ext.volume);
  react->f_react_other = mkarr(app->use_gpu, app->neut_basis.num_basis, s->local_ext.volume);
  react->proj_max = gkyl_proj_maxwellian_on_basis_new(&s->grid,
    &app->confBasis, &app->neut_basis, app->neut_basis.poly_order+1, s->vel_map, app->use_gpu);  

  int vdim = app->vdim+1; // neutral species are 3v otherwise
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

    if (gk_find_neut_species(app, react->react_type[i].donor_nm)) {
      react->donor_idx[i] = gk_find_neut_species_idx(app, react->react_type[i].donor_nm);
      gk_neut_species_moment_init(app, &app->neut_species[react->donor_idx[i]], &react->moms_donor[i], "FiveMoments");   
    }
    if (gk_find_neut_species(app, react->react_type[i].partner_nm)) {
      react->partner_idx[i] = gk_find_neut_species_idx(app, react->react_type[i].partner_nm);

      // Compute a minimum representable temperature based on the smallest dv in the grid.
      react->neut_vtsq_min[i] = gk_neut_species_react_get_vtsq_min(app, &app->neut_species[react->partner_idx[i]]);

      gk_neut_species_moment_init(app, &app->neut_species[react->partner_idx[i]], &react->moms_partner[i], "FiveMoments");

      react->prim_vars_cxi[i] = mkarr(app->use_gpu, (2+app->vdim)*app->confBasis.num_basis, app->local_ext.volume);
      react->prim_vars_cxn[i] = mkarr(app->use_gpu, (2+app->vdim)*app->confBasis.num_basis, app->local_ext.volume);
    }
    react->coeff_react[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->vt_sq_iz1[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->vt_sq_iz2[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_elc[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_ion[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_partner[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_mod[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->prim_vars[i] = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
    react->prim_vars_donor[i] = mkarr(app->use_gpu, (vdim+1)*app->confBasis.num_basis, app->local_ext.volume);
    if (react->react_id[i] == GKYL_REACT_IZ) {
      struct gkyl_dg_iz_inp iz_inp = {
        .grid = &s->grid, 
        .cbasis = &app->confBasis, 
        .pbasis = &app->neut_basis, 
        .conf_rng = &app->local, 
        .conf_rng_ext = &app->local_ext,
        .phase_rng = &s->local, 
        .mass_ion = react->react_type[i].ion_mass, 
        .type_ion = react->react_type[i].ion_id, 
        .charge_state = react->react_type[i].charge_state, 
        .type_self = react->type_self[i], 
        .all_gk = false, 
      };
      react->iz[i] = gkyl_dg_iz_new(&iz_inp, app->use_gpu); 
    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      struct gkyl_dg_recomb_inp recomb_inp = {
        .grid = &s->grid,
        .cbasis = &app->confBasis, 
        .pbasis = &app->neut_basis, 
        .conf_rng = &app->local, 
        .conf_rng_ext = &app->local_ext,
        .phase_rng = &s->local, 
        .mass_self = s->info.mass, 
        .type_ion = react->react_type[i].ion_id, 
        .charge_state = react->react_type[i].charge_state, 
        .type_self = react->type_self[i], 
        .all_gk = false, 
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
gk_neut_species_react_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
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
      
      gk_neut_species_moment_calc(&react->moms_donor[i], app->neut_species[react->donor_idx[i]].local,
        app->local, fin_neut[react->donor_idx[i]]);
      for (int j=0; j<react->moms_donor[i].num_mom; ++j) {
        gkyl_dg_div_op_range(react->moms_donor[i].mem_geo, app->confBasis, j, react->moms_donor[i].marr, j,
          react->moms_donor[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      }

      // compute ionization reaction rate
      gkyl_dg_iz_coll(react->iz[i], react->moms_elc[i].marr, react->moms_donor[i].marr,
    	app->gk_geom->b_i, react->prim_vars[i], react->prim_vars_donor[i],
        react->vt_sq_iz1[i], react->vt_sq_iz2[i], react->coeff_react[i], 0);
    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      // compute needed moments
      gk_species_moment_calc(&react->moms_elc[i], app->species[react->elc_idx[i]].local,
        app->local, fin[react->elc_idx[i]]);
      gk_species_moment_calc(&react->moms_ion[i], app->species[react->ion_idx[i]].local,
        app->local, fin[react->ion_idx[i]]);

      for (int j=0; j<react->moms_elc[i].num_mom; ++j) {
        gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, j, react->moms_elc[i].marr, j,
          react->moms_elc[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
        gkyl_dg_div_op_range(react->moms_ion[i].mem_geo, app->confBasis, j, react->moms_ion[i].marr, j,
          react->moms_ion[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      }

      gkyl_array_set_range(react->m0_elc[i], 1.0, react->moms_elc[i].marr, &app->local);
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
gk_neut_species_react_rhs(gkyl_gyrokinetic_app *app, const struct gk_neut_species *s,
  struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  for (int i=0; i<react->num_react; ++i) {
    gkyl_array_clear(react->f_react, 0.0);
    // neutral species can only be the donor or the receiver in reaction
    if (react->react_id[i] == GKYL_REACT_IZ) {
      // donor update is -n_elc*coeff_react*f_donor
      gkyl_array_set(react->f_react, 1.0, fin);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, react->f_react, 
          react->coeff_react[i], react->f_react, &app->local, &s->local);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, react->f_react, 
          react->m0_elc[i], react->f_react, &app->local, &s->local);  
      gkyl_array_accumulate(rhs, -1.0, react->f_react);    

    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      gkyl_proj_maxwellian_on_basis_prim_mom(react->proj_max, &s->local, &app->local,
        react->moms_ion[i].marr, react->prim_vars[i], react->f_react);
      gk_neut_species_moment_calc(&s->m0, s->local_ext, app->local_ext, react->f_react);
      gkyl_dg_div_op_range(s->m0.mem_geo, app->confBasis, 0, react->m0_mod[i], 0,
        react->m0_ion[i], 0, s->m0.marr, &app->local);
      gkyl_dg_mul_op_range(app->confBasis, 0, react->m0_mod[i], 0, react->m0_mod[i], 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
        react->m0_mod[i], react->f_react, &app->local_ext, &s->local_ext);

      // receiver update is n_elc*coeff_react*fmax(n_ion, upar_ion, vt_ion^2)
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, react->f_react,
          react->coeff_react[i], react->f_react, &app->local, &s->local);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, react->f_react,
          react->m0_elc[i], react->f_react, &app->local, &s->local);
      // multiply by configuration space Jacobian before final accumulation
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, react->f_react,
          app->gk_geom->jacobgeo, react->f_react, &app->local_ext, &s->local_ext);
      gkyl_array_accumulate(rhs, 1.0, react->f_react);
    }
    else if (react->react_id[i] == GKYL_REACT_CX) {
      // here prim_vars[i] is prim_vars_neut_gk
      gkyl_proj_maxwellian_on_basis_prim_mom(react->proj_max, &s->local, &app->local,
        react->moms_ion[i].marr, react->prim_vars_cxi[i], react->f_react);

      // scale to correct m0 and multiply f
      gk_neut_species_moment_calc(&s->m0, s->local_ext, app->local_ext, react->f_react);
      gkyl_dg_div_op_range(s->m0.mem_geo, app->confBasis, 0, react->m0_mod[i], 0,
        react->m0_ion[i], 0, s->m0.marr, &app->local);
      gkyl_dg_mul_op_range(app->confBasis, 0, react->m0_mod[i], 0, react->m0_mod[i], 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, react->f_react,
        react->m0_mod[i], react->f_react, &app->local_ext, &s->local_ext);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, react->f_react,
        react->m0_partner[i], react->f_react, &app->local, &s->local);
      
      // neut update is coeff_react*(n_neut*f_ion - n_ion*f_neut)
      gkyl_array_set(react->f_react_other, 1.0, fin);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, react->f_react_other,
        react->m0_ion[i], react->f_react_other, &app->local, &s->local);
      gkyl_array_accumulate(react->f_react, -1.0, react->f_react_other);
      gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->neut_basis, react->f_react,
        react->coeff_react[i], react->f_react, &app->local, &s->local);
      gkyl_array_accumulate(rhs, 1.0, react->f_react);
    }
  }
}

void 
gk_neut_species_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react)
{
  gkyl_array_release(react->f_react);
  gkyl_proj_maxwellian_on_basis_release(react->proj_max);
  for (int i=0; i<react->num_react; ++i) {
    gk_species_moment_release(app, &react->moms_elc[i]);
    gk_species_moment_release(app, &react->moms_ion[i]);
    if (gk_find_neut_species(app, react->react_type[i].donor_nm))
      gk_neut_species_moment_release(app, &react->moms_donor[i]);
    if (gk_find_neut_species(app, react->react_type[i].partner_nm))
      gk_neut_species_moment_release(app, &react->moms_partner[i]);

    gkyl_array_release(react->coeff_react[i]);
    gkyl_array_release(react->vt_sq_iz1[i]);
    gkyl_array_release(react->vt_sq_iz2[i]);
    gkyl_array_release(react->m0_elc[i]);
    gkyl_array_release(react->m0_ion[i]);
    gkyl_array_release(react->m0_partner[i]);
    gkyl_array_release(react->m0_mod[i]);
    gkyl_array_release(react->prim_vars[i]);
    gkyl_array_release(react->prim_vars_donor[i]); 
    if (react->react_id[i] == GKYL_REACT_IZ) 
      gkyl_dg_iz_release(react->iz[i]);
    else if (react->react_id[i] == GKYL_REACT_RECOMB)  
      gkyl_dg_recomb_release(react->recomb[i]);
    else if (react->react_id[i] == GKYL_REACT_CX)
      gkyl_dg_cx_release(react->cx[i]);
  }
}
