#include <assert.h>
#include <gkyl_array_rio.h>
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

void 
gk_species_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_species *s, struct gk_react *react)
{
  // distribution function which holds update for each reaction
  // form depend on react->type_self, e.g., for ionization and react->type_self == GKYL_SELF_ELC
  // react->f_react = n_elc*coeff_react*(2*fmax(n_elc, upar_donor, vtiz^2) - f_elc)
  react->f_react = mkarr(app->use_gpu, app->basis.num_basis, s->local_ext.volume);
  react->proj_max = gkyl_proj_maxwellian_on_basis_new(&s->grid,
    &app->confBasis, &app->basis, app->basis.poly_order+1, app->use_gpu);

  for (int i=0; i<react->num_react; ++i) {
    react->react_id[i] = react->react_type[i].react_id;
    react->type_self[i] = react->react_type[i].type_self;
    // Fetch pointers to species objects
    react->species_elc[i] = gk_find_species(app, react->react_type[i].elc_nm);
    react->species_ion[i] = gk_find_species(app, react->react_type[i].ion_nm);
    // Fetch index of species for indexing arrays
    react->elc_idx[i] = gk_find_species_idx(app, react->react_type[i].elc_nm);
    react->ion_idx[i] = gk_find_species_idx(app, react->react_type[i].ion_nm);

    gk_species_moment_init(app, &app->species[react->elc_idx[i]], &react->moms_elc[i], "ThreeMoments");
    gk_species_moment_init(app, &app->species[react->ion_idx[i]], &react->moms_ion[i], "ThreeMoments");

    // If all the reacting species are gyrokinetic species, need to use 
    // gk methods to fetch pointers and indices, otherwise use gk_neut methods
    // to get the necessary neutral species information
    if (react->all_gk && gk_find_species(app, react->react_type[i].donor_nm)) {
      react->donor_idx[i] = gk_find_species_idx(app, react->react_type[i].donor_nm);
      gk_species_moment_init(app, &app->species[react->donor_idx[i]], &react->moms_donor[i], "ThreeMoments");
    }
    else if (gk_find_neut_species(app, react->react_type[i].donor_nm)) {
      react->donor_idx[i] = gk_find_neut_species_idx(app, react->react_type[i].donor_nm);
      gk_neut_species_moment_init(app, &app->neut_species[react->donor_idx[i]], &react->moms_donor[i], "FiveMoments");   
    }

    react->coeff_react[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->coeff_react_host[i] = react->coeff_react[i];
    react->vt_sq_iz[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_elc[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_ion[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->m0_donor[i] = mkarr(app->use_gpu, app->confBasis.num_basis, app->local_ext.volume);
    react->prim_vars[i] = mkarr(app->use_gpu, 2*app->confBasis.num_basis, app->local_ext.volume);
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
        .base = ".",
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
        .base = ".",
      };
      react->recomb[i] = gkyl_dg_recomb_new(&recomb_inp, app->use_gpu);
    }
  }
}

// computes reaction coefficients
void
gk_species_react_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_species *species,
  struct gk_react *react, const struct gkyl_array *f_self, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  for (int i=0; i<react->num_react; ++i) {
    if (react->react_id[i] == GKYL_REACT_IZ) {
      // compute needed moments
      gk_species_moment_calc(&react->moms_elc[i], app->species[react->elc_idx[i]].local,
        app->local, fin[react->elc_idx[i]]);
      if (react->all_gk)
        gk_species_moment_calc(&react->moms_donor[i], app->species[react->donor_idx[i]].local,
          app->local, fin[react->donor_idx[i]]);
      else
        gk_neut_species_moment_calc(&react->moms_donor[i], app->neut_species[react->donor_idx[i]].local,
          app->local, fin_neut[react->donor_idx[i]]);

      gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, 0, react->moms_elc[i].marr, 0, react->moms_elc[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, 1, react->moms_elc[i].marr, 1, react->moms_elc[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, 2, react->moms_elc[i].marr, 2, react->moms_elc[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      
      gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, 0, react->moms_donor[i].marr, 0, react->moms_donor[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, 1, react->moms_donor[i].marr, 1, react->moms_donor[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, 2, react->moms_donor[i].marr, 2, react->moms_donor[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      
      gkyl_array_set_range(react->m0_elc[i], 1.0, react->moms_elc[i].marr, &app->local);
      gkyl_array_set_range(react->m0_donor[i], 1.0, react->moms_donor[i].marr, &app->local);

      // When do moments become infinitely large? 
      //gkyl_grid_sub_array_write(&species->grid, &species->local, react->moms_elc[i].marr, "moms_elc.gkyl");

      // compute ionization reaction rate
      gkyl_dg_iz_coll(react->iz[i], react->moms_elc[i].marr, react->moms_donor[i].marr,
    	app->gk_geom->b_i, react->vt_sq_iz[i], react->prim_vars[i],
        react->coeff_react[i], 0);
    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      // compute needed moments
      gk_species_moment_calc(&react->moms_elc[i], app->species[react->elc_idx[i]].local,
        app->local, fin[react->elc_idx[i]]);
      gk_species_moment_calc(&react->moms_ion[i], app->species[react->ion_idx[i]].local,
        app->local, fin[react->ion_idx[i]]);

      gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, 0, react->moms_elc[i].marr, 0, react->moms_elc[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, 1, react->moms_elc[i].marr, 1, react->moms_elc[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_div_op_range(react->moms_elc[i].mem_geo, app->confBasis, 2, react->moms_elc[i].marr, 2, react->moms_elc[i].marr, 0, app->gk_geom->jacobgeo, &app->local);

      gkyl_dg_div_op_range(react->moms_ion[i].mem_geo, app->confBasis, 0, react->moms_ion[i].marr, 0, react->moms_ion[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_div_op_range(react->moms_ion[i].mem_geo, app->confBasis, 1, react->moms_ion[i].marr, 1, react->moms_ion[i].marr, 0, app->gk_geom->jacobgeo, &app->local);
      gkyl_dg_div_op_range(react->moms_ion[i].mem_geo, app->confBasis, 2, react->moms_ion[i].marr, 2, react->moms_ion[i].marr, 0, app->gk_geom->jacobgeo, &app->local);

      gkyl_array_set_range(react->m0_elc[i], 1.0, react->moms_elc[i].marr, &app->local);
      
      // compute recombination reaction rate
      gkyl_dg_recomb_coll(react->recomb[i], react->moms_elc[i].marr, react->moms_ion[i].marr,
        app->gk_geom->b_i, react->prim_vars[i], react->coeff_react[i], 0);
    }
  }
}

// updates the reaction terms in the rhs
void
gk_species_react_rhs(gkyl_gyrokinetic_app *app, const struct gk_species *s,
  struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  //fprintf(stdout, "\nspecies react");
  for (int i=0; i<react->num_react; ++i) {
    gkyl_array_clear(react->f_react, 0.0);

    if (react->react_id[i] == GKYL_REACT_IZ) {
      if (react->type_self[i] == GKYL_SELF_ELC) {
        gkyl_array_set_offset(react->prim_vars[i], 1.0, react->vt_sq_iz[i], 1*app->confBasis.num_basis);
        gkyl_proj_gkmaxwellian_on_basis_prim_mom(react->proj_max, &s->local, &app->local,
          react->moms_elc[i].marr, react->prim_vars[i],
          app->gk_geom->bmag, app->gk_geom->jacobtot, s->info.mass, react->f_react);

        // electron update is n_elc*coeff_react*(2*fmax(n_elc, upar_donor, vtiz^2) - f_elc)
        gkyl_array_scale(react->f_react, 2.0);
        gkyl_array_accumulate(react->f_react, -1.0, fin);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_donor[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, 1.0, react->f_react);

	//fprintf(stdout, "\nApplying iz elc rhs");
      }
      else if (react->type_self[i] == GKYL_SELF_ION) {
        gkyl_proj_gkmaxwellian_on_basis_prim_mom(react->proj_max, &s->local, &app->local,
          react->moms_donor[i].marr, react->prim_vars[i], //react->moms_donor[i].marr, react->prim_vars[i],
          app->gk_geom->bmag, app->gk_geom->jacobtot, s->info.mass, react->f_react);

        // ion update is n_elc*coeff_react*fmax(n_donor, upar_donor, vt_donor^2)
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_elc[i], react->f_react, &app->local, &s->local);
	gkyl_array_accumulate(rhs, 1.0, react->f_react);

	//fprintf(stdout, "\nApplying iz ion rhs");
      }
      else {
        // donor update is -n_elc*coeff_react*f_donor
        gkyl_array_set(react->f_react, 1.0, fin);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_elc[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, -1.0, react->f_react);
	//fprintf(stdout, "\nApplying iz donor rhs");
      }
    }
    else if (react->react_id[i] == GKYL_REACT_RECOMB) {
      if (react->type_self[i] == GKYL_SELF_ELC || react->type_self[i] == GKYL_SELF_ION) {
        // update is -n_elc*coeff_react*f_s where s = (electron, ion)
        gkyl_array_set(react->f_react, 1.0, fin);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
	//gkyl_grid_sub_array_write(&app->grid, &app->local, react->coeff_react[i], "coef_recomb.gkyl");
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_elc[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, -1.0, react->f_react);
	//gkyl_grid_sub_array_write(&s->grid, &s->local, react->f_react, "f_recomb.gkyl");
      }
      else {
        gkyl_proj_gkmaxwellian_on_basis_lab_mom(react->proj_max, &s->local, &app->local,
          react->moms_ion[i].marr, app->gk_geom->bmag, app->gk_geom->jacobtot, s->info.mass, react->f_react);

        // receiver update is n_elc*coeff_react*fmax(n_ion, upar_ion, vt_ion^2)
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->coeff_react[i], react->f_react, &app->local, &s->local);
        gkyl_dg_mul_conf_phase_op_range(&app->confBasis, &app->basis, react->f_react,
            react->m0_elc[i], react->f_react, &app->local, &s->local);
        gkyl_array_accumulate(rhs, 1.0, react->f_react);
      }
    }
  }
}

void 
gk_species_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react)
{
  gkyl_array_release(react->f_react);
  gkyl_proj_maxwellian_on_basis_release(react->proj_max);
  for (int i=0; i<react->num_react; ++i) {
    gk_species_moment_release(app, &react->moms_elc[i]);
    gk_species_moment_release(app, &react->moms_ion[i]);
    if (react->all_gk && gk_find_species(app, react->react_type[i].donor_nm))
      gk_species_moment_release(app, &react->moms_donor[i]);
    else if (gk_find_neut_species(app, react->react_type[i].donor_nm))
      gk_neut_species_moment_release(app, &react->moms_donor[i]);

    gkyl_array_release(react->coeff_react[i]);
    gkyl_array_release(react->vt_sq_iz[i]);
    gkyl_array_release(react->m0_elc[i]);
    gkyl_array_release(react->m0_ion[i]);
    gkyl_array_release(react->m0_donor[i]);
    gkyl_array_release(react->prim_vars[i]); 
    if (react->react_id[i] == GKYL_REACT_IZ) 
      gkyl_dg_iz_release(react->iz[i]);
    else if (react->react_id[i] == GKYL_REACT_RECOMB)  
      gkyl_dg_recomb_release(react->recomb[i]);
  }
}