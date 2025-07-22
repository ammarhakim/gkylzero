#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

void 
gk_neut_species_recycle_react_scale_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns, 
  struct gkyl_gyrokinetic_recycling_reaction_scaling_inp *inp, struct gk_recycle_react_scale *rrs)
{
  rrs->num_react = inp->num_react; 
  // Initialize information about reactions from input struct.
  for (int i=0; i<rrs->num_react; ++i) 
    rrs->react_type[i] = inp->react_type[i];
  rrs->write_diagnostics = inp->write_diagnostics;

  int num_mom = 1; // M0.
  if (app->use_gpu){
    rrs->red_integ_mom = gkyl_cu_malloc(sizeof(double[num_mom]));
    rrs->red_integ_mom_global = gkyl_cu_malloc(sizeof(double[num_mom]));
  }
  else {
    rrs->red_integ_mom = gkyl_malloc(sizeof(double[num_mom]));
    rrs->red_integ_mom_global = gkyl_malloc(sizeof(double[num_mom]));
  }

  rrs->num_boundaries = ns->info.recycling_reaction_scaling.num_boundaries;
  for (int j=0; j < rrs->num_boundaries; ++j) {
    int dir  = inp->boundaries_dir[j];
    int edge = inp->boundaries_edge[j];

    // Source adaptation on periodic, zero flux, or reflect boundary is not allowed.
    assert(ns->bc_is_np[dir]);
    if (edge == GKYL_LOWER_EDGE) {
      assert(ns->lower_bc[dir].type != GKYL_SPECIES_ZERO_FLUX);
      assert(ns->lower_bc[dir].type != GKYL_SPECIES_REFLECT);
    } else {
      assert(ns->upper_bc[dir].type != GKYL_SPECIES_ZERO_FLUX);
      assert(ns->upper_bc[dir].type != GKYL_SPECIES_REFLECT);
    }

    // Default scenario: we set the ranges to the full range of the ghost cells.
    rrs->boundaries_conf_ghost[j] = edge == GKYL_LOWER_EDGE ? app->lower_ghost[dir] : app->upper_ghost[dir];
    rrs->boundaries_dir[j]  = dir;
    rrs->boundaries_edge[j] = edge;

    // Specific scenario if we are in a inner wall limited case. We select only SOL range in parallel direction.
    if (edge == GKYL_LOWER_EDGE? ns->lower_bc[dir].type == GKYL_SPECIES_GK_IWL
                               : ns->upper_bc[dir].type == GKYL_SPECIES_GK_IWL)
    {
      rrs->boundaries_conf_ghost[j] = edge == GKYL_LOWER_EDGE ? app->lower_ghost_par_sol : app->upper_ghost_par_sol;
    }
  }
}

void 
gk_neut_species_scale_by_react_cross_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns,
  struct gk_react *react)
{
  // react->write_diagnostics = ns->info.react_neut.write_diagnostics;
  // distribution function which holds update for each reaction
  // form depend on react->type_self, e.g., for recombination and react->type_self == GKYL_SELF_RECVR
  // react->f_react = n_elc*coeff_react*fmax(n_ion, upar_ion b_i, vt_ion^2)
  react->f_react = mkarr(app->use_gpu, ns->basis.num_basis, ns->local_ext.volume);

  for (int i=0; i<react->num_react; ++i) {
    // Fetch index of species for indexing arrays
    react->elc_idx[i] = gk_find_species_idx(app, react->react_type[i].elc_nm);
    react->ion_idx[i] = gk_find_species_idx(app, react->react_type[i].ion_nm);

    react->coeff_react[i] = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    if (react->write_diagnostics) {
      react->coeff_react_host[i] = react->coeff_react[i];
      if (app->use_gpu) {
        react->coeff_react_host[i] = mkarr(false, react->coeff_react[i]->ncomp, react->coeff_react[i]->size);
      }
    }

    struct gkyl_dg_iz_inp iz_inp = {
      .cbasis = &app->basis, 
      .conf_rng = &app->local, 
      .mass_ion = react->react_type[i].ion_mass, 
      .type_ion = react->react_type[i].ion_id, 
      .charge_state = react->react_type[i].charge_state, 
      .type_self = GKYL_SELF_ION, // Could be GKYL_SELF_DONOR. It just can't be
                                  // GKYL_SELF_ELC because we don't need to
                                  // compute the ionization temperatures.
    };
  }
}

// Computes reaction coefficients.
void
gk_neut_species_scale_by_react_cross_moms(gkyl_gyrokinetic_app *app, const struct gk_neut_species *species,
  struct gk_react *react, const struct gkyl_array *fin[], const struct gkyl_array *fin_neut[])
{
  struct timespec wst = gkyl_wall_clock();    
  for (int i=0; i<react->num_react; ++i) {
    struct gk_species *gks_elc = &app->species[react->elc_idx[i]]; 
    struct gk_species *gks_ion = &app->species[react->ion_idx[i]]; 

    // Compute needed electron Maxwellian moments (J*n, u_par, T/m).
    gk_species_moment_calc(&gks_elc->lte.moms, 
      gks_elc->local, app->local, fin[react->elc_idx[i]]);

    // Divide out the Jacobian from the electron density for computing reaction rates.
    gkyl_dg_div_op_range(gks_elc->lte.moms.mem_geo, app->basis, 
      0, gks_elc->lte.moms.marr, 0, gks_elc->lte.moms.marr, 0, 
      app->gk_geom->jacobgeo, &app->local); 

    // Compute ionization reaction rate from input electron primitive moments.
    gkyl_dg_iz_coll(react->iz[i], gks_elc->lte.moms.marr, 
      react->coeff_react[i], react->coeff_react[i], react->coeff_react[i], 0);
  }
  app->stat.neut_species_react_mom_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_neut_species_scale_by_react_rhs(gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_react *react, const struct gkyl_array *fin, struct gkyl_array *rhs)
{
  struct timespec wst = gkyl_wall_clock();  
  for (int i=0; i<react->num_react; ++i) {
    gkyl_array_clear(react->f_react, 0.0);

    struct gk_species *gks_elc = &app->species[react->elc_idx[i]];
    struct gk_species *gks_ion = &app->species[react->ion_idx[i]];

    // Neutral species can only be the donor or the receiver in reaction
    // donor update is -n_elc*coeff_react*f_donor.

    // Accumulate -n_elc*(J*f_donor) (*note* Jacobian factor already included in fin)
    gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &ns->basis, react->f_react,
      -1.0, gks_elc->lte.moms.marr, fin, &app->local, &ns->local);  

    // Accumulate reaction update to rhs 
    gkyl_dg_mul_conf_phase_op_accumulate_range(&app->basis, &ns->basis, rhs,
      1.0, react->coeff_react[i], react->f_react, &app->local, &ns->local);  
  }
  app->stat.neut_species_react_tm += gkyl_time_diff_now_sec(wst);
}

void
gk_neut_species_scale_by_react_write(gkyl_gyrokinetic_app* app, struct gk_neut_species *gkns, struct gk_react *gkr,
  int ridx, double tm, int frame)
{
  // React diagnostics usually written from gk_species.
  // In the case of static gk_species, write_diagnostics flag
  // can be used to check reaction rates from gk_neut_species.
  if (gkr->write_diagnostics) {
    struct timespec wst = gkyl_wall_clock();
    // Compute reaction rate
    const struct gkyl_array *fin[app->num_species];
    const struct gkyl_array *fin_neut[app->num_neut_species];
    for (int i=0; i<app->num_species; ++i) {
      fin[i] = app->species[i].f;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      fin_neut[i] = app->neut_species[i].f;
    }
    gk_neut_species_scale_by_react_cross_moms(app, gkns, gkr, fin, fin_neut);
    app->stat.neut_species_diag_calc_tm += gkyl_time_diff_now_sec(wst);
    
    struct timespec wtm = gkyl_wall_clock();
    struct gkyl_msgpack_data *mt = gk_array_meta_new( (struct gyrokinetic_output_meta) {
        .frame = frame,
        .stime = tm,
        .poly_order = app->poly_order,
        .basis_type = app->basis.id
      }
    );

    if (app->use_gpu)
      gkyl_array_copy(gkr->coeff_react_host[ridx], gkr->coeff_react[ridx]);
    
    const char *fmt = "%s-%s_%s_%s_iz_react_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, gkns->info.name,
      gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].ion_nm, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, gkns->info.name,
     gkr->react_type[ridx].elc_nm, gkr->react_type[ridx].ion_nm, frame);
  
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, 
      gkr->coeff_react_host[ridx], fileNm);
    app->stat.n_neut_diag_io += 1;
    
    gk_array_meta_release(mt); 
    app->stat.neut_species_diag_io_tm += gkyl_time_diff_now_sec(wtm);
  }
}

void 
gk_neut_species_scale_by_react_release(const struct gkyl_gyrokinetic_app *app, const struct gk_react *react)
{
  gkyl_array_release(react->f_react);

  for (int i=0; i<react->num_react; ++i) {
    gkyl_array_release(react->coeff_react[i]);
    if (react->write_diagnostics) {
      if (app->use_gpu) {
        gkyl_array_release(react->coeff_react_host[i]);
      }
    }

    gkyl_dg_iz_release(react->iz[i]);
  }
}
