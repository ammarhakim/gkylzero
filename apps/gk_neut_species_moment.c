#include <assert.h>
#include <gkyl_gyrokinetic_priv.h>

// initialize neutral species moment object
void
gk_neut_species_moment_init(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s,
  struct gk_species_moment *sm, enum gkyl_distribution_moments mom_type, bool is_integrated, bool is_external)
{
  sm->is_integrated = is_integrated;
  sm->is_maxwellian_moms = mom_type == GKYL_F_MOMENT_LTE;
  sm->is_external = is_external;

  if (sm->is_integrated) {
    // Create moment operator.
    struct gkyl_mom_canonical_pb_auxfields can_pb_inp = {.hamil = s->hamil};
    sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->basis, 
      &s->basis, &app->local, &s->local_vel, &s->local, s->model_id, &can_pb_inp, 
      mom_type, sm->is_integrated, app->use_gpu);

    sm->num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);

    // Allocate arrays to hold moments.
    sm->marr = mkarr(app->use_gpu, sm->num_mom, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu) {
      sm->marr_host = mkarr(false, sm->num_mom, app->local_ext.volume); 
    }
  }
  else {
    // Create moment operator.
    if (sm->is_maxwellian_moms) {
      struct gkyl_vlasov_lte_moments_inp inp_mom = {
        .phase_grid = &s->grid,
        .vel_grid = &s->grid_vel, 
        .conf_basis = &app->basis,
        .phase_basis = &s->basis,
        .conf_range =  &app->local,
        .conf_range_ext = &app->local_ext,
        .vel_range = &s->local_vel,
        .phase_range = &s->local,
        .h_ij = s->g_ij,
        .h_ij_inv = s->gij,
        .det_h = app->gk_geom->geo_int.jacobgeo,
        .hamil = s->hamil,
        .model_id = s->model_id,
        .use_gpu = app->use_gpu,
      };
      sm->vlasov_lte_moms = gkyl_vlasov_lte_moments_inew(&inp_mom);
      sm->num_mom = 5; // (n, ux, uy, uz, T/m).
    }
    else {
      struct gkyl_mom_canonical_pb_auxfields can_pb_inp = {.hamil = s->hamil};
      sm->mcalc = gkyl_dg_updater_moment_new(&s->grid, &app->basis, 
        &s->basis, &app->local, &s->local_vel, &s->local, s->model_id, &can_pb_inp, 
        mom_type, sm->is_integrated, app->use_gpu);

      sm->num_mom = gkyl_dg_updater_moment_num_mom(sm->mcalc);
    }

    // Allocate arrays to hold moments.
    sm->marr = mkarr(app->use_gpu, sm->num_mom*app->basis.num_basis, app->local_ext.volume);
    sm->marr_host = sm->marr;
    if (app->use_gpu)
      sm->marr_host = mkarr(false, sm->num_mom*app->basis.num_basis, app->local_ext.volume);

    // Bin Op memory for rescaling moment by inverse of Jacobian
    if (app->use_gpu) {
      sm->mem_geo = gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis);
    }
    else {
      sm->mem_geo = gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);
    }
  }
}

static void gk_neut_species_read_moment(struct gkyl_gyrokinetic_app *app, struct gkyl_array *out, cstr fileNm)
{
  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &app->local, app->poly_order);
  struct gkyl_array* nnodal = mkarr(false, 1, nrange.volume);

  FILE *ptr = fopen(fileNm.str,"r");
  size_t status;

  int nr = gkyl_range_shape(&nrange, 0);
  int nz = gkyl_range_shape(&nrange, 1);
  int idx[2];

  for(int ir = 0; ir < nr; ir++){
    idx[0] = ir;
    for(int iz = 0; iz < nz; iz++){
      idx[1] = iz;
      // set psi
      double *nnodal_n = gkyl_array_fetch(nnodal, gkyl_range_idx(&nrange, idx));
      status = fscanf(ptr,"%lf", nnodal_n);
    }
  }

  fclose(ptr);

  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&app->basis, &app->grid, false);
  gkyl_nodal_ops_n2m(n2m, &app->basis, &app->grid, &nrange, &app->local, 1, nnodal, out, false);
  gkyl_array_release(nnodal);
}

void gk_neut_species_read_lte(struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
    cstr fileNm = cstr_from_fmt("gkeyll_text_input/%s-%s_ux.txt", app->name, s->info.name);
    gk_neut_species_read_moment(app, s->m0.marr_host, fileNm);
    gkyl_array_set_offset(s->lte.moms.marr_host, 1.0, s->m0.marr_host, 1*app->basis.num_basis);
    fileNm = cstr_from_fmt("gkeyll_text_input/%s-%s_uy.txt", app->name, s->info.name);
    gk_neut_species_read_moment(app, s->m0.marr_host, fileNm);
    gkyl_array_set_offset(s->lte.moms.marr_host, 1.0, s->m0.marr_host, 2*app->basis.num_basis);
    fileNm = cstr_from_fmt("gkeyll_text_input/%s-%s_uz.txt", app->name, s->info.name);
    gk_neut_species_read_moment(app, s->m0.marr_host, fileNm);
    gkyl_array_set_offset(s->lte.moms.marr_host, 1.0, s->m0.marr_host, 3*app->basis.num_basis);

    fileNm = cstr_from_fmt("gkeyll_text_input/%s-%s_Temp.txt", app->name, s->info.name);
    gk_neut_species_read_moment(app, s->m0.marr_host, fileNm);
    gkyl_array_scale(s->m0.marr_host, 1/s->info.mass);
    gkyl_array_set_offset(s->lte.moms.marr_host, 1.0, s->m0.marr_host, 4*app->basis.num_basis);

    fileNm = cstr_from_fmt("gkeyll_text_input/%s-%s_M0.txt", app->name, s->info.name);
    gk_neut_species_read_moment(app, s->m0.marr_host, fileNm);
    gkyl_array_set_offset(s->lte.moms.marr_host, 1.0, s->m0.marr_host, 0*app->basis.num_basis);

    if (app->use_gpu) {
      gkyl_array_copy(s->m0.marr, s->m0.marr_host);
      gkyl_array_copy(s->lte.moms.marr, s->lte.moms.marr_host);
    }
}


void
gk_neut_species_moment_calc(const struct gk_species_moment *sm,
  const struct gkyl_range phase_rng, const struct gkyl_range conf_rng,
  const struct gkyl_array *fin)
{
  if (sm->is_external)
    return;
  if (sm->is_maxwellian_moms) {
    gkyl_vlasov_lte_moments_advance(sm->vlasov_lte_moms, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }
  else {
    gkyl_dg_updater_moment_advance(sm->mcalc, 
      &phase_rng, &conf_rng, fin, sm->marr);
  }  
}

void
gk_neut_species_moment_release(const struct gkyl_gyrokinetic_app *app, const struct gk_species_moment *sm)
{
  gkyl_array_release(sm->marr);
  if (app->use_gpu)
    gkyl_array_release(sm->marr_host);

  if (sm->is_integrated) {
    gkyl_dg_updater_moment_release(sm->mcalc);
  }
  else {
    if (sm->is_maxwellian_moms) {
      gkyl_vlasov_lte_moments_release(sm->vlasov_lte_moms);
    }
    else {
      gkyl_dg_updater_moment_release(sm->mcalc);
    }

    // Free the weak division memory.
    gkyl_dg_bin_op_mem_release(sm->mem_geo);
  }
}
