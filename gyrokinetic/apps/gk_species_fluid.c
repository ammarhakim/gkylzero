#include <gkyl_gk_species_priv.h>

void
gk_species_fluid_file_import_init(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_gyrokinetic_ic_import inp
)
{
  if (inp.type != GKYL_IC_IMPORT_F || inp.enforce_positivity || inp.conf_scale ||
      inp.jacobtot_inv_file_name[0] || inp.jacobvel_file_name[0]) {
    gkyl_exit("Fluid initial-state import supports unmodified configuration-space moments only.");
  }
  enum gkyl_array_rio_status status = gk_species_read(app, gks, inp.file_name);
  if (status != GKYL_ARRAY_RIO_SUCCESS) {
    fprintf(
      stderr, "Unable to import fluid species '%s' from '%s': %s\n", gks->info.name, inp.file_name,
      gkyl_array_rio_status_msg(status)
    );
    gkyl_exit("Fluid initial-state import failed.");
  }
}

static void
gk_species_fluid_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  gk_species_release_integrated_mom(app, gks);
  if (gks->info.init_from_file.type == GKYL_IC_IMPORT_NONE) {
    gk_species_projection_release(app, &gks->proj_init);
  }
  gk_species_moment_release(app, &gks->m0);
  gk_species_moment_release(app, &gks->lte.moms);
  for (int m = 0; m < gks->info.num_diag_moments; ++m) {
    gk_species_moment_release(app, &gks->moms[m]);
  }
  gkyl_free(gks->moms);
  gk_species_scaling_release(app, &gks->sca);
  gk_species_lbo_release(app, &gks->lbo);
  gk_species_bgk_release(app, &gks->bgk);
  gkyl_array_release(gks->m0_gyroavg);
  gkyl_array_release(gks->gyro_phi);
  gkyl_array_release(gks->f);
  if (app->use_gpu) {
    gkyl_array_release(gks->f_host);
  }
  gkyl_comm_release(gks->comm);
  gkyl_msgpack_map_elem_release(gks->io_meta_basic_len, gks->io_meta_basic);
  gkyl_msgpack_map_elem_release(gks->io_meta_conf_len, gks->io_meta_conf);
  gkyl_msgpack_map_elem_release(gks->io_meta_phase_len, gks->io_meta_phase);
}

void
gk_species_fluid_init(struct gkyl_gk *inp, gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Fluid species currently have algebraic density closures, not transport equations.
  if (gks->info.collisionless.type || gks->info.flr.type || gks->info.source.source_id ||
      gks->info.source_bgk.source_bgk_id || gks->info.radiation.radiation_id ||
      gks->info.react.num_react || gks->info.react_neut.num_react || gks->info.positivity.type ||
      gks->info.positivity.quasineutrality_rescale || gks->info.damping.type ||
      gks->info.anomalous_diffusion.anomalous_diff_id ||
      gks->info.time_rate_multiplier.num_multipliers || gks->info.scale_with_polarization ||
      gks->info.time_rate_diagnostics || gks->info.write_omega_cfl || gks->info.mapc2p.mapping ||
      gks->info.boundary_flux_diagnostics.num_diag_moments ||
      gks->info.boundary_flux_diagnostics.num_integrated_diag_moments) {
    gkyl_exit("Fluid gyrokinetic species: requested kinetic-only operator or diagnostic.");
  }
  if (gks->info.collisions.collision_id && !gks->info.collisions.not_in_dfdt) {
    gkyl_exit("Fluid gyrokinetic collisions require not_in_dfdt=true.");
  }
  if (gks->info.scaling.type == GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE) {
    gkyl_exit("Recycling/ionization balance is not supported for charged fluid species.");
  }
  for (int d = 0; d < GKYL_MAX_VDIM; ++d) {
    if (gks->info.cells[d] || gks->info.lower[d] || gks->info.upper[d]) {
      gkyl_exit("Fluid gyrokinetic species must omit the velocity grid entirely.");
    }
  }
  // Fluid moments always represent parallel and perpendicular velocity coordinates.
  gks->info.vdim = 2;
  for (int b = 0; b < 2 * app->cdim; ++b) {
    if (gks->info.bcs[b].type != GKYL_BC_GK_SKIP &&
        gks->info.bcs[b].type != GKYL_BC_GK_SPECIES_PERIODIC) {
      gkyl_exit("Kinetic boundary conditions are not supported for fluid gyrokinetic species.");
    }
  }
  if (gks->info.mass <= 0.0 || gks->info.charge == 0.0) {
    gkyl_exit("Invalid fluid gyrokinetic species mass or charge.");
  }
  gks->info.is_static = true;
  gks->release_species_func = gk_species_fluid_release;
  gks->grid = app->grid;
  gks->basis = app->basis;
  gks->basis_on_dev = app->basis_on_dev;
  gks->local = app->local;
  gks->local_ext = app->local_ext;
  gks->global = app->global;
  gks->global_ext = app->global_ext;
  gks->comm = gkyl_comm_acquire(app->comm);
  int ghost[GKYL_MAX_CDIM] = {1, 1, 1};
  for (int d = 0; d < app->cdim; ++d) {
    gkyl_skin_ghost_ranges(
      &gks->local_lower_skin[d], &gks->local_lower_ghost[d], d, GKYL_LOWER_EDGE, &gks->local_ext,
      ghost
    );
    gkyl_skin_ghost_ranges(
      &gks->local_upper_skin[d], &gks->local_upper_ghost[d], d, GKYL_UPPER_EDGE, &gks->local_ext,
      ghost
    );
    gkyl_skin_ghost_ranges(
      &gks->global_lower_skin[d], &gks->global_lower_ghost[d], d, GKYL_LOWER_EDGE, &gks->global_ext,
      ghost
    );
    gkyl_skin_ghost_ranges(
      &gks->global_upper_skin[d], &gks->global_upper_ghost[d], d, GKYL_UPPER_EDGE, &gks->global_ext,
      ghost
    );
  }
  gks->f = mkarr(app->use_gpu, 4 * app->basis.num_basis, app->local_ext.volume);
  gks->f_host = app->use_gpu ? mkarr(false, gks->f->ncomp, gks->f->size) : gks->f;
  gk_species_init_io(app, gks);
  gk_species_moment_init(app, gks, &gks->m0, GKYL_F_MOMENT_M0, false);
  gk_species_moment_init(app, gks, &gks->lte.moms, GKYL_F_MOMENT_MAXWELLIAN, false);
  gks->m0_gyroavg = gkyl_array_acquire(gks->m0.marr);
  gks->gyro_phi = gkyl_array_acquire(app->field->phi_smooth);
  gks->gyroaverage = gk_species_gyroaverage_disabled;
  gks->moms = gkyl_malloc(sizeof(struct gk_species_moment) * gks->info.num_diag_moments);
  for (int m = 0; m < gks->info.num_diag_moments; ++m) {
    gk_species_moment_init(app, gks, &gks->moms[m], gks->info.diag_moments[m], false);
  }
  if (gks->info.init_from_file.type == GKYL_IC_IMPORT_NONE) {
    gk_species_projection_init(app, gks, gks->info.projection, &gks->proj_init);
  } else {
    gk_species_file_import_init(app, gks, gks->info.init_from_file);
  }
  gk_species_collisionless_init(app, gks, &gks->collisionless);
  gk_species_anomalous_diff_init(app, gks, &gks->anom_diff);
  gk_species_damping_init(app, gks, &gks->damping);
  gk_species_fdot_multiplier_init(app, gks, &gks->fdot_mult);
  gk_species_source_bgk_init(app, gks, &gks->bgk_src);
  gk_species_positivity_init(app, gks, &gks->positivity);
  gk_species_react_init(app, gks, gks->info.react, &gks->react, true);
  gk_species_react_init(app, gks, gks->info.react_neut, &gks->react_neut, false);
  gk_species_bflux_init(
    app, gks, &gks->bflux, GK_SPECIES_BFLUX_NONE, (struct gkyl_phase_diagnostics_inp){}
  );
  gk_species_scaling_init(app, gks, &gks->sca);
  gk_species_lbo_init(app, gks, &gks->lbo);
  gk_species_bgk_init(app, gks, &gks->bgk);
  gk_species_init_static(inp, app, gks);
  gk_species_init_integrated_mom(app, gks);
  gks->write_func = gk_species_write_dynamic;
  gks->write_mom_func = gk_species_write_mom_dynamic;
  gks->calc_integrated_mom_func = gk_species_calc_integrated_mom_dynamic;
  gks->write_integrated_mom_func = gk_species_write_integrated_mom_dynamic;
}
