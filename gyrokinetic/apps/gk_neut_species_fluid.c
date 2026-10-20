#include <gkyl_gk_neut_species_priv.h>

static void
gk_neut_species_fluid_diffusion_dirichlet_ghost(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, int dir,
  enum gkyl_edge_loc edge, const struct gkyl_range *ghost, double number_density,
  const struct gkyl_array *number_density_field, struct gkyl_array *density,
  struct gkyl_array *tensor
)
{
  // Communication has already populated internal-rank ghosts. Only a rank
  // owning the physical edge supplies the exterior Dirichlet state.
  if ((edge == GKYL_LOWER_EDGE && gkns->local.lower[dir] != gkns->global.lower[dir]) ||
      (edge == GKYL_UPPER_EDGE && gkns->local.upper[dir] != gkns->global.upper[dir])) {
    return;
  }

  const int skin_idx_dir =
    edge == GKYL_LOWER_EDGE ? gkns->local.lower[dir] : gkns->local.upper[dir];
  const int num_basis = app->basis.num_basis;
  const int num_tensor = app->cdim * app->cdim;
  const double rho_boundary_modal =
    pow(sqrt(2.0), app->cdim) * gkns->info.mass * number_density;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, ghost);
  while (gkyl_range_iter_next(&iter)) {
    int skin_idx[GKYL_MAX_CDIM];
    for (int d = 0; d < app->cdim; ++d) {
      skin_idx[d] = iter.idx[d];
    }
    skin_idx[dir] = skin_idx_dir;

    const double *rho_skin =
      gkyl_array_cfetch(density, gkyl_range_idx(&gkns->local_ext, skin_idx));
    double *rho_ghost =
      gkyl_array_fetch(density, gkyl_range_idx(&gkns->local_ext, iter.idx));

    // Mirror the skin polynomial and choose its exterior value so the
    // arithmetic trace at the physical face is rho_b:
    //   rho_g = 2 rho_b - R(rho_s).
    app->basis.flip_odd_sign(dir, rho_skin, rho_ghost);
    for (int k = 0; k < num_basis; ++k) {
      rho_ghost[k] = -rho_ghost[k];
    }
    if (number_density_field) {
      const double *n_boundary =
        gkyl_array_cfetch(number_density_field, gkyl_range_idx(&gkns->local_ext, iter.idx));
      for (int k = 0; k < num_basis; ++k) {
        rho_ghost[k] += 2.0 * gkns->info.mass * n_boundary[k];
      }
    } else {
      rho_ghost[0] += 2.0 * rho_boundary_modal;
    }

    const double *tensor_skin =
      gkyl_array_cfetch(tensor, gkyl_range_idx(&gkns->local_ext, skin_idx));
    double *tensor_ghost =
      gkyl_array_fetch(tensor, gkyl_range_idx(&gkns->local_ext, iter.idx));
    for (int k = 0; k < num_tensor; ++k) {
      app->basis.flip_odd_sign(
        dir, &tensor_skin[k * num_basis], &tensor_ghost[k * num_basis]
      );
    }
  }
}

static void
gk_neut_species_fluid_recycling_density(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, int dir,
  enum gkyl_edge_loc edge, const struct gkyl_range *ghost
)
{
  const int edge_idx = edge == GKYL_LOWER_EDGE ? 0 : 1;
  const struct gkyl_gyrokinetic_bc *bc =
    edge == GKYL_LOWER_EDGE ? &gkns->lower_bc[dir] : &gkns->upper_bc[dir];
  gkyl_array_clear_range(gkns->diffusion_bc_density, 0.0, ghost);

  // A unit-density isotropic Maxwellian emits the one-sided particle flux
  // Gamma_M,1=sqrt(T_n/(2*pi*m_n)). Match the kinetic recycling convention:
  // n_b*Gamma_M,1=R*Gamma_i, locally along the target.
  const double unit_flux = sqrt(bc->emission.emission_temp / (2.0 * M_PI * gkns->info.mass));
  for (int i = 0; i < gkns->diffusion_recycling_num_species[edge_idx]; ++i) {
    struct gk_species *ion = &app->species[gkns->diffusion_recycling_ion_idx[edge_idx][i]];
    assert(ion->bflux.current_moms);
    gkyl_array_clear_range(gkns->diffusion_bc_density_tmp, 0.0, ghost);
    gk_species_bflux_get_flux_mom(
      &ion->bflux, dir, edge, GKYL_F_MOMENT_M0, ion->bflux.current_moms,
      gkns->diffusion_bc_density_tmp, ghost
    );
    gkyl_array_accumulate_range(
      gkns->diffusion_bc_density, bc->emission.recycling_frac / unit_flux,
      gkns->diffusion_bc_density_tmp, ghost
    );
  }
}

static void
gk_neut_species_fluid_apply_bc_dynamic(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, struct gkyl_array *f
)
{
  struct timespec wst = gkyl_wall_clock();

  gkyl_comm_array_per_sync(
    gkns->comm, &gkns->local, &gkns->local_ext, app->num_periodic_dir, app->periodic_dirs, f
  );
  gkyl_comm_array_sync(gkns->comm, &gkns->local, &gkns->local_ext, f);

  app->stat.neut_species_bc_tm += gkyl_time_diff_now_sec(wst);
}

static double
gk_neut_species_fluid_rhs_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  double omega_cfl = 1 / DBL_MAX;
  gkyl_array_clear(gkns->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Collisionless terms.
  struct timespec wst = gkyl_wall_clock();
  if (gkns->has_diffusion) {
    if (gkns->use_reaction_rate_diffusion) {
      struct gk_react *react = &gkns->react_neut;
      int ionization_idx = gkns->ionization_react_idx;
      int cx_idx = gkns->diffusion_cx_react_idx;
      struct gk_species *electron = &app->species[react->elc_idx[ionization_idx]];

      // Evaluate <sigma v>_cx for the diffusion closure independently of the
      // neutral-density guard used by the reaction source.
      gkyl_dg_cx_coll_rate(
        react->cx[cx_idx], app->species[react->ion_idx[cx_idx]].lte.moms.marr,
        gkns->lte.moms.marr, react->upar_ion[cx_idx], gkns->diffusion_coeff
      );

      // denominator = ne*<sigma v>_cx.
      gkyl_dg_mul_op_range(
        &app->basis, 0, gkns->diffusion_moment_ratio, 0, electron->lte.moms.marr, 0,
        gkns->diffusion_coeff, &gkns->local
      );

      // Apply a user-selected floor to the cell-average collision frequency.
      // When activated, a cell below the floor is replaced by the constant DG
      // representation of the floor before weak division.
      const double minimum_frequency = gkns->info.diffusion.min_collision_frequency;
      if (minimum_frequency > 0.0) {
        const double dg_norm = pow(sqrt(2.0), app->cdim);
        struct gkyl_range_iter iter;
        gkyl_range_iter_init(&iter, &gkns->local);
        while (gkyl_range_iter_next(&iter)) {
          long loc = gkyl_range_idx(&gkns->local, iter.idx);
          double *frequency = gkyl_array_fetch(gkns->diffusion_moment_ratio, loc);
          if (frequency[0] / dg_norm < minimum_frequency) {
            for (int k = 0; k < app->basis.num_basis; ++k) {
              frequency[k] = 0.0;
            }
            frequency[0] = dg_norm * minimum_frequency;
          }
        }
      }

      // Isotropic scalar DG coefficient: D = vti^2/(ne*<sigma v>_cx).
      gkyl_dg_div_op_range(
        gkns->diffusion_div_mem, &app->basis, 0, gkns->diffusion_coeff, 0,
        react->vt_sq_ion[cx_idx], 0, gkns->diffusion_moment_ratio, &gkns->local
      );

      // The variable-coefficient surface terms read neighboring cells.
      // Populate coefficient ghosts independently of the state BC update.
      gkyl_comm_array_per_sync(
        gkns->comm, &gkns->local, &gkns->local_ext, gkns->num_periodic_dir,
        gkns->periodic_dirs, gkns->diffusion_coeff
      );
      gkyl_comm_array_sync(
        gkns->comm, &gkns->local, &gkns->local_ext, gkns->diffusion_coeff
      );
    }

    // The evolved density is N=J*rho. Form the physical density rho=N/J and
    // the mapped tensor K^{ij}=J*D*g_neut^{ij} used by
    //   d_t N = d_i(K^{ij} d_j rho).
    // Use weak division by J instead of multiplication by its separately
    // projected inverse. This preserves a spatially constant physical density
    // exactly when the stored state is N=J*rho.
    gkyl_dg_div_op_range(
      gkns->diffusion_div_mem, &app->basis, 0, gkns->diffusion_density, 0, fin, 0,
      app->gk_geom->geo_int.jacobgeo, &gkns->local
    );
    gkyl_dg_mul_op_range(
      &app->basis, 0, gkns->diffusion_geom_factor, 0, gkns->diffusion_coeff, 0,
      app->gk_geom->geo_int.jacobgeo, &gkns->local
    );

    const int metric_idx_1x[1][1] = {{5}};
    const int metric_idx_2x[2][2] = {{0, 2}, {2, 5}};
    const int metric_idx_3x[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};
    for (int i = 0; i < app->cdim; ++i) {
      for (int j = 0; j < app->cdim; ++j) {
        int gij_idx = app->cdim == 1 ? metric_idx_1x[i][j]
                                     : app->cdim == 2 ? metric_idx_2x[i][j]
                                                      : metric_idx_3x[i][j];
        gkyl_dg_mul_op_range(
          &app->basis, i * app->cdim + j, gkns->diffusion_tensor, 0,
          gkns->diffusion_geom_factor, gij_idx, app->gk_geom->geo_int.gij_neut, &gkns->local
        );
      }
    }
    gkyl_comm_array_per_sync(
      gkns->comm, &gkns->local, &gkns->local_ext, gkns->num_periodic_dir, gkns->periodic_dirs,
      gkns->diffusion_density
    );
    gkyl_comm_array_sync(
      gkns->comm, &gkns->local, &gkns->local_ext, gkns->diffusion_density
    );
    gkyl_comm_array_per_sync(
      gkns->comm, &gkns->local, &gkns->local_ext, gkns->num_periodic_dir, gkns->periodic_dirs,
      gkns->diffusion_tensor
    );
    gkyl_comm_array_sync(
      gkns->comm, &gkns->local, &gkns->local_ext, gkns->diffusion_tensor
    );

    // Fixed-density or locally recycled open-field-line boundaries.
    const int parallel_dir = app->cdim - 1;
    const struct gkyl_range *lower_ghost = app->gk_geom->has_LCFS
      ? &app->local_lower_ghost_par_sol
      : &gkns->local_lower_ghost[parallel_dir];
    const struct gkyl_range *upper_ghost = app->gk_geom->has_LCFS
      ? &app->local_upper_ghost_par_sol
      : &gkns->local_upper_ghost[parallel_dir];
    if (gkns->lower_bc[parallel_dir].type == GKYL_BC_GK_SPECIES_RECYCLE) {
      gk_neut_species_fluid_recycling_density(
        app, gkns, parallel_dir, GKYL_LOWER_EDGE, lower_ghost
      );
      gk_neut_species_fluid_diffusion_dirichlet_ghost(
        app, gkns, parallel_dir, GKYL_LOWER_EDGE, lower_ghost, 0.0,
        gkns->diffusion_bc_density,
        gkns->diffusion_density, gkns->diffusion_tensor
      );
    } else if (gkns->info.diffusion.lower_bc_type == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET) {
      gk_neut_species_fluid_diffusion_dirichlet_ghost(
        app, gkns, parallel_dir, GKYL_LOWER_EDGE, lower_ghost,
        gkns->info.diffusion.lower_bc_density, 0, gkns->diffusion_density,
        gkns->diffusion_tensor
      );
    }
    if (gkns->upper_bc[parallel_dir].type == GKYL_BC_GK_SPECIES_RECYCLE) {
      gk_neut_species_fluid_recycling_density(
        app, gkns, parallel_dir, GKYL_UPPER_EDGE, upper_ghost
      );
      gk_neut_species_fluid_diffusion_dirichlet_ghost(
        app, gkns, parallel_dir, GKYL_UPPER_EDGE, upper_ghost, 0.0,
        gkns->diffusion_bc_density,
        gkns->diffusion_density, gkns->diffusion_tensor
      );
    } else if (gkns->info.diffusion.upper_bc_type == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET) {
      gk_neut_species_fluid_diffusion_dirichlet_ghost(
        app, gkns, parallel_dir, GKYL_UPPER_EDGE, upper_ghost,
        gkns->info.diffusion.upper_bc_density, 0, gkns->diffusion_density,
        gkns->diffusion_tensor
      );
    }

    gkyl_dg_updater_conf_diffusion_advance(
      gkns->diffusion_slvr, &gkns->local, gkns->diffusion_density, gkns->cflrate, rhs
    );

    // Density is the only independently diffused quantity. Propagate its
    // change to the dependent conserved moments while holding their specific
    // values (moment/rho) fixed.
    for (int m = 1; m < gkns->num_moments; ++m) {
      gkyl_dg_div_op_range(
        gkns->diffusion_div_mem, &app->basis, m, gkns->diffusion_moment_ratio, m, fin, 0, fin,
        &gkns->local
      );
      gkyl_dg_mul_op_range(
        &app->basis, m, rhs, m, gkns->diffusion_moment_ratio, 0, rhs, &gkns->local
      );
    }
  }
  app->stat.neut_species_collisionless_tm += gkyl_time_diff_now_sec(wst);

  // Ionization loss and other enabled reactions with charged species.
  gk_neut_species_react_rhs(app, gkns, &gkns->react_neut, fin, rhs);

  // Add the explicit ionization frequency ne*<sigma v>_iz to the CFL rate.
  if (gkns->ionization_react_idx >= 0) {
    struct gk_react *react = &gkns->react_neut;
    for (int i = 0; i < react->num_react; ++i) {
      if (react->react_id[i] != GKYL_REACT_IZ) {
        continue;
      }
      struct gk_species *electron = &app->species[react->elc_idx[i]];
      gkyl_dg_mul_op_range(
        &app->basis, 0, gkns->ionization_rate, 0, electron->lte.moms.marr, 0,
        react->coeff_react[i], &gkns->local
      );
      gkyl_dg_calc_average_range(
        &app->basis, 0, gkns->ionization_cflrate, 0, gkns->ionization_rate, gkns->local
      );
      gkyl_array_accumulate_range(
        gkns->cflrate, 1.0, gkns->ionization_cflrate, &gkns->local
      );
    }
  }

  // Compute volume-integrated reactions in sca.
  gk_neut_species_scaling_rhs(app, gkns, &gkns->sca, fin, rhs);

  app->stat.n_neut_species_omega_cfl += 1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(gkns->omega_cfl, gkns->cflrate, GKYL_MAX, &gkns->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, gkns->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_cfl_ho[0] = gkns->omega_cfl[0];
  }

  omega_cfl = omega_cfl_ho[0];

  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl / omega_cfl;
}

static double
gk_neut_species_fluid_rhs_implicit_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
)
{
  double omega_cfl = 1 / DBL_MAX;
  gkyl_array_clear(gkns->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // No implicit terms yet.

  gkyl_array_accumulate(gkyl_array_scale(rhs, dt), 1.0, fin);

  app->stat.n_neut_species_omega_cfl += 1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(gkns->omega_cfl, gkns->cflrate, GKYL_MAX, &gkns->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, gkns->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_cfl_ho[0] = gkns->omega_cfl[0];
  }
  omega_cfl = omega_cfl_ho[0];

  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl / omega_cfl;
}

static void
gk_neut_species_fluid_release_dynamic(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns
)
{
  // Release memory allocated for dynamic neutrals.
  gkyl_array_release(gkns->cflrate);

  if (gkns->has_diffusion) {
    gkyl_array_release(gkns->diffusion_coeff);
    gkyl_array_release(gkns->diffusion_tensor);
    gkyl_array_release(gkns->diffusion_geom_factor);
    gkyl_array_release(gkns->diffusion_density);
    if (gkns->diffusion_bc_density) {
      gkyl_array_release(gkns->diffusion_bc_density);
    }
    if (gkns->diffusion_bc_density_tmp) {
      gkyl_array_release(gkns->diffusion_bc_density_tmp);
    }
    gkyl_array_release(gkns->diffusion_moment_ratio);
    gkyl_dg_bin_op_mem_release(gkns->diffusion_div_mem);
    gkyl_dg_updater_conf_diffusion_release(gkns->diffusion_slvr);
  }
  if (gkns->ionization_react_idx >= 0) {
    gkyl_array_release(gkns->ionization_rate);
    gkyl_array_release(gkns->ionization_cflrate);
  }

  if (app->use_gpu) {
    gkyl_cu_free(gkns->omega_cfl);
  } else {
    gkyl_free(gkns->omega_cfl);
  }

  // Release integrated mom data.
  gk_neut_species_moment_release(app, &gkns->integ_moms);

  // Release integrated mom diag data.
  gkyl_dynvec_release(gkns->integ_diag);

  if (app->use_gpu) {
    gkyl_cu_free(gkns->red_integ_diag);
    gkyl_cu_free(gkns->red_integ_diag_global);
  } else {
    gkyl_free(gkns->red_integ_diag);
    gkyl_free(gkns->red_integ_diag_global);
  }
}

static void
gk_neut_species_fluid_release(const gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns)
{
  // Release resources for fluid neutral species.
  gkyl_msgpack_map_elem_release(gkns->io_meta_basic_len, gkns->io_meta_basic);
  gkyl_msgpack_map_elem_release(gkns->io_meta_conf_len, gkns->io_meta_conf);

  gkyl_array_release(gkns->f);
  gkyl_array_release(gkns->f1);
  gkyl_array_release(gkns->fnew);
  gkyl_array_release(gkns->f_host);

  if (gkns->info.init_from_file.type == 0) {
    gk_neut_species_projection_release(app, &gkns->proj_init);
  }
  gkyl_comm_release(gkns->comm);

  for (int i = 0; i < gkns->info.num_diag_moments; ++i) {
    gk_neut_species_moment_release(app, &gkns->moms[i]);
  }
  gkyl_free(gkns->moms);

  gk_neut_species_bgk_release(app, &gkns->bgk);

  gk_neut_species_positivity_release(app, &gkns->positivity);

  gk_neut_species_react_release(app, &gkns->react_neut);

  // Free boundary flux memory.
  gk_neut_species_bflux_release(app, gkns, &gkns->bflux);

  gk_neut_species_lte_release(app, &gkns->lte);

  // Free memory for the object that scales the species according to a balance
  // between recycling and reactions.
  gk_neut_species_scaling_release(app, &gkns->sca);

  gkns->release_is_static_func(app, gkns);
}

static void
gk_neut_species_fluid_init_dynamic(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns
)
{
  int cdim = app->cdim;

  // Allocate additional moment arrays for time stepping.
  gkns->f1 = mkarr(app->use_gpu, gkns->f->ncomp, gkns->f->size);
  gkns->fnew = mkarr(app->use_gpu, gkns->f->ncomp, gkns->f->size);

  // Allocate cflrate (scalar array).
  gkns->cflrate = mkarr(app->use_gpu, 1, gkns->local_ext.volume);

  gkns->omega_cfl = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));

  // Initialize diffusion of neutral density (component 0) only. The same
  // scalar coefficient is used in every configuration-space direction.
  gkns->use_reaction_rate_diffusion = gkns->info.diffusion.use_reaction_rates;
  gkns->has_diffusion =
    gkns->info.diffusion.D > 0.0 || gkns->use_reaction_rate_diffusion;
  gkns->ionization_react_idx = -1;
  gkns->diffusion_cx_react_idx = -1;
  gkns->diffusion_coeff = 0;
  gkns->diffusion_tensor = 0;
  gkns->diffusion_geom_factor = 0;
  gkns->diffusion_density = 0;
  gkns->diffusion_bc_density = 0;
  gkns->diffusion_bc_density_tmp = 0;
  gkns->diffusion_moment_ratio = 0;
  gkns->diffusion_div_mem = 0;
  gkns->diffusion_slvr = 0;
  gkns->diffusion_recycling_num_species[0] = 0;
  gkns->diffusion_recycling_num_species[1] = 0;
  gkns->ionization_rate = 0;
  gkns->ionization_cflrate = 0;

  for (int i = 0; i < gkns->info.react_neut.num_react; ++i) {
    enum gkyl_react_id react_id = gkns->info.react_neut.react_type[i].react_id;
    if (react_id == GKYL_REACT_IZ && gkns->ionization_react_idx < 0) {
      gkns->ionization_react_idx = i;
    } else if (react_id == GKYL_REACT_CX && gkns->diffusion_cx_react_idx < 0) {
      gkns->diffusion_cx_react_idx = i;
    }
  }
  if (gkns->ionization_react_idx >= 0) {
    gkns->ionization_rate =
      mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkns->ionization_cflrate = mkarr(app->use_gpu, 1, app->local_ext.volume);
  }
  if (gkns->has_diffusion) {
    // Configuration-space restriction of the single-integration-by-parts plus
    // recovery-DG anomalous-diffusion scheme. Generated CPU kernels currently
    // cover 1x-3x serendipity p1.
    assert(cdim > 0 && cdim <= 3);
    assert(app->basis.poly_order == 1);
    assert(!app->use_gpu);

    gkns->diffusion_coeff =
      mkarr(false, app->basis.num_basis, app->local_ext.volume);
    gkns->diffusion_tensor =
      mkarr(false, cdim * cdim * app->basis.num_basis, app->local_ext.volume);
    gkns->diffusion_geom_factor =
      mkarr(false, app->basis.num_basis, app->local_ext.volume);
    gkns->diffusion_density =
      mkarr(false, app->basis.num_basis, app->local_ext.volume);

    const int parallel_dir = cdim - 1;
    const bool recycle_lower =
      gkns->lower_bc[parallel_dir].type == GKYL_BC_GK_SPECIES_RECYCLE;
    const bool recycle_upper =
      gkns->upper_bc[parallel_dir].type == GKYL_BC_GK_SPECIES_RECYCLE;
    if (recycle_lower || recycle_upper) {
      gkns->diffusion_bc_density =
        mkarr(false, app->basis.num_basis, app->local_ext.volume);
      gkns->diffusion_bc_density_tmp =
        mkarr(false, app->basis.num_basis, app->local_ext.volume);
      gkyl_array_clear(gkns->diffusion_bc_density, 0.0);
      gkyl_array_clear(gkns->diffusion_bc_density_tmp, 0.0);
    }
    for (int edge_idx = 0; edge_idx < 2; ++edge_idx) {
      const struct gkyl_gyrokinetic_bc *bc = edge_idx == 0
        ? &gkns->lower_bc[parallel_dir]
        : &gkns->upper_bc[parallel_dir];
      if (bc->type != GKYL_BC_GK_SPECIES_RECYCLE) {
        continue;
      }
      assert(bc->emission.emission_temp > 0.0);
      assert(bc->emission.num_species > 0);
      gkns->diffusion_recycling_num_species[edge_idx] = bc->emission.num_species;
      for (int i = 0; i < bc->emission.num_species; ++i) {
        gkns->diffusion_recycling_ion_idx[edge_idx][i] =
          gk_find_species_idx(app, bc->emission.in_species[i]);
        assert(gkns->diffusion_recycling_ion_idx[edge_idx][i] >= 0);
      }
    }
    gkyl_array_clear(gkns->diffusion_coeff, 0.0);
    gkyl_array_clear(gkns->diffusion_tensor, 0.0);
    gkyl_array_clear(gkns->diffusion_geom_factor, 0.0);
    gkyl_array_clear(gkns->diffusion_density, 0.0);
    if (!gkns->use_reaction_rate_diffusion) {
      // Modal constant coefficient in cdim dimensions.
      gkyl_array_shiftc(
        gkns->diffusion_coeff, pow(sqrt(2.0), cdim) * gkns->info.diffusion.D, 0
      );
    }

    if (gkns->use_reaction_rate_diffusion) {
      assert(gkns->ionization_react_idx >= 0);
      assert(gkns->diffusion_cx_react_idx >= 0);
    }

    gkns->diffusion_moment_ratio =
      mkarr(app->use_gpu, gkns->f->ncomp, gkns->f->size);
    gkns->diffusion_div_mem = app->use_gpu
      ? gkyl_dg_bin_op_mem_cu_dev_new(app->local.volume, app->basis.num_basis)
      : gkyl_dg_bin_op_mem_new(app->local.volume, app->basis.num_basis);

    int diffusion_zero_flux_flags[2 * GKYL_MAX_CDIM] = {0};
    for (int d = 0; d < cdim; ++d) {
      diffusion_zero_flux_flags[d] =
        diffusion_zero_flux_flags[d + cdim] = gkns->bc_is_np[d];
    }
    if (recycle_lower ||
        gkns->info.diffusion.lower_bc_type == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET) {
      diffusion_zero_flux_flags[parallel_dir] = 0;
    }
    if (recycle_upper ||
        gkns->info.diffusion.upper_bc_type == GKYL_NEUT_FLUID_DIFFUSION_DIRICHLET) {
      diffusion_zero_flux_flags[parallel_dir + cdim] = 0;
    }
    gkns->diffusion_slvr = gkyl_dg_updater_conf_diffusion_new(
      &app->grid, &app->basis, &app->local_ext, gkns->diffusion_tensor,
      app->gk_geom->geo_int.jacobgeo_inv, diffusion_zero_flux_flags, app->use_gpu
    );
  }

  // Allocate data for integrated moments.
  gk_neut_species_moment_init(app, gkns, &gkns->integ_moms, GKYL_F_MOMENT_M0M1M2, true);

  // Allocate data for integrated diagnostics.
  if (app->use_gpu) {
    gkns->red_integ_diag = gkyl_cu_malloc(sizeof(double[gkns->integ_moms.num_mom]));
    gkns->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[gkns->integ_moms.num_mom]));
  } else {
    gkns->red_integ_diag = gkyl_malloc(sizeof(double[gkns->integ_moms.num_mom]));
    gkns->red_integ_diag_global = gkyl_malloc(sizeof(double[gkns->integ_moms.num_mom]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments.
  gkns->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, gkns->integ_moms.num_mom);
  gkns->is_first_integ_write_call = true;

  // Set function pointers
  gkns->rhs_func = gk_neut_species_fluid_rhs_dynamic;
  gkns->rhs_implicit_func = gk_neut_species_fluid_rhs_implicit_dynamic;
  gkns->bc_func = gk_neut_species_fluid_apply_bc_dynamic;
  gkns->release_func = gk_neut_species_fluid_release;
  gkns->release_is_static_func = gk_neut_species_fluid_release_dynamic;
  gkns->step_f_func = gk_neut_species_step_f_dynamic;
  gkns->combine_func = gk_neut_species_combine_dynamic;
  gkns->copy_func = gk_neut_species_copy_range_dynamic;
  gkns->write_func = gk_neut_species_write_dynamic;
  gkns->write_mom_func =
    gk_neut_species_write_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  gkns->calc_integrated_mom_func =
    gk_neut_species_calc_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  gkns->write_integrated_mom_func =
    gk_neut_species_write_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  gkns->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

static void
gk_neut_species_fluid_init_static(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns
)
{
  // Set pointers for RK methods.
  gkns->f1 = gkyl_array_acquire(gkns->f);
  gkns->fnew = gkyl_array_acquire(gkns->f);

  // Set function pointers
  gkns->rhs_func = gk_neut_species_rhs_static;
  gkns->rhs_implicit_func = gk_neut_species_rhs_implicit_static;
  gkns->bc_func = gk_neut_species_apply_bc_static;
  gkns->release_func = gk_neut_species_fluid_release;
  gkns->release_is_static_func = gk_neut_species_release_static;
  gkns->step_f_func = gk_neut_species_step_f_static;
  gkns->combine_func = gk_neut_species_combine_static;
  gkns->copy_func = gk_neut_species_copy_range_static;
  gkns->write_func = gk_neut_species_write_init_only;
  gkns->write_mom_func = gk_neut_species_write_mom_init_only;
  gkns->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_static;
  gkns->write_integrated_mom_func = gk_neut_species_write_integrated_mom_static;
  gkns->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

void
gk_neut_species_fluid_init(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns
)
{
  gkns->is_fluid = true; // Fluid neutrals.
  assert(
    gkns->info.vdim == 0
  ); // Ensure user provided vdim=0 in input file, or didn't provide it at all.

  gkns->model_id = GKYL_MODEL_DEFAULT;
  gkns->field_id = GKYL_FIELD_NULL;

  int cdim = app->cdim;

  // Number of moments depends on eqn_type.
  gkns->num_moments = 5; // rho, rho*ux, rho*uy, rho*uy, totE

  // Use the same basis as conf-space.
  gkns->basis = app->basis;
  gkns->basis_on_dev = app->basis_on_dev;

  // Use the same grid as conf-space.
  gkns->grid = app->grid;
  gkns->global_ext = app->global_ext;
  gkns->global = app->global;

  // Use the same communicator as conf-space.
  gkns->comm = gkyl_comm_acquire(app->comm);

  // Use the same local range as conf-space.
  gkns->local_ext = app->local_ext;
  gkns->local = app->local;

  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // add the parallel direction in case TS BCs are needed.
  gkns->num_periodic_dir = app->num_periodic_dir;
  for (int d = 0; d < gkns->num_periodic_dir; ++d) {
    gkns->periodic_dirs[d] = app->periodic_dirs[d];
  }

  // In an LCFS-spanning domain, populate parallel ghosts periodically first.
  // The diffusion RHS later overwrites only the SOL subset with recycled or
  // fixed Dirichlet data, leaving the confined subset periodic.
  if (app->gk_geom->has_LCFS) {
    gkns->periodic_dirs[gkns->num_periodic_dir] = app->cdim - 1;
    gkns->num_periodic_dir += 1;
  }

  for (int d = 0; d < app->cdim; ++d) {
    gkns->bc_is_np[d] = true;
  }
  for (int d = 0; d < gkns->num_periodic_dir; ++d) {
    gkns->bc_is_np[gkns->periodic_dirs[d]] = false;
  }

  // Store the BCs from the input file.
  for (int d = 0; d < app->cdim; ++d) {
    struct gkyl_gyrokinetic_bc *bc_lo =
      gk_fetch_bc_with_dir_edge(gkns->info.bcs, 2 * app->cdim, d, GKYL_LOWER_EDGE);
    if (bc_lo != 0) {
      gkns->lower_bc[d] = *bc_lo;
    } else {
      gkns->lower_bc[d].type = GKYL_BC_GK_SKIP;
    }

    struct gkyl_gyrokinetic_bc *bc_up =
      gk_fetch_bc_with_dir_edge(gkns->info.bcs, 2 * app->cdim, d, GKYL_UPPER_EDGE);
    if (bc_up != 0) {
      gkns->upper_bc[d] = *bc_up;
    } else {
      gkns->upper_bc[d].type = GKYL_BC_GK_SKIP;
    }
  }

  // Species properties metadata.
  struct gkyl_msgpack_map_elem io_meta_sprop[] = {
    {.key = "mass", .elem_type = GKYL_MP_DOUBLE, .dval = gkns->info.mass},
    {.key = "charge", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "gas_gamma", .elem_type = GKYL_MP_DOUBLE, .dval = gkns->info.gas_gamma},
    {.key = "vdim", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = gkns->info.vdim}
  };

  // Metadata for integrated quantities.
  const struct gkyl_msgpack_map_elem *io_meta_basic_union[] = {app->io_meta_basic, io_meta_sprop};
  int io_meta_basic_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0])
  };
  gkns->io_meta_basic = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_basic_union) / sizeof(io_meta_basic_union[0]), io_meta_basic_union_len,
    io_meta_basic_union, &gkns->io_meta_basic_len
  );

  // Metadata for conf-space quantities.
  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    {.key = "value_form", .elem_type = GKYL_MP_STRING, .cval = "modal"},
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = app->basis.poly_order},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = app->basis.id},
    {.key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0}
  };
  const struct gkyl_msgpack_map_elem *io_meta_conf_union[] = {
    app->io_meta_basic, io_meta_sprop, io_meta_conf
  };
  int io_meta_conf_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0]),
    sizeof(io_meta_conf) / sizeof(io_meta_conf[0])
  };
  gkns->io_meta_conf = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_conf_union) / sizeof(io_meta_conf_union[0]), io_meta_conf_union_len,
    io_meta_conf_union, &gkns->io_meta_conf_len
  );

  // Metadata for phase-space quantities.
  gkns->io_meta_phase = gkns->io_meta_conf;
  gkns->io_meta_phase_len = gkns->io_meta_conf_len;

  // Allocate distribution function array for initialization and I/O.
  gkns->f = mkarr(app->use_gpu, gkns->num_moments * gkns->basis.num_basis, gkns->local_ext.volume);
  gkns->f_host = app->use_gpu ? mkarr(false, gkns->f->ncomp, gkns->f->size) :
                                gkyl_array_acquire(gkns->f);

  // Create skin/ghost ranges fir applying BCs. Only used for dynamic neutrals but included here to avoid
  // code duplication since the "ghost" array is needed.
  int ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    ghost[d] = 1;
  }

  // Create skin/ghost ranges.
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &gkns->local_lower_skin[dir], &gkns->local_lower_ghost[dir], dir, GKYL_LOWER_EDGE,
      &gkns->local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &gkns->local_upper_skin[dir], &gkns->local_upper_ghost[dir], dir, GKYL_UPPER_EDGE,
      &gkns->local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &gkns->global_lower_skin[dir], &gkns->global_lower_ghost[dir], dir, GKYL_LOWER_EDGE,
      &gkns->global_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &gkns->global_upper_skin[dir], &gkns->global_upper_ghost[dir], dir, GKYL_UPPER_EDGE,
      &gkns->local_ext, ghost
    );
  }

  // Initialize projection routine for initial conditions.
  if (gkns->info.init_from_file.type == 0) {
    gk_neut_species_projection_init(app, gkns, gkns->info.projection, &gkns->proj_init);
  }

  // Allocate objects for computing diagnostic moments.
  int ndm = gkns->info.num_diag_moments;
  gkns->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m = 0; m < ndm; ++m) {
    gk_neut_species_moment_init(app, gkns, &gkns->moms[m], gkns->info.diag_moments[m], false);
  }

  // Initialize boundary fluxes.
  gkns->bflux = (struct gk_boundary_fluxes){};
  // Additional bflux moments to step in time.
  struct gkyl_phase_diagnostics_inp add_bflux_moms_inp = (struct gkyl_phase_diagnostics_inp){};
  // Set the operation type for the bflux app.
  enum gkyl_species_bflux_type bflux_type = GK_SPECIES_BFLUX_NONE;
  gk_neut_species_bflux_init(app, gkns, &gkns->bflux, bflux_type, add_bflux_moms_inp);

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine
  // Projection routine optionally corrects all the Maxwellian/LTE moments
  // This routine is utilized by both reactions and BGK collisions
  gkns->lte = (struct gk_lte){};
  struct correct_all_moms_inp corr_inp = {
    .correct_all_moms = false,
    .max_iter = 0,
    .iter_eps = 10,
    .use_last_converged = false,
  };
  gk_neut_species_lte_init(app, gkns, &gkns->lte, corr_inp);

  // Initialize the object that scales the species according to a balance
  // between recycling and reactions.
  gkns->sca = (struct gk_scaling){};
  gk_neut_species_scaling_init(app, gkns, &gkns->sca);

  // Initialize BGK collisions with null type (not applicable to fluids).
  gkns->bgk = (struct gk_bgk_collisions){};
  gkns->info.collisions.collision_id = 0;
  gk_neut_species_bgk_init(app, gkns, &gkns->bgk);

  // Initialize positivity enforcing operator with null type (NYI for fluids).
  gkns->positivity = (struct gk_positivity){};
  gk_neut_species_positivity_init(app, gkns, &gkns->positivity);

  // Initialize reactions with charged species. Fluid neutrals use the
  // ionization reaction as a density loss; CX supplies the diffusion rate.
  gkns->react_neut = (struct gk_react){};
  gk_neut_species_react_init(app, gkns, gkns->info.react_neut, &gkns->react_neut);

  gkns->src = (struct gk_source){};
  if (!gkns->info.is_static) {
    gk_neut_species_fluid_init_dynamic(gk, app, gkns);
  } else {
    gk_neut_species_fluid_init_static(gk, app, gkns);
  }
}
