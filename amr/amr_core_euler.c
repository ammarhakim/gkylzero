#include <gkyl_amr_core.h>
#include <gkyl_amr_block_priv.h>
#include <gkyl_amr_block_coupled_priv.h>
#include <gkyl_amr_patch_priv.h>

void
euler1d_run_single(int argc, char **argv, struct euler1d_single_init* init)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int base_Nx = init->base_Nx;
  int ref_factor = init->ref_factor;

  double coarse_x1 = init->coarse_x1;
  double coarse_x2 = init->coarse_x2;

  double refined_x1 = init->refined_x1;
  double refined_x2 = init->refined_x2;

  evalf_t eval = init->eval;
  double gas_gamma = init->gas_gamma;

  char euler_output[32];
  strcpy(euler_output, init->euler_output);

  bool low_order_flux = init->low_order_flux;
  int num_frames = init->num_frames;
  double cfl_frac = init->cfl_frac;
  double t_end = init->t_end;

  int ndim = 1;
  int num_patches = 3;
  int Nx = base_Nx;

  struct euler_patch_data coarse_pdata[num_patches];
  struct gkyl_job_pool *coarse_job_pool = gkyl_thread_pool_new(app_args.num_threads);

#ifdef AMR_DEBUG
  gkyl_rect_grid_init(&coarse_pdata[0].grid, 1, (double []) { refined_x1 }, (double []) { refined_x2 }, (int []) { Nx } );
#else
  gkyl_rect_grid_init(&coarse_pdata[0].grid, 1, (double []) { refined_x1 }, (double []) { refined_x2 }, (int []) { Nx * ref_factor } );
#endif
  
  gkyl_rect_grid_init(&coarse_pdata[1].grid, 1, (double []) { coarse_x1 }, (double []) { refined_x1 }, (int []) { Nx } );
  gkyl_rect_grid_init(&coarse_pdata[2].grid, 1, (double []) { refined_x2 }, (double []) { coarse_x2 }, (int []) { Nx } );

#ifdef AMR_DEBUG
  struct euler_patch_data fine_pdata[num_patches];
  struct gkyl_job_pool *fine_job_pool = gkyl_thread_pool_new(app_args.num_threads);

  gkyl_rect_grid_init(&fine_pdata[0].grid, 1, (double []) { refined_x1 }, (double []) { refined_x2 }, (int []) { Nx * ref_factor } );
  gkyl_rect_grid_init(&fine_pdata[1].grid, 1, (double []) { coarse_x1 }, (double []) { refined_x1 }, (int []) { Nx * ref_factor } );
  gkyl_rect_grid_init(&fine_pdata[2].grid, 1, (double []) { refined_x2 }, (double []) { coarse_x2 }, (int []) { Nx * ref_factor } );
#endif
  
  for (int i = 0; i < num_patches; i++) {
    coarse_pdata[i].fv_proj = gkyl_fv_proj_new(&coarse_pdata[i].grid, 1, 5, eval, 0);

#ifdef AMR_DEBUG
    fine_pdata[i].fv_proj = gkyl_fv_proj_new(&fine_pdata[i].grid, 1, 5, eval, 0);
#endif
  }

  for (int i = 0; i < num_patches; i++) {
    gkyl_create_grid_ranges(&coarse_pdata[i].grid, (int []) { 2 }, &coarse_pdata[i].ext_range, &coarse_pdata[i].range );
    coarse_pdata[i].geom = gkyl_wave_geom_new(&coarse_pdata[i].grid, &coarse_pdata[i].ext_range, 0, 0, false);

#ifdef AMR_DEBUG
    gkyl_create_grid_ranges(&fine_pdata[i].grid, (int []) { 2 }, &fine_pdata[i].ext_range, &fine_pdata[i].range);
    fine_pdata[i].geom = gkyl_wave_geom_new(&fine_pdata[i].grid, &fine_pdata[i].ext_range, 0, 0, false);
#endif
  }

  for (int i = 0; i < num_patches; i++) {
    if (low_order_flux) {
      struct gkyl_wv_euler_inp inp = {
        .gas_gamma = gas_gamma,
        .rp_type = WV_EULER_RP_HLL,
        .use_gpu = app_args.use_gpu,
      };
      coarse_pdata[i].euler = gkyl_wv_euler_inew(&inp);
    }
    else {
      coarse_pdata[i].euler = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    }

    coarse_pdata[i].slvr[0] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
        .grid = &coarse_pdata[i].grid,
        .equation = coarse_pdata[i].euler,
        .limiter = GKYL_MONOTONIZED_CENTERED,
        .num_up_dirs = 1,
        .update_dirs = { 0 },
        .cfl = cfl_frac,
        .geom = coarse_pdata[i].geom,
      }
    );

#ifdef AMR_DEBUG
    if (low_order_flux) {
      struct gkyl_wv_euler_inp inp = {
        .gas_gamma = gas_gamma,
        .rp_type = WV_EULER_RP_HLL,
        .use_gpu = app_args.use_gpu,
      };
      fine_pdata[i].euler = gkyl_wv_euler_inew(&inp);
    }
    else {
      fine_pdata[i].euler = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    }

    fine_pdata[i].slvr[0] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
        .grid = &fine_pdata[i].grid,
        .equation = fine_pdata[i].euler,
        .limiter = GKYL_MONOTONIZED_CENTERED,
        .num_up_dirs = 1,
        .update_dirs = { 0 },
        .cfl = cfl_frac,
        .geom = fine_pdata[i].geom,
      }
    );
#endif
  }

  struct gkyl_block_topo *ptopo = create_patch_topo();

  for (int i = 0; i < num_patches; i++) {
    euler_patch_bc_updaters_init(&coarse_pdata[i], &ptopo->conn[i]);

#ifdef AMR_DEBUG
    euler_patch_bc_updaters_init(&fine_pdata[i], &ptopo->conn[i]);
#endif
  }

  for (int i = 0; i < num_patches; i++) {
    coarse_pdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 5, coarse_pdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      coarse_pdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 5, coarse_pdata[i].ext_range.volume);
    }
#ifdef AMR_DEBUG
    fine_pdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 5, fine_pdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      fine_pdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 5, fine_pdata[i].ext_range.volume);
    }
#endif
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_patches; i++) {
    gkyl_job_pool_add_work(coarse_job_pool, euler_init_job_func_patch, &coarse_pdata[i]);

#ifdef AMR_DEBUG
    gkyl_job_pool_add_work(fine_job_pool, euler_init_job_func_patch, &fine_pdata[i]);
#endif
  }
  gkyl_job_pool_wait(coarse_job_pool);

#ifdef AMR_DEBUG
  gkyl_job_pool_wait(fine_job_pool);
#endif
#else
  for (int i = 0; i < num_patches; i++) {
    euler_init_job_func_patch(&coarse_pdata[i]);

#ifdef AMR_DEBUG
    euler_init_job_func_patch(&fine_pdata[i]);
#endif
  }
#endif

#ifdef AMR_DEBUG
  char coarse0[64];
  snprintf(coarse0, 64, "%s_coarse_0", euler_output);
  euler_write_sol_patch(coarse0, num_patches, coarse_pdata);

  char fine0[64];
  snprintf(fine0, 64, "%s_fine_0", euler_output);
  euler_write_sol_patch(fine0, num_patches, fine_pdata);

  char fine0p0[64];
  char coarse0p0[64];
  char p0[64];
  snprintf(fine0p0, 64, "%s_fine_0_p0.gkyl", euler_output);
  snprintf(coarse0p0, 64, "%s_coarse_0_p0.gkyl", euler_output);
  snprintf(p0, 64, "%s_0_p0.gkyl", euler_output);
  rename(fine0p0, p0);
  remove(coarse0p0);

  for (int i = 1; i < 3; i++) {
    char buf_old[64];
    char buf_new[64];
    char buf_del[64];

    snprintf(buf_old, 64, "%s_coarse_0_p%d.gkyl", euler_output, i);
    snprintf(buf_new, 64, "%s_0_p%d.gkyl", euler_output, i);
    snprintf(buf_del, 64, "%s_fine_0_p%d.gkyl", euler_output, i);

    rename(buf_old, buf_new);
    remove(buf_del);
  }
#else
  char amr0[64];
  snprintf(amr0, 64, "%s_0", euler_output);
  euler_write_sol_patch(amr0, num_patches, coarse_pdata);
#endif

  double coarse_t_curr = 0.0;
  double fine_t_curr = 0.0;
  double coarse_dt = euler_max_dt_patch(num_patches, coarse_pdata);

#ifdef AMR_DEBUG
  double fine_dt = euler_max_dt_patch(num_patches, fine_pdata);
#else
  double fine_dt = (1.0 / ref_factor) * coarse_dt;
#endif

  struct sim_stats stats = { };

  struct timespec tm_start = gkyl_wall_clock();

  long coarse_step = 1;
  long num_steps = app_args.num_steps;

  double io_trigger = t_end / num_frames;

  while ((coarse_t_curr < t_end) && (coarse_step <= num_steps)) {
    printf("Taking coarse (level 0) time-step %ld at t = %g; ", coarse_step, coarse_t_curr);
    struct gkyl_update_status coarse_status = euler_update_patch(coarse_job_pool, ptopo, coarse_pdata, coarse_t_curr, coarse_dt, &stats);
    printf(" dt = %g\n", coarse_status.dt_actual);

    if (!coarse_status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }

    for (long fine_step = 1; fine_step < ref_factor + 1; fine_step++) {
#ifdef AMR_DEBUG
      printf("   Taking fine (level 1) time-step %ld at t = %g; ", fine_step, fine_t_curr);
      struct gkyl_update_status fine_status = euler_update_patch(fine_job_pool, ptopo, fine_pdata, fine_t_curr, fine_dt, &stats);
      printf(" dt = %g\n", fine_status.dt_actual);

      if (!fine_status.success) {
        printf("** Update method failed! Aborting simulation ....\n");
        break;
      }

      fine_t_curr += fine_status.dt_actual;
      fine_dt = fine_status.dt_suggested;
#else
      printf("   Taking fine (level 1) time-step %ld at t = %g", fine_step, fine_t_curr);
      printf(" dt = %g\n", (1.0 / ref_factor) * coarse_status.dt_actual);

      fine_t_curr += (1.0 / ref_factor) * coarse_status.dt_actual;
      fine_dt = (1.0 / ref_factor) * coarse_status.dt_suggested;
#endif
    }

    for (int i = 1; i < num_frames; i++) {
      if (coarse_t_curr < (i * io_trigger) && (coarse_t_curr + coarse_status.dt_actual) > (i * io_trigger)) {
#ifdef AMR_DEBUG
      char buf_coarse[64];
      char buf_fine[64];

      snprintf(buf_coarse, 64, "%s_coarse_%d", euler_output, i);
      snprintf(buf_fine, 64, "%s_fine_%d", euler_output, i);

      euler_write_sol_patch(buf_coarse, num_patches, coarse_pdata);
      euler_write_sol_patch(buf_fine, num_patches, fine_pdata);

      char buf_fine_old[64];
      char buf_fine_new[64];
      char buf_coarse_old[64];

      snprintf(buf_fine_old, 64, "%s_fine_%d_p0.gkyl", euler_output, i);
      snprintf(buf_fine_new, 64, "%s_%d_p0.gkyl", euler_output, i);
      snprintf(buf_coarse_old, 64, "%s_coarse_%d_p0.gkyl", euler_output, i);

      rename(buf_fine_old, buf_fine_new);
      remove(buf_coarse_old);

      for (int j = 1; j < 3; j++) {
        char buf_old[64];
        char buf_new[64];
        char buf_del[64];

        snprintf(buf_old, 64, "%s_coarse_%d_p%d.gkyl", euler_output, i, j);
        snprintf(buf_new, 64, "%s_%d_p%d.gkyl", euler_output, i, j);
        snprintf(buf_del, 64, "%s_fine_%d_p%d.gkyl", euler_output, i, j);

        rename(buf_old, buf_new);
        remove(buf_del);
      }
#else
      char buf[64];
      snprintf(buf, 64, "%s_%d", euler_output, i);

      euler_write_sol_patch(buf, num_patches, coarse_pdata);
#endif
      }
    }

    coarse_t_curr += coarse_status.dt_actual;
    coarse_dt = coarse_status.dt_suggested;

    coarse_step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

#ifdef AMR_DEBUG
  char buf_coarse[64];
  char buf_fine[64];

  snprintf(buf_coarse, 64, "%s_coarse_%d", euler_output, num_frames);
  snprintf(buf_fine, 64, "%s_fine_%d", euler_output, num_frames);

  euler_write_sol_patch(buf_coarse, num_patches, coarse_pdata);
  euler_write_sol_patch(buf_fine, num_patches, fine_pdata);

  char buf_fine_old[64];
  char buf_fine_new[64];
  char buf_coarse_old[64];

  snprintf(buf_fine_old, 64, "%s_fine_%d_p0.gkyl", euler_output, num_frames);
  snprintf(buf_fine_new, 64, "%s_%d_p0.gkyl", euler_output, num_frames);
  snprintf(buf_coarse_old, 64, "%s_coarse_%d_p0.gkyl", euler_output,num_frames);

  rename(buf_fine_old, buf_fine_new);
  remove(buf_coarse_old);

  for (int i = 1; i < 3; i++) {
    char buf_old[64];
    char buf_new[64];
    char buf_del[64];

    snprintf(buf_old, 64, "%s_coarse_%d_p%d.gkyl", euler_output, num_frames, i);
    snprintf(buf_new, 64, "%s_%d_p%d.gkyl", euler_output, num_frames, i);
    snprintf(buf_del, 64, "%s_fine_%d_p%d.gkyl", euler_output, num_frames, i);

    rename(buf_old, buf_new);
    remove(buf_del);
  }
#else
  char buf[64];
  snprintf(buf, 64, "%s_%d", euler_output, num_frames);

  euler_write_sol_patch(buf, num_patches, coarse_pdata);
#endif

  printf("Total run-time: %g. Failed steps: %d\n", tm_total_sec, stats.nfail);

  for (int i = 0; i < num_patches; i++) {
    gkyl_fv_proj_release(coarse_pdata[i].fv_proj);
    gkyl_wv_eqn_release(coarse_pdata[i].euler);
    euler_patch_bc_updaters_release(&coarse_pdata[i]);
    gkyl_wave_geom_release(coarse_pdata[i].geom);

#ifdef AMR_DEBUG
    gkyl_fv_proj_release(fine_pdata[i].fv_proj);
    gkyl_wv_eqn_release(fine_pdata[i].euler);
    euler_patch_bc_updaters_release(&fine_pdata[i]);
    gkyl_wave_geom_release(fine_pdata[i].geom);
#endif

    gkyl_wave_prop_release(coarse_pdata[i].slvr[0]);

#ifdef AMR_DEBUG
    gkyl_wave_prop_release(fine_pdata[i].slvr[0]);
#endif
    
    gkyl_array_release(coarse_pdata[i].fdup);
#ifdef AMR_DEBUG
    gkyl_array_release(fine_pdata[i].fdup);
#endif
    
    gkyl_array_release(coarse_pdata[i].f[0]);

#ifdef AMR_DEBUG
    gkyl_array_release(fine_pdata[i].f[0]);
#endif
  }

  gkyl_block_topo_release(ptopo);
  gkyl_job_pool_release(coarse_job_pool);
#ifdef AMR_DEBUG
  gkyl_job_pool_release(fine_job_pool);
#endif
}

void
euler2d_run_single(int argc, char **argv, struct euler2d_single_init* init)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int base_Nx = init->base_Nx;
  int base_Ny = init->base_Ny;
  int ref_factor = init->ref_factor;

  double coarse_x1 = init->coarse_x1;
  double coarse_y1 = init->coarse_y1;
  double coarse_x2 = init->coarse_x2;
  double coarse_y2 = init->coarse_y2;

  double refined_x1 = init->refined_x1;
  double refined_y1 = init->refined_y1;
  double refined_x2 = init->refined_x2;
  double refined_y2 = init->refined_y2;

  evalf_t eval = init->eval;
  double gas_gamma = init->gas_gamma;

  char euler_output[32];
  strcpy(euler_output, init->euler_output);
  
  bool low_order_flux = init->low_order_flux;
  int num_frames = init->num_frames;
  double cfl_frac = init->cfl_frac;
  double t_end = init->t_end;

  int ndim = 2;
  int num_blocks = 9;
  int Nx = base_Nx;
  int Ny = base_Ny;

  struct euler_block_data coarse_bdata[num_blocks];
  struct gkyl_job_pool *coarse_job_pool = gkyl_thread_pool_new(app_args.num_threads);

#ifdef AMR_DEBUG
  gkyl_rect_grid_init(&coarse_bdata[0].grid, 2, (double []) { refined_x1, refined_y1 }, (double []) { refined_x2, refined_y2 },
    (int []) { Nx, Ny });
#else
  gkyl_rect_grid_init(&coarse_bdata[0].grid, 2, (double []) { refined_x1, refined_y1 }, (double []) { refined_x2, refined_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
#endif

  gkyl_rect_grid_init(&coarse_bdata[1].grid, 2, (double []) { coarse_x1, refined_y2 }, (double []) { refined_x1, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&coarse_bdata[2].grid, 2, (double []) { refined_x1, refined_y2 }, (double []) { refined_x2, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&coarse_bdata[3].grid, 2, (double []) { refined_x2, refined_y2 }, (double []) { coarse_x2, coarse_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&coarse_bdata[4].grid, 2, (double []) { coarse_x1, refined_y1 }, (double []) { refined_x1, refined_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&coarse_bdata[5].grid, 2, (double []) { refined_x2, refined_y1 }, (double []) { coarse_x2, refined_y2 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&coarse_bdata[6].grid, 2, (double []) { coarse_x1, coarse_y1 }, (double []) { refined_x1, refined_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&coarse_bdata[7].grid, 2, (double []) { refined_x1, coarse_y1 }, (double []) { refined_x2, refined_y1 },
    (int []) { Nx, Ny });
  gkyl_rect_grid_init(&coarse_bdata[8].grid, 2, (double []) { refined_x2, coarse_y1 }, (double []) { coarse_x2, refined_y1 },
    (int []) { Nx, Ny });

#ifdef AMR_DEBUG
  struct euler_block_data fine_bdata[num_blocks];
  struct gkyl_job_pool *fine_job_pool = gkyl_thread_pool_new(app_args.num_threads);

  gkyl_rect_grid_init(&fine_bdata[0].grid, 2, (double []) { refined_x1, refined_y1 }, (double []) { refined_x2, refined_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
  gkyl_rect_grid_init(&fine_bdata[1].grid, 2, (double []) { coarse_x1, refined_y2 }, (double []) { refined_x1, coarse_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
  gkyl_rect_grid_init(&fine_bdata[2].grid, 2, (double []) { refined_x1, refined_y2 }, (double []) { refined_x2, coarse_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
  gkyl_rect_grid_init(&fine_bdata[3].grid, 2, (double []) { refined_x2, refined_y2 }, (double []) { coarse_x2, coarse_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
  gkyl_rect_grid_init(&fine_bdata[4].grid, 2, (double []) { coarse_x1, refined_y1 }, (double []) { refined_x1, refined_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
  gkyl_rect_grid_init(&fine_bdata[5].grid, 2, (double []) { refined_x2, refined_y1 }, (double []) { coarse_x2, refined_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
  gkyl_rect_grid_init(&fine_bdata[6].grid, 2, (double []) { coarse_x1, coarse_y1 }, (double []) { refined_x1, refined_y1 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
  gkyl_rect_grid_init(&fine_bdata[7].grid, 2, (double []) { refined_x1, coarse_y1 }, (double []) { refined_x2, refined_y1 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
  gkyl_rect_grid_init(&fine_bdata[8].grid, 2, (double []) { refined_x2, coarse_y1 }, (double []) { coarse_x2, refined_y1 },
    (int []) { Nx * ref_factor, Ny * ref_factor });
#endif

  for (int i = 0; i < num_blocks; i++) {
    coarse_bdata[i].fv_proj = gkyl_fv_proj_new(&coarse_bdata[i].grid, 2, 5, eval, 0);

#ifdef AMR_DEBUG
    fine_bdata[i].fv_proj = gkyl_fv_proj_new(&fine_bdata[i].grid, 2, 5, eval, 0);
#endif
  }
  
  for (int i = 0; i < num_blocks; i++) {
    gkyl_create_grid_ranges(&coarse_bdata[i].grid, (int []) { 2, 2 }, &coarse_bdata[i].ext_range, &coarse_bdata[i].range);
    coarse_bdata[i].geom = gkyl_wave_geom_new(&coarse_bdata[i].grid, &coarse_bdata[i].ext_range, 0, 0, false);

#ifdef AMR_DEBUG
    gkyl_create_grid_ranges(&fine_bdata[i].grid, (int []) { 2, 2 }, &fine_bdata[i].ext_range, &fine_bdata[i].range);
    fine_bdata[i].geom = gkyl_wave_geom_new(&fine_bdata[i].grid, &fine_bdata[i].ext_range, 0, 0, false);
#endif
  }

  for (int i = 0; i < num_blocks; i++) {
    if (low_order_flux) {
      struct gkyl_wv_euler_inp inp = {
        .gas_gamma = gas_gamma,
        .rp_type = WV_EULER_RP_HLL,
        .use_gpu = app_args.use_gpu,
      };
      coarse_bdata[i].euler = gkyl_wv_euler_inew(&inp);
    }
    else {
      coarse_bdata[i].euler = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    }

    for (int d = 0; d < ndim; d++) {
      coarse_bdata[i].slvr[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &coarse_bdata[i].grid,
          .equation = coarse_bdata[i].euler,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = coarse_bdata[i].geom,
        }
      );
    }

#ifdef AMR_DEBUG
    if (low_order_flux) {
      struct gkyl_wv_euler_inp inp = {
        .gas_gamma = gas_gamma,
        .rp_type = WV_EULER_RP_HLL,
        .use_gpu = app_args.use_gpu,
      };
      fine_bdata[i].euler = gkyl_wv_euler_inew(&inp);
    }
    else {
      fine_bdata[i].euler = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    }

    for (int d = 0; d < ndim; d++) {
      fine_bdata[i].slvr[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &fine_bdata[i].grid,
          .equation = fine_bdata[i].euler,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = fine_bdata[i].geom,
        }
      );
    }
#endif
  }

  struct gkyl_block_topo *btopo = create_block_topo();

  for (int i = 0; i < num_blocks; i++) {
    euler_block_bc_updaters_init(&coarse_bdata[i], &btopo->conn[i]);

#ifdef AMR_DEBUG
    euler_block_bc_updaters_init(&fine_bdata[i], &btopo->conn[i]);
#endif
  }

  for (int i = 0; i < num_blocks; i++) {
    coarse_bdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      coarse_bdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);
    }

#ifdef AMR_DEBUG
    fine_bdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      fine_bdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);
    }
#endif
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_blocks; i++) {
    gkyl_job_pool_add_work(coarse_job_pool, euler_init_job_func_block, &coarse_bdata[i]);

#ifdef AMR_DEBUG
    gkyl_job_pool_add_work(fine_job_pool, euler_init_job_func_block, &fine_bdata[i]);
#endif
  }
  gkyl_job_pool_wait(coarse_job_pool);

#ifdef AMR_DEBUG
  gkyl_job_pool_wait(fine_job_pool);
#endif
#else
  for (int i = 0; i < num_blocks; i++) {
    euler_init_job_func_block(&coarse_bdata[i]);

#ifdef AMR_DEBUG
    euler_init_job_func_block(&fine_bdata[i]);
#endif
  }
#endif

#ifdef AMR_DEBUG
  char coarse0[64];
  snprintf(coarse0, 64, "%s_coarse_0", euler_output);
  euler_write_sol_block(coarse0, num_blocks, coarse_bdata);

  char fine0[64];
  snprintf(fine0, 64, "%s_fine_0", euler_output);
  euler_write_sol_block(fine0, num_blocks, fine_bdata);

  char fine0b0[64];
  char coarse0b0[64];
  char b0[64];
  snprintf(fine0b0, 64, "%s_fine_0_b0.gkyl", euler_output);
  snprintf(coarse0b0, 64, "%s_coarse_0_b0.gkyl", euler_output);
  snprintf(b0, 64, "%s_0_b0.gkyl", euler_output);
  rename(fine0b0, b0);
  remove(coarse0b0);

  for (int i = 1; i < 9; i++) {
    char buf_old[64];
    char buf_new[64];
    char buf_del[64];

    snprintf(buf_old, 64, "%s_coarse_0_b%d.gkyl", euler_output, i);
    snprintf(buf_new, 64, "%s_0_b%d.gkyl", euler_output, i);
    snprintf(buf_del, 64, "%s_fine_0_b%d.gkyl", euler_output, i);

    rename(buf_old, buf_new);
    remove(buf_del);
  }
#else
  char amr0[64];
  snprintf(amr0, 64, "%s_0", euler_output);
  euler_write_sol_block(amr0, num_blocks, coarse_bdata);
#endif

  double coarse_t_curr = 0.0;
  double fine_t_curr = 0.0;
  double coarse_dt = euler_max_dt_block(num_blocks, coarse_bdata);

#ifdef AMR_DEBUG
  double fine_dt = euler_max_dt_block(num_blocks, fine_bdata);
#else
  double fine_dt = (1.0 / ref_factor) * coarse_dt;
#endif

  struct sim_stats stats = { };

  struct timespec tm_start = gkyl_wall_clock();

  long coarse_step = 1;
  long num_steps = app_args.num_steps;

  double io_trigger = t_end / num_frames;

  while ((coarse_t_curr < t_end) && (coarse_step <= num_steps)) {
    printf("Taking coarse (level 0) time-step %ld at t = %g; ", coarse_step, coarse_t_curr);
    struct gkyl_update_status coarse_status = euler_update_block(coarse_job_pool, btopo, coarse_bdata, coarse_t_curr, coarse_dt, &stats);
    printf(" dt = %g\n", coarse_status.dt_actual);

    if (!coarse_status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }

    for (long fine_step = 1; fine_step < ref_factor + 1; fine_step++) {
#ifdef AMR_DEBUG
      printf("   Taking fine (level 1) time-step %ld at t = %g; ", fine_step, fine_t_curr);
      struct gkyl_update_status fine_status = euler_update_block(fine_job_pool, btopo, fine_bdata, fine_t_curr, fine_dt, &stats);
      printf(" dt = %g\n", fine_status.dt_actual);

      if (!fine_status.success) {
        printf("** Update method failed! Aborting simulation ....\n");
        break;
      }

      fine_t_curr += fine_status.dt_actual;
      fine_dt = fine_status.dt_suggested;
#else
      printf("   Taking fine (level 1) time-step %ld at t = %g", fine_step, fine_t_curr);
      printf(" dt = %g\n", (1.0 / ref_factor) * coarse_status.dt_actual);

      fine_t_curr += (1.0 / ref_factor) * coarse_status.dt_actual;
      fine_dt = (1.0 / ref_factor) * coarse_status.dt_suggested;
#endif
    }

    for (int i = 1; i < num_frames; i++) {
      if (coarse_t_curr < (i * io_trigger) && (coarse_t_curr + coarse_status.dt_actual) > (i * io_trigger)) {
#ifdef AMR_DEBUG
      char buf_coarse[64];
      char buf_fine[64];

      snprintf(buf_coarse, 64, "%s_coarse_%d", euler_output, i);
      snprintf(buf_fine, 64, "%s_fine_%d", euler_output, i);

      euler_write_sol_block(buf_coarse, num_blocks, coarse_bdata);
      euler_write_sol_block(buf_fine, num_blocks, fine_bdata);

      char buf_fine_old[64];
      char buf_fine_new[64];
      char buf_coarse_old[64];

      snprintf(buf_fine_old, 64, "%s_fine_%d_b0.gkyl", euler_output, i);
      snprintf(buf_fine_new, 64, "%s_%d_b0.gkyl", euler_output, i);
      snprintf(buf_coarse_old, 64, "%s_coarse_%d_b0.gkyl", euler_output, i);

      rename(buf_fine_old, buf_fine_new);
      remove(buf_coarse_old);

      for (int j = 1; j < 9; j++) {
        char buf_old[64];
        char buf_new[64];
        char buf_del[64];

        snprintf(buf_old, 64, "%s_coarse_%d_b%d.gkyl", euler_output, i, j);
        snprintf(buf_new, 64, "%s_%d_b%d.gkyl", euler_output, i, j);
        snprintf(buf_del, 64, "%s_fine_%d_b%d.gkyl", euler_output, i, j);

        rename(buf_old, buf_new);
        remove(buf_del);
      }
#else
      char buf[64];
      snprintf(buf, 64, "%s_%d", euler_output, i);

      euler_write_sol_block(buf, num_blocks, coarse_bdata);
#endif
      }
    }

    coarse_t_curr += coarse_status.dt_actual;
    coarse_dt = coarse_status.dt_suggested;

    coarse_step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

#ifdef AMR_DEBUG
  char buf_coarse[64];
  char buf_fine[64];

  snprintf(buf_coarse, 64, "%s_coarse_%d", euler_output, num_frames);
  snprintf(buf_fine, 64, "%s_fine_%d", euler_output, num_frames);

  euler_write_sol_block(buf_coarse, num_blocks, coarse_bdata);
  euler_write_sol_block(buf_fine, num_blocks, fine_bdata);

  char buf_fine_old[64];
  char buf_fine_new[64];
  char buf_coarse_old[64];

  snprintf(buf_fine_old, 64, "%s_fine_%d_b0.gkyl", euler_output, num_frames);
  snprintf(buf_fine_new, 64, "%s_%d_b0.gkyl", euler_output, num_frames);
  snprintf(buf_coarse_old, 64, "%s_coarse_%d_b0.gkyl", euler_output,num_frames);

  rename(buf_fine_old, buf_fine_new);
  remove(buf_coarse_old);

  for (int i = 1; i < 9; i++) {
    char buf_old[64];
    char buf_new[64];
    char buf_del[64];

    snprintf(buf_old, 64, "%s_coarse_%d_b%d.gkyl", euler_output, num_frames, i);
    snprintf(buf_new, 64, "%s_%d_b%d.gkyl", euler_output, num_frames, i);
    snprintf(buf_del, 64, "%s_fine_%d_b%d.gkyl", euler_output, num_frames, i);

    rename(buf_old, buf_new);
    remove(buf_del);
  }
#else
  char buf[64];
  snprintf(buf, 64, "%s_%d", euler_output, num_frames);

  euler_write_sol_block(buf, num_blocks, coarse_bdata);
#endif

  printf("Total run-time: %g. Failed steps: %d\n", tm_total_sec, stats.nfail);

  for (int i = 0; i < num_blocks; i++) {
    gkyl_fv_proj_release(coarse_bdata[i].fv_proj);
    gkyl_wv_eqn_release(coarse_bdata[i].euler);
    euler_block_bc_updaters_release(&coarse_bdata[i]);
    gkyl_wave_geom_release(coarse_bdata[i].geom);

#ifdef AMR_DEBUG
    gkyl_fv_proj_release(fine_bdata[i].fv_proj);
    gkyl_wv_eqn_release(fine_bdata[i].euler);
    euler_block_bc_updaters_release(&fine_bdata[i]);
    gkyl_wave_geom_release(fine_bdata[i].geom);
#endif

    for (int d = 0; d < ndim; d++) {
      gkyl_wave_prop_release(coarse_bdata[i].slvr[d]);

#ifdef AMR_DEBUG
      gkyl_wave_prop_release(fine_bdata[i].slvr[d]);
#endif
    }
    
    gkyl_array_release(coarse_bdata[i].fdup);
#ifdef AMR_DEBUG
    gkyl_array_release(fine_bdata[i].fdup);
#endif
    
    for (int d = 0; d < ndim; d++) {
      gkyl_array_release(coarse_bdata[i].f[d]);

#ifdef AMR_DEBUG
      gkyl_array_release(fine_bdata[i].f[d]);
#endif
    }
  }

  gkyl_block_topo_release(btopo);
  gkyl_job_pool_release(coarse_job_pool);
#ifdef AMR_DEBUG
  gkyl_job_pool_release(fine_job_pool);
#endif
}