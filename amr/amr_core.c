#include <gkyl_amr_core.h>
#include <gkyl_amr_block_coupled_priv.h>

void
euler2d_run_single(int argc, char **argv, struct euler2d_single_init* init)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int base_Nx = init -> base_Nx;
  int base_Ny = init -> base_Ny;
  int ref_factor = init -> ref_factor;

  double coarse_x1 = init -> coarse_x1;
  double coarse_y1 = init -> coarse_y1;
  double coarse_x2 = init -> coarse_x2;
  double coarse_y2 = init -> coarse_y2;

  double refined_x1 = init -> refined_x1;
  double refined_y1 = init -> refined_y1;
  double refined_x2 = init -> refined_x2;
  double refined_y2 = init -> refined_y2;

  evalf_t eval = init -> eval;
  double gas_gamma = init -> gas_gamma;
  
  int num_frames = init -> num_frames;
  double cfl_frac = init -> cfl_frac;
  double t_end = init -> t_end;

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
    //coarse_bdata[i].euler = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    struct gkyl_wv_euler_inp inp = {
      .gas_gamma = gas_gamma,
      .rp_type = WV_EULER_RP_HLLC,
      .use_gpu = app_args.use_gpu,
    };
    coarse_bdata[i].euler = gkyl_wv_euler_inew(&inp);

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
    fine_bdata[i].euler = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);

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
    euler_block_bc_updaters_init(&coarse_bdata[i], &btopo -> conn[i]);

#ifdef AMR_DEBUG
    euler_block_bc_updaters_init(&fine_bdata[i], &btopo -> conn[i]);
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
    gkyl_job_pool_add_work(coarse_job_pool, euler_init_job_func, &coarse_bdata[i]);

#ifdef AMR_DEBUG
    gkyl_job_pool_add_work(fine_job_pool, euler_init_job_func, &fine_bdata[i]);
#endif
  }
  gkyl_job_pool_wait(coarse_job_pool);

#ifdef AMR_DEBUG
  gkyl_job_pool_wait(fine_job_pool);
#endif
#else
  for (int i = 0; i < num_blocks; i++) {
    euler_init_job_func(&coarse_bdata[i]);

#ifdef AMR_DEBUG
    euler_init_job_func(&fine_bdata[i]);
#endif
  }
#endif

#ifdef AMR_DEBUG
  euler_write_sol("euler_amr_coarse_0", num_blocks, coarse_bdata);
  euler_write_sol("euler_amr_fine_0", num_blocks, fine_bdata);

  rename("euler_amr_fine_0_b0.gkyl", "euler_amr_0_b0.gkyl");
  remove("euler_amr_coarse_0_b0.gkyl");

  for (int i = 1; i < 9; i++) {
    char buf_old[32];
    char buf_new[32];
    char buf_del[32];

    snprintf(buf_old, 32, "euler_amr_coarse_0_b%d.gkyl", i);
    snprintf(buf_new, 32, "euler_amr_0_b%d.gkyl", i);
    snprintf(buf_del, 32, "euler_amr_fine_0_b%d.gkyl", i);

    rename(buf_old, buf_new);
    remove(buf_del);
  }
#else
  euler_write_sol("euler_amr_0", num_blocks, coarse_bdata);
#endif

  double coarse_t_curr = 0.0;
  double fine_t_curr = 0.0;
  double coarse_dt = euler_max_dt(num_blocks, coarse_bdata);

#ifdef AMR_DEBUG
  double fine_dt = euler_max_dt(num_blocks, fine_bdata);
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
    struct gkyl_update_status coarse_status = euler_update(coarse_job_pool, btopo, coarse_bdata, coarse_t_curr, coarse_dt, &stats);
    printf(" dt = %g\n", coarse_status.dt_actual);

    if (!coarse_status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }

    for (long fine_step = 1; fine_step < ref_factor + 1; fine_step++) {
#ifdef AMR_DEBUG
      printf("   Taking fine (level 1) time-step %ld at t = %g; ", fine_step, fine_t_curr);
      struct gkyl_update_status fine_status = euler_update(fine_job_pool, btopo, fine_bdata, fine_t_curr, fine_dt, &stats);
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
      char buf_coarse[32];
      char buf_fine[32];

      snprintf(buf_coarse, 32, "euler_amr_coarse_%d", i);
      snprintf(buf_fine, 32, "euler_amr_fine_%d", i);

      euler_write_sol(buf_coarse, num_blocks, coarse_bdata);
      euler_write_sol(buf_fine, num_blocks, fine_bdata);

      char buf_fine_old[32];
      char buf_fine_new[32];
      char buf_coarse_old[32];

      snprintf(buf_fine_old, 32, "euler_amr_fine_%d_b0.gkyl", i);
      snprintf(buf_fine_new, 32, "euler_amr_%d_b0.gkyl", i);
      snprintf(buf_coarse_old, 32, "euler_amr_coarse_%d_b0.gkyl", i);

      rename(buf_fine_old, buf_fine_new);
      remove(buf_coarse_old);

      for (int j = 1; j < 9; j++) {
        char buf_old[32];
        char buf_new[32];
        char buf_del[32];

        snprintf(buf_old, 32, "euler_amr_coarse_%d_b%d.gkyl", i, j);
        snprintf(buf_new, 32, "euler_amr_%d_b%d.gkyl", i, j);
        snprintf(buf_del, 32, "euler_amr_fine_%d_b%d.gkyl", i, j);

        rename(buf_old, buf_new);
        remove(buf_del);
      }
#else
      char buf[32];
      snprintf(buf, 32, "euler_amr_%d", i);

      euler_write_sol(buf, num_blocks, coarse_bdata);
#endif
      }
    }

    coarse_t_curr += coarse_status.dt_actual;
    coarse_dt = coarse_status.dt_suggested;

    coarse_step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

#ifdef AMR_DEBUG
  char buf_coarse[32];
  char buf_fine[32];

  snprintf(buf_coarse, 32, "euler_amr_coarse_%d", num_frames);
  snprintf(buf_fine, 32, "euler_amr_fine_%d", num_frames);

  euler_write_sol(buf_coarse, num_blocks, coarse_bdata);
  euler_write_sol(buf_fine, num_blocks, fine_bdata);

  char buf_fine_old[32];
  char buf_fine_new[32];
  char buf_coarse_old[32];

  snprintf(buf_fine_old, 32, "euler_amr_fine_%d_b0.gkyl", num_frames);
  snprintf(buf_fine_new, 32, "euler_amr_%d_b0.gkyl", num_frames);
  snprintf(buf_coarse_old, 32, "euler_amr_coarse_%d_b0.gkyl", num_frames);

  rename(buf_fine_old, buf_fine_new);
  remove(buf_coarse_old);

  for (int i = 1; i < 9; i++) {
    char buf_old[32];
    char buf_new[32];
    char buf_del[32];

    snprintf(buf_old, 32, "euler_amr_coarse_%d_b%d.gkyl", num_frames, i);
    snprintf(buf_new, 32, "euler_amr_%d_b%d.gkyl", num_frames, i);
    snprintf(buf_del, 32, "euler_amr_fine_%d_b%d.gkyl", num_frames, i);

    rename(buf_old, buf_new);
    remove(buf_del);
  }
#else
  char buf[32];
  snprintf(buf, 32, "euler_amr_%d", num_frames);

  euler_write_sol(buf, num_blocks, coarse_bdata);
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

void
five_moment_2d_run_single(int argc, char **argv, struct five_moment_2d_single_init* init)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int base_Nx = init -> base_Nx;
  int base_Ny = init -> base_Ny;
  int ref_factor = init -> ref_factor;

  double coarse_x1 = init -> coarse_x1;
  double coarse_y1 = init -> coarse_y1;
  double coarse_x2 = init -> coarse_x2;
  double coarse_y2 = init -> coarse_y2;

  double refined_x1 = init -> refined_x1;
  double refined_y1 = init -> refined_y1;
  double refined_x2 = init -> refined_x2;
  double refined_y2 = init -> refined_y2;

  evalf_t eval_elc = init -> eval_elc;
  evalf_t eval_ion = init -> eval_ion;
  evalf_t eval_field = init -> eval_field;

  double gas_gamma = init -> gas_gamma;
  double k0_elc = init -> k0_elc;
  double k0_ion = init -> k0_ion;

  double light_speed = init -> light_speed;
  double e_fact = init -> e_fact;
  double b_fact = init -> b_fact;

  double epsilon0 = init -> epsilon0;
  double mass_elc = init -> mass_elc;
  double charge_elc = init -> charge_elc;
  double mass_ion = init -> mass_ion;
  double charge_ion = init -> charge_ion;

  bool periodic_x = init -> periodic_x;
  bool periodic_y = init -> periodic_y;

  bool wall_x = init -> wall_x;
  bool wall_y = init -> wall_y;

  double cfl_frac = init -> cfl_frac;
  double t_end = init -> t_end;
  int num_frames = init -> num_frames;

  int ndim = 2;
  int num_blocks = 9;
  int Nx = base_Nx;
  int Ny = base_Ny;

  struct five_moment_block_data coarse_bdata[num_blocks];
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
  struct five_moment_block_data fine_bdata[num_blocks];
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
    coarse_bdata[i].fv_proj_elc = gkyl_fv_proj_new(&coarse_bdata[i].grid, 2, 5, eval_elc, 0);
    coarse_bdata[i].fv_proj_ion = gkyl_fv_proj_new(&coarse_bdata[i].grid, 2, 5, eval_ion, 0);
    coarse_bdata[i].fv_proj_maxwell = gkyl_fv_proj_new(&coarse_bdata[i].grid, 2, 8, eval_field, 0);

#ifdef AMR_DEBUG
    fine_bdata[i].fv_proj_elc = gkyl_fv_proj_new(&fine_bdata[i].grid, 2, 5, eval_elc, 0);
    fine_bdata[i].fv_proj_ion = gkyl_fv_proj_new(&fine_bdata[i].grid, 2, 5, eval_ion, 0);
    fine_bdata[i].fv_proj_maxwell = gkyl_fv_proj_new(&fine_bdata[i].grid, 2, 8, eval_field, 0);
#endif
  }

  for (int i = 0; i < num_blocks; i++) {
    gkyl_create_grid_ranges(&coarse_bdata[i].grid, (int []) { 2, 2 }, &coarse_bdata[i].ext_range, &coarse_bdata[i].range);
    coarse_bdata[i].geom = gkyl_wave_geom_new(&coarse_bdata[i].grid, &coarse_bdata[i].ext_range, 0, 0, false);
    
    coarse_bdata[i].periodic_x = periodic_x;
    coarse_bdata[i].periodic_y = periodic_y;

    coarse_bdata[i].wall_x = wall_x;
    coarse_bdata[i].wall_y = wall_y;

#ifdef AMR_DEBUG
    gkyl_create_grid_ranges(&fine_bdata[i].grid, (int []) { 2, 2 }, &fine_bdata[i].ext_range, &fine_bdata[i].range);
    fine_bdata[i].geom = gkyl_wave_geom_new(&fine_bdata[i].grid, &fine_bdata[i].ext_range, 0, 0, false);

    fine_bdata[i].periodic_x = periodic_x;
    fine_bdata[i].periodic_y = periodic_y;

    fine_bdata[i].wall_x = wall_x;
    fine_bdata[i].wall_y = wall_y;
#endif
  }

  for (int i = 0; i < num_blocks; i++) {
    coarse_bdata[i].euler_elc = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    coarse_bdata[i].euler_ion = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    coarse_bdata[i].maxwell = gkyl_wv_maxwell_new(light_speed, e_fact, b_fact);

    for (int d = 0; d < ndim; d++) {
      coarse_bdata[i].slvr_elc[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &coarse_bdata[i].grid,
          .equation = coarse_bdata[i].euler_elc,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = coarse_bdata[i].geom,
        }
      );
      coarse_bdata[i].slvr_ion[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &coarse_bdata[i].grid,
          .equation = coarse_bdata[i].euler_ion,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = coarse_bdata[i].geom,
        }
      );
      coarse_bdata[i].slvr_maxwell[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &coarse_bdata[i].grid,
          .equation = coarse_bdata[i].maxwell,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = coarse_bdata[i].geom,
        }
      );
    }

    struct gkyl_moment_em_coupling_inp coarse_src_inp = {
      .grid = &coarse_bdata[i].grid,
      .nfluids = 2,
      .epsilon0 = epsilon0,
    };

    coarse_src_inp.param[0] = (struct gkyl_moment_em_coupling_data) {
      .type = coarse_bdata[i].euler_elc -> type,
      .charge = charge_elc,
      .mass = mass_elc,
      .k0 = k0_elc,
    };
    coarse_src_inp.param[1] = (struct gkyl_moment_em_coupling_data) {
      .type = coarse_bdata[i].euler_ion -> type,
      .charge = charge_ion,
      .mass = mass_ion,
      .k0 = k0_ion,
    };

    coarse_bdata[i].src_slvr = gkyl_moment_em_coupling_new(coarse_src_inp);

#ifdef AMR_DEBUG
    fine_bdata[i].euler_elc = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    fine_bdata[i].euler_ion = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    fine_bdata[i].maxwell = gkyl_wv_maxwell_new(light_speed, e_fact, b_fact);

    for (int d = 0; d < ndim; d++) {
      fine_bdata[i].slvr_elc[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &fine_bdata[i].grid,
          .equation = fine_bdata[i].euler_elc,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = fine_bdata[i].geom,
        }
      );
      fine_bdata[i].slvr_ion[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &fine_bdata[i].grid,
          .equation = fine_bdata[i].euler_ion,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = fine_bdata[i].geom,
        }
      );
      fine_bdata[i].slvr_maxwell[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &fine_bdata[i].grid,
          .equation = fine_bdata[i].maxwell,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = fine_bdata[i].geom, 
        }
      );
    }

    struct gkyl_moment_em_coupling_inp fine_src_inp = {
      .grid = &fine_bdata[i].grid,
      .nfluids = 2,
      .epsilon0 = epsilon0,
    };

    fine_src_inp.param[0] = (struct gkyl_moment_em_coupling_data) {
      .type = fine_bdata[i].euler_elc -> type,
      .charge = charge_elc,
      .mass = mass_elc,
      .k0 = k0_elc,
    };
    fine_src_inp.param[1] = (struct gkyl_moment_em_coupling_data) {
      .type = fine_bdata[i].euler_ion -> type,
      .charge = charge_ion,
      .mass = mass_ion,
      .k0 = k0_ion,
    };

    fine_bdata[i].src_slvr = gkyl_moment_em_coupling_new(fine_src_inp);
#endif
  }

  struct gkyl_block_topo *btopo = create_block_topo();

  for (int i = 0; i < num_blocks; i++) {
    five_moment_block_bc_updaters_init(&coarse_bdata[i], &btopo -> conn[i]);

#ifdef AMR_DEBUG
    five_moment_block_bc_updaters_init(&fine_bdata[i], &btopo -> conn[i]);
#endif
  }

  for (int i = 0; i < num_blocks; i++) {
    coarse_bdata[i].fdup_elc = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);
    coarse_bdata[i].fdup_ion = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);
    coarse_bdata[i].fdup_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, coarse_bdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      coarse_bdata[i].f_elc[d] = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);
      coarse_bdata[i].f_ion[d] = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);
      coarse_bdata[i].f_maxwell[d] = gkyl_array_new(GKYL_DOUBLE, 8, coarse_bdata[i].ext_range.volume);
    }

    coarse_bdata[i].app_accel_elc = gkyl_array_new(GKYL_DOUBLE, 3, coarse_bdata[i].ext_range.volume);
    coarse_bdata[i].app_accel_ion = gkyl_array_new(GKYL_DOUBLE, 3, coarse_bdata[i].ext_range.volume);
    coarse_bdata[i].rhs_source_elc = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);
    coarse_bdata[i].rhs_source_ion = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);
    coarse_bdata[i].ext_em = gkyl_array_new(GKYL_DOUBLE, 6, coarse_bdata[i].ext_range.volume);
    coarse_bdata[i].app_current = gkyl_array_new(GKYL_DOUBLE, 3, coarse_bdata[i].ext_range.volume);
    coarse_bdata[i].nT_source_elc = gkyl_array_new(GKYL_DOUBLE, 2, coarse_bdata[i].ext_range.volume);
    coarse_bdata[i].nT_source_ion = gkyl_array_new(GKYL_DOUBLE, 2, coarse_bdata[i].ext_range.volume);

#ifdef AMR_DEBUG
    fine_bdata[i].fdup_elc = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);
    fine_bdata[i].fdup_ion = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);
    fine_bdata[i].fdup_maxwell = gkyl_array_new(GKYL_DOUBLE, 8, fine_bdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++) {
      fine_bdata[i].f_elc[d] = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);
      fine_bdata[i].f_ion[d] = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);
      fine_bdata[i].f_maxwell[d] = gkyl_array_new(GKYL_DOUBLE, 8, fine_bdata[i].ext_range.volume);
    }

    fine_bdata[i].app_accel_elc = gkyl_array_new(GKYL_DOUBLE, 3, fine_bdata[i].ext_range.volume);
    fine_bdata[i].app_accel_ion = gkyl_array_new(GKYL_DOUBLE, 3, fine_bdata[i].ext_range.volume);
    fine_bdata[i].rhs_source_elc = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);
    fine_bdata[i].rhs_source_ion = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);
    fine_bdata[i].ext_em = gkyl_array_new(GKYL_DOUBLE, 6, fine_bdata[i].ext_range.volume);
    fine_bdata[i].app_current = gkyl_array_new(GKYL_DOUBLE, 3, fine_bdata[i].ext_range.volume);
    fine_bdata[i].nT_source_elc = gkyl_array_new(GKYL_DOUBLE, 2, fine_bdata[i].ext_range.volume);
    fine_bdata[i].nT_source_ion = gkyl_array_new(GKYL_DOUBLE, 2, fine_bdata[i].ext_range.volume);
#endif
  }

#ifdef AMR_USETHREADS
  for (int i = 0; i < num_blocks; i++) {
    gkyl_job_pool_add_work(coarse_job_pool, five_moment_init_job_func, &coarse_bdata[i]);

#ifdef AMR_DEBUG
    gkyl_job_pool_add_work(fine_job_pool, five_moment_init_job_func, &fine_bdata[i]);
#endif
  }
  gkyl_job_pool_wait(coarse_job_pool);

#ifdef AMR_DEBUG
  gkyl_job_pool_wait(fine_job_pool);
#endif
#else
  for (int i = 0; i < num_blocks; i++) {
    five_moment_init_job_func(&coarse_bdata[i]);

#ifdef AMR_DEBUG
    five_moment_init_job_func(&fine_bdata[i]);
#endif
  }
#endif

#ifdef AMR_DEBUG
  five_moment_write_sol("5m_amr_coarse_0", num_blocks, coarse_bdata);
  five_moment_write_sol("5m_amr_fine_0", num_blocks, fine_bdata);

  rename("5m_amr_fine_0_elc_b0.gkyl", "5m_amr_0_elc_b0.gkyl");
  rename("5m_amr_fine_0_ion_b0.gkyl", "5m_amr_0_ion_b0.gkyl");
  rename("5m_amr_fine_0_field_b0.gkyl", "5m_amr_0_field_b0.gkyl");

  remove("5m_amr_coarse_0_elc_b0.gkyl");
  remove("5m_amr_coarse_0_ion_b0.gkyl");
  remove("5m_amr_coarse_0_field_b0.gkyl");

  for (int i = 1; i < 9; i++) {
    char buf_old_elc[32];
    char buf_old_ion[32];
    char buf_old_field[32];

    char buf_new_elc[32];
    char buf_new_ion[32];
    char buf_new_field[32];

    char buf_del_elc[32];
    char buf_del_ion[32];
    char buf_del_field[32];

    snprintf(buf_old_elc, 32, "5m_amr_coarse_0_elc_b%d.gkyl", i);
    snprintf(buf_old_ion, 32, "5m_amr_coarse_0_ion_b%d.gkyl", i);
    snprintf(buf_old_field, 32, "5m_amr_coarse_0_field_b%d.gkyl", i);

    snprintf(buf_new_elc, 32, "5m_amr_0_elc_b%d.gkyl", i);
    snprintf(buf_new_ion, 32, "5m_amr_0_ion_b%d.gkyl", i);
    snprintf(buf_new_field, 32, "5m_amr_0_field_b%d.gkyl", i);
    
    snprintf(buf_del_elc, 32, "5m_amr_fine_0_elc_b%d.gkyl", i);
    snprintf(buf_del_ion, 32, "5m_amr_fine_0_ion_b%d.gkyl", i);
    snprintf(buf_del_field, 32, "5m_amr_fine_0_field_b%d.gkyl", i);

    rename(buf_old_elc, buf_new_elc);
    rename(buf_old_ion, buf_new_ion);
    rename(buf_old_field, buf_new_field);

    remove(buf_del_elc);
    remove(buf_del_ion);
    remove(buf_del_field);
  }
#else
  five_moment_write_sol("5m_amr_0", num_blocks, coarse_bdata);
#endif

  double coarse_t_curr = 0.0;
  double fine_t_curr = 0.0;
  double coarse_dt = five_moment_max_dt(num_blocks, coarse_bdata);

#ifdef AMR_DEBUG
  double fine_dt = five_moment_max_dt(num_blocks, fine_bdata);
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
    struct gkyl_update_status coarse_status = five_moment_update(coarse_job_pool, btopo, coarse_bdata, coarse_t_curr, coarse_dt, &stats);
    printf(" dt = %g\n", coarse_status.dt_actual);

    if (!coarse_status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }

    for (long fine_step = 1; fine_step < ref_factor + 1; fine_step++) {
#ifdef AMR_DEBUG
      printf("   Taking fine (level 1) time-step %ld at t = %g; ", fine_step, fine_t_curr);
      struct gkyl_update_status fine_status = five_moment_update(fine_job_pool, btopo, fine_bdata, fine_t_curr, fine_dt, &stats);
      printf(" dt = %g\n", fine_status.dt_actual);

      if (!fine_status.success) {
        printf("** Update method failed! Aborting simulation ....\n");
        break;
      }

      fine_t_curr += fine_status.dt_actual;
      fine_dt += fine_status.dt_suggested;
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
      char buf_coarse[32];
      char buf_fine[32];

      snprintf(buf_coarse, 32, "5m_amr_coarse_%d", i);
      snprintf(buf_fine, 32, "5m_amr_fine_%d", i);

      five_moment_write_sol(buf_coarse, num_blocks, coarse_bdata);
      five_moment_write_sol(buf_fine, num_blocks, fine_bdata);

      char buf_fine_old_elc[32];
      char buf_fine_old_ion[32];
      char buf_fine_old_field[32];

      char buf_fine_new_elc[32];
      char buf_fine_new_ion[32];
      char buf_fine_new_field[32];

      char buf_coarse_old_elc[32];
      char buf_coarse_old_ion[32];
      char buf_coarse_old_field[32];

      snprintf(buf_fine_old_elc, 32, "5m_amr_fine_%d_elc_b0.gkyl", i);
      snprintf(buf_fine_old_ion, 32, "5m_amr_fine_%d_ion_b0.gkyl", i);
      snprintf(buf_fine_old_field, 32, "5m_amr_fine_%d_field_b0.gkyl", i);

      snprintf(buf_fine_new_elc, 32, "5m_amr_%d_elc_b0.gkyl", i);
      snprintf(buf_fine_new_ion, 32, "5m_amr_%d_ion_b0.gkyl", i);
      snprintf(buf_fine_new_field, 32, "5m_amr_%d_field_b0.gkyl", i);

      snprintf(buf_coarse_old_elc, 32, "5m_amr_coarse_%d_elc_b0.gkyl", i);
      snprintf(buf_coarse_old_ion, 32, "5m_amr_coarse_%d_ion_b0.gkyl", i);
      snprintf(buf_coarse_old_field, 32, "5m_amr_coarse_%d_field_b0.gkyl", i);

      rename(buf_fine_old_elc, buf_fine_new_elc);
      rename(buf_fine_old_ion, buf_fine_new_ion);
      rename(buf_fine_old_field, buf_fine_new_field);

      remove(buf_coarse_old_elc);
      remove(buf_coarse_old_ion);
      remove(buf_coarse_old_field);

      for (int j = 1; j < 9; j++) {
        char buf_old_elc[32];
        char buf_old_ion[32];
        char buf_old_field[32];

        char buf_new_elc[32];
        char buf_new_ion[32];
        char buf_new_field[32];

        char buf_del_elc[32];
        char buf_del_ion[32];
        char buf_del_field[32];

        snprintf(buf_old_elc, 32, "5m_amr_coarse_%d_elc_b%d.gkyl", i, j);
        snprintf(buf_old_ion, 32, "5m_amr_coarse_%d_ion_b%d.gkyl", i, j);
        snprintf(buf_old_field, 32, "5m_amr_coarse_%d_field_b%d.gkyl", i, j);

        snprintf(buf_new_elc, 32, "5m_amr_%d_elc_b%d.gkyl", i, j);
        snprintf(buf_new_ion, 32, "5m_amr_%d_ion_b%d.gkyl", i, j);
        snprintf(buf_new_field, 32, "5m_amr_%d_field_b%d.gkyl", i, j);

        snprintf(buf_del_elc, 32, "5m_amr_fine_%d_elc_b%d.gkyl", i, j);
        snprintf(buf_del_ion, 32, "5m_amr_fine_%d_ion_b%d.gkyl", i, j);
        snprintf(buf_del_field, 32, "5m_amr_fine_%d_field_b%d.gkyl", i, j);

        rename(buf_old_elc, buf_new_elc);
        rename(buf_old_ion, buf_new_ion);
        rename(buf_old_field, buf_new_field);

        remove(buf_del_elc);
        remove(buf_del_ion);
        remove(buf_del_field);
      }
#else
      char buf[32];
      snprintf(buf, 32, "5m_amr_%d", i);

      five_moment_write_sol(buf, num_blocks, coarse_bdata);
#endif
      }
    }

    coarse_t_curr += coarse_status.dt_actual;
    coarse_dt = coarse_status.dt_suggested;

    coarse_step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

#ifdef AMR_DEBUG
  char buf_coarse[32];
  char buf_fine[32];

  snprintf(buf_coarse, 32, "5m_amr_coarse_%d", num_frames);
  snprintf(buf_fine, 32, "5m_amr_fine_%d", num_frames);

  five_moment_write_sol(buf_coarse, num_blocks, coarse_bdata);
  five_moment_write_sol(buf_fine, num_blocks, fine_bdata);

  char buf_fine_old_elc[32];
  char buf_fine_old_ion[32];
  char buf_fine_old_field[32];

  char buf_fine_new_elc[32];
  char buf_fine_new_ion[32];
  char buf_fine_new_field[32];

  char buf_coarse_old_elc[32];
  char buf_coarse_old_ion[32];
  char buf_coarse_old_field[32];

  snprintf(buf_fine_old_elc, 32, "5m_amr_fine_%d_elc_b0.gkyl", num_frames);
  snprintf(buf_fine_old_ion, 32, "5m_amr_fine_%d_ion_b0.gkyl", num_frames);
  snprintf(buf_fine_old_field, 32, "5m_amr_fine_%d_field_b0.gkyl", num_frames);

  snprintf(buf_fine_new_elc, 32, "5m_amr_%d_elc_b0.gkyl", num_frames);
  snprintf(buf_fine_new_ion, 32, "5m_amr_%d_ion_b0.gkyl", num_frames);
  snprintf(buf_fine_new_field, 32, "5m_amr_%d_field_b0.gkyl", num_frames);

  snprintf(buf_coarse_old_elc, 32, "5m_amr_coarse_%d_elc_b0.gkyl", num_frames);
  snprintf(buf_coarse_old_ion, 32, "5m_amr_coarse_%d_ion_b0.gkyl", num_frames);
  snprintf(buf_coarse_old_field, 32, "5m_amr_coarse_%d_field_b0.gkyl", num_frames);

  rename(buf_fine_old_elc, buf_fine_new_elc);
  rename(buf_fine_old_ion, buf_fine_new_ion);
  rename(buf_fine_old_field, buf_fine_new_field);

  remove(buf_coarse_old_elc);
  remove(buf_coarse_old_ion);
  remove(buf_coarse_old_field);

  for (int i = 1; i < 9; i++) {
    char buf_old_elc[32];
    char buf_old_ion[32];
    char buf_old_field[32];

    char buf_new_elc[32];
    char buf_new_ion[32];
    char buf_new_field[32];

    char buf_del_elc[32];
    char buf_del_ion[32];
    char buf_del_field[32];

    snprintf(buf_old_elc, 32, "5m_amr_coarse_%d_elc_b%d.gkyl", num_frames, i);
    snprintf(buf_old_ion, 32, "5m_amr_coarse_%d_ion_b%d.gkyl", num_frames, i);
    snprintf(buf_old_field, 32, "5m_amr_coarse_%d_field_b%d.gkyl", num_frames, i);

    snprintf(buf_new_elc, 32, "5m_amr_%d_elc_b%d.gkyl", num_frames, i);
    snprintf(buf_new_ion, 32, "5m_amr_%d_ion_b%d.gkyl", num_frames, i);
    snprintf(buf_new_field, 32, "5m_amr_%d_field_b%d.gkyl", num_frames, i);

    snprintf(buf_del_elc, 32, "5m_amr_fine_%d_elc_b%d.gkyl", num_frames, i);
    snprintf(buf_del_ion, 32, "5m_amr_fine_%d_ion_b%d.gkyl", num_frames, i);
    snprintf(buf_del_field, 32, "5m_amr_fine_%d_field_b%d.gkyl", num_frames, i);

    rename(buf_old_elc, buf_new_elc);
    rename(buf_old_ion, buf_new_ion);
    rename(buf_old_field, buf_new_field);

    remove(buf_del_elc);
    remove(buf_del_ion);
    remove(buf_del_field);
  }
#else
  char buf[32];
  snprintf(buf, 32, "5m_amr_%d", num_frames);

  five_moment_write_sol(buf, num_blocks, coarse_bdata);
#endif

  printf("Total run-time: %g. Failed steps: %d\n", tm_total_sec, stats.nfail);

  for (int i = 0; i < num_blocks; i++) {
    gkyl_fv_proj_release(coarse_bdata[i].fv_proj_elc);
    gkyl_fv_proj_release(coarse_bdata[i].fv_proj_ion);
    gkyl_fv_proj_release(coarse_bdata[i].fv_proj_maxwell);

    gkyl_wv_eqn_release(coarse_bdata[i].euler_elc);
    gkyl_wv_eqn_release(coarse_bdata[i].euler_ion);
    gkyl_wv_eqn_release(coarse_bdata[i].maxwell);

    five_moment_block_bc_updaters_release(&coarse_bdata[i]);
    gkyl_wave_geom_release(coarse_bdata[i].geom);

#ifdef AMR_DEBUG
    gkyl_fv_proj_release(fine_bdata[i].fv_proj_elc);
    gkyl_fv_proj_release(fine_bdata[i].fv_proj_ion);
    gkyl_fv_proj_release(fine_bdata[i].fv_proj_maxwell);

    gkyl_wv_eqn_release(fine_bdata[i].euler_elc);
    gkyl_wv_eqn_release(fine_bdata[i].euler_ion);
    gkyl_wv_eqn_release(fine_bdata[i].maxwell);

    five_moment_block_bc_updaters_release(&fine_bdata[i]);
    gkyl_wave_geom_release(fine_bdata[i].geom);
#endif

    for (int d = 0; d < ndim; d++) {
      gkyl_wave_prop_release(coarse_bdata[i].slvr_elc[d]);
      gkyl_wave_prop_release(coarse_bdata[i].slvr_ion[d]);
      gkyl_wave_prop_release(coarse_bdata[i].slvr_maxwell[d]);

#ifdef AMR_DEBUG
      gkyl_wave_prop_release(fine_bdata[i].slvr_elc[d]);
      gkyl_wave_prop_release(fine_bdata[i].slvr_ion[d]);
      gkyl_wave_prop_release(fine_bdata[i].slvr_maxwell[d]);
#endif
    }

    gkyl_array_release(coarse_bdata[i].fdup_elc);
    gkyl_array_release(coarse_bdata[i].fdup_ion);
    gkyl_array_release(coarse_bdata[i].fdup_maxwell);

#ifdef AMR_DEBUG
    gkyl_array_release(fine_bdata[i].fdup_elc);
    gkyl_array_release(fine_bdata[i].fdup_ion);
    gkyl_array_release(fine_bdata[i].fdup_maxwell);
#endif

    for(int d = 0; d < ndim; d++) {
      gkyl_array_release(coarse_bdata[i].f_elc[d]);
      gkyl_array_release(coarse_bdata[i].f_ion[d]);
      gkyl_array_release(coarse_bdata[i].f_maxwell[d]);

#ifdef AMR_DEBUG
      gkyl_array_release(fine_bdata[i].f_elc[d]);
      gkyl_array_release(fine_bdata[i].f_ion[d]);
      gkyl_array_release(fine_bdata[i].f_maxwell[d]);
#endif
    }

    gkyl_array_release(coarse_bdata[i].app_accel_elc);
    gkyl_array_release(coarse_bdata[i].app_accel_ion);
    gkyl_array_release(coarse_bdata[i].rhs_source_elc);
    gkyl_array_release(coarse_bdata[i].rhs_source_ion);
    gkyl_array_release(coarse_bdata[i].ext_em);
    gkyl_array_release(coarse_bdata[i].app_current);
    gkyl_array_release(coarse_bdata[i].nT_source_elc);
    gkyl_array_release(coarse_bdata[i].nT_source_ion);

#ifdef AMR_DEBUG
    gkyl_array_release(fine_bdata[i].app_accel_elc);
    gkyl_array_release(fine_bdata[i].app_accel_ion);
    gkyl_array_release(fine_bdata[i].rhs_source_elc);
    gkyl_array_release(fine_bdata[i].rhs_source_ion);
    gkyl_array_release(fine_bdata[i].ext_em);
    gkyl_array_release(fine_bdata[i].app_current);
    gkyl_array_release(fine_bdata[i].nT_source_elc);
    gkyl_array_release(fine_bdata[i].nT_source_ion);
#endif
  }

  gkyl_block_topo_release(btopo);
  gkyl_job_pool_release(coarse_job_pool);
#ifdef AMR_DEBUG
  gkyl_job_pool_release(fine_job_pool);
#endif
}