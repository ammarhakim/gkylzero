#include <gkyl_amr_core.h>
#include <gkyl_amr_block_priv.h>

void
euler2d_run_level1(int argc, char **argv, struct euler2d_level1_init* init)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem)
  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  long baseNx = init -> baseNx;
  long baseNy = init -> baseNy;
  long ref_factor = init -> ref_factor;

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

  double cfl_frac = init -> cfl_frac;
  double t_end = init -> t_end;

  int ndim = 2;
  int num_blocks = 9;
  long Nx = baseNx;
  long Ny = baseNy;

  struct euler_block_data coarse_bdata[num_blocks];
  struct euler_block_data fine_bdata[num_blocks];

  struct gkyl_job_pool *coarse_job_pool = gkyl_thread_pool_new(app_args.num_threads);
  struct gkyl_job_pool *fine_job_pool = gkyl_thread_pool_new(app_args.num_threads);

  gkyl_rect_grid_init(&coarse_bdata[0].grid, 2, (double []) { refined_x1, refined_y1 }, (double []) { refined_x2, refined_y2 }, (int []) { Nx, Ny });
  gkyl_rect_grid_init(&fine_bdata[0].grid, 2, (double []) { refined_x1, refined_y1 }, (double []) { refined_x2, refined_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  gkyl_rect_grid_init(&coarse_bdata[1].grid, 2, (double []) { coarse_x1, refined_y2 }, (double []) { refined_x1, coarse_y2 }, (int []) { Nx, Ny });
  gkyl_rect_grid_init(&fine_bdata[1].grid, 2, (double []) { coarse_x1, refined_y2 }, (double []) { refined_x1, coarse_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  gkyl_rect_grid_init(&coarse_bdata[2].grid, 2, (double []) { refined_x1, refined_y2 }, (double []) { refined_x2, coarse_y2 }, (int []) { Nx, Ny });
  gkyl_rect_grid_init(&fine_bdata[2].grid, 2, (double []) { refined_x1, refined_y2 }, (double []) { refined_x2, coarse_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  gkyl_rect_grid_init(&coarse_bdata[3].grid, 2, (double []) { refined_x2, refined_y2 }, (double []) { coarse_x2, coarse_y2 }, (int []) { Nx, Ny });
  gkyl_rect_grid_init(&fine_bdata[3].grid, 2, (double []) { refined_x2, refined_y2 }, (double []) { coarse_x2, coarse_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  gkyl_rect_grid_init(&coarse_bdata[4].grid, 2, (double []) { coarse_x1, refined_y1 }, (double []) { refined_x1, refined_y2 }, (int []) { Nx, Ny });
  gkyl_rect_grid_init(&fine_bdata[4].grid, 2, (double []) { coarse_x1, refined_y1 }, (double []) { refined_x1, refined_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  gkyl_rect_grid_init(&coarse_bdata[5].grid, 2, (double []) { refined_x2, refined_y1 }, (double []) { coarse_x2, refined_y2 }, (int []) { Nx, Ny });
  gkyl_rect_grid_init(&fine_bdata[5].grid, 2, (double []) { refined_x2, refined_y1 }, (double []) { coarse_x2, refined_y2 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  gkyl_rect_grid_init(&coarse_bdata[6].grid, 2, (double []) { coarse_x1, coarse_y1 }, (double []) { refined_x1, refined_y1 }, (int []) { Nx, Ny });
  gkyl_rect_grid_init(&fine_bdata[6].grid, 2, (double []) { coarse_x1, coarse_y1 }, (double []) { refined_x1, refined_y1 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  gkyl_rect_grid_init(&coarse_bdata[7].grid, 2, (double []) { refined_x1, coarse_y1 }, (double []) { refined_x2, refined_y1 }, (int []) { Nx, Ny });
  gkyl_rect_grid_init(&fine_bdata[7].grid, 2, (double []) { refined_x1, coarse_y1 }, (double []) { refined_x2, refined_y1 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  gkyl_rect_grid_init(&coarse_bdata[8].grid, 2, (double []) { refined_x2, coarse_y1 }, (double []) { coarse_x2, refined_y1 }, (int []) { Nx, Ny });
  gkyl_rect_grid_init(&fine_bdata[8].grid, 2, (double []) { refined_x2, coarse_y1 }, (double []) { coarse_x2, refined_y1 },
    (int []) { Nx * ref_factor, Ny * ref_factor });

  for (int i = 0; i < num_blocks; i++)
  {
    coarse_bdata[i].fv_proj = gkyl_fv_proj_new(&coarse_bdata[i].grid, 2, 5, eval, 0);
    fine_bdata[i].fv_proj = gkyl_fv_proj_new(&fine_bdata[i].grid, 2, 5, eval, 0);
  }

  for (int i = 0; i < num_blocks; i++)
  {
    gkyl_create_grid_ranges(&coarse_bdata[i].grid, (int []) { 2, 2 }, &coarse_bdata[i].ext_range, &coarse_bdata[i].range);
    gkyl_create_grid_ranges(&fine_bdata[i].grid, (int []) { 2, 2 }, &fine_bdata[i].ext_range, &fine_bdata[i].range);

    coarse_bdata[i].geom = gkyl_wave_geom_new(&coarse_bdata[i].grid, &coarse_bdata[i].ext_range, 0, 0, false);
    fine_bdata[i].geom = gkyl_wave_geom_new(&fine_bdata[i].grid, &fine_bdata[i].ext_range, 0, 0, false);
  }

  for (int i = 0; i < num_blocks; i++)
  {
    coarse_bdata[i].euler = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);
    fine_bdata[i].euler = gkyl_wv_euler_new(gas_gamma, app_args.use_gpu);

    for (int d = 0; d < ndim; d++)
    {
      coarse_bdata[i].slvr[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &coarse_bdata[i].grid,
          .equation = coarse_bdata[i].euler,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = coarse_bdata[i].geom
        }
      );

      fine_bdata[i].slvr[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &fine_bdata[i].grid,
          .equation = fine_bdata[i].euler,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = fine_bdata[i].geom
        }
      );
    }
  }

  struct gkyl_block_topo *btopo = create_block_topo();

  for (int i = 0; i < num_blocks; i++)
  {
    euler_block_bc_updaters_init(&coarse_bdata[i], &btopo -> conn[i]);
    euler_block_bc_updaters_init(&fine_bdata[i], &btopo -> conn[i]);
  }

  for (int i = 0; i < num_blocks; i++)
  {
    coarse_bdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);
    fine_bdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++)
    {
      coarse_bdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 5, coarse_bdata[i].ext_range.volume);
      fine_bdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 5, fine_bdata[i].ext_range.volume);
    }
  }

  for (int i = 0; i < num_blocks; i++)
  {
    gkyl_job_pool_add_work(coarse_job_pool, euler_init_job_func, &coarse_bdata[i]);
    gkyl_job_pool_add_work(fine_job_pool, euler_init_job_func, &fine_bdata[i]);
  }
  gkyl_job_pool_wait(coarse_job_pool);
  gkyl_job_pool_wait(fine_job_pool);

  euler_write_sol("euler_amr_coarse_0", num_blocks, coarse_bdata);
  euler_write_sol("euler_amr_fine_0", num_blocks, fine_bdata);

  double coarse_t_curr = 0.0;
  double fine_t_curr = 0.0;
  double coarse_dt = euler_max_dt(num_blocks, coarse_bdata);
  double fine_dt = euler_max_dt(num_blocks, fine_bdata);

  struct sim_stats stats = { };

  struct timespec tm_start = gkyl_wall_clock();

  long coarse_step = 1;
  long num_steps = app_args.num_steps;

  while ((coarse_t_curr < t_end) && (coarse_step <= num_steps))
  {
    printf("Taking coarse (level 0) time-step %ld at t = %g; ", coarse_step, coarse_t_curr);
    struct gkyl_update_status coarse_status = euler_update(coarse_job_pool, btopo, coarse_bdata, coarse_t_curr, coarse_dt, &stats);
    printf(" dt = %g\n", coarse_status.dt_actual);

    if (!coarse_status.success)
    {
      printf("** Update method failed! Aborting simulaton ....\n");
      break;
    }

    for (long fine_step = 1; fine_step < ref_factor + 1; fine_step++)
    {
      printf("   Taking fine (level 1) time-step %ld at t = %g; ", fine_step, fine_t_curr);
      struct gkyl_update_status fine_status = euler_update(fine_job_pool, btopo, fine_bdata, fine_t_curr, fine_dt, &stats);
      printf(" dt = %g\n", fine_status.dt_actual);
      
      if (!fine_status.success)
      {
        printf("   ** Update method failed! Aborting simulation ....\n");
        break;
      }
      
      fine_t_curr += fine_status.dt_actual;
      
      fine_dt = fine_status.dt_suggested;
    }

    coarse_t_curr += coarse_status.dt_actual;
    coarse_dt = coarse_status.dt_suggested;

    coarse_step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

  euler_write_sol("euler_amr_coarse_1", num_blocks, coarse_bdata);
  euler_write_sol("euler_amr_fine_1", num_blocks, fine_bdata);

  printf("Total run-time: %g. Failed steps: %d\n", tm_total_sec, stats.nfail);

  for (int i = 0; i < num_blocks; i++)
  {
    gkyl_fv_proj_release(coarse_bdata[i].fv_proj);
    gkyl_fv_proj_release(fine_bdata[i].fv_proj);

    gkyl_wv_eqn_release(coarse_bdata[i].euler);
    gkyl_wv_eqn_release(fine_bdata[i].euler);

    euler_block_bc_updaters_release(&coarse_bdata[i]);
    euler_block_bc_updaters_release(&fine_bdata[i]);

    gkyl_wave_geom_release(coarse_bdata[i].geom);
    gkyl_wave_geom_release(fine_bdata[i].geom);

    for (int d = 0; d < ndim; d++)
    {
      gkyl_wave_prop_release(coarse_bdata[i].slvr[d]);
      gkyl_wave_prop_release(fine_bdata[i].slvr[d]);
    }

    gkyl_array_release(coarse_bdata[i].fdup);
    gkyl_array_release(fine_bdata[i].fdup);

    for (int d = 0; d < ndim + 1; d++)
    {
      gkyl_array_release(coarse_bdata[i].f[d]);
      gkyl_array_release(fine_bdata[i].f[d]);
    }
  }

  gkyl_block_topo_release(btopo);
  gkyl_job_pool_release(coarse_job_pool);
  gkyl_job_pool_release(fine_job_pool);
}