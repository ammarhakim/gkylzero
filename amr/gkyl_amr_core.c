#include <gkyl_amr_core.h>
#include <gkyl_amr_block_priv.h>

void
euler2d_run_level1(int argc, char **argv, evalf_t eval, int baseNx, int baseNy, int ref_factor, double refined_x1, double refined_y1, double refined_x2, double refined_y2,
    double coarse_x1, double coarse_y1, double coarse_x2, double coarse_y2, double cfl_frac, double t_end)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem)
  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  int ndim = 2;
  int num_blocks = 9;
  int Nx = baseNx;
  int Ny = baseNy;

  struct euler_block_data bdata[num_blocks];
  struct gkyl_job_pool *job_pool = gkyl_thread_pool_new(app_args.num_threads);

  gkyl_rect_grid_init(&bdata[0].grid, 2,
    (double []) { refined_x1, refined_y1 },
    (double []) { refined_x2, refined_y2 },
    (int []) { Nx, Ny }
  );

  gkyl_rect_grid_init(&bdata[1].grid, 2,
    (double []) { coarse_x1, refined_y2 },
    (double []) { refined_x1, coarse_y2 },
    (int []) { Nx, Ny }
  );

  gkyl_rect_grid_init(&bdata[2].grid, 2,
    (double []) { refined_x1, refined_y2 },
    (double []) { refined_x2, coarse_y2 },
    (int []) { Nx, Ny }
  );

  gkyl_rect_grid_init(&bdata[3].grid, 2,
    (double []) { refined_x2, refined_y2 },
    (double []) { coarse_x2, coarse_y2 },
    (int []) { Nx, Ny }
  );

  gkyl_rect_grid_init(&bdata[4].grid, 2,
    (double []) { coarse_x1, refined_y1 },
    (double []) { refined_x1, refined_y2 },
    (int []) { Nx, Ny }
  );

  gkyl_rect_grid_init(&bdata[5].grid, 2,
    (double []) { refined_x2, refined_y1 },
    (double []) { coarse_x2, refined_y2 },
    (int []) { Nx, Ny }
  );

  gkyl_rect_grid_init(&bdata[6].grid, 2,
    (double []) { coarse_x1, coarse_y1 },
    (double []) { refined_x1, refined_y1 },
    (int []) { Nx, Ny }
  );

  gkyl_rect_grid_init(&bdata[7].grid, 2,
    (double []) { refined_x1, coarse_y1 },
    (double []) { refined_x2, refined_y1 },
    (int []) { Nx, Ny }
  );

  gkyl_rect_grid_init(&bdata[8].grid, 2,
    (double []) { refined_x2, coarse_y1 },
    (double []) { coarse_x2, refined_y1 },
    (int []) { Nx, Ny }
  );

  for (int i = 0; i < num_blocks; i++)
  {
    bdata[i].fv_proj = gkyl_fv_proj_new(&bdata[i].grid, 2, 5, eval, 0);
  }

  for (int i = 0; i < num_blocks; i++)
  {
    gkyl_create_grid_ranges(&bdata[i].grid, (int []) { 2, 2 }, &bdata[i].ext_range, &bdata[i].range);
    bdata[i].geom = gkyl_wave_geom_new(&bdata[i].grid, &bdata[i].ext_range, 0, 0, false);
  }

  for (int i = 0; i < num_blocks; i++)
  {
    bdata[i].euler = gkyl_wv_euler_new(1.4, false);

    for (int d = 0; d < ndim; d++)
    {
      bdata[i].slvr[d] = gkyl_wave_prop_new(& (struct gkyl_wave_prop_inp) {
          .grid = &bdata[i].grid,
          .equation = bdata[i].euler,
          .limiter = GKYL_MONOTONIZED_CENTERED,
          .num_up_dirs = 1,
          .update_dirs = { d },
          .cfl = cfl_frac,
          .geom = bdata[i].geom
        }
      );
    }
  }

  struct gkyl_block_topo *btopo = create_block_topo();

  for (int i = 0; i < num_blocks; i++)
  {
    euler_block_bc_updaters_init(&bdata[i], &btopo -> conn[i]);
  }

  for (int i = 0; i < num_blocks; i++)
  {
    bdata[i].fdup = gkyl_array_new(GKYL_DOUBLE, 5, bdata[i].ext_range.volume);

    for (int d = 0; d < ndim + 1; d++)
    {
      bdata[i].f[d] = gkyl_array_new(GKYL_DOUBLE, 5, bdata[i].ext_range.volume);
    }
  }

  for (int i = 0; i < num_blocks; i++)
  {
    gkyl_job_pool_add_work(job_pool, euler_init_job_func, &bdata[i]);
  }
  gkyl_job_pool_wait(job_pool);

  euler_write_sol("euler_amr_0", num_blocks, bdata);

  double t_curr = 0.0;
  double dt = euler_max_dt(num_blocks, bdata);

  struct sim_stats stats = { };

  struct timespec tm_start = gkyl_wall_clock();

  long step = 1;
  long num_steps = app_args.num_steps;

  while ((t_curr < t_end) && (step <= num_steps))
  {
    printf("Taking time-step %ld at t = %g ...", step, t_curr);
    struct gkyl_update_status status = euler_update(job_pool, btopo, bdata, t_curr, dt, &stats);
    printf(" dt = %g\n", status.dt_actual);

    if (!status.success)
    {
      printf ("** Update method failed! Aborting simulaton ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    step += 1;
  }

  double tm_total_sec = gkyl_time_diff_now_sec(tm_start);

  euler_write_sol("euler_amr_1", num_blocks, bdata);

  printf("Total run-time: %g. Failed steps: %d\n", tm_total_sec, stats.nfail);

  for (int i = 0; i < num_blocks; i++)
  {
    gkyl_fv_proj_release(bdata[i].fv_proj);

    gkyl_wv_eqn_release(bdata[i].euler);
    euler_block_bc_updaters_release(&bdata[i]);

    gkyl_wave_geom_release(bdata[i].geom);

    for (int d = 0; d < ndim; d++)
    {
      gkyl_wave_prop_release(bdata[i].slvr[d]);
    }

    gkyl_array_release(bdata[i].fdup);

    for (int d = 0; d < ndim + 1; d++)
    {
      gkyl_array_release(bdata[i].f[d]);
    }
  }

  gkyl_block_topo_release(btopo);
  gkyl_job_pool_release(job_pool);
}