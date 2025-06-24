#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_null_pool.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_thread_pool.h>
#include <gkyl_util.h>

void evalFunc(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = x*x + sin(3*y)*cos(z) + z*z;
}

struct thread_data {
  struct gkyl_range range; // thread-specific range
    
  gkyl_proj_on_basis *proj; // shared updater
  struct gkyl_array *f; // shared field
};

void
thread_worker(void *ctx)
{
  struct thread_data *td = ctx;
  gkyl_proj_on_basis_advance(td->proj, 0.0, &td->range, td->f);
}

int
parse_args(int argc, char **argv)
{
  int c, nthread = 1;
  while ((c = getopt(argc, argv, "+hn:")) != -1)
    switch (c)
    {
      case 'h':
        nthread = 0;
        break;

      case 'n':
        nthread = atoi(optarg);
        break;

      case '?':
        return 0;
    }
  return nthread;
}

void
proj_with_job_pool(const struct gkyl_job_pool *job_pool, const char *pname)
{
  int poly_order = 2;
  double lower[] = {-2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {100, 100, 100};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc, 0);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  int max_thread = job_pool->pool_size;
  
  struct thread_data td[max_thread];
  for (int tid=0; tid<max_thread; ++tid)
    td[tid]  = (struct thread_data) {
      .range = gkyl_range_split(&arr_range, max_thread, tid),
      .proj = projDistf,
      .f = distf
    };

  struct timespec tstart = gkyl_wall_clock();

  // run projection updater on threads
  for (int tid=0; tid<max_thread; ++tid)
    gkyl_job_pool_add_work(job_pool, thread_worker, &td[tid]);
  gkyl_job_pool_wait(job_pool);

  double tm = gkyl_time_sec(gkyl_time_diff(tstart, gkyl_wall_clock()));
  printf("%d %s took %g to update\n", job_pool->pool_size, pname, tm);

  // construct file name and write data out
  const char *fmt = "%s-%d.gkyl";
  int sz = snprintf(0, 0, fmt, pname, max_thread);
  char fileNm[sz+1]; // ensures no buffer overflow  
  snprintf(fileNm, sizeof fileNm, fmt, pname, max_thread);
  
  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, fileNm);
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

int
main(int argc, char **argv)
{
  int max_thread = parse_args(argc, argv);
  if (max_thread < 1) {
    printf("Usage: rt_job_pool_proj -n <num-threads>\n");
    exit(1);
  }

  // run with thread-pool
  struct gkyl_job_pool *thread_pool = gkyl_thread_pool_new(max_thread);
  proj_with_job_pool(thread_pool, "thread_pool");

  // run with null pool
  struct gkyl_job_pool *null_pool = gkyl_null_pool_new(max_thread);
  proj_with_job_pool(null_pool, "null_pool");
  
  gkyl_job_pool_release(thread_pool);
  gkyl_job_pool_release(null_pool);
  
  return 0;
}
