#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <gkyl_util.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>

void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = x*x + sin(3*y)*cos(z) + z*z;
}

struct thread_data {
    struct gkyl_range range; // thread-specific range
    
    gkyl_proj_on_basis *proj; // shared updater
    struct gkyl_array *f; // shared field
};

void*
thread_worker(void *ctx)
{
  struct thread_data *td = ctx;
  gkyl_proj_on_basis_advance(td->proj, 0.0, &td->range, td->f);
  return NULL;
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

int
main(int argc, char **argv)
{
  int max_thread = parse_args(argc, argv);
  if (max_thread < 1) {
    printf("Usage: app_pthread_proj -n <num-threads>\n");
    exit(1);
  }
  
  int polyOrder = 2;
  double lower[] = {-2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {100, 100, 100};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, polyOrder);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    polyOrder+1, 1, evalFunc, 0);

  // create array range: no ghost-cells 
  struct gkyl_range arr_range;
  gkyl_range_init_from_shape(&arr_range, 3, cells);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.numBasis, arr_range.volume);

  pthread_t thread[max_thread];  
  struct thread_data td[max_thread];
  for (int tid=0; tid<max_thread; ++tid)
    td[tid]  = (struct thread_data) {
      .range = gkyl_range_split(&arr_range, max_thread, tid),
      .proj = projDistf,
      .f = distf
    };

  struct timespec tstart = gkyl_wall_clock();
  
  for (int tid=0; tid<max_thread; ++tid) {
    int rc = pthread_create(&thread[tid], 0, thread_worker, &td[tid]);
    if (rc) {
      printf("Thread creation failed with %d\n", rc);
      exit(-1);
    }    
  }
  for (int i=0; i<max_thread; ++i)
    pthread_join(thread[i], 0);

  double tm = gkyl_time_sec(gkyl_time_diff(tstart, gkyl_wall_clock()));
  printf("%d threads took %g to update\n", max_thread, tm);

  gkyl_grid_array_write(&grid, &arr_range, distf, "pthread.gkyl");
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  
  return 0;
}
