#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include <gkyl_util.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>

void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = x*x;
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
  printf("thread_worker: %ld %ld \n", td->range.th_start, td->range.th_len);
  
  gkyl_proj_on_basis_advance(td->proj, 0.0, &td->range, td->f);
  return NULL;
}

int
main(void)
{
  int polyOrder = 1;
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {50000};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, polyOrder);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    polyOrder+1, 1, evalFunc, NULL);

  // create array range: no ghost-cells 
  struct gkyl_range arr_range;
  gkyl_range_init_from_shape(&arr_range, 1, cells);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(sizeof(double[basis.numBasis]),
    arr_range.volume);

  int max_thread = 2;
  pthread_t thread[max_thread];  
  struct thread_data td[max_thread];

  struct timespec tstart = gkyl_wall_clock();
  
  for (int tid=0; tid<max_thread; ++tid) {

    gkyl_range_thread(&arr_range, max_thread, tid);
    td[tid]  = (struct thread_data) { .range = arr_range, .proj = projDistf, .f = distf };
    
    int rc = pthread_create(&thread[tid], NULL, thread_worker, &td[tid]);
    
    if (rc) {
      printf("Thread creation failed with %d\n", rc);
      exit(-1);
    }    
  }

  for (int i=0; i<max_thread; ++i)
    pthread_join(thread[i], NULL);

  double tm = gkyl_time_sec(gkyl_time_diff(tstart, gkyl_wall_clock()));
  printf("%d threads took %g to update\n", max_thread, tm);

  gkyl_range_thread(&arr_range, 1, 0);
  gkyl_grid_array_write(&grid, &arr_range, distf, "pthread.gkyl");
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
  
  return 0;
}
