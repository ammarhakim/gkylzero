extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_wave_prop.h>
}

// CPU interface to create and track a GPU object
gkyl_wave_prop*
gkyl_wave_prop_cu_dev_new(gkyl_wave_prop_inp winp)
{
  // STEP: CREATE HOST OBJECT
  gkyl_wave_prop *up = (gkyl_wave_prop *)gkyl_malloc(sizeof(gkyl_wave_prop));

  // STEP: SET HOST OR COMMON HOST/DEVICE DATA IN HOST OBJECT
  up->grid = *(winp.grid);
  up->ndim = up->grid.ndim;
  up->num_up_dirs = winp.num_up_dirs;
  for (int i=0; i<winp.num_up_dirs; ++i)
    up->update_dirs[i] = winp.update_dirs[i];
  up->limiter = winp.limiter == 0 ? GKYL_MONOTONIZED_CENTERED : winp.limiter;
  up->cfl = winp.cfl;
  int nghost[3] = { 2, 2, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&up->grid, nghost, &ext_range, &range);

  up->equation = winp.equation->on_dev;
  up->geom = winp.geom->on_dev;

  // force null pointers that can be handled by gkyl_array_release
  up->waves = NULL;
  up->speeds = NULL;
  up->flux2 = NULL;

  // STEP: COPY HOST OBJECT TO DEVICE OBJECT
  gkyl_wave_prop *up_dev =
    (gkyl_wave_prop*) gkyl_cu_malloc(sizeof(gkyl_wave_prop));
  gkyl_cu_memcpy(up_dev, up, sizeof(gkyl_wave_prop), GKYL_CU_MEMCPY_H2D);

  up->equation = gkyl_wv_eqn_acquire(winp.equation);
  up->geom = gkyl_wave_geom_acquire(winp.geom);

  // STEP: SET DEVICE DATA

  // STEP: KEEP POINTER TO THE DEVICE OBJECT
  up->on_dev = up_dev;

  return up;
}

__global__ void
do_gkyl_wave_prop_cu_dev_advance(
    const gkyl_wave_prop *wv,
    double tm,
    double dt,
    const struct gkyl_range *update_range,
    const struct gkyl_array *qin,
    struct gkyl_array *qout,
    struct gkyl_wave_prop_status *status_dev)
{
}

struct gkyl_wave_prop_status
gkyl_wave_prop_cu_dev_advance(
    const gkyl_wave_prop *wv,
    double tm,
    double dt,
    const struct gkyl_range *update_range,
    const struct gkyl_array *qin,
    struct gkyl_array *qout)
{
  int nthreads = update_range->nthreads;
  int nblocks = update_range->nblocks;

  struct gkyl_wave_prop_status *status_dev =
    (struct gkyl_wave_prop_status *) gkyl_cu_malloc(
       sizeof(struct gkyl_wave_prop_status));

  gkyl_array_copy(qout, qin);
  do_gkyl_wave_prop_cu_dev_advance<<<nthreads, nblocks>>>(
      wv, tm, dt, update_range, qin, qout, status_dev);
  checkCuda(cudaGetLastError());

  struct gkyl_wave_prop_status status;
  gkyl_cu_memcpy(&status, status_dev, sizeof(struct gkyl_wave_prop_status),
                 GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(status_dev);

  return status;
}

