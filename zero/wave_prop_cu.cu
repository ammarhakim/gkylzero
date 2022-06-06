extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wave_prop_priv.h>
}

GKYL_CU_D static void limit_waves_cu() {}

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
    struct gkyl_wave_prop_status *status)
{
  return;

  int ndim = update_range->ndim;
  // int meqn = wv->equation->num_equations, mwaves = wv->equation->num_waves;
  const int meqn = 8;
  const int mwaves = 6;

  double cfla = 0.0, cfl = wv->cfl, cflm = 1.1*cfl;

  double ql_local[meqn], qr_local[meqn];
  double waves_local[meqn*mwaves];
  double delta[meqn], amdq[meqn], apdq[meqn];  

  int idxl[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  for (int d=0; d<wv->num_up_dirs; ++d) {
    int dir = wv->update_dirs[d];

    double dtdx = dt/wv->grid.dx[dir];

    // upper/lower bounds in direction 'd'. Note these are edge indices
    // int loidx = update_range->lower[dir]-1;
    // int upidx = update_range->upper[dir]+2;

    struct gkyl_range slice_range = *update_range;
    struct gkyl_range_iter iter;
    // gkyl_range_init(&slice_range, 1, (int[]) { loidx }, (int[]) { upidx } );
    // while (gkyl_range_iter_next(&iter)) {
    {
      gkyl_copy_int_arr(ndim, iter.idx, idxl);
      gkyl_copy_int_arr(ndim, iter.idx, idxr);

      // for (int i=loidx; i<upidx; ++i) {
      {
        int i = 0;
        idxl[dir] = i-1; idxr[dir] = i;

        // geometry in cell
        const struct gkyl_wave_cell_geom *cg =
          gkyl_wave_geom_get(wv->geom, idxr);

        long sidx = gkyl_ridx(slice_range, i);
        long lidx = gkyl_range_idx(update_range, idxl);
        long ridx = gkyl_range_idx(update_range, idxr);        

        const double *qinl = (const double*)gkyl_array_cfetch(qin, lidx);
        const double *qinr = (const double*)gkyl_array_cfetch(qin, ridx);

        wv->equation->rotate_to_local_func(
            cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinl, ql_local);
        wv->equation->rotate_to_local_func(
            cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinr, qr_local);

        calc_jump(meqn, ql_local, qr_local, delta);
        double *s = (double *)gkyl_array_fetch(wv->speeds, sidx);
        wv->equation->waves_func(
            wv->equation, delta, ql_local, qr_local, waves_local, s);

        double lenr = cg->lenr[dir];
        double *waves = (double *)gkyl_array_fetch(wv->waves, sidx);
        for (int mw=0; mw<mwaves; ++mw) {
          // rotate waves back
          wv->equation->rotate_to_global_func(
            cg->tau1[dir], cg->tau2[dir], cg->norm[dir], &waves_local[mw*meqn],
            &waves[mw*meqn]
          );

          s[mw] *= lenr; // rescale speeds
        }

        wv->equation->qfluct_func(
            wv->equation, qinl, qinr, waves, s, amdq, apdq);

        double *qoutl = (double *)gkyl_array_fetch(qout, lidx);
        double *qoutr = (double *)gkyl_array_fetch(qout, ridx);

        calc_first_order_update(meqn, dtdx/cg->kappa, qoutl, qoutr, amdq, apdq);
        cfla = calc_cfla(mwaves, cfla, dtdx/cg->kappa, s);
      }

      if (cfla > cflm) {
        status->success = 0;
        status->dt_suggested = dt*cfl/cfla;
        return;
      }

      // apply limiters to waves for all edges in update range,
      // including edges that are on the range boundary
      // limit_waves_cu(wv, &slice_range,
      //   update_range->lower[dir], update_range->upper[dir]+1, wv->waves,
      //   wv->speeds);
      limit_waves_cu();

      // get the kappa in the first ghost cell on left (needed in the
      // second order flux calculation)
      idxl[dir] = update_range->lower[dir]-1;
      const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);
      double kappal = cg->kappa;

      // FIXME gkyl_array_clear(wv->flux2, 0.0);
      // compute second-order correction fluxes at each interface:
      // note that there is one extra edge than cell
      for (int i=update_range->lower[dir]; i<=update_range->upper[dir]+1; ++i) {
        long sidx = gkyl_ridx(slice_range, i);

        const double *waves = (const double*)gkyl_array_cfetch(wv->waves, sidx);
        const double *s = (const double*)gkyl_array_cfetch(wv->speeds, sidx);
        double *flux2 = (double *)gkyl_array_fetch(wv->flux2, sidx);

        idxl[dir] = i;
        const struct gkyl_wave_cell_geom *cg =
          gkyl_wave_geom_get(wv->geom, idxl);
        double kappar = cg->kappa;

        for (int mw=0; mw<mwaves; ++mw)
          calc_second_order_flux(
              meqn, dtdx/(0.5*(kappal+kappar)), s[mw], &waves[mw*meqn], flux2);

        kappal = kappar;
      }

      // add second correction flux to solution in each interior cell
      for (int i=update_range->lower[dir]; i<=update_range->upper[dir]; ++i) {
        idxl[dir] = i;
        const struct gkyl_wave_cell_geom *cg =
          gkyl_wave_geom_get(wv->geom, idxl);

        calc_second_order_update(meqn, dtdx/cg->kappa,
          (double *)gkyl_array_fetch(
            qout, gkyl_range_idx(update_range, idxl)),
          (const double *)gkyl_array_cfetch(
            wv->flux2, gkyl_ridx(slice_range, i)),
          (const double *)gkyl_array_cfetch(
            wv->flux2, gkyl_ridx(slice_range, i+1))
        );
      }
    }
  }

  double dt_suggested = dt*cfl/fmax(cfla, DBL_MIN);

  status->success = 1;
  status->dt_suggested = dt_suggested > dt ? dt_suggested : dt;
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

