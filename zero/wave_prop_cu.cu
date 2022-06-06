extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wave_prop_priv.h>
}

GKYL_CU_D static void
limit_waves_cu(const gkyl_wave_prop *wv, const struct gkyl_range *slice_range,
  int lower, int upper, struct gkyl_array *waves, const struct gkyl_array *speed)
{
  int meqn = wv->equation->num_equations, mwave = wv->equation->num_waves;

  for (int mw=0; mw<mwave; ++mw) {
    const double *wl = (const double*)gkyl_array_cfetch(
        waves, gkyl_ridx(*slice_range, lower-1));
    const double *wr = (const double*)gkyl_array_cfetch(
        waves, gkyl_ridx(*slice_range, lower));

    double dotr = wave_dot_prod(meqn, &wl[mw*meqn], &wr[mw*meqn]);

    for (int i=lower; i<=upper; ++i) {
      double dotl = dotr;
      
      double * GKYL_RESTRICT wi =
        (double *)gkyl_array_fetch(waves, gkyl_ridx(*slice_range, i));
      const double * GKYL_RESTRICT wi1 =
        (double *)gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, i+1));
      
      double wnorm2 = wave_dot_prod(meqn, &wi[mw*meqn], &wi[mw*meqn]);
      dotr = wave_dot_prod(meqn, &wi[mw*meqn], &wi1[mw*meqn]);

      if (wnorm2 > 0) {
        const double *s = (const double *)gkyl_array_cfetch(
            speed, gkyl_ridx(*slice_range, i));
        double r = s[mw] > 0 ? dotl/wnorm2 : dotr/wnorm2;
        double theta = limiter_function(r, wv->limiter);
        wave_rescale(meqn, theta, &wi[mw*meqn]);
      }
    }
  }
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
    const struct gkyl_range update_range,
    const struct gkyl_array *qin,
    struct gkyl_array *qout,
    struct gkyl_wave_prop_status *status)
{
  int ndim = update_range.ndim;
  int idxl[3], idxc[3], idxr[3];
  double xcl[3], xcc[3], xcr[3];

  // int meqn = wv->equation->num_equations, mwave = wv->equation->num_waves;
  // FIXME
  const int meqn = 8;
  const int mwave = 6;

  double cfla = 0.0, cfl = wv->cfl, cflm = 1.1*cfl;

  double ql_local[meqn], qr_local[meqn];
  double waves_local[meqn*mwave];
  double delta[meqn], amdq[meqn], apdq[meqn];

  // assign buffer for each thread; note that each thread does RP on four edges
  // to update one cell
  extern __shared__ double dummy[];
  int base = 0;

  double *waves = dummy + base + (meqn * mwave * 4) * (threadIdx.x);
  base += (meqn * mwave * 5) * (blockDim.x);

  double *speeds = dummy + base + (mwave * 4) * (threadIdx.x);
  base += (mwave * 4) * (blockDim.x);

  double *flux2 = dummy + base + (meqn * 4) * (threadIdx.x);
  base += (meqn * 4) * (blockDim.x);

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < update_range.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&update_range, linc1, idxc);
    gkyl_rect_grid_cell_center(&wv->grid, idxc, xcc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&update_range, idxc);

    gkyl_copy_int_arr(ndim, idxc, idxl);
    gkyl_copy_int_arr(ndim, idxc, idxr);

    for (int d=0; d<wv->num_up_dirs; ++d) {
      int dir = wv->update_dirs[d];

      double dtdx = dt/wv->grid.dx[dir];

      // solve RP on four edges to update one cell
      for (int i=-1; i<3; ++i) {
        idxl[dir] = i-1; idxr[dir] = i;
        int sidx = i + 1;

        // geometry in cell
        const struct gkyl_wave_cell_geom *cg =
          gkyl_wave_geom_get(wv->geom, idxr);

        long lidx = gkyl_range_idx(&update_range, idxl);
        long ridx = gkyl_range_idx(&update_range, idxr);        

        const double *qinl = (const double*)gkyl_array_cfetch(qin, lidx);
        const double *qinr = (const double*)gkyl_array_cfetch(qin, ridx);

        wv->equation->rotate_to_local_func(
            cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinl, ql_local);
        wv->equation->rotate_to_local_func(
            cg->tau1[dir], cg->tau2[dir], cg->norm[dir], qinr, qr_local);

        calc_jump(meqn, ql_local, qr_local, delta);
        double *s = speeds + mwave * sidx;
        wv->equation->waves_func(
            wv->equation, delta, ql_local, qr_local, waves_local, s);

        double lenr = cg->lenr[dir];
        double *my_waves = waves + mwave * meqn * sidx;
        for (int mw=0; mw<mwave; ++mw) {
          // rotate waves back
          wv->equation->rotate_to_global_func(cg->tau1[dir], cg->tau2[dir],
              cg->norm[dir], &waves_local[mw*meqn], &my_waves[mw*meqn]
          );

          s[mw] *= lenr; // rescale speeds
        }

        wv->equation->qfluct_func(
            wv->equation, qinl, qinr, my_waves, s, amdq, apdq);

        double *qoutl = (double *)gkyl_array_fetch(qout, lidx);
        double *qoutr = (double *)gkyl_array_fetch(qout, ridx);

        calc_first_order_update(meqn, dtdx/cg->kappa, qoutl, qoutr, amdq, apdq);
        cfla = calc_cfla(mwave, cfla, dtdx/cg->kappa, s);
    }

      if (cfla > cflm) {
        status->success = 0;
        status->dt_suggested = dt*cfl/cfla;
        status->dt_suggested = 1.0; // for testing
        return;
      }
    //
    //   // apply limiters to waves for all edges in update range,
    //   // including edges that are on the range boundary
    //   limit_waves_cu(wv, &slice_range, update_range.lower[dir],
    //      update_range.upper[dir]+1, waves, wv->speeds);

      // get the kappa in the first ghost cell on left (needed in the
      // second order flux calculation)
      idxl[dir] = -1;
      const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);
      double kappal = cg->kappa;

      // FIXME gkyl_array_clear(wv->flux2, 0.0);
      // compute second-order correction fluxes at each interface:
      // note that there is one extra edge than cell
      for (int i=0; i<2; ++i) {
        int sidx = i + 1;

        const double *my_waves = waves + mwave * meqn * sidx;
        const double *s = speeds + mwave * sidx;
        double *flux2 = flux2 + meqn * sidx;

        idxl[dir] = i;
        const struct gkyl_wave_cell_geom *cg =
          gkyl_wave_geom_get(wv->geom, idxl);
        double kappar = cg->kappa;

        for (int mw=0; mw<mwave; ++mw)
          calc_second_order_flux(
              meqn, dtdx/(0.5*(kappal+kappar)), s[mw], &my_waves[mw*meqn], flux2);

        kappal = kappar;
      }

      // add second correction flux to solution in each interior cell
      idxl[dir] = 0;
      cg = gkyl_wave_geom_get(wv->geom, idxl);

      calc_second_order_update(meqn, dtdx/cg->kappa,
        (double *)gkyl_array_fetch( qout, gkyl_range_idx(&update_range, idxl)),
        flux2 + meqn * 0, flux2 + meqn * 1);
    }
  }

  double dt_suggested = dt*cfl/fmax(cfla, DBL_MIN);
  dt_suggested = 1.0; // for testing

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

  int meqn = wv->equation->num_equations, mwave = wv->equation->num_waves;
  int shared_mem_size = 0;
  shared_mem_size += (meqn * mwave * 4) * (nthreads);
  shared_mem_size += (mwave * 4) * (nthreads);
  shared_mem_size += (meqn * 4) * (nthreads);
  shared_mem_size *= sizeof(double);

  struct gkyl_wave_prop_status *status_dev =
    (struct gkyl_wave_prop_status *) gkyl_cu_malloc(
       sizeof(struct gkyl_wave_prop_status));

  gkyl_array_copy(qout, qin);
  do_gkyl_wave_prop_cu_dev_advance<<<nthreads, nblocks, shared_mem_size>>>(
      wv->on_dev, tm, dt, *update_range, qin->on_dev, qout->on_dev, status_dev);
  checkCuda(cudaGetLastError());

  struct gkyl_wave_prop_status status;
  gkyl_cu_memcpy(&status, status_dev, sizeof(struct gkyl_wave_prop_status),
                 GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(status_dev);

  return status;
}

