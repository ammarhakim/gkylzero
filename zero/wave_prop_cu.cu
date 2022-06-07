extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>
#include <gkyl_wave_prop.h>
#include <gkyl_wave_prop_priv.h>
}

GKYL_CU_D static void limit_waves_cu(const gkyl_wave_prop *wv,
                                     const struct gkyl_range *slice_range,
                                     int lower, int upper,
                                     struct gkyl_array *waves,
                                     const struct gkyl_array *speed) {
  int meqn = wv->equation->num_equations, mwave = wv->equation->num_waves;

  for (int mw = 0; mw < mwave; ++mw) {
    const double *wl = (const double *)gkyl_array_cfetch(
        waves, gkyl_ridx(*slice_range, lower - 1));
    const double *wr = (const double *)gkyl_array_cfetch(
        waves, gkyl_ridx(*slice_range, lower));

    double dotr = wave_dot_prod(meqn, &wl[mw * meqn], &wr[mw * meqn]);

    for (int i = lower; i <= upper; ++i) {
      double dotl = dotr;

      double *GKYL_RESTRICT wi =
          (double *)gkyl_array_fetch(waves, gkyl_ridx(*slice_range, i));
      const double *GKYL_RESTRICT wi1 =
          (double *)gkyl_array_cfetch(waves, gkyl_ridx(*slice_range, i + 1));

      double wnorm2 = wave_dot_prod(meqn, &wi[mw * meqn], &wi[mw * meqn]);
      dotr = wave_dot_prod(meqn, &wi[mw * meqn], &wi1[mw * meqn]);

      if (wnorm2 > 0) {
        const double *s = (const double *)gkyl_array_cfetch(
            speed, gkyl_ridx(*slice_range, i));
        double r = s[mw] > 0 ? dotl / wnorm2 : dotr / wnorm2;
        double theta = limiter_function(r, wv->limiter);
        wave_rescale(meqn, theta, &wi[mw * meqn]);
      }
    }
  }
}

// CPU interface to create and track a GPU object
gkyl_wave_prop *gkyl_wave_prop_cu_dev_new(gkyl_wave_prop_inp winp) {
  // STEP: CREATE HOST OBJECT
  gkyl_wave_prop *up = (gkyl_wave_prop *)gkyl_malloc(sizeof(gkyl_wave_prop));

  // STEP: SET HOST OR COMMON HOST/DEVICE DATA IN HOST OBJECT
  up->grid = *(winp.grid);
  up->ndim = up->grid.ndim;
  up->num_up_dirs = winp.num_up_dirs;
  for (int i = 0; i < winp.num_up_dirs; ++i)
    up->update_dirs[i] = winp.update_dirs[i];
  up->limiter = winp.limiter == 0 ? GKYL_MONOTONIZED_CENTERED : winp.limiter;
  up->cfl = winp.cfl;
  int nghost[3] = {2, 2, 2};
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
      (gkyl_wave_prop *)gkyl_cu_malloc(sizeof(gkyl_wave_prop));
  gkyl_cu_memcpy(up_dev, up, sizeof(gkyl_wave_prop), GKYL_CU_MEMCPY_H2D);

  up->equation = gkyl_wv_eqn_acquire(winp.equation);
  up->geom = gkyl_wave_geom_acquire(winp.geom);

  // STEP: SET DEVICE DATA

  // STEP: KEEP POINTER TO THE DEVICE OBJECT
  up->on_dev = up_dev;

  return up;
}

__global__ void do_gkyl_wave_prop_cu_dev_advance(
    const gkyl_wave_prop *wv, double tm, double dt,
    const struct gkyl_range update_range, const struct gkyl_array *qin,
    struct gkyl_array *qout, struct gkyl_wave_prop_status *status) {
  int ndim = update_range.ndim;
  int idxl[3], idxc[3], idxr[3];

  // int meqn = wv->equation->num_equations, mwave = wv->equation->num_waves;
  // FIXME
  const int meqn = 8;
  const int mwave = 6;

  double cfla = 0.0, cfl = wv->cfl, cflm = 1.1 * cfl;

  double ql_local[meqn], qr_local[meqn];
  double waves_local[meqn * mwave];
  double delta[meqn], amdq[meqn], apdq[meqn];

  // assign buffers for each thread to solve RP on four edges, and computing
  // fluxes on two edges for updating one cell
  extern __shared__ double dummy[];
  int base = 0;

  double *waves = dummy + base + (meqn * mwave * 4) * (threadIdx.x);
  base += (meqn * mwave * 4) * (blockDim.x);

  double *speeds = dummy + base + (mwave * 4) * (threadIdx.x);
  base += (mwave * 4) * (blockDim.x);

  double *flux2 = dummy + base + (meqn * 2) * (threadIdx.x);
  base += (meqn * 2) * (blockDim.x);

  // iterate over all cells, one thread updates one cell a time */
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x;
       linc1 < update_range.volume; linc1 += blockDim.x * gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&update_range, linc1, idxc);

    gkyl_copy_int_arr(ndim, idxc, idxl);
    gkyl_copy_int_arr(ndim, idxc, idxr);

    for (int d = 0; d < wv->num_up_dirs; ++d) {
      int dir = wv->update_dirs[d];

      double dtdx = dt / wv->grid.dx[dir];

      /****************************************/
      /* SOLVE RIEMANN PROBLEMS ON FOUR EDGES */
      /****************************************/
      for (int i = 0; i <= 3; ++i) {
        idxl[dir] += i - 2;
        idxr[dir] += i - 1; // left and right cells of the edge

        const struct gkyl_wave_cell_geom *cg =
            gkyl_wave_geom_get(wv->geom, idxr);

        long lidx = gkyl_range_idx(&update_range, idxl);
        long ridx = gkyl_range_idx(&update_range, idxr);

        const double *qinl = (const double *)gkyl_array_cfetch(qin, lidx);
        const double *qinr = (const double *)gkyl_array_cfetch(qin, ridx);

        wv->equation->rotate_to_local_func(cg->tau1[dir], cg->tau2[dir],
                                           cg->norm[dir], qinl, ql_local);
        wv->equation->rotate_to_local_func(cg->tau1[dir], cg->tau2[dir],
                                           cg->norm[dir], qinr, qr_local);

        calc_jump(meqn, ql_local, qr_local, delta);
        double *s = speeds + mwave * i;
        wv->equation->waves_func(wv->equation, delta, ql_local, qr_local,
                                 waves_local, s);

        double lenr = cg->lenr[dir];
        double *my_waves = waves + mwave * meqn * i;
        for (int mw = 0; mw < mwave; ++mw) {
          wv->equation->rotate_to_global_func(
              cg->tau1[dir], cg->tau2[dir], cg->norm[dir],
              &waves_local[mw * meqn], &my_waves[mw * meqn]);
          s[mw] *= lenr;
        }

        wv->equation->qfluct_func(wv->equation, qinl, qinr, my_waves, s, amdq,
                                  apdq);

        double *qoutl = (double *)gkyl_array_fetch(qout, lidx);
        double *qoutr = (double *)gkyl_array_fetch(qout, ridx);

        calc_first_order_update(meqn, dtdx / cg->kappa, qoutl, qoutr, amdq,
                                apdq);
        cfla = calc_cfla(mwave, cfla, dtdx / cg->kappa, s);
      }

      if (cfla > cflm) {
        status->success = 0;
        status->dt_suggested = dt * cfl / cfla;
        return;
      }

      /****************************/
      /* LIMIT WAVES ON TWO EDGES */
      /****************************/

      //   limit_waves_cu(wv, &slice_range, update_range.lower[dir],
      //      update_range.upper[dir]+1, waves, wv->speeds);

      /*******************************************************************/
      /* COMPUTE 2ND-ORDER FLUXES ON LEFT AND RIGHT EDGES OF TARGET CELL */
      /*******************************************************************/

      for (int i = 0; i < meqn * 2; ++i)
        flux2[i] = 0.;

      idxl[dir] = idxc[dir] - 1;
      const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wv->geom, idxl);
      double kappal = cg->kappa;

      for (int i = 1; i <= 2; ++i) {
        const double *my_waves = waves + mwave * meqn * i;
        const double *s = speeds + mwave * i;

        // we stored flux2 on two edges only, thus the shift is i-1 not i
        double *my_flux2 = flux2 + meqn * (i - 1);

        idxl[dir] = idxc[dir] + i - 1;

        const struct gkyl_wave_cell_geom *cg =
            gkyl_wave_geom_get(wv->geom, idxl);
        double kappar = cg->kappa;

        for (int mw = 0; mw < mwave; ++mw)
          calc_second_order_flux(meqn, dtdx / (0.5 * (kappal + kappar)), s[mw],
                                 &my_waves[mw * meqn], my_flux2);

        kappal = kappar;
      }

      /*********************************************/
      /* ADD 2ND-ORDER CORRECTION ONTO TARGET CELL */
      /*********************************************/
      long linc = gkyl_range_idx(&update_range, idxc);
      double *qc = (double *)gkyl_array_fetch(qout, linc);
      cg = gkyl_wave_geom_get(wv->geom, idxc);
      calc_second_order_update(meqn, dtdx / cg->kappa, qc, flux2, flux2 + meqn);
    }
  }

  double dt_suggested = dt * cfl / fmax(cfla, DBL_MIN);
  status->dt_suggested = dt_suggested > dt ? dt_suggested : dt;
  status->success = 1;
}

struct gkyl_wave_prop_status
gkyl_wave_prop_cu_dev_advance(const gkyl_wave_prop *wv, double tm, double dt,
                              const struct gkyl_range *update_range,
                              const struct gkyl_array *qin,
                              struct gkyl_array *qout) {
  int nthreads = update_range->nthreads;
  int nblocks = update_range->nblocks;

  int meqn = wv->equation->num_equations, mwave = wv->equation->num_waves;
  int shared_mem_size = 0;
  shared_mem_size += (meqn * mwave * 4) * (nthreads);
  shared_mem_size += (mwave * 4) * (nthreads);
  shared_mem_size += (meqn * 2) * (nthreads);
  shared_mem_size *= sizeof(double);

  struct gkyl_wave_prop_status *status_dev =
      (struct gkyl_wave_prop_status *)gkyl_cu_malloc(
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
