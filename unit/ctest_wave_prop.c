#include <math.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wave_prop.h>
#include <acutest.h>

#ifdef GKYL_HAVE_CUDA

// XXX following is for dev testing only, not for real unit test

static struct gkyl_array*
mkarr1(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static void
my_mapc2p(double t, const double *xc, double *xp, void *ctx)
{
  for (int i=0; i<1; ++i) xp[i] = xc[i] * 1.1;
}

void
do_test_wave_prop_maxwel_1d(bool use_gpu)
{
  int nthreads = 16;
  int nblocks = 16;
  int ncells_per_block = nthreads;
  int ncells = nblocks * ncells_per_block;

  int ndim = 1;
  int cells[] = {ncells};
  int ghost[] = {2};
  double lower[] = {1.};
  double upper[] = {2.};

  struct gkyl_rect_grid grid;
  struct gkyl_range range, range_ext;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, &range_ext, &range);

  if (use_gpu) {
    range.nthreads = nthreads;
    range.nblocks = nblocks;
  }

  // create geometry, equation, and solver
  struct gkyl_wave_geom *geom;
  struct gkyl_wv_eqn *eqn;
  struct gkyl_wave_prop *wv;

  double c = 1.0, e_fact = 0.5, b_fact = 0.5;
  if (use_gpu) {
    geom = gkyl_wave_geom_cu_dev_new(&grid, &range_ext, my_mapc2p, &ndim);
    eqn = gkyl_wv_maxwell_cu_dev_new(c, e_fact, b_fact);
    wv = gkyl_wave_prop_cu_dev_new( (struct gkyl_wave_prop_inp) {
        .grid = &grid,
        .equation = eqn,
        .limiter = GKYL_MONOTONIZED_CENTERED,
        .num_up_dirs = 1,
        .update_dirs = { 0 },
        .cfl = 1.0,
        .geom = geom,
        }
      );
  } else {
    geom = gkyl_wave_geom_new(&grid, &range_ext, my_mapc2p, &ndim);
    eqn = gkyl_wv_maxwell_new(c, e_fact, b_fact);
    wv = gkyl_wave_prop_new( (struct gkyl_wave_prop_inp) {
        .grid = &grid,
        .equation = eqn,
        .limiter = GKYL_MONOTONIZED_CENTERED,
        .num_up_dirs = 1,
        .update_dirs = { 0 },
        .cfl = 1.0,
        .geom = geom,
        }
      );
  }

  // allocate and set data
  struct gkyl_array *qin_h = mkarr1(false, 8, range_ext.volume);
  struct gkyl_array *qou_h = mkarr1(false, 8, range_ext.volume);
  struct gkyl_array *qin_d, *qou_d;
  if (use_gpu) {
    qin_d = mkarr1(true, 8, range_ext.volume);
    qou_d = mkarr1(true, 8, range_ext.volume);
  } else {
    qin_d = qin_h;
    qou_d = qou_h;
  }

  // apply initial condition
  double *ptr = qin_h->data;
  for(int i=0; i< 8 * (ncells+4); i++) {
    ptr[i] = 1. + sin(i + 1.);
  }
  gkyl_array_copy(qin_d, qin_h);

  // advance solution
  double t = 0.0;
  double dt = 1e-4;
  double nsteps = 100;
  double step = 0;

  struct gkyl_wave_prop_status status;
  if (use_gpu) {
    while(step < nsteps) {
      status = gkyl_wave_prop_cu_dev_advance(wv, t, dt, &range, qin_d, qou_d);
      t += dt;
      step += 1;
    }
  } else {
    while(step < nsteps) {
      status = gkyl_wave_prop_advance(wv, t, dt, &range, qin_h, qou_h);
      t += dt;
      step += 1;
    }
  }

  printf("status OK? %d, dt %g\n", status.success, status.dt_suggested);

  // release data
  gkyl_wave_prop_release(wv);
  gkyl_wv_eqn_release(eqn);
  gkyl_wave_geom_release(geom);
  gkyl_array_release(qin_h);
  gkyl_array_release(qou_h);
  gkyl_array_release(qin_d);
  gkyl_array_release(qou_d);
}

void
test_wave_prop_maxwel_1d_cpu() {
  do_test_wave_prop_maxwel_1d(false); 
}

void
test_wave_prop_maxwel_1d_gpu() {
  do_test_wave_prop_maxwel_1d(true); 
  checkCuda(cudaGetLastError());
}
#endif

// map (r,theta) -> (x,y)
void
mapc2p_polar(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  double r = xc[0], th = xc[1];
  xp[0] = r*cos(th); xp[1] = r*sin(th);
}

void
test_wave_prop_2d()
{
  /* create grid and wave_geom */
  int ndim = 2;
  double lower[] = {0.25, 0}, upper[] = {1.25, M_PI / 2.};
  int cells[] = {10, 10};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  int nghost[GKYL_MAX_DIM] = { 2, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);
  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(
      &grid, &ext_range, mapc2p_polar, &ndim);

  /* create gkyl_wv_eqn */
  double c = 299792458.0, e_fact = 2.0, b_fact = 2.5;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_new(c, e_fact, b_fact);

  /* create gkyl_wave_prop */
  gkyl_wave_prop *slvr[ndim];
  for (int d=0; d<ndim; ++d) {
    slvr[d] = gkyl_wave_prop_new( (struct gkyl_wave_prop_inp) {
        .grid = &grid,
        .equation = eqn,
        .limiter = GKYL_MONOTONIZED_CENTERED,
        .num_up_dirs = 1,
        .update_dirs = { d },
        .cfl = 1.0,
        .geom = wg,
      }
    );
  }

  gkyl_wv_eqn_release(eqn);
  gkyl_wave_geom_release(wg);
  for (int d=0; d<ndim; ++d) {
    gkyl_wave_prop_release(slvr[d]);
  }
}

#ifdef GKYL_HAVE_CUDA
int cu_wave_prop_test(gkyl_wave_prop **slvr, const int ndim);

void
test_wave_prop_2d_cu()
{
  /* create grid and wave_geom */
  // in order to conveniently check validity of the wave_geom object on device,
  // use only one cell and no ghost cell
  int ndim = 2;
  double lower[] = {0.25, 0}, upper[] = {1.25, M_PI / 2.};
  int cells[] = {1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  int nghost[GKYL_MAX_DIM] = { 0, 0 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);
  struct gkyl_wave_geom *wg = gkyl_wave_geom_cu_dev_new(
      &grid, &range, mapc2p_polar, &ndim);

  /* create gkyl_wv_eqn */
  double c = 299792458.0, e_fact = 2.0, b_fact = 2.5;
  struct gkyl_wv_eqn *eqn = gkyl_wv_maxwell_cu_dev_new(c, e_fact, b_fact);

  /* create gkyl_wave_prop */
  gkyl_wave_prop *slvr[ndim];
  for (int d=0; d<ndim; ++d) {
    slvr[d] = gkyl_wave_prop_cu_dev_new( (struct gkyl_wave_prop_inp) {
        .grid = &grid,
        .equation = eqn,
        .limiter = GKYL_MONOTONIZED_CENTERED,
        .num_up_dirs = 1,
        .update_dirs = { d },
        .cfl = 1.0,
        .geom = wg,
      }
    );
  }

  cu_wave_prop_test(slvr, ndim);

  gkyl_wv_eqn_release(eqn);
  gkyl_wave_geom_release(wg);
  for (int d=0; d<ndim; ++d) {
    gkyl_wave_prop_release(slvr[d]);
  }
}
#endif

TEST_LIST = {
  { "wave_prop_2d", test_wave_prop_2d },
#ifdef GKYL_HAVE_CUDA
  { "wave_prop_maxwell_1d_cpu", test_wave_prop_maxwel_1d_cpu },
  { "wave_prop_maxwell_1d_gpu", test_wave_prop_maxwel_1d_gpu },
  { "wave_prop_2d_cu", test_wave_prop_2d_cu },
#endif
  { NULL, NULL },
};
