#include <math.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_maxwell.h>
#include <gkyl_wave_prop.h>
#include <acutest.h>

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
  { "wave_prop_2d_cu", test_wave_prop_2d_cu },
#endif
  { NULL, NULL },
};
