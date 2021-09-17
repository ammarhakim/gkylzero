#include <acutest.h>

#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>

#include <math.h>

static void
nomapc2p(double t, const double *xc, double *xp, void *ctx)
{
  int *ndim = ctx;
  for (int i=0; i<(*ndim); ++i) xp[i] = xc[i];
}

void
test_wv_geom_1d_1()
{
  int ndim = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {10};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, ndim, cells);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, nomapc2p, &ndim);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, 1.0, 1e-15) );
    TEST_CHECK( cg->lenr[0] == 1.0 );

    TEST_CHECK( cg->norm[0][0] == 1.0 );
    TEST_CHECK( cg->tau1[0][1] == 1.0 );
    TEST_CHECK( cg->tau2[0][2] == 1.0 );
  }

  gkyl_wave_geom_release(wg);
}

static void
mapc2p(double t, const double *xc, double *xp, void *ctx)
{
  // quadratic mapping
  int *ndim = ctx;
  for (int i=0; i<(*ndim); ++i) xp[i] = xc[i]*xc[i];
}

void
test_wv_geom_1d_2()
{
  int ndim = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, ndim, cells);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, mapc2p, &ndim);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  // cell 1
  do {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, (int[]) { range.lower[0]+0 } );
    TEST_CHECK( gkyl_compare_double( cg->kappa, 0.5*0.5/0.5, 1e-15) );
    TEST_CHECK( cg->lenr[0] == 1.0 );
    TEST_CHECK( cg->norm[0][0] == 1.0 );
    TEST_CHECK( cg->tau1[0][1] == 1.0 );
    TEST_CHECK( cg->tau2[0][2] == 1.0 );
  } while (0);

  // cell 2
  do {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, (int[]) { range.lower[0]+1 } );
    TEST_CHECK( gkyl_compare_double( cg->kappa, (1-0.5*0.5)/0.5, 1e-15) );
    TEST_CHECK( cg->lenr[0] == 1.0 );
    TEST_CHECK( cg->norm[0][0] == 1.0 );
    TEST_CHECK( cg->tau1[0][1] == 1.0 );
    TEST_CHECK( cg->tau2[0][2] == 1.0 );
  } while (0);  

  gkyl_wave_geom_release(wg);
}

void
test_wv_geom_2d_1()
{
  int ndim = 2;
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int cells[] = {2, 2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, ndim, cells);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, nomapc2p, &ndim);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[1], 1.0, 1e-15) );

    // normal to left face is ex
    TEST_CHECK( gkyl_compare_double( cg->norm[0][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm[0][2], 0.0, 1e-15) );

    // tangent1 to left face is ey
    TEST_CHECK( gkyl_compare_double( cg->tau1[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1[0][1], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1[0][2], 0.0, 1e-15) );

    // tangent2 to left face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][2], 1.0, 1e-15) );

    // normal to bottom face is ey
    TEST_CHECK( gkyl_compare_double( cg->norm[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm[1][1], -1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm[1][0], 0.0, 1e-15) );

    // tangent1 to bottom face is ex
    TEST_CHECK( gkyl_compare_double( cg->tau1[1][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1[1][2], 0.0, 1e-15) );

    // tangent2 to bottom face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][2], 1.0, 1e-15) );    
  }

  gkyl_wave_geom_release(wg);
}

static void
mapc2p_2d(double t, const double *xc, double *xp, void *ctx)
{
  double x = xc[0], y = xc[1];

  xp[0] = 0.375*((x+1.0)*y+x+1.0)-0.25*((x+1.0)*y-1.0*x-1.0)-0.125*((x-1.0)*y+x-1.0);
  xp[1] = 0.25*((x+1.0)*y+x+1.0)-0.25*((x-1.0)*y+x-1.0);
}

void
test_wv_geom_2d_2()
{
  int ndim = 2;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {1, 1};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, ndim, cells);

  struct gkyl_wave_geom *wg = gkyl_wave_geom_new(&grid, &range, mapc2p_2d, &ndim);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  
  while (gkyl_range_iter_next(&iter)) {
    const struct gkyl_wave_cell_geom *cg = gkyl_wave_geom_get(wg, iter.idx);

    TEST_CHECK( gkyl_compare_double( cg->kappa, 1.0/4.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[0], 1.118033988749895/2.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->lenr[1], 1.0/2.0, 1e-15) );

    // normal to left face is ex
    TEST_CHECK( gkyl_compare_double( cg->norm[0][0], 0.8944271909999159, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm[0][1], -0.4472135954999579, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm[0][2], 0.0, 1e-15) );

    // tangent1 to left face
    TEST_CHECK( gkyl_compare_double( cg->tau1[0][0], 0.4472135954999579, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1[0][1], 0.8944271909999159, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1[0][2], 0.0, 1e-15) );

    // tangent2 to left face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][2], 1.0, 1e-15) );

    // normal to bottom face is ey
    TEST_CHECK( gkyl_compare_double( cg->norm[1][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm[1][1], -1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->norm[1][0], 0.0, 1e-15) );

    // tangent1 to bottom face is ex
    TEST_CHECK( gkyl_compare_double( cg->tau1[1][0], 1.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1[1][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau1[1][2], 0.0, 1e-15) );

    // tangent2 to bottom face is ez
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][0], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][1], 0.0, 1e-15) );
    TEST_CHECK( gkyl_compare_double( cg->tau2[0][2], 1.0, 1e-15) );    
  }

  gkyl_wave_geom_release(wg);
}

TEST_LIST = {
  { "wv_geom_1d_1", test_wv_geom_1d_1 },
  { "wv_geom_1d_2", test_wv_geom_1d_2 },
  { "wv_geom_2d_1", test_wv_geom_2d_1 },
  { "wv_geom_2d_2", test_wv_geom_2d_2 },
  { NULL, NULL },
};
