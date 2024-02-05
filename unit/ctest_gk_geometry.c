#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>

#include <gkyl_comm.h>

#include <gkyl_calc_bmag.h>
#include <gkyl_mirror_geo_priv.h>

// Helper Functions

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Functions for this test
void mapc2p(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], z = xn[2];
  fout[0] = r*cos(theta); fout[1] = r*sin(theta); fout[2] = z;
}

void exact_gij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  fout[0] = 1.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
  fout[3] = r*r;
  fout[4] = 0.0;
  fout[5] = 1.0;
}

void bmag_func(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx){
  //printf("callig func\n");
  fout[0] = 0.0398;
}

void
test_3x_p1()
{
  struct gkyl_basis basis;
  int poly_order = 1;
  gkyl_cart_modal_serendip(&basis, 3, poly_order);
  
  
  double Lz = 1.8049e+01;
  double Lx = 1.2534e+00;
  double Rmax = Lx/2;
  double Rmin = 0.05*Rmax;
  int Nz = 10;

  double lower[3] = {Rmin, -M_PI, -Lz/2};
  double upper[3] = {Rmax,  M_PI,  Lz/2};
  int cells[3] = { 18, 18, Nz };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  
  struct gkyl_range ext_range, range;
  int nghost[3] = { 1,1,1};
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);
  struct gk_geometry* gkgeom = gkyl_gk_geometry_mapc2p_new(&grid, &range, &ext_range, &basis, mapc2p, 0, bmag_func, 0, false);
  gkyl_gk_geometry_release(gkgeom);
}

void
test_uniform_to_nonuniform_fraczero()
{
  // Tests the function map_theta_to_z when using a zero mapping_fraction
  // This should provide a completely uniform mapping
  struct arc_length_ctx arc_app;
  arc_app.mapping_order_expander = 10;
  arc_app.mapping_order_center = 10;
  arc_app.theta_min = -1.0;
  arc_app.theta_max = 1.0;
  arc_app.theta_throat = 0.5;
  arc_app.mapping_frac = 0.0;
  double dtheta = 0.1;
  for (int i=0; i<20; i++){
    double theta = arc_app.theta_min + i*dtheta;
    double nonuniform_coordinate = map_theta_to_z(theta, &arc_app);
    TEST_CHECK( gkyl_compare_double(nonuniform_coordinate, theta, 1e-12) );
  }
}

void
test_uniform_to_nonuniform_linear()
{
  // Tests the function map_theta_to_z when using a linear polynomial
  // This should provide a completely uniform mapping
  struct arc_length_ctx arc_app;
  arc_app.mapping_order_expander = 1;
  arc_app.mapping_order_center = 1;
  arc_app.theta_min = -1.0;
  arc_app.theta_max = 1.0;
  arc_app.theta_throat = 0.5;
  arc_app.mapping_frac = 1.0;
  double dtheta = 0.1;
  for (int i=0; i<20; i++){
    double theta = arc_app.theta_min + i*dtheta;
    double nonuniform_coordinate = map_theta_to_z(theta, &arc_app);
    TEST_CHECK( gkyl_compare_double(nonuniform_coordinate, theta, 1e-12) );
  }
}

void
test_uniform_to_nonuniform_fixed_points()
{
  // Tests the function map_theta_to_z such that
  // the fixed points at theta_min, \pm z_m, 0, and theta_max are preserved
  struct arc_length_ctx arc_app;
  arc_app.mapping_order_expander = 10;
  arc_app.mapping_order_center = 10;
  arc_app.theta_min = -1.0;
  arc_app.theta_max = 1.0;
  arc_app.theta_throat = 0.5;
  arc_app.mapping_frac = 1.0;
  
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(arc_app.theta_min, &arc_app), arc_app.theta_min, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-arc_app.theta_throat, &arc_app), -arc_app.theta_throat, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(0.0, &arc_app), 0.0, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(arc_app.theta_throat, &arc_app), arc_app.theta_throat, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(arc_app.theta_max, &arc_app), arc_app.theta_max, 1e-12) );
}

void
test_uniform_to_nonuniform_arbitrary_points()
{
  // Tests the function map_theta_to_z such that
  // the fixed points at theta_min, \pm z_m, 0, and theta_max are preserved
  struct arc_length_ctx arc_app;
  arc_app.mapping_order_expander = 2;
  arc_app.mapping_order_center = 2;
  arc_app.theta_min = -1.0;
  arc_app.theta_max = 1.0;
  arc_app.theta_throat = 0.5;
  arc_app.mapping_frac = 0.5;
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-1.0, &arc_app), -1.0,  1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-0.9, &arc_app), -0.86, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-0.8, &arc_app), -0.74, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-0.7, &arc_app), -0.64, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-0.6, &arc_app), -0.56, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-0.5, &arc_app), -0.50, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-0.4, &arc_app), -0.44, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-0.3, &arc_app), -0.36, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-0.2, &arc_app), -0.26, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z(-0.1, &arc_app), -0.14, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.0, &arc_app),   0.0, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.1, &arc_app),  0.14, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.2, &arc_app),  0.26, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.3, &arc_app),  0.36, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.4, &arc_app),  0.44, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.5, &arc_app),  0.50, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.6, &arc_app),  0.56, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.7, &arc_app),  0.64, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.8, &arc_app),  0.74, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 0.9, &arc_app),  0.86, 1e-12) );
  TEST_CHECK( gkyl_compare_double(map_theta_to_z( 1.0, &arc_app),  1.0,  1e-12) );
}


TEST_LIST = {
  { "test_3x_p1", test_3x_p1},
  //{ "test_3x_p2", test_3x_p2},
  {"test_uniform_to_nonuniform_fraczero", test_uniform_to_nonuniform_fraczero},
  {"test_uniform_to_nonuniform_linear", test_uniform_to_nonuniform_linear},
  {"test_uniform_to_nonuniform_fixed_points", test_uniform_to_nonuniform_fixed_points},
  {"test_uniform_to_nonuniform_arbitrary_points", test_uniform_to_nonuniform_arbitrary_points},
  { NULL, NULL },
};
