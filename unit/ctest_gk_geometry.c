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

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
      .geometry_id = GKYL_MAPC2P,
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = 0,
      .bmag_func = bmag_func, // magnetic field magnitude
      .bmag_ctx =0 ,
      .grid = grid,
      .local = range,
      .local_ext = ext_range,
      .global = range,
      .global_ext = ext_range,
      .basis = basis,
      .geo_grid = grid,
      .geo_local = range,
      .geo_local_ext = ext_range,
      .geo_global = range,
      .geo_global_ext = ext_range,
      .geo_basis = basis,
  };


  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);

  gkyl_gk_geometry_release(gk_geom);
}

TEST_LIST = {
  { "test_3x_p1", test_3x_p1},
  //{ "test_3x_p2", test_3x_p2},
  { NULL, NULL },
};
