#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


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

#include <gkyl_efit.h>
#include <gkyl_calc_bmag.h>
#include <gkyl_tok_geo.h>


#include <gkyl_calc_metric.h>
#include <gkyl_calc_derived_geo.h>

#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>





// Z is constant at 8.429
// R goes from 4.9 to 5.9
void horizontal_pfunc_upper(double s, double* RZ){
  RZ[0] = 4.9 + s;
  RZ[1] = 8.429;
}

// Z is constant at -8.429
// R goes from 4.9 to 5.3
void horizontal_pfunc_lower(double s, double* RZ){
  RZ[0] = 4.9 + s;
  RZ[1] = -8.429;
}

// R is constant at 4.9
// R goes from 8.0 to 8.5
void vertical_pfunc_upper(double s, double* RZ){
  RZ[0] = 4.9;
  RZ[1] = 8.0 + s/2;
}

void vertical_pfunc_lower(double s, double* RZ){
  RZ[0] = 4.9;
  RZ[1] = -8.0 - s/2;
}

// Actual Plate info:
// p1 = [5.7413,8.5471]
// p2 = [5.8544, 8.5229]
// p3 = [5.8549, 8.4258]
// Try a Different (slanted) plate instead
// p1 [5.151,8.516]
// p2 [5.852, 8.434]
void shaped_pfunc_upper(double s, double* RZ){
    RZ[0] = 0.8 + (0.916 - 0.8)*s;
    RZ[1] = -1.2 + (-1.329 + 1.2)*s;
}

void shaped_pfunc_lower(double s, double* RZ){
    RZ[0] = 1.6 + (1.8 - 1.6)*s;
    RZ[1] = -1.26 + (-1.1 + 1.26)*s;
}




void
test_fixed_z()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };

  double clower[] = { 0.16, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.17501, 0.01, M_PI-1e-14 };

  int ccells[] = { 1, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .rmin = 0.0,
    .rmax = 5.0,
    .ftype = GKYL_SOL_SN_LO,
    .rclose = 2.5,
    .rright = 2.5,
    .rleft = 0.7,
    .zmin = -1.3,
    .zmax = 1.0,
    .zmin_left = -1.3,
    .zmin_right = -1.3,
  }; 
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}

void
test_shaped_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_efit_inp efit_inp = {
      // psiRZ and related inputs
      .filepath = "gyrokinetic/data/eqdsk/asdex.geqdsk",
      .rz_poly_order = 2,
      .flux_poly_order = 1,
    };

  double clower[] = { 0.16, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.17501, 0.01, M_PI-1e-14 };

  int ccells[] = { 1, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_SN_LO,
    .rmin = 0.0,
    .rmax = 5.0,
    .rclose = 2.5,
    .rright = 2.5,
    .rleft = 0.7,
    .zmin = -1.3,
    .zmax = 1.0,
    .zmin_left = -1.2,
    .zmin_right = -1.0,
  }; 
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .efit_info = efit_inp,
    .tok_grid_info = ginp,
    .grid = cgrid,
    .local = clocal,
    .local_ext = clocal_ext,
    .global = clocal,
    .global_ext = clocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = clocal,
    .geo_local_ext = clocal_ext,
    .geo_global = clocal,
    .geo_global_ext = clocal_ext,
    .geo_basis = cbasis,
  };

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&geometry_inp); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
}

TEST_LIST = {
  //{ "test_fixed_z", test_fixed_z},
  { "test_shaped_plate", test_shaped_plate},
  { NULL, NULL },
};
