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
    RZ[0] = 5.151 + (5.852 - 5.151)*s;
    RZ[1] = 8.516 + (8.434 - 8.516)*s;
}

void shaped_pfunc_lower(double s, double* RZ){
    RZ[0] = 5.151 + (5.852 - 5.151)*s;
    RZ[1] = -(8.516 + (8.434 - 8.516)*s);
}




void
test_fixed_z()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_tok_geo_inp inp = {
      // psiRZ and related inputs
      .filepath = "./efit_data/input.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = false,
      .quad_param = {  .eps = 1e-10 }
    };

  double clower[] = { 0.934, -0.01, -3.14 };
  double cupper[] = {1.0, 0.01, 3.14 };

  int ccells[] = { 2, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_geo_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
    .write_node_coord_array = true,
    .node_file_nm = "step_outboard_fixed_z_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_horizontal_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_tok_geo_inp inp = {
      // psiRZ and related inputs
      .filepath = "./efit_data/input.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = true,
      // can set plate func if you dont want a fixed zmin and zmax
      .plate_func_lower = horizontal_pfunc_lower,
      .plate_func_upper = horizontal_pfunc_upper,
      .quad_param = {  .eps = 1e-10 }
    };

  double clower[] = { 0.934, -0.01, -3.14 };
  double cupper[] = {1.3, 0.01, 3.14 };

  int ccells[] = { 2, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_geo_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
    .write_node_coord_array = true,
    .node_file_nm = "step_outboard_horizontal_plate_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_vertical_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_tok_geo_inp inp = {
      // psiRZ and related inputs
      .filepath = "./efit_data/input.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = true,
      // can set plate func if you dont want a fixed zmin and zmax
      .plate_func_lower = vertical_pfunc_lower,
      .plate_func_upper = vertical_pfunc_upper,
      .quad_param = {  .eps = 1e-10 }
    };

  double clower[] = { 0.934, -0.01, -3.14 };
  double cupper[] = {1.4688, 0.01, 3.14 };

  int ccells[] = { 2, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_geo_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
    .write_node_coord_array = true,
    .node_file_nm = "step_outboard_vertical_plate_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_shaped_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_tok_geo_inp inp = {
      // psiRZ and related inputs
      .filepath = "./efit_data/input.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = true,
      // can set plate func if you dont want a fixed zmin and zmax
      .plate_func_lower = shaped_pfunc_lower,
      .plate_func_upper = shaped_pfunc_upper,
      .quad_param = {  .eps = 1e-10 }
    };

  double clower[] = { 0.934, -0.01, -3.14 };
  double cupper[] = {1.4688, 0.01, 3.14 };

  int ccells[] = { 2, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_geo_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
    .write_node_coord_array = true,
    .node_file_nm = "step_outboard_shaped_plate_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

TEST_LIST = {
  //{ "test_fixed_z", test_fixed_z},
  //{ "test_horizontal_plate", test_horizontal_plate},
  //{ "test_vertical_plate", test_vertical_plate},
  { "test_shaped_plate", test_shaped_plate},
  { NULL, NULL },
};
