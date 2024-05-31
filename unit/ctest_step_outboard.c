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

void
write_geometry(gk_geometry *up, struct gkyl_rect_grid grid, struct gkyl_range local, const char *name)
{
  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, name, "jacobtot_inv");
  char fileNm[sz+1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "mapc2p");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->mc2p, fileNm);
  sprintf(fileNm, fmt, name, "bmag");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->bmag, fileNm);
  sprintf(fileNm, fmt, name, "g_ij");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->g_ij, fileNm);
  sprintf(fileNm, fmt, name, "dxdz");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->dxdz, fileNm);
  sprintf(fileNm, fmt, name, "dzdx");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->dzdx, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobgeo, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobgeo_inv, fileNm);
  sprintf(fileNm, fmt, name, "gij");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gij, fileNm);
  sprintf(fileNm, fmt, name, "b_i");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->b_i, fileNm);
  sprintf(fileNm, fmt, name, "cmag");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->cmag, fileNm);
  sprintf(fileNm, fmt, name, "jacobtot");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobtot, fileNm);
  sprintf(fileNm, fmt, name, "jacobtot_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobtot_inv, fileNm);
  sprintf(fileNm, fmt, name, "bmag_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->bmag_inv, fileNm);
  sprintf(fileNm, fmt, name, "bmag_inv_sq");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->bmag_inv_sq, fileNm);
  sprintf(fileNm, fmt, name, "gxxj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gxxj, fileNm);
  sprintf(fileNm, fmt, name, "gxyj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gxyj, fileNm);
  sprintf(fileNm, fmt, name, "gyyj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gyyj, fileNm);
  sprintf(fileNm, fmt, name, "gxzj");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gxzj, fileNm);
  sprintf(fileNm, fmt, name, "eps2");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->eps2, fileNm);
}






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
// Reasonable Diagonal plate:
void shaped_pfunc_upper(double s, double* RZ){
    RZ[0] = 5.151 + (5.852 - 5.151)*s;
    RZ[1] = 8.516 + (8.434 - 8.516)*s;
}

void shaped_pfunc_lower(double s, double* RZ){
    RZ[0] = 5.151 + (5.852 - 5.151)*s;
    RZ[1] = -(8.516 + (8.434 - 8.516)*s);
}

//X-pt info
//upper xpt R, Z = 2.50498, 6.14226
//lower xpt R, Z = 2.50504, -6.14214





void
test_fixed_z()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_tok_geo_efit_inp inp = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/step.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = false,
      .quad_param = {  .eps = 1e-10 }
    };

  //double clower[] = { 0.934, -0.01, -3.14 };
  //double cupper[] = {1.0, 0.01, 3.14 };

  double psisep = 1.5098198350000001;
  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };

  int ccells[] = { 1, 1, 8 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rright = 6.2,
    .rleft = 1.1,
    .zmin = -8.3,
    .zmax = 8.3,
    .rmin = 1.1,
    .rmax = 6.2,
  }; 
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .tok_efit_info = &inp,
    .tok_grid_info = &ginp,
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
  //write_geometry(up, cgrid, clocal, "step_outboard_fixed_z");

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



  struct gkyl_tok_geo_efit_inp inp = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/step.geqdsk",
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

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
  }; 

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .tok_efit_info = &inp,
    .tok_grid_info = &ginp,
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
  printf("total time = %g\n", cpu_time_used);
}

void
test_vertical_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_tok_geo_efit_inp inp = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/step.geqdsk",
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

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
  }; 
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .tok_efit_info = &inp,
    .tok_grid_info = &ginp,
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
  printf("total time = %g\n", cpu_time_used);
}

void
test_shaped_plate()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_tok_geo_efit_inp inp = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/step.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = true,
      // can set plate func if you dont want a fixed zmin and zmax
      .plate_func_lower = shaped_pfunc_lower,
      .plate_func_upper = shaped_pfunc_upper,
      .quad_param = {  .eps = 1e-10 }
    };

  //double clower[] = { 0.934, -0.01, -3.14 };
  //double cupper[] = {1.4688, 0.01, 3.14 };

  double psisep = 1.5098198350000001;
  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };

  int ccells[] = { 2, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
  }; 

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .tok_efit_info = &inp,
    .tok_grid_info = &ginp,
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
  printf("total time = %g\n", cpu_time_used);
}

TEST_LIST = {
  { "test_fixed_z", test_fixed_z},
  //{ "test_horizontal_plate", test_horizontal_plate},
  //{ "test_vertical_plate", test_vertical_plate},
  //{ "test_shaped_plate", test_shaped_plate},
  { NULL, NULL },
};
