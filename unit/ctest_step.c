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
#include <gkyl_tok_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>



//// Coords that worked for outer
//.zmin = -6.14213,
//.zmax = 6.14226,
//
//Coords that worked for core
//.zxpt_lo = -6.14214,
//.zxpt_up = 6.1423,

double psisep = 1.5098198350000001;
struct gkyl_tok_geo_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./efit_data/step.geqdsk",
    .rzpoly_order = 2,
    .fluxpoly_order = 1,
    .plate_spec = false,
    .quad_param = {  .eps = 1e-10 }
  };

int cpoly_order = 1;
struct gkyl_basis cbasis;
int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
struct gkyl_rect_grid cgrid;
struct gkyl_range clocal, clocal_ext;

void
test_core()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { psisep, -0.01, -M_PI+1e-14 };
  double cupper[] = {1.8, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 16 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_CORE,
    .rclose = 6.2,
    .rleft= 1.1,
    .rright= 6.2,
    //.zxpt_lo = -6.2,
    //.zxpt_up = 6.2,
    //.zxpt_lo = -6.14214,
    //.zxpt_up = 6.1423,
    .zxpt_lo = -6.142,
    .zxpt_up = 6.142,

    .write_node_coord_array = true,
    .node_file_nm = "stepcore_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_outer()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 16 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rright = 6.2,
    .zmin = -8.3,
    .zmax = 8.3,
    .write_node_coord_array = true,
    .node_file_nm = "stepouter_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_2()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 24 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT_LO,
    .rright = 6.2,
    .zmin = -8.3,
    .zxpt_lo = -6.5,
    .write_node_coord_array = true,
    .node_file_nm = "step2_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_3()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 8 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT_MID,
    .rright= 6.2,
    .zxpt_lo = -6.142,
    .zxpt_up = 6.142,
    .write_node_coord_array = true,
    .node_file_nm = "step3_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_4()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 8 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT_UP,
    .rright = 6.2,
    .zxpt_up = 6.142,
    .zmax = 8.3,
    .write_node_coord_array = true,
    .node_file_nm = "step4_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

TEST_LIST = {
  //{ "test_core", test_core},
  //{ "test_outer", test_outer}, //Works
  //{"test_3", test_3}, //Works
  {"test_2", test_2}, //Does not work cmag not 1
  //{"test_4", test_4}, // same as 2
  { NULL, NULL },
};
