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
//.zmin = -1.093941629526259813,
//.zmax = 1.094605624589544826,
//
//Coords that worked for core
//.zxpt_lo = -1.093941629526259814,
//.zxpt_up = 1.09460562458954483,

double psisep = - .12501253600000001 ;
struct gkyl_tok_geo_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./efit_data/mast.geqdsk",
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
  double cupper[] = {-0.11, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 16 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_CORE,
    .rclose = 2.0,
    .rleft= 0.5,
    .rright= 2.0,
    .zxpt_lo = -1.0939416295262598,
    .zxpt_up = 1.0946056245895448,

    .write_node_coord_array = true,
    .node_file_nm = "mastcore_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_11()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { psisep, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.11, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 8 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_CORE_R,
    .rclose = 2.0,
    .rleft= 0.5,
    .rright= 2.0,
    .zxpt_lo = -1.0939416295262598,
    .zxpt_up = 1.0946056245895448,

    .write_node_coord_array = true,
    .node_file_nm = "mast11_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_12()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { psisep, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.11, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 8 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_CORE_L,
    .rclose = 2.0,
    .rleft = 0.5,
    .rright= 2.0,
    .zxpt_lo = -1.0939416295262598,
    .zxpt_up = 1.0946056245895448,

    .write_node_coord_array = true,
    .node_file_nm = "mast12_nodes.gkyl"
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

  double clower[] = { -0.13, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 16 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rright = 2.0,
    .zmin = -1.7,
    .zmax = 1.7,
    .write_node_coord_array = true,
    .node_file_nm = "mastouter_nodes.gkyl"
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

  double clower[] = { -0.13, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 8 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT_LO,
    .rright = 2.0,
    .zmin = -1.7,
    //.zxpt_lo = -1.0939416295262598,
    .zxpt_lo = -1.0925,
    .write_node_coord_array = true,
    .node_file_nm = "mast2_nodes.gkyl"
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

  double clower[] = { -0.13, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 8 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT_MID,
    .rright= 2.0,
    .zxpt_lo = -1.0925,
    .zxpt_up = 1.092,
    .write_node_coord_array = true,
    .node_file_nm = "mast3_nodes.gkyl"
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

  double clower[] = { -0.13, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 8 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT_UP,
    .rright = 2.0,
    //.zxpt_up = 1.0946056245895448,
    .zxpt_up = 1.092,
    .zmax = 1.7,
    .write_node_coord_array = true,
    .node_file_nm = "mast4_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_8()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { -0.13, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 8 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_IN_MID,
    .rleft = 0.5,
    .zxpt_lo = -1.0939416295262598,
    .zxpt_up = 1.0946056245895448,
    .write_node_coord_array = true,
    .node_file_nm = "mast8_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_7()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { -0.13, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 4 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_IN_UP,
    .rleft = 0.5,
    .zxpt_up = 1.0946056245895448,
    .zmax = 1.3,
    .write_node_coord_array = true,
    .node_file_nm = "mast7_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_9()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { -0.13, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 4 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_IN_LO,
    .rleft = 2.0,
    .zxpt_lo = -1.0939416295262598,
    .zmin = -1.3,
    .write_node_coord_array = true,
    .node_file_nm = "mast9_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_1()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { psisep, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.11, 0.01, M_PI-1e-14 };
  int ccells[] = { 1, 1, 4 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_PF_LO_R,
    .rright = 2.0,
    .zxpt_lo = -1.0939416295262598,
    .zmin = -8.4,
    .write_node_coord_array = true,
    .node_file_nm = "mast1_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

TEST_LIST = {
  //{ "test_core", test_core},
  //{ "test_outer", test_outer}, // Works
  {"test_3", test_3}, // Works. Good cmag
  {"test_2", test_2}, // cmag -> 0 at xpt  but ok
  {"test_4", test_4}, // cmag ->0 at xpt but ok
  //{"test_11", test_11}, // Works. Good cmag. Even with nup hack for PF
  //{"test_12", test_12}, // Works. Good cmag. Even with nup and nlo hack for PF
  //{"test_8", test_8},  // Works. Good cmag. Verified correct orientation
  //{"test_7", test_7}, // Works. cmag->0 at xpt but ok
  //{"test_9", test_9}, // Works. cmag-<0 at xpt but ok
  //{"test_1", test_1},  // Nodes work but cmag looks bad. Great if zxpt_lo=2.

  { NULL, NULL },
};
