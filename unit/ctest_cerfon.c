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


  struct gkyl_tok_geo_efit_inp inp = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/cerfon.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = false,
      .quad_param = {  .eps = 1e-10 }
    };


void
test_11()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_CORE_R,
    .rclose = 6.0,
    .rright= 6.0,
    .zxpt_lo = -4.3,
    .zxpt_up = 4.3,

    .write_node_coord_array = true,
    .node_file_nm = "cerfon11_nodes.gkyl"
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

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_CORE_L,
    .rclose = 6.0,
    .rleft= 0.25,
    .zxpt_lo = -4.3,
    .zxpt_up = 4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon12_nodes.gkyl"
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

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_PF_LO_R,
    .rright = 6.0,
    .zmin = -5.8,
    .zxpt_lo = -4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon1_nodes.gkyl"
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

  double clower[] = { -0.2, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_SOL_DN_OUT_LO,
    .rright = 6.0,
    .zmin = -5.8,
    .zxpt_lo = -4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon2_nodes.gkyl"
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

  double clower[] = { -0.2, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_SOL_DN_OUT_MID,
    .rright = 6.0,
    .zxpt_lo = -4.3,
    .zxpt_up = 4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon3_nodes.gkyl"
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

  double clower[] = { -0.2, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_SOL_DN_OUT_UP,
    .rright = 6.0,
    .zmax = 5.8,
    .zxpt_up = 4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon4_nodes.gkyl"
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

  double clower[] = { -0.01, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_SOL_DN_IN_UP,
    .rleft = 0.25,
    .zxpt_up = 4.3,
    .zmax = 5.8,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon7_nodes.gkyl"
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

  double clower[] = { -0.01, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_SOL_DN_IN_MID,
    .rleft = 0.25,
    .zxpt_lo = -4.3,
    .zxpt_up = 4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon8_nodes.gkyl"
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

  double clower[] = { -0.01, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.0, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_SOL_DN_IN_LO,
    .rleft = 0.25,
    .zmin = -5.8,
    .zxpt_lo = -4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon9_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_10()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_PF_LO_L,
    .rleft= 0.25,
    .zmin = -5.8,
    .zxpt_lo = -4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon10_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_5()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_PF_UP_R,
    .rright = 6.0,
    .zmax = 5.8,
    .zxpt_up = 4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon5_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_6()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_PF_UP_L,
    .rleft = 0.25,
    .zmax = 5.8,
    .zxpt_up = 4.3,
    .write_node_coord_array = true,
    .node_file_nm = "cerfon6_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}









TEST_LIST = {
  { "test_1", test_1}, //cmag not so great
  { "test_2", test_2},
  { "test_3", test_3},
  { "test_11", test_11},
  { "test_12", test_12},
  { "test_4", test_4},
  {"test_7", test_7},
  {"test_8", test_8},
  {"test_9", test_9},
  {"test_10", test_10}, //cmag not so great
  {"test_5", test_5}, //cmag not so great
  {"test_6", test_6}, //cmag pretty good but not great
  { NULL, NULL },
};