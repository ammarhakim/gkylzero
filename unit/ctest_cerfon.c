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



void
test_deep_core()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();


  struct gkyl_tok_geo_efit_inp inp = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/cerfon.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = false,
      .quad_param = {  .eps = 1e-10 }
    };


  double clower[] = {0.006+3e-10, -0.01, -M_PI + 1e-14 };
  double cupper[] = {0.2, 0.01, M_PI - 1e-14 };

  int ccells[] = { 1, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  printf("CGRID INFO:\n cgrid.lower = %g,%g,%g\n cgrid.upper = %g,%g,%g\n cgrid.dx= %g,%g,%g\n", cgrid.lower[0],cgrid.lower[1], cgrid.lower[2],cgrid.upper[0],cgrid.upper[1], cgrid.upper[2], cgrid.dx[0], cgrid.dx[1], cgrid.dx[2]);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);


  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_CORE,
    .rclose = 6.0,
    .rleft= 0.25,
    .rright= 6.0,
    .zmin = -4.4,
    .zmax = 4.4,

    .write_node_coord_array = true,
    .node_file_nm = "cerfoncore_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 


  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_boundary()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();


  struct gkyl_tok_geo_efit_inp inp = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/cerfon.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = false,
      .quad_param = {  .eps = 1e-10 }
    };


  double clower[] = { 0.0, -0.01, -M_PI+1e-14 };
  double cupper[] = {0.1, 0.01, M_PI-1e-14 };

  int ccells[] = { 1, 1, 32 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  printf("CGRID INFO:\n cgrid.lower = %g,%g,%g\n cgrid.upper = %g,%g,%g\n cgrid.dx= %g,%g,%g\n", cgrid.lower[0],cgrid.lower[1], cgrid.lower[2],cgrid.upper[0],cgrid.upper[1], cgrid.upper[2], cgrid.dx[0], cgrid.dx[1], cgrid.dx[2]);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);


  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_CORE,
    .rclose = 6.0,
    .rleft= 0.25,
    .rright= 6.0,
    .zmin = -4.3,
    .zmax = 4.3,

    .write_node_coord_array = true,
    .node_file_nm = "cerfonbry_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 


  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_outboard()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();



  struct gkyl_tok_geo_efit_inp inp = {
      // psiRZ and related inputs
      .filepath = "./data/eqdsk/cerfon.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = false,
      .quad_param = {  .eps = 1e-10 }
    };

  double clower[] = { -0.2, -0.01, -M_PI+1e-14 };
  double cupper[] = {-0.1, 0.01, M_PI-1e-14 };

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
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.0,
    .zmin = -5.8,
    .zmax = 5.8,
    .write_node_coord_array = true,
    .node_file_nm = "cerfonoutboard_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}


TEST_LIST = {
  //{ "test_deep_core", test_deep_core},
  //{ "test_boundary", test_boundary},
  { "test_outboard", test_outboard},
  { NULL, NULL },
};
