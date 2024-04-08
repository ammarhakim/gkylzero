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
      .filepath = "./data/eqdsk/step.geqdsk",
      .rzpoly_order = 2,
      .fluxpoly_order = 1,
      .plate_spec = false,
      .quad_param = {  .eps = 1e-10 }
    };


  double clower[] = { 2.01, -0.01, -M_PI + 1e-14 };
  double cupper[] = {2.1, 0.01, M_PI - 1e-14 };

  int ccells[] = { 1, 1, 16 };



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
    .rclose = 6.2,
    .rleft= 1.1,
    .rright= 6.2,
    .zxpt_lo = -6.2,
    .zxpt_up = 6.2,
    .rmin = 1.5,
    .rmax = 6.2,

    .write_node_coord_array = true,
    .node_file_nm = "stepcore_nodes.gkyl"
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
test_boundary()
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

  double psisep = 1.5098198350000001;
  double clower[] = { psisep, -0.01, -M_PI+1e-14 };
  double cupper[] = {1.8, 0.01, M_PI-1e-14 };
  //double cupper[] = {1.8, 0.01, 0.197};
  //double cupper[] = {1.8, 0.01, 0.197109};

  //double clower[] = { 1.5098198350000001, -0.01, -2.8 };
  //double cupper[] = {1.8, 0.01, 2.8 };

  int ccells[] = { 1, 1, 8 };



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
    .rclose = 6.2,
    .rleft= 1.1,
    .rright= 6.2,
    .zxpt_lo = -6.2,
    .zxpt_up = 6.2,
    .rmin = 1.5,
    .rmax = 6.2,

    .write_node_coord_array = true,
    .node_file_nm = "stepbry_nodes.gkyl"
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
  //{ "test_deep_core", test_deep_core},
  { "test_boundary", test_boundary},
  { NULL, NULL },
};
