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
test_1()
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


  double clower[] = { 2.01, -0.01, -3.14 };
  double cupper[] = {2.1, 0.01, 3.14 };

  int ccells[] = { 4, 1, 64 };



  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  printf("CGRID INFO:\n cgrid.lower = %g,%g,%g\n cgrid.upper = %g,%g,%g\n cgrid.dx= %g,%g,%g\n", cgrid.lower[0],cgrid.lower[1], cgrid.lower[2],cgrid.upper[0],cgrid.upper[1], cgrid.upper[2], cgrid.dx[0], cgrid.dx[1], cgrid.dx[2]);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);


  struct gkyl_tok_geo_geo_inp ginp = {
    .ftype = GKYL_CORE,
    .rclose = 6.2,
    .rleft= 1.1,
    .rright= 6.2,
    .zmin = -6.2,
    .zmax = 6.2,

    .write_node_coord_array = true,
    .node_file_nm = "stepccore_nodes.gkyl"
  }; 

  struct gk_geometry* up = gkyl_gk_geometry_tok_new(&cgrid, &clocal, &clocal_ext, &cbasis, &inp, &ginp, false); 


  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

TEST_LIST = {
  { "test_1", test_1},
  { NULL, NULL },
};
