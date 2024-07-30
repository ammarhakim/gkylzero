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
#include <gkyl_mirror_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mirror.h>


void
test_lores()
{
  struct gkyl_mirror_geo_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/wham.geqdsk",
    .rzpoly_order = 2,
    .fluxpoly_order = 1,
    .plate_spec = false,
    .quad_param = {  .eps = 1e-10 }
  };

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 1e-3, -0.01, -M_PI+1e-14 };
  double cupper[] = { 3e-3,  0.01,  M_PI-1e-14 };

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);


struct gkyl_mirror_geo_grid_inp ginp = {
  .rclose = 0.2,
  .zmin = -2.48,
  .zmax =  2.48,
};

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_MIRROR,
    .mirror_efit_info = &inp,
    .mirror_grid_info = &ginp,
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

  struct gk_geometry* up = gkyl_gk_geometry_mirror_new(&geometry_inp); 


  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_hires()
{
  struct gkyl_mirror_geo_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/wham_hires.geqdsk",
    .rzpoly_order = 2,
    .fluxpoly_order = 1,
    .plate_spec = false,
    .quad_param = {  .eps = 1e-10 }
  };

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 1e-3, -0.01, -M_PI+1e-14 };
  double cupper[] = { 3e-3,  0.01,  M_PI-1e-14 };

  int ccells[] = { 1, 1, 8 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);

  struct gkyl_range clocal, clocal_ext;
  int cnghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);


struct gkyl_mirror_geo_grid_inp ginp = {
  .rclose = 0.2,
  .zmin = -2.48,
  .zmax =  2.48,
};

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_MIRROR,
    .mirror_efit_info = &inp,
    .mirror_grid_info = &ginp,
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

  struct gk_geometry* up = gkyl_gk_geometry_mirror_new(&geometry_inp); 


  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

TEST_LIST = {
  // { "test_lores", test_lores },
  { "test_hires", test_hires },
  { NULL, NULL },
};
