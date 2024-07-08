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
  sprintf(fileNm, fmt, name, "normals");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->normals, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobgeo, fileNm);
  sprintf(fileNm, fmt, name, "jacobgeo_inv");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->jacobgeo_inv, fileNm);
  sprintf(fileNm, fmt, name, "gij");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->gij, fileNm);
  sprintf(fileNm, fmt, name, "b_i");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->b_i, fileNm);
  sprintf(fileNm, fmt, name, "bcart");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->bcart, fileNm);
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





//// Coords that worked for outer
//.zmin = -6.14213,
//.zmax = 6.14226,
//
//Coords that worked for core
//.zxpt_lo = -6.14214,
//.zxpt_up = 6.1423,

//double psisep = 1.5098198350000001;
double psisep = 1.5093121789263824;

//from python
//double Zxpt_up = 6.167215740711948;
double Zxpt_up = 6.1672157229985665;
double Rxpt_up = 2.490007067421558;

double Zxpt_lo = -6.167210302704182;
double Rxpt_lo = 2.4900075684870013;

struct gkyl_tok_geo_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/step.geqdsk",
    .rzpoly_order = 2,
    //.rz_basis_type = GKYL_BASIS_MODAL_TENSOR,
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
  int ccells[] = { 4, 1, 8 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_CORE_R,
    .rclose = 6.2,
    .rleft= 1.1,
    .rright= 6.2,
    .rmin=1.5,
    .rmax=6.2,
    //.zxpt_lo = -6.1,
    //.zxpt_up = 6.1,
    .zxpt_lo = Zxpt_lo,
    .zxpt_up = Zxpt_up,

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
  write_geometry(up, cgrid, clocal, "stepcore");
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_core_l()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { psisep, -0.01, -M_PI+1e-14 };
  double cupper[] = {1.8, 0.01, M_PI-1e-14 };
  int ccells[] = { 4, 1, 8 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_CORE_L,
    .rclose = 1.1,
    .rleft= 1.1,
    .rright= 6.2,
    .rmin=1.5,
    .rmax=6.2,
    //.zxpt_lo = -6.1,
    //.zxpt_up = 6.1,
    .zxpt_lo = Zxpt_lo,
    .zxpt_up = Zxpt_up,

    .write_node_coord_array = true,
    .node_file_nm = "stepcorel_nodes.gkyl"
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
  write_geometry(up, cgrid, clocal, "stepcorel");
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
  int ccells[] = { 4, 1, 8 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rright = 6.2,
    .rleft = 1.1,
    .rmin = 1.1,
    .rmax = 6.2,
    //.zmin = -6.1,
    //.zmax = 6.1,
    .zmin = 0.0,
    .zmax = Zxpt_up,
    .write_node_coord_array = true,
    .node_file_nm = "stepouter_nodes.gkyl"
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
  write_geometry(up, cgrid, clocal, "stepouter");
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_upper()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep, 0.01, M_PI-1e-14 };
  int ccells[] = { 4, 1, 8 };



  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT_UP,
    .rright = 6.2,
    .rleft = 1.1,
    .rmin = 1.1,
    .rmax = 6.2,
    .zxpt_up = Zxpt_up,
    .zmax = 8.3,
    .write_node_coord_array = true,
    .node_file_nm = "stepupper_nodes.gkyl"
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
  write_geometry(up, cgrid, clocal, "stepupper");
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}



TEST_LIST = {
  //{ "test_core", test_core},
  //{ "test_core_l", test_core_l},
  { "test_outer", test_outer},
  //{ "test_upper", test_upper},
  { NULL, NULL },
};
