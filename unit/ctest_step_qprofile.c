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
#include <gkyl_nodal_ops.h>




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
  sprintf(fileNm, fmt, name, "qprofile");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->qprofile, fileNm);
}





//// Coords that worked for outer
//.zmin = -6.14213,
//.zmax = 6.14226,
//
//Coords that worked for core

//double psisep = 1.5098198350000001;
//double psisep = 1.5093121789263824;
//double psisep = 1.5093073363672997;
double psisep = 1.5093065418975686;

// After dg reflect
//double psisep = 1.5093065418975689;


//from python
//double Zxpt_up = 6.167215740711948;
//double Zxpt_up = 6.1672157229985665;
//double Rxpt_up = 2.490007067421558;

//double Zxpt_lo = -6.16726668549029;
//double Rxpt_lo = 2.4901129589826647;
//
//double Zxpt_up = 6.16726668549029;
//double Rxpt_up = 2.4901129589826647;

// After Dg reflect:
double Zxpt_lo = -6.1672666854902927;
double Rxpt_lo = 2.4901129589826714;

double Zxpt_up = 6.1672666854902927;
double Rxpt_up = 2.4901129589826714;


// Somer vertical and horizontal plates
//void shaped_pfunc_lower_inner(double s, double* RZ){
//  RZ[0] = 1.71;
//  RZ[1] = -6.2 - s;
//}
//
//void shaped_pfunc_upper_inner(double s, double* RZ){
//  RZ[0] = 2.235;
//  RZ[1] = 6.267 + s;
//}

void shaped_pfunc_lower_outer(double s, double* RZ){
  RZ[0] = 3.5+2.0*s;
  RZ[1] = -8.29;
}

void shaped_pfunc_upper_outer(double s, double* RZ){
  RZ[0] = 3.5+2.0*s;
  RZ[1] = 8.29;
}

// Actual Plate info (Inboard Upper Plate):
// p1 = [1.651,6.331]
// p2 = [1.8, 6.777]
// Try a Different (slanted) plate instead
// p1 [5.151,8.516]
// p2 [5.852, 8.434]
// Reasonable Diagonal plate:
void shaped_pfunc_upper_inner(double s, double* RZ){
    RZ[0] = 1.651 + (1.8 - 1.651)*s;
    RZ[1] = 6.331 + (6.777 - 6.331)*s;
}

void shaped_pfunc_lower_inner(double s, double* RZ){
    RZ[0] = 1.65 + (1.8 - 1.65)*s;
    RZ[1] = -(6.33 + (6.777 - 6.33)*s);
}

// Actual Plate info:
// p1 = [5.7413,8.5471]
// p2 = [5.8544, 8.5229]
// p3 = [5.8549, 8.4258]
// Try a Different (slanted) plate instead
// p1 [5.151,8.516]
// p2 [5.852, 8.434]
// Reasonable Diagonal plate:
//void shaped_pfunc_upper_outer(double s, double* RZ){
//    RZ[0] = 5.151 + (5.852 - 5.151)*s;
//    RZ[1] = 8.516 + (8.434 - 8.516)*s;
//}
//
//void shaped_pfunc_lower_outer(double s, double* RZ){
//    RZ[0] = 5.151 + (5.852 - 5.151)*(2*s-1);
//    RZ[1] = -(8.516 + (8.434 - 8.516)*(2*s-1));
//}

struct gkyl_efit_inp inp_inner= {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

struct gkyl_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

struct gkyl_efit_inp inp_outer = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/step.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
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

  //double clower[] = { 1.6, -0.01, -M_PI+3e-16 };
  //double cupper[] = {6.0, 0.01, M_PI-3e-16 };
  double clower[] = { psisep, -0.01, -M_PI+3e-16 };
  double cupper[] = {6.0, 0.01, M_PI-3e-16 };
  int ccells[] = { 16,1,4 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_CORE,
    .rclose = 6.2,
    .rleft= 1.1,
    .rright= 6.2,
    .rmin=1.5,
    .rmax=6.2,
    .plate_spec = true,
    .plate_func_lower = shaped_pfunc_lower_outer,
    .plate_func_upper = shaped_pfunc_upper_outer,
  }; 
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .efit_info = inp_outer,
    .tok_grid_info = ginp,
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
  // Create Nodal Range and Grid and Write Nodal Coordinates
  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &clocal, cpoly_order);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&cbasis, &cgrid, false);
  gkyl_nodal_ops_m2n(n2m, &cbasis, &cgrid, &nrange, &clocal, 3, mc2p_nodal, up->mc2p);
  gkyl_nodal_ops_release(n2m);
  struct gkyl_rect_grid ngrid;
  gkyl_gk_geometry_init_nodal_grid(&ngrid, &cgrid, &nrange);
  gkyl_grid_sub_array_write(&ngrid, &nrange, 0,  mc2p_nodal, "stepcore_nodes.gkyl");
  gkyl_array_release(mc2p_nodal);
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}

void
test_sol()
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 0.934, -0.01, -M_PI+1e-14 };
  double cupper[] = {psisep-0.01, 0.01, M_PI-1e-14 };
  int ccells[] = { 16,1,4 };

  gkyl_rect_grid_init(&cgrid, 3, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cnghost, &clocal_ext, &clocal);
  gkyl_cart_modal_serendip(&cbasis, 3, cpoly_order);

  struct gkyl_tok_geo_grid_inp ginp = {
    .ftype = GKYL_SOL_DN_OUT,
    .rclose = 6.2,
    .rleft= 1.1,
    .rright= 6.2,
    .rmin=1.5,
    .rmax=6.2,
    .plate_spec = true,
    .plate_func_lower = shaped_pfunc_lower_outer,
    .plate_func_upper = shaped_pfunc_upper_outer,
  }; 
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_TOKAMAK,
    .efit_info = inp_outer,
    .tok_grid_info = ginp,
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
  write_geometry(up, cgrid, clocal, "stepsol");
  // Create Nodal Range and Grid and Write Nodal Coordinates
  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &clocal, cpoly_order);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&cbasis, &cgrid, false);
  gkyl_nodal_ops_m2n(n2m, &cbasis, &cgrid, &nrange, &clocal, 3, mc2p_nodal, up->mc2p);
  gkyl_nodal_ops_release(n2m);
  struct gkyl_rect_grid ngrid;
  gkyl_gk_geometry_init_nodal_grid(&ngrid, &cgrid, &nrange);
  gkyl_grid_sub_array_write(&ngrid, &nrange, 0,  mc2p_nodal, "stepsol_nodes.gkyl");
  gkyl_array_release(mc2p_nodal);
  gkyl_gk_geometry_release(up);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("total time = %g\n", cpu_time_used);
}


TEST_LIST = {
  { "test_core", test_core},
  { "test_sol", test_sol},
  { NULL, NULL },
};
