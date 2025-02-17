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
#include <gkyl_nodal_ops.h>

void
write_geometry(gk_geometry *up, struct gkyl_rect_grid grid, struct gkyl_range local, const char *name)
{
  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, name, "jacobtot_inv");
  char fileNm[sz+1]; // ensure no buffer overflow

  sprintf(fileNm, fmt, name, "mapc2p");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->mc2p, fileNm);
  sprintf(fileNm, fmt, name, "mapc2nu");
  gkyl_grid_sub_array_write(&grid, &local, 0,  up->mc2nu_pos, fileNm);
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


  // Create Nodal Range and Grid and Write Nodal Coordinates
  struct gkyl_range nrange;
  gkyl_gk_geometry_init_nodal_range(&nrange, &local, up->basis.poly_order);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->basis, &grid, false);
  gkyl_nodal_ops_m2n(n2m, &up->basis, &grid, &nrange, &local, 3, mc2p_nodal, up->mc2p);
  gkyl_nodal_ops_release(n2m);
  struct gkyl_rect_grid ngrid;
  gkyl_gk_geometry_init_nodal_grid(&ngrid, &grid, &nrange);
  sprintf(fileNm, fmt, name, "nodes");
  gkyl_grid_sub_array_write(&ngrid, &nrange, 0,  mc2p_nodal, fileNm);
  gkyl_array_release(mc2p_nodal);
}

void
test_lores()
{
  struct gkyl_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/wham.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };

  clock_t start, end;
  double cpu_time_used;
  start = clock();

  double clower[] = { 1e-10, -0.01, -M_PI+1e-14 };
  double cupper[] = { 2e-3,  0.01,  M_PI-1e-14 };

  int ccells[] = { 4, 1, 8 };

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
    .zmin = -2.4,
    .zmax =  2.4,
  };

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_MIRROR,
    .efit_info = inp,
    .mirror_grid_info = ginp,
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
  struct gkyl_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/wham_hires.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
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
  .zmin = -2.0,
  .zmax =  2.0,
};

  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = GKYL_MIRROR,
    .efit_info = inp,
    .mirror_grid_info = ginp,
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

// def psi_f(R, Z):
//     Bmag = 0.5
//     return Bmag/2 * R**2

// Functions for this test
void mapc2p(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double psi = xn[0], alpha = xn[1], zeta = xn[2];
  fout[0] = sqrt(psi * 4 ); // Function fed is psi = 0.5/2 * R^2 from the efit file
  fout[1] = zeta * 1.0 / M_PI; // Note that this does not have pi-1e-2 in it because the coordinate zeta is always defined -pi to pi
  fout[2] = -alpha; // There is a minus due to conventions
}

void exact_gij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r*r/4;
  fout[0] = 1/psi; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = r*r; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = 1/(M_PI*M_PI); // g_33
}

void exact_g_contra_ij(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double r = xn[0], theta = xn[1], phi = xn[2];
  double psi = r*r/4;
  fout[0] = r*r/4; // g_11
  fout[1] = 0.0; // g_12
  fout[2] = 0.0; // g_13
  fout[3] = 1/psi/4; // g_22
  fout[4] = 0.0; // g_23
  fout[5] = (M_PI*M_PI); // g_33
}

void bmag_func(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx){
  fout[0] = 0.5;
}

void
test_3x_p1_straight_cylinder()
{
  // Very similar to the unit test in ctest_gk_geometry.c
  // The geometry is created to extend from Z = -1 to 1, R = (0.001, 1) in units meters
  // Magnetic field = 0.5 uniform everywhere
  // Psi = 0.5/2 * R^2
  // Check write_efit_straight_cylinder.py for how the efit file is written
  // There are some important differences between how the numerical geometry is calculated compared to mapc2p geometries

  struct gkyl_basis basis;
  int poly_order = 1;
  int cdim = 3;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);
  
  double psiMax = 0.2;
  double psiMin = 1e-4;
  int Nz = 10;

  double lower[3] = {psiMin, -M_PI, -M_PI+1e-2};
  double upper[3] = {psiMax,  M_PI,  M_PI-1e-2};
  // int cells[3] = { 18, 18, Nz };
  int cells[3] = { 8, 8, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  
  struct gkyl_range ext_range, range;
  int nghost[3] = { 1,1,1};
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_efit_inp inp = {
    // psiRZ and related inputs
    .filepath = "./data/eqdsk/straight_cylinder.geqdsk",
    .rz_poly_order = 2,
    .flux_poly_order = 1,
    .reflect = true,
  };
  struct gkyl_mirror_geo_grid_inp ginp = {
    .rclose = 0.5,
    .zmin = -1.,
    .zmax =  1.,
  };
  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_MIRROR,
    .efit_info = inp,
    .mirror_grid_info = ginp,
    .grid = grid,
    .local = range,
    .local_ext = ext_range,
    .global = range,
    .global_ext = ext_range,
    .basis = basis,
    .geo_grid = grid,
    .geo_local = range,
    .geo_local_ext = ext_range,
    .geo_global = range,
    .geo_global_ext = ext_range,
    .geo_basis = basis,
  };

  struct gk_geometry *gk_geom = gkyl_gk_geometry_mirror_new(&geometry_input);

  write_geometry(gk_geom, grid, range, "straight_cylinder");

  // Define nodal operations
  enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
  int cidx[3];
  int nodes[] = { 1, 1, 1 };
  for (int d=0; d<grid.ndim; ++d)
    nodes[d] = grid.cells[d] + 1;
  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);

  // Check that |bhat|=1 at nodes
  struct gkyl_array* bhat_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, bhat_nodal, gk_geom->bcart);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *bhat_n = gkyl_array_fetch(bhat_nodal, gkyl_range_idx(&nrange, cidx));
        double bhat_mag = sqrt(bhat_n[0]*bhat_n[0] + bhat_n[1]*bhat_n[1] + bhat_n[2]*bhat_n[2]);
        TEST_CHECK( gkyl_compare( bhat_mag, 1.0, 1e-12) );
      }
    }
  }

  // Check that the duals are what they should be
  // Broken
  struct gkyl_array* dualmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, dualmag_nodal, gk_geom->dualmag);
  struct gkyl_array* mapc2p_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, mapc2p_nodal, gk_geom->mc2p);
  // for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
  //   for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
  //     for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
  //       cidx[PSI_IDX] = ip;
  //       cidx[AL_IDX] = ia;
  //       cidx[TH_IDX] = it;
  //       double *dualmag_n = gkyl_array_fetch(dualmag_nodal, gkyl_range_idx(&nrange, cidx));
  //       double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
  //       double e2mag = 1/mapc2p_n[0]; // 1/R
  //       TEST_CHECK( gkyl_compare( dualmag_n[0], 1.0, 1e-8) );
  //       TEST_CHECK( gkyl_compare( dualmag_n[1], e2mag, 1e-8) );
  //       TEST_CHECK( gkyl_compare( dualmag_n[2], 1.0, 1e-8) );
  //     }
  //   }
  // }


  // Check bmag is what it should be
  struct gkyl_array* bmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, bmag_nodal, gk_geom->bmag);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
        double xn[3] = {psi, alpha, theta};
        double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(&nrange, cidx));
        double bmag_anal[1];
        bmag_func(0, xn, bmag_anal, 0);
        TEST_CHECK( gkyl_compare( bmag_n[0], bmag_anal[0], 1e-8) );
      }
    }
  }

  // Check bmag_inv_sq = 1/B^2
  struct gkyl_array* bmag_inv_sq_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, bmag_inv_sq_nodal, gk_geom->bmag_inv_sq);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *bmag_inv_sq_n = gkyl_array_fetch(bmag_inv_sq_nodal, gkyl_range_idx(&nrange, cidx));
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
        double xn[3] = {psi, alpha, theta};
        double bmag_anal[1];
        bmag_func(0, xn, bmag_anal, 0);
        TEST_CHECK( gkyl_compare( bmag_inv_sq_n[0], 1/(bmag_anal[0]*bmag_anal[0]), 1e-6) );
      }
    }
  }

  // Check bmag_inv = 1/B
  struct gkyl_array* bmag_inv_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, bmag_inv_nodal, gk_geom->bmag_inv);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *bmag_inv_n = gkyl_array_fetch(bmag_inv_nodal, gkyl_range_idx(&nrange, cidx));
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
        double xn[3] = {psi, alpha, theta};
        double bmag_anal[1];
        bmag_func(0, xn, bmag_anal, 0);
        TEST_CHECK( gkyl_compare( bmag_inv_n[0], 1/bmag_anal[0], 1e-8) );
      }
    }
  }

  // Check cmag = 1
  struct gkyl_array* cmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, cmag_nodal, gk_geom->cmag);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *cmag_n = gkyl_array_fetch(cmag_nodal, gkyl_range_idx(&nrange, cidx));
        TEST_CHECK( gkyl_compare( cmag_n[0], 1.0, 1e-8) );
      }
    }
  }

  // Check eps2 
  // What is this?
  struct gkyl_array* eps2_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, eps2_nodal, gk_geom->eps2);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *eps2_n = gkyl_array_fetch(eps2_nodal, gkyl_range_idx(&nrange, cidx));
        TEST_CHECK( gkyl_compare( eps2_n[0], 0.0, 1e-8) );
      }
    }
  }

  // Check g_ij
  struct gkyl_array* gij_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 6, gij_nodal, gk_geom->g_ij);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *gij_n = gkyl_array_fetch(gij_nodal, gkyl_range_idx(&nrange, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double r = mapc2p_n[0];
        double xn[3] = {r, 0.0, 0.0};
        double fout[6];
        exact_gij(0.0, xn, fout, 0);
        for (int i=0; i<6; ++i)
          TEST_CHECK( gkyl_compare( gij_n[i], fout[i], 1e-6) );
      }
    }
  }

  // Check g^ij
  struct gkyl_array* gij_contra_nodal = gkyl_array_new(GKYL_DOUBLE, 6, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 6, gij_contra_nodal, gk_geom->gij);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *gij_contra_n = gkyl_array_fetch(gij_contra_nodal, gkyl_range_idx(&nrange, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double r = mapc2p_n[0];
        double xn[3] = {r, 0.0, 0.0};
        double fout[6];
        exact_g_contra_ij(0.0, xn, fout, 0);
        for (int i=0; i<6; ++i)
          TEST_CHECK( gkyl_compare( gij_contra_n[i], fout[i], 1e-6) );
      }
    }
  }

  // Check that Jacobgeo is what it should be. This is the Jacobian for the problem
  struct gkyl_array* jacobgeo_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, jacobgeo_nodal, gk_geom->jacobgeo);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobgeo_n = gkyl_array_fetch(jacobgeo_nodal, gkyl_range_idx(&nrange, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double jacobian_analytic = 2/M_PI;
        TEST_CHECK( gkyl_compare( jacobgeo_n[0], jacobian_analytic, 1e-6) );
      }
    }
  }

  // Check jacobgeo_inv
  struct gkyl_array* jacobgeo_inv_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, jacobgeo_inv_nodal, gk_geom->jacobgeo_inv);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double *jacobgeo_inv_n = gkyl_array_fetch(jacobgeo_inv_nodal, gkyl_range_idx(&nrange, cidx));
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double jacobian_analytic = 2/M_PI;
        TEST_CHECK( gkyl_compare( jacobgeo_inv_n[0], 1/jacobian_analytic, 1e-6) );
      }
    }
  }

  // Check jacobtot
  struct gkyl_array* jacobtot_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, jacobtot_nodal, gk_geom->jacobtot);
  // for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
  //   for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
  //     for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
  //       cidx[PSI_IDX] = ip;
  //       cidx[AL_IDX] = ia;
  //       cidx[TH_IDX] = it;
  //       double *jacobtot_n = gkyl_array_fetch(jacobtot_nodal, gkyl_range_idx(&nrange, cidx));
  //       // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
  //       double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
  //       double radius = mapc2p_n[0];
  //       TEST_CHECK( gkyl_compare( jacobtot_n[0], radius, 1e-8) );
  //     }
  //   }
  // }

  // Check jacobtot_inv
  struct gkyl_array* jacobtot_inv_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, jacobtot_inv_nodal, gk_geom->jacobtot_inv);
  // for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
  //   for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
  //     for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
  //       cidx[PSI_IDX] = ip;
  //       cidx[AL_IDX] = ia;
  //       cidx[TH_IDX] = it;
  //       double *jacobtot_inv_n = gkyl_array_fetch(jacobtot_inv_nodal, gkyl_range_idx(&nrange, cidx));
  //       // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
  //       double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
  //       double radius = mapc2p_n[0];
  //       TEST_CHECK( gkyl_compare( jacobtot_inv_n[0], 1/radius, 1e-8) );
  //     }
  //   }
  // }

  // Check mapc2p
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX]; 
        // mapc2p_n[0] = R, mapc2p_n[1] = Theta, mapc2p_n[2] = Z_cylindrical
        double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
        double xn[3] = {psi, alpha, theta};
        double fout[3];
        mapc2p(0.0, xn, fout, 0);
        for (int i=0; i<3; ++i)
          TEST_CHECK( gkyl_compare( mapc2p_n[i], fout[i], 1e-8) );
      }
    }
  }

  // Check mc2nu_pos
  struct gkyl_array* mc2nu_pos_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, mc2nu_pos_nodal, gk_geom->mc2nu_pos);
  for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
    for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
      for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
        cidx[PSI_IDX] = ip;
        cidx[AL_IDX] = ia;
        cidx[TH_IDX] = it;
        double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
        double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
        double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
        double xn[3] = {psi, alpha, theta};
        double *mc2nu_pos_n = gkyl_array_fetch(mc2nu_pos_nodal, gkyl_range_idx(&nrange, cidx));
        for (int i=0; i<3; ++i)
          TEST_CHECK( gkyl_compare( mc2nu_pos_n[i], xn[i], 1e-8) );
      }
    }
  }

  // Check normals
  struct gkyl_array* normals_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
  gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, normals_nodal, gk_geom->normals);
  // for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
  //   for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
  //     for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
  //       cidx[PSI_IDX] = ip;
  //       cidx[AL_IDX] = ia;
  //       cidx[TH_IDX] = it;
  //       double *normals_n = gkyl_array_fetch(normals_nodal, gkyl_range_idx(&nrange, cidx));
  //       double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
  //       double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
  //       double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
  //       double xn[3] = {psi, alpha, theta};
  //       double fout[3];
  //       normal(0.0, xn, fout, 0);
  //       for (int i=0; i<3; ++i)
  //         TEST_CHECK( gkyl_compare( normals_n[i], fout[i], 1e-8) );
  //     }
  //   }
  // }


  gkyl_array_release(bhat_nodal);
  gkyl_array_release(dualmag_nodal);
  gkyl_array_release(bmag_nodal);
  gkyl_array_release(bmag_inv_nodal);
  gkyl_array_release(bmag_inv_sq_nodal);
  gkyl_array_release(cmag_nodal);
  gkyl_array_release(eps2_nodal);
  gkyl_array_release(gij_nodal);
  gkyl_array_release(gij_contra_nodal);
  gkyl_array_release(jacobgeo_nodal);
  gkyl_array_release(jacobgeo_inv_nodal);
  gkyl_array_release(jacobtot_nodal);
  gkyl_array_release(jacobtot_inv_nodal);
  gkyl_array_release(mapc2p_nodal);
  gkyl_array_release(mc2nu_pos_nodal);
  gkyl_array_release(normals_nodal);
  gkyl_nodal_ops_release(n2m);
  gkyl_gk_geometry_release(gk_geom);
}

// void
// mapz(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
// {
//   double Lz = 1.8;
//   double a = -Lz/2;
//   fout[0] = -1/(2*a) * pow(a - xn[0], 2) + a;
// }

// void
// jacob_mapz(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
// {
//   double Lz = 1.8;
//   double a = -Lz/2;
//   fout[0] = 1 - xn[0]/a;
// }

// void
// test_3x_p1_geometry_quantities_position_map()
// {  
//   // Very similar to the unit test in ctest_gk_geometry.c
//   struct gkyl_basis basis;
//   int poly_order = 1;
//   int cdim = 3;
//   gkyl_cart_modal_serendip(&basis, 3, poly_order);
  
//   double Lz = 1.8;
//   double Rmax = 0.9;
//   double Rmin = 0.1;
//   int Nz = 10;

//   double lower[3] = {Rmin, -M_PI, -Lz/2};
//   double upper[3] = {Rmax,  M_PI,  Lz/2};
//   int cells[3] = { 18, 18, Nz };
//   struct gkyl_rect_grid grid;
//   gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  
//   struct gkyl_range ext_range, range;
//   int nghost[3] = { 1,1,1};
//   gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

//   struct gkyl_position_map_inp pos_map_inp = {
//     .maps = {0, 0, mapz},
//     .ctxs = {0, 0, 0},
//   };

//   // Configuration space geometry initialization
//   struct gkyl_position_map *pos_map = gkyl_position_map_new(pos_map_inp, grid, range, 
//     ext_range, range, ext_range, basis);

//   struct gkyl_efit_inp inp = {
//     // psiRZ and related inputs
//     .filepath = "./data/eqdsk/straight_cylinder.geqdsk",
//     .rz_poly_order = 2,
//     .flux_poly_order = 1,
//     .reflect = true,
//   };
//   struct gkyl_mirror_geo_grid_inp ginp = {
//     .rclose = 0.1,
//     .zmin = -1.0,
//     .zmax =  1.0,
//   };
//   // Initialize geometry
//   struct gkyl_gk_geometry_inp geometry_input = {
//     .geometry_id = GKYL_MIRROR,
//     .efit_info = inp,
//     .mirror_grid_info = ginp,
//     .position_map = pos_map,
//     .grid = grid,
//     .local = range,
//     .local_ext = ext_range,
//     .global = range,
//     .global_ext = ext_range,
//     .basis = basis,
//     .geo_grid = grid,
//     .geo_local = range,
//     .geo_local_ext = ext_range,
//     .geo_global = range,
//     .geo_global_ext = ext_range,
//     .geo_basis = basis,
//   };

//   struct gk_geometry *gk_geom = gkyl_gk_geometry_mirror_new(&geometry_input);

//   gkyl_position_map_set(pos_map, gk_geom->mc2nu_pos);

//   // Define nodal operations
//   enum { PSI_IDX, AL_IDX, TH_IDX }; // arrangement of computational coordinates
//   int cidx[3];
//   int nodes[] = { 1, 1, 1 };
//   for (int d=0; d<grid.ndim; ++d)
//     nodes[d] = grid.cells[d] + 1;
//   struct gkyl_range nrange;
//   gkyl_range_init_from_shape(&nrange, grid.ndim, nodes);
//   struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&basis, &grid, false);

//   // Check that Jacobgeo is what it should be. J = R in cylindrical coordinates
//   // We have a contribution from the position map too, given as dZ/dz
//   struct gkyl_array* jacobgeo_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
//   gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, jacobgeo_nodal, gk_geom->jacobgeo);
//   struct gkyl_array* mapc2p_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
//   gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 3, mapc2p_nodal, gk_geom->mc2p);
//   for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
//     for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
//       for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
//         cidx[PSI_IDX] = ip;
//         cidx[AL_IDX] = ia;
//         cidx[TH_IDX] = it;
//         double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
//         double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
//         double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
//         double *jacobgeo_n = gkyl_array_fetch(jacobgeo_nodal, gkyl_range_idx(&nrange, cidx));
//         // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
//         double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
//         double radius = sqrt(mapc2p_n[0]*mapc2p_n[0] + mapc2p_n[1]*mapc2p_n[1]);
//         double fout[1];
//         jacob_mapz(0.0, &theta, fout, 0);
//         double jacob_anal = radius * fout[0];
//         TEST_CHECK( gkyl_compare( jacobgeo_n[0], jacob_anal, 1e-8) );
//       }
//     }
//   }

//   // Check bmag is what it should be
//   struct gkyl_array* bmag_nodal = gkyl_array_new(GKYL_DOUBLE, grid.ndim, nrange.volume);
//   gkyl_nodal_ops_m2n(n2m, &basis, &grid, &nrange, &range, 1, bmag_nodal, gk_geom->bmag);
//   for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
//     for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
//       for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
//         cidx[PSI_IDX] = ip;
//         cidx[AL_IDX] = ia;
//         cidx[TH_IDX] = it;
//         double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
//         double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
//         double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX];
//         mapz(0.0, &theta, &theta, 0);
//         double xn[3] = {psi, alpha, theta};
//         double *bmag_n = gkyl_array_fetch(bmag_nodal, gkyl_range_idx(&nrange, cidx));
//         double bmag_anal[1];
//         bmag_func(0, xn, bmag_anal, 0);
//         TEST_CHECK( gkyl_compare( bmag_n[0], bmag_anal[0], 1e-8) );
//       }
//     }
//   }

//   // Check mapc2p
//   for (int ia=nrange.lower[AL_IDX]; ia<=nrange.upper[AL_IDX]; ++ia){
//     for (int ip=nrange.lower[PSI_IDX]; ip<=nrange.upper[PSI_IDX]; ++ip) {
//       for (int it=nrange.lower[TH_IDX]; it<=nrange.upper[TH_IDX]; ++it) {
//         cidx[PSI_IDX] = ip;
//         cidx[AL_IDX] = ia;
//         cidx[TH_IDX] = it;
//         double psi = grid.lower[PSI_IDX] + ip*(grid.upper[PSI_IDX]-grid.lower[PSI_IDX])/grid.cells[PSI_IDX];
//         double alpha = grid.lower[AL_IDX] + ia*(grid.upper[AL_IDX]-grid.lower[AL_IDX])/grid.cells[AL_IDX];
//         double theta = grid.lower[TH_IDX] + it*(grid.upper[TH_IDX]-grid.lower[TH_IDX])/grid.cells[TH_IDX]; 
//         // mapc2p_n[0] = x, mapc2p_n[1] = y, mapc2p_n[2] = z
//         double *mapc2p_n = gkyl_array_fetch(mapc2p_nodal, gkyl_range_idx(&nrange, cidx));
//         mapz(0.0, &theta, &theta, 0);
//         double xn[3] = {psi, alpha, theta};
//         double fout[3];
//         mapc2p(0.0, xn, fout, 0);
//         for (int i=0; i<3; ++i)
//           TEST_CHECK( gkyl_compare( mapc2p_n[i], fout[i], 1e-8) );
//       }
//     }
//   }
//   gkyl_array_release(bmag_nodal);
//   gkyl_array_release(jacobgeo_nodal);
//   gkyl_array_release(mapc2p_nodal);
//   gkyl_nodal_ops_release(n2m);
//   gkyl_gk_geometry_release(gk_geom);
// }





















TEST_LIST = {
  // { "test_lores", test_lores },
  //{ "test_hires", test_hires },
  { "test_3x_p1_straight_cylinder", test_3x_p1_straight_cylinder },
  // { "test_3x_p1_straight_cylinder_position_map", test_3x_p1_straight_cylinder_position_map },
  { NULL, NULL },
};
