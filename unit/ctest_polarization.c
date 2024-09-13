#include "gkyl_array.h"
#include <stdio.h>
#include <stdlib.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <acutest.h>
#include <gkyl_gyrokinetic_pol_density.h>

void evalFunc1x_1(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

void evalFunc1x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = -pow(xn[0]-0.5, 2) + 1.0;
}

void evalFunc2x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = (-pow(xn[0]-0.5, 2) + 1.0) * (-pow(xn[1]-0.5, 2) + 1.0);
}

void evalFunc3x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = (-pow(xn[0]-0.5, 2) + 1.0) * (-pow(xn[1]-0.5, 2) + 1.0) * (-pow(xn[2]-0.5, 2) + 1.0);
}

void test_1x_flat() 
{  
  int cells[] = {8};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, false);
  struct gkyl_array *npol = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, 0.0, &localRange, epsilon);
  gkyl_eval_on_nodes_release(epsilon_proj);

  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, 1, poly_order+1);
  struct gkyl_array *phi_pol = gkyl_array_new(GKYL_DOUBLE, phi_pol_basis.num_basis, localRange_ext.volume);
  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&grid, &phi_pol_basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &localRange, phi_pol);
  gkyl_eval_on_nodes_release(phi_pol_proj);

  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol, linidx);
    for (int i=0; i<basis.num_basis; ++i) {
      TEST_CHECK( npol_d[i] == 0.0 );
    }
  }
  gkyl_array_release(npol);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void test_1x_quad() 
{  
  int cells[] = {8};
  int poly_order = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, false);
  struct gkyl_array *npol = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, 0.0, &localRange, epsilon);
  gkyl_eval_on_nodes_release(epsilon_proj);

  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, 1, poly_order+1);
  struct gkyl_array *phi_pol = gkyl_array_new(GKYL_DOUBLE, phi_pol_basis.num_basis, localRange_ext.volume);
  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&grid, &phi_pol_basis,
    1, evalFunc1x_quad, NULL);
  gkyl_eval_on_nodes_advance(phi_pol_proj, 0.0, &localRange, phi_pol);
  gkyl_eval_on_nodes_release(phi_pol_proj);

  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol, linidx);
    TEST_CHECK( gkyl_compare(npol_d[0], 2.8284271247462, 1e-14) );
    TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-12) );
  }

  gkyl_array_release(npol);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void
test_2x_quad() 
{
  int cells[] = {7, 7};
  int poly_order = 1;
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  double time = 0.0;
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, false);
  struct gkyl_array *npol = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, time, &localRange, epsilon);
  gkyl_eval_on_nodes_release(epsilon_proj);

  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, dim, poly_order+1);
  struct gkyl_array *phi_pol = gkyl_array_new(GKYL_DOUBLE, phi_pol_basis.num_basis, localRange_ext.volume);
  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&grid, &phi_pol_basis,
    1, evalFunc2x_quad, NULL);
  gkyl_eval_on_nodes_advance(phi_pol_proj, time, &localRange, phi_pol);
  gkyl_eval_on_nodes_release(phi_pol_proj);

  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol, linidx);
    if (conf_iter.idx[1] == 1) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.2585034013604615, 1e-13) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-14) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.1413919026586975, 1e-13) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, 1e-14) );
    } else if (conf_iter.idx[1] == 2) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.6666666666666150, 1e-13) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-14) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.0942612684391317, 1e-12) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, 1e-14) );
    } else if (conf_iter.idx[1] == 3) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.9115646258503931, 1e-12) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-14) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.0471306342195649, 1e-12) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, 1e-14) );
    } else if (conf_iter.idx[1] == 4) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.9931972789115648, 1e-12) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-14) );
      TEST_CHECK( gkyl_compare(npol_d[2], 0.0, 1e-12) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, 1e-14) );
    } else if (conf_iter.idx[1] == 5) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.9115646258503931, 1e-12) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-14) );
      TEST_CHECK( gkyl_compare(npol_d[2], -0.0471306342195649, 1e-12) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, 1e-14) );
    } else if (conf_iter.idx[1] == 6) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.6666666666666150, 1e-13) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-14) );
      TEST_CHECK( gkyl_compare(npol_d[2], -0.0942612684391317, 1e-12) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, 1e-14) );
    } else if (conf_iter.idx[1] == 7) {
      TEST_CHECK( gkyl_compare(npol_d[0], 3.2585034013604615, 1e-13) );
      TEST_CHECK( gkyl_compare(npol_d[1], 0.0, 1e-14) );
      TEST_CHECK( gkyl_compare(npol_d[2], -0.1413919026586975, 1e-13) );
      TEST_CHECK( gkyl_compare(npol_d[3], 0.0, 1e-14) );
    }
  }

  gkyl_array_release(npol);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

void
test_3x_flat() 
{
  int cells[] = {7, 7, 7};
  int poly_order = 1;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  double time = 0.0;
  int dim = sizeof(lower)/sizeof(lower[0]);
  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, dim, lower, upper, cells);

  // Basis functions.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, dim, poly_order);

  int ghost[] = { 1, 1 };
  struct gkyl_range localRange, localRange_ext; // local, local-ext ranges.
  gkyl_create_grid_ranges(&grid, ghost, &localRange_ext, &localRange);

  struct gkyl_gyrokinetic_pol_density* npol_op = gkyl_gyrokinetic_pol_density_new(basis, grid, false);
  struct gkyl_array *npol = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);
  struct gkyl_array *epsilon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, localRange_ext.volume);

  struct gkyl_eval_on_nodes *epsilon_proj = gkyl_eval_on_nodes_new(&grid, &basis,
    1, evalFunc1x_1, NULL);
  gkyl_eval_on_nodes_advance(epsilon_proj, time, &localRange, epsilon);
  gkyl_eval_on_nodes_release(epsilon_proj);

  struct gkyl_basis phi_pol_basis;
  gkyl_cart_modal_tensor(&phi_pol_basis, dim, poly_order+1);
  struct gkyl_array *phi_pol = gkyl_array_new(GKYL_DOUBLE, phi_pol_basis.num_basis, localRange_ext.volume);
  struct gkyl_eval_on_nodes *phi_pol_proj = gkyl_eval_on_nodes_new(&grid, &phi_pol_basis,
    1, evalFunc1x_quad, NULL);
  gkyl_eval_on_nodes_advance(phi_pol_proj, time, &localRange, phi_pol);
  gkyl_eval_on_nodes_release(phi_pol_proj);

  gkyl_gyrokinetic_pol_density_advance(npol_op, &localRange, epsilon, phi_pol, npol);
  // // Save the arrays to look at
  const char *fmt = "%s-%s.gkyl";
  char name[32] = "pol";
  int sz = gkyl_calc_strlen(fmt, name, "jacobtot_inv");
  char fileNm[sz+1]; // ensure no buffer overflow
  sprintf(fileNm, fmt, name, "epsilon");
  gkyl_grid_sub_array_write(&grid, &localRange, 0, epsilon, fileNm);
  sprintf(fileNm, fmt, name, "phi");
  gkyl_grid_sub_array_write(&grid, &localRange, 0, phi_pol, fileNm);
  sprintf(fileNm, fmt, name, "npol");
  gkyl_grid_sub_array_write(&grid, &localRange, 0, npol, fileNm);

  // read the components of npol
  struct gkyl_range_iter conf_iter;
  gkyl_range_iter_init(&conf_iter, &localRange);
  while (gkyl_range_iter_next(&conf_iter)) {
    long linidx = gkyl_range_idx(&localRange, conf_iter.idx);
    double *npol_d = gkyl_array_fetch(npol, linidx);
    double tol = 1e-12;
    TEST_CHECK( gkyl_compare(npol_d[0], 5.6568542494927714, tol) );
    TEST_CHECK( gkyl_compare(npol_d[1], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[2], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[3], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[4], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[5], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[6], 0.0, tol) );
    TEST_CHECK( gkyl_compare(npol_d[7], 0.0, tol) );
  }

  gkyl_array_release(npol);
  gkyl_array_release(epsilon);
  gkyl_array_release(phi_pol);
  gkyl_gyrokinetic_pol_density_release(npol_op);
}

TEST_LIST = {
  { "test_1x_flat", test_1x_flat }, 
  { "test_1x_quad", test_1x_quad },
  { "test_2x_quad", test_2x_quad },
  { "test_3x_flat", test_3x_flat},
  { NULL, NULL },
};