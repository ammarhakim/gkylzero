#include "math.h"

#include <acutest.h>
// #include <gkyl_array.h>
#include <gkyl_array_rio.h>
// #include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_dg_fpo_vlasov_diff_coeff.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff.h>
#include <gkyl_proj_on_basis.h>
// #include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
// #include <gkyl_util.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <time.h>

// Allocate array filled with zeros
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Project a Gaussian in 1x3v
void eval_1x3v_fn(double t, const double *xn, double* restrict fout, void *ctx)
{
  double vx = xn[1], vy = xn[2], vz = xn[3];

  fout[0] = exp(-vx*vx)*exp(-vy*vy)*exp(-vz*vz);
}

void eval_1x3v_fn_1st_derivatives(double t, const double *xn, double* restrict fout, void *ctx)
{
  double vx = xn[1], vy = xn[2], vz = xn[3];

  fout[0] = -2.0*vx*exp(-vx*vx)*exp(-vy*vy)*exp(-vz*vz);
  fout[1] = -2.0*vy*exp(-vx*vx)*exp(-vy*vy)*exp(-vz*vz);
  fout[2] = -2.0*vz*exp(-vx*vx)*exp(-vy*vy)*exp(-vz*vz);

}

// Project the diagonal terms of the diffusion tensor for a Gaussian
void eval_1x3v_fn_2nd_derivatives(double t, const double *xn, double* restrict fout, void *ctx)
{
  double vx = xn[1], vy = xn[2], vz = xn[3];

  fout[0] = (-2.0 + 4.0*vx*vx)*exp(-vx*vx)*exp(-vy*vy)*exp(-vz*vz);
  fout[1] = (-2.0 + 4.0*vy*vy)*exp(-vx*vx)*exp(-vy*vy)*exp(-vz*vz);
  fout[2] = (-2.0 + 4.0*vz*vz)*exp(-vx*vx)*exp(-vy*vy)*exp(-vz*vz);
}

void
test(int Nv, int poly_order, double eps)
{
  double L = 5.0;
  int ndim = 4;
  int ncells[] = {2, Nv, Nv, Nv};

  double lower[ndim], upper[ndim];
  int cells[ndim], ghost[ndim];
  for (int n=0; n<ndim; ++n) {
    lower[n] = -L;
    upper[n] = L;
    cells[n] = ncells[n];
    ghost[n] = 0;
  }

  // Grid
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Local, local-ext phase-space ranges
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project the function and analytic derivatives
  int num_quads = 8;
  struct gkyl_proj_on_basis *proj_fn, *proj_1st_derivatives, *proj_2nd_derivatives;
  proj_fn = gkyl_proj_on_basis_new(&grid, &basis,
    num_quads, 1, eval_1x3v_fn, NULL);

  proj_1st_derivatives = gkyl_proj_on_basis_new(&grid, &basis,
    num_quads, 3, eval_1x3v_fn_1st_derivatives, NULL);

  proj_2nd_derivatives = gkyl_proj_on_basis_new(&grid, &basis,
    num_quads, 3, eval_1x3v_fn_2nd_derivatives, NULL);

  // Create arrays for Gaussian and derivative
  struct gkyl_array *fn, *fn_1st_deriv, *fn_2nd_deriv;
  struct gkyl_array *drag_coeff, *diff_coeff;
  struct gkyl_array *drag_coeff_diff, *diff_coeff_diff;
  fn = mkarr(basis.num_basis, local.volume);
  fn_1st_deriv = mkarr(3*basis.num_basis, local.volume);
  fn_2nd_deriv = mkarr(3*basis.num_basis, local.volume);

  drag_coeff = mkarr(3*basis.num_basis, local.volume);
  diff_coeff = mkarr(3*basis.num_basis, local.volume);
 
  drag_coeff_diff = mkarr(3*basis.num_basis, local.volume);
  diff_coeff_diff = mkarr(3*basis.num_basis, local.volume); 

  // Project Gaussian and derivatives onto grid
  gkyl_proj_on_basis_advance(proj_fn, 0.0, &local, fn);
  gkyl_proj_on_basis_advance(proj_1st_derivatives, 0.0, &local, fn_1st_deriv);
  gkyl_proj_on_basis_advance(proj_2nd_derivatives, 0.0, &local, fn_2nd_deriv);
  
  // Call updaters to calculate drag and diffusion coefficients
  struct timespec tm1 = gkyl_wall_clock();
  gkyl_calc_fpo_drag_coeff_recovery(&grid, basis, &local, fn, drag_coeff);
  double drag_tm = gkyl_time_diff_now_sec(tm1);

  printf("\nDrag coefficient on (%d)^3 took %g sec\n", cells[1], drag_tm);

  struct timespec tm2 = gkyl_wall_clock();
  gkyl_calc_fpo_diff_coeff_recovery(&grid, basis, &local, fn, diff_coeff);
  double diff_tm = gkyl_time_diff_now_sec(tm2);

  printf("Diffusion coefficient diagonals on (%d)^3 took %g sec\n", cells[1], diff_tm);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);

  // Compare derivatives from recovery with analytic solutions
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);

    // Check drag coefficient components against 1st derivatives
    const double *drag_coeff_d = gkyl_array_cfetch(drag_coeff, linidx);
    const double *fn_1st_deriv_d = gkyl_array_cfetch(fn_1st_deriv, linidx);

    double *drag_coeff_diff_d = gkyl_array_fetch(drag_coeff_diff, linidx);

    for (int m=0; m<basis.num_basis; ++m) {
      TEST_CHECK(gkyl_compare(drag_coeff_d[m], fn_1st_deriv_d[m], eps));
      TEST_MSG("Expected %.13e, coefficient (%d) in cell (%d, %d, %d, %d)", fn_1st_deriv_d[m], m, iter.idx[0], iter.idx[1], iter.idx[2], iter.idx[3]);
      TEST_MSG("Produced: %.13e, coefficient (%d)", drag_coeff_d[m], m);
      drag_coeff_diff_d[m] = drag_coeff_d[m] - fn_1st_deriv_d[m];
    }

    // Check diffusion coefficient components against 2nd derivatives
    const double *diff_coeff_d = gkyl_array_cfetch(diff_coeff, linidx);
    const double *fn_2nd_deriv_d = gkyl_array_cfetch(fn_2nd_deriv, linidx);

    double *diff_coeff_diff_d = gkyl_array_fetch(diff_coeff_diff, linidx);

    for (int m=0; m<basis.num_basis; ++m) {
      TEST_CHECK(gkyl_compare(diff_coeff_d[m], fn_2nd_deriv_d[m], eps));
      TEST_MSG("Expected %.13e, coefficient (%d) in cell (%d, %d, %d, %d)", fn_2nd_deriv_d[m], m, iter.idx[0], iter.idx[1], iter.idx[2], iter.idx[3]);
      TEST_MSG("Produced: %.13e, coefficient (%d)", diff_coeff_d[m], m);
      diff_coeff_diff_d[m] = diff_coeff_d[m] - fn_2nd_deriv_d[m];
    }
  }

  // Write output files
  char fname1[64], fname2[64], fname3[64], fname4[64];
  sprintf(fname1, "unit-out/ctest_drag_coeff_difference_%dc_p%d.gkyl", Nv, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, drag_coeff_diff, fname1);

  sprintf(fname2, "unit-out/ctest_diff_coeff_difference_%dc_p%d.gkyl", Nv, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, diff_coeff_diff, fname2);
  //
  // sprintf(fname3, "unit-out/ctest_drag_coeff_%dc_p%d.gkyl", Nv, poly_order);
  // gkyl_grid_sub_array_write(&grid, &local, drag_coeff, fname3);
  //
  // sprintf(fname4, "unit-out/ctest_diff_coeff_%dc_p%d.gkyl", Nv, poly_order);
  // gkyl_grid_sub_array_write(&grid, &local, diff_coeff, fnam43);

  gkyl_array_release(fn);
  gkyl_array_release(fn_1st_deriv);
  gkyl_array_release(fn_2nd_deriv);
  gkyl_array_release(drag_coeff);
  gkyl_array_release(diff_coeff);
  gkyl_array_release(drag_coeff_diff);
  gkyl_array_release(diff_coeff_diff);

  gkyl_proj_on_basis_release(proj_fn);
  gkyl_proj_on_basis_release(proj_1st_derivatives); 
  gkyl_proj_on_basis_release(proj_2nd_derivatives);
}

void test_1x3v_8_p1() { test(8, 1, 1.0); }
void test_1x3v_16_p1() { test(16, 1, 1.0); }
void test_1x3v_32_p1() { test(32, 1, 1.0); }
void test_1x3v_64_p1() { test(64, 1, 1.0); }

void test_1x3v_8_p2() { test(8, 2, 1.0); }
void test_1x3v_16_p2() { test(16, 2, 1.0); }
void test_1x3v_32_p2() { test(32, 2, 1.0); }
void test_1x3v_64_p2() { test(64, 2, 1.0); }

TEST_LIST = {
  { "test_1x3v_8_p1", test_1x3v_8_p1 },
  { "test_1x3v_16_p1", test_1x3v_16_p1 },
  { "test_1x3v_32_p1", test_1x3v_32_p1 },
  { "test_1x3v_64_p1", test_1x3v_64_p1 },

  { "test_1x3v_8_p2", test_1x3v_8_p2 },
  { "test_1x3v_16_p2", test_1x3v_16_p2 },
  { "test_1x3v_32_p2", test_1x3v_32_p2 },
  { "test_1x3v_64_p2", test_1x3v_64_p2 },

  { NULL, NULL },
};
