#include "gkyl_util.h"
#include <acutest.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>

void f_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2 + x;
}

void g_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2*x*x + 8;
}

void
test_1d(int poly_order)
{
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_1d, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, g_1d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);  

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  struct gkyl_array *h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h, 0, distf, 0, distg);
  // f_bar = h/g = f
  gkyl_dg_div_op(basis, 0, f_bar, 0, h, 0, distg);
  // g_bar = h/f = g
  gkyl_dg_div_op(basis, 0, g_bar, 0, h, 0, distf);

  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
  }

  // Test range methods
  gkyl_array_clear(h, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(basis, 0, f_bar, 0, h, 0, distg, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(basis, 0, g_bar, 0, h, 0, distf, arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
  }

  // mean ops
  struct gkyl_array *mvals = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);
  gkyl_array_clear(mvals, 0.0);

  // means are stored in h[0]
  gkyl_dg_calc_average_range(basis, 0, mvals, 0, distf, arr_range);
  // L2 are stored in h[1]
  gkyl_dg_calc_l2_range(basis, 1, mvals, 0, distf, arr_range);

  double al2[2];
  gkyl_array_reduce_range(al2, mvals, GKYL_SUM, arr_range);

  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 2.5, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 19.0/3.0, 1e-14) );
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(h);
  gkyl_array_release(mvals);
}

void
test_mean_ops_1d(int poly_order)
{
}

void f_2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  fout[0] = 2 + x + y;
}

void g_2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  fout[0] = 2*x*y + 8;
}

void
test_2d(int poly_order)
{
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int cells[] = {2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_2d, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, g_2d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  struct gkyl_array *h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h, 0, distf, 0, distg);
  // f_bar = h/g = f
  gkyl_dg_div_op(basis, 0, f_bar, 0, h, 0, distg);
  // g_bar = h/f = g
  gkyl_dg_div_op(basis, 0, g_bar, 0, h, 0, distf);

  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
  }

  // Test range methods
  gkyl_array_clear(h, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(basis, 0, f_bar, 0, h, 0, distg, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(basis, 0, g_bar, 0, h, 0, distf, arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
  }

  // mean ops
  struct gkyl_array *mvals = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);
  gkyl_array_clear(mvals, 0.0);

  // means are stored in h[0]
  gkyl_dg_calc_average_range(basis, 0, mvals, 0, distf, arr_range);
  // L2 are stored in h[1]
  gkyl_dg_calc_l2_range(basis, 1, mvals, 0, distf, arr_range);

  double al2[2];
  gkyl_array_reduce_range(al2, mvals, GKYL_SUM, arr_range);

  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 3.0, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 55.0/6.0, 1e-14) );  
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(h);
  gkyl_array_release(mvals);
}

void f_3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = 5 + x + y + z;
}

void g_3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = ((x*y*z + 8)*(x*y*z + 8) + 100);
}

void fg_3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];

  double f[1], g[1];
  f_3d(t, xn, f, ctx);
  g_3d(t, xn, g, ctx);
  
  fout[0] = f[0]*g[0];
}


void
test_3d(int poly_order)
{
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, f_3d, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, g_3d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  struct gkyl_array *h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h, 0, distf, 0, distg);
  // f_bar = h/g = f
  gkyl_dg_div_op(basis, 0, f_bar, 0, h, 0, distg);
  // g_bar = h/f = g
  gkyl_dg_div_op(basis, 0, g_bar, 0, h, 0, distf);

  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }

  // Test range methods
  gkyl_array_clear(h, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(basis, 0, f_bar, 0, h, 0, distg, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(basis, 0, g_bar, 0, h, 0, distf, arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }

  // mean ops
  struct gkyl_array *mvals = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);
  gkyl_array_clear(mvals, 0.0);

  // means are stored in h[0]
  gkyl_dg_calc_average_range(basis, 0, mvals, 0, distf, arr_range);
  // L2 are stored in h[1]
  gkyl_dg_calc_l2_range(basis, 1, mvals, 0, distf, arr_range);

  double al2[2];
  gkyl_array_reduce_range(al2, mvals, GKYL_SUM, arr_range);

  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 6.5, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 85.0/2.0, 1e-14) );
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(h);
  gkyl_array_release(mvals);
}

void
test_1d_p1()
{
  test_1d(1);
}

void
test_1d_p2()
{
  test_1d(2);
}

void
test_1d_p3()
{
  test_1d(3);
}

void
test_2d_p1()
{
  test_2d(1);
}

void
test_2d_p2()
{
  test_2d(2);
}

void
test_2d_p3()
{
  test_2d(3);
}

void
test_3d_p1()
{
  test_3d(1);
}

void
test_3d_p2()
{
  test_3d(2);
}

void f_3d_p3(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = ((5 + x)*(5 + y)*(5 + z)*(5 + x)*(5 + y)*(5 + z)*(5 + x)*(5 + y)*(5 + z) + 100);
}

void g_3d_p3(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = ((x*y*z + 8)*(x*y*z + 8)*(x*y*z + 8) + 100);
}

void
test_3d_p3()
{
  int poly_order = 3;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, f_3d_p3, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, g_3d_p3, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  
  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  struct gkyl_array *h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h, 0, distf, 0, distg);
  // f_bar = h/g = f
  gkyl_dg_div_op(basis, 0, f_bar, 0, h, 0, distg);
  // g_bar = h/f = g
  gkyl_dg_div_op(basis, 0, g_bar, 0, h, 0, distf);

  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }

  // Test range methods
  gkyl_array_clear(h, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(basis, 0, f_bar, 0, h, 0, distg, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(basis, 0, g_bar, 0, h, 0, distf, arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(h);
}

// Cuda specific tests
#ifdef GKYL_HAVE_CUDA

void
test_1d_cu(int poly_order)
{
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_1d, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, g_1d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // make device copies of arrays
  struct gkyl_array *distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);  

  // project distribution function on basis on CPU
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  // copy host array to device
  gkyl_array_copy(distf_cu, distf);
  gkyl_array_copy(distg_cu, distg);

  // Product array only needs to be initialized on GPU
  struct gkyl_array *h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h_cu, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
  // f_bar = h/g = f
  gkyl_dg_div_op(basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
  // g_bar = h/f = g
  gkyl_dg_div_op(basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
  }

  // Test range methods
  gkyl_array_clear(h_cu, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);
  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
  }
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(distf_cu);
  gkyl_array_release(distg_cu);
  gkyl_array_release(f_bar_cu);
  gkyl_array_release(g_bar_cu);
  gkyl_array_release(h_cu);
}

void
test_2d_cu(int poly_order)
{
  double lower[] = {0.0, 0.0}, upper[] = {1.0, 1.0};
  int cells[] = {2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_2d, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, g_2d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // make device copies of arrays
  struct gkyl_array *distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);  

  // project distribution function on basis on CPU
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  // copy host array to device
  gkyl_array_copy(distf_cu, distf);
  gkyl_array_copy(distg_cu, distg);

  // Product array only needs to be initialized on GPU
  struct gkyl_array *h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h_cu, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
  // f_bar = h/g = f
  gkyl_dg_div_op(basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
  // g_bar = h/f = g
  gkyl_dg_div_op(basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
  }

  // Test range methods
  gkyl_array_clear(h_cu, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);
  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-12) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-12) );
    }
  }
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(distf_cu);
  gkyl_array_release(distg_cu);
  gkyl_array_release(f_bar_cu);
  gkyl_array_release(g_bar_cu);
  gkyl_array_release(h_cu);
}

void
test_3d_cu(int poly_order)
{
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, f_3d, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, g_3d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // make device copies of arrays
  struct gkyl_array *distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);  

  // project distribution function on basis on CPU
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  // copy host array to device
  gkyl_array_copy(distf_cu, distf);
  gkyl_array_copy(distg_cu, distg);

  // Product array only needs to be initialized on GPU
  struct gkyl_array *h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h_cu, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
  // f_bar = h/g = f
  gkyl_dg_div_op(basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
  // g_bar = h/f = g
  gkyl_dg_div_op(basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }

  // Test range methods
  gkyl_array_clear(h_cu, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);
  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(distf_cu);
  gkyl_array_release(distg_cu);
  gkyl_array_release(f_bar_cu);
  gkyl_array_release(g_bar_cu);
  gkyl_array_release(h_cu);
}

void
test_1d_p1_cu()
{
  test_1d_cu(1);
}

void
test_1d_p2_cu()
{
  test_1d_cu(2);
}

void
test_1d_p3_cu()
{
  test_1d_cu(3);
}

void
test_2d_p1_cu()
{
  test_2d_cu(1);
}

void
test_2d_p2_cu()
{
  test_2d_cu(2);
}

void
test_2d_p3_cu()
{
  test_2d_cu(3);
}

void
test_3d_p1_cu()
{
  test_3d_cu(1);
}

void
test_3d_p2_cu()
{
  test_3d_cu(2);
}

void
test_3d_p3_cu()
{
  int poly_order = 3;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {1.0, 1.0, 1.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, f_3d_p3, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, 5, 1, g_3d_p3, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // make device copies of arrays
  struct gkyl_array *distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);  

  // project distribution function on basis on CPU
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  // copy host array to device
  gkyl_array_copy(distf_cu, distf);
  gkyl_array_copy(distg_cu, distg);

  // Product array only needs to be initialized on GPU
  struct gkyl_array *h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(h_cu, 0.0);

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);

  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(g_bar, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
  // f_bar = h/g = f
  gkyl_dg_div_op(basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
  // g_bar = h/f = g
  gkyl_dg_div_op(basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  for (size_t i=0; i<arr_range.volume; ++i) {
    const double *f_d = gkyl_array_cfetch(distf, i);
    const double *fbar_d = gkyl_array_cfetch(f_bar, i);
    const double *g_d = gkyl_array_cfetch(distg, i);
    const double *gbar_d = gkyl_array_cfetch(g_bar, i);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }

  // Test range methods
  gkyl_array_clear(h_cu, 0.0);
  gkyl_array_clear(f_bar_cu, 0.0);
  gkyl_array_clear(g_bar_cu, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);
  // copy from device and check if things are ok
  gkyl_array_copy(f_bar, f_bar_cu);
  gkyl_array_copy(g_bar, g_bar_cu);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&arr_range, iter.idx);
    const double *f_d = gkyl_array_cfetch(distf, loc);
    const double *fbar_d = gkyl_array_cfetch(f_bar, loc);
    const double *g_d = gkyl_array_cfetch(distg, loc);
    const double *gbar_d = gkyl_array_cfetch(g_bar, loc);
    for (int k=0; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare(f_d[k], fbar_d[k], 1e-10) );
      TEST_CHECK( gkyl_compare(g_d[k], gbar_d[k], 1e-10) );
    }
  }
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(distf_cu);
  gkyl_array_release(distg_cu);
  gkyl_array_release(f_bar_cu);
  gkyl_array_release(g_bar_cu);
  gkyl_array_release(h_cu);
}

#endif

TEST_LIST = {
  { "test_1d_p1", test_1d_p1 },
  { "test_1d_p2", test_1d_p2 },
  { "test_1d_p3", test_1d_p3 },
  { "test_2d_p1", test_2d_p1 },
  { "test_2d_p2", test_2d_p2 },
  { "test_2d_p3", test_2d_p3 },
  { "test_3d_p1", test_3d_p1 },
  { "test_3d_p2", test_3d_p2 },
  { "test_3d_p3", test_3d_p3 },
#ifdef GKYL_HAVE_CUDA
  { "test_1d_p1_cu", test_1d_p1_cu },
  { "test_1d_p2_cu", test_1d_p2_cu },
  { "test_1d_p3_cu", test_1d_p3_cu },
  { "test_2d_p1_cu", test_2d_p1_cu },
  { "test_2d_p2_cu", test_2d_p2_cu },
  { "test_2d_p3_cu", test_2d_p3_cu },
  { "test_3d_p1_cu", test_3d_p1_cu },
  { "test_3d_p2_cu", test_3d_p2_cu },
  { "test_3d_p3_cu", test_3d_p3_cu },
#endif
  { NULL, NULL },
};
