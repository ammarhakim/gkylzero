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
test_1d(int poly_order, bool use_gpu)
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
  struct gkyl_array *distf_cu, *distg_cu;
  if (use_gpu) {
    distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  }

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(distf_cu, distf);
    gkyl_array_copy(distg_cu, distg);
  }

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);

  struct gkyl_array *h, *f_bar_cu, *g_bar_cu, *h_cu;
  if (use_gpu) {
    f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);

    // Product array only needs to be initialized on GPU
    h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h_cu, 0.0);
  } else {
    h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h, 0.0);
  }

  gkyl_dg_bin_op_mem *mem;
  if (use_gpu) {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);
  
    // h = f*g
    gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);
  
    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
  } else {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(basis, 0, h, 0, distf, 0, distg);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, basis, 0, f_bar, 0, h, 0, distg);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, basis, 0, g_bar, 0, h, 0, distf);
  }

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
<<<<<<< HEAD
  gkyl_array_clear(h, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, basis, 0, f_bar, 0, h, 0, distg, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, basis, 0, g_bar, 0, h, 0, distf, arr_range);
=======
  if (use_gpu) {
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(h_cu, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);

    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
  } else {
    gkyl_array_clear(f_bar, 0.0);
    gkyl_array_clear(g_bar, 0.0);
    gkyl_array_clear(h, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, basis, 0, f_bar, 0, h, 0, distg, arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, basis, 0, g_bar, 0, h, 0, distf, arr_range);
  }
>>>>>>> main

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
  double al2[2];
  if (use_gpu) {
    struct gkyl_array *mvals_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals_cu, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(basis, 0, mvals_cu, 0, distf_cu, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(basis, 1, mvals_cu, 0, distf_cu, arr_range);

    double* al2_cu = (double*) gkyl_cu_malloc(sizeof(double[2]));
    gkyl_array_reduce_range(al2_cu, mvals_cu, GKYL_SUM, arr_range);
  
    gkyl_cu_memcpy(al2, al2_cu, sizeof(double[2]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(al2_cu);
    gkyl_array_release(mvals_cu);
  } else {
    struct gkyl_array *mvals = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(basis, 0, mvals, 0, distf, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(basis, 1, mvals, 0, distf, arr_range);

    gkyl_array_reduce_range(al2, mvals, GKYL_SUM, arr_range);
    gkyl_array_release(mvals);
  }


  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 2.5, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 19.0/3.0, 1e-14) );
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_dg_bin_op_mem_release(mem);
  if (use_gpu) {
    gkyl_array_release(distf_cu);
    gkyl_array_release(distg_cu);
    gkyl_array_release(f_bar_cu);
    gkyl_array_release(g_bar_cu);
    gkyl_array_release(h_cu);
  } else {
    gkyl_array_release(h);
  }
}

void f_1d2d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1 + x;
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
test_2d(int poly_order, bool use_gpu)
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
  struct gkyl_array *distf_cu, *distg_cu;
  if (use_gpu) {
    distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  }

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(distf_cu, distf);
    gkyl_array_copy(distg_cu, distg);
  }

  // create conf-space grid:
  double clower[] = {lower[0]}, cupper[] = {upper[0]};
  int ccells[] = {cells[0]};
  int cdim = sizeof(clower)/sizeof(clower[0]);
  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, cdim, clower, cupper, ccells);
  struct gkyl_range arr_crange, arr_ext_crange;
  gkyl_create_grid_ranges(&cgrid, nghost, &arr_ext_crange, &arr_crange);
  // conf-space basis functions
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  // create a conf-space factor
  struct gkyl_array *cfield = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
  // project conf-space function onto basis.
  gkyl_proj_on_basis *proj_cfield = gkyl_proj_on_basis_new(&cgrid, &cbasis, poly_order+1, 1, f_1d2d, NULL);
  gkyl_proj_on_basis_advance(proj_cfield, 0.0, &arr_crange, cfield);
  struct gkyl_array *cfield_cu;
  if (use_gpu) {
    cfield_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
    gkyl_array_copy(cfield_cu, cfield);
  }

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *w_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(w_bar, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);

  struct gkyl_array *h, *f_bar_cu, *g_bar_cu, *w_bar_cu, *h_cu;
  if (use_gpu) {
    f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    w_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(w_bar_cu, 0.0);

    // Product array only needs to be initialized on GPU
    h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h_cu, 0.0);
  } else {
    h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h, 0.0);
  }
  
  gkyl_dg_bin_op_mem *mem;
  if (use_gpu) {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
  } else {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(basis, 0, h, 0, distf, 0, distg);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, basis, 0, f_bar, 0, h, 0, distg);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, basis, 0, g_bar, 0, h, 0, distf);
  }

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
<<<<<<< HEAD
  gkyl_array_clear(h, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, basis, 0, f_bar, 0, h, 0, distg, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, basis, 0, g_bar, 0, h, 0, distf, arr_range);
=======
  if (use_gpu) {
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(w_bar_cu, 0.0);
    gkyl_array_clear(h_cu, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(cbasis, basis, w_bar_cu, cfield_cu, distf_cu, arr_crange, arr_range);

    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
    gkyl_array_copy(w_bar, w_bar_cu);
  } else {
    gkyl_array_clear(f_bar, 0.0);
    gkyl_array_clear(g_bar, 0.0);
    gkyl_array_clear(w_bar, 0.0);
    gkyl_array_clear(h, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, basis, 0, f_bar, 0, h, 0, distg, arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, basis, 0, g_bar, 0, h, 0, distf, arr_range);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(cbasis, basis, w_bar, cfield, distf, arr_crange, arr_range);
  }
>>>>>>> main

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

    const double *wbar_d = gkyl_array_cfetch(w_bar, loc);
    int cidx[cdim];
    for (int d=0; d<cdim; d++) cidx[d] = iter.idx[d];
    long cloc = gkyl_range_idx(&arr_crange, cidx);
    const double *cf_d = gkyl_array_cfetch(cfield, cloc);
    if (poly_order==1) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[0]+cf_d[1]*f_d[1]), wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[1]+cf_d[1]*f_d[0]), wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[1]*f_d[3]+cf_d[0]*f_d[2]), wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[3]+cf_d[1]*f_d[2]), wbar_d[3], 1e-12) );
    } else if (poly_order==2) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[0]+cf_d[1]*f_d[1]+cf_d[2]*f_d[4]), wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*(cf_d[1]*f_d[4]+cf_d[2]*f_d[1])
                              +0.7071067811865475*(cf_d[0]*f_d[1]+cf_d[1]*f_d[0]), wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[0]*f_d[2]+cf_d[1]*f_d[3]+cf_d[2]*f_d[6]), wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.632455532033676*(cf_d[1]*f_d[6]+cf_d[2]*f_d[3])
                              +0.7071067811865475*(cf_d[0]*f_d[3]+cf_d[1]*f_d[2]), wbar_d[3], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[4]
                              +0.7071067811865475*(cf_d[0]*f_d[4]+cf_d[2]*f_d[0])
                              +0.6324555320336759*cf_d[1]*f_d[1], wbar_d[4], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*(cf_d[1]*f_d[7]+cf_d[0]*f_d[5]), wbar_d[5], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[6]
                              +0.7071067811865475*(cf_d[0]*f_d[6]+cf_d[2]*f_d[2])
                              +0.632455532033676*cf_d[1]*f_d[3], wbar_d[6], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[7]
                              +0.7071067811865475*(cf_d[0]*f_d[7]+cf_d[1]*f_d[5]), wbar_d[7], 1e-12) );
    } else if (poly_order==3) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[3]*f_d[8]+0.7071067811865475*cf_d[2]*f_d[4]+0.7071067811865475*cf_d[1]*f_d[1]+0.7071067811865475*cf_d[0]*f_d[0], wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081186*cf_d[2]*f_d[8]+0.6210590034081186*cf_d[3]*f_d[4]+0.6324555320336759*cf_d[1]*f_d[4]+0.6324555320336759*f_d[1]*cf_d[2]+0.7071067811865475*cf_d[0]*f_d[1]+0.7071067811865475*f_d[0]*cf_d[1], wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[3]*f_d[10]+0.7071067811865475*cf_d[2]*f_d[6]+0.7071067811865475*cf_d[1]*f_d[3]+0.7071067811865475*cf_d[0]*f_d[2], wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081187*cf_d[2]*f_d[10]+0.6210590034081187*cf_d[3]*f_d[6]+0.632455532033676*cf_d[1]*f_d[6]+0.6324555320336759*cf_d[2]*f_d[3]+0.7071067811865475*cf_d[0]*f_d[3]+0.7071067811865475*cf_d[1]*f_d[2], wbar_d[3], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[3]*f_d[8]+0.6210590034081186*cf_d[1]*f_d[8]+0.4517539514526256*cf_d[2]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[4]+0.6210590034081186*f_d[1]*cf_d[3]+0.7071067811865475*f_d[0]*cf_d[2]+0.6324555320336759*cf_d[1]*f_d[1], wbar_d[4], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[5], wbar_d[5], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4216370213557839*cf_d[3]*f_d[10]+0.6210590034081187*cf_d[1]*f_d[10]+0.4517539514526256*cf_d[2]*f_d[6]+0.7071067811865475*cf_d[0]*f_d[6]+0.6210590034081187*cf_d[3]*f_d[3]+0.632455532033676*cf_d[1]*f_d[3]+0.7071067811865475*cf_d[2]*f_d[2], wbar_d[6], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[7]+0.7071067811865475*cf_d[1]*f_d[5], wbar_d[7], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[8]+0.7071067811865475*cf_d[0]*f_d[8]+0.421637021355784*cf_d[3]*f_d[4]+0.6210590034081186*cf_d[1]*f_d[4]+0.7071067811865475*f_d[0]*cf_d[3]+0.6210590034081186*f_d[1]*cf_d[2], wbar_d[8], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[11]+0.7071067811865475*cf_d[0]*f_d[9], wbar_d[9], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[10]+0.4216370213557839*cf_d[3]*f_d[6]+0.6210590034081187*cf_d[1]*f_d[6]+0.6210590034081187*cf_d[2]*f_d[3]+0.7071067811865474*f_d[2]*cf_d[3], wbar_d[10], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[0]*f_d[11]+0.7071067811865474*cf_d[1]*f_d[9], wbar_d[11], 1e-12) );
    }
  }

  // mean ops
  double al2[2];
  if (use_gpu) {
    struct gkyl_array *mvals_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals_cu, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(basis, 0, mvals_cu, 0, distf_cu, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(basis, 1, mvals_cu, 0, distf_cu, arr_range);

    double* al2_cu = (double*) gkyl_cu_malloc(sizeof(double[2]));
    gkyl_array_reduce_range(al2_cu, mvals_cu, GKYL_SUM, arr_range);

    gkyl_cu_memcpy(al2, al2_cu, sizeof(double[2]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(al2_cu);
    gkyl_array_release(mvals_cu);
  } else {
    struct gkyl_array *mvals = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(basis, 0, mvals, 0, distf, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(basis, 1, mvals, 0, distf, arr_range);

    gkyl_array_reduce_range(al2, mvals, GKYL_SUM, arr_range);
    gkyl_array_release(mvals);
  }

  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 3.0, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 55.0/6.0, 1e-14) );  
  
  gkyl_proj_on_basis_release(proj_cfield);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(cfield);
  gkyl_array_release(w_bar);
  gkyl_dg_bin_op_mem_release(mem);  
  if (use_gpu) {
    gkyl_array_release(distf_cu);
    gkyl_array_release(distg_cu);
    gkyl_array_release(cfield_cu);
    gkyl_array_release(f_bar_cu);
    gkyl_array_release(g_bar_cu);
    gkyl_array_release(w_bar_cu);
    gkyl_array_release(h_cu);
  } else {
    gkyl_array_release(h);
  }
}

void f_1d3d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 5 + x;
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
test_3d(int poly_order, bool use_gpu)
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
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_3d, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, g_3d, NULL);

  // create array range: no ghost-cells in velocity space
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *distf_cu, *distg_cu;
  if (use_gpu) {
    distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    distg_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  }

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);
  gkyl_proj_on_basis_advance(projDistg, 0.0, &arr_range, distg);

  if (use_gpu) {
    // copy host array to device
    gkyl_array_copy(distf_cu, distf);
    gkyl_array_copy(distg_cu, distg);
  }

  // create conf-space grid:
  double clower[] = {lower[0]}, cupper[] = {upper[0]};
  int ccells[] = {cells[0]};
  int cdim = sizeof(clower)/sizeof(clower[0]);
  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, cdim, clower, cupper, ccells);
  struct gkyl_range arr_crange, arr_ext_crange;
  gkyl_create_grid_ranges(&cgrid, nghost, &arr_ext_crange, &arr_crange);
  // conf-space basis functions
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  // create a conf-space factor
  struct gkyl_array *cfield = gkyl_array_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
  // project conf-space function onto basis.
  gkyl_proj_on_basis *proj_cfield = gkyl_proj_on_basis_new(&cgrid, &cbasis, poly_order+1, 1, f_1d3d, NULL);
  gkyl_proj_on_basis_advance(proj_cfield, 0.0, &arr_crange, cfield);
  struct gkyl_array *cfield_cu;
  if (use_gpu) {
    cfield_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, cbasis.num_basis, arr_crange.volume);
    gkyl_array_copy(cfield_cu, cfield);
  }

  struct gkyl_array *f_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *g_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *w_bar = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  gkyl_array_clear(w_bar, 0.0);

  struct gkyl_array *h, *f_bar_cu, *g_bar_cu, *w_bar_cu, *h_cu;
  if (use_gpu) {
    f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    w_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(w_bar_cu, 0.0);

    // Product array only needs to be initialized on GPU
    h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h_cu, 0.0);
  } else {
    h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h, 0.0);
  }

  gkyl_dg_bin_op_mem *mem;
  if (use_gpu) {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
  } else {
    // allocate memory
    mem = gkyl_dg_bin_op_mem_new(f_bar->size, basis.num_basis);

    // h = f*g
    gkyl_dg_mul_op(basis, 0, h, 0, distf, 0, distg);
    // f_bar = h/g = f
    gkyl_dg_div_op(mem, basis, 0, f_bar, 0, h, 0, distg);
    // g_bar = h/f = g
    gkyl_dg_div_op(mem, basis, 0, g_bar, 0, h, 0, distf);
  }

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
<<<<<<< HEAD
  gkyl_array_clear(h, 0.0);
  gkyl_array_clear(f_bar, 0.0);
  gkyl_array_clear(g_bar, 0.0);
  // h = f*g
  gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, basis, 0, f_bar, 0, h, 0, distg, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, basis, 0, g_bar, 0, h, 0, distf, arr_range);
=======
  if (use_gpu) {
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(w_bar_cu, 0.0);
    gkyl_array_clear(h_cu, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(cbasis, basis, w_bar_cu, cfield_cu, distf_cu, arr_crange, arr_range);

    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
    gkyl_array_copy(w_bar, w_bar_cu);
  } else {
    gkyl_array_clear(f_bar, 0.0);
    gkyl_array_clear(g_bar, 0.0);
    gkyl_array_clear(w_bar, 0.0);
    gkyl_array_clear(h, 0.0);
    // h = f*g
    gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, arr_range);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, basis, 0, f_bar, 0, h, 0, distg, arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, basis, 0, g_bar, 0, h, 0, distf, arr_range);
    // w = cfield*f
    gkyl_dg_mul_conf_phase_op_range(cbasis, basis, w_bar, cfield, distf, arr_crange, arr_range);
  }
>>>>>>> main

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

    const double *wbar_d = gkyl_array_cfetch(w_bar, loc);
    int cidx[cdim];
    for (int d=0; d<cdim; d++) cidx[d] = iter.idx[d];
    long cloc = gkyl_range_idx(&arr_crange, cidx);
    const double *cf_d = gkyl_array_cfetch(cfield, cloc);
    if (poly_order==1) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[1]+0.7071067811865475*cf_d[0]*f_d[0], wbar_d[0], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[0]*f_d[1]+0.7071067811865475*f_d[0]*cf_d[1], wbar_d[1], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[2], wbar_d[2], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[3], wbar_d[3], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[0]*f_d[4]+0.7071067811865475*cf_d[1]*f_d[2], wbar_d[4], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[0]*f_d[5]+0.7071067811865475*cf_d[1]*f_d[3], wbar_d[5], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[6], wbar_d[6], 1e-8) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[0]*f_d[7]+0.7071067811865475*cf_d[1]*f_d[6], wbar_d[7], 1e-8) );
    } else if (poly_order==2) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[1]*f_d[1]+0.7071067811865475*cf_d[0]*f_d[0], wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[1]*f_d[7]+0.6324555320336759*f_d[1]*cf_d[2]+0.7071067811865475*cf_d[0]*f_d[1]+0.7071067811865475*f_d[0]*cf_d[1], wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[2], wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[2]*f_d[13]+0.7071067811865475*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[3], wbar_d[3], 1e-12) );
      TEST_CHECK( gkyl_compare(0.632455532033676*cf_d[1]*f_d[11]+0.6324555320336759*cf_d[2]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[4]+0.7071067811865475*cf_d[1]*f_d[2], wbar_d[4], 1e-12) );
      TEST_CHECK( gkyl_compare(0.632455532033676*cf_d[1]*f_d[13]+0.6324555320336759*cf_d[2]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[5]+0.7071067811865475*cf_d[1]*f_d[3], wbar_d[5], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[2]*f_d[17]+0.7071067811865475*cf_d[1]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[6], wbar_d[6], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[7]+0.7071067811865475*f_d[0]*cf_d[2]+0.6324555320336759*cf_d[1]*f_d[1], wbar_d[7], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[12]+0.7071067811865475*cf_d[0]*f_d[8], wbar_d[8], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[15]+0.7071067811865475*cf_d[0]*f_d[9], wbar_d[9], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[1]*f_d[17]+0.6324555320336759*cf_d[2]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[10]+0.7071067811865475*cf_d[1]*f_d[6], wbar_d[10], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[0]*f_d[11]+0.632455532033676*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[2]*f_d[2], wbar_d[11], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[12]+0.7071067811865475*cf_d[0]*f_d[12]+0.7071067811865475*cf_d[1]*f_d[8], wbar_d[12], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[13]+0.7071067811865475*cf_d[0]*f_d[13]+0.632455532033676*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[2]*f_d[3], wbar_d[13], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[18]+0.7071067811865475*cf_d[0]*f_d[14], wbar_d[14], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[15]+0.7071067811865475*cf_d[0]*f_d[15]+0.7071067811865475*cf_d[1]*f_d[9], wbar_d[15], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[19]+0.7071067811865475*cf_d[0]*f_d[16], wbar_d[16], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4517539514526256*cf_d[2]*f_d[17]+0.7071067811865475*cf_d[0]*f_d[17]+0.6324555320336759*cf_d[1]*f_d[10]+0.7071067811865475*cf_d[2]*f_d[6], wbar_d[17], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[18]+0.7071067811865475*cf_d[0]*f_d[18]+0.7071067811865475*cf_d[1]*f_d[14], wbar_d[18], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[19]+0.7071067811865475*cf_d[0]*f_d[19]+0.7071067811865475*cf_d[1]*f_d[16], wbar_d[19], 1e-12) );
    } else if (poly_order==3) {
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[3]*f_d[17]+0.7071067811865475*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[1]*f_d[1]+0.7071067811865475*cf_d[0]*f_d[0], wbar_d[0], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081186*cf_d[2]*f_d[17]+0.6210590034081186*cf_d[3]*f_d[7]+0.6324555320336759*cf_d[1]*f_d[7]+0.6324555320336759*f_d[1]*cf_d[2]+0.7071067811865475*cf_d[0]*f_d[1]+0.7071067811865475*f_d[0]*cf_d[1], wbar_d[1], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[3]*f_d[23]+0.7071067811865475*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[2], wbar_d[2], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[3]*f_d[25]+0.7071067811865475*cf_d[2]*f_d[13]+0.7071067811865475*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[3], wbar_d[3], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081187*cf_d[2]*f_d[23]+0.6210590034081187*cf_d[3]*f_d[11]+0.632455532033676*cf_d[1]*f_d[11]+0.6324555320336759*cf_d[2]*f_d[4]+0.7071067811865475*cf_d[0]*f_d[4]+0.7071067811865475*cf_d[1]*f_d[2], wbar_d[4], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081187*cf_d[2]*f_d[25]+0.6210590034081187*cf_d[3]*f_d[13]+0.632455532033676*cf_d[1]*f_d[13]+0.6324555320336759*cf_d[2]*f_d[5]+0.7071067811865475*cf_d[0]*f_d[5]+0.7071067811865475*cf_d[1]*f_d[3], wbar_d[5], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[3]*f_d[29]+0.7071067811865475*cf_d[2]*f_d[20]+0.7071067811865475*cf_d[1]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[6], wbar_d[6], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[3]*f_d[17]+0.6210590034081186*cf_d[1]*f_d[17]+0.4517539514526256*cf_d[2]*f_d[7]+0.7071067811865475*cf_d[0]*f_d[7]+0.6210590034081186*f_d[1]*cf_d[3]+0.7071067811865475*f_d[0]*cf_d[2]+0.6324555320336759*cf_d[1]*f_d[1], wbar_d[7], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[12]+0.7071067811865475*cf_d[0]*f_d[8], wbar_d[8], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[15]+0.7071067811865475*cf_d[0]*f_d[9], wbar_d[9], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6210590034081186*cf_d[2]*f_d[29]+0.6210590034081186*cf_d[3]*f_d[20]+0.6324555320336759*cf_d[1]*f_d[20]+0.6324555320336759*cf_d[2]*f_d[10]+0.7071067811865475*cf_d[0]*f_d[10]+0.7071067811865475*cf_d[1]*f_d[6], wbar_d[10], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4216370213557839*cf_d[3]*f_d[23]+0.6210590034081187*cf_d[1]*f_d[23]+0.4517539514526256*cf_d[2]*f_d[11]+0.7071067811865475*cf_d[0]*f_d[11]+0.6210590034081187*cf_d[3]*f_d[4]+0.632455532033676*cf_d[1]*f_d[4]+0.7071067811865475*cf_d[2]*f_d[2], wbar_d[11], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[12]+0.7071067811865475*cf_d[0]*f_d[12]+0.7071067811865475*cf_d[1]*f_d[8], wbar_d[12], 1e-12) );
      TEST_CHECK( gkyl_compare(0.4216370213557839*cf_d[3]*f_d[25]+0.6210590034081187*cf_d[1]*f_d[25]+0.4517539514526256*cf_d[2]*f_d[13]+0.7071067811865475*cf_d[0]*f_d[13]+0.6210590034081187*cf_d[3]*f_d[5]+0.632455532033676*cf_d[1]*f_d[5]+0.7071067811865475*cf_d[2]*f_d[3], wbar_d[13], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[21]+0.7071067811865475*cf_d[0]*f_d[14], wbar_d[14], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[15]+0.7071067811865475*cf_d[0]*f_d[15]+0.7071067811865475*cf_d[1]*f_d[9], wbar_d[15], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865475*cf_d[1]*f_d[22]+0.7071067811865475*cf_d[0]*f_d[16], wbar_d[16], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[17]+0.7071067811865475*cf_d[0]*f_d[17]+0.421637021355784*cf_d[3]*f_d[7]+0.6210590034081186*cf_d[1]*f_d[7]+0.7071067811865475*f_d[0]*cf_d[3]+0.6210590034081186*f_d[1]*cf_d[2], wbar_d[17], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[24]+0.7071067811865475*cf_d[0]*f_d[18], wbar_d[18], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[27]+0.7071067811865475*cf_d[0]*f_d[19], wbar_d[19], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[3]*f_d[29]+0.6210590034081186*cf_d[1]*f_d[29]+0.4517539514526256*cf_d[2]*f_d[20]+0.7071067811865475*cf_d[0]*f_d[20]+0.6210590034081186*cf_d[3]*f_d[10]+0.6324555320336759*cf_d[1]*f_d[10]+0.7071067811865475*cf_d[2]*f_d[6], wbar_d[20], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[21]+0.7071067811865475*cf_d[0]*f_d[21]+0.7071067811865475*cf_d[1]*f_d[14], wbar_d[21], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[22]+0.7071067811865475*cf_d[0]*f_d[22]+0.7071067811865475*cf_d[1]*f_d[16], wbar_d[22], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[23]+0.7071067811865475*cf_d[0]*f_d[23]+0.4216370213557839*cf_d[3]*f_d[11]+0.6210590034081187*cf_d[1]*f_d[11]+0.6210590034081187*cf_d[2]*f_d[4]+0.7071067811865474*f_d[2]*cf_d[3], wbar_d[23], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[24]+0.7071067811865475*cf_d[0]*f_d[24]+0.7071067811865474*cf_d[1]*f_d[18], wbar_d[24], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[25]+0.7071067811865475*cf_d[0]*f_d[25]+0.4216370213557839*cf_d[3]*f_d[13]+0.6210590034081187*cf_d[1]*f_d[13]+0.6210590034081187*cf_d[2]*f_d[5]+0.7071067811865474*cf_d[3]*f_d[3], wbar_d[25], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[30]+0.7071067811865475*cf_d[0]*f_d[26], wbar_d[26], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[27]+0.7071067811865475*cf_d[0]*f_d[27]+0.7071067811865474*cf_d[1]*f_d[19], wbar_d[27], 1e-12) );
      TEST_CHECK( gkyl_compare(0.7071067811865474*cf_d[1]*f_d[31]+0.7071067811865475*cf_d[0]*f_d[28], wbar_d[28], 1e-12) );
      TEST_CHECK( gkyl_compare(0.421637021355784*cf_d[2]*f_d[29]+0.7071067811865475*cf_d[0]*f_d[29]+0.421637021355784*cf_d[3]*f_d[20]+0.6210590034081186*cf_d[1]*f_d[20]+0.6210590034081186*cf_d[2]*f_d[10]+0.7071067811865475*cf_d[3]*f_d[6], wbar_d[29], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[30]+0.7071067811865475*cf_d[0]*f_d[30]+0.7071067811865474*cf_d[1]*f_d[26], wbar_d[30], 1e-12) );
      TEST_CHECK( gkyl_compare(0.6324555320336759*cf_d[2]*f_d[31]+0.7071067811865475*cf_d[0]*f_d[31]+0.7071067811865474*cf_d[1]*f_d[28], wbar_d[31], 1e-12) );
    }
  }

  // mean ops
  double al2[2];
  if (use_gpu) {
    struct gkyl_array *mvals_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals_cu, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(basis, 0, mvals_cu, 0, distf_cu, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(basis, 1, mvals_cu, 0, distf_cu, arr_range);

    double* al2_cu = (double*) gkyl_cu_malloc(sizeof(double[2]));
    gkyl_array_reduce_range(al2_cu, mvals_cu, GKYL_SUM, arr_range);

    gkyl_cu_memcpy(al2, al2_cu, sizeof(double[2]), GKYL_CU_MEMCPY_D2H);
    gkyl_cu_free(al2_cu);
    gkyl_array_release(mvals_cu);
  } else {
    struct gkyl_array *mvals = gkyl_array_new(GKYL_DOUBLE, 2, arr_range.volume);
    gkyl_array_clear(mvals, 0.0);

    // means are stored in h[0]
    gkyl_dg_calc_average_range(basis, 0, mvals, 0, distf, arr_range);
    // L2 are stored in h[1]
    gkyl_dg_calc_l2_range(basis, 1, mvals, 0, distf, arr_range);

    gkyl_array_reduce_range(al2, mvals, GKYL_SUM, arr_range);
    gkyl_array_release(mvals);
  }

  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 6.5, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 85.0/2.0, 1e-14) );

  gkyl_proj_on_basis_release(proj_cfield);
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_array_release(distf);
  gkyl_array_release(distg);
  gkyl_array_release(f_bar);
  gkyl_array_release(g_bar);
  gkyl_array_release(cfield);
  gkyl_array_release(w_bar);
  gkyl_dg_bin_op_mem_release(mem);
  if (use_gpu) {
    gkyl_array_release(distf_cu);
    gkyl_array_release(distg_cu);
    gkyl_array_release(cfield_cu);
    gkyl_array_release(f_bar_cu);
    gkyl_array_release(g_bar_cu);
    gkyl_array_release(w_bar_cu);
    gkyl_array_release(h_cu);
  } else {
    gkyl_array_release(h);
  }
}

void test_1d_p1(){ test_1d(1, false); }
void test_1d_p2(){ test_1d(2, false); }
void test_1d_p3(){ test_1d(3, false); }

void test_2d_p1(){ test_2d(1, false); }
void test_2d_p2(){ test_2d(2, false); }
void test_2d_p3(){ test_2d(3, false); }

void test_3d_p1(){ test_3d(1, false); }
void test_3d_p2(){ test_3d(2, false); }

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

  // allocate memory
  gkyl_dg_bin_op_mem *mem = gkyl_dg_bin_op_mem_new(f_bar->size, basis.num_basis);  

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h, 0, distf, 0, distg);
  // f_bar = h/g = f
  gkyl_dg_div_op(mem, basis, 0, f_bar, 0, h, 0, distg);
  // g_bar = h/f = g
  gkyl_dg_div_op(mem, basis, 0, g_bar, 0, h, 0, distf);

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
  gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, basis, 0, f_bar, 0, h, 0, distg, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, basis, 0, g_bar, 0, h, 0, distf, arr_range);

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
  gkyl_dg_bin_op_mem_release(mem);  
}

// Cuda specific tests
#ifdef GKYL_HAVE_CUDA

<<<<<<< HEAD
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

  // allocate memory
  gkyl_dg_bin_op_mem *mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);
  
  // h = f*g
  gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
  // f_bar = h/g = f
  gkyl_dg_div_op(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
  // g_bar = h/f = g
  gkyl_dg_div_op(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

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
  gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);

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

  // mean ops
  struct gkyl_array *mvals_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, arr_range.volume);
  gkyl_array_clear(mvals_cu, 0.0);

  // means are stored in h[0]
  gkyl_dg_calc_average_range(basis, 0, mvals_cu, 0, distf_cu, arr_range);
  // L2 are stored in h[1]
  gkyl_dg_calc_l2_range(basis, 1, mvals_cu, 0, distf_cu, arr_range);
=======
void test_1d_p1_cu(){ test_1d(1, true); }
void test_1d_p2_cu(){ test_1d(2, true); }
void test_1d_p3_cu(){ test_1d(3, true); }
>>>>>>> main

void test_2d_p1_cu(){ test_2d(1, true); }
void test_2d_p2_cu(){ test_2d(2, true); }
void test_2d_p3_cu(){ test_2d(3, true); }

<<<<<<< HEAD
  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 2.5, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 19.0/3.0, 1e-14) );
  
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
  gkyl_array_release(mvals_cu);
  gkyl_cu_free(al2_cu);
  gkyl_dg_bin_op_mem_release(mem);  
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

  gkyl_dg_bin_op_mem *mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);  

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
  // f_bar = h/g = f
  gkyl_dg_div_op(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
  // g_bar = h/f = g
  gkyl_dg_div_op(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

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
  gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);

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

  // mean ops
  struct gkyl_array *mvals_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, arr_range.volume);
  gkyl_array_clear(mvals_cu, 0.0);

  // means are stored in h[0]
  gkyl_dg_calc_average_range(basis, 0, mvals_cu, 0, distf_cu, arr_range);
  // L2 are stored in h[1]
  gkyl_dg_calc_l2_range(basis, 1, mvals_cu, 0, distf_cu, arr_range);

  double* al2_cu = (double*) gkyl_cu_malloc(sizeof(double[2]));
  gkyl_array_reduce_range(al2_cu, mvals_cu, GKYL_SUM, arr_range);

  double al2[2];
  gkyl_cu_memcpy(al2, al2_cu, sizeof(double[2]), GKYL_CU_MEMCPY_D2H);

  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 3.0, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 55.0/6.0, 1e-14) );
  
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
  gkyl_array_release(mvals_cu);
  gkyl_cu_free(al2_cu);
  gkyl_dg_bin_op_mem_release(mem);
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

  gkyl_dg_bin_op_mem *mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
  // f_bar = h/g = f
  gkyl_dg_div_op(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
  // g_bar = h/f = g
  gkyl_dg_div_op(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

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
  gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);

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

  // mean ops
  struct gkyl_array *mvals_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 2, arr_range.volume);
  gkyl_array_clear(mvals_cu, 0.0);

  // means are stored in h[0]
  gkyl_dg_calc_average_range(basis, 0, mvals_cu, 0, distf_cu, arr_range);
  // L2 are stored in h[1]
  gkyl_dg_calc_l2_range(basis, 1, mvals_cu, 0, distf_cu, arr_range);

  double* al2_cu = (double*) gkyl_cu_malloc(sizeof(double[2]));
  gkyl_array_reduce_range(al2_cu, mvals_cu, GKYL_SUM, arr_range);

  double al2[2];
  gkyl_cu_memcpy(al2, al2_cu, sizeof(double[2]), GKYL_CU_MEMCPY_D2H);

  double vol = grid.cellVolume;
  TEST_CHECK( gkyl_compare(al2[0]*vol, 6.5, 1e-14) );
  TEST_CHECK( gkyl_compare(al2[1]*vol, 85.0/2.0, 1e-14) );
  
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
  gkyl_array_release(mvals_cu);
  gkyl_cu_free(al2_cu);
  gkyl_dg_bin_op_mem_release(mem);
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
=======
void test_3d_p1_cu(){ test_3d(1, true); }
void test_3d_p2_cu(){ test_3d(2, true); }
>>>>>>> main

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

  gkyl_dg_bin_op_mem *mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);

  // h = f*g
  gkyl_dg_mul_op(basis, 0, h_cu, 0, distf_cu, 0, distg_cu);
  // f_bar = h/g = f
  gkyl_dg_div_op(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu);
  // g_bar = h/f = g
  gkyl_dg_div_op(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu);

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
  gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, &arr_range);
  // f_bar = h/g = f
  gkyl_dg_div_op_range(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, arr_range);
  // g_bar = h/f = g
  gkyl_dg_div_op_range(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, arr_range);

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
  gkyl_dg_bin_op_mem_release(mem);
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
