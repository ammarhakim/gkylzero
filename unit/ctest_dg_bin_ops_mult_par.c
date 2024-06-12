#include <acutest.h>
#include <time.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

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

void fv2_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2. + x;
  fout[1] = 1. - x;
}
void gv2_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2*x*x + 8;
  fout[1] = x*x - 8;
}

void fv3_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2. + x;
  fout[1] = 1. - x;
  fout[2] = 0.5 + x;
}
void gv3_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 2*x*x + 8;
  fout[1] = x*x - 8;
  fout[2] = 0.5*x*x + 4;
}

void check_dot_product_1d(const double *fv1_d, const double *gv1_d, const double *fvdgv1_d,
                          const double *fv2_d, const double *gv2_d, const double *fvdgv2_d,
                          const double *fv3_d, const double *gv3_d, const double *fvdgv3_d,
                          int poly_order) {
  if (poly_order == 1) {
    TEST_CHECK( gkyl_compare(fvdgv1_d[0], 0.7071067811865475*(fv1_d[1]*gv1_d[1]+fv1_d[0]*gv1_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv1_d[1], 0.7071067811865475*(fv1_d[0]*gv1_d[1]+gv1_d[0]*fv1_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[0], 0.7071067811865475*(fv2_d[3]*gv2_d[3]+fv2_d[2]*gv2_d[2]+fv2_d[1]*gv2_d[1]+fv2_d[0]*gv2_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[1], 0.7071067811865475*(fv2_d[2]*gv2_d[3]+gv2_d[2]*fv2_d[3]+fv2_d[0]*gv2_d[1]+gv2_d[0]*fv2_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[0], 0.7071067811865475*(fv3_d[5]*gv3_d[5]+fv3_d[4]*gv3_d[4]+fv3_d[3]*gv3_d[3]+fv3_d[2]*gv3_d[2]+fv3_d[1]*gv3_d[1]+fv3_d[0]*gv3_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[1], 0.7071067811865475*(fv3_d[4]*gv3_d[5]+gv3_d[4]*fv3_d[5]+fv3_d[2]*gv3_d[3]+gv3_d[2]*fv3_d[3]+fv3_d[0]*gv3_d[1]+gv3_d[0]*fv3_d[1]), 1e-12) );
  } else if (poly_order == 2) {
    TEST_CHECK( gkyl_compare(fvdgv1_d[0], 0.7071067811865475*(fv1_d[2]*gv1_d[2]+fv1_d[1]*gv1_d[1]+fv1_d[0]*gv1_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv1_d[1], 0.1414213562373095*(4.47213595499958*fv1_d[1]*gv1_d[2]+4.47213595499958*gv1_d[1]*fv1_d[2]+5.0*fv1_d[0]*gv1_d[1]+5.0*gv1_d[0]*fv1_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv1_d[2], 0.02020305089104421*((22.3606797749979*fv1_d[2]+35.0*fv1_d[0])*gv1_d[2]+35.0*gv1_d[0]*fv1_d[2]+31.30495168499706*fv1_d[1]*gv1_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[0], 0.7071067811865475*(fv2_d[5]*gv2_d[5]+fv2_d[4]*gv2_d[4]+fv2_d[3]*gv2_d[3]+fv2_d[2]*gv2_d[2]+fv2_d[1]*gv2_d[1]+fv2_d[0]*gv2_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[1], 0.1414213562373095*(4.47213595499958*fv2_d[4]*gv2_d[5]+4.47213595499958*gv2_d[4]*fv2_d[5]+5.0*fv2_d[3]*gv2_d[4]+5.0*gv2_d[3]*fv2_d[4]+4.47213595499958*fv2_d[1]*gv2_d[2]+4.47213595499958*gv2_d[1]*fv2_d[2]+5.0*fv2_d[0]*gv2_d[1]+5.0*gv2_d[0]*fv2_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv2_d[2], 0.02020305089104421*((22.3606797749979*fv2_d[5]+35.0*fv2_d[3])*gv2_d[5]+35.0*gv2_d[3]*fv2_d[5]+31.30495168499706*fv2_d[4]*gv2_d[4]+(22.3606797749979*fv2_d[2]+35.0*fv2_d[0])*gv2_d[2]+35.0*gv2_d[0]*fv2_d[2]+31.30495168499706*fv2_d[1]*gv2_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[0], 0.7071067811865475*(fv3_d[8]*gv3_d[8]+fv3_d[7]*gv3_d[7]+fv3_d[6]*gv3_d[6]+fv3_d[5]*gv3_d[5]+fv3_d[4]*gv3_d[4]+fv3_d[3]*gv3_d[3]+fv3_d[2]*gv3_d[2]+fv3_d[1]*gv3_d[1]+fv3_d[0]*gv3_d[0]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[1], 0.1414213562373095*(4.47213595499958*fv3_d[7]*gv3_d[8]+4.47213595499958*gv3_d[7]*fv3_d[8]+5.0*fv3_d[6]*gv3_d[7]+5.0*gv3_d[6]*fv3_d[7]+4.47213595499958*fv3_d[4]*gv3_d[5]+4.47213595499958*gv3_d[4]*fv3_d[5]+5.0*fv3_d[3]*gv3_d[4]+5.0*gv3_d[3]*fv3_d[4]+4.47213595499958*fv3_d[1]*gv3_d[2]+4.47213595499958*gv3_d[1]*fv3_d[2]+5.0*fv3_d[0]*gv3_d[1]+5.0*gv3_d[0]*fv3_d[1]), 1e-12) );
    TEST_CHECK( gkyl_compare(fvdgv3_d[2], 0.02020305089104421*((22.3606797749979*fv3_d[8]+35.0*fv3_d[6])*gv3_d[8]+35.0*gv3_d[6]*fv3_d[8]+31.30495168499706*fv3_d[7]*gv3_d[7]+(22.3606797749979*fv3_d[5]+35.0*fv3_d[3])*gv3_d[5]+35.0*gv3_d[3]*fv3_d[5]+31.30495168499706*fv3_d[4]*gv3_d[4]+(22.3606797749979*fv3_d[2]+35.0*fv3_d[0])*gv3_d[2]+35.0*gv3_d[0]*fv3_d[2]+31.30495168499706*fv3_d[1]*gv3_d[1]), 1e-12) );
  }
  return;
}

void
test_1d(int poly_order, bool use_gpu)
{
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {10000};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, f_1d, NULL);
  gkyl_proj_on_basis *projDistg = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, g_1d, NULL);
  // projection updaters for vector fields.
  gkyl_proj_on_basis *projfv2 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 2, fv2_1d, NULL);
  gkyl_proj_on_basis *projgv2 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 2, gv2_1d, NULL);
  gkyl_proj_on_basis *projfv3 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 3, fv3_1d, NULL);
  gkyl_proj_on_basis *projgv3 = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 3, gv3_1d, NULL);

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
  gkyl_dg_bin_op_mem *mem;
  if (use_gpu) {
    f_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    g_bar_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);

    // Product array only needs to be initialized on GPU.
    h_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h_cu, 0.0);
    // allocate memory
    mem = gkyl_dg_bin_op_mem_cu_dev_new(f_bar->size, basis.num_basis);
  } 
  else {
    h = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
    gkyl_array_clear(h, 0.0);
    // allocate memory
    mem = gkyl_dg_bin_op_mem_new(f_bar->size, basis.num_basis);
  }

  // Test range methods
  if (use_gpu) {
    gkyl_array_clear(f_bar_cu, 0.0);
    gkyl_array_clear(g_bar_cu, 0.0);
    gkyl_array_clear(h_cu, 0.0);
    // h = f*g
    clock_t start = clock();
    gkyl_dg_mul_op_range(basis, 0, h_cu, 0, distf_cu, 0, distg_cu, &arr_range);
    clock_t end = clock();
    double time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for multiplication on GPU: %f\n", time_taken);

    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, basis, 0, f_bar_cu, 0, h_cu, 0, distg_cu, &arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, basis, 0, g_bar_cu, 0, h_cu, 0, distf_cu, &arr_range);
    // copy from device and check if things are ok
    gkyl_array_copy(f_bar, f_bar_cu);
    gkyl_array_copy(g_bar, g_bar_cu);
  } 
  else {
    gkyl_array_clear(f_bar, 0.0);
    gkyl_array_clear(g_bar, 0.0);
    gkyl_array_clear(h, 0.0);
    // h = f*g

    clock_t start = clock();    
    gkyl_dg_mul_op_range(basis, 0, h, 0, distf, 0, distg, &arr_range);
    clock_t end = clock();
    double time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time taken for multiplication on CPU: %f\n", time_taken);
    // f_bar = h/g = f
    gkyl_dg_div_op_range(mem, basis, 0, f_bar, 0, h, 0, distg, &arr_range);
    // g_bar = h/f = g
    gkyl_dg_div_op_range(mem, basis, 0, g_bar, 0, h, 0, distf, &arr_range);
  }

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
  
  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistg);
  gkyl_proj_on_basis_release(projfv2);
  gkyl_proj_on_basis_release(projgv2);
  gkyl_proj_on_basis_release(projfv3);
  gkyl_proj_on_basis_release(projgv3);
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

void test_1d_p1(){ test_1d(1, false); }
void test_1d_p2(){ test_1d(2, false); }
void test_1d_p3(){ test_1d(3, false); }

// Cuda specific tests
#ifdef GKYL_HAVE_CUDA

void test_1d_p1_cu(){ test_1d(1, true); }
void test_1d_p2_cu(){ test_1d(2, true); }
void test_1d_p3_cu(){ test_1d(3, true); }

#endif

TEST_LIST = {
  { "test_1d_p1", test_1d_p1 },
  { "test_1d_p2", test_1d_p2 },
  { "test_1d_p3", test_1d_p3 },
#ifdef GKYL_HAVE_CUDA
  { "test_1d_p1_cu", test_1d_p1_cu },
  { "test_1d_p2_cu", test_1d_p2_cu },
  { "test_1d_p3_cu", test_1d_p3_cu },
#endif
  { NULL, NULL },
};
