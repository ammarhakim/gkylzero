// Test integration of a gkyl_array over a range.
//

#include <acutest.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_average.h>
#include <math.h>
#include <assert.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size, bool use_gpu)
{
  struct gkyl_array *a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void evalFunc_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {-6.0}, upper[] = {6.0}; // Has to match the test below.
  fout[0] = 1./(upper[0]-lower[0]);
}

void test_1x(int poly_order, bool use_gpu)
{
  double lower[] = {-6.0}, upper[] = {6.0};
  int cells[] = {16};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int nc = 1;

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // local, local-ext ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // projection updater for dist-function
  gkyl_proj_on_basis *projf;
  projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_1x, NULL);

  // create distribution function array
  struct gkyl_array *distf = mkarr(nc*basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *distf_ho = use_gpu? mkarr(nc*basis.num_basis, local_ext.volume, false) : distf;

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projf, 0.0, &local, distf_ho);
  if (use_gpu) gkyl_array_copy(distf, distf_ho);

  // integrate distribution function.
  struct gkyl_array_integrate *integ_up = gkyl_array_integrate_new(&grid, &basis, nc, GKYL_ARRAY_INTEGRATE_OP_NONE, use_gpu);

  double *fint = use_gpu? gkyl_cu_malloc(nc*sizeof(double)) : gkyl_malloc(nc*sizeof(double));
  struct gkyl_array *weight = mkarr(nc*basis.num_basis, local_ext.volume, use_gpu);
  gkyl_array_integrate_advance(integ_up, distf, 1., weight, &local, fint);

  gkyl_array_integrate_release(integ_up);

  double *fint_ho = gkyl_malloc(nc*sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(fint_ho, fint, nc*sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(fint_ho, fint, nc*sizeof(double));

  printf("\t declare the scalar structures to perform the x average\n");
  // integrate the 1D array
  // scalar gkyl_array (1D)
  struct gkyl_array *avg_res = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *avg_res_ho = use_gpu? mkarr(avg_res->ncomp, avg_res->size, false)
                                        : gkyl_array_acquire(avg_res);

  printf("\t create 1D integral updater and advance it\n");
  // Create an array average updater
  struct gkyl_array_average *avg_x;
  gkyl_array_average_new(&grid, &basis, GKYL_ARRAY_AVERAGE_OP, use_gpu);
  gkyl_array_average_advance(avg_x, &local, &local, distf, avg_res);
  gkyl_array_average_release(avg_x);

  // check output
  struct gkyl_range_iter iter_;
  gkyl_range_iter_init(&iter_, &local);
  while (gkyl_range_iter_next(&iter_)) {
    long lidx = gkyl_range_idx(&local, iter_.idx);
    const double *m_i = gkyl_array_cfetch(avg_res, lidx);
    printf("\t\tm_[%ld]=%g\n",lidx,m_i[0]);
  }

  TEST_CHECK( gkyl_compare( 1.0, fint_ho[0], 1e-12) );

  gkyl_array_release(avg_res);
  gkyl_array_release(avg_res_ho);
  gkyl_array_release(distf);
  gkyl_array_release(weight);
  gkyl_proj_on_basis_release(projf);
  gkyl_free(fint_ho);
  if (use_gpu)
    gkyl_cu_free(fint);
  else
    gkyl_free(fint);
}

void evalFunc_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {0., -6.0}, upper[] = {2., 6.0}; // Has to match the test below.
  fout[0] = 1./((upper[0]-lower[0])*(upper[1]-lower[1]));
}

void test_2x(int poly_order, bool use_gpu)
{
  double lower[] = {0., -6.0}, upper[] = {2., 6.0};
  int cells[] = {6, 16};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int nc = 1;

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1, 0 };
  struct gkyl_range local, local_ext; // local, local-ext ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // projection updater for dist-function
  gkyl_proj_on_basis *projf;
  projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_2x, NULL);

  // create distribution function array
  struct gkyl_array *distf = mkarr(nc*basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *distf_ho = use_gpu? mkarr(nc*basis.num_basis, local_ext.volume, false) : distf;

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projf, 0.0, &local, distf_ho);
  if (use_gpu) gkyl_array_copy(distf, distf_ho);

  // integrate distribution function.
  struct gkyl_array_integrate *integ_up = gkyl_array_integrate_new(&grid, &basis, nc, GKYL_ARRAY_INTEGRATE_OP_NONE, use_gpu);

  double *fint = use_gpu? gkyl_cu_malloc(nc*sizeof(double)) : gkyl_malloc(nc*sizeof(double));
  struct gkyl_array *weight = mkarr(nc*basis.num_basis, local_ext.volume, use_gpu);
  gkyl_array_integrate_advance(integ_up, distf, 1., weight, &local, fint);

  gkyl_array_integrate_release(integ_up);

  double *fint_ho = gkyl_malloc(nc*sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(fint_ho, fint, nc*sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(fint_ho, fint, nc*sizeof(double));

  TEST_CHECK( gkyl_compare( 1.0, fint_ho[0], 1e-12) );

  gkyl_array_release(distf);
  gkyl_array_release(weight);
  if (use_gpu) gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(projf);
  gkyl_free(fint_ho);
  if (use_gpu)
    gkyl_cu_free(fint);
  else
    gkyl_free(fint);
}

void test_1x_cpu()
{
  // p=1
  test_1x(1, false);

  // p=2
  // test_1x(2, false);

}

void test_2x_cpu()
{
  // p=1
  test_2x(1, false);

  // p=2
  // test_2x(2, false);
}

#ifdef GKYL_HAVE_CUDA
void test_1x_gpu()
{
  // p=1
  test_1x(1, true);

  // p=2
  // test_1x(2, true);
}

void test_2x_gpu()
{
  // p=1
  test_2x_nc1_op(1, true);
  
  // p=2
  // test_2x_nc1_op(2, true);
}
#endif

TEST_LIST = {
  { "test_1x_cpu", test_1x_cpu },
  { "test_2x_cpu", test_2x_cpu },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_gpu", test_1x_gpu },
  { "test_2x_gpu", test_2x_gpu },
#endif
  { NULL, NULL },
};
