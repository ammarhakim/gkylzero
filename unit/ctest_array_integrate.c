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
#include <math.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size, bool use_gpu)
{
  struct gkyl_array *a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void evalFunc_1x_nc1_op_none(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {-6.0}, upper[] = {6.0}; // Has to match the test below.
  fout[0] = 1./(upper[0]-lower[0]);
}

void evalFunc_1x_nc1_op_sq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {-6.0}, upper[] = {6.0}; // Has to match the test below.
  fout[0] = 1./sqrt(upper[0]-lower[0]);
}

void test_1x_nc1_op(enum gkyl_array_integrate_op integ_op, int poly_order, bool use_gpu)
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
  if (integ_op == GKYL_ARRAY_INTEGRATE_OP_SQ)
    projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_1x_nc1_op_sq, NULL);
  else
    projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_1x_nc1_op_none, NULL);

  // create distribution function array
  struct gkyl_array *distf = mkarr(nc*basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *distf_ho = use_gpu? mkarr(nc*basis.num_basis, local_ext.volume, false) : distf;

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projf, 0.0, &local, distf_ho);
  if (use_gpu) gkyl_array_copy(distf, distf_ho);

  if (integ_op == GKYL_ARRAY_INTEGRATE_OP_ABS)
    gkyl_array_scale(distf, -1.);

  // integrate distribution function.
  struct gkyl_array_integrate *integ_up = gkyl_array_integrate_new(&grid, &basis, nc, integ_op, use_gpu);

  double *fint = use_gpu? gkyl_cu_malloc(nc*sizeof(double)) : gkyl_malloc(nc*sizeof(double));
  gkyl_array_integrate_advance(integ_up, distf, 1., &local, fint);

  gkyl_array_integrate_release(integ_up);

  double *fint_ho = gkyl_malloc(nc*sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(fint_ho, fint, nc*sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(fint_ho, fint, nc*sizeof(double));

  TEST_CHECK( gkyl_compare( 1.0, fint_ho[0], 1e-12) );

  gkyl_array_release(distf);
  gkyl_proj_on_basis_release(projf);
  gkyl_free(fint_ho);
  if (use_gpu)
    gkyl_cu_free(fint);
  else
    gkyl_free(fint);
}

void evalFunc_1x_nc3_op_none(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {-6.0}, upper[] = {6.0}; // Has to match the test below.
  fout[0] = 1./(upper[0]-lower[0]);
  fout[1] = 1.5/(upper[0]-lower[0]);
  fout[2] = 2.5/(upper[0]-lower[0]);
}

void evalFunc_1x_nc3_op_sq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {-6.0}, upper[] = {6.0}; // Has to match the test below.
  fout[0] = 1./sqrt(upper[0]-lower[0]);
  fout[1] = 1.5/sqrt(upper[0]-lower[0]);
  fout[2] = 2.5/sqrt(upper[0]-lower[0]);
}

void test_1x_nc3_op(enum gkyl_array_integrate_op integ_op, int poly_order, bool use_gpu)
{
  double lower[] = {-6.0}, upper[] = {6.0};
  int cells[] = {16};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int nc = 3;

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
  if (integ_op == GKYL_ARRAY_INTEGRATE_OP_SQ)
    projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_1x_nc3_op_sq, NULL);
  else
    projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_1x_nc3_op_none, NULL);

  // create distribution function array
  struct gkyl_array *distf = mkarr(nc*basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *distf_ho = use_gpu? mkarr(nc*basis.num_basis, local_ext.volume, false) : distf;

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projf, 0.0, &local, distf_ho);
  if (use_gpu) gkyl_array_copy(distf, distf_ho);

  if (integ_op == GKYL_ARRAY_INTEGRATE_OP_ABS)
    gkyl_array_scale(distf, -1.);

  // integrate distribution function.
  struct gkyl_array_integrate *integ_up = gkyl_array_integrate_new(&grid, &basis, nc, integ_op, use_gpu);

  double *fint = use_gpu? gkyl_cu_malloc(nc*sizeof(double)) : gkyl_malloc(nc*sizeof(double));
  gkyl_array_integrate_advance(integ_up, distf, 1., &local, fint);

  gkyl_array_integrate_release(integ_up);

  double *fint_ho = gkyl_malloc(nc*sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(fint_ho, fint, nc*sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(fint_ho, fint, nc*sizeof(double));

  TEST_CHECK( gkyl_compare( 1.0, fint_ho[0], 1e-12) );
  TEST_CHECK( gkyl_compare( integ_op == GKYL_ARRAY_INTEGRATE_OP_SQ? 1.5*1.5 : 1.5, fint_ho[1], 1e-12) );
  TEST_CHECK( gkyl_compare( integ_op == GKYL_ARRAY_INTEGRATE_OP_SQ? 2.5*2.5 : 2.5, fint_ho[2], 1e-12) );

  gkyl_array_release(distf);
  gkyl_proj_on_basis_release(projf);
  gkyl_free(fint_ho);
  if (use_gpu)
    gkyl_cu_free(fint);
  else
    gkyl_free(fint);
}

void evalFunc_2x_nc1_op_none(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {0., -6.0}, upper[] = {2., 6.0}; // Has to match the test below.
  fout[0] = 1./((upper[0]-lower[0])*(upper[1]-lower[1]));
}

void evalFunc_2x_nc1_op_sq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {0., -6.0}, upper[] = {2., 6.0}; // Has to match the test below.
  fout[0] = 1./sqrt((upper[0]-lower[0])*(upper[1]-lower[1]));
}

void test_2x_nc1_op(enum gkyl_array_integrate_op integ_op, int poly_order, bool use_gpu)
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
  if (integ_op == GKYL_ARRAY_INTEGRATE_OP_SQ)
    projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_2x_nc1_op_sq, NULL);
  else
    projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_2x_nc1_op_none, NULL);

  // create distribution function array
  struct gkyl_array *distf = mkarr(nc*basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *distf_ho = use_gpu? mkarr(nc*basis.num_basis, local_ext.volume, false) : distf;

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projf, 0.0, &local, distf_ho);
  if (use_gpu) gkyl_array_copy(distf, distf_ho);

  if (integ_op == GKYL_ARRAY_INTEGRATE_OP_ABS)
    gkyl_array_scale(distf, -1.);

  // integrate distribution function.
  struct gkyl_array_integrate *integ_up = gkyl_array_integrate_new(&grid, &basis, nc, integ_op, use_gpu);

  double *fint = use_gpu? gkyl_cu_malloc(nc*sizeof(double)) : gkyl_malloc(nc*sizeof(double));
  gkyl_array_integrate_advance(integ_up, distf, 1., &local, fint);

  gkyl_array_integrate_release(integ_up);

  double *fint_ho = gkyl_malloc(nc*sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(fint_ho, fint, nc*sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(fint_ho, fint, nc*sizeof(double));

  TEST_CHECK( gkyl_compare( 1.0, fint_ho[0], 1e-12) );

  gkyl_array_release(distf);
  gkyl_proj_on_basis_release(projf);
  gkyl_free(fint_ho);
  if (use_gpu)
    gkyl_cu_free(fint);
  else
    gkyl_free(fint);
}

void evalFunc_2x_nc3_op_none(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {0., -6.0}, upper[] = {2., 6.0}; // Has to match the test below.
  fout[0] = 1./((upper[0]-lower[0])*(upper[1]-lower[1]));
  fout[1] = 1.5/((upper[0]-lower[0])*(upper[1]-lower[1]));
  fout[2] = 2.5/((upper[0]-lower[0])*(upper[1]-lower[1]));
}

void evalFunc_2x_nc3_op_sq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {0., -6.0}, upper[] = {2., 6.0}; // Has to match the test below.
  fout[0] = 1./sqrt((upper[0]-lower[0])*(upper[1]-lower[1]));
  fout[1] = 1.5/sqrt((upper[0]-lower[0])*(upper[1]-lower[1]));
  fout[2] = 2.5/sqrt((upper[0]-lower[0])*(upper[1]-lower[1]));
}

void test_2x_nc3_op(enum gkyl_array_integrate_op integ_op, int poly_order, bool use_gpu)
{
  double lower[] = {0., -6.0}, upper[] = {2., 6.0};
  int cells[] = {6, 16};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  int nc = 3;

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
  if (integ_op == GKYL_ARRAY_INTEGRATE_OP_SQ)
    projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_2x_nc3_op_sq, NULL);
  else
    projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_2x_nc3_op_none, NULL);

  // create distribution function array
  struct gkyl_array *distf = mkarr(nc*basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *distf_ho = use_gpu? mkarr(nc*basis.num_basis, local_ext.volume, false) : distf;

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projf, 0.0, &local, distf_ho);
  if (use_gpu) gkyl_array_copy(distf, distf_ho);

  if (integ_op == GKYL_ARRAY_INTEGRATE_OP_ABS)
    gkyl_array_scale(distf, -1.);

  // integrate distribution function.
  struct gkyl_array_integrate *integ_up = gkyl_array_integrate_new(&grid, &basis, nc, integ_op, use_gpu);

  double *fint = use_gpu? gkyl_cu_malloc(nc*sizeof(double)) : gkyl_malloc(nc*sizeof(double));
  gkyl_array_integrate_advance(integ_up, distf, 1., &local, fint);

  gkyl_array_integrate_release(integ_up);

  double *fint_ho = gkyl_malloc(nc*sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(fint_ho, fint, nc*sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(fint_ho, fint, nc*sizeof(double));

  TEST_CHECK( gkyl_compare( 1.0, fint_ho[0], 1e-12) );
  TEST_CHECK( gkyl_compare( integ_op == GKYL_ARRAY_INTEGRATE_OP_SQ? 1.5*1.5 : 1.5, fint_ho[1], 1e-12) );
  TEST_CHECK( gkyl_compare( integ_op == GKYL_ARRAY_INTEGRATE_OP_SQ? 2.5*2.5 : 2.5, fint_ho[2], 1e-12) );

  gkyl_array_release(distf);
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
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 1, false);
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 1, false);
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 1, false);

  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 1, false);
  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 1, false);
  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 1, false);

  // p=2
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 2, false);
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 2, false);
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 2, false);

  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 2, false);
  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 2, false);
  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 2, false);
}

void test_2x_cpu()
{
  // p=1
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 1, false);
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 1, false);
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 1, false);

  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 1, false);
  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 1, false);
  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 1, false);

  // p=2
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 2, false);
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 2, false);
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 2, false);

  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 2, false);
  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 2, false);
  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 2, false);
}

#ifdef GKYL_HAVE_CUDA
void test_1x_gpu()
{
  // p=1
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 1, true);
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 1, true);
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 1, true);

  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 1, true);
  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 1, true);
  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 1, true);

  // p=2
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 2, true);
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 2, true);
  test_1x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 2, true);

  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 2, true);
  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 2, true);
  test_1x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 2, true);
}

void test_2x_gpu()
{
  // p=1
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 1, true);
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 1, true);
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 1, true);

  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 1, true);
  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 1, true);
  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 1, true);

  // p=2
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 2, true);
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 2, true);
  test_2x_nc1_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 2, true);

  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_NONE, 2, true);
  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_ABS, 2, true);
  test_2x_nc3_op(GKYL_ARRAY_INTEGRATE_OP_SQ, 2, true);
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
