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

// to compute the volume
void evalFunc_one(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 1.0;
}

// to weight the integral
void evalFunc_weight(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 10.0;
}
// function to evaluate in the 1x testcase
void evalFunc_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {-6.0}, upper[] = {6.0}; // Has to match the test below.
  double Lx = upper[0]-lower[0];
  double k_x = 2.*M_PI/Lx;
  double phi = 0.0;

  fout[0] = 1./Lx;
  // fout[0] = sin(k_x*x + phi);
  // fout[0] = 2.0*pow(x,2) + 1.3*x - 1.2;
}

void evalFunc_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double lower[] = {0., -6.0}, upper[] = {2., 6.0}; // Has to match the test below.
  double Lx = upper[0]-lower[0];
  double Ly = upper[1]-lower[1];
  double k_x = 2.*M_PI/Lx;
  double k_y = 2.*M_PI/Lx;
  double phi = 0.5;

  // fout[0] = 1./(Lx*Ly);
  fout[0] = sin(k_x*x) * cos(k_y*y + phi);
  // fout[0] = 2.0*pow(x,2) - 3.0*pow(y,2) + 1.3*x*y - 1.2;
}

//----------------- TESTS ------------------

void test_1x(int poly_order, bool use_gpu)
{
  printf("\n");
  double lower[] = {-6.0}, upper[] = {6.0};
  int cells[] = {16};
  int ndim = sizeof(lower)/sizeof(lower[0]);

  //------------------ 1. Define grid and basis ------------------
  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = { 1 };
  struct gkyl_range local, local_ext; // local, local-ext ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  //------------------ 2. Project the const functions to evaluate the volume of the domain ------------------
  //.projection updater for the volume evaluation
  gkyl_proj_on_basis *proj_one;
  proj_one = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, evalFunc_one, NULL);
  //.create and project the const function
  struct gkyl_array *one_ga = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  //.project distribution function on basis
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local, one_ga);
  gkyl_proj_on_basis_release(proj_one);

  //------------------ 3. Project the target function ------------------
  //.projection updater for dist-function
  gkyl_proj_on_basis *projf;
  projf = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, evalFunc_1x, NULL);
  //.create distribution function array
  struct gkyl_array *distf = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  //.project distribution function on basis
  gkyl_proj_on_basis_advance(projf, 0.0, &local, distf);
  gkyl_proj_on_basis_release(projf);

  //------------------ 4. Average through array integrate of target function and volume function ------------------
  //.integrate updater
  struct gkyl_array_integrate *integ_up = gkyl_array_integrate_new(&grid, &basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, use_gpu);
  double *fint = use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  double *vint = use_gpu? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));
  struct gkyl_array *unused_weights = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  // Integrate the input function
  gkyl_array_integrate_advance(integ_up, distf, 1., unused_weights, &local, fint);
  // Integrate the space to get the volume
  gkyl_array_integrate_advance(integ_up, one_ga,  1., unused_weights, &local, vint);
  // release the updater
  gkyl_array_integrate_release(integ_up);
  gkyl_array_release(unused_weights);

  // Get the integration results back to the host
  double *fint_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(fint_ho, fint, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(fint_ho, fint, sizeof(double));
    
  double *vint_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(vint_ho, vint, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(vint_ho, vint, sizeof(double));
  
  //------------------ 5. Average through the array average routine ------------------
  // declare the lower skin (ls) range to perform the x average
  struct gkyl_range ls_rng;
  gkyl_range_lower_skin(&ls_rng, &local, 0, 1);
  // const basis
  struct gkyl_basis const_basis;
  gkyl_cart_modal_serendip(&const_basis, 1, 0);
  // declare a single cell gkyl array (integral will be first coeff)
  struct gkyl_array *avg_res = mkarr(const_basis.num_basis, ls_rng.volume, use_gpu);

  //.projection updater for the weight evaluation
  gkyl_proj_on_basis *proj_weight;
  proj_weight = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, evalFunc_weight, NULL);
  //.create and project the const function
  struct gkyl_array *weight = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  //.project distribution function on basis
  gkyl_proj_on_basis_advance(proj_weight, 0.0, &local, weight);
  gkyl_proj_on_basis_release(proj_weight);

  // create full average updater and advance it
  struct gkyl_array_average *avg_full;
  avg_full = gkyl_array_average_new(&grid, basis, const_basis, local, ls_rng, weight, GKYL_ARRAY_AVERAGE_OP, use_gpu);
  // run the updater (this will integrate)
  gkyl_array_average_advance(avg_full, distf, avg_res);
  // release the updater
  gkyl_array_average_release(avg_full);
  // Fetch the integration result and send it back to the host
  const double *avg = gkyl_array_cfetch(avg_res, 0);
  double *avg_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(avg_ho, avg, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(avg_ho, avg, sizeof(double));

  //------------------ 6. Checks ------------------
  printf("Volume integral : %g\n", vint_ho[0]);
  printf("Average array results: %g\n",avg_ho[0]);
  printf("Ground truth: %g\n", fint_ho[0]);
  TEST_CHECK( gkyl_compare( avg_ho[0], fint_ho[0], 1e-12) );

  gkyl_array_release(avg_res);
  gkyl_array_release(distf);
  gkyl_array_release(one_ga);
  gkyl_array_release(weight);

  gkyl_free(fint_ho);
  if (use_gpu)
    gkyl_cu_free(fint);
  else
    gkyl_free(fint);
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

  // projection updater for constant one function
  gkyl_proj_on_basis *proj_one;
  proj_one = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, nc, evalFunc_one, NULL);

  // create distribution function array
  struct gkyl_array *one_ga = mkarr(nc*basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *one_ga_ho = use_gpu? mkarr(nc*basis.num_basis, local_ext.volume, false) : one_ga;

  // project distribution function on basis
  gkyl_proj_on_basis_advance(proj_one, 0.0, &local, one_ga_ho);
  if (use_gpu) gkyl_array_copy(one_ga, one_ga_ho);


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

  // 1D basis functions
  struct gkyl_basis basis_x;
  gkyl_cart_modal_serendip(&basis_x, 1, poly_order);
  // 1D grid
  struct gkyl_rect_grid grid_x;
  gkyl_rect_grid_init(&grid_x, 1, &lower[0], &upper[0], &cells[0]);
  // 1D Ranges
  int ghost_x[] = { ghost[0] };
  struct gkyl_range local_x, local_x_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&grid_x, ghost_x, &local_x_ext, &local_x);
  // 1D gkyl array
  struct gkyl_array *favg_x = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);
  struct gkyl_array *favg_x_ho = use_gpu? mkarr(favg_x->ncomp, favg_x->size, false)
                                        : gkyl_array_acquire(favg_x);

  // Create an array average updater
  struct gkyl_array_average *avg_yz;
  gkyl_array_average_new(&grid, basis, basis_x, local, local_x, one_ga, GKYL_ARRAY_AVERAGE_OP_X, use_gpu);
  gkyl_array_average_advance(avg_yz, distf, favg_x);
  gkyl_array_average_release(avg_yz);

  // declare the lower skin (ls) structure to perform the x average
  struct gkyl_range ls_rng;
  gkyl_range_lower_skin(&ls_rng, &local_x_ext, 0, 1);
  // const basis
  struct gkyl_basis const_basis;
  gkyl_cart_modal_serendip(&const_basis, 1, 0);
  printf("const_basis.num_basis = %d",const_basis.num_basis);
  // declare a single cell gkyl array (integral will be first coeff)
  struct gkyl_array *favg_arr = mkarr(const_basis.num_basis, ls_rng.volume, use_gpu);

  // projection updater for constant one function
  gkyl_proj_on_basis *proj_one_x;
  proj_one_x = gkyl_proj_on_basis_new(&grid_x, &basis_x, poly_order+1, nc, evalFunc_one, NULL);

  // create distribution function array
  struct gkyl_array *one_ga_x = mkarr(nc*basis_x.num_basis, local_x_ext.volume, use_gpu);
  struct gkyl_array *one_ga_x_ho = use_gpu? mkarr(nc*basis.num_basis, local_ext.volume, false) : one_ga_x;

  // project distribution function on basis
  gkyl_proj_on_basis_advance(proj_one_x, 0.0, &local_x, one_ga_x_ho);
  if (use_gpu) gkyl_array_copy(one_ga_x, one_ga_x_ho);

  // create full average updater and advance it
  struct gkyl_array_average *avg_x;
  gkyl_array_average_new(&grid_x, basis_x, const_basis,  local_x, ls_rng, NULL, GKYL_ARRAY_AVERAGE_OP, use_gpu);
  gkyl_array_average_advance(avg_x, favg_x, favg_arr);
  gkyl_array_average_release(avg_x);


  // check output
  const double *favg = gkyl_array_cfetch(favg_arr, 0);
  double *favg_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(favg_ho, favg, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(favg_ho, favg, nc*sizeof(double));

  printf("Average array results: %g\n",favg_ho[0]);
  printf("Ground truth: %g\n", fint_ho[0]);
  TEST_CHECK( gkyl_compare( favg_ho[0], fint_ho[0], 1e-12) );

  gkyl_array_release(favg_x);
  gkyl_array_release(favg_x_ho);
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
  // { "test_2x_cpu", test_2x_cpu },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_gpu", test_1x_gpu },
  { "test_2x_gpu", test_2x_gpu },
#endif
  { NULL, NULL },
};
