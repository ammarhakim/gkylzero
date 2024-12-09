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
#include <gkyl_dg_bin_ops.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size, bool use_gpu)
{
  struct gkyl_array *a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}
//----------------- TEST 1x ------------------
// function to average
void evalFunc_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {-4.0}, upper[] = {6.0}; // Has to match the test below.
  double Lx = upper[0]-lower[0];
  double k_x = 2.*M_PI/Lx;
  double phi = 0.5;

  fout[0] = x*sin(k_x*x + phi);
}
// to weight the integral
void evalWeight_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1+x*x;
}
double solution_1x() { 
  // Solution from a trapz integration with Python (see code at the end)
  return -0.4189328208844751;
  }

void test_1x(int poly_order, bool use_gpu)
{
  //------------- 1. Define grid and basis ----------------
  double lower[] = {-4.0}, upper[] = {6.0};
  int cells[] = {32};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Define the reduced range and basis for averaging
  struct gkyl_range red_local;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);
  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  //------------- 2. Project the target function ----------------
  gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalFunc_1x, NULL);

  struct gkyl_array *fx_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(projf, 0.0, &local, fx_c);

  gkyl_proj_on_basis_release(projf);

  //------------- 3. Project the weight function ----------------
  gkyl_proj_on_basis *proj_weight = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalWeight_1x, NULL);

  struct gkyl_array *wx_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(proj_weight, 0.0, &local_ext, wx_c);

  gkyl_proj_on_basis_release(proj_weight);

  //------------- 4. Compute weighted average ----------------
    struct gkyl_array_average_inp inp_avg_full = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = red_basis,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &red_local,
    .weights = wx_c,
    .op = GKYL_ARRAY_AVERAGE_OP,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *avg_full = gkyl_array_average_new(&inp_avg_full);

  struct gkyl_array *avgf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(avg_full, fx_c, avgf_c);

  gkyl_array_average_release(avg_full);

  //------------- 5. Fetch and transfer results ----------------
  const double *avg_c0 = gkyl_array_cfetch(avgf_c, 0);
  double *avg_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(avg_c0_ho, avg_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(avg_c0_ho, avg_c0, sizeof(double));

  //------------- 6. Check results ----------------
  double result = avg_c0_ho[0]*0.5*sqrt(2);
  double rel_err = fabs(result - solution_1x()) / fabs(solution_1x());
  printf("\n Checking one-step average (X to scalar)\n");
  printf("Result: %g, solution: %g\n",result,solution_1x());
  printf("Relative error: %e\n", fabs(result - solution_1x()) / fabs(solution_1x()));
  TEST_CHECK(gkyl_compare(result, solution_1x(), 1e-4));

  //------------- 7. Clean up ----------------
  gkyl_array_release(avgf_c);
  gkyl_array_release(fx_c);
  gkyl_array_release(wx_c);
}
//------------------------------------------------------
//----------------- TEST 2x ----------------------------
//------------------------------------------------------
void evalFunc_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double lower[] = {-4., -3.}, upper[] = {6., 5.};
  double Lx = upper[0]-lower[0], Ly = upper[1]-lower[1];
  double k_x = 2.*M_PI/Lx, k_y = 2.*M_PI/Ly;
  double phi = 0.5;
  fout[0] = 1 + sin(k_x*x + k_y*y);
  fout[0] = x * y * sin(1.5*k_x*x + 0.75*k_y*y + phi) * cos(1.42*k_y*y);
}
void evalWeight_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = 1 + x*x + y*y;
}
double solution_2x(){
  // Solution from a trapz integration with Python (see code at the end)
  return -0.6715118302909872;
}

// One step averaging
void test_2x_1step(int poly_order, bool use_gpu)
{
  printf("\n");

  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {64, 32};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Define the reduced range and basis for averaging
  struct gkyl_range red_local;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // Project the target and weight functions
  gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalFunc_2x, NULL);

  struct gkyl_array *fxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxy_c);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *proj_weight = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalWeight_2x, NULL);

  struct gkyl_array *wxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(proj_weight, 0.0, &local_ext, wxy_c);

  gkyl_proj_on_basis_release(proj_weight);

  // Perform the one step average
  struct gkyl_array_average_inp inp_avg_xy = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = red_basis,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &red_local,
    .weights = wxy_c,
    .op = GKYL_ARRAY_AVERAGE_OP,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *avg_xy = gkyl_array_average_new(&inp_avg_xy);

  struct gkyl_array *avgf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(avg_xy, fxy_c, avgf_c);

  gkyl_array_average_release(avg_xy);

  //------------- Check results ----------------
  const double *avgf_c0 = gkyl_array_cfetch(avgf_c, 0);
  double *avgf_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu){
    gkyl_cu_memcpy(avgf_c0_ho, avgf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else{
    memcpy(avgf_c0_ho, avgf_c0, sizeof(double));
  }

  double result = avgf_c0_ho[0] * sqrt(2)/2;
  double rel_err = fabs(result - solution_2x()) / fabs(solution_2x());
  printf("\nChecking one-step average (XY to scalar)\n");
  printf("\tResult: %g, solution: %g\n",result,solution_2x());
  printf("\tRelative error: %e\n", rel_err);
  TEST_CHECK(gkyl_compare(rel_err, 0, 1e-4));

  //------------------ Clean up ------------------
  gkyl_array_release(avgf_c);
  gkyl_array_release(fxy_c);
  gkyl_array_release(wxy_c);
}


void test_2x_intx_inty(int poly_order, bool use_gpu)
{
  // Define grids and basis
  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {32, 24};
  int ghost[] = {0, 0};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_rect_grid grid_y;
  gkyl_rect_grid_init(&grid_y, 1, &lower[1], &upper[1], &cells[1]);
  
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  struct gkyl_basis basis_y;
  gkyl_cart_modal_serendip(&basis_y, 1, poly_order);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);
  
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_range local_y, local_y_ext;
  int remove_dir [GKYL_MAX_CDIM] = {0}; 
  remove_dir[0] = 1;
  int tmp_idx[GKYL_MAX_CDIM] = {0}; 
  tmp_idx[0] = local.lower[0];
  gkyl_range_deflate(&local_y_ext, &local_ext, remove_dir, tmp_idx);
  gkyl_range_deflate(&local_y, &local, remove_dir, tmp_idx);

  struct gkyl_range red_local;
  int local_lower_x[] = {local.lower[0]};
  gkyl_range_init(&red_local, 1, local_lower_x, local_lower_x);

  gkyl_proj_on_basis *proj_weight = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalWeight_2x, NULL);

  struct gkyl_array *wxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(proj_weight, 0.0, &local_ext, wxy_c);

  gkyl_proj_on_basis_release(proj_weight);

  //  Integration over x only
  struct gkyl_array_average_inp inp_int_x = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = basis_y,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &local_y,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP_Y,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_x = gkyl_array_average_new(&inp_int_x);

  struct gkyl_array *wy_c = mkarr(basis_y.num_basis, local_y_ext.volume, use_gpu);
  gkyl_array_average_advance(int_x, wxy_c, wy_c);

  gkyl_array_average_release(int_x);

  //  Integration over remaining dimensions (x)
      struct gkyl_array_average_inp inp_int_y = {
    .grid = &grid_y,
    .tot_basis = basis_y,
    .sub_basis = red_basis,
    .tot_rng = &local_y,
    .tot_rng_ext = &local_y_ext,
    .sub_rng = &red_local,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_y = gkyl_array_average_new(&inp_int_y);

  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);

  gkyl_array_average_advance(int_y, wy_c, intw_c);

  gkyl_array_average_release(int_y);

  // Retrieve the computed average from the device (if applicable)
  const double *intw_c0 = gkyl_array_cfetch(intw_c, 0);
  double *intw_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intw_c0_ho, intw_c0, sizeof(double));

  const double result = intw_c0_ho[0];

  //------------------ Check results ------------------
  double solution = 1333.33542278013;
  double rel_err = fabs(result - solution) / fabs(solution);
  printf("\nChecking two-step integration (XY to Y then Y to scalar)\n");   
  printf("\tResult: %g, solution: %g\n",result, solution);
  printf("\tRelative error: %e\n", rel_err );
  TEST_CHECK(gkyl_compare(rel_err, 0.0, 1e-4));

  //------------------ Clean up ------------------
  gkyl_array_release(wxy_c);
  gkyl_array_release(wy_c);
  gkyl_array_release(intw_c);
}

void test_2x_avgx_avgy(int poly_order, bool use_gpu)
{
  // Define grids and basis
  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {32, 24};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {0, 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_range local_y, local_y_ext;
  gkyl_range_init(&local_y, 1, &local.lower[1], &local.upper[1]);
  gkyl_range_init(&local_y_ext, 1, &local_ext.lower[1], &local_ext.upper[1]);

  struct gkyl_basis basis_y;
  gkyl_cart_modal_serendip(&basis_y, 1, poly_order);

  struct gkyl_rect_grid grid_y;
  gkyl_rect_grid_init(&grid_y, 1, &lower[1], &upper[1], &cells[1]);

  struct gkyl_range red_local;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // Project the target function and weights
  gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalFunc_2x, NULL);

  struct gkyl_array *fxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxy_c);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *proj_weight = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalWeight_2x, NULL);

  struct gkyl_array *wxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(proj_weight, 0.0, &local_ext, wxy_c);

  gkyl_proj_on_basis_release(proj_weight);


  //Average over x only

  // Create and run the array average updater to average on x only
  struct gkyl_array_average_inp inp_avg_x = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = basis_y,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &local_y,
    .weights = wxy_c,
    .op = GKYL_ARRAY_AVERAGE_OP_Y,
    .use_gpu = use_gpu
  };
  struct gkyl_array *fy_c = mkarr(basis_y.num_basis, local_y_ext.volume, use_gpu);

  struct gkyl_array_average *avg_x = gkyl_array_average_new(&inp_avg_x);
  gkyl_array_average_advance(avg_x, fxy_c, fy_c);
  
  gkyl_array_average_release(avg_x);
  /*
  fy_c is DG coeff of int[w(x,y) f(x,y)]dx / int[w(x,y)]dx
  */

  // obtain x integral of the weights too
  struct gkyl_array_average_inp inp_int_x = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = basis_y,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &local_y,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP_Y,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_x = gkyl_array_average_new(&inp_int_x);

  struct gkyl_array *wy_c = mkarr(basis_y.num_basis, local_y_ext.volume, use_gpu);
  gkyl_array_average_advance(int_x, wxy_c, wy_c);

  gkyl_array_average_release(int_x);
  /*
  wy_c is DG coeff of int[w(x,y)]dx
  */

  // we now remove manually the denominator
  gkyl_dg_mul_op_range(basis_y, 0, fy_c, 0, fy_c, 0, wy_c, &local_y);
  /*
  fy_c is DG coeff of int[w(x,y) f(x,y)]dx
  */

  //Average over y now
  struct gkyl_array_average_inp inp_int_y = {
    .grid = &grid_y,
    .tot_basis = basis_y,
    .sub_basis = red_basis,
    .tot_rng = &local_y,
    .tot_rng_ext = &local_y_ext,
    .sub_rng = &red_local,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_y = gkyl_array_average_new(&inp_int_y);

  struct gkyl_array *intf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_y, fy_c, intf_c);
  /*
  here intf_c is DG coeff of int[int[w(x,y) f(x,y)]dx]dy
  */

  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_y, wy_c, intw_c);
  /*
  here intw_c is DG coeff of int[int[w(x,y)]dx]dy
  */

  gkyl_array_average_release(int_y);

  // Retrieve the computed average from the device (if applicable)
  const double *intf_c0  = gkyl_array_cfetch(intf_c, 0);
  double *intf_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(intf_c0_ho, intf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intf_c0_ho, intf_c0, sizeof(double));

  const double *intw_c0 = gkyl_array_cfetch(intw_c, 0);
  double *intw_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intw_c0_ho, intw_c0, sizeof(double));

  double integral_wf = intf_c0_ho[0];
  double integral_w = intw_c0_ho[0];
  double result = integral_wf/integral_w;
  double rel_err = fabs(result - solution_2x()) / fabs(solution_2x());
  //------------------ 5. Check results two step avg ------------------
  printf("\nChecking two-step average (XY to Y then Y to scalar)\n");   
  printf("\tResult: %g, solution: %g\n",result, solution_2x());
  printf("\tRelative error: %e\n", rel_err);
  TEST_CHECK(gkyl_compare(rel_err, 0, 1e-2));


  //------------------ Clean up ------------------
  gkyl_array_release(fxy_c);
  gkyl_array_release(wxy_c);
  gkyl_array_release(fy_c);
  gkyl_array_release(wy_c);
  gkyl_array_release(intf_c);
  gkyl_array_release(intw_c);
}


void test_2x_avgy_avgx(int poly_order, bool use_gpu)
{
  printf("\n");

  //------------------ 1. Initialization ------------------
  // 1.1 Define grid and basis
  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {32, 24};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  // Initialize the grid
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Initialize the polynomial basis
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Create ranges (local and extended, no ghost cells in this case)
  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // 1.2 Project the target function
  // Create a projection updater for the target function
  gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalFunc_2x, NULL);

  // Create an array to store the projected target function
  struct gkyl_array *fxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  // Project the target function onto the basis
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxy_c);
  gkyl_proj_on_basis_release(projf);

  // 1.3 Project the weight function
  // Create a projection updater for the weight function
  gkyl_proj_on_basis *proj_weight = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalWeight_2x, NULL);

  // Create an array to store the projected weight function
  struct gkyl_array *wxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  // Project the weight function onto the basis
  gkyl_proj_on_basis_advance(proj_weight, 0.0, &local_ext, wxy_c);
  gkyl_proj_on_basis_release(proj_weight);

  //------------------ 2. Two steps average ------------------
  //  5.1 Average over y only
  // Define the reduced range and basis for averaging along x only
  struct gkyl_range local_x, local_x_ext;
  gkyl_range_init(&local_x, 1, &local.lower[0], &local.upper[0]);
  gkyl_range_init(&local_x_ext, 1, &local_ext.lower[0], &local_ext.upper[0]);

  struct gkyl_basis basis_x;
  gkyl_cart_modal_serendip(&basis_x, 1, poly_order);
  // Create an array to store the averaged result
  struct gkyl_array *fx_c = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);
  struct gkyl_array *wx_c = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);

  // Create and run the array average updater to average on x only
  struct gkyl_array_average_inp inp_avg_x = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = basis_x,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &local_x,
    .weights = wxy_c,
    .op = GKYL_ARRAY_AVERAGE_OP_X,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *avg_x = gkyl_array_average_new(&inp_avg_x);
  gkyl_array_average_advance(avg_x, fxy_c, fx_c);
  gkyl_array_average_release(avg_x);
  /*
  fx_c is DG coeff of int[w(x,y) f(x,y)]dx / int[w(x,y)]dx
  */

  // obtain x integral of the weights too
  struct gkyl_array_average_inp inp_int_x = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = basis_x,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &local_x,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP_X,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_x = gkyl_array_average_new(&inp_int_x);
  gkyl_array_average_advance(int_x, wxy_c, wx_c);
  gkyl_array_average_release(int_x);
  /*
  wx_c is DG coeff of int[w(x,y)]dx
  */

  // we now remove manually the denominator
  gkyl_dg_mul_op_range(basis_x, 0, fx_c, 0, fx_c, 0, wx_c, &local_x);
  /*
  fx_c is DG coeff of int[w(x,y) f(x,y)]dx
  */

  //  5.2 Average over y now
  struct gkyl_rect_grid grid_x;
  gkyl_rect_grid_init(&grid_x, 1, &lower[0], &upper[0], &cells[0]);

  // Define the reduced range and basis for scalar result
  struct gkyl_range red_local;
  gkyl_range_init(&red_local, 1, &local.lower[1], &local.lower[1]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // Create an array to store the averaged result
  struct gkyl_array *intf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);

  // Create and run the array average updater to average on x only
      struct gkyl_array_average_inp inp_int_y = {
    .grid = &grid_x,
    .tot_basis = basis_x,
    .sub_basis = red_basis,
    .tot_rng = &local_x,
    .tot_rng_ext = &local_x_ext,
    .sub_rng = &red_local,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_y = gkyl_array_average_new(&inp_int_y);
  // printf("-compute int(wf)dy\n");
  gkyl_array_average_advance(int_y, fx_c, intf_c);
  /*
  here intf_c is DG coeff of int[int[w(x,y) f(x,y)]dx]dy
  */

  // obtain full integral of weights too
  gkyl_array_average_advance(int_y, wx_c, intw_c);
  gkyl_array_average_release(int_y);
  /*
  here intw_c is DG coeff of int[int[w(x,y)]dx]dy
  */

  // Retrieve the computed average from the device (if applicable)
  const double *intf_c0  = gkyl_array_cfetch(intf_c, 0);
  double *intf_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(intf_c0_ho, intf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intf_c0_ho, intf_c0, sizeof(double));

  const double *intw_c0 = gkyl_array_cfetch(intw_c, 0);
  double *intw_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intw_c0_ho, intw_c0, sizeof(double));

  double integral_wf = intf_c0_ho[0];
  double integral_w = intw_c0_ho[0];
  double result = integral_wf/integral_w;
  double rel_err = fabs(result - solution_2x()) / fabs(solution_2x());
  //------------------ 5. Check results two step avg ------------------
  printf("Checking two-step average (XY to X then X to scalar)\n");   
  printf("\tResult: %g, solution: %g\n",result,solution_2x());
  printf("\tRelative error: %e\n", rel_err);
  TEST_CHECK(gkyl_compare(rel_err, 0.0, 1e-3));


  //------------------ Clean up ------------------
  gkyl_array_release(fxy_c);
  gkyl_array_release(wxy_c);
  gkyl_array_release(fx_c);
  gkyl_array_release(wx_c);
  gkyl_array_release(intf_c);
  gkyl_array_release(intw_c);
}

//----------------- TEST 3x ------------------
// function to average
void evalFunc_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  double lower[] = {-4., -3., -2.}, upper[] = {6., 5., 4.};
  double Lx = upper[0]-lower[0];
  double Ly = upper[1]-lower[1];
  double Lz = upper[2]-lower[2];
  double k_x = 2.*M_PI/Lx;
  double k_y = 2.*M_PI/Ly;
  double k_z = 2.*M_PI/Lz;
  double phi = 0.5;
  fout[0] = x * y * z * sin(1.5*k_x*x + 0.75*k_y*y + 0.5*k_z*z + phi) * cos(1.42*k_y*y) * cos(4.20*k_z*z);
}
// to weight the integral
void evalWeight_3x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = 1 + x*x + y*y + z*z;
}
// Solution from a trapz integration with Python (see code at the end)
double solution_3x(){
  return 0.06652603651728471;
}

// two setps first avg on x then on y and z
void test_3x_avgx_avgyz(int poly_order, bool use_gpu)
{
  // Define grids and basis
  double lower[] = {-4.0, -3.0, -2.0}, upper[] = {6.0, 5.0, 4.0};
  int cells[] = {64, 48, 32};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // 2dim subgrid
  struct gkyl_rect_grid grid_yz;
  double yz_grid_lower[] = {grid.lower[1], grid.lower[2]};
  double yz_grid_upper[] = {grid.upper[1], grid.upper[2]};
  int    yz_grid_cells[] = {grid.cells[1], grid.cells[2]};
  gkyl_rect_grid_init(&grid_yz, 2, yz_grid_lower, yz_grid_upper, yz_grid_cells);

  // Define the reduced range and basis for averaging along x only
  int ghost_yz[] = {ghost[1], ghost[2]};
  struct gkyl_range local_yz, local_yz_ext;
  gkyl_create_grid_ranges(&grid_yz, ghost_yz, &local_yz_ext, &local_yz);
  struct gkyl_basis basis_yz;
  gkyl_cart_modal_serendip(&basis_yz, 2, poly_order);

  // Define the reduced range and basis for full reduction
  struct gkyl_range red_local;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // Project the target function and weights
  gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalFunc_3x, NULL);

  struct gkyl_array *fxyz_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxyz_c);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *proj_weight = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalWeight_3x, NULL);

  struct gkyl_array *wxyz_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(proj_weight, 0.0, &local_ext, wxyz_c);

  gkyl_proj_on_basis_release(proj_weight);

  // Create and run the array average updater to average on x only
  struct gkyl_array_average_inp inp_avg_xyz_to_yz = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = basis_yz,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &local_yz,
    .weights = wxyz_c,
    .op = GKYL_ARRAY_AVERAGE_OP_YZ,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *avg_xyz_to_yz = gkyl_array_average_new(&inp_avg_xyz_to_yz);

  struct gkyl_array *fyz_c = mkarr(basis_yz.num_basis, local_yz_ext.volume, use_gpu);
  gkyl_array_average_advance(avg_xyz_to_yz, fxyz_c, fyz_c);
  // fy_c is DG coeff of int[w(x,y,z) f(x,y,z)]dx / int[w(x,y,z)]dx
  gkyl_array_average_release(avg_xyz_to_yz);

  // obtain x integral of the weights too
  struct gkyl_array_average_inp inp_int_xyz_to_yz = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = basis_yz,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &local_yz,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP_YZ,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_xyz_to_yz = gkyl_array_average_new(&inp_int_xyz_to_yz);

  struct gkyl_array *wyz_c = mkarr(basis_yz.num_basis, local_yz_ext.volume, use_gpu);
  gkyl_array_average_advance(int_xyz_to_yz, wxyz_c, wyz_c);
  // wy_c is DG coeff of int[w(x,y)]dy
  gkyl_array_average_release(int_xyz_to_yz);

  // we now remove manually the denominator
  gkyl_dg_mul_op_range(basis_yz, 0, fyz_c, 0, fyz_c, 0, wyz_c, &local_yz);
  // fy_c is DG coeff of int[w(x,y) f(x,y)]dy

  // Create and run the array average updater to average on x only
      struct gkyl_array_average_inp inp_int_y = {
    .grid = &grid_yz,
    .tot_basis = basis_yz,
    .sub_basis = red_basis,
    .tot_rng = &local_yz,
    .tot_rng_ext = &local_yz_ext,
    .sub_rng = &red_local,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_y = gkyl_array_average_new(&inp_int_y);
  struct gkyl_array *intf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_y, fyz_c, intf_c);
  // here intf_c is DG coeff of int[int[w(x,y) f(x,y)]dy]dx

  // obtain full integral of weights too
  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_y, wyz_c, intw_c);
  // here intw_c is DG coeff of int[int[w(x,y)]dy]dx

  gkyl_array_average_release(int_y);

  // Retrieve the computed average from the device (if applicable)
  const double *intf_c0  = gkyl_array_cfetch(intf_c, 0);
  double *intf_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(intf_c0_ho, intf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intf_c0_ho, intf_c0, sizeof(double));

  const double *intw_c0 = gkyl_array_cfetch(intw_c, 0);
  double *intw_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intw_c0_ho, intw_c0, sizeof(double));

  const double integral_wf = intf_c0_ho[0];
  const double integral_w = intw_c0_ho[0];
  const double result = integral_wf/integral_w;
  const double solution = solution_3x();
  const double rel_err = fabs(result - solution) / fabs(solution);

  // Check results two step avg
  printf("\nChecking two-step average (XYZ to YZ then YZ to scalar)\n");   
  printf("\tResult: %g, solution: %g\n",result,solution);
  printf("\tRelative error: %e\n", rel_err);
  TEST_CHECK(gkyl_compare(rel_err, 0, 1e-3));

  //------------------ Clean up ------------------
  gkyl_array_release(fxyz_c);
  gkyl_array_release(wxyz_c);
  gkyl_array_release(fyz_c);
  gkyl_array_release(wyz_c);
  gkyl_array_release(intf_c);
  gkyl_array_release(intw_c);
}

void test_3x_avgyz_avgx(int poly_order, bool use_gpu)
{
  // Define grids and basis
  double lower[] = {-4.0, -3.0, -2.0}, upper[] = {6.0, 5.0, 4.0};
  int cells[] = {128, 64, 48};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // 2dim subgrid
  struct gkyl_rect_grid grid_x;
  gkyl_rect_grid_init(&grid_x, 1, &grid.lower[0], &grid.upper[0], &grid.cells[0]);

  // Define the reduced range and basis for averaging along x only
  struct gkyl_range local_x, local_x_ext;
  gkyl_create_grid_ranges(&grid_x, &ghost[0], &local_x_ext, &local_x);
  struct gkyl_basis basis_x;
  gkyl_cart_modal_serendip(&basis_x, 1, poly_order);

  // Define the reduced range and basis for full reduction
  struct gkyl_range red_local;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // Project the target function and weights
  gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalFunc_3x, NULL);

  struct gkyl_array *fxyz_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxyz_c);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *proj_weight = gkyl_proj_on_basis_new(
      &grid, &basis, poly_order + 1, 1, evalWeight_3x, NULL);

  struct gkyl_array *wxyz_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_on_basis_advance(proj_weight, 0.0, &local_ext, wxyz_c);

  gkyl_proj_on_basis_release(proj_weight);


  // Create and run the array average updater to average on x only
  struct gkyl_array_average_inp inp_avg_xyz_to_x = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = basis_x,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &local_x,
    .weights = wxyz_c,
    .op = GKYL_ARRAY_AVERAGE_OP_X,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *avg_xyz_to_x = gkyl_array_average_new(&inp_avg_xyz_to_x);

  struct gkyl_array *fx_c = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);
  gkyl_array_average_advance(avg_xyz_to_x, fxyz_c, fx_c);
  // fy_c is DG coeff of int[w(x,y,z) f(x,y,z)]dydz / int[w(x,y,z)]dydz

  gkyl_array_average_release(avg_xyz_to_x);

  // obtain x integral of the weights too
  struct gkyl_array_average_inp inp_int_xyz_to_x = {
    .grid = &grid,
    .tot_basis = basis,
    .sub_basis = basis_x,
    .tot_rng = &local,
    .tot_rng_ext = &local_ext,
    .sub_rng = &local_x,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP_X,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_xyz_to_x = gkyl_array_average_new(&inp_int_xyz_to_x);

  struct gkyl_array *wx_c = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);
  gkyl_array_average_advance(int_xyz_to_x, wxyz_c, wx_c);
  // wy_c is DG coeff of int[w(x,y)]dydz

  gkyl_array_average_release(int_xyz_to_x);

  // we now remove manually the denominator
  gkyl_dg_mul_op_range(basis_x, 0, fx_c, 0, fx_c, 0, wx_c, &local_x);
  // fy_c is DG coeff of int[w(x,y) f(x,y)]dy

  // Create and run the array average updater to average on x only
      struct gkyl_array_average_inp inp_int_x = {
    .grid = &grid_x,
    .tot_basis = basis_x,
    .sub_basis = red_basis,
    .tot_rng = &local_x,
    .tot_rng_ext = &local_x_ext,
    .sub_rng = &red_local,
    .weights = NULL,
    .op = GKYL_ARRAY_AVERAGE_OP,
    .use_gpu = use_gpu
  };
  struct gkyl_array_average *int_x = gkyl_array_average_new(&inp_int_x);

  struct gkyl_array *intf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_x, fx_c, intf_c);
  // intf_c is DG coeff of int[int[int[w(x,y) f(x,y)]dy]dz]dx

  // obtain full integral of weights too
  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_x, wx_c, intw_c);
  // intw_c is DG coeff of int[int[int[w(x,y)]dy]dz]dx

  gkyl_array_average_release(int_x);

  // Retrieve the computed average from the device (if applicable)
  const double *intf_c0  = gkyl_array_cfetch(intf_c, 0);
  double *intf_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(intf_c0_ho, intf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intf_c0_ho, intf_c0, sizeof(double));

  const double *intw_c0 = gkyl_array_cfetch(intw_c, 0);
  double *intw_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu)
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(intw_c0_ho, intw_c0, sizeof(double));

  const double integral_wf = intf_c0_ho[0];
  const double integral_w = intw_c0_ho[0];
  const double result = integral_wf/integral_w;
  const double solution = solution_3x();
  const double rel_err = fabs(result - solution) / fabs(solution);

  //------------------ 5. Check results two step avg ------------------
  printf("\nChecking two-step average (XYZ to X then X to scalar)\n");   
  printf("\tResult: %g, solution: %g\n",result,solution);
  printf("\tRelative error: %e\n", rel_err);
  TEST_CHECK(gkyl_compare(rel_err, 0, 1e-3));

  //------------------ Clean up ------------------
  gkyl_array_release(fxyz_c);
  gkyl_array_release(wxyz_c);
  gkyl_array_release(fx_c);
  gkyl_array_release(wx_c);
  gkyl_array_release(intf_c);
  gkyl_array_release(intw_c);
}

void test_1x_cpu()
{
  // p=1
  test_1x(1, false);

  // p=2
  // test_1x(2, false);

}

void test_2x_cpu_1step()
{
  // p=1
  test_2x_1step(1, false);

  // p=2
  // test_2x(2, false);
}

void test_2x_cpu_2steps()
{
  // p=1
  test_2x_intx_inty(1, false);
  test_2x_avgx_avgy(1, false);
  test_2x_avgy_avgx(1, false);

  // p=2
  // test_2x_2steps(2, false);
}

void test_3x_cpu()
{
  // p=1
  test_3x_avgx_avgyz(1, false);
  test_3x_avgyz_avgx(1, false);

  // p=2
  // test_3x(2, false);
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

void test_3x_gpu()
{
  // p=1
  test_3x(1, true);

  // p=2
  // test_3x(2, true);
}
#endif

TEST_LIST = {
  { "test_1x_cpu", test_1x_cpu },
  { "test_2x_cpu_1step", test_2x_cpu_1step },
  { "test_2x_cpu_2steps", test_2x_cpu_2steps },
  { "test_3x_cpu", test_3x_cpu },
#ifdef GKYL_HAVE_CUDA
  { "test_1x_gpu", test_1x_gpu },
  { "test_2x_gpu_1step", test_2x_gpu_1step },
  { "test_2x_gpu_2steps", test_2x_gpu_2steps },
  { "test_3x_gpu", test_3x_gpu },
#endif
  { NULL, NULL },
};


//-------- PYTHON CODE FOR THE SOLUTION OF TEST_1X
/*
import numpy as np
from scipy.integrate import quad

# Define the function to integrate
def f(x, k_x, phi):
    return x * np.sin(k_x * x + phi)

def w(x):
    return 1+x**2

def fw(x, k_x, phi):
    return w(x)*x * np.sin(k_x * x + phi)


# Parameters
a, b = -4, 6
k_x = 2.0 * np.pi / (b-a)  # Wave number
phi = 0.5  # Phase

# Perform numerical integration
int_fw, efw = quad(fw, a, b, args=(k_x, phi))
int_w,  ew = quad(w, a, b)

# Result
normalized_result = int_fw / int_w
total_error = efw + ew

print("Normalized Result:", normalized_result)
print("Total Error:", total_error)
*/
/* OUTPUT
int(wf): -43.2897248247291
int(w): 103.33333333333334
Normalized Result: -0.4189328208844751
*/

//-------- PYTHON CODE FOR THE SOLUTION OF TEST_2X
/*
import numpy as np

# Define the function to integrate
def f(x, y, k_x, k_y, phi):
    return x * y * np.sin(1.5 * k_x * x + 0.75*k_y * y + phi) * np.cos(1.42*k_y*y)

def w(x, y):
    return 1 + x**2 + y**2

def fw(x, y, k_x, k_y, phi):
    return w(x, y) * f(x, y, k_x, k_y, phi)

# Parameters
ax, bx = -4, 6  # Bounds for x
ay, by = -3, 5  # Bounds for y
k_x = 2.0 * np.pi / (bx - ax)  # Wave number for x
k_y = 2.0 * np.pi / (by - ay)  # Wave number for y
phi = 0.5  # Phase

# Create grid points
nx, ny = 1024, 1024  # Number of grid points in x and y
x = np.linspace(ax, bx, nx); y = np.linspace(ay, by, ny)
dx = (bx - ax) / (nx - 1); dy = (by - ay) / (ny - 1)
# Create 2D grids
X, Y = np.meshgrid(x, y, indexing="ij")

# Evaluate the functions on the grid
f_values = f(X, Y, k_x, k_y, phi)
w_values = w(X, Y)
fw_values = fw(X, Y, k_x, k_y, phi)

# Perform integration using the trapezoidal rule
int_fw = np.trapz(np.trapz(fw_values, x=y, axis=1), x=x)
int_w = np.trapz(np.trapz(w_values, x=y, axis=1), x=x)

# Normalize the result
normalized_result = int_fw / int_w

# Output
print("Normalized Result:", normalized_result)
*/
/* OUTPUT
int(wf): -895.3505101428923
int(w): 1333.33542278013
Normalized Result: -0.6715118302909872
*/

//-------- PYTHON CODE FOR THE SOLUTION OF TEST_3X
/*
import numpy as np

# Define the function to integrate
def f(x, y, z, k_x, k_y, k_z, phi):
    # return 1 + np.sin(k_x*x)
    return x * y * z * np.sin(1.5*k_x*x + 0.75*k_y*y + 0.5*k_z*z + phi) * np.cos(1.42*k_y*y) * np.cos(4.20*k_z*z)

def w(x, y, z):
    return np.ones(np.shape(x)) + x**2 + y**2 + z**2

def fw(x, y, z, k_x, k_y, k_z, phi):
    return w(x, y, z) * f(x, y, z, k_x, k_y, k_z, phi)

# Parameters
ax, bx = -4, 6  # Bounds for x
ay, by = -3, 5  # Bounds for y
az, bz = -2, 4  # Bounds for z
k_x = 2.0 * np.pi / (bx - ax)  # Wave number for x
k_y = 2.0 * np.pi / (by - ay)  # Wave number for y
k_z = 2.0 * np.pi / (bz - az)  # Wave number for y
phi = 0.5  # Phase

# Create grid points
nx, ny, nz = 256,256,256 # Number of grid points in x and y
x = np.linspace(ax, bx, nx); y = np.linspace(ay, by, ny); z = np.linspace(az, bz, nz)
dx = (bx - ax) / (nx - 1); dy = (by - ay) / (ny - 1)
# Create 2D grids
X, Y, Z = np.meshgrid(x, y, z)

# Evaluate the functions on the grid
f_values = f(X, Y, Z, k_x, k_y, k_z, phi)
w_values = w(X, Y, Z)
fw_values = fw(X, Y, Z, k_x, k_y, k_z, phi)

# Perform integration using the trapezoidal rule
int_fw = np.trapz(np.trapz(np.trapz(fw_values, x=z, axis=2), x=y, axis=1), x=x)
int_w  = np.trapz(np.trapz(np.trapz( w_values, x=z, axis=2), x=y, axis=1), x=x)

# Normalize the result
normalized_result = int_fw / int_w

# Output
print("int(wf):", int_fw)
print("int(w):", int_w)
print("Normalized Result:", normalized_result)
*/
/* OUTPUT
int(wf): 659.9546515953209
int(w): 9920.246059207997
Normalized Result: 0.06652603651728471
*/