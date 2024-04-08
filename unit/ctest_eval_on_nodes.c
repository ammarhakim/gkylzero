#include <acutest.h>

#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <math.h>

void evalFunc_1x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = x*x;
}

void evalFunc_1x_trig(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = sin((2.*M_PI/(4.))*x);
}

void test_1x(int poly_order, int test_func_op)
{
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evup;
  if (test_func_op==0)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_1x_quad, NULL);
  else if (test_func_op==1)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_1x_trig, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_eval_on_nodes_advance(evup, 0.0, &arr_range, distf);

      gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, "ctest_eval_on_nodes_distf.gkyl");

  double *dfl = gkyl_array_fetch(distf, 0);  // left cell
  double *dfr = gkyl_array_fetch(distf, 1);  // right cell
  if (poly_order == 1) {
    if (test_func_op==0) {
      TEST_CHECK( gkyl_compare( 2.82842712474619e+00, dfl[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.63299316185545e+00, dfl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.82842712474619e+00, dfr[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.63299316185545e+00, dfr[1], 1e-12) );
    } else if (test_func_op==1) {
      TEST_CHECK( gkyl_compare(-8.65956056235493e-17, dfl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.99959962173949e-17, dfl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 8.65956056235493e-17, dfr[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.99959962173949e-17, dfr[1], 1e-12) );
    }
  } else if (poly_order == 2) {
    if (test_func_op==0) {
      TEST_CHECK( gkyl_compare( 1.88561808316413e+00, dfl[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.63299316185545e+00, dfl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.21637021355784e-01, dfl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.88561808316413e+00, dfr[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.63299316185545e+00, dfr[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.21637021355784e-01, dfr[2], 1e-12) );
    } else if (test_func_op==1) {
      TEST_CHECK( gkyl_compare(-9.42809041582064e-01, dfl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.99959962173949e-17, dfl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.21637021355784e-01, dfl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 9.42809041582064e-01, dfr[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.99959962173949e-17, dfr[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-4.21637021355784e-01, dfr[2], 1e-12) );
    }
  }

  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(distf);
}

void evalFunc_2x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = x*x + x*y + y*y;
}

void evalFunc_2x_trig(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = sin((2.*M_PI/(4.))*x) + sin((2.*M_PI/(4.))*x)*cos((2.*M_PI/(4.))*y) + cos((2.*M_PI/(4.))*y);
}

void test_2x(int poly_order, int test_func_op)
{
  double lower[] = {-2.0, -2.0}, upper[] = {2.0, 2.0};
  int cells[] = {2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evup;
  if (test_func_op==0)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_2x_quad, NULL);
  else if (test_func_op==1)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_2x_trig, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_eval_on_nodes_advance(evup, 0.0, &arr_range, distf);

  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, "ctest_eval_on_nodes_distf.gkyl");

  double *dfll = gkyl_array_fetch(distf, 0);  // left, low cell
  double *dflu = gkyl_array_fetch(distf, 1);  // left, up cell
  double *dfrl = gkyl_array_fetch(distf, 2);  // right, low cell
  double *dfru = gkyl_array_fetch(distf, 3);  // right, up cell
  if (poly_order == 1) {
    if (test_func_op==0) {
      TEST_CHECK( gkyl_compare( 1.00000000000000e+01, dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513776e+00, dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513776e+00, dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.00000000000000e+00, dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00, dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00, dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.00000000000000e+00, dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00, dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00, dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.00000000000000e+01, dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513776e+00, dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513776e+00, dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfru[3], 1e-12) );
    } else if (test_func_op==1) {
      TEST_CHECK( gkyl_compare(-1.11022302462516e-16 , dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17 , dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00 , dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.77555756156289e-17 , dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.11022302462516e-16 , dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17 , dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00 , dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17 , dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.11022302462516e-16 , dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17 , dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00 , dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.77555756156289e-17 , dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00 , dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17 , dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00 , dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17 , dfru[3], 1e-12) );
    }
  } else if (poly_order == 2) {
    if (test_func_op==0) {
      TEST_CHECK( gkyl_compare( 7.33333333333333e+00, dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513775e+00, dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513775e+00, dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999945e-01, dfll[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999945e-01, dfll[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.11022302462516e-16, dfll[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.11022302462516e-16, dfll[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.33333333333333e+00, dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00, dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00, dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dflu[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dflu[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dflu[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dflu[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.33333333333333e+00, dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00, dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00, dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfrl[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfrl[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dfrl[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dfrl[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 7.33333333333333e+00, dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513775e+00, dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513775e+00, dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfru[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfru[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfru[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfru[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfru[7], 1e-12) );
    } else if (test_func_op==1) {
      TEST_CHECK( gkyl_compare(-1.33333333333333e+00, dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.93889390390723e-17, dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.84900179459750e-01, dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.77555756156289e-17, dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfll[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfll[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.44265186329548e-01, dfll[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfll[7], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.33333333333333e+00, dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.93889390390723e-17, dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.84900179459750e-01, dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17, dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dflu[4], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17, dflu[5], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.44265186329548e-01, dflu[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dflu[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.33333333333333e+00, dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.92450089729875e+00, dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.77555756156289e-17, dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-5.96284793999944e-01, dfrl[4], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.38777878078145e-16, dfrl[5], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.44265186329548e-01, dfrl[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfrl[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.33333333333333e+00, dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.93889390390723e-17, dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.92450089729875e+00, dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17, dfru[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-5.96284793999944e-01, dfru[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 8.32667268468867e-17, dfru[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.44265186329548e-01, dfru[6], 1e-12) );
      TEST_CHECK( gkyl_compare(-6.93889390390723e-17, dfru[7], 1e-12) );
    }
  }

  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(distf);
}

void evalFunc_3x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = x*x + x*y + y*y + x*z + y*z + z*z;
}

void evalFunc_3x_trig(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = sin((2.*M_PI/(4.))*x) + sin((2.*M_PI/(4.))*x)*cos((2.*M_PI/(4.))*y) + cos((2.*M_PI/(4.))*y)
           + sin((2.*M_PI/(4.))*x)*cos((2.*M_PI/(4.))*z)
           + cos((2.*M_PI/(4.))*y)*cos((2.*M_PI/(4.))*z) + cos((2.*M_PI/(4.))*z);
}

void test_3x(int poly_order, int test_func_op)
{
  double lower[] = {-2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evup;
  if (test_func_op==0)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_3x_quad, NULL);
  else if (test_func_op==1)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_3x_trig, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_eval_on_nodes_advance(evup, 0.0, &arr_range, distf);

  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, "ctest_eval_on_nodes_distf.gkyl");

  double *dfllf = gkyl_array_fetch(distf, 0);  // left, low, front cell
  double *dfllb = gkyl_array_fetch(distf, 1);  // left, low, back cell
  double *dfluf = gkyl_array_fetch(distf, 2);  // left, up, front cell
  double *dflub = gkyl_array_fetch(distf, 3);  // left, up, back cell
  double *dfrlf = gkyl_array_fetch(distf, 4);  // right, low, front cell
  double *dfrlb = gkyl_array_fetch(distf, 5);  // right, low, back cell
  double *dfruf = gkyl_array_fetch(distf, 6);  // right, up, front cell
  double *dfrub = gkyl_array_fetch(distf, 7);  // right, up, back cell
  if (poly_order == 1) {
    if (test_func_op==0) {
      double fref[] = {
         2.54558441227157e+01, -6.53197264742181e+00, -6.53197264742181e+00, -6.53197264742181e+00, 9.42809041582064e-01, 9.42809041582062e-01,
         9.42809041582062e-01, -5.55111512312578e-17, 1.41421356237310e+01, -3.26598632371091e+00, -3.26598632371090e+00, 7.77156117237610e-16,
         9.42809041582064e-01, 9.42809041582063e-01, 9.42809041582063e-01, -5.55111512312578e-17, 1.41421356237309e+01, -3.26598632371091e+00,
         3.33066907387547e-16, -3.26598632371090e+00, 9.42809041582064e-01, 9.42809041582063e-01, 9.42809041582064e-01, 5.55111512312578e-17,
         1.41421356237310e+01, 0.00000000000000e+00, 3.26598632371091e+00, 3.26598632371091e+00, 9.42809041582064e-01, 9.42809041582064e-01,
         9.42809041582064e-01, 0.00000000000000e+00, 1.41421356237309e+01, 2.22044604925031e-16, -3.26598632371090e+00, -3.26598632371090e+00,
         9.42809041582064e-01, 9.42809041582064e-01, 9.42809041582063e-01, 0.00000000000000e+00, 1.41421356237310e+01, 3.26598632371091e+00,
         8.88178419700125e-16, 3.26598632371091e+00, 9.42809041582064e-01, 9.42809041582064e-01, 9.42809041582064e-01, 0.00000000000000e+00,
         1.41421356237310e+01, 3.26598632371091e+00, 3.26598632371091e+00, 8.88178419700125e-16, 9.42809041582064e-01, 9.42809041582064e-01,
         9.42809041582064e-01, 0.00000000000000e+00, 2.54558441227157e+01, 6.53197264742181e+00, 6.53197264742181e+00, 6.53197264742181e+00,
         9.42809041582064e-01, 9.42809041582065e-01, 9.42809041582065e-01, 0.00000000000000e+00,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    } else if (test_func_op==1) {
      double fref[] = {
         8.88178419700125e-16, 2.22044604925031e-16, 1.63299316185545e+00, 1.63299316185545e+00, 5.55111512312578e-17, 5.55111512312578e-17,
         9.42809041582064e-01, 2.77555756156289e-17, -1.22124532708767e-15, 8.32667268468867e-17, 1.63299316185545e+00, -1.63299316185545e+00,
         6.93889390390723e-17, -6.93889390390723e-17, -9.42809041582063e-01, 2.77555756156289e-17, 7.77156117237610e-16, -2.22044604925031e-16,
         -1.63299316185545e+00, 1.63299316185545e+00, -5.55111512312578e-17, 5.55111512312578e-17, -9.42809041582064e-01, -5.55111512312578e-17,
         -1.16573417585641e-15, 8.32667268468867e-17, -1.63299316185545e+00, -1.63299316185545e+00, -6.93889390390723e-17, -6.93889390390723e-17,
         9.42809041582063e-01, 1.38777878078145e-17, 6.66133814775094e-16, 2.22044604925031e-16, 1.63299316185545e+00, 1.63299316185545e+00,
         1.11022302462516e-16, 1.11022302462516e-16, 9.42809041582064e-01, 0.00000000000000e+00, -7.77156117237610e-16, 1.66533453693773e-16,
         1.63299316185545e+00, -1.63299316185545e+00, 8.32667268468867e-17, -8.32667268468867e-17, -9.42809041582063e-01, 2.77555756156289e-17,
         9.99200722162641e-16, -2.22044604925031e-16, -1.63299316185545e+00, 1.63299316185545e+00, -1.11022302462516e-16, 1.11022302462516e-16,
         -9.42809041582064e-01, -8.32667268468867e-17, -1.27675647831893e-15, 1.66533453693773e-16, -1.63299316185545e+00, -1.63299316185545e+00,
         -8.32667268468867e-17, -8.32667268468867e-17, 9.42809041582063e-01, 1.38777878078145e-17,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    }
  } else if (poly_order == 2) {
    if (test_func_op==0) {
      double fref[] = {
         1.97989898732233e+01, -6.53197264742181e+00, -6.53197264742181e+00, -6.53197264742181e+00, 9.42809041582064e-01, 9.42809041582063e-01,
         9.42809041582064e-01, 8.43274042711568e-01, 8.43274042711567e-01, 8.43274042711568e-01, 0.00000000000000e+00, 1.52655665885959e-16,
         5.55111512312578e-17, 1.52655665885959e-16, 2.77555756156289e-16, 2.22044604925031e-16, 2.22044604925031e-16, 0.00000000000000e+00,
         0.00000000000000e+00, 0.00000000000000e+00, 8.48528137423857e+00, -3.26598632371090e+00, -3.26598632371090e+00, 2.22044604925031e-16,
         9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711568e-01, 8.43274042711568e-01, 8.43274042711568e-01,
         0.00000000000000e+00, 5.55111512312578e-17, 5.55111512312578e-17, -5.55111512312578e-17, -5.55111512312578e-17, 1.66533453693773e-16,
         1.66533453693773e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 8.48528137423858e+00, -3.26598632371090e+00,
         2.22044604925031e-16, -3.26598632371090e+00, 9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711568e-01,
         8.43274042711568e-01, 8.43274042711568e-01, 0.00000000000000e+00, -5.55111512312578e-17, 1.66533453693773e-16, 5.55111512312578e-17,
         1.66533453693773e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00,
         8.48528137423857e+00, 0.00000000000000e+00, 3.26598632371090e+00, 3.26598632371090e+00, 9.42809041582063e-01, 9.42809041582063e-01,
         9.42809041582063e-01, 8.43274042711568e-01, 8.43274042711568e-01, 8.43274042711568e-01, 0.00000000000000e+00, -1.11022302462516e-16,
         1.11022302462516e-16, -1.11022302462516e-16, -1.11022302462516e-16, 0.00000000000000e+00, 1.11022302462516e-16, 0.00000000000000e+00,
         0.00000000000000e+00, 0.00000000000000e+00, 8.48528137423857e+00, 0.00000000000000e+00, -3.26598632371090e+00, -3.26598632371090e+00,
         9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711568e-01, 8.43274042711568e-01, 8.43274042711568e-01,
         0.00000000000000e+00, 1.66533453693773e-16, -5.55111512312578e-17, 1.66533453693773e-16, 1.66533453693773e-16, -5.55111512312578e-17,
         1.11022302462516e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 8.48528137423857e+00, 3.26598632371090e+00,
         0.00000000000000e+00, 3.26598632371090e+00, 9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711568e-01,
         8.43274042711568e-01, 8.43274042711568e-01, 0.00000000000000e+00, 1.11022302462516e-16, -2.22044604925031e-16, 0.00000000000000e+00,
         -2.22044604925031e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00,
         8.48528137423857e+00, 3.26598632371090e+00, 3.26598632371090e+00, 4.44089209850063e-16, 9.42809041582063e-01, 9.42809041582063e-01,
         9.42809041582063e-01, 8.43274042711568e-01, 8.43274042711568e-01, 8.43274042711568e-01, 0.00000000000000e+00, 0.00000000000000e+00,
         -1.11022302462516e-16, 1.11022302462516e-16, 0.00000000000000e+00, -1.11022302462516e-16, -1.11022302462516e-16, 0.00000000000000e+00,
         0.00000000000000e+00, 0.00000000000000e+00, 1.97989898732233e+01, 6.53197264742181e+00, 6.53197264742181e+00, 6.53197264742181e+00,
         9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711567e-01, 8.43274042711567e-01, 8.43274042711568e-01,
         0.00000000000000e+00, -2.22044604925031e-16, -2.22044604925031e-16, -2.22044604925031e-16, -2.22044604925031e-16, -2.22044604925031e-16,
         -2.22044604925031e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    } else if (test_func_op==1) {
      double fref[] = {
        -1.88561808316413e+00, 8.32667268468867e-17, 5.44331053951817e-01, 5.44331053951817e-01, 1.38777878078145e-17, 6.93889390390723e-17,
        9.42809041582063e-01, 8.43274042711568e-01, -5.55111512312578e-17, -5.55111512312578e-17, 0.00000000000000e+00, 4.86864495560148e-01,
        -2.77555756156289e-17, 4.86864495560148e-01, -8.32667268468867e-17, 2.77555756156289e-17, -8.32667268468867e-17, -1.38777878078145e-17,
        -1.38777878078145e-17, 2.77555756156289e-17, -1.88561808316413e+00, 5.55111512312578e-17, 5.44331053951817e-01, -5.44331053951817e-01,
        3.46944695195361e-17, -8.32667268468867e-17, -9.42809041582063e-01, 8.43274042711568e-01, -2.77555756156289e-17, -1.38777878078145e-16,
        -1.38777878078145e-17, 4.86864495560148e-01, -4.16333634234434e-17, -4.86864495560147e-01, 1.11022302462516e-16, 2.08166817117217e-17,
        -2.77555756156289e-17, 0.00000000000000e+00, 0.00000000000000e+00, 2.08166817117217e-17, -1.88561808316413e+00, 5.55111512312578e-17,
        -5.44331053951817e-01, 5.44331053951817e-01, -1.38777878078145e-17, 8.32667268468867e-17, -9.42809041582063e-01, 8.43274042711568e-01,
        -5.55111512312578e-17, -8.32667268468867e-17, 0.00000000000000e+00, -4.86864495560148e-01, -1.38777878078145e-17, 4.86864495560148e-01,
        -9.71445146547012e-17, 2.77555756156289e-17, 8.32667268468867e-17, 1.38777878078145e-17, 0.00000000000000e+00, -1.38777878078145e-17,
        -1.88561808316413e+00, 6.93889390390723e-17, -5.44331053951817e-01, -5.44331053951817e-01, -6.93889390390723e-18, -9.02056207507940e-17,
        9.42809041582063e-01, 8.43274042711568e-01, -4.16333634234434e-17, -6.93889390390723e-17, 1.38777878078145e-17, -4.86864495560147e-01,
        -2.08166817117217e-17, -4.86864495560147e-01, 1.31838984174237e-16, 2.08166817117217e-17, 1.17961196366423e-16, 6.93889390390723e-18,
        6.93889390390723e-18, -2.08166817117217e-17, 1.88561808316413e+00, 8.32667268468867e-17, 2.72165526975909e+00, 2.72165526975909e+00,
        4.16333634234434e-17, 6.93889390390723e-17, 9.42809041582063e-01, -8.43274042711568e-01, -1.66533453693773e-16, -5.55111512312578e-17,
        2.77555756156289e-17, -4.86864495560148e-01, 0.00000000000000e+00, -4.86864495560148e-01, -8.32667268468867e-17, 0.00000000000000e+00,
        -2.77555756156289e-17, 1.38777878078145e-17, 0.00000000000000e+00, 1.38777878078145e-17, 1.88561808316413e+00, 1.11022302462516e-16,
        2.72165526975909e+00, -2.72165526975909e+00, 6.24500451351651e-17, -8.32667268468867e-17, -9.42809041582063e-01, -8.43274042711568e-01,
        -1.94289029309402e-16, -1.94289029309402e-16, -1.38777878078145e-17, -4.86864495560148e-01, -1.38777878078145e-17, 4.86864495560148e-01,
        9.71445146547012e-17, -1.38777878078145e-17, -4.16333634234434e-17, 0.00000000000000e+00, 0.00000000000000e+00, 6.93889390390723e-18,
        1.88561808316413e+00, 1.11022302462516e-16, -2.72165526975909e+00, 2.72165526975909e+00, -4.16333634234434e-17, 8.32667268468867e-17,
        -9.42809041582063e-01, -8.43274042711568e-01, -1.38777878078145e-16, 2.77555756156289e-17, -2.77555756156289e-17, 4.86864495560148e-01,
        -2.77555756156289e-17, -4.86864495560148e-01, -1.24900090270330e-16, 0.00000000000000e+00, 5.55111512312578e-17, -1.38777878078145e-17,
        0.00000000000000e+00, 0.00000000000000e+00, 1.88561808316413e+00, 6.93889390390723e-17, -2.72165526975909e+00, -2.72165526975909e+00,
        -3.46944695195361e-17, -9.02056207507940e-17, 9.42809041582063e-01, -8.43274042711568e-01, -1.66533453693773e-16, -1.66533453693773e-16,
        1.38777878078145e-17, 4.86864495560148e-01, -1.38777878078145e-17, 4.86864495560148e-01, 8.32667268468867e-17, -1.38777878078145e-17,
        4.16333634234434e-17, -6.93889390390723e-18, -6.93889390390723e-18, -6.93889390390723e-18,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    }
  }

  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(distf);
}

void test_1x_p1_quad() { test_1x(1, 0); };
void test_1x_p1_trig() { test_1x(1, 1); };
void test_2x_p1_quad() { test_2x(1, 0); };
void test_2x_p1_trig() { test_2x(1, 1); };
void test_3x_p1_quad() { test_3x(1, 0); };
void test_3x_p1_trig() { test_3x(1, 1); };

void test_1x_p2_quad() { test_1x(2, 0); };
void test_1x_p2_trig() { test_1x(2, 1); };
void test_2x_p2_quad() { test_2x(2, 0); };
void test_2x_p2_trig() { test_2x(2, 1); };
void test_3x_p2_quad() { test_3x(2, 0); };
void test_3x_p2_trig() { test_3x(2, 1); };

TEST_LIST = {
  { "test_1x_p1_quad", test_1x_p1_quad },
  { "test_1x_p1_trig", test_1x_p1_trig },
  { "test_2x_p1_quad", test_2x_p1_quad },
  { "test_2x_p1_trig", test_2x_p1_trig },
  { "test_3x_p1_quad", test_3x_p1_quad },
  { "test_3x_p1_trig", test_3x_p1_trig },
//
  { "test_1x_p2_quad", test_1x_p2_quad },
  { "test_1x_p2_trig", test_1x_p2_trig },
  { "test_2x_p2_quad", test_2x_p2_quad },
  { "test_2x_p2_trig", test_2x_p2_trig },
  { "test_3x_p2_quad", test_3x_p2_quad },
  { "test_3x_p2_trig", test_3x_p2_trig },
  { NULL, NULL },
};
