#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_iz_react_rate.h>
#include <gkyl_array_rio.h>
#include <stdio.h>

void eval_dens(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_vtSq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40*1.602e-19/9.1e-31*x*x; //fabs(x);
}

void
test_iz_react_rate()
{

  int poly_order = 1;
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projDens = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_dens, NULL);
  gkyl_proj_on_basis *projVtSq = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_vtSq, NULL);
  
  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create moment arrays
  struct gkyl_array *m0Elc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *m0Neut = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *vtSqElc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *vtSqNeut = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  struct gkyl_array *coefIz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projDens, 0.0, &arr_range, m0Elc);
  gkyl_proj_on_basis_advance(projDens, 0.0, &arr_range, m0Neut);
  gkyl_proj_on_basis_advance(projVtSq, 0.0, &arr_range, vtSqElc);
  gkyl_proj_on_basis_advance(projVtSq, 0.0, &arr_range, vtSqNeut);

  double echarge = 1.602176634e-19;
  double emass = 9.109383701528e-31;
  
  // Parameters for H
  double A = 0.291e-13;
  double E = 13.6*echarge/emass; //same units as vtSq
  double K = 0.39;
  double P = 0.;
  double X = 0.232;
  
  struct gkyl_iz_react_rate *reactRate = gkyl_iz_react_rate_new(&basis, emass, A, E, K, P, X);

  gkyl_iz_react_rate_coef(reactRate, arr_range, m0Elc, m0Neut, vtSqElc, vtSqNeut, coefIz);
  gkyl_grid_sub_array_write(&grid, &arr_range, coefIz, "ctest_react_rate_1x.gkyl");

  // left cell
  double *cl = gkyl_array_fetch(coefIz, 0);
  TEST_CHECK( gkyl_compare(4.5659992900909018e+05, cl[0], 1e-12) );
  TEST_CHECK( gkyl_compare(-4.4853312014181618e+01, cl[1], 1e-12) );

  gkyl_iz_react_rate_release(reactRate);
}

#ifdef GKYL_HAVE_CUDA


#endif

TEST_LIST = {
  { "iz_react_rate", test_iz_react_rate },
#ifdef GKYL_HAVE_CUDA
/*  { "cu_dg_gyrokinetic", test_cu_dg_gyrokinetic }, */
#endif  
  { NULL, NULL },
};
