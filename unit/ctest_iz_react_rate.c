#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_iz.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_const.h>
#include <stdio.h>

// Global variables
double echarge = GKYL_ELEMENTARY_CHARGE;
double emass = GKYL_ELECTRON_MASS;

void eval_m0(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_m2(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40*echarge/emass*1.0e19;  //fabs(x);
}

void
test_iz_react_rate()
{

  int poly_order = 1;
  int cdim = 1, vdim = 1;
  int pdim = cdim + vdim; 
  double lower[] = {-2.0,-1.0}, upper[] = {2.0,1.0};
  int ghost[] = {0, 0};
  int cells[] = {1,2};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange, &phaseRange);

  // basis functions
  struct gkyl_basis phaseBasis, basis; // phase-space, conf-space basis

  /* Force hybrid basis (p=2 in velocity space). */
  gkyl_cart_modal_hybrid(&phaseBasis, cdim, vdim);
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  // projection updater for moments
  gkyl_proj_on_basis *projM0 = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projM2 = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_m2, NULL);

   // create moment arrays
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *m2 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_elc = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *cflRate = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
  struct gkyl_array *vtSqIz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *coefIz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);

  // project moments on basis
  gkyl_proj_on_basis_advance(projM0, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projM2, 0.0, &confRange, m2);

  gkyl_array_set_offset(moms_neut, 1.0, m0, 0);
  gkyl_array_set_offset(moms_elc, 1.0, m0, 0);
  gkyl_array_set_offset(moms_neut, 1.0, m2, 2*basis.num_basis);
  gkyl_array_set_offset(moms_elc, 1.0, m2, 2*basis.num_basis);

  gkyl_grid_sub_array_write(&confGrid, &confRange, moms_neut, "ctest_moms_neut_1x.gkyl");
  gkyl_grid_sub_array_write(&confGrid, &confRange, moms_elc, "ctest_moms_elc_1x.gkyl");

  struct gkyl_dg_iz *reactRate = gkyl_dg_iz_new(&basis, &phaseBasis, &confRange, &phaseRange,
  						echarge, emass, GKYL_H, true, false);

  gkyl_dg_iz_react_rate(reactRate, moms_elc, moms_neut, vtSqIz, cflRate, coefIz);
  gkyl_grid_sub_array_write(&confGrid, &confRange, vtSqIz, "ctest_vtSqIz_1x.gkyl");
  gkyl_grid_sub_array_write(&confGrid, &confRange, coefIz, "ctest_react_rate_1x.gkyl");
    
  // left cell
  double *cl_vt = gkyl_array_fetch(vtSqIz, 0);
  TEST_CHECK( gkyl_compare(3.8470971703792085e+12, cl_vt[0], 1e-12) );
  TEST_CHECK( gkyl_compare(0.0, cl_vt[1], 1e-12) );

  double *cl_coef = gkyl_array_fetch(coefIz, 0);
  TEST_CHECK( gkyl_compare(3.362239235468358e-14, cl_coef[0], 1e-16) );
  TEST_CHECK( gkyl_compare(0.0, cl_coef[1], 1e-16) );
  
  //gkyl_dg_iz_release(reactRate);
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
