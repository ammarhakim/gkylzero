#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_iz.h>
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
  int cdim = 1, vdim = 1;
  int pdim = cdim + vdim; 
  double lower[] = {-2.0,-1.0}, upper[] = {2.0,1.0};
  int ghost[] = {0, 0};
  int cells[] = {16,8};

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
  gkyl_proj_on_basis *projDens = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_dens, NULL);
  gkyl_proj_on_basis *projVtSq = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_vtSq, NULL);

  // create moment arrays
  struct gkyl_array *m0Neut = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *vtSqElc = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *vtSqNeut = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *cflRate = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *coefIz = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projDens, 0.0, &confRange, m0Neut);
  gkyl_proj_on_basis_advance(projVtSq, 0.0, &confRange, vtSqElc);
  gkyl_proj_on_basis_advance(projVtSq, 0.0, &confRange, vtSqNeut);

  gkyl_grid_sub_array_write(&confGrid, &confRange, m0Neut, "ctest_iz_m0_1x.gkyl");
  gkyl_grid_sub_array_write(&confGrid, &confRange, vtSqElc, "ctest_iz_vtSqElc_1x.gkyl");

  double echarge = 1.602176634e-19;
  double emass = 9.109383701528e-31;
  
  struct gkyl_dg_iz *reactRate = gkyl_dg_iz_new(&basis, echarge, emass, GKYL_H, false);

  //gkyl_dg_iz_temp(reactRate, &confRange, vtSqElc, coefIz); 

  gkyl_dg_iz_react_rate(reactRate, &confRange, &phaseRange, m0Neut, vtSqNeut, vtSqElc, cflRate, coefIz);
  gkyl_grid_sub_array_write(&confGrid, &confRange, coefIz, "ctest_react_rate_1x.gkyl");

  // left cell
  double *cl = gkyl_array_fetch(coefIz, 0);
  TEST_CHECK( gkyl_compare(4.5659992900909018e+05, cl[0], 1e-12) );
  TEST_CHECK( gkyl_compare(-4.4853312014181618e+01, cl[1], 1e-12) );

  gkyl_dg_iz_release(reactRate);
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
