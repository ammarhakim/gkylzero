#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_cx.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_array_ops.h>
#include <stdio.h>
#include <math.h>
#include <gkyl_const.h>

// Global variables
double echarge = GKYL_ELEMENTARY_CHARGE;
double imass = GKYL_PROTON_MASS; // hyrdogen ion

void eval_m0(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19;
}
void eval_ion_u(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = sqrt(40*echarge/imass); //fabs(x);
}
void eval_neut_u(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.; //fabs(x);
}
void eval_ion_vtSq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 40*echarge/imass; //fabs(x);
}
void eval_neut_vtSq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 10*echarge/imass; //fabs(x);
}
void eval_ones(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.;
}
  
void
test_cx_react_rate()
{

  int poly_order = 1;
  int cdim = 1, vdim = 1;
  int pdim = cdim + vdim; 
  double lower[] = {-2.0,-1.0}, upper[] = {2.0,1.0};
  int ghost[] = {0, 0};
  int cells[] = {1,1};

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
    poly_order+1, 1, eval_m0, NULL);
  gkyl_proj_on_basis *projIonU = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_ion_u, NULL);
  gkyl_proj_on_basis *projNeutU = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_neut_u, NULL);
  gkyl_proj_on_basis *projIonVtSq = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_ion_vtSq, NULL);
  gkyl_proj_on_basis *projNeutVtSq = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_neut_vtSq, NULL);
  gkyl_proj_on_basis *projOnes = gkyl_proj_on_basis_new(&confGrid, &basis,
    poly_order+1, 1, eval_ones, NULL);
 
  // create moment arrays
  struct gkyl_array *m0 = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *uIon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *uNeut = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *vtSqIon = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *vtSqNeut = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *moms_neut = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *moms_ion = gkyl_array_new(GKYL_DOUBLE, 3*basis.num_basis, confRange.volume);
  struct gkyl_array *m2_temp = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *b_hat = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);

  struct gkyl_array *coefCx = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, confRange.volume);
  struct gkyl_array *cflRate = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
  
  // project moments on basis
  gkyl_proj_on_basis_advance(projDens, 0.0, &confRange, m0);
  gkyl_proj_on_basis_advance(projIonU, 0.0, &confRange, uIon);
  gkyl_proj_on_basis_advance(projNeutU, 0.0, &confRange, uNeut);
  gkyl_proj_on_basis_advance(projIonVtSq, 0.0, &confRange, vtSqIon);
  gkyl_proj_on_basis_advance(projNeutVtSq, 0.0, &confRange, vtSqNeut);
  gkyl_proj_on_basis_advance(projOnes, 0.0, &confRange, b_hat);

  gkyl_dg_mul_op_range(basis, 0, m2_temp, 0, uIon, 0, uIon, &confRange);
  gkyl_array_accumulate_range(m2_temp, 1.0, vtSqIon, confRange);
  gkyl_dg_mul_op_range(basis, 0, m2_temp, 0, m2_temp, 0, m0, &confRange);

  gkyl_array_set_offset(moms_ion, 1.0, m0, 0);
  gkyl_array_set_offset(moms_ion, 1.0, uIon, 1*basis.num_basis);
  gkyl_dg_mul_op_range(basis, 1, moms_ion, 1, moms_ion, 0, m0, &confRange);
  gkyl_array_set_offset(moms_ion, 1.0, m2_temp, 2*basis.num_basis);

  gkyl_dg_mul_op_range(basis, 0, m2_temp, 0, uNeut, 0, uNeut, &confRange);
  gkyl_array_accumulate_range(m2_temp, 1.0, vtSqNeut, confRange);
  gkyl_dg_mul_op_range(basis, 0, m2_temp, 0, m2_temp, 0, m0, &confRange); 

  gkyl_array_set_offset(moms_neut, 1.0, m0, 0);
  gkyl_array_set_offset(moms_neut, 1.0, uNeut, 1*basis.num_basis);
  gkyl_dg_mul_op_range(basis, 1, moms_neut, 1, moms_neut, 0, m0, &confRange);
  gkyl_array_set_offset(moms_neut, 1.0, m2_temp, 2*basis.num_basis);

  gkyl_grid_sub_array_write(&confGrid, &confRange, moms_neut, "ctest_cx_moms_neut_1x.gkyl");
  gkyl_grid_sub_array_write(&confGrid, &confRange, moms_ion, "ctest_cx_moms_ion_1x.gkyl");
 
  struct gkyl_dg_cx *reactRate = gkyl_dg_cx_new(&confGrid, &basis, &phaseBasis, &confRange, &phaseRange, imass, GKYL_CX_H, false, false);

  gkyl_dg_cx_react_rate(reactRate, moms_ion, moms_neut, b_hat, cflRate, coefCx);
  gkyl_grid_sub_array_write(&confGrid, &confRange, coefCx, "ctest_cx_react_rate_1x.gkyl");

  // Calculate reaction rate analytically to compare
  double m0_n = 1e19;
  double u_i = sqrt(40*echarge/imass);
  double u_n = 0.0;
  double vt2_i = 40*echarge/imass;
  double vt2_n = 10*echarge/imass;

  double v_cx_sq = 4./GKYL_PI*(vt2_i + vt2_n) + pow((u_i - u_n),2.0);
  double v_cx = sqrt(v_cx_sq); 

  double sig_cx = 1.12e-18 - 7.15e-20*log(v_cx);
  printf("\n%e", v_cx*sig_cx*sqrt(2));
  
  // left cell
  double *cl = gkyl_array_fetch(coefCx, 0);
  printf("\n%e", cl[0]); 
  TEST_CHECK( gkyl_compare(4.186509619868415e-14, cl[0], 1e-16) ); 
  //TEST_CHECK( gkyl_compare(1., cl[0], 1e-12) ); 
  
  gkyl_dg_cx_release(reactRate);
}

#ifdef GKYL_HAVE_CUDA


#endif

TEST_LIST = {
  { "cx_react_rate", test_cx_react_rate },
#ifdef GKYL_HAVE_CUDA
/*  { "cu_dg_gyrokinetic", test_cu_dg_gyrokinetic }, */
#endif  
  { NULL, NULL },
};
