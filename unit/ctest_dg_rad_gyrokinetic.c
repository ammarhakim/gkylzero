// A test for the vlasov LBO collision updater
//
#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_updater_rad_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_const.h>
#include <gkyl_array_rio.h>
#include <math.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// allocate cu_dev array
static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

void ni_prof(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vx = xn[1];
  double vy  = xn[2];
  fout[0] = 1.0;
}

// May need some velocity transformation of vperp
void maxwellian1x2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vx = xn[1];
  double vy = xn[2];
  double *te = ctx;
  fout[0] = 1.0/(2*M_PI)*exp(-(pow(vx, 2) + pow(vy, 2))/(2*te[0]));
}

void maxwellian1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vx = xn[1];
  double *te = ctx;
  fout[0] = 1.0/sqrt(2*M_PI)*exp(-(pow(vx, 2))/(2*te[0]));
}

void
test_1x1v_p2()
{

  // rad specific variables
  //double *b = {1.0, 1.0};
  //double *data = {0.16, 8000.1, 0.9, -3.9, 3.1}; // H fit params
  struct gkyl_array *fit_params = gkyl_array_new(GKYL_DOUBLE, 1, 5);
  struct gkyl_array *bmag = gkyl_array_new(GKYL_DOUBLE, 1, 2);
  double *te = (double*) malloc(sizeof(double)*1);
  double *data = fit_params->data;
  double *b = bmag->data;
  te[0]=30;
  b[0]=1;
  b[1]=1;
  data[0]=0.16;
  data[1]=8000.1;
  data[2]=0.9;
  data[3]=-3.9;
  data[4]=3.1;
  double *temp = (double *)gkyl_array_fetch(bmag,0);
  double *temp2 = (double *)gkyl_array_fetch(fit_params,0);
  printf("\n");
  for (int i=0; i<2;i++) {
    printf("b[i]=%f, bmag[%d]=%f ",i,b[i],i,temp[i]);
    printf("data[i]=%f, fit_params[%d]=%f\n",i,data[i],i,temp2[i]);
    }


   // initialize grid and ranges
  int cdim = 1, vdim = 1;
  int pdim = cdim+vdim;

  int cells[] = {2, 4};
  int ghost[] = {0, 0};
  double lower[] = {0., -1.};
  double upper[] = {1., 5.e6};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize hyper_dg slvr
  int up_dirs[] = {1};
  int zero_flux_flags[] = {1};

  gkyl_dg_updater_collisions *slvr;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;
  printf("Before dg updater\n");
  //  slvr = gkyl_dg_updater_lbo_vlasov_new(&phaseGrid, &confBasis, &basis, &confRange, model_id, false);
 
  //Initilize vnu and vsqnu
  struct gkyl_array *vnu = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
  struct gkyl_array *vsqnu = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
  
  printf("int %d\n",phaseRange.ndim);
  slvr = gkyl_dg_updater_rad_gyrokinetic_new(&phaseGrid, &confBasis, &basis, &confRange,  &phaseRange, bmag, fit_params, vnu, vsqnu, false);
   // printf("Created vlasov updater\n");
  
  
  printf("After dg updater (ctest)\n");
  gkyl_proj_on_basis *projF = gkyl_proj_on_basis_new(&phaseGrid, &basis, poly_order+1, 1, maxwellian1x1v, te);
  gkyl_proj_on_basis *projNi = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, ni_prof, NULL);
  printf("After projection (ctest)\n");
  struct gkyl_array *cflrate, *rhs, *fin, *nI;
  cflrate = mkarr(1, phaseRange_ext.volume);
  rhs = mkarr(basis.num_basis, phaseRange_ext.volume);
  fin = mkarr(basis.num_basis, phaseRange_ext.volume);
  nI = mkarr(confBasis.num_basis, confRange_ext.volume);
  //nuUSum = mkarr(vdim*confBasis.num_basis, confRange_ext.volume);
  //nuVtSqSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  //nuPrimMomsSum = mkarr((vdim+1)*confBasis.num_basis, confRange_ext.volume);

  gkyl_proj_on_basis_advance(projF, 0.0, &phaseRange_ext, fin);
  gkyl_proj_on_basis_advance(projNi, 0.0, &confRange_ext, nI);
  gkyl_proj_on_basis_release(projF);
  gkyl_proj_on_basis_release(projNi);
  printf("After projection release (ctest)\n");

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_dg_updater_rad_gyrokinetic_advance(slvr, &phaseRange, nI, vnu, vsqnu, fin, cflrate, rhs);
    printf("After updater advance, n=%i\n",n);
  }

  // Take 2nd moment of f to find energy
  //struct gkyl_mom_type *m2 = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, GKYL_ELECTRON_MASS, "M2", false);
  
  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d1, *rhs_d2;
  int i = 0;
  while(val==0) {
    rhs_d1 = gkyl_array_fetch(rhs, i);
    val = rhs_d1[0];
    if(val==0) i++;
    }
  TEST_CHECK(i == linl);

  // get linear index of some other cell
  int idx1[] = {1, 1};
  int idx2[] = {2, 3};
  int linl1 = gkyl_range_idx(&phaseRange, idx1);
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d1 = gkyl_array_fetch(rhs, linl1);
  rhs_d2 = gkyl_array_fetch(rhs, linl2);
  TEST_CHECK( gkyl_compare_double(rhs_d1[0], -0.20752566554314, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[1], 0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[2], 0.52331159647886, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[3], -1.3061806047047, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[4], 0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[5], -0.59363481492976, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[6], -1.3061806047047, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d1[7], 1.5729635651345, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[0], 0.20752566554313, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[1], -0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[2], -0.021201659682202, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[3], 0.05635142305445, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[4], -0.49773088011319, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[5], -0.00090746168862665, 1e-10) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[6], 0.056351423054451, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d2[7], 0.015080609728081, 1e-12) );

  // release memory for moment data object
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(fin);
  gkyl_array_release(nI);
  gkyl_dg_updater_rad_gyrokinetic_release(slvr);
}


#ifdef GKYL_HAVE_CUDA

#endif

TEST_LIST = {
  { "test_1x1v_p2", test_1x1v_p2 },
  // { "test_1x2v_p2", test_1x2v_p2 },
  //  #ifdef GKYL_HAVE_CUDA
  // { "test_1x1v_p2_cu", test_1x1v_p2_cu },
  // { "test_1x2v_p2_cu", test_1x2v_p2_cu },
  // #endif
  { NULL, NULL },
};
