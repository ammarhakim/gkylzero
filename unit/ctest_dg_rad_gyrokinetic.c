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
#include <gkyl_mom_calc.h>
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
mkarr_cu(long nc, long size, bool use_gpu)
{
  if (use_gpu) {
    struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
    return a;
  } else {
    struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
    return a;
  }
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
  fout[0] = 1.0/(2*M_PI*te[0])*exp(-(pow(vx, 2) + pow(vy, 2))/(2*te[0]));
}

void maxwellian1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vx = xn[1];
  double *te = ctx;
  fout[0] = 1.0/sqrt(2*M_PI*te[0])*exp(-(pow(vx, 2))*GKYL_ELECTRON_MASS/(2*te[0]));
}

void
test_1x1v_p2()
{

  // rad specific variables
  //double *b = {1.0, 1.0};
  //double *data = {0.16, 8000.1, 0.9, -3.9, 3.1}; // H fit params
  struct gkyl_array *fit_params = mkarr(1, 5);
  struct gkyl_array *bmag = mkarr(1, 2);
  double *te = (double*) malloc(sizeof(double)*1);
  double *data = fit_params->data;
  double *b = bmag->data;
  bool use_gpu = false;
  te[0]=30;
  b[0]=1;
  b[1]=1;
  data[0]=0.16;
  data[1]=8000.1;
  data[2]=0.9;
  data[3]=-5.9;
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

  int cells[] = {2, 5};
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
  int poly_order = 1;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize hyper_dg slvr
  int up_dirs[] = {1};
  int zero_flux_flags[] = {1};

  gkyl_dg_updater_collisions *slvr;
  enum gkyl_model_id model_id = GKYL_MODEL_DEFAULT;
 
  //Initilize vnu and vsqnu
  struct gkyl_array *vnu = mkarr(basis.num_basis, phaseRange.volume);
  struct gkyl_array *vsqnu = mkarr(basis.num_basis, phaseRange.volume);
  
  slvr = gkyl_dg_updater_rad_gyrokinetic_new(&phaseGrid, &confBasis, &basis, &confRange,  &phaseRange, bmag, fit_params, vnu, vsqnu, false);
  
  gkyl_proj_on_basis *projF = gkyl_proj_on_basis_new(&phaseGrid, &basis, poly_order+1, 1, maxwellian1x1v, te);
  gkyl_proj_on_basis *projNi = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, ni_prof, NULL);
  
  struct gkyl_array *cflrate, *rhs, *fin, *nI, *fmax;
  cflrate = mkarr(1, phaseRange_ext.volume);
  rhs = mkarr(basis.num_basis, phaseRange_ext.volume);
  printf("Size of rhs %d,%d\n",basis.num_basis, phaseRange_ext.volume);
  fin = mkarr(basis.num_basis, phaseRange_ext.volume);
  nI = mkarr(confBasis.num_basis, confRange_ext.volume);
  fmax = mkarr(basis.num_basis, phaseRange_ext.volume);
  //nuUSum = mkarr(vdim*confBasis.num_basis, confRange_ext.volume);
  //nuVtSqSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  //nuPrimMomsSum = mkarr((vdim+1)*confBasis.num_basis, confRange_ext.volume);

  gkyl_proj_on_basis_advance(projF, 0.0, &phaseRange_ext, fin);
  gkyl_proj_on_basis_advance(projNi, 0.0, &confRange_ext, nI);
  gkyl_array_copy(fmax, fin);
  gkyl_proj_on_basis_release(projF);
  gkyl_proj_on_basis_release(projNi);

  // run hyper_dg_advance
  int nrep = 1;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_dg_updater_rad_gyrokinetic_advance(slvr, &phaseRange, nI, vnu, vsqnu, fin, cflrate, rhs);
    printf("After updater advance, n=%i\n",n);
  }

  // Take 2nd moment of f to find energy
  struct gkyl_mom_type *m0 = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confRange, GKYL_ELECTRON_MASS, "M0", use_gpu);
  struct gkyl_mom_type *m2 = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confRange, GKYL_ELECTRON_MASS, "M2", use_gpu);
  printf("After moment new, cdim=%d\n",m2->cdim);
  gkyl_gyrokinetic_set_bmag(m0, bmag);
  gkyl_gyrokinetic_set_bmag(m2, bmag);
  struct gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&phaseGrid, m0, use_gpu);
  struct gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&phaseGrid, m2, use_gpu);
  printf("After moment calc new\n");
  struct gkyl_array *m0_ho = mkarr(confBasis.num_basis, confRange_ext.volume);
  struct gkyl_array *m0final = mkarr_cu(confBasis.num_basis, confRange_ext.volume, use_gpu);
  struct gkyl_array *m2_ho = mkarr(confBasis.num_basis, confRange_ext.volume);
  struct gkyl_array *m2final = mkarr_cu(confBasis.num_basis, confRange_ext.volume, use_gpu);
  printf("After array calculations\n");
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(m0calc, &phaseRange, &confRange, fmax, m0final);
    gkyl_mom_calc_advance_cu(m2calc, &phaseRange, &confRange, fmax, m2final);
  } else {
    gkyl_mom_calc_advance(m0calc, &phaseRange, &confRange, rhs, m0final);
    gkyl_mom_calc_advance(m2calc, &phaseRange, &confRange, rhs, m2final);
  }
  printf("After mom calc advance\n");
  double *m00 = gkyl_array_fetch(m0final, 0+ghost[0]);
  double *m20 = gkyl_array_fetch(m2final, 0+ghost[0]);
  double *m21 = gkyl_array_fetch(rhs, 0+ghost[0]);

  //  double cell_avg0 = m20[0]/pow(sqrt(2),cdim);
  double cell_avg0 = m20[0]/m00[0];

  double correct = 4.419192399328895e-32;
  for (int i=0; i<30; i++){
    printf("cell_avg=%e, correct energy=%e, density=%e, m2=%e\n",cell_avg0, correct, m00[0], m20[0]);
  }
  
  TEST_CHECK( gkyl_compare( cell_avg0, correct, 1e-12));
  TEST_CHECK( cell_avg0>0);


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
