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
#include <gkyl_dg_vlasov_lbo.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_array_rio.h>
#include <math.h>

static inline double sq(double x) { return x*x; }

void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vx = xn[1];
  if(vx>-1.0 && vx<1.0) {
    fout[0] = 1.0;
  } else {
    fout[0] = 0.0;
  }
}

void
test_1x1v_p2()
{
  int poly_order = 2;
  double lower[] = {0.0, -4.0}, upper[] = {1.0, 4.0};
  int cells[] = {4, 24};
  int ghost[] = {0, 0};
  int pdim = sizeof(lower)/sizeof(lower[0]);
  int vdim = 1, cdim = 1;

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // basis functions
  struct gkyl_basis basis, confBasis;
  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDist = gkyl_proj_on_basis_new(&phaseGrid, &basis, poly_order+1, 1, evalFunc, NULL);

  // create array range: no ghost-cells in velocity space
  struct gkyl_range arr_range;
  gkyl_range_init_from_shape(&arr_range, pdim, cells);

  // create distribution function
  //struct gkyl_array *dist = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange.volume);
  struct gkyl_array *dist = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDist, 0.0, &arr_range, dist);

  int up_dirs[] = {0, 1};
  int zero_flux_flags[] = {0, 1};

  struct gkyl_dg_eqn *eqn;
  gkyl_hyper_dg *slvr;

  eqn = gkyl_dg_vlasov_lbo_new(&confBasis, &basis);
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1);

  double nuSum = 1.0;
  double nuUSum[vdim];
  double nuVtSqSum[vdim];
  
  int nnu = vdim*confBasis.num_basis;
  for(int i=0; i< nnu; i++) {
   nuUSum[i] = 0.0;
   if (i==0) {nuVtSqSum[i] = 1.0;} else {nuVtSqSum[i] = 0.0;}
  }

  struct gkyl_array *cflrate, *rhs, *fIn;
  double *cfl_ptr;
  cflrate = gkyl_array_new(GKYL_DOUBLE, 1, phaseRange_ext.volume);
  rhs = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ext.volume);
  fIn = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, phaseRange_ext.volume);
  gkyl_array_copy_range(fIn, dist, phaseRange);
  cfl_ptr = gkyl_malloc(sizeof(double));

  // Write the initial distribution array to file.
  const char *fmt = "%s-%s_%d.gkyl";
  char name[] = "distf";
  char momName[] = "elc";
  int frame = 0;
  int sz = snprintf(0, 0, fmt, name, momName, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, name, momName, frame);
  gkyl_grid_sub_array_write(&phaseGrid, &phaseRange, dist, fileNm);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_vlasov_lbo_set_nuSum(eqn, nuSum);
    gkyl_vlasov_lbo_set_nuUSum(eqn, nuUSum);
    gkyl_vlasov_lbo_set_nuVtSqSum(eqn, nuVtSqSum);
    gkyl_hyper_dg_advance(slvr, phaseRange, dist, cflrate, rhs);
    gkyl_array_reduce(cfl_ptr, cflrate, GKYL_MAX);
    gkyl_array_copy_range(fIn, rhs, phaseRange);
  }

  // Write the final distribution array to file.
  frame = 1;
  sz = snprintf(0, 0, fmt, name, momName, frame);
  fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, name, momName, frame);
  gkyl_grid_sub_array_write(&phaseGrid, &phaseRange, rhs, fileNm);

  // release memory for moment data object
  gkyl_proj_on_basis_release(projDist);
  gkyl_array_release(dist);
  gkyl_array_release(rhs);
  gkyl_array_release(fIn);
  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "test_1x1v_p2", test_1x1v_p2 },
  { NULL, NULL },
};
