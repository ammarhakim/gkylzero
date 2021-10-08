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

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
test_1x2v_p2()
{
  // initialize grid and ranges
  int cdim = 1, vdim = 2;
  int pdim = cdim+vdim;

  int cells[] = {24, 12, 12};
  int ghost[] = {1, 0, 0};
  double lower[] = {0., -1., -1.};
  double upper[] = {1., 1., 1.};

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

  // initialize eqn
  struct gkyl_dg_eqn *eqn;
  eqn = gkyl_dg_vlasov_lbo_new(&confBasis, &basis, &confRange);

  // initialize hyper_dg slvr
  int up_dirs[] = {1, 2};
  int zero_flux_flags[] = {1, 1};

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, vdim, up_dirs, zero_flux_flags, 1);

  struct gkyl_array *cflrate, *rhs, *fin, *nuSum, *nuUSum, *nuVtSqSum;
  double *cfl_ptr;
  cflrate = mkarr(1, phaseRange_ext.volume);
  rhs = mkarr(basis.num_basis, phaseRange_ext.volume);
  fin = mkarr(basis.num_basis, phaseRange_ext.volume);
  nuSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  nuUSum = mkarr(vdim*confBasis.num_basis, confRange_ext.volume);
  nuVtSqSum = mkarr(confBasis.num_basis, confRange_ext.volume);
  
  cfl_ptr = gkyl_malloc(sizeof(double));

  // set initial condition
  int nf = phaseRange_ext.volume;
  double *fin_d;
  fin_d = fin->data;
  for(int i=0; i<nf; i++) {
    for(int j=0; j<basis.num_basis; j++) {
      fin_d[i*basis.num_basis + j] = (double)(2*i+(11 + j)% nf) / nf  * ((i%2 == 0) ? 1 : -1);
      //printf("Values: %d, %d, \%f\n", nf, i, fin_d[i*basis.num_basis + j]);
    }
  }

  int nem = confRange_ext.volume*confBasis.num_basis;
  gkyl_array_clear(nuSum, 1.0);
  gkyl_array_clear(nuUSum, 0.0);
  gkyl_array_clear(nuVtSqSum, 1.0);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_vlasov_lbo_set_nuSum(eqn, nuSum);
    gkyl_vlasov_lbo_set_nuUSum(eqn, nuUSum);
    gkyl_vlasov_lbo_set_nuVtSqSum(eqn, nuVtSqSum);
    gkyl_hyper_dg_advance(slvr, phaseRange, fin, cflrate, rhs);

    gkyl_array_reduce(cfl_ptr, cflrate, GKYL_MAX);
  }

  printf("CFL: %f\n", cfl_ptr[0]);
  //TEST_CHECK( gkyl_compare_double(cfl_ptr[0], 2.5178875733842702e+01, 1e-12) );

  // get linear index of first non-ghost cell
  int idx[] = {0, 0, 0, 0, 0};
  int linl = gkyl_range_idx(&phaseRange, idx);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d;
  int i = 0;
  while(val==0) {
    rhs_d = gkyl_array_fetch(rhs, i);
    val = rhs_d[0];
    if(val==0) i++;
    }
  TEST_CHECK(i == linl);

  // get linear index of some other cell
  int idx2[] = {5, 2, 4};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs, linl2);
  for(int i=0; i<basis.num_basis; i++) {
    printf("[%d, %d, %d, %d] = %f\n", 1, 2, 2, i, rhs_d[i]);
  }

  // release memory for moment data object
  //gkyl_free(cfl_ptr);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(fin);
  gkyl_array_release(nuSum);
  gkyl_array_release(nuUSum);
  gkyl_array_release(nuVtSqSum);
  gkyl_dg_eqn_release(eqn);
  gkyl_hyper_dg_release(slvr);
}

TEST_LIST = {
  { "test_1x2v_p2", test_1x2v_p2 },
  { NULL, NULL },
};
