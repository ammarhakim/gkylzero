#include <gkylzero.h>
#include <gkyl_dg_vlasov_priv.h>

#define TEST_NO_MAIN
#include <acutest.h>

extern "C" {
void test_vlasov_2x3v_p1_cu();
}

__global__
void ker_cu_hyper_dg_test(const gkyl_hyper_dg *slvr, int *nfail)
{
  *nfail = 0;
  
  GKYL_CU_CHECK( slvr->ndim == 5, nfail );
  GKYL_CU_CHECK( slvr->numBasis == 32, nfail );
  GKYL_CU_CHECK( slvr->num_up_dirs == 5, nfail );

  GKYL_CU_CHECK( slvr->grid.ndim == 5, nfail );

  GKYL_CU_CHECK( slvr->equation->num_equations == 1, nfail );

  // DO NOT DO THIS IN PRODUCTION! ONLY FOR TESTING
  struct dg_vlasov *vlasov = container_of(slvr->equation, struct dg_vlasov, eqn);

  GKYL_CU_CHECK( vlasov->cdim == 2, nfail );
  GKYL_CU_CHECK( vlasov->pdim == 5, nfail );
  GKYL_CU_CHECK( vlasov->conf_range.volume == 8*8, nfail );
}

static int
hyper_dg_test(const gkyl_hyper_dg *slvr)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_cu_hyper_dg_test<<<1,1>>>(slvr, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;  
}

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

void test_vlasov_2x3v_p1_cu()
{
  // initialize grid and ranges on host
  int cdim = 2, vdim = 3;
  int pdim = cdim+vdim;

  int cells[] = {8, 8, 8, 8, 8};
  int ghost[] = {1, 1, 0, 0, 0};
  double lower[] = {0., 0., -1., -1., -1.};
  double upper[] = {1., 1., 1., 1., 1.};

  struct gkyl_rect_grid confGrid;
  struct gkyl_range confRange, confRange_ext;
  gkyl_rect_grid_init(&confGrid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&confGrid, ghost, &confRange_ext, &confRange);

  struct gkyl_rect_grid phaseGrid;
  struct gkyl_range phaseRange, phaseRange_ext;
  gkyl_rect_grid_init(&phaseGrid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phaseGrid, ghost, &phaseRange_ext, &phaseRange);

  // initialize basis (note: basis has no device implementation)
  int poly_order = 1;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn on device
  struct gkyl_dg_eqn *eqn_cu;
  eqn_cu = gkyl_dg_vlasov_cu_dev_new(&confBasis, &basis, &confRange);

  // initialize hyper_dg slvr
  int up_dirs[] = {0, 1, 2, 3, 4};
  int zero_flux_flags[] = {0, 0, 1, 1, 1};
  double maxs_init[] = {0., 0., 0., 0., 0.};

  gkyl_hyper_dg *slvr_cu;
  slvr_cu = gkyl_hyper_dg_cu_dev_new(&phaseGrid, &basis, eqn_cu, pdim, up_dirs, zero_flux_flags, 1, maxs_init);

  // basic checks
  int nfail = hyper_dg_test(slvr_cu);

  TEST_CHECK( nfail == 0 );

  // initialize host arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  fin = mkarr(basis.numBasis, phaseRange_ext.volume);
  rhs = mkarr(basis.numBasis, phaseRange_ext.volume);
  cflrate = mkarr(1, phaseRange_ext.volume);
  qmem = mkarr(8*confBasis.numBasis, confRange_ext.volume);

  // set initial condition
  int nf = phaseRange_ext.volume*basis.numBasis;
  double *fin_d = (double*) fin->data;
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }

  int nem = confRange_ext.volume*confBasis.numBasis;
  double *qmem_d = (double*) qmem->data;
  for(int i=0; i< nem; i++) {
    qmem_d[i] = (double)(-i+27 % nem) / nem  * ((i%2 == 0) ? 1 : -1);
  }

  // initialize device arrays 
  struct gkyl_array *fin_cu = mkarr_cu(basis.numBasis, phaseRange_ext.volume);
  struct gkyl_array *rhs_cu = mkarr_cu(basis.numBasis, phaseRange_ext.volume);
  struct gkyl_array *cflrate_cu = mkarr_cu(1, phaseRange_ext.volume);
  struct gkyl_array *qmem_cu = mkarr_cu(8*confBasis.numBasis, confRange_ext.volume);
  struct gkyl_array *maxs_by_cell_cu = mkarr_cu(GKYL_MAX_DIM, phaseRange_ext.volume);

  // also set up maxs arrays on host and device
  double* maxs = (double*) gkyl_malloc(GKYL_MAX_DIM*sizeof(double));
  double* maxs_cu = (double*) gkyl_cu_malloc(GKYL_MAX_DIM*sizeof(double));
  for(int i=0; i<GKYL_MAX_DIM; i++) 
    maxs[i] = 0.;

  // copy initial conditions to device
  gkyl_array_copy(fin_cu, fin);
  gkyl_array_copy(qmem_cu, qmem);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    // zero out array struct device data
    gkyl_array_clear_cu(rhs_cu, 0.0);
    gkyl_array_clear_cu(cflrate_cu, 0.0);
    gkyl_array_clear_cu(maxs_by_cell_cu, 0.0);

    // zero out maxs values in slvr_cu
    cudaMemset(slvr_cu->maxs, 0., GKYL_MAX_DIM*sizeof(double));

    // set pointer to EM fields in vlasov equation object (on device)
    gkyl_vlasov_set_qmem(eqn_cu, qmem_cu); // must set EM fields to use

    // advance hyper_dg
    gkyl_hyper_dg_advance_cu(slvr_cu, phaseRange, fin_cu, cflrate_cu, rhs_cu, maxs_by_cell_cu);
  
    // reduction to get maxs (maximum over cells of maxs_by_cell)
    gkyl_array_reduce(maxs_by_cell_cu, GKYL_MAX, slvr_cu->maxs);
  }

  // check results of maxs
  gkyl_cu_memcpy(maxs, slvr_cu->maxs, GKYL_MAX_DIM*sizeof(double), GKYL_CU_MEMCPY_D2H);
  TEST_CHECK( gkyl_compare_double(maxs[0], 0.0000000000000000e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[1], 0.0000000000000000e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[2], 1.1565625000000002e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[3], 1.1571875000000000e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[4], 1.1578125000000001e+00, 1e-12) );

  // copy result from device to host
  gkyl_array_copy(rhs, rhs_cu);

  // get linear index of first non-ghost cell
  int idx[] = {0, 0, 0, 0, 0};
  int linl = gkyl_range_idx(&phaseRange, idx);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d;
  int i = 0;
  while(val==0) {
    rhs_d = (double*) gkyl_array_fetch(rhs, i);
    val = rhs_d[1];
    if(val==0) i++;
  }
  TEST_CHECK(i == linl);

  // check data in first non-ghost cell
  rhs_d = (double*) gkyl_array_fetch(rhs, linl);
  TEST_CHECK( gkyl_compare_double(rhs_d[0],  4.894191610931403e+00 , 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  1.331341236610166e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -8.324741199843084e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  -2.673619244471518e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  4.243853346589546e+00 , 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  -6.618097617712298e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -2.204812073157034e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  1.233289696343498e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -1.462768984400967e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  1.711468109826621e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10], -4.311578661561084e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], -1.055566810350123e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12], -1.462989088297130e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 3.475894627821628e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14], 1.049594917536618e+01 , 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 9.761290631219577e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16], -1.621943928500327e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], 2.722448306561996e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18], -1.069595454463273e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[19], -4.546812665639302e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[20], -2.315401840681718e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[21], 1.084202513678910e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[22], 1.965414256580261e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[23], 9.307169707410683e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[24], -1.159443470189338e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[25], 8.959151842026718e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[26], -1.922720351907200e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[27], 1.939005961012885e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[28], -2.187959418345774e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[29], 9.872685698935690e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[30], -4.888428376000754e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[31], 2.508827706127977e+01 , 1e-12) );

  // get linear index of some other cell
  int idx2[] = {5, 2, 4, 7, 1};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = (double*) gkyl_array_fetch(rhs, linl2);

  TEST_CHECK( gkyl_compare_double(rhs_d[0],  4.7846096908261693e-02 , 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  1.1169641822719267e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -5.7985839866517281e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  5.6000000000000183e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  -3.5215390309173478e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  3.1215390309173235e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -5.9451852405668802e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  1.0813746062398621e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -4.9765990416585353e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  3.4930985546332849e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10], -5.7062455762946534e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], 3.1215390309173224e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12], -3.5100138306029467e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 4.8317490964220532e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14], -3.2905989232414845e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 7.9094010767584866e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16], -6.8157586550966215e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], 5.1498741574033772e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18], -3.5101085195914536e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 4.8317866239642733e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[20], -6.0284065563081157e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[21], 3.5332118170962823e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[22], -5.7103189399232193e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[23], 1.0413981651809296e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[24], -4.9366836100923372e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[25], 7.2784609690826585e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[26], -6.0284408668301865e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[27], 6.0244518005727826e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[28], -6.0036637472496942e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[29], 1.0962068462170658e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[30], -5.7779480060160040e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[31], 6.0644701637051824e+01 , 1e-12) );

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  gkyl_array_release(fin_cu);
  gkyl_array_release(rhs_cu);
  gkyl_array_release(cflrate_cu);
  gkyl_array_release(qmem_cu);
  gkyl_array_release(maxs_by_cell_cu);

  gkyl_free(maxs);
  gkyl_cu_free(maxs_cu);

//  gkyl_hyper_dg_release(slvr_cu);
//  gkyl_dg_eqn_release(eqn_cu);
}
