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

void test_vlasov_2x3v_p1_cu()
{
  cudaError_t err;

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

  // clone grid and ranges to device
  struct gkyl_rect_grid *confGrid_cu = gkyl_rect_grid_clone_on_cu_dev(&confGrid);    
  struct gkyl_rect_grid *phaseGrid_cu = gkyl_rect_grid_clone_on_cu_dev(&phaseGrid);    
  struct gkyl_range *confRange_cu = gkyl_range_clone_on_cu_dev(&confRange);
  struct gkyl_range *confRange_ext_cu = gkyl_range_clone_on_cu_dev(&confRange_ext);
  struct gkyl_range *phaseRange_cu = gkyl_range_clone_on_cu_dev(&phaseRange);
  struct gkyl_range *phaseRange_ext_cu = gkyl_range_clone_on_cu_dev(&phaseRange_ext);

  err = cudaGetLastError();
  if (err != cudaSuccess) {
    char str[100];
    sprintf(str, "\nafter range clone: CUDA error: %s\n", cudaGetErrorString(err));
    gkyl_exit(str);
  }

  // initialize basis (note: basis has no device implementation)
  int poly_order = 1;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn on device
  struct gkyl_dg_eqn *eqn_cu;
  eqn_cu = gkyl_dg_vlasov_cu_dev_new(&confBasis, &basis, confRange_cu);

  err = cudaGetLastError();
  if (err != cudaSuccess) {
    char str[100];
    sprintf(str, "\nafter vlasov init: CUDA error: %s\n", cudaGetErrorString(err));
    gkyl_exit(str);
  }

  // initialize hyper_dg slvr
  int up_dirs[] = {0, 1, 2, 3, 4};
  int zero_flux_flags[] = {0, 0, 1, 1, 1};

  gkyl_hyper_dg *slvr_cu;
  slvr_cu = gkyl_hyper_dg_cu_dev_new(phaseGrid_cu, &basis, eqn_cu, pdim, up_dirs, zero_flux_flags, 1);
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    char str[100];
    sprintf(str, "\nafter hyper_dg init: CUDA error: %s\n", cudaGetErrorString(err));
    gkyl_exit(str);
  }

  // basic checks
  int nfail = hyper_dg_test(slvr_cu);
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    char str[100];
    sprintf(str, "\nafter hyper_dg_test: CUDA error: %s\n", cudaGetErrorString(err));
    gkyl_exit(str);
  }

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

  // initialize device arrays as clones of host arrays
  struct gkyl_array *fin_cu = gkyl_array_clone_on_cu_dev(fin);
  struct gkyl_array *rhs_cu = gkyl_array_clone_on_cu_dev(rhs);
  struct gkyl_array *cflrate_cu = gkyl_array_clone_on_cu_dev(cflrate);
  struct gkyl_array *qmem_cu = gkyl_array_clone_on_cu_dev(qmem);

  // get host pointers to device struct data (need better way of handling these...)
  void *rhs_data_cu;
  gkyl_cu_memcpy(&rhs_data_cu, &(rhs_cu->data), sizeof(void*), GKYL_CU_MEMCPY_D2H);
  void *cflrate_data_cu;
  gkyl_cu_memcpy(&cflrate_data_cu, &(cflrate_cu->data), sizeof(void*), GKYL_CU_MEMCPY_D2H);

  double *maxs_cu;
  cudaMalloc((void**) &maxs_cu, sizeof(double)*5);

  err = cudaGetLastError();
  if (err != cudaSuccess) {
    char str[100];
    sprintf(str, "\nBefore hyper_dg advance: CUDA error: %s\n", cudaGetErrorString(err));
    gkyl_exit(str);
  }

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    cudaMemset(maxs_cu, 0., sizeof(double)*5);
    cudaMemset(rhs_data_cu, 0.0, rhs->size*rhs->esznc);
    cudaMemset(cflrate_data_cu, 0.0, sizeof(double)*phaseRange_ext.volume);
    gkyl_vlasov_set_qmem_cu(eqn_cu, qmem_cu); // must set EM fields to use

    int dB = 256;
    int dG = phaseRange.volume/dB + 1;
    gkyl_hyper_dg_advance_cu<<<dG,dB>>>(slvr_cu, phaseRange_cu, fin_cu, cflrate_cu, rhs_cu, maxs_cu);
    err = cudaGetLastError();
    if (err != cudaSuccess) {
      char str[100];
      sprintf(str, "\nCUDA error: %s\n", cudaGetErrorString(err));
      gkyl_exit(str);
    }
  }

  gkyl_array_clear(rhs, 0.0);
  gkyl_cu_memcpy(rhs->data, rhs_data_cu, rhs->size*rhs->esznc, GKYL_CU_MEMCPY_D2H);

  // get linear index of first non-ghost cell
  int idx[] = {0, 0, 0, 0, 0};
  int linl = gkyl_range_idx(&phaseRange, idx);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d;
  int i = 0;
  while(val==0) {
    rhs_d = (double*) gkyl_array_fetch(rhs, i);
    val = rhs_d[0];
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

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

//  gkyl_hyper_dg_release(slvr_cu);
//  gkyl_dg_eqn_release(eqn_cu);
}
