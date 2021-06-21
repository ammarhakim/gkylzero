#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_hyper_dg.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static struct gkyl_array*
mkarr1(bool use_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (use_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

int hyper_dg_kernel_test(const gkyl_hyper_dg *slvr);

void
test_vlasov_1x2v_p2_(bool use_gpu)
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
  if (use_gpu)
    eqn = gkyl_dg_vlasov_cu_dev_new(&confBasis, &basis, &confRange);
  else
    eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange);

  // initialize hyper_dg slvr
  int up_dirs[] = {0, 1, 2};
  int zero_flux_flags[] = {0, 1, 1};
  double maxs_init[] = {0., 0., 0.};

  gkyl_hyper_dg *slvr;
  if (use_gpu)
    slvr = gkyl_hyper_dg_cu_dev_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, maxs_init);
  else
    slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, maxs_init);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *maxs_by_cell, *qmem;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  fin = mkarr1(use_gpu, basis.numBasis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.numBasis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  maxs_by_cell = mkarr1(use_gpu, pdim, phaseRange_ext.volume);
  qmem = mkarr1(use_gpu, 8*confBasis.numBasis, confRange_ext.volume);

  double *cfl_ptr;
  if (use_gpu)
    cfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    cfl_ptr = gkyl_malloc(sizeof(double));

  // set initial condition
  int nf = phaseRange_ext.volume*basis.numBasis;
  double *fin_d;
  if (use_gpu) {
    fin_h = mkarr1(false, basis.numBasis, phaseRange_ext.volume);
    fin_d = fin_h->data;
  } else {
    fin_d = fin->data;
  }
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(fin, fin_h);

  int nem = confRange_ext.volume*confBasis.numBasis;
  double *qmem_d;
  if (use_gpu) {
    qmem_h = mkarr1(false, 8*confBasis.numBasis, confRange_ext.volume);
    qmem_d = qmem_h->data;
  } else {
    qmem_d = qmem->data;
  }
  for(int i=0; i< nem; i++) {
    qmem_d[i] = (double)(-i+27 % nem) / nem  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(qmem, qmem_h);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    if (use_gpu)
      gkyl_cu_memset(slvr->maxs, 0, GKYL_MAX_DIM*sizeof(double));
    else
      for (int i=0; i<pdim; i++) slvr->maxs[i] = 0.;
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_array_clear(maxs_by_cell, 0.0);
    gkyl_vlasov_set_qmem(eqn, qmem); // must set EM fields to use
    if (use_gpu) 
      gkyl_hyper_dg_advance_cu(slvr, phaseRange, fin, cflrate, rhs, maxs_by_cell);
    else
      gkyl_hyper_dg_advance(slvr, phaseRange, fin, cflrate, rhs, maxs_by_cell);

    gkyl_array_reduce(cfl_ptr, cflrate, GKYL_MAX);
    gkyl_array_reduce_range(slvr->maxs, maxs_by_cell, GKYL_MAX, phaseRange);
  }

  double *cfl_ptr_h;
  if (use_gpu) {
    cfl_ptr_h = (double*) gkyl_malloc(sizeof(double));
    gkyl_cu_memcpy(cfl_ptr_h, cfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
    cfl_ptr = cfl_ptr_h;
  }
  TEST_CHECK( gkyl_compare_double(cfl_ptr[0], 2.5178875733842702e+01, 1e-12) );

  // check results of maxs
  double *maxs;
  if (use_gpu) {
    maxs = (double*) gkyl_malloc(GKYL_MAX_DIM*sizeof(double));
    gkyl_cu_memcpy(maxs, slvr->maxs, GKYL_MAX_DIM*sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    maxs = slvr->maxs;
  }
    
  //printf("maxs\n");
  //for(int i=0; i<pdim; i++) printf("%.16e\n", maxs[i]);
  TEST_CHECK( gkyl_compare_double(maxs[0], 0.0000000000000000e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[1], 9.6634593835823734e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[2], 9.9844695137960127e-02, 1e-12) );

  // get linear index of first non-ghost cell
  int idx[] = {0, 0, 0, 0, 0};
  int linl = gkyl_range_idx(&phaseRange, idx);

  if (use_gpu) {
    rhs_h = mkarr1(false, basis.numBasis, phaseRange_ext.volume);
    gkyl_array_copy(rhs_h, rhs);
    rhs = rhs_h;
  }

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

  // check data in first non-ghost cell
  rhs_d = gkyl_array_fetch(rhs, linl);

  //printf("first cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0],  1.0546061550871522e+00 , 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  3.5130683492927017e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -5.2594564188908226e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  -2.2234588402505593e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  -3.0877761745598523e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  -5.7612934659242940e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -5.1953592457808213e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  1.6308176265956352e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -5.6284864625968856e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  4.9034564530949115e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10], -3.0404733449958361e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], -2.4690047322837483e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12], -1.1449365614068894e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 8.8675676590624661e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14], -5.3388252825551130e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 1.1838756032288330e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16], 3.6532177298741773e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], -2.4342828318438372e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18], -1.1414305139756932e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 1.9010564666002338e+01 , 1e-12) );

  // get linear index of some other cell
  int idx2[] = {5, 2, 4};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs, linl2);

  //printf("second cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0],  7.9777292154304535e-01 , 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  -2.7182122357080729e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -3.1823319852967531e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  -1.3560732323031068e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  -9.3361523823747120e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  -6.4546439524194739e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -3.0046857486230216e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  5.7739758711920025e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -5.5255780047650815e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  5.4569313596807967e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10], -9.3082985569032928e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], -5.9264062443889756e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12], -3.4252162191461373e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 5.2936631523369080e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14], -5.5255780047650616e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 3.7711009940704315e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16], 2.8771855264879633e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], -5.8883885687489297e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18], -3.4263348605987986e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 4.1758987770172112e+01 , 1e-12) );

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(maxs_by_cell);
  gkyl_array_release(qmem);

  if (!use_gpu) {
    gkyl_hyper_dg_release(slvr);
    gkyl_dg_eqn_release(eqn);
  } else {
    // need to figure out how to release on device
  }
}

void
test_vlasov_2x3v_p1_(bool use_gpu)
{
  // initialize grid and ranges
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

  // initialize basis
  int poly_order = 1;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn
  struct gkyl_dg_eqn *eqn;
  if (use_gpu)
    eqn = gkyl_dg_vlasov_cu_dev_new(&confBasis, &basis, &confRange);
  else
    eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange);

  // initialize hyper_dg slvr
  int up_dirs[] = {0, 1, 2, 3, 4};
  int zero_flux_flags[] = {0, 0, 1, 1, 1};
  double maxs_init[] = {0., 0., 0., 0., 0.};

  gkyl_hyper_dg *slvr;
  if (use_gpu) {
    slvr = gkyl_hyper_dg_cu_dev_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, maxs_init);
    // basic check
    int nfail = hyper_dg_kernel_test(slvr);

    TEST_CHECK( nfail == 0 );
  } else {
    slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, maxs_init);
  }

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *maxs_by_cell, *qmem;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  fin = mkarr1(use_gpu, basis.numBasis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.numBasis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  maxs_by_cell = mkarr1(use_gpu, pdim, phaseRange_ext.volume);
  qmem = mkarr1(use_gpu, 8*confBasis.numBasis, confRange_ext.volume);

  // set initial condition
  int nf = phaseRange_ext.volume*basis.numBasis;
  double *fin_d;
  if (use_gpu) {
    fin_h = mkarr1(false, basis.numBasis, phaseRange_ext.volume);
    fin_d = fin_h->data;
  } else {
    fin_d = fin->data;
  }
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(fin, fin_h);

  int nem = confRange_ext.volume*confBasis.numBasis;
  double *qmem_d;
  if (use_gpu) {
    qmem_h = mkarr1(false, 8*confBasis.numBasis, confRange_ext.volume);
    qmem_d = qmem_h->data;
  } else {
    qmem_d = qmem->data;
  }
  for(int i=0; i< nem; i++) {
    qmem_d[i] = (double)(-i+27 % nem) / nem  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(qmem, qmem_h);

  // run hyper_dg_advance
  int nrep = 10;
  for(int n=0; n<nrep; n++) {
    if (use_gpu)
      gkyl_cu_memset(slvr->maxs, 0, GKYL_MAX_DIM*sizeof(double));
    else
      for (int i=0; i<pdim; i++) slvr->maxs[i] = 0.;
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_array_clear(maxs_by_cell, 0.0);
    gkyl_vlasov_set_qmem(eqn, qmem); // must set EM fields to use
    if (use_gpu) 
      gkyl_hyper_dg_advance_cu(slvr, phaseRange, fin, cflrate, rhs, maxs_by_cell);
    else
      gkyl_hyper_dg_advance(slvr, phaseRange, fin, cflrate, rhs, maxs_by_cell);

    gkyl_array_reduce_range(slvr->maxs, maxs_by_cell, GKYL_MAX, phaseRange);
  }

  // check results of maxs
  double *maxs;
  if (use_gpu) {
    maxs = (double*) gkyl_malloc(GKYL_MAX_DIM*sizeof(double));
    gkyl_cu_memcpy(maxs, slvr->maxs, GKYL_MAX_DIM*sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    maxs = slvr->maxs;
  }
    
  TEST_CHECK( gkyl_compare_double(maxs[0], 0.0000000000000000e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[1], 0.0000000000000000e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[2], 1.1565625000000002e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[3], 1.1571875000000000e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[4], 1.1578125000000001e+00, 1e-12) );

  // get linear index of first non-ghost cell
  int idx[] = {0, 0, 0, 0, 0};
  int linl = gkyl_range_idx(&phaseRange, idx);

  if (use_gpu) {
    rhs_h = mkarr1(false, basis.numBasis, phaseRange_ext.volume);
    gkyl_array_copy(rhs_h, rhs);
    rhs = rhs_h;
  }

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

  // check data in first non-ghost cell
  rhs_d = gkyl_array_fetch(rhs, linl);

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
  rhs_d = gkyl_array_fetch(rhs, linl2);

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
  gkyl_array_release(maxs_by_cell);
  gkyl_array_release(qmem);

  if (!use_gpu) {
    gkyl_hyper_dg_release(slvr);
    gkyl_dg_eqn_release(eqn);
  } else {
    // need to figure out how to release on device
  }
}


void
test_vlasov_1x2v_p2()
{
  test_vlasov_1x2v_p2_(false);
}

void
test_vlasov_1x2v_p2_cu()
{
  test_vlasov_1x2v_p2_(true);
}

void
test_vlasov_2x3v_p1()
{
  test_vlasov_2x3v_p1_(false);
}

void
test_vlasov_2x3v_p1_cu()
{
  test_vlasov_2x3v_p1_(true);
}

#ifndef GKYL_HAVE_CUDA
int hyper_dg_kernel_test(const gkyl_hyper_dg *slvr) {
  return 0;
}
#endif

TEST_LIST = {
  { "test_vlasov_1x2v_p2", test_vlasov_1x2v_p2 },
  { "test_vlasov_2x3v_p1", test_vlasov_2x3v_p1 },
#ifdef GKYL_HAVE_CUDA
  { "test_vlasov_1x2v_p2_cu", test_vlasov_1x2v_p2_cu },
  { "test_vlasov_2x3v_p1_cu", test_vlasov_2x3v_p1_cu },
#endif
  { NULL, NULL },
};

