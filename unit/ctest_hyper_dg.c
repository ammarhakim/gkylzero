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

  int cells[] = {8, 8, 8};
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
  bool no_iter = false;
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
  TEST_CHECK( gkyl_compare_double(maxs[1], 4.4358364283918261e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(maxs[2], 5.2704627669472981e-02, 1e-12) );

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

  TEST_CHECK( gkyl_compare_double(rhs_d[0],  1.0129572842068597e+00 , 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  9.6300495720038637e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -3.8613141902236516e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  -1.8709837234164641e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  -2.5957207692607465e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  -3.9764150953959811e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -4.2125270511714819e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  1.4135512227445936e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -7.0080241696921286e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  5.0632236079517057e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10],  -2.5203869273562059e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11],  -2.0442213321820187e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12],  -7.7956265730269116e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13],  8.0299504751403443e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14],  -7.5378638043775570e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[15],  9.2306475566436959e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16],  3.2601578306660519e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17],  -2.0208869486349624e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18],  -7.7839511018276157e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[19],  1.4657430686821954e+01 , 1e-12) );

  // get linear index of some other cell
  int idx2[] = {5, 2, 4};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs, linl2);

  TEST_CHECK( gkyl_compare_double(rhs_d[0],  6.4031853395068872e-01 , 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  4.5052817963633860e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -1.7625365528308432e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  -5.5968146604931235e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  -5.0004235571209051e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  2.4237255186803299e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -1.5315964451549977e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  3.6445955578949061e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -1.3268449456435782e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  1.2376043070340126e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10], -4.9637963801073468e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], -3.3903170306672529e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12], -1.4027410188568949e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 3.3781390320955396e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14], -1.3268449456435838e-01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 1.8610326300838032e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16], 1.5547005383792492e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], -3.3407603143472507e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18], -1.4037630451059291e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 2.1102202333430139e+01 , 1e-12) );

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
  bool no_iter = false;
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

