#include <acutest.h>
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

void
test_vlasov_2x3v_p1()
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
  eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange);

  // initialize hyper_dg slvr
  int up_dirs[] = {0, 1, 2, 3, 4};
  int zero_flux_flags[] = {0, 0, 1, 1, 1};

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  double maxs[] = {0., 0., 0., 0., 0.};
  fin = mkarr(basis.numBasis, phaseRange_ext.volume);
  rhs = mkarr(basis.numBasis, phaseRange_ext.volume);
  cflrate = mkarr(1, phaseRange_ext.volume);
  qmem = mkarr(8*confBasis.numBasis, confRange_ext.volume);

  // set initial condition
  int nf = phaseRange_ext.volume*basis.numBasis;
  double *fin_d = fin->data;
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(rand() % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }

  int nem = confRange_ext.volume*confBasis.numBasis;
  double *qmem_d = qmem->data;
  for(int i=0; i< nem; i++) {
    qmem_d[i] = (double)(rand() % nem) / nem  * ((i%2 == 0) ? 1 : -1);
  }

  // run hyper_dg_advance
  int nrep = 10;
  bool no_iter = true;
  for(int n=0; n<nrep; n++) {
    maxs[0] = 0., maxs[1] = 0., maxs[2] = 0., maxs[3] = 0., maxs[4] = 0.;
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_vlasov_set_qmem(eqn, qmem); // must set EM fields to use
    if(no_iter) 
      gkyl_hyper_dg_advance_no_iter(slvr, &phaseRange, fin, cflrate, rhs, maxs);
    else
      gkyl_hyper_dg_advance(slvr, &phaseRange, fin, cflrate, rhs, maxs);
  }

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

  // check data in first non-ghost cell
  rhs_d = gkyl_array_fetch(rhs, linl);
  TEST_CHECK( gkyl_compare(rhs_d[0],  1.0304792379738015e+00 , 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[1],  2.2627411769164027e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[2],  -1.6021470304561980e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[3],  9.1255327714994259e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[4],  7.5287043354735150e-01 , 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[5],  1.2757696533129009e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[6],  -2.3849313633072203e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[7],  1.5255581104717683e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[8],  -1.7549070431710209e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[9],  2.3405492704104105e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[10], -5.9634768140956034e+00, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[11], 1.9029908016788081e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[12], -1.6488906565561070e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[13], 4.3643099820838202e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[14], -1.9919208272924148e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[15], 5.5953966545084100e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[16], -5.6423559056750186e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[17], 3.6758785139811039e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[18], -4.1536565352187928e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[19], 1.4459466599159571e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[20], -3.8816095010640602e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[21], 2.2188443035540793e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[22], -4.0054647851542654e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[23], 1.3107709274756042e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[24], -3.6563382180307023e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[25], 1.2227888246147115e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[26], -7.9066068660090735e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[27], 4.5900671713860326e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[28], -4.3585940109181344e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[29], 5.8017670830532770e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[30], -1.0205587194140962e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[31], 2.7201257739631469e+01 , 1e-12) );

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);
}

#ifdef GKYL_HAVE_CUDA

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

  // clone grid and ranges to device
  struct gkyl_rect_grid *confGrid_cu = gkyl_rect_grid_clone_on_cu_dev(&confGrid);    
  struct gkyl_rect_grid *phaseGrid_cu = gkyl_rect_grid_clone_on_cu_dev(&phaseGrid);    
  struct gkyl_range *confRange_cu = gkyl_range_clone_on_cu_dev(&confRange);
  struct gkyl_range *confRange_ext_cu = gkyl_range_clone_on_cu_dev(&confRange_ext);
  struct gkyl_range *phaseRange_cu = gkyl_range_clone_on_cu_dev(&phaseRange);
  struct gkyl_range *phaseRange_ext_cu = gkyl_range_clone_on_cu_dev(&phaseRange_ext);

  // initialize basis (note: basis has no device implementation)
  int poly_order = 1;
  struct gkyl_basis basis, confBasis; // phase-space, conf-space basis

  gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn on device
  struct gkyl_dg_eqn *eqn_cu;
  eqn_cu = gkyl_dg_vlasov_cu_dev_new(&confBasis, &basis, &confRange_cu);

  // initialize hyper_dg slvr
  int up_dirs[] = {0, 1, 2, 3, 4};
  int zero_flux_flags[] = {0, 0, 1, 1, 1};

  gkyl_hyper_dg *slvr_cu;
  slvr_cu = gkyl_hyper_dg_cu_dev_new(&phaseGrid_cu, &basis, eqn_cu, pdim, up_dirs, zero_flux_flags, 1);

  // initialize host arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  fin = mkarr(basis.numBasis, phaseRange_ext.volume);
  rhs = mkarr(basis.numBasis, phaseRange_ext.volume);
  cflrate = mkarr(1, phaseRange_ext.volume);
  qmem = mkarr(8*confBasis.numBasis, confRange_ext.volume);

  // set initial condition
  int nf = phaseRange_ext.volume*basis.numBasis;
  double *fin_d = fin->data;
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(rand() % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }

  int nem = confRange_ext.volume*confBasis.numBasis;
  double *qmem_d = qmem->data;
  for(int i=0; i< nem; i++) {
    qmem_d[i] = (double)(rand() % nem) / nem  * ((i%2 == 0) ? 1 : -1);
  }

  // initialize device arrays as clones of host arrays
  struct gkyl_array *fin_cu, *rhs_cu, *cflrate_cu, *qmem_cu;
  struct gkyl_array *fin_cu = gkyl_array_clone_on_cu_dev(fin);
  struct gkyl_array *rhs_cu = gkyl_array_clone_on_cu_dev(rhs);
  struct gkyl_array *cflrate_cu = gkyl_array_clone_on_cu_dev(cflrate);
  struct gkyl_array *qmem_cu = gkyl_array_clone_on_cu_dev(qmem);

  double *maxs_cu;
  cudaMalloc((void**) &maxs_cu, sizeof(int)*5);

  // run hyper_dg_advance
  int nrep = 1;
  for(int n=0; n<nrep; n++) {
    cudaMemset(maxs_cu, 0, sizeof(int)*5);
    cudaMemset(rhs_cu->data, 0.0);
    cudaMemset(cflrate_cu->data, 0.0);
    gkyl_vlasov_set_qmem(eqn_cu, qmem_cu); // must set EM fields to use
    gkyl_hyper_dg_advance_cu<<<1,1>>>(slvr_cu, &phaseRange_cu, fin_cu, cflrate_cu, rhs_cu, maxs_cu);
  }

  gkyl_array_clear(rhs, 0.0);
  gkyl_array_copy(rhs, rhs_cu);

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

  // check data in first non-ghost cell
  rhs_d = gkyl_array_fetch(rhs, linl);
  TEST_CHECK( gkyl_compare(rhs_d[0],  1.0304792379738015e+00 , 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[1],  2.2627411769164027e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[2],  -1.6021470304561980e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[3],  9.1255327714994259e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[4],  7.5287043354735150e-01 , 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[5],  1.2757696533129009e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[6],  -2.3849313633072203e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[7],  1.5255581104717683e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[8],  -1.7549070431710209e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[9],  2.3405492704104105e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[10], -5.9634768140956034e+00, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[11], 1.9029908016788081e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[12], -1.6488906565561070e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[13], 4.3643099820838202e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[14], -1.9919208272924148e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[15], 5.5953966545084100e-01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[16], -5.6423559056750186e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[17], 3.6758785139811039e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[18], -4.1536565352187928e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[19], 1.4459466599159571e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[20], -3.8816095010640602e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[21], 2.2188443035540793e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[22], -4.0054647851542654e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[23], 1.3107709274756042e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[24], -3.6563382180307023e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[25], 1.2227888246147115e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[26], -7.9066068660090735e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[27], 4.5900671713860326e+01 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[28], -4.3585940109181344e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare(rhs_d[29], 5.8017670830532770e+00 , 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[30], -1.0205587194140962e+01, 1e-12) );
  TEST_CHECK( gkyl_compare(rhs_d[31], 2.7201257739631469e+01 , 1e-12) );

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  //gkyl_hyper_dg_release(slvr_cu);
  //gkyl_dg_eqn_release(eqn_cu);


}

#endif

TEST_LIST = {
  { "test_vlasov_2x3v_p1", test_vlasov_2x3v_p1 },
#ifdef GKYL_HAVE_CUDA
  { "test_vlasov_2x3v_p1_cu", test_vlasov_2x3v_p1_cu },
#endif
  { NULL, NULL },
};

