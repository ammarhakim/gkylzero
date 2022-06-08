#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_ghost_surf_calc.h>

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

//#define mkarr1(use_gpu, nc, size) (fprintf(stderr, "mkarr1: %d\n", __LINE__), mkarr1_(use_gpu, nc, size))

int ghost_surf_calc_kernel_test(const gkyl_ghost_surf_calc *slvr);

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
  enum gkyl_field_id field_id = GKYL_FIELD_E_B;
  eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange, field_id, use_gpu);

  // initialize hyper_dg slvr
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2};
  int zero_flux_flags[GKYL_MAX_DIM] = {0, 1, 1};

  gkyl_ghost_surf_calc *slvr;
  slvr = gkyl_ghost_surf_calc_new(&phaseGrid, eqn);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  qmem = mkarr1(use_gpu, 8*confBasis.num_basis, confRange_ext.volume);

  double *cfl_ptr;
  if (use_gpu)
    cfl_ptr = gkyl_cu_malloc(sizeof(double));
  else
    cfl_ptr = gkyl_malloc(sizeof(double));

  // set initial condition
  int nf = phaseRange_ext.volume*basis.num_basis;
  double *fin_d;
  if (use_gpu) {
    fin_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
    fin_d = fin_h->data;
  } else {
    fin_d = fin->data;
  }
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(fin, fin_h);

  int nem = confRange_ext.volume*confBasis.num_basis;
  double *qmem_d;
  if (use_gpu) {
    qmem_h = mkarr1(false, 8*confBasis.num_basis, confRange_ext.volume);
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
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_vlasov_set_auxfields(eqn,
      (struct gkyl_dg_vlasov_auxfields) { .qmem = qmem }); // Must set EM fields to use.
    if (use_gpu)
      gkyl_ghost_surf_calc_advance_cu(slvr, &phaseRange_ext, &confRange_ext, fin, rhs);
    else
      gkyl_ghost_surf_calc_advance(slvr, &phaseRange_ext, &confRange_ext, fin, rhs);
  }

  double cfl_ptr_h[1];
  if (use_gpu)
    gkyl_cu_memcpy(cfl_ptr_h, cfl_ptr, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    cfl_ptr_h[0] = cfl_ptr[0];
  TEST_CHECK( gkyl_compare_double(cfl_ptr_h[0], 2.5178875733842702e+01, 1e-12) );

  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  rhs_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
  gkyl_array_copy(rhs_h, rhs);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d;
  rhs_d = gkyl_array_fetch(rhs_h, 0);

  //printf("first cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0],  1.1873679155168162e+00, 1e-12) );
  printf("rhs_d[0] = %f, %f\n", rhs_d[0], 1.1873679155168162e+00);
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  2.2131609273296690e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -5.1421989058620037e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  -2.1015559414749911e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  -3.1274299616328875e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  -5.7980678014735298e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -5.1794065426655207e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  1.6555476342162777e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -3.7839306076374341e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  6.0821819985513692e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10],  -3.0765279980185458e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11],  -2.4545628123190717e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12],  -1.1912775540624752e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13],  9.1266993445490581e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14],  -4.5599700127748999e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[15],  1.1971791673432133e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16],  3.6053906660129060e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17],  -2.4385757910567278e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18],  -1.1795732611711490e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[19],  1.8857755500762241e+01, 1e-12) );

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_h);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  gkyl_ghost_surf_calc_release(slvr);
  gkyl_dg_eqn_release(eqn);

  if (use_gpu) {
    gkyl_cu_free(cfl_ptr);
    
    gkyl_array_release(fin_h);
    gkyl_array_release(qmem_h);
  }
  else {
    gkyl_free(cfl_ptr);
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
  enum gkyl_field_id field_id = GKYL_FIELD_E_B;
  eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange, field_id, use_gpu);

  // initialize hyper_dg slvr
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2, 3, 4};
  int zero_flux_flags[GKYL_MAX_DIM] = {0, 0, 1, 1, 1};

  gkyl_ghost_surf_calc *slvr;
  slvr = gkyl_ghost_surf_calc_new(&phaseGrid, eqn);

  // initialize arrays
  struct gkyl_array *fin, *rhs, *cflrate, *qmem;
  struct gkyl_array *fin_h, *qmem_h, *rhs_h;
  
  fin = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  rhs = mkarr1(use_gpu, basis.num_basis, phaseRange_ext.volume);
  cflrate = mkarr1(use_gpu, 1, phaseRange_ext.volume);
  qmem = mkarr1(use_gpu, 8*confBasis.num_basis, confRange_ext.volume);

  // set initial condition
  int nf = phaseRange_ext.volume*basis.num_basis;
  double *fin_d;
  if (use_gpu) {
    fin_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
    fin_d = fin_h->data;
  } else {
    fin_d = fin->data;
  }
  for(int i=0; i< nf; i++) {
    fin_d[i] = (double)(2*i+11 % nf) / nf  * ((i%2 == 0) ? 1 : -1);
  }
  if (use_gpu) gkyl_array_copy(fin, fin_h);

  int nem = confRange_ext.volume*confBasis.num_basis;
  double *qmem_d;
  if (use_gpu) {
    qmem_h = mkarr1(false, 8*confBasis.num_basis, confRange_ext.volume);
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
    gkyl_array_clear(rhs, 0.0);
    gkyl_array_clear(cflrate, 0.0);
    gkyl_vlasov_set_auxfields(eqn,
      (struct gkyl_dg_vlasov_auxfields) { .qmem = qmem }); // must set EM fields to use
    if (use_gpu)
      gkyl_ghost_surf_calc_advance_cu(slvr, &phaseRange_ext, &confRange_ext, fin, rhs);
    else
      gkyl_ghost_surf_calc_advance(slvr, &phaseRange_ext, &confRange_ext, fin, rhs);
  }

  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  rhs_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
  gkyl_array_copy(rhs_h, rhs);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d;
  rhs_d = gkyl_array_fetch(rhs_h, 0);

  TEST_CHECK( gkyl_compare_double(rhs_d[0],  5.7844432740288285e+00, 1e-12) );
  printf("rhs_d[0] = %f, %f\n", rhs_d[0], 5.7844432740288285e+00);
  TEST_CHECK( gkyl_compare_double(rhs_d[1],  1.2420962193260433e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[2],  -7.4300925874823767e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[3],  4.2941817918870306e-02, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[4],  3.2534445627054236e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[5],  -6.7449532498993019e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[6],  -2.2944967853674846e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[7],  9.6125249593607940e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[8],  -1.1903507003208647e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[9],  1.8104933822676330e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[10],  -5.3016754408285633e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[11],  6.5658328154814427e-02, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[12],  -1.4502876530701883e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[13],  3.3487215685532625e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[14],  -1.5449712957466879e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[15],  -1.5127753981281999e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[16],  -1.8947433067488099e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[17],  2.8214423785413441e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[18],  -2.1313744728598749e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[19],  1.0159572731623626e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[20],  -2.3026686627466411e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[21],  2.2887034800363711e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[22],  -1.0083684729265356e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[23],  1.1799884697653724e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[24],  -1.4090960447465381e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[25],  -6.0458094032640641e+00, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[26],  -2.9837921269569524e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[27],  3.1443247788436210e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[28],  -1.9379257908892853e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[29],  1.6834219017489403e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[30],  -1.1869770681550872e+01, 1e-12) ); 
  TEST_CHECK( gkyl_compare_double(rhs_d[31],  3.2089428097916851e+01, 1e-12) );
  
  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_h);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  gkyl_ghost_surf_calc_release(slvr);
  gkyl_dg_eqn_release(eqn);

  if (use_gpu) {
    gkyl_array_release(fin_h);
    gkyl_array_release(qmem_h);
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
int ghost_surf_calc_kernel_test(const gkyl_ghost_surf_calc *slvr) {
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

