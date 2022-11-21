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

int hyper_dg_kernel_test(const gkyl_hyper_dg *slvr);

void
test_vlasov_2x3v_p1_(bool use_gpu)
{
  // initialize grid and ranges
  int cdim = 2, vdim = 3;
  int pdim = cdim+vdim;

  int cells[] = {8, 8, 8, 8, 8, 8};
  int ghost[] = {1, 1, 1, 0, 0, 0};
  double lower[] = {0., 0., 0., -1., -1., -1.};
  double upper[] = {1., 1., 1., 1., 1., 1.};

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

  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in velocity space). */
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // initialize eqn
  struct gkyl_dg_eqn *eqn;
  enum gkyl_field_id field_id = GKYL_FIELD_NULL;
  eqn = gkyl_dg_vlasov_new(&confBasis, &basis, &confRange, field_id, use_gpu);

  // initialize hyper_dg slvr // CHECK THIS (TNB)
  int up_dirs[GKYL_MAX_DIM] = {0, 1, 2, 3, 4, 5};
  int zero_flux_flags[GKYL_MAX_DIM] = {0, 0, 0, 1, 1, 1};

  gkyl_hyper_dg *slvr;
  slvr = gkyl_hyper_dg_new(&phaseGrid, &basis, eqn, pdim, up_dirs, zero_flux_flags, 1, use_gpu);

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
    qmem_d[i] = 0.; //(double)(-i+27 % nem) / nem  * ((i%2 == 0) ? 1 : -1);
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
      gkyl_hyper_dg_advance_cu(slvr, &phaseRange, fin, cflrate, rhs);
    else
      gkyl_hyper_dg_advance(slvr, &phaseRange, fin, cflrate, rhs);
  }

  // get linear index of first non-ghost cell
  // 1-indexed for interfacing with G2 Lua layer
  int idx[] = {1, 1, 1, 1, 1, 1};
  int linl = gkyl_range_idx(&phaseRange, idx);

  rhs_h = mkarr1(false, basis.num_basis, phaseRange_ext.volume);
  gkyl_array_copy(rhs_h, rhs);

  // check that ghost cells are empty
  double val = 0;
  double *rhs_d;
  int i = 0;
  while (val==0) {
    rhs_d = gkyl_array_fetch(rhs_h, i);
    val = rhs_d[0];
    if(val==0) i++;
  }
  TEST_CHECK(i == linl);

  // check data in first non-ghost cell
  rhs_d = gkyl_array_fetch(rhs_h, linl);

  printf("\nfirst cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  //for(int i=0; i<rhs->ncomp; i++) printf("  TEST_CHECK( gkyl_compare_double(rhs_d[%d], %.16e, 1e-12) );\n", i, rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0], 1.3690598923241531e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[1], 2.5959124607127443e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2], -9.1488834927265028e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3], -3.7667112687365285e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4], -5.5533422537472843e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[5], 1.0309401076758480e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6], -1.1074510946701672e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7], 4.0855765035618780e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8], -1.1163550491206880e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[9], 3.5436950494916047e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10], -7.9861410154013406e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], 5.7730481386913730e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12], -2.9883756828964786e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 1.0294587787959923e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14], -1.8108249704942395e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 1.1783500590114784e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16], -9.2287937163858977e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], 9.8524611325354616e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18], -6.0553984851628695e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 9.7347366830467656e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[20], -8.9061171167058877e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[21], 2.2236181248744577e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[22], -8.6851072802314935e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[23], -6.2311671918347092e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[24], -1.1738397710451357e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[25], 1.0226951861308613e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[26], -9.5434496888054916e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[27], 8.6677665786627554e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[28], -8.3838134849355939e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[29], 6.6911322691834907e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[30], -8.5490688922265754e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[31], 8.7237671583388074e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[32], 1.0334747013046115e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[33], 7.0975673823576679e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[34], -1.0353770140648715e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[35], -6.7914231596387657e-05, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[36], 1.5553796806952253e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[37], 9.1181559558118011e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[38], -4.4368046638159059e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[39], 1.0521370129389039e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[40], -2.5328004264881254e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[41], 1.0735589126771989e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[42], 8.0253459362875812e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[43], 9.1179421184190090e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[44], -1.0546866726854645e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[45], 7.4978411411601920e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[46], -1.0215084283778118e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[47], 1.1103845770127094e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[48], -6.3986417153681030e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[49], 3.2284501760088669e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[50], -8.5705274011403283e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[51], 6.3986417153680406e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[52], 1.7976350479367634e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[53], 9.4521357289442633e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[54], 1.6211996623166575e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[55], 9.1170087011384435e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[56], 4.4217588044088992e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[57], 1.0693322352950108e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[58], -1.5880055981183322e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[59], 1.1192740889871851e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[60], -8.9937579923000968e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[61], 5.2285772693619348e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[62], -8.7349229823544192e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[63], 9.0295862131402451e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[64], -1.0309401076758483e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[65], 3.3886652512965904e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[66], -1.0095523624279904e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[67], 7.9999999999999544e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[68], 1.5999999999999998e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[69], 1.0937174059046534e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[70], 1.4609845870289315e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[71], 1.0602022582867146e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[72], 7.5401693355437813e-03, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[73], 1.0740603221996334e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[74], -1.9547005383792518e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[75], 9.7082717581449387e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[76], -1.0514555600634669e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[77], 4.9887923446496557e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[78], -1.0255717792950517e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[79], 8.9070574748053115e+00, 1e-12) );

/*   // get linear index of some other cell */
/*   // 1-indexed for interfacing with G2 Lua layer */
  int idx2[] = {6, 3, 5, 8, 2, 1};
  int linl2 = gkyl_range_idx(&phaseRange, idx2);
  rhs_d = gkyl_array_fetch(rhs_h, linl2);

  printf("second cell rhs\n");
  //for(int i=0; i<rhs->ncomp; i++) printf("%.16e\n", rhs_d[i]);
  //for(int i=0; i<rhs->ncomp; i++) printf("  TEST_CHECK( gkyl_compare_double(rhs_d[%d], %.16e, 1e-12) );\n", i, rhs_d[i]);
  TEST_CHECK( gkyl_compare_double(rhs_d[0], 4.0254663528750839e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[1], 3.5055661041459551e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[2], -5.9719988627674560e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[3], 2.4679735333677346e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[4], 1.7433013587762800e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[5], -2.4254663528751841e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[6], -8.0895487891781030e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[7], 3.3979716522998558e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[8], -5.1465327291404527e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[9], 3.1469742464311096e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[10], -6.0378612387375462e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[11], -2.5213307702818638e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[12], -3.1458461565011586e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[13], 4.7745069764400398e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[14], 8.4558149005685354e-02, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[15], 2.5532012342852122e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[16], -8.7692916515462798e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[17], 8.3320979782297044e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[18], -3.3355757126962551e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[19], 5.1777943243601726e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[20], -8.8211329316249746e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[21], 2.9755756064845059e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[22], -5.6564985249657092e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[23], 3.1846905993719712e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[24], -4.7375950449271848e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[25], 2.4921330770281909e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[26], -8.2558126402449730e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[27], 8.6266738402830455e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[28], -8.9430491127651436e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[29], 3.4356255437200701e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[30], -6.0941623438163560e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[31], 9.4168094543678635e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[32], 4.4898811618556794e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[33], 3.5314166191155195e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[34], -4.7884211522148370e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[35], -3.4570856028669938e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[36], -2.3501026779855829e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[37], 8.8615528467765316e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[38], -3.4968136794301024e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[39], 5.0823106334645843e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[40], -2.8525893551345447e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[41], 4.8266011563965989e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[42], 4.2589410541797218e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[43], 8.5764073350286012e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[44], -7.9694811442965104e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[45], 3.5274306067035454e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[46], -4.7665520252619011e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[47], 8.3504454968954363e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[48], 2.0587714126142315e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[49], 3.1478592976962837e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[50], -5.5828271610272871e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[51], -2.0587714126142381e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[52], -2.3216542289014268e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[53], 9.1772094823822343e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[54], -3.0993674868130277e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[55], 5.9146055987912519e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[56], -2.8340314704207628e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[57], 4.8694333787823062e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[58], 1.1106299860327151e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[59], 8.1556274558128763e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[60], -8.7748814590552357e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[61], 3.1678718186298134e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[62], -5.6064269779058783e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[63], 8.7713124022956279e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[64], 2.4254663528751652e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[65], 3.1418987396851609e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[66], -4.7545872990580406e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[67], -2.6564064605509718e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[68], -2.3456406460551036e+00, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[69], 8.3529863247802268e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[70], -3.0934069288019053e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[71], 5.0903639709630930e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[72], -2.8301365035276184e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[73], 4.8270971532602601e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[74], 1.5017059221717682e-01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[75], 8.9799405237043672e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[76], -7.9470780435029781e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[77], 3.1659112606186909e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[78], -4.7786072571909806e+01, 1e-12) );
  TEST_CHECK( gkyl_compare_double(rhs_d[79], 8.8212949996830361e+01, 1e-12) );

  // clean up
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_h);
  gkyl_array_release(cflrate);
  gkyl_array_release(qmem);

  gkyl_hyper_dg_release(slvr);
  gkyl_dg_eqn_release(eqn);

  if (use_gpu) {
    gkyl_array_release(fin_h);
    gkyl_array_release(qmem_h);
  }
}

void
test_vlasov_3x3v_p1()
{
  test_vlasov_2x3v_p1_(false);
}

TEST_LIST = {
  { "test_vlasov_3x3v_p1", test_vlasov_3x3v_p1 },
  { NULL, NULL },
};

