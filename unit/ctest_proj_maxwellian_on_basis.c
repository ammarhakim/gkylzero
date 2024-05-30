#include "gkyl_array.h"
#include "gkyl_util.h"
#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void eval_M0(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void eval_M1i_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}

void eval_M2_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double n = 1.0, vth2 = 1.0, ux = 0.5;
  double x = xn[0];
  fout[0] = n*vth2 + n*ux*ux;
}

void eval_udrift_1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}

void eval_vtsq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}

void
test_1x1v(int poly_order, bool use_gpu)
{
  double lower[] = {0.1, -6.0}, upper[] = {1.0, 6.0};
  int cells[] = {2, 32};
  const int vdim = 1;
  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim-vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  int vGhost[vdim] = { 0 };
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M2_1v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // proj_maxwellian expects the moments as a single array.
  struct gkyl_array *moms_ho = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(moms_ho, 1., m0 , 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m1i, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m2 , (vdim+1)*confBasis.num_basis);
  struct gkyl_array *moms;
  if (use_gpu) { // copy host array to device
    moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, (vdim+2)*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(moms, moms_ho);
  } else {
    moms = moms_ho;
  }

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu;
  if (use_gpu)  // create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // projection updater to compute Maxwellian
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, gvm, use_gpu);

  if (use_gpu) {
    gkyl_proj_maxwellian_on_basis_lab_mom(proj_max, &local, &confLocal, moms, distf_cu);
    gkyl_array_copy(distf, distf_cu);
  } else {
    gkyl_proj_maxwellian_on_basis_lab_mom(proj_max, &local, &confLocal, moms, distf);
  }

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p1_vals[] = {  7.5586260555876306e-01, 0.0000000000000000e+00, 2.5444480361615229e-02, 0.0000000000000000e+00, -3.5764626824659473e-03, 0.0000000000000000e+00 };
  double p2_vals[] = {  7.5586260555876306e-01,  3.3461741853476639e-17,  2.5444480361615243e-02,
    9.2374888136351991e-18, 3.6409663532636887e-16, -3.5764626824658884e-03,
    5.2417618568471655e-17, -3.0935326627861718e-18 };
  double p3_vals[] = { 7.5586259435651881e-01, -3.7627349011120229e-17,  2.5444733034629394e-02,
    -9.8131366622101836e-18, -1.1779064054377974e-16, -3.5692341122334492e-03,
    -5.9602614442099326e-18, -5.9602614442099326e-18,  2.3674858872649741e-17,
    -1.1364096470021106e-04,  5.7229637382275112e-19,  4.0417433257763653e-18
  };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]) { 1, 17 }));
  
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );    
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 3) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p3_vals[i], fv[i], 1e-12) );
  }

//  // write distribution function to file
//  char fname[1024];
//  sprintf(fname, "ctest_proj_maxwellian_on_basis_test_1x1v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);

  // now perform the maxwellian projection using primitive moments.
  struct gkyl_array *udrift, *vtsq;
  udrift = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift_1v, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq, NULL);

  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);

  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr((vdim+1)*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1., udrift, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., vtsq  , vdim*confBasis.num_basis);
  struct gkyl_array *prim_moms;
  if (use_gpu) { // copy host array to device
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, (vdim+1)*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms, prim_moms_ho);
  } else {
    prim_moms = prim_moms_ho;
  }

  if (use_gpu) {
    gkyl_proj_maxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, moms, prim_moms, distf_cu);
    gkyl_array_copy(distf, distf_cu);
  } else {
    gkyl_proj_maxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, moms, prim_moms, distf);
  }

  fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]) { 1, 17 }));
  
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );    
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 3) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p3_vals[i], fv[i], 1e-12) );
  }

//  // write distribution function to file
//  char fname[1024];
//  sprintf(fname, "ctest_proj_maxwellian_on_basis_test_1x1v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  gkyl_array_release(udrift); gkyl_array_release(vtsq);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);

  gkyl_array_release(moms_ho); gkyl_array_release(prim_moms_ho);
  if (use_gpu) {
    gkyl_array_release(moms); gkyl_array_release(prim_moms);
  }

  gkyl_array_release(distf);
  if (use_gpu)
    gkyl_array_release(distf_cu);
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);
}

void eval_M1i_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; fout[1] = 0.25;
}

void eval_M2_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double n = 1.0, vth2 = 1.0, ux = 0.5, uy = 0.25;
  double x = xn[0];
  fout[0] = 2*n*vth2 + n*(ux*ux+uy*uy);
}

void eval_udrift_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; fout[1] = 0.25;
}

void
test_1x2v(int poly_order, bool use_gpu)
{
  double lower[] = {0.1, -6.0, -6.0}, upper[] = {1.0, 6.0, 6.0};
  int cells[] = {2, 16, 16};
  const int vdim = 2;
  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim-vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  int vGhost[vdim] = {0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_M1i_2v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M2_2v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // proj_maxwellian expects the moments as a single array.
  struct gkyl_array *moms_ho = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(moms_ho, 1., m0 , 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m1i, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m2 , (vdim+1)*confBasis.num_basis);
  struct gkyl_array *moms;
  if (use_gpu) { // copy host array to device
    moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, (vdim+2)*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(moms, moms_ho);
  } else {
    moms = moms_ho;
  }

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu;
  if (use_gpu)  // create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // projection updater to compute Maxwellian
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, gvm, use_gpu);

  if (use_gpu) {
    gkyl_proj_maxwellian_on_basis_lab_mom(proj_max, &local, &confLocal, moms, distf_cu);
    gkyl_array_copy(distf, distf_cu);
  } else {
    gkyl_proj_maxwellian_on_basis_lab_mom(proj_max, &local, &confLocal, moms, distf);
  }

  // values to compare  at index (1, 9, 9) [remember, lower-left index is (1,1,1)]
  double p1_vals[] = {  
    4.2337474137655012e-01, 0.0000000000000000e+00, 1.1241037221784697e-02, -1.1241037221784659e-02, 0.0000000000000000e+00, 0.0000000000000000e+00,
    -2.9846116329636935e-04, 0.0000000000000000e+00, -8.7555592875305493e-03, 0.0000000000000000e+00, 2.3246915375411265e-04, 0.0000000000000000e+00,
    -8.7555592875305354e-03, 0.0000000000000000e+00, -2.3246915375410571e-04, 0.0000000000000000e+00 };
  double p2_vals[] = { 4.2337474137655023e-01,  5.0502880544733958e-17,  1.1241037221784692e-02,
    -1.1241037221784697e-02, -4.7391077427032355e-18,  2.1997861612039929e-18,
    -2.9846116329638447e-04,  1.7051554382338028e-16, -8.7555592875305649e-03,
    -8.7555592875305441e-03,  2.7340310249380901e-20, -7.7853344310736426e-18,
    -8.4644052716641422e-19, -1.2104293181211046e-18,  2.3246915375412010e-04,
    -8.4644052716641422e-19, -2.3246915375411175e-04, -6.1682793769044501e-19,
    2.8526190142631692e-18,  1.1178955382863621e-18  };
  double p3_vals[] = { 4.2337367481234789e-01, -6.0016526586247659e-18,  1.1242923855289584e-02,
    -1.1242923855289581e-02, -3.9968671012250635e-19, -3.9968671012250635e-19,
    -2.9856210798149093e-04, -7.3108077460892585e-17, -8.6780466947367265e-03,
    -8.6780466947367196e-03, -4.8918503715040388e-19, -1.0074751858744896e-17,
    3.8030359490695606e-18,  1.2009507351043501e-18,  2.3045036573145683e-04,
    1.2009507351043501e-18, -2.3045036573145770e-04, -5.3277715857024823e-18,
    -2.0839365136087289e-04,  2.0839365136087100e-04,  6.2223560245147039e-18,
    1.8855473345726862e-18,  1.8855473345726862e-18, -1.3367114158171800e-17,
    5.1067364964265611e-19,  6.0263289067347938e-18,  5.5340095371337203e-06,
    3.8847760981017097e-19,  5.5340095371246130e-06, -2.6567508608833698e-18,
    2.5841550939375440e-18,  4.1575074896653521e-19 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
  
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 3) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p3_vals[i], fv[i], 1e-10) );
  }

//  // write distribution function to file
//  char fname[1024];
//  sprintf(fname, "ctest_proj_maxwellian_on_basis_test_1x2v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1i); gkyl_array_release(m2);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);

  // now perform the maxwellian projection using primitive moments.
  struct gkyl_array *udrift, *vtsq;
  udrift = mkarr(vdim*confBasis.num_basis, confLocal_ext.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift_2v, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq, NULL);

  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);

  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr((vdim+1)*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1., udrift, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., vtsq  , vdim*confBasis.num_basis);
  struct gkyl_array *prim_moms;
  if (use_gpu) { // copy host array to device
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, (vdim+1)*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms, prim_moms_ho);
  } else {
    prim_moms = prim_moms_ho;
  }

  if (use_gpu) {
    gkyl_proj_maxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, moms, prim_moms, distf_cu);
    gkyl_array_copy(distf, distf_cu);
  } else {
    gkyl_proj_maxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, moms, prim_moms, distf);
  }

  fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
  
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 3) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p3_vals[i], fv[i], 1e-10) );
  }

//  // write distribution function to file
//  char fname[1024];
//  sprintf(fname, "ctest_proj_maxwellian_on_basis_test_1x2v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  gkyl_array_release(udrift); gkyl_array_release(vtsq);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);

  gkyl_array_release(moms_ho); gkyl_array_release(prim_moms_ho);
  if (use_gpu) {
    gkyl_array_release(moms); gkyl_array_release(prim_moms);
  }

  gkyl_array_release(distf);
  if (use_gpu)
    gkyl_array_release(distf_cu);
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);

}

void eval_udrift_2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}
void eval_M1_2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1], udrift[1];
  eval_M0(t, xn, den, ctx);
  eval_udrift_2v_gk(t, xn, udrift, ctx);
  fout[0] = den[0]*udrift[0];
}
void eval_M2_2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1], udrift[1], vtsq[1];
  eval_M0(t, xn, den, ctx);
  eval_udrift_2v_gk(t, xn, udrift, ctx);
  eval_vtsq(t, xn, vtsq, ctx);
  fout[0] = 3.*den[0]*vtsq[0] + den[0]*(udrift[0]*udrift[0]);
}
void eval_bmag(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}
void eval_jacob_tot(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}


void
test_1x2v_gk(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double lower[] = {0.1, -6.0, 0.0}, upper[] = {1.0, 6.0, 6.0};
  int cells[] = {2, 16, 16};
  const int vdim = 2;
  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim-vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  int vGhost[vdim] = {0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  // create moment arrays
  struct gkyl_array *m0, *m1, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M1_2v_gk, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M2_2v_gk, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1, 0.0, &confLocal, m1);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // proj_maxwellian expects the moments as a single array.
  struct gkyl_array *moms_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(moms_ho, 1., m0, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m1, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_ho, 1., m2, 2*confBasis.num_basis);
  struct gkyl_array *moms;
  if (use_gpu) { // copy host array to device
    moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(moms, moms_ho);
  } else {
    moms = moms_ho;
  }

  // create bmag and jacob_tot arrays
  struct gkyl_array *bmag, *jacob_tot;
  bmag = mkarr(confBasis.num_basis, confLocal_ext.volume);
  jacob_tot = mkarr(confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *bmag_cu, *jacob_tot_cu;
  if (use_gpu) { // create device copies
    bmag_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
    jacob_tot_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  }
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bmag, NULL);
  gkyl_proj_on_basis *proj_jac = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_jacob_tot, NULL);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag);
  gkyl_proj_on_basis_advance(proj_jac, 0.0, &confLocal, jacob_tot);
  if (use_gpu) {  // copy host array to device
    gkyl_array_copy(bmag_cu, bmag);
    gkyl_array_copy(jacob_tot_cu, jacob_tot);
  }

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu;
  if (use_gpu)  // create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // projection updater to compute Maxwellian
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, gvm, use_gpu);

  if (use_gpu) {
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max, &local, &confLocal, moms,
                                            bmag_cu, jacob_tot_cu, mass, distf_cu);
    gkyl_array_copy(distf, distf_cu);
  } else {
    gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_max, &local, &confLocal, moms,
                                            bmag, jacob_tot, mass, distf);
  }

  // values to compare  at index (1, 9, 9) [remember, lower-left index is (1,1,1)]
  double p1_vals[] = {  
     7.2307139183122714e-03, 0.0000000000000000e+00, 1.9198293226362615e-04, -7.7970439910196674e-04, 0.0000000000000000e+00, 0.0000000000000000e+00,
    -2.0701958137127286e-05, 0.0000000000000000e+00, -1.4953406100022537e-04, 0.0000000000000000e+00, 1.6124599381836546e-05, 0.0000000000000000e+00,
    -8.2719200283232917e-19, 0.0000000000000000e+00, -3.4806248503322844e-20, 0.0000000000000000e+00, };
  double p2_vals[] = { 
    7.2307468609012666e-03, 0.0000000000000000e+00, 1.9198380692343289e-04, -7.8092230706225602e-04, 0.0000000000000000e+00, 0.0000000000000000e+00,
    -2.0734294852987710e-05, 3.6591823321385775e-18, -1.4953474226616330e-04, 3.7739922227981074e-05, 0.0000000000000000e+00, 7.0473141211557788e-19,
    0.0000000000000000e+00, -4.8789097761847700e-19, 1.6149786206441256e-05, 0.0000000000000000e+00, 1.0020339643610290e-06, 5.4210108624275222e-20,
    0.0000000000000000e+00, 0.0000000000000000e+00 };

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
  
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
  }

//  // write distribution function to file
//  char fname[1024];
//  sprintf(fname, "ctest_proj_gkmaxwellian_on_basis_test_1x2v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0); gkyl_array_release(m1); gkyl_array_release(m2);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);

  // now perform the maxwellian projection using primitive moments.
  struct gkyl_array *udrift, *vtsq;
  udrift = mkarr(confBasis.num_basis, confLocal_ext.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift_2v_gk, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq, NULL);

  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);

  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set(prim_moms_ho, 1., moms);
  gkyl_array_set_offset(prim_moms_ho, 1., udrift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., vtsq  , 2*confBasis.num_basis);
  struct gkyl_array *prim_moms;
  if (use_gpu) { // copy host array to device
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms, prim_moms_ho);
  } else {
    prim_moms = prim_moms_ho;
  }

  if (use_gpu) {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
                                             bmag_cu, jacob_tot_cu, mass, distf_cu);
    gkyl_array_copy(distf, distf_cu);
  } else {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
                                             bmag, jacob_tot, mass, distf);
  }

  fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
  
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
  }

//  // write distribution function to file
//  char fname[1024];
//  sprintf(fname, "ctest_proj_maxwellian_on_basis_test_1x2v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, distf, fname);

  gkyl_array_release(udrift); gkyl_array_release(vtsq);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);

  gkyl_array_release(moms_ho); gkyl_array_release(prim_moms_ho);
  gkyl_array_release(bmag); gkyl_array_release(jacob_tot);
  if (use_gpu) {
    gkyl_array_release(moms); gkyl_array_release(prim_moms);
    gkyl_array_release(bmag_cu); gkyl_array_release(jacob_tot_cu);
  }

  gkyl_array_release(distf);
  if (use_gpu)
    gkyl_array_release(distf_cu);
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);

}

void test_1x1v_p0() { test_1x1v(0, false); }
void test_1x1v_p1() { test_1x1v(1, false); }
void test_1x1v_p2() { test_1x1v(2, false); }
void test_1x1v_p3() { test_1x1v(3, false); }

void test_1x2v_p0() { test_1x2v(0, false); }
void test_1x2v_p1() { test_1x2v(1, false); }
void test_1x2v_p2() { test_1x2v(2, false); }
void test_1x2v_p3() { test_1x2v(3, false); }

void test_1x2v_p1_gk() { test_1x2v_gk(1, false); }
void test_1x2v_p2_gk() { test_1x2v_gk(2, false); }

#ifdef GKYL_HAVE_CUDA
void test_1x1v_p0_gpu() { test_1x1v(0, true); }
void test_1x1v_p1_gpu() { test_1x1v(1, true); }
void test_1x1v_p2_gpu() { test_1x1v(2, true); }
void test_1x1v_p3_gpu() { test_1x1v(3, true); }

void test_1x2v_p0_gpu() { test_1x2v(0, true); }
void test_1x2v_p1_gpu() { test_1x2v(1, true); }
void test_1x2v_p2_gpu() { test_1x2v(2, true); }
void test_1x2v_p3_gpu() { test_1x2v(3, true); }

void test_1x2v_p1_gk_gpu() { test_1x2v_gk(1, true); }
void test_1x2v_p2_gk_gpu() { test_1x2v_gk(2, true); }
#endif

TEST_LIST = {
  { "test_1x1v_p0", test_1x1v_p0 },
  { "test_1x1v_p1", test_1x1v_p1 },
  { "test_1x1v_p2", test_1x1v_p2 },
  { "test_1x1v_p3", test_1x1v_p3 },
  
  { "test_1x2v_p0", test_1x2v_p0 },  
  { "test_1x2v_p1", test_1x2v_p1 },
  { "test_1x2v_p2", test_1x2v_p2 },
  { "test_1x2v_p3", test_1x2v_p3 },

  { "test_1x2v_p1_gk", test_1x2v_p1_gk },
  { "test_1x2v_p2_gk", test_1x2v_p2_gk },

#ifdef GKYL_HAVE_CUDA
  { "test_1x1v_p0_gpu", test_1x1v_p0_gpu },
  { "test_1x1v_p1_gpu", test_1x1v_p1_gpu },
  { "test_1x1v_p2_gpu", test_1x1v_p2_gpu },
  { "test_1x1v_p3_gpu", test_1x1v_p3_gpu },
  
  { "test_1x2v_p0_gpu", test_1x2v_p0_gpu },  
  { "test_1x2v_p1_gpu", test_1x2v_p1_gpu },
  { "test_1x2v_p2_gpu", test_1x2v_p2_gpu },
  { "test_1x2v_p3_gpu", test_1x2v_p3_gpu },

  { "test_1x2v_p1_gk_gpu", test_1x2v_p1_gk_gpu },
  { "test_1x2v_p2_gk_gpu", test_1x2v_p2_gk_gpu },
#endif
  { NULL, NULL },
};
