#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_eqn_type.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>
#include <gkyl_util.h>
#include <math.h>

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void 
eval_M0(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void 
eval_M1i_1v_no_drift(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0;
}

void 
eval_M2_1v_no_drift(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
eval_M1i_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; 
}

void 
eval_M2_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
eval_M1i_2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
  fout[1] = 0.25;
}

void 
eval_M2_2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
eval_M1i_3v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
  fout[1] = 0.25;
  fout[2] = -0.5;
}

void 
eval_M2_3v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void 
test_1x1v_no_drift(int poly_order)
{
  double lower[] = {0.1, -15.0}, upper[] = {1.0, 15.0};
  int cells[] = {2, 32};
  int vdim = 1, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  moms = mkarr((vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v_no_drift, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_1v_no_drift, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);
  gkyl_array_set_offset_range(moms, 1.0, m0, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m1i, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m2, (vdim+1)*confBasis.num_basis, &confLocal);

  // build the p_over_gamma
  struct gkyl_array *p_over_gamma;
  p_over_gamma = mkarr(vdim * velBasis.num_basis, velLocal.volume);

  // build gamma
  struct gkyl_array *gamma;
  gamma = mkarr(velBasis.num_basis, velLocal.volume);

  // build gamma_inv
  struct gkyl_array *gamma_inv;
  gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);

  // Make GammaV2, GammaV, GammaV_inv
  gkyl_calc_sr_vars_init_p_vars(&vel_grid, &velBasis, &velLocal,
    p_over_gamma, gamma, gamma_inv);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .p_over_gamma = p_over_gamma,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .h_ij_inv = 0,
    .det_h = 0,
    .model_id = GKYL_MODEL_SR,
    .mass = 1.0,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p1_vals[] = {5.3918752026566863e-01, -1.0910243387206232e-17, -6.0196985297046972e-02,
    5.0006050167249552e-18};
  double p2_vals[] = {5.3922143701031633e-01, -9.6625223288531320e-18, -5.7898881215132203e-02,
    7.8842251929957589e-18, 1.9166441863144966e-17, -1.0173903909543560e-02,
    1.7916734900988946e-17, 1.4245174569363429e-18};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]){1, 17}));

  if (poly_order == 1)
    for (int i = 0; i < basis.num_basis; ++i)
      TEST_CHECK(gkyl_compare_double(p1_vals[i], fv[i], 1e-12));

  if (poly_order == 2)
    for (int i = 0; i < basis.num_basis; ++i)
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_mj_on_basis_test_1x1v_p%d_no_drift.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

void test_1x1v_no_drift_p1() { test_1x1v_no_drift(1); }
void test_1x1v_no_drift_p2() { test_1x1v_no_drift(2); }

void 
test_1x1v(int poly_order)
{
  double lower[] = {0.1, -15.0}, upper[] = {1.0, 15.0}; // +/- 15 on velocity
  int cells[] = {2, 32};                                // default {2, 32}
  int vdim = 1, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1]}, velUpper[] = {upper[1]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  moms = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_1v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);
  gkyl_array_set_offset_range(moms, 1.0, m0, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m1i, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m2, (vdim+1)*confBasis.num_basis, &confLocal);

  // build the p_over_gamma
  struct gkyl_array *p_over_gamma;
  p_over_gamma = mkarr(vdim * velBasis.num_basis, velLocal.volume);

  // build gamma
  struct gkyl_array *gamma;
  gamma = mkarr(velBasis.num_basis, velLocal.volume);

  // build gamma_inv
  struct gkyl_array *gamma_inv;
  gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);

  // Make GammaV2, GammaV, GammaV_inv
  gkyl_calc_sr_vars_init_p_vars(&vel_grid, &velBasis, &velLocal,
    p_over_gamma, gamma, gamma_inv);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .p_over_gamma = p_over_gamma,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .h_ij_inv = 0,
    .det_h = 0,
    .model_id = GKYL_MODEL_SR,
    .mass = 1.0,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf);

  // test accuracy of the projection:
  struct gkyl_vlasov_lte_moments_inp inp_mom = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .p_over_gamma = p_over_gamma,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .h_ij_inv = 0,
    .det_h = 0,
    .model_id = GKYL_MODEL_SR,
    .mass = 1.0,
    .use_gpu = false,
  };
  gkyl_vlasov_lte_moments *lte_moms = gkyl_vlasov_lte_moments_inew( &inp_mom );
  gkyl_vlasov_lte_moments_advance(lte_moms, &local, &confLocal, distf, moms);
  gkyl_array_set_offset_range(m0, 1.0, moms, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m1i, 1.0, moms, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(m2, 1.0, moms, (vdim+1)*confBasis.num_basis, &confLocal);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {5.9020018022791720e-01, 1.8856465819367569e-17, 1.6811060851198739e-02,
    -3.1835607256007792e-18, -2.5559407924922751e-17, -1.6328957375267440e-02,
    3.5362637935871408e-17, -1.1626136441969583e-17};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]){1, 17}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
    }
  }
  
  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_mj_on_basis_test_1x1v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // release memory for moment data object
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

// special note, the p1 basis does not function
void test_1x1v_p2() { test_1x1v(2); }

void 
test_1x2v(int poly_order)
{
  double lower[] = {0.1, -15.0, -15.0}, upper[] = {1.0, 15.0, 15.0};
  int cells[] = {2, 16, 16};
  int vdim = 2, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2]}, velUpper[] = {upper[1], upper[2]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  moms = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_2v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_2v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);
  gkyl_array_set_offset_range(moms, 1.0, m0, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m1i, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m2, (vdim+1)*confBasis.num_basis, &confLocal);

  // build the p_over_gamma
  struct gkyl_array *p_over_gamma;
  p_over_gamma = mkarr(vdim * velBasis.num_basis, velLocal.volume);

  // build gamma
  struct gkyl_array *gamma;
  gamma = mkarr(velBasis.num_basis, velLocal.volume);

  // build gamma_inv
  struct gkyl_array *gamma_inv;
  gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);

  // Make GammaV2, GammaV, GammaV_inv
  gkyl_calc_sr_vars_init_p_vars(&vel_grid, &velBasis, &velLocal,
    p_over_gamma, gamma, gamma_inv);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .p_over_gamma = p_over_gamma,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .h_ij_inv = 0,
    .det_h = 0,
    .model_id = GKYL_MODEL_SR,
    .mass = 1.0,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf);

  // values to compare  at index (1, 9, 9) [remember, lower-left index is (1,1,1)]
  double p2_vals[] = {1.7020667884226476e-01, -7.7674914557148726e-18, -3.9516229859383111e-03,
    -2.8737491097032285e-02, 5.5726088991130445e-19, 1.1528499648312598e-18,
    8.1874847179781978e-03, -3.2664700781515490e-17, -1.0895575012022617e-02,
    -8.4007005283266260e-03, -5.6362918661943106e-19, 1.2361294855201663e-17,
    1.5787813110263945e-18, 1.6049196559163938e-17, 2.1063605116438087e-03,
    9.9876074849903799e-19, -6.4748916982029943e-04, -1.1513826111652373e-17,
    -9.2775152713167211e-19, 7.4400681776181471e-19};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]){1, 9, 9}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
    }
  }

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_mj_on_basis_test_1x2v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

void test_1x2v_p2() { test_1x2v(2); }

void 
test_1x3v(int poly_order)
{
  double lower[] = {0.1, -15.0, -15.0, -15.0}, upper[] = {1.0, 15.0, 15.0, 15.0};
  int cells[] = {2, 16, 16, 16};
  int vdim = 3, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  double velLower[] = {lower[1], lower[2], lower[3]}, velUpper[] = {upper[1], upper[2], upper[3]};
  int confCells[] = {cells[0]};
  int velCells[] = {cells[1], cells[2], cells[3]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vel_grid;
  gkyl_rect_grid_init(&vel_grid, vdim, velLower, velUpper, velCells);

  // velocity range
  int velGhost[] = {0, 0, 0};
  struct gkyl_range velLocal, velLocal_ext;
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = {confGhost[0], 0, 0, 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  moms = mkarr((vdim+2)*confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_3v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_3v, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);
  gkyl_array_set_offset_range(moms, 1.0, m0, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m1i, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms, 1.0, m2, (vdim+1)*confBasis.num_basis, &confLocal);

  // build the p_over_gamma
  struct gkyl_array *p_over_gamma;
  p_over_gamma = mkarr(vdim * velBasis.num_basis, velLocal.volume);

  // build gamma
  struct gkyl_array *gamma;
  gamma = mkarr(velBasis.num_basis, velLocal.volume);

  // build gamma_inv
  struct gkyl_array *gamma_inv;
  gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);

  // Make GammaV2, GammaV, GammaV_inv
  gkyl_calc_sr_vars_init_p_vars(&vel_grid, &velBasis, &velLocal,
    p_over_gamma, gamma, gamma_inv);

  // create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute LTE distribution
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_lte = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .vel_range = &velLocal,
    .p_over_gamma = p_over_gamma,
    .gamma = gamma,
    .gamma_inv = gamma_inv,
    .h_ij_inv = 0,
    .det_h = 0,
    .model_id = GKYL_MODEL_SR,
    .mass = 1.0,
    .use_gpu = false,
  };  
  gkyl_vlasov_lte_proj_on_basis *proj_lte = gkyl_vlasov_lte_proj_on_basis_inew(&inp_lte);
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms, distf);

  // values to compare  at index (1, 9, 9, 9) [remember, lower-left index is (1,1,1,1)]
  double p2_vals[] = {1.6326923473662415e-02, -2.7779798362812092e-19, -5.7251397678571113e-06,
    -2.9413756182605218e-03, -1.0882706908197464e-02, -9.1691545498145794e-20, -1.6291112189383686e-19,
    7.3891738935669528e-04, 2.6215943616563620e-19, 4.7510775105676502e-04, 2.3920085423587236e-03,
    -9.9802375797330017e-18, -1.2803270453557807e-03, -9.9127347429027999e-04, 2.5692430697696867e-03,
    6.1658843868991359e-20, -1.3264339369366439e-19, 1.0168790481577535e-19, -6.9285728763006371e-04,
    9.6703166952672561e-19, -1.2973734048803087e-19, 2.1959726150563245e-18, 2.8815110668801877e-04,
    5.0624860779664815e-20, -8.0565695192437930e-05, 6.1758274994024305e-18, 8.9364068681703317e-04,
    6.1884319642738651e-04, -1.4541027076081463e-19, -2.9745137044417144e-04, -7.3247363003543797e-04,
    5.0969448116179235e-21, -9.1037608148881686e-19, -2.0316355770111445e-19, -1.1045702861983222e-20,
    -1.6121338502579294e-18, -1.1349043120445892e-19, -1.1128896927091706e-18, -2.0910403428929313e-04,
    -1.0857925819912504e-19, 7.0221942561665685e-05, -3.5211619192312874e-20, 9.6248525403692007e-20,
    2.6937616824160315e-04, 5.4200306973849285e-19, -7.5794415473298147e-20, -1.8943801946011106e-20,
    6.8793916231015913e-22};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[4]){1, 9, 9, 9}));

  if (poly_order == 2) {
    for (int i = 0; i < basis.num_basis; ++i) {
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
    }
  }

  // write distribution function to file
  char fname[1024];
  sprintf(fname, "ctest_proj_mj_on_basis_test_1x3v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

void test_1x3v_p2() { test_1x3v(2); }

TEST_LIST = {
  {"test_1x1v_no_drift_p1", test_1x1v_no_drift_p1},
  {"test_1x1v_no_drift_p2", test_1x1v_no_drift_p2},
  {"test_1x1v_p2", test_1x1v_p2},
  {"test_1x2v_p2", test_1x2v_p2},
  {"test_1x3v_p2", test_1x3v_p2},
  {NULL, NULL},
};
