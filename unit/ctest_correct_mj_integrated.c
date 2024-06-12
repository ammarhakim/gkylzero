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
eval_M1i_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; // 0.5;
}

void 
eval_M2(double t, const double *xn, double *restrict fout, void *ctx)
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
eval_M1i_3v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; fout[1] = 0.5; fout[2] = 0.5;
}

void
eval_M0_x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.1 + sin(xn[0]);
}

void 
eval_M1i_1v_x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5*(sin(xn[0])); // 0.5;
}

void 
eval_M2_x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T*(1.1 + sin(xn[0]));
}

void 
test_1x1v(int poly_order)
{
  double lower[] = {0.1, -10.0}, upper[] = {1.0, 10.0};
  int cells[] = {2, 32}; // 1001
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

  // basis functions
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

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, NULL);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

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
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
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
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m) 
  gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Correct the distribution function
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
  double p2_vals[] = {0.4106556323526475, -8.940762710879627e-17,
    0.06572788982671821, 8.645045365809577e-18, -5.979556483302724e-17,
    -0.001036545017544019, 2.229425706102836e-17, 2.764128755933108e-17};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]){1, 16}));

  if (poly_order == 2)
    for (int i = 0; i < basis.num_basis; ++i)
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}


void 
test_1x1v_spatially_varied(int poly_order)
{
  double lower[] = {0.0, -10.0}, upper[] = {10.0, 10.0};
  int cells[] = {10, 32}; // 1001
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

  // basis functions
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

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0_x, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_1v_x, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2_x, NULL);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

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
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
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
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m)
  gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Correct the distribution function
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

  // Write the output (moments)
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_n_corr.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m0_corr,fname);
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_vb_corr.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m1i_corr,fname);
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_T_corr.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m2_corr,fname);


  // Write the output (moments)
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_n.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m0,fname);
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_vb.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m1i,fname);
  //sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_x_T.gkyl", poly_order);
  //gkyl_grid_sub_array_write(&confGrid,&confLocal,0,m2,fname);


  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double m0_vals[] = {2.2057459659974716e+00,3.4938172762332936e-01,-2.4953670224477906e-02};
  double m1i_vals[] = {3.2505552369353341e-01,1.7469086381166479e-01,-1.2476835112239079e-02};
  double m2_vals[] = {2.2057459659974712e+00,3.4938172762332981e-01,-2.4953670224477753e-02};


  //const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[2]){1, 16}));
  const double *m0_fixed = gkyl_array_cfetch(m0, gkyl_range_idx(&confLocal_ext, (int[1]){1}));
  const double *m1i_fixed = gkyl_array_cfetch(m1i, gkyl_range_idx(&confLocal_ext, (int[1]){1}));
  const double *m2_fixed = gkyl_array_cfetch(m2, gkyl_range_idx(&confLocal_ext, (int[1]){1}));

  if (poly_order == 2)
    for (int i = 0; i < confBasis.num_basis; ++i){
      TEST_CHECK(gkyl_compare_double(m0_vals[i], m0_fixed[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(m1i_vals[i], m1i_fixed[i], 1e-12));
      TEST_CHECK(gkyl_compare_double(m2_vals[i], m2_fixed[i], 1e-12));
    }

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}


void 
test_1x2v(int poly_order)
{
  double lower[] = {0.1, -10.0, -10.0}, upper[] = {1.0, 10.0, 10.0};
  int cells[] = {2, 32, 32};
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

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_2v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, NULL);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

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
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
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
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m)
  gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Correct the distribution function
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

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x2v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {0.1196584827807841, -3.488028281807569e-18,
    0.01964687504797331, 0.01333312935386793, 3.040820909981071e-18,
    1.049969351462672e-18, 0.0024642991346404, -5.334926505468131e-18,
    -0.0002280262498167821, -0.001033628149770621, -1.153462471614103e-18,
    -8.381126915873462e-18, -1.442233011966233e-18, 8.098746105906206e-18,
    -4.291622780592622e-05, 2.924904640105737e-19, -0.0002003167531638971,
    2.378358894760202e-18, -2.237263192050087e-19, -2.237263192050087e-19};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[3]){1, 16, 16}));

  if (poly_order == 2)
    for (int i = 0; i < basis.num_basis; ++i)
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_moments_release(lte_moms);
  gkyl_vlasov_lte_proj_on_basis_release(proj_lte);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

void 
test_1x3v(int poly_order)
{
  double lower[] = {0.1, -10.0, -10.0, -10.0}, upper[] = {1.0, 10.0, 10.0, 10.0};
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

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr, *moms_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms_corr = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2, *moms;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  moms = gkyl_array_new(GKYL_DOUBLE, (vdim+2) * confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, vdim, eval_M1i_3v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order + 1, 1, eval_M2, NULL);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);
  gkyl_array_set_offset_range(moms_corr, 1.0, m0_corr, 0*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m1i_corr, 1*confBasis.num_basis, &confLocal);
  gkyl_array_set_offset_range(moms_corr, 1.0, m2_corr, (vdim+1)*confBasis.num_basis, &confLocal);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

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
  // Project LTE distribution function (and correct its density internally)
  gkyl_vlasov_lte_proj_on_basis_advance(proj_lte, &local, &confLocal, moms_corr, distf);

  // Create a MJ with corrected moments
  struct gkyl_vlasov_lte_correct_inp inp_corr = {
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
    .max_iter = 100,
    .eps = 1e-12,
  };
  gkyl_vlasov_lte_correct *corr_mj = gkyl_vlasov_lte_correct_inew( &inp_corr );
  // Correct the other moments (V_drift, T/m)
  gkyl_vlasov_lte_correct_all_moments(corr_mj, distf, moms_corr, &local, &confLocal);
  gkyl_vlasov_lte_correct_release(corr_mj);

  // Correct the distribution function
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

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x3v_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {0.002127623222951445, -2.162483920741371e-19, 0.001142561434831357,
    0.001142561434831355, 0.001142561434831355, 1.764595896820067e-20,
    -8.294481910314426e-20, 0.0006425084967732312, 9.516906111084089e-21,
    0.000642508496773231, 0.0006425084967732302, -3.192214167141901e-18,
    0.000186472915491083, 0.0001864729154910819, 0.0001864729154910819,
    -5.325792575489166e-20, 8.137013563789262e-20, 4.733285231645305e-20,
    0.000380214364112031, -1.642788966842025e-18, 1.704766833543138e-19,
    -1.642788966842023e-18, 0.0001121226604555277, 1.95905162472965e-20,
    0.0001121226604555275, -1.906839759279303e-18, 0.0001121226604555277,
    0.0001121226604555272, 6.544111088116724e-20, 0.0001121226604555273,
    0.000112122660455527, 1.685634497069759e-20, -5.335890207604426e-19,
    1.154778075974092e-19, -3.540835950960802e-20, -7.517892185638519e-19,
    1.657731966330816e-19, -7.688078602245709e-19, 7.168010462744494e-05,
    -1.026066499177183e-20, 7.168010462744473e-05, 7.331147141885321e-20,
    7.33114714188532e-20, 7.168010462744468e-05, -3.388438170490261e-19,
    3.636875421167469e-20, -3.907432934183393e-20, 3.636875421167467e-20};

  const double *fv = gkyl_array_cfetch(distf, gkyl_range_idx(&local_ext, (int[4]){1, 8, 8, 8}));

  if (poly_order == 2)
    for (int i = 0; i < basis.num_basis; ++i)
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(moms);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(moms_corr);
  gkyl_array_release(distf);
  gkyl_vlasov_lte_moments_release(lte_moms);
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
void test_1x1v_p2_spatially_varied() { test_1x1v_spatially_varied(2); }
void test_1x2v_p2() { test_1x2v(2); }
void test_1x3v_p2() { test_1x3v(2); }

TEST_LIST = {
  {"test_1x1v_p2", test_1x1v_p2},
  {"test_1x1v_p2_spatially_varied", test_1x1v_p2_spatially_varied},
  {"test_1x2v_p2", test_1x2v_p2},
  {"test_1x3v_p2", test_1x3v_p2},
  {NULL, NULL},
};
