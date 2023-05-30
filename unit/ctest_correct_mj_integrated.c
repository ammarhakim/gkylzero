#include "gkyl_array.h"
#include "gkyl_util.h"
#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_correct_mj.h>
#include <gkyl_mj_moments.h>
#include <gkyl_proj_mj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <math.h>

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges
{
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

  for (int d = 0; d < ndim; ++d)
  {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
                           d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
                           d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void eval_M0(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void eval_M0_null(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0;
}

void eval_M1i_1v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5; // 0.5;
}

void eval_M1i_1v_null(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0; // 0.5;
}

void eval_M2(double t, const double *xn, double *restrict fout, void *ctx)
{
  double T = 1.0;
  double x = xn[0];
  fout[0] = T;
}

void eval_M2_null(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0;
}

void eval_M1i_2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
  fout[1] = 0.25;
}

void eval_M1i_2v_null(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0;
  fout[1] = 0.0;
}

void eval_M1i_3v(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  // case : 1
  fout[0] = 0.5;
  fout[1] = 0.5;
  fout[2] = 0.5;
  // fout[0] = 0.5; fout[1] = 0.25; fout[2] = -0.5;
}

void eval_M1i_3v_null(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = 0.0;
}

void test_1x1v(int poly_order)
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
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);

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
  struct gkyl_array *distf_mj;
  distf_mj = mkarr(basis.num_basis, local_ext.volume);

  // Create a MJ with corrected moments
  gkyl_correct_mj *corr_mj = gkyl_correct_mj_new(&grid, &confBasis, &basis, &confLocal, &velLocal, confLocal.volume, confLocal_ext.volume, false);
  // gkyl_proj_mj_on_basis *proj_mj = gkyl_proj_mj_on_basis_new(&grid,
  //                                                            &confBasis, &basis, poly_order + 1);
  // gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(proj_mj, &local, &confLocal, m0_corr, m1i_corr, m2_corr, distf_mj);
  // gkyl_correct_mj_fix_m0(corr_mj, p_over_gamma, distf_mj, m0_corr, m1i_corr, &local, &confLocal);
  gkyl_correct_mj_fix(corr_mj, distf_mj, m0_corr, m1i_corr, m2_corr, &local, &confLocal, poly_order, &confLocal_ext, &velLocal, &velBasis, &vel_grid);
  gkyl_correct_mj_release(corr_mj);

  // Write the output
  char fname[1024];
  sprintf(fname, "ctest_correct_mj_integrated_1x1v_p%d_with_fix_x10.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, distf_mj, fname);

  // Correct the distribution function
  gkyl_mj_moments *mj_moms = gkyl_mj_moments_new(&grid, &confBasis, &basis, &confLocal, &velLocal, confLocal.volume, confLocal_ext.volume, false);
  gkyl_mj_moments_advance(mj_moms, p_over_gamma, gamma, gamma_inv, distf_mj, m0, m1i, m2, &local, &confLocal);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {0.4106556323526475, -8.940762710879627e-17,
                      0.06572788982671821, 8.645045365809577e-18, -5.979556483302724e-17,
                      -0.001036545017544019, 2.229425706102836e-17, 2.764128755933108e-17};

  const double *fv = gkyl_array_cfetch(distf_mj, gkyl_range_idx(&local_ext, (int[2]){1, 16}));

  if (poly_order == 2)
  {
    for (int i = 0; i < basis.num_basis; ++i)
    {
      // printf("fv[%d] = %1.16g\n", i, fv[i]);
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
    }
  }

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(distf_mj);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
}

void test_1x2v(int poly_order)
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
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = {confGhost[0], 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);

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
  struct gkyl_array *distf_mj;
  distf_mj = mkarr(basis.num_basis, local_ext.volume);

  // Create a MJ with corrected moments
  gkyl_correct_mj *corr_mj = gkyl_correct_mj_new(&grid, &confBasis, &basis, &confLocal, &velLocal, confLocal.volume, confLocal_ext.volume, false);
  gkyl_correct_mj_fix(corr_mj, distf_mj, m0_corr, m1i_corr, m2_corr, &local, &confLocal, poly_order, &confLocal_ext, &velLocal, &velBasis, &vel_grid);
  gkyl_correct_mj_release(corr_mj);

  // Correct the distribution function
  gkyl_mj_moments *mj_moms = gkyl_mj_moments_new(&grid, &confBasis, &basis, &confLocal, &velLocal, confLocal.volume, confLocal_ext.volume, false);
  gkyl_mj_moments_advance(mj_moms, p_over_gamma, gamma, gamma_inv, distf_mj, m0, m1i, m2, &local, &confLocal);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {0.1196584827807841, -3.488028281807569e-18,
                      0.01964687504797331, 0.01333312935386793, 3.040820909981071e-18,
                      1.049969351462672e-18, 0.0024642991346404, -5.334926505468131e-18,
                      -0.0002280262498167821, -0.001033628149770621, -1.153462471614103e-18,
                      -8.381126915873462e-18, -1.442233011966233e-18, 8.098746105906206e-18,
                      -4.291622780592622e-05, 2.924904640105737e-19, -0.0002003167531638971,
                      2.378358894760202e-18, -2.237263192050087e-19, -2.237263192050087e-19};

  const double *fv = gkyl_array_cfetch(distf_mj, gkyl_range_idx(&local_ext, (int[3]){1, 16, 16}));

  if (poly_order == 2)
  {
    for (int i = 0; i < basis.num_basis; ++i)
    {
      // printf("%1.16g,\n", fv[i]);
      // printf("fv[%d] = %1.16g\n", i, fv[i]);
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
    }
  }

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(distf_mj);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
}

void test_1x3v(int poly_order)
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
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&vel_grid, velGhost, &velLocal_ext, &velLocal);

  // basis functions
  struct gkyl_basis basis, confBasis, velBasis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order);

  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = {confGhost[0], 0, 0, 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr;
  m0_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *m0, *m1i, *m2;
  m0 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  m1i = gkyl_array_new(GKYL_DOUBLE, vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);

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
  struct gkyl_array *distf_mj;
  distf_mj = mkarr(basis.num_basis, local_ext.volume);

  // Create a MJ with corrected moments
  gkyl_correct_mj *corr_mj = gkyl_correct_mj_new(&grid, &confBasis, &basis, &confLocal, &velLocal, confLocal.volume, confLocal_ext.volume, false);
  gkyl_correct_mj_fix(corr_mj, distf_mj, m0_corr, m1i_corr, m2_corr, &local, &confLocal, poly_order, &confLocal_ext, &velLocal, &velBasis, &vel_grid);
  gkyl_correct_mj_release(corr_mj);

  // Correct the distribution function
  gkyl_mj_moments *mj_moms = gkyl_mj_moments_new(&grid, &confBasis, &basis, &confLocal, &velLocal, confLocal.volume, confLocal_ext.volume, false);
  gkyl_mj_moments_advance(mj_moms, p_over_gamma, gamma, gamma_inv, distf_mj, m0, m1i, m2, &local, &confLocal);

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {0.002127623222952959, -1.631506586025315e-19, 0.001142561434832162,
                      0.00114256143483216, 0.00114256143483216, -2.51371116629004e-19,
                      -2.01075727593333e-19, 0.0006425084967736795, -2.595001694861082e-19,
                      0.0006425084967736792, 0.0006425084967736784, -3.510815947551172e-18,
                      0.0001864729154912117, 0.000186472915491211, 0.0001864729154912109,
                      4.028479705065517e-20, -1.268594757705666e-19, -2.111921481277014e-19,
                      0.000380214364112294, -1.816029978707947e-18, 1.480556426471607e-19,
                      -1.941768451297124e-18, 0.0001121226604556045, -2.797821897768886e-20,
                      0.0001121226604556044, -1.816029978707946e-18, 0.0001121226604556046,
                      0.0001121226604556042, -5.757070789733834e-20, 0.0001121226604556041,
                      0.0001121226604556039, -8.190440447282892e-20, -8.088564395679194e-19,
                      8.497928576027658e-20, -1.5611492311066e-20, -7.504319976751439e-19,
                      9.536202206769661e-21, -9.434845009607244e-19, 7.168010462749378e-05,
                      -9.105457586457297e-20, 7.168010462749364e-05, -5.777782848963335e-20,
                      -1.332209120431403e-19, 7.168010462749348e-05, -4.780821736467332e-19,
                      6.948823307748629e-20, -5.954850476020656e-21, -8.13979340295276e-20};

  const double *fv = gkyl_array_cfetch(distf_mj, gkyl_range_idx(&local_ext, (int[4]){1, 8, 8, 8}));

  if (poly_order == 2)
  {
    for (int i = 0; i < basis.num_basis; ++i)
    {
      // printf("%1.16g,\n", fv[i]);
      //  printf("fv[%d] = %1.16g\n", i, fv[i]);
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
    }
  }

  // release memory for moment data object
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(distf_mj);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_array_release(p_over_gamma);
}

// special note, the p1 basis does not function
void test_1x1v_p2() { test_1x1v(2); }
void test_1x2v_p2() { test_1x2v(2); }
void test_1x3v_p2() { test_1x3v(2); }

TEST_LIST = {
    {"test_1x1v_p2", test_1x1v_p2},
    {"test_1x2v_p2", test_1x2v_p2},
    {"test_1x3v_p2", test_1x3v_p2},
    {NULL, NULL},
};
