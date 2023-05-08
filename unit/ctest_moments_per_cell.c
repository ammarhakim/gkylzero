#include "gkyl_array.h"
#include "gkyl_util.h"
#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_mom_type.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_mj_moments.h>
#include <gkyl_proj_mj_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <math.h>

#include <gkyl_correct_mj.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size)
{
  struct gkyl_array *a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Projection functions for p/(gamma) = v in special relativistic systems
// Simplifies to p/sqrt(1 + p^2) where c = 1
static void
ev_p_over_gamma_1p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0] / sqrt(1.0 + xn[0] * xn[0]);
}
static void
ev_p_over_gamma_2p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0] / sqrt(1.0 + xn[0] * xn[0] + xn[1] * xn[1]);
  out[1] = xn[1] / sqrt(1.0 + xn[0] * xn[0] + xn[1] * xn[1]);
}
static void
ev_p_over_gamma_3p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = xn[0] / sqrt(1.0 + xn[0] * xn[0] + xn[1] * xn[1] + xn[2] * xn[2]);
  out[1] = xn[1] / sqrt(1.0 + xn[0] * xn[0] + xn[1] * xn[1] + xn[2] * xn[2]);
  out[2] = xn[2] / sqrt(1.0 + xn[0] * xn[0] + xn[1] * xn[1] + xn[2] * xn[2]);
}

static const evalf_t p_over_gamma_func[3] = {ev_p_over_gamma_1p, ev_p_over_gamma_2p, ev_p_over_gamma_3p};

// Projection functions for gamma = sqrt(1 + p^2) in special relativistic systems
static void
ev_gamma_1p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = sqrt(1.0 + xn[0] * xn[0]);
}
static void
ev_gamma_2p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = sqrt(1.0 + xn[0] * xn[0] + xn[1] * xn[1]);
}
static void
ev_gamma_3p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = sqrt(1.0 + xn[0] * xn[0] + xn[1] * xn[1] + xn[2] * xn[2]);
}

static const evalf_t gamma_func[3] = {ev_gamma_1p, ev_gamma_2p, ev_gamma_3p};

// Projection functions for gamma_inv = 1/sqrt(1 + p^2) in special relativistic systems
static void
ev_gamma_inv_1p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 1.0 / sqrt(1.0 + xn[0] * xn[0]);
}
static void
ev_gamma_inv_2p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 1.0 / sqrt(1.0 + xn[0] * xn[0] + xn[1] * xn[1]);
}
static void
ev_gamma_inv_3p(double t, const double *xn, double *out, void *ctx)
{
  out[0] = 1.0 / sqrt(1.0 + xn[0] * xn[0] + xn[1] * xn[1] + xn[2] * xn[2]);
}

static const evalf_t gamma_inv_func[3] = {ev_gamma_inv_1p, ev_gamma_inv_2p, ev_gamma_inv_3p};

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
  int cells[] = {2, 256}; // {16, 32, 64, 128}
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
  int poly_order_fixed = 2;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order_fixed);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order_fixed);
  gkyl_cart_modal_serendip(&velBasis, vdim, poly_order_fixed);

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

  // create moment arrays
  struct gkyl_array *m0, *m1i, *m2;
  m0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // Create a copy for comparison
  struct gkyl_array *m0_corr, *m1i_corr, *m2_corr;
  m0_corr = mkarr(confBasis.num_basis, confLocal_ext.volume);
  m1i_corr = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  m2_corr = mkarr(confBasis.num_basis, confLocal_ext.volume);

  // Create a copy for differences (d) and differences of differences (dd)
  struct gkyl_array *dm0, *dm1i, *dm2;
  dm0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  dm1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  dm2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  struct gkyl_array *ddm0, *ddm1i, *ddm2;
  ddm0 = mkarr(confBasis.num_basis, confLocal_ext.volume);
  ddm1i = mkarr(vdim * confBasis.num_basis, confLocal_ext.volume);
  ddm2 = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
                                                       poly_order + 1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1i = gkyl_proj_on_basis_new(&confGrid, &confBasis,
                                                        poly_order + 1, vdim, eval_M1i_1v, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
                                                       poly_order + 1, 1, eval_M2, NULL);

  gkyl_proj_on_basis *proj_m0_null = gkyl_proj_on_basis_new(&confGrid, &confBasis,
                                                            poly_order + 1, 1, eval_M0_null, NULL);
  gkyl_proj_on_basis *proj_m1i_null = gkyl_proj_on_basis_new(&confGrid, &confBasis,
                                                             poly_order + 1, vdim, eval_M1i_1v_null, NULL);
  gkyl_proj_on_basis *proj_m2_null = gkyl_proj_on_basis_new(&confGrid, &confBasis,
                                                            poly_order + 1, 1, eval_M2_null, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2);

  // create a copy for the correct intial value
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_corr);
  gkyl_proj_on_basis_advance(proj_m1i, 0.0, &confLocal, m1i_corr);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_corr);

  // create a copy for the differences & differences of differences
  gkyl_proj_on_basis_advance(proj_m0_null, 0.0, &confLocal, dm0);
  gkyl_proj_on_basis_advance(proj_m1i_null, 0.0, &confLocal, dm1i);
  gkyl_proj_on_basis_advance(proj_m2_null, 0.0, &confLocal, dm2);

  gkyl_proj_on_basis_advance(proj_m0_null, 0.0, &confLocal, ddm0);
  gkyl_proj_on_basis_advance(proj_m1i_null, 0.0, &confLocal, ddm1i);
  gkyl_proj_on_basis_advance(proj_m2_null, 0.0, &confLocal, ddm2);

  // build the p_over_gamma
  struct gkyl_array *p_over_gamma;
  p_over_gamma = mkarr(vdim * velBasis.num_basis, velLocal.volume);
  gkyl_proj_on_basis *p_over_gamma_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &vel_grid,
      .basis = &velBasis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = 8,
      .num_ret_vals = vdim,
      .eval = p_over_gamma_func[vdim - 1], // ev_p_over_gamma_1p,
      .ctx = 0});
  gkyl_proj_on_basis_advance(p_over_gamma_proj, 0.0, &velLocal, p_over_gamma);

  // build gamma
  struct gkyl_array *gamma;
  gamma = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_proj_on_basis *gamma_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &vel_grid,
      .basis = &velBasis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = 8,
      .num_ret_vals = 1,
      .eval = gamma_func[vdim - 1],
      .ctx = 0});
  gkyl_proj_on_basis_advance(gamma_proj, 0.0, &velLocal, gamma);

  // build gamma_inv
  struct gkyl_array *gamma_inv;
  gamma_inv = mkarr(velBasis.num_basis, velLocal.volume);
  gkyl_proj_on_basis *gamma_inv_proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
      .grid = &vel_grid,
      .basis = &velBasis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = 8,
      .num_ret_vals = 1,
      .eval = gamma_inv_func[vdim - 1],
      .ctx = 0});
  gkyl_proj_on_basis_advance(gamma_inv_proj, 0.0, &velLocal, gamma_inv);

  // create distribution function array
  struct gkyl_array *distf_mj;
  distf_mj = mkarr(basis.num_basis, local_ext.volume);

  // projection updater to compute mj
  gkyl_proj_mj_on_basis *proj_mj = gkyl_proj_mj_on_basis_new(&grid,
                                                             &confBasis, &basis, poly_order + 1);

  // 1. Project the MJ with the intially correct moments
  gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(proj_mj, &local, &confLocal, m0_corr, m1i_corr, m2_corr, distf_mj);

  // 1. (NEW) For moments in a cell, and moment object (all moments)
  struct gkyl_array *m;
  m = mkarr(vdim + 2, local_ext.volume);
  struct gkyl_mom_type *sr_int_mom = gkyl_int_mom_vlasov_sr_new(&confBasis, &basis, &confLocal, &velLocal, 0);
  struct gkyl_array *GammaV2;
  struct gkyl_array *Gamma_inv;
  GammaV2 = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  Gamma_inv = gkyl_array_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);

  // Convergence scheme
  for (int i = 0; i < 1; ++i)
  { // 3000

    // printf("\n----------- ************************* ---------\n");
    // printf("----------- Begining iterative Loop: T = %d ---------\n", i);
    // printf("----------- ************************* ---------\n\n");

    // write distribution function to file
    char fname[1024];
    if ((i % 100) == 0)
    {
      sprintf(fname, "ctest_moments_per_cell_f_1x1v_p%d_iteration_%03d.gkyl", poly_order, i);
      gkyl_grid_sub_array_write(&grid, &local, distf_mj, fname);
    }

    // 2. Calculate the new moments
    // calculate the moments of the dist (n, vb, T -> m0, m1i, m2)
    gkyl_mj_moments *mj_moms = gkyl_mj_moments_new(&grid, &confBasis, &basis, &confLocal, &velLocal, confLocal.volume, confLocal_ext.volume, false);
    gkyl_mj_moments_advance(mj_moms, p_over_gamma, gamma, gamma_inv, distf_mj, m0, m1i, m2, &local, &confLocal);

    // 2. (NEW) Calculate the moments / set aux fields ( already calculated in mj_moms)
    // Set Aux Fields
    // (Gamma^2 = 1/(1-vb^2)) calculate
    gkyl_calc_sr_vars_Gamma2(&confBasis, &basis,
                             &confLocal, m1i, GammaV2);

    // (Gamma_inv = sqrt(1-vb^2))
    gkyl_calc_sr_vars_Gamma_inv(&confBasis, &basis,
                                &confLocal, m1i, Gamma_inv);
    gkyl_mom_vlasov_sr_set_auxfields(sr_int_mom,
                                     (struct gkyl_mom_vlasov_sr_auxfields){.p_over_gamma = p_over_gamma,
                                                                           .gamma = gamma,
                                                                           .gamma_inv = gamma_inv,
                                                                           .V_drift = m1i,
                                                                           .GammaV2 = GammaV2,
                                                                           .GammaV_inv = Gamma_inv});

    // 3. (NEW) Calc the moments in all cells
    struct gkyl_range_iter phase_iter;
    double xc[ndim];
    gkyl_range_iter_init(&phase_iter, &local);
    while (gkyl_range_iter_next(&phase_iter))
    {
      long pidx = gkyl_range_idx(&local, phase_iter.idx);
      gkyl_rect_grid_cell_center(&grid, phase_iter.idx, xc);
      gkyl_mom_type_calc(sr_int_mom, xc, grid.dx, phase_iter.idx,
                         gkyl_array_cfetch(distf_mj, pidx), gkyl_array_fetch(m, pidx), 0);
    }

    // 4. (NEW) Print the moments per cell
    sprintf(fname, "ctest_moments_per_cell_moments_1x1v_nquad%d_NV_%03d.gkyl", (poly_order + 1), cells[1]);
    gkyl_grid_sub_array_write(&grid, &local, m, fname);

    // a. Calculate  ddMi^(k+1) =  Mi_corr - Mi_new
    // ddm0 = m0_corr - m0;
    //  Compute out = out + a*inp. Returns out.
    gkyl_array_clear_range(ddm0, 0.0, confLocal);
    gkyl_array_accumulate_range(ddm0, -1.0, m0, confLocal);
    gkyl_array_accumulate_range(ddm0, 1.0, m0_corr, confLocal);
    gkyl_array_clear_range(ddm1i, 0.0, confLocal);
    gkyl_array_accumulate_range(ddm1i, -1.0, m1i, confLocal);
    gkyl_array_accumulate_range(ddm1i, 1.0, m1i_corr, confLocal);
    gkyl_array_clear_range(ddm2, 0.0, confLocal);
    gkyl_array_accumulate_range(ddm2, -1.0, m2, confLocal);
    gkyl_array_accumulate_range(ddm2, 1.0, m2_corr, confLocal);

    // b. Calculate  dMi^(k+1) = dM0^k + ddMi^(k+1) | where dM0^0 = 0
    // dm_new = dm_old + ddm0;
    gkyl_array_accumulate_range(dm0, 1.0, ddm0, confLocal);
    gkyl_array_accumulate_range(dm1i, 1.0, ddm1i, confLocal);
    gkyl_array_accumulate_range(dm2, 1.0, ddm2, confLocal);

    // 4. Diagnostic ouputs
    if ((i % 100) == 0)
    {
      struct gkyl_range_iter biter;
      gkyl_range_iter_init(&biter, &confLocal);
      while (gkyl_range_iter_next(&biter))
      {
        long midx = gkyl_range_idx(&confLocal, biter.idx);
        const double *m0_corr_local = gkyl_array_cfetch(m0_corr, midx);
        const double *m0_local = gkyl_array_cfetch(m0, midx);
        const double *dm0_local = gkyl_array_cfetch(dm0, midx);
        const double *ddm0_local = gkyl_array_cfetch(ddm0, midx);
        const double *m1i_corr_local = gkyl_array_cfetch(m1i_corr, midx);
        const double *m1i_local = gkyl_array_cfetch(m1i, midx);
        const double *dm1i_local = gkyl_array_cfetch(dm1i, midx);
        const double *ddm1i_local = gkyl_array_cfetch(ddm1i, midx);
        const double *m2_corr_local = gkyl_array_cfetch(m2_corr, midx);
        const double *m2_local = gkyl_array_cfetch(m2, midx);
        const double *dm2_local = gkyl_array_cfetch(dm2, midx);
        const double *ddm2_local = gkyl_array_cfetch(ddm2, midx);
        printf("\n------- n interation : %d ------\n", i);
        printf("n_corr: %g\n", m0_corr_local[0]);
        printf("n: %g\n", m0_local[0]);
        printf("dn: %g\n", dm0_local[0]);
        printf("ddn: %g\n", ddm0_local[0]);
        printf("Diff (n - n_corr): %g\n", (m0_local[0] - m0_corr_local[0]));
        printf("------- vb interation : %d ------\n", i);
        printf("vb_corr: %g\n", m1i_corr_local[0]);
        printf("vb: %g\n", m1i_local[0]);
        printf("dvb: %g\n", dm1i_local[0]);
        printf("ddvb: %g\n", ddm1i_local[0]);
        printf("Diff (vb - vb_corr): %g\n", (m1i_local[0] - m1i_corr_local[0]));
        printf("------- T interation : %d ------\n", i);
        printf("m2_corr: %g\n", m2_corr_local[0]);
        printf("m2: %g\n", m2_local[0]);
        printf("dm2: %g\n", dm2_local[0]);
        printf("ddm2: %g\n", ddm2_local[0]);
        printf("Diff (vb - vb_corr): %g\n", (m2_local[0] - m2_corr_local[0]));
      }
    }

    // c. Calculate  M0^(k+1) = m0 + dm^(k+1)
    // m0 = m0_corr + dm_new;
    gkyl_array_clear_range(m0, 0.0, confLocal);
    gkyl_array_accumulate_range(m0, 1.0, m0_corr, confLocal);
    gkyl_array_accumulate_range(m0, 1.0, dm0, confLocal);
    gkyl_array_clear_range(m1i, 0.0, confLocal);
    gkyl_array_accumulate_range(m1i, 1.0, m1i_corr, confLocal);
    gkyl_array_accumulate_range(m1i, 1.0, dm1i, confLocal);
    gkyl_array_clear_range(m2, 0.0, confLocal);
    gkyl_array_accumulate_range(m2, 1.0, m2_corr, confLocal);
    gkyl_array_accumulate_range(m2, 1.0, dm2, confLocal);

    // 3. Update the dist_mj using the corrected moments
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(proj_mj, &local, &confLocal, m0, m1i, m2, distf_mj);

    // Release the memory
    gkyl_mj_moments_release(mj_moms);
  }
  // end iteration loop

  // values to compare  at index (1, 17) [remember, lower-left index is (1,1)]
  double p2_vals[] = {0.4106556323526475, -8.940762710879627e-17,
                      0.06572788982671821, 8.645045365809577e-18, -5.979556483302724e-17,
                      -0.001036545017544019, 2.229425706102836e-17, 2.764128755933108e-17};

  const double *fv = gkyl_array_cfetch(distf_mj, gkyl_range_idx(&local_ext, (int[2]){1, 16}));

  if (poly_order == 2)
  {
    for (int i = 0; i < basis.num_basis; ++i)
    {
      // printf("fv[%d] = %1.16g\n",i,fv[i]);
      TEST_CHECK(gkyl_compare_double(p2_vals[i], fv[i], 1e-12));
    }
  }

  // release memory for moment data object
  gkyl_proj_on_basis_release(p_over_gamma_proj);
  gkyl_proj_on_basis_release(gamma_proj);
  gkyl_proj_on_basis_release(gamma_inv_proj);
  gkyl_array_release(m0);
  gkyl_array_release(m1i);
  gkyl_array_release(m2);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1i_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(dm0);
  gkyl_array_release(dm1i);
  gkyl_array_release(dm2);
  gkyl_array_release(ddm0);
  gkyl_array_release(ddm1i);
  gkyl_array_release(ddm2);
  gkyl_array_release(distf_mj);
  gkyl_proj_mj_on_basis_release(proj_mj);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1i);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_proj_on_basis_release(proj_m0_null);
  gkyl_proj_on_basis_release(proj_m1i_null);
  gkyl_proj_on_basis_release(proj_m2_null);
  gkyl_array_release(p_over_gamma);
  // 5. (NEW) Release the new data
  // gkyl_int_mom_vlasov_sr_release(sr_int_mom);
  gkyl_array_release(m);
  gkyl_array_release(GammaV2);
  gkyl_array_release(Gamma_inv);
}

// special note, the p1 basis does not function
void test_1x1v_p2() { test_1x1v(2); }

TEST_LIST = {
    {"test_1x1v_p2", test_1x1v_p2},
    {NULL, NULL},
};
