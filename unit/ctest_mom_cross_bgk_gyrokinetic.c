#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_cross_bgk.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <math.h>
#include <stdio.h>

// Allocate cu_dev array
static struct gkyl_array*
mkarr(long nc, long size, bool use_gpu)
{
  struct gkyl_array* a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// Create ghost ranges
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

  for (int d = 0; d < ndim; ++d){
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
                           d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
                           d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void eval_M0_e(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 10.0;
}
void eval_M0_i(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 10.0;
}
void eval_M1_e(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1];
  double udrift = 0.0;
  eval_M0_e(t, xn, den, ctx);
  fout[0] = den[0]*udrift;
}
void eval_M1_i(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1];
  double udrift = 1.0;
  eval_M0_i(t, xn, den, ctx);
  fout[0] = den[0]*udrift;
}
void eval_M2_e(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1];
  double udrift = 0.0, vtsq = 100.0;
  eval_M0_e(t, xn, den, ctx);
  fout[0] = den[0]*vtsq + den[0]*udrift*udrift;
}
void eval_M2_i(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1];
  double udrift = 1.0, vtsq = 36.0;
  eval_M0_i(t, xn, den, ctx);
  fout[0] = den[0]*vtsq + den[0]*udrift*udrift;
}

void test_1x1v(int poly_order, bool use_gpu)
{
  double me = 1.0; //9.10938215e-31;
  double mi = 2.0; //1.672621637e-27;
  double beta = 0.0;
  double nu_ei = 10.0; 
  double nu_ie = 5.0;

  double lower[] = {-0.5, -5.0}, upper[] = {0.5, 5.0};
  int cells[] = {2, 32};
  int vdim = 1, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // Grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // Basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order==1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Configuration space range
  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost;
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  // Phase space range
  int ghost[] = {confGhost[0], 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost;
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create moment arrays
  // (1) electron
  struct gkyl_array *m0_e_ho, *m1_e_ho, *m2_e_ho;
  m0_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0_e, NULL);
  gkyl_proj_on_basis *proj_m1_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1_e, NULL);
  gkyl_proj_on_basis *proj_m2_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2_e, NULL);
  gkyl_proj_on_basis_advance(proj_m0_e, 0.0, &confLocal, m0_e_ho);
  gkyl_proj_on_basis_advance(proj_m1_e, 0.0, &confLocal, m1_e_ho);
  gkyl_proj_on_basis_advance(proj_m2_e, 0.0, &confLocal, m2_e_ho);
  struct gkyl_array *m0_e, *m1_e, *m2_e;
  if (use_gpu) {
    m0_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_e, m0_e_ho);
    gkyl_array_copy(m1_e, m1_e_ho);
    gkyl_array_copy(m2_e, m2_e_ho);
  } else {
    m0_e = m0_e_ho;
    m1_e = m1_e_ho;
    m2_e = m2_e_ho;
  }
  struct gkyl_array *moms_e = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_e, 1., m0_e, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_e, 1., m1_e, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_e, 1., m2_e, 2*confBasis.num_basis);
  // (2) ion
  struct gkyl_array *m0_i_ho, *m1_i_ho, *m2_i_ho;
  m0_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0_i, NULL);
  gkyl_proj_on_basis *proj_m1_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1_i, NULL);
  gkyl_proj_on_basis *proj_m2_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2_i, NULL);
  gkyl_proj_on_basis_advance(proj_m0_i, 0.0, &confLocal, m0_i_ho);
  gkyl_proj_on_basis_advance(proj_m1_i, 0.0, &confLocal, m1_i_ho);
  gkyl_proj_on_basis_advance(proj_m2_i, 0.0, &confLocal, m2_i_ho);
  struct gkyl_array *m0_i, *m1_i, *m2_i;
  if (use_gpu) {
    m0_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_i, m0_i_ho);
    gkyl_array_copy(m1_i, m1_i_ho);
    gkyl_array_copy(m2_i, m2_i_ho);
  } else {
    m0_i = m0_i_ho;
    m1_i = m1_i_ho;
    m2_i = m2_i_ho;
  }
  struct gkyl_array *moms_i = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_i, 1., m0_i, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_i, 1., m1_i, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_i, 1., m2_i, 2*confBasis.num_basis);
  
  // Calculate the cross moments
  gkyl_mom_cross_bgk_gyrokinetic *momCrossCalc = gkyl_mom_cross_bgk_gyrokinetic_new(&basis, &confBasis);
  struct gkyl_array *moms_cross = mkarr(6*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_mom_cross_bgk_gyrokinetic_advance(momCrossCalc, &confLocal_ext, beta, me, moms_e, mi, moms_i, nu_ei, nu_ie, moms_cross);
  gkyl_mom_cross_bgk_gyrokinetic_release(momCrossCalc);

  // Write out on host
  char fname[1024];
  sprintf(fname, "ctest_mom_cross_bgk_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, moms_cross, fname);

  // Compare with the expected cross moments
  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&confLocal, idx);
    const double *momsCross = gkyl_array_cfetch(moms_cross, linidx);
    TEST_CHECK( gkyl_compare(10.0, momsCross[0*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(10.0, momsCross[1*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(5.0, momsCross[2*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(5.0, momsCross[3*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(913.3333333333334, momsCross[4*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(456.6666666666667, momsCross[5*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_MSG("Produced: %.13e, \t%.13e, \t%.13e, %.13e, \t%.13e, \t%.13e", momsCross[0*confBasis.num_basis]/sqrt(2), momsCross[1*confBasis.num_basis]/sqrt(2), momsCross[2*confBasis.num_basis]/sqrt(2), momsCross[3*confBasis.num_basis]/sqrt(2), momsCross[4*confBasis.num_basis]/sqrt(2), momsCross[5*confBasis.num_basis]/sqrt(2));
  } // The hard coded numbers are expected values. The basis is looked up in maxima with: load("basis-precalc/basisSer1x"); polyOrder:1$ basis:basisC[polyOrder];

  // Release memory for moment data object
  gkyl_array_release(m0_e_ho);
  gkyl_array_release(m1_e_ho);
  gkyl_array_release(m2_e_ho);
  gkyl_array_release(m0_i_ho);
  gkyl_array_release(m1_i_ho);
  gkyl_array_release(m2_i_ho);
  if (use_gpu) {
    gkyl_array_release(m0_e);
    gkyl_array_release(m1_e);
    gkyl_array_release(m2_e);
    gkyl_array_release(m0_i);
    gkyl_array_release(m1_i);
    gkyl_array_release(m2_i);
  }  
  gkyl_array_release(moms_e);
  gkyl_array_release(moms_i);
  gkyl_array_release(moms_cross);
  gkyl_proj_on_basis_release(proj_m0_e);
  gkyl_proj_on_basis_release(proj_m1_e);
  gkyl_proj_on_basis_release(proj_m2_e);
  gkyl_proj_on_basis_release(proj_m0_i);
  gkyl_proj_on_basis_release(proj_m1_i);
  gkyl_proj_on_basis_release(proj_m2_i);
}

void eval_M2_e_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1];
  double udrift = 0.0, vtsq = 100.0;
  eval_M0_e(t, xn, den, ctx);
  fout[0] = 3*den[0]*vtsq + den[0]*udrift*udrift;
}
void eval_M2_i_2v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den[1];
  double udrift = 1.0, vtsq = 36.0;
  eval_M0_i(t, xn, den, ctx);
  fout[0] = 3*den[0]*vtsq + den[0]*udrift*udrift;
}

void test_1x2v(int poly_order, bool use_gpu)
{
  double me = 1.0; //9.10938215e-31;
  double mi = 2.0; //1.672621637e-27;
  double beta = 0.0;
  double nu_ei = 10.0; 
  double nu_ie = 5.0;

  double lower[] = {-0.5, -5.0, 0.0}, upper[] = {0.5, 5.0, 5.0};
  int cells[] = {2, 32, 32};
  int vdim = 2, cdim = 1;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // Grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // Basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order==1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Configuration space range
  int confGhost[] = {1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost;
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  // Phase space range
  int ghost[] = {confGhost[0], 0, 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost;
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create moment arrays
  // (1) electron
  struct gkyl_array *m0_e_ho, *m1_e_ho, *m2_e_ho;
  m0_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0_e, NULL);
  gkyl_proj_on_basis *proj_m1_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1_e, NULL);
  gkyl_proj_on_basis *proj_m2_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2_e_2v, NULL);
  gkyl_proj_on_basis_advance(proj_m0_e, 0.0, &confLocal, m0_e_ho);
  gkyl_proj_on_basis_advance(proj_m1_e, 0.0, &confLocal, m1_e_ho);
  gkyl_proj_on_basis_advance(proj_m2_e, 0.0, &confLocal, m2_e_ho);
  struct gkyl_array *m0_e, *m1_e, *m2_e;
  if (use_gpu) {
    m0_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_e, m0_e_ho);
    gkyl_array_copy(m1_e, m1_e_ho);
    gkyl_array_copy(m2_e, m2_e_ho);
  } else {
    m0_e = m0_e_ho;
    m1_e = m1_e_ho;
    m2_e = m2_e_ho;
  }
  struct gkyl_array *moms_e = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_e, 1., m0_e, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_e, 1., m1_e, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_e, 1., m2_e, 2*confBasis.num_basis);
  // (2) ion
  struct gkyl_array *m0_i_ho, *m1_i_ho, *m2_i_ho;
  m0_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0_i, NULL);
  gkyl_proj_on_basis *proj_m1_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1_i, NULL);
  gkyl_proj_on_basis *proj_m2_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2_i_2v, NULL);
  gkyl_proj_on_basis_advance(proj_m0_i, 0.0, &confLocal, m0_i_ho);
  gkyl_proj_on_basis_advance(proj_m1_i, 0.0, &confLocal, m1_i_ho);
  gkyl_proj_on_basis_advance(proj_m2_i, 0.0, &confLocal, m2_i_ho);
  struct gkyl_array *m0_i, *m1_i, *m2_i;
  if (use_gpu) {
    m0_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_i, m0_i_ho);
    gkyl_array_copy(m1_i, m1_i_ho);
    gkyl_array_copy(m2_i, m2_i_ho);
  } else {
    m0_i = m0_i_ho;
    m1_i = m1_i_ho;
    m2_i = m2_i_ho;
  }
  struct gkyl_array *moms_i = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_i, 1., m0_i, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_i, 1., m1_i, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_i, 1., m2_i, 2*confBasis.num_basis);
  
  // Calculate the cross moments
  gkyl_mom_cross_bgk_gyrokinetic *momCrossCalc = gkyl_mom_cross_bgk_gyrokinetic_new(&basis, &confBasis);
  struct gkyl_array *moms_cross = mkarr(6*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_mom_cross_bgk_gyrokinetic_advance(momCrossCalc, &confLocal_ext, beta, me, moms_e, mi, moms_i, nu_ei, nu_ie, moms_cross);
  gkyl_mom_cross_bgk_gyrokinetic_release(momCrossCalc);

  // Write out on host
  char fname[1024];
  sprintf(fname, "ctest_mom_cross_bgk_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, moms_cross, fname);

  // Compare with the expected cross moments
  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&confLocal, idx);
    const double *momsCross = gkyl_array_cfetch(moms_cross, linidx);
    TEST_CHECK( gkyl_compare(10.0, momsCross[0*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(10.0, momsCross[1*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(5.0, momsCross[2*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(5.0, momsCross[3*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(2726.6666666666665, momsCross[4*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_CHECK( gkyl_compare(1363.3333333333333, momsCross[5*confBasis.num_basis]/sqrt(2), 1e-12) );
    TEST_MSG("Produced: %.13e, \t%.13e, \t%.13e, %.13e, \t%.13e, \t%.13e", momsCross[0*confBasis.num_basis]/sqrt(2), momsCross[1*confBasis.num_basis]/sqrt(2), momsCross[2*confBasis.num_basis]/sqrt(2), momsCross[3*confBasis.num_basis]/sqrt(2), momsCross[4*confBasis.num_basis]/sqrt(2), momsCross[5*confBasis.num_basis]/sqrt(2));
  }

  // Release memory for moment data object
  gkyl_array_release(m0_e_ho);
  gkyl_array_release(m1_e_ho);
  gkyl_array_release(m2_e_ho);
  gkyl_array_release(m0_i_ho);
  gkyl_array_release(m1_i_ho);
  gkyl_array_release(m2_i_ho);
  if (use_gpu) {
    gkyl_array_release(m0_e);
    gkyl_array_release(m1_e);
    gkyl_array_release(m2_e);
    gkyl_array_release(m0_i);
    gkyl_array_release(m1_i);
    gkyl_array_release(m2_i);
  }  
  gkyl_array_release(moms_e);
  gkyl_array_release(moms_i);
  gkyl_array_release(moms_cross);
  gkyl_proj_on_basis_release(proj_m0_e);
  gkyl_proj_on_basis_release(proj_m1_e);
  gkyl_proj_on_basis_release(proj_m2_e);
  gkyl_proj_on_basis_release(proj_m0_i);
  gkyl_proj_on_basis_release(proj_m1_i);
  gkyl_proj_on_basis_release(proj_m2_i);
}

void test_2x2v(int poly_order, bool use_gpu)
{
  double me = 1.0; //9.10938215e-31;
  double mi = 2.0; //1.672621637e-27;
  double beta = 0.0;
  double nu_ei = 10.0; 
  double nu_ie = 5.0;

  double lower[] = {-0.5, -0.5, -5.0, 0.0}, upper[] = {0.5, 0.5, 5.0, 5.0};
  int cells[] = {2, 2, 32, 32};
  int vdim = 2, cdim = 2;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0], lower[1]}, confUpper[] = {upper[0], upper[1]};
  int confCells[] = {cells[0], cells[1]};

  // Grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // Basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order==1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Configuration space range
  int confGhost[] = {1, 1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost;
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  // Phase space range
  int ghost[] = {confGhost[0], confGhost[1], 0, 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost;
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create moment arrays
  // (1) electron
  struct gkyl_array *m0_e_ho, *m1_e_ho, *m2_e_ho;
  m0_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0_e, NULL);
  gkyl_proj_on_basis *proj_m1_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1_e, NULL);
  gkyl_proj_on_basis *proj_m2_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2_e_2v, NULL);
  gkyl_proj_on_basis_advance(proj_m0_e, 0.0, &confLocal, m0_e_ho);
  gkyl_proj_on_basis_advance(proj_m1_e, 0.0, &confLocal, m1_e_ho);
  gkyl_proj_on_basis_advance(proj_m2_e, 0.0, &confLocal, m2_e_ho);
  struct gkyl_array *m0_e, *m1_e, *m2_e;
  if (use_gpu) {
    m0_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_e, m0_e_ho);
    gkyl_array_copy(m1_e, m1_e_ho);
    gkyl_array_copy(m2_e, m2_e_ho);
  } else {
    m0_e = m0_e_ho;
    m1_e = m1_e_ho;
    m2_e = m2_e_ho;
  }
  struct gkyl_array *moms_e = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_e, 1., m0_e, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_e, 1., m1_e, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_e, 1., m2_e, 2*confBasis.num_basis);
  // (2) ion
  struct gkyl_array *m0_i_ho, *m1_i_ho, *m2_i_ho;
  m0_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0_i, NULL);
  gkyl_proj_on_basis *proj_m1_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1_i, NULL);
  gkyl_proj_on_basis *proj_m2_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2_i_2v, NULL);
  gkyl_proj_on_basis_advance(proj_m0_i, 0.0, &confLocal, m0_i_ho);
  gkyl_proj_on_basis_advance(proj_m1_i, 0.0, &confLocal, m1_i_ho);
  gkyl_proj_on_basis_advance(proj_m2_i, 0.0, &confLocal, m2_i_ho);
  struct gkyl_array *m0_i, *m1_i, *m2_i;
  if (use_gpu) {
    m0_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_i, m0_i_ho);
    gkyl_array_copy(m1_i, m1_i_ho);
    gkyl_array_copy(m2_i, m2_i_ho);
  } else {
    m0_i = m0_i_ho;
    m1_i = m1_i_ho;
    m2_i = m2_i_ho;
  }
  struct gkyl_array *moms_i = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_i, 1., m0_i, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_i, 1., m1_i, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_i, 1., m2_i, 2*confBasis.num_basis);
  
  // Calculate the cross moments
  gkyl_mom_cross_bgk_gyrokinetic *momCrossCalc = gkyl_mom_cross_bgk_gyrokinetic_new(&basis, &confBasis);
  struct gkyl_array *moms_cross = mkarr(6*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_mom_cross_bgk_gyrokinetic_advance(momCrossCalc, &confLocal_ext, beta, me, moms_e, mi, moms_i, nu_ei, nu_ie, moms_cross);
  gkyl_mom_cross_bgk_gyrokinetic_release(momCrossCalc);

  // Write out on host
  char fname[1024];
  sprintf(fname, "ctest_mom_cross_bgk_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, moms_cross, fname);

  // Compare with the expected cross moments
  for (int k=0; k<cells[0]; k++) {
    for (int l=0; l<cells[1]; l++) {
      int idx[] = {k+1, l+1};
      long linidx = gkyl_range_idx(&confLocal, idx);
      const double *momsCross = gkyl_array_cfetch(moms_cross, linidx);
      TEST_CHECK( gkyl_compare(10.0, momsCross[0*confBasis.num_basis]/2, 1e-12) );
      TEST_CHECK( gkyl_compare(10.0, momsCross[1*confBasis.num_basis]/2, 1e-12) );
      TEST_CHECK( gkyl_compare(5.0, momsCross[2*confBasis.num_basis]/2, 1e-12) );
      TEST_CHECK( gkyl_compare(5.0, momsCross[3*confBasis.num_basis]/2, 1e-12) );
      TEST_CHECK( gkyl_compare(2726.6666666666665, momsCross[4*confBasis.num_basis]/2, 1e-12) );
      TEST_CHECK( gkyl_compare(1363.3333333333333, momsCross[5*confBasis.num_basis]/2, 1e-12) );
      TEST_MSG("Produced: %.13e, \t%.13e, \t%.13e, %.13e, \t%.13e, \t%.13e", momsCross[0*confBasis.num_basis]/2, momsCross[1*confBasis.num_basis]/2, momsCross[2*confBasis.num_basis]/2, momsCross[3*confBasis.num_basis]/2, momsCross[4*confBasis.num_basis]/2, momsCross[5*confBasis.num_basis]/2);
    }
  }

  // Release memory for moment data object
  gkyl_array_release(m0_e_ho);
  gkyl_array_release(m1_e_ho);
  gkyl_array_release(m2_e_ho);
  gkyl_array_release(m0_i_ho);
  gkyl_array_release(m1_i_ho);
  gkyl_array_release(m2_i_ho);
  if (use_gpu) {
    gkyl_array_release(m0_e);
    gkyl_array_release(m1_e);
    gkyl_array_release(m2_e);
    gkyl_array_release(m0_i);
    gkyl_array_release(m1_i);
    gkyl_array_release(m2_i);
  }  
  gkyl_array_release(moms_e);
  gkyl_array_release(moms_i);
  gkyl_array_release(moms_cross);
  gkyl_proj_on_basis_release(proj_m0_e);
  gkyl_proj_on_basis_release(proj_m1_e);
  gkyl_proj_on_basis_release(proj_m2_e);
  gkyl_proj_on_basis_release(proj_m0_i);
  gkyl_proj_on_basis_release(proj_m1_i);
  gkyl_proj_on_basis_release(proj_m2_i);
}

void test_3x2v(int poly_order, bool use_gpu)
{
  double me = 1.0; //9.10938215e-31;
  double mi = 2.0; //1.672621637e-27;
  double beta = 0.0;
  double nu_ei = 10.0; 
  double nu_ie = 5.0;

  double lower[] = {-0.5, -0.5, -0.5, -5.0, 0.0}, upper[] = {0.5, 0.5, 0.5, 5.0, 5.0};
  int cells[] = {2, 2, 2, 32, 32};
  int vdim = 2, cdim = 3;
  int ndim = cdim + vdim;

  double confLower[] = {lower[0], lower[1], lower[2]}, confUpper[] = {upper[0], upper[1], upper[2]};
  int confCells[] = {cells[0], cells[1], cells[2]};

  // Grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // Basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order==1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Configuration space range
  int confGhost[] = {1, 1, 1};
  struct gkyl_range confLocal, confLocal_ext;
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost;
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  // Phase space range
  int ghost[] = {confGhost[0], confGhost[1], confGhost[2], 0, 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost;
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  // Create moment arrays
  // (1) electron
  struct gkyl_array *m0_e_ho, *m1_e_ho, *m2_e_ho;
  m0_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_e_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0_e, NULL);
  gkyl_proj_on_basis *proj_m1_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1_e, NULL);
  gkyl_proj_on_basis *proj_m2_e = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2_e_2v, NULL);
  gkyl_proj_on_basis_advance(proj_m0_e, 0.0, &confLocal, m0_e_ho);
  gkyl_proj_on_basis_advance(proj_m1_e, 0.0, &confLocal, m1_e_ho);
  gkyl_proj_on_basis_advance(proj_m2_e, 0.0, &confLocal, m2_e_ho);
  struct gkyl_array *m0_e, *m1_e, *m2_e;
  if (use_gpu) {
    m0_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_e = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_e, m0_e_ho);
    gkyl_array_copy(m1_e, m1_e_ho);
    gkyl_array_copy(m2_e, m2_e_ho);
  } else {
    m0_e = m0_e_ho;
    m1_e = m1_e_ho;
    m2_e = m2_e_ho;
  }
  struct gkyl_array *moms_e = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_e, 1., m0_e, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_e, 1., m1_e, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_e, 1., m2_e, 2*confBasis.num_basis);
  // (2) ion
  struct gkyl_array *m0_i_ho, *m1_i_ho, *m2_i_ho;
  m0_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_i_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0_i, NULL);
  gkyl_proj_on_basis *proj_m1_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1_i, NULL);
  gkyl_proj_on_basis *proj_m2_i = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2_i_2v, NULL);
  gkyl_proj_on_basis_advance(proj_m0_i, 0.0, &confLocal, m0_i_ho);
  gkyl_proj_on_basis_advance(proj_m1_i, 0.0, &confLocal, m1_i_ho);
  gkyl_proj_on_basis_advance(proj_m2_i, 0.0, &confLocal, m2_i_ho);
  struct gkyl_array *m0_i, *m1_i, *m2_i;
  if (use_gpu) {
    m0_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_i = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_i, m0_i_ho);
    gkyl_array_copy(m1_i, m1_i_ho);
    gkyl_array_copy(m2_i, m2_i_ho);
  } else {
    m0_i = m0_i_ho;
    m1_i = m1_i_ho;
    m2_i = m2_i_ho;
  }
  struct gkyl_array *moms_i = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_i, 1., m0_i, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_i, 1., m1_i, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_i, 1., m2_i, 2*confBasis.num_basis);
  
  // Calculate the cross moments
  gkyl_mom_cross_bgk_gyrokinetic *momCrossCalc = gkyl_mom_cross_bgk_gyrokinetic_new(&basis, &confBasis);
  struct gkyl_array *moms_cross = mkarr(6*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_mom_cross_bgk_gyrokinetic_advance(momCrossCalc, &confLocal_ext, beta, me, moms_e, mi, moms_i, nu_ei, nu_ie, moms_cross);
  gkyl_mom_cross_bgk_gyrokinetic_release(momCrossCalc);

  // Write out on host
  char fname[1024];
  sprintf(fname, "ctest_mom_cross_bgk_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, moms_cross, fname);

  // Compare with the expected cross moments
  for (int k=0; k<cells[0]; k++) {
    for (int l=0; l<cells[1]; l++) {
      for (int m=0; m<cells[1]; m++) {
        int idx[] = {k+1, l+1, m+1};
        long linidx = gkyl_range_idx(&confLocal, idx);
        const double *momsCross = gkyl_array_cfetch(moms_cross, linidx);
        TEST_CHECK( gkyl_compare(10.0, momsCross[0*confBasis.num_basis]/sqrt(8), 1e-12) );
        TEST_CHECK( gkyl_compare(10.0, momsCross[1*confBasis.num_basis]/sqrt(8), 1e-12) );
        TEST_CHECK( gkyl_compare(5.0, momsCross[2*confBasis.num_basis]/sqrt(8), 1e-12) );
        TEST_CHECK( gkyl_compare(5.0, momsCross[3*confBasis.num_basis]/sqrt(8), 1e-12) );
        TEST_CHECK( gkyl_compare(2726.6666666666665, momsCross[4*confBasis.num_basis]/sqrt(8), 1e-12) );
        TEST_CHECK( gkyl_compare(1363.3333333333333, momsCross[5*confBasis.num_basis]/sqrt(8), 1e-12) );
        TEST_MSG("Produced: %.13e, \t%.13e, \t%.13e, %.13e, \t%.13e, \t%.13e", momsCross[0*confBasis.num_basis]/sqrt(8), momsCross[1*confBasis.num_basis]/sqrt(8), momsCross[2*confBasis.num_basis]/sqrt(8), momsCross[3*confBasis.num_basis]/sqrt(8), momsCross[4*confBasis.num_basis]/sqrt(8), momsCross[5*confBasis.num_basis]/sqrt(8));
      }
    }
  }

  // Release memory for moment data object
  gkyl_array_release(m0_e_ho);
  gkyl_array_release(m1_e_ho);
  gkyl_array_release(m2_e_ho);
  gkyl_array_release(m0_i_ho);
  gkyl_array_release(m1_i_ho);
  gkyl_array_release(m2_i_ho);
  if (use_gpu) {
    gkyl_array_release(m0_e);
    gkyl_array_release(m1_e);
    gkyl_array_release(m2_e);
    gkyl_array_release(m0_i);
    gkyl_array_release(m1_i);
    gkyl_array_release(m2_i);
  }  
  gkyl_array_release(moms_e);
  gkyl_array_release(moms_i);
  gkyl_array_release(moms_cross);
  gkyl_proj_on_basis_release(proj_m0_e);
  gkyl_proj_on_basis_release(proj_m1_e);
  gkyl_proj_on_basis_release(proj_m2_e);
  gkyl_proj_on_basis_release(proj_m0_i);
  gkyl_proj_on_basis_release(proj_m1_i);
  gkyl_proj_on_basis_release(proj_m2_i);
}
// Run the test
void test_1x1v_p1() {test_1x1v(1, false);}
//void test_1x1v_p2() {test_1x1v(2, true);}
void test_1x2v_p1() {test_1x2v(1, false);}
//void test_1x2v_p2() {test_1x2v(2, true);}
void test_2x2v_p1() {test_2x2v(1, false);}
void test_3x2v_p1() {test_3x2v(1, false);}

TEST_LIST = {
  {"test_1x1v_p1", test_1x1v_p1},
  //{"test_1x1v_p2", test_1x1v_p2},
  {"test_1x2v_p1", test_1x2v_p1},
  //{"test_1x2v_p2", test_1x2v_p2},
  {"test_2x2v_p1", test_2x2v_p1},
  {"test_3x2v_p1", test_3x2v_p1},
  {NULL, NULL},
};
