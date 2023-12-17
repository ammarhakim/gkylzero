#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_correct_maxwellian_gyrokinetic.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <math.h>

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

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; 
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0];
  fout[0] = cos((2.*M_PI/(2.*2.*M_PI))*x);
}

void eval_M0(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}
void eval_M1(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.5;
}
void eval_M2(double t, const double *xn, double* restrict fout, void *ctx)
{ 
  double n = 1.0, vtsq = 1.0, ux = 0.5;
  double x = xn[0];
  fout[0] = n*vtsq + n*ux*ux;
}
void eval_vtsq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}

void test_1x1v(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double err_max = 1e-14, iter_max = 50;
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
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
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

  // Initialize geometry
  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&confGrid, &confLocal, &confLocal_ext, &confBasis, 
    mapc2p, 0, bmag_func, 0, use_gpu);

  // Create correct moment arrays
  struct gkyl_array *m0_in_ho, *m1_in_ho, *m2_in_ho;
  m0_in_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_in_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_in_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2, NULL);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_in_ho);
  gkyl_proj_on_basis_advance(proj_m1, 0.0, &confLocal, m1_in_ho);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_in_ho);
  struct gkyl_array *m0_in, *m1_in, *m2_in;
  if (use_gpu) {
    m0_in = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_in = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_in = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_in, m0_in_ho);
    gkyl_array_copy(m1_in, m1_in_ho);
    gkyl_array_copy(m2_in, m2_in_ho);
  } else {
    m0_in = m0_in_ho;
    m1_in = m1_in_ho;
    m2_in = m2_in_ho;
  }

  // Project the Maxwellian on basis
  // (1) proj_maxwellian expects the moments as a single array
  struct gkyl_array *moms_in = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_in, 1., m0_in, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m1_in, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m2_in, 2*confBasis.num_basis);
  // (2) create distribution function array
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&grid, &confBasis, &basis, poly_order+1, use_gpu);
  struct gkyl_array *fM = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_maxwellian, &local, &confLocal, moms_in, gk_geom->bmag, gk_geom->jacobtot, mass, fM);
  // (3) copy from device to host
  struct gkyl_array *fM_ho;
  if (use_gpu) {
    fM_ho = mkarr(basis.num_basis, local_ext.volume, false);
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }
  // (4) write out on host
  char fname_fM_ic[1024];
  sprintf(fname_fM_ic, "ctest_correct_maxwellian_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fM_ho, fname_fM_ic);
 
  // Create a Maxwellian with corrected moments
  gkyl_correct_maxwellian_gyrokinetic *corr_max = gkyl_correct_maxwellian_gyrokinetic_new(&grid, &confBasis, &basis, &confLocal, &confLocal_ext, 
    mass, gk_geom, use_gpu);
  gkyl_correct_maxwellian_gyrokinetic_fix(corr_max, fM, moms_in, err_max, iter_max, &confLocal, &local, &confLocal_ext);
  gkyl_correct_maxwellian_gyrokinetic_release(corr_max);
  gkyl_array_clear(fM_ho, 0.0);
  if (use_gpu) {
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }
  char fname_fM_corr[1024];
  sprintf(fname_fM_corr, "ctest_correct_maxwellian_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fM_ho, fname_fM_corr);

  // Calculate the corrected moments
  // (1) create the calculators
  struct gkyl_mom_type *MOMS_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, gk_geom, "ThreeMoments", use_gpu);
  gkyl_mom_calc *momsCalc = gkyl_mom_calc_new(&grid, MOMS_t, use_gpu);
  // (2) calculate the moments and copy from host to device
  struct gkyl_array *moms_corr = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  struct gkyl_array *moms_corr_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume, false);
  if (use_gpu) { 
    gkyl_mom_calc_advance_cu(momsCalc, &local, &confLocal, fM, moms_corr); 
    gkyl_array_copy(moms_corr_ho, moms_corr); 
  } else {
    gkyl_mom_calc_advance(momsCalc, &local, &confLocal, fM, moms_corr);
    moms_corr_ho = moms_corr;
  }
  // (4) compare the correct moments with the input moments
  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&confLocal, idx);
    const double *m0in = gkyl_array_cfetch(m0_in_ho, linidx);
    const double *m1in = gkyl_array_cfetch(m1_in_ho, linidx);
    const double *m2in = gkyl_array_cfetch(m2_in_ho, linidx);
    const double *momsCorr = gkyl_array_cfetch(moms_corr_ho, linidx);
    for (int m=0; m<confBasis.num_basis; m++) {
      TEST_CHECK( gkyl_compare(m0in[m], momsCorr[m+0*confBasis.num_basis], 1e-12) );
      TEST_CHECK( gkyl_compare(m1in[m], momsCorr[m+1*confBasis.num_basis], 1e-12) );
      TEST_CHECK( gkyl_compare(m2in[m], momsCorr[m+2*confBasis.num_basis], 1e-12) );
      TEST_MSG("Expected: %.13e, \t%.13e, \t%.13e, \tin cell (%d)", m0in[m], m1in[m], m2in[m], idx[0]);
      TEST_MSG("Produced: %.13e, \t%.13e, \t%.13e", momsCorr[m+0*confBasis.num_basis], momsCorr[m+1*confBasis.num_basis], momsCorr[m+2*confBasis.num_basis]);
    }
  }
  
  // Release memory for moment data object
  gkyl_gk_geometry_release(gk_geom);  
  gkyl_array_release(m0_in);
  gkyl_array_release(m1_in);
  gkyl_array_release(m2_in);
  gkyl_array_release(m0_in_ho);
  gkyl_array_release(m1_in_ho);
  gkyl_array_release(m2_in_ho);
  gkyl_array_release(moms_in);
  gkyl_array_release(moms_corr);
  gkyl_array_release(moms_corr_ho);
  gkyl_array_release(fM);
  gkyl_array_release(fM_ho);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_mom_calc_release(momsCalc);
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

void test_1x2v(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double err_max = 1e-14, iter_max = 50;
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

  // Initialize geometry
  struct gk_geometry *gk_geom = gkyl_gk_geometry_mapc2p_new(&confGrid, &confLocal, &confLocal_ext, &confBasis, 
    mapc2p, 0, bmag_func, 0, use_gpu);

  // Create correct moment arrays
  struct gkyl_array *m0_in_ho, *m1_in_ho, *m2_in_ho;
  m0_in_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_in_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_in_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_m1 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M1_2v_gk, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_M2_2v_gk, NULL);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_in_ho);
  gkyl_proj_on_basis_advance(proj_m1, 0.0, &confLocal, m1_in_ho);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_in_ho);
  struct gkyl_array *m0_in, *m1_in, *m2_in;
  if (use_gpu) {
    m0_in = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m1_in = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    m2_in = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    gkyl_array_copy(m0_in, m0_in_ho);
    gkyl_array_copy(m1_in, m1_in_ho);
    gkyl_array_copy(m2_in, m2_in_ho);
  } else {
    m0_in = m0_in_ho;
    m1_in = m1_in_ho;
    m2_in = m2_in_ho;
  }

  // Project the Maxwellian on basis
  // (1) proj_maxwellian expects the moments as a single array
  struct gkyl_array *moms_in = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms_in, 1., m0_in, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m1_in, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m2_in, 2*confBasis.num_basis);
  // (2) create distribution function array
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&grid, &confBasis, &basis, poly_order+1, use_gpu);
  struct gkyl_array *fM = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_maxwellian, &local, &confLocal, moms_in, gk_geom->bmag, gk_geom->jacobtot, mass, fM);
  // (3) copy from device to host
  struct gkyl_array *fM_ho;
  if (use_gpu) {
    fM_ho = mkarr(basis.num_basis, local_ext.volume, false);
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }
  // (4) write out on host
  char fname_fM_ic[1024];
  sprintf(fname_fM_ic, "ctest_correct_maxwellian_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fM_ho, fname_fM_ic);
 
  // Create a Maxwellian with corrected moments
  gkyl_correct_maxwellian_gyrokinetic *corr_max = gkyl_correct_maxwellian_gyrokinetic_new(&grid, &confBasis, &basis, &confLocal, &confLocal_ext, 
    mass, gk_geom, use_gpu);
  gkyl_correct_maxwellian_gyrokinetic_fix(corr_max, fM, moms_in, err_max, iter_max, &confLocal, &local, &confLocal_ext);
  gkyl_correct_maxwellian_gyrokinetic_release(corr_max);
  gkyl_array_clear(fM_ho, 0.0);
  if (use_gpu) {
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }
  char fname_fM_corr[1024];
  sprintf(fname_fM_corr, "ctest_correct_maxwellian_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fM_ho, fname_fM_corr);

  // Calculate the corrected moments
  // (1) create the calculators
  struct gkyl_mom_type *MOMS_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, gk_geom, "ThreeMoments", use_gpu);
  gkyl_mom_calc *momsCalc = gkyl_mom_calc_new(&grid, MOMS_t, use_gpu);
  // (2) calculate the moments and copy from host to device
  struct gkyl_array *moms_corr = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  struct gkyl_array *moms_corr_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume, false);
  if (use_gpu) { 
    gkyl_mom_calc_advance_cu(momsCalc, &local, &confLocal, fM, moms_corr); 
    gkyl_array_copy(moms_corr_ho, moms_corr); 
  } else {
    gkyl_mom_calc_advance(momsCalc, &local, &confLocal, fM, moms_corr);
    moms_corr_ho = moms_corr;
  }
  // (4) compare the correct moments with the input moments
  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&confLocal, idx);
    const double *m0in = gkyl_array_cfetch(m0_in_ho, linidx);
    const double *m1in = gkyl_array_cfetch(m1_in_ho, linidx);
    const double *m2in = gkyl_array_cfetch(m2_in_ho, linidx);
    const double *momsCorr = gkyl_array_cfetch(moms_corr_ho, linidx);
    for (int m=0; m<confBasis.num_basis; m++) {
      TEST_CHECK( gkyl_compare(m0in[m], momsCorr[m+0*confBasis.num_basis], 1e-12) );
      TEST_CHECK( gkyl_compare(m1in[m], momsCorr[m+1*confBasis.num_basis], 1e-12) );
      TEST_CHECK( gkyl_compare(m2in[m], momsCorr[m+2*confBasis.num_basis], 1e-12) );
      TEST_MSG("Expected: %.13e, \t%.13e, \t%.13e, \tin cell (%d)", m0in[m], m1in[m], m2in[m], idx[0]);
      TEST_MSG("Produced: %.13e, \t%.13e, \t%.13e", momsCorr[m+0*confBasis.num_basis], momsCorr[m+1*confBasis.num_basis], momsCorr[m+2*confBasis.num_basis]);
    }
  }

  // Release memory for moment data object
  gkyl_gk_geometry_release(gk_geom);  
  gkyl_array_release(m0_in);
  gkyl_array_release(m1_in);
  gkyl_array_release(m2_in);
  gkyl_array_release(m0_in_ho);
  gkyl_array_release(m1_in_ho);
  gkyl_array_release(m2_in_ho);
  gkyl_array_release(moms_in);
  gkyl_array_release(moms_corr);
  gkyl_array_release(moms_corr_ho);
  gkyl_array_release(fM);
  gkyl_array_release(fM_ho);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_mom_calc_release(momsCalc);
}

// Run the test
void test_1x1v_p1() {test_1x1v(1, false);}
void test_1x1v_p2() {test_1x1v(2, false);}
void test_1x2v_p1() {test_1x2v(1, false);}
void test_1x2v_p2() {test_1x2v(2, false);}

#ifdef GKYL_HAVE_CUDA
void test_1x1v_p1_gpu() {test_1x1v(1, true);}
void test_1x1v_p2_gpu() {test_1x1v(2, true);}
void test_1x2v_p1_gpu() {test_1x2v(1, true);}
void test_1x2v_p2_gpu() {test_1x2v(2, true);}
#endif

TEST_LIST = {
  {"test_1x1v_p1", test_1x1v_p1},
  {"test_1x1v_p2", test_1x1v_p2},
  {"test_1x2v_p1", test_1x2v_p1},
  {"test_1x2v_p2", test_1x2v_p2},
#ifdef GKYL_HAVE_CUDA
  {"test_1x1v_p1_gpu", test_1x1v_p1_gpu},
  {"test_1x1v_p2_gpu", test_1x1v_p2_gpu},
  {"test_1x2v_p1_gpu", test_1x2v_p1_gpu},
  {"test_1x2v_p2_gpu", test_1x2v_p2_gpu},
#endif
  {NULL, NULL},
};
