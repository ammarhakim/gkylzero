#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_correct_maxwellian_gyrokinetic.h>
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

void eval_jacob_tot(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void eval_bmag_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
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

  // Create bmag and jacob_tot arrays
  struct gkyl_array *bmag_ho, *jacob_tot_ho;
  bmag_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  jacob_tot_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_bmag_1x, NULL);
  gkyl_proj_on_basis *proj_jac = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_jacob_tot, NULL);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_proj_on_basis_advance(proj_jac, 0.0, &confLocal, jacob_tot_ho);
  struct gkyl_array *bmag, *jacob_tot;
  if (use_gpu) { 
    // create device copies
    bmag = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    jacob_tot = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    // copy host array to device
    gkyl_array_copy(bmag, bmag_ho);
    gkyl_array_copy(jacob_tot, jacob_tot_ho);
  } else {
    bmag = bmag_ho;
    jacob_tot = jacob_tot_ho;
  }

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
  struct gkyl_array *moms = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms, 1., m0_in, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms, 1., m1_in, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms, 1., m2_in, 2*confBasis.num_basis);
  // (2) create distribution function array
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&grid, &confBasis, &basis, poly_order+1, use_gpu);
  struct gkyl_array *fM = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_maxwellian, &local, &confLocal, moms, bmag, jacob_tot, mass, fM);
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
  gkyl_correct_maxwellian_gyrokinetic *corr_max = gkyl_correct_maxwellian_gyrokinetic_new(&grid, &confBasis, &basis, &confLocal, &confLocal_ext, bmag, mass, poly_order, use_gpu);
  gkyl_correct_maxwellian_gyrokinetic_fix(corr_max, fM, m0_in, m1_in, m2_in, jacob_tot, bmag, mass, err_max, iter_max, &confLocal, &local, &confLocal_ext, use_gpu);
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
  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M0", use_gpu);
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M1", use_gpu);
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2", use_gpu);
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, M0_t, use_gpu);
  gkyl_mom_calc *m1calc = gkyl_mom_calc_new(&grid, M1_t, use_gpu);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, M2_t, use_gpu);
  // (2) create moment arrays
  struct gkyl_array *m0_corr, *m1_corr, *m2_corr;
  m0_corr = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  m1_corr = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  m2_corr = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(m0calc, &local, &confLocal, fM, m0_corr);
    gkyl_mom_calc_advance_cu(m1calc, &local, &confLocal, fM, m1_corr);
    gkyl_mom_calc_advance_cu(m2calc, &local, &confLocal, fM, m2_corr);
  } else {
    gkyl_mom_calc_advance(m0calc, &local, &confLocal, fM, m0_corr);
    gkyl_mom_calc_advance(m1calc, &local, &confLocal, fM, m1_corr);
    gkyl_mom_calc_advance(m2calc, &local, &confLocal, fM, m2_corr);
  }    
  // (3) copy from device to host
  struct gkyl_array *m0_corr_ho, *m1_corr_ho, *m2_corr_ho; 
  m0_corr_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_corr_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_corr_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  if (use_gpu) {
    gkyl_array_copy(m0_corr_ho, m0_corr);
    gkyl_array_copy(m1_corr_ho, m1_corr);
    gkyl_array_copy(m2_corr_ho, m2_corr);
  } else {
    m0_corr_ho = m0_corr;
    m1_corr_ho = m1_corr;
    m2_corr_ho = m2_corr;
  }
  // (4) compare the corrected moments with the input moments 
  for (int k=0; k<cells[0]; k++) {
     int idx[] = {k+1};
     long linidx = gkyl_range_idx(&confLocal, idx);
     const double *m0in = gkyl_array_cfetch(m0_in_ho, linidx);
     const double *m1in = gkyl_array_cfetch(m1_in_ho, linidx);
     const double *m2in = gkyl_array_cfetch(m2_in_ho, linidx);
     const double *m0corr = gkyl_array_cfetch(m0_corr_ho, linidx);
     const double *m1corr = gkyl_array_cfetch(m1_corr_ho, linidx);
     const double *m2corr = gkyl_array_cfetch(m2_corr_ho, linidx);
     for (int m=0; m<basis.num_basis; m++) {
       TEST_CHECK( gkyl_compare(m0in[m], m0corr[m], 1e-12) );
       TEST_CHECK( gkyl_compare(m1in[m], m1corr[m], 1e-12) );
       TEST_CHECK( gkyl_compare(m2in[m], m2corr[m], 1e-12) );
       TEST_MSG("Expected: %.13e, \t%.13e, \t%.13e, \tin cell (%d)", m0in[m], m1in[m], m2in[m], idx[0]);
       TEST_MSG("Produced: %.13e, \t%.13e, \t%.13e", m0corr[m], m1corr[m], m2corr[m]);
     }
   }

  // Release memory for moment data object
  gkyl_array_release(bmag);
  gkyl_array_release(jacob_tot);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(jacob_tot_ho);
  gkyl_array_release(m0_in);
  gkyl_array_release(m1_in);
  gkyl_array_release(m2_in);
  gkyl_array_release(m0_in_ho);
  gkyl_array_release(m1_in_ho);
  gkyl_array_release(m2_in_ho);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(m0_corr_ho);
  gkyl_array_release(m1_corr_ho);
  gkyl_array_release(m2_corr_ho);
  gkyl_array_release(moms);
  gkyl_array_release(fM);
  gkyl_array_release(fM_ho);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_mom_calc_release(m0calc);
  gkyl_mom_calc_release(m1calc);
  gkyl_mom_calc_release(m2calc);
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

  // Create bmag and jacob_tot arrays
  struct gkyl_array *bmag_ho, *jacob_tot_ho;
  bmag_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  jacob_tot_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_bmag_1x, NULL);
  gkyl_proj_on_basis *proj_jac = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_jacob_tot, NULL);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_proj_on_basis_advance(proj_jac, 0.0, &confLocal, jacob_tot_ho);
  struct gkyl_array *bmag, *jacob_tot;
  if (use_gpu) { 
    // create device copies
    bmag = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    jacob_tot = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
    // copy host array to device
    gkyl_array_copy(bmag, bmag_ho);
    gkyl_array_copy(jacob_tot, jacob_tot_ho);
  } else {
    bmag = bmag_ho;
    jacob_tot = jacob_tot_ho;
  }

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
  struct gkyl_array *moms = mkarr(3*confBasis.num_basis, confLocal_ext.volume, use_gpu);
  gkyl_array_set_offset(moms, 1., m0_in, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms, 1., m1_in, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms, 1., m2_in, 2*confBasis.num_basis);
  // (2) create distribution function array
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&grid, &confBasis, &basis, poly_order+1, use_gpu);
  struct gkyl_array *fM = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_proj_gkmaxwellian_on_basis_lab_mom(proj_maxwellian, &local, &confLocal, moms, bmag, jacob_tot, mass, fM);
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
  gkyl_correct_maxwellian_gyrokinetic *corr_max = gkyl_correct_maxwellian_gyrokinetic_new(&grid, &confBasis, &basis, &confLocal, &confLocal_ext, bmag, mass, poly_order, use_gpu);
  gkyl_correct_maxwellian_gyrokinetic_fix(corr_max, fM, m0_in, m1_in, m2_in, jacob_tot, bmag, mass, err_max, iter_max, &confLocal, &local, &confLocal_ext, use_gpu);
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
  struct gkyl_mom_type *M0_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M0", use_gpu);
  struct gkyl_mom_type *M1_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M1", use_gpu);
  struct gkyl_mom_type *M2_t = gkyl_mom_gyrokinetic_new(&confBasis, &basis, &confLocal, mass, "M2", use_gpu);
  gkyl_gyrokinetic_set_bmag(M0_t, bmag);
  gkyl_gyrokinetic_set_bmag(M1_t, bmag);
  gkyl_gyrokinetic_set_bmag(M2_t, bmag);
  gkyl_mom_calc *m0calc = gkyl_mom_calc_new(&grid, M0_t, use_gpu);
  gkyl_mom_calc *m1calc = gkyl_mom_calc_new(&grid, M1_t, use_gpu);
  gkyl_mom_calc *m2calc = gkyl_mom_calc_new(&grid, M2_t, use_gpu);
  // (2) create moment arrays
  struct gkyl_array *m0_corr, *m1_corr, *m2_corr;
  m0_corr = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  m1_corr = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  m2_corr = mkarr(confBasis.num_basis, confLocal_ext.volume, use_gpu);
  if (use_gpu) {
    gkyl_mom_calc_advance_cu(m0calc, &local, &confLocal, fM, m0_corr);
    gkyl_mom_calc_advance_cu(m1calc, &local, &confLocal, fM, m1_corr);
    gkyl_mom_calc_advance_cu(m2calc, &local, &confLocal, fM, m2_corr);
  } else {
    gkyl_mom_calc_advance(m0calc, &local, &confLocal, fM, m0_corr);
    gkyl_mom_calc_advance(m1calc, &local, &confLocal, fM, m1_corr);
    gkyl_mom_calc_advance(m2calc, &local, &confLocal, fM, m2_corr);
  }    
  // (3) copy from device to host
  struct gkyl_array *m0_corr_ho, *m1_corr_ho, *m2_corr_ho; 
  m0_corr_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m1_corr_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  m2_corr_ho = mkarr(confBasis.num_basis, confLocal_ext.volume, false);
  if (use_gpu) {
    gkyl_array_copy(m0_corr_ho, m0_corr);
    gkyl_array_copy(m1_corr_ho, m1_corr);
    gkyl_array_copy(m2_corr_ho, m2_corr);
  } else {
    m0_corr_ho = m0_corr;
    m1_corr_ho = m1_corr;
    m2_corr_ho = m2_corr;
  }
  // (4) compare the corrected moments with the input moments 
  for (int k=0; k<cells[0]; k++) {
     int idx[] = {k+1};
     long linidx = gkyl_range_idx(&confLocal, idx);
     const double *m0in = gkyl_array_cfetch(m0_in_ho, linidx);
     const double *m1in = gkyl_array_cfetch(m1_in_ho, linidx);
     const double *m2in = gkyl_array_cfetch(m2_in_ho, linidx);
     const double *m0corr = gkyl_array_cfetch(m0_corr_ho, linidx);
     const double *m1corr = gkyl_array_cfetch(m1_corr_ho, linidx);
     const double *m2corr = gkyl_array_cfetch(m2_corr_ho, linidx);
     for (int m=0; m<basis.num_basis; m++) {
       TEST_CHECK( gkyl_compare(m0in[m], m0corr[m], 1e-12) );
       TEST_CHECK( gkyl_compare(m1in[m], m1corr[m], 1e-12) );
       TEST_CHECK( gkyl_compare(m2in[m], m2corr[m], 1e-12) );
       TEST_MSG("Expected: %.13e, \t%.13e, \t%.13e, \tin cell (%d)", m0in[m], m1in[m], m2in[m], idx[0]);
       TEST_MSG("Produced: %.13e, \t%.13e, \t%.13e", m0corr[m], m1corr[m], m2corr[m]);
     }
   }

  // Release memory for moment data object
  gkyl_array_release(bmag);
  gkyl_array_release(jacob_tot);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(jacob_tot_ho);
  gkyl_array_release(m0_in);
  gkyl_array_release(m1_in);
  gkyl_array_release(m2_in);
  gkyl_array_release(m0_in_ho);
  gkyl_array_release(m1_in_ho);
  gkyl_array_release(m2_in_ho);
  gkyl_array_release(m0_corr);
  gkyl_array_release(m1_corr);
  gkyl_array_release(m2_corr);
  gkyl_array_release(m0_corr_ho);
  gkyl_array_release(m1_corr_ho);
  gkyl_array_release(m2_corr_ho);
  gkyl_array_release(moms);
  gkyl_array_release(fM);
  gkyl_array_release(fM_ho);
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_mom_calc_release(m0calc);
  gkyl_mom_calc_release(m1calc);
  gkyl_mom_calc_release(m2calc);
}

// Run the test
void test_1x1v_p1() {test_1x1v(1, true);}
void test_1x2v_p1() {test_1x2v(1, true);}

TEST_LIST = {
  {"test_1x1v_p1", test_1x1v_p1},
  {"test_1x2v_p1", test_1x2v_p1},
  {NULL, NULL},
};
