// Tests for the Vlasov FPO routines to compute moments and boundary corrections
// and resulting corrections to drag and diffusion coefficients.
//
// Included are tests for p=1 hybrid and p=2 serendipity.
// As shown in Rodman 2025 PhD Thesis, the p=2 computation enforces
// conservation of momentum and energy, so this result is taken as the
// "correct" result for the hybrid basis.

#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>

#include <gkyl_dg_fpo_vlasov_diff_coeff.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff.h>
#include <gkyl_fpo_vlasov_coeff_correct.h>
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_fpo_vlasov.h>
#include <gkyl_dg_updater_moment.h>

struct fpo_ctx {
  double n0;
  double ux0;
  double uy0;
  double uz0;
  double vth0;
  double gamma0;
};

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

// allocate cu_dev array
static struct gkyl_array*
mkarr_cu(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  return a;
}

static inline double
bump_maxwellian(double n, double vx, double vy, double vz, double ux, double uy, double uz, double vt, double bA, double bUx, double bUy, double bUz, double bS, double bVt)
{
  double v2 = (vx - ux)*(vx - ux) + (vy - uy)*(vy - uy) + (vz - uz)*(vz - uz);
  double bv2 = (vx - bUx)*(vx - bUx) + (vy - bUy)*(vy - bUy) + (vz - bUz)*(vz - bUz);
  return n/pow(sqrt(2*M_PI*vt*vt), 3)*exp(-v2/(2*vt*vt)) + n/pow(sqrt(2*M_PI*bVt*bVt), 3)*exp(-bv2/(2*bVt*bVt))*(bA*bA)/(bv2 + bS*bS);
}

void
eval_distf_square(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  // double v2 = vx*vx+vy*vy+vz*vz;
  // fout[0] = 2.5/pow(sqrt(2*M_PI*1.0*1.0), 3)*exp(-v2/(2.0*1.0*1.0));
  double width = 2.0; // corresponds to a final vth of 2/3
  if(vx>-width && vx<width && vy>-width && vy<width && vz>-width && vz<width) {
    fout[0] = 0.5;
  } else {
    fout[0] = 0.0;
  }
}

void
eval_distf_bump(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = bump_maxwellian(1.0, vx, vy, vz, 0.0, 0.0, 0.0, app->vth0, 
    sqrt(0.15), 4.0*app->vth0, 0.0, 0.0, 0.14, 3.0*app->vth0);
}

void
eval_drag_coeff(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = vy*vz;
  fout[1] = vx*vz;
  fout[2] = vy*vz;
}

void
eval_diff_coeff(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  fout[0] = vx*vx;
  fout[1] = vx*vy;
  fout[2] = vx*vz;
  fout[3] = vy*vx;
  fout[4] = vy*vy;
  fout[5] = vz*vz;
  fout[6] = vz*vx;
  fout[7] = vz*vy;
  fout[8] = vz*vz;
}

struct fpo_ctx
create_ctx() {
  struct fpo_ctx ctx = {
    .n0 = 1.0,
    .ux0 = 0.0,
    .uy0 = 0.0,
    .uz0 = 0.0,
    .vth0 = 1.0,
    .gamma0 = 1.0,
  };
return ctx;
}

void test_1x3v_square_p1(int NV)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {2, NV, NV, NV};
  int cells_vel[] = {NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 4.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};
  double lower_vel[] = {-L, -L, -L};
  double upper_vel[] = {L, L, L};

  // Configuration space grid
  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  // Phase space grid
  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // Velocity space grid
  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range, vel_range_ext;
  gkyl_rect_grid_init(&vel_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vel_grid, ghost, &vel_range_ext, &vel_range);

  // initialize basis
  int poly_order = 1;
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_hybrid(&phase_basis, cdim, vdim);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_hybrid(&surf_basis, cdim, vdim-1);

  int num_quad = poly_order+3;
  gkyl_proj_on_basis *proj_distf_square = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, 1, eval_distf_square, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim*vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume);
  moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume);
  fpo_moms = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume);
  boundary_corrections = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume);
  drag_diff_coeff_corrs = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume);

  // Initialize updater to compute Five Moments
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_new(&phase_grid,
    &conf_basis, &phase_basis, &conf_range, &vel_range, &phase_range, 0, 0,
    GKYL_F_MOMENT_M0M1M2, 0, 0);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&conf_basis,
    &phase_basis, &phase_range, 0);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {
    .a = drag_coeff, .D = diff_coeff };
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, 0);

  // Initialize updater to compute boundary corrections
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = -L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&phase_grid,
    &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, 0);

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc = gkyl_fpo_coeff_correct_new(&phase_grid,
    &conf_basis, &conf_range, 0);

  // Project distribution function and coefficients
  gkyl_proj_on_basis_advance(proj_distf_square, 0.0, &phase_range, distf);
  gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
  gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  gkyl_array_scale(drag_coeff, 10.0);

  // Compute Five Moments
  gkyl_dg_updater_moment_advance(mcalc, &phase_range, &conf_range, distf, moms);

  // Compute moments, boundary corrections, and coefficient corrections
  gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, 
    distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(coeff_correct_calc,
    &conf_range, &phase_range, fpo_moms, boundary_corrections,
    moms, drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf,
    diff_coeff, diff_coeff_surf, 0);

  const double *fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
  const double *bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
  const double *coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);

  if (NV == 4) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -5.8049951547205602e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], -7.1054273576010019e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -1.2414965746617999e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], -8.4706424896707739e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 8.4475690125752760e+02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -1.4356862112509038e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 3.7613110601685491e-32, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.8140609858501767e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], 2.2204460492503151e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 3.8796767958181285e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 4.9742589423178196e-34, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], 2.6470757780221199e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -8.7995510547659190e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -3.7443003455994221e-17, 1e-12) );
  }
  if (NV == 8) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], -6.0396132539608516e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], 3.1245918384568393e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], -3.5527136788005009e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -1.6293560982028782e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -1.7763568394002505e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], 1.1098291389402597e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 8.4475690125752794e+02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -5.3596093236574169e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 1.8873791418627693e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], -9.7643494951776407e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], 1.1102230246251585e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 5.0917378068840032e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 5.5511151231257926e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], -3.4682160591883174e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -8.7995510547659315e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 3.7129898741968476e-16, 1e-12) );
  }
}

void test_1x3v_bump_p1(int NV)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {2, NV, NV, NV};
  int cells_vel[] = {NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 4.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};
  double lower_vel[] = {-L, -L, -L};
  double upper_vel[] = {L, L, L};

  // Configuration space grid
  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  // Phase space grid
  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // Velocity space grid
  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range, vel_range_ext;
  gkyl_rect_grid_init(&vel_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vel_grid, ghost, &vel_range_ext, &vel_range);

  // initialize basis
  int poly_order = 1;
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_hybrid(&phase_basis, cdim, vdim);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_hybrid(&surf_basis, cdim, vdim-1);

  int num_quad = poly_order+3;
  gkyl_proj_on_basis *proj_distf_bump = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, 1, eval_distf_bump, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim*vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume);
  moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume);
  fpo_moms = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume);
  boundary_corrections = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume);
  drag_diff_coeff_corrs = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume);

  // Initialize updater to compute Five Moments
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_new(&phase_grid,
    &conf_basis, &phase_basis, &conf_range, &vel_range, &phase_range, 0, 0,
    GKYL_F_MOMENT_M0M1M2, 0, 0);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&conf_basis,
    &phase_basis, &phase_range, 0);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {
    .a = drag_coeff, .D = diff_coeff };
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, 0);

  // Initialize updater to compute boundary corrections
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = -L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&phase_grid,
    &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, 0);

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc = gkyl_fpo_coeff_correct_new(&phase_grid,
    &conf_basis, &conf_range, 0);

  // Project distribution function and coefficients
  gkyl_proj_on_basis_advance(proj_distf_bump, 0.0, &phase_range, distf);
  gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
  gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  gkyl_array_scale(drag_coeff, 10.0);

  // Compute Five Moments
  gkyl_dg_updater_moment_advance(mcalc, &phase_range, &conf_range, distf, moms);

  // Compute moments, boundary corrections, and coefficient corrections
  gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, 
    distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(coeff_correct_calc,
    &conf_range, &phase_range, fpo_moms, boundary_corrections,
    moms, drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf,
    diff_coeff, diff_coeff_surf, 0);

  const double *fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
  const double *bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
  const double *coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);

  if (NV == 4) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.1685375343585043e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -4.1710971760487405e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], -2.5777990853015353e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -2.4714311593325799e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -4.0939474033052647e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], 2.3189744101340079e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 2.0443723328499349e+01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 2.2091690949038620e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.7187609281138045e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], -4.9965538495170426e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], -5.2922618544448685e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 5.8209395583300142e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 1.0096632731271260e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], 1.3085323203632530e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], -3.8875043712455232e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], -8.5316721167048962e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 1.6228429243719039e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 1.2205121280494322e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -1.2348316747910406e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.8586360148702548e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], 2.5544463610232715e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 3.2846815926326794e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 3.9897666484364568e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], -2.2420438408924226e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -6.7917415582715881e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -8.3754623025259556e-17, 1e-12) );
  }
  if (NV == 8) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.1500173245089075e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -1.3814606337380109e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[2], -2.4654757402320371e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -2.6887040032578782e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -4.6403852982379590e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[5], 1.9031362467424888e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 2.0436859018474077e+01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 9.3620119931091115e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.9720462549479656e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 2.9401345535727655e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[2], -3.7608262858090935e-19, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 1.0130431394984520e-19, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], -3.0408482806429382e-19, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[5], -5.0549075583401925e-20, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], -3.9888185019791869e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], -1.1963559167756827e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[8], 1.6589403275025838e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 1.0962162989424440e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 2.4172048297609167e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 2.2824915248939008e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[2], 2.2831057173690739e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 2.7319520996301737e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 4.5082533140980421e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[5], -1.8857979951717035e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], -6.8722195535544301e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -5.1851480741430756e-16, 1e-12) );
  }
}

void test_1x3v_square_p2(int NV)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {2, NV, NV, NV};
  int cells_vel[] = {NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 4.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};
  double lower_vel[] = {-L, -L, -L};
  double upper_vel[] = {L, L, L};

  // Configuration space grid
  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  // Phase space grid
  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // Velocity space grid
  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range, vel_range_ext;
  gkyl_rect_grid_init(&vel_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vel_grid, ghost, &vel_range_ext, &vel_range);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);

  int num_quad = poly_order+2;
  gkyl_proj_on_basis *proj_distf_square = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, 1, eval_distf_square, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim*vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume);
  moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume);
  fpo_moms = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume);
  boundary_corrections = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume);
  drag_diff_coeff_corrs = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume);

  // Initialize updater to compute Five Moments
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_new(&phase_grid,
    &conf_basis, &phase_basis, &conf_range, &vel_range, &phase_range, 0, 0,
    GKYL_F_MOMENT_M0M1M2, 0, 0);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&conf_basis,
    &phase_basis, &phase_range, 0);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {
    .a = drag_coeff, .D = diff_coeff };
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, 0);

  // Initialize updater to compute boundary corrections
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = -L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&phase_grid,
    &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, 0);

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc = gkyl_fpo_coeff_correct_new(&phase_grid,
    &conf_basis, &conf_range, 0);

  // Project distribution function and coefficients
  gkyl_proj_on_basis_advance(proj_distf_square, 0.0, &phase_range, distf);
  gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
  gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  gkyl_array_scale(drag_coeff, 10.0);

  // Compute Five Moments
  gkyl_dg_updater_moment_advance(mcalc, &phase_range, &conf_range, distf, moms);

  // Compute moments, boundary corrections, and coefficient corrections
  gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, 
    distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(coeff_correct_calc,
    &conf_range, &phase_range, fpo_moms, boundary_corrections,
    moms, drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf,
    diff_coeff, diff_coeff_surf, 0);

  const double *fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
  const double *bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
  const double *coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);

  if (NV == 4) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -5.8049951547205602e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -7.1054273576010019e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -1.2414965746617999e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], -8.4706424896707739e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[9], 8.4475690125752760e+02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[10], -1.4356862112509038e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[16], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[18], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[19], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[21], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[22], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 3.7613110601685491e-32, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.8140609858501767e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 2.2204460492503151e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 3.8796767958181285e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 4.9742589423178196e-34, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], 2.6470757780221199e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[9], -8.7995510547659190e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[10], -3.7443003455994221e-17, 1e-12) );
  }
  if (NV == 8) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], -6.0396132539608516e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], 3.1245918384568393e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -3.5527136788005009e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -1.6293560982028782e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -1.7763568394002505e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 1.1098291389402597e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[9], 8.4475690125752794e+02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[10], -5.3596093236574169e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[16], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[18], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[19], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[21], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[22], 0.0000000000000000e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 1.8873791418627693e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], -9.7643494951776407e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 1.1102230246251585e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 5.0917378068840032e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 5.5511151231257926e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -3.4682160591883174e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[9], -8.7995510547659315e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[10], 3.7129898741968476e-16, 1e-12) );
  }
}

void test_1x3v_bump_p2(int NV)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {2, NV, NV, NV};
  int cells_vel[] = {NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 4.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};
  double lower_vel[] = {-L, -L, -L};
  double upper_vel[] = {L, L, L};

  // Configuration space grid
  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  // Phase space grid
  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // Velocity space grid
  struct gkyl_rect_grid vel_grid;
  struct gkyl_range vel_range, vel_range_ext;
  gkyl_rect_grid_init(&vel_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&vel_grid, ghost, &vel_range_ext, &vel_range);

  // initialize basis
  int poly_order = 2;
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
  gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
  gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);

  int num_quad = poly_order+2;
  gkyl_proj_on_basis *proj_distf_bump = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, 1, eval_distf_bump, &ctx);
  gkyl_proj_on_basis *proj_drag_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim, eval_drag_coeff, &ctx);
  gkyl_proj_on_basis *proj_diff_coeff = gkyl_proj_on_basis_new(&phase_grid, &phase_basis,
    num_quad, vdim*vdim, eval_diff_coeff, &ctx);

  struct gkyl_array *distf, *moms;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  struct gkyl_array *fpo_moms, *boundary_corrections, *drag_diff_coeff_corrs;
  distf = mkarr(phase_basis.num_basis, phase_range_ext.volume);
  moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume);
  fpo_moms = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume);
  boundary_corrections = mkarr(2*(vdim+1)*conf_basis.num_basis, conf_range_ext.volume);
  drag_diff_coeff_corrs = mkarr((vdim+1)*conf_basis.num_basis, conf_range_ext.volume);

  // Initialize updater to compute Five Moments
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_new(&phase_grid,
    &conf_basis, &phase_basis, &conf_range, &vel_range, &phase_range, 0, 0,
    GKYL_F_MOMENT_M0M1M2, 0, 0);

  // Initialize updater to compute correction moments
  const struct gkyl_mom_type* fpo_mom_type = gkyl_mom_fpo_vlasov_new(&conf_basis,
    &phase_basis, &phase_range, 0);
  struct gkyl_mom_fpo_vlasov_auxfields fpo_mom_auxfields = {
    .a = drag_coeff, .D = diff_coeff };
  gkyl_mom_fpo_vlasov_set_auxfields(fpo_mom_type, fpo_mom_auxfields);
  struct gkyl_mom_calc *fpo_mom_calc = gkyl_mom_calc_new(&phase_grid, fpo_mom_type, 0);

  // Initialize updater to compute boundary corrections
  double v_bounds[2*GKYL_MAX_DIM];
  for (int d=0; d<vdim; ++d) {
    v_bounds[d] = -L;
    v_bounds[d + vdim] = -L;
  }
  struct gkyl_mom_calc_bcorr *bcorr_calc = gkyl_mom_calc_bcorr_fpo_vlasov_new(&phase_grid,
    &conf_basis, &phase_basis, &phase_range, v_bounds, diff_coeff, 0);

  // Initialize updater to compute corrections to coefficients
  gkyl_fpo_coeff_correct *coeff_correct_calc = gkyl_fpo_coeff_correct_new(&phase_grid,
    &conf_basis, &conf_range, 0);

  // Project distribution function and coefficients
  gkyl_proj_on_basis_advance(proj_distf_bump, 0.0, &phase_range, distf);
  gkyl_proj_on_basis_advance(proj_drag_coeff, 0.0, &phase_range, drag_coeff);
  gkyl_proj_on_basis_advance(proj_diff_coeff, 0.0, &phase_range, diff_coeff);
  gkyl_array_scale(drag_coeff, 10.0);

  // Compute Five Moments
  gkyl_dg_updater_moment_advance(mcalc, &phase_range, &conf_range, distf, moms);

  // Compute moments, boundary corrections, and coefficient corrections
  gkyl_mom_calc_advance(fpo_mom_calc, &phase_range, &conf_range, distf, fpo_moms);
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &phase_range, &conf_range, 
    distf, boundary_corrections);
  gkyl_fpo_coeff_correct_advance(coeff_correct_calc,
    &conf_range, &phase_range, fpo_moms, boundary_corrections,
    moms, drag_diff_coeff_corrs, drag_coeff, drag_coeff_surf,
    diff_coeff, diff_coeff_surf, 0);

  const double *fpo_moms_c = gkyl_array_cfetch(fpo_moms, 1);
  const double *bcorr_c = gkyl_array_cfetch(boundary_corrections, 1);
  const double *coeff_corrs_c = gkyl_array_cfetch(drag_diff_coeff_corrs, 1);

  if (NV == 4) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.1685375343585043e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -4.1710971760487405e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -2.5777990853015353e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -2.4714311593325799e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -4.0939474033052647e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 2.3189744101340079e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[9], 2.0443723328499349e+01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[10], 2.2091690949038620e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.7187609281138045e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], -4.9965538495170426e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], -5.2922618544448685e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 5.8209395583300142e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], 1.0096632731271260e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], 1.3085323203632530e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], -3.8875043712455232e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], -8.5316721167048962e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 1.6228429243719039e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 1.2205121280494322e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -3.0466081046842675e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[16], 1.1163350799398796e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[18], -2.7105054312137611e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[19], 1.2370836834802779e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[21], -2.6421262074837537e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[22], -2.8346836004093335e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], -1.2348316747910406e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 1.8586360148702548e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 2.5544463610232715e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 3.2846815926326794e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 3.9897666484364568e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -2.2420438408924226e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[9], -6.7917415582715881e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[10], -8.3754623025259556e-17, 1e-12) );
  }
  if (NV == 8) {
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[0], 1.1500173245089075e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[1], -1.3814606337380109e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[3], -2.4654757402320371e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[4], -2.6887040032578782e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[6], -4.6403852982379590e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[7], 1.9031362467424888e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[9], 2.0436859018474077e+01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(fpo_moms_c[10], 9.3620119931091115e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[0], 9.9720462549479656e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[1], 2.9401345535727655e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[3], -3.7608262858090935e-19, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[4], 1.0130431394984520e-19, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[6], -3.0408482806429382e-19, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[7], -5.0549075583401925e-20, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[9], -3.9888185019791869e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[10], -1.1963559167756827e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[12], 1.6589403275025838e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[13], 1.0962162989424440e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[15], -1.8458541986565713e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[16], 6.7388104875646603e-18, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[18], -1.1519648082658485e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[19], 1.6166989952224706e-19, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[21], -5.1860559890122271e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(bcorr_c[22], -4.3842122124302845e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[0], 2.4172048297609167e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[1], 2.2824915248939008e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[3], 2.2831057173690739e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[4], 2.7319520996301737e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[6], 4.5082533140980421e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[7], -1.8857979951717035e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[9], -6.8722195535544301e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(coeff_corrs_c[10], -5.1851480741430756e-16, 1e-12) );
  }
}

void test_1x3v_square_p1_4() { test_1x3v_square_p1(4); }
void test_1x3v_square_p1_8() { test_1x3v_square_p1(8); }
void test_1x3v_bump_p1_4() { test_1x3v_bump_p1(4); }
void test_1x3v_bump_p1_8() { test_1x3v_bump_p1(8); }
void test_1x3v_square_p2_4() { test_1x3v_square_p2(4); }
void test_1x3v_square_p2_8() { test_1x3v_square_p2(8); }
void test_1x3v_bump_p2_4() { test_1x3v_bump_p2(4); }
void test_1x3v_bump_p2_8() { test_1x3v_bump_p2(8); }

TEST_LIST = {
  { "test_1x3v_square_p1_4", test_1x3v_square_p1_4 },
  { "test_1x3v_square_p1_8", test_1x3v_square_p1_8 },
  { "test_1x3v_bump_p1_4", test_1x3v_bump_p1_4 },
  { "test_1x3v_bump_p1_8", test_1x3v_bump_p1_8 },
  { "test_1x3v_square_p2_4", test_1x3v_square_p2_4 },
  { "test_1x3v_square_p2_8", test_1x3v_square_p2_8 },
  { "test_1x3v_bump_p2_4", test_1x3v_bump_p2_4 },
  { "test_1x3v_bump_p2_8", test_1x3v_bump_p2_8 },
  { NULL, NULL }
};
