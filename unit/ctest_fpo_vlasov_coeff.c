// Tests for the Vlasov FPO Maxwellian Rosenbluth potential projection and
// drag and diffusion coefficient computation

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
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>

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

void eval_gamma(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  fout[0] = app->gamma0;
}

void eval_lte_moms(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct fpo_ctx *app = ctx;
  fout[0] = app->n0;
  fout[1] = app->ux0;
  fout[2] = app->uy0;
  fout[3] = app->uz0;
  fout[4] = app->vth0*app->vth0;
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

void test_1x3v(int poly_order, int NV)
{
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {1, NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 5.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};

  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // initialize basis
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  if (poly_order == 1) {
    gkyl_cart_modal_hybrid(&phase_basis, cdim, vdim);
    gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
    gkyl_cart_modal_hybrid(&surf_basis, cdim, vdim-1);
  }
  else {
    gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
    gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
  }

  gkyl_proj_on_basis *proj_gamma = gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order+1, 1, eval_gamma, &ctx);
  gkyl_proj_on_basis *proj_lte = gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order+1, 5, eval_lte_moms, &ctx);

  struct gkyl_array *lte_moms, *gamma, *h, *g, *h_surf, *g_surf, *dhdv_surf, *dgdv_surf, *d2gdv2_surf;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  lte_moms = mkarr(5*conf_basis.num_basis, conf_range_ext.volume);
  gamma = mkarr(conf_basis.num_basis, conf_range_ext.volume);
  h = mkarr(phase_basis.num_basis, phase_range_ext.volume);
  g = mkarr(phase_basis.num_basis, phase_range_ext.volume);
  h_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  g_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  dhdv_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  dgdv_surf = mkarr(2*vdim*surf_basis.num_basis, phase_range_ext.volume);
  d2gdv2_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  drag_coeff = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume);
  drag_coeff_surf = mkarr(vdim*surf_basis.num_basis, phase_range_ext.volume);
  diff_coeff = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume);
  diff_coeff_surf = mkarr(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume);

  // Array of relative offsets for recovery stencils
  long offsets[36];

  // Initialize updaters for potentials and coeffs
  gkyl_proj_maxwellian_pots_on_basis *pot_slvr = gkyl_proj_maxwellian_pots_on_basis_new(
    &phase_grid, &conf_range, &phase_range, &conf_basis, &phase_basis, poly_order+1, 0);

  gkyl_fpo_vlasov_coeff_recovery *coeff_recovery = gkyl_fpo_vlasov_coeff_recovery_new(&phase_grid,
    &phase_basis, &phase_range_ext, offsets, 0);

  // Project moments and compute potentials
  gkyl_proj_on_basis_advance(proj_gamma, 0.0, &conf_range, gamma);
  gkyl_proj_on_basis_advance(proj_lte, 0.0, &conf_range, lte_moms);
  gkyl_proj_maxwellian_pots_on_basis_advance(pot_slvr, &phase_range, &conf_range, 
    lte_moms, h, g, h_surf, g_surf,
    dhdv_surf, dgdv_surf, d2gdv2_surf);

  // Compute drag and diffusion coefficients
  gkyl_calc_fpo_drag_coeff_recovery(coeff_recovery, &phase_grid, phase_basis, &phase_range,
    &conf_range, gamma, h, dhdv_surf, drag_coeff, drag_coeff_surf, 0);

  gkyl_calc_fpo_diff_coeff_recovery(coeff_recovery, &phase_grid, phase_basis, 
    &phase_range, &conf_range, gamma,
    g, g_surf, dgdv_surf, d2gdv2_surf, 
    diff_coeff, diff_coeff_surf, 0);

  const char *fmt_drag = "ctest_fpo_drag_coeff_p%d_%d.gkyl";
  const char *fmt_diff = "ctest_fpo_diff_coeff_p%d_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt_drag, NV);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof(fileNm), fmt_drag, poly_order, NV);
  gkyl_grid_sub_array_write(&phase_grid, &phase_range, drag_coeff, fileNm);
  snprintf(fileNm, sizeof(fileNm), fmt_diff, poly_order, NV);
  gkyl_grid_sub_array_write(&phase_grid, &phase_range, diff_coeff, fileNm);

  // Compare values in corner cell and interior cell (of velocity space)
  // Checking first five components of a_x, a_y, D_xx, D_xy, and D_yx
  int idx_corner[] = {1, 1, 1, 1};
  int idx_int[] = {1, 2, 2, 2};
  long lin_corner = gkyl_range_idx(&phase_range, idx_corner);
  long lin_int = gkyl_range_idx(&phase_range, idx_int);

  double *drag_corner, *drag_int;
  double *diff_corner, *diff_int;
  drag_corner = gkyl_array_fetch(drag_coeff, lin_corner);
  drag_int = gkyl_array_fetch(drag_coeff, lin_int);
  diff_corner = gkyl_array_fetch(diff_coeff, lin_corner);
  diff_int = gkyl_array_fetch(diff_coeff, lin_int);

  if ((poly_order == 1) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 1.0957802247535949e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -3.3267081287092361e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 6.2671381068314137e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 2.1289992906787691e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 2.1289992906787663e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[40], 1.0957802247535935e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[41], -9.1780687842145553e-20, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[42], 2.1289992906787722e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[43], 6.2671381068322811e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[44], 2.1289992906787684e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 5.9970334620789123e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], -1.0681454292690527e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], -1.7386662731286925e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 1.8927389073864298e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 1.8927389073864298e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[40], 5.9970334620789101e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[41], -9.4689458638985962e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[42], 1.8927389073864304e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[43], -1.7386662731286939e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[44], 1.8927389073864301e-01, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 4.1009998892550925e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], -6.0484480028533050e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 7.5275719022215526e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 2.7322422345082420e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 2.7322422345070763e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[40], -1.8839022423843843e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[41], 2.3629602838952619e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[42], 1.5462754276315768e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[43], -6.9468342856394199e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[44], -3.3461567239070754e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[120], -1.8839022423843754e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[121], -1.5696515912891607e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[122], -6.9468342856394199e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[123], 1.5462754276330340e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[124], -3.3461567239070303e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.1361360918450480e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 2.2178533940871595e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 3.2279185521764703e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 7.7070123529148868e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 7.7070123529149312e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[40], -1.8881645182685905e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[41], 1.1525237590321227e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[42], 7.3518832495257314e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[43], 7.3518832495257244e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[44], -4.1452406621216631e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[120], -1.8881645182685905e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[121], 1.1525237590321225e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[122], 7.3518832495257244e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[123], 7.3518832495257980e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[124], -4.1452406621216596e-02, 1e-12) );
  }
  if ((poly_order == 1) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 8.0437908204001179e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 8.0952120429861365e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 3.4791321315376378e-05, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 6.6446451201719334e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 6.6446451201719638e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[40], 8.0437908204000944e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[41], 7.4614544409679005e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[42], 6.6446451201719473e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[43], 3.4791321315482196e-05, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[44], 6.6446451201719768e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 1.5767743298028405e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 7.7056929110321426e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], 1.7245670963886850e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 1.8263418949811702e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 1.8263418949811674e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[40], 1.5767743298028425e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[41], 8.5564459289356172e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[42], 1.8263418949811708e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[43], 1.7245670963899064e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[44], 1.8263418949811702e-02, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 3.5183645245625783e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 7.8518510325525888e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 2.8049423283107444e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 5.6189623383584575e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 5.6189623384617082e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[40], -1.6733413784066911e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[41], 9.4931002501411847e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[42], 5.1421074331812729e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[43], -2.5349166383959246e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[44], -1.3205693830535516e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[120], -1.6733413784067111e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[121], -1.9208246714452471e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[122], -2.5349166383959246e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[123], 5.1421074331651953e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[124], -1.3205693830535881e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 4.9271652157063811e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], -1.5535631812761674e-14, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 5.2728454298921860e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 2.2095876646950788e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 2.2095876646944478e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[40], -2.1782351553934687e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[41], -6.6317305234838796e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[42], 1.9949327795609818e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[43], 1.9949327795785372e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[44], -2.3045715258555029e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[120], -2.1782351553934809e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[121], -6.6317305234838796e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[122], 1.9949327795687906e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[123], 1.9949327795698245e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[124], -2.3045715258555054e-02, 1e-12) );
  }
  if ((poly_order == 2) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 1.0957802247535955e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -8.7261341959348119e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 6.2671381068308239e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 2.1289992906787701e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 2.1289992906787698e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[48], 1.0957802247535939e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[49], -9.6252648582924661e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[50], 2.1289992906787708e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[51], 6.2671381068341372e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[52], 2.1289992906787684e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 5.9970334620789090e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 2.3329030243168463e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], -1.7386662731286934e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 1.8927389073864292e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 1.8927389073864287e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 5.9970334620789056e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], -3.7342583274017147e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 1.8927389073864290e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], -1.7386662731286939e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 1.8927389073864287e-01, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 4.1009998892551014e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 6.1261464547349370e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 7.5275719022222631e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 2.7322422345077424e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 2.7322422345064101e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -1.8839022423843735e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], 1.2575542081729871e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 1.5462754276341141e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -6.9468342856310655e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -3.3461567239071108e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -1.8839022423843774e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 6.2951722469947693e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -6.9468342856310655e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 1.5462754276338881e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -3.3461567239071038e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.1361360918450383e+00, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], -8.7971390732944644e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 3.2279185521765386e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 7.7070123529148257e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 7.7070123529149104e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -1.8881645182685813e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -1.5235998965191322e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 7.3518832495254510e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 7.3518832495256550e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -4.1452406621216797e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -1.8881645182685813e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -1.5235998965191346e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 7.3518832495255509e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 7.3518832495255718e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -4.1452406621216763e-02, 1e-12) );
  }
  if ((poly_order == 2) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 8.0437908204000333e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 1.3923949659298653e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 3.4791321315259284e-05, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 6.6446451201719629e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 6.6446451201719560e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[48], 8.0437908204000944e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[49], 6.7501580973303227e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[50], 6.6446451201719221e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[51], 3.4791321315502145e-05, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[52], 6.6446451201719664e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 1.5767743298028425e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 1.5657358927561106e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], 1.7245670963813022e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 1.8263418949811636e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 1.8263418949811653e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 1.5767743298028414e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], 4.9931040519873304e-17, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 1.8263418949811653e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], 1.7245670963879635e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 1.8263418949811663e-02, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 3.5183645245622230e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 3.0193105866998722e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 2.8049423283334818e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 5.6189623383440246e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 5.6189623383851028e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -1.6733413784067072e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], -3.2279542438494002e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 5.1421074332042591e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -2.5349166383736468e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -1.3205693830536195e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -1.6733413784066992e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 2.9522573920776913e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -2.5349166383847855e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 5.1421074332355361e-04, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -1.3205693830536233e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 4.9271652157067025e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 9.6478421124524351e-15, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 5.2728454298930957e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 2.2095876646910642e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 2.2095876646944391e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -2.1782351553935089e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -2.9026776924397315e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 1.9949327795625894e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 1.9949327795660060e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -2.3045715258555935e-02, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -2.1782351553935128e-01, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -2.9026776924397320e-16, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 1.9949327795618288e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 1.9949327795633935e-03, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -2.3045715258555859e-02, 1e-12) );
  }

  // Release memory
  gkyl_array_release(lte_moms);
  gkyl_array_release(gamma);
  gkyl_array_release(h);
  gkyl_array_release(g);
  gkyl_array_release(h_surf);
  gkyl_array_release(g_surf);
  gkyl_array_release(dhdv_surf);
  gkyl_array_release(dgdv_surf);
  gkyl_array_release(d2gdv2_surf);
  gkyl_array_release(drag_coeff);
  gkyl_array_release(drag_coeff_surf);
  gkyl_array_release(diff_coeff);
  gkyl_array_release(diff_coeff_surf);

  gkyl_proj_maxwellian_pots_on_basis_release(pot_slvr);
  gkyl_fpo_vlasov_coeff_recovery_release(coeff_recovery);
}

void test_1x3v_cu(int poly_order, int NV)
{
  bool use_gpu = true;
  int cdim = 1, vdim = 3;
  int pdim = cdim+vdim;

  struct fpo_ctx ctx = create_ctx();

  int cells[] = {1, NV, NV, NV};
  int ghost[] = {0, 0, 0, 0};

  double L = 5.0;
  double lower[] = {0.0, -L, -L, -L};
  double upper[] = {1.0, L, L, L};

  struct gkyl_rect_grid conf_grid;
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &conf_range_ext, &conf_range);

  struct gkyl_rect_grid phase_grid;
  struct gkyl_range phase_range, phase_range_ext;
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&phase_grid, ghost, &phase_range_ext, &phase_range);

  // initialize basis
  struct gkyl_basis phase_basis, conf_basis, surf_basis;

  if (poly_order == 1) {
    gkyl_cart_modal_hybrid(&phase_basis, cdim, vdim);
    gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
    gkyl_cart_modal_hybrid(&surf_basis, cdim, vdim-1);
  }
  else {
    gkyl_cart_modal_serendip(&phase_basis, pdim, poly_order);
    gkyl_cart_modal_serendip(&conf_basis, cdim, poly_order);
    gkyl_cart_modal_serendip(&surf_basis, pdim-1, poly_order);
  }

  gkyl_proj_on_basis *proj_gamma = gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order+1, 1, eval_gamma, &ctx);
  gkyl_proj_on_basis *proj_lte = gkyl_proj_on_basis_new(&conf_grid, &conf_basis, poly_order+1, 5, eval_lte_moms, &ctx);

  struct gkyl_array *lte_moms, *gamma, *lte_moms_ho, *gamma_ho;
  struct gkyl_array *h, *g, *h_surf, *g_surf, *dhdv_surf, *dgdv_surf, *d2gdv2_surf;
  struct gkyl_array *drag_coeff, *drag_coeff_surf, *diff_coeff, *diff_coeff_surf;
  lte_moms = mkarr_cu(5*conf_basis.num_basis, conf_range_ext.volume);
  gamma = mkarr_cu(conf_basis.num_basis, conf_range_ext.volume);
  lte_moms_ho = mkarr(5*conf_basis.num_basis, conf_range_ext.volume);
  gamma_ho = mkarr(conf_basis.num_basis, conf_range_ext.volume);
  h = mkarr_cu(phase_basis.num_basis, phase_range_ext.volume);
  g = mkarr_cu(phase_basis.num_basis, phase_range_ext.volume);
  h_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  g_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  dhdv_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  dgdv_surf = mkarr_cu(2*vdim*surf_basis.num_basis, phase_range_ext.volume);
  d2gdv2_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  drag_coeff = mkarr_cu(vdim*phase_basis.num_basis, phase_range_ext.volume);
  drag_coeff_surf = mkarr_cu(vdim*surf_basis.num_basis, phase_range_ext.volume);
  diff_coeff = mkarr_cu(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume);
  diff_coeff_surf = mkarr_cu(2*vdim*vdim*surf_basis.num_basis, phase_range_ext.volume);

  // Array of relative offsets for recovery stencils
  long offsets[36];

  // Initialize updaters for potentials and coeffs
  gkyl_proj_maxwellian_pots_on_basis *pot_slvr = gkyl_proj_maxwellian_pots_on_basis_new(
    &phase_grid, &conf_range, &phase_range, &conf_basis, &phase_basis, poly_order+1, use_gpu);

  gkyl_fpo_vlasov_coeff_recovery *coeff_recovery = gkyl_fpo_vlasov_coeff_recovery_new(&phase_grid,
    &phase_basis, &phase_range_ext, offsets, use_gpu);

  // Project moments and compute potentials
  gkyl_proj_on_basis_advance(proj_gamma, 0.0, &conf_range, gamma_ho);
  gkyl_proj_on_basis_advance(proj_lte, 0.0, &conf_range, lte_moms_ho);
  gkyl_array_copy(gamma, gamma_ho);
  gkyl_array_copy(lte_moms, lte_moms_ho);
  gkyl_proj_maxwellian_pots_on_basis_advance(pot_slvr, &phase_range, &conf_range, 
    lte_moms, h, g, h_surf, g_surf,
    dhdv_surf, dgdv_surf, d2gdv2_surf);

  // Compute drag and diffusion coefficients
  gkyl_calc_fpo_drag_coeff_recovery(coeff_recovery, &phase_grid, phase_basis, &phase_range,
    &conf_range, gamma, h, dhdv_surf, drag_coeff, drag_coeff_surf, use_gpu);

  gkyl_calc_fpo_diff_coeff_recovery(coeff_recovery, &phase_grid, phase_basis, 
    &phase_range, &conf_range, gamma,
    g, g_surf, dgdv_surf, d2gdv2_surf, 
    diff_coeff, diff_coeff_surf, use_gpu);

  struct gkyl_array *drag_coeff_ho = mkarr(vdim*phase_basis.num_basis, phase_range_ext.volume);
  struct gkyl_array *diff_coeff_ho = mkarr(vdim*vdim*phase_basis.num_basis, phase_range_ext.volume);
  gkyl_array_copy(drag_coeff_ho, drag_coeff);
  gkyl_array_copy(diff_coeff_ho, diff_coeff);

  // Compare values in corner cell and interior cell (of velocity space)
  // Checking first five components of a_x, a_y, D_xx, D_xy, and D_yx
  int idx_corner[] = {1, 1, 1, 1};
  int idx_int[] = {1, 2, 2, 2};
  long lin_corner = gkyl_range_idx(&phase_range, idx_corner);
  long lin_int = gkyl_range_idx(&phase_range, idx_int);

  double *drag_corner, *drag_int;
  double *diff_corner, *diff_int;
  drag_corner = gkyl_array_fetch(drag_coeff_ho, lin_corner);
  drag_int = gkyl_array_fetch(drag_coeff_ho, lin_int);
  diff_corner = gkyl_array_fetch(diff_coeff_ho, lin_corner);
  diff_int = gkyl_array_fetch(diff_coeff_ho, lin_int);

  if ((poly_order == 1) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 1.0957802247535949e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -3.3267081287092361e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 6.2671381068314137e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 2.1289992906787691e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 2.1289992906787663e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[40], 1.0957802247535935e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[41], -9.1780687842145553e-20, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[42], 2.1289992906787722e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[43], 6.2671381068322811e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[44], 2.1289992906787684e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 5.9970334620789123e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], -1.0681454292690527e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], -1.7386662731286925e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 1.8927389073864298e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 1.8927389073864298e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[40], 5.9970334620789101e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[41], -9.4689458638985962e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[42], 1.8927389073864304e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[43], -1.7386662731286939e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[44], 1.8927389073864301e-01, 1e-10) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 4.1009998892550925e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], -6.0484480028533050e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 7.5275719022215526e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 2.7322422345082420e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 2.7322422345070763e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[40], -1.8839022423843843e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[41], 2.3629602838952619e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[42], 1.5462754276315768e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[43], -6.9468342856394199e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[44], -3.3461567239070754e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[120], -1.8839022423843754e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[121], -1.5696515912891607e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[122], -6.9468342856394199e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[123], 1.5462754276330340e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[124], -3.3461567239070303e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.1361360918450480e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 2.2178533940871595e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 3.2279185521764703e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 7.7070123529148868e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 7.7070123529149312e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[40], -1.8881645182685905e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[41], 1.1525237590321227e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[42], 7.3518832495257314e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[43], 7.3518832495257244e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[44], -4.1452406621216631e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[120], -1.8881645182685905e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[121], 1.1525237590321225e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[122], 7.3518832495257244e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[123], 7.3518832495257980e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[124], -4.1452406621216596e-02, 1e-10) );
  }
  if ((poly_order == 1) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 8.0437908204001179e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 8.0952120429861365e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 3.4791321315376378e-05, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 6.6446451201719334e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 6.6446451201719638e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[40], 8.0437908204000944e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[41], 7.4614544409679005e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[42], 6.6446451201719473e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[43], 3.4791321315482196e-05, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[44], 6.6446451201719768e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 1.5767743298028405e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 7.7056929110321426e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], 1.7245670963886850e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 1.8263418949811702e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 1.8263418949811674e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[40], 1.5767743298028425e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[41], 8.5564459289356172e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[42], 1.8263418949811708e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[43], 1.7245670963899064e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[44], 1.8263418949811702e-02, 1e-10) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 3.5183645245625783e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 7.8518510325525888e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 2.8049423283107444e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 5.6189623383584575e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 5.6189623384617082e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[40], -1.6733413784066911e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[41], 9.4931002501411847e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[42], 5.1421074331812729e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[43], -2.5349166383959246e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[44], -1.3205693830535516e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[120], -1.6733413784067111e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[121], -1.9208246714452471e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[122], -2.5349166383959246e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[123], 5.1421074331651953e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[124], -1.3205693830535881e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 4.9271652157063811e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], -1.5535631812761674e-14, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 5.2728454298921860e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 2.2095876646950788e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 2.2095876646944478e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[40], -2.1782351553934687e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[41], -6.6317305234838796e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[42], 1.9949327795609818e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[43], 1.9949327795785372e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[44], -2.3045715258555029e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[120], -2.1782351553934809e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[121], -6.6317305234838796e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[122], 1.9949327795687906e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[123], 1.9949327795698245e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[124], -2.3045715258555054e-02, 1e-10) );
  }
  if ((poly_order == 2) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 1.0957802247535955e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -8.7261341959348119e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 6.2671381068308239e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 2.1289992906787701e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 2.1289992906787698e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[48], 1.0957802247535939e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[49], -9.6252648582924661e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[50], 2.1289992906787708e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[51], 6.2671381068341372e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[52], 2.1289992906787684e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 5.9970334620789090e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 2.3329030243168463e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], -1.7386662731286934e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 1.8927389073864292e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 1.8927389073864287e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 5.9970334620789056e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], -3.7342583274017147e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 1.8927389073864290e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], -1.7386662731286939e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 1.8927389073864287e-01, 1e-10) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 4.1009998892551014e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 6.1261464547349370e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 7.5275719022222631e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 2.7322422345077424e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 2.7322422345064101e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -1.8839022423843735e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], 1.2575542081729871e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 1.5462754276341141e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -6.9468342856310655e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -3.3461567239071108e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -1.8839022423843774e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 6.2951722469947693e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -6.9468342856310655e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 1.5462754276338881e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -3.3461567239071038e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.1361360918450383e+00, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], -8.7971390732944644e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 3.2279185521765386e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 7.7070123529148257e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 7.7070123529149104e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -1.8881645182685813e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -1.5235998965191322e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 7.3518832495254510e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 7.3518832495256550e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -4.1452406621216797e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -1.8881645182685813e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -1.5235998965191346e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 7.3518832495255509e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 7.3518832495255718e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -4.1452406621216763e-02, 1e-10) );
  }
  if ((poly_order == 2) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 8.0437908204000333e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 1.3923949659298653e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 3.4791321315259284e-05, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 6.6446451201719629e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 6.6446451201719560e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[48], 8.0437908204000944e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[49], 6.7501580973303227e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[50], 6.6446451201719221e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[51], 3.4791321315502145e-05, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[52], 6.6446451201719664e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[0], 1.5767743298028425e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[1], 1.5657358927561106e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[2], 1.7245670963813022e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[3], 1.8263418949811636e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[4], 1.8263418949811653e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 1.5767743298028414e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], 4.9931040519873304e-17, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 1.8263418949811653e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], 1.7245670963879635e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 1.8263418949811663e-02, 1e-10) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 3.5183645245622230e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 3.0193105866998722e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 2.8049423283334818e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 5.6189623383440246e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 5.6189623383851028e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -1.6733413784067072e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], -3.2279542438494002e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 5.1421074332042591e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -2.5349166383736468e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -1.3205693830536195e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -1.6733413784066992e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 2.9522573920776913e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -2.5349166383847855e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 5.1421074332355361e-04, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -1.3205693830536233e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 4.9271652157067025e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 9.6478421124524351e-15, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 5.2728454298930957e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 2.2095876646910642e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 2.2095876646944391e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -2.1782351553935089e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -2.9026776924397315e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 1.9949327795625894e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 1.9949327795660060e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -2.3045715258555935e-02, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -2.1782351553935128e-01, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -2.9026776924397320e-16, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 1.9949327795618288e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 1.9949327795633935e-03, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -2.3045715258555859e-02, 1e-10) );
  }
  // Release memory
  gkyl_array_release(lte_moms);
  gkyl_array_release(gamma);
  gkyl_array_release(h);
  gkyl_array_release(g);
  gkyl_array_release(h_surf);
  gkyl_array_release(g_surf);
  gkyl_array_release(dhdv_surf);
  gkyl_array_release(dgdv_surf);
  gkyl_array_release(d2gdv2_surf);
  gkyl_array_release(drag_coeff);
  gkyl_array_release(drag_coeff_surf);
  gkyl_array_release(diff_coeff);
  gkyl_array_release(diff_coeff_surf);

  gkyl_array_release(drag_coeff_ho);
  gkyl_array_release(diff_coeff_ho);

  gkyl_proj_maxwellian_pots_on_basis_release(pot_slvr);
  gkyl_fpo_vlasov_coeff_recovery_release(coeff_recovery);
}

void test_1x3v_p1_4() { test_1x3v(1, 4); }
void test_1x3v_p1_8() { test_1x3v(1, 8); }
void test_1x3v_p2_4() { test_1x3v(2, 4); }
void test_1x3v_p2_8() { test_1x3v(2, 8); }
void test_1x3v_p1_4_cu() { test_1x3v_cu(1, 4); }
void test_1x3v_p1_8_cu() { test_1x3v_cu(1, 8); }
void test_1x3v_p2_4_cu() { test_1x3v_cu(2, 4); }
void test_1x3v_p2_8_cu() { test_1x3v_cu(2, 8); }

TEST_LIST = {
  { "test_1x3v_p1_4", test_1x3v_p1_4 },
  { "test_1x3v_p1_8", test_1x3v_p1_8 },
  { "test_1x3v_p2_4", test_1x3v_p2_4 },
  { "test_1x3v_p2_8", test_1x3v_p2_8 },
  #ifdef GKYL_HAVE_CUDA
  { "test_1x3v_p1_4_cu", test_1x3v_p1_4_cu },
  { "test_1x3v_p1_8_cu", test_1x3v_p1_8_cu },
  { "test_1x3v_p2_4_cu", test_1x3v_p2_4_cu },
  { "test_1x3v_p2_8_cu", test_1x3v_p2_8_cu },
  #endif
  { NULL, NULL }
};
