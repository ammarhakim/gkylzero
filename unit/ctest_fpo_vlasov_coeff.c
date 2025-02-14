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
  gkyl_grid_sub_array_write(&phase_grid, &phase_range, 0, drag_coeff, fileNm);
  snprintf(fileNm, sizeof(fileNm), fmt_diff, poly_order, NV);
  gkyl_grid_sub_array_write(&phase_grid, &phase_range, 0, diff_coeff, fileNm);

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
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.109578022475359, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.000626713810683, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.021289992906788, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.021289992906788, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[40], 0.599703346207891, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[41], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[42], 0.189273890738643, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[43], -0.173866627312870, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[44], 0.189273890738643, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.410099988925513, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], -0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.075275719022201, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.002732242234505, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.002732242234505, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[40], -0.188390224238438, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[41], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[42], 0.001546275427632, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[43], -0.006946834285642, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[44], -0.033461567239071, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[120], -0.188390224238438, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[121], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[122], -0.006946834285639, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[123], 0.001546275427633, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[124], -0.033461567239070, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.136136091845050, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 0.000000000000002, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.322791855217647, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.077070123529149, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.077070123529150, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[40], -0.188816451826859, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[41], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[42], 0.073518832495259, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[43], 0.073518832495259, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[44], -0.041452406621217, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[120], -0.188816451826859, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[121], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[122], 0.073518832495257, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[123], 0.073518832495258, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[124], -0.041452406621217, 1e-12) );
  }
  if ((poly_order == 1) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.080437908204001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.000034791321315, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.006644645120172, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.006644645120172, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], 0.003504550358724, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.001437056889635, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.000000000000000, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.351836452456254, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 0.000000000000008, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.028049423283142, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.000561896233837, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.000561896233847, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], -0.000651059231579, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], -0.000818114814429, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], 0.000000000000034, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.003269344998163, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.006818263469524, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], -0.000171167653601, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.000000000000018, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 0.492716521570660, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], -0.000000000000012, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.052728454298929, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.002209587664696, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.002209587664694, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -0.001393838600309, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], -0.001393838600310, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.000000000000002, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], 0.001805116419436, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], -0.000009080754942, 1e-8) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.000189185210018, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], 0.000000000000003, 1e-12) );
  }
  if ((poly_order == 2) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.109578022475359, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.000626713810683, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.021289992906788, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.021289992906788, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.599703346207891, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.189273890738643, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], -0.173866627312870, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.189273890738643, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.410099988925508, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.075275719022243, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.002732242234508, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.002732242234506, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -0.188390224238438, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 0.001546275427634, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -0.006946834285631, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -0.033461567239071, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.188390224238438, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.006946834285631, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 0.001546275427634, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.033461567239071, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.136136091845041, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.322791855217654, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.077070123529149, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.077070123529149, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.188816451826858, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 0.073518832495255, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.073518832495256, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.041452406621217, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -0.188816451826858, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 0.073518832495256, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.073518832495257, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -0.041452406621217, 1e-12) );
  }
  if ((poly_order == 2) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.080437908204001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.000034791321315, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.006644645120172, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.006644645120172, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.157677432980284, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.018263418949812, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], 0.000172456709639, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.018263418949812, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.351836452456224, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.028049423283322, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.000561896233836, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.000561896233841, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -0.167334137840670, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], 0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 0.000514210743318, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -0.002534916638363, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -0.013205693830536, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.167334137840670, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.002534916638374, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 0.000514210743324, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.013205693830536, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 0.492716521570680, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 0.000000000000009, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.052728454298940, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.002209587664691, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.002209587664695, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.217823515539350, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 0.001994932779555, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.001994932779561, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.023045715258556, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -0.217823515539349, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 0.001994932779561, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.001994932779563, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -0.023045715258556, 1e-12) );
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
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.109578022475359, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.000626713810683, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.021289992906788, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.021289992906788, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[40], 0.599703346207891, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[41], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[42], 0.189273890738643, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[43], -0.173866627312870, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[44], 0.189273890738643, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.410099988925513, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], -0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.075275719022201, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.002732242234505, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.002732242234505, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[40], -0.188390224238438, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[41], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[42], 0.001546275427632, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[43], -0.006946834285642, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[44], -0.033461567239071, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[120], -0.188390224238438, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[121], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[122], -0.006946834285639, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[123], 0.001546275427633, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[124], -0.033461567239070, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.136136091845050, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 0.000000000000002, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.322791855217647, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.077070123529149, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.077070123529150, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[40], -0.188816451826859, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[41], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[42], 0.073518832495259, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[43], 0.073518832495259, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[44], -0.041452406621217, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[120], -0.188816451826859, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[121], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[122], 0.073518832495257, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[123], 0.073518832495258, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[124], -0.041452406621217, 1e-12) );
  }
  if ((poly_order == 1) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.080437908204001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.000034791321315, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.006644645120172, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.006644645120172, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], 0.003504550358724, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.001437056889635, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.000000000000000, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.351836452456254, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 0.000000000000008, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.028049423283142, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.000561896233837, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.000561896233847, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], -0.000651059231579, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], -0.000818114814429, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], 0.000000000000034, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.003269344998163, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.006818263469524, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], -0.000171167653601, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.000000000000018, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 0.492716521570660, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], -0.000000000000012, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.052728454298929, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.002209587664696, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.002209587664694, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -0.001393838600309, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], -0.001393838600310, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.000000000000002, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], 0.001805116419436, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], -0.000009080754942, 1e-8) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.000189185210018, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], 0.000000000000003, 1e-12) );
  }
  if ((poly_order == 2) && (NV == 4)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.109578022475359, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.000626713810683, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.021289992906788, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.021289992906788, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.599703346207891, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.189273890738643, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], -0.173866627312870, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.189273890738643, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.410099988925508, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.075275719022243, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.002732242234508, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.002732242234506, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -0.188390224238438, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 0.001546275427634, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -0.006946834285631, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -0.033461567239071, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.188390224238438, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.006946834285631, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 0.001546275427634, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.033461567239071, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 1.136136091845041, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.322791855217654, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.077070123529149, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.077070123529149, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.188816451826858, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 0.073518832495255, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.073518832495256, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.041452406621217, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -0.188816451826858, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 0.073518832495256, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.073518832495257, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -0.041452406621217, 1e-12) );
  }
  if ((poly_order == 2) && (NV == 8)) {
    TEST_CHECK( gkyl_compare_double(drag_corner[0], 0.080437908204001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[1], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[2], 0.000034791321315, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_corner[3], 0.006644645120172, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_corner[4], 0.006644645120172, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[48], 0.157677432980284, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[49], 0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[50], 0.018263418949812, 1e-12) );
    TEST_CHECK( gkyl_compare_double(drag_int[51], 0.000172456709639, 1e-10) );
    TEST_CHECK( gkyl_compare_double(drag_int[52], 0.018263418949812, 1e-12) );

    TEST_CHECK( gkyl_compare_double(diff_corner[0], 0.351836452456224, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[1], 0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[2], 0.028049423283322, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[3], 0.000561896233836, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[4], 0.000561896233841, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[48], -0.167334137840670, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[49], 0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[50], 0.000514210743318, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[51], -0.002534916638363, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[52], -0.013205693830536, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[144], -0.167334137840670, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[145], 0.000000000000001, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[146], -0.002534916638374, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_corner[147], 0.000514210743324, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_corner[148], -0.013205693830536, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[0], 0.492716521570680, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[1], 0.000000000000009, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[2], 0.052728454298940, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[3], 0.002209587664691, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[4], 0.002209587664695, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[48], -0.217823515539350, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[49], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[50], 0.001994932779555, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[51], 0.001994932779561, 1e-10) );
    TEST_CHECK( gkyl_compare_double(diff_int[52], -0.023045715258556, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[144], -0.217823515539349, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[145], -0.000000000000000, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[146], 0.001994932779561, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[147], 0.001994932779563, 1e-12) );
    TEST_CHECK( gkyl_compare_double(diff_int[148], -0.023045715258556, 1e-12) );
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
