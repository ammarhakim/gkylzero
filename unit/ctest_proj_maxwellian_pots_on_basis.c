#include <acutest.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <gkyl_fpo_proj_maxwellian_pots_on_basis.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_proj_on_basis.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void eval_M0(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void eval_udrift(double t, const double *xn, double* restrict fout, void *ctx)
{
  fout[0] = 0.5;
  fout[1] = 0.5;
  fout[2] = 0.0;
}

void eval_vtsq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = 1.0;
  fout[0] = vtsq;
}

void eval_gamma(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

void
test_1x3v(int poly_order)
{
  double lower[] = {-1.0, -6.0, -6.0, -6.0};
  double upper[] = {1.0, 6.0, 6.0, 6.0};
  int cells[] = {2, 16, 16, 16};
  int cdim = 1;
  int vdim = 3;
  int pdim = cdim + vdim;

  double confLower[] = {lower[0]};
  double confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);

  // basis functions
  struct gkyl_basis basis, confBasis, surfBasis;
  if (poly_order == 1) {
    gkyl_cart_modal_hybrid(&basis, cdim, vdim);
    gkyl_cart_modal_hybrid(&surfBasis, cdim, vdim-1);
  } else {
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);
    gkyl_cart_modal_serendip(&surfBasis, pdim-1, poly_order);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 0 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0, 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // create moment and potential arrays
  struct gkyl_array *m0, *u_drift, *vtsq, *prim_moms, *gamma;
  struct gkyl_array *fpo_h, *fpo_g, *fpo_h_surf, *fpo_g_surf;
  struct gkyl_array *fpo_dhdv_surf, *fpo_dgdv_surf, *fpo_d2gdv2_surf;
  m0 = mkarr(confBasis.num_basis, confLocal.volume);
  u_drift = mkarr(vdim*confBasis.num_basis, confLocal.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal.volume);
  prim_moms = mkarr((vdim+1)*confBasis.num_basis, confLocal.volume);
  gamma = mkarr(confBasis.num_basis, confLocal.volume);

  fpo_h = mkarr(basis.num_basis, local.volume);
  fpo_g = mkarr(basis.num_basis, local.volume);
 
  fpo_h_surf = mkarr(3*surfBasis.num_basis, local.volume);
  fpo_g_surf = mkarr(3*surfBasis.num_basis, local.volume); 
  fpo_dhdv_surf = mkarr(3*surfBasis.num_basis, local.volume);
  fpo_dgdv_surf = mkarr(9*surfBasis.num_basis, local.volume); 
  fpo_d2gdv2_surf = mkarr(9*surfBasis.num_basis, local.volume); 

  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_M0, NULL);
  gkyl_proj_on_basis *proj_u_drift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq, NULL);
  gkyl_proj_on_basis *proj_gamma = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_gamma, NULL);

  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0);
  gkyl_proj_on_basis_advance(proj_u_drift, 0.0, &confLocal, u_drift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);
  gkyl_proj_on_basis_advance(proj_gamma, 0.0, &confLocal, gamma);

  gkyl_array_set_offset(prim_moms, 1., u_drift, 0);
  gkyl_array_set_offset(prim_moms, 1., vtsq, vdim*confBasis.num_basis);

  // initialize potential projection updater
  gkyl_proj_maxwellian_pots_on_basis *mpob;
  int num_quad = poly_order+1;
  mpob = gkyl_proj_maxwellian_pots_on_basis_new(&grid, &confBasis, &basis, num_quad);
  gkyl_proj_maxwellian_pots_on_basis_lab_mom(mpob, &local, &confLocal, m0, prim_moms, fpo_h, 
    fpo_g, fpo_h_surf, fpo_g_surf, fpo_dhdv_surf, fpo_dgdv_surf, fpo_d2gdv2_surf);

  // Write potentials to files
  char fname_H[1024], fname_G[1024];
  sprintf(fname_H, "ctest_proj_max_pots_H_p%d.gkyl", poly_order);
  sprintf(fname_G, "ctest_proj_max_pots_G_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fpo_h, fname_H);
  gkyl_grid_sub_array_write(&grid, &local, fpo_g, fname_G);

  int idxc[GKYL_MAX_DIM] = {1, 4, 4, 4};
  const double *h_d = gkyl_array_cfetch(fpo_h, gkyl_range_idx(&local, idxc));

  // values to compare at index (1, 4, 4, 4)
  double p1_vals[] = {
    6.2150679410041698e-01, 7.7260531932436161e-17, 1.2598270479977382e-02,
    1.2598270479977382e-02, 1.0975947323022133e-02, 7.3258164081454418e-18,
    -6.5519713996690150e-18, 7.6568765677525823e-04, 3.8692250423821360e-19,
    6.6734911584838039e-04
  };
  double p2_vals[] = {
    6.2150679410041709e-01, -3.9131161402363641e-17, 1.2598270479977382e-02,
    1.2598270479977391e-02, 1.0975947323022123e-02, -1.1444419625373463e-17,
    -1.0360787695126197e-18, 7.6568765677526604e-04, 7.2242128683151019e-19,
    6.6734911584839080e-04
  };

  if (poly_order == 1) {
    for (int i=0; i<10; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], h_d[i], 1e-12) );
  }

  if (poly_order == 2) {
    for (int i=0; i<10; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], h_d[i], 1e-12) );
  }

  gkyl_array_release(m0);
  gkyl_array_release(u_drift);
  gkyl_array_release(vtsq);
  gkyl_array_release(prim_moms);
  gkyl_array_release(gamma);
  gkyl_array_release(fpo_h);
  gkyl_array_release(fpo_g);
  gkyl_array_release(fpo_h_surf);
  gkyl_array_release(fpo_g_surf);
  gkyl_array_release(fpo_dhdv_surf);
  gkyl_array_release(fpo_dgdv_surf);
  gkyl_array_release(fpo_d2gdv2_surf);

  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_u_drift);
  gkyl_proj_on_basis_release(proj_vtsq);
  gkyl_proj_on_basis_release(proj_gamma);

  gkyl_proj_maxwellian_pots_on_basis_release(mpob);
}

void test_1x3v_p1() { test_1x3v(1); }
void test_1x3v_p2() { test_1x3v(2); }

TEST_LIST = {
  { "test_1x3v_p1", test_1x3v_p1 },
  { "test_1x3v_p2", test_1x3v_p2 },
  { NULL, NULL },
};

