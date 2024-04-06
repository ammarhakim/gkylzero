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
  struct gkyl_array *fpo_h, *fpo_g, *fpo_dhdv, *fpo_d2gdv2, *fpo_h_surf, *fpo_g_surf;
  struct gkyl_array *fpo_dhdv_surf, *fpo_dgdv_surf, *fpo_d2gdv2_surf;
  m0 = mkarr(confBasis.num_basis, confLocal.volume);
  u_drift = mkarr(vdim*confBasis.num_basis, confLocal.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal.volume);
  prim_moms = mkarr((vdim+1)*confBasis.num_basis, confLocal.volume);
  gamma = mkarr(confBasis.num_basis, confLocal.volume);

  fpo_h = mkarr(basis.num_basis, local.volume);
  fpo_g = mkarr(basis.num_basis, local.volume);
  fpo_dhdv = mkarr(vdim*basis.num_basis, local.volume);
  fpo_d2gdv2 = mkarr(vdim*vdim*basis.num_basis, local.volume);
 
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

  gkyl_array_set_offset(prim_moms, 1.0, m0, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, u_drift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, vtsq, (vdim+1)*confBasis.num_basis);

  // initialize potential projection updater
  gkyl_proj_maxwellian_pots_on_basis *mpob;
  int num_quad = poly_order+1;
  mpob = gkyl_proj_maxwellian_pots_on_basis_new(&grid, &confBasis, &basis, num_quad);
  gkyl_proj_maxwellian_pots_on_basis_advance(mpob, &local, &confLocal, prim_moms, 
    fpo_h, fpo_g, fpo_h_surf, fpo_g_surf, fpo_dhdv_surf, fpo_dgdv_surf, fpo_d2gdv2_surf);
  gkyl_proj_maxwellian_pots_deriv_on_basis_advance(mpob, &local, &confLocal, prim_moms,
    fpo_dhdv, fpo_d2gdv2);

  // Write potentials to files
  char fname_H[1024], fname_G[1024], fname_dHdv[1024], fname_d2Gdv2[1024];
  sprintf(fname_H, "ctest_proj_max_pots_H_p%d.gkyl", poly_order);
  sprintf(fname_G, "ctest_proj_max_pots_G_p%d.gkyl", poly_order);
  sprintf(fname_dHdv, "ctest_proj_max_pots_dHdv_p%d.gkyl", poly_order);
  sprintf(fname_d2Gdv2, "ctest_proj_max_pots_d2Gdv2_p%d.gkyl", poly_order);
  gkyl_grid_sub_array_write(&grid, &local, fpo_h, fname_H);
  gkyl_grid_sub_array_write(&grid, &local, fpo_g, fname_G);
  gkyl_grid_sub_array_write(&grid, &local, fpo_dhdv, fname_dHdv);
  gkyl_grid_sub_array_write(&grid, &local, fpo_d2gdv2, fname_d2Gdv2);

  int idxc[GKYL_MAX_DIM] = {1, 4, 4, 4};
  const double *h_d = gkyl_array_cfetch(fpo_h, gkyl_range_idx(&local, idxc));
  const double *g_d = gkyl_array_cfetch(fpo_g, gkyl_range_idx(&local, idxc));
  const double *dhdv_d = gkyl_array_cfetch(fpo_dhdv, gkyl_range_idx(&local, idxc));
  const double *d2gdv2_d = gkyl_array_cfetch(fpo_d2gdv2, gkyl_range_idx(&local, idxc));

  // values to compare at index (1, 4, 4, 4)
  double p1_H_vals[] = {
    6.2150679410041698e-01, 7.7260531932436161e-17, 1.2598270479977382e-02,
    1.2598270479977382e-02, 1.0975947323022133e-02, 7.3258164081454418e-18,
    -6.5519713996690150e-18, 7.6568765677525823e-04, 3.8692250423821360e-19,
    6.6734911584838039e-04
  };
  double p1_G_vals[] = {
    2.6394564743713552e+01, 5.2445164890279720e-15, -5.0845721303132352e-01,
    -5.0845721303132241e-01, -4.4287347013728501e-01, -2.6731370225071298e-16,
    1.7677550759934962e-16, -9.7763833130525158e-03, 1.2126435636809179e-16,
    -8.5171671066874167e-03
  };
  double p1_dHdv_vals[] = {
    5.8142162630751108e-02, -3.1174110011968858e-18, 2.9155662245575957e-04,
    3.5347334739749213e-03, 3.0807671103840447e-03, -9.5372448325446511e-19,
    -5.2004361426026333e-19, 1.6101648225076210e-04, -8.6362745266061537e-20,
    1.4050794069816070e-04
  };
  double p1_d2Gdv2_vals[] = {
    3.9758218609470125e-01, -3.8675439530490767e-17, 2.2955951238790299e-02,
    -3.2433108729885128e-04, -2.8617646838544433e-04, 1.1949780380790645e-17, 
    1.5414395249298035e-18, 8.6387382233785550e-04, 1.5414395249298035e-18, 
    7.5238867378917028e-04
  };
  double p2_H_vals[] = {
    6.2150679410041709e-01, -3.9131161402363641e-17, 1.2598270479977382e-02,
    1.2598270479977391e-02, 1.0975947323022123e-02, -1.1444419625373463e-17,
    -1.0360787695126197e-18, 7.6568765677526604e-04, 7.2242128683151019e-19,
    6.6734911584839080e-04
  };
  double p2_G_vals[] = {
    2.6394564743713545e+01, 6.0307699679378142e-16, -5.0845721303132241e-01,
    -5.0845721303132196e-01, -4.4287347013728540e-01, 1.4523185095403825e-16,
    3.6727645587906959e-16, -9.7763833130526043e-03, 1.0968261639733243e-17,
    -8.5171671066874496e-03
  };
  double p2_dHdv_vals[] = {
    5.8142162630751094e-02, -2.8464736213248744e-19, 2.9155662245576125e-04,
    3.5347334739749218e-03, 3.0807671103840447e-03, -4.6974766245459400e-19,
    3.9761407553380955e-19, 1.6101648225076224e-04, -1.6219065213807818e-19,
    1.4050794069816159e-04
  };
  double p2_d2Gdv2_vals[] = {
    3.9758218609470114e-01, 5.9068367118966121e-17, 2.2955951238790302e-02,
    -3.2433108729885095e-04, -2.8617646838543836e-04, 4.1035357856357695e-18,
    6.3408883368215515e-19, 8.6387382233785929e-04, 2.0698095261161238e-18,
    7.5238867378917668e-04
  };

  if (poly_order == 1) {
    for (int i=0; i<10; ++i) {
      TEST_CHECK( gkyl_compare_double(p1_H_vals[i], h_d[i], 1e-12) );
      TEST_CHECK( gkyl_compare_double(p1_G_vals[i], g_d[i], 1e-12) );
      TEST_CHECK( gkyl_compare_double(p1_dHdv_vals[i], dhdv_d[i], 1e-12) );
      TEST_CHECK( gkyl_compare_double(p1_d2Gdv2_vals[i], d2gdv2_d[i], 1e-12) );
    }
  }

  if (poly_order == 2) {
    for (int i=0; i<10; ++i) {
      TEST_CHECK( gkyl_compare_double(p2_H_vals[i], h_d[i], 1e-12) );
      TEST_CHECK( gkyl_compare_double(p2_G_vals[i], g_d[i], 1e-12) );
      TEST_CHECK( gkyl_compare_double(p2_dHdv_vals[i], dhdv_d[i], 1e-12) );
      TEST_CHECK( gkyl_compare_double(p2_d2Gdv2_vals[i], d2gdv2_d[i], 1e-12) );
    }
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

