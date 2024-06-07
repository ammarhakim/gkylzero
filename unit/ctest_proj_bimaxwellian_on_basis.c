#include "gkyl_array.h"
#include "gkyl_util.h"
#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_proj_bimaxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>
#include <gkyl_array_ops.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void eval_bmag(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}
void eval_jacob_tot(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}
void eval_moms_1x2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den = 1.0;
  double upar = 0.5;
  double tpar = 1.3;
  double tperp = 0.6;
  double mass = 1.0;
  fout[0] = den;  // Density.
  fout[1] = den*upar;  // Parallel momentum.
  fout[2] = den*(tpar/mass+upar*upar);  // Parallel kinetic energy.
  fout[3] = den*tperp/mass;  // Perpendicular kinetic energy.
}
void eval_prim_moms_1x2v_gk(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double den = 1.0;
  double upar = 0.5;
  double tpar = 1.3;
  double tperp = 0.6;
  double mass = 1.0;
  fout[0] = den;  // Density.
  fout[1] = upar;  // Parallel drift speed.
  fout[2] = tpar/mass;  // Parallel temperature divided by mass (vtpar^2).
  fout[3] = tperp/mass;  // Perpendicular temperature divided by mass (vtperp^2).
}

void
test_1x2v_gk(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double lower[] = {0.1, -6.0, 0.0}, upper[] = {1.0, 6.0, 6.0};
  int cells[] = {2, 16, 16};
  int vdim = 2, cdim = 1;
  int ndim = cdim+vdim;

  double confLower[] = {lower[0]}, confUpper[] = {upper[0]};
  int confCells[] = {cells[0]};

  double velLower[vdim], velUpper[vdim];
  int velCells[vdim];
  for (int d=0; d<vdim; d++) {
    velLower[d] = lower[cdim+d];
    velUpper[d] = upper[cdim+d];
    velCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[] = { confGhost[0], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // create bmag and jacob_tot arrays
  struct gkyl_array *bmag_ho, *jacob_tot_ho, *bmag, *jacob_tot;
  bmag_ho = mkarr(confBasis.num_basis, confLocal_ext.volume);
  jacob_tot_ho = mkarr(confBasis.num_basis, confLocal_ext.volume);
  if (use_gpu) { // create device copies
    bmag = gkyl_array_cu_dev_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
    jacob_tot = gkyl_array_cu_dev_new(GKYL_DOUBLE, confBasis.num_basis, confLocal_ext.volume);
  } else {
    bmag = bmag_ho;
    jacob_tot = jacob_tot_ho;
  }
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bmag, NULL);
  gkyl_proj_on_basis *proj_jac = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_jacob_tot, NULL);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_proj_on_basis_advance(proj_jac, 0.0, &confLocal, jacob_tot_ho);
  gkyl_array_copy(bmag, bmag_ho);
  gkyl_array_copy(jacob_tot, jacob_tot_ho);

  // create moment arrays
  struct gkyl_array *moms, *moms_ho;
  moms_ho = mkarr(4*confBasis.num_basis, confLocal_ext.volume);
  if (use_gpu) { // create device copies
    moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 4*confBasis.num_basis, confLocal_ext.volume);
  } else {
    moms = moms_ho;
  }

  gkyl_proj_on_basis *proj_moms = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 4, eval_moms_1x2v_gk, NULL);

  gkyl_proj_on_basis_advance(proj_moms, 0.0, &confLocal, moms_ho);
  gkyl_array_copy(moms, moms_ho);

  // create distribution function array
  struct gkyl_array *distf_ho = mkarr(basis.num_basis, local_ext.volume);

  struct gkyl_array *distf;
  if (use_gpu)
    distf = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
  else
    distf = distf_ho;

  // projection updater to compute Maxwellian
  gkyl_proj_bimaxwellian_on_basis *proj_max = gkyl_proj_bimaxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, gvm, use_gpu);

  gkyl_proj_bimaxwellian_on_basis_gyrokinetic_lab_mom(proj_max, &local, &confLocal, moms,
                                                      bmag, jacob_tot, mass, distf);
  gkyl_array_copy(distf_ho, distf);

  // values to compare  at index (1, 9, 9) [remember, lower-left index is (1,1,1)]
  double p1_vals[] = {
    1.2844077246026704e-03, 8.5474047609025796e-20, 2.6350176904782896e-05, -2.2925318340133869e-04, -1.7211242035879901e-20, -9.3936424834558415e-21,
    -4.7032276611997613e-06, 9.8938122762577098e-21, -2.0497552553122845e-05, -1.1378998408032301e-20, 3.6585961643864254e-06, -1.1378998408032301e-20,
  };
  double p2_vals[] = {
    1.2844524987967697e-03, 3.9389697868298960e-20, 2.6351095466632584e-05, -2.3024734563371324e-04, -6.2767367023478419e-21, -6.2661441748961667e-21,
    -4.7236233269944659e-06, -1.3380402930807540e-18, -2.0498267093668484e-05, 1.8487818299271491e-05, -1.2472007912769041e-20, -7.6464635967265530e-20,
    -2.2254527342990302e-20, 2.2846722504428259e-19, 3.6744617592548295e-06, -5.3138683979042951e-21, 3.7928554417563121e-07, -1.5490008912633003e-21,
    -1.5490008912633003e-21, -1.5490008912633003e-21,
  };

  const double *fv = gkyl_array_cfetch(distf_ho, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
  
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
  }

//  // write distribution function to file
//  char fname[1024];
//  sprintf(fname, "ctest_proj_bimaxwellian_on_basis_gyrokinetic_test1_1x2v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, distf_ho, fname);

  // release memory for moment data object
  gkyl_array_release(moms);
  if (use_gpu) {
    gkyl_array_release(moms_ho);
  }
  gkyl_proj_on_basis_release(proj_moms);

  // now perform the maxwellian projection using primitive moments.
  struct gkyl_array *prim_moms, *prim_moms_ho;
  prim_moms_ho = mkarr(4*confBasis.num_basis, confLocal_ext.volume);
  if (use_gpu) {
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 4*confBasis.num_basis, confLocal_ext.volume);
  } else {
    prim_moms = prim_moms_ho;
  }

  gkyl_proj_on_basis *proj_prim_moms = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 4, eval_prim_moms_1x2v_gk, NULL);

  gkyl_proj_on_basis_advance(proj_prim_moms, 0.0, &confLocal, prim_moms_ho);
  gkyl_array_copy(prim_moms, prim_moms_ho);

  gkyl_proj_bimaxwellian_on_basis_gyrokinetic_prim_mom(proj_max, &local, &confLocal, prim_moms,
                                                       bmag, jacob_tot, mass, distf);
  gkyl_array_copy(distf_ho, distf);

  fv = gkyl_array_cfetch(distf_ho, gkyl_range_idx(&local_ext, (int[3]) { 1, 9, 9 }));
  
  if (poly_order == 1) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p1_vals[i], fv[i], 1e-12) );
  }

  if (poly_order == 2) {
    for (int i=0; i<basis.num_basis; ++i)
      TEST_CHECK( gkyl_compare_double(p2_vals[i], fv[i], 1e-12) );
  }

//  // write distribution function to file
//  char fname[1024];
//  sprintf(fname, "ctest_proj_bimaxwellian_on_basis_gyrokinetic_test2_1x2v_p%d.gkyl", poly_order);
//  gkyl_grid_sub_array_write(&grid, &local, distf_ho, fname);

  gkyl_velocity_map_release(gvm);
  gkyl_array_release(prim_moms);
  if (use_gpu) {
    gkyl_array_release(prim_moms_ho);
  }
  gkyl_proj_on_basis_release(proj_prim_moms);

  gkyl_array_release(bmag); gkyl_array_release(jacob_tot);
  if (use_gpu) {
    gkyl_array_release(bmag_ho); gkyl_array_release(jacob_tot_ho);
  }

  gkyl_array_release(distf);
  if (use_gpu) {
    gkyl_array_release(distf_ho);
  }
  gkyl_proj_bimaxwellian_on_basis_release(proj_max);

}

void test_1x2v_p1_gk() { test_1x2v_gk(1, false); }
void test_1x2v_p2_gk() { test_1x2v_gk(2, false); }

#ifdef GKYL_HAVE_CUDA
void test_1x2v_p1_gk_gpu() { test_1x2v_gk(1, true); }
void test_1x2v_p2_gk_gpu() { test_1x2v_gk(2, true); }
#endif

TEST_LIST = {
  { "test_1x2v_p1_gk", test_1x2v_p1_gk },
  { "test_1x2v_p2_gk", test_1x2v_p2_gk },

#ifdef GKYL_HAVE_CUDA
  { "test_1x2v_p1_gk_gpu", test_1x2v_p1_gk_gpu },
  { "test_1x2v_p2_gk_gpu", test_1x2v_p2_gk_gpu },
#endif
  { NULL, NULL },
};
