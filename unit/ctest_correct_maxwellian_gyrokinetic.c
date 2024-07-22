#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gyrokinetic_maxwellian_correct.h>
#include <gkyl_gyrokinetic_maxwellian_moments.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <math.h>


static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

void
mapc2p_3x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}


void
bmag_func_3x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  fout[0] = 0.5;
}

void eval_n(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0e19*(1.0 + 0.5*cos(x));
}
void eval_vtsq(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double vtsq = (10.0*1.602e-19/9.1e-31)*exp(-x*x/(M_PI));
  fout[0] = vtsq;
}
void eval_upar(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];

  double vtsq[1];
  eval_vtsq(t, xn, vtsq, ctx);
  double vt = sqrt(vtsq[0]); 

  fout[0] = vt;
}

void eval_n_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], z = xn[1];
  fout[0] = 1.0e19*(1.0 + 0.5*cos(x)*sin(z));
}
void eval_vtsq_2x(double t, const double *xn, double* restrict fout, void *ctx)
{
  // in 2D, just initialize a constant temperature
  double x = xn[0], z = xn[1];
  double vtsq = (10.0*1.602e-19/9.1e-31);//*exp(-x*x/(M_PI))*exp(-z*z/(M_PI));
  fout[0] = vtsq;
}
void eval_upar_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double vtsq[1];
  eval_vtsq_2x(t, xn, vtsq, ctx);
  double vt = sqrt(vtsq[0]); 

  fout[0] = vt;
}

void test_1x1v(int poly_order, bool use_gpu)
{
  double mass = 9.1e-31;
  double err_max = 1.0e-10, iter_max = 50;
  double vt = sqrt(10.0*1.602e-19/9.1e-31); // reference temperature
  double lower[] = {-M_PI, -4.0*vt}, upper[] = {M_PI, 4.0*vt};
  int cells[] = {4, 16};
  const int vdim = 1;

  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim - vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
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
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  } 
  else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1, 1 }; // 3 elements because it's used by geo.
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[ndim] = { 0 };
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p_3x, // mapping of computational to physical space
      .c2p_ctx = 0,
      .bmag_func = bmag_func_3x, // magnetic field magnitude
      .bmag_ctx = 0,
      .grid = confGrid,
      .local = confLocal,
      .local_ext = confLocal_ext,
      .global = confLocal,
      .global_ext = confLocal_ext,
      .basis = confBasis,
  };
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, confGhost, &geometry_input.geo_local_ext, &geometry_input.geo_local);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  struct gk_geometry* gk_geom_3d;
  gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  // deflate geometry if necessary
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // Create correct moment arrays
  struct gkyl_array *m0_in_ho, *m1_in_ho, *m2_in_ho;
  m0_in_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  m1_in_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  m2_in_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_n, NULL);
  gkyl_proj_on_basis *proj_m1 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_upar, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_vtsq, NULL);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_in_ho);
  gkyl_proj_on_basis_advance(proj_m1, 0.0, &confLocal, m1_in_ho);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_in_ho);
  struct gkyl_array *m0_in, *m1_in, *m2_in;
  if (use_gpu) {
    m0_in = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    m1_in = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    m2_in = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
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
  struct gkyl_array *moms_in = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(moms_in, 1., m0_in, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m1_in, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m2_in, 2*confBasis.num_basis);
  // (2) create distribution function array
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, gvm, use_gpu);
  struct gkyl_array *fM = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_maxwellian, &local, &confLocal, 
    moms_in, gk_geom->bmag, gk_geom->jacobtot, mass, fM);
  // (3) copy from device to host
  struct gkyl_array *fM_ho;
  if (use_gpu) {
    fM_ho = mkarr(basis.num_basis, local_ext.volume, false);
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }
  // Optionally write out the inital distribution function and desired initial moments  
  // char fname_fM_ic[1024];
  // sprintf(fname_fM_ic, "ctest_correct_maxwellian_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&grid, &local, 0, fM_ho, fname_fM_ic);
  // char moms_ic[1024];
  // sprintf(moms_ic, "ctest_correct_maxwellian_moms_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&confGrid, &confLocal, 0, moms_in, moms_ic);
  // Create a Maxwellian with corrected moments
  struct gkyl_gyrokinetic_maxwellian_correct_inp inp = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .mass = mass, 
    .max_iter = iter_max, 
    .eps = err_max, 
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .divide_jacobgeo = true, 
    .use_gpu = use_gpu
  };
  gkyl_gyrokinetic_maxwellian_correct *corr_max = gkyl_gyrokinetic_maxwellian_correct_inew(&inp);
  // First correct the density
  gkyl_gyrokinetic_maxwellian_correct_density_moment(corr_max, 
    fM, moms_in, &local, &confLocal);
  // Now correct all the moments
  struct gkyl_gyrokinetic_maxwellian_correct_status status_corr;
  status_corr = gkyl_gyrokinetic_maxwellian_correct_all_moments(corr_max, 
    fM, moms_in, &local, &confLocal);
  gkyl_gyrokinetic_maxwellian_correct_release(corr_max);
  if (use_gpu) {
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }

  // Compute the moments of our corrected distribution function
  struct gkyl_gyrokinetic_maxwellian_moments_inp inp_mom = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .divide_jacobgeo = true, 
    .mass = mass, 
    .use_gpu = use_gpu,
  };
  gkyl_gyrokinetic_maxwellian_moments *max_moms = gkyl_gyrokinetic_maxwellian_moments_inew( &inp_mom );
  // (2) calculate the moments and copy from host to device
  struct gkyl_array *moms_corr = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_corr_ho;
  gkyl_gyrokinetic_maxwellian_moments_advance(max_moms, &local, &confLocal, fM, moms_corr);
  gkyl_gyrokinetic_maxwellian_moments_release(max_moms);
  if (use_gpu) { 
    moms_corr_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume, false);
    gkyl_array_copy(moms_corr_ho, moms_corr); 
  } 
  else {
    moms_corr_ho = moms_corr;
  }
  // Optionally write out the corrected distribution function and moments  
  // char fname_fM_corr[1024];
  // sprintf(fname_fM_corr, "ctest_correct_maxwellian_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&grid, &local, 0, fM_ho, fname_fM_corr);
  // char moms_corr_ic[1024];
  // sprintf(moms_corr_ic, "ctest_correct_maxwellian_moms_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&confGrid, &confLocal, 0, moms_corr_ho, moms_corr_ic);

  // (4) compare the correct moments with the input moments
  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&confLocal, idx);
    const double *m0in = gkyl_array_cfetch(m0_in_ho, linidx);
    const double *m1in = gkyl_array_cfetch(m1_in_ho, linidx);
    const double *m2in = gkyl_array_cfetch(m2_in_ho, linidx);
    const double *momsCorr = gkyl_array_cfetch(moms_corr_ho, linidx);
    // Check the cell averages
    TEST_CHECK( gkyl_compare(m0in[0], momsCorr[0*confBasis.num_basis], 1.0e-10) );
    TEST_CHECK( gkyl_compare(m1in[0], momsCorr[1*confBasis.num_basis], 1.0e-10) );
    TEST_CHECK( gkyl_compare(m2in[0], momsCorr[2*confBasis.num_basis], 1.0e-10) );
    TEST_MSG("Expected cell average: %.13e, \t%.13e, \t%.13e, \tin cell (%d)", m0in[0], m1in[0], m2in[0], idx[0]);
    TEST_MSG("Produced cell average: %.13e, \t%.13e, \t%.13e", momsCorr[0*confBasis.num_basis], momsCorr[1*confBasis.num_basis], momsCorr[2*confBasis.num_basis]);
  }
 
  // Release memory for moment data object
  gkyl_gk_geometry_release(gk_geom);  
  gkyl_velocity_map_release(gvm);
  gkyl_array_release(m0_in);
  gkyl_array_release(m1_in);
  gkyl_array_release(m2_in);
  gkyl_array_release(moms_in);
  gkyl_array_release(moms_corr);
  gkyl_array_release(fM);
  if (use_gpu) {
    gkyl_array_release(m0_in_ho);
    gkyl_array_release(m1_in_ho);
    gkyl_array_release(m2_in_ho);
    gkyl_array_release(moms_corr_ho);
    gkyl_array_release(fM_ho);
  }
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_proj_maxwellian_on_basis_release(proj_maxwellian);
}

void test_1x2v(int poly_order, bool use_gpu)
{
  double mass = 9.1e-31;
  double err_max = 1.0e-10, iter_max = 50;
  double vt = sqrt(10.0*1.602e-19/9.1e-31); // reference temperature
  double lower[] = {-M_PI, -4.0*vt, 0.0}, upper[] = {M_PI, 4.0*vt, 4.0*vt*vt*mass};
  int cells[] = {4, 16, 16};
  const int vdim = 2;

  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim - vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
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
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1, 1 }; // 3 elements because it's used by geo.
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[ndim] = { 0 };
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p_3x, // mapping of computational to physical space
      .c2p_ctx = 0,
      .bmag_func = bmag_func_3x, // magnetic field magnitude
      .bmag_ctx =0 ,
      .grid = confGrid,
      .local = confLocal,
      .local_ext = confLocal_ext,
      .global = confLocal,
      .global_ext = confLocal_ext,
      .basis = confBasis,
  };
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, confGhost, &geometry_input.geo_local_ext, &geometry_input.geo_local);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  struct gk_geometry* gk_geom_3d;
  gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  // deflate geometry if necessary
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // Create correct moment arrays
  struct gkyl_array *m0_in_ho, *m1_in_ho, *m2_in_ho;
  m0_in_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  m1_in_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  m2_in_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_n, NULL);
  gkyl_proj_on_basis *proj_m1 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_upar, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_vtsq, NULL);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_in_ho);
  gkyl_proj_on_basis_advance(proj_m1, 0.0, &confLocal, m1_in_ho);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_in_ho);
  struct gkyl_array *m0_in, *m1_in, *m2_in;
  if (use_gpu) {
    m0_in = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    m1_in = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    m2_in = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
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
  struct gkyl_array *moms_in = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(moms_in, 1., m0_in, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m1_in, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m2_in, 2*confBasis.num_basis);
  // (2) create distribution function array
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, gvm, use_gpu);
  struct gkyl_array *fM = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_maxwellian, &local, &confLocal, 
    moms_in, gk_geom->bmag, gk_geom->jacobtot, mass, fM);
  // (3) copy from device to host
  struct gkyl_array *fM_ho;
  if (use_gpu) {
    fM_ho = mkarr(basis.num_basis, local_ext.volume, false);
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }
  // Optionally write out the inital distribution function and desired initial moments  
  // char fname_fM_ic[1024];
  // sprintf(fname_fM_ic, "ctest_correct_maxwellian_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&grid, &local, 0, fM_ho, fname_fM_ic);
  // char moms_ic[1024];
  // sprintf(moms_ic, "ctest_correct_maxwellian_moms_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&confGrid, &confLocal, 0, moms_in, moms_ic);
  // Create a Maxwellian with corrected moments
  struct gkyl_gyrokinetic_maxwellian_correct_inp inp = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .divide_jacobgeo = true, 
    .mass = mass, 
    .max_iter = iter_max, 
    .eps = err_max, 
    .use_gpu = use_gpu
  };
  gkyl_gyrokinetic_maxwellian_correct *corr_max = gkyl_gyrokinetic_maxwellian_correct_inew(&inp);
  struct gkyl_gyrokinetic_maxwellian_correct_status status_corr;
  status_corr = gkyl_gyrokinetic_maxwellian_correct_all_moments(corr_max, 
    fM, moms_in, &local, &confLocal);
  gkyl_gyrokinetic_maxwellian_correct_release(corr_max);
  if (use_gpu) {
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }
  
  // Compute the moments of our corrected distribution function
  struct gkyl_gyrokinetic_maxwellian_moments_inp inp_mom = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .divide_jacobgeo = true, 
    .mass = mass, 
    .use_gpu = use_gpu,
  };
  gkyl_gyrokinetic_maxwellian_moments *max_moms = gkyl_gyrokinetic_maxwellian_moments_inew( &inp_mom );
  // (2) calculate the moments and copy from host to device
  struct gkyl_array *moms_corr = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_corr_ho;
  gkyl_gyrokinetic_maxwellian_moments_advance(max_moms, &local, &confLocal, fM, moms_corr);
  gkyl_gyrokinetic_maxwellian_moments_release(max_moms);
  if (use_gpu) { 
    moms_corr_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume, false);
    gkyl_array_copy(moms_corr_ho, moms_corr); 
  } 
  else {
    moms_corr_ho = moms_corr;
  }
  // Optionally write out the corrected distribution function and moments  
  // char fname_fM_corr[1024];
  // sprintf(fname_fM_corr, "ctest_correct_maxwellian_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&grid, &local, 0, fM_ho, fname_fM_corr);
  // char moms_corr_ic[1024];
  // sprintf(moms_corr_ic, "ctest_correct_maxwellian_moms_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&confGrid, &confLocal, 0, moms_corr_ho, moms_corr_ic);

  // (4) compare the correct moments with the input moments
  for (int k=0; k<cells[0]; k++) {
    int idx[] = {k+1};
    long linidx = gkyl_range_idx(&confLocal, idx);
    const double *m0in = gkyl_array_cfetch(m0_in_ho, linidx);
    const double *m1in = gkyl_array_cfetch(m1_in_ho, linidx);
    const double *m2in = gkyl_array_cfetch(m2_in_ho, linidx);
    const double *momsCorr = gkyl_array_cfetch(moms_corr_ho, linidx);
    // Check the cell averages
    TEST_CHECK( gkyl_compare(m0in[0], momsCorr[0*confBasis.num_basis], 1.0e-10) );
    TEST_CHECK( gkyl_compare(m1in[0], momsCorr[1*confBasis.num_basis], 1.0e-10) );
    TEST_CHECK( gkyl_compare(m2in[0], momsCorr[2*confBasis.num_basis], 1.0e-10) );
    TEST_MSG("Expected cell average: %.13e, \t%.13e, \t%.13e, \tin cell (%d)", m0in[0], m1in[0], m2in[0], idx[0]);
    TEST_MSG("Produced cell average: %.13e, \t%.13e, \t%.13e", momsCorr[0*confBasis.num_basis], momsCorr[1*confBasis.num_basis], momsCorr[2*confBasis.num_basis]);
  }

  // Release memory for moment data object
  gkyl_gk_geometry_release(gk_geom);  
  gkyl_velocity_map_release(gvm);
  gkyl_array_release(m0_in);
  gkyl_array_release(m1_in);
  gkyl_array_release(m2_in);
  gkyl_array_release(moms_in);
  gkyl_array_release(moms_corr);
  gkyl_array_release(fM);
  if (use_gpu) {
    gkyl_array_release(m0_in_ho);
    gkyl_array_release(m1_in_ho);
    gkyl_array_release(m2_in_ho);
    gkyl_array_release(moms_corr_ho);
    gkyl_array_release(fM_ho);
  }
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_proj_maxwellian_on_basis_release(proj_maxwellian);
}

void test_2x2v(int poly_order, bool use_gpu)
{
  double mass = 9.1e-31;
  double err_max = 1.0e-10, iter_max = 50;
  double vt = sqrt(10.0*1.602e-19/9.1e-31); // reference temperature
  double lower[] = {-M_PI, -M_PI, -4.0*vt, 0.0}, upper[] = {M_PI, M_PI, 4.0*vt, 4.0*vt*vt*mass};
  int cells[] = {4, 4, 16, 16};
  const int vdim = 2;

  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim - vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
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
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1, 1, 1 }; // 3 elements because it's used by geo.
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[] = { 0, 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[ndim] = { 0 };
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0},
      .mapc2p = mapc2p_3x, // mapping of computational to physical space
      .c2p_ctx = 0,
      .bmag_func = bmag_func_3x, // magnetic field magnitude
      .bmag_ctx =0 ,
      .grid = confGrid,
      .local = confLocal,
      .local_ext = confLocal_ext,
      .global = confLocal,
      .global_ext = confLocal_ext,
      .basis = confBasis,
  };
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, confGhost, &geometry_input.geo_local_ext, &geometry_input.geo_local);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  struct gk_geometry* gk_geom_3d;
  gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  // deflate geometry if necessary
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);

  // Initialize velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // Create correct moment arrays
  struct gkyl_array *m0_in_ho, *m1_in_ho, *m2_in_ho;
  m0_in_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  m1_in_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  m2_in_ho = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_n_2x, NULL);
  gkyl_proj_on_basis *proj_m1 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_upar_2x, NULL);
  gkyl_proj_on_basis *proj_m2 = gkyl_proj_on_basis_new(&confGrid, &confBasis, poly_order+1, 1, eval_vtsq_2x, NULL);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0_in_ho);
  gkyl_proj_on_basis_advance(proj_m1, 0.0, &confLocal, m1_in_ho);
  gkyl_proj_on_basis_advance(proj_m2, 0.0, &confLocal, m2_in_ho);
  struct gkyl_array *m0_in, *m1_in, *m2_in;
  if (use_gpu) {
    m0_in = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    m1_in = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    m2_in = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
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
  struct gkyl_array *moms_in = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(moms_in, 1., m0_in, 0*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m1_in, 1*confBasis.num_basis);
  gkyl_array_set_offset(moms_in, 1., m2_in, 2*confBasis.num_basis);
  // (2) create distribution function array
  gkyl_proj_maxwellian_on_basis *proj_maxwellian = gkyl_proj_maxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, gvm, use_gpu);
  struct gkyl_array *fM = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_maxwellian, &local, &confLocal, 
    moms_in, gk_geom->bmag, gk_geom->jacobtot, mass, fM);
  // (3) copy from device to host
  struct gkyl_array *fM_ho;
  if (use_gpu) {
    fM_ho = mkarr(basis.num_basis, local_ext.volume, false);
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }
  // Optionally write out the inital distribution function and desired initial moments  
  // char fname_fM_ic[1024];
  // sprintf(fname_fM_ic, "ctest_correct_maxwellian_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&grid, &local, fM_ho, fname_fM_ic);
  // char moms_ic[1024];
  // sprintf(moms_ic, "ctest_correct_maxwellian_moms_%dx%dv_p%d.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&confGrid, &confLocal, moms_in, moms_ic);
  // Create a Maxwellian with corrected moments
  struct gkyl_gyrokinetic_maxwellian_correct_inp inp = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .divide_jacobgeo = true, 
    .mass = mass, 
    .max_iter = iter_max, 
    .eps = err_max, 
    .use_gpu = use_gpu
  };
  gkyl_gyrokinetic_maxwellian_correct *corr_max = gkyl_gyrokinetic_maxwellian_correct_inew(&inp);
  struct gkyl_gyrokinetic_maxwellian_correct_status status_corr;
  status_corr = gkyl_gyrokinetic_maxwellian_correct_all_moments(corr_max, 
    fM, moms_in, &local, &confLocal);
  gkyl_gyrokinetic_maxwellian_correct_release(corr_max);
  if (use_gpu) {
    gkyl_array_copy(fM_ho, fM);
  } else {
    fM_ho = fM;
  }
  
  // Compute the moments of our corrected distribution function
  struct gkyl_gyrokinetic_maxwellian_moments_inp inp_mom = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .divide_jacobgeo = true, 
    .mass = mass, 
    .use_gpu = use_gpu,
  };
  gkyl_gyrokinetic_maxwellian_moments *max_moms = gkyl_gyrokinetic_maxwellian_moments_inew( &inp_mom );
  // (2) calculate the moments and copy from host to device
  struct gkyl_array *moms_corr = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_corr_ho;
  gkyl_gyrokinetic_maxwellian_moments_advance(max_moms, &local, &confLocal, fM, moms_corr);
  gkyl_gyrokinetic_maxwellian_moments_release(max_moms);
  if (use_gpu) { 
    moms_corr_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume, false);
    gkyl_array_copy(moms_corr_ho, moms_corr); 
  } 
  else {
    moms_corr_ho = moms_corr;
  }
  // Optionally write out the corrected distribution function and moments  
  // char fname_fM_corr[1024];
  // sprintf(fname_fM_corr, "ctest_correct_maxwellian_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&grid, &local, fM_ho, fname_fM_corr);
  // char moms_corr_ic[1024];
  // sprintf(moms_corr_ic, "ctest_correct_maxwellian_moms_%dx%dv_p%d_corr.gkyl", cdim, vdim, poly_order);
  // gkyl_grid_sub_array_write(&confGrid, &confLocal,moms_corr_ho, moms_corr_ic);

  // (4) compare the correct moments with the input moments
  for (int i=0; i<cells[0]; i++) {
    for (int j=0; j<cells[1]; j++) {
      int idx[] = {i+confGhost[0], j+confGhost[1]};
      long linidx = gkyl_range_idx(&confLocal, idx);
      const double *m0in = gkyl_array_cfetch(m0_in_ho, linidx);
      const double *m1in = gkyl_array_cfetch(m1_in_ho, linidx);
      const double *m2in = gkyl_array_cfetch(m2_in_ho, linidx);
      const double *momsCorr = gkyl_array_cfetch(moms_corr_ho, linidx);
      // Check the cell averages
      TEST_CHECK( gkyl_compare(m0in[0], momsCorr[0*confBasis.num_basis], 1.0e-10) );
      TEST_CHECK( gkyl_compare(m1in[0], momsCorr[1*confBasis.num_basis], 1.0e-10) );
      TEST_CHECK( gkyl_compare(m2in[0], momsCorr[2*confBasis.num_basis], 1.0e-10) );
      TEST_MSG("Expected cell average: %.13e, \t%.13e, \t%.13e, \tin cell (%d)", m0in[0], m1in[0], m2in[0], idx[0]);
      TEST_MSG("Produced cell average: %.13e, \t%.13e, \t%.13e", momsCorr[0*confBasis.num_basis], momsCorr[1*confBasis.num_basis], momsCorr[2*confBasis.num_basis]);
    }
  }

  // Release memory for moment data object
  gkyl_gk_geometry_release(gk_geom);  
  gkyl_velocity_map_release(gvm);
  gkyl_array_release(m0_in);
  gkyl_array_release(m1_in);
  gkyl_array_release(m2_in);
  gkyl_array_release(moms_in);
  gkyl_array_release(moms_corr);
  gkyl_array_release(fM);
  if (use_gpu) {
    gkyl_array_release(m0_in_ho);
    gkyl_array_release(m1_in_ho);
    gkyl_array_release(m2_in_ho);
    gkyl_array_release(moms_corr_ho);
    gkyl_array_release(fM_ho);
  }
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_m1);
  gkyl_proj_on_basis_release(proj_m2);
  gkyl_proj_maxwellian_on_basis_release(proj_maxwellian);
}

// Run the test
void test_1x1v_p1() {test_1x1v(1, false);}
void test_1x2v_p1() {test_1x2v(1, false);}
void test_2x2v_p1() {test_2x2v(1, false);}

#ifdef GKYL_HAVE_CUDA
void test_1x1v_p1_gpu() {test_1x1v(1, true);}
void test_1x2v_p1_gpu() {test_1x2v(1, true);}
void test_2x2v_p1_gpu() {test_2x2v(1, true);}
#endif

TEST_LIST = {
  {"test_1x1v_p1", test_1x1v_p1},
  {"test_1x2v_p1", test_1x2v_p1},
  {"test_2x2v_p1", test_2x2v_p1},
#ifdef GKYL_HAVE_CUDA
  {"test_1x1v_p1_gpu", test_1x1v_p1_gpu},
  {"test_1x2v_p1_gpu", test_1x2v_p1_gpu},
  {"test_2x2v_p1_gpu", test_2x2v_p1_gpu},
#endif
  {NULL, NULL},
};
