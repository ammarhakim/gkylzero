// A test for the vlasov LBO collision updater
//
#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_updater_rad_gyrokinetic.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <math.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_app_priv.h>
#include <gkyl_read_radiation.h>
#include <gkyl_dg_bin_ops.h>

void
mapc2p_3x(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bmag_func_3x(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0;
}

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0e19;
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_vthsq(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double *te=ctx;
  fout[0] = te[0]*GKYL_ELEMENTARY_CHARGE/GKYL_ELECTRON_MASS;
}

void
test_1x(int poly_order, bool use_gpu, double te)
{
  double mass = GKYL_ELECTRON_MASS;
  double charge = -1.0*GKYL_ELEMENTARY_CHARGE;
  int vdim = 2;
  int cdim = 1;
  int ndim = cdim + vdim;
  double confLower[cdim], confUpper[cdim], vLower[vdim], vUpper[vdim];
  int confCells[cdim], vCells[vdim];

  double vth = sqrt(te*GKYL_ELEMENTARY_CHARGE/GKYL_ELECTRON_MASS);
  double vth10eV = sqrt(10.0*GKYL_ELEMENTARY_CHARGE/GKYL_ELECTRON_MASS);
  // Phase space and Configuration space extents and resolution
  double lower[] = {-1.0, -4*vth, 0.0};
  int scale = 1; //(int)vth/vth10eV;
  double upper[] = {1.0, 4*vth, 9*vth*vth*GKYL_ELECTRON_MASS};
  int cells[] = {2, 1028*scale, 512*scale*scale};
  confLower[0] = lower[0]; 
  confUpper[0] = upper[0];
  confCells[0] = cells[0];
  vLower[0] = lower[1];
  vLower[1] = lower[2];
  vUpper[0] = upper[1];
  vUpper[1] = upper[2];
  vCells[0] = cells[1];
  vCells[1] = cells[2];

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis, surf_vpar_basis, surf_mu_basis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
    gkyl_cart_modal_serendip(&surf_mu_basis, ndim-1, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
    // constant mu surface
    gkyl_cart_modal_gkhybrid(&surf_mu_basis, cdim, poly_order);
  }
  // constant vpar surface
  gkyl_cart_modal_serendip(&surf_vpar_basis, ndim-1, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  int confGhost[] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], 0 , 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int vGhost[] = {0,0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  // Initialize geometry
  struct gkyl_gyrokinetic_geometry geometry_input = {
      .geometry_id = GKYL_MAPC2P,
      .world = {0.0, 0.0},
      .mapc2p = mapc2p_3x, // mapping of computational to physical space
      .c2p_ctx = 0,
      .bmag_func = bmag_func_3x, // magnetic field magnitude
      .bmag_ctx = 0 
  };
  struct gkyl_rect_grid geo_grid;
  struct gkyl_range geo_local;
  struct gkyl_range geo_local_ext;
  struct gkyl_basis geo_basis;
  bool geo_3d_use_gpu = use_gpu;
  geo_grid = agument_grid(confGrid, geometry_input);
  gkyl_create_grid_ranges(&geo_grid, ghost, &geo_local_ext, &geo_local);
  geo_3d_use_gpu = false;
  gkyl_cart_modal_serendip(&geo_basis, 3, poly_order);
  struct gk_geometry* gk_geom_3d;
  gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geo_grid, &geo_local, &geo_local_ext, &geo_basis, 
      geometry_input.mapc2p, geometry_input.c2p_ctx, geometry_input.bmag_func,  geometry_input.bmag_ctx, geo_3d_use_gpu);
  // deflate geometry if necessary
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &confGrid, &confLocal, &confLocal_ext, 
      &confBasis, use_gpu);
  gkyl_gk_geometry_release(gk_geom_3d);

  // allocate drag coefficients in vparallel and mu for each collision
  // vnu = v_par*nu(v)
  // vsqnu = 2*mu*nu(v)
  // Note that through the spatial variation of B, both these drag coefficients depend on the full phase space
  struct gkyl_array *vnu, *vsqnu, *vnu_surf, *vsqnu_surf, *nvnu, *nvsqnu, *nvnu_surf, *nvsqnu_surf,
    *nvnu_host, *nvsqnu_host, *nvnu_surf_host, *nvsqnu_surf_host;
  vnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vsqnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  vnu_surf = mkarr(use_gpu, surf_vpar_basis.num_basis, local_ext.volume);
  vsqnu_surf = mkarr(use_gpu, surf_mu_basis.num_basis, local_ext.volume);


  struct all_radiation_states *rad_data=gkyl_read_rad_fit_params();
  double a[1], alpha[1], beta[1], gamma[1], v0[1];
  int atomic_z = 5;
  int charge_state = 0;
  int num_ne = 1;
  int status = gkyl_get_fit_params(*rad_data, atomic_z, charge_state, a, alpha, beta, gamma, v0, num_ne);
  if (status == 1) {
    printf("No radiation fits exist for z=%d, charge state=%d\n",atomic_z, charge_state);
  }
  printf("a=%e,alpha=%e,beta=%e,gamma=%e,V0=%e\n",a[0], alpha[0], beta[0], gamma[0], v0[0]);
  /*  a = 0.153650876536253;
  alpha = 8000.006932403581;
  beta = 0.892102642790662;
  gamma = -3.923194017288736;
  v0 = 3.066473173090881;*/
  /*a = 1.1e-5;
  alpha = 3;
  beta = 1;
  gamma = -1.5;
  v0 = 1;*/
  
  struct gkyl_dg_calc_gk_rad_vars *calc_gk_rad_vars = gkyl_dg_calc_gk_rad_vars_new(&grid, &confBasis, &basis, 
		  charge, mass, gk_geom, a[0], alpha[0], beta[0], gamma[0], v0[0], use_gpu);

  gkyl_dg_calc_gk_rad_vars_nu_advance(calc_gk_rad_vars, &confLocal, &local, vnu_surf, vnu, vsqnu_surf, vsqnu);

  nvnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nvsqnu = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  nvnu_surf = mkarr(use_gpu, surf_vpar_basis.num_basis, local_ext.volume);
  nvsqnu_surf = mkarr(use_gpu, surf_mu_basis.num_basis, local_ext.volume);
  
  nvnu_host = nvnu;
  nvsqnu_host = nvsqnu;
  nvnu_surf_host = nvnu_surf;
  nvsqnu_surf_host = nvsqnu_surf;
  if (use_gpu) {
    nvnu_host = mkarr(false, basis.num_basis, local_ext.volume);
    nvsqnu_host = mkarr(false, basis.num_basis, local_ext.volume);
    nvnu_surf_host = mkarr(false, surf_vpar_basis.num_basis, local_ext.volume);
    nvsqnu_surf_host = mkarr(false, surf_mu_basis.num_basis, local_ext.volume);
  }

  if (use_gpu) {
    gkyl_array_copy(nvnu, nvnu_host);
    gkyl_array_copy(nvsqnu, nvsqnu_host);
    gkyl_array_copy(nvnu_surf, nvnu_surf_host);
    gkyl_array_copy(nvsqnu_surf, nvsqnu_surf_host);
  }

  // Initialize distribution function with proj_gkmaxwellian_on_basis
  struct gkyl_array *f = mkarr(use_gpu, basis.num_basis, local_ext.volume);

  struct gkyl_array *nI = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  // Project n, udrift, and vt^2 based on input functions
  struct gkyl_array *m0 = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *udrift = mkarr(false, vdim*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_density, NULL);

  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_upar, 0);
  
  double ctx[1], Lz[1];
  ctx[0]=te;
  gkyl_get_fit_lz(*rad_data, atomic_z, charge_state, log10(1e19), ctx, Lz);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vthsq, ctx);
  printf("Te = %f\n",ctx[0]);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0); 
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);
  //gkyl_array_scale(vtsq, 1.0/mass);

  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms = mkarr(false, 2*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms, 1.0, udrift, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, vtsq, 1*confBasis.num_basis);

  // Initialize Maxwellian projection object
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, use_gpu);

  // If on GPUs, need to copy n, udrift, and vt^2 onto device
  struct gkyl_array *prim_moms_dev, *m0_dev;
  if (use_gpu) {
    prim_moms_dev = mkarr(use_gpu, 2*confBasis.num_basis, confLocal_ext.volume);

    gkyl_array_copy(prim_moms_dev, prim_moms);
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local_ext, &confLocal_ext, m0_dev, prim_moms_dev,
      gk_geom->bmag, gk_geom->bmag, mass, f);
  }
  else {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local_ext, &confLocal_ext, m0, prim_moms,
      gk_geom->bmag, gk_geom->bmag, mass, f);
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local_ext, &confLocal_ext, m0, prim_moms,
      gk_geom->bmag, gk_geom->bmag, GKYL_PROTON_MASS, nI);
  }

  // initialize solver 
  struct gkyl_dg_updater_collisions *slvr;
  struct gkyl_dg_rad_gyrokinetic_auxfields drag_inp = { .nvnu_surf = nvnu_surf, .nvnu = nvnu, 
    .nvsqnu_surf = nvsqnu_surf, .nvsqnu = nvsqnu};
  slvr = gkyl_dg_updater_rad_gyrokinetic_new(&grid, &confBasis, &basis, &local, &drag_inp, use_gpu);

  struct gkyl_array *cflrate, *rhs, *fmax;
  cflrate = mkarr(use_gpu, 1, local_ext.volume);
  rhs = mkarr(use_gpu, basis.num_basis, local_ext.volume);

  gkyl_grid_sub_array_write(&grid, &local, vnu, "ctest_dg_rad_gyrokinetic_vnu.gkyl");
  gkyl_grid_sub_array_write(&confGrid, &confLocal, vtsq, "ctest_dg_rad_gyrokinetic_vthsq.gkyl");
  gkyl_grid_sub_array_write(&confGrid, &confLocal, m0, "ctest_dg_rad_gyrokinetic_M0.gkyl");
  gkyl_grid_sub_array_write(&confGrid, &confLocal, gk_geom->bmag, "ctest_dg_rad_gyrokinetic_bmag.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, vsqnu, "ctest_dg_rad_gyrokinetic_vsqnu.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, vnu_surf, "ctest_dg_rad_gyrokinetic_vnu_surf.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, vsqnu_surf, "ctest_dg_rad_gyrokinetic_vsqnu_surf.gkyl");
  printf("Advancing hyper dg...\n");
  // run hyper_dg_advance
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(cflrate, 0.0);
  gkyl_array_clear(nvnu_surf, 0.0);
  gkyl_array_clear(nvnu, 0.0);
  gkyl_array_clear(nvsqnu_surf, 0.0);
  gkyl_array_clear(nvsqnu, 0.0);

  // Assumed electron and ion density are the same and uniform
  gkyl_dg_calc_gk_rad_vars_nI_nu_advance(calc_gk_rad_vars, 
    &confLocal, &local, vnu_surf, vnu, vsqnu_surf, vsqnu, 
    m0, nvnu_surf, nvnu, nvsqnu_surf, nvsqnu);

  gkyl_grid_sub_array_write(&grid, &local, f, "ctest_dg_rad_gyrokinetic_f.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, nvnu, "ctest_dg_rad_gyrokinetic_nvnu.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, nvsqnu, "ctest_dg_rad_gyrokinetic_nvsqnu.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, nvnu_surf, "ctest_dg_rad_gyrokinetic_nvnu_surf.gkyl");
  gkyl_grid_sub_array_write(&grid, &local, nvsqnu_surf, "ctest_dg_rad_gyrokinetic_nvsqnu_surf.gkyl");

  gkyl_dg_updater_rad_gyrokinetic_advance(slvr, &local, f, cflrate, rhs);

  gkyl_grid_sub_array_write(&grid, &local, rhs, "ctest_dg_rad_gyrokinetic_rhs.gkyl");
  // Take 2nd moment of rhs to find energy loss
  struct gkyl_dg_updater_moment *m2_calc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &confBasis, &basis,
    &confLocal, &vLocal, GKYL_ELECTRON_MASS, gk_geom, "M2", false, use_gpu);
  struct gkyl_array *m2_final = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  gkyl_dg_updater_moment_gyrokinetic_advance(m2_calc, &local, &confLocal, rhs, m2_final);

  double *m00 = gkyl_array_fetch(m0, 0+ghost[0]);
  double *m20 = gkyl_array_fetch(m2_final, 0+ghost[0]);
  //double *m21 = gkyl_array_fetch(rhs, 0+ghost[0]);
  
  double cell_avg_m2 = m20[0]/pow(sqrt(2.0),cdim);
  double cell_avg_m0 = m00[0]/pow(sqrt(2.0),cdim);
  // two factors of density, one for the electrons and one for the ions
  double cell_avg0 = 1.0/2.0*GKYL_ELECTRON_MASS*cell_avg_m2/(cell_avg_m0*cell_avg_m0);

  //double correct = 4.419192427285379e-32;
  double correct = Lz[0];;
  //  for (int i=0; i<30; i++){
printf("cell_avg=%e, correct energy=%e, ratio=%e, percent error = %e, density=%.10e, m2=%e\n", cell_avg0, correct, cell_avg0/correct, (fabs(cell_avg0)-correct)/correct, m00[0], m20[0]);
    //}
  struct gkyl_array *dem = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *emissivity = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  //gkyl_dg_mul_op(confBasis, 0, dem, 0, m0, 0, m0);
  //gkyl_dg_bin_op_mem *mem = gkyl_dg_bin_op_mem_new(emissivity->size, confBasis.num_basis);
  // gkyl_dg_div_op(mem, confBasis, 0, emissivity, 0, m2_final, 0, dem);
  //double *emissivity_a = gkyl_array_fetch(emissivity, 0+ghost[0]);
  //printf("Emissivity from bin_ops=%e",0.5*GKYL_ELECTRON_MASS*emissivity_a[0]);

  
  TEST_CHECK( gkyl_compare( -correct*1e30, cell_avg0*1e30, 1e-12));

  // Release memory
  gkyl_array_release(vnu);
  gkyl_array_release(vsqnu);
  if (use_gpu) {
    gkyl_array_release(nvnu_host);
    gkyl_array_release(nvsqnu_host);      
  }
  gkyl_dg_calc_gk_rad_vars_release(calc_gk_rad_vars);

  gkyl_array_release(m0);
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_array_release(prim_moms);
  if (use_gpu) {
    gkyl_array_release(m0_dev);
    gkyl_array_release(prim_moms_dev);      
  }
  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);
  gkyl_proj_maxwellian_on_basis_release(proj_max);  

  gkyl_array_release(f);
  gkyl_array_release(rhs);
  gkyl_array_release(cflrate);
  gkyl_dg_updater_rad_gyrokinetic_release(slvr);
  //gkyl_dg_updater_moment_gyrokinetic_release(m0_calc);
  gkyl_dg_updater_moment_gyrokinetic_release(m2_calc);
  //gkyl_array_release(m0_final);
  gkyl_array_release(m2_final);
}

void test_1x2v_p1_10eV() { test_1x(1, false, 10.0); }
void test_1x2v_p1_30eV() { test_1x(1, false, 30.0); }
void test_1x2v_p1_50eV() { test_1x(1, false, 50.0); }
void test_1x2v_p1_100eV() { test_1x(1, false, 100.0); }
void test_1x2v_p1_500eV() { test_1x(1, false, 500.0); }
void test_1x2v_p1_1000eV() { test_1x(1, false, 1000.0); }
void test_1x2v_p1_5000eV() { test_1x(1, false, 5000.0); }
void test_1x2v_p1_10000eV() { test_1x(1, false, 10000.0); }
void test_1x2v_p1_50000eV() { test_1x(1, false, 50000.0); }
void test_1x2v_p1_100000eV() { test_1x(1, false, 100000.0); }
void test_1x2v_p2() { test_1x(2, false, 30.0); }


// #ifdef GKYL_HAVE_CUDA

// void test_1x1v_p1_gpu() { test_1x(1, 1, true); }
// void test_1x2v_p1_gpu() { test_1x(1, 2, true); }
// void test_1x1v_p2_gpu() { test_1x(2, 1, true); }
// void test_1x2v_p2_gpu() { test_1x(2, 2, true); }

// #endif

TEST_LIST = {
	     //  { "test_1x2v_p1_10eV", test_1x2v_p1_10eV },
  { "test_1x2v_p1_30eV", test_1x2v_p1_30eV },
  { "test_1x2v_p1_50eV", test_1x2v_p1_50eV },
  { "test_1x2v_p1_100eV", test_1x2v_p1_100eV },
  { "test_1x2v_p1_500eV", test_1x2v_p1_500eV },
  { "test_1x2v_p1_1000eV", test_1x2v_p1_1000eV },
  { "test_1x2v_p1_5000eV", test_1x2v_p1_5000eV },
  { "test_1x2v_p1_10000eV", test_1x2v_p1_10000eV },
  { "test_1x2v_p1_50000eV", test_1x2v_p1_50000eV },
  { "test_1x2v_p1_100000eV", test_1x2v_p1_100000eV },
  // { "test_1x1v_p2", test_1x1v_p2 },
  // { "test_1x2v_p2", test_1x2v_p2 },

// #ifdef GKYL_HAVE_CUDA
//   { "test_1x1v_p1_gpu", test_1x1v_p1_gpu },
//   { "test_1x2v_p1_gpu", test_1x2v_p1_gpu },
//   { "test_1x1v_p2_gpu", test_1x1v_p2_gpu },
//   { "test_1x2v_p2_gpu", test_1x2v_p2_gpu },

// #endif
  { NULL, NULL },
};
