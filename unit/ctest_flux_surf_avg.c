// A test for the integrated moments
//
#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_const.h>

#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_radiation_read.h>
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
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

struct test_ctx {
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM];
  double B0;
  double vt;
  double mass;
};

void eval_bmag_3x(double t, const double *xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];

  struct test_ctx *pars = ctx;
  double B0 = pars->B0;

  fout[0] = B0;
}

void
eval_density(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 3.0e19;
}

void
eval_upar(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_vthsq(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // D ion mass
  fout[0] =150*eV/mi;
}

void
eval_weights(double t, const double * GKYL_RESTRICT xn, double* GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 1.0;
}

void
test_3x_option(bool use_gpu)
{

  /* 0. Prologue
  We set up a simulation 3x2v simulation domain very similar to TCV turb simulations.
  */
  printf("\n0. Prologue\n");
  int poly_order=1;
  double eps0 = GKYL_EPSILON0;
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mi = 2.014*GKYL_PROTON_MASS; // D ion mass
  double qi = eV; // ion charge
  double Ti = 150*eV;
  double B0 = 1.0; // Magnetic field magnitude in Tesla
  double n0 = 3.0e19; // Particle density in 1/m^3
  double vtIon = sqrt(Ti/mi);
  double vpar_max_ion = 6.0*vtIon;
  double mu_max_ion = 12.*mi*vtIon*vtIon/(2.0*B0);

  // Phase space and Configuration space extents and resolution
  double lower[] = {0.0, 0.0, 0.0, -vpar_max_ion, 0.0};
  double upper[] = {1.0, 1.0, 1.0, vpar_max_ion, mu_max_ion};
  int cells[] = {12, 4, 16, 12, 6};
  const int vdim = 2;
  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim - vdim;

  double confLower[cdim], confUpper[cdim], vLower[vdim], vUpper[vdim];
  int confCells[cdim], vCells[vdim];
  for (int i = 0; i < cdim; i++){
    confLower[i] = lower[i]; 
    confUpper[i] = upper[i];
    confCells[i] = cells[i];
  }
  for (int i = 0; i < vdim; i++){
    vLower[i] = lower[cdim+i]; 
    vUpper[i] = upper[cdim+i];
    vCells[i] = cells[cdim+i];
  }

  printf("\t 0.1 set grid\n");
  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  printf("\t 0.2 basis functions\n");
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
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[] = {1, 1, 1};
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int ghost[] = { confGhost[0], confGhost[1], confGhost[2], 0 , 0};
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  int vGhost[] = {0, 0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  printf("\t 0.3 initialize geometry\n");
  // Initialize geometry

    struct test_ctx proj_ctx = {
    .lower = {lower[0], lower[1], lower[2], lower[3], lower[4]},
    .upper = {upper[0], upper[1], upper[2], upper[3], upper[4]},
    .cells = {cells[0], cells[1], cells[2], cells[3], cells[4]},
    .B0 = B0,
    .vt = vtIon,
    .mass = mi,
  };

  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_MAPC2P,
    .c2p_ctx = 0,
    .mapc2p = mapc2p,
    .bmag_ctx = &proj_ctx,
    .bmag_func = eval_bmag_3x,
    .grid = confGrid,
    .local = confLocal,
    .local_ext = confLocal_ext,
    .global = confLocal,
    .global_ext = confLocal_ext,
    .basis = confBasis,
    .geo_grid = confGrid,
    .geo_local = confLocal,
    .geo_local_ext = confLocal_ext,
    .geo_global = confLocal,
    .geo_global_ext = confLocal_ext,
    .geo_basis = confBasis,
  };

  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  struct gk_geometry* gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
  gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
  if (use_gpu) {  // If we are on the gpu, copy from host
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  printf("\t 0.4 project moment functions into DG\n");
  // Project n, udrift, and vt^2 based on input functions
  struct gkyl_array *m0 = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *udrift = mkarr(false, vdim*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_density, NULL);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_upar, 0);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vthsq, 0);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &confLocal, m0); 
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);
  
  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms = mkarr(false, 3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms, 1.0, m0, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, udrift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, vtsq, 2*confBasis.num_basis);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  printf("\t 0.5 initialize f with maxwellian projection\n");
  // Initialize Maxwellian projection object
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&grid,
    &confBasis, &basis, poly_order+1, gvm, use_gpu);

  // Initialize distribution function with proj_gkmaxwellian_on_basis
  struct gkyl_array *f = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  // If on GPUs, need to copy n, udrift, and vt^2 onto device
  struct gkyl_array *prim_moms_dev, *m0_dev;
  if (use_gpu) {
    prim_moms_dev = mkarr(use_gpu, 3*confBasis.num_basis, confLocal_ext.volume);
    m0_dev = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms_dev, prim_moms);
    gkyl_array_copy(m0_dev, m0);
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local_ext, &confLocal_ext, prim_moms_dev,
      gk_geom->bmag, gk_geom->bmag, mi, f);
  }
  else {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local_ext, &confLocal_ext, prim_moms,
      gk_geom->bmag, gk_geom->bmag, mi, f);
  }

  /* 1. Integrating with int_moment 
  
  */
  // Initialize integrated moment calculator
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &confBasis, &basis,
    &confLocal, mi, gvm, gk_geom, "Integrated", true, use_gpu);  
  int num_mom = 4;

  struct gkyl_array *marr = mkarr(use_gpu, num_mom, confLocal_ext.volume);
  struct gkyl_array *marr_host = marr;
  if (use_gpu)
    marr_host = mkarr(false, num_mom, local_ext.volume);  

  double *red_integ_diag_global;

  if (use_gpu)
    red_integ_diag_global = gkyl_cu_malloc(4*sizeof(double));
  else
    red_integ_diag_global = gkyl_malloc(4*sizeof(double));

  // Now calculate the integrated moments
  double avals_global[2+vdim]; // M0, M1, M2par, M2perp
  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc, &local, &confLocal, f, marr);
  gkyl_array_reduce_range(red_integ_diag_global, marr, GKYL_SUM, &confLocal);


  if (use_gpu)
    gkyl_cu_memcpy(avals_global, red_integ_diag_global, sizeof(double[2+vdim]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(avals_global, red_integ_diag_global, sizeof(double[2+vdim]));

  /* 2. Integrating manually the energy moment
    We want now to compute the energy moment M2 and integrating over the whole volume manually
    (this could also be done by calling integrated moments with the dg_updater_moment)
  */
  printf("2. Integrating manually the energy moment\n");
  // Set space to store the M2 moment
  struct gkyl_array *m2arr = mkarr(use_gpu, num_mom, confLocal_ext.volume);
  struct gkyl_array *m2arr_host = m2arr;
  if (use_gpu)
    m2arr_host = mkarr(false, num_mom, local_ext.volume);  

  // Initialize M2 moment calculator
  struct gkyl_dg_updater_moment *m2calc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &confBasis, &basis,
    &confLocal, mi, gvm, gk_geom, "M0", false, use_gpu);    

  // Compute the moment
  gkyl_dg_updater_moment_gyrokinetic_advance(m2calc, &local, &confLocal, f, m2arr);

  // Create a integration updater
  gkyl_array_integrate *local_int_calc = gkyl_array_integrate_new(&grid, &basis, cdim, GKYL_ARRAY_INTEGRATE_OP_NONE, use_gpu);
  // Create weights and initialize by proj on basis
  struct gkyl_array *weights = mkarr(false, confBasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_weights = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_weights, NULL);
  gkyl_proj_on_basis_advance(proj_weights, 0.0, &confLocal, weights); 

  // compute the integral
  double local_int;
  gkyl_array_integrate_advance(local_int_calc, m2arr, 1.0, weights, &confLocal, &local_int);
  // gkyl_array_integrate_advance(local_int_calc, m2arr, 1.0, gk_geom->jacobgeo, &confLocal, &local_int);

  double integral_global;
  if (use_gpu)
    gkyl_cu_memcpy(&integral_global, &local_int, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(&integral_global, &local_int, sizeof(double));

  // Check the integrated moments are correct.
  printf("integrated moment result of M2 is: %g\n", avals_global[0]);
  printf("array integrate of M2 is: %g\n", integral_global);
  // Check of intM1 really just checks the drift velocity is close to zero, this will not be perfect.
  // TEST_CHECK( gkyl_compare( avals_global[0]/3e19, 1.0, 1e-2));
  // TEST_CHECK( gkyl_compare( avals_global[1]/3e19, 0.0, 1e-2));
  // TEST_CHECK( gkyl_compare( avals_global[2]/2.14e+29, 1.0, 1e-2));
  // TEST_CHECK( gkyl_compare( avals_global[3]/4.21e+29, 1.0, 1e-2));


  gkyl_array_release(m0);
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_array_release(prim_moms);

  gkyl_array_release(marr);
  gkyl_array_release(m2arr);

  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);
  gkyl_proj_on_basis_release(proj_weights);
  gkyl_proj_maxwellian_on_basis_release(proj_max);  
  gkyl_velocity_map_release(gvm);

  gkyl_array_release(f);

  if (use_gpu) {
    gkyl_array_release(m0_dev);
    gkyl_array_release(prim_moms_dev);   
  }

  gkyl_dg_updater_moment_gyrokinetic_release(mcalc);
  gkyl_dg_updater_moment_gyrokinetic_release(m2calc);
  gkyl_array_integrate_release(local_int_calc);
  gkyl_gk_geometry_release(gk_geom);
}



void test_3x() { test_3x_option(false); }

#ifdef GKYL_HAVE_CUDA

void test_3x_gpu() { test_3x_option(true); }
#endif

TEST_LIST = {
  { "test_3x", test_3x },
#ifdef GKYL_HAVE_CUDA
  { "test_3x_gpu", test_3x_gpu },
#endif
  { NULL, NULL },
};