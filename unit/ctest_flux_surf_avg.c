/* A test for the flux surface average
description of the test and the function tested here:...
*/ 
#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_average.h>
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

  fout[0] = B0*(1.0+0.5*sin(z))/(x+0.6); // the x dependence will modify the grid shape
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
  double lower[] = {-0.5, -0.5, -0.5, -vpar_max_ion, 0.0};
  double upper[] = { 0.5,  0.5,  0.5, vpar_max_ion, mu_max_ion};
  int cells[] = {4, 6, 8, 12, 6};
  const int vdim = 2;
  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim - vdim;

  double clower[cdim], cupper[cdim], vLower[vdim], vUpper[vdim];
  int ccells[cdim], vCells[vdim];
  for (int i = 0; i < cdim; i++){
    clower[i] = lower[i]; 
    cupper[i] = upper[i];
    ccells[i] = cells[i];
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
  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, cdim, clower, cupper, ccells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  printf("\t 0.2 basis functions\n");
  // basis functions
  struct gkyl_basis basis, cbasis;
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  } else if (poly_order == 1) {
    /* Force hybrid basis (p=2 in vpar). */
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  }
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);

  // Ranges
  int confGhost[] = {1, 1, 1};
  struct gkyl_range crng_xyz, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&cgrid, confGhost, &confLocal_ext, &crng_xyz);

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
    .grid = cgrid,
    .local = crng_xyz,
    .local_ext = confLocal_ext,
    .global = crng_xyz,
    .global_ext = confLocal_ext,
    .basis = cbasis,
    .geo_grid = cgrid,
    .geo_local = crng_xyz,
    .geo_local_ext = confLocal_ext,
    .geo_global = crng_xyz,
    .geo_global_ext = confLocal_ext,
    .geo_basis = cbasis,
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
  struct gkyl_array *m0 = mkarr(false, cbasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *udrift = mkarr(false, vdim*cbasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *vtsq = mkarr(false, cbasis.num_basis, confLocal_ext.volume);
  gkyl_proj_on_basis *proj_m0 = gkyl_proj_on_basis_new(&cgrid, &cbasis,
    poly_order+1, 1, eval_density, NULL);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&cgrid, &cbasis,
    poly_order+1, vdim, eval_upar, 0);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&cgrid, &cbasis,
    poly_order+1, 1, eval_vthsq, 0);
  gkyl_proj_on_basis_advance(proj_m0, 0.0, &crng_xyz, m0); 
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &crng_xyz, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &crng_xyz, vtsq);
  
  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms = mkarr(false, 3*cbasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms, 1.0, m0, 0*cbasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, udrift, 1*cbasis.num_basis);
  gkyl_array_set_offset(prim_moms, 1.0, vtsq, 2*cbasis.num_basis);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  printf("\t 0.5 initialize f with maxwellian projection\n");
  // Initialize Maxwellian projection object
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_new(&grid,
    &cbasis, &basis, poly_order+1, gvm, use_gpu);

  // Initialize distribution function with proj_gkmaxwellian_on_basis
  struct gkyl_array *f = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  // If on GPUs, need to copy n, udrift, and vt^2 onto device
  struct gkyl_array *prim_moms_dev, *m0_dev;
  if (use_gpu) {
    prim_moms_dev = mkarr(use_gpu, 3*cbasis.num_basis, confLocal_ext.volume);
    m0_dev = mkarr(use_gpu, cbasis.num_basis, confLocal_ext.volume);
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
  printf("1. Integration with int_moment\n");
  struct gkyl_dg_updater_moment *mcalc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &cbasis, &basis,
    &crng_xyz, mi, gvm, gk_geom, "Integrated", true, use_gpu);  
  int num_mom = 4;

  struct gkyl_array *marr = mkarr(use_gpu, num_mom, confLocal_ext.volume);
  struct gkyl_array *marr_host = use_gpu? mkarr(false, marr->ncomp, marr->size)
                                        : gkyl_array_acquire(marr);
  double *red_integ_diag_global;

  if (use_gpu)
    red_integ_diag_global = gkyl_cu_malloc(4*sizeof(double));
  else
    red_integ_diag_global = gkyl_malloc(4*sizeof(double));

  // Now calculate the integrated moments
  double avals_global[2+vdim]; // M0, M1, M2par, M2perp
  gkyl_dg_updater_moment_gyrokinetic_advance(mcalc, &local, &crng_xyz, f, marr);
  gkyl_array_reduce_range(red_integ_diag_global, marr, GKYL_SUM, &crng_xyz);


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
  struct gkyl_array *mom_xyz = mkarr(use_gpu, num_mom, confLocal_ext.volume);
  struct gkyl_array *mom_xyz_ho = use_gpu? mkarr(false, mom_xyz->ncomp, mom_xyz->size)
                                        : gkyl_array_acquire(mom_xyz);

  // Initialize M2 moment calculator
  struct gkyl_dg_updater_moment *mom_calc = gkyl_dg_updater_moment_gyrokinetic_new(&grid, &cbasis, &basis,
    &crng_xyz, mi, gvm, gk_geom, "M0", false, use_gpu);    

  // Compute the moment
  gkyl_dg_updater_moment_gyrokinetic_advance(mom_calc, &local, &crng_xyz, f, mom_xyz);

  // Create a integration updater
  gkyl_array_integrate *local_int_calc = gkyl_array_integrate_new(&cgrid, &cbasis, cdim, GKYL_ARRAY_INTEGRATE_OP_NONE, use_gpu);
  // Create weights but useless (jacobian is already in the volume of the cell)
  struct gkyl_array *weights = mkarr(use_gpu, cbasis.num_basis, confLocal_ext.volume);

  // compute the local integral
  double local_int;
  gkyl_array_integrate_advance(local_int_calc, mom_xyz, 1.0, weights, &crng_xyz, &local_int);

  double integral_global;
  if (use_gpu)
    gkyl_cu_memcpy(&integral_global, &local_int, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(&integral_global, &local_int, sizeof(double));

  /* 3. Integrate in 1D the surface integral
    We compute the 1D DG representation of the yz surface integral and integrate it
    over x.
  */
   printf("3. 1D integration of the surface integral DG representation\n");
   printf("\t 3.1 declare 1D structures to perform the yz\n");
  // 1D basis functions
  struct gkyl_basis cbasis_x;
  gkyl_cart_modal_serendip(&cbasis_x, 1, poly_order);
  // 1D grid
  struct gkyl_rect_grid cgrid_x;
  gkyl_rect_grid_init(&cgrid_x, 1, &clower[0], &cupper[0], &ccells[0]);
  // 1D Ranges
  int cghost_x[] = { confGhost[0] };
  struct gkyl_range crng_x, crng_x_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&cgrid_x, cghost_x, &crng_x_ext, &crng_x);
  // 1D gkyl array
  struct gkyl_array *mom_x = mkarr(use_gpu, cbasis_x.num_basis, crng_x_ext.volume);
  struct gkyl_array *mom_x_ho = mom_x;
  if (use_gpu)
    mom_x_ho = mkarr(false, num_mom, crng_x_ext.volume);  
  printf("\t 3.2 create the 2D integral updater and advance it\n");
  // Create an array average updater
  struct gkyl_array_average *avg_yz;
  gkyl_array_average_new(&cgrid, &cbasis, GKYL_ARRAY_AVERAGE_OP_X, use_gpu);
  gkyl_array_average_advance(avg_yz, &crng_xyz, &crng_x, mom_xyz, mom_x);
  gkyl_array_average_release(avg_yz);

  // check output
  struct gkyl_range_iter iter_x;
  gkyl_range_iter_init(&iter_x, &crng_x);
  while (gkyl_range_iter_next(&iter_x)) {
    long lidx = gkyl_range_idx(&crng_x, iter_x.idx);
    const double *m_i = gkyl_array_cfetch(mom_x, lidx);
    printf("\t\tm_x[%ld]=%g\n",lidx,m_i[0]);
  }

  printf("\t 3.3 declare the scalar structures to perform the x average\n");
  // integrate the 1D array
  // scalar gkyl_array (1D)
  struct gkyl_array *mom_int = mkarr(use_gpu, cbasis_x.num_basis, crng_x_ext.volume);
  struct gkyl_array *mom_int_ho = use_gpu? mkarr(false, mom_int->ncomp, mom_int->size)
                                        : gkyl_array_acquire(mom_int);
  // declare a gkyl range for a scalar
  struct gkyl_range rng_0D;
  gkyl_range_lower_skin(&rng_0D, &crng_x, 0, 1);

  printf("\t 3.4 create 1D integral updater and advance it\n");
  // Create an array average updater
  struct gkyl_array_average *avg_x;
  gkyl_array_average_new(&cgrid_x, &cbasis_x, GKYL_ARRAY_AVERAGE_OP, use_gpu);
  gkyl_array_average_advance(avg_x, &crng_x, &crng_x, mom_x, mom_int);
  gkyl_array_average_release(avg_x);

  // check output
  struct gkyl_range_iter iter_;
  gkyl_range_iter_init(&iter_, &rng_0D);
  while (gkyl_range_iter_next(&iter_)) {
    long lidx = gkyl_range_idx(&rng_0D, iter_.idx);
    const double *m_i = gkyl_array_cfetch(mom_int, lidx);
    printf("\t\tm_[%ld]=%g\n",lidx,m_i[0]);
  }
  /* 4. Check the results

  */
  printf("4. Checks\n");
  // Check the integrated moments are correct.
  printf("\tintegrated moment result of M0 is: %g\n", avals_global[0]);
  printf("\tarray integrate of M0 is: %g\n", integral_global);
  // Check of intM1 really just checks the drift velocity is close to zero, this will not be perfect.
  TEST_CHECK( gkyl_compare( avals_global[0], integral_global, 1e-2));

  gkyl_array_release(m0);
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_array_release(prim_moms);

  gkyl_array_release(marr);
  gkyl_array_release(marr_host);
  gkyl_array_release(mom_xyz);
  gkyl_array_release(mom_xyz_ho);
  gkyl_array_release(mom_x);
  gkyl_array_release(mom_x_ho);

  gkyl_proj_on_basis_release(proj_m0);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);
  gkyl_proj_maxwellian_on_basis_release(proj_max);  
  gkyl_velocity_map_release(gvm);

  gkyl_array_release(f);

  if (use_gpu) {
    gkyl_array_release(m0_dev);
    gkyl_array_release(prim_moms_dev);   
  }

  gkyl_dg_updater_moment_gyrokinetic_release(mcalc);
  gkyl_dg_updater_moment_gyrokinetic_release(mom_calc);
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