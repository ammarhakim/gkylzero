#include "gkyl_array.h"
#include "gkyl_util.h"
#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_const.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_gyrokinetic_maxwellian_moments.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_velocity_map.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(long nc, long size)
{
  struct gkyl_array* a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
};

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
  int ndim = parent->ndim;

  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);
  }
}

void eval_den(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double r = sqrt(x*x + y*y);
  double n0 = 2.0e18; 

  double mi = 3.973 * GKYL_PROTON_MASS;
  double qi = GKYL_ELEMENTARY_CHARGE;
  double Te = 6.0 * GKYL_ELEMENTARY_CHARGE;
  double B0 = 0.0398;

  double c_s = sqrt(Te / mi);
  double omega_ci = fabs(qi * B0 / mi);
  double rho_s = c_s / omega_ci;
  double L_perp = 100*rho_s;
  
  //pcg32_random_t rng = gkyl_pcg32_init(0);
  //double perturb = 2.0e-3 * (1.0 - 0.5 * gkyl_pcg32_rand_double(&rng));
  double perturb = 0.0; 

  double n = 0.0;

  if (r < 0.5 * L_perp) {
    n = ((1.0 - (1.0 / 20.0)) * pow(1.0 - (r / (0.5 * L_perp)) * (r / (0.5 * L_perp)), 3.0) + (1.0 / 20.0)) * n0 * (1.0 + perturb);
  }
  else {
    n = (1.0 / 20.0) * n0 * (1.0 + perturb);
  }

  // Set number density.
  fout[0] = n;
}

void eval_udrift(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 0.0;
}

void eval_vtsq_elc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double r = sqrt(x*x + y*y);

  double me = GKYL_ELECTRON_MASS;
  double mi = 3.973 * GKYL_PROTON_MASS;
  double qi = GKYL_ELEMENTARY_CHARGE;
  double Te = 6.0 * GKYL_ELEMENTARY_CHARGE;
  double B0 = 0.0398;

  double c_s = sqrt(Te / mi);
  double omega_ci = fabs(qi * B0 / mi);
  double rho_s = c_s / omega_ci;
  double L_perp = 100*rho_s;
  
  double T = 0.0;

  if (r < 0.5 * L_perp) {
    T = ((1.0 - (1.0 / 5.0)) * pow(1.0 - (r / (0.5 * L_perp)) * (r / (0.5 * L_perp)), 3.0) + (1.0 / 5.0)) * Te;
  }
  else {
    T = (1.0 / 5.0) * Te;
  }

  fout[0] = T / me;
}

void eval_vtsq_ion(double t, const double *xn, double* restrict fout, void *ctx)
{
  double m = 3.973 * GKYL_PROTON_MASS;
  double T = 1.0 * GKYL_ELEMENTARY_CHARGE;

  fout[0] = T / m;
}

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void
bmag_func(double t, const double *xc, double* GKYL_RESTRICT fout, void *ctx)
{
  double x = xc[0], y = xc[1], z = xc[2];
  fout[0] = 0.0398;
}

void
test_elc(int poly_order, bool use_gpu)
{
  double me = GKYL_ELECTRON_MASS;
  double mi = 3.973 * GKYL_PROTON_MASS;
  double qi = GKYL_ELEMENTARY_CHARGE;
  double Te = 6.0 * GKYL_ELEMENTARY_CHARGE;
  double Ti = 1.0 * GKYL_ELEMENTARY_CHARGE;

  double B0 = 0.0398;
 
  double c_s = sqrt(Te / mi);
  double vte = sqrt(Te / me);
  double omega_ci = fabs(qi * B0 / mi);
  double rho_s = c_s / omega_ci;

  double Lx = 100.0* rho_s;
  double Ly = 100.0* rho_s;
  double Lz = 36.0 * 40.0 * rho_s;

  double lower[] = {-0.5*Lx, -0.5*Ly, -0.5*Lz, -4.0*vte, 0.0}; 
  double upper[] = {0.5*Lx, 0.5*Ly, 0.5*Lz, 4.0*vte, (3.0/2.0)*0.5*me*pow(4.0*vte,2)/(2.0*B0)}; 
  int cells[] = {8, 8, 4, 8, 4};
  const int vdim = 2;
  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim-vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_basis *basis_on_dev, *conf_basis_on_dev;
  if (use_gpu) {
#ifdef GKYL_HAVE_CUDA
    basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    conf_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    if (poly_order == 1) 
      gkyl_cart_modal_gkhybrid_cu_dev(basis_on_dev, cdim, vdim);
    else
      gkyl_cart_modal_serendip_cu_dev(basis_on_dev, ndim, poly_order);
    gkyl_cart_modal_serendip_cu_dev(conf_basis_on_dev, cdim, poly_order);
#endif
  }
  else { 
    basis_on_dev = &basis;
    conf_basis_on_dev = &confBasis;
  }

  int confGhost[] = { 1, 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], confGhost[2], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  int vGhost[] = {0, 0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  // Create primitive moment arrays
  struct gkyl_array *den, *udrift, *vtsq;
  den = mkarr(confBasis.num_basis, confLocal_ext.volume);
  udrift = mkarr(confBasis.num_basis, confLocal_ext.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_den = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_den, NULL);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq_elc, NULL);

  gkyl_proj_on_basis_advance(proj_den, 0.0, &confLocal, den);
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);

  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1., den, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., udrift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., vtsq  , 2*confBasis.num_basis);
  struct gkyl_array *prim_moms;
  if (use_gpu) { // copy host array to device
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms, prim_moms_ho);
  } else {
    prim_moms = prim_moms_ho;
  }

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
      .geometry_id = GKYL_MAPC2P,
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = 0,
      .bmag_func = bmag_func, // magnetic field magnitude
      .bmag_ctx =0 ,
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
  struct gk_geometry* gk_geom;
  gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);

  // If we are on the gpu, copy from host
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu, *distf_ho;
  if (use_gpu) { // create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    distf_ho = mkarr(basis.num_basis, local_ext.volume);
  }

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // Projection updater to compute Maxwellian
  struct gkyl_proj_maxwellian_on_basis_inp inp_proj = {
    .grid = &grid,
    .phase_basis = &basis,
    .conf_basis = &confBasis,
    .phase_basis_on_dev = basis_on_dev, 
    .conf_basis_on_dev = conf_basis_on_dev, 
    .phase_range_ext = &local_ext, 
    .conf_range_ext = &confLocal_ext, 
    .vel_map = gvm,
    .use_gpu = use_gpu,
  };
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_inew(&inp_proj);
  
  if (use_gpu) {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
      gk_geom->bmag, gk_geom->jacobtot, me, distf_cu);
    gkyl_array_copy(distf, distf_cu);
    gkyl_grid_sub_array_read(&grid, &local, distf_ho, "ctest_rt_gk_bgk_IC_3x2v_p1-elc_cpu.gkyl"); // load cpu data from file.
  } else {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
      gk_geom->bmag, gk_geom->jacobtot, me, distf);
  }

  // Write distribution function to file
  char fname[1024];
  if (use_gpu) {
    sprintf(fname, "ctest_rt_gk_bgk_IC_3x2v_p%d-elc_gpu.gkyl", poly_order);
  } 
  else {
    sprintf(fname, "ctest_rt_gk_bgk_IC_3x2v_p%d-elc_cpu.gkyl", poly_order);
  }   
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Calculate the moments and copy from device to host.
  struct gkyl_gyrokinetic_maxwellian_moments_inp inp_calc = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .divide_jacobgeo = true,
    .mass = me,
    .use_gpu = use_gpu,
  };
  gkyl_gyrokinetic_maxwellian_moments *calc_moms = gkyl_gyrokinetic_maxwellian_moments_inew( &inp_calc );

  struct gkyl_array *moms;
  moms = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_cu;
  if (use_gpu) {
    moms_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_gyrokinetic_maxwellian_moments_advance(calc_moms, &local, &confLocal, distf_cu, moms_cu);
    gkyl_array_copy(moms, moms_cu);
  }
  else {
    gkyl_gyrokinetic_maxwellian_moments_advance(calc_moms, &local, &confLocal, distf, moms);
  }
  gkyl_gyrokinetic_maxwellian_moments_release(calc_moms);

  // Write moments to file
  char fname_moms[1024];
  if (use_gpu) {
    sprintf(fname_moms, "ctest_rt_gk_bgk_IC_3x2v_p%d-elc_moms_gpu.gkyl", poly_order);
  } 
  else {
    sprintf(fname_moms, "ctest_rt_gk_bgk_IC_3x2v_p%d-elc_moms_cpu.gkyl", poly_order);
  }   
  gkyl_grid_sub_array_write(&confGrid, &confLocal, 0, moms, fname_moms);

  // Compare the GPU results against the CPU results
  if (use_gpu) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&local_ext, iter.idx);
      const double *fv_cpu = gkyl_array_fetch(distf_ho, loc);
      const double *fv_gpu = gkyl_array_fetch(distf, loc);
      for (int i=0; i<basis.num_basis; ++i)
        TEST_CHECK( gkyl_compare_double(fv_cpu[i], fv_gpu[i], 1e-8) );
    }
  }

  gkyl_array_release(den); 
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_proj_on_basis_release(proj_den);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);

  gkyl_array_release(prim_moms_ho);
  if (use_gpu)
    gkyl_array_release(prim_moms);

  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(distf);
  if (use_gpu) {
    gkyl_array_release(distf_cu);
    gkyl_array_release(distf_ho);
  }
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);

  gkyl_array_release(moms);
  if (use_gpu) 
    gkyl_array_release(moms_cu);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) 
    gkyl_cu_free(basis_on_dev);
    gkyl_cu_free(conf_basis_on_dev);
#endif  
}

void
test_ion(int poly_order, bool use_gpu)
{
  double mi = 3.973 * GKYL_PROTON_MASS;
  double qi = GKYL_ELEMENTARY_CHARGE;
  double Te = 6.0 * GKYL_ELEMENTARY_CHARGE;
  double Ti = 1.0 * GKYL_ELEMENTARY_CHARGE;

  double B0 = 0.0398;
 
  double c_s = sqrt(Te / mi);
  double vti = sqrt(Ti / mi);
  double omega_ci = fabs(qi * B0 / mi);
  double rho_s = c_s / omega_ci;

  double Lx = 100.0* rho_s;
  double Ly = 100.0* rho_s;
  double Lz = 36.0 * 40.0 * rho_s;

  double lower[] = {-0.5*Lx, -0.5*Ly, -0.5*Lz, -4.0*vti, 0.0}; 
  double upper[] = {0.5*Lx, 0.5*Ly, 0.5*Lz, 4.0*vti, (3.0/2.0)*0.5*mi*pow(4.0*vti,2)/(2.0*B0)}; 
  int cells[] = {8, 8, 4, 8, 4};
  const int vdim = 2;
  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim-vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_basis *basis_on_dev, *conf_basis_on_dev;
  if (use_gpu) {
#ifdef GKYL_HAVE_CUDA
    basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    conf_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    if (poly_order == 1) 
      gkyl_cart_modal_gkhybrid_cu_dev(basis_on_dev, cdim, vdim);
    else
      gkyl_cart_modal_serendip_cu_dev(basis_on_dev, ndim, poly_order);
    gkyl_cart_modal_serendip_cu_dev(conf_basis_on_dev, cdim, poly_order);
#endif
  }
  else { 
    basis_on_dev = &basis;
    conf_basis_on_dev = &confBasis;
  }

  int confGhost[] = { 1, 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], confGhost[2], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  int vGhost[] = {0, 0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  // Create primitive moment arrays
  struct gkyl_array *den, *udrift, *vtsq;
  den = mkarr(confBasis.num_basis, confLocal_ext.volume);
  udrift = mkarr(confBasis.num_basis, confLocal_ext.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_den = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_den, NULL);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq_elc, NULL);

  gkyl_proj_on_basis_advance(proj_den, 0.0, &confLocal, den);
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);

  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1., den, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., udrift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., vtsq  , 2*confBasis.num_basis);
  struct gkyl_array *prim_moms;
  if (use_gpu) { // copy host array to device
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms, prim_moms_ho);
  } else {
    prim_moms = prim_moms_ho;
  }

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
      .geometry_id = GKYL_MAPC2P,
      .mapc2p = mapc2p, // mapping of computational to physical space
      .c2p_ctx = 0,
      .bmag_func = bmag_func, // magnetic field magnitude
      .bmag_ctx =0 ,
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
  struct gk_geometry* gk_geom;
  gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);

  // If we are on the gpu, copy from host
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu, *distf_ho;
  if (use_gpu)  {// create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    distf_ho = mkarr(basis.num_basis, local_ext.volume);
  }

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // Projection updater to compute Maxwellian
  struct gkyl_proj_maxwellian_on_basis_inp inp_proj = {
    .grid = &grid,
    .phase_basis = &basis,
    .conf_basis = &confBasis,
    .phase_basis_on_dev = basis_on_dev, 
    .conf_basis_on_dev = conf_basis_on_dev, 
    .phase_range_ext = &local_ext, 
    .conf_range_ext = &confLocal_ext, 
    .vel_map = gvm,
    .use_gpu = use_gpu,
  };
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_inew(&inp_proj);
  
  if (use_gpu) {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
      gk_geom->bmag, gk_geom->jacobtot, mi, distf_cu);
    gkyl_array_copy(distf, distf_cu);
    gkyl_grid_sub_array_read(&grid, &local, distf_ho, "ctest_rt_gk_bgk_IC_3x2v_p1-ion_cpu.gkyl"); // load cpu data from file.
  } else {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
      gk_geom->bmag, gk_geom->jacobtot, mi, distf);
  }

  // Write distribution function to file
  char fname[1024];
  if (use_gpu) {
    sprintf(fname, "ctest_rt_gk_bgk_IC_3x2v_p%d-ion_gpu.gkyl", poly_order);
  } 
  else {
    sprintf(fname, "ctest_rt_gk_bgk_IC_3x2v_p%d-ion_cpu.gkyl", poly_order);
  }   
  gkyl_grid_sub_array_write(&grid, &local, 0, distf, fname);

  // Compare the GPU results against the CPU results
  if (use_gpu) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&local_ext, iter.idx);
      const double *fv_cpu = gkyl_array_fetch(distf_ho, loc);
      const double *fv_gpu = gkyl_array_fetch(distf, loc);
      for (int i=0; i<basis.num_basis; ++i)
        TEST_CHECK( gkyl_compare_double(fv_cpu[i], fv_gpu[i], 1e-8) );
    }
  }

  // Calculate the moments and copy from device to host.
  struct gkyl_gyrokinetic_maxwellian_moments_inp inp_calc = {
    .phase_grid = &grid,
    .conf_basis = &confBasis,
    .phase_basis = &basis,
    .conf_range =  &confLocal,
    .conf_range_ext = &confLocal_ext,
    .gk_geom = gk_geom,
    .vel_map = gvm,
    .divide_jacobgeo = true,
    .mass = mi,
    .use_gpu = use_gpu,
  };
  gkyl_gyrokinetic_maxwellian_moments *calc_moms = gkyl_gyrokinetic_maxwellian_moments_inew( &inp_calc );

  struct gkyl_array *moms;
  moms = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *moms_cu;
  if (use_gpu) {
    moms_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_gyrokinetic_maxwellian_moments_advance(calc_moms, &local, &confLocal, distf_cu, moms_cu);
    gkyl_array_copy(moms, moms_cu);
  }
  else {
    gkyl_gyrokinetic_maxwellian_moments_advance(calc_moms, &local, &confLocal, distf, moms);
  }
  gkyl_gyrokinetic_maxwellian_moments_release(calc_moms);

  // Write moments to file
  char fname_moms[1024];
  if (use_gpu) {
    sprintf(fname_moms, "ctest_rt_gk_bgk_IC_3x2v_p%d-ion_moms_gpu.gkyl", poly_order);
  } 
  else {
    sprintf(fname_moms, "ctest_rt_gk_bgk_IC_3x2v_p%d-ion_moms_cpu.gkyl", poly_order);
  }   
  gkyl_grid_sub_array_write(&confGrid, &confLocal, 0, moms, fname_moms);

  gkyl_array_release(den); 
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_proj_on_basis_release(proj_den);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);

  gkyl_array_release(prim_moms_ho);
  if (use_gpu)
    gkyl_array_release(prim_moms);

  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(distf);
  if (use_gpu) {
    gkyl_array_release(distf_cu);
    gkyl_array_release(distf_ho);
  }
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);

  gkyl_array_release(moms);
  if (use_gpu) 
    gkyl_array_release(moms_cu);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) 
    gkyl_cu_free(basis_on_dev);
    gkyl_cu_free(conf_basis_on_dev);
#endif  
}

/*
void
test_3x2v_gk_compare(int poly_order, bool use_gpu)
{
  double mass = 1.0;
  double lower[] = {0.1, 0.1, 0.1, -6.0, 0.0}, upper[] = {1.0, 1.0, 1.0, 6.0, 6.0};
  int cells[] = {2, 2, 2, 16, 16};
  const int vdim = 2;
  const int ndim = sizeof(cells)/sizeof(cells[0]);
  const int cdim = ndim-vdim;

  double confLower[cdim], confUpper[cdim];
  int confCells[cdim];
  for (int d=0; d<cdim; d++) {
    confLower[d] = lower[d];
    confUpper[d] = upper[d];
    confCells[d] = cells[d];
  }
  double vLower[vdim], vUpper[vdim];
  int vCells[vdim];
  for (int d=0; d<vdim; d++) {
    vLower[d] = lower[cdim+d];
    vUpper[d] = upper[cdim+d];
    vCells[d] = cells[cdim+d];
  }

  // grids
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid vGrid;
  gkyl_rect_grid_init(&vGrid, vdim, vLower, vUpper, vCells);

  // basis functions
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1) 
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  struct gkyl_basis *basis_on_dev, *conf_basis_on_dev;
  if (use_gpu) {
#ifdef GKYL_HAVE_CUDA
    basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    conf_basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    if (poly_order == 1) 
      gkyl_cart_modal_gkhybrid_cu_dev(basis_on_dev, cdim, vdim);
    else
      gkyl_cart_modal_serendip_cu_dev(basis_on_dev, ndim, poly_order);
    gkyl_cart_modal_serendip_cu_dev(conf_basis_on_dev, cdim, poly_order);
#endif
  }
  else { 
    basis_on_dev = &basis;
    conf_basis_on_dev = &confBasis;
  }

  int confGhost[] = { 1, 1, 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);
  struct skin_ghost_ranges confSkin_ghost; // conf-space skin/ghost
  skin_ghost_ranges_init(&confSkin_ghost, &confLocal_ext, confGhost);

  int ghost[] = { confGhost[0], confGhost[1], confGhost[2], 0, 0 };
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  struct skin_ghost_ranges skin_ghost; // phase-space skin/ghost
  skin_ghost_ranges_init(&skin_ghost, &local_ext, ghost);

  int vGhost[] = {0, 0};
  struct gkyl_range vLocal, vLocal_ext;
  gkyl_create_grid_ranges(&vGrid, vGhost, &vLocal_ext, &vLocal);

  // Create primitive moment arrays
  struct gkyl_array *den, *udrift, *vtsq;
  den = mkarr(confBasis.num_basis, confLocal_ext.volume);
  udrift = mkarr(confBasis.num_basis, confLocal_ext.volume);
  vtsq = mkarr(confBasis.num_basis, confLocal_ext.volume);

  gkyl_proj_on_basis *proj_den = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_den_3x, NULL);
  gkyl_proj_on_basis *proj_udrift = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, vdim, eval_udrift_2v_gk, NULL);
  gkyl_proj_on_basis *proj_vtsq = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_vtsq, NULL);

  gkyl_proj_on_basis_advance(proj_den, 0.0, &confLocal, den);
  gkyl_proj_on_basis_advance(proj_udrift, 0.0, &confLocal, udrift);
  gkyl_proj_on_basis_advance(proj_vtsq, 0.0, &confLocal, vtsq);

  // proj_maxwellian expects the primitive moments as a single array.
  struct gkyl_array *prim_moms_ho = mkarr(3*confBasis.num_basis, confLocal_ext.volume);
  gkyl_array_set_offset(prim_moms_ho, 1., den, 0*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., udrift, 1*confBasis.num_basis);
  gkyl_array_set_offset(prim_moms_ho, 1., vtsq  , 2*confBasis.num_basis);
  struct gkyl_array *prim_moms;
  if (use_gpu) { // copy host array to device
    prim_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*confBasis.num_basis, confLocal_ext.volume);
    gkyl_array_copy(prim_moms, prim_moms_ho);
  } else {
    prim_moms = prim_moms_ho;
  }

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
      .geometry_id = GKYL_MAPC2P,
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
      .geo_grid = confGrid,
      .geo_local = confLocal,
      .geo_local_ext = confLocal_ext,
      .geo_global = confLocal,
      .geo_global_ext = confLocal_ext,
      .geo_basis = confBasis,
  };
  struct gk_geometry* gk_geom;
  gk_geom = gkyl_gk_geometry_mapc2p_new(&geometry_input);

  // If we are on the gpu, copy from host
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }

  // Create distribution function array
  struct gkyl_array *distf;
  distf = mkarr(basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_cu, *distf_ho;
  if (use_gpu)  { // create device copy.
    distf_cu  = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, local_ext.volume);
    distf_ho = mkarr(basis.num_basis, local_ext.volume); // array to load cpu data.
  }

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, vGrid,
    local, local_ext, vLocal, vLocal_ext, use_gpu);

  // Projection updater to compute Maxwellian
  struct gkyl_proj_maxwellian_on_basis_inp inp_proj = {
    .grid = &grid,
    .phase_basis = &basis,
    .conf_basis = &confBasis,
    .phase_basis_on_dev = basis_on_dev, 
    .conf_basis_on_dev = conf_basis_on_dev, 
    .phase_range_ext = &local_ext, 
    .conf_range_ext = &confLocal_ext, 
    .vel_map = gvm,
    .use_gpu = use_gpu,
  };
  gkyl_proj_maxwellian_on_basis *proj_max = gkyl_proj_maxwellian_on_basis_inew(&inp_proj);
  
  if (use_gpu) {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
      gk_geom->bmag, gk_geom->jacobtot, mass, distf_cu);
    gkyl_array_copy(distf, distf_cu);
    gkyl_grid_sub_array_read(&grid, &local, distf_ho, "ctest_proj_gkmaxwellian_on_basis_prim_mom_test_3x2v_p1_cpu.gkyl"); // load cpu data from file.
  } else {
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(proj_max, &local, &confLocal, prim_moms,
      gk_geom->bmag, gk_geom->jacobtot, mass, distf);
  }

  // Compare the GPU results against the CPU results
  if (use_gpu) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&local_ext, iter.idx);
      const double *fv_cpu = gkyl_array_fetch(distf_ho, loc);
      const double *fv_gpu = gkyl_array_fetch(distf, loc);
      for (int i=0; i<basis.num_basis; ++i)
        TEST_CHECK( gkyl_compare_double(fv_cpu[i], fv_gpu[i], 1e-8) );
    }
  }

  gkyl_array_release(den); 
  gkyl_array_release(udrift); 
  gkyl_array_release(vtsq);
  gkyl_proj_on_basis_release(proj_den);
  gkyl_proj_on_basis_release(proj_udrift);
  gkyl_proj_on_basis_release(proj_vtsq);

  gkyl_array_release(prim_moms_ho);
  if (use_gpu) {
    gkyl_array_release(prim_moms);
  }

  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(distf);
  if (use_gpu)
    gkyl_array_release(distf_cu);
    gkyl_array_release(distf_ho);
  gkyl_proj_maxwellian_on_basis_release(proj_max);
  gkyl_velocity_map_release(gvm);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) 
    gkyl_cu_free(basis_on_dev);
    gkyl_cu_free(conf_basis_on_dev);
#endif  
}
*/
void test_IC_elc() { test_elc(1, false); }
void test_IC_ion() { test_ion(1, false); }
// No p2 geometry in 3D

#ifdef GKYL_HAVE_CUDA
void test_IC_elc_gpu() { test_elc(1, true); }
void test_IC_ion_gpu() { test_ion(1, true); }
#endif

TEST_LIST = {
  { "test_IC_elc", test_IC_elc },
  { "test_IC_ion", test_IC_ion },

#ifdef GKYL_HAVE_CUDA
  { "test_IC_elc_gpu", test_IC_elc_gpu },
  { "test_IC_ion_gpu", test_IC_ion_gpu },
#endif
  { NULL, NULL },
};
