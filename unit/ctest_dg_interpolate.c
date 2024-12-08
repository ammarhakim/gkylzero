#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_velocity_map.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_dg_updater_moment_gyrokinetic.h>
#include <gkyl_array_integrate.h>
#include <gkyl_util.h>
#include <gkyl_array_rio.h>
#include <acutest.h>

// Allocate array (filled with zeros).
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

struct test_ctx {
  double n0; // Density.
  double upar; // Parallel flow speed.
  double T0; // Temperature.
  double mass; // Particle mass.
  double B0; // Magnetic field.
  int vdim; // Number of velocity space dimensions.
  double vpar_max; // Maximum vpar of the grid.
};

void
mapc2p(double t, const double *xc, double* GKYL_RESTRICT xp, void *ctx)
{
  xp[0] = xc[0]; xp[1] = xc[1]; xp[2] = xc[2];
}

void eval_bmag_1x(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];

  struct test_ctx *tctx = ctx;
  double B0 = tctx->B0;

  fout[0] = B0;
}

void eval_distf_1x1v(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], vpar = xn[1];

  struct test_ctx *tctx = ctx;
  double n0 = tctx->n0;
  double upar = tctx->upar;
  double T0 = tctx->T0;
  double B0 = tctx->B0;
  double mass = tctx->mass;
  int vdim = tctx->vdim;

  double vtsq = T0/mass;

  fout[0] = (n0/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vpar-upar,2))/(2.0*vtsq));
//  double a = 2.0,  b = 5.0;
//  double c[] = {a/12.0 - 0.5, 0.0};
//  double d[] = {0.0, b/12.0 - 0.5};
//
//  fout[0] = (pow(x,2)/2.0+c[0]*x+c[1])*(pow(vpar,2)/2.0+d[0]*vpar+d[1]);
}

static void calc_moms(struct gkyl_rect_grid *grid, struct gkyl_basis *confBasis, struct gkyl_basis *basis,
  struct gkyl_range *confLocal, struct gkyl_range *local, double mass, struct gkyl_velocity_map *gvm,
  struct gk_geometry *gk_geom, bool use_gpu, struct gkyl_array *distf, struct gkyl_array *moms)
{
  struct gkyl_dg_updater_moment* mom_op = gkyl_dg_updater_moment_gyrokinetic_new(grid, confBasis,
    basis, confLocal, mass, gvm, gk_geom, "ThreeMoments", false, use_gpu);
  gkyl_dg_updater_moment_gyrokinetic_advance(mom_op, local, confLocal, distf, moms);
  gkyl_dg_updater_moment_gyrokinetic_release(mom_op);
}

static void calc_int_moms(int num_mom, struct gkyl_rect_grid *confGrid, struct gkyl_basis *confBasis,
  struct gkyl_range *confLocal, bool use_gpu, struct gkyl_array *moms, double *int_moms)
{
  // Compute the volume integral of the moments.
  double *integrated_moms = use_gpu? gkyl_cu_malloc(num_mom*sizeof(double)) : gkyl_malloc(num_mom*sizeof(double));

  struct gkyl_array_integrate* integ_op = gkyl_array_integrate_new(confGrid, confBasis,
    num_mom, GKYL_ARRAY_INTEGRATE_OP_NONE, use_gpu);

  gkyl_array_integrate_advance(integ_op, moms, 1.0, moms, confLocal, integrated_moms); 

  if (use_gpu)
    gkyl_cu_memcpy(int_moms, integrated_moms, sizeof(double[num_mom]), GKYL_CU_MEMCPY_D2H);
  else
    memcpy(int_moms, integrated_moms, sizeof(double[num_mom]));

  if (use_gpu)
    gkyl_cu_free(integrated_moms);
  else
    gkyl_free(integrated_moms);
  gkyl_array_integrate_release(integ_op);
}

static struct gk_geometry* init_gk_geo(int poly_order, struct gkyl_rect_grid confGrid, struct gkyl_basis confBasis,
  struct gkyl_range confLocal, struct gkyl_range confLocal_ext, void *bmag_ctx, bool use_gpu)
{
  // Initialize GK geometry.
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_MAPC2P,
    .world = {0.0},  .mapc2p = mapc2p,  .c2p_ctx = 0,
    .bmag_func = eval_bmag_1x,  .bmag_ctx = bmag_ctx,
    .basis = confBasis,  .grid = confGrid,
    .local = confLocal,  .local_ext = confLocal_ext,
    .global = confLocal, .global_ext = confLocal_ext,
  };
  int geo_ghost[3] = {1, 1, 1};
  geometry_input.geo_grid = gkyl_gk_geometry_augment_grid(confGrid, geometry_input);
  gkyl_cart_modal_serendip(&geometry_input.geo_basis, 3, poly_order);
  gkyl_create_grid_ranges(&geometry_input.geo_grid, geo_ghost, &geometry_input.geo_global_ext, &geometry_input.geo_global);
  memcpy(&geometry_input.geo_local, &geometry_input.geo_global, sizeof(struct gkyl_range));
  memcpy(&geometry_input.geo_local_ext, &geometry_input.geo_global_ext, sizeof(struct gkyl_range));
  // Deflate geometry.
  struct gk_geometry* gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_input);
  struct gk_geometry *gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_input);
  gkyl_gk_geometry_release(gk_geom_3d);
  // If we are on the gpu, copy from host.
  if (use_gpu) {
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(gk_geom, &geometry_input, use_gpu);
    gkyl_gk_geometry_release(gk_geom);
    gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);
  }
  return gk_geom;
}

void
test_1x1v_gk(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  const int cdim = 1,  vdim = 1;
  double x_min = 0.0;
  double x_max = 1.0;
  double vpar_min = -6.0;
  double vpar_max =  6.0;
  double lower[] = {x_min, vpar_min}, upper[] = {x_max, vpar_max};
  double mass = 1.0;

  const int ndim = cdim+vdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar = 0, // Parallel flow speed.
    .T0 = 1.0, // Temperature.
    .B0 = 1.0, // Magnetic field.
    .vdim = vdim, // Number of velocity space dimensions.
    .vpar_max = vpar_max, // Maximum vpar of the grid.
    .mass = mass, // Particle mass.
  };

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

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // Basis functions.
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create bmag arrays.
  struct gkyl_array *bmag = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size)
                                      : gkyl_array_acquire(bmag);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bmag_1x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_array_copy(bmag, bmag_ho);

  // Create distribution function arrays.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_1x1v, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize geometry.
  struct gk_geometry *gk_geom = init_gk_geo(poly_order, confGrid, confBasis, confLocal, confLocal_ext, &proj_ctx, use_gpu);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  int num_mom = 3;
  // Calculate the moments.
  struct gkyl_array *moms = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_ext.volume);
  calc_moms(&grid, &confBasis, &basis, &confLocal, &local, mass, gvm, gk_geom, use_gpu, distf, moms);

  // Calculate the integrated moments.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &confGrid, &confBasis, &confLocal, use_gpu, moms, int_moms);

  // Write donor distribution function to file.
  char fname0[1024];
  sprintf(fname0, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);

  // Write target moments to file.
  char fname0m[1024];
  sprintf(fname0m, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_do_mom.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, moms, fname0m);

  printf("\n  Donor int_moms  = %g %g %g\n", int_moms[0],int_moms[1],int_moms[2]);

  // Target grids.
  int confCells_tar[cdim], velCells_tar[vdim];
  for (int d=0; d<cdim; d++)
    confCells_tar[d] = cells_tar[d];
  for (int d=0; d<vdim; d++)
    velCells_tar[d] = cells_tar[cdim+d];

  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);
  struct gkyl_rect_grid confGrid_tar;
  gkyl_rect_grid_init(&confGrid_tar, cdim, confLower, confUpper, confCells_tar);
  struct gkyl_rect_grid velGrid_tar;
  gkyl_rect_grid_init(&velGrid_tar, vdim, velLower, velUpper, velCells_tar);

  // Target ranges.
  struct gkyl_range confLocal_tar, confLocal_tar_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid_tar, confGhost, &confLocal_tar_ext, &confLocal_tar);

  struct gkyl_range velLocal_tar, velLocal_tar_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid_tar, velGhost, &velLocal_tar_ext, &velLocal_tar);

  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target geometry.
  struct gk_geometry *gk_geom_tar = init_gk_geo(poly_order, confGrid_tar, confBasis,
    confLocal_tar, confLocal_tar_ext, &proj_ctx, use_gpu);

  // Target velocity space mapping.
  struct gkyl_velocity_map *gvm_tar = gkyl_velocity_map_new(c2p_in, grid_tar, velGrid_tar,
    local_tar, local_tar_ext, velLocal_tar, velLocal_tar_ext, use_gpu);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(cdim, &basis,
    &grid, &grid_tar, use_gpu);

  gkyl_dg_interpolate_advance(interp, &local, &local_tar, distf, distf_tar);

  // Calculate the moments.
  struct gkyl_array *moms_tar = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_tar_ext.volume);
  calc_moms(&grid_tar, &confBasis, &basis, &confLocal_tar, &local_tar, mass, gvm_tar, gk_geom_tar, use_gpu, distf_tar, moms_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &confGrid_tar, &confBasis, &confLocal_tar, use_gpu, moms_tar, int_moms_tar);

  // Write target distribution function to file.
  char fname1[1024];
  sprintf(fname1, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);

  // Write target moments to file.
  char fname1m[1024];
  sprintf(fname1m, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_tar_mom.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
  gkyl_grid_sub_array_write(&confGrid_tar, &confLocal_tar, NULL, moms_tar, fname1m);

  printf("\n  Target int_moms = %g %g %g\n", int_moms_tar[0],int_moms_tar[1],int_moms_tar[2]);

  // Check the results against accepted results.
  gkyl_array_copy(distf_tar_ho, distf_tar);

//  if (cells[0]     == 3 && cells[1]     == 3 &&
//      cells_tar[0] == 3 && cells_tar[1] == 4) {
//    // From 3x3 to 3x4.
//    const double sol[72] =  {
//       5.3155006858710552e-03,  2.3016793138989707e-03,
//      -8.9097263763831017e-04, -3.8580246913580163e-04,
//      -1.9515639104739080e-18, -1.3010426069826053e-18,
//       8.5733882030178117e-04,  3.7123859901596213e-04,
//      -2.2109321008061795e-03, -9.5736168267032413e-04,
//      -3.4081206790120182e-04, -1.4757595435937272e-04,
//      -8.4019204389574834e-03, -3.6381382703564403e-03,
//      -3.1349037250236887e-03, -1.3574531321444889e-03,
//      -3.4081206790120399e-04, -1.4757595435937066e-04,
//      -2.2462277091906721e-02, -9.7264512942182284e-03,
//      -4.4548631881915352e-03, -1.9290123456790012e-03,
//       3.3827107781547738e-17,  1.3877787807814457e-17,
//      
//       5.3155006858710560e-03, -2.3016793138989672e-03,
//      -8.9097263763830996e-04,  3.8580246913580066e-04,
//      -2.6020852139652106e-18,  1.9515639104739080e-18,
//       8.5733882030178084e-04, -3.7123859901596159e-04,
//      -2.2109321008061790e-03,  9.5736168267032337e-04,
//      -3.4081206790119987e-04,  1.4757595435937933e-04,
//      -8.4019204389574765e-03,  3.6381382703564308e-03,
//      -3.1349037250236861e-03,  1.3574531321444885e-03,
//      -3.4081206790120225e-04,  1.4757595435936459e-04,
//      -2.2462277091906714e-02,  9.7264512942182128e-03,
//      -4.4548631881915352e-03,  1.9290123456790020e-03,
//       3.5561831257524545e-17,  0.0000000000000000e+00,
//      
//      -1.0631001371742102e-02, -6.9050379416969168e-03,
//       1.7819452752766182e-03,  1.1574074074074058e-03,
//       1.7347234759768071e-18,  3.9031278209478160e-18,
//      -1.7146776406035597e-03, -1.1137157970478892e-03,
//       4.4218642016123546e-03,  2.8720850480109753e-03,
//       6.8162413580241882e-04,  4.4272786307813410e-04,
//       1.6803840877914932e-02,  1.0914414811069327e-02,
//       6.2698074500473714e-03,  4.0723593964334722e-03,
//       6.8162413580239584e-04,  4.4272786307807642e-04,
//       4.4924554183813387e-02,  2.9179353882654718e-02,
//       8.9097263763830599e-03,  5.7870370370370072e-03,
//      -1.1449174941446927e-16, -5.7245874707234634e-17,
//    };
//    long i = 0;
//    struct gkyl_range_iter iter;
//    gkyl_range_iter_init(&iter, &local_tar);
//    while (gkyl_range_iter_next(&iter)) {
//      long linidx = gkyl_range_idx(&local_tar, iter.idx);
//      const double *f_c = gkyl_array_cfetch(distf_tar_ho, linidx);
//      for (int m=0; m<basis.num_basis; m++) {
//        TEST_CHECK( gkyl_compare(sol[i], f_c[m], 1e-12) );
//        TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[i], iter.idx[0], iter.idx[1]);
//        TEST_MSG("Produced: %.13e", f_c[m]);
//        i += 1;
//      }
//    }
//  }
//  else if (cells[0]     == 6 && cells[1]      == 12 &&
//           cells_tar[0] == 12 && cells_tar[1] == 12) {
//    // From 6x12 to 12x12.
//    // First 2 cells along vpar only.
//    const double sol[144] =  {
//       2.0924425582990402e-03,  9.8842276988000142e-04,
//      -2.5522653682347228e-05, -1.2056327160493826e-05,
//       0.0000000000000000e+00, -2.1684043449710089e-19,
//       1.9156164266117964e-03,  9.0489408510140948e-04,
//      -7.6567961047041644e-05, -3.6168981481481215e-05,
//      -9.7578195523695399e-19, -6.5052130349130266e-19,
//      
//       5.5164394718792854e-03,  9.8842276988000142e-04,
//      -6.7286996071642870e-05, -1.2056327160493826e-05,
//      -6.5052130349130266e-19, -5.4210108624275222e-19,
//       5.0502614883401898e-03,  9.0489408510140948e-04,
//      -2.0186098821492793e-04, -3.6168981481481215e-05,
//      -3.3610267347050637e-18, -9.7578195523695399e-19,
//      
//       7.7991040809327833e-03,  3.2947425662666952e-04,
//      -9.5129890997840399e-05, -4.0187757201647585e-06,
//      -4.3368086899420177e-19,  4.3368086899420177e-19,
//       7.1400248628257848e-03,  3.0163136170047259e-04,
//      -2.8538967299351870e-04, -1.2056327160493864e-05,
//      -4.3368086899420177e-18, -4.3368086899420177e-19,
//      
//       8.9404363854595319e-03,  3.2947425662666660e-04,
//      -1.0905133846093873e-04, -4.0187757201646043e-06,
//      -4.3368086899420177e-19,  0.0000000000000000e+00,
//       8.1849065500685849e-03,  3.0163136170046972e-04,
//      -3.2715401538281413e-04, -1.2056327160493699e-05,
//      -5.2041704279304213e-18, -4.3368086899420177e-19,
//      
//       8.9404363854595319e-03, -3.2947425662666725e-04,
//      -1.0905133846093873e-04,  4.0187757201644849e-06,
//      -8.6736173798840355e-19,  0.0000000000000000e+00,
//       8.1849065500685849e-03, -3.0163136170046988e-04,
//      -3.2715401538281456e-04,  1.2056327160493606e-05,
//      -4.7704895589362195e-18, -8.6736173798840355e-19,
//      
//       7.7991040809327825e-03, -3.2947425662667007e-04,
//      -9.5129890997840399e-05,  4.0187757201642968e-06,
//      -4.3368086899420177e-19,  0.0000000000000000e+00,
//       7.1400248628257865e-03, -3.0163136170047275e-04,
//      -2.8538967299351914e-04,  1.2056327160493569e-05,
//      -4.3368086899420177e-18, -6.5052130349130266e-19,
//      
//       5.5164394718792862e-03, -9.8842276988000142e-04,
//      -6.7286996071642654e-05,  1.2056327160493798e-05,
//      -1.9515639104739080e-18,  1.0842021724855044e-19,
//       5.0502614883401907e-03, -9.0489408510140948e-04,
//      -2.0186098821492793e-04,  3.6168981481481120e-05,
//      -3.6862873864507151e-18, -1.1926223897340549e-18,
//      
//       2.0924425582990449e-03, -9.8842276987999838e-04,
//      -2.5522653682347445e-05,  1.2056327160494414e-05,
//       1.7347234759768071e-18,  2.9273458657108620e-18,
//       1.9156164266118014e-03, -9.0489408510140644e-04,
//      -7.6567961047041861e-05,  3.6168981481481635e-05,
//      -5.4210108624275222e-19,  2.4936649967166602e-18,
//      
//      -2.4728866598079449e-03, -1.6473712831333374e-03,
//       3.0163136170047178e-05,  2.0093878600823011e-05,
//       1.9515639104739080e-18, -8.2399365108898337e-18,
//      -2.2639103223593862e-03, -1.5081568085023513e-03,
//       9.0489408510140097e-05,  6.0281635802468706e-05,
//       1.9515639104739080e-18,  4.3368086899420177e-19,
//      
//      -8.1795481824416833e-03, -1.6473712831333337e-03,
//       9.9770373485539264e-05,  2.0093878600825159e-05,
//       8.6736173798840355e-19,  2.1684043449710089e-19,
//      -7.4883187585733707e-03, -1.5081568085023476e-03,
//       2.9931112045661676e-04,  6.0281635802470603e-05,
//       7.1557343384043293e-18, -3.4694469519536142e-18,
//      
//      -1.5027542009602192e-02, -2.3063197963866674e-03,
//       1.8329905826413098e-04,  2.8131430041155208e-05,
//       3.4694469519536142e-18,  8.6736173798840355e-19,
//      -1.3757608882030174e-02, -2.1114195319032864e-03,
//       5.4989717479239056e-04,  8.4394290123459062e-05,
//       1.0408340855860843e-17,  1.7347234759768071e-18,
//      
//      -2.3016868141289417e-02, -2.3063197963866674e-03,
//       2.8074919050581973e-04,  2.8131430041149062e-05,
//       1.7347234759768071e-18,  8.6736173798840355e-19,
//      -2.1071780692729743e-02, -2.1114195319032864e-03,
//       8.4224757151745767e-04,  8.4394290123453438e-05,
//       1.3010426069826053e-17,  4.3368086899420177e-18,
//    };
//    long i = 0;
//    struct gkyl_range_iter iter;
//    gkyl_range_iter_init(&iter, &local_tar);
//    while (gkyl_range_iter_next(&iter)) {
//      if (iter.idx[1] < 3) {
//        long linidx = gkyl_range_idx(&local_tar, iter.idx);
//        const double *f_c = gkyl_array_cfetch(distf_tar_ho, linidx);
//        for (int m=0; m<basis.num_basis; m++) {
//          TEST_CHECK( gkyl_compare(sol[i], f_c[m], 1e-12) );
//          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[i], iter.idx[0], iter.idx[1]);
//          TEST_MSG("Produced: %.13e", f_c[m]);
//          i += 1;
//        }
//      }
//    }
//  }

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(moms_tar);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_velocity_map_release(gvm_tar);
  gkyl_gk_geometry_release(gk_geom_tar);
  gkyl_array_release(moms);
  gkyl_velocity_map_release(gvm);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(bmag);
  gkyl_array_release(distf);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_bmag);
  gkyl_proj_on_basis_release(proj_distf);
}

void
test_1x2v_gk(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  const int cdim = 1,  vdim = 2;
  double x_min = 0.0;
  double x_max = 1.0;
  double vpar_min = -6.0;
  double vpar_max =  6.0;
  double mu_max =  0.5;
  double lower[] = {x_min, vpar_min, 0.0}, upper[] = {x_max, vpar_max, mu_max};
  double mass = 1.0;

  const int ndim = cdim+vdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar = 0, // Parallel flow speed.
    .T0 = 1.0, // Temperature.
    .B0 = 1.0, // Magnetic field.
    .vdim = vdim, // Number of velocity space dimensions.
    .vpar_max = vpar_max, // Maximum vpar of the grid.
    .mass = mass, // Particle mass.
  };

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

  // Grids.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid confGrid;
  gkyl_rect_grid_init(&confGrid, cdim, confLower, confUpper, confCells);
  struct gkyl_rect_grid velGrid;
  gkyl_rect_grid_init(&velGrid, vdim, velLower, velUpper, velCells);

  // Basis functions.
  struct gkyl_basis basis, confBasis;
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&confBasis, cdim, poly_order);

  // Ranges
  int confGhost[GKYL_MAX_CDIM] = { 1 };
  struct gkyl_range confLocal, confLocal_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid, confGhost, &confLocal_ext, &confLocal);

  int velGhost[3] = { 0 };
  struct gkyl_range velLocal, velLocal_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid, velGhost, &velLocal_ext, &velLocal);

  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; d++) ghost[d] = confGhost[d];
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Create bmag arrays.
  struct gkyl_array *bmag = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  struct gkyl_array *bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size)
                                      : gkyl_array_acquire(bmag);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bmag_1x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_array_copy(bmag, bmag_ho);

  // Create distribution function arrays.
  struct gkyl_array *distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                                       : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_1x1v, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Initialize geometry.
  struct gk_geometry *gk_geom = init_gk_geo(poly_order, confGrid, confBasis, confLocal, confLocal_ext, &proj_ctx, use_gpu);

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  int num_mom = 3;
  // Calculate the moments.
  struct gkyl_array *moms = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_ext.volume);
  calc_moms(&grid, &confBasis, &basis, &confLocal, &local, mass, gvm, gk_geom, use_gpu, distf, moms);

  // Calculate the integrated moments.
  double int_moms[num_mom];
  for (int i=0; i<num_mom; i++) int_moms[i] = 0.0;
  calc_int_moms(num_mom, &confGrid, &confBasis, &confLocal, use_gpu, moms, int_moms);

  // Write donor distribution function to file.
  char fname0[1024];
  sprintf(fname0, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);

  // Write target moments to file.
  char fname0m[1024];
  sprintf(fname0m, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_do_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
  gkyl_grid_sub_array_write(&confGrid, &confLocal, NULL, moms, fname0m);

  printf("\n  Donor int_moms  = %g %g %g\n", int_moms[0],int_moms[1],int_moms[2]);

  // Target grids.
  int confCells_tar[cdim], velCells_tar[vdim];
  for (int d=0; d<cdim; d++)
    confCells_tar[d] = cells_tar[d];
  for (int d=0; d<vdim; d++)
    velCells_tar[d] = cells_tar[cdim+d];

  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);
  struct gkyl_rect_grid confGrid_tar;
  gkyl_rect_grid_init(&confGrid_tar, cdim, confLower, confUpper, confCells_tar);
  struct gkyl_rect_grid velGrid_tar;
  gkyl_rect_grid_init(&velGrid_tar, vdim, velLower, velUpper, velCells_tar);

  // Target ranges.
  struct gkyl_range confLocal_tar, confLocal_tar_ext; // local, local-ext conf-space ranges
  gkyl_create_grid_ranges(&confGrid_tar, confGhost, &confLocal_tar_ext, &confLocal_tar);

  struct gkyl_range velLocal_tar, velLocal_tar_ext; // local, local-ext vel-space ranges
  gkyl_create_grid_ranges(&velGrid_tar, velGhost, &velLocal_tar_ext, &velLocal_tar);

  struct gkyl_range local_tar, local_tar_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target geometry.
  struct gk_geometry *gk_geom_tar = init_gk_geo(poly_order, confGrid_tar, confBasis,
    confLocal_tar, confLocal_tar_ext, &proj_ctx, use_gpu);

  // Target velocity space mapping.
  struct gkyl_velocity_map *gvm_tar = gkyl_velocity_map_new(c2p_in, grid_tar, velGrid_tar,
    local_tar, local_tar_ext, velLocal_tar, velLocal_tar_ext, use_gpu);

  // Target field.
  struct gkyl_array *distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  struct gkyl_array *distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                                           : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(cdim, &basis,
    &grid, &grid_tar, use_gpu);

  gkyl_dg_interpolate_advance(interp, &local, &local_tar, distf, distf_tar);

  // Calculate the moments.
  struct gkyl_array *moms_tar = mkarr(use_gpu, num_mom*confBasis.num_basis, confLocal_tar_ext.volume);
  calc_moms(&grid_tar, &confBasis, &basis, &confLocal_tar, &local_tar, mass, gvm_tar, gk_geom_tar, use_gpu, distf_tar, moms_tar);

  // Calculate the integrated moments of the target.
  double int_moms_tar[num_mom];
  for (int i=0; i<num_mom; i++) int_moms_tar[i] = 0.0;
  calc_int_moms(num_mom, &confGrid_tar, &confBasis, &confLocal_tar, use_gpu, moms_tar, int_moms_tar);

  // Write target distribution function to file.
  char fname1[1024];
  sprintf(fname1, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);

  // Write target moments to file.
  char fname1m[1024];
  sprintf(fname1m, "ctest_dg_interp_1x2v_p%d_N%dx%dx%d-N%dx%dx%d_tar_mom.gkyl", poly_order, cells[0], cells[1], cells[2], cells_tar[0], cells_tar[1], cells_tar[2]);
  gkyl_grid_sub_array_write(&confGrid_tar, &confLocal_tar, NULL, moms_tar, fname1m);

  printf("\n  Target int_moms = %g %g %g\n", int_moms_tar[0],int_moms_tar[1],int_moms_tar[2]);

  // Check the results against accepted results.
  gkyl_array_copy(distf_tar_ho, distf_tar);

//  if (cells[0]     == 3 && cells[1]     == 3 &&
//      cells_tar[0] == 3 && cells_tar[1] == 4) {
//    // From 3x3 to 3x4.
//    const double sol[72] =  {
//       5.3155006858710552e-03,  2.3016793138989707e-03,
//      -8.9097263763831017e-04, -3.8580246913580163e-04,
//      -1.9515639104739080e-18, -1.3010426069826053e-18,
//       8.5733882030178117e-04,  3.7123859901596213e-04,
//      -2.2109321008061795e-03, -9.5736168267032413e-04,
//      -3.4081206790120182e-04, -1.4757595435937272e-04,
//      -8.4019204389574834e-03, -3.6381382703564403e-03,
//      -3.1349037250236887e-03, -1.3574531321444889e-03,
//      -3.4081206790120399e-04, -1.4757595435937066e-04,
//      -2.2462277091906721e-02, -9.7264512942182284e-03,
//      -4.4548631881915352e-03, -1.9290123456790012e-03,
//       3.3827107781547738e-17,  1.3877787807814457e-17,
//      
//       5.3155006858710560e-03, -2.3016793138989672e-03,
//      -8.9097263763830996e-04,  3.8580246913580066e-04,
//      -2.6020852139652106e-18,  1.9515639104739080e-18,
//       8.5733882030178084e-04, -3.7123859901596159e-04,
//      -2.2109321008061790e-03,  9.5736168267032337e-04,
//      -3.4081206790119987e-04,  1.4757595435937933e-04,
//      -8.4019204389574765e-03,  3.6381382703564308e-03,
//      -3.1349037250236861e-03,  1.3574531321444885e-03,
//      -3.4081206790120225e-04,  1.4757595435936459e-04,
//      -2.2462277091906714e-02,  9.7264512942182128e-03,
//      -4.4548631881915352e-03,  1.9290123456790020e-03,
//       3.5561831257524545e-17,  0.0000000000000000e+00,
//      
//      -1.0631001371742102e-02, -6.9050379416969168e-03,
//       1.7819452752766182e-03,  1.1574074074074058e-03,
//       1.7347234759768071e-18,  3.9031278209478160e-18,
//      -1.7146776406035597e-03, -1.1137157970478892e-03,
//       4.4218642016123546e-03,  2.8720850480109753e-03,
//       6.8162413580241882e-04,  4.4272786307813410e-04,
//       1.6803840877914932e-02,  1.0914414811069327e-02,
//       6.2698074500473714e-03,  4.0723593964334722e-03,
//       6.8162413580239584e-04,  4.4272786307807642e-04,
//       4.4924554183813387e-02,  2.9179353882654718e-02,
//       8.9097263763830599e-03,  5.7870370370370072e-03,
//      -1.1449174941446927e-16, -5.7245874707234634e-17,
//    };
//    long i = 0;
//    struct gkyl_range_iter iter;
//    gkyl_range_iter_init(&iter, &local_tar);
//    while (gkyl_range_iter_next(&iter)) {
//      long linidx = gkyl_range_idx(&local_tar, iter.idx);
//      const double *f_c = gkyl_array_cfetch(distf_tar_ho, linidx);
//      for (int m=0; m<basis.num_basis; m++) {
//        TEST_CHECK( gkyl_compare(sol[i], f_c[m], 1e-12) );
//        TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[i], iter.idx[0], iter.idx[1]);
//        TEST_MSG("Produced: %.13e", f_c[m]);
//        i += 1;
//      }
//    }
//  }
//  else if (cells[0]     == 6 && cells[1]      == 12 &&
//           cells_tar[0] == 12 && cells_tar[1] == 12) {
//    // From 6x12 to 12x12.
//    // First 2 cells along vpar only.
//    const double sol[144] =  {
//       2.0924425582990402e-03,  9.8842276988000142e-04,
//      -2.5522653682347228e-05, -1.2056327160493826e-05,
//       0.0000000000000000e+00, -2.1684043449710089e-19,
//       1.9156164266117964e-03,  9.0489408510140948e-04,
//      -7.6567961047041644e-05, -3.6168981481481215e-05,
//      -9.7578195523695399e-19, -6.5052130349130266e-19,
//      
//       5.5164394718792854e-03,  9.8842276988000142e-04,
//      -6.7286996071642870e-05, -1.2056327160493826e-05,
//      -6.5052130349130266e-19, -5.4210108624275222e-19,
//       5.0502614883401898e-03,  9.0489408510140948e-04,
//      -2.0186098821492793e-04, -3.6168981481481215e-05,
//      -3.3610267347050637e-18, -9.7578195523695399e-19,
//      
//       7.7991040809327833e-03,  3.2947425662666952e-04,
//      -9.5129890997840399e-05, -4.0187757201647585e-06,
//      -4.3368086899420177e-19,  4.3368086899420177e-19,
//       7.1400248628257848e-03,  3.0163136170047259e-04,
//      -2.8538967299351870e-04, -1.2056327160493864e-05,
//      -4.3368086899420177e-18, -4.3368086899420177e-19,
//      
//       8.9404363854595319e-03,  3.2947425662666660e-04,
//      -1.0905133846093873e-04, -4.0187757201646043e-06,
//      -4.3368086899420177e-19,  0.0000000000000000e+00,
//       8.1849065500685849e-03,  3.0163136170046972e-04,
//      -3.2715401538281413e-04, -1.2056327160493699e-05,
//      -5.2041704279304213e-18, -4.3368086899420177e-19,
//      
//       8.9404363854595319e-03, -3.2947425662666725e-04,
//      -1.0905133846093873e-04,  4.0187757201644849e-06,
//      -8.6736173798840355e-19,  0.0000000000000000e+00,
//       8.1849065500685849e-03, -3.0163136170046988e-04,
//      -3.2715401538281456e-04,  1.2056327160493606e-05,
//      -4.7704895589362195e-18, -8.6736173798840355e-19,
//      
//       7.7991040809327825e-03, -3.2947425662667007e-04,
//      -9.5129890997840399e-05,  4.0187757201642968e-06,
//      -4.3368086899420177e-19,  0.0000000000000000e+00,
//       7.1400248628257865e-03, -3.0163136170047275e-04,
//      -2.8538967299351914e-04,  1.2056327160493569e-05,
//      -4.3368086899420177e-18, -6.5052130349130266e-19,
//      
//       5.5164394718792862e-03, -9.8842276988000142e-04,
//      -6.7286996071642654e-05,  1.2056327160493798e-05,
//      -1.9515639104739080e-18,  1.0842021724855044e-19,
//       5.0502614883401907e-03, -9.0489408510140948e-04,
//      -2.0186098821492793e-04,  3.6168981481481120e-05,
//      -3.6862873864507151e-18, -1.1926223897340549e-18,
//      
//       2.0924425582990449e-03, -9.8842276987999838e-04,
//      -2.5522653682347445e-05,  1.2056327160494414e-05,
//       1.7347234759768071e-18,  2.9273458657108620e-18,
//       1.9156164266118014e-03, -9.0489408510140644e-04,
//      -7.6567961047041861e-05,  3.6168981481481635e-05,
//      -5.4210108624275222e-19,  2.4936649967166602e-18,
//      
//      -2.4728866598079449e-03, -1.6473712831333374e-03,
//       3.0163136170047178e-05,  2.0093878600823011e-05,
//       1.9515639104739080e-18, -8.2399365108898337e-18,
//      -2.2639103223593862e-03, -1.5081568085023513e-03,
//       9.0489408510140097e-05,  6.0281635802468706e-05,
//       1.9515639104739080e-18,  4.3368086899420177e-19,
//      
//      -8.1795481824416833e-03, -1.6473712831333337e-03,
//       9.9770373485539264e-05,  2.0093878600825159e-05,
//       8.6736173798840355e-19,  2.1684043449710089e-19,
//      -7.4883187585733707e-03, -1.5081568085023476e-03,
//       2.9931112045661676e-04,  6.0281635802470603e-05,
//       7.1557343384043293e-18, -3.4694469519536142e-18,
//      
//      -1.5027542009602192e-02, -2.3063197963866674e-03,
//       1.8329905826413098e-04,  2.8131430041155208e-05,
//       3.4694469519536142e-18,  8.6736173798840355e-19,
//      -1.3757608882030174e-02, -2.1114195319032864e-03,
//       5.4989717479239056e-04,  8.4394290123459062e-05,
//       1.0408340855860843e-17,  1.7347234759768071e-18,
//      
//      -2.3016868141289417e-02, -2.3063197963866674e-03,
//       2.8074919050581973e-04,  2.8131430041149062e-05,
//       1.7347234759768071e-18,  8.6736173798840355e-19,
//      -2.1071780692729743e-02, -2.1114195319032864e-03,
//       8.4224757151745767e-04,  8.4394290123453438e-05,
//       1.3010426069826053e-17,  4.3368086899420177e-18,
//    };
//    long i = 0;
//    struct gkyl_range_iter iter;
//    gkyl_range_iter_init(&iter, &local_tar);
//    while (gkyl_range_iter_next(&iter)) {
//      if (iter.idx[1] < 3) {
//        long linidx = gkyl_range_idx(&local_tar, iter.idx);
//        const double *f_c = gkyl_array_cfetch(distf_tar_ho, linidx);
//        for (int m=0; m<basis.num_basis; m++) {
//          TEST_CHECK( gkyl_compare(sol[i], f_c[m], 1e-12) );
//          TEST_MSG("Expected: %.13e in cell (%d,%d)", sol[i], iter.idx[0], iter.idx[1]);
//          TEST_MSG("Produced: %.13e", f_c[m]);
//          i += 1;
//        }
//      }
//    }
//  }

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(moms_tar);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_velocity_map_release(gvm_tar);
  gkyl_gk_geometry_release(gk_geom_tar);
  gkyl_array_release(moms);
  gkyl_velocity_map_release(gvm);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(bmag);
  gkyl_array_release(distf);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_bmag);
  gkyl_proj_on_basis_release(proj_distf);
}

void test_1x1v_gk_ho()
{
//  int cells_do0[] = {8, 8};
//  int cells_tar0[] = {8, 4};
//  test_1x1v_gk(cells_do0, cells_tar0, 1, false);

  int cells_do1[] = {2, 4};
  int cells_tar1[] = {2, 2};
  test_1x1v_gk(cells_do1, cells_tar1, 1, false);

//  int cells_do1[] = {6, 12};
//  int cells_tar1[] = {12, 12};
//  test_1x1v_gk(cells_do1, cells_tar1, 1, false);
}

void test_1x1v_gk_dev()
{
  int cells_do0[] = {3, 3};
  int cells_tar0[] = {3, 4};
  test_1x1v_gk(cells_do0, cells_tar0, 1, true);

  int cells_do1[] = {6, 12};
  int cells_tar1[] = {12, 12};
  test_1x1v_gk(cells_do1, cells_tar1, 1, true);
}

TEST_LIST = {
  { "test_1x1v_gk_ho", test_1x1v_gk_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_1x1v_gk_dev", test_1x1v_gk_dev },
#endif
  { NULL, NULL },
};
