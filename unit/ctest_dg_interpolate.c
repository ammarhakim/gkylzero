#include <gkyl_array.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>
#include <gkyl_velocity_map.h>
#include <gkyl_dg_interpolate.h>
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
  double B0 = tctx->B0;
  int vdim = tctx->vdim;

  double a = 2.0,  b = 5.0;
  double c[] = {a/12.0 - 0.5, 0.0};
  double d[] = {0.0, b/12.0 - 0.5};

//  fout[0] = (n0/pow(2.0*M_PI*vtsq,vdim/2.0)) * exp(-(pow(vpar-upar,2))/(2.0*vtsq));
  fout[0] = (pow(x,2)/2.0+c[0]*x+c[1])*(pow(vpar,2)/2.0+d[0]*vpar+d[1]);
}

void
test_1x1v(const int *cells, const int *cells_tar, int poly_order, bool use_gpu)
{
  const int cdim = 1,  vdim = 1;
  double x_min = 0.0;
  double x_max = 1.0;
  double vpar_min = 0.0;
  double vpar_max = 1.0;
  double lower[] = {x_min, vpar_min}, upper[] = {x_max, vpar_max};

  const int ndim = cdim+vdim;

  struct test_ctx proj_ctx = {
    .n0 = 1.0, // Density.
    .upar = 0, // Parallel flow speed.
    .B0 = 1.0, // Magnetic field.
    .vdim = vdim, // Number of velocity space dimensions.
    .vpar_max = vpar_max, // Maximum vpar of the grid.
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
  struct gkyl_array *bmag_ho, *bmag;
  bmag = mkarr(use_gpu, confBasis.num_basis, confLocal_ext.volume);
  bmag_ho = use_gpu? mkarr(false, bmag->ncomp, bmag->size)
                   : gkyl_array_acquire(bmag);
  gkyl_proj_on_basis *proj_bmag = gkyl_proj_on_basis_new(&confGrid, &confBasis,
    poly_order+1, 1, eval_bmag_1x, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_bmag, 0.0, &confLocal, bmag_ho);
  gkyl_array_copy(bmag, bmag_ho);

  // Create distribution function arrays.
  struct gkyl_array *distf_ho, *distf;
  distf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  distf_ho = use_gpu? mkarr(false, distf->ncomp, distf->size)
                    : gkyl_array_acquire(distf);
  gkyl_proj_on_basis *proj_distf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, eval_distf_1x1v, &proj_ctx);
  gkyl_proj_on_basis_advance(proj_distf, 0.0, &local, distf_ho);
  gkyl_array_copy(distf, distf_ho);

  // Write donor distribution function to file.
  char fname0[1024];
  sprintf(fname0, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_do.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
  gkyl_grid_sub_array_write(&grid, &local, NULL, distf_ho, fname0);

  // Initialize geometry
  struct gkyl_gk_geometry_inp geometry_input = {
    .geometry_id = GKYL_MAPC2P,
    .world = {0.0},  .mapc2p = mapc2p,  .c2p_ctx = 0,
    .bmag_func = eval_bmag_1x,  .bmag_ctx = &proj_ctx,
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

  // Velocity space mapping.
  struct gkyl_mapc2p_inp c2p_in = { };
  struct gkyl_velocity_map *gvm = gkyl_velocity_map_new(c2p_in, grid, velGrid,
    local, local_ext, velLocal, velLocal_ext, use_gpu);

  // Target grid and range.
  struct gkyl_rect_grid grid_tar;
  gkyl_rect_grid_init(&grid_tar, ndim, lower, upper, cells_tar);
  struct gkyl_range local_tar, local_tar_ext; // Target local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid_tar, ghost, &local_tar_ext, &local_tar);

  // Target field.
  struct gkyl_array *distf_tar_ho, *distf_tar;
  distf_tar = mkarr(use_gpu, basis.num_basis, local_tar_ext.volume);
  distf_tar_ho = use_gpu? mkarr(false, distf_tar->ncomp, distf_tar->size)
                    : gkyl_array_acquire(distf_tar);

  // Create the interpolation operator and interpolate onto the target grid.
  struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(cdim, &basis,
    &grid, &grid_tar, use_gpu);

  gkyl_dg_interpolate_advance(interp, &local, &local_tar, distf, distf_tar);
  gkyl_array_copy(distf_tar_ho, distf_tar);

  // Write donor distribution function to file.
  char fname1[1024];
  sprintf(fname1, "ctest_dg_interp_1x1v_p%d_N%dx%d-N%dx%d_tar.gkyl", poly_order, cells[0], cells[1], cells_tar[0], cells_tar[1]);
  gkyl_grid_sub_array_write(&grid_tar, &local_tar, NULL, distf_tar_ho, fname1);

  gkyl_dg_interpolate_release(interp);
  gkyl_array_release(distf_tar);
  gkyl_array_release(distf_tar_ho);
  gkyl_velocity_map_release(gvm);
  gkyl_gk_geometry_release(gk_geom);
  gkyl_array_release(bmag);
  gkyl_array_release(distf);
  gkyl_array_release(bmag_ho);
  gkyl_array_release(distf_ho);
  gkyl_proj_on_basis_release(proj_bmag);
  gkyl_proj_on_basis_release(proj_distf);
}

void test_1x1v_ho()
{
  int cells_do[] = {6, 12};
  int cells_tar[] = {12, 12};
//  int cells_do[] = {3, 3};
//  int cells_tar[] = {3, 4};
  test_1x1v(cells_do, cells_tar, 1, false);
}

void test_1x1v_dev()
{
  int cells_do[] = {6, 12};
  int cells_tar[] = {12, 12};
  test_1x1v(cells_do, cells_tar, 1, true);
}

TEST_LIST = {
  { "test_1x1v_ho", test_1x1v_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_1x1v_dev", test_1x1v_dev },
#endif
  { NULL, NULL },
};
