#include <gkyl_dg_interpolate.h>
#include <gkyl_dg_interpolate_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>

struct gkyl_dg_interpolate*
gkyl_dg_interpolate_new(int cdim, const struct gkyl_basis *pbasis,
  const struct gkyl_rect_grid *grid_do, const struct gkyl_rect_grid *grid_tar,
  bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_dg_interpolate *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;
  up->ndim = pbasis->ndim;
  up->grid_do = *grid_do;
  up->grid_tar = *grid_tar;

  // Perform some basic checks:
  assert(grid_do->ndim == grid_tar->ndim);
  for (int d=0; d<up->ndim; d++) {
    assert(fabs(grid_do->lower[d] - grid_tar->lower[d]) < 1e-14);
    assert(fabs(grid_do->upper[d] - grid_tar->upper[d]) < 1e-14);

    const int num_prime_facs_max = 32;
    int prime_factors[num_prime_facs_max];
    int num_prime_facs_do = dg_interp_prime_factors(grid_do->cells[d], prime_factors, num_prime_facs_max);
    for (int k=0; k<num_prime_facs_do; k++)
      assert(prime_factors[k] == 2 || prime_factors[k] == 3); // Only (2^a)*(3^b) grids supported.
 
    int num_prime_facs_tar = dg_interp_prime_factors(grid_tar->cells[d], prime_factors, num_prime_facs_max);
    for (int k=0; k<num_prime_facs_tar; k++)
      assert(prime_factors[k] == 2 || prime_factors[k] == 3); // Only (2^a)*(3^b) grids supported.
  }

  // Choose kernels that translates the DG coefficients.
  up->interp = dg_interp_choose_gk_interp_kernel(*pbasis);

  for (int d=0; d<up->ndim; d++) {
    // Ratio of cell-lengths in each direction.
    up->dxRat[d] = grid_do->dx[d]/grid_tar->dx[d];

    // Map from grid to stencil index in each direction.
    if (up->dxRat[d] > 1)
      up->grid2stencil[d] = dg_interp_index_stencil_map_refine;
    else
      up->grid2stencil[d] = dg_interp_index_stencil_map_coarsen;
  }

  // Interior (away from boundaries) stencil size, in each direction.
  int intStencilSize[GKYL_MAX_DIM] = {-1};
  for (int d=0; d<up->ndim; d++) {
    if (up->dxRat[d] > 1) {
      // Mesh refinement.
      // Brute force search. Start with the size of the boundary stencil.
      int maxSize = floor(up->dxRat[d]) + ceil( up->dxRat[d]-floor(up->dxRat[d]) );
      for (int i=2; i<grid_do->cells[d]; i++) {
         double decimalL = 1-((i-1)*up->dxRat[d]-floor((i-1)*up->dxRat[d]));
         double decimalU = 1-(ceil(i*up->dxRat[d])-i*up->dxRat[d]);
         int currSize = floor(up->dxRat[d]-decimalL-decimalU) + ceil(decimalL) + ceil(decimalU);
         maxSize = GKYL_MAX2(maxSize, currSize);
      }
      intStencilSize[d] = maxSize;
    }
    else {
      // Mesh coarsening, or no change in resolution (dxRat=1).
      intStencilSize[d] = 1+ceil(1*(1/up->dxRat[d]-floor(1/up->dxRat[d])));
    }
  }

  // This updater will loop through the donor grid. At each cell it will give
  // the interpolation kernel the target grid cells, one at a time, and the
  // kernel will add the corresponding contributions to the target field DG
  // coefficients in that each target grid cell.

  // There are 3^ndim stencils or regions. For example, in 2D the regions are
  //   .-----.-----.-----.
  //   |     |     |     |
  //   |  6  |  4  |  8  |
  //   |     |     |     |
  //   .-----.-----.-----.
  //   |     |     |     |
  //   |  1  |  0  |  2  |
  //   |     |     |     |
  //   .-----.-----.-----.
  //   |     |     |     |
  //   |  5  |  3  |  7  |
  //   |     |     |     |
  //   .-----.-----.-----.
  // For each region we will make a list of the number of target
  // cells to fetch contributions from in each direction.
  up->offset_upper = gkyl_malloc(pow(3,up->ndim) * up->ndim * sizeof(int));
  int stencil_sizes[((int) round(pow(3,up->ndim))) * up->ndim];
  int sc = 0; // Current number of stencils.
  // Interior stencil (region 0).
  for (int d=0; d<up->ndim; d++) {
    stencil_sizes[sc*up->ndim+d] = intStencilSize[d];
    up->offset_upper[sc*up->ndim+d] = stencil_sizes[sc*up->ndim+d];
  }
  sc++;
  int num_stencil_curr = 1; // Like sc but updated outside of the psI loop.
  // Other stencils.
  for (int dI=0; dI<up->ndim; dI++) {
    int stencils_added = 0;
    for (int psI=0; psI<num_stencil_curr; psI++) {
      for (int mp=-1; mp<2; mp += 2) {
        for (int d=0; d<up->ndim; d++) {
          stencil_sizes[sc*up->ndim+d] = stencil_sizes[psI*up->ndim+d];
        }
        stencil_sizes[sc*up->ndim+dI] = floor(up->dxRat[dI]) + ceil(up->dxRat[dI] - floor(up->dxRat[dI]));

        for (int d=0; d<up->ndim; d++)
          up->offset_upper[sc*up->ndim+d] = stencil_sizes[sc*up->ndim+d];

        sc++;
        stencils_added++;
      }
    }
    num_stencil_curr += stencils_added;
  }

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    // Allocate a device copy of the updater.
    gkyl_dg_interpolate_new_cu(up, *pbasis);
  }
#endif

  return up;
}

void
gkyl_dg_interpolate_advance(gkyl_dg_interpolate* up,
  const struct gkyl_range *range_do, const struct gkyl_range *range_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{

  gkyl_array_clear_range(ftar, 0.0, range_tar);

#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_dg_interpolate_advance_cu(up, range_do, range_tar, fdo, ftar);
    return;
  }
#endif

  int idx_tar[GKYL_MAX_DIM] = {-1};
  int idx_tar_ll[GKYL_MAX_DIM] = {-1};
  double xc_do[GKYL_MAX_DIM];
  double xc_tar[GKYL_MAX_DIM];

  // Loop over the donor range.
  struct gkyl_range_iter iter_do;
  gkyl_range_iter_init(&iter_do, range_do);
  while (gkyl_range_iter_next(&iter_do)) {

    int *idx_do = iter_do.idx;

    gkyl_rect_grid_cell_center(&up->grid_do, idx_do, xc_do);

    long linidx_do = gkyl_range_idx(range_do, idx_do);
    const double *fdo_c = gkyl_array_cfetch(fdo, linidx_do);

    // Compute the target-grid index of the "lower-left" corner this cell contributes to.
    for (int d=0; d<up->ndim; d++) {
      double eveOI = up->dxRat[d]*(idx_do[d]-1);
      idx_tar_ll[d] = ceil(eveOI)+((int) ceil(eveOI-floor(eveOI))+1) % 2;
    }

    // Get the index to the stencil for this donor cell.
    int idx_sten = 1;
    for (int d=0; d<up->ndim; d++) {
      idx_sten = up->grid2stencil[d](d, idx_do[d], up->grid_do.cells[d], up->dxRat[d], idx_sten);
    }
    idx_sten -= 1;

    // Loop over the target-grid cells this donor cell contributes to.
    // Since we can't use range_iter's in CUDA kernels we have to use 4 nested
    // loops instead.
    for (int dI=0; dI<up->ndim; dI++) {

      int offset[GKYL_MAX_DIM] = {0};

      for (int oI=0; oI<up->offset_upper[idx_sten*up->ndim+dI]; oI++) {
        offset[dI] = oI;

        for (int eI=0; eI<dI; eI++) {

          for (int pI=0; pI<up->offset_upper[idx_sten*up->ndim+eI]; pI++) {
            offset[eI] = pI;

            for (int d=0; d<up->ndim; d++)
              idx_tar[d] = idx_tar_ll[d] + offset[d];

            gkyl_rect_grid_cell_center(&up->grid_tar, idx_tar, xc_tar);

            long linidx_tar = gkyl_range_idx(range_tar, idx_tar);
            double *ftar_c = gkyl_array_fetch(ftar, linidx_tar);

            up->interp(xc_do, xc_tar, up->grid_do.dx, up->grid_tar.dx, fdo_c, ftar_c);

          }
        }
      }
    }

  }
}

void
gkyl_dg_interpolate_release(gkyl_dg_interpolate* up)
{
  // Release memory associated with this updater.
  if (!up->use_gpu)
    gkyl_free(up->offset_upper);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->offset_upper);
    if (GKYL_IS_CU_ALLOC(up->flags))
      gkyl_cu_free(up->on_dev);
  }
#endif
  gkyl_free(up);
}
