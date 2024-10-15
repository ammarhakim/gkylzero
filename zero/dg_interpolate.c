#include <gkyl_dg_interpolate.h>
#include <gkyl_dg_interpolate_priv.h>
#include <gkyl_alloc.h>

static int dg_interp_prime_factors(int n, int *pfs, int pfs_size)
{
  // Find the prime factos of number `n`, and put them into `pfs`. We assume
  // there are fewer than `pfs_size` prime factors.
  int pf_count = 0;
  int c = 2;
  while (n > 1) {
    if (n % c == 0) {
      pfs[pf_count] = c;
      pf_count++;
      assert(pf_count < pfs_size);
      n /= c;
    }
    else
      c++;
  }
  return pf_count;
}

static int dg_interp_index_stencil_map_refine(int dir, int idx,
  int num_cells, double dx_rat, int stencilIn)
{
  // Given an index 'idx' to a cell in the coarse grid with 'num_cells'
  // cells, return the index of the refinement stencil needed, within
  // the table that holds stencils. Here 'dx_rat' is the ratio of
  // donor to target cell length in the 'dir' direction.
  double remDecL = (idx-1)*dx_rat-floor((idx-1)*dx_rat);
  double remDecU = ceil(idx*dx_rat)-idx*dx_rat;
  int stencilOut = stencilIn;
  if ((idx == 1) ||   // First cell.
      (remDecL == 0) ||    // Interior cell with a left-boundary-like stencil.
      ((remDecL <= 0.5) && (remDecU <= 0.5))) {
    stencilOut = 2*stencilIn + pow(3,dir-1) - 1;
  }
  else if ((idx == num_cells) || // Last cell.
          (remDecU == 0)) { // Interior cell with a right-boundary-like stencil.
    stencilOut = 2*stencilIn + pow(3,dir-1);
  }
  return stencilOut;
}

static int dg_interp_index_stencil_map_coarsen(int dir, int idx,
  int num_cells, double dx_rat, int stencilIn)
{
  // Given an index 'idx' to a cell in the fine grid with 'num_cells'
  // cells, return the index of the coarsening stencil needed, within
  // the table that holds stencils. Here 'dx_rat' is the ratio of
  // donor to target cell length in the 'dir' direction.
  double remDecL = (idx-1)*dx_rat-floor(idx*dx_rat);
  double remDecU = ceil(idx*dx_rat)-idx*dx_rat;
  int stencilOut = stencilIn;
  if ((idx == 1) || // First cell.
      (remDecL == 0) || // Interior cell with a left-boundary-like stencil.
      ((remDecL > 0) && (remDecU > 0))) {
     stencilOut = 2*stencilIn + pow(3,dir-1) - 1;
  }
  else if ((idx == num_cells) || // Last cell.
          (remDecU == 0)) { // Interior cell with a right-boundary-like stencil.
     stencilOut = 2*stencilIn + pow(3,dir-1);
  }
  return stencilOut;
}


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
  up->interp = dg_interp_choose_gk_interp_kernel(pbasis);

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
  up->stencils = gkyl_malloc(pow(3,up->ndim)*sizeof(struct dg_interp_stencils));
  // Interior stencil (region 0).
  int sc = 0; // Current number of stencils.
  for (int d=0; d<up->ndim; d++) {
    up->stencils[sc].sz[d] = intStencilSize[d];
  }
  int offsets_lo[GKYL_MAX_DIM] = {0};
  int offsets_up[GKYL_MAX_DIM] = {-1};
  for (int d=0; d<up->ndim; d++)
    offsets_up[d] = up->stencils[sc].sz[d]-1;
  printf("\nsc=%d | offsets=%d:%d,%d:%d\n", sc, offsets_lo[0], offsets_up[0], offsets_lo[1], offsets_up[1]);
  gkyl_range_init(&up->stencils[sc].offsets, up->ndim, offsets_lo, offsets_up);
  sc++;
  int num_stencil_curr = 1; // Like sc but updated outside of the psI loop.
  // Other stencils.
  for (int dI=0; dI<up->ndim; dI++) {
    int stencils_added = 0;
    for (int psI=0; psI<num_stencil_curr; psI++) {
      for (int mp=-1; mp<2; mp += 2) {
        for (int d=0; d<up->ndim; d++) {
          up->stencils[sc].sz[d] = up->stencils[psI].sz[d];
        }
        up->stencils[sc].sz[dI] = floor(up->dxRat[dI]) + ceil(up->dxRat[dI] - floor(up->dxRat[dI]));

        for (int d=0; d<up->ndim; d++)
          offsets_up[d] = up->stencils[sc].sz[d]-1;
  printf("sc=%d | offsets=%d:%d,%d:%d\n", sc, offsets_lo[0], offsets_up[0], offsets_lo[1], offsets_up[1]);
        gkyl_range_init(&up->stencils[sc].offsets, up->ndim, offsets_lo, offsets_up);

        sc++;
        stencils_added++;
      }
    }
    num_stencil_curr += stencils_added;
  }

  return up;
}

void
gkyl_dg_interpolate_advance(gkyl_dg_interpolate* up,
  const struct gkyl_range *range_do, const struct gkyl_range *range_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
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
    struct gkyl_range_iter iter_tar_off;
    gkyl_range_iter_init(&iter_tar_off, &up->stencils[idx_sten].offsets);
    while (gkyl_range_iter_next(&iter_tar_off)) {
      for (int d=0; d<up->ndim; d++) {
        idx_tar[d] = idx_tar_ll[d] + iter_tar_off.idx[d];
      }

      gkyl_rect_grid_cell_center(&up->grid_tar, idx_tar, xc_tar);

      long linidx_tar = gkyl_range_idx(range_tar, idx_tar);
      double *ftar_c = gkyl_array_fetch(ftar, linidx_tar);

      up->interp(xc_do, xc_tar, up->grid_do.dx, up->grid_tar.dx, fdo_c, ftar_c);
    }
  }
}

void
gkyl_dg_interpolate_release(gkyl_dg_interpolate* up)
{
  // Release memory associated with this updater.
  gkyl_free(up->stencils);
  gkyl_free(up);
}
