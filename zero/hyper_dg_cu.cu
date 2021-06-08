#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_hyper_dg_set_update_vol_cu_kernel(gkyl_hyper_dg *hdg, int update_vol_term)
{
  gkyl_hyper_dg_set_update_vol(hdg, update_vol_term);
}

__global__ void
gkyl_hyper_dg_advance_cu_kernel(const gkyl_hyper_dg* hdg, const struct gkyl_range update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs, double* GKYL_RESTRICT maxs)
{
  int ndim = hdg->ndim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM];
  // integer used for selecting between left-edge zero-flux BCs and right-edge zero-flux BCs
  int edge;

  double maxs_old[GKYL_MAX_DIM];
  for (int i=0; i<hdg->ndim; ++i)
    maxs_old[i] = maxs[i];

  for(unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x; 
      linc1 < update_range.volume;
      linc1 += blockDim.x*gridDim.x)
  {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={0,0,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&update_range, linc1, idxc);
    gkyl_rect_grid_cell_center(&hdg->grid, idxc, xcc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&update_range, idxc);

    if (hdg->update_vol_term) {
      double cflr = hdg->equation->vol_term(
        hdg->equation, xcc, hdg->grid.dx, idxc,
        (double*) gkyl_array_cfetch(fIn, linc), (double*) gkyl_array_fetch(rhs, linc)
      );
      double *cflrate_d = (double*) gkyl_array_fetch(cflrate, linc);
      cflrate_d[0] += cflr; // frequencies are additive
    }
    
    for (int d=0; d<hdg->num_up_dirs; ++d) {
      int dir = hdg->update_dirs[d];
      gkyl_copy_int_arr(ndim, idxc, idxl);
      gkyl_copy_int_arr(ndim, idxc, idxr);
      // TODO: fix for arbitrary subrange
      if (hdg->zero_flux_flags[d] && (idxc[dir] == update_range.lower[dir] || idxc[dir] == update_range.upper[dir])) {
        edge = (idxc[dir] == update_range.lower[dir]) ? -1 : 1;
        // use idxl to store interior edge index (first index away from skin cell)
        idxl[dir] = idxl[dir]-edge;

        gkyl_rect_grid_cell_center(&hdg->grid, idxl, xcl);
        long linl = gkyl_range_idx(&update_range, idxl);

        double mdir = hdg->equation->boundary_surf_term(hdg->equation,
          dir, xcl, xcc, hdg->grid.dx, hdg->grid.dx,
          maxs_old[dir], idxl, idxc, edge,
          (double*) gkyl_array_cfetch(fIn, linl), (double*) gkyl_array_cfetch(fIn, linc),
          (double*) gkyl_array_fetch(rhs, linc)
        );
        maxs[dir] = fmax(maxs[dir], mdir);         
      }
      else {
        idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;
        gkyl_rect_grid_cell_center(&hdg->grid, idxl, xcl);
        gkyl_rect_grid_cell_center(&hdg->grid, idxr, xcr);
        long linl = gkyl_range_idx(&update_range, idxl); 
        long linr = gkyl_range_idx(&update_range, idxr);

        double mdir = hdg->equation->surf_term(hdg->equation,
          dir, xcl, xcc, xcr, hdg->grid.dx, hdg->grid.dx, hdg->grid.dx,
          maxs_old[dir], idxl, idxc, idxr,
          (double*) gkyl_array_cfetch(fIn, linl), (double*) gkyl_array_cfetch(fIn, linc), (double*) gkyl_array_cfetch(fIn, linr),
          (double*) gkyl_array_fetch(rhs, linc)
        );
        maxs[dir] = fmax(maxs[dir], mdir);
      }
    }
  }
}

// wrapper to call advance kernel on device
void
gkyl_hyper_dg_advance_cu(const gkyl_hyper_dg* hdg, const struct gkyl_range update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs, double* GKYL_RESTRICT maxs)
{
  gkyl_hyper_dg_advance_cu_kernel<<<update_range.nblocks, update_range.nthreads>>>(hdg, update_range,
    fIn->on_device, cflrate->on_device, rhs->on_device, maxs);
}

void
gkyl_hyper_dg_set_update_vol_cu(gkyl_hyper_dg *hdg, int update_vol_term)
{
  gkyl_hyper_dg_set_update_vol_cu_kernel<<<1,1>>>(hdg, update_vol_term);
}

gkyl_hyper_dg*
gkyl_hyper_dg_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation_cu,
  int num_up_dirs, int update_dirs[], int zero_flux_flags[], int update_vol_term)
{
  gkyl_hyper_dg *up = (gkyl_hyper_dg*) gkyl_malloc(sizeof(gkyl_hyper_dg));

  up->ndim = basis->ndim;
  up->numBasis = basis->numBasis;
  up->num_up_dirs = num_up_dirs;
  up->grid = *grid;

  for (int i=0; i<num_up_dirs; ++i) {
    up->update_dirs[i] = update_dirs[i];
    up->zero_flux_flags[i] = zero_flux_flags[i];
  }
  up->update_vol_term = update_vol_term;

  up->equation = equation_cu; // NB: RHS a pointer to device

  // copy host struct to device struct
  gkyl_hyper_dg *up_cu = (gkyl_hyper_dg*) gkyl_cu_malloc(sizeof(gkyl_hyper_dg));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_hyper_dg), GKYL_CU_MEMCPY_H2D);

  gkyl_free(up);

  return up_cu;
}
