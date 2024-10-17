/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_dg_interpolate_priv.h>
#include <gkyl_util.h>
}

// CUDA kernel to set device pointer to interpolating kernel.
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_interp_set_cu_dev_ptrs(struct gkyl_dg_interpolate *up, struct gkyl_basis pbasis, const double *dxRat)
{
  up->interp = dg_interp_choose_gk_interp_kernel(pbasis);

  for (int d=0; d<up->ndim; d++) {
    // Map from grid to stencil index in each direction.
    if (dxRat[d] > 1)
      up->grid2stencil[d] = dg_interp_index_stencil_map_refine;
    else
      up->grid2stencil[d] = dg_interp_index_stencil_map_coarsen;
  }
}

void
gkyl_dg_interpolate_new_cu(struct gkyl_dg_interpolate *up, struct gkyl_basis pbasis)
{
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  size_t sz = pow(3,up->ndim) * up->ndim * sizeof(int);
  int *offset_upper_cu = (int *) gkyl_cu_malloc(sz);
  gkyl_cu_memcpy(offset_upper_cu, up->offset_upper, sz, GKYL_CU_MEMCPY_H2D);
  gkyl_free(up->offset_upper);
  up->offset_upper = offset_upper_cu;

  struct gkyl_dg_interpolate *up_cu = (struct gkyl_dg_interpolate*) gkyl_cu_malloc(sizeof(gkyl_dg_interpolate));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_interpolate), GKYL_CU_MEMCPY_H2D);

  dg_interp_set_cu_dev_ptrs<<<1,1>>>(up_cu, pbasis, up_cu->dxRat); 

  // Set parent on_dev pointer.
  up->on_dev = up_cu;
}

__global__ static void
gkyl_dg_interpolate_advance_cu_ker(struct gkyl_dg_interpolate *up,
  const struct gkyl_range range_do, const struct gkyl_range range_tar,
  const struct gkyl_array* GKYL_RESTRICT fdo, struct gkyl_array* GKYL_RESTRICT ftar)
{
  int idx_do[GKYL_MAX_DIM] = {-1};
  int idx_tar[GKYL_MAX_DIM] = {-1};
  int idx_tar_ll[GKYL_MAX_DIM] = {-1};
  double xc_do[GKYL_MAX_DIM];
  double xc_tar[GKYL_MAX_DIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < range_do.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&range_do, tid, idx_do);

    gkyl_rect_grid_cell_center(&up->grid_do, idx_do, xc_do);

    long linidx_do = gkyl_range_idx(&range_do, idx_do);
    const double *fdo_c = (const double *) gkyl_array_cfetch(fdo, linidx_do);

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

            long linidx_tar = gkyl_range_idx(&range_tar, idx_tar);
            double *ftar_c = (double *) gkyl_array_fetch(ftar, linidx_tar);

            up->interp(xc_do, xc_tar, up->grid_do.dx, up->grid_tar.dx, fdo_c, ftar_c);

          }
        }
      }
    }
  }
}


void
gkyl_dg_interpolate_advance_cu(gkyl_dg_interpolate* up,
  const struct gkyl_range *range_do, const struct gkyl_range *range_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  int nblocks = range_do->nblocks, nthreads = range_do->nthreads;

  gkyl_dg_interpolate_advance_cu_ker<<<nblocks, nthreads>>>
    (up->on_dev, *range_do, *range_tar, fdo->on_dev, ftar->on_dev);
}
