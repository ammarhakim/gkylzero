/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_parproj_priv.h>
}

// CUDA kernel to set device pointers to l2g, RHS src and solution
// kernels. Doing function pointer stuff in here avoids troublesome
// cudaMemcpyFromSymbol.
__global__ static void
fem_parproj_set_cu_ker_ptrs(struct gkyl_fem_parproj_kernels* kers, enum gkyl_basis_type b_type,
  int dim, int poly_order)
{

  // Set l2g kernels.
  const local2global_kern_list *local2global_kernels;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      local2global_kernels = ser_loc2glob_list;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      local2global_kernels = ten_loc2glob_list;
      break;
    default:
      assert(false);
      break;
  }
  kers->l2g = local2global_kernels[dim].kernels[poly_order];

  // Set RHS stencil kernels.
  const srcstencil_kern_list *srcstencil_kernels;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      srcstencil_kernels = ser_srcstencil_list;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      srcstencil_kernels = ten_srcstencil_list;
      break;
    default:
      assert(false);
  }
  kers->srcker = srcstencil_kernels[dim].kernels[poly_order];

  // Set the get solution stencil kernel.
  const solstencil_kern_list *solstencil_kernels;
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      solstencil_kernels = ser_solstencil_list;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      solstencil_kernels = ten_solstencil_list;
      break;
    default:
      assert(false);
  }
  kers->solker = solstencil_kernels[dim].kernels[poly_order];

}

__global__ void
gkyl_fem_parproj_set_rhs_kernel(double *rhs_global, struct gkyl_array *rhs_local, struct gkyl_range range, struct gkyl_range perp_range2d, struct gkyl_range par_range1d, struct gkyl_fem_parproj_kernels *kers, long numnodes_global)
{
  int idx[PARPROJ_MAX_DIM];
  long globalidx[32];
  int parnum_cells = range.upper[range.ndim-1]-range.lower[range.ndim-1]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume;
       linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);
    const double *local_d = (const double*) gkyl_array_cfetch(rhs_local, start);

    int idx1d[] = {idx[range.ndim-1]};
    long paridx = gkyl_range_idx(&par_range1d, idx1d);
    kers->l2g(parnum_cells, paridx, globalidx);

    int idx2d[] = {perp_range2d.lower[0], perp_range2d.lower[0]};
    for (int d=0; d<range.ndim-1; d++) idx2d[d] = idx[d];
    long perpidx2d = gkyl_range_idx(&perp_range2d, idx2d);
    long perpProbOff = perpidx2d*numnodes_global;

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    kers->srcker(local_d, perpProbOff, globalidx, rhs_global);
  }
}

__global__ void
gkyl_fem_parproj_get_sol_kernel(struct gkyl_array *x_local, const double *x_global, struct gkyl_range range, struct gkyl_range perp_range2d, struct gkyl_range par_range1d, struct gkyl_fem_parproj_kernels *kers, long numnodes_global)
{
  int idx[GKYL_MAX_DIM];
  long globalidx[32];
  int parnum_cells = range.upper[range.ndim-1]-range.lower[range.ndim-1]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume;
       linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);
    double *local_d = (double*) gkyl_array_cfetch(x_local, start);

    int idx1d[] = {idx[range.ndim-1]};
    long paridx = gkyl_range_idx(&par_range1d, idx1d);
    kers->l2g(parnum_cells, paridx, globalidx);

    int idx2d[] = {perp_range2d.lower[0], perp_range2d.lower[0]};
    for (int d=0; d<range.ndim-1; d++) idx2d[d] = idx[d];
    long perpidx2d = gkyl_range_idx(&perp_range2d, idx2d);
    long perpProbOff = perpidx2d*numnodes_global;

    if (range.ndim==3)
    printf("tId=%d | bId=%d | blockDim=%d | gridDim=%d | linc1=%d | idx=(%d,%d,%d) | paridx=%d | perpidx2d=%d | perpProbOff=%d\n",threadIdx.x, blockIdx.x, blockDim.x, gridDim.x, linc1, idx[0], idx[1], idx[2], paridx, perpidx2d, perpProbOff);

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    kers->solker(x_global, perpProbOff, globalidx, local_d);
  }
}

void
fem_parproj_choose_kernels_cu(const struct gkyl_basis* basis, bool isperiodic, struct gkyl_fem_parproj_kernels *kers)
{
  fem_parproj_set_cu_ker_ptrs<<<1,1>>>(kers, basis->b_type, basis->ndim, basis->poly_order);
}

void
gkyl_fem_parproj_set_rhs_cu(gkyl_fem_parproj *up, const struct gkyl_array *rhsin)
{
  double *rhs_cu = gkyl_cusolver_get_rhs_ptr(up->prob_cu, 0);
  gkyl_cu_memset(rhs_cu, 0, sizeof(double)*up->numnodes_global*up->perp_range.volume);
  gkyl_fem_parproj_set_rhs_kernel<<<rhsin->nblocks, rhsin->nthreads>>>(rhs_cu, rhsin->on_dev, up->solve_range, up->perp_range2d, up->par_range1d, up->kernels_cu, up->numnodes_global);
}

void
gkyl_fem_parproj_solve_cu(gkyl_fem_parproj *up, struct gkyl_array *phiout)
{
  // do linear solve with cusolver
  gkyl_cusolver_solve(up->prob_cu);

  double *x_cu = gkyl_cusolver_get_sol_ptr(up->prob_cu, 0);

  gkyl_fem_parproj_get_sol_kernel<<<phiout->nblocks, phiout->nthreads>>>(phiout->on_dev, x_cu, up->solve_range, up->perp_range2d, up->par_range1d, up->kernels_cu, up->numnodes_global);
}
