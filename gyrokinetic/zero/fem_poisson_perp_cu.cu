/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_fem_poisson_perp.h>
#include <gkyl_fem_poisson_perp_priv.h>
}

// CUDA kernel to set device pointers to l2g kernel function.
// Doing function pointer stuff in here avoids troublesome
// cudaMemcpyFromSymbol.
__global__ static void
fem_poisson_perp_set_cu_l2gker_ptrs(struct gkyl_fem_poisson_perp_kernels* kers, int ndim, enum gkyl_basis_type b_type,
  int poly_order, const int *bckey)
{
  int ndim_perp = ndim-1;

  // Set l2g kernels.
  const local2global_kern_bcx_list_2x *local2global_2x_kernels;
  const local2global_kern_bcx_list_3x *local2global_3x_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      local2global_2x_kernels = ser_loc2glob_list_2x;
      local2global_3x_kernels = ser_loc2glob_list_3x;
      break;
    default:
      assert(false);
      break;
  }

  for (int k=0; k<GKYL_IPOW(2,ndim_perp); k++)
    kers->l2g[k] = ndim == 2? CK2x(local2global_2x_kernels, poly_order, k, bckey[0])
                            : CK3x(local2global_3x_kernels, poly_order, k, bckey[0], bckey[1]);
}

// CUDA kernel to set device pointers to RHS src and solution kernels.
__global__ static void
fem_poisson_perp_set_cu_ker_ptrs(struct gkyl_fem_poisson_perp_kernels* kers, enum gkyl_basis_type b_type,
  int ndim, int poly_order, const int *bckey)
{
  int ndim_perp = ndim-1;

  // Set RHS stencil kernels.
  const srcstencil_kern_bcx_list_2x *srcstencil_2x_kernels;
  const srcstencil_kern_bcx_list_3x *srcstencil_3x_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
        srcstencil_2x_kernels = ser_srcstencil_list_2x;
        srcstencil_3x_kernels = ser_srcstencil_list_3x;
      break;
//    case GKYL_BASIS_MODAL_TENSOR:
//      break;
    default:
      assert(false);
  }

  for (int k=0; k<GKYL_IPOW(3,ndim_perp); k++)
    kers->srcker[k] = ndim == 2? CK2x(srcstencil_2x_kernels, poly_order, k, bckey[0])
                               : CK3x(srcstencil_3x_kernels, poly_order, k, bckey[0], bckey[1]);

  // Set the get solution stencil kernel.
  const solstencil_kern_list *solstencil_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
        solstencil_kernels = ser_solstencil_list;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      break;
    default:
      assert(false);
  }

  kers->solker = solstencil_kernels[ndim].kernels[poly_order];
}

void
fem_poisson_perp_choose_kernels_cu(const struct gkyl_basis* basis, const struct gkyl_poisson_bc *bcs,
  const bool *isdirperiodic, struct gkyl_fem_poisson_perp_kernels *kers)
{

  int poly_order = basis->poly_order;
  int ndim = basis->ndim;
  int ndim_perp = ndim-1;

  int bckey[GKYL_MAX_CDIM] = {-1,-1,-1};
  for (int d=0; d<ndim_perp; d++) bckey[d] = isdirperiodic[d] ? 0 : 1;
  int *bckey_d = (int *) gkyl_cu_malloc(GKYL_MAX_CDIM*sizeof(int));
  gkyl_cu_memcpy(bckey_d, bckey, GKYL_MAX_CDIM*sizeof(int), GKYL_CU_MEMCPY_H2D);

  fem_poisson_perp_set_cu_l2gker_ptrs<<<1,1>>>(kers, ndim, basis->b_type, poly_order, bckey_d);

  for (int d=0; d<ndim_perp; d++) {
         if (bcs->lo_type[d]==GKYL_POISSON_PERIODIC  && bcs->up_type[d]==GKYL_POISSON_PERIODIC ) { bckey[d] = 0; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 1; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_NEUMANN  ) { bckey[d] = 2; }
    else if (bcs->lo_type[d]==GKYL_POISSON_NEUMANN   && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 3; }
    else { assert(false); }
  };
  gkyl_cu_memcpy(bckey_d, bckey, GKYL_MAX_CDIM*sizeof(int), GKYL_CU_MEMCPY_H2D);

  fem_poisson_perp_set_cu_ker_ptrs<<<1,1>>>(kers, basis->b_type, ndim, poly_order, bckey_d);

  gkyl_cu_free(bckey_d);
}

__global__ void
gkyl_fem_poisson_perp_set_rhs_kernel(struct gkyl_array *epsilon, const double *dx, double *rhs_global,
  struct gkyl_array *rhs_local, const struct gkyl_range range, struct gkyl_range par_range1d,
  const double *bcvals, struct gkyl_fem_poisson_perp_kernels *kers, long numnodes_global)
{
  int idx[GKYL_MAX_CDIM];  int idx0[GKYL_MAX_CDIM];  int num_cells[GKYL_MAX_CDIM];
  long globalidx[32];

  int ndim = range.ndim;
  int ndim_perp = ndim-1;

  for (int d=0; d<ndim; d++) num_cells[d] = range.upper[d]-range.lower[d]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume; linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);

    const double *local_d = (const double*) gkyl_array_cfetch(rhs_local, start);
    const double *epsilon_d = (const double*) gkyl_array_cfetch(epsilon, start);

    int keri = idx_to_inup_ker(ndim_perp, num_cells, idx);
    for (size_t d=0; d<ndim; d++) idx0[d] = idx[d]-1;
    kers->l2g[keri](num_cells, idx0, globalidx);

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    keri = idx_to_inloup_ker(ndim_perp, num_cells, idx);
    int idx1d[] = {idx[ndim-1]};
    long paridx = gkyl_range_idx(&par_range1d, idx1d);
    long parProbOff = paridx*numnodes_global;
    kers->srcker[keri](epsilon_d, dx, local_d, bcvals, parProbOff, globalidx, rhs_global);
  }
}

void
gkyl_fem_poisson_perp_set_rhs_cu(gkyl_fem_poisson_perp *up, struct gkyl_array *rhsin)
{
  gkyl_culinsolver_clear_rhs(up->prob_cu, 0);
  double *rhs_cu = gkyl_culinsolver_get_rhs_ptr(up->prob_cu, 0);

  gkyl_fem_poisson_perp_set_rhs_kernel<<<rhsin->nblocks, rhsin->nthreads>>>(up->epsilon->on_dev, up->dx_cu,
    rhs_cu, rhsin->on_dev, *up->solve_range, up->par_range1d, up->bcvals_cu, up->kernels_cu, up->numnodes_global);
}

__global__ void
gkyl_fem_poisson_perp_get_sol_kernel(struct gkyl_array *x_local, const double *x_global, struct gkyl_range range,
  struct gkyl_range par_range1d, struct gkyl_fem_poisson_perp_kernels *kers, long numnodes_global)
{
  int idx[GKYL_MAX_CDIM];  int idx0[GKYL_MAX_CDIM];  int num_cells[GKYL_MAX_CDIM];
  long globalidx[32];

  int ndim = range.ndim;
  int ndim_perp = ndim-1;

  for (int d=0; d<ndim; d++) num_cells[d] = range.upper[d]-range.lower[d]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume; linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);

    double *local_d = (double*) gkyl_array_cfetch(x_local, start);

    int keri = idx_to_inup_ker(ndim_perp, num_cells, idx);
    for (size_t d=0; d<ndim; d++) idx0[d] = idx[d]-1;
    kers->l2g[keri](num_cells, idx0, globalidx);

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    int idx1d[] = {idx[ndim-1]};
    long paridx = gkyl_range_idx(&par_range1d, idx1d);
    long parProbOff = paridx*numnodes_global;
    kers->solker(x_global, parProbOff, globalidx, local_d);
  }
}

void
gkyl_fem_poisson_perp_solve_cu(struct gkyl_fem_poisson_perp *up, struct gkyl_array *phiout)
{
  gkyl_culinsolver_solve(up->prob_cu);
  double *x_cu = gkyl_culinsolver_get_sol_ptr(up->prob_cu, 0);

  gkyl_fem_poisson_perp_get_sol_kernel<<<phiout->nblocks, phiout->nthreads>>>(phiout->on_dev,
    x_cu, *up->solve_range, up->par_range1d, up->kernels_cu, up->numnodes_global);
}
