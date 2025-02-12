/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_fem_poisson.h>
#include <gkyl_fem_poisson_priv.h>
}

// CUDA kernel to set device pointers to l2g kernel function.
// Doing function pointer stuff in here avoids troublesome
// cudaMemcpyFromSymbol.
__global__ static void
fem_poisson_set_cu_l2gker_ptrs(struct gkyl_fem_poisson_kernels* kers, enum gkyl_basis_type b_type,
  int dim, int poly_order, const int *bckey)
{

  // Set l2g kernels.
  const local2global_kern_bcx_list_1x *local2global_1x_kernels;
  const local2global_kern_bcx_list_2x *local2global_2x_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      local2global_1x_kernels = ser_loc2glob_list_1x;
      local2global_2x_kernels = ser_loc2glob_list_2x;
      break;
    default:
      assert(false);
      break;
  }

  for (int k=0; k<(int)(pow(2,dim)+0.5); k++) {
    if (dim == 1) {
      kers->l2g[k] = CK1(local2global_1x_kernels, poly_order, k, bckey[0]);
    } else if ( dim == 2) {
      kers->l2g[k] = CK2(local2global_2x_kernels, poly_order, k, bckey[0], bckey[1]);
//    } else if (dim == 3) {
//      kers->l2g[k] = CK3(ser_loc2glob_list_3x, poly_order, k, bckey[0], bckey[1], bckey[2]);
    }
  }

}

// CUDA kernel to set device pointers to RHS src and solution kernels.
__global__ static void
fem_poisson_set_cu_ker_ptrs(struct gkyl_fem_poisson_kernels* kers, enum gkyl_basis_type b_type,
  int dim, int poly_order, const int *bckey, bool isvareps)
{

  // Set RHS stencil kernels.
  const srcstencil_kern_bcx_list_1x *srcstencil_1x_kernels;
  const srcstencil_kern_bcx_list_2x *srcstencil_2x_kernels;

  if (isvareps) {
    switch (b_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
          srcstencil_1x_kernels = ser_srcstencil_vareps_list_1x;
          srcstencil_2x_kernels = ser_srcstencil_vareps_list_2x;
        break;
//      case GKYL_BASIS_MODAL_TENSOR:
//        break;
      default:
        assert(false);
    }
  } else {
    switch (b_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
          srcstencil_1x_kernels = ser_srcstencil_consteps_list_1x;
          srcstencil_2x_kernels = ser_srcstencil_consteps_list_2x;
        break;
//      case GKYL_BASIS_MODAL_TENSOR:
//        break;
      default:
        assert(false);
    }
  }

  for (int k=0; k<(int)(pow(3,dim)+0.5); k++) {
    if (dim == 1) {
      kers->srcker[k] = CK1(srcstencil_1x_kernels, poly_order, k, bckey[0]);
    } else if (dim == 2) {
      kers->srcker[k] = CK2(srcstencil_2x_kernels, poly_order, k, bckey[0], bckey[1]);
//  } else if (dim == 3) {
//    kers->srcker[k] = CK3(srcstencil_3x_kernels, poly_order, k, bckey[0], bckey[1], bckey[2]);
    }
  }

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

  kers->solker = solstencil_kernels[dim].kernels[poly_order];

}

void
fem_poisson_choose_kernels_cu(const struct gkyl_basis* basis, const struct gkyl_poisson_bc *bcs, bool isvareps, const bool *isdirperiodic, struct gkyl_fem_poisson_kernels *kers)
{

  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  int bckey[GKYL_MAX_CDIM] = {-1};
  for (int d=0; d<basis->ndim; d++) bckey[d] = isdirperiodic[d] ? 0 : 1;
  int *bckey_d = (int *) gkyl_cu_malloc(sizeof(int[GKYL_MAX_CDIM]));
  gkyl_cu_memcpy(bckey_d, bckey, sizeof(int[GKYL_MAX_CDIM]), GKYL_CU_MEMCPY_H2D);

  fem_poisson_set_cu_l2gker_ptrs<<<1,1>>>(kers, basis->b_type, dim, poly_order, bckey_d);
  
  for (int d=0; d<basis->ndim; d++) {
    if (bcs->lo_type[d]==GKYL_POISSON_PERIODIC && bcs->up_type[d]==GKYL_POISSON_PERIODIC) { bckey[d] = 0; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 1; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_NEUMANN) { bckey[d] = 2; }
    else if (bcs->lo_type[d]==GKYL_POISSON_NEUMANN && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 3; }
    else if (bcs->lo_type[d]==GKYL_POISSON_DIRICHLET && bcs->up_type[d]==GKYL_POISSON_ROBIN) { bckey[d] = 4; }
    else if (bcs->lo_type[d]==GKYL_POISSON_ROBIN && bcs->up_type[d]==GKYL_POISSON_DIRICHLET) { bckey[d] = 5; }
    else { assert(false); }
  };
  gkyl_cu_memcpy(bckey_d, bckey, sizeof(int[GKYL_MAX_CDIM]), GKYL_CU_MEMCPY_H2D);

  fem_poisson_set_cu_ker_ptrs<<<1,1>>>(kers, basis->b_type, dim, poly_order, bckey_d, isvareps);

  gkyl_cu_free(bckey_d);
}

__global__ void
gkyl_fem_poisson_set_rhs_kernel(struct gkyl_array *epsilon, bool isvareps, const double *dx, double *rhs_global, struct gkyl_array *rhs_local, struct gkyl_range range, 
  const double *bcvals, struct gkyl_fem_poisson_kernels *kers)
{
  int idx[GKYL_MAX_CDIM];  int idx0[GKYL_MAX_CDIM];  int num_cells[GKYL_MAX_CDIM];
  long globalidx[32];
  for (int d=0; d<GKYL_MAX_CDIM; d++) num_cells[d] = range.upper[d]-range.lower[d]+1;

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
    const double *epsilon_d = isvareps? (const double*) gkyl_array_cfetch(epsilon, start)
                                      : (const double*) gkyl_array_cfetch(epsilon, 0);

    int keri = idx_to_inup_ker(range.ndim, num_cells, idx);
    for (size_t d=0; d<range.ndim; d++) idx0[d] = idx[d]-1;
    kers->l2g[keri](num_cells, idx0, globalidx);

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    keri = idx_to_inloup_ker(range.ndim, num_cells, idx);
    kers->srcker[keri](epsilon_d, dx, local_d, bcvals, globalidx, rhs_global);
  }
}
__global__ void
gkyl_fem_poisson_apply_target_bias_kernel(double *rhs_global, struct gkyl_range range, bool is_z_edge, 
  int idxLCFS_m, double target_corner_bias, struct gkyl_fem_poisson_kernels *kers)
{  
  int idx[GKYL_MAX_CDIM];  int idx0[GKYL_MAX_CDIM];  
  long globalidx[32];

  int num_cells[GKYL_MAX_CDIM];
  for (int d=0; d<GKYL_MAX_CDIM; d++) num_cells[d] = range.upper[d]-range.lower[d]+1;

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
       linc1 < range.volume;
       linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    int keri = idx_to_inup_ker(range.ndim, num_cells, idx);
    for (size_t d=0; d<range.ndim; d++) idx0[d] = idx[d]-1;
    kers->l2g[keri](num_cells, idx0, globalidx);

    // modify the RHS for x=LCFS corresponding lines
    // we overwrite the left nodes (globalidx[0]) of the cells to the right of the LCFS
    // i.e. cells with x coordinate (idx[0]) equal to idxLCFS_m+1
    if (is_z_edge && idx[0] == idxLCFS_m+1) 
        atomicExch((unsigned long long int*) &rhs_global[globalidx[0]],__double_as_longlong(target_corner_bias));
  }
}

__global__ void
gkyl_fem_poisson_get_sol_kernel(struct gkyl_array *x_local, const double *x_global, struct gkyl_range range, struct gkyl_fem_poisson_kernels *kers)
{
  int idx[GKYL_MAX_CDIM];  int idx0[GKYL_MAX_CDIM];  int num_cells[GKYL_MAX_CDIM];
  long globalidx[32];
  for (int d=0; d<GKYL_MAX_CDIM; d++) num_cells[d] = range.upper[d]-range.lower[d]+1;

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

    int keri = idx_to_inup_ker(range.ndim, num_cells, idx);
    for (size_t d=0; d<range.ndim; d++) idx0[d] = idx[d]-1;
    kers->l2g[keri](num_cells, idx0, globalidx);

    // Apply the RHS source stencil. It's mostly the mass matrix times a
    // modal-to-nodal operator times the source, modified by BCs in skin cells.
    kers->solker(x_global, globalidx, local_d);
  }
}

void 
gkyl_fem_poisson_set_rhs_cu(gkyl_fem_poisson *up, struct gkyl_array *rhsin)
{
  gkyl_culinsolver_clear_rhs(up->prob_cu, 0);
  double *rhs_cu = gkyl_culinsolver_get_rhs_ptr(up->prob_cu, 0);

  gkyl_fem_poisson_set_rhs_kernel<<<rhsin->nblocks, rhsin->nthreads>>>(up->epsilon->on_dev, up->isvareps, up->dx_cu, 
    rhs_cu, rhsin->on_dev, *up->solve_range, up->bcvals_cu, up->kernels_cu); 

  gkyl_fem_poisson_apply_target_bias_kernel<<<rhsin->nblocks, rhsin->nthreads>>>(rhs_cu, *up->solve_range, 
   up->is_z_edge, up->idxLCFS_m, up->target_corner_bias, up->kernels_cu);
}	

void
gkyl_fem_poisson_solve_cu(gkyl_fem_poisson *up, struct gkyl_array *phiout)
{
  // Do linear solve with cusolver.
  gkyl_culinsolver_solve(up->prob_cu);
  double *x_cu = gkyl_culinsolver_get_sol_ptr(up->prob_cu, 0);

  gkyl_fem_poisson_get_sol_kernel<<<phiout->nblocks, phiout->nthreads>>>(phiout->on_dev, x_cu, *up->solve_range, up->kernels_cu); 
}

