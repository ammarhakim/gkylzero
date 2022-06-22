// Private header for fem_parproj updater.
#pragma once
#include <gkyl_fem_parproj_kernels.h>
#ifdef GKYL_HAVE_CUDA
#include <gkyl_cusolver_ops.h>
#endif

// Function pointer type for local-to-global mapping.
typedef void (*local2global_t)(const int numCellsPar, const int parIdx,
  long *globalIdxs);

// For use in kernel tables.
typedef struct { local2global_t kernels[3]; } local2global_kern_list;

// Serendipity local-to-global kernels.
GKYL_CU_D
static const local2global_kern_list ser_loc2glob_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_local_to_global_1x_ser_p1, fem_parproj_local_to_global_1x_ser_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_local_to_global_3x_ser_p1, fem_parproj_local_to_global_3x_ser_p2 }
};

// Tensor product local-to-global kernels.
GKYL_CU_D
static const local2global_kern_list ten_loc2glob_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_local_to_global_1x_tensor_p1, fem_parproj_local_to_global_1x_tensor_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_local_to_global_3x_tensor_p1, fem_parproj_local_to_global_3x_tensor_p2 }
};

// Function pointer type for rhs source kernels.
typedef void (*srcstencil_t)(const double *rho, long nodeOff, const long *globalIdxs,
  double *bsrc);

// For use in kernel tables.
typedef struct { srcstencil_t kernels[3]; } srcstencil_kern_list;

// Serendipity src kernels.
GKYL_CU_D
static const srcstencil_kern_list ser_srcstencil_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_src_stencil_1x_ser_p1, fem_parproj_src_stencil_1x_ser_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_src_stencil_3x_ser_p1, fem_parproj_src_stencil_3x_ser_p2 }
};

// Tensor product src kernels.
GKYL_CU_D
static const srcstencil_kern_list ten_srcstencil_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_src_stencil_1x_tensor_p1, fem_parproj_src_stencil_1x_tensor_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_src_stencil_3x_tensor_p1, fem_parproj_src_stencil_3x_tensor_p2 }
};

// Struct containing pointers to the various kernels. Needed to create a similar struct on the GPU.
struct gkyl_fem_parproj_kernels {
  local2global_t l2g;  // Pointer to local-to-global kernel.

  srcstencil_t srcker;  // RHS source kernels.

//  solstencil_t solker;  // Kernel that takes the solution and converts it to modal.
};

struct gkyl_fem_parproj {
  void *ctx; // evaluation context.
  struct gkyl_rect_grid grid;
  int ndim; // grid's number of dimensions.
  int num_basis; // number of basis functions.
  enum gkyl_basis_type basis_type;
  int poly_order;
  int pardir; // parallel (z) direction.
  int parnum_cells; // number of cells in parallel (z) direction.
  bool isperiodic; // =true if parallel direction is periodic.

  struct gkyl_range local_range, local_range_ext;
  struct gkyl_range perp_range; // range of perpendicular cells.
  struct gkyl_range_iter perp_iter;
  struct gkyl_range par_range; // range of parallel cells.
  struct gkyl_range_iter par_iter;
  struct gkyl_range par_range1d; // 1D range of parallel cells.
  struct gkyl_range_iter par_iter1d;
  struct gkyl_range perp_range2d; // 2D range of perpendicular cells.
  struct gkyl_range_iter perp_iter2d;

  int numnodes_local;
  long numnodes_global;

  struct gkyl_array *brhs;

  struct gkyl_superlu_prob* prob;
  struct gkyl_mat *local_mass; // local mass matrix.
  struct gkyl_mat *local_nodtomod; // local nodal-to-modal matrix.

  long *globalidx;

  struct gkyl_fem_parproj_kernels *kernels;
  struct gkyl_fem_parproj_kernels *kernels_cu;
  bool use_gpu;  
};

void
choose_kernels_cu(const struct gkyl_basis* basis, bool isperiodic, struct gkyl_fem_parproj_kernels *kers);

GKYL_CU_D
static local2global_t
choose_local2global_kernel(const int dim, const int basis_type, const int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_loc2glob_list[dim].kernels[poly_order];
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_loc2glob_list[dim].kernels[poly_order];
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static srcstencil_t
choose_srcstencil_kernel(const int dim, const int basis_type, const int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_srcstencil_list[dim].kernels[poly_order];
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_srcstencil_list[dim].kernels[poly_order];
    default:
      assert(false);
      break;
  }
}
