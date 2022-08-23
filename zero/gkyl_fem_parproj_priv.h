// Private header for fem_parproj updater.
#pragma once
#include <gkyl_fem_parproj_kernels.h>
#include <gkyl_basis.h>
#include <gkyl_superlu_ops.h>
#ifdef GKYL_HAVE_CUDA
#include <gkyl_cusolver_ops.h>
#endif

// Function pointer type for local-to-global mapping.
typedef void (*local2global_t)(const int numCellsPar, const int parIdx,
  long *globalIdxs);

typedef struct { local2global_t kernels[3]; } local2global_kern_list;  // For use in kernel tables.

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

// Function pointer type for lhs kernels.
typedef void (*lhsstencil_t)(const double *weight, const long *globalIdxs, gkyl_mat_triples *tri);

typedef struct { lhsstencil_t kernels[3]; } lhsstencil_kern_list;  // For use in kernel tables.

// Serendipity unweighted lhs kernels.
static const lhsstencil_kern_list ser_lhsstencil_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_lhs_stencil_1x_ser_p1, fem_parproj_lhs_stencil_1x_ser_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_lhs_stencil_3x_ser_p1, fem_parproj_lhs_stencil_3x_ser_p2 }
};

// Tensor product unweighted lhs kernels.
static const lhsstencil_kern_list ten_lhsstencil_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_lhs_stencil_1x_tensor_p1, fem_parproj_lhs_stencil_1x_tensor_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_lhs_stencil_3x_tensor_p1, fem_parproj_lhs_stencil_3x_tensor_p2 }
};

// Serendipity weighted lhs kernels.
static const lhsstencil_kern_list ser_wlhsstencil_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_weighted_lhs_stencil_1x_ser_p1, fem_parproj_weighted_lhs_stencil_1x_ser_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_weighted_lhs_stencil_3x_ser_p1, fem_parproj_weighted_lhs_stencil_3x_ser_p2 }
};

// Tensor product weighted lhs kernels.
static const lhsstencil_kern_list ten_wlhsstencil_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_weighted_lhs_stencil_1x_tensor_p1, fem_parproj_weighted_lhs_stencil_1x_tensor_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_weighted_lhs_stencil_3x_tensor_p1, fem_parproj_weighted_lhs_stencil_3x_tensor_p2 }
};

// Function pointer type for rhs source kernels.
typedef void (*srcstencil_t)(const double *rho, long nodeOff, const long *globalIdxs,
  double *bsrc);

typedef struct { srcstencil_t kernels[3]; } srcstencil_kern_list;  // For use in kernel tables.

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

// Function pointer type for kernels that convert the solution from nodal to
// modal.
typedef void (*solstencil_t)(const double *sol_nodal_global, long nodeOff,
  const long *globalIdxs, double *sol_modal_local);

typedef struct { solstencil_t kernels[3]; } solstencil_kern_list;  // For use in kernel tables.

// Serendipity sol kernels.
GKYL_CU_D
static const solstencil_kern_list ser_solstencil_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_sol_stencil_1x_ser_p1, fem_parproj_sol_stencil_1x_ser_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_sol_stencil_3x_ser_p1, fem_parproj_sol_stencil_3x_ser_p2 }
};

// Tensor product sol kernels.
GKYL_CU_D
static const solstencil_kern_list ten_solstencil_list[] = {
  { NULL, NULL, NULL }, // No 0D basis functions
  { NULL, fem_parproj_sol_stencil_1x_tensor_p1, fem_parproj_sol_stencil_1x_tensor_p2 },
  { NULL, NULL, NULL }, // No 2D basis functions
  { NULL, fem_parproj_sol_stencil_3x_tensor_p1, fem_parproj_sol_stencil_3x_tensor_p2 }
};

// Struct containing pointers to the various kernels. Needed to create a similar struct on the GPU.
struct gkyl_fem_parproj_kernels {
  local2global_t l2g;  // Pointer to local-to-global kernel.

  lhsstencil_t lhsker;  // Weighted LHS kernel.

  srcstencil_t srcker;  // RHS source kernel.

  solstencil_t solker;  // Kernel that takes the solution and converts it to modal.
};

struct gkyl_fem_parproj {
  struct gkyl_rect_grid grid;
  int ndim; // grid's number of dimensions.
  int num_basis; // number of basis functions.
  enum gkyl_basis_type basis_type;
  int poly_order;
  int pardir; // parallel (z) direction.
  int parnum_cells; // number of cells in parallel (z) direction.
  bool isperiodic; // =true if parallel direction is periodic.

  struct gkyl_range local_range, local_range_ext;
  struct gkyl_range solve_range;
  struct gkyl_range_iter solve_iter;
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
#ifdef GKYL_HAVE_CUDA
  struct gkyl_cusolver_prob* prob_cu;
  struct gkyl_array *brhs_cu;
#endif

  long *globalidx;

  struct gkyl_fem_parproj_kernels *kernels;
  struct gkyl_fem_parproj_kernels *kernels_cu;
  bool use_gpu;  
};

void
fem_parproj_choose_kernels_cu(const struct gkyl_basis* basis, bool isperiodic, struct gkyl_fem_parproj_kernels *kers);

GKYL_CU_D
static local2global_t
fem_parproj_choose_local2global_kernel(int dim, int basis_type, int poly_order)
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
static lhsstencil_t
fem_parproj_choose_lhs_kernel(int dim, int basis_type, int poly_order, bool isweighted)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return isweighted? ser_wlhsstencil_list[dim].kernels[poly_order] : ser_lhsstencil_list[dim].kernels[poly_order]; 
    case GKYL_BASIS_MODAL_TENSOR:
      return isweighted? ten_wlhsstencil_list[dim].kernels[poly_order] : ten_lhsstencil_list[dim].kernels[poly_order];
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static srcstencil_t
fem_parproj_choose_srcstencil_kernel(int dim, int basis_type, int poly_order)
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

GKYL_CU_D
static solstencil_t
fem_parproj_choose_solstencil_kernel(int dim, int basis_type, int poly_order)
{
  switch (basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_solstencil_list[dim].kernels[poly_order];
    case GKYL_BASIS_MODAL_TENSOR:
      return ten_solstencil_list[dim].kernels[poly_order];
    default:
      assert(false);
      break;
  }
}
