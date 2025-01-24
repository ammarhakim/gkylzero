// Private header for fem_parproj_multib updater.
#pragma once
#include <gkyl_fem_parproj_priv.h>
#include <gkyl_basis.h>
#include <gkyl_superlu_ops.h>
#ifdef GKYL_HAVE_CUDA
#include <gkyl_culinsolver_ops.h>
#endif

struct gkyl_fem_parproj_multib {
  int ndim; // Grid's number of dimensions.
  int num_basis; // Number of basis functions.
  enum gkyl_basis_type basis_type; // Type of DG basis.
  int poly_order; // Polynomial order of basis function.
  int pardir; // Parallel (z) direction.
  int *parnum_cells; // Number of cells in parallel (z) direction.
  bool isdirichlet; // =true if parallel direction has periodic BCs.
  struct gkyl_array *weight_rhs; // The RHS weight.

  int num_blocks; // Number of blocks connected along the magnetic field.
  const struct gkyl_range *mbz_range; // Multiblock range along the magnetic field.
  struct gkyl_range *solve_range; // Range to perform the projection in.
  struct gkyl_range perp_range2d; // 2D range of perpendicular cells.
  struct gkyl_range *par_range1d; // 1D range of parallel cells.
  struct gkyl_range_iter perp_iter2d; // 2D iterator.
  struct gkyl_range_iter par_iter1d; // 1D iterator.

  int numnodes_perp; // Number of local perp nodes.
  long *numnodes_global; // Total number of nodes in linear problem.
  long *globalidx; // Global indices (in linear problem) of each node within a cell.

  struct gkyl_superlu_prob* prob; // A SuperLU linear problem.
  struct gkyl_array *brhs; // RHS vector/matrix of linear problem (flat).
  struct gkyl_fem_parproj_kernels *kernels;
#ifdef GKYL_HAVE_CUDA
  struct gkyl_culinsolver_prob* prob_cu; // A CUDA linear problem.
  struct gkyl_array *brhs_cu; // Device RHS vector/matrix of linear problem (flat).
  struct gkyl_fem_parproj_kernels *kernels_cu;
#endif

  bool use_gpu; // Whether to run on the GPU.
};

#ifdef GKYL_HAVE_CUDA
/**
 * Assign the right-side vector with the discontinuous (DG) source field
 * on the NVIDIA GPU.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 * @param phibc Potential to use for Dirichlet BCs (only use ghost cells).
 */
void gkyl_fem_parproj_multib_set_rhs_cu(struct gkyl_fem_parproj_multib *up, const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Solve the linear problem
 * on the NVIDIA GPU.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_multib_solve_cu(struct gkyl_fem_parproj_multib* up, struct gkyl_array *phiout);
#endif
