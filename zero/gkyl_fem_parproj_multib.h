#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_array_ops.h>
#include <gkyl_mat.h>
#include <gkyl_mat_triples.h>
#include <gkyl_superlu_ops.h>
#include <gkyl_fem_parproj.h>

// Object type
typedef struct gkyl_fem_parproj_multib gkyl_fem_parproj_multib;

/**
 * Create new updater to project a multiblock DG field onto the weighted
 * FEM (nodal) basis in order to make the field continuous. That is, we solve
 *    wgtL*phi_{fem} \doteq wgtR*rho_{dg}
 * where wgt is the weight field, phi_{fem} is the (continuous field)
 * we wish to compute, rho_{dg} is the (discontinuous) input source field,
 * and \doteq implies weak equality with respect to the FEM basis.
 * Free using gkyl_fem_parproj_multib_release method.
 *
 * @param solve_range Range in which to perform the projection operation.
 * @param basis Basis functions of the DG field.
 * @param bctype Type of boundary condition (see gkyl_fem_parproj_bc_type).
 * @param weight_left Weight on left-side of the operator (time-independent,
                      global across blocks connected in the parallel direction).
 * @param weight_right Weight on right-side of the operator (time-independent,
                       global across blocks connected in the parallel direction).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_fem_parproj_multib* gkyl_fem_parproj_multib_new(int num_blocks, const struct gkyl_range *mbz_range,
  struct gkyl_range *solve_ranges, const struct gkyl_basis *basis, enum gkyl_fem_parproj_bc_type bctype,
  const struct gkyl_array *weight_left, const struct gkyl_array *weight_right,  bool use_gpu);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 * @param phibc Potential to use for Dirichlet BCs (only use ghost cells).
 */
void gkyl_fem_parproj_multib_set_rhs(struct gkyl_fem_parproj_multib* up,
  const struct gkyl_array *rhsin, const struct gkyl_array *phibc);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_multib_solve(struct gkyl_fem_parproj_multib* up, struct gkyl_array *phiout);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_parproj_multib_release(struct gkyl_fem_parproj_multib *up);
