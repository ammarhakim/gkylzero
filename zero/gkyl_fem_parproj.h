#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_mat.h>
#include <gkyl_mat_triples.h>
#include <gkyl_superlu_ops.h>

// Object type
typedef struct gkyl_fem_parproj gkyl_fem_parproj;

/**
 * Create new updater to project a DG field onto the FEM (nodal) basis
 * in order to make the field continuous or, thanks to the option to pass
 * a multiplicative weight, solve 1D algebraic equations in which the output
 * field is continuous (but the input may not be). That is, we solve
 *    wgt*phi_{fem} \doteq rho_{dg}
 * where wgt is the weight field, phi_{fem} is the (continuous field)
 * we wish to compute, rho_{dg} is the (discontinuous) input source field,
 * and \doteq implies weak equality with respect to the FEM basis.
 * Free using gkyl_fem_parproj_release method.
 *
 * @param grid Grid object
 * @param basis Basis functions of the DG field.
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
gkyl_fem_parproj* gkyl_fem_parproj_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  void *ctx);

/**
 * Set the multiplicative weight.
 *
 * @param up FEM project updater to run.
 * @param tm Time at which projection must be computed
 * @param inw Input array weight.
 */
void gkyl_fem_parproj_set_weight(const gkyl_fem_parproj *up,
  double tm, const struct gkyl_array *inw);

/**
 * Begin assembling the right-side source vector and, if necessary,
 * the stiffness matrix.
 * In parallel simulations there will be an option to use
 * non-blocking MPI here and later checking that assembly is complete.
 *
 * @param up FEM project updater to run.
 * @param tm Time at which projection must be computed
 * @param src Input source field.
 */
void gkyl_fem_parproj_begin_assembly(const gkyl_fem_parproj *up,
  double tm, const struct gkyl_array *src);

/**
 * Assign the right-side vector with the discontinuous (DG) source field.
 *
 * @param up FEM project updater to run.
 * @param rhsin DG field to set as RHS source.
 */
void gkyl_fem_parproj_set_rhs(gkyl_fem_parproj* up, const struct gkyl_array *rhsin);

/**
 * Solve the linear problem.
 *
 * @param up FEM project updater to run.
 */
void gkyl_fem_parproj_solve(gkyl_fem_parproj* up, struct gkyl_array *phiout);

/**
 * Compute the projection onto the FEM basis. 
 *
 * @param up FEM project updater to run
 * @param tm Time at which projection must be computed
 * @param out Output array
 */
void gkyl_fem_parproj_advance(gkyl_fem_parproj *up,
  double tm, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_parproj_release(gkyl_fem_parproj *up);
