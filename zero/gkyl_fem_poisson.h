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

// Object type
typedef struct gkyl_fem_poisson gkyl_fem_poisson;

/**
 * Create new updater to solve the Poisson problem
 *   - nabla^2 phi = rho
 * using a FEM to ensure phi is continuous. The input is the
 * DG field rho, which is translated to FEM. The output is the
 * DG field phi, after we've translted the FEM solution to DG.
 * Free using gkyl_fem_poisson_release method.
 *
 * @param grid Grid object
 * @param basis Basis functions of the DG field.
 * @param isdirperiodic boolean array indicating periodic directions.
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
gkyl_fem_poisson* gkyl_fem_poisson_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const bool *isdirperiodic, void *ctx);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_fem_poisson_release(gkyl_fem_poisson *up);
