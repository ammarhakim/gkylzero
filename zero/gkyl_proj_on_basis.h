#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_real_type.h>
#include <gkyl_rect_grid.h>

/**
 * Type of function to project.
 *
 * @param t Time to evaluate function
 * @param xn Coordinates for evaluation
 * @param fout Output vector of 'num_ret_vals'
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*evalf_t)(gkyl_real t, const gkyl_real *xn, gkyl_real *fout, void *ctx);

// Object type
typedef struct gkyl_proj_on_basis gkyl_proj_on_basis;

/**
 * Create new updater to project function on basis functions on a
 * grid. Free using gkyl_proj_on_basis_release method.
 *
 * @param grid Grid object
 * @param basis Basis functions to project on
 * @param num_quad Number of quadrature nodes
 * @param num_ret_vals Number of values 'eval' sets
 * @param eval Function to project.
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
gkyl_proj_on_basis* gkyl_proj_on_basis_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx);

/**
 * Compute projection on basis. The update_rng MUST be a sub-range of
 * the range on which the array is defined. That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param pob Project on basis updater to run
 * @param tm Time at which projection must be computed
 * @param update_rng Range on which to run projection.
 * @param out Output array
 */
void gkyl_proj_on_basis_advance(const gkyl_proj_on_basis *pob,
  gkyl_real tm, const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_on_basis_release(gkyl_proj_on_basis* pob);
