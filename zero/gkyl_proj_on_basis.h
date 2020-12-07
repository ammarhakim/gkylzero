#pragma once

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Type of function to project
typedef void (*evalf_t)(double t, const double *xn, double *fout);

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
 * @return New updater pointer. Free using
 */
gkyl_proj_on_basis* gkyl_proj_on_basis_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_quad, int num_ret_vals, evalf_t eval);

/**
 * Compute projection on basis.
 *
 * @param pob Project on basis updater to run
 * @param tm Time at which projection must be computed
 * @param up_rng Range on which to run projection
 * @param arr_rng Range on which output array is defined
 * @param out Output array
 */
void gkyl_proj_on_basis_advance(const gkyl_proj_on_basis* pob,
  double tm,
  const struct gkyl_range *up_rng, const struct gkyl_range *arr_rng,
  struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_on_basis_release(gkyl_proj_on_basis* pob);
