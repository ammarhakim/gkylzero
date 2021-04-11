#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

/**
 * Type of function to apply BC
 *
 * @param t Time at which BC is applied
 * @param dir Direction in which BC is applied
 * @param skin Pointer to data in skin-cell
 * @param ghost Pointer to data in ghost-cell
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*rect_bc_func_t)(double t, int dir, const double *skin, double *ghost, void *ctx);

// Object type
typedef struct gkyl_rect_apply_bc gkyl_rect_apply_bc;

/**
 * Create new updater to apply set of boundary conditions.
 *
 * @param grid Grid object
 * @param dir Direction to apply BC 
 * @param edge Edge to apply BC (0: left-edge; 1: right-edge)
 * @param Boundary condition function to apply
 * @param ctx Context to pass to bcfunc.
 * @return New updater pointer.
 */
gkyl_rect_apply_bc* gkyl_rect_apply_bc_new(const struct gkyl_rect_grid *grid,
  int dir, int edge, rect_bc_func_t bcfunc, void *ctx);

/**
 * Apply boundary condition on specified field.
 *
 * @param bc BC updater to run
 * @param tm Time at which BC is applied
 * @param out Output array
 */
void gkyl_rect_apply_bc_advance(const gkyl_rect_apply_bc *bc, double tm, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param bc Updater to delete.
 */
void gkyl_rect_apply_bc_release(gkyl_rect_apply_bc* bc);
