#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type.
typedef struct gkyl_dg_fv_interpolate gkyl_dg_fv_interpolate;

/**
 * Create a new updater to interpolate a field from cell-average finite
 * volume representation to p=1 modal DG representation on a 2x coarsened grid
 * and vice-versa. 
 * 
*/
struct gkyl_dg_fv_interpolate*
gkyl_dg_fv_interpolate_new(struct gkyl_basis *basis, 
  struct gkyl_rect_grid *grid_dg, struct gkyl_rect_grid *grid_fv);

/**
 * Interpolate from finite volume to p=1 modal DG.
 * 
 */
void
gkyl_dg_fv_interpolate_fv_to_modal_advance(gkyl_dg_fv_interpolate *up, 
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *fv_in, int comp_in,
  struct gkyl_array *dg_out, int comp_out);

/**
 * Interpolate from finite volume to p=1 modal DG using recovery-like definition.
 * 
 */
void
gkyl_dg_fv_interpolate_fv_to_modal_recov_advance(gkyl_dg_fv_interpolate *up, 
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *fv_in, int comp_in,
  struct gkyl_array *dg_out, int comp_out);

/**
 * Interpolate from p=1 modal DG to finite volume.
 * 
 */
void
gkyl_dg_fv_interpolate_modal_to_fv_advance(gkyl_dg_fv_interpolate *up, 
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *dg_in, int comp_in,
  struct gkyl_array *fv_out, int comp_out);

void
gkyl_dg_fv_interpolate_release(gkyl_dg_fv_interpolate *up);