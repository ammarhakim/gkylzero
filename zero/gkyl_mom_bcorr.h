#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

/** Flags for indicating acting edge of velocity space */
enum gkyl_vel_edge { GKYL_VX_LOWER, GKYL_VY_LOWER, GKYL_VZ_LOWER, GKYL_VX_UPPER, GKYL_VY_UPPER, GKYL_VZ_UPPER };

// Object type
typedef struct gkyl_mom_bcorr gkyl_mom_bcorr;

/**
 * Create new updater to update boundary corrections.
 *
 * @param grid Grid object
 * @param momt Moment type object for boundary correction
 */
gkyl_mom_bcorr* gkyl_mom_bcorr_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

gkyl_mom_bcorr* gkyl_mom_bcorr_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt);

/**
 * Compute boundary correction moments.
 *
 * @param bcorr Boundary correction updater object
 * @param phase_rng Phase space range on which to compute.
 * @param conf_rng Configuration space range on which to compute.
 * @param fIn Input to updater
 * @param out Output
 */
void gkyl_mom_bcorr_advance(gkyl_mom_bcorr *bcorr, struct gkyl_range phase_rng, struct gkyl_range conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *out);

void gkyl_mom_bcorr_advance_cu(gkyl_mom_bcorr *bcorr, struct gkyl_range phase_rng, struct gkyl_range conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *out);
  
/**
 * Delete updater.
 *
 * @param bcorr Updater to delete.
 */
void gkyl_mom_bcorr_release(gkyl_mom_bcorr* up);
