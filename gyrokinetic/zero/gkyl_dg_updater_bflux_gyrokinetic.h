#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_updater_gyrokinetic.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_bflux_gyrokinetic gkyl_dg_updater_bflux_gyrokinetic;

// return type for drag and diffusion timers
struct gkyl_dg_updater_bflux_gyrokinetic_tm {
  double bflux_tm; // time for bflux updates
};

/**
 * Create new updater to compute the boundary flux contributions from
 * gyrokinetic equations. These are DG surface terms using hyper dg.
 *
 * @param grid Grid object.
 * @param cdim Number of configuration space dimensions.
 * @param gyrokinetic objects which performs the gyrokinetic DG update.
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * 
 * @return New boundary flux updater object.
 */
struct gkyl_dg_updater_bflux_gyrokinetic* gkyl_dg_updater_bflux_gyrokinetic_new(const struct gkyl_rect_grid *grid, 
  int cdim, const gkyl_dg_updater_gyrokinetic *gyrokinetic, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Boundary flux updater.
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater.
 * @param rhs RHS output.
 */
void gkyl_dg_updater_bflux_gyrokinetic_advance(struct gkyl_dg_updater_bflux_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs);

void gkyl_dg_updater_bflux_gyrokinetic_advance_cu(struct gkyl_dg_updater_bflux_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent computing boundary fluxes.
 *
 * @param up Boundary flux updater.
 * @return timers
 */
struct gkyl_dg_updater_bflux_gyrokinetic_tm gkyl_dg_updater_bflux_gyrokinetic_get_tm(const struct gkyl_dg_updater_bflux_gyrokinetic *up);

/**
 * Delete updater.
 *
 * @param up Boundary flux updater.
 */
void gkyl_dg_updater_bflux_gyrokinetic_release(struct gkyl_dg_updater_bflux_gyrokinetic *up);
