#pragma once

#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_positivity_shift_gyrokinetic gkyl_positivity_shift_gyrokinetic;

/**
 * Create a new upater which shifts f up to a floor wherever f<0, and computes
 * contribution to the integrated moments of the shift in f in each cell (the
 * final reduction needs to be performed outside of this updater, as it may
 * involve a communicator).
 *
 * @param cdim Configuration-space dimensions.
 * @param pbasis Phase-space basis.
 * @param pgrid Phase-space grid.
 * @param mass Mass of species.
 * @param gk_geom Geometry object.
 * @param use_gpu bool to determine if on GPU.
 * @return New positivity shift updater pointer.
 */
struct gkyl_positivity_shift_gyrokinetic*
gkyl_positivity_shift_gyrokinetic_new(int cdim, struct gkyl_basis pbasis, struct gkyl_rect_grid pgrid,
  double mass, const struct gk_geometry *gk_geom, bool use_gpu);

/**
 * Run the positivity shift updater in the indicated range.
 *
 * @param up Positivity shift updater.
 * @param phase_rng Phase-space range.
 * @param conf_rng Config-space range.
 * @param distf Distribution function array.
 * @param mom Output moment array.
 */
void
gkyl_positivity_shift_gyrokinetic_advance(gkyl_positivity_shift_gyrokinetic* up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  struct gkyl_array *GKYL_RESTRICT distf, struct gkyl_array *GKYL_RESTRICT mom);

/**
 * Release the memory associated with this positivity shift updater.
 *
 * @param up Positivity shift updater.
 */
void
gkyl_positivity_shift_gyrokinetic_release(gkyl_positivity_shift_gyrokinetic* up);
