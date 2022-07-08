#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_diffusion gkyl_dg_updater_diffusion;

// return type for drag and diffusion timers
struct gkyl_dg_updater_diffusion_tm {
  double diffusion_tm; // time for diffusion updates
};

/**
 * Create new updater to update diffusion equations using hyper dg or hyper dg gen stencil.
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param conf_range Config space range
 * @return New diff updater object
 */
gkyl_dg_updater_diffusion* gkyl_dg_updater_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  enum gkyl_diffusion_id diffusion_id, bool use_gpu);

gkyl_dg_updater_diffusion* gkyl_dg_updater_diffusion_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  enum gkyl_diffusion_id diffusion_id);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param diff diffusion updater object
 * @param update_rng Range on which to compute.
 * @param D Diffusion tensor
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_diffusion_advance(gkyl_dg_updater_diffusion *diff, enum gkyl_diffusion_id diffusion_id,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *D, 
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

void gkyl_dg_updater_diffusion_advance_cu(gkyl_dg_updater_diffusion *diff, enum gkyl_diffusion_id diffusion_id,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *D, 
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in diffusion terms
 *
 * @param diffusion Updater object
 * @return timers
 */
struct gkyl_dg_updater_diffusion_tm gkyl_dg_updater_diffusion_get_tm(const gkyl_dg_updater_diffusion *diff);

/**
 * Delete updater.
 *
 * @param diffusion Updater to delete.
 */
void gkyl_dg_updater_diffusion_release(gkyl_dg_updater_diffusion* diff);
