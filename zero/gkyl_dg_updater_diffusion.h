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
 * @param grid Grid object.
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param diffusion_id Diffusion type (constant/varying & fluid/vlasov/etc).
 * @param diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @param conf_range Configuration space range (to index diff coefficient).
 * @param use_gpu Whether to run on host or device.
 * @return New diff updater object
 */
struct gkyl_dg_updater_diffusion* gkyl_dg_updater_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  enum gkyl_diffusion_id diffusion_id, bool *diff_in_dir, int diff_order,
  const struct gkyl_range *conf_range, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up diffusion updater object.
 * @param update_rng Range on which to compute.
 * @param coeff Diffusion coefficient/tensor.
 * @param fIn Input to updater.
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_diffusion_advance(struct gkyl_dg_updater_diffusion *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *coeff,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in diffusion terms
 *
 * @param diffusion Updater object
 * @return timers
 */
struct gkyl_dg_updater_diffusion_tm gkyl_dg_updater_diffusion_get_tm(const struct gkyl_dg_updater_diffusion *up);

/**
 * Delete updater.
 *
 * @param diffusion Updater to delete.
 */
void gkyl_dg_updater_diffusion_release(struct gkyl_dg_updater_diffusion *up);
