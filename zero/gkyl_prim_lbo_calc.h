#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_prim_lbo_calc gkyl_prim_lbo_calc;

/**
 * Create new updater to compute primitive moments of distribution
 * function. Free using gkyl_prim_vlasov_calc_release.
 *
 * @param grid Grid object
 * @param prim Pointer to primitive moment type object
 * @return New updater pointer.
 */
gkyl_prim_lbo_calc* gkyl_prim_lbo_calc_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim);

gkyl_prim_lbo_calc* gkyl_prim_lbo_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim);

/**
 * Compute primitive moments of distribution function. The phase_rng and conf_rng
 * MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Primitive moment calculator updater to run
 * @param cbasis_rng Config-space basis functions
 * @param conf_rng Config-space range
 * @param moms Moments of distribution function (Zeroth, First, and Second)
 * @param boundary_corrections Momentum and Energy boundary corrections
 * @param uout Output drift velocity primitive moment array
 * @param vtSqout Output thermal velocity primitive moment array
 */
void gkyl_prim_lbo_calc_advance(gkyl_prim_lbo_calc* calc, struct gkyl_basis cbasis,
  struct gkyl_range conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections,
  struct gkyl_array *uout, struct gkyl_array *vtSqout);

void gkyl_prim_lbo_calc_advance_cu(gkyl_prim_lbo_calc* calc, struct gkyl_basis cbasis,
  struct gkyl_range conf_rng, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections,
  struct gkyl_array* uout, struct gkyl_array* vtSqout);

/**
 * Delete pointer to primitive moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_prim_lbo_calc_release(gkyl_prim_lbo_calc* calc);
