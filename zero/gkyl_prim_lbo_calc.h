#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_prim_lbo.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_prim_lbo_calc gkyl_prim_lbo_calc;

struct gkyl_prim_lbo_calc {
  struct gkyl_rect_grid grid;
  struct gkyl_prim_lbo *prim;
};

/**
 * Create new updater to compute primitive moments of distribution
 * function. Free using gkyl_prim_vlasov_calc_release.
 *
 * @param grid Grid object
 * @param prim Pointer to primitive moment type object
 * @return New updater pointer.
 */
gkyl_prim_lbo_calc* gkyl_prim_lbo_calc_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo *prim);

/**
 * Compute primitive moments of distribution function. The phase_rng and conf_rng
 * MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Primitive moment calculator updater to run
 * @param conf_rng Config-space basis functions
 * @param conf_rng Config-space range
 * @param m0 Zeroeth moment of distribution function
 * @param m1 First moment of distribution function
 * @param m2 Second moment of distribution function
 * @param cM Momentum boundary correction
 * @param cE Energy boundary correction
 * @param uout Output drift velocity primitive moment array
 * @param vtSqout Output thermal velocity primitive moment array
 */
void gkyl_prim_lbo_calc_advance(gkyl_prim_lbo_calc* calc, const struct gkyl_basis cbasis,
  const struct gkyl_range conf_rng, const struct gkyl_array *m0, const struct gkyl_array *m1,
  const struct gkyl_array *m2, const struct gkyl_array *cM, const struct gkyl_array *cE,
  struct gkyl_array *uout, struct gkyl_array *vtSqout);

/**
 * Delete pointer to primitive moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_prim_lbo_calc_release(gkyl_prim_lbo_calc* calc);
