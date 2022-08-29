#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_prim_lbo_cross_calc gkyl_prim_lbo_cross_calc;

/**
 * Create new updater to compute cross-primitive moments of distribution
 * function. Free using gkyl_prim_vlasov_cross_calc_release.
 *
 * @param grid Grid object
 * @param prim Pointer to primitive moment type object
 * @return New updater pointer.
 */
gkyl_prim_lbo_cross_calc* gkyl_prim_lbo_cross_calc_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim);

gkyl_prim_lbo_cross_calc* gkyl_prim_lbo_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim);

/**
 * Compute cross-primitive moments of distribution function. The conf_rng
 * MUST be a sub-range of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Primitive moment calculator updater to run
 * @param cbasis_rng Config-space basis functions
 * @param conf_rng Config-space range
 * @param greene Greene's factor
 * @param self_m Mass of the species
 * @param self_u Drift velocity of the species
 * @param self_vtsq Thermal velocity of the species
 * @param cross_m Mass of the colliding species
 * @param cross_u Drift velocity of the colliding species
 * @param cross_vtsq Thermal velocity of the colliding species
 * @param moms Moments of distribution function (Zeroth, First, and Second)
 * @param boundary_corrections Momentum and Energy boundary corrections
 * @param u_out Output drift velocity primitive moment array
 * @param vtsq_out Output thermal velocity primitive moment array
 */
void gkyl_prim_lbo_cross_calc_advance(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_basis cbasis, const struct gkyl_range *conf_rng,
  const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  double cross_m, const struct gkyl_array *cross_u, const struct gkyl_array *cross_vtsq, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections, 
  struct gkyl_array *u_out, struct gkyl_array *vtsq_out);

void gkyl_prim_lbo_cross_calc_advance_cu(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_basis cbasis, const struct gkyl_range *conf_rng,
  const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_u, const struct gkyl_array *self_vtsq,
  double cross_m, const struct gkyl_array *cross_u, const struct gkyl_array *cross_vtsq, 
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections, 
  struct gkyl_array *u_out, struct gkyl_array *vtsq_out);

/**
 * Delete pointer to primitive moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_prim_lbo_cross_calc_release(gkyl_prim_lbo_cross_calc* calc);

/**
 * Return pointer to primitive moment type structure.

 * @param calc Updater pointer.
 */
const struct gkyl_prim_lbo_type* gkyl_prim_lbo_cross_calc_get_prim(gkyl_prim_lbo_cross_calc* calc);

// "derived" class constructors
gkyl_prim_lbo_cross_calc* 
gkyl_prim_lbo_vlasov_cross_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_cross_calc* 
gkyl_prim_lbo_gyrokinetic_cross_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_cross_calc* 
gkyl_prim_lbo_vlasov_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_cross_calc* 
gkyl_prim_lbo_gyrokinetic_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);
