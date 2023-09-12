#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Identifiers for different charge exchange types
enum gkyl_dg_cx_type
{
  GKYL_CX_H, // Hydrogen plasma
  GKYL_CX_D, // Deuterium plasma
  GKYL_CX_NE, // Neon plasma
  GKYL_CX_HE, // Helium plasma
};

// Object type
typedef struct gkyl_dg_cx gkyl_dg_cx;

/**
 * Create new updater to calculate charge exchange reaction rate
 * @param grid Phase space grid
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param conf_rng Configuration-space range
 * @param phase_rng Phase
 * @param mass_ion Ion mass
 * @param is_gk Boolen for whether gk or vlasov
 * @param type_ion Enum for type of ion for charge exchange (support H^+, D^+, Ar^+, Ne^+)
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_cx* gkyl_dg_cx_new(const struct gkyl_rect_grid *grid,
  struct gkyl_basis *cbasis, struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  double mass_ion, enum gkyl_dg_cx_type type_ion, 
  bool is_gk, bool use_gpu);

/**
 * Compute CX reaction rate coefficient for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param cx charge exchange object.
 * @param moms_neut Input neutral moments 
 * @param moms_ion Input ion moments
 * @param b_hat Input unit B vector
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param coef_cx Output reaction rate coefficient
 */

void gkyl_dg_cx_react_rate(const struct gkyl_dg_cx *cx, const struct gkyl_array *moms_ion,
  const struct gkyl_array *moms_neut, const struct gkyl_array *b_hat,
  struct gkyl_array *cflrate, struct gkyl_array *coef_cx);

/**
 * Delete updater.
 *
 * @param cx Updater to delete.
 */
void gkyl_dg_cx_release(gkyl_dg_cx *cx);
