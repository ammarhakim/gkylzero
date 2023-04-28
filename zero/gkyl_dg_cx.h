#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Identifiers for different charge exchange types
enum gkyl_dg_cx_type
{
  GKYL_H, // Hydrogen plasma
  GKYL_D, // Deuterium plasma
  GKYL_NE, // Neon plasma
  GKYL_HE, // Helium plasma
};

// Object type
typedef struct gkyl_dg_cx gkyl_dg_cx;

/**
 * Create new updater to calculate charge exchange reaction rate
 * @param cbasis Configuration-space basis-functions
 * @param type_ion Enum for type of ion for charge exchange (support H^+, D^+, Ar^+, Ne^+)
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_cx* gkyl_dg_cx_new(const struct gkyl_basis* cbasis, 
  enum gkyl_dg_cx_type type_ion, 
  bool use_gpu);

/**
 * Compute CX reaction rate coefficient for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param cx charge exchange object.
 * @param update_rng Update range (Configuration space)
 * @param phase_rng Phase-space range (for indexing cflrate for stable timestep)
 * @param m0_neut Input neutral density 
 * @param u_neut Input neutral fluid flow
 * @param u_ion Input ion fluid flow
 * @param vth_sq_neut Input neutral vth^2 = T/m
 * @param vth_sq_ion Input ion vth^2 = T/m
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param coef_cx Output reaction rate coefficient
 */

void gkyl_dg_cx_react_rate(const struct gkyl_dg_cx *cx,
  const struct gkyl_range *update_rng, const struct gkyl_range *phase_rng, 
  const struct gkyl_array *n_neut, const struct gkyl_array *u_neut, const struct gkyl_array *vth_sq_neut,
  const struct gkyl_array *u_ion, const struct gkyl_array *vth_sq_ion,
  struct gkyl_array *cflrate, struct gkyl_array *coef_cx);

/**
 * Delete updater.
 *
 * @param cx Updater to delete.
 */
void gkyl_dg_cx_release(gkyl_dg_cx *cx);
