#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Identifiers for different ionization types
enum gkyl_dg_iz_type
{
  GKYL_H, // Hydrogen plasma
  GKYL_AR, // Argon plasma
};

// Object type
typedef struct gkyl_dg_iz gkyl_dg_iz;

/**
 * Create new updater to calculate ionization temperature or reaction rate
 * @param cbasis Configuration-space basis-functions
 * @param elem_charge Elementary charge value
 * @param mass_elc Mass of the electron value
 * @param type_ion Enum for type of ion for ionization (support H^+ and Ar^+)
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_iz* gkyl_dg_iz_new(const struct gkyl_basis* cbasis, 
  double elem_charge, double mass_elc, enum gkyl_dg_iz_type type_ion, 
  bool use_gpu);

/**
 * Compute ionization temperature for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param iz Ionization object.
 * @param update_rng Update range (Configuration space)
 * @param vth_sq_elc Input electron vth^2 = T/m
 * @param vth_sq_iz Output ionization vth^2 = T_iz/m
 */

void gkyl_dg_iz_temp(const struct gkyl_dg_iz *iz,
  const struct gkyl_range *update_rng, const struct gkyl_array *vth_sq_elc,
  struct gkyl_array *vth_sq_iz);

/**
 * Compute reaction rate coefficient for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param iz Ionization object.
 * @param update_rng Update range (Configuration space)
 * @param phase_rng Phase-space range (for indexing cflrate for stable timestep)
 * @param n_neut Input neutral density 
 * @param vth_sq_neut Input neutral vth^2 = T/m
 * @param vth_sq_elc Input electron vth^2 = T/m
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param coef_iz Output reaction rate coefficient
 */

void gkyl_dg_iz_react_rate(const struct gkyl_dg_iz *viz,
  const struct gkyl_range *update_rng, const struct gkyl_range *phase_rng, 
  const struct gkyl_array *n_neut, const struct gkyl_array *vth_sq_neut, const struct gkyl_array *vth_sq_elc,
  struct gkyl_array *cflrate, struct gkyl_array *coef_iz);

/**
 * Delete updater.
 *
 * @param iz Updater to delete.
 */
void gkyl_dg_iz_release(gkyl_dg_iz *iz);
