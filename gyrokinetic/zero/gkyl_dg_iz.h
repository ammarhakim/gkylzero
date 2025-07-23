#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_dg_iz_inp {
  const struct gkyl_rect_grid* grid; // Grid object needed for fmax
  struct gkyl_basis* cbasis; // Configuration-space basis-functions
  struct gkyl_basis* pbasis; // Phase-space basis-functions
  const struct gkyl_range *conf_rng; // Configuration-space range
  const struct gkyl_range *conf_rng_ext; // Configuration-space extended range
  const struct gkyl_range *phase_rng; // Phase-space range
  double mass_ion; // Mass of the ion 
  enum gkyl_ion_type type_ion; // Enum for type of ion for ionization (H thru O, Ar)
  int charge_state; // Ion charge state
  enum gkyl_react_self_type type_self; // Species type (ion, electron or donor)
};

// Object type
typedef struct gkyl_dg_iz gkyl_dg_iz;

/**
 * Create new updater to calculate ionization temperature or reaction rate
 * @param gkyl_dg_iz_inp
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_iz* gkyl_dg_iz_new(struct gkyl_dg_iz_inp *inp, bool use_gpu); 

/**
 * Create new ionization updater type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_iz* gkyl_dg_iz_cu_dev_new(struct gkyl_dg_iz_inp *inp);

/**
 * Compute ionization collision term for use in neutral reactions. 
 * 
 *
 * @param iz Ionization object.
 * @param prim_vars_elc (n, upar, T/m) for the electrons
 * @param vtSq_iz1 First thermal speed for ionization fmax (primary electrons)
 * @param vtSq_iz2 Second thermal speed for ionization fmax (secondary electrons)
 * @param coef_iz Output reaction rate coefficient
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T]) 
 */

void gkyl_dg_iz_coll(const struct gkyl_dg_iz *up, 
  const struct gkyl_array *prim_vars_elc, 
  struct gkyl_array *vtSq_iz1, struct gkyl_array *vtSq_iz2,
  struct gkyl_array *coef_iz, struct gkyl_array *cflrate);

void gkyl_dg_iz_coll_cu(const struct gkyl_dg_iz *up, 
  const struct gkyl_array *prim_vars_elc, 
  struct gkyl_array *vtSq_iz1, struct gkyl_array *vtSq_iz2,
  struct gkyl_array *coef_iz, struct gkyl_array *cflrate);
  
/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_iz_release(struct gkyl_dg_iz *up);
