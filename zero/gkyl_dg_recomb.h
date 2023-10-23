#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Identifiers for different ionrecombation types
enum gkyl_dg_recomb_type
{
  GKYL_RECOMB_H, // Hydrogen plasma
  GKYL_RECOMB_HE, // Helium plasma
  GKYL_RECOMB_LI, // Lithium plasma
  GKYL_RECOMB_BE, // Beryllium ions
  GKYL_RECOMB_B,  // Boron ions
  GKYL_RECOMB_C,  // Carbon ions
  GKYL_RECOMB_N,  // Nitrogen ions
  GKYL_RECOMB_O,  // Oxygen ions
};

// Identifiers for self species to determine form of collision operator
enum gkyl_dg_recomb_self
{
  GKYL_RECOMB_ELC, // Electron species
  GKYL_RECOMB_ION, // Reacting ion species (increases charge state)
  GKYL_RECOMB_RECVR, // Resulting species (receives electron)
};

// Object type
typedef struct gkyl_dg_recomb gkyl_dg_recomb;

/**
 * Create new updater to calculate ionrecombation temperature or reaction rate
 * @param grid Grid object needed for fmax
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions (GK)
 * @param conf_rng Configuration range
 * @param phase_rng Phase range
 * @param elem_charge Elementary charge value
 * @param mass_elc Mass of the electron value
 * @param type_ion Enum for type of ion for ionrecombation (support H, He, Li)
 * @param charge_state Int for ion charge state
 * @param type_self Enum for species type (electron, ion or neutral)
 * @param all_gk Boolean to indicate if all 3 species are gyrokinetic
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_recomb* gkyl_dg_recomb_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis,
  struct gkyl_basis* pbasis, const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  double elem_charge, double mass_elc, double mass_self, enum gkyl_dg_recomb_type type_ion,
  int charge_state, enum gkyl_dg_recomb_self type_self, bool all_gk, const char *base, bool use_gpu); 

/**
 * Create new ionrecombation updater type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_recomb* gkyl_dg_recomb_cu_dev_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis,
  struct gkyl_basis* pbasis, const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  double elem_charge, double mass_elc, double mass_self, enum gkyl_dg_recomb_type type_ion,
  int charge_state, enum gkyl_dg_recomb_self type_self, bool all_gk, const char *base); 

/**
 * Compute ionrecombation collision term for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param recomb Ionrecombation object.
 * @param moms_elc Input electron moments
 * @param moms_ion Input ion moments
 * @param bmag Magnetic field used for GK fmax 
 * @param jacob_tot Total Jacobian used for GK fmax
 * @param b_i Unit bmag vector in Cartesian (X,Y,Z) components
 * @param distf_self Species self distribution function
 * @param coll_recomb Output reaction rate coefficient
 */
void gkyl_dg_recomb_coll(const struct gkyl_dg_recomb *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_ion,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i, 
  const struct gkyl_array *f_self, struct gkyl_array *coll_recomb, struct gkyl_array *cflrate);

void gkyl_dg_recomb_coll_cu(const struct gkyl_dg_recomb *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_ion,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *f_self, struct gkyl_array *coll_recomb, struct gkyl_array *cflrate);
/**
 * Delete updater.
 *
 * @param recomb Updater to delete.
 */
void gkyl_dg_recomb_release(gkyl_dg_recomb *recomb);
