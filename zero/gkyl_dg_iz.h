#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Identifiers for different ionization types
enum gkyl_dg_iz_type
{
  GKYL_IZ_H, // Hydrogen plasma
  GKYL_IZ_HE, // Helium plasma
  GKYL_IZ_LI, // Lithium plasma
};

// Identifiers for self species to determine form of collision operator
enum gkyl_dg_iz_self
{
  GKYL_IZ_ELC, // Electron species
  GKYL_IZ_ION, // Resulting ion species (increases charge state)
  GKYL_IZ_DONOR, // Reacting species (donates electron)
};

// Object type
typedef struct gkyl_dg_iz gkyl_dg_iz;

/**
 * Create new updater to calculate ionization temperature or reaction rate
 * @param grid Grid object needed for fmax
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions (GK)
 * @param conf_rng Configuration range
 * @param phase_rng Phase range
 * @param elem_charge Elementary charge value
 * @param mass_elc Mass of the electron value
 * @param type_ion Enum for type of ion for ionization (support H, He, Li)
 * @param charge_state Int for ion charge state
 * @param type_self Enum for species type (electron, ion or neutral)
 * @param all_gk Boolean to indicate if all 3 species are gyrokinetic
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_iz* gkyl_dg_iz_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, double elem_charge,
  double mass_elc, enum gkyl_dg_iz_type type_ion, int charge_state, enum gkyl_dg_iz_self type_self, bool all_gk, bool use_gpu); 

/**
 * Create new ionization updater type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_iz* gkyl_dg_iz_cu_dev_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, double elem_charge,
  double mass_elc, enum gkyl_dg_iz_type type_ion, int charge_state, enum gkyl_dg_iz_self type_self, bool all_gk); 

/**
 * Compute ionization collision term for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param iz Ionization object.
 * @param moms_elc Input electron moments
 * @param moms_donor Input neutral moments
 * @param bmag Magnetic field used for GK fmax 
 * @param jacob_tot Total Jacobian used for GK fmax
 * @param bhat_vec Unit bmag vector in Cartesian (X,Y,Z) components
 * @param distf_self Species self distribution function
 * @param coll_iz Output reaction rate coefficient
 */

void gkyl_dg_iz_coll(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_donor,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *distf_self, struct gkyl_array *coll_iz, struct gkyl_array *cflrate);

void gkyl_dg_iz_coll_cu(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_donor,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *distf_self, struct gkyl_array *coll_iz, struct gkyl_array *cflrate);

/**
 * Delete updater.
 *
 * @param iz Updater to delete.
 */
void gkyl_dg_iz_release(gkyl_dg_iz *iz);
