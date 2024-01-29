#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Identifiers for different ionization types
enum gkyl_dg_iz_type
{
  GKYL_IZ_H = 0,  // Hydrogen ions
  GKYL_IZ_HE = 1, // Helium ions
  GKYL_IZ_LI = 2, // Lithium ions
  GKYL_IZ_BE = 3, // Beryllium ions
  GKYL_IZ_B = 4,  // Boron ions
  GKYL_IZ_C = 5,  // Carbon ions
  GKYL_IZ_N = 6,  // Nitrogen ions
  GKYL_IZ_O = 7,  // Oxygen ions
  GKYL_IZ_AR = 8, // Argon ions
};

// Identifiers for self species to determine form of collision operator
enum gkyl_dg_iz_self
{
  GKYL_IZ_ELC = 0, // Electron species
  GKYL_IZ_ION = 1, // Resulting ion species (increases charge state)
  GKYL_IZ_DONOR = 2, // Reacting species (donates electron)
};

struct gkyl_dg_iz_inp {
  const struct gkyl_rect_grid* grid; // Grid object needed for fmax
  struct gkyl_basis* cbasis; // Configuration-space basis-functions
  struct gkyl_basis* pbasis; // Phase-space basis-functions
  const struct gkyl_range *conf_rng; // Configuration range
  const struct gkyl_range *phase_rng; // Phase range
  double mass_ion; // Mass of the ion 
  enum gkyl_dg_iz_type type_ion; // Enum for type of ion for ionization (H thru 0)
  int charge_state; // Ion charge state
  enum gkyl_dg_iz_self type_self; // Species type (ion, electron or donor)
  bool all_gk; // To indicate if all 3 interacting species are GK or not
  const char* base; // File path to locate adas data
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
 * @param up Updater to delete.
 */
void gkyl_dg_iz_release(struct gkyl_dg_iz *up);
