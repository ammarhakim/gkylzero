#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_sheath_rarefaction_pot gkyl_sheath_rarefaction_pot;

/**
 * Create a new updater to modify the electrostatic potential at the
 * boundary to account for the rarefaction wave.
 *
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param is_elc_boltz Whether using the Boltzmann electron model.
 * @param elem_charge elementary charge (i.e. charge of singly ionized ion).
 * @param mass_e electron mass.
 * @param mass_i ion mass.
 * @param temp_boltz_elc Electron temperature (for Boltzmann electrons).
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_sheath_rarefaction_pot*
gkyl_sheath_rarefaction_pot_new(enum gkyl_edge_loc edge, const struct gkyl_basis *basis,
  bool is_elc_boltz, double elem_charge, double mass_e, double mass_i, double temp_boltz_elc,
  bool use_gpu);

/**
 * Modify the electrostatic potential at the boundary.
 *
 * @param up updater object.
 * @param skin_range Range of the skin where potential is modified.
 * @param surf_range Surface range the wall potential lives on.
 * @param moms_e first four moments (m0, m1, m2par, m2perp) of the electron distribution.
 * @param moms_i first four moments (m0, m1, m2par, m2perp) of the ion distribution.
 * @param phi_wall Wall potential.
 * @param phi Electrostatic potential.
 */
void gkyl_sheath_rarefaction_pot_advance(const struct gkyl_sheath_rarefaction_pot *up,
  const struct gkyl_range *skin_range, const struct gkyl_range *surf_range,
  const struct gkyl_array *moms_e, const struct gkyl_array *moms_i,
  const struct gkyl_array *phi_wall, struct gkyl_array *phi);

/**
 * Free memory associated with sheath_rarefaction_pot updater.
 *
 * @param up BC updater.
 */
void gkyl_sheath_rarefaction_pot_release(struct gkyl_sheath_rarefaction_pot *up);
