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
 * @param local_range_ext Local extended range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param basis Basis on which coefficients in array are expanded (a device pointer if use_gpu=true).
 * @param grid cartesian grid dynamic field is defined on.
 * @param elem_charge elementary charge (i.e. charge of singly ionized ion).
 * @param mass_e electron mass.
 * @param mass_i ion mass.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_sheath_rarefaction_pot*
gkyl_sheath_rarefaction_pot_new(enum gkyl_edge_loc edge, const struct gkyl_range *local_range_ext,
  const int *num_ghosts, const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid,
  double elem_charge, double mass_e, double mass_i, bool use_gpu);

/**
 * Modify the electrostatic potential at the boundary.
 *
 * @param up updater object.
 * @param moms_e first threee moments of the electron distribution.
 * @param m2par_e v_par^2 moment of the electron distribution.
 * @param moms_i first threee moments of the ion distribution.
 * @param m2par_i v_par^2 moment of the ion distribution.
 * @param phi_wall Wall potential.
 * @param phi Electrostatic potential.
 */
void gkyl_sheath_rarefaction_pot_advance(const struct gkyl_sheath_rarefaction_pot *up,
  const struct gkyl_array *moms_e, const struct gkyl_array *m2par_e,
  const struct gkyl_array *moms_i, const struct gkyl_array *m2par_i,
  const struct gkyl_array *phi_wall, struct gkyl_array *phi);

/**
 * Free memory associated with sheath_rarefaction_pot updater.
 *
 * @param up BC updater.
 */
void gkyl_sheath_rarefaction_pot_release(struct gkyl_sheath_rarefaction_pot *up);
