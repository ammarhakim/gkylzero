#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_pkpm_kinetic_vel_auxfields {
  const struct gkyl_array *c_0; // Normalization closure (c_0 F source).
  const struct gkyl_array *alpha_0; // Pressure force a_P/xi closing the first moment.
  const struct gkyl_array *c_h; // Heating rate closing the scalar second moment.
  const struct gkyl_array *vel_flux_surf; // Modal expansion of the fluxes on the v_par faces (each cell's lower face; the velocity ghost cell's slot is the upper domain edge).
};

// Input packaged as a struct
struct gkyl_dg_pkpm_kinetic_vel_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_basis *phase_basis; // Phase-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range for indexing auxiliary fields.
  const struct gkyl_range *phase_range; // Phase-space range for indexing the stored fluxes.
  struct gkyl_dg_pkpm_kinetic_vel_auxfields auxfields; // Auxiliary fields.
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new equation object for the velocity-space pass (stage B) of the PKPM
 * kinetic update: the pressure-force volume term integral(d_v phi a_P [F_0, G]) and
 * the surface term scattering the stored v_par fluxes (characteristic outflow fluxes at
 * the velocity domain edges). Run after the configuration-space pass and the a_P closure, accumulating
 * into the same right-hand side.
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_kinetic_vel_inp struct.
 * @return Pointer to the equation object.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_kinetic_vel_inew(const struct gkyl_dg_pkpm_kinetic_vel_inp *inp);

/**
 * Create a new equation object for the velocity-space pass on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_kinetic_vel_cu_dev_inew(
  const struct gkyl_dg_pkpm_kinetic_vel_inp *inp
);

/**
 * Set the auxiliary fields needed in the update.
 *
 * @param eqn Equation pointer.
 * @param auxin Struct of auxiliary fields.
 */
void gkyl_dg_pkpm_kinetic_vel_set_auxfields(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_vel_auxfields auxin
);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set the auxiliary fields needed in the update.
 *
 * @param eqn Equation pointer.
 * @param auxin Struct of auxiliary fields.
 */
void gkyl_dg_pkpm_kinetic_vel_set_auxfields_cu(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_vel_auxfields auxin
);
#endif
