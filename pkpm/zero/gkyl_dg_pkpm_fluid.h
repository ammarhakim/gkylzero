#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_pkpm_fluid_auxfields {
  const struct gkyl_array *moms; // Kinetic moments [rho, p_parallel, p_perp, M1].
  const struct gkyl_array *bvar; // Magnetic field unit vector [bx, by, bz].
  const struct gkyl_array *u; // Flow velocity [ux, uy, uz].
  const struct gkyl_array *fluid_flux_surf; // Modal expansion of the momentum fluxes on the faces.
};

// Input packaged as a struct
struct gkyl_dg_pkpm_fluid_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range for indexing auxiliary fields.
  struct gkyl_dg_pkpm_fluid_auxfields auxfields; // Auxiliary fields.
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new equation object for the momentum equation of the refactored PKPM
 * solver: the exact volume term integral(d_d phi (M_i u_d + P_id)) with the pressure
 * tensor formed from its primitive factors, and the surface term scattering the stored
 * momentum fluxes.
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_fluid_inp struct.
 * @return Pointer to the equation object.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_fluid_inew(const struct gkyl_dg_pkpm_fluid_inp *inp);

/**
 * Create a new equation object for the momentum equation on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_fluid_cu_dev_inew(const struct gkyl_dg_pkpm_fluid_inp *inp);

/**
 * Set the auxiliary fields needed in the update.
 *
 * @param eqn Equation pointer.
 * @param auxin Struct of auxiliary fields.
 */
void gkyl_dg_pkpm_fluid_set_auxfields(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_fluid_auxfields auxin
);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set the auxiliary fields needed in the update.
 *
 * @param eqn Equation pointer.
 * @param auxin Struct of auxiliary fields.
 */
void gkyl_dg_pkpm_fluid_set_auxfields_cu(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_fluid_auxfields auxin
);
#endif
