#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_pkpm_kinetic_conf_auxfields {
  const struct gkyl_array *bvar; // Magnetic field unit vector [bx, by, bz].
  const struct gkyl_array *grad_b; // Weak gradient of b (gkyl_dg_calc_pkpm_bvar_grad).
  const struct gkyl_array *u; // Flow velocity [ux, uy, uz].
  const struct gkyl_array *grad_u; // Weak gradient of u (gkyl_dg_calc_pkpm_vars_grad_u).
  const struct gkyl_array *theta; // T_perp/m.
  const struct gkyl_array *nu_sum; // Collisionality (sum over self and cross collisions).
  const struct gkyl_array *nu_vthsq_sum; // Sum of nu*vth^2 over self and cross collisions.
  const struct gkyl_array
    *conf_flux_surf; // Modal expansion of the fluxes on the configuration-space faces.
};

// Input packaged as a struct
struct gkyl_dg_pkpm_kinetic_conf_inp {
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions.
  const struct gkyl_basis *phase_basis; // Phase-space basis functions.
  const struct gkyl_range *conf_range; // Configuration-space range for indexing auxiliary fields.
  const struct gkyl_range *phase_range; // Phase-space range for indexing the stored fluxes.
  struct gkyl_dg_pkpm_kinetic_conf_auxfields auxfields; // Auxiliary fields.
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new equation object for the configuration-space pass (stage A) of the
 * PKPM kinetic update: the volume terms of the streaming, of the -v_par (bb:grad u)
 * and mirror forces in v_par, and of the G perpendicular source, together with the
 * surface term scattering the stored configuration-space fluxes. Everything except
 * the pressure-force volume term and the v_par surface term, which belong to the
 * velocity-space pass (stage B) once the force closure is known.
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_kinetic_conf_inp struct.
 * @return Pointer to the equation object.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_kinetic_conf_inew(const struct gkyl_dg_pkpm_kinetic_conf_inp *inp);

/**
 * Create a new equation object for the configuration-space pass on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_eqn *gkyl_dg_pkpm_kinetic_conf_cu_dev_inew(
  const struct gkyl_dg_pkpm_kinetic_conf_inp *inp
);

/**
 * Set the auxiliary fields needed in the update.
 *
 * @param eqn Equation pointer.
 * @param auxin Struct of auxiliary fields.
 */
void gkyl_dg_pkpm_kinetic_conf_set_auxfields(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_conf_auxfields auxin
);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set the auxiliary fields needed in the update.
 *
 * @param eqn Equation pointer.
 * @param auxin Struct of auxiliary fields.
 */
void gkyl_dg_pkpm_kinetic_conf_set_auxfields_cu(
  const struct gkyl_dg_eqn *eqn, struct gkyl_dg_pkpm_kinetic_conf_auxfields auxin
);
#endif
