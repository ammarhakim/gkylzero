#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_boltzmann_photon_auxfields { 
  const struct gkyl_array *kpar_abs; // |kpar| absolute value of the parallel momentum
  const struct gkyl_array *jacob_vel_inv; // inverse velocity space Jacobian for mapped velocity grids (default is linear mapping)
};

/**
 * Create a new photon Boltzmann equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param vel_range Velocity space range for use in indexing |kpar|
 * @param light_speed Speed of light
 * @param rho_curv Radius of curvature of the magnetic field
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to photon Boltzmann equation object
 */
struct gkyl_dg_eqn* gkyl_dg_boltzmann_photon_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* vel_range, 
  double light_speed, double rho_curv, bool use_gpu);

/**
 * Create new photon Boltzmann equation object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_eqn* gkyl_dg_boltzmann_photon_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* vel_range, 
  double light_speed, double rho_curv);

/**
 * Set the auxiliary fields needed in updating the equation system (such as |kpar| for advection in kperp).
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_boltzmann_photon_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_boltzmann_photon_auxfields auxin);


#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set the auxiliary fields needed in updating the equation system (such as |kpar| for advection in kperp).
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_boltzmann_photon_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_boltzmann_photon_auxfields auxin);

#endif
