#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_dg_calc_pkpm_em_coupling gkyl_dg_calc_pkpm_em_coupling;

/**
 * Create new updater to compute the source update for the fluid-EM coupling in the PKPM system.
 * Solves the linear system for num_species number of momentum equations and Ampere's Law:
 * 
 * J_s^{n+1} - 0.5*dt*(q_s^2/m_s^2*rho_s^n*E^{n+1} + q_s/m_s*J_s^{n+1} x B^n) = J_s^n + ext_forces
 * epsilon0*E^{n+1} + 0.5*dt*sum_s J_s^{n+1} = epsilon0*E^{n} - ext_currents
 * 
 * We use a time-centered scheme to solve for J_bar = (J^{n+1} + J^{n})/2.0, E_bar = (E^{n+1} + E^{n})/2.0.
 * Note that there are no source terms for rho^n and B^n so we can use the mass density and the 
 * magnetic field at the known time step as part of the implicit update.
 * rho^n is given by the moment of the PKPM kinetic equation
 * 
 * Updater also stores function pointer to reconstructing the energy in Euler/5-moment 
 * at the new time from the old pressure and new kinetic energy (computed from the updated momentum).
 * 
 * @param cbasis           Configuration space basis functions
 * @param mem_range        Configuration space range that sets the size of the bin_op memory
 * @param num_species       Number of fluid species in the update
 * @param qbym[num_species] Array of charge/mass for each fluid species
 * @param epsilon0         Permitivvity of free space
 * @param use_gpu          bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_pkpm_em_coupling* 
gkyl_dg_calc_pkpm_em_coupling_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range *mem_range, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0,
  bool use_gpu);

/**
 * Create new updater to compute fluid variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_pkpm_em_coupling* 
gkyl_dg_calc_pkpm_em_coupling_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range *mem_range, 
  int num_species, double qbym[GKYL_MAX_SPECIES], double epsilon0);

/**
 * Compute the updated fluid momentum and electric field implicitly from time-centered source solve.
 *
 * @param up                            Updater for computing fluid-EM coupling in the PKPM system.
 * @param app_accel[num_species]        Input array of applied accelerations (external forces) on each species 
 * @param ext_em                        Input array of external electromagnetic fields
 * @param app_current                   Input array of applied currents
 * @param vlasov_pkpm_moms[num_species] Input array of velocity moments of PKPM kinetic equation
 * @param euler_pkpm[num_species]       Input and output array of fluid momentum in PKPM update
 *                                      (update is done in place with momentum modified to new time)
 * @param em                            Input and output array of EM fields
 *                                      (update is done in place with electric field modified to new time)
 * 
 */
void gkyl_dg_calc_pkpm_em_coupling_advance(struct gkyl_dg_calc_pkpm_em_coupling *up, double dt, 
  const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], 
  const struct gkyl_array* ext_em, const struct gkyl_array* app_current, 
  const struct gkyl_array* vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  struct gkyl_array* euler_pkpm[GKYL_MAX_SPECIES], struct gkyl_array* em);

/**
 * Delete pointer to updater to compute fluid variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_pkpm_em_coupling_release(struct gkyl_dg_calc_pkpm_em_coupling *up);

/**
 * Host-side wrappers for fluid vars operations on device
 */

void gkyl_dg_calc_pkpm_em_coupling_advance_cu(struct gkyl_dg_calc_pkpm_em_coupling *up, double dt, 
  const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], 
  const struct gkyl_array* ext_em, const struct gkyl_array* app_current, 
  const struct gkyl_array* vlasov_pkpm_moms[GKYL_MAX_SPECIES], 
  struct gkyl_array* euler_pkpm[GKYL_MAX_SPECIES], struct gkyl_array* em);
