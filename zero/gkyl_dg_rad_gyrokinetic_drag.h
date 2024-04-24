#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Struct containing the pointers to auxiliary fields. These are all summed over radiating species
struct gkyl_dg_rad_gyrokinetic_auxfields {
  // A, alpha, beta, gamma, and V0 are fitting parameters.
  // C = 8/sqrt(pi)*(2*charge/mass)^(gamma/2)
  // D = A*(alpha+beta)/C
  // vmag = sqrt(vpar^2 + 2*B*mu/mass)
  // nu(vpar,mu) = D*vmag^(gamma)/(beta*(vmag/V0)^-alpha + alpha*(vmag/V0)^beta)
  const struct gkyl_array *nvnu_surf;  // surface drag for vpar direction: n*vpar*nu(vpar,mu)
  const struct gkyl_array *nvnu;  // volume drag for vpar direction: n*vpar*nu(vpar,mu)
  const struct gkyl_array *nvsqnu_surf;  // surface drag for mu direction: 2*n*mu*nu(vpar,mu)
  const struct gkyl_array *nvsqnu;  // volume drag for mu direction: 2*n*mu*nu(vpar,mu)
  const struct gkyl_array *vtsq;  // Vth squared of all configuration space
  double vtsq_min;  // Vth squared of lowest temperature at which radiation is evaluated
};

/**
 * Create a new gyrokinetic RAD drag term equation object.
 *
 * @param conf_basis Configuration-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param phase_range Phase-space range for use in indexing drag coefficients
 * @param conf_range Configuration-space range for use in indexing temperature
 * @return Pointer to RAD equation object
 */
struct gkyl_dg_eqn* gkyl_dg_rad_gyrokinetic_drag_new(const struct gkyl_basis* conf_basis, 
  const struct gkyl_basis* phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, bool use_gpu);

/**
 * Create a new RAD equation object that lives on NV-GPU
 *
 * @param conf_basis Configuration-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param phase_range Phase-space range for use in indexing drag coefficients
 * @param conf_range Configuration-space range for use in indexing temperature
 * @return Pointer to RAD equation object
 */
struct gkyl_dg_eqn* gkyl_dg_rad_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* conf_basis, 
  const struct gkyl_basis* phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range);

/**
 * Set auxiliary fields needed in updating the drag flux term.
 * These are nvnu_sum, nvsqnu_sum
 * 
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_rad_gyrokinetic_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_rad_gyrokinetic_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields needed in updating the drag flux term.
 */
void gkyl_rad_gyrokinetic_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_rad_gyrokinetic_auxfields auxin);

#endif
