#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>


// Struct containing the pointers to auxiliary fields.
struct gkyl_mom_vlasov_sr_auxfields { 
  const struct gkyl_array *gamma; // gamma = sqrt(1 + p^2)
  const struct gkyl_array *gamma_inv; // gamma_inv = 1/gamma = 1/sqrt(1 + p^2)
  const struct gkyl_array *V_drift; // bulk fluid velocity (computed from M0*V_drift = M1i with weak division)
  const struct gkyl_array *GammaV2; // Gamma^2 = 1/(1 - V_drift^2/c^2), Lorentz boost factor squared from bulk fluid velocity
};

/**
 * Create new special relativistic Vlasov moment type object. 
 * Valid 'mom' strings are "M0", "M1i", "Ni" = (M0, M1i) (1+vdim components)
 * "Energy" = gamma*mc^2 moment (m = c = 1)
 * "Pressure" = n*T where n is rest frame density and T is temperature
 * "Pressure" moment is 1/vdim*gamma_inv*(Gamma2*(gamma - V_drift . p)^2 - 1)
 * "Tij" = stress energy tensor (Energy, Energy flux (vdim components), Stress tensor (vdim*(vdim+1))/2 components))
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vel_range Velocity space range
 * @param mom Name of moment to compute.
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_mom_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  const char *mom, bool use_gpu);

/**
 * Create new special relativistic Vlasov moment type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  const char *mom);

/**
 * Create new special relativistic Vlasov integrated moment type
 * object.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vel_range Velocity space range
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_int_mom_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  bool use_gpu);

/**
 * Create new special relativistic Vlasov integrated moment type
 * object on NV-GPU: see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_int_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range);

/**
 * Set the auxiliary fields needed in computing moments.
 * Not every auxiliary field is used by every moment, but setting pointers for 
 * all of them for flexibility.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_vlasov_sr_set_auxfields(const struct gkyl_mom_type *momt,
  struct gkyl_mom_vlasov_sr_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields needed in computing moments.
 * Not every auxiliary field is used by every moment, but setting pointers for 
 * all of them for flexibility.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_vlasov_sr_set_auxfields_cu(const struct gkyl_mom_type *momt,
  struct gkyl_mom_vlasov_sr_auxfields auxin);

#endif
