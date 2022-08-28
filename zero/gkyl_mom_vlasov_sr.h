#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>


// Struct containing the pointers to auxiliary fields.
struct gkyl_mom_vlasov_sr_auxfields { 
  const struct gkyl_array *p_over_gamma; // p/gamma (velocity)
};

/**
 * Create new special relativistic Vlasov moment type object. 
 * Valid 'mom' strings are "M0", "M1i"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vel_range Velocity space range
 * @param mom Name of moment to compute.
 */
struct gkyl_mom_type* gkyl_mom_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* vel_range,
  const char *mom);

/**
 * Create new special relativistic Vlasov moment type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* vel_range,
  const char *mom);

/**
 * Create new special relativistic Vlasov integrated moment type
 * object.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vel_range Velocity space range
 */
struct gkyl_mom_type* gkyl_int_mom_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* vel_range);

/**
 * Create new special relativistic Vlasov integrated moment type
 * object on NV-GPU: see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_int_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* vel_range);

/**
 * Set the auxiliary fields p/gamma needed in computing moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_vlasov_sr_set_auxfields(const struct gkyl_mom_type *momt,
  struct gkyl_mom_vlasov_sr_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields p/gamma needed in computing moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_vlasov_sr_set_auxfields_cu(const struct gkyl_mom_type *momt,
  struct gkyl_mom_vlasov_sr_auxfields auxin);

#endif
