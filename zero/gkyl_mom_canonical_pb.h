#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_mom_canonical_pb_auxfields { 
  const struct gkyl_array *hamil; // hamiltonian function
};

/**
 * Create new canonical-pb moment type object. 
 * Valid 'mom' strings are "MEnergy" which is the integral(hamil*f) 
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param phase_range Phase space range
 * @param mom Name of moment to compute.
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_mom_canonical_pb_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, 
  const char *mom, bool use_gpu);

/**
 * Create new canonical-pb moment type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range,
  const char *mom);

/**
 * Create new canonical-pb integrated moment type
 * object. Integrates (MEnergy).
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vel_range Velocity space range
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* gkyl_int_mom_canonical_pb_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, 
  bool use_gpu);

/**
 * Create new canonical-pb integrated moment type
 * object on NV-GPU: see new() method above for documentation.
 */
struct gkyl_mom_type* gkyl_int_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range);

/**
 * Set the auxiliary fields needed in computing moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_canonical_pb_set_auxfields(const struct gkyl_mom_type *momt,
  struct gkyl_mom_canonical_pb_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields needed in computing moments.
 * 
 * @param momt moment type.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_mom_canonical_pb_set_auxfields_cu(const struct gkyl_mom_type *momt,
  struct gkyl_mom_canonical_pb_auxfields auxin);

#endif