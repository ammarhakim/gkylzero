#pragma once

#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_prim_lbo_vlasov_with_fluid_auxfields { 
  const struct gkyl_array *fluid;
};

/**
 * Create a new Vlasov primitive moment object.
 * Unique object for when a Vlasov solver is coupled to a fluid solver
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration-space range
 * @return Pointer to Vlasov (with fluid coupling) primitive moment object
 */
struct gkyl_prim_lbo_type* 
gkyl_prim_lbo_vlasov_with_fluid_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  bool use_gpu);

/**
 * Create a new Vlasov primitive moment object for fluid coupling on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_prim_lbo_type* 
gkyl_prim_lbo_vlasov_with_fluid_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Set the auxiliary fields (e.g. nT_perp or nT_z) 
 * needed in calculating the primitive moments.
 * 
 * @param prim prim_lbo_type pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_prim_lbo_vlasov_with_fluid_set_auxfields(const struct gkyl_prim_lbo_type *prim,
  struct gkyl_prim_lbo_vlasov_with_fluid_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set the auxiliary fields (e.g. nT_perp or nT_z) 
 * needed in calculating the primitive moments.
 * 
 * @param prim prim_lbo_type pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_prim_lbo_vlasov_with_fluid_set_auxfields_cu(const struct gkyl_prim_lbo_type *prim,
  struct gkyl_prim_lbo_vlasov_with_fluid_auxfields auxin);

#endif
