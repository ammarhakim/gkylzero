#pragma once

#include <gkyl_basis.h>
#include <gkyl_dg_prim_vars_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_prim_vars_auxfields { 
  const struct gkyl_array *b_i; // covariant components of the field aligned unit vector.
};

/**
 * Create new primitive variables type object for transforming between Vlasov and GK
 * Contains kernels for solving a linear system for u_par b_i from GK or 
 * upar (u_i . b_i) and vth_GK^2 (M2/M0 - upar^2)
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param prim_nm Name of primitive variable (u_i, vth2, or prim)
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_dg_prim_vars_type* 
gkyl_dg_prim_vars_transform_vlasov_gk_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const char *prim_nm, bool use_gpu);

/**
 * Create new Vlasov primitive variables type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_prim_vars_type* 
gkyl_dg_prim_vars_transform_vlasov_gk_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const char *prim_nm);
