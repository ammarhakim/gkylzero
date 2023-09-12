#pragma once

#include <gkyl_basis.h>
#include <gkyl_dg_prim_vars_type.h>

/**
 * Create new Vlasov primitive variables type object.
 * Contains kernels for solving a linear system for u_i, vth^2, or both u_i and vth^2
 * Since each solve is independent, all the solves are done at the same time.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param prim_nm Name of primitive variable (u_i, vth2, or prim)
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_dg_prim_vars_type* 
gkyl_dg_prim_vars_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *prim_nm, bool use_gpu);

/**
 * Create new Vlasov primitive variables type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_prim_vars_type* 
gkyl_dg_prim_vars_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *prim_nm);
