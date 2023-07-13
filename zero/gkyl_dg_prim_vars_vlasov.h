#pragma once

#include <gkyl_basis.h>
#include <gkyl_dg_prim_vars_type.h>

/**
 * Create new Vlasov primitive variables type object.
 * Solves a (vdim+1) linear system for u_i, vth^2
 * Since each solve is independent, all the solves are done at the same time.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_dg_prim_vars_type* 
gkyl_dg_prim_vars_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, bool use_gpu);

/**
 * Create new Vlasov primitive variables type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_prim_vars_type* 
gkyl_dg_prim_vars_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);
