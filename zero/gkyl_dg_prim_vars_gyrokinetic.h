#pragma once

#include <gkyl_basis.h>
#include <gkyl_dg_prim_vars_type.h>

/**
 * Create new gyrokinetic primitive variables type object.
 * Solves a 2 component linear system for u_par, vth^2
 * Since each solve is independent, both the solves are done at the same time.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_dg_prim_vars_type* 
gkyl_dg_prim_vars_gyrokinetic_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, bool use_gpu);

/**
 * Create new gyrokinetic primitive variables type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_prim_vars_type* 
gkyl_dg_prim_vars_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);
