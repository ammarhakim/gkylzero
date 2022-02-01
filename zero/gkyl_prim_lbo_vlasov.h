#pragma once

#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>

/**
 * Create a new Vlasov primitive moment object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @return Pointer to Vlasov primitive moment object
 */
struct gkyl_prim_lbo_type* gkyl_prim_lbo_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);

/**
 * Create a new Vlasov primitive moment object that lives on NV-GPU.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @return Pointer to Vlasov primitive moment object
 */
struct gkyl_prim_lbo_type* gkyl_prim_lbo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);
