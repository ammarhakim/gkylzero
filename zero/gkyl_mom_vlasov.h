#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>

/**
 * Create new Vlasov moment type object. Valid 'mom' strings are "M0",
 * "M1i", "M2", "M2ij", "M3i", "M3ijk", "FiveMoments"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param mom Name of moment to compute.
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* 
gkyl_mom_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom, bool use_gpu);

/**
 * Create new Vlasov moment type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_mom_type* 
gkyl_mom_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom);

/**
 * Create new integrated Vlasov moment type object. Lab-frame
 * integrated moments (M0, M1i and M2) are computed.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* 
gkyl_int_mom_vlasov_new(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, bool use_gpu);

/**
 * Create new integrated Vlasov moment type object on NV-GPU:
 * see new() method above for documentation.
 */
struct gkyl_mom_type* 
gkyl_int_mom_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis);
