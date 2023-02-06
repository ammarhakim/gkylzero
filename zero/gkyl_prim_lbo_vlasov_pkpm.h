#pragma once

#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_prim_lbo_vlasov_pkpm_auxfields { 
  const struct gkyl_array *pvar;
};

/**
 * Create a new Vlasov primitive moment object.
 * Unique object for Vlasov solver in parallel-kinetic-perpendicular-moment model
 * No momentum corrections, only energy corrections for self-primitive moments
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration-space range
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to Vlasov (with fluid coupling) primitive moment object
 */
struct gkyl_prim_lbo_type* 
gkyl_prim_lbo_vlasov_pkpm_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  bool use_gpu);

/**
 * Create a new Vlasov PKPM primitive moment object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_prim_lbo_type* 
gkyl_prim_lbo_vlasov_pkpm_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Set the auxiliary fields (e.g. P = p_perp I + (p_parallel - p_perp) bb) needed in calculating the primitive moments.
 * 
 * @param prim prim_lbo_type pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_prim_lbo_vlasov_pkpm_set_auxfields(const struct gkyl_prim_lbo_type *prim,
  struct gkyl_prim_lbo_vlasov_pkpm_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set the auxiliary fields (e.g. P = p_perp I + (p_parallel - p_perp) bb) needed in calculating the primitive moments.
 * 
 * @param prim prim_lbo_type pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_prim_lbo_vlasov_pkpm_set_auxfields_cu(const struct gkyl_prim_lbo_type *prim,
  struct gkyl_prim_lbo_vlasov_pkpm_auxfields auxin);

#endif
