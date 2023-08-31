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
 * @param conf_range Configuration-space range
 * @param prim_nm Name of primitive variable (u_par_i, u_par, or prim)
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

/**
 * Set the auxiliary fields(e.g. b_i) needed in computing the primitive variables.
 * 
 * @param pvt Primitive variables type pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(const struct gkyl_dg_prim_vars_type *pvt, struct gkyl_dg_prim_vars_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. b_i) needed in computing the primitive variables.
 * 
 * @param pvt Primitive variables type pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields_cu(const struct gkyl_dg_prim_vars_type *pvt, struct gkyl_dg_prim_vars_auxfields auxin);


#endif