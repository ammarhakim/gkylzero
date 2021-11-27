#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

/**
 * Create a new Vlasov-Poisson equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param field_id enum to determine what type of Vlasov-Poisson system (phi only vs. phi and A)
 * @return Pointer to Vlasov-Poisson equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_poisson_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, enum gkyl_field_id field_id);

/**
 * Create a new Vlasov-Poisson equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param field_id enum to determine what type of Vlasov-Poisson system (phi only vs. phi and A)
 * @return Pointer to Vlasov-Poisson equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_poisson_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, enum gkyl_field_id field_id);

/**
 * Set the fac_phi = factor*phi array needed in updating the force terms.
 * This factor is fac = q/m for plasmas and fac = G*m for self-gravitating systems
 * 
 * @param eqn Equation pointer
 * @param fac_phi Pointer to potential scaled by fac,
 */
void gkyl_vlasov_poisson_set_fac_phi(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *fac_phi);

/**
 * Set the vecA = q/m*A array needed in updating the force terms when running with external fields.
 * 
 * @param eqn Equation pointer
 * @param vecA Pointer to vector potential scaled by q/m,
 */
void gkyl_vlasov_poisson_set_vecA(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *vecA);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device functions to set the fac_phi = factor*phi array needed in updating the force terms.
 * 
 * @param eqn Equation pointer
 * @param fac_phi Pointer to potential scaled by fac,
 */
void gkyl_vlasov_poisson_set_fac_phi_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *fac_phi);

/**
 * CUDA device functions to set the vecA = q/m*A array needed in updating the force terms when running with external fields.
 * 
 * @param eqn Equation pointer
 * @param vecA Pointer to vector potential scaled by q/m,
 */
void gkyl_vlasov_poisson_set_vecA_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *vecA);

#endif
