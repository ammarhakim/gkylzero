#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_vlasov_poisson_auxfields { 
  const struct gkyl_array *field; // (q/m)*(phi,A_ext).
};

/**
 * Create a new Vlasov Poisson equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing field.
 * @param phase_range Phase space range.
 * @param model_id enum to determine what type of Vlasov-Poisson model.
 * @param field_id enum to determine what type of fields (e.g. phi, or phi and A_ext).
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_poisson_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const struct gkyl_range* phase_range,
  enum gkyl_vpmodel_id model_id, enum gkyl_vpfield_id field_id, bool use_gpu);

/**
 * Set the auxiliary fields (e.g. q/m*EM) needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_poisson_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_poisson_auxfields auxin);
