#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_vlasov_poisson_auxfields { 
  const struct gkyl_array *field; // q/m*(phi_tot,A_ext)
};

/**
 * Create a new Vlasov Poisson equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param vel_range Velocity space range for use in indexing p/gamma (velocity)
 * @param vmap Velocity mapping projected onto DG basis.
 * @param field_id type of EM fields (phi only, or (phi_tot,A_ext)).
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to Vlasov equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_poisson_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range,
  const struct gkyl_range* vel_range, const struct gkyl_array *vmap,
  enum gkyl_field_id field_id, bool use_gpu);

/**
 * Set the auxiliary fields (e.g. q/m*EM) needed in updating the force terms.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_dg_vlasov_poisson_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_poisson_auxfields auxin);
