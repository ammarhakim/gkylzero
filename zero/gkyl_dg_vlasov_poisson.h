#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_vlasov_poisson_auxfields {
  const struct gkyl_array *fac_phi; // Pointer to fac*phi, where phi is the potential.
  const struct gkyl_array *vecA; // Pointer to q/m*A, where A is the vector potential.
};

/**
 * Create a new Vlasov-Poisson equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param field_id enum to determine what type of Vlasov-Poisson system (phi only vs. phi and A)
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to Vlasov-Poisson equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_poisson_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  enum gkyl_field_id field_id, bool use_gpu);

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
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, 
  enum gkyl_field_id field_id);

/**
 * Set the auxiliary fields (e.g. q/m, external force) needed in updating the force terms.
 * 
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_poisson_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_poisson_auxfields auxin);

/**
 * Set up function to apply wall boundary conditions.
 * 
 * @param eqn Equation pointer.
 * @param dir Direction to apply wall boundary conditions.
 * @param pbasis Phase space basis
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods
 */

struct gkyl_array_copy_func* gkyl_vlasov_poisson_wall_bc_create(const struct gkyl_dg_eqn *eqn, 
  int dir, const struct gkyl_basis* pbasis);

/**
 * Set up function to apply absorbing boundary conditions.
 * 
 * @param eqn Equation pointer.
 * @param dir Direction to apply absorbing boundary conditions.
 * @param pbasis Phase space basis
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods
 */

struct gkyl_array_copy_func* gkyl_vlasov_poisson_absorb_bc_create(const struct gkyl_dg_eqn *eqn, 
  int dir, const struct gkyl_basis* pbasis);

/**
 * Release boundary conditions function.
 * 
 * @param bc Pointer to array_copy_func.
 */

void gkyl_vlasov_poisson_bc_release(struct gkyl_array_copy_func* bc);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. q/m, external force) needed in
 * updating the force terms.
 * 
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_vlasov_poisson_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_poisson_auxfields auxin);

/**
 * CUDA device function to set up function to apply wall boundary conditions.
 * 
 * @param eqn Equation pointer.
 * @param dir Direction to apply wall boundary conditions.
 * @param pbasis Phase space basis
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods
 */

struct gkyl_array_copy_func* gkyl_vlasov_poisson_wall_bc_create_cu(const struct gkyl_dg_eqn *eqn, 
  int dir, const struct gkyl_basis* pbasis);

/**
 * CUDA device function to set up function to apply absorbing boundary conditions.
 * 
 * @param eqn Equation pointer.
 * @param dir Direction to apply absorbing boundary conditions.
 * @param pbasis Phase space basis
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods
 */

struct gkyl_array_copy_func* gkyl_vlasov_poisson_absorb_bc_create_cu(const struct gkyl_dg_eqn *eqn, 
  int dir, const struct gkyl_basis* pbasis);

#endif
