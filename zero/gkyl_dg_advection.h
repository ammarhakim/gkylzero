#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_advection_auxfields { 
  const struct gkyl_array *u;
};

/**
 * Create a new advection equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing advection velocity
 * @return Pointer to advection equation object
 */
struct gkyl_dg_eqn* gkyl_dg_advection_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu);

/**
 * Create a new advection equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing advection velocity
 * @return Pointer to advection equation object
 */
struct gkyl_dg_eqn* gkyl_dg_advection_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range);

/**
 * Set the auxiliary fields (e.g. advection velocity u) needed in updating advection equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_advection_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_advection_auxfields auxin);

/**
 * Set up function to apply absorbing boundary conditions.
 * 
 * @param eqn Equation pointer.
 * @param dir Direction to apply absorbing boundary conditions.
 * @param cbasis Configuration space basis
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods
 */

struct gkyl_array_copy_func* gkyl_advection_absorb_bc_create(const struct gkyl_dg_eqn *eqn, 
  int dir, const struct gkyl_basis* cbasis);

/**
 * Release boundary conditions function.
 * 
 * @param bc Pointer to array_copy_func.
 */

void gkyl_advection_bc_release(struct gkyl_array_copy_func* bc);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. advection velocity u) needed in updating advection equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_advection_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_advection_auxfields auxin);

/**
 * CUDA device function to set up function to apply absorbing boundary conditions.
 * 
 * @param eqn Equation pointer.
 * @param dir Direction to apply absorbing boundary conditions.
 * @param cbasis Configuration space basis
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods
 */

struct gkyl_array_copy_func* gkyl_advection_absorb_bc_create_cu(const struct gkyl_dg_eqn *eqn, 
  int dir, const struct gkyl_basis* cbasis);

#endif
