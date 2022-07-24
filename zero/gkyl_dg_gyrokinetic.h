#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h> 
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_gyrokinetic_auxfields {
  const struct gkyl_array *bmag; // Pointer to magnetic field magnitude.
  const struct gkyl_array *jacobtot_inv; // Pointer to 1/(conf Jacobian * gyro-center coords Jacobian).
  const struct gkyl_array *cmag; // Pointer to parallel gradient coefficient.
  const struct gkyl_array *b_i; // Pointer to covariant components of B-field unit vector.
  const struct gkyl_array *phi; // Pointer to electrostatic potential.
  const struct gkyl_array *apar; // Pointer to A_\parallel.
  const struct gkyl_array *apardot; // Pointer to d(A_parallel)/dt.
};

/**
 * Create a new Gyrokinetic equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @return Pointer to Gyrokinetic equation object
 */
struct gkyl_dg_eqn* gkyl_dg_gyrokinetic_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const double charge, const double mass, bool use_gpu);

/**
 * Create a new Gyrokinetic equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param field_id enum to determine what type of EM fields (Gyrokinetic-Maxwell vs. neutrals)
 * @return Pointer to Gyrokinetic equation object
 */
struct gkyl_dg_eqn* gkyl_dg_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const double charge, const double mass);

/**
 * Set the auxiliary fields (e.g. geometry & EM fields) needed in computing
 * gyrokinetic updates.
 *
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gyrokinetic_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_gyrokinetic_auxfields auxin);

/**
 * Release boundary conditions function.
 * 
 * @param bc Pointer to array_copy_func.
 */

void gkyl_gyrokinetic_bc_release(struct gkyl_array_copy_func* bc);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set the auxiliary fields (e.g. geometry & EM fields)
 * needed in computing gyrokinetic updates.
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_gyrokinetic_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_gyrokinetic_auxfields auxin);

#endif


