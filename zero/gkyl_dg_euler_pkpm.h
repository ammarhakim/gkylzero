#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_euler_pkpm_auxfields { 
  const struct gkyl_array *moms;
  const struct gkyl_array *u_i;
  const struct gkyl_array *div_p;
  const struct gkyl_array *vth_sq;
};

/**
 * Create a new Euler equation object for parallel-kinetic-perpendicular-moment (pkpm) model.
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing velocity
 * @return Pointer to Euler equation object for parallel-kinetic-perpendicular-moment (pkpm) model
 */
struct gkyl_dg_eqn* gkyl_dg_euler_pkpm_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu);

/**
 * Create a new Euler equation object for parallel-kinetic-perpendicular-moment (pkpm) model that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing velocity
 * @return Pointer to Euler equation object for parallel-kinetic-perpendicular-moment (pkpm) model
 */
struct gkyl_dg_eqn* gkyl_dg_euler_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range);

/**
 * Set the auxiliary fields (e.g. velocity u = rho*u/rho) needed in updating Euler equation for parallel-kinetic-perpendicular-moment (pkpm) model.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_pkpm_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_pkpm_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. velocity u = rho*u/rho) needed in updating Euler equation for parallel-kinetic-perpendicular-moment (pkpm) model.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_pkpm_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_pkpm_auxfields auxin);

#endif
