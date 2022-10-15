#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_euler_iso_auxfields {
  const struct gkyl_array *u_i;
};

/**
 * Create a new isothermal euler equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing velocity
 * @param vth Constant thermal velocity
 * @param use_gpu Boolean to determine whether equation object is on host or device
 * @return Pointer to euler equation object
 */
struct gkyl_dg_eqn* gkyl_dg_euler_iso_new(const struct gkyl_basis* cbasis,
  const struct gkyl_range* conf_range, double vth, bool use_gpu);

struct gkyl_dg_eqn* gkyl_dg_euler_iso_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_range* conf_range, double vth);

/**
 * Set the auxiliary fields (e.g. velocity u = rho*u/rho) needed in updating isothermal euler equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_iso_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_iso_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields (e.g. velocity u = rho*u/rho) needed in updating isothermal euler equation.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_iso_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_auxfields auxin);

#endif