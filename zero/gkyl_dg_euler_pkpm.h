#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_euler_pkpm_auxfields { 
  const struct gkyl_array *vlasov_pkpm_moms;
  const struct gkyl_array *pkpm_u;
  const struct gkyl_array *pkpm_u_surf;
  const struct gkyl_array *pkpm_p_ij;
  const struct gkyl_array *pkpm_lax;
};

/**
 * Create a new Euler equation object for parallel-kinetic-perpendicular-moment (pkpm) model.
 * Euler equations are solved at order p (kinetic equations are solved at order 2*p)
 *
 * @param cbasis Configuration space basis functions
 * @param conf_range Configuration space range for use in indexing auxiliary variables
 * @param wv_eqn Wave equation object which contains information and functions for the specific fluid equation
 * @param geom Wave geometry object for computing fluctuations local to surfaces
 * @return Pointer to Euler equation object for parallel-kinetic-perpendicular-moment (pkpm) model
 */
struct gkyl_dg_eqn* gkyl_dg_euler_pkpm_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, bool use_gpu);

/**
 * Create new Euler equation object parallel-kinetic-perpendicular-moment (pkpm) model the lives on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_eqn* gkyl_dg_euler_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom);

/**
 * Set the auxiliary fields needed in updating Euler equation for parallel-kinetic-perpendicular-moment (pkpm) model.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_pkpm_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_pkpm_auxfields auxin);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device function to set auxiliary fields needed in updating Euler equation for parallel-kinetic-perpendicular-moment (pkpm) model.
 * 
 * @param eqn Equation pointer.
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_euler_pkpm_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_pkpm_auxfields auxin);

#endif
