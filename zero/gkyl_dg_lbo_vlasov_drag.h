#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_lbo_vlasov_drag_auxfields {
  const struct gkyl_array *nuSum; // total collisionality 
  const struct gkyl_array *nuPrimMomsSum; // sum of nu*prim_moms (first vdim components are sum_r nu_sr u_sr)
  const struct gkyl_array *vmap; // velocity space mapping for mapped velocity grids
  const struct gkyl_array *jacob_vel_inv; // inverse velocity space Jacobian for mapped velocity grids
};

/**
 * Create a new Vlasov LBO drag equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @param vel_range Velocity space range for use in indexing nonuniform velocity map and inverse Jacobian (if present)
 * @param pgrid Phase-space grid object.
 * @param use_vmap bool to determine if we are using mapped velocity grid kernels
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to LBO equation object
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_vlasov_drag_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const struct gkyl_range* vel_range,
  const struct gkyl_rect_grid *pgrid, bool use_vmap, bool use_gpu);

/**
 * Create a new Vlasov LBO drag equation object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_vlasov_drag_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const struct gkyl_range* vel_range,
  const struct gkyl_rect_grid *pgrid, bool use_vmap);

/**
 * Set auxiliary fields needed in updating the drag flux term.
 * These are nu, nu*prim_moms (nu*u, nu*vt^2).
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_vlasov_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_vlasov_drag_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields needed in updating the drag flux term.
 * These are nu, nu*prim_moms (nu*u, nu*vt^2).
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_vlasov_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_vlasov_drag_auxfields auxin);

#endif
