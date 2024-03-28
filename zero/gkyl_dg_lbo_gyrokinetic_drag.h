#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_lbo_gyrokinetic_drag_auxfields { 
  const struct gkyl_array *nuSum;
  const struct gkyl_array *nuPrimMomsSum;
  const struct gkyl_array *m2self;
};

/**
 * Create a new gyrokinetic LBO drag term equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @param vel_range Velocity space range for use in indexing velocity mappping.
 * @param pgrid Phase-space grid object.
 * @param mass Species mass
 * @param gk_geom Gyrokinetic geometry object.
 * @param vmap Velocity space mapping.
 * @param vmap_prime Derivative of the velocity space mapping.
 * @param jacobvel Velocity space mapping Jacobian.
 * @param bounds_vel Velocity at the boundaries.
 * @return Pointer to LBO equation object
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_gyrokinetic_drag_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, const struct gkyl_rect_grid *pgrid,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_array *vmap, 
  const struct gkyl_array *vmap_prime, const struct gkyl_array *jacobvel, double *bounds_vel, bool use_gpu);

/**
 * Set auxiliary fields needed in updating the drag flux term.
 * These are bmag, nu, nu*u, and nu*vt^2.
 * 
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_gyrokinetic_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_gyrokinetic_drag_auxfields auxin);
