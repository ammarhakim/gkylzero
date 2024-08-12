#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

/**
 * Construct a velocity map using a cubic C^1 representation based on input
 * mapping function. Returns vdim 1D mappings and vdim 1D quadratic 
 * representations of the inverse Jacobian. 
 *
 * @param vgrid Velocity-space grid object
 * @param vrange Velocity-space range
 * @param eval_vmap Velocity mapping used by dg_basis_ops to construct C^1 mapping.
 * @param ctx Context for function evaluation. Can be NULL.
 * @param vmap C^1 cubic representation of mapping in each velocity-space dimension
 * @param jacob_vel_inv Inverse of the Jacobian (derivative of vmap) in each velocity space dimension.
 * @param jacob_vel_gauss Total velocity space Jacobian evaluated at Gauss-Legendre quadrature points.
 */
void gkyl_velocity_map_cubic_new(const struct gkyl_rect_grid *vgrid, 
  const struct gkyl_range *vrange, 
  evalf_t eval_vmap, void *ctx, 
  struct gkyl_array *vmap, struct gkyl_array *jacob_vel_inv, struct gkyl_array *jacob_vel_gauss);