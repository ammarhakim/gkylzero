#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

/**
 * Compute the momentum grid variables for special relativistic simulations
 * Uses special kernels which convert between a Gauss-Lobatto nodal basis and
 * our modal basis to insure continuity of the momentum grid variables.
 *
 * @param vgrid Momentum-space grid
 * @param vbasis Momentum-space basis
 * @param vrange Momentum-space range
 * @param gamma Output array of particle Lorentz boost factor, gamma = sqrt(1 + p^2) 
 * @param gamma_inv Output array of inverse particle Lorentz boost factor, 1/gamma = 1/sqrt(1 + p^2) 
 */
void gkyl_calc_sr_vars_init_p_vars(const struct gkyl_rect_grid *vgrid, 
  const struct gkyl_basis *vbasis, const struct gkyl_range *vrange,
  struct gkyl_array* gamma, struct gkyl_array* gamma_inv);

/**
 * Compute the square of the Lorentz boost factor for a given bulk velocity, V.
 * GammaV2 = 1/(1 - V^2/c^2)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute (1 - V^2/c^2) using basis_exp_sq (see gkyl_basis_*_exp_sq.h in kernels/basis/)
 * 2. Compute 1/(1 - V^2/c^2) using basis_inv (see gkyl_basis_*_inv.h in kernels/basis/)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param V Input array which contain bulk velocity
 * @param GammaV2 Output array of the square of the Lorentz boost factor
 */
void gkyl_calc_sr_vars_GammaV2(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV2);

/**
 * Compute the inverse of the Lorentz boost factor for a given bulk velocity, V.
 * GammaV_inv = sqrt(1 - V^2/c^2)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute GammV2_inv = 1 - V^2/c^2 using basis_exp_sq 
 *    (see gkyl_basis_*_exp_sq.h in kernels/basis/)
 * 2. Project onto quadrature points, evaluate square root point wise, 
 *    and project back onto modal basis using basis_sqrt (see gkyl_basis_*_sqrt.h in kernels/basis/)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param V Input array which contain bulk velocity
 * @param GammaV_inv Output array of inverse Lorentz boost factor
 */
void gkyl_calc_sr_vars_GammaV_inv(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV_inv);

/**
 * Host-side wrappers for sr vars operations on device
 */

void gkyl_calc_sr_vars_init_p_vars_cu(const struct gkyl_rect_grid *vgrid, 
  const struct gkyl_basis *vbasis, const struct gkyl_range *vrange,
  struct gkyl_array* gamma, struct gkyl_array* gamma_inv);

void gkyl_calc_sr_vars_GammaV2_cu(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV2);

void gkyl_calc_sr_vars_GammaV_inv_cu(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV_inv);
