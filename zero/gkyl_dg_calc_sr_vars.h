#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

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
 * @param Gamma2V Output array of the square of the Lorentz boost factor
 */
void gkyl_calc_sr_vars_Gamma2(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV2);

/**
 * Compute the Lorentz boost factor for a given bulk velocity, V.
 * GammaV = 1/sqrt(1 - V^2/c^2)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute 1/(1 - V^2/c^2) using basis_exp_sq and basis_inv
 *    (see gkyl_basis_*_exp_sq.h and gkyl_basis_*_inv.h in kernels/basis/)
 * 2. Project onto quadrature points, evaluate square root point wise, 
 *    and project back onto modal basis using basis_sqrt (see gkyl_basis_*_sqrt.h in kernels/basis/)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param V Input array which contain bulk velocity
 * @param Gamma Output array of Lorentz boost factor
 */
void gkyl_calc_sr_vars_Gamma(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV);

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
 * @param Gamma Output array of Lorentz boost factor
 */
void gkyl_calc_sr_vars_Gamma_inv(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV_inv);
