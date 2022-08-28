#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>

/**
 * Compute the Lorentz boost factor for a given bulk velocity, V.
 * GammaV = 1/sqrt(1 - V^2/c^2)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute 1/(1 - V^2/c^2) using weak multiplication and division
 * 2. Project onto quadrature points and evaluate square root point wise
 * 3. Project back onto modal basis
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param V Input array which contain bulk velocity
 * @param Gamma Output array of Lorentz boost factor
 */
void gkyl_calc_sr_vars_Gamma(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, struct gkyl_range range,
  const struct gkyl_array* V, struct gkyl_array* GammaV);

/**
 * Compute the integration factor required to compute the plasma rest frame pressure
 * for a relativistic plasma, P = integral(f(x,p) p_fac(x,p) dp)
 * p_fac = [GammaV^2*(gamma_p - V . p)^2 - 1]/gamma_p
 * where GammaV = 1/sqrt(1 - V^2/c^2), gamma_p = sqrt(1 + p^2/c^2)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Pre-compute gamma_p (geometric factor)
 * 2. Compute GammaV^2 = 1/(1 - V^2/c^2) using weak multiplication and division
 * 3. Compute p_fac using weak multiplication and division
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param gamma_p Input array which contains particle Lorentz boost factor
 * @param V Input array which contains bulk velocity
 * @param Gamma Output array of Lorentz boost factor
 */
void gkyl_calc_sr_vars_p_fac(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, struct gkyl_range range,
  const struct gkyl_array* gamma_p, const struct gkyl_array* V, struct gkyl_array* p_fac);
