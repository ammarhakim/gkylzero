#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>

/**
 * Compute the magnetic field unit vector (b_i = B_i/|B|, three components) 
 * and unit tensor (b_i b_j = B_i B_j/|B|^2, 6 components)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute unit tensor (b_i b_j = B_i B_j/|B|^2, 6 components) first using weak multiplication and division
 * 2. Project onto quadrature points and evaluate square root of diagonal components point wise
 * 3. Project diagonal components back onto modal basis to get b_i
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param em Input array which contain EM fields (Ex, Ey, Ez, Bx, By, Bz)
 * @param bvar Output array of magnetic field unit vector and unit tensor
 */
void gkyl_calc_em_vars_bvar(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, struct gkyl_range range,
  const struct gkyl_array* em, struct gkyl_array* bvar);

/**
 * Compute the ExB velocity (E x B/|B|^2, 3 components)
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param em Input array which contain EM fields (Ex, Ey, Ez, Bx, By, Bz)
 * @param ExB Output array of E x B velocity
 */
void gkyl_calc_em_vars_ExB(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, struct gkyl_range range,
  const struct gkyl_array* em, struct gkyl_array* ExB);

/**
 * Compute the ExB velocity inverse Lorentz boost factor and Lorentz boost factor
 * kappa_inv = sqrt(1 - |E x B|^2/(c^2 |B|^4)) = sqrt(1 - |E_perp|^2/(c^2 |B|^2))
 * kappa = (1/sqrt(1 - |E_perp|^2/(c^2 |B|^2)))
 *
 * @param mem Pre-allocated space for use in the division 
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param ExB Input array which contain the E x B velocity
 * @param kappa_inv Output array of E x B velocity *inverse* Lorentz boost factor
 * @param kappa Output array of E x B velocity Lorentz boost factor
 */
void gkyl_calc_em_vars_ExB_lorentz(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, struct gkyl_range range,
  const struct gkyl_array* ExB, struct gkyl_array* kappa_inv, struct gkyl_array* kappa);
