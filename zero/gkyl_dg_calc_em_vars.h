#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

/**
 * Compute the magnetic field unit vector (b_i = B_i/|B|, three components) 
 * and unit tensor (b_i b_j = B_i B_j/|B|^2, 6 components)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute unit tensor (b_i b_j = B_i B_j/|B|^2, 6 components) first using basis_exp_sq and basis_inv
 *    (see gkyl_basis_*_exp_sq.h and gkyl_basis_*_inv.h in kernels/basis/)
 * 2. Project diagonal components onto quadrature points, evaluate square root point wise, 
 *    and project back onto modal basis using basis_sqrt to obtain b_i (see gkyl_basis_*_sqrt.h in kernels/basis/)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param em Input array which contain EM fields (Ex, Ey, Ez, Bx, By, Bz)
 * @param bvar Output array of magnetic field unit vector and unit tensor
 */
void gkyl_calc_em_vars_bvar(const struct gkyl_basis* cbasis, const struct gkyl_range *range, 
  const struct gkyl_array* em, struct gkyl_array* bvar);

/**
 * Compute the ExB velocity (E x B/|B|^2, 3 components)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param em Input array which contain EM fields (Ex, Ey, Ez, Bx, By, Bz)
 * @param ExB Output array of E x B velocity
 */
void gkyl_calc_em_vars_ExB(const struct gkyl_basis* cbasis, const struct gkyl_range *range, 
  const struct gkyl_array* em, struct gkyl_array* ExB);

/**
 * Compute b_hat/kappa, the magnetic field unit vector divided by the ExB velocity Lorentz boost factor
 * for use in the special relativistic parallel-kinetic-perpendicular-moment (pkpm) model
 * Used in both collisionless advection div ([p_parallel/gamma b_hat/kappa] f) 
 * and moments (J_parallel = int [p_parallel/gamma b_hat/kappa] f)
 * Note: 1/kappa = 1 - |E x B|^2/(c^2 |B|^4)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute b_i^2 * (1 - |E x B|^2/(c^2 |B|^4)) using weak multiplication
 *    Note we have b_i^2 because we are storing the magnetic field unit tensor and we can
 *    obtain (E x B/|B|^2)^2 using basis_exp_sq (see gkyl_basis_*_exp_sq.h in kernels/basis/)
 * 2. Project onto quadrature points, evaluate square root point wise, 
 *    and project back onto modal basis using basis_sqrt (see gkyl_basis_*_sqrt.h in kernels/basis/)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param bvar Input array which contains magnetic field unit vector and unit tensor
 * @param ExB Input array which contain the E x B velocity
 * @param kappa_inv_b Output array of magnetic field unit vector 
 *                    divided by E x B velocity Lorentz boost factor b/kappa
 */
void gkyl_calc_em_vars_pkpm_kappa_inv_b(const struct gkyl_basis* cbasis, const struct gkyl_range *range, 
  const struct gkyl_array* bvar, const struct gkyl_array* ExB, struct gkyl_array* kappa_inv_b);

/**
 * Host-side wrappers for em vars operations on device
 */

void gkyl_calc_em_vars_bvar_cu(const struct gkyl_basis* cbasis, const struct gkyl_range *range, 
  const struct gkyl_array* em, struct gkyl_array* bvar);

void gkyl_calc_em_vars_ExB_cu(const struct gkyl_basis* cbasis, const struct gkyl_range *range, 
  const struct gkyl_array* em, struct gkyl_array* ExB);

void gkyl_calc_em_vars_pkpm_kappa_inv_b_cu(const struct gkyl_basis* cbasis, const struct gkyl_range *range, 
  const struct gkyl_array* bvar, const struct gkyl_array* ExB, struct gkyl_array* kappa_inv_b);

