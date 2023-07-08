#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>

void gkyl_calc_em_vars_bvar_basis_inv(struct gkyl_basis basis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* bvar);

/**
 * Compute the magnetic field unit vector (b_i = B_i/|B|, three components) 
 * and unit tensor (b_i b_j = B_i B_j/|B|^2, 6 components)
 * Note order of operations is designed to minimize aliasing errors
 * 1. Compute unit tensor (b_i b_j = B_i B_j/|B|^2, 6 components) first using basis_exp_sq and basis_inv
 *    (see gkyl_basis_*_exp_sq.h and gkyl_basis_*_inv.h in kernels/basis/)
 * 2. Project diagonal components onto quadrature points, evaluate square root point wise, 
 *    and project back onto modal basis using basis_sqrt to obtain b_i (see gkyl_basis_*_sqrt.h in kernels/basis/)
 *
 * @param mem   Pre-allocated matrix memory for use in inversion operations
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param em    Input array which contain EM fields (Ex, Ey, Ez, Bx, By, Bz)
 * @param bvar  Output array of magnetic field unit vector and unit tensor
 */
void gkyl_calc_em_vars_bvar(gkyl_dg_bin_op_mem *mem, 
  struct gkyl_basis basis, const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* bvar);

void gkyl_calc_em_vars_ExB_basis_inv(struct gkyl_basis basis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* ExB);

/**
 * Compute the ExB velocity (E x B/|B|^2, 3 components)
 *
 * @param basis Basis functions used in expansions
 * @param range Range to apply division operator
 * @param em Input array which contain EM fields (Ex, Ey, Ez, Bx, By, Bz)
 * @param ExB Output array of E x B velocity
 */
void gkyl_calc_em_vars_ExB(gkyl_dg_bin_op_mem *mem, 
  struct gkyl_basis basis, const struct gkyl_range *range, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* ExB);

/**
 * Host-side wrappers for em vars operations on device
 */

void gkyl_calc_em_vars_bvar_basis_inv_cu(struct gkyl_basis basis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* bvar);

void gkyl_calc_em_vars_ExB_basis_inv_cu(struct gkyl_basis basis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* ExB);

void gkyl_calc_em_vars_bvar_cu(gkyl_dg_bin_op_mem *mem, 
  struct gkyl_basis basis, const struct gkyl_range *range, 
  const struct gkyl_array* em, struct gkyl_array* bvar);

void gkyl_calc_em_vars_ExB_cu(gkyl_dg_bin_op_mem *mem, 
  struct gkyl_basis basis, const struct gkyl_range *range, 
  const struct gkyl_array* em, struct gkyl_array* ExB);

