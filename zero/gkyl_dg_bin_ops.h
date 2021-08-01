#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

/**
 * Compute out = lop*rop. The c_oop, c_lop and c_rop are the
 * components into the DG fields to multiply (in case the field is a
 * vector field). For scalar fields c_oop = c_rop = c_lop = 0, for
 * example.
 *
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field in which to store product
 * @param out Output DG field
 * @param c_lop Component of left operand to use in product
 * @param lop Left operand DG field
 * @param c_rop Component of right operand to use in product
 * @param rop Right operand DG field
 */
void gkyl_dg_mul_op(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop);

/**
 * Same as gkyl_dg_mul_op, except operator is applied only on
 * specified range (sub-range of range containing the DG fields).
 *
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field in which to store product
 * @param out Output DG field
 * @param c_lop Component of left operand to use in product
 * @param lop Left operand DG field
 * @param c_rop Component of right operand to use in product
 * @param rop Right operand DG field
 * @param range Range to apply multiplication operator
 */
void gkyl_dg_mul_op_range(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range);

/**
 * Return FLOP count for multiplication operation.
 *
 * @param basis Basis functions used in expansions
 * @return FLOP count
 */
struct gkyl_kern_op_count gkyl_dg_mul_op_count(struct gkyl_basis basis);
