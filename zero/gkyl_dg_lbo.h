#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

/**
 * Create a new LBO equation object.
 *
 * @param Basis functions
 * @param range Range for use in indexing LBO tensor
 * @return Pointer to LBO equation object
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double nuUSum,
  const double* nuUSum_l, const double* nuUSum_r, const double* nuVtSqSum_l, const double* nuVtSqSum_r);

/**
 * Create a new LBO equation object that lives on NV-GPU
 *
 * @param basis Basis functions
 * @param range Range for use in indexing LBO tensor
 * @return Pointer to LBO equation object
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double nuUSum,
  const double* nuUSum_l, const double* nuUSum_r, const double* nuVtSqSum_l, const double* nuVtSqSum_r);
