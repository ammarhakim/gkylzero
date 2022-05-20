#pragma once

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

/**
 * Create a new diffusion equation object.
 *
 * @param Basis functions
 * @param range Range for use in indexing diffusion tensor
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn* gkyl_dg_const_diffusion_new(const struct gkyl_basis* basis, const double* D);

/**
 * Create a new diffusion equation object that lives on NV-GPU
 *
 * @param basis Basis functions
 * @param range Range for use in indexing diffusion tensor
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn* gkyl_dg_const_diffusion_cu_dev_new(const struct gkyl_basis* basis, const double* D);

