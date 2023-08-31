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
 * @param cbasis Configuration space basis functions
 * @param D Constant diffusion coefficient
 * @param order Integer order of diffusion (default is grad^2, supports grad^4 and grad^6)
 * @param diffusion_id Enum for identifying type of diffusion ()
 * @param use_gpu bool to determine if on GPU
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn* gkyl_dg_diffusion_new(const struct gkyl_basis* cbasis, 
  double D, int order, enum gkyl_diffusion_id diffusion_id, bool use_gpu);

/**
 * Create a new diffusion equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param D Constant diffusion coefficient
 * @param diffusion_id Enum for identifying type of diffusion (default isotropic grad^2, also support grad^4 and grad^6)
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn* gkyl_dg_diffusion_cu_dev_new(const struct gkyl_basis* cbasis, 
  double D, int order, enum gkyl_diffusion_id diffusion_id);
