#pragma once

#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>

/**
 * Create a new Maxwell equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param lightSpeed Speed of light
 * @param elcErrorSpeedFactor Factor multiplying lightSpeed for div E correction
 * @param mgnErrorSpeedFactor Factor multiplying lightSpeed for div B correction
 * @return Pointer to Maxwell equation object
 */
struct gkyl_dg_eqn* gkyl_dg_maxwell_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor, bool use_gpu);

/*
 * Create a new Maxwell equation object that lives on NV-GPU.
 *
 * @param cbasis Configuration space basis functions
 * @param lightSpeed Speed of light
 * @param elcErrorSpeedFactor Factor multiplying lightSpeed for div E correction
 * @param mgnErrorSpeedFactor Factor multiplying lightSpeed for div B correction
 * @return Pointer to Maxwell equation object
 */
struct gkyl_dg_eqn* gkyl_dg_maxwell_cu_dev_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor);

/**
 * Set up function to apply wall boundary conditions.
 * 
 * @param eqn Equation pointer.
 * @param dir Direction to apply wall boundary conditions.
 * @param cbasis basis
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods
 */

struct gkyl_array_copy_func* gkyl_maxwell_wall_bc_create(const struct gkyl_dg_eqn *eqn, 
  int dir, const struct gkyl_basis* cbasis);

/**
 * Release wall boundary conditions function.
 * 
 * @param bc Pointer to array_copy_func.
 */

void gkyl_maxwell_bc_release(struct gkyl_array_copy_func* bc);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set up function to apply wall boundary conditions.
 * 
 * @param eqn Equation pointer.
 * @param dir Direction to apply wall boundary conditions.
 * @param cbasis basis
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods
 */

struct gkyl_array_copy_func* gkyl_maxwell_wall_bc_create_cu(const struct gkyl_dg_eqn *eqn, 
  int dir, const struct gkyl_basis* cbasis);

#endif