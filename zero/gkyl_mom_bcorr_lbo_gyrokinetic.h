#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

/**
 * Create new LBO Gyrokinetic boundary correction moment type object.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vBoundary Values at the edges of velocity space.
 * @param mass Mass of species
 * @param vel_range Velocity space range for use in indexing velocity mappping.
 * @param vmap_prime Derivative of the velocity space mapping.
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_mom_type* 
gkyl_mom_bcorr_lbo_gyrokinetic_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_basis* pbasis, const double* vBoundary, double mass,
  const struct gkyl_range *vel_range, const struct gkyl_array *vmap_prime, bool use_gpu);
