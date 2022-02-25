#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>

/**
 * Create new lbo gyrokinetic boundary correction moment type object.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vBoundary Values at the edges of velocity space.
 * @param mass Mass of species
 */
struct gkyl_mom_type* gkyl_mom_bcorr_lbo_gyrokinetic_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double* vBoundary, double mass);

struct gkyl_mom_type* gkyl_mom_bcorr_lbo_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double* vBoundary, double mass);
