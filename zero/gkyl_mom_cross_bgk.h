#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>

// Type for storing preallocated memory 
typedef struct gkyl_mom_cross_bgk_gyrokinetic gkyl_mom_cross_bgk_gyrokinetic;

/**
 * Allocate memory for use in cross moments calculation.
 *
 * @ param
 */
gkyl_mom_cross_bgk_gyrokinetic* gkyl_mom_cross_bgk_gyrokinetic_new(
  const struct gkyl_basis *phase_basis, const struct gkyl_basis *conf_basis);

/**
 * Compute the cross moments with moments of each species.
 *
 * @param
 */
void gkyl_mom_cross_bgk_gyrokinetic_advance(
  gkyl_mom_cross_bgk_gyrokinetic *up,
  const struct gkyl_range *conf_rng, const double beta,
  const double m_self, const struct gkyl_array *moms_self,
  const double m_other, const struct gkyl_array *moms_other,
  const double nu_sr, const double nu_rs, struct gkyl_array *moms_cross);

/**
 * Release memory needed in the cross moments calculation.
 *
 * @param up Memory to release.
 */
void gkyl_mom_cross_bgk_gyrokinetic_release(gkyl_mom_cross_bgk_gyrokinetic *up);
