#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>

/**
 * Create new Vlasov moment type object. Valid 'mom' strings are "M0",
 * "M1i", "M2", "M2ij", "M3i", "M3ijk"
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param mom Name of moment to compute.
 */
struct gkyl_mom_type* gkyl_vlasov_lbo_mom_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom);

/**
 * Set the value of v at the boundary.
 * 
 * @param momt Moment type pointer
 * @param qmem Boundary value of v,
 */
void gkyl_lbo_mom_set_vBoundary(const struct gkyl_mom_type *momt, const double vBoundary);

/**
 * Set which boundary to evaluate.
 * 
 * @param momt Moment type pointer
 * @param atLower True if evaluate at lower boundary, false if at upper,
 */
void gkyl_lbo_mom_set_atLower(const struct gkyl_mom_type *momt, const bool atLower);
