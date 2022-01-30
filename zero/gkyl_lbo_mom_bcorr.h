#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>

/**
 * Create new lbo boundary correction moment type object. Valid 'mom' strings are "f", "vf".
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param mom Name of moment to compute.
 * @param vBoundary Values at the edges of velocity space.
 */
struct gkyl_mom_type* gkyl_vlasov_lbo_mom_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const char *mom, const double* vBoundary);

struct gkyl_mom_type* gkyl_vlasov_lbo_mom_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const char *mom, const double* vBoundary);
