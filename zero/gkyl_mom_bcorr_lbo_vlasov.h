#pragma once

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>

/**
 * Create new lbo vlasov boundary correction moment type object.
 *
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions
 * @param vBoundary Values at the edges of velocity space.
 */
struct gkyl_mom_type* gkyl_mom_bcorr_lbo_vlasov_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary);

struct gkyl_mom_type* gkyl_mom_bcorr_lbo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary);
