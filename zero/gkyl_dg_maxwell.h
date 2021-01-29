#pragma once

#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_real_type.h>

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
  gkyl_real lightSpeed, gkyl_real elcErrorSpeedFactor, gkyl_real mgnErrorSpeedFactor);
