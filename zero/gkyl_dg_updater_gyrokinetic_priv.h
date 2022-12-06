#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_updater_gyrokinetic.h>

struct gkyl_dg_updater_gyrokinetic {
  struct gkyl_dg_eqn *dgeqn; // Equation object
  struct gkyl_hyper_dg *hyperdg; // solvers for specific gyrokinetic equation

  enum gkyl_gkeqn_id eqn_id;  // ID to distinguish the equation type.

  bool use_gpu;
};
