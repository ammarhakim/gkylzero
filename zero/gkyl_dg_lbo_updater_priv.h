#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_lbo_updater {
  struct gkyl_dg_eqn *coll_drag; // Collision drag equation
  struct gkyl_dg_eqn *coll_diff; // Collision diffusion equation
  struct gkyl_hyper_dg *drag; // solvers for drag terms
  struct gkyl_hyper_dg *diff; // solvers for diffusion terms
};
