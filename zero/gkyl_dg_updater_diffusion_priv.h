#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_diffusion {
  struct gkyl_dg_eqn *eqn_diffusion; // Equation object
  struct gkyl_hyper_dg *up_diffusion; // solvers for specific diffusion equation

  double diffusion_tm; // total time spent in computing diffusion equation
};
