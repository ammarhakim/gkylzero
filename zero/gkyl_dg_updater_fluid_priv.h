#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_fluid {
  struct gkyl_dg_eqn *eqn_fluid; // Equation object
  struct gkyl_hyper_dg *up_fluid; // solvers for specific fluid equation

  double fluid_tm; // total time spent in computing fluid equation
};
