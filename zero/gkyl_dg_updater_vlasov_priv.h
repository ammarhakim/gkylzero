#pragma once

#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_vlasov {
  struct gkyl_dg_eqn *eqn_vlasov; // Equation object
  struct gkyl_hyper_dg *up_vlasov; // solvers for specific Vlasov equation

  double vlasov_tm; // total time spent in computing vlasov equation
};
