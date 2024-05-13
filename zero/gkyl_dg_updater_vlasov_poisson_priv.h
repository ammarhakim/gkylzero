#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_vlasov_poisson {
  bool use_gpu; // Boolean for if the update is performed on GPUs.
  struct gkyl_dg_eqn *eqn_vlasov; // Equation object.
  struct gkyl_hyper_dg *hdg_vlasov; // solvers for specific Vlasov-Poisson equation.

  double vlasov_poisson_tm; // total time spent in computing Vlasov-Poisson equation.
};
