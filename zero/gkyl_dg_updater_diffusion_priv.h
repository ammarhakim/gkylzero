#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_dg_eqn.h>

struct gkyl_dg_updater_diffusion {
  enum gkyl_diffusion_id diffusion_id; // Type of diffusion (e.g. isotropic vs. general diffusion tensor)
  bool use_gpu; // Boolean for if the update is performed on GPUs
  struct gkyl_dg_eqn *eqn_diffusion; // Equation object
  struct gkyl_hyper_dg *up_diffusion; // solvers for specific diffusion equation

  double diffusion_tm; // total time spent in computing diffusion equation
};
