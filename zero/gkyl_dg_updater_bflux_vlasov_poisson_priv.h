#pragma once

#include <gkyl_ghost_surf_calc.h>

struct gkyl_dg_updater_bflux_vlasov_poisson {
  struct gkyl_ghost_surf_calc *slvr;

  double bflux_tm; // total time spent in computing boundary fluxes.
};
