#pragma once

#include <gkyl_proj_on_basis.h>
#include <gkyl_util.h>

// Initialization data for a 2D simulation using the Euler equations, run with static, block-structured mesh refinement with a single refinement patch.
struct euler2d_single_init {
  int base_Nx;
  int base_Ny;
  int ref_factor;

  double coarse_x1;
  double coarse_y1;
  double coarse_x2;
  double coarse_y2;

  double refined_x1;
  double refined_y1;
  double refined_x2;
  double refined_y2;

  evalf_t eval;
  double gas_gamma;

  double cfl_frac;
  double t_end;
};

/**
* Run a 2D simulation using the Euler equations, with static, block-structured mesh refinement with a single refinement patch.
*
* @param argc Nubmer of command line arguments passed to the function.
* @param argv Array of command line arguments passed to the function/
* @param init Initialization data for the 2D Euler equations.
*/
void euler2d_run_single(int argc, char **argv, struct euler2d_single_init* init);