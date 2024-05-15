#pragma once

#include <gkyl_proj_on_basis.h>
#include <gkyl_util.h>

// Initialization data for a 1D simulation using the Euler equations, run with static, patch-structured mesh refinement with a single refinement patch.
struct euler1d_single_init {
  int base_Nx;
  int ref_factor;

  double coarse_x1;
  double coarse_x2;
  
  double refined_x1;
  double refined_x2;

  evalf_t eval;
  double gas_gamma;

  bool low_order_flux;
  double cfl_frac;
  double t_end;
  int num_frames;
};

/**
* Run a 1D simulation using the Euler equations, with static, patch-structured mesh refinement with a single refinement patch.
*
* @param argc Number of command line arguments passed to the function.
* @param argv Array of command line arguments passed to the function.
* @param init Initialization data for the 1D Euler equations.
*/
void euler1d_run_single(int argc, char **argv, struct euler1d_single_init* init);

// Initialization data for a 1D simulation using the general relativistic Euler equations, run with static, patch-structured mesh refinement with a single refinement patch.
struct gr_euler1d_single_init {
  int base_Nx;
  int ref_factor;

  double coarse_x1;
  double coarse_x2;
  
  double refined_x1;
  double refined_x2;

  evalf_t eval;
  double gas_gamma;
  struct gkyl_gr_spacetime *spacetime;

  bool low_order_flux;
  double cfl_frac;
  double t_end;
  int num_frames;
};

/**
* Run a 1D simulation using the general relativistic Euler equations, with static, patch-structured mesh refinement with a single refinement patch.
*
* @param argc Number of command line arguments passed to the function.
* @param argv Array of command line arguments passed to the function.
* @param init Initialization data for the 1D general relativistic Euler equations.
*/
void gr_euler1d_run_single(int argc, char **argv, struct gr_euler1d_single_init* init);

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

  bool low_order_flux;
  double cfl_frac;
  double t_end;
  int num_frames;
};

/**
* Run a 2D simulation using the Euler equations, with static, block-structured mesh refinement with a single refinement patch.
*
* @param argc Number of command line arguments passed to the function.
* @param argv Array of command line arguments passed to the function.
* @param init Initialization data for the 2D Euler equations.
*/
void euler2d_run_single(int argc, char **argv, struct euler2d_single_init* init);

// Initialization data for a 2D simulation using the general relativistic Euler equations, run with static, block-structured mesh refinement with a single refinement patch.
struct gr_euler2d_single_init {
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
  struct gkyl_gr_spacetime *spacetime;

  bool low_order_flux;
  double cfl_frac;
  double t_end;
  int num_frames;
};

/**
* Run a 2D simulation using the general relativistic Euler equations, with static, block-structured mesh refinement with a single refinement patch.
*
* @param argc Number of command line arguments passed to the function.
* @param argv Array of command line arguments passed to the function.
* @param init Initialization data for the 2D general relativistic Euler equations.
*/
void gr_euler2d_run_single(int argc, char **argv, struct gr_euler2d_single_init* init);

// Initialization data for a 2D simulation using the coupled five-moment equations, run with static, block-structured mesh refinement with a single refinement patch.
struct five_moment_2d_single_init {
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

  evalf_t eval_elc;
  evalf_t eval_ion;
  evalf_t eval_field;

  double gas_gamma;
  double k0_elc;
  double k0_ion;

  double light_speed;
  double e_fact;
  double b_fact;

  double epsilon0;
  double mass_elc;
  double charge_elc;
  double mass_ion;
  double charge_ion;

  bool periodic_x;
  bool periodic_y;

  bool wall_x;
  bool wall_y;

  double cfl_frac;
  double t_end;
  int num_frames;
};

/**
* Run a 2D simulation using the coupled five-moment equations, with static, block-structured mesh refinement with a single refinement patch.
*
* @param argc Number of command line arguments passed to the function.
* @param argv Array of command line arguments passed to the function.
* @param init Initialization data for the 2D coupled five-moment equations.
*/
void five_moment_2d_run_single(int argc, char **argv, struct five_moment_2d_single_init* init);