#pragma once

// Private header, not for direct use in user code

// header files for ADAS data
#include "bilinear_interp.h"
#include "adf11.h"

// Primary struct in this updater.
struct gkyl_dg_iz {
  struct gkyl_rect_grid grid; // conf grid object
  int cdim; // number of configuration space dimensions
  int poly_order; // polynomial order of DG basis

  const struct gkyl_range *conf_rng;
  const struct gkyl_range *phase_rng; 
  bool use_gpu; 
  
  double elem_charge; // elementary charge value
  double mass_elc; // mass of the electron

  int vdim_gk;
  int vdim_vl;

  int resM0;
  int resTe;
  double dlogM0;
  double dlogTe;
  double *M0q;
  double *Teq;
  double *recomb_data;
  double *ioniz_data;
  double E; 

  struct gkyl_array *udrift_neut;
  struct gkyl_array *vtSq_elc;

  struct gkyl_dg_prim_vars_type *prim_vars_neut_udrift;
  struct gkyl_dg_prim_vars_type *prim_vars_elc_vtSq;

  struct gkyl_iz_kernels *kernels;  // iz_react_rate kernel.
  struct gkyl_iz_kernels *kernels_cu;  // device copy.
  
  //dg_iz_react_ratef_t react_rate; // pointer to reaction rate kernel
};


