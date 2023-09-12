#pragma once

// Private header, not for direct use in user code

// header files for ADAS data
#include "bilinear_interp.h"
#include "adf11.h"

// Primary struct in this updater.
struct gkyl_dg_iz {
  struct gkyl_rect_grid *grid; // conf grid object
  int cdim; // number of configuration space dimensions
  int poly_order; // polynomial order of DG basis

  const struct gkyl_range *conf_rng;
  const struct gkyl_range *phase_rng; 
  bool use_gpu; 
  
  double elem_charge; // elementary charge value
  double mass_elc; // mass of the electron

  int resM0;
  int resTe;
  double dlogM0;
  double dlogTe;
  double E;

  double minLogM0, minLogTe, maxLogM0, maxLogTe;

  struct gkyl_basis *cbasis;
  struct gkyl_basis *pbasis;

  struct gkyl_array *M0q;
  struct gkyl_array *Teq;
  struct gkyl_array *ioniz_data;
  struct gkyl_array *prim_vars_neut;
  struct gkyl_array *vtSq_elc;
  struct gkyl_array *vtSq_iz;
  struct gkyl_array *prim_vars_fmax; 
  struct gkyl_array *coef_iz;
  struct gkyl_array *fmax_iz; 

  struct gkyl_dg_prim_vars_type *calc_prim_vars_neut_gk;
  struct gkyl_dg_prim_vars_type *calc_prim_vars_elc_vtSq;

  struct gkyl_proj_maxwellian_on_basis *proj_max;
  
  struct gkyl_dg_iz *on_dev; // pointer to itself or device data
  //dg_iz_react_ratef_t react_rate; // pointer to reaction rate kernel
};


