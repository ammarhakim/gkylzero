#pragma once

#include <read_adas.h>

// Private header, not for direct use in user code

// Primary struct in this updater.
struct gkyl_dg_recomb {
  const struct gkyl_rect_grid *grid; // conf grid object
  int cdim; // number of configuration space dimensions
  int poly_order; // polynomial order of DG basis

  const struct gkyl_range *conf_rng;
  const struct gkyl_range *phase_rng; 
  bool use_gpu;
  bool all_gk; 
  
  double elem_charge; // elementary charge value
  double mass_elc; // mass of the electron
  double mass_self; // mass of species, used for gk maxwellian projection

  int vdim;
  int resM0;
  int resTe;
  double dlogM0;
  double dlogTe;
  double minLogM0, minLogTe, maxLogM0, maxLogTe;

  enum gkyl_dg_recomb_self type_self;

  struct gkyl_basis *cbasis;
  struct gkyl_basis *pbasis;

  struct gkyl_range adas_rng;
  struct gkyl_basis adas_basis;
  struct gkyl_basis *basis_on_dev;
  
  struct gkyl_array *recomb_data;
  struct gkyl_array *vtSq_elc;
  struct gkyl_array *coef_m0;
  struct gkyl_array *coef_recomb;
  struct gkyl_array *prim_vars_ion;
  struct gkyl_array *udrift_ion;
  struct gkyl_array *vtSq_ion;

  struct gkyl_dg_prim_vars_type *calc_prim_vars_elc_vtSq;
  struct gkyl_dg_prim_vars_type *calc_prim_vars_ion_udrift;
  struct gkyl_dg_prim_vars_type *calc_prim_vars_ion_vtSq;

  struct gkyl_proj_maxwellian_on_basis *proj_max;
  
  struct gkyl_dg_recomb *on_dev; // pointer to itself or device data
  //dg_recomb_react_ratef_t react_rate; // pointer to reaction rate kernel
};


