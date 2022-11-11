#pragma once

// Private header, not for direct use in user code

#include <gkyl_dg_iz_kernels.h>

typedef void (*dg_iz_tempf_t)(double elem_charge, double mass, double E, 
  const double *vt_sq, double* GKYL_RESTRICT vt_sq_iz);
typedef double (*dg_iz_react_ratef_t)(double elem_charge, double mass, double E, double A, double K, double P, double X, 
  const double *n_neut, const double *vt_sq_neut, const double *vt_sq_elc, double* GKYL_RESTRICT coef_iz);

struct gkyl_dg_iz {
  struct gkyl_rect_grid grid; // grid object
  int cdim; // number of configuration space dimensions
  int poly_order; // polynomial order of DG basis
  
  double elem_charge; // elementary charge value
  double mass_elc; // mass of the electron

  double E; // Voronov ionization energy. 
  double A; // Voronov constant. 
  double K; // Voronov constant. 
  double P; // Voronov constant. 
  double X; // Voronov constant.
  
  dg_iz_tempf_t iz_temp; // pointer to ionization temperature kernel
  dg_iz_react_ratef_t react_rate; // pointer to reaction rate kernel
};

// for use in kernel tables
typedef struct { dg_iz_tempf_t kernels[3]; } gkyl_iz_temp_kern_list;
typedef struct { dg_iz_react_ratef_t kernels[3]; } gkyl_iz_react_rate_kern_list;

//
// Serendipity basis kernels
// 

// Ionization temperature kernel list
GKYL_CU_D
static const gkyl_iz_temp_kern_list ser_iz_temp_kernels[] = {
  { NULL, iz_temp_1x, iz_temp_1x }, // 0
  { NULL, iz_temp_2x, iz_temp_2x }, // 1
  { NULL, iz_temp_3x, iz_temp_3x }, // 2
};

// Voronov reaction rate kernel list
GKYL_CU_D
static const gkyl_iz_react_rate_kern_list ser_iz_react_rate_kernels[] = {
  { NULL, iz_react_rate_1x_ser_p1, iz_react_rate_1x_ser_p2 }, // 0
  { NULL, iz_react_rate_2x_ser_p1, iz_react_rate_2x_ser_p2 }, // 1
  { NULL, iz_react_rate_3x_ser_p1, iz_react_rate_3x_ser_p2 }, // 2
};
