#pragma once

// Private header, not for direct use in user code

#include <gkyl_dg_cx_kernels.h>

typedef double (*dg_cx_react_ratef_t)(double a, double b, const double *m0_neut, const double *u_ion, const double *u_neut, const double *vt_sq_ion, double vt_sq_ion_min, const double *vt_sq_neut, double vt_sq_neut_min, double* GKYL_RESTRICT v_sigma_cx)  ;

struct gkyl_dg_cx {
  struct gkyl_rect_grid grid; // grid object
  int cdim; // number of configuration space dimensions
  int vdim_vl; // number of vlasov velocity space dimensions
  int poly_order; // polynomial order of DG basis
  struct gkyl_basis *basis;
  const struct gkyl_range *conf_rng;
  const struct gkyl_range *phase_rng; 
  
  double a; // Fitting function coefficient.
  double b; // Fitting function coefficient.
  double mass_ion;

  struct gkyl_array *m2_temp;
  struct gkyl_array *udrift_neut;
  struct gkyl_array *udrift_ion; 
  struct gkyl_array *vth_sq_neut;
  struct gkyl_array *vth_sq_ion;
  struct gkyl_dg_bin_op_mem *mem; 

  dg_cx_react_ratef_t react_rate; // pointer to reaction rate kernel
};

// for use in kernel tables
typedef struct { dg_cx_react_ratef_t kernels[5]; } gkyl_cx_react_rate_kern_list;

//
// Serendipity basis kernels
// 

// CX reaction rate kernel list
GKYL_CU_D
static const gkyl_cx_react_rate_kern_list ser_cx_react_rate_kernels[] = {
  { NULL, sigma_cx_1x1v_ser_p1, sigma_cx_1x1v_ser_p2 }, // 0
  { NULL, sigma_cx_1x2v_ser_p1, sigma_cx_1x2v_ser_p2 }, // 1
  { NULL, sigma_cx_1x3v_ser_p1, sigma_cx_1x3v_ser_p2 }, // 2
  { NULL, sigma_cx_2x3v_ser_p1, sigma_cx_2x3v_ser_p2 }, // 3
  { NULL, sigma_cx_3x3v_ser_p1, sigma_cx_3x3v_ser_p2 }, // 4
};
