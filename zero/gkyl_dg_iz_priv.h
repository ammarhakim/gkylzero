#pragma once

// Private header, not for direct use in user code

#include <gkyl_dg_iz_kernels.h>

typedef void (*dg_iz_react_ratef_t)(double elem_charge, double mass, double E, double A, double K, double P, double X, 
  const double *n_neut, const double *vt_sq_neut, const double *vt_sq_elc, double* GKYL_RESTRICT coef_iz);

// for use in kernel tables
typedef struct { dg_iz_react_ratef_t kernels[3]; } gkyl_iz_react_rate_kern_list;

//
// Serendipity basis kernels
// 

// Voronov reaction rate kernel list
GKYL_CU_D
static const gkyl_iz_react_rate_kern_list ser_iz_react_rate_list[] = {
  { NULL, iz_react_rate_1x_ser_p1, iz_react_rate_1x_ser_p2 }, // 0
  { NULL, iz_react_rate_2x_ser_p1, iz_react_rate_2x_ser_p2 }, // 1
  { NULL, iz_react_rate_3x_ser_p1, iz_react_rate_3x_ser_p2 }, // 2
};

struct gkyl_iz_kernels {
  dg_iz_react_ratef_t react_ratef;  // react_ratef kernel.
};

// Primary struct in this updater.
struct gkyl_dg_iz {
  struct gkyl_rect_grid grid; // conf grid object
  int cdim; // number of configuration space dimensions
  int poly_order; // polynomial order of DG basis
  struct gkyl_basis *basis;
  const struct gkyl_range *conf_rng;
  const struct gkyl_range *phase_rng; 
  
  double elem_charge; // elementary charge value
  double mass_elc; // mass of the electron

  double E; // Voronov ionization energy. 
  double A; // Voronov constant. 
  double K; // Voronov constant. 
  double P; // Voronov constant. 
  double X; // Voronov constant.

  int vdim_gk;
  int vdim_vl;

  struct gkyl_array *m2_temp;
  struct gkyl_array *u_sq_temp;
  struct gkyl_array *vth_sq_neut;
  struct gkyl_array *vth_sq_elc;
  struct gkyl_dg_bin_op_mem *mem; 

  struct gkyl_iz_kernels *kernels;  // iz_react_rate kernel.
  struct gkyl_iz_kernels *kernels_cu;  // device copy.
  
  //dg_iz_react_ratef_t react_rate; // pointer to reaction rate kernel
};

void
gkyl_iz_choose_react_ratef_kernel_cu(const struct gkyl_basis *basis, struct gkyl_iz_kernels *kers);


GKYL_CU_D
static dg_iz_react_ratef_t
gkyl_iz_choose_react_ratef_kernel(const struct gkyl_basis *basis)
{
  int cdim = basis->ndim;
  int poly_order = basis->poly_order;
  return ser_iz_react_rate_list[cdim-1].kernels[poly_order];
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to apply the sheath BC.

 * @param up Ionization updater.
 * @param moms_neut Neutral moments.
 * @param coef_iz Ionization reaction rate parameter.
 */
gkyl_iz_react_rate_cu(const struct gkyl_dg_iz *up, const struct gkyl_array *moms_neut,
  const struct gkyl_array *coef_iz);

#endif
