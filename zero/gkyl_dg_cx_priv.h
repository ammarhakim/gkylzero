#pragma once

// Private header, not for direct use in user code

#include <gkyl_dg_cx_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

typedef double (*dg_cx_react_ratef_t)(const double a, const double b, double vt_sq_ion_min, double vt_sq_neut_min, const double *m0, const double *prim_vars_ion, const double *prim_vars_neut, double* GKYL_RESTRICT v_sigma_cx) ;

// for use in kernel tables
typedef struct { dg_cx_react_ratef_t kernels[3]; } gkyl_cx_react_rate_kern_list;

//
// Serendipity basis kernels
// 

// CX reaction rate kernel list
GKYL_CU_D
static const gkyl_cx_react_rate_kern_list ser_cx_react_rate_kernels[] = {
  { NULL, sigma_cx_1x1v_ser_p1, sigma_cx_1x1v_ser_p2 }, // 0
  { NULL, sigma_cx_1x2v_ser_p1, sigma_cx_1x2v_ser_p2 }, // 1
  { NULL, sigma_cx_1x3v_ser_p1, sigma_cx_1x3v_ser_p2 }, // 2
  { NULL, sigma_cx_2x2v_ser_p1, sigma_cx_2x2v_ser_p2 }, // 3
  { NULL, sigma_cx_2x3v_ser_p1, sigma_cx_2x3v_ser_p2 }, // 4
  { NULL, sigma_cx_3x3v_ser_p1, sigma_cx_3x3v_ser_p2 }, // 5
};

struct gkyl_dg_cx_kernels {
  dg_cx_react_ratef_t react_rate;
};

struct gkyl_dg_cx {
  const struct gkyl_rect_grid *grid; // grid object
  int cdim; // number of configuration space dimensions
  int vdim_vl;
  int vdim_gk;
  int poly_order; // polynomial order of DG basis
  struct gkyl_basis *cbasis;
  struct gkyl_basis *pbasis_gk;
  struct gkyl_basis *pbasis_vl;
  const struct gkyl_range *conf_rng;
  const struct gkyl_range *conf_rng_ext;
  const struct gkyl_range *phase_rng; 
  
  double a; // Fitting function coefficient.
  double b; // Fitting function coefficient.
  double mass_ion;
  double mass_neut;
  double vt_sq_ion_min;
  double vt_sq_neut_min; 
  
  enum gkyl_ion_type type_ion;
  enum gkyl_react_self_type type_self;

  struct gkyl_dg_prim_vars_type *calc_prim_vars_ion;
  struct gkyl_dg_prim_vars_type *calc_prim_vars_neut;
  struct gkyl_dg_prim_vars_type *calc_prim_vars_neut_gk;

  bool use_gpu;
  struct gkyl_dg_cx *on_dev; // pointer to itself or device data

  struct gkyl_dg_cx_kernels *kernels;
  //dg_cx_react_ratef_t react_rate; // pointer to reaction rate kernel
};

void
dg_cx_choose_kernel_cu(struct gkyl_dg_cx_kernels *kernels,
  struct gkyl_basis pbasis_vl, struct gkyl_basis cbasis);

GKYL_CU_D
static void dg_cx_choose_kernel(struct gkyl_dg_cx_kernels *kernels,
  struct gkyl_basis pbasis_vl, struct gkyl_basis cbasis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    dg_cx_choose_kernel_cu(kernels, pbasis_vl, cbasis);
    return;
  }
#endif

  enum gkyl_basis_type basis_type = pbasis_vl.b_type;
  int pdim = pbasis_vl.ndim;
  int cdim = cbasis.ndim;
  int vdim = pdim - cdim;
  int poly_order = pbasis_vl.poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_HYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->react_rate = ser_cx_react_rate_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
}
