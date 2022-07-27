#pragma once

#include <gkyl_euler_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in euler DG equation object creation
// functions

// Types for various kernels
typedef double (*euler_vol_t)(const double *w, const double *dx, double gas_gamma, 
  const double *u, const double *q, double* GKYL_RESTRICT out);

typedef void (*euler_surf_t)(const double *w, const double *dx, double gas_gamma, 
  const double *ul, const double *uc, const double *ur,
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { euler_vol_t kernels[3]; } gkyl_dg_euler_vol_kern_list;
typedef struct { euler_surf_t kernels[3]; } gkyl_dg_euler_surf_kern_list;

typedef void (*bc_funcf_t)(size_t nc, double *out, const double *inp, void *ctx);

//
// Serendipity basis kernels
// 

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_euler_vol_kern_list ser_vol_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // { NULL, euler_vol_1x_ser_p1, euler_vol_1x_ser_p2 }, // 0
  // { NULL, euler_vol_2x_ser_p1, euler_vol_2x_ser_p2 }, // 1
  // { NULL, euler_vol_3x_ser_p1, euler_vol_3x_ser_p2 }, // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // { NULL, euler_surfx_1x_ser_p1, euler_surfx_1x_ser_p2 }, // 0
  // { NULL, euler_surfx_2x_ser_p1, euler_surfx_2x_ser_p2 }, // 1
  // { NULL, euler_surfx_3x_ser_p1, euler_surfx_3x_ser_p2 }, // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // { NULL, NULL, NULL }, // 0
  // { NULL, euler_surfy_2x_ser_p1, euler_surfy_2x_ser_p2 }, // 1
  // { NULL, euler_surfy_3x_ser_p1, euler_surfy_3x_ser_p2 }, // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // { NULL, NULL, NULL },                 // 0
  // { NULL, NULL, NULL },                 // 1
  // { NULL, euler_surfz_3x_ser_p1, euler_surfz_3x_ser_p2 }, // 2
};

struct dg_euler {
  struct gkyl_dg_eqn eqn; // Base object    
  euler_vol_t vol; // pointer to volume kernel
  euler_surf_t surf[3]; // pointers to surface kernels
  bc_funcf_t absorb_bc; // Absorbing BCs function
  double gas_gamma; // adiabatic index
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_euler_auxfields auxfields; // Auxiliary fields.
};


/**
 * Free euler equation object
 *
 * @param ref Reference counter for euler equation
 */
void gkyl_euler_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double*  dx, 
  const int*  idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);
  return euler->vol(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx), qIn, qRhsOut);
}

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);

  long cidx_l = gkyl_range_idx(&euler->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&euler->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&euler->conf_range, idxR);

  euler->surf[dir](xcC, dxC, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx_l),
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx_c),
    (const double*) gkyl_array_cfetch(euler->auxfields.u, cidx_r), 
    qInL, qInC, qInR, qRhsOut);
}

GKYL_CU_D
static void
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  
}
