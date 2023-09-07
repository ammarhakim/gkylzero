#pragma once

#include <gkyl_euler_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in euler DG equation object creation
// functions

// Types for various kernels
typedef double (*euler_surf_t)(const double *w, const double *dx, const double gas_gamma, 
  const double *ul, const double *uc, const double *ur,
  const double *pl, const double *pc, const double *pr,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_euler_vol_kern_list;
typedef struct { euler_surf_t kernels[3]; } gkyl_dg_euler_surf_kern_list;

struct dg_euler {
  struct gkyl_dg_eqn eqn; // Base object  
  euler_surf_t surf[3]; // pointers to surface kernels
  double gas_gamma; // adiabatic index
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_euler_auxfields auxfields; // Auxiliary fields.
};

GKYL_CU_DH
static double
kernel_euler_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_1x_ser_p1(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_ij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_1x_ser_p2(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_ij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_2x_ser_p1(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_ij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_2x_ser_p2(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_ij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_3x_ser_p1(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_ij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  long cidx = gkyl_range_idx(&euler->conf_range, idx);

  return euler_vol_3x_ser_p2(xc, dx, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_ij, cidx),
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_euler_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_euler_vol_1x_ser_p1, kernel_euler_vol_1x_ser_p2 }, // 0
  { NULL, kernel_euler_vol_2x_ser_p1, kernel_euler_vol_2x_ser_p2 }, // 1
  { NULL, kernel_euler_vol_3x_ser_p1, kernel_euler_vol_3x_ser_p2 }, // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, euler_surfx_1x_ser_p1, euler_surfx_1x_ser_p2 }, // 0
  { NULL, euler_surfx_2x_ser_p1, euler_surfx_2x_ser_p2 }, // 1
  { NULL, euler_surfx_3x_ser_p1, euler_surfx_3x_ser_p2 }, // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, euler_surfy_2x_ser_p1, euler_surfy_2x_ser_p2 }, // 1
  { NULL, euler_surfy_3x_ser_p1, euler_surfy_3x_ser_p2 }, // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_euler_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, euler_surfz_3x_ser_p1, euler_surfz_3x_ser_p2 }, // 2
};

/**
 * Free euler equation object
 *
 * @param ref Reference counter for euler equation
 */
void gkyl_euler_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
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

  return euler->surf[dir](xcC, dxC, euler->gas_gamma, 
    (const double*) gkyl_array_cfetch(euler->auxfields.u_i, cidx_l),
    (const double*) gkyl_array_cfetch(euler->auxfields.u_i, cidx_c),
    (const double*) gkyl_array_cfetch(euler->auxfields.u_i, cidx_r), 
    (const double*) gkyl_array_cfetch(euler->auxfields.p_ij, cidx_l),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_ij, cidx_c),
    (const double*) gkyl_array_cfetch(euler->auxfields.p_ij, cidx_r), 
    qInL, qInC, qInR, qRhsOut);
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  return 0.;
}
