#pragma once

#include <gkyl_euler_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in Euler DG equation object (for parallel-kinetic-perpendicular-moment (pkpm) model) creation
// functions

// Types for various kernels
typedef void (*euler_pkpm_surf_t)(const double *w, const double *dx, 
  const double *u_il, const double *u_ic, const double *u_ir,
  const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir,
  const double *p_perpl, const double *p_perpc, const double *p_perpr,
  const double *statevecl, const double *statevecc, const double *statevecr, 
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_euler_pkpm_vol_kern_list;
typedef struct { euler_pkpm_surf_t kernels[3]; } gkyl_dg_euler_pkpm_surf_kern_list;

struct dg_euler_pkpm {
  struct gkyl_dg_eqn eqn; // Base object  
  euler_pkpm_surf_t surf[3]; // pointers to surface kernels
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_euler_pkpm_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_euler_pkpm_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);
  long cidx = gkyl_range_idx(&euler_pkpm->conf_range, idx);

  return euler_pkpm_vol_1x_ser_p1(xc, dx, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.div_p, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_perp_i, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.p_perp, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_euler_pkpm_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);
  long cidx = gkyl_range_idx(&euler_pkpm->conf_range, idx);

  return euler_pkpm_vol_1x_ser_p2(xc, dx, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_i, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.div_p, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_perp_i, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.p_perp, cidx),
    qIn, qRhsOut);
}

// PKPM Volume kernel list
GKYL_CU_D
static const gkyl_dg_euler_pkpm_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_euler_pkpm_vol_1x_ser_p1, kernel_euler_pkpm_vol_1x_ser_p2 }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
};

// PKPM Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, euler_pkpm_surfx_1x_ser_p1, euler_pkpm_surfx_1x_ser_p2 }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
};

// PKPM Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
};

// PKPM Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
};

/**
 * Free euler equation object
 *
 * @param ref Reference counter for euler equation
 */
void gkyl_euler_pkpm_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);

  long cidx_l = gkyl_range_idx(&euler_pkpm->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&euler_pkpm->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&euler_pkpm->conf_range, idxR);

  // Note for surface moments from Vlasov equation, center index owns *left* edge
  euler_pkpm->surf[dir](xcC, dxC, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_i, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_i, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_i, cidx_r),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_perp_i, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_perp_i, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u_perp_i, cidx_r),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.p_perp, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.p_perp, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.p_perp, cidx_r),
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
