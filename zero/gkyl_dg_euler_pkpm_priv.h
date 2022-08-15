#pragma once

#include <gkyl_euler_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in Euler DG equation object (for parallel-kinetic-perpendicular-moment (pkpm) model) creation
// functions

// Types for various kernels
typedef void (*euler_pkpm_perp_pressure_t)(const double *uvar, const double *ppar, const double *statevec, double* GKYL_RESTRICT out);

typedef double (*euler_pkpm_vol_t)(const double *w, const double *dx, 
  const double *uvar, const double *ppar, const double *pperp, const double *qpar, 
  const double *statevec, double* GKYL_RESTRICT out);

typedef void (*euler_pkpm_surf_t)(const double *w, const double *dx, 
  const double *ul, const double *uc, const double *ur,
  const double *pparl, const double *pparc, const double *pparr,
  const double *pperpl, const double *pperpc, const double *pperpr,
  const double *qparl, const double *qparc, const double *qparr,
  const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { euler_pkpm_perp_pressure_t kernels[3]; } gkyl_dg_euler_pkpm_perp_pressure_kern_list;
typedef struct { euler_pkpm_vol_t kernels[3]; } gkyl_dg_euler_pkpm_vol_kern_list;
typedef struct { euler_pkpm_surf_t kernels[3]; } gkyl_dg_euler_pkpm_surf_kern_list;

//
// Serendipity basis kernels
// 

// Perpendicular pressure from state variables kernel list
GKYL_CU_D
static const gkyl_dg_euler_pkpm_perp_pressure_kern_list ser_perp_pressure_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  // { NULL, euler_pkpm_perp_pressure_1x_ser_p1, euler_pkpm_perp_pressure_1x_ser_p2 }, // 0
  // { NULL, euler_pkpm_perp_pressure_2x_ser_p1, euler_pkpm_perp_pressure_2x_ser_p2 }, // 1
  // { NULL, euler_pkpm_perp_pressure_3x_ser_p1, euler_pkpm_perp_pressure_3x_ser_p2 }, // 2
};

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_euler_pkpm_vol_kern_list ser_vol_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  // { NULL, euler_pkpm_vol_1x_ser_p1, euler_pkpm_vol_1x_ser_p2 }, // 0
  // { NULL, euler_pkpm_vol_2x_ser_p1, euler_pkpm_vol_2x_ser_p2 }, // 1
  // { NULL, euler_pkpm_vol_3x_ser_p1, euler_pkpm_vol_3x_ser_p2 }, // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  // { NULL, euler_pkpm_surfx_1x_ser_p1, euler_pkpm_surfx_1x_ser_p2 }, // 0
  // { NULL, euler_pkpm_surfx_2x_ser_p1, euler_pkpm_surfx_2x_ser_p2 }, // 1
  // { NULL, euler_pkpm_surfx_3x_ser_p1, euler_pkpm_surfx_3x_ser_p2 }, // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  // { NULL, euler_pkpm_surfy_2x_ser_p1, euler_pkpm_surfy_2x_ser_p2 }, // 1
  // { NULL, euler_pkpm_surfy_3x_ser_p1, euler_pkpm_surfy_3x_ser_p2 }, // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_euler_pkpm_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  // { NULL, NULL, NULL }, // 1
  // { NULL, euler_pkpm_surfz_3x_ser_p1, euler_pkpm_surfz_3x_ser_p2 }, // 2
};

struct dg_euler_pkpm {
  struct gkyl_dg_eqn eqn; // Base object  
  euler_pkpm_perp_pressure_t perp_pressure; // pointer to pressure computation kernel
  euler_pkpm_vol_t vol; // pointer to volume kernel
  euler_pkpm_surf_t surf[3]; // pointers to surface kernels
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_euler_pkpm_auxfields auxfields; // Auxiliary fields.
};


/**
 * Free euler equation object
 *
 * @param ref Reference counter for euler equation
 */
void gkyl_euler_pkpm_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double*  dx, 
  const int*  idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);
  long cidx = gkyl_range_idx(&euler_pkpm->conf_range, idx);

  // Compute perpendicular pressure from state variables and parallel pressure in cell
  euler_pkpm->perp_pressure(
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.ppar, cidx),
    qIn, (double*) gkyl_array_fetch(euler_pkpm->auxfields.pperp, cidx)
  );

  return euler_pkpm->vol(xc, dx, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.ppar, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pperp, cidx),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.qpar, cidx),
    qIn, qRhsOut);
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
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);

  long cidx_l = gkyl_range_idx(&euler_pkpm->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&euler_pkpm->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&euler_pkpm->conf_range, idxR);

  // Compute perpendicular pressure from state variables 
  // and parallel pressure in cell in left, center, and right cells
  // TO DO: Inefficient, computes the perpendicular pressure in each cell 3 times
  euler_pkpm->perp_pressure(
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.ppar, cidx_l),
    qInL, (double*) gkyl_array_fetch(euler_pkpm->auxfields.pperp, cidx_l)
  );
  euler_pkpm->perp_pressure(
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.ppar, cidx_c),
    qInC, (double*) gkyl_array_fetch(euler_pkpm->auxfields.pperp, cidx_c)
  );
  euler_pkpm->perp_pressure(
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u, cidx_r),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.ppar, cidx_r),
    qInR, (double*) gkyl_array_fetch(euler_pkpm->auxfields.pperp, cidx_r)
  );

  euler_pkpm->surf[dir](xcC, dxC, 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.u, cidx_r), 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.ppar, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.ppar, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.ppar, cidx_r), 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pperp, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pperp, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.pperp, cidx_r), 
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.qpar, cidx_l),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.qpar, cidx_c),
    (const double*) gkyl_array_cfetch(euler_pkpm->auxfields.qpar, cidx_r), 
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
