#pragma once

// Private header, not for direct use in user code.

#include <gkyl_euler_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// Types for various kernels.
typedef double (*euler_iso_vol_t)(const double *w, const double *dxv, const double vth,
  const double* statevec, const double *uvar, double* GKYL_RESTRICT out);

typedef void (*euler_iso_surf_t)(const double *w, const double *dxv, double vth_,
  const double *statevecl, const double *statevecc, const double *statevecr,
  const double *uvarl, const double *uvarc, const double *uvar, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { euler_iso_vol_t kernels[3]; } gkyl_dg_euler_iso_vol_kern_list;
typedef struct { euler_iso_surf_t kernels[3]; } gkyl_dg_euler_iso_surf_kern_list;

//
// Serendipity basis kernels.
//

// Volume kernel list.
GKYL_CU_D
static const gkyl_dg_euler_iso_vol_kern_list ser_vol_kernels[] = { //TODO: rename all kernels using snake case and remove pOrder=3 cases
  { NULL, euler_iso_vol_1x_ser_p1, euler_iso_vol_1x_ser_p2 }, // 0
  { NULL, euler_iso_vol_2x_ser_p1, euler_iso_vol_2x_ser_p2 }, // 1
  { NULL, euler_iso_vol_3x_ser_p1, euler_iso_vol_3x_ser_p2 }, // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_euler_iso_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, euler_iso_surfx_1x_ser_p1, euler_iso_surfx_1x_ser_p2 }, // 0
  { NULL, euler_iso_surfx_2x_ser_p1, euler_iso_surfx_2x_ser_p2 }, // 1
  { NULL, euler_iso_surfx_3x_ser_p1, euler_iso_surfx_3x_ser_p1 }, // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_euler_iso_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, euler_iso_surfy_2x_ser_p1, euler_iso_surfy_2x_ser_p2 }, // 1
  { NULL, euler_iso_surfy_3x_ser_p1, euler_iso_surfy_3x_ser_p2 }, // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_euler_iso_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, euler_iso_surfz_3x_ser_p1, euler_iso_surfz_3x_ser_p2 }, // 2
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

struct dg_euler_iso {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  euler_iso_vol_t vol; // Volume kernel.
  euler_iso_surf_t surf[3]; // Surface terms.
  struct gkyl_range conf_range; // Configuration space range.
  double vth;
  struct gkyl_dg_euler_iso_auxfields auxfields; // Auxiliary fields.
};

/**
 * Free gyrokinetic eqn object.
 *
 * @param ref Reference counter for gyrokinetic eqn
 */
void gkyl_euler_iso_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double*  dx,
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_euler_iso *euler_iso = container_of(eqn, struct dg_euler_iso, eqn);
  long cidx = gkyl_range_idx(&euler_iso->conf_range, idx);
  return euler_iso->vol(xc, dx, euler_iso->vth,
    (const double*) gkyl_array_cfetch(euler_iso->auxfields.u_i, cidx), qIn, qRhsOut);
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
  struct dg_euler_iso *euler_iso = container_of(eqn, struct dg_euler_iso, eqn);

  long cidx_l = gkyl_range_idx(&euler_iso->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&euler_iso->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&euler_iso->conf_range, idxR);

  euler_iso->surf[dir](xcC, dxC, euler_iso->vth,
    (const double*) gkyl_array_cfetch(euler_iso->auxfields.u_i, cidx_l),
    (const double*) gkyl_array_cfetch(euler_iso->auxfields.u_i, cidx_c),
    (const double*) gkyl_array_cfetch(euler_iso->auxfields.u_i, cidx_r),
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
