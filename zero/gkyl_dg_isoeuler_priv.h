#pragma once

// Private header, not for direct use in user code.

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_isoeuler_kernels.h>
#include <gkyl_dg_isoeuler.h>


// Types for various kernels.
typedef double (*isoeuler_vol_t)(const double *w, const double *dxv, const double vth_,
  const double* statevec, const double *uvar, double* GKYL_RESTRICT out);

typedef void (*isoeuler_surf_t)(const double *w, const double *dxv, double vth_,
  const double *statevecl, const double *statevecc, const double *statevecr,
  const double *uvarl, const double *uvarc, const double *uvar, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { isoeuler_vol_t kernels[3]; } gkyl_dg_isoeuler_vol_kern_list;
typedef struct { isoeuler_surf_t kernels[3]; } gkyl_dg_isoeuler_surf_kern_list;

//
// Serendipity basis kernels.
//

// Volume kernel list.
GKYL_CU_D
static const gkyl_dg_isoeuler_vol_kern_list ser_vol_kernels[] = { //TODO: rename all kernels using snake case and remove pOrder=3 cases
  // { NULL, NULL, NULL }, // 0
  // { NULL, NULL, NULL }, // 1
  // { NULL, NULL, NULL }, // 2
  { NULL, isoeuler_vol_1x1v_ser_p1, isoeuler_vol_1x1v_ser_p2 }, // 0
  { NULL, isoeuler_vol_2x2v_ser_p1, isoeuler_vol_2x2v_ser_p2 }, // 1
  { NULL, isoeuler_vol_3x3v_ser_p1, isoeuler_vol_3x3v_ser_p2 }, // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_isoeuler_surf_kern_list ser_surf_x_kernels[] = {
  // { NULL, NULL, NULL }, // 0
  // { NULL, NULL, NULL }, // 1
  // { NULL, NULL, NULL }, //
  { NULL, isoeuler_surfx_1x1v_ser_p1, isoeuler_surfx_1x1v_ser_p2 }, // 0
  { NULL, isoeuler_surfx_2x2v_ser_p1, isoeuler_surfx_2x2v_ser_p2 }, // 1
  { NULL, isoeuler_surfx_3x3v_ser_p1, isoeuler_surfx_3x3v_ser_p1 }, // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_isoeuler_surf_kern_list ser_surf_y_kernels[] = {
  // { NULL, NULL, NULL }, // 0
  // { NULL, NULL, NULL }, // 1
  // { NULL, NULL, NULL }, // 2
  { NULL, NULL, NULL }, // 0
  { NULL, isoeuler_surfy_2x2v_ser_p1, isoeuler_surfy_2x2v_ser_p2 }, // 1
  { NULL, isoeuler_surfy_3x3v_ser_p1, isoeuler_surfy_3x3v_ser_p2 }, // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_isoeuler_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, isoeuler_surfz_3x3v_ser_p1, isoeuler_surfz_3x3v_ser_p2 }, // 2
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim].kernels[poly_order]

struct dg_isoeuler {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  isoeuler_vol_t vol; // Volume kernel.
  isoeuler_surf_t surf[3]; // Surface terms.
  struct gkyl_range conf_range; // Configuration space range.
  double vth;
  struct gkyl_dg_isoeuler_auxfields auxfields; // Auxiliary fields.
};

/**
 * Free gyrokinetic eqn object.
 *
 * @param ref Reference counter for gyrokinetic eqn
 */
void gkyl_isoeuler_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double*  dx,
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_isoeuler *isoeuler = container_of(eqn, struct dg_isoeuler, eqn);
  long cidx = gkyl_range_idx(&isoeuler->conf_range, idx);
  return isoeuler->vol(xc, dx, isoeuler->vth,
    (const double*) gkyl_array_cfetch(isoeuler->auxfields.uvar, cidx), qIn, qRhsOut);
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
  struct dg_isoeuler *isoeuler = container_of(eqn, struct dg_isoeuler, eqn);

  long cidx_l = gkyl_range_idx(&isoeuler->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&isoeuler->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&isoeuler->conf_range, idxR);

  isoeuler->surf[dir](xcC, dxC, isoeuler->vth,
    (const double*) gkyl_array_cfetch(isoeuler->auxfields.uvar, cidx_l),
    (const double*) gkyl_array_cfetch(isoeuler->auxfields.uvar, cidx_c),
    (const double*) gkyl_array_cfetch(isoeuler->auxfields.uvar, cidx_r),
    qInL, qInC, qInR, qRhsOut);
}
