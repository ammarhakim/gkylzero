#pragma once

// Private header, not for direct use in user code

#include <gkyl_lbo_vlasov_kernels.h>

// Types for various kernels
typedef double (*lbo_vlasov_drag_vol_t)(const double *w, const double *dxv,
  const double *nuSum, const double *nuUSum, const double *nuVtSqSum,
  const double *f, double* GKYL_RESTRICT out);

typedef void (*lbo_vlasov_drag_surf_t)(const double *w, const double *dxv,
  const double *nuSum, const double *nuUSum, const double *nuVtSqSum, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef void (*lbo_vlasov_drag_boundary_surf_t)(const double *w, const double *dxv,
  const double *nuSum, const double *nuUSum, const double *nuVtSqSum, 
  const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { lbo_vlasov_drag_vol_t kernels[3]; } gkyl_dg_lbo_vlasov_drag_vol_kern_list;
typedef struct { lbo_vlasov_drag_surf_t kernels[3]; } gkyl_dg_lbo_vlasov_drag_surf_kern_list;
typedef struct { lbo_vlasov_drag_boundary_surf_t kernels[3]; } gkyl_dg_lbo_vlasov_drag_boundary_surf_kern_list;

//
// Serendipity basis kernels
//

GKYL_CU_D
static const gkyl_dg_lbo_vlasov_drag_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, lbo_vlasov_drag_vol_1x1v_ser_p1, lbo_vlasov_drag_vol_1x1v_ser_p2 }, // 0
  { NULL, lbo_vlasov_drag_vol_1x2v_ser_p1, lbo_vlasov_drag_vol_1x2v_ser_p2 }, // 1
  { NULL, lbo_vlasov_drag_vol_1x3v_ser_p1, lbo_vlasov_drag_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, lbo_vlasov_drag_vol_2x2v_ser_p1, lbo_vlasov_drag_vol_2x2v_ser_p2 }, // 3
  { NULL, lbo_vlasov_drag_vol_2x3v_ser_p1, NULL               }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_drag_vol_3x3v_ser_p1, NULL               }, // 5
};

// Constant nu surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_drag_surf_kern_list ser_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, lbo_vlasov_drag_surfvx_1x1v_ser_p1, lbo_vlasov_drag_surfvx_1x1v_ser_p2 }, // 0
  { NULL, lbo_vlasov_drag_surfvx_1x2v_ser_p1, lbo_vlasov_drag_surfvx_1x2v_ser_p2 }, // 1
  { NULL, lbo_vlasov_drag_surfvx_1x3v_ser_p1, lbo_vlasov_drag_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, lbo_vlasov_drag_surfvx_2x2v_ser_p1, lbo_vlasov_drag_surfvx_2x2v_ser_p2 }, // 3
  { NULL, lbo_vlasov_drag_surfvx_2x3v_ser_p1, NULL }, // 
  // 3x kernels
  { NULL, lbo_vlasov_drag_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Constant nu surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_drag_surf_kern_list ser_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, lbo_vlasov_drag_surfvy_1x2v_ser_p1, lbo_vlasov_drag_surfvy_1x2v_ser_p2 }, // 1
  { NULL, lbo_vlasov_drag_surfvy_1x3v_ser_p1, lbo_vlasov_drag_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, lbo_vlasov_drag_surfvy_2x2v_ser_p1, lbo_vlasov_drag_surfvy_2x2v_ser_p2 }, // 3
  { NULL, lbo_vlasov_drag_surfvy_2x3v_ser_p1, NULL }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_drag_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Constant nu surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_drag_surf_kern_list ser_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, lbo_vlasov_drag_surfvz_1x3v_ser_p1, lbo_vlasov_drag_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, lbo_vlasov_drag_surfvz_2x3v_ser_p1, NULL }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_drag_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_drag_boundary_surf_kern_list ser_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, lbo_vlasov_drag_boundary_surfvx_1x1v_ser_p1, lbo_vlasov_drag_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, lbo_vlasov_drag_boundary_surfvx_1x2v_ser_p1, lbo_vlasov_drag_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, lbo_vlasov_drag_boundary_surfvx_1x3v_ser_p1, lbo_vlasov_drag_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, lbo_vlasov_drag_boundary_surfvx_2x2v_ser_p1, lbo_vlasov_drag_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, lbo_vlasov_drag_boundary_surfvx_2x3v_ser_p1, NULL }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_drag_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_drag_boundary_surf_kern_list ser_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, lbo_vlasov_drag_boundary_surfvy_1x2v_ser_p1, lbo_vlasov_drag_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, lbo_vlasov_drag_boundary_surfvy_1x3v_ser_p1, lbo_vlasov_drag_boundary_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, lbo_vlasov_drag_boundary_surfvy_2x2v_ser_p1, lbo_vlasov_drag_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, lbo_vlasov_drag_boundary_surfvy_2x3v_ser_p1, NULL }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_drag_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_drag_boundary_surf_kern_list ser_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, lbo_vlasov_drag_boundary_surfvz_1x3v_ser_p1, lbo_vlasov_drag_boundary_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, lbo_vlasov_drag_boundary_surfvz_2x3v_ser_p1, NULL }, // 4
  // 3x kernels
  { NULL, lbo_vlasov_drag_boundary_surfvz_3x3v_ser_p1, NULL }, // 5
};

//
// Tensor-product basis kernels
//

/* GKYL_CU_D */
/* static const gkyl_dg_lbo_vlasov_drag_vol_kern_list ten_vol_kernels[] = { */
/*   // 1x kernels */
/*   { NULL, lbo_vlasov_drag_vol_1x1v_ser_p1, lbo_vlasov_drag_vol_1x1v_tensor_p2 }, // 0 */
/*   { NULL, lbo_vlasov_drag_vol_1x2v_ser_p1, lbo_vlasov_drag_vol_1x2v_tensor_p2 }, // 1 */
/*   { NULL, lbo_vlasov_drag_vol_1x3v_ser_p1, lbo_vlasov_drag_vol_1x3v_tensor_p2 }, // 2 */
/*   // 2x kernels */
/*   { NULL, lbo_vlasov_drag_vol_2x2v_ser_p1, lbo_vlasov_drag_vol_2x2v_tensor_p2 }, // 3 */
/*   { NULL, lbo_vlasov_drag_vol_2x3v_ser_p1, NULL               }, // 4 */
/*   // 3x kernels */
/*   { NULL, lbo_vlasov_drag_vol_3x3v_ser_p1, NULL               }, // 5 */
/* }; */

/* // Constant nu surface kernel list: vx-direction */
/* GKYL_CU_D */
/* static const gkyl_dg_lbo_vlasov_drag_surf_kern_list ten_surf_vx_kernels[] = { */
/*   // 1x kernels */
/*   { NULL, lbo_vlasov_drag_surfvx_1x1v_ser_p1, lbo_vlasov_drag_surfvx_1x1v_tensor_p2 }, // 0 */
/*   { NULL, lbo_vlasov_drag_surfvx_1x2v_ser_p1, lbo_vlasov_drag_surfvx_1x2v_tensor_p2 }, // 1 */
/*   { NULL, lbo_vlasov_drag_surfvx_1x3v_ser_p1, lbo_vlasov_drag_surfvx_1x3v_tensor_p2 }, // 2 */
/*   // 2x kernels */
/*   { NULL, lbo_vlasov_drag_surfvx_2x2v_ser_p1, lbo_vlasov_drag_surfvx_2x2v_tensor_p2 }, // 3 */
/*   { NULL, lbo_vlasov_drag_surfvx_2x3v_ser_p1, NULL }, //  */
/*   // 3x kernels */
/*   { NULL, lbo_vlasov_drag_surfvx_3x3v_ser_p1, NULL                   }, // 5 */
/* }; */

/* // Constant nu surface kernel list: vy-direction */
/* GKYL_CU_D */
/* static const gkyl_dg_lbo_vlasov_drag_surf_kern_list ten_surf_vy_kernels[] = { */
/*   // 1x kernels */
/*   { NULL, NULL, NULL }, // 0 */
/*   { NULL, lbo_vlasov_drag_surfvy_1x2v_ser_p1, lbo_vlasov_drag_surfvy_1x2v_tensor_p2 }, // 1 */
/*   { NULL, lbo_vlasov_drag_surfvy_1x3v_ser_p1, lbo_vlasov_drag_surfvy_1x3v_tensor_p2 }, // 2 */
/*   // 2x kernels */
/*   { NULL, lbo_vlasov_drag_surfvy_2x2v_ser_p1, lbo_vlasov_drag_surfvy_2x2v_tensor_p2 }, // 3 */
/*   { NULL, lbo_vlasov_drag_surfvy_2x3v_ser_p1, NULL }, // 4 */
/*   // 3x kernels */
/*   { NULL, lbo_vlasov_drag_surfvy_3x3v_ser_p1, NULL                   }, // 5 */
/* }; */

/* // Constant nu surface kernel list: vz-direction */
/* GKYL_CU_D */
/* static const gkyl_dg_lbo_vlasov_drag_surf_kern_list ten_surf_vz_kernels[] = { */
/*   // 1x kernels */
/*   { NULL, NULL, NULL }, // 0 */
/*   { NULL, NULL, NULL }, // 1 */
/*   { NULL, lbo_vlasov_drag_surfvz_1x3v_ser_p1, lbo_vlasov_drag_surfvz_1x3v_tensor_p2 }, // 2 */
/*   // 2x kernels */
/*   { NULL, NULL, NULL }, // 3 */
/*   { NULL, lbo_vlasov_drag_surfvz_2x3v_ser_p1, NULL }, // 4 */
/*   // 3x kernels */
/*   { NULL, lbo_vlasov_drag_surfvz_3x3v_ser_p1, NULL }, // 5 */
/* }; */

/* // Constant nu boundary surface kernel (zero-flux BCs) list: vx-direction */
/* GKYL_CU_D */
/* static const gkyl_dg_lbo_vlasov_drag_boundary_surf_kern_list ten_boundary_surf_vx_kernels[] = { */
/*   // 1x kernels */
/*   { NULL, lbo_vlasov_drag_boundary_surfvx_1x1v_ser_p1, lbo_vlasov_drag_boundary_surfvx_1x1v_tensor_p2 }, // 0 */
/*   { NULL, lbo_vlasov_drag_boundary_surfvx_1x2v_ser_p1, lbo_vlasov_drag_boundary_surfvx_1x2v_tensor_p2 }, // 1 */
/*   { NULL, lbo_vlasov_drag_boundary_surfvx_1x3v_ser_p1, lbo_vlasov_drag_boundary_surfvx_1x3v_tensor_p2 }, // 2 */
/*   // 2x kernels */
/*   { NULL, lbo_vlasov_drag_boundary_surfvx_2x2v_ser_p1, lbo_vlasov_drag_boundary_surfvx_2x2v_tensor_p2 }, // 3 */
/*   { NULL, lbo_vlasov_drag_boundary_surfvx_2x3v_ser_p1, NULL }, // 4 */
/*   // 3x kernels */
/*   { NULL, lbo_vlasov_drag_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5 */
/* }; */

/* // Constant nu boundary surface kernel (zero-flux BCs) list: vy-direction */
/* GKYL_CU_D */
/* static const gkyl_dg_lbo_vlasov_drag_boundary_surf_kern_list ten_boundary_surf_vy_kernels[] = { */
/*   // 1x kernels */
/*   { NULL, NULL, NULL }, // 0 */
/*   { NULL, lbo_vlasov_drag_boundary_surfvy_1x2v_ser_p1, lbo_vlasov_drag_boundary_surfvy_1x2v_tensor_p2 }, // 1 */
/*   { NULL, lbo_vlasov_drag_boundary_surfvy_1x3v_ser_p1, lbo_vlasov_drag_boundary_surfvy_1x3v_tensor_p2 }, // 2 */
/*   // 2x kernels */
/*   { NULL, lbo_vlasov_drag_boundary_surfvy_2x2v_ser_p1, lbo_vlasov_drag_boundary_surfvy_2x2v_tensor_p2 }, // 3 */
/*   { NULL, lbo_vlasov_drag_boundary_surfvy_2x3v_ser_p1, NULL }, // 4 */
/*   // 3x kernels */
/*   { NULL, lbo_vlasov_drag_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5 */
/* }; */

/* // Constant nu boundary surface kernel (zero-flux BCs) list: vz-direction */
/* GKYL_CU_D */
/* static const gkyl_dg_lbo_vlasov_drag_boundary_surf_kern_list ten_boundary_surf_vz_kernels[] = { */
/*   // 1x kernels */
/*   { NULL, NULL, NULL }, // 0 */
/*   { NULL, NULL, NULL }, // 1 */
/*   { NULL, lbo_vlasov_drag_boundary_surfvz_1x3v_ser_p1, lbo_vlasov_drag_boundary_surfvz_1x3v_tensor_p2 }, // 2 */
/*   // 2x kernels */
/*   { NULL, NULL, NULL }, // 3 */
/*   { NULL, lbo_vlasov_drag_boundary_surfvz_2x3v_ser_p1, NULL }, // 4 */
/*   // 3x kernels */
/*   { NULL, lbo_vlasov_drag_boundary_surfvz_3x3v_ser_p1, NULL }, // 5 */
/* }; */

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst, cdim, vd, poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

struct dg_lbo_vlasov_drag {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  lbo_vlasov_drag_vol_t vol; // Volume kernel
  lbo_vlasov_drag_surf_t surf[3]; // Surface terms for acceleration
  lbo_vlasov_drag_boundary_surf_t boundary_surf[3]; // Surface terms for acceleration
  struct gkyl_range conf_range; // Configuration space range.
  struct gkyl_dg_lbo_vlasov_drag_auxfields auxfields; // Auxiliary fields.
};

void gkyl_lbo_vlasov_drag_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_drag *lbo_vlasov_drag = container_of(eqn, struct dg_lbo_vlasov_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_drag->conf_range, idx);
  return lbo_vlasov_drag->vol(xc, dx, 
    (const double*) gkyl_array_cfetch(lbo_vlasov_drag->auxfields.nuSum, cidx), 
    (const double*) gkyl_array_cfetch(lbo_vlasov_drag->auxfields.nuUSum, cidx), 
    (const double*) gkyl_array_cfetch(lbo_vlasov_drag->auxfields.nuVtSqSum, cidx), 
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
  struct dg_lbo_vlasov_drag *lbo_vlasov_drag = container_of(eqn, struct dg_lbo_vlasov_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_drag->conf_range, idxC);
  if (dir >= lbo_vlasov_drag->cdim) {
    lbo_vlasov_drag->surf[dir-lbo_vlasov_drag->cdim](xcC, dxC, 
      (const double*) gkyl_array_cfetch(lbo_vlasov_drag->auxfields.nuSum, cidx), 
      (const double*) gkyl_array_cfetch(lbo_vlasov_drag->auxfields.nuUSum, cidx), 
      (const double*) gkyl_array_cfetch(lbo_vlasov_drag->auxfields.nuVtSqSum, cidx), 
      qInL, qInC, qInR, qRhsOut);
  }
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
  struct dg_lbo_vlasov_drag *lbo_vlasov_drag = container_of(eqn, struct dg_lbo_vlasov_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_drag->conf_range, idxSkin);
  if (dir >= lbo_vlasov_drag->cdim) {
    lbo_vlasov_drag->boundary_surf[dir-lbo_vlasov_drag->cdim](xcSkin, dxSkin, 
      (const double*) gkyl_array_cfetch(lbo_vlasov_drag->auxfields.nuSum, cidx), 
      (const double*) gkyl_array_cfetch(lbo_vlasov_drag->auxfields.nuUSum, cidx), 
      (const double*) gkyl_array_cfetch(lbo_vlasov_drag->auxfields.nuVtSqSum, cidx),  
      edge, qInSkin, qInEdge, qRhsOut);
  }
}

