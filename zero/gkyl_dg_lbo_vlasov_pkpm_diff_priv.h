#pragma once

// Private header, not for direct use in user code

#include <gkyl_lbo_vlasov_kernels.h>

// Types for various kernels
typedef double (*lbo_vlasov_pkpm_diff_vol_t)(const double *w, const double *dxv,
  const double *nuVtSq, const double *f, double* GKYL_RESTRICT out);

typedef void (*lbo_vlasov_pkpm_diff_surf_t)(const double *w, const double *dxv,
  const double *nuVtSq, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef void (*lbo_vlasov_pkpm_diff_boundary_surf_t)(const double *w, const double *dxv,
  const double *nuVtSqSum,  const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);


// for use in kernel tables
typedef struct { lbo_vlasov_pkpm_diff_vol_t kernels[3]; } gkyl_dg_lbo_vlasov_pkpm_diff_vol_kern_list;
typedef struct { lbo_vlasov_pkpm_diff_surf_t kernels[3]; } gkyl_dg_lbo_vlasov_pkpm_diff_surf_kern_list;
typedef struct { lbo_vlasov_pkpm_diff_boundary_surf_t kernels[3]; } gkyl_dg_lbo_vlasov_pkpm_diff_boundary_surf_kern_list;

//
// Serendipity basis kernels
//

GKYL_CU_D
static const gkyl_dg_lbo_vlasov_pkpm_diff_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, lbo_vlasov_pkpm_diff_vol_1x1v_ser_p2 }, // 0
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  // // 2x kernels
  // { NULL, lbo_vlasov_pkpm_diff_vol_2x1v_ser_p1, lbo_vlasov_pkpm_diff_vol_2x1v_ser_p2 }, // 1
  // // 3x kernels
  // { NULL, lbo_vlasov_pkpm_diff_vol_3x1v_ser_p1, lbo_vlasov_pkpm_diff_vol_3x1v_ser_p2 }, // 2
};

// Constant nu surface kernel list: vpar-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_pkpm_diff_surf_kern_list ser_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, NULL, lbo_vlasov_pkpm_diff_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  // // 2x kernels
  // { NULL, lbo_vlasov_pkpm_diff_surfvpar_2x1v_ser_p1, lbo_vlasov_pkpm_diff_surfvpar_2x1v_ser_p2 }, // 1
  // // 3x kernels
  // { NULL, lbo_vlasov_pkpm_diff_surfvpar_3x1v_ser_p1, lbo_vlasov_pkpm_diff_surfvpar_3x1v_ser_p2 }, // 2
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vpar-direction
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_pkpm_diff_boundary_surf_kern_list ser_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, NULL, lbo_vlasov_pkpm_diff_boundary_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  // // 2x kernels
  // { NULL, lbo_vlasov_pkpm_diff_boundary_surfvpar_2x1v_ser_p1, lbo_vlasov_pkpm_diff_boundary_surfvpar_2x1v_ser_p2 }, // 1
  // // 3x kernels
  // { NULL, lbo_vlasov_pkpm_diff_boundary_surfvpar_3x1v_ser_p1, lbo_vlasov_pkpm_diff_boundary_surfvpar_3x1v_ser_p2 }, // 2
};

struct dg_lbo_vlasov_pkpm_diff {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  lbo_vlasov_pkpm_diff_vol_t vol; // Volume kernel
  lbo_vlasov_pkpm_diff_surf_t surf; // Surface terms for acceleration
  lbo_vlasov_pkpm_diff_boundary_surf_t boundary_surf; // Surface terms for acceleration
  struct gkyl_range conf_range; // Configuration space range.
  struct gkyl_dg_lbo_vlasov_pkpm_diff_auxfields auxfields; // Auxiliary fields.
};

void gkyl_lbo_vlasov_pkpm_diff_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_pkpm_diff *lbo_vlasov_pkpm_diff = container_of(eqn, struct dg_lbo_vlasov_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_diff->conf_range, idx);
  return lbo_vlasov_pkpm_diff->vol(xc, dx, 
    (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_diff->auxfields.nuVtSq, cidx), 
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
  struct dg_lbo_vlasov_pkpm_diff *lbo_vlasov_pkpm_diff = container_of(eqn, struct dg_lbo_vlasov_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_diff->conf_range, idxC);
  if (dir >= lbo_vlasov_pkpm_diff->cdim) {
    lbo_vlasov_pkpm_diff->surf(xcC, dxC, 
      (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_diff->auxfields.nuVtSq, cidx), 
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
  struct dg_lbo_vlasov_pkpm_diff *lbo_vlasov_pkpm_diff = container_of(eqn, struct dg_lbo_vlasov_pkpm_diff, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_diff->conf_range, idxSkin);
  if (dir >= lbo_vlasov_pkpm_diff->cdim) {
    lbo_vlasov_pkpm_diff->boundary_surf(xcSkin, dxSkin, 
      (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_diff->auxfields.nuVtSq, cidx),  
      edge, qInSkin, qInEdge, qRhsOut);
  }
}

