#pragma once

// Private header, not for direct use in user code

#include <gkyl_lbo_vlasov_pkpm_kernels.h>

// Types for various kernels
typedef void (*lbo_vlasov_pkpm_drag_surf_t)(const double *w, const double *dxv,
  const double *nu, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef void (*lbo_vlasov_pkpm_drag_boundary_surf_t)(const double *w, const double *dxv,
  const double *nu, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_lbo_vlasov_pkpm_drag_vol_kern_list;
typedef struct { lbo_vlasov_pkpm_drag_surf_t kernels[3]; } gkyl_dg_lbo_vlasov_pkpm_drag_surf_kern_list;
typedef struct { lbo_vlasov_pkpm_drag_boundary_surf_t kernels[3]; } gkyl_dg_lbo_vlasov_pkpm_drag_boundary_surf_kern_list;

struct dg_lbo_vlasov_pkpm_drag {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  lbo_vlasov_pkpm_drag_surf_t surf; // Surface terms for acceleration
  lbo_vlasov_pkpm_drag_boundary_surf_t boundary_surf; // Surface terms for acceleration
  struct gkyl_range conf_range; // Configuration space range.
  struct gkyl_dg_lbo_vlasov_pkpm_drag_auxfields auxfields; // Auxiliary fields.
  double vMaxSq;
};

GKYL_CU_DH
static inline bool
checkPrimMomCross(struct dg_lbo_vlasov_pkpm_drag *lbo_vlasov_pkpm_drag,
  const double* nu_p, const double* nuVtSq_p) {
  bool noPrimMomCross = true;
  noPrimMomCross = noPrimMomCross && ((nuVtSq_p[0]>0.)
    && (nuVtSq_p[0]/nu_p[0] < lbo_vlasov_pkpm_drag->vMaxSq));
  return noPrimMomCross;
}

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_lbo_vlasov_pkpm_drag_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_pkpm_drag *lbo_vlasov_pkpm_drag = container_of(eqn, struct dg_lbo_vlasov_pkpm_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_drag->conf_range, idx);

  const double* nu_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nu, cidx);
  const double* nuVtSq_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nuVtSq, cidx);
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_pkpm_drag, nu_p, nuVtSq_p);
  if (noPrimMomCross) {
    return lbo_vlasov_pkpm_drag_vol_1x1v_ser_p1(xc, dx, nu_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_pkpm_drag_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_pkpm_drag *lbo_vlasov_pkpm_drag = container_of(eqn, struct dg_lbo_vlasov_pkpm_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_drag->conf_range, idx);

  const double* nu_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nu, cidx);
  const double* nuVtSq_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nuVtSq, cidx);
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_pkpm_drag, nu_p, nuVtSq_p);
  if (noPrimMomCross) {
    return lbo_vlasov_pkpm_drag_vol_1x1v_ser_p2(xc, dx, nu_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_pkpm_drag_vol_1x1v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_pkpm_drag *lbo_vlasov_pkpm_drag = container_of(eqn, struct dg_lbo_vlasov_pkpm_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_drag->conf_range, idx);

  const double* nu_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nu, cidx);
  const double* nuVtSq_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nuVtSq, cidx);
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_pkpm_drag, nu_p, nuVtSq_p);
  if (noPrimMomCross) {
    return lbo_vlasov_pkpm_drag_vol_1x1v_tensor_p2(xc, dx, nu_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_pkpm_drag_vol_2x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_pkpm_drag *lbo_vlasov_pkpm_drag = container_of(eqn, struct dg_lbo_vlasov_pkpm_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_drag->conf_range, idx);

  const double* nu_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nu, cidx);
  const double* nuVtSq_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nuVtSq, cidx);
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_pkpm_drag, nu_p, nuVtSq_p);
  if (noPrimMomCross) {
    return lbo_vlasov_pkpm_drag_vol_2x1v_ser_p1(xc, dx, nu_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_pkpm_drag_vol_2x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_pkpm_drag *lbo_vlasov_pkpm_drag = container_of(eqn, struct dg_lbo_vlasov_pkpm_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_drag->conf_range, idx);

  const double* nu_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nu, cidx);
  const double* nuVtSq_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nuVtSq, cidx);
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_pkpm_drag, nu_p, nuVtSq_p);
  if (noPrimMomCross) {
    return lbo_vlasov_pkpm_drag_vol_2x1v_ser_p2(xc, dx, nu_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_vlasov_pkpm_drag_vol_3x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_pkpm_drag *lbo_vlasov_pkpm_drag = container_of(eqn, struct dg_lbo_vlasov_pkpm_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_drag->conf_range, idx);

  const double* nu_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nu, cidx);
  const double* nuVtSq_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nuVtSq, cidx);
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_pkpm_drag, nu_p, nuVtSq_p);
  if (noPrimMomCross) {
    return lbo_vlasov_pkpm_drag_vol_3x1v_ser_p1(xc, dx, nu_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

// Volume kernel list (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_pkpm_drag_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_lbo_vlasov_pkpm_drag_vol_1x1v_ser_p1, kernel_lbo_vlasov_pkpm_drag_vol_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, kernel_lbo_vlasov_pkpm_drag_vol_2x1v_ser_p1, kernel_lbo_vlasov_pkpm_drag_vol_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, kernel_lbo_vlasov_pkpm_drag_vol_3x1v_ser_p1, NULL }, // 2
};

// Volume kernel list (Tensor basis)
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_pkpm_drag_vol_kern_list ten_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_lbo_vlasov_pkpm_drag_vol_1x1v_ser_p1, kernel_lbo_vlasov_pkpm_drag_vol_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, kernel_lbo_vlasov_pkpm_drag_vol_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, kernel_lbo_vlasov_pkpm_drag_vol_3x1v_ser_p1, NULL }, // 2
};


// Constant nu surface kernel list: vpar-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_pkpm_drag_surf_kern_list ser_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_vlasov_pkpm_drag_surfvpar_1x1v_ser_p1, lbo_vlasov_pkpm_drag_surfvpar_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, lbo_vlasov_pkpm_drag_surfvpar_2x1v_ser_p1, lbo_vlasov_pkpm_drag_surfvpar_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, lbo_vlasov_pkpm_drag_surfvpar_3x1v_ser_p1, NULL }, // 2
};

// Constant nu surface kernel list: vpar-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_pkpm_drag_surf_kern_list ten_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_vlasov_pkpm_drag_surfvpar_1x1v_ser_p1, lbo_vlasov_pkpm_drag_surfvpar_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, lbo_vlasov_pkpm_drag_surfvpar_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, lbo_vlasov_pkpm_drag_surfvpar_3x1v_ser_p1, NULL }, // 2
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vpar-direction (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_pkpm_drag_boundary_surf_kern_list ser_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_vlasov_pkpm_drag_boundary_surfvpar_1x1v_ser_p1, lbo_vlasov_pkpm_drag_boundary_surfvpar_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, lbo_vlasov_pkpm_drag_boundary_surfvpar_2x1v_ser_p1, lbo_vlasov_pkpm_drag_boundary_surfvpar_2x1v_ser_p2 }, // 1
  // 3x kernels
  { NULL, lbo_vlasov_pkpm_drag_boundary_surfvpar_3x1v_ser_p1, NULL }, // 2
};

// Constant nu boundary surface kernel (zero-flux BCs) list: vpar-direction (Tensor basis)
GKYL_CU_D
static const gkyl_dg_lbo_vlasov_pkpm_drag_boundary_surf_kern_list ten_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_vlasov_pkpm_drag_boundary_surfvpar_1x1v_ser_p1, lbo_vlasov_pkpm_drag_boundary_surfvpar_1x1v_tensor_p2 }, // 0
  // 2x kernels
  { NULL, lbo_vlasov_pkpm_drag_boundary_surfvpar_2x1v_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, lbo_vlasov_pkpm_drag_boundary_surfvpar_3x1v_ser_p1, NULL }, // 2
};

void gkyl_lbo_vlasov_pkpm_drag_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_vlasov_pkpm_drag *lbo_vlasov_pkpm_drag = container_of(eqn, struct dg_lbo_vlasov_pkpm_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_drag->conf_range, idxC);

  const double* nu_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nu, cidx);
  const double* nuVtSq_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nuVtSq, cidx);
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_pkpm_drag, nu_p, nuVtSq_p);
  if ((dir >= lbo_vlasov_pkpm_drag->cdim) && (noPrimMomCross)) {
    lbo_vlasov_pkpm_drag->surf(xcC, dxC, 
      nu_p, qInL, qInC, qInR, qRhsOut);
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
  struct dg_lbo_vlasov_pkpm_drag *lbo_vlasov_pkpm_drag = container_of(eqn, struct dg_lbo_vlasov_pkpm_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_vlasov_pkpm_drag->conf_range, idxSkin);

  const double* nu_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nu, cidx);
  const double* nuVtSq_p = (const double*) gkyl_array_cfetch(lbo_vlasov_pkpm_drag->auxfields.nuVtSq, cidx);
  bool noPrimMomCross = checkPrimMomCross(lbo_vlasov_pkpm_drag, nu_p, nuVtSq_p);
  if ((dir >= lbo_vlasov_pkpm_drag->cdim) && (noPrimMomCross)) {
    lbo_vlasov_pkpm_drag->boundary_surf(xcSkin, dxSkin, 
      nu_p, edge, qInSkin, qInEdge, qRhsOut);
  }
}

