#pragma once

// Private header for the gyrokinetic LBO drag equation object.
// Not for direct use in user code.

#include <gkyl_lbo_gyrokinetic_kernels.h>

// Types for various kernels
typedef double (*lbo_gyrokinetic_drag_surf_t)(const double *w, const double *dxv, const double m_,
  const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*lbo_gyrokinetic_drag_boundary_surf_t)(const double *w, const double *dxv, const double m_,
  const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum,
  const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_lbo_gyrokinetic_drag_vol_kern_list;
typedef struct { lbo_gyrokinetic_drag_surf_t kernels[3]; } gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list;
typedef struct { lbo_gyrokinetic_drag_boundary_surf_t kernels[3]; } gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list;

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst, cdim, vd, poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

struct dg_lbo_gyrokinetic_drag {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  lbo_gyrokinetic_drag_surf_t surf[2]; // Surface terms for acceleration.
  lbo_gyrokinetic_drag_boundary_surf_t boundary_surf[2]; // Surface terms for acceleration.
  struct gkyl_range conf_range; // Configuration space range.
  double mass; // Species mass.
  struct gkyl_dg_lbo_gyrokinetic_drag_auxfields auxfields; // Auxiliary fields.
  double vparMax, vparMaxSq;
  int num_cbasis;
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.)) {
    return lbo_gyrokinetic_drag_vol_1x1v_ser_p1(xc, dx, lbo_gyrokinetic_drag->mass, 
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.)) {
    return lbo_gyrokinetic_drag_vol_1x1v_ser_p2(xc, dx, lbo_gyrokinetic_drag->mass, 
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.)) {
  return lbo_gyrokinetic_drag_vol_1x2v_ser_p1(xc, dx, lbo_gyrokinetic_drag->mass, 
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.)) {
  return lbo_gyrokinetic_drag_vol_1x2v_ser_p2(xc, dx, lbo_gyrokinetic_drag->mass, 
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.)) {
  return lbo_gyrokinetic_drag_vol_2x2v_ser_p1(xc, dx, lbo_gyrokinetic_drag->mass, 
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.)) {
  return lbo_gyrokinetic_drag_vol_2x2v_ser_p2(xc, dx, lbo_gyrokinetic_drag->mass, 
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_3x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.)) {
  return lbo_gyrokinetic_drag_vol_3x2v_ser_p1(xc, dx, lbo_gyrokinetic_drag->mass, 
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

GKYL_CU_DH
static double
kernel_lbo_gyrokinetic_drag_vol_3x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idx);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.)) {
  return lbo_gyrokinetic_drag_vol_3x2v_ser_p2(xc, dx, lbo_gyrokinetic_drag->mass, 
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qIn, qRhsOut);
  } else {
    return 0.;
  }
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_lbo_gyrokinetic_drag_vol_1x1v_ser_p1, kernel_lbo_gyrokinetic_drag_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_lbo_gyrokinetic_drag_vol_1x2v_ser_p1, kernel_lbo_gyrokinetic_drag_vol_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_lbo_gyrokinetic_drag_vol_2x2v_ser_p1, kernel_lbo_gyrokinetic_drag_vol_2x2v_ser_p2 }, // 3
  // 3x kernels
  { NULL, kernel_lbo_gyrokinetic_drag_vol_3x2v_ser_p1, kernel_lbo_gyrokinetic_drag_vol_3x2v_ser_p2 }, // 4
};

// Surface kernel list: vpar-direction
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list ser_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_gyrokinetic_drag_surfvpar_1x1v_ser_p1, lbo_gyrokinetic_drag_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, lbo_gyrokinetic_drag_surfvpar_1x2v_ser_p1, lbo_gyrokinetic_drag_surfvpar_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, lbo_gyrokinetic_drag_surfvpar_2x2v_ser_p1, lbo_gyrokinetic_drag_surfvpar_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, lbo_gyrokinetic_drag_surfvpar_3x2v_ser_p1, lbo_gyrokinetic_drag_surfvpar_3x2v_ser_p2 }, // 3
};

// Surface kernel list: mu-direction
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list ser_surf_mu_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, lbo_gyrokinetic_drag_surfmu_1x2v_ser_p1, lbo_gyrokinetic_drag_surfmu_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, lbo_gyrokinetic_drag_surfmu_2x2v_ser_p1, lbo_gyrokinetic_drag_surfmu_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, lbo_gyrokinetic_drag_surfmu_3x2v_ser_p1, lbo_gyrokinetic_drag_surfmu_3x2v_ser_p2 }, // 3
};

// Boundary surface kernel (zero-flux BCs) list: vpar-direction
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list ser_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p1, lbo_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, lbo_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p1, lbo_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfvpar_2x2v_ser_p1, lbo_gyrokinetic_drag_boundary_surfvpar_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p1, lbo_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p2 }, // 3
};

// Constant nu boundary surface kernel (zero-flux BCs) list: mu-direction
GKYL_CU_D
static const gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list ser_boundary_surf_mu_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, lbo_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p1, lbo_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfmu_2x2v_ser_p1, lbo_gyrokinetic_drag_boundary_surfmu_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, lbo_gyrokinetic_drag_boundary_surfmu_3x2v_ser_p1, lbo_gyrokinetic_drag_boundary_surfmu_3x2v_ser_p2 }, // 3
};

void gkyl_lbo_gyrokinetic_drag_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idxC);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx);
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((dir >= lbo_gyrokinetic_drag->cdim) &&
      (fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.))
  {
    return lbo_gyrokinetic_drag->surf[dir-lbo_gyrokinetic_drag->cdim](xcC, dxC, lbo_gyrokinetic_drag->mass,
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, qInL, qInC, qInR, qRhsOut);
  }
  return 0.;
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&lbo_gyrokinetic_drag->conf_range, idxSkin);
  const double* nuSum_p     = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuSum, cidx);
  const double* nuPrimMomsSum_p = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum, cidx); 
  const double* m2self_p    = (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.m2self, cidx);
  const double* nuUSum_p    = nuPrimMomsSum_p;
  const double* nuVtSqSum_p = &nuPrimMomsSum_p[lbo_gyrokinetic_drag->num_cbasis];
  if ((dir >= lbo_gyrokinetic_drag->cdim) &&
      (fabs(nuUSum_p[0]/nuSum_p[0]) < lbo_gyrokinetic_drag->vparMax) &&
      (nuVtSqSum_p[0]>0.) && (nuVtSqSum_p[0]/nuSum_p[0] < lbo_gyrokinetic_drag->vparMaxSq) &&
      (m2self_p[0]>0.))
  {
    return lbo_gyrokinetic_drag->boundary_surf[dir-lbo_gyrokinetic_drag->cdim](xcSkin, dxSkin, 
      lbo_gyrokinetic_drag->mass,
      (const double*) gkyl_array_cfetch(lbo_gyrokinetic_drag->auxfields.bmag_inv, cidx), 
      nuSum_p, nuPrimMomsSum_p, edge, qInSkin, qInEdge, qRhsOut);
  }
  return 0.;
}

