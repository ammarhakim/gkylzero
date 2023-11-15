#pragma once

// Private header for the gyrokinetic RAD drag equation object.
// Not for direct use in user code.

#include <gkyl_rad_gyrokinetic_kernels.h>

// Types for various kernels
typedef double (*rad_gyrokinetic_drag_surf_t)(const double *w, const double *dxv, const double *nI,
  const double *nuField, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*rad_gyrokinetic_drag_boundary_surf_t)(const double *w, const double *dxv,
  const double *nI, const double *nuField, const int edge, const double *fEdge, const double *fSkin,
  double* GKYL_RESTRICT out);

// Struct containing data for vnu calculation
struct gkyl_dg_rad_vnu_params {
  const struct gkyl_array *bmag;
  const struct gkyl_array *fit_params;
  const struct gkyl_rect_grid *pgrid;
  const struct gkyl_range *conf_range;
  int cdim, vdim;
};
/*
struct gkyl_dg_rad_vnu_params* gkyl_dg_rad_vnu_params_new(struct gkyl_array *bmag, struct gkyl_array *fit_params, int cdim, int vdim) {
  struct gkyl_dg_rad_vnu_params* vnu_params = gkyl_malloc(sizeof(struct gkyl_dg_rad_vnu_params));

  vnu_params->cdim = cdim;
  vnu_params->vdim = vdim;
  vnu_params->bmag = bmag;
  vnu_params->fit_params = fit_params;

  return vnu_params;
  }*/

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_rad_gyrokinetic_drag_vol_kern_list;
typedef struct { rad_gyrokinetic_drag_surf_t kernels[3]; } gkyl_dg_rad_gyrokinetic_drag_surf_kern_list;
typedef struct { rad_gyrokinetic_drag_boundary_surf_t kernels[3]; } gkyl_dg_rad_gyrokinetic_drag_boundary_surf_kern_list;

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst, cdim, vd, poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

struct dg_rad_gyrokinetic_drag {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  rad_gyrokinetic_drag_surf_t surf[2]; // Surface terms for acceleration.
  rad_gyrokinetic_drag_boundary_surf_t boundary_surf[2]; // Surface terms for acceleration.
  struct gkyl_range conf_range; // Configuration space range.
  struct gkyl_dg_rad_gyrokinetic_drag_auxfields auxfields; // Auxiliary fields.
  double vparMax, vparMaxSq;
  int num_cbasis;
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_drag_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idx);
  return rad_gyrokinetic_drag_vol_1x1v_ser_p1(xc, dx,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_drag_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idx);
  return rad_gyrokinetic_drag_vol_1x1v_ser_p1(xc, dx,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_drag_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idx);
  return rad_gyrokinetic_drag_vol_1x1v_ser_p1(xc, dx,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_drag_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idx);
  return rad_gyrokinetic_drag_vol_1x1v_ser_p1(xc, dx,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_drag_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idx);
  return rad_gyrokinetic_drag_vol_1x1v_ser_p1(xc, dx,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_drag_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
    struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idx);
  return rad_gyrokinetic_drag_vol_1x1v_ser_p1(xc, dx,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_drag_vol_3x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idx);
  return rad_gyrokinetic_drag_vol_1x1v_ser_p1(xc, dx,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx), qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_rad_gyrokinetic_drag_vol_3x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idx);
  return rad_gyrokinetic_drag_vol_1x1v_ser_p1(xc, dx,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx), 
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx), qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_drag_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_rad_gyrokinetic_drag_vol_1x1v_ser_p1, kernel_rad_gyrokinetic_drag_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_rad_gyrokinetic_drag_vol_1x2v_ser_p1, kernel_rad_gyrokinetic_drag_vol_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_rad_gyrokinetic_drag_vol_2x2v_ser_p1, kernel_rad_gyrokinetic_drag_vol_2x2v_ser_p2 }, // 3
  // 3x kernels
  { NULL, kernel_rad_gyrokinetic_drag_vol_3x2v_ser_p1, kernel_rad_gyrokinetic_drag_vol_3x2v_ser_p2 }, // 4
};

// Surface kernel list: vpar-direction
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_drag_surf_kern_list ser_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, rad_gyrokinetic_drag_surfvpar_1x1v_ser_p1, rad_gyrokinetic_drag_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, rad_gyrokinetic_drag_surfvpar_1x2v_ser_p1, rad_gyrokinetic_drag_surfvpar_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, rad_gyrokinetic_drag_surfvpar_2x2v_ser_p1, rad_gyrokinetic_drag_surfvpar_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, rad_gyrokinetic_drag_surfvpar_3x2v_ser_p1, rad_gyrokinetic_drag_surfvpar_3x2v_ser_p2 }, // 3
};

// Surface kernel list: mu-direction
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_drag_surf_kern_list ser_surf_mu_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, rad_gyrokinetic_drag_surfmu_1x2v_ser_p1, rad_gyrokinetic_drag_surfmu_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, rad_gyrokinetic_drag_surfmu_2x2v_ser_p1, rad_gyrokinetic_drag_surfmu_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, rad_gyrokinetic_drag_surfmu_3x2v_ser_p1, rad_gyrokinetic_drag_surfmu_3x2v_ser_p2 }, // 3
};

// Boundary surface kernel (zero-flux BCs) list: vpar-direction
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_drag_boundary_surf_kern_list ser_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, rad_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p1, rad_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, rad_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p1, rad_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, rad_gyrokinetic_drag_boundary_surfvpar_2x2v_ser_p1, rad_gyrokinetic_drag_boundary_surfvpar_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, rad_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p1, rad_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p2 }, // 3
};

// Constant nu boundary surface kernel (zero-flux BCs) list: mu-direction
GKYL_CU_D
static const gkyl_dg_rad_gyrokinetic_drag_boundary_surf_kern_list ser_boundary_surf_mu_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, rad_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p1, rad_gyrokinetic_drag_boundary_surfmu_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, rad_gyrokinetic_drag_boundary_surfmu_2x2v_ser_p1, rad_gyrokinetic_drag_boundary_surfmu_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, rad_gyrokinetic_drag_boundary_surfmu_3x2v_ser_p1, rad_gyrokinetic_drag_boundary_surfmu_3x2v_ser_p2 }, // 3
};

void gkyl_rad_gyrokinetic_drag_free(const struct gkyl_ref_count* ref);

// Function called by eval_on_nodes to calculate vnu
void vnu_calc(double t, const double *xn, double *fout, void *ctx);

// Function called by eval_on_nodes to calculate vsqnu
void vsqnu_calc(double t, const double *xn, double *fout, void *ctx);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idxC);
  if (dir >= rad_gyrokinetic_drag->cdim+1) {
    return rad_gyrokinetic_drag->surf[dir-rad_gyrokinetic_drag->cdim](xcC, dxC,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx),
								      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx),
      qInL, qInC, qInR, qRhsOut);
  } else if (dir >= rad_gyrokinetic_drag->cdim) {
    return rad_gyrokinetic_drag->surf[dir-rad_gyrokinetic_drag->cdim](xcC, dxC,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx),
								      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx),
      qInL, qInC, qInR, qRhsOut);
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
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idxSkin);
  if (dir >= rad_gyrokinetic_drag->cdim+1) {
    return rad_gyrokinetic_drag->boundary_surf[dir-rad_gyrokinetic_drag->cdim](xcSkin, dxSkin,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx),
									       (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx),
      edge, qInSkin, qInEdge, qRhsOut);
  } else if (dir >= rad_gyrokinetic_drag->cdim) {
    return rad_gyrokinetic_drag->boundary_surf[dir-rad_gyrokinetic_drag->cdim](xcSkin, dxSkin,
      (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx),
									       (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, cidx),
      edge, qInSkin, qInEdge, qRhsOut);
  }
  return 0.;
}

