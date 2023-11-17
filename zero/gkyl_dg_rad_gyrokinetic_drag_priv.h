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
  struct gkyl_range prange; // phase space range.
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
  printf("In kernel_rad_gyro (drag_priv.h)\n");
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  long cidx = gkyl_range_idx(&rad_gyrokinetic_drag->conf_range, idx);
  printf("Before call to kernel, cidx=%d\n", cidx);
  printf("cdim=%d, pdim=%d\n",rad_gyrokinetic_drag->cdim,rad_gyrokinetic_drag->pdim);
  printf("Before call to kernel, idx[0]=%d, idx[1]=%d\n", idx[0],idx[1]);
  const double *a = (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.nI, cidx);
  printf("After a, a[0]=%f, a[1]=%f\n",a[0],a[1]);
  printf("conf int %d\n",rad_gyrokinetic_drag->conf_range.ndim);
  printf("int %d\n",rad_gyrokinetic_drag->prange.ndim);
  long pidx = gkyl_range_idx(&rad_gyrokinetic_drag->prange, idx);
  printf("pidx=%d\n",pidx);
  const double *ap = (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.bmag, cidx);
  printf("After ap, ap[0]=%f, ap[1]=%f\n",ap[0],ap[1]);
  const double *b = (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vnu, pidx);
  printf("After bb\n");
  const double *c = (const double*) gkyl_array_cfetch(rad_gyrokinetic_drag->auxfields.vsqnu, cidx);
  printf("After c\n");

  for (int i=0;i<2;i++){
    printf("w=%f,dxv=%f,a=%f,b=%f,c=%f,f=%f\n",xc[i],dx[i],a[i],b[i],c[i],qIn[i]);
  }
  
  // For testing
  const double w[2]={1.,2.};
  const double dxv[2] ={1.,2.};
  const double nI[8]={1.,2.,3.,4.,5.,6.,7.,8.};
  const double vnu[8]={1.,2.,3.,4.,5.,6.,7.,8.};
  const double vsqnu[8]={1.2,2.2,3.2,4.2,5.2,6.2,7.2,8.2};
  const double f[8]={1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5};
  double *rhs;
  //printf("before test call\n"); this works
  //rad_gyrokinetic_drag_vol_1x1v_ser_p2(w,dxv,nI,vnu,vsqnu,f,rhs);
  

  
  return rad_gyrokinetic_drag_vol_1x1v_ser_p2(xc, dx,
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
  return rad_gyrokinetic_drag_vol_1x2v_ser_p1(xc, dx,
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
  return rad_gyrokinetic_drag_vol_1x2v_ser_p2(xc, dx,
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
  return rad_gyrokinetic_drag_vol_2x2v_ser_p1(xc, dx,
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
  return rad_gyrokinetic_drag_vol_2x2v_ser_p2(xc, dx,
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
  return rad_gyrokinetic_drag_vol_3x2v_ser_p1(xc, dx,
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
  return rad_gyrokinetic_drag_vol_3x2v_ser_p2(xc, dx,
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

