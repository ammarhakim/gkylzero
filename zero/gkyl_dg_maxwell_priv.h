#pragma once

#include <gkyl_maxwell_kernels.h>

// private header for use in Maxwell DG equation object creation
// functions

// Types for various kernels
typedef double (*maxwell_vol_t)(const gkyl_maxwell_inp *meq,
  const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out);

typedef double (*maxwell_surf_t)(const gkyl_maxwell_inp *meq,
  const double *w, const double *dx, const double tau,
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);

// Volume kernel list
static struct { maxwell_vol_t kernels[3]; } vol_kernels[] = {
  { NULL, maxwell_vol_1x_ser_p1, maxwell_vol_1x_ser_p2 }, // 0
  { NULL, maxwell_vol_2x_ser_p1, maxwell_vol_2x_ser_p2 }, // 1
  { NULL, maxwell_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
static struct { maxwell_surf_t kernels[3]; } surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_ser_p2 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_ser_p2 }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
static struct { maxwell_surf_t kernels[3]; } surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_ser_p2 }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
static struct { maxwell_surf_t kernels[3]; } surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, maxwell_surfz_3x_ser_p1, NULL }, // 2
};


struct dg_maxwell {
  struct gkyl_dg_eqn eqn; // Base object    
  gkyl_maxwell_inp maxwell_data; // Parameters needed by kernels
  maxwell_vol_t vol; // pointer to volume kernel
  maxwell_surf_t surf[3]; // pointers to surface kernels
};

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double*  dx, 
  const int*  idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  return maxwell->vol(&maxwell->maxwell_data, xc, dx, qIn, qRhsOut);
}

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  double maxsOld, const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  double maxs = maxwell->maxwell_data.c;
  return maxwell->surf[dir](&maxwell->maxwell_data, xcC, dxC,
    maxs, qInL, qInC, qInR, qRhsOut);
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  double maxsOld, const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  return 0;
}
