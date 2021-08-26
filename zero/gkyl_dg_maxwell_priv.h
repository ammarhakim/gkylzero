#pragma once

#include <gkyl_maxwell_kernels.h>

// private header for use in Maxwell DG equation object creation
// functions

// Types for various kernels
typedef double (*maxwell_vol_t)(const gkyl_maxwell_inp *meq,
  const double *w, const double *dx, const double *q, double* GKYL_RESTRICT out);

typedef void (*maxwell_surf_t)(const gkyl_maxwell_inp *meq, const double *w, const double *dx,
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { maxwell_vol_t kernels[3]; } gkyl_dg_maxwell_vol_kern_list;
typedef struct { maxwell_surf_t kernels[3]; } gkyl_dg_maxwell_surf_kern_list;

//
// Serendipity basis kernels
// 

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_maxwell_vol_kern_list ser_vol_kernels[] = {
  { NULL, maxwell_vol_1x_ser_p1, maxwell_vol_1x_ser_p2 }, // 0
  { NULL, maxwell_vol_2x_ser_p1, maxwell_vol_2x_ser_p2 }, // 1
  { NULL, maxwell_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_ser_p2 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_ser_p2 }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_ser_p2 }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, maxwell_surfz_3x_ser_p1, NULL }, // 2
};

//
// Tensor-product basis kernels
// 

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_maxwell_vol_kern_list ten_vol_kernels[] = {
  { NULL, maxwell_vol_1x_ser_p1, maxwell_vol_1x_tensor_p2 }, // 0
  { NULL, maxwell_vol_2x_ser_p1, maxwell_vol_2x_tensor_p2 }, // 1
  { NULL, maxwell_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_x_kernels[] = {
  { NULL, maxwell_surfx_1x_ser_p1, maxwell_surfx_1x_tensor_p2 }, // 0
  { NULL, maxwell_surfx_2x_ser_p1, maxwell_surfx_2x_tensor_p2 }, // 1
  { NULL, maxwell_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, maxwell_surfy_2x_ser_p1, maxwell_surfy_2x_tensor_p2 }, // 1
  { NULL, maxwell_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_maxwell_surf_kern_list ten_surf_z_kernels[] = {
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
static void
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  maxwell->surf[dir](&maxwell->maxwell_data, xcC, dxC,
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
