#pragma once

#include <gkyl_advection_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in advection DG equation object creation
// functions

// Types for various kernels
typedef double (*advection_vol_t)(const double *w, const double *dx, 
  const double *u, const double *q, double* GKYL_RESTRICT out);

typedef void (*advection_surf_t)(const double *w, const double *dx,
  const double *ul, const double *uc, const double *ur,
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { advection_vol_t kernels[3]; } gkyl_dg_advection_vol_kern_list;
typedef struct { advection_surf_t kernels[3]; } gkyl_dg_advection_surf_kern_list;

//
// Serendipity basis kernels
// 

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_advection_vol_kern_list ser_vol_kernels[] = {
  { NULL, advection_vol_1x_ser_p1, advection_vol_1x_ser_p2 }, // 0
  { NULL, advection_vol_2x_ser_p1, advection_vol_2x_ser_p2 }, // 1
  { NULL, advection_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_advection_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, advection_surfx_1x_ser_p1, advection_surfx_1x_ser_p2 }, // 0
  { NULL, advection_surfx_2x_ser_p1, advection_surfx_2x_ser_p2 }, // 1
  { NULL, advection_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_advection_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, advection_surfy_2x_ser_p1, advection_surfy_2x_ser_p2 }, // 1
  { NULL, advection_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_advection_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, advection_surfz_3x_ser_p1, NULL }, // 2
};

//
// Tensor-product basis kernels
// 

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_advection_vol_kern_list ten_vol_kernels[] = {
  { NULL, advection_vol_1x_ser_p1, advection_vol_1x_tensor_p2 }, // 0
  { NULL, advection_vol_2x_ser_p1, advection_vol_2x_tensor_p2 }, // 1
  { NULL, advection_vol_3x_ser_p1, NULL },              // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_advection_surf_kern_list ten_surf_x_kernels[] = {
  { NULL, advection_surfx_1x_ser_p1, advection_surfx_1x_tensor_p2 }, // 0
  { NULL, advection_surfx_2x_ser_p1, advection_surfx_2x_tensor_p2 }, // 1
  { NULL, advection_surfx_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_advection_surf_kern_list ten_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, advection_surfy_2x_ser_p1, advection_surfy_2x_tensor_p2 }, // 1
  { NULL, advection_surfy_3x_ser_p1, NULL },                 // 2
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_advection_surf_kern_list ten_surf_z_kernels[] = {
  { NULL, NULL, NULL },                 // 0
  { NULL, NULL, NULL },                 // 1
  { NULL, advection_surfz_3x_ser_p1, NULL }, // 2
};

struct dg_advection {
  struct gkyl_dg_eqn eqn; // Base object    
  advection_vol_t vol; // pointer to volume kernel
  advection_surf_t surf[3]; // pointers to surface kernels
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_advection_auxfields auxfields; // Auxiliary fields.
};


/**
 * Free advection equation object
 *
 * @param ref Reference counter for advection equation
 */
void gkyl_advection_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double*  dx, 
  const int*  idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  long cidx = gkyl_range_idx(&advection->conf_range, idx);
  return advection->vol(xc, dx, 
    (const double*) gkyl_array_cfetch(advection->auxfields.u, cidx), qIn, qRhsOut);
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
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);

  long cidx_l = gkyl_range_idx(&advection->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&advection->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&advection->conf_range, idxR);

  advection->surf[dir](xcC, dxC, 
    (const double*) gkyl_array_cfetch(advection->auxfields.u, cidx_l),
    (const double*) gkyl_array_cfetch(advection->auxfields.u, cidx_c),
    (const double*) gkyl_array_cfetch(advection->auxfields.u, cidx_r), 
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
