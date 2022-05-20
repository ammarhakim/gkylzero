#pragma once

#include <gkyl_const_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in advection DG equation object creation
// functions

// Types for various kernels
typedef double (*const_diffusion_vol_t)(const double* w, const double* dx,
  const double* D, const double* f, double* GKYL_RESTRICT out);

typedef void (*const_diffusion_surf_t)(const double *w, const double *dx,
  const double* D, const double *fl, const double *fc, const double *fr,
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { const_diffusion_vol_t kernels[3]; } gkyl_dg_const_diffusion_vol_kern_list;
typedef struct { const_diffusion_surf_t kernels[3]; } gkyl_dg_const_diffusion_surf_kern_list;

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_const_diffusion_vol_kern_list vol_kernels[] = {
  { NULL, const_diffusion_vol_1x_ser_p1, const_diffusion_vol_1x_ser_p2 }, // 1
  { NULL, const_diffusion_vol_2x_ser_p1, NULL }, // 2
  { NULL, NULL, NULL },
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_const_diffusion_surf_kern_list surf_x_kernels[] = {
  { NULL, const_diffusion_surfx_1x_ser_p1, const_diffusion_surfx_1x_ser_p2 },
  { NULL, const_diffusion_surfx_2x_ser_p1, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_const_diffusion_surf_kern_list surf_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, const_diffusion_surfy_2x_ser_p1, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_const_diffusion_surf_kern_list surf_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, NULL, NULL },
};


struct dg_const_diffusion {
  struct gkyl_dg_eqn eqn; // Base object
  const_diffusion_vol_t vol; // Volume kernel
  const_diffusion_surf_t surf[3]; // Surface terms
  double D[GKYL_MAX_DIM]; // pointer to diffusion tensor
};

/**
 * Free constant diffusion equation object
 *
 * @param ref Reference counter for constant diffusion equation
 */
void gkyl_const_diffusion_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn* eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_const_diffusion* const_diffusion = container_of(eqn, struct dg_const_diffusion, eqn);

  return const_diffusion->vol(xc, dx, const_diffusion->D, qIn, qRhsOut);
}

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_const_diffusion* const_diffusion = container_of(eqn, struct dg_const_diffusion, eqn);

  const_diffusion->surf[dir](xcC, dxC, const_diffusion->D,
    qInL, qInC, qInR, qRhsOut);
}

GKYL_CU_D
static void
boundary_surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  
}
