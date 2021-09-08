#pragma once

// Private header, not for direct use in user code

#include <gkyl_const_diffusion_kernels.h>

// Types for various kernels
typedef double (*const_diffusion_vol_t)(const double* w, const double* dx,
  const double* D, const double* f, double* GKYL_RESTRICT out);

typedef void (*const_diffusion_surf_t)(const double *w, const double *dx,
  const double* D, const double *fl, const double *fc, const double *fr,
  double* GKYL_RESTRICT out);

// Volume kernel list
GKYL_CU_D
static struct { const_diffusion_vol_t kernels[3]; } vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, const_diffusion_vol_1x_ser_p1, const_diffusion_vol_1x_ser_p2 }, // 1
  { NULL, const_diffusion_vol_2x_ser_p1, NULL }, // 2
};


// Surface kernel list: x-direction
GKYL_CU_D
static struct { const_diffusion_surf_t kernels[3]; } surf_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, const_diffusion_surfx_1x_ser_p1, const_diffusion_surfx_1x_ser_p2 }, // 1
  { NULL, const_diffusion_surfx_2x_ser_p1, NULL }, // 2
};

// Surface kernel list: x-direction
GKYL_CU_D
static struct { const_diffusion_surf_t kernels[3]; } surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, const_diffusion_surfy_2x_ser_p1, NULL }, // 2
};


// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,dim,poly_order) lst[dim].kernels[poly_order]

struct dg_const_diffusion {
  struct gkyl_dg_eqn eqn; // Base object
  int dim; // Dimensions
  const_diffusion_vol_t vol; // Volume kernel
  const_diffusion_surf_t surf[3]; // Surface terms
  //diffusion_const_boundary_surf_t boundary_surf[3]; // Surface terms for acceleration
  const double* D; // pointer to diffusion tensor
};

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_const_diffusion* const_diffusion = container_of(eqn, struct dg_const_diffusion, eqn);

  return const_diffusion->vol(xc, dx, const_diffusion->D, qIn, qRhsOut);
}

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  double maxsOld, const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_const_diffusion* const_diffusion = container_of(eqn, struct dg_const_diffusion, eqn);

  const_diffusion->surf[dir](xcC, dxC, const_diffusion->D, qInL, qInC, qInR, qRhsOut);

  return 0.0;
}

// GKYL_CU_D
// static double
// boundary_surf(const struct gkyl_dg_eqn *eqn,
//   int dir,
//   const double* xcEdge, const double* xcSkin,
//   const double* dxEdge, const double* dxSkin,
//   double maxsOld, const int* idxEdge, const int* idxSkin, const int edge,
//   const double* qInEdge, const double* qInSkin,
//   double* GKYL_RESTRICT qRhsOut)
// {
//   struct dg_diffusion *diffusion = container_of(eqn, struct dg_vlasov, eqn);

//   assert(false);
//   return 0.0;
// }

