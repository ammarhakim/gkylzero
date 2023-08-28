#pragma once

#include <gkyl_dg_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in diffusion DG equation object creation
// functions

// Types for various kernels
typedef double (*diffusion_surf_t)(const double *w, const double *dx,
  const double* D, const double *ql, const double *qc, const double *qr,
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_diffusion_vol_kern_list;
typedef struct { diffusion_surf_t kernels[3]; } gkyl_dg_diffusion_surf_kern_list;

struct dg_diffusion {
  struct gkyl_dg_eqn eqn;
  diffusion_surf_t surf[3];
  struct gkyl_range conf_range;
  struct gkyl_dg_diffusion_auxfields auxfields;
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  
  return dg_diffusion_vol_1x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  
  return dg_diffusion_vol_1x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  
  return dg_diffusion_vol_2x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  
  return dg_diffusion_vol_2x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  
  return dg_diffusion_vol_3x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  
  return dg_diffusion_vol_3x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_diffusion_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_dg_diffusion_vol_1x_ser_p1, kernel_dg_diffusion_vol_1x_ser_p2 },
  { NULL, kernel_dg_diffusion_vol_2x_ser_p1, kernel_dg_diffusion_vol_2x_ser_p2 },
  { NULL, kernel_dg_diffusion_vol_3x_ser_p1, kernel_dg_diffusion_vol_3x_ser_p2 },
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, dg_diffusion_surfx_1x_ser_p1, dg_diffusion_surfx_1x_ser_p2 },
  { NULL, dg_diffusion_surfx_2x_ser_p1, dg_diffusion_surfx_2x_ser_p2 },
  { NULL, dg_diffusion_surfx_3x_ser_p1, dg_diffusion_surfx_3x_ser_p2 },
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_diffusion_surfy_2x_ser_p1, dg_diffusion_surfy_2x_ser_p2 },
  { NULL, dg_diffusion_surfy_3x_ser_p1, dg_diffusion_surfy_3x_ser_p2 },
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_diffusion_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_diffusion_surfz_3x_ser_p1, dg_diffusion_surfz_3x_ser_p2 },
};

/**
 * Free diffusion equation object
 *
 * @param ref Reference counter for constant diffusion equation
 */
void gkyl_diffusion_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&diffusion->conf_range, idxC);
  
  return diffusion->surf[dir](xcC, dxC,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx), 
    qInL, qInC, qInR, qRhsOut);
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{ 
  return 0.;
}
