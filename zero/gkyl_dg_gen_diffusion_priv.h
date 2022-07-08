#pragma once

#include <gkyl_dg_gen_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in gen_diffusion DG equation object creation
// functions

// Types for various kernels
typedef double (*gen_diffusion_vol_t)(const double* w, const double* dx,
  const double* D, const double* q, double* GKYL_RESTRICT out);

typedef void (*gen_diffusion_surf_t)(const double *w, const double *dx,
  const double* D, const double *q[],
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { gen_diffusion_vol_t kernels[3]; } gkyl_dg_gen_diffusion_vol_kern_list;
typedef struct { gen_diffusion_surf_t kernels[3]; } gkyl_dg_gen_diffusion_surf_kern_list;

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_gen_diffusion_vol_kern_list ser_vol_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, dg_gen_diffusion_vol_2x_ser_p1, dg_gen_diffusion_vol_2x_ser_p2 },
  { NULL, dg_gen_diffusion_vol_3x_ser_p1, dg_gen_diffusion_vol_3x_ser_p2 },
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, dg_gen_diffusion_surfx_2x_ser_p1, dg_gen_diffusion_surfx_2x_ser_p2 },
  { NULL, dg_gen_diffusion_surfx_3x_ser_p1, dg_gen_diffusion_surfx_3x_ser_p2 },
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL }, // no y-direction in 1D
  { NULL, dg_gen_diffusion_surfy_2x_ser_p1, dg_gen_diffusion_surfy_2x_ser_p2 },
  { NULL, dg_gen_diffusion_surfy_3x_ser_p1, dg_gen_diffusion_surfy_3x_ser_p2 },
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_gen_diffusion_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL }, // no z-direction in 1D
  { NULL, NULL, NULL }, // no z-direction in 2D
  { NULL, dg_gen_diffusion_surfz_3x_ser_p1, dg_gen_diffusion_surfz_3x_ser_p2 },
};

struct dg_gen_diffusion {
  struct gkyl_dg_eqn eqn;
  gen_diffusion_vol_t vol;
  gen_diffusion_surf_t surf[3];
  struct gkyl_range conf_range;
  struct gkyl_dg_gen_diffusion_auxfields auxfields;
};

/**
 * Free gen_diffusion equation object
 *
 * @param ref Reference counter for constant gen_diffusion equation
 */
void gkyl_gen_diffusion_free(const struct gkyl_ref_count* ref);


GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn* eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gen_diffusion* gen_diffusion = container_of(eqn, struct dg_gen_diffusion, eqn);
  
  long cidx = gkyl_range_idx(&gen_diffusion->conf_range, idx);
  
  return gen_diffusion->vol(xc, dx,
    (const double*) gkyl_array_cfetch(gen_diffusion->auxfields.Dij, cidx),
    qIn, qRhsOut);
}

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xc, const double* dx, const int* idx,
  const double* qIn[],
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gen_diffusion* gen_diffusion = container_of(eqn, struct dg_gen_diffusion, eqn);
  long cidx = gkyl_range_idx(&gen_diffusion->conf_range, idx);
  
  gen_diffusion->surf[dir](xc, dx,
    (const double*) gkyl_array_cfetch(gen_diffusion->auxfields.Dij, cidx), 
    qIn, qRhsOut);
}
