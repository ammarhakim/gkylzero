#pragma once

#include <dg_euleriso_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// Types for various kernels
typedef double (*diffusion_euler_iso_surf_t)(const double *w, const double *dx,
  const double* D, const double *uvarl, const double *uvarc, const double *uvarr,
  const double *statevecl, const double *statevecc, const double *statevecr,
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_diffusion_euler_iso_vol_kern_list;
typedef struct { diffusion_euler_iso_surf_t kernels[3]; } gkyl_dg_diffusion_euler_iso_surf_kern_list;

struct dg_diffusion_euler_iso {
  struct gkyl_dg_eqn eqn;
  diffusion_euler_iso_surf_t surf[3];
  struct gkyl_range conf_range;
  struct gkyl_dg_diffusion_euler_iso_auxfields auxfields;
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_dg_diffusion_euler_iso_vol_1x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_euler_iso* diffusion = container_of(eqn, struct dg_diffusion_euler_iso, eqn);

  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);

  return dg_euleriso_diffusion_vol_1x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    (const double*) gkyl_array_cfetch(diffusion->auxfields.u_i, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_euler_iso_vol_1x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_euler_iso* diffusion = container_of(eqn, struct dg_diffusion_euler_iso, eqn);

  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);

  return dg_euleriso_diffusion_vol_1x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    (const double*) gkyl_array_cfetch(diffusion->auxfields.u_i, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_euler_iso_vol_2x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_euler_iso* diffusion = container_of(eqn, struct dg_diffusion_euler_iso, eqn);

  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);

  return dg_euleriso_diffusion_vol_2x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    (const double*) gkyl_array_cfetch(diffusion->auxfields.u_i, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_euler_iso_vol_2x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_euler_iso* diffusion = container_of(eqn, struct dg_diffusion_euler_iso, eqn);

  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);

  return dg_euleriso_diffusion_vol_2x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    (const double*) gkyl_array_cfetch(diffusion->auxfields.u_i, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_euler_iso_vol_3x_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_euler_iso* diffusion = container_of(eqn, struct dg_diffusion_euler_iso, eqn);

  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);

  return dg_euleriso_diffusion_vol_3x_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    (const double*) gkyl_array_cfetch(diffusion->auxfields.u_i, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_dg_diffusion_euler_iso_vol_3x_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_euler_iso* diffusion = container_of(eqn, struct dg_diffusion_euler_iso, eqn);

  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);

  return dg_euleriso_diffusion_vol_3x_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx),
    (const double*) gkyl_array_cfetch(diffusion->auxfields.u_i, cidx),
    qIn, qRhsOut);
}

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_diffusion_euler_iso_vol_kern_list ser_vol_kernels[] = {
  { NULL, kernel_dg_diffusion_euler_iso_vol_1x_ser_p1, kernel_dg_diffusion_euler_iso_vol_1x_ser_p2 },
  { NULL, kernel_dg_diffusion_euler_iso_vol_2x_ser_p1, kernel_dg_diffusion_euler_iso_vol_2x_ser_p2 },
  { NULL, kernel_dg_diffusion_euler_iso_vol_3x_ser_p1, kernel_dg_diffusion_euler_iso_vol_3x_ser_p2 },
};

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_diffusion_euler_iso_surf_kern_list ser_surf_x_kernels[] = {
  { NULL, dg_euleriso_diffusion_surfx_1x_ser_p1, dg_euleriso_diffusion_surfx_1x_ser_p2 },
  { NULL, dg_euleriso_diffusion_surfx_2x_ser_p1, dg_euleriso_diffusion_surfx_2x_ser_p2 },
  { NULL, dg_euleriso_diffusion_surfx_3x_ser_p1, dg_euleriso_diffusion_surfx_3x_ser_p2 },
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_diffusion_euler_iso_surf_kern_list ser_surf_y_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, dg_euleriso_diffusion_surfy_2x_ser_p1, dg_euleriso_diffusion_surfy_2x_ser_p2 },
  { NULL, dg_euleriso_diffusion_surfy_3x_ser_p1, dg_euleriso_diffusion_surfy_3x_ser_p2 },
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_diffusion_euler_iso_surf_kern_list ser_surf_z_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, dg_euleriso_diffusion_surfz_3x_ser_p1, dg_euleriso_diffusion_surfz_3x_ser_p2 },
};

/**
 * Free diffusion equation object
 *
 * @param ref Reference counter for constant diffusion equation
 */
void gkyl_diffusion_euler_iso_free(const struct gkyl_ref_count* ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR,
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_euler_iso* diffusion = container_of(eqn, struct dg_diffusion_euler_iso, eqn);

  long cidx_l = gkyl_range_idx(&diffusion->conf_range, idxL);
  long cidx_c = gkyl_range_idx(&diffusion->conf_range, idxC);
  long cidx_r = gkyl_range_idx(&diffusion->conf_range, idxR);

  return diffusion->surf[dir](xcC, dxC,
    (const double*) gkyl_array_cfetch(diffusion->auxfields.D, cidx_c),
    (const double*) gkyl_array_cfetch(diffusion->auxfields.u_i, cidx_l),
    (const double*) gkyl_array_cfetch(diffusion->auxfields.u_i, cidx_c),
    (const double*) gkyl_array_cfetch(diffusion->auxfields.u_i, cidx_r),
    qInL, qInC, qInR, qRhsOut);
}
