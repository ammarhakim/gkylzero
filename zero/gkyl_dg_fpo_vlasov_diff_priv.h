#pragma once

#include <gkyl_fpo_vlasov_kernels.h>

// private header for use in fpo_vlasov_diff DG equation object creation
// functions

// Types for various kernels
typedef double (*fpo_vlasov_diff_vol_t)(const double* w, const double* dx,
  const double* g, const double* f, double* GKYL_RESTRICT out);

typedef void (*fpo_vlasov_diff_surf_t)(const double *w, const double *dx,
  const double* g, const double *f[],
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { fpo_vlasov_diff_vol_t kernels[3]; } gkyl_dg_fpo_vlasov_diff_vol_kern_list;
typedef struct { fpo_vlasov_diff_surf_t kernels[3]; } gkyl_dg_fpo_vlasov_diff_surf_kern_list;
typedef struct { fpo_vlasov_diff_surf_t kernels[3]; } gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list;

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_vol_kern_list ser_vol_kernels[] = {
  // { NULL, fpo_vlasov_diff_vol_1x3v_ser_p1, fpo_vlasov_diff_vol_1x3v_ser_p2 }, // 0
  // { NULL, fpo_vlasov_diff_vol_2x3v_ser_p1, fpo_vlasov_diff_vol_2x3v_ser_p2 }, // 1
  // { NULL, fpo_vlasov_diff_vol_3x3v_ser_p1, NULL }, // 2
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: xx-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_surf_kern_list ser_surf_xx_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: xy-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_surf_kern_list ser_surf_xy_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: xz-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_surf_kern_list ser_surf_xz_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: yx-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_surf_kern_list ser_surf_yx_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: yy-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_surf_kern_list ser_surf_yy_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: yz-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_surf_kern_list ser_surf_yz_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: zx-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_surf_kern_list ser_surf_zx_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: zy-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_surf_kern_list ser_surf_zy_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Surface kernel list: zz-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_surf_kern_list ser_surf_zz_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: xx-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_xx_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: xy-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_xy_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: xz-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_xz_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: yx-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_yx_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: yy-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_yy_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: yz-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_yz_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: zx-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_zx_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: zy-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_zy_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

// Boundary Surface kernel list: zz-direction
GKYL_CU_D
static const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list ser_boundary_surf_zz_kernels[] = {
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
  { NULL, NULL, NULL },
};

struct dg_fpo_vlasov_diff {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  fpo_vlasov_diff_vol_t vol; // Volume kernel
  fpo_vlasov_diff_surf_t surf[3][3]; // Surface terms for acceleration
  fpo_vlasov_diff_surf_t boundary_surf[3][3]; // Surface terms for acceleration
  struct gkyl_range phase_range; // Configuration space range.
  struct gkyl_dg_fpo_vlasov_diff_auxfields auxfields; // Auxiliary fields.
};

/**
 * Free fpo_vlasov_diff equation object
 *
 * @param ref Reference counter for constant fpo_vlasov_diff equation
 */
void gkyl_fpo_vlasov_diff_free(const struct gkyl_ref_count* ref);


GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn* eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);
  
  long pidx = gkyl_range_idx(&fpo_vlasov_diff->phase_range, idx);
  
  return fpo_vlasov_diff->vol(xc, dx,
    (const double*) gkyl_array_cfetch(fpo_vlasov_diff->auxfields.g, pidx),
    qIn, qRhsOut);
}

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn* eqn, int dir1, int dir2,
  const double* xc, const double* dxc, const int* idxc,
  const int idx[][GKYL_MAX_DIM], const double* qIn[],
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);
  int cdim = fpo_vlasov_diff->cdim;
  long pidx = gkyl_range_idx(&fpo_vlasov_diff->phase_range, idxc);

  if (dir1 >= cdim && dir2 >= cdim) {
    fpo_vlasov_diff->surf[dir1-cdim][dir2-cdim](xc, dxc,
      (const double*) gkyl_array_cfetch(fpo_vlasov_diff->auxfields.g, pidx), 
      qIn, qRhsOut);
  }
}

GKYL_CU_D
static void
boundary_surf(const struct gkyl_dg_eqn* eqn, int dir1, int dir2,
  const double* xc, const double* dxc, const int* idxc,
  const int idx[][GKYL_MAX_DIM], const double* qIn[],
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);
  int cdim = fpo_vlasov_diff->cdim;
  long pidx = gkyl_range_idx(&fpo_vlasov_diff->phase_range, idxc);

  if (dir1 >= cdim && dir2 >= cdim) {
    fpo_vlasov_diff->surf[dir1-cdim][dir2-cdim](xc, dxc,
      (const double*) gkyl_array_cfetch(fpo_vlasov_diff->auxfields.g, pidx), 
      qIn, qRhsOut);
  }
}
