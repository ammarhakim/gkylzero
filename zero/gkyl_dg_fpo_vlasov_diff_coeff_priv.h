#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_fpo_vlasov_kernels.h>

// Kernel function pointers
typedef void (*fpo_diff_coeff_t)(const double *dxv, const double *G_l, const double *G_c, const double *G_r, double *diff_coeff);

typedef void (*fpo_diff_coeff_surf_t)(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff);

// For use in kernel tables
typedef struct { fpo_diff_coeff_t kernels[3]; } gkyl_dg_fpo_diff_coeff_kern_list;
typedef struct { fpo_diff_coeff_surf_t kernels[3]; } gkyl_dg_fpo_diff_coeff_surf_kern_list;

// diffusion coefficient diagonal vx term volume kernels
GKYL_CU_D
static const gkyl_dg_fpo_diff_coeff_kern_list ser_fpo_diff_coeff_vx_kernels[] = {
  NULL, fpo_diff_coeff_recov_vx_1x3v_ser_p1, fpo_diff_coeff_recov_vx_1x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_vx_2x3v_ser_p1, fpo_diff_coeff_recov_vx_2x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_vx_3x3v_ser_p1, NULL,
};

// diffusion coefficient diagonal vy term volume kernels
GKYL_CU_D
static const gkyl_dg_fpo_diff_coeff_kern_list ser_fpo_diff_coeff_vy_kernels[] = {
  NULL, fpo_diff_coeff_recov_vy_1x3v_ser_p1, fpo_diff_coeff_recov_vy_1x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_vy_2x3v_ser_p1, fpo_diff_coeff_recov_vy_2x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_vy_3x3v_ser_p1, NULL,
};

// diffusion coefficient diagonal vz term volume kernels
GKYL_CU_D
static const gkyl_dg_fpo_diff_coeff_kern_list ser_fpo_diff_coeff_vz_kernels[] = {
  NULL, fpo_diff_coeff_recov_vz_1x3v_ser_p1, fpo_diff_coeff_recov_vz_1x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_vz_2x3v_ser_p1, fpo_diff_coeff_recov_vz_2x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_vz_3x3v_ser_p1, NULL,
};

// diffusion coefficient diagonal vx term surface kernels
GKYL_CU_D
static const gkyl_dg_fpo_diff_coeff_surf_kern_list ser_fpo_diff_coeff_vx_surf_kernels[] = {
  NULL, fpo_diff_coeff_recov_surf_vx_1x3v_ser_p1, fpo_diff_coeff_recov_surf_vx_1x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_surf_vx_2x3v_ser_p1, fpo_diff_coeff_recov_surf_vx_2x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_surf_vx_3x3v_ser_p1, NULL,
};

// diffusion coefficient diagonal vy term surface kernels
GKYL_CU_D
static const gkyl_dg_fpo_diff_coeff_surf_kern_list ser_fpo_diff_coeff_vy_surf_kernels[] = {
  NULL, fpo_diff_coeff_recov_surf_vy_1x3v_ser_p1, fpo_diff_coeff_recov_surf_vy_1x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_surf_vy_2x3v_ser_p1, fpo_diff_coeff_recov_surf_vy_2x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_surf_vy_3x3v_ser_p1, NULL,
};

// diffusion coefficient diagonal vz term surface kernels
GKYL_CU_D
static const gkyl_dg_fpo_diff_coeff_surf_kern_list ser_fpo_diff_coeff_vz_surf_kernels[] = {
  NULL, fpo_diff_coeff_recov_surf_vz_1x3v_ser_p1, fpo_diff_coeff_recov_surf_vz_1x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_surf_vz_2x3v_ser_p1, fpo_diff_coeff_recov_surf_vz_2x3v_ser_p2,
  NULL, fpo_diff_coeff_recov_surf_vz_3x3v_ser_p1, NULL,
};

GKYL_CU_D
static const fpo_diff_coeff_t
choose_ser_fpo_diff_coeff_recovery_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ser_fpo_diff_coeff_vx_kernels[cdim-1].kernels[poly_order];
  else if (dir == 1)
    return ser_fpo_diff_coeff_vy_kernels[cdim-1].kernels[poly_order];
  else if (dir == 2)
    return ser_fpo_diff_coeff_vz_kernels[cdim-1].kernels[poly_order];
  else 
    return NULL;
};

GKYL_CU_D
static const fpo_diff_coeff_surf_t
choose_ser_fpo_diff_coeff_recovery_surf_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ser_fpo_diff_coeff_vx_surf_kernels[cdim-1].kernels[poly_order];
  else if (dir == 1)
    return ser_fpo_diff_coeff_vy_surf_kernels[cdim-1].kernels[poly_order];
  else if (dir == 2)
    return ser_fpo_diff_coeff_vz_surf_kernels[cdim-1].kernels[poly_order];
  else 
    return NULL;
};
