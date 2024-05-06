// Private header: not for direct use
#pragma once


#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_fpo_vlasov_kernels.h>

static void
create_offsets(const int num_up_dirs, const bool is_edge_lower[], const bool is_edge_upper[], const int update_dirs[], const struct gkyl_range *range, long offsets[])
{
  // Construct the offsets *only* in the directions being updated.
  // No need to load the neighbors that are not needed for the update.
  int lower_offset[GKYL_MAX_DIM] = {0};
  int upper_offset[GKYL_MAX_DIM] = {0};
  for (int d=0; d<num_up_dirs; ++d) {
    int dir = update_dirs[d];
    lower_offset[dir] = -1 + is_edge_lower[d];
    upper_offset[dir] = 1 - is_edge_upper[d];
  }  

  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, lower_offset, upper_offset);
  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);
  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

static int idx_to_inloup_ker(int dim, const int *idx, const int *dirs, const int *num_cells) {
  int iout = 0;

  for (int d=0; d<dim; ++d) {
    if (idx[dirs[d]] == 1) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[dirs[d]] == num_cells[dirs[d]]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}

// Kernel function pointers
typedef int (*fpo_drag_coeff_t)(const double *dxv, const double *gamma,
  const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff,
  double *drag_coeff_surf, double *sgn_drag_coeff_surf);

// For use in kernel tables
typedef struct { fpo_drag_coeff_t kernels[3]; } gkyl_dg_fpo_drag_coeff_kern_list;
typedef struct { gkyl_dg_fpo_drag_coeff_kern_list list[3]; } gkyl_dg_fpo_drag_coeff_stencil_list;

// drag coefficient kernel lists
GKYL_CU_D
static const gkyl_dg_fpo_drag_coeff_stencil_list ser_fpo_drag_coeff_1x3v_vx_kernels = {
  {
    {NULL, NULL, NULL},
    {fpo_drag_coeff_1x3v_vx_ser_p1_invx, fpo_drag_coeff_1x3v_vx_ser_p1_lovx, fpo_drag_coeff_1x3v_vx_ser_p1_upvx},
    {fpo_drag_coeff_1x3v_vx_ser_p2_invx, fpo_drag_coeff_1x3v_vx_ser_p2_lovx, fpo_drag_coeff_1x3v_vx_ser_p2_upvx}
  }
};

GKYL_CU_D
static const gkyl_dg_fpo_drag_coeff_stencil_list ser_fpo_drag_coeff_1x3v_vy_kernels = {
  {
    {NULL, NULL, NULL},
    {fpo_drag_coeff_1x3v_vy_ser_p1_invy, fpo_drag_coeff_1x3v_vy_ser_p1_lovy, fpo_drag_coeff_1x3v_vy_ser_p1_upvy},
    {fpo_drag_coeff_1x3v_vy_ser_p2_invy, fpo_drag_coeff_1x3v_vy_ser_p2_lovy, fpo_drag_coeff_1x3v_vy_ser_p2_upvy}
  }
};

GKYL_CU_D
static const gkyl_dg_fpo_drag_coeff_stencil_list ser_fpo_drag_coeff_1x3v_vz_kernels = {
  {
    {NULL, NULL, NULL},
    {fpo_drag_coeff_1x3v_vz_ser_p1_invz, fpo_drag_coeff_1x3v_vz_ser_p1_lovz, fpo_drag_coeff_1x3v_vz_ser_p1_upvz},
    {fpo_drag_coeff_1x3v_vz_ser_p2_invz, fpo_drag_coeff_1x3v_vz_ser_p2_lovz, fpo_drag_coeff_1x3v_vz_ser_p2_upvz}
  }
};

GKYL_CU_D
static fpo_drag_coeff_t
choose_ser_fpo_drag_coeff_recovery_kern(int dir, int cdim, int poly_order, int stencil_idx)
{
  if (dir == 0)
    return ser_fpo_drag_coeff_1x3v_vx_kernels.list[poly_order].kernels[stencil_idx];
  else if (dir == 1)
    return ser_fpo_drag_coeff_1x3v_vy_kernels.list[poly_order].kernels[stencil_idx];
  else if (dir == 2)
    return ser_fpo_drag_coeff_1x3v_vz_kernels.list[poly_order].kernels[stencil_idx];
  else
    return NULL;
};

