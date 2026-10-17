// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_fluid_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

typedef double (*pkpm_fluid_flux_t)(
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, const double *w,
  const double *dxv, const double *moms_l, const double *moms_r, const double *bvar_l,
  const double *bvar_r, const double *u_l, const double *u_r, const double *euler_l,
  const double *euler_r, double *GKYL_RESTRICT flux, double *GKYL_RESTRICT mass_flux
);

// for use in kernel tables
typedef struct {
  pkpm_fluid_flux_t kernels[3];
} gkyl_dg_pkpm_fluid_flux_kern_list;

struct gkyl_dg_pkpm_fluid_flux_surf {
  struct gkyl_rect_grid conf_grid; // Configuration-space grid.
  int cdim; // Configuration-space dimensions.
  int num_surf_basis; // Surface modal basis size per component and direction.
  bool use_gpu; // bool to determine if on GPU.
  const struct gkyl_wv_eqn *wv_eqn; // Ten-moment wave equation (device pointer on device).
  const struct gkyl_wave_geom *geom; // Wave geometry (device pointer on device).
  pkpm_fluid_flux_t fluid_flux[3]; // Flux kernel in each direction.

  uint32_t flags;
  struct gkyl_dg_pkpm_fluid_flux_surf *on_dev; // pointer to itself or device data
};

// Momentum flux kernels (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_fluid_flux_kern_list ser_pkpm_fluid_flux_kernels[] = {
  {pkpm_fluid_flux_x_1x_ser_p1, NULL, NULL}, // 0
  {pkpm_fluid_flux_x_2x_ser_p1, pkpm_fluid_flux_y_2x_ser_p1, NULL}, // 1
  {pkpm_fluid_flux_x_3x_ser_p1, pkpm_fluid_flux_y_3x_ser_p1, pkpm_fluid_flux_z_3x_ser_p1}, // 2
};

GKYL_CU_D
static pkpm_fluid_flux_t
choose_pkpm_fluid_flux_kern(enum gkyl_basis_type b_type, int cdim, int poly_order, int dir)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_fluid_flux_kernels[cdim - 1].kernels[dir];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}
