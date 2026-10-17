// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

typedef double (*pkpm_vel_flux_t)(
  const double *w, const double *dxv, const double *bvar, const double *grad_b, const double *u,
  const double *grad_u, const double *theta, const double *a_P, const double *f_l,
  const double *f_r, double *GKYL_RESTRICT flux
);

struct gkyl_dg_pkpm_vel_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  int num_surf_basis; // Surface modal basis size per component.
  bool use_gpu; // bool to determine if on GPU.
  pkpm_vel_flux_t vel_flux; // Flux kernel on the v_par face.

  uint32_t flags;
  struct gkyl_dg_pkpm_vel_flux_surf *on_dev; // pointer to itself or device data
};

// v_par flux kernels (p=1 serendipity only; indexed by cdim-1)
GKYL_CU_D static const pkpm_vel_flux_t ser_pkpm_vel_flux_kernels[] = {
  pkpm_vel_flux_vpar_1x1v_ser_p1,
  pkpm_vel_flux_vpar_2x1v_ser_p1,
  pkpm_vel_flux_vpar_3x1v_ser_p1,
};

GKYL_CU_D
static pkpm_vel_flux_t
choose_pkpm_vel_flux_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_vel_flux_kernels[cdim - 1];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}
