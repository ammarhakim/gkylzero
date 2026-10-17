// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

typedef double (*pkpm_conf_flux_t)(
  const double *w, const double *dxv, const double *b_l, const double *b_r, const double *u_l,
  const double *u_r, const double *theta_l, const double *theta_r, const double *f_l,
  const double *f_r, double *GKYL_RESTRICT flux, double *GKYL_RESTRICT mass_flux
);

typedef void (*pkpm_conf_flux_corr_t)(
  double mass, const double *moms_l, const double *moms_r, const double *roe_mass_flux,
  const double *kin_mass_flux, const double *f_l, const double *f_r, double *GKYL_RESTRICT flux
);

// for use in kernel tables
typedef struct {
  pkpm_conf_flux_t kernels[3];
} gkyl_dg_pkpm_conf_flux_kern_list;

typedef struct {
  pkpm_conf_flux_corr_t kernels[3];
} gkyl_dg_pkpm_conf_flux_corr_kern_list;

struct gkyl_dg_pkpm_conf_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  int num_surf_basis; // Phase-space surface basis size per component and direction.
  int num_face_basis; // Configuration face basis size (mass fluxes) per direction.
  bool acoustic_floor; // Add the perpendicular acoustic speed to the HLL speed bounds.
  bool skip_domain_faces[GKYL_MAX_CDIM]; // Leave the domain-boundary faces unconstrained.
  bool use_gpu; // bool to determine if on GPU.
  pkpm_conf_flux_t conf_flux[3]; // Flux kernel in each configuration direction.
  pkpm_conf_flux_corr_t conf_flux_corr[3]; // Mass-flux constraint kernel in each direction.

  uint32_t flags;
  struct gkyl_dg_pkpm_conf_flux_surf *on_dev; // pointer to itself or device data
};

// Configuration-space HLL flux kernels (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_conf_flux_kern_list ser_pkpm_conf_flux_kernels[] = {
  {pkpm_conf_flux_x_1x1v_ser_p1, NULL, NULL}, // 0
  {pkpm_conf_flux_x_2x1v_ser_p1, pkpm_conf_flux_y_2x1v_ser_p1, NULL}, // 1
  {pkpm_conf_flux_x_3x1v_ser_p1, pkpm_conf_flux_y_3x1v_ser_p1, pkpm_conf_flux_z_3x1v_ser_p1}, // 2
};

// Mass-flux constraint kernels (same indexing)
GKYL_CU_D static const gkyl_dg_pkpm_conf_flux_corr_kern_list ser_pkpm_conf_flux_corr_kernels[] = {
  {pkpm_conf_flux_corr_x_1x1v_ser_p1, NULL, NULL}, // 0
  {pkpm_conf_flux_corr_x_2x1v_ser_p1, pkpm_conf_flux_corr_y_2x1v_ser_p1, NULL}, // 1
  {pkpm_conf_flux_corr_x_3x1v_ser_p1, pkpm_conf_flux_corr_y_3x1v_ser_p1,
    pkpm_conf_flux_corr_z_3x1v_ser_p1}, // 2
};

GKYL_CU_D
static pkpm_conf_flux_t
choose_pkpm_conf_flux_kern(enum gkyl_basis_type b_type, int cdim, int poly_order, int dir)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_conf_flux_kernels[cdim - 1].kernels[dir];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

GKYL_CU_D
static pkpm_conf_flux_corr_t
choose_pkpm_conf_flux_corr_kern(enum gkyl_basis_type b_type, int cdim, int poly_order, int dir)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_conf_flux_corr_kernels[cdim - 1].kernels[dir];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}
