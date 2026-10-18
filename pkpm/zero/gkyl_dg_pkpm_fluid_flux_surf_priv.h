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
  const double *bvar_r, const double *u_l, const double *u_r, const double *fluid_l,
  const double *fluid_r, double *GKYL_RESTRICT flux
);

typedef void (*pkpm_fluid_flux_absorb_accum_t)(
  int edge, const double *w, const double *dxv, const double *bvar, const double *u,
  const double *xi, const double *f, double *GKYL_RESTRICT I
);
typedef void (*pkpm_fluid_flux_absorb_t)(
  int edge, const double *bvar, const double *u, const double *xi, const double *fluid,
  const double *I, double *GKYL_RESTRICT flux
);

// for use in kernel tables
typedef struct {
  pkpm_fluid_flux_t kernels[3];
} gkyl_dg_pkpm_fluid_flux_kern_list;
typedef struct {
  pkpm_fluid_flux_absorb_accum_t kernels[3];
} gkyl_dg_pkpm_fluid_flux_absorb_accum_kern_list;
typedef struct {
  pkpm_fluid_flux_absorb_t kernels[3];
} gkyl_dg_pkpm_fluid_flux_absorb_kern_list;

struct gkyl_dg_pkpm_fluid_flux_surf {
  struct gkyl_rect_grid conf_grid; // Configuration-space grid.
  int cdim; // Configuration-space dimensions.
  int num_surf_basis; // Surface modal basis size per component and direction.
  bool use_gpu; // bool to determine if on GPU.
  const struct gkyl_wv_eqn *wv_eqn; // Ten-moment wave equation (device pointer on device).
  const struct gkyl_wave_geom *geom; // Wave geometry (device pointer on device).
  pkpm_fluid_flux_t fluid_flux[3]; // Flux kernel in each direction.
  struct gkyl_rect_grid phase_grid; // Phase-space grid (absorbing fluxes only).
  bool has_phase_grid; // Whether the phase-space grid was given (absorbing fluxes available).
  int num_face_nodes; // Transverse nodes of a face (3^(cdim-1)) of the absorbing-flux kernels.
  pkpm_fluid_flux_absorb_accum_t fluid_flux_absorb_accum[3]; // Upwind-flux moment kernel per direction.
  pkpm_fluid_flux_absorb_t fluid_flux_absorb[3]; // Absorbing-face flux kernel per direction.

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

// Absorbing-boundary flux kernels (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_fluid_flux_absorb_accum_kern_list ser_pkpm_fluid_flux_absorb_accum_kernels[] = {
  {pkpm_fluid_flux_absorb_accum_x_1x1v_ser_p1, NULL, NULL}, // 0
  {pkpm_fluid_flux_absorb_accum_x_2x1v_ser_p1, pkpm_fluid_flux_absorb_accum_y_2x1v_ser_p1, NULL}, // 1
  {pkpm_fluid_flux_absorb_accum_x_3x1v_ser_p1, pkpm_fluid_flux_absorb_accum_y_3x1v_ser_p1,
    pkpm_fluid_flux_absorb_accum_z_3x1v_ser_p1}, // 2
};
GKYL_CU_D static const gkyl_dg_pkpm_fluid_flux_absorb_kern_list ser_pkpm_fluid_flux_absorb_kernels[] = {
  {pkpm_fluid_flux_absorb_x_1x_ser_p1, NULL, NULL}, // 0
  {pkpm_fluid_flux_absorb_x_2x_ser_p1, pkpm_fluid_flux_absorb_y_2x_ser_p1, NULL}, // 1
  {pkpm_fluid_flux_absorb_x_3x_ser_p1, pkpm_fluid_flux_absorb_y_3x_ser_p1, pkpm_fluid_flux_absorb_z_3x_ser_p1}, // 2
};

GKYL_CU_D
static pkpm_fluid_flux_absorb_accum_t
choose_pkpm_fluid_flux_absorb_accum_kern(enum gkyl_basis_type b_type, int cdim, int poly_order, int dir)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_fluid_flux_absorb_accum_kernels[cdim - 1].kernels[dir];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

GKYL_CU_D
static pkpm_fluid_flux_absorb_t
choose_pkpm_fluid_flux_absorb_kern(enum gkyl_basis_type b_type, int cdim, int poly_order, int dir)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_fluid_flux_absorb_kernels[cdim - 1].kernels[dir];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

// Moments of the kinetic upwind flux through the absorbing face (dir, edge) of the skin
// cell idx, summed over the velocity cells of the phase-space range: I[4 * num_face_nodes].
GKYL_CU_DH static inline void
pkpm_fluid_flux_absorb_moments(
  const struct gkyl_dg_pkpm_fluid_flux_surf *up, int dir, int edge, const int *idx,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const double *b_c, const double *u_c, const double *xi_c, const struct gkyl_array *f, double *I
)
{
  int cdim = up->cdim;
  int pidx[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  for (int k = 0; k < 4 * up->num_face_nodes; ++k) {
    I[k] = 0.0;
  }
  gkyl_copy_int_arr(cdim, idx, pidx);
  for (int iv = phase_range->lower[cdim]; iv <= phase_range->upper[cdim]; ++iv) {
    pidx[cdim] = iv;
    long ploc = gkyl_range_idx(phase_range_ext, pidx);
    gkyl_rect_grid_cell_center(&up->phase_grid, pidx, xc);
    up->fluid_flux_absorb_accum[dir](
      edge, xc, up->phase_grid.dx, b_c, u_c, xi_c, (const double *)gkyl_array_cfetch(f, ploc), I
    );
  }
}
