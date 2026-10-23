// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_fluid_kernels.h>
#include <gkyl_pkpm_harm_coeffs.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

typedef double (*pkpm_fluid_flux_accum_t)(
  const double *w, const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *fc, const double *f_l, const double *f_r, double *GKYL_RESTRICT I
);
typedef void (*pkpm_fluid_flux_t)(
  const double *bvar_l, const double *bvar_r, const double *u_l, const double *u_r, const double *xi_l,
  const double *xi_r, const double *fluid_l, const double *fluid_r, const double *I,
  double *GKYL_RESTRICT flux
);

// Nodal moments per side and transverse face node: [I_0, I_1, I_2, I_G, I_M (3)].
#define PKPM_FLUID_FLUX_NUM_MOMS 7
// Largest number of transverse face nodes (5^(cdim-1), cdim <= 3).
#define PKPM_FLUID_FLUX_MAX_NODES 25

// for use in kernel tables
typedef struct {
  pkpm_fluid_flux_accum_t kernels[3];
} gkyl_dg_pkpm_fluid_flux_accum_kern_list;
typedef struct {
  pkpm_fluid_flux_t kernels[3];
} gkyl_dg_pkpm_fluid_flux_kern_list;

struct gkyl_dg_pkpm_fluid_flux_surf {
  struct gkyl_rect_grid conf_grid; // Configuration-space grid.
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int num_surf_basis; // Surface modal basis size per component and direction.
  int num_face_nodes; // Transverse nodes of a face (5^(cdim-1)).
  int num_face_coef; // Face coefficients per direction (PKPM_FC_NUM times num_face_nodes).
  bool use_gpu; // bool to determine if on GPU.
  const struct gkyl_vlasov_position_map *pos_map; // Position map (acquired; host pointer).
  const struct gkyl_array *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device).
  pkpm_fluid_flux_accum_t fluid_flux_accum[3]; // Kinetic-flux moment kernel per direction.
  pkpm_fluid_flux_t fluid_flux[3]; // Flux assembly kernel per direction.
  uint32_t flags;
  struct gkyl_dg_pkpm_fluid_flux_surf *on_dev; // pointer to itself or device data
};

// Kinetic-flux moment kernels (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_fluid_flux_accum_kern_list ser_pkpm_fluid_flux_accum_kernels[] = {
  {pkpm_fluid_flux_accum_x_1x1v_ser_p1, NULL, NULL}, // 0
  {pkpm_fluid_flux_accum_x_2x1v_ser_p1, pkpm_fluid_flux_accum_y_2x1v_ser_p1, NULL}, // 1
  {pkpm_fluid_flux_accum_x_3x1v_ser_p1, pkpm_fluid_flux_accum_y_3x1v_ser_p1,
    pkpm_fluid_flux_accum_z_3x1v_ser_p1}, // 2
};

// Flux assembly kernels (p=1 serendipity only; indexed [cdim-1].kernels[dir])
GKYL_CU_D static const gkyl_dg_pkpm_fluid_flux_kern_list ser_pkpm_fluid_flux_kernels[] = {
  {pkpm_fluid_flux_x_1x_ser_p1, NULL, NULL}, // 0
  {pkpm_fluid_flux_x_2x_ser_p1, pkpm_fluid_flux_y_2x_ser_p1, NULL}, // 1
  {pkpm_fluid_flux_x_3x_ser_p1, pkpm_fluid_flux_y_3x_ser_p1, pkpm_fluid_flux_z_3x_ser_p1}, // 2
};

GKYL_CU_D
static pkpm_fluid_flux_accum_t
choose_pkpm_fluid_flux_accum_kern(enum gkyl_basis_type b_type, int cdim, int poly_order, int dir)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_fluid_flux_accum_kernels[cdim - 1].kernels[dir];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

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

// Flux of the face between the configuration cells idx_l (lower) and idx_r (upper) in
// direction dir: the one-sided moments of the kinetic HLL flux of the eight blocks (from the face
// coefficients owned by idx_r) summed over the velocity cells of phase_range, then the moment map
// with the physical traces of both cells.
// Either cell may be a ghost cell (indexed through the extended ranges). Returns the CFL rate of
// the r cell from the widest speed bounds of the face.
GKYL_CU_DH static inline double
pkpm_fluid_flux_face(
  const struct gkyl_dg_pkpm_fluid_flux_surf *up, int dir, const int *idx_l, const int *idx_r,
  const struct gkyl_range *conf_range_ext, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_range_ext, const struct gkyl_array *bvar, const struct gkyl_array *u,
  const struct gkyl_array *xi, const struct gkyl_array *face_coef, const struct gkyl_array *fluid,
  const struct gkyl_array *f, double *flux
)
{
  int cdim = up->cdim;
  int pidx_l[GKYL_MAX_DIM], pidx_r[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  double I[2 * PKPM_FLUID_FLUX_NUM_MOMS * PKPM_FLUID_FLUX_MAX_NODES];
  for (int k = 0; k < 2 * PKPM_FLUID_FLUX_NUM_MOMS * up->num_face_nodes; ++k) {
    I[k] = 0.0;
  }
  long cidx_l = gkyl_range_idx(conf_range_ext, idx_l), cidx_r = gkyl_range_idx(conf_range_ext, idx_r);
  const double *jac_l = (const double *)gkyl_array_cfetch(up->jacob_pos, cidx_l);
  const double *jac_r = (const double *)gkyl_array_cfetch(up->jacob_pos, cidx_r);
  const double *b_l = (const double *)gkyl_array_cfetch(bvar, cidx_l), *b_r = (const double *)gkyl_array_cfetch(bvar, cidx_r);
  const double *u_l = (const double *)gkyl_array_cfetch(u, cidx_l), *u_r = (const double *)gkyl_array_cfetch(u, cidx_r);
  const double *xi_l = (const double *)gkyl_array_cfetch(xi, cidx_l), *xi_r = (const double *)gkyl_array_cfetch(xi, cidx_r);
  const double *fc = &((const double *)gkyl_array_cfetch(face_coef, cidx_r))[dir * up->num_face_coef];
  double amax = 0.0;
  gkyl_copy_int_arr(cdim, idx_l, pidx_l);
  gkyl_copy_int_arr(cdim, idx_r, pidx_r);
  for (int jv = phase_range->lower[cdim]; jv <= phase_range->upper[cdim]; ++jv) {
    pidx_l[cdim] = jv;
    pidx_r[cdim] = jv;
    gkyl_rect_grid_cell_center(&up->phase_grid, pidx_r, xc);
    amax = fmax(amax, up->fluid_flux_accum[dir](
      xc, up->phase_grid.dx, jac_l, jac_r, fc,
      (const double *)gkyl_array_cfetch(f, gkyl_range_idx(phase_range_ext, pidx_l)),
      (const double *)gkyl_array_cfetch(f, gkyl_range_idx(phase_range_ext, pidx_r)), I
    ));
  }
  up->fluid_flux[dir](
    b_l, b_r, u_l, u_r, xi_l, xi_r, (const double *)gkyl_array_cfetch(fluid, cidx_l),
    (const double *)gkyl_array_cfetch(fluid, cidx_r), I, flux
  );
  // CFL rate of the r cell: (2p+1)/2 dx1 / J_d, as the kinetic flux
  return 1.5 * (2.0 / up->conf_grid.dx[dir]) / jac_r[2 * dir] * amax;
}
