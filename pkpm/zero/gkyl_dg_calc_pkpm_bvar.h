#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

// Object type
typedef struct gkyl_dg_calc_pkpm_bvar gkyl_dg_calc_pkpm_bvar;

/**
 * Create new updater to compute the magnetic field unit vector b = B/|B| for the
 * PKPM model by direct normalization at the 2^cdim Gauss-Legendre nodes of each cell,
 * interpolated to the p=1 modal basis. No unit tensor b_i b_j is formed: consumers
 * build b_i b_j from the two p=1 components wherever the product is needed.
 *
 * Policy for a vanishing field: at a node where |B| = 0 the unit vector is set to
 * zero; the number of such nodes encountered in the last advance is available through
 * gkyl_dg_calc_pkpm_bvar_num_null_b. Only p=1 serendipity bases are supported.
 *
 * @param cbasis  Configuration space basis functions.
 * @param use_gpu bool to determine if on GPU.
 * @return New updater pointer.
 */
struct gkyl_dg_calc_pkpm_bvar *gkyl_dg_calc_pkpm_bvar_new(
  const struct gkyl_rect_grid *conf_grid, const struct gkyl_basis *cbasis,
  const struct gkyl_vlasov_position_map *pos_map, bool use_gpu
);

/**
 * Create new updater to compute the magnetic field unit vector on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_calc_pkpm_bvar *gkyl_dg_calc_pkpm_bvar_cu_dev_new(
  const struct gkyl_rect_grid *conf_grid, const struct gkyl_basis *cbasis, const struct gkyl_vlasov_position_map *pos_map
);

/**
 * Compute the magnetic field unit vector b = B/|B| (three p=1 components).
 *
 * @param up         Updater for computing the magnetic field unit vector.
 * @param conf_range Configuration space range to compute b over (the caller passes
 *                   the extended range when b is needed in ghost cells).
 * @param em         Input total electromagnetic field [Ex, Ey, Ez, Bx, By, Bz, phi, psi].
 * @param bvar       Output magnetic field unit vector [bx, by, bz].
 */
void gkyl_dg_calc_pkpm_bvar_advance(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, struct gkyl_array *bvar
);

/**
 * Number of Gauss-Legendre nodes at which |B| vanished in the last advance
 * (b was set to zero there).
 *
 * @param up Updater for computing the magnetic field unit vector.
 * @return Number of null-field nodes in the last advance.
 */
int gkyl_dg_calc_pkpm_bvar_num_null_b(const struct gkyl_dg_calc_pkpm_bvar *up);

/**
 * Compute the weak (surface-averaged) gradient of b, the p=1 expansions of d_i b_j for
 * every configuration direction i and component j, stored per cell as
 * [direction][component][mode] (cdim * 3 * num_basis entries), with the same kernels
 * and the same face-averaged form as the gradient of u (see
 * gkyl_dg_calc_pkpm_vars_grad_u). It needs b in the neighbors of every cell of
 * conf_range. The kinetic mirror terms form div b from it.
 *
 * @param up         Updater for computing the magnetic field unit vector.
 * @param conf_range Configuration space range (interior cells).
 * @param bvar       Input magnetic field unit vector, defined on the extended range.
 * @param grad_b     Output weak gradient of b.
 */
void gkyl_dg_calc_pkpm_bvar_grad(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, struct gkyl_array *grad_b
);

/**
 * Host-side wrapper for computing the weak gradient of b on device.
 */
void gkyl_dg_calc_pkpm_bvar_grad_cu(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, struct gkyl_array *grad_b
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_pkpm_bvar_release(struct gkyl_dg_calc_pkpm_bvar *up);

/**
 * Host-side wrapper for computing the magnetic field unit vector on device.
 */
void gkyl_dg_calc_pkpm_bvar_advance_cu(
  struct gkyl_dg_calc_pkpm_bvar *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *em, struct gkyl_array *bvar
);
