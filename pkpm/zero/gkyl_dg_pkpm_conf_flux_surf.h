#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

// Object type
typedef struct gkyl_dg_pkpm_conf_flux_surf gkyl_dg_pkpm_conf_flux_surf;

// Input packaged as a struct
struct gkyl_dg_pkpm_conf_flux_surf_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis functions (p=1 hybrid, one velocity dimension).
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space position map (REQUIRED; per-direction Jacobians of the cells).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater computing the configuration-space surface fluxes of the PKPM
 * kinetic equations of the eight blocks [F, G, M0_perp, N_perp] (plan 17.14). On the lower
 * face of every cell in every configuration direction, and on the upper ghost-owned face of
 * the last interior cell, the traces of the blocks are evaluated at the face Gauss-Legendre nodes
 * (five per transverse configuration direction, three in w) from the modal fields of the two
 * cells, the coefficients u_n, xi, b, theta at the same transverse nodes once per configuration
 * face (gkyl_dg_pkpm_conf_flux_surf_face_coef), and the HLL flux of
 * J(X) = a X + T(X), a = {u_n} + {xi} b_n w, with every coefficient at its face average
 * {q} = (q_l + q_r)/2 and the field direction b = {b}/|{b}| (T the first-harmonic coupling
 * {xi} M0_n, {xi} N_n, {xi} G (n - b b_n), {xi} H (n - b b_n), H with {theta}) and the bounds
 * a -+ c of this frozen-coefficient flux, c = {xi} |b x n| sqrt((2 + sqrt 2) {theta}), is formed
 * (pkpm_harm_conf_flux, plans 17.41, 17.42); the time step sees the envelope of both sides'
 * speeds u_n + xi w b_n -+ c and of a -+ c. The nodal flux is projected onto the surface modal
 * basis. Each cell owns its *lower* face; the
 * stored layout per cell is [direction][block][surface mode] with 8 * num_surf_basis entries
 * per direction.
 *
 * In the factored update the distributions are normalized: the flux moves only the shape;
 * the density is moved by the fluid system. The stored
 * distributions carry the position-map Jacobian J of their cell: the kernels divide
 * the traces by J_d of the owning cell, which equals forming the physical flux of the
 * physical traces and multiplying by the face Jacobian J/J_d.
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_conf_flux_surf_inp struct.
 * @return Pointer to configuration-space surface flux updater.
 */
struct gkyl_dg_pkpm_conf_flux_surf *gkyl_dg_pkpm_conf_flux_surf_inew(
  const struct gkyl_dg_pkpm_conf_flux_surf_inp *inp
);

/**
 * Create a new updater computing the configuration-space surface fluxes on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_pkpm_conf_flux_surf *gkyl_dg_pkpm_conf_flux_surf_cu_dev_inew(
  const struct gkyl_dg_pkpm_conf_flux_surf_inp *inp
);

/**
 * Compute the face coefficients of the configuration-space flux once per configuration face
 * (plan 17.42): at the transverse face nodes of the lower face of every cell of conf_range in every
 * direction, and of the upper ghost-owned face of the last interior cell, the face averages of
 * u_n, xi, b and theta that the flux function uses and the per-side speeds of the bounds
 * (pkpm_harm_face_coef, gkyl_pkpm_harm_coeffs.h). Each cell owns its *lower* face; the layout per
 * cell is [direction][coefficient][node], with gkyl_dg_pkpm_conf_flux_surf_face_coef_ncomp
 * components. The kinetic flux of every velocity cell of the face and the fluid flux read them.
 *
 * @param up Configuration-space surface flux updater.
 * @param conf_range Configuration-space range of the interior cells (a sub-range of the extended
 *                   range, so ghost neighbors resolve through the same indexer).
 * @param bvar Input magnetic field unit vector [bx, by, bz].
 * @param u Input flow velocity [ux, uy, uz].
 * @param xi Input thermal speed [xi, xi^2, p] (the streaming speed is xi w b_n + u_n).
 * @param theta Input T_perp/T (the closure H and the speed bounds).
 * @param face_coef Output face coefficients.
 */
void gkyl_dg_pkpm_conf_flux_surf_face_coef(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, struct gkyl_array *face_coef
);

/**
 * Host-side wrapper for computing the face coefficients on device.
 */
void gkyl_dg_pkpm_conf_flux_surf_face_coef_cu(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, struct gkyl_array *face_coef
);

/**
 * Compute the configuration-space surface fluxes.
 *
 * @param up Configuration-space surface flux updater.
 * @param conf_range Configuration-space range for indexing the face coefficients (a sub-range of
 *                   the extended range, so ghost neighbors resolve through the same indexer).
 * @param phase_range Phase-space range to compute fluxes over (interior cells).
 * @param phase_range_ext Extended phase-space range for indexing ghost cells.
 * @param face_coef Input face coefficients (gkyl_dg_pkpm_conf_flux_surf_face_coef).
 * @param fin Input normalized distributions [F, G, M0_perp, N_perp].
 * @param cflrate Input cflrate. Accumulated to from the maximum speed bound at the face nodes.
 * @param conf_flux_surf Output modal configuration-space fluxes.
 */
void gkyl_dg_pkpm_conf_flux_surf_advance(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *face_coef, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *conf_flux_surf
);

/**
 * Host-side wrapper for computing configuration-space fluxes on device.
 */
void gkyl_dg_pkpm_conf_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *face_coef, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *conf_flux_surf
);

/**
 * Number of components per cell of the face-coefficient array (cdim directions times
 * PKPM_FC_NUM coefficients times 5^(cdim-1) transverse face nodes).
 *
 * @param up Updater.
 * @return Number of components.
 */
int gkyl_dg_pkpm_conf_flux_surf_face_coef_ncomp(const struct gkyl_dg_pkpm_conf_flux_surf *up);

/**
 * Number of surface modal basis functions per component and direction in the flux
 * array this updater fills (the array needs 8 * cdim times this many components per cell).
 *
 * @param up Updater.
 * @return Surface modal basis size.
 */
int gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_conf_flux_surf *up);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_conf_flux_surf_release(struct gkyl_dg_pkpm_conf_flux_surf *up);
