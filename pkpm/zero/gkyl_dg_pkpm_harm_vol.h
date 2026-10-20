#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_position_map.h>

// Object type
typedef struct gkyl_dg_pkpm_harm_vol gkyl_dg_pkpm_harm_vol;

// Input packaged as a struct
struct gkyl_dg_pkpm_harm_vol_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis (p=1 hybrid, one velocity dimension).
  double qbym; // Charge to mass ratio of the species (gyration of the first harmonic).
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space position map (REQUIRED; the flux divergence carries 1/J_d).
  bool use_gpu; // Whether to run on the GPU.
};

/**
 * Create a new updater for the volume integrals of the known first-harmonic terms of the
 * configuration-space pass (stage A) of the normalized PKPM kinetic update, on the eight blocks
 * [F, G, M0_perp (x, y, z), N_perp (x, y, z)]: the w-flux volume terms int d_w phi Phi, the
 * sources int phi S and the configuration-space flux couplings int d_x phi . T of the
 * first-harmonic terms (the two-field terms of F and G and the streaming of all eight blocks stay
 * in the configuration-pass equation object), with the definitions of gkyl_pkpm_harm_coeffs.h. The
 * integrals are exact and modal: per configuration cell every coefficient is built as a tensor
 * Legendre expansion one p=1 factor at a time (1D recurrences, direction by direction), keeping
 * only its projection onto degree <= 2 per direction, which is all that the integral against
 * (test function) x (distribution) sees; per phase cell the projections are contracted with the
 * distribution through exact 1D tables. p=1 serendipity only.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_dg_pkpm_harm_vol *gkyl_dg_pkpm_harm_vol_inew(const struct gkyl_dg_pkpm_harm_vol_inp *inp);

/**
 * Create the updater on NV-GPU. See the inew method for documentation.
 */
struct gkyl_dg_pkpm_harm_vol *gkyl_dg_pkpm_harm_vol_cu_dev_inew(
  const struct gkyl_dg_pkpm_harm_vol_inp *inp
);

/**
 * Accumulate the first-harmonic volume terms into rhs and their source rates into cflrate.
 *
 * @param up Updater.
 * @param conf_range Configuration-space range indexing the configuration fields.
 * @param phase_range Phase-space range of the cells to update.
 * @param bvar Magnetic field unit vector [bx, by, bz].
 * @param grad_b Weak gradient of b, [d_i bx, d_i by, d_i bz] per direction i.
 * @param u Flow velocity [ux, uy, uz].
 * @param grad_u Weak gradient of u.
 * @param xi Thermal speed [xi, xi^2, p] (the first expansion is used).
 * @param grad_xi Weak gradient of [xi, 0, 0].
 * @param grad_ln_n grad(ln n), cdim expansions.
 * @param theta T_perp/T.
 * @param dtB_over_magB d_t B/|B| of the stage.
 * @param em Total physical electromagnetic field of the stage (B at component 3).
 * @param fin Distributions [F, G, M0_perp, N_perp].
 * @param cflrate Phase-space time-step rates, incremented.
 * @param rhs Right-hand side, incremented.
 */
void gkyl_dg_pkpm_harm_vol_advance(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *xi,
  const struct gkyl_array *grad_xi, const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
);

/**
 * Host-side wrapper of the advance method on NV-GPU.
 */
void gkyl_dg_pkpm_harm_vol_advance_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *xi,
  const struct gkyl_array *grad_xi, const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_harm_vol_release(struct gkyl_dg_pkpm_harm_vol *up);
