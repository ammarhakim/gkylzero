#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

// Object type
typedef struct gkyl_dg_pkpm_fluid_flux_surf gkyl_dg_pkpm_fluid_flux_surf;

// Input packaged as a struct
struct gkyl_dg_pkpm_fluid_flux_surf_inp {
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid (the velocity cells of the moments).
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space position map (per-direction Jacobians).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater computing the fluid face fluxes of the factored PKPM model as the
 * moments of the kinetic configuration-space flux (pkpm-fluid-flux-kinetic.mac): on the
 * lower face of every cell in every direction, and on the upper ghost-owned face of the
 * last interior cell, the HLL flux of the two-valued speed xi w b_n + u_n is formed per
 * normalized velocity and face node from the traces of [F, G] of the two cells, its
 * one-sided moments are accumulated over the velocity cells, and the fluxes of
 * [rho, M, E] are assembled from them with the physical traces (rho, u, xi, b) of each
 * side. The fluid system and the distributions thus see one discrete transport; with a
 * zero trace on one side (absorbing wall) the flux is the one-sided upwind flux, and with
 * mirrored traces (reflecting wall) the mass and energy fluxes vanish exactly. Each cell
 * owns its *lower* face; the stored layout per cell is [direction][component][surface
 * mode] with 5 * num_surf_basis entries per direction.
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_fluid_flux_surf_inp struct.
 * @return Pointer to the fluid surface flux updater.
 */
struct gkyl_dg_pkpm_fluid_flux_surf *gkyl_dg_pkpm_fluid_flux_surf_inew(
  const struct gkyl_dg_pkpm_fluid_flux_surf_inp *inp
);

/**
 * Create a new updater computing the fluid surface fluxes on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_pkpm_fluid_flux_surf *gkyl_dg_pkpm_fluid_flux_surf_cu_dev_inew(
  const struct gkyl_dg_pkpm_fluid_flux_surf_inp *inp
);

/**
 * Compute the fluid surface fluxes.
 *
 * @param up Fluid surface flux updater.
 * @param conf_range Configuration-space range to compute fluxes over (interior cells).
 * @param conf_range_ext Extended configuration-space range for indexing ghost cells.
 * @param phase_range Phase-space range (its velocity extent is summed over).
 * @param phase_range_ext Extended phase-space range for indexing ghost cells.
 * @param bvar Input magnetic field unit vector [bx, by, bz] (ghost cells filled).
 * @param u Input flow velocity [ux, uy, uz] (ghost cells filled).
 * @param xi Input thermal speed [xi, xi^2, p] (ghost cells filled).
 * @param theta Input T_perp/T = M0(G) (ghost cells filled).
 * @param fluid Input PHYSICAL fluid state [rho, rho ux, rho uy, rho uz, E] (Jacobian divided out; ghost cells filled).
 * @param f Input Jacobian-weighted normalized distributions J [F, G] (ghost cells filled).
 * @param cflrate Input cflrate. Accumulated to from the speed bounds (kinetic streaming and perpendicular sound) of the faces.
 * @param fluid_flux_surf Output modal fluxes [F_rho, F_x, F_y, F_z, F_E] per direction, for the evolved J-weighted fluid.
 */
void gkyl_dg_pkpm_fluid_flux_surf_advance(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_range_ext, const struct gkyl_array *bvar, const struct gkyl_array *u,
  const struct gkyl_array *xi, const struct gkyl_array *theta, const struct gkyl_array *fluid,
  const struct gkyl_array *f, struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf
);

/**
 * Host-side wrapper for computing the fluid fluxes on device.
 */
void gkyl_dg_pkpm_fluid_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_range_ext, const struct gkyl_array *bvar, const struct gkyl_array *u,
  const struct gkyl_array *xi, const struct gkyl_array *theta, const struct gkyl_array *fluid,
  const struct gkyl_array *f, struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf
);

/**
 * Number of surface basis functions per component and direction of the stored fluxes.
 *
 * @param up Fluid surface flux updater.
 * @return Surface basis size.
 */
int gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_fluid_flux_surf *up);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_fluid_flux_surf_release(struct gkyl_dg_pkpm_fluid_flux_surf *up);
