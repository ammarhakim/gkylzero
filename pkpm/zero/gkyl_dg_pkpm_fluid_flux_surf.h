#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Object type
typedef struct gkyl_dg_pkpm_fluid_flux_surf gkyl_dg_pkpm_fluid_flux_surf;

// Input packaged as a struct
struct gkyl_dg_pkpm_fluid_flux_surf_inp {
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_wv_eqn *wv_eqn; // Ten-moment wave equation for the Roe fluctuations.
  const struct gkyl_wave_geom *geom; // Wave geometry (face normals and tangents).
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid (optional; needed for absorbing fluxes).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater computing the fluid surface fluxes of the factored PKPM model. On
 * the lower face of every cell in every direction, and on the upper ghost-owned face of
 * the last interior cell, the two traces of the fluid state [M, rho], the normalized
 * kinetic moments, b and u are evaluated at three Gauss-Legendre nodes per transverse
 * direction, the pressure tensor is formed pointwise from its factors, and the ten-moment
 * Roe momentum and mass flux
 * of the nodal states (the HLL flux of the momentum with the normal speed bounds at a
 * node where a trace has rho <= 0 or P_nn <= 0; no Lax penalization anywhere) is
 * projected onto the surface basis. Each cell owns its *lower* face; the stored layout
 * per cell is [direction][component][surface mode] with 5 * num_surf_basis entries
 * per direction.
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_fluid_flux_surf_inp struct.
 * @return Pointer to momentum surface flux updater.
 */
struct gkyl_dg_pkpm_fluid_flux_surf *gkyl_dg_pkpm_fluid_flux_surf_inew(
  const struct gkyl_dg_pkpm_fluid_flux_surf_inp *inp
);

/**
 * Create a new updater computing the momentum surface fluxes on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_pkpm_fluid_flux_surf *gkyl_dg_pkpm_fluid_flux_surf_cu_dev_inew(
  const struct gkyl_dg_pkpm_fluid_flux_surf_inp *inp
);

/**
 * Compute the momentum surface fluxes.
 *
 * @param up Momentum surface flux updater.
 * @param conf_range Configuration-space range to compute fluxes over (interior cells).
 * @param conf_range_ext Extended configuration-space range for indexing ghost cells.
 * @param vlasov_pkpm_moms Input kinetic coupling moments [p_par, p_perp, q_par, q_perp].
 * @param bvar Input magnetic field unit vector [bx, by, bz].
 * @param u Input flow velocity [ux, uy, uz].
 * @param fluid Input fluid state [rho, rho ux, rho uy, rho uz, E].
 * @param cflrate Input cflrate. Accumulated to from the normal speed bound at the face nodes.
 * @param fluid_flux_surf Output modal fluxes [F_rho, F_x, F_y, F_z, F_E] per direction.
 */
void gkyl_dg_pkpm_fluid_flux_surf_advance(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *fluid,
  struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf
);

/**
 * Host-side wrapper for computing momentum fluxes on device.
 */
void gkyl_dg_pkpm_fluid_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *fluid,
  struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf
);

/**
 * Overwrite the stored fluxes on the absorbing boundary faces (dir, edge) with the
 * fluxes assembled from the moments of the kinetic upwind flux of the skin cells
 * (pkpm-fluid-flux-absorb.mac): the fluid loses exactly the mass, momentum and energy
 * the distribution functions lose through the wall. The flux of a lower boundary is
 * stored in the skin cell (its lower face), that of an upper boundary in the ghost cell.
 * The updater must have been created with the phase-space grid.
 *
 * @param up Fluid surface flux updater.
 * @param dir Boundary direction.
 * @param edge Boundary edge.
 * @param conf_skin Skin cells of the boundary (configuration space).
 * @param conf_range_ext Extended configuration-space range for indexing.
 * @param phase_range Phase-space range (interior; its velocity extent is used).
 * @param phase_range_ext Extended phase-space range for indexing.
 * @param bvar Input magnetic field unit vector [bx, by, bz].
 * @param u Input flow velocity [ux, uy, uz].
 * @param xi Input thermal speed [xi, xi^2, p].
 * @param fluid Input fluid state [rho, rho ux, rho uy, rho uz, E].
 * @param f Input normalized distribution functions [F, G].
 * @param fluid_flux_surf Modal fluxes, overwritten on the boundary faces.
 */
void gkyl_dg_pkpm_fluid_flux_surf_absorb(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *conf_skin, const struct gkyl_range *conf_range_ext,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *fluid, const struct gkyl_array *f, struct gkyl_array *fluid_flux_surf
);

/**
 * Host-side wrapper for the absorbing-boundary fluxes on device.
 */
void gkyl_dg_pkpm_fluid_flux_surf_absorb_cu(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *conf_skin, const struct gkyl_range *conf_range_ext,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *fluid, const struct gkyl_array *f, struct gkyl_array *fluid_flux_surf
);

/**
 * Number of surface modal basis functions per component and direction in the flux
 * array this updater fills (the array needs 3 * cdim times this many components per cell).
 *
 * @param up Updater.
 * @return Surface modal basis size.
 */
int gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_fluid_flux_surf *up);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_fluid_flux_surf_release(struct gkyl_dg_pkpm_fluid_flux_surf *up);
