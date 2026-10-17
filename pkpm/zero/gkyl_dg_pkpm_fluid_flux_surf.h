#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
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
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater computing the momentum surface fluxes of the PKPM model. On the
 * lower face of every cell in every direction, and on the upper ghost-owned face of the
 * last interior cell, the two traces of rho, p_par, p_perp, b, u and the momentum M are
 * evaluated at three Gauss-Legendre nodes per transverse direction, the pressure tensor
 * is formed pointwise from its primitive factors, and the ten-moment Roe momentum flux
 * of the nodal states (the HLL flux of the momentum with the normal speed bounds at a
 * node where a trace has rho <= 0 or P_nn <= 0; no Lax penalization anywhere) is
 * projected onto the surface basis, together with the mass flux of the same solve
 * (to which the kinetic configuration-space flux is constrained). Each cell owns its
 * *lower* face; the stored layout
 * per cell is [direction][component][surface mode] with 3 * num_surf_basis entries
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
 * @param moms Input kinetic moments [rho, p_parallel, p_perp, M1].
 * @param bvar Input magnetic field unit vector [bx, by, bz].
 * @param u Input flow velocity [ux, uy, uz].
 * @param euler Input momentum [rho ux, rho uy, rho uz].
 * @param cflrate Input cflrate. Accumulated to from the normal speed bound at the face nodes.
 * @param fluid_flux_surf Output modal momentum fluxes.
 * @param mass_flux_surf Output modal mass fluxes of the same Riemann solve,
 *                       [direction][surface mode] (cdim * num_surf_basis per cell).
 */
void gkyl_dg_pkpm_fluid_flux_surf_advance(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_array *moms,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *euler,
  struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf,
  struct gkyl_array *mass_flux_surf
);

/**
 * Host-side wrapper for computing momentum fluxes on device.
 */
void gkyl_dg_pkpm_fluid_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_array *moms,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *euler,
  struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf,
  struct gkyl_array *mass_flux_surf
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
