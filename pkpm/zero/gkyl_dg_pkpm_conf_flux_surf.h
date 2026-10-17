#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_pkpm_conf_flux_surf gkyl_dg_pkpm_conf_flux_surf;

// Input packaged as a struct
struct gkyl_dg_pkpm_conf_flux_surf_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis functions (p=1 hybrid, one velocity dimension).
  bool acoustic_floor; // Widen the HLL speed bounds by the perpendicular acoustic speed
                       // sqrt(3 theta (1 - b_n^2)) (default false: plain upwind HLL).
  bool skip_domain_faces[GKYL_MAX_CDIM]; // Per direction: leave the two domain-boundary faces
                                         // out of the mass-flux constraint (non-periodic boundary
                                         // conditions, where the ghost state is a boundary
                                         // construct and the kinetic flux is the physical one).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater computing the configuration-space surface fluxes of the PKPM
 * kinetic equations. On the lower face of every cell in every configuration direction,
 * and on the upper ghost-owned face of the last interior cell, the two traces of the
 * streaming speed a = v_par b_n + u_n are evaluated at the face Gauss-Legendre nodes
 * (three per transverse configuration direction and in v_par) from the modal b and u
 * of the two cells, the HLL flux with the speed bounds S_L = min(a_l - c_l, a_r - c_r, 0),
 * S_R = max(a_l + c_l, a_r + c_r, 0) is formed for F_0 and G (c = 0 without the acoustic
 * floor, which gives max(a_l, 0) f_l + min(a_r, 0) f_r), and the nodal flux is projected
 * onto the surface modal basis. Each cell owns its *lower* face; the stored layout per
 * cell is [direction][component][surface mode] with 2 * num_surf_basis entries per
 * direction.
 *
 * The velocity integral of the F_0 flux (the kinetic mass flux through the face) is
 * also produced, per phase-space cell on the configuration face basis, and the
 * constrain() method corrects the stored fluxes so that this integral equals the mass
 * flux of the momentum solve: the density then evolves with the same Riemann solver as
 * the momentum, which removes the drift of u = M/rho from mismatched numerical fluxes.
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
 * Compute the configuration-space surface fluxes.
 *
 * @param up Configuration-space surface flux updater.
 * @param conf_range Configuration-space range for indexing b, u and theta (a sub-range of
 *                   the extended range, so ghost neighbors resolve through the same indexer).
 * @param phase_range Phase-space range to compute fluxes over (interior cells).
 * @param phase_range_ext Extended phase-space range for indexing ghost cells.
 * @param bvar Input magnetic field unit vector [bx, by, bz].
 * @param u Input flow velocity [ux, uy, uz].
 * @param theta Input p_perp/rho (used only with the acoustic floor).
 * @param fin Input distribution functions [F_0, G].
 * @param cflrate Input cflrate. Accumulated to from the maximum speed bound at the face nodes.
 * @param conf_flux_surf Output modal configuration-space fluxes.
 * @param mass_flux Output velocity integral of the F_0 flux of each phase-space cell on
 *                  the configuration face basis, [direction][face mode]
 *                  (cdim * num_face_basis entries per cell).
 */
void gkyl_dg_pkpm_conf_flux_surf_advance(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf,
  struct gkyl_array *mass_flux
);

/**
 * Host-side wrapper for computing configuration-space fluxes on device.
 */
void gkyl_dg_pkpm_conf_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf,
  struct gkyl_array *mass_flux
);

/**
 * Constrain the stored configuration-space fluxes to the mass flux of the momentum
 * solve. The per-cell mass fluxes of advance() are first summed over the velocity
 * column into mass_flux_face (cleared here); then, on every face, the defect
 * roe_mass_flux/mass - mass_flux_face is added to the F_0 and G fluxes with the shape
 * of the face-averaged distribution normalized by the face-averaged density, so that
 * the velocity integral of the corrected F_0 flux equals roe_mass_flux/mass. The
 * constraint is one-sided: the correction is applied only at face nodes where it adds
 * dissipation to the density (defect times density jump <= 0), i.e. where the momentum
 * solve dissipates the density more than the per-node kinetic upwinding does; where the
 * kinetic flux is the more dissipative one it is left alone (removing that dissipation
 * drives velocity cells of F_0 negative at steep gradients). No correction is applied at
 * a face node where the face-averaged density is not positive, nor on the
 * domain-boundary faces of a direction flagged in skip_domain_faces.
 *
 * @param up Configuration-space surface flux updater.
 * @param conf_range Configuration-space range (interior cells; sub-range of the extended range).
 * @param phase_range Phase-space range (interior cells).
 * @param phase_range_ext Extended phase-space range for indexing ghost cells.
 * @param mass Species mass (moms and roe_mass_flux carry the mass; the fluxes do not).
 * @param moms Input kinetic moments [rho, p_parallel, p_perp, M1] (valid in the ghost cells).
 * @param roe_mass_flux Input mass flux of the momentum solve, [direction][face mode] per
 *                      configuration cell (the lower faces; ghost cells hold the upper
 *                      domain faces), as filled by gkyl_dg_pkpm_fluid_flux_surf_advance.
 * @param mass_flux Input per-cell kinetic mass fluxes from advance().
 * @param mass_flux_face Scratch/output kinetic mass flux per configuration cell, same layout
 *                       as roe_mass_flux.
 * @param fin Input distribution functions [F_0, G].
 * @param conf_flux_surf Stored modal configuration-space fluxes, corrected in place.
 */
void gkyl_dg_pkpm_conf_flux_surf_constrain(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext, double mass,
  const struct gkyl_array *moms, const struct gkyl_array *roe_mass_flux,
  const struct gkyl_array *mass_flux, struct gkyl_array *mass_flux_face,
  const struct gkyl_array *fin, struct gkyl_array *conf_flux_surf
);

/**
 * Host-side wrapper for the mass-flux constraint on device.
 */
void gkyl_dg_pkpm_conf_flux_surf_constrain_cu(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext, double mass,
  const struct gkyl_array *moms, const struct gkyl_array *roe_mass_flux,
  const struct gkyl_array *mass_flux, struct gkyl_array *mass_flux_face,
  const struct gkyl_array *fin, struct gkyl_array *conf_flux_surf
);

/**
 * Number of surface modal basis functions per component and direction in the flux
 * array this updater fills (the array needs 2 * cdim times this many components per cell).
 *
 * @param up Updater.
 * @return Surface modal basis size.
 */
int gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_conf_flux_surf *up);

/**
 * Number of configuration face basis functions per direction in the mass-flux arrays
 * (the arrays need cdim times this many components per cell).
 *
 * @param up Updater.
 * @return Configuration face basis size.
 */
int gkyl_dg_pkpm_conf_flux_surf_num_face_basis(const struct gkyl_dg_pkpm_conf_flux_surf *up);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_conf_flux_surf_release(struct gkyl_dg_pkpm_conf_flux_surf *up);
