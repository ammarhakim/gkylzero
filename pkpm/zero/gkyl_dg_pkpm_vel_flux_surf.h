#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_pkpm_vel_flux_surf gkyl_dg_pkpm_vel_flux_surf;

// Input packaged as a struct
struct gkyl_dg_pkpm_vel_flux_surf_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid.
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions (p=1 serendipity).
  const struct gkyl_basis *phase_basis; // Phase-space basis functions (p=1 hybrid, one velocity dimension).
  bool use_gpu; // bool to determine if on GPU.
};

/**
 * Create a new updater computing the v_par surface fluxes of the PKPM kinetic
 * equations (velocity-space pass, stage B). On the lower v_par face of every cell the
 * single-valued configuration-space coefficients alpha = a_P/xi + (c - bb:grad(u)) w_face - (b.grad xi) w_face^2,
 * beta = xi div(b) and theta are evaluated at the 3^cdim face nodes, the split upwind
 * fluxes alpha [F_0, G] + beta [G, H] with H = 4 theta G - 2 theta^2 F_0 are formed,
 * and the nodal fluxes are projected onto the surface modal basis. Each cell owns its
 * *lower* face, stored as [component][surface mode]; the faces at the two velocity
 * domain edges carry zero flux (written explicitly).
 *
 * @param inp Input parameters defined in gkyl_dg_pkpm_vel_flux_surf_inp struct.
 * @return Pointer to v_par surface flux updater.
 */
struct gkyl_dg_pkpm_vel_flux_surf *gkyl_dg_pkpm_vel_flux_surf_inew(
  const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp
);

/**
 * Create a new updater computing the v_par surface fluxes on NV-GPU.
 * See new() method for documentation.
 */
struct gkyl_dg_pkpm_vel_flux_surf *gkyl_dg_pkpm_vel_flux_surf_cu_dev_inew(
  const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp
);

/**
 * Compute the v_par surface fluxes on the INTERIOR faces (the lower face of every cell but
 * the lowest and of the velocity ghost cell; the domain-edge faces are set by the closure updater gkyl_dg_calc_pkpm_vars_closure and left
 * untouched here).
 *
.
 *
 * @param up v_par surface flux updater.
 * @param conf_range Configuration-space range for indexing the configuration fields.
 * @param phase_range Phase-space range to compute fluxes over.
 * @param bvar Input magnetic field unit vector [bx, by, bz].
 * @param grad_b Input weak gradient of b (gkyl_dg_calc_pkpm_bvar_grad).
 * @param u Input flow velocity [ux, uy, uz].
 * @param grad_u Input weak gradient of u (gkyl_dg_calc_pkpm_vars_grad_u).
 * @param xi Input thermal speed [xi, xi^2].
 * @param grad_xi Input weak gradient of [xi, 0, 0].
 * @param theta Input T_perp/T.
 * @param alpha_0 Input pressure force a_P/xi closing the first moment.
 * @param c_h Input heating rate closing the scalar second moment.
 * @param fin Input distribution functions [F_0, G].
 * @param cflrate Input cflrate. Accumulated to from the maximum speed at the face nodes.
 * @param vel_flux_surf Output modal v_par fluxes.
 */
void gkyl_dg_pkpm_vel_flux_surf_advance(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *theta, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf
);

/**
 * Host-side wrapper for computing v_par fluxes on device.
 */
void gkyl_dg_pkpm_vel_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *theta, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf
);

/**
 * Number of surface modal basis functions per component in the flux array this
 * updater fills (the array needs 2 times this many components per cell).
 *
 * @param up Updater.
 * @return Surface modal basis size.
 */
int gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_vel_flux_surf *up);




/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_pkpm_vel_flux_surf_release(struct gkyl_dg_pkpm_vel_flux_surf *up);
