#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type
typedef struct gkyl_proj_bimaxwellian_on_basis gkyl_proj_bimaxwellian_on_basis;

// Object inputs packaged as a struct.
struct gkyl_proj_bimaxwellian_on_basis_inp {
  const struct gkyl_rect_grid *grid; // grid on which to project
  const struct gkyl_basis *conf_basis, *phase_basis; // basis functions
  int num_quad; // number of quadrature points
  bool use_gpu; // whether to use the GPU.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
};

/**
 * Create new updater to project a biMaxwellian on basis functions.
 * Free using gkyl_proj_bimaxwellian_on_basis_release method.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
gkyl_proj_bimaxwellian_on_basis* gkyl_proj_bimaxwellian_on_basis_inew(
  const struct gkyl_proj_bimaxwellian_on_basis_inp *inp);

/**
 * Create new updater to project a biMaxwellian on basis functions,
 * passing inputs separately. Free using the
 * gkyl_proj_bimaxwellian_on_basis_release method.
 *
 * @param grid Grid object
 * @param conf_basis Conf-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param num_quad Number of quadrature nodes (in 1D).
 * @param vel_map Velocity space mapping object.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_proj_bimaxwellian_on_basis* gkyl_proj_bimaxwellian_on_basis_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad, const struct gkyl_velocity_map *vel_map, bool use_gpu);

/**
 * Compute projection of a gyrokinetic Maxwellian on basis.
 * This method takes lab-frame moments (m0, m1, m2par, m2perp)
 * to compute the projection of Maxwellian on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param moms velocity moments (m0, m1i, m2par, m2perp)
 * @param bmag Magnetic field magnitude.
 * @param jacob_tot Total jacobian (conf * guiding center jacobian). 
 * @param mass Species mass.
 * @param fmax Output Maxwellian
 */
void gkyl_proj_bimaxwellian_on_basis_gyrokinetic_lab_mom(const gkyl_proj_bimaxwellian_on_basis *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *bmag,
  const struct gkyl_array *jacob_tot, double mass, struct gkyl_array *fmax);

/**
 * Compute projection of a gyrokinetic Maxwellian on basis. This
 * method takes primitive (fluid-frame) moments (m0, upar, Tpar,
 * Tperp) to compute the projection of Maxwellian on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param prim_moms primitive velocity moments (m0, upar, tpar, tperp)
 * @param bmag Magnetic field magnitude.
 * @param jacob_tot Total jacobian (conf * guiding center jacobian). 
 * @param mass Species mass.
 * @param fmax Output Maxwellian
 */
void gkyl_proj_bimaxwellian_on_basis_gyrokinetic_prim_mom(const gkyl_proj_bimaxwellian_on_basis *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *prim_moms, const struct gkyl_array *bmag,
  const struct gkyl_array *jacob_tot, double mass, struct gkyl_array *fmax);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_bimaxwellian_on_basis_release(gkyl_proj_bimaxwellian_on_basis *up);
