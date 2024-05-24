#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_proj_maxwellian_on_basis gkyl_proj_maxwellian_on_basis;

/**
 * Create new updater to project Maxwellian on basis functions. Free
 * using gkyl_proj_maxwellian_on_basis_release method.
 *
 * @param grid Grid object
 * @param conf_basis Conf-space basis functions
 * @param phase_basis Phase-space basis functions
 * @param num_quad Number of quadrature nodes (in 1D).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_proj_maxwellian_on_basis* gkyl_proj_maxwellian_on_basis_new(
  const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  int num_quad, bool use_gpu);

/**
 * Compute projection of Maxwellian on basis. This method takes
 * lab-frame moments to compute the projection of Maxwellian on basis
 * functions.
 *
 * @param mob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param moms velocity moments (m0, m1i, m2)
 * @param fmax Output Maxwellian
 */
void gkyl_proj_maxwellian_on_basis_lab_mom(const gkyl_proj_maxwellian_on_basis *mob,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms, struct gkyl_array *fmax);

/**
 * Compute projection of Maxwellian on basis. This method takes
 * primitive (fluid-frame) moments to compute the projection of
 * Maxwellian on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param moms velocity moments (m0, m1i, m2)
 * @param prim_moms (primitive moments udrift, vtsq=T/m)
 * @param fmax Output Maxwellian
 */
void gkyl_proj_maxwellian_on_basis_prim_mom(const gkyl_proj_maxwellian_on_basis *mob,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *moms, const struct gkyl_array *prim_moms, struct gkyl_array *fmax);

/**
 * Compute projection of a gyrokinetic Maxwellian on basis.
 * This method takes lab-frame moments to compute the projection
 * of Maxwellian on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param moms velocity moments (m0, m1i, m2)
 * @param bmag Magnetic field magnitude.
 * @param jacob_tot Total jacobian (conf * guiding center jacobian). 
 * @param mass Species mass.
 * @param fmax Output Maxwellian
 */
void gkyl_proj_gkmaxwellian_on_basis_lab_mom(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *bmag,
  const struct gkyl_array *jacob_tot, double mass, struct gkyl_array *fmax);

/**
 * Compute projection of a gyrokinetic Maxwellian on basis. This
 * method takes primitive (fluid-frame) moments to compute the
 * projection of Maxwellian on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param prim_moms (primitive moments n, upar, vtsq=T/m)
 * @param bmag Magnetic field magnitude.
 * @param jacob_tot Total jacobian (conf * guiding center jacobian). 
 * @param mass Species mass.
 * @param fmax Output Maxwellian
 */
void gkyl_proj_gkmaxwellian_on_basis_prim_mom(const gkyl_proj_maxwellian_on_basis *up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *prim_moms,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, double mass,
  struct gkyl_array *fmax);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_maxwellian_on_basis_release(gkyl_proj_maxwellian_on_basis* mob);
