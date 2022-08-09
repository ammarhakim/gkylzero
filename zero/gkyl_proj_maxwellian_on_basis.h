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
 * @param num_quad Number of quadrature nodes
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
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param M0 Number density moment
 * @param M1i Momentum in lab-frame
 * @param M2 Energy in lab-frame
 * @param fmax Output Maxwellian
 */
void gkyl_proj_maxwellian_on_basis_lab_mom(const gkyl_proj_maxwellian_on_basis *pob,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *M0, const struct gkyl_array *M1i, const struct gkyl_array *M2,  
  struct gkyl_array *fmax);

/**
 * Compute projection of Maxwellian on basis. This method takes
 * primitive (fluid-frame) moments to compute the projection of
 * Maxwellian on basis functions.
 *
 * @param pob Project on basis updater to run
 * @param phase_rng Phase-space range
 * @param conf_rng Config-space range
 * @param n Number density moment
 * @param vel Velocity vector
 * @param vth2 Square of thermal velocity (vth2 = T/m)
 * @param fmax Output Maxwellian
 */
void gkyl_proj_maxwellian_on_basis_prim_mom(const gkyl_proj_maxwellian_on_basis *pob,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array *n, const struct gkyl_array *vel, const struct gkyl_array *vth2,  
  struct gkyl_array *fmax);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_proj_maxwellian_on_basis_release(gkyl_proj_maxwellian_on_basis* pob);
