#pragma once

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_basis.h>

// Object type.
typedef struct gkyl_dg_diffusion_gyrokinetic_proj_coeff gkyl_dg_diffusion_gyrokinetic_proj_coeff;

/**
 * Create an updater to project the diffusion coefficient
 *   g^{ij}*J*nu*(1+xi*v^2/v_t^2)
 * onto the phase basis, where g^{ij} is the contravariant metric, J is the
 * Jacobian, nu and xi are coefficients based on D and chi (see below), v^2
 * is the squared velocity, and v_t^2 is the squared thermal speed.
 *
 * @param cdim Configuration space dimensions.
 * @param pbasis Phase basis object.
 * @param grid Grid object.
 * @param D Particle diffusion coefficient for each direction.
 * @param chi Heat diffusivity for each direction.
 * @param mass Particle mass.
 * @param vtsq_min Minimum temperature supported by the grid.
 * @param diff_in_dir Whether diffusion is applied along a direction.
 * @param use_gpu Whether to run on the GPU.
 *
 * @return New updater to project the GK diffusion coefficient.
 */
struct gkyl_dg_diffusion_gyrokinetic_proj_coeff*
gkyl_dg_diffusion_gyrokinetic_proj_coeff_new(int cdim, struct gkyl_basis pbasis, struct gkyl_rect_grid *grid,
  const double *D, const double *chi, double mass, double vtsq_min, const bool *diff_in_dir, bool use_gpu);

/**
 * Project the diffusion coefficient onto the basis on the specified range.
 *
 * @param up Updater to run.
 * @param conf_rng Configuration space range to operate in.
 * @param vel_rng Velocity space range to operate in..
 * @param phase_rng Phase space range to operate in.
 * @param gijJ Contravariant metric multiplied by the Jacobian.
 * @param vmap Velocity space mapping.
 * @param vmapSq Square velocity space mapping.
 * @param bmag Magnetic field amplitude.
 * @param vtsq Thermal speed squared.
 * @param out Output diffusion coefficient.
 */
void
gkyl_dg_diffusion_gyrokinetic_proj_coeff_advance(gkyl_dg_diffusion_gyrokinetic_proj_coeff* up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *vel_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_array *GKYL_RESTRICT gijJ, struct gkyl_array *GKYL_RESTRICT vmap,
  struct gkyl_array *GKYL_RESTRICT vmapSq, struct gkyl_array *GKYL_RESTRICT bmag, struct gkyl_array *GKYL_RESTRICT vtsq,
  struct gkyl_array *GKYL_RESTRICT out);

/**
 * Free memory allocated for this updater.
 *
 * @param up Updater to free.
 */
void
gkyl_dg_diffusion_gyrokinetic_proj_coeff_release(gkyl_dg_diffusion_gyrokinetic_proj_coeff* up);
