#pragma once

#include <gkyl_dg_vlasov_pkpm.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_flash gkyl_dg_updater_flash;

struct gkyl_dg_flash_auxfields {};

/**
 * Initialize updater for advancing the kinetic portion of the hybrid system.
 * 
 * @param conf_grid 
 * @param phase_grid 
 * @param conf_basis 
 * @param phase_basis 
 * @param conf_range 
 * @param phase_range 
 * @param is_zero_flux_dir 
 * @param aux 
 * @return gkyl_dg_updater_flash* 
 */
gkyl_dg_updater_flash*
gkyl_dg_updater_flash_new(const struct gkyl_rect_grid *conf_grid,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_basis *conf_basis, 
  const struct gkyl_basis *phase_basis, const struct gkyl_range *conf_range, 
  const struct gkyl_range *phase_range, const bool *is_zero_flux_dir,
  struct gkyl_dg_vlasov_pkpm_auxfields *aux);

/**
 * Compute right-hand-side of PKPM kinetic update.
 * 
 * @param up 
 * @param update_phase_rng 
 * @param update_conf_rng 
 * @param fIn 
 * @param cflrate_f 
 * @param rhs_f 
 */
void
gkyl_dg_updater_flash_advance(gkyl_dg_updater_flash *up, 
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng, 
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate_f,
  struct gkyl_array* GKYL_RESTRICT rhs_f);

/**
 * Compute maximum stable timestep for kinetic update without performing the update.
 * Useful for external solvers that set timestep before update call.
 * 
 * @param up 
 * @param phase_grid 
 * @param update_phase_rng 
 * @param update_conf_rng 
 * @param pkpm_accel 
 * @param pkpm_lax 
 * @param div_b 
 * @param max_b 
 * @param cflrate 
 * @param cfl 
 * @return double 
 */

void
gkyl_dg_updater_flash_compute_cflrate(gkyl_dg_updater_flash *up,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_range *update_phase_rng,
  const struct gkyl_range *update_conf_rng, const struct gkyl_array *pkpm_accel, 
  const struct gkyl_array *pkpm_lax, const struct gkyl_array *div_b, 
  const struct gkyl_array *max_b, struct gkyl_array *cflrate);

void
gkyl_dg_updater_flash_compute_cflrate_lbo(gkyl_dg_updater_flash *up,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_range *update_phase_rng,
  const struct gkyl_range *update_conf_rng, const struct gkyl_array *nuSum, 
  const struct gkyl_array *nuPrimMomsSum, struct gkyl_array *cflrate);

void
gkyl_dg_updater_flash_release(gkyl_dg_updater_flash *up);