#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_mom_pkpm_surf_calc gkyl_mom_pkpm_surf_calc;

/**
 * Create new updater to compute moments of distribution
 * function. Free using gkyl_mom_pkpm_surf_calc_new_release.
 *
 * @param vel_grid Velocity space grid object
 * @param mass Mass of species
 * @return New updater pointer.
 */
struct gkyl_mom_pkpm_surf_calc* gkyl_mom_pkpm_surf_calc_new(const struct gkyl_rect_grid *vel_grid,
  double mass);

/**
 * Create new updater to compute moments of distribution function on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_mom_pkpm_surf_calc* gkyl_mom_pkpm_surf_calc_cu_dev_new(const struct gkyl_rect_grid *vel_grid,
  double mass);

/**
 * Compute moment of distribution function. The phase_rng and conf_rng
 * MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Surface Moment calculator updater to run
 * @param conf_edge_rng Range of configuration space edges (where the output is stored)
 * @param vel_rng Velocity space range
 * @param conf_cell_rng Config-space range (range over cells)
 * @param phase_cell_rng Phase-space range (range over cells)
 * @param uin Input velocity vector
 * @param bin Input magnetic field unit vector and tensor
 * @param fin Input distribution function array
 * @param mout Output moment array
 */
void gkyl_mom_pkpm_surf_calc_advance(const struct gkyl_mom_pkpm_surf_calc* calc,
  const struct gkyl_range *conf_edge_rng, const struct gkyl_range *vel_rng, 
  const struct gkyl_range *conf_cell_rng, const struct gkyl_range *phase_cell_rng,
  const struct gkyl_array *GKYL_RESTRICT uin, const struct gkyl_array *GKYL_RESTRICT bin, 
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout);

void gkyl_mom_pkpm_surf_calc_advance_cu(const struct gkyl_mom_pkpm_surf_calc* calc,
  const struct gkyl_range *conf_edge_rng, const struct gkyl_range *vel_rng, 
  const struct gkyl_range *conf_cell_rng, const struct gkyl_range *phase_cell_rng,
  const struct gkyl_array *GKYL_RESTRICT uin, const struct gkyl_array *GKYL_RESTRICT bin, 
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout);

/**
 * Delete pointer to moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_mom_pkpm_surf_calc_release(struct gkyl_mom_pkpm_surf_calc* calc);
