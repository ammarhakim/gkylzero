#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_calc_gyrokinetic_vars gkyl_dg_calc_gyrokinetic_vars;

/**
 * Create new updater to compute gyrokinetic variables needed in 
 * updates and used for diagnostics. Methods compute:
 * Bstar_Bmag : Volume expansion of time-independent component of Bstar/Bmag
 * alpha_surf : Surface expansion of phase space flux alpha
 * 
 * 
 * @param phase_grid Phase space grid (for getting cell spacing and cell center)
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis functions
 * @param charge Species charge
 * @param mass Species mass
 * @param gk_geom Geometry struct
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_dg_calc_gyrokinetic_vars* 
gkyl_dg_calc_gyrokinetic_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const double charge, const double mass, const struct gk_geometry *gk_geom, bool use_gpu);

/**
 * Create new updater to compute gyrokinetic variables on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_calc_gyrokinetic_vars* 
gkyl_dg_calc_gyrokinetic_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  double charge, double mass, const struct gk_geometry *gk_geom);

/**
 * Compute volume expansion of the time-independent component of Bstar/Bmag
 * 
 * @param up Updater for computing gyrokinetic variables 
 * @param conf_range Configuration space range (should only be local range because geometry only defined on local range)
 * @param phase_range Phase space range 
 * @param Bstar_Bmag Output volume expansion of Bstar/Bmag, time-independent component
 */
void gkyl_dg_calc_gyrokinetic_vars_Bstar_Bmag(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, struct gkyl_array* Bstar_Bmag);

/**
 * Compute surface expansion of phase space flux alpha
 * Computes the Poisson bracket of alpha = dz/dt = {z, H} and then evaluates the resulting polynomial
 * expansion at a surface and projects the evaluated quantity onto the surface basis.
 * 
 * Note: Each cell stores the surface expansion on the *lower* edge of the cell
 * @param up Updater for computing gyrokinetic variables 
 * @param conf_range Configuration space range (should only be local range because geometry only defined on local range)
 * @param phase_range Phase space range 
 * @param phase_ext_range Extended Phase space range (so we obtain alpha_surf at all the needed surfaces)
 * @param phi Electrostatic potential
 * @param phi Bstar/Bmag time-independent component pre-computed
 * @param alpha_surf Output surface expansion in a cell on the *lower* edge in each direction 
 */
void gkyl_dg_calc_gyrokinetic_vars_alpha_surf(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, const struct gkyl_range *phase_ext_range, 
  const struct gkyl_array *phi, const struct gkyl_array *Bstar_Bmag, struct gkyl_array* alpha_surf);

/**
 * Delete pointer to updater to compute gyrokinetic variables.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_calc_gyrokinetic_vars_release(struct gkyl_dg_calc_gyrokinetic_vars *up);

/**
 * Host-side wrappers for gyrokinetic vars operations on device
 */

void gkyl_dg_calc_gyrokinetic_vars_Bstar_Bmag_cu(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, struct gkyl_array* Bstar_Bmag);

void gkyl_dg_calc_gyrokinetic_vars_alpha_surf_cu(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, const struct gkyl_range *phase_ext_range, 
  const struct gkyl_array *phi, const struct gkyl_array *Bstar_Bmag, struct gkyl_array* alpha_surf);
