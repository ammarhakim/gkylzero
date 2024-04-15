#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_gyrokinetic_maxwellian_correct gkyl_gyrokinetic_maxwellian_correct;

// input packaged as a struct
struct gkyl_gyrokinetic_maxwellian_correct_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  const struct gkyl_range *conf_range; // Configuration-space range
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gkyl_range *vel_range; // velocity space range
  const struct gk_geometry *gk_geom; // Geometry object
  bool divide_jacobgeo; // Boolean for if we are dividing out the configuration-space Jacobian from density
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
                           // *even if* the scheme fails to converge. 
  double mass; // Mass factor 
  bool use_gpu; // bool for gpu useage
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
};

// Correction status
struct gkyl_gyrokinetic_maxwellian_correct_status {
  bool iter_converged; // true if iterations converged
  int num_iter; // number of iterations for the correction
  double error[3]; // error in each moment (n, u_par, T/m)
};  

/**
 * Create new updater to correct the gyrokinetic Maxwellian distribution function
 * so that its moments match desired input moments.
 *
 * @param inp Input parameters defined in gkyl_gyrokinetic_maxwellian_correct_inp struct.
 * @return New updater pointer.
 */
struct gkyl_gyrokinetic_maxwellian_correct* 
gkyl_gyrokinetic_maxwellian_correct_inew(const struct gkyl_gyrokinetic_maxwellian_correct_inp *inp);

/**
 * Fix the gyrokinetic Maxwellian distribution function so that its density matches a target density.
 *
 * @param up gyrokinetic Maxwellian distribution function moment correction updater
 * @param f_max gyrokinetic Maxwellian distribution function to fix (modified in-place)
 * @param moms_target Target density, n
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 */
void
gkyl_gyrokinetic_maxwellian_correct_density_moment(gkyl_gyrokinetic_maxwellian_correct *up,
  struct gkyl_array *f_max, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Fix the gyrokinetic Maxwellian distribution function
 * so that *all* its primitive moments (n, u_par, T/m) match target moments.
 * NOTE: If this algorithm fails, the returns the original distribution function
 * with only the desired density moment corrected.
 *
 * @param up gyrokinetic Maxwellian distribution function moment correction updater
 * @param f_max gyrokinetic Maxwellian distribution function to fix (modified in-place)
 * @param moms_target Target primitive moments (n, u_par, T/m)
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @return Status of correction
 */
struct gkyl_gyrokinetic_maxwellian_correct_status 
gkyl_gyrokinetic_maxwellian_correct_all_moments(gkyl_gyrokinetic_maxwellian_correct *up,
  struct gkyl_array *f_max, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Host-side wrapper for computing the absolute value of the 
 * difference in cell averages between the target moments and iterative moments.
 */
void gkyl_gyrokinetic_maxwellian_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int nc, const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_gyrokinetic_maxwellian_correct_release(gkyl_gyrokinetic_maxwellian_correct *up);
