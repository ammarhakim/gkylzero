#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_correct_maxwellian_gyrokinetic gkyl_correct_maxwellian_gyrokinetic;

// Input struct
struct gkyl_correct_maxwellian_gyrokinetic_inp
{
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions

  const struct gkyl_range *conf_local; // local Configuration-space range
  const struct gkyl_range *conf_local_ext; // extended Configuration-space range
  double mass; // species mass
  const struct gk_geometry *gk_geom; // geometry struct
  long max_iter; // maximum allowed iterations
  double eps_err; // desired error tolerance
  bool use_gpu; // Boolean for if updater is on GPU

  struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
};


/**
 * Create new updater to correct a Maxwellian to match specified
 * moments.
 *
 * @param inp Input struct 
 * @return Pointer to gyrokinetic correct maxwellian struct
 */
struct gkyl_correct_maxwellian_gyrokinetic* 
gkyl_correct_maxwellian_gyrokinetic_new(const struct gkyl_correct_maxwellian_gyrokinetic_inp *inp);

/**
 * Fix the Maxwellian so that it's moments match desired moments.
 *
 * @param cmax Maxwellian-fix updater
 * @param fM Distribution function to fix (modified in-place)
 * @param moms_tar Target moments to match
 * @param conf_local Local configuration space range
 * @param phase_local Local phase-space range
 */
void gkyl_correct_maxwellian_gyrokinetic_advance(gkyl_correct_maxwellian_gyrokinetic *cmax,
  struct gkyl_array *fM, const struct gkyl_array *moms_tar, 
  const struct gkyl_range *conf_local, const struct gkyl_range *phase_local);

/**
 * Delete updater.
 *
 * @param cmax Updater to delete.
 */
void gkyl_correct_maxwellian_gyrokinetic_release(gkyl_correct_maxwellian_gyrokinetic* cmax);
