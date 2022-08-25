#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_wv_mhd.h>

struct gkyl_mhd_src_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  enum gkyl_wv_mhd_div_constraint divergence_constraint;
  double glm_ch;
  double glm_cp;
  double dxyz_min;
  double cfl; // cfl used in the hyperbolic part
};

// Object type
typedef struct gkyl_mhd_src gkyl_mhd_src;

/**
 * Create new updater.
 *
 * @param inp Input parameters to updater
 */
gkyl_mhd_src *gkyl_mhd_src_new(struct gkyl_mhd_src_inp inp);

/**
 * @param mes
 * @param dt
 */
void gkyl_mhd_src_advance(const gkyl_mhd_src *mes,
                          double dt,
                          const struct gkyl_range *update_rng,
                          struct gkyl_array *q,
                          const struct gkyl_array *app_accel);

/**
 * Delete updater.
 *
 * @param mes Updater to delete.
 */
void gkyl_mhd_src_release(gkyl_mhd_src *mes);
