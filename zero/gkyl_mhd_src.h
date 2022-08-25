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
  double glm_alpha; // Mignone & Tzeferacos, JCP (2010) 229, 2117, Eq (27).
  double dxyz_min;
};

// Object type
struct gkyl_mhd_src {
  struct gkyl_rect_grid grid;
  int ndim;
  enum gkyl_wv_mhd_div_constraint divergence_constraint;
  double glm_ch;
  double glm_alpha; // Mignone & Tzeferacos, JCP (2010) 229, 2117, Eq (27).
  double dxyz_min;
};
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
