#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Limiters
enum gkyl_wave_limiter {
  GKYL_NO_LIMITER,
  GKYL_MIN_MOD,
  GKYL_SUPERBEE,
  GKYL_VAN_LEER,
  GKYL_MONOTONIZED_CENTERED,
  GKYL_BEAM_WARMING, GKYL_ZERO
};

struct gkyl_wave_prop_advance_status {
    int success; // 1 if step worked, 0 otherwise
    double dt_suggested; // suggested time-step
};

// Object type
typedef struct gkyl_wave_prop gkyl_wave_prop;

struct gkyl_wave_prop_inp {
    const struct gkyl_rect_grid *grid; // grid on which to solve equations
    const struct gkyl_wv_eqn *equation; // equation solver
    enum gkyl_wave_limiter limiter; // limiter to use
    int num_up_dirs; // number of update directions
    int update_dirs[GKYL_MAX_DIM]; // directions to update
    double cfl, cflm; // CFL and maximum CFL number
};

/**
 * Create new updater to update equations using wave-propagation
 * algorithm.
 *
 * @param winp Input for creating updater. See gkyl_wave_prop_inp above.
 */
gkyl_wave_prop* gkyl_wave_prop_new(struct gkyl_wave_prop_inp winp);

/**
 * Compute wave-propagation update. The update_rng MUST be a sub-range
 * of the range on which the array is defined. That is, it must be
 * either the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param wv Updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
struct gkyl_wave_prop_advance_status gkyl_wave_prop_advance(const gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *qin, struct gkyl_array *qout);
  
/**
 * Delete updater.
 *
 * @param wv Updater to delete.
 */
void gkyl_wave_prop_release(gkyl_wave_prop* wv);
