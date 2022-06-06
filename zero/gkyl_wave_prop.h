#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Limiters
enum gkyl_wave_limiter {
  GKYL_NO_LIMITER = 1, // to allow default to be 0
  GKYL_MIN_MOD,
  GKYL_SUPERBEE,
  GKYL_VAN_LEER,
  GKYL_MONOTONIZED_CENTERED,
  GKYL_BEAM_WARMING,
  GKYL_ZERO
};

struct gkyl_wave_prop_status {
  int success; // 1 if step worked, 0 otherwise
  double dt_suggested; // suggested time-step
};

// Object type for updater
struct gkyl_wave_prop {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  enum gkyl_wave_limiter limiter; // limiter to use
  double cfl; // CFL number
  const struct gkyl_wv_eqn *equation; // equation object
  struct gkyl_wave_geom *geom; // geometry object
  struct gkyl_array *waves, *speeds, *flux2; // data for 1D slice update
  struct gkyl_wave_prop *on_dev; // to itself or device copy
};
typedef struct gkyl_wave_prop gkyl_wave_prop;

// Parameters for constructor
struct gkyl_wave_prop_inp {
  const struct gkyl_rect_grid *grid; // grid on which to solve equations
  const struct gkyl_wv_eqn *equation; // equation solver
  enum gkyl_wave_limiter limiter; // limiter to use
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  double cfl; // CFL number to use

  const struct gkyl_wave_geom *geom; // geometry
};

/**
 * Create new updater to update equations using wave-propagation
 * algorithm.
 *
 * @param winp Input for creating updater. See gkyl_wave_prop_inp above.
 */
gkyl_wave_prop* gkyl_wave_prop_new(struct gkyl_wave_prop_inp winp);

/**
 * Create new updater on CUDA device to update equations using wave-propagation
 * algorithm.
 *
 * @param winp Input for creating updater. See gkyl_wave_prop_inp above.
 */
gkyl_wave_prop* gkyl_wave_prop_cu_dev_new(struct gkyl_wave_prop_inp winp);


/**
 * Compute wave-propagation update. The update_rng MUST be a sub-range
 * of the range on which the array is defined. That is, it must be
 * either the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param wv Updater object
 * @param tm Current time
 * @param dt time-step
 * @param update_rng Range on which to compute.
 * @param qin Input to updater
 * @param qout Solution at tm+dt
 */
struct gkyl_wave_prop_status gkyl_wave_prop_advance(const gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *qin, struct gkyl_array *qout);

/**
 * Compute an estimate of maximum stable time-step for given input
 * state 'qin'
 *
 * @param wv Updater object
 * @param qin Input to compute dt for
 * @return maximum stable time-step
 */
double gkyl_wave_prop_max_dt(const gkyl_wave_prop *wv, const struct gkyl_range *update_range,
  const struct gkyl_array *qin);
  
/**
 * Delete updater.
 *
 * @param wv Updater to delete.
 */
void gkyl_wave_prop_release(gkyl_wave_prop* wv);

/**
 * Compute wave-propagation update on GPU. The update_rng MUST be a sub-range
 * of the range on which the array is defined. That is, it must be
 * either the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param wv Updater object
 * @param tm Current time
 * @param dt time-step
 * @param update_rng Range on which to compute.
 * @param qin Input to updater
 * @param qout Solution at tm+dt
 */
struct gkyl_wave_prop_status gkyl_wave_prop_cu_dev_advance(const gkyl_wave_prop *wv,
  double tm, double dt, const struct gkyl_range *update_range,
  const struct gkyl_array *qin, struct gkyl_array *qout);

