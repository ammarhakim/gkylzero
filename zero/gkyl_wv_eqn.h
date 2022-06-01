#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <gkyl_evalf_def.h>

// Forward declare for use in function pointers
struct gkyl_wv_eqn;

/**
 * Function pointer for computing waves and speeds from left/right conserved variables. The
 * 'waves' array has size num_equations X num_waves in length. The 'm'
 * wave (m = 0 ... num_waves-1) is stored starting at location
 * waves[m*num_equations].
 *
 * @param eqn Equation object
 * @param delta Jump across interface to split
 * @param ql Conserved variables on left of interface
 * @param qr Conserved variables on right of interface
 * @param waves On output, waves 
 * @param speeds On output wave speeds[num_wave]
 * @return Maximum wave speed.
 */
typedef double (*wv_waves_t)(const struct gkyl_wv_eqn *eqn, 
  const double *delta, const double *ql, const double *qr, double *waves, double *speeds);

// Function pointer for computing q-fluctuations from waves
typedef void (*wv_qfluct_t)(const struct gkyl_wv_eqn *eqn, 
  const double *ql, const double *qr, const double *waves, const double *speeds,
  double *amdq, double *apdq);

// Function pointer for computing maximum speed given local state
typedef double (*wv_max_speed_t)(const struct gkyl_wv_eqn *eqn, const double *q);

/**
 * Function pointer to rotate state (conserved/primitive) vector to local tangent-normal coordinate frame.
 *
 * @param eqn Equation object
 * @param tau1 Tangent vector
 * @param tau2 Tangent vector
 * @param norm Normal vector such that norm = tau1 x tau2
 * @param qglobal State vector in global coordinates
 * @param qlocal State vector in local coordinates
 */
typedef void (*wv_rotate_to_local)(const double *tau1, const double *tau2, const double *norm,
  const double *qglobal, double *qlocal);

/**
 * Function pointer to rotate state (conserved/primitive) vector to local global coordinate frame.
 *
 * @param eqn Equation object
 * @param tau1 Tangent vector
 * @param tau2 Tangent vector
 * @param norm Normal vector such that norm = tau1 x tau2
 * @param qglobal State vector in global coordinates
 * @param qlocal State vector in local coordinates
 */
typedef void (*wv_rotate_to_global)(const double *tau1, const double *tau2, const double *norm,
  const double *qlocal, double *qglobal);

struct gkyl_wv_eqn {
  enum gkyl_eqn_type type; // Equation type
  int num_equations; // number of equations in system
  int num_waves; // number of waves in system
  wv_waves_t waves_func; // function to compute waves and speeds
  wv_qfluct_t qfluct_func; // function to compute q-fluctuations
  wv_max_speed_t max_speed_func; // function to compute max-speed
  wv_rotate_to_local rotate_to_local_func; // function to rotate to local frame
  wv_rotate_to_global rotate_to_global_func; // function to rotate to global frame

  wv_bc_func_t wall_bc_func; // function to apply wall BC

  uint32_t flags;  
  struct gkyl_ref_count ref_count; // reference count
  struct gkyl_wv_eqn *on_dev; // pointer to itself or device data
};

/**
 * Acquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_eqn_acquire(const struct gkyl_wv_eqn* eqn);

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn);
