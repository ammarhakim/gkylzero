#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <gkyl_evalf_def.h>

// Forward declare for use in function pointers
struct gkyl_wv_eqn;

// Function pointer for compute waves from RP solver
typedef double (*wv_waves_t)(const struct gkyl_wv_eqn *eqn, 
  const double *delta, const double *ql, const double *qr, double *waves, double *speeds);

// Function pointer for compute q-fluctuations from waves
typedef void (*wv_qfluct_t)(const struct gkyl_wv_eqn *eqn, 
  const double *ql, const double *qr, const double *waves, const double *speeds,
  double *amdq, double *apdq);

// Function pointer to compute maximum speed given local state
typedef double (*wv_max_speed_t)(const struct gkyl_wv_eqn *eqn, const double *q);

// Function pointer to rotate conserved variables to local
// tangent-normal frame: tau1 X tau2 = norm
typedef void (*wv_rotate_to_local)(const double *tau1, const double *tau2, const double *norm,
  const double *qglobal, double *qlocal);

// Function pointer to rotate conserved variables to local
// tangent-normal frame: tau1 X tau2 = norm
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
