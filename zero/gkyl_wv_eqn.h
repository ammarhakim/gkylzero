#pragma once

#include <gkyl_eqn_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_util.h>
#include <gkyl_evalf_def.h>

// Flux type for use in wave/qfluct methods
enum gkyl_wv_flux_type { GKYL_WV_HIGH_ORDER_FLUX, GKYL_WV_LOW_ORDER_FLUX };

// Forward declare for use in function pointers
struct gkyl_wv_eqn;

// Function pointer for compute waves from RP solver
typedef double (*wv_waves_t)(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *speeds);

// Function pointer for compute q-fluctuations from waves
typedef void (*wv_qfluct_t)(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *speeds,
  double *amdq, double *apdq);

// Function pointer to check if invariant domain is preserved
typedef bool (*wv_check_inv)(const struct gkyl_wv_eqn *eqn, const double *q);

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

// Function pointer to convert conserved variables to Riemann
// variables given an input state 'qstate'
typedef void (*wv_cons_to_riem)(const struct gkyl_wv_eqn *eqn,
  const double *qstate, const double *qin, double *wout);

struct gkyl_wv_eqn {
  enum gkyl_eqn_type type; // Equation type
  int num_equations; // number of equations in system
  int num_waves; // number of waves in system
  wv_waves_t waves_func; // function to compute waves and speeds
  wv_qfluct_t qfluct_func; // function to compute q-fluctuations
  wv_check_inv check_inv_func; // function to check invariant domains
  wv_max_speed_t max_speed_func; // function to compute max-speed
  wv_rotate_to_local rotate_to_local_func; // function to rotate to local frame
  wv_rotate_to_global rotate_to_global_func; // function to rotate to global frame

  wv_cons_to_riem cons_to_riem; // function to convert cons to Riemann vars

  wv_bc_func_t wall_bc_func; // function to apply wall BC
  wv_bc_func_t no_slip_bc_func; // function to apply no-slip BC
  
  struct gkyl_ref_count ref_count; // reference count
};

/**
 * Acquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_eqn_acquire(const struct gkyl_wv_eqn* eqn);

/**
 * Compute waves and speeds from left/right conserved variables. The
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
inline double
gkyl_wv_eqn_waves(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, double *waves, double *speeds)
{
  return eqn->waves_func(eqn, type, delta, ql, qr, waves, speeds);
}

/**
 * Compute waves and speeds from left/right conserved variables. The
 * 'waves' array has size num_equations X num_waves in length. The 'm'
 * wave (m = 0 ... num_waves-1) is stored starting at location
 * waves[m*num_equations].
 *
 * @param eqn Equation object
 * @param ql Conserved variables on left of interface
 * @param qr Conserved variables on right of interface
 * @param waves Waves computed from waves() method
 * @param speeds Wave speeds[num_wave]
 * @param amdq On output, the left-going fluctuations.
 * @param apdq On output, the right-going fluctuations.
 */
inline void
gkyl_wv_eqn_qfluct(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double *waves, const double *speeds,
  double *amdq, double *apdq)
{
  eqn->qfluct_func(eqn, type, ql, qr, waves, speeds, amdq, apdq);  
}

/**
 * Compute waves and speeds from left/right conserved variables. The
 * 'waves' array has size num_equations X num_waves in length. The 'm'
 * wave (m = 0 ... num_waves-1) is stored starting at location
 * waves[m*num_equations].
 *
 * @param eqn Equation object
 * @param q Conserved variables
 * @return maximum wave-speed in direction 'dir'
 */
inline double
gkyl_wv_eqn_max_speed(const struct gkyl_wv_eqn *eqn, const double *q)
{
  return eqn->max_speed_func(eqn, q);
}

/**
 * Rotate state (conserved/primitive) vector to local tangent-normal coordinate frame.
 *
 * @param eqn Equation object
 * @param tau1 Tangent vector
 * @param tau2 Tangent vector
 * @param norm Normal vector such that norm = tau1 x tau2
 * @param qglobal State vector in global coordinates
 * @param qlocal State vector in local coordinates
 */
inline void
gkyl_wv_eqn_rotate_to_local(const struct gkyl_wv_eqn* eqn,
  const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal)
{
  return eqn->rotate_to_local_func(tau1, tau2, norm, qglobal, qlocal);
}

/**
 * Rotate state (conserved/primitive) vector to global coordinate frame.
 *
 * @param eqn Equation object
 * @param tau1 Tangent vector
 * @param tau2 Tangent vector
 * @param norm Normal vector such that norm = tau1 x tau2
 * @param qlocal State vector in local coordinates
 * @param qglobal State vector in local coordinates
 */
inline void
gkyl_wv_eqn_rotate_to_global(const struct gkyl_wv_eqn* eqn,
  const double *tau1, const double *tau2, const double *norm,
  const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal)
{
  return eqn->rotate_to_global_func(tau1, tau2, norm, qlocal, qglobal);
}

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn);
