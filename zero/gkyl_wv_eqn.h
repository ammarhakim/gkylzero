#pragma once

#include <gkyl_ref_count.h>

// Forward declare for use in function pointers
struct gkyl_wv_eqn;

// Function pointer for compute waves from RP solver
typedef double (*wv_waves_t)(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *ql, const double *qr, double *waves, double *speeds);

struct gkyl_wv_eqn {
    int num_equations; // number of equations in system
    int num_waves; // number of waves in system
    wv_waves_t wave_func; // function to compute waves and speeds
    struct gkyl_ref_count ref_count; // reference count
};

/**
 * Aquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 */
struct gkyl_wv_eqn* gkyl_wv_eqn_aquire(const struct gkyl_wv_eqn* eqn);

/**
 * Compute waves and speeds from left/right conserved variables. The
 * 'waves' array has size num_equations X num_waves in length. The 'm'
 * wave (m = 0 ... num_waves-1) is stored starting at location
 * waves[m*num_equations].
 *
 * @param eqn Equation object
 * @param dir Direction to compute waves/speeds
 * @param ql Conserved variables on left of interface
 * @param qr Conserved variables on right of interface
 * @param waves On output, waves 
 * @param speed On output wave speeds[num_wave]
 * @return Maximum wave speed.
 */
double gkyl_wv_waves(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *ql, const double *qr, double *waves, double *speeds);

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn);
