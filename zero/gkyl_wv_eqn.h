#pragma once

#include <gkyl_ref_count.h>

// Forward declare for use in function pointers
struct gkyl_wv_eqn;

// Function pointer for compute waves from RP solver
typedef double (*wv_rp_t)(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *ql, const double *qr, double *waves, double *speeds);

struct gkyl_wv_eqn {
    int num_equations; // number of equations in system
    int num_wave; // number of waves in system
    wv_rp_t rp; // Riemann solver
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
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_wv_eqn_release(const struct gkyl_wv_eqn* eqn);
