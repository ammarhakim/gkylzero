#pragma once

#include <gkyl_ref_count.h>

// Forward declare for use in function pointers
struct gkyl_dg_eqn;

// Function pointer type for volume kernel
typedef double (*vol_termf_t)(const struct gkyl_dg_eqn *eqn, 
  const double*  xc, const double*  dx, const int*  idx,
  const double* qIn, double *qRhsOut);


// Function pointer type for surface kernel
typedef double (*surf_termf_t)(const struct gkyl_dg_eqn *eqn, int dir,
  const double*  xcL, const double*  xcR, const double*  dxL, const double* dxR,
  double maxsOld, const int*  idxL, const int*  idxR,
  const double* qInL, const double*  qInR, double *qRhsOutL, double *qRhsOutR);

struct gkyl_dg_eqn {
    int num_equations; // Number of equations in system
    vol_termf_t vol_term; // Volume term kernel
    surf_termf_t surf_term; // Surface term kernel
    surf_termf_t boundary_surf_term; // Boundary surface term kernel
    struct gkyl_ref_count ref_count; // Reference count     
};

/**
 * Aquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 */
struct gkyl_dg_eqn* gkyl_dg_eqn_aquire(const struct gkyl_dg_eqn* eqn);

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_dg_eqn_release(const struct gkyl_dg_eqn* eqn);
