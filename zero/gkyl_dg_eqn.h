#pragma once

#include <gkyl_real_type.h>
#include <gkyl_ref_count.h>

// Forward declare for use in function pointers
struct gkyl_dg_eqn;

// Function pointer type for volume kernel
typedef gkyl_real (*vol_termf_t)(const struct gkyl_dg_eqn *eqn, 
  const gkyl_real*  xc, const gkyl_real*  dx, const int*  idx,
  const gkyl_real* qIn, gkyl_real* restrict qRhsOut);

// Function pointer type for surface kernel
typedef gkyl_real (*surf_termf_t)(const struct gkyl_dg_eqn *eqn, int dir,
  const gkyl_real*  xcL, const gkyl_real*  xcR, const gkyl_real*  dxL, const gkyl_real* dxR,
  gkyl_real maxsOld, const int*  idxL, const int*  idxR,
  const gkyl_real* qInL, const gkyl_real*  qInR, gkyl_real* restrict qRhsOutL, gkyl_real* restrict qRhsOutR);

struct gkyl_dg_eqn {
    int num_equations; // number of equations in system
    vol_termf_t vol_term; // volume term kernel
    surf_termf_t surf_term; // surface term kernel
    surf_termf_t boundary_surf_term; // boundary surface term kernel
    struct gkyl_ref_count ref_count; // reference count     
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
