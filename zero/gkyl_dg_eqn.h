#pragma once

#include <gkyl_ref_count.h>

// Forward declare for use in function pointers
struct gkyl_dg_eqn;

// Function pointer type for volume kernel
typedef double (*vol_termf_t)(const struct gkyl_dg_eqn *eqn, 
  const double*  xc, const double*  dx, const int*  idx,
  const double* qIn, double* GKYL_RESTRICT qRhsOut);

// Function pointer type for surface kernel
typedef double (*surf_termf_t)(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  double maxsOld, const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut);

// Function pointer type for surface kernel
typedef double (*boundary_surf_termf_t)(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  double maxsOld, const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut);

struct gkyl_dg_eqn {
  int num_equations; // number of equations in system
  vol_termf_t vol_term; // volume term kernel
  surf_termf_t surf_term; // surface term kernel
  boundary_surf_termf_t boundary_surf_term; // boundary surface term kernel
  struct gkyl_ref_count ref_count; // reference count     
};

/**
 * Aquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 */
static inline
struct gkyl_dg_eqn* gkyl_dg_eqn_aquire(const struct gkyl_dg_eqn* eqn) 
{
  gkyl_ref_count_inc(&eqn->ref_count);
  return (struct gkyl_dg_eqn*) eqn;
}

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
static inline
void gkyl_dg_eqn_release(const struct gkyl_dg_eqn* eqn)
{
  gkyl_ref_count_dec(&eqn->ref_count);
}
