#pragma once

#include <gkyl_ref_count.h>
#include <gkyl_util.h>

#include <stdbool.h>

// Forward declare for use in function pointers
struct gkyl_dg_eqn;

// Function pointer type for volume kernel
typedef double (*vol_termf_t)(const struct gkyl_dg_eqn *eqn, 
  const double*  xc, const double*  dx, const int*  idx,
  const double* qIn, double* GKYL_RESTRICT qRhsOut);

// Function pointer type for surface kernel
typedef void (*surf_termf_t)(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut);

// Function pointer type for surface kernel
typedef void (*boundary_surf_termf_t)(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut);

struct gkyl_dg_eqn {
  int num_equations; // number of equations in system
  vol_termf_t vol_term; // volume term kernel
  surf_termf_t surf_term; // surface term kernel
  boundary_surf_termf_t boundary_surf_term; // boundary surface term kernel

  uint32_t flags;
  struct gkyl_ref_count ref_count; // reference count
  struct  gkyl_dg_eqn *on_dev; // pointer to itself or device data
};

/**
 * Check if equation is on device.
 *
 * @param eqn Equation to check
 * @return true if eqn on device, false otherwise
 */
bool gkyl_dg_eqn_is_cu_dev(const struct gkyl_dg_eqn *eqn);

/**
 * Acquire pointer to equation object. Delete using the release()
 * method
 *
 * @param eqn Equation object.
 */
struct gkyl_dg_eqn* gkyl_dg_eqn_acquire(const struct gkyl_dg_eqn* eqn);

/**
 * Delete equation object
 *
 * @param eqn Equation object to delete.
 */
void gkyl_dg_eqn_release(const struct gkyl_dg_eqn *eqn);
