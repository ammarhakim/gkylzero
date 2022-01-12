#pragma once

#include <gkyl_mat.h>
#include <gkyl_ref_count.h>

// Forward declare for use in function pointers
struct gkyl_prim_lbo;

// Self-primitive moment kernel pointer type
typedef void (*self_prim_t)(const struct gkyl_prim_lbo *prim, struct gkyl_mat *A,
  struct gkyl_mat *rhs, const double *m0, const double *m1, const double *m2,
  const double *cM, const double *cE);

struct gkyl_prim_lbo {
  int cdim; // config-space dim
  int pdim; // phase-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_phase; // number of basis functions in phase-space
  self_prim_t self_prim; // moment calculation kernel
  struct gkyl_ref_count ref_count; // reference count     
};

/**
 * Acquire primitive moment object pointer. Delete using the release()
 * method
 *
 * @param prim Primitive moment object.
 */
struct gkyl_prim_lbo* gkyl_prim_lbo_acquire(const struct gkyl_prim_lbo* prim);

/**
 * Delete primitive moment object
 *
 * @param prim Primitive moment object to delete.
 */
void gkyl_prim_lbo_release(const struct gkyl_prim_lbo* prim);
