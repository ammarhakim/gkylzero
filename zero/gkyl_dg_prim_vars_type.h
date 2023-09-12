#pragma once

#include <gkyl_ref_count.h>
#include <stdint.h>

// Forward declare for use in function pointers
struct gkyl_dg_prim_vars_type;

/**
 * Function pointer type to compute the needed primitive variables.
 */
typedef void (*pvf_t)(const struct gkyl_dg_prim_vars_type *pvt,
  const int *idx, const double *in, double* out);

struct gkyl_dg_prim_vars_type {
  int cdim; // config-space dim
  int vdim; // velocity-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_mom; // number of components in primitive variables
  pvf_t kernel; // primitive variable calculation kernel
  struct gkyl_ref_count ref_count; // reference count

  uint32_t flags;
  struct gkyl_dg_prim_vars_type *on_dev; // pointer to itself or device data
};

/**
 * Check if primitive variables type is on device.
 *
 * @param pvt Primitive variable type to check
 * @return true if pvt on device, false otherwise
 */
bool gkyl_dg_prim_vars_type_is_cu_dev(const struct gkyl_dg_prim_vars_type *pvt);

/**
 * Acquire pointer to primitive variables object. Delete using the release() method
 *
 * @param pvt Primitive variable object to get pointer from.
 * @return acquired object
 */
struct gkyl_dg_prim_vars_type* gkyl_dg_prim_vars_type_acquire(const struct gkyl_dg_prim_vars_type* pvt);

/**
 * Delete Primitive variable object
 *
 * @param pvt Primitive variable object to delete.
 */
void gkyl_dg_prim_vars_type_release(const struct gkyl_dg_prim_vars_type* pvt);


