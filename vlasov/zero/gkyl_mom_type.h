#pragma once

#include <gkyl_ref_count.h>
#include <stdint.h>

// Forward declare for use in function pointers
struct gkyl_mom_type;

/**
 * Function pointer type to compute the needed moment.
 */
typedef void (*momf_t)(const struct gkyl_mom_type *momt,
  const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param);

struct gkyl_mom_type {
  int cdim; // config-space dim
  int pdim; // phase-space dim
  int poly_order; // polynomal order
  int num_config; // number of basis functions in config-space
  int num_phase; // number of basis functions in phase-space
  int num_mom; // number of components in moment
  momf_t kernel; // moment calculation kernel
  struct gkyl_ref_count ref_count; // reference count

  uint32_t flags;
  struct gkyl_mom_type *on_dev; // pointer to itself or device data
};

/**
 * Check if moment type is on device.
 *
 * @param momt Moment type to check
 * @return true if momt on device, false otherwise
 */
bool gkyl_mom_type_is_cu_dev(const struct gkyl_mom_type *momt);

/**
 * Acquire pointer to moment object. Delete using the release() method
 *
 * @param momt Moment object to get pointer from.
 * @return acquired object
 */
struct gkyl_mom_type* gkyl_mom_type_acquire(const struct gkyl_mom_type* momt);

/**
 * Delete moment object
 *
 * @param momt Moment object to delete.
 */
void gkyl_mom_type_release(const struct gkyl_mom_type* momt);

/**
 * Calculate moment specified by mom_type object.
 *
 * @param momt Moment type object
 * @param xc Cell center coordinates
 * @param dx Cell size in each direction
 * @param idx Index into phase-space cell
 * @param f Input pointer to distribution function in cell
 * @param out On output, contribution to moment from phase-space cell
 */
void gkyl_mom_type_calc(const struct gkyl_mom_type* momt,
  const double *xc, const double *dx, const int *idx,
  const double *f, double* GKYL_RESTRICT out, void *param);

/**
 * Get number of moments specified by mom_type object
 *
 * @param momt Moment type object
 * returns int Number of moments
 */
int gkyl_mom_type_num_mom(const struct gkyl_mom_type* momt);
