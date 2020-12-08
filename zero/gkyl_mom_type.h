#pragma once

#include <gkyl_ref_count.h>

/**
 * Function pointer type to compute the needed moment.
 */
typedef void (*momf_t)(const double *xc, const double *dx,
  const long *idx, const double *f, double* out);

struct gkyl_mom_type {
    int cdim; // config-space dim
    int pdim; // phase-space dim
    int polyOrder; // polynomal order
    int num_config; // number of basis functions in config-space
    int num_phase; // number of basis functions in phase-space
    int num_mom; // number of components in moment
    momf_t kernel; // pointer to moment calculation kernel
    struct gkyl_ref_count ref_count; // reference count    
};

/**
 * Aquire pointer to moment object. Delete using the release() method
 *
 * @param momt Moment object to delete.
 */
struct gkyl_mom_type* gkyl_mom_type_aquire(const struct gkyl_mom_type* momt);

/**
 * Delete moment object
 *
 * @param momt Moment object to delete.
 */
void gkyl_mom_type_release(struct gkyl_mom_type* momt);

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
  const double *xc, const double *dx, const long *idx,
  const double *f, double* out);
