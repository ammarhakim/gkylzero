#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>

// Object type
typedef struct gkyl_vtsq gkyl_vtsq;

/**
 * Create a new updater that computes the thermal speed squared given the
 * velocity moments (M_0, M_1, M_2) using:
 *   v_t^2 = (M_2 - (M_1/M_0) . M_1)/(d_v*M_0*M_0)
 * where d_v is the number of physical velocity dimensions accounted for.
 *
 * @param basis Basis object (configuration space).
 * @param range Range in which we'll compute m0_s*delta_s.
 * @param m1comps Number of vector components in M_1.
 * @param vdim_phys Number of physical velocity dimensions accounted for.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_vtsq* gkyl_vtsq_new(
  const struct gkyl_basis *basis, const struct gkyl_range *range,
  int m1comps, int vdim_phys, bool use_gpu);

/**
 * Compute m0_s*delta_s.
 *
 * @param up Struct defining this updater..
 * @param basis Basis object (configuration space).
 * @param moms Velocity moments of the distribution function.
 * @param range Range in which we'll compute m0_s*delta_s.
 * @param out Output array.
 */
void gkyl_vtsq_advance(struct gkyl_vtsq *up, struct gkyl_basis basis,
  const struct gkyl_array *moms, const struct gkyl_range *range,
  struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_vtsq_release(struct gkyl_vtsq *up);
