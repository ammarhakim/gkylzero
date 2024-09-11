#pragma once

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_array.h>

// Object type.
typedef struct gkyl_translate_dim_gyrokinetic gkyl_translate_dim_gyrokinetic;

/**
 * Create a new updater that translates the DG coefficients of a lower
 * dimensional donor field to those of a higher dimensional target
 * field (for gyrokinetics).
 *
 * For now this is meant for translating:
 *   x,vpar,mu -> x,y,vpar,mu
 *   x,vpar,mu -> x,y,z,vpar,mu
 *   x,y,vpar,mu -> x,y,z,vpar,mu
 *
 * @param cdim_do Configuration space dimension of the donor field.
 * @param pbasis_do Phase basis of the donor field.
 * @param cdim_tar Configuration space dimension of the target field.
 * @param pbasis_tar Phase basis of the target field.
 * @param use_gpu Whether to run it on the GPU or not.
 */
struct gkyl_translate_dim_gyrokinetic*
gkyl_translate_dim_gyrokinetic_new(int cdim_do, struct gkyl_basis pbasis_do,
  int cdim_tar, struct gkyl_basis pbasis_tar, bool use_gpu);

/**
 * Run the updater that translates the DG coefficients of a lower dimensional
 * donor field to those of a higher dimensional target field.
 *
 * @param up Updater object. 
 * @param phase_rng_do Phase range of the donor field.
 * @param phase_rng_tar Phase range of the target field.
 * @param fdo Donor field.
 * @param ftar target field.
 */
void
gkyl_translate_dim_gyrokinetic_advance(gkyl_translate_dim_gyrokinetic* up,
  const struct gkyl_range *phase_rng_do, const struct gkyl_range *phase_rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar);

/**
 * Release the memory associated with the translate_dim updater. 
 *
 * @param up translate_dim updater.
 */
void
gkyl_translate_dim_gyrokinetic_release(gkyl_translate_dim_gyrokinetic* up);
