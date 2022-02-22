#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

/**
 * Create a new Gyrokinetic equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @return Pointer to Gyrokinetic equation object
 */
struct gkyl_dg_eqn* gkyl_dg_gyrokinetic_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Create a new Gyrokinetic equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing EM field
 * @param field_id enum to determine what type of EM fields (Gyrokinetic-Maxwell vs. neutrals)
 * @return Pointer to Gyrokinetic equation object
 */
struct gkyl_dg_eqn* gkyl_dg_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Set the geometry-related fields needed in computing gyrokinetic updates.
 *
 * @param eqn Equation pointer.
 * @param bmag Pointer to magnetic field magnitude.
 * @param jacobtot_inv Pointer to 1/(total jacobian field).
 * @param cmag Pointer to field with parallel gradient coefficient.
 * @param b_i Pointer to field with covariant components of B-field unit vector.
 */
void gkyl_gyrokinetic_set_geo_fields(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *bmag,
  const struct gkyl_array *jacobtot_inv, const struct gkyl_array *cmag, const struct gkyl_array *b_i);

/**
 * Set the electromanetic field pointers needed in computing gyrokinetic updates.
 *
 * @param eqn Equation pointer.
 * @param phi Pointer to electrostatic potential field.
 * @param apar Pointer to parallel vector potential field.
 * @param apardot Pointer to d(apar)/dt field.
 */
void gkyl_gyrokinetic_set_em_fields(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *phi,
  const struct gkyl_array *apar, const struct gkyl_array *apardot);

#ifdef GKYL_HAVE_CUDA
/**
 * CUDA device functions to set the geometry-related fields needed in computing gyrokinetic updates.
 *
 * @param eqn Equation pointer
 * @param bmag Pointer to magnetic field magnitude field.
 * @param jacobtot_inv Pointer to 1/(total jacobian field).
 * @param cmag Pointer to field with parallel gradient coefficient.
 * @param b_i Pointer to field with covariant components of B-field unit vector.
 */
void gkyl_gyrokinetic_set_geo_fields_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *bmag,
  const struct gkyl_array *jacobtot_inv, const struct gkyl_array *cmag, const struct gkyl_array *b_i);

/**
 * Set the electromanetic field pointers needed in computing gyrokinetic updates.
 *
 * @param eqn Equation pointer.
 * @param phi Pointer to electrostatic potential field.
 * @param apar Pointer to parallel vector potential field.
 * @param apardot Pointer to d(apar)/dt field.
 */
void gkyl_gyrokinetic_set_em_fields_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *phi,
  const struct gkyl_array *apar, const struct gkyl_array *apardot);
#endif

