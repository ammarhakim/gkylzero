#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_lbo_vlasov_drag_auxfields {
  const struct gkyl_array *nuSum;
  const struct gkyl_array *nuUSum;
  const struct gkyl_array *nuVtSqSum;
};

/**
 * Create a new LBO equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @return Pointer to LBO equation object
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_vlasov_drag_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Create a new LBO equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @return Pointer to LBO equation object
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_vlasov_drag_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Set auxiliary fields needed in updating the drag flux term.
 * These are nu, nu*u, and nu*vt^2.
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_vlasov_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_vlasov_drag_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * Set auxiliary fields needed in updating the drag flux term.
 * These are nu, nu*u, and nu*vt^2.
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_vlasov_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_vlasov_drag_auxfields auxin);

#endif
