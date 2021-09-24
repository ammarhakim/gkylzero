#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

/**
 * Create a new LBO equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @return Pointer to LBO equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_lbo_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Create a new LBO equation object that lives on NV-GPU
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @return Pointer to LBO equation object
 */
struct gkyl_dg_eqn* gkyl_dg_vlasov_lbo_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Set the nu needed in updating the lbo terms.
 * 
 * @param eqn Equation pointer
 * @param nuSum Value of collision frequency in the cell.
 */
void gkyl_vlasov_lbo_set_nuSum(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuSum);

/**
 * Set the nu*u needed in updating the drag flux term.
 * 
 * @param eqn Equation pointer
 * @param nuUSum Pointer to nu multiplied by drift.
 */
void gkyl_vlasov_lbo_set_nuUSum(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuUSum);

/**
 * Set the nu*vt^2 needed in updating the diffusion flux term.
 * 
 * @param eqn Equation pointer
 * @param nuVtSqSum Pointer to nu multiplied by thermal velocity.
 */
void gkyl_vlasov_lbo_set_nuVtSqSum(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuVtSqSum);
