#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_lbo_vlasov_pkpm_drag_auxfields {
  const struct gkyl_array *nu;
};

/**
 * Create a new LBO drag term equation object for Vlasov equation in 
 * parallel-kinetic-perpendicular-moment (pkpm) model.
 * Note that pkpm model is in local rest frame so drag term is only div_v(v f).
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @param use_gpu Bool to determine if equation object is on host or device
 * @return Pointer to LBO drag term equation object
 */
struct gkyl_dg_eqn* gkyl_dg_lbo_vlasov_pkpm_drag_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, bool use_gpu);

struct gkyl_dg_eqn* gkyl_dg_lbo_vlasov_pkpm_drag_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range);

/**
 * Set auxiliary fields needed in updating the drag flux term (nu = the collision frequency).
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_vlasov_pkpm_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_vlasov_pkpm_drag_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields needed in updating the drag flux term  (nu = the collision frequency).
 *
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_lbo_vlasov_pkpm_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_vlasov_pkpm_drag_auxfields auxin);

#endif
