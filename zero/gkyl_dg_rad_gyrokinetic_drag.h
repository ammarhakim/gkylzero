#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Struct containing the pointers to auxiliary fields.
// Okay if these are not const?
struct gkyl_dg_rad_gyrokinetic_drag_auxfields { 
  const struct gkyl_array *nI;
  const struct gkyl_array *vnu;
  const struct gkyl_array *vsqnu;
  const struct gkyl_array *bmag;
  const struct gkyl_array *fit_params;
};

/**
 * Create a new gyrokinetic RAD drag term equation object.
 *
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis functions
 * @param conf_range Configuration space range for use in indexing primitive moments
 * @param pgrid Phase-space grid object.
 * @param bmag Magnetic field
 * @param fit_params Fit parameters in nu
 * @return Pointer to RAD equation object
 */
struct gkyl_dg_eqn* gkyl_dg_rad_gyrokinetic_drag_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
						     const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid, const struct gkyl_array *bmag, const struct gkyl_array *fit_params, bool use_gpu);

/**
 * TO DO: Create a new RAD equation object that lives on NV-GPU
 */
struct gkyl_dg_eqn* gkyl_dg_rad_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid, double mass);

/**
 * Set auxiliary fields needed in updating the drag flux term.
 * These are bmag, vnu, vsqnu, and nI
 * 
 * @param eqn Equation pointer
 * @param auxfields Pointer to struct of aux fields.
 */
void gkyl_rad_gyrokinetic_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_rad_gyrokinetic_drag_auxfields auxin);

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set auxiliary fields needed in updating the drag flux term.
 */
void gkyl_rad_gyrokinetic_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_rad_gyrokinetic_drag_auxfields auxin);

#endif
