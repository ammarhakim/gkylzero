#pragma once

#include <gkyl_gk_geometry.h>
/**
 * Create a new geometry object using mirror input (efit)
 *
 * @param grid Grid on which geometry lives
 * @param range Range on which geometry should be constructed
 * @param range_ext 
 * @param basis configuration space basis
 * @param mirror_rz_ctx RZ Context for use in mapping (efit file and associated quantities)
 * @param mirror_comp_ctx computation domain context for calculating mapping
 */
struct gk_geometry* gkyl_gk_geometry_mirror_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, void* mirror_rz_ctx, void* mirror_comp_ctx, bool use_gpu);

/**
 * Create a new gk geometry object that lives on NV-GPU: see new() method
 * above for documentation.
 */

struct gk_geometry* gkyl_gk_geometry_mirror_cu_dev_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, void* mirror_rz_ctx, void* mirror_comp_ctx);

