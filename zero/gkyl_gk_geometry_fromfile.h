#pragma once

#include <gkyl_gk_geometry.h>
/**
 * Create a new geometry object by reading geometric quantities from file. 
 *
 * @param grid Grid on which geometry lives
 * @param range Range on which geometry should be constructed
 * @param range_ext 
 * @param basis configuration space basis
 */
struct gk_geometry* gkyl_gk_geometry_fromfile_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, const struct gkyl_basis* basis, bool use_gpu);

/**
 * Create a new gk geometry object that lives on NV-GPU: see new() method
 * above for documentation.
 */

struct gk_geometry* gkyl_gk_geometry_fromfile_cu_dev_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, const struct gkyl_basis* basis);


