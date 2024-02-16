#pragma once

#include <gkyl_gk_geometry.h>
/**
 * Create a new geometry object using mirror input (efit)
 *
 * @param grid Grid on which geometry lives
 * @param local, _ext local config-space range and extended range on which geometry should be constructed
 * @param global, _ext global config-space range and extended range.
 * @param basis configuration space basis
 * @param mirror_rz_ctx RZ Context for use in mapping (efit file and associated quantities)
 * @param mirror_comp_ctx computation domain context for calculating mapping
 */
struct gk_geometry* gkyl_gk_geometry_mirror_new(const struct gkyl_rect_grid* grid, 
    const struct gkyl_range *local, const struct gkyl_range* local_ext, 
    const struct gkyl_range *global, const struct gkyl_range* global_ext,
    const struct gkyl_basis* basis, void* mirror_rz_ctx, void* mirror_comp_ctx, bool use_gpu);

/**
 * Create a new gk geometry object that lives on NV-GPU: see new() method
 * above for documentation.
 */

struct gk_geometry* gkyl_gk_geometry_mirror_cu_dev_new(const struct gkyl_rect_grid* grid, 
    const struct gkyl_range *local, const struct gkyl_range* local_ext, 
    const struct gkyl_range *global, const struct gkyl_range* global_ext,
    const struct gkyl_basis* basis, void* mirror_rz_ctx, void* mirror_comp_ctx);


