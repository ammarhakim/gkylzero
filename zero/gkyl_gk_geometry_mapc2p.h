#pragma once

#include <gkyl_gk_geometry.h>
/**
 * Create a new geometry object using mac2p. 
 *
 * @param grid Grid on which geometry lives
 * @param range Range on which geometry should be constructed
 * @param range_ext 
 * @param basis configuration space basis
 * @param mapc2p Mapping from computational to physical space
 * @param mapc2p_ctx Context for use in mapping
 * @param bmag function which gives |B| in computational space
 * @param bmag_ctx Context for calculating |B|
 */
struct gk_geometry* gkyl_gk_geometry_mapc2p_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range,
    const struct gkyl_range* range_ext, const struct gkyl_range *global, const struct gkyl_range* global_ext,
    const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx, bool use_gpu);

/**
 * Create a new gk geometry object that lives on NV-GPU: see new() method
 * above for documentation.
 */

struct gk_geometry* gkyl_gk_geometry_mapc2p_cu_dev_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range,
    const struct gkyl_range* range_ext, const struct gkyl_range *global, const struct gkyl_range* global_ext,
    const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx);


