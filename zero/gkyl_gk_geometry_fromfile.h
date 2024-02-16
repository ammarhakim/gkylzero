#pragma once

#include <gkyl_gk_geometry.h>
/**
 * Create a new geometry object by reading geometric quantities from file. 
 *
 * @ param geo_host gk_geometry object on the host
 *   Unused if calling with use_gpu=false
 *   If use_gpu=true, geo_host must already be initialized.
 * @param grid Grid on which geometry lives
 * @param local Range on which geometry should be constructed
 * @param local_ext 
 * @param basis configuration space basis
 * @param use_gpu whether or not to use gpu
 */
struct gk_geometry* gkyl_gk_geometry_fromfile_new(struct gk_geometry* geo_host, const struct gkyl_rect_grid* grid, const struct gkyl_range *local, const struct gkyl_range* local_ext, const struct gkyl_range *global, const struct gkyl_range* global_ext, const struct gkyl_basis* basis, bool use_gpu);

/**
 * Create a new gk geometry object that lives on NV-GPU: see new() method
 * above for documentation.
 */

struct gk_geometry* gkyl_gk_geometry_fromfile_cu_dev_new(struct gk_geometry* geo_host, const struct gkyl_rect_grid* grid, const struct gkyl_range *local, const struct gkyl_range* local_ext, const struct gkyl_range *global, const struct gkyl_range* global_ext, const struct gkyl_basis* basis);



