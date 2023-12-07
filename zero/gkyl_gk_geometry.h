#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

typedef struct gk_geometry gk_geometry;

struct gk_geometry {
  // stuff for mapc2p and finite differences array
  struct gkyl_range range;
  struct gkyl_range range_ext;
  struct gkyl_basis basis;
  struct gkyl_rect_grid grid;

  struct gkyl_array* bmag;
  struct gkyl_array* g_ij;
  struct gkyl_array* jacobgeo;
  struct gkyl_array* jacobgeo_inv;
  struct gkyl_array* gij;
  struct gkyl_array* b_i;
  struct gkyl_array* cmag;
  struct gkyl_array* jacobtot;
  struct gkyl_array* jacobtot_inv;
  struct gkyl_array* bmag_inv;
  struct gkyl_array* bmag_inv_sq;
  struct gkyl_array* gxxj;
  struct gkyl_array* gxyj;
  struct gkyl_array* gyyj;

  bool tokamak;

  uint32_t flags;
  struct gkyl_ref_count ref_count;  
  struct gk_geometry *on_dev; // pointer to itself or device object
};


/**
 * Create a new wave geometry object. 
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
struct gk_geometry* gkyl_gk_geometry_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx, bool tokamak, bool use_gpu);

/**
 * Create a new wave geometry object that lives on NV-GPU: see new() method
 * above for documentation.
 */

struct gk_geometry* gkyl_gk_geometry_cu_dev_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx, bool tokamak);

/**
 * Acquire pointer to geometry object. The pointer must be released
 * using gkyl_wave_geom_release method.
 *
 * @param up Geometry to which a pointer is needed
 * @return Pointer to acquired geometry
 */
struct gk_geometry* gkyl_gk_geometry_acquire(const struct gk_geometry* up);



/**
 * Release geometry object.
 *
 * @param wg Wave geometry object to release.
 */
void gkyl_gk_geometry_release(const struct gk_geometry *up);
