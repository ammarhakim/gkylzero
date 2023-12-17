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
 * Acquire pointer to gk geometry object. The pointer must be released
 * using gkyl_gk_geometry_release method.
 *
 * @param up Geometry to which a pointer is needed
 * @return Pointer to acquired geometry
 */
struct gk_geometry* gkyl_gk_geometry_acquire(const struct gk_geometry* up);



void gkyl_gk_geometry_free(const struct gkyl_ref_count *ref);

/**
 * Release gk geometry object.
 *
 * @param up gk geometry object to release.
 */
void gkyl_gk_geometry_release(const struct gk_geometry *up);
