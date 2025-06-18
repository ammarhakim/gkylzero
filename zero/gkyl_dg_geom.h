#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_math.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Geometry information for lower surfaces (left, bottom, front) for a cell
struct gkyl_dg_surf_geom {
  double area_elem; // area element: Jc*norm(e^d), where e^d is dual normal to face 'd'
  struct gkyl_vec3 norm; // norm[d] is the normal to face perp to direction 'd'
  // tau1[d] X tau2[d] = norm[d] are tangents to face perp to direction 'd'
  struct gkyl_vec3 tau1, tau2;
};

// Geometry information for a single cell
struct gkyl_dg_cell_geom {
  struct gkyl_vec3 tang[3]; // tangent vectors, e_i
  struct gkyl_vec3 dual[3]; // dual vectors, e^i
  double Jc; // Jacobian = e_1*(e_2 X e_3)  = 1/e^1*(e^2 X e^3)
};  

// geometry information over a range of cells: 
struct gkyl_dg_geom {
  struct gkyl_range range; // range over which geometry is defined
  struct gkyl_array *surf_geom[GKYL_MAX_CDIM]; // surface geometry in dir 'd' in each cell
  struct gkyl_array *cell_geom; // cell geometry
  
  struct gkyl_range surf_quad_range; // range for indexing surface nodes
  struct gkyl_range vol_quad_range; // range for indexing volume nodes
  
  uint32_t flags;
  struct gkyl_ref_count ref_count;  
  struct gkyl_dg_geom *on_dev; // pointer to itself or device object
};

// Input for constructing geometry with mapc2p
struct gkyl_dg_geom_inp {
  const struct gkyl_rect_grid *grid; // grid over which geometry needs to be defined
  struct gkyl_range *range; // range of grid: this needs to be extended grid
  bool skip_geometry_creation; // true if only allocation is needed (geom must be computed eleswhere)
  int nquad; // number of quadrature points in each direction
  evalf_t mapc2p; // mapping function: set mapc2p = 0 to use identity map
  void *ctx; // context for use in mapc2p
  bool use_gpu; // does this live on the GPU?
};

/**
 * Create a new DG geometry object from mapc2p. 
 *
 * @param inp Inputs for use in constructing geometry
 */
struct gkyl_dg_geom* gkyl_dg_geom_new(const struct gkyl_dg_geom_inp *inp);

/**
 * Acquire pointer to geometry object. The pointer must be released
 * using gkyl_dg_geom_release method.
 *
 * @param dgg Geometry to which a pointer is needed
 * @return Pointer to acquired geometry
 */
struct gkyl_dg_geom* gkyl_dg_geom_acquire(const struct gkyl_dg_geom* dgg);

/**
 * Get pointer to geometry on the surface normal to 'd' given by idx
 * into the range over which the geometry was constructed. The
 * returned array needs to be further indexed to fetch the specific
 * struct at a quadrature point
 *
 * @param dgg DG geometry object
 * @param d Direction to which surface is perpendicular
 * @param idx Index into grid
 * @return Pointer to surface geometry at all quadrature nodes in cell @a idx
 */
GKYL_CU_DH
static inline const struct gkyl_dg_surf_geom*
gkyl_dg_geom_get_surf(const struct gkyl_dg_geom *dgg, int d, const int *idx)
{
  return (const struct gkyl_dg_surf_geom*) gkyl_array_cfetch(dgg->surf_geom[d], gkyl_range_idx(&dgg->range, idx));
}

/**
 * Get pointer to geometry on the surface normal to 'd' given by idx
 * into the range over which the geometry was constructed. The
 * returned array needs to be further indexed to fetch the specific
 * struct at a quadrature point
 *
 * @param dgg DG geometry object
 * @param idx Index (ndim-1) of surface quadrature node
 * @return Linear index for indexing surface quadrature array
 */
GKYL_CU_DH
static inline long
gkyl_dg_geom_surf_quad_idx(const struct gkyl_dg_geom *dgg, const int *idx)
{
  return gkyl_range_idx(&dgg->surf_quad_range, idx);
}

/**
 * Release geometry object.
 *
 * @param dgg Dg geometry object to release.
 */
void gkyl_dg_geom_release(const struct gkyl_dg_geom *dgg);
