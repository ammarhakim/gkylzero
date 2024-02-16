#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_eqn_type.h>


typedef struct gk_geometry gk_geometry;

struct gk_geometry {
  // stuff for mapc2p and finite differences array
  struct gkyl_range local;
  struct gkyl_range local_ext;
  struct gkyl_range global;
  struct gkyl_range global_ext;
  struct gkyl_basis basis;
  struct gkyl_rect_grid grid;

  struct gkyl_array* mc2p;
  struct gkyl_array* bmag;
  struct gkyl_array* g_ij;
  struct gkyl_array* dxdz;
  struct gkyl_array* dzdx;
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
  struct gkyl_array* gxzj;
  struct gkyl_array* eps2; // eps2 = Jg^33 - J/g_33

  uint32_t flags;
  struct gkyl_ref_count ref_count;  
  struct gk_geometry *on_dev; // pointer to itself or device object
};


// Input struct gor geometry creation
struct gkyl_gyrokinetic_geometry_inp {
  enum gkyl_geometry_id geometry_id;

  void *c2p_ctx; // context for mapc2p function
  // pointer to mapc2p function: xc are the computational space
  // coordinates and on output xp are the corresponding physical space
  // coordinates.
  void (*mapc2p)(double t, const double *xc, double *xp, void *ctx);

  void *bmag_ctx; // context for bmag function
  // pointer to bmag function
  void (*bmag_func)(double t, const double *xc, double *xp, void *ctx);

  struct gkyl_tok_geo_efit_inp *tok_efit_info; // context with RZ data such as efit file for a tokamak
  struct gkyl_tok_geo_grid_inp *tok_grid_info; // context for tokamak geometry with computational domain info

  struct gkyl_mirror_geo_efit_inp *mirror_efit_info; // context with RZ data such as efit file for a mirror
  struct gkyl_mirror_geo_grid_inp *mirror_grid_info; // context for mirror geometry with computational domain info

  double world[3]; // extra computational coordinates for cases with reduced dimensionality

  // 3D grid ranges and basis
  struct gkyl_rect_grid geo_grid;
  struct gkyl_range geo_local;
  struct gkyl_range geo_local_ext;
  struct gkyl_range geo_global;
  struct gkyl_range geo_global_ext;
  struct gkyl_basis geo_basis;

  // Grid ranges and basis with cdim of simulation
  struct gkyl_rect_grid grid;
  struct gkyl_range local;
  struct gkyl_range local_ext;
  struct gkyl_range global;
  struct gkyl_range global_ext;
  struct gkyl_basis basis;

};


struct gkyl_rect_grid 
augment_grid(struct gkyl_rect_grid grid, struct gkyl_gyrokinetic_geometry_inp geometry);

void 
augment_local(const struct gkyl_range *inrange, const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);

/**
 * deflate geometry to lower dimensionality
 * param up_3d 3d geometry object to deflate
 * param grid deflated grid
 * param local deflated local range
 * param local_ext deflated local extended range
 * param basis deflated basis
 * param use_gpu whether or not to use gpu
 */
struct gk_geometry* gkyl_gk_geometry_deflate(const struct gk_geometry* up_3d, struct gkyl_gyrokinetic_geometry_inp *geometry_inp);

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
