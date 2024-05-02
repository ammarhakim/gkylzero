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

  // These 20 DG fields contain the geometric quantities needed to solve the
  // GK Equation and Poisson Equation and to apply certain BC's
  // The first 19 are defined on the configuration space domain. The last is a single element.
  struct gkyl_array* mc2p; // 3 components. Cartesian X,Y, and Z
  struct gkyl_array* bmag; // 1 component. B Magnitude of magnetic field
  struct gkyl_array* g_ij; // 6 components. 
                           // Metric coefficients g_{ij} Stored in order g_11, g12, g_13, g_22, g_23, g_33
  struct gkyl_array* dxdz; // 9 components.
                           // Cartesian components of tangent Vectors stored in order e_1, e_2, e_3
  struct gkyl_array* dzdx; // 9 components.
                           // Cartesian components of dual vectors stroed in order e^1, e^2, e^3
  struct gkyl_array* jacobgeo; // 1 component. Configuration space jacobian J
  struct gkyl_array* jacobgeo_inv; // 1 component. 1/J
  struct gkyl_array* gij; // Matric coefficients g^{ij}. See g_ij for order.
  struct gkyl_array* b_i; // 3 components. Contravariant components of magnetic field vector b_1, b_2, b_3.
  struct gkyl_array* cmag; // 1 component. C = JB/sqrt(g_33)
  struct gkyl_array* jacobtot; // 1 component. Phase space Jacobian = JB
  struct gkyl_array* jacobtot_inv; // 1 component. 1/(JB)
  struct gkyl_array* bmag_inv; // 1 component. 1/B.
  struct gkyl_array* bmag_inv_sq; // 1 component. 1/B^2.
  struct gkyl_array* gxxj; // 1 component. g^{xx} * J. For poisson solve.
  struct gkyl_array* gxyj; // 1 component. g^{xy} * J. For poisson solve.
  struct gkyl_array* gyyj; // 1 component. g^{yy} * J. For poisson solve.
  struct gkyl_array* gxzj; // 1 component. g^{xz} * J. For poisson solve if z derivatives are kept.
  struct gkyl_array* eps2; // 1 component. eps2 = Jg^33 - J/g_33. For poisson if z derivatives are kept.
  struct gkyl_array* bmag_mid; // 1 component. B at center of domain.

  uint32_t flags;
  struct gkyl_ref_count ref_count;  
  struct gk_geometry *on_dev; // pointer to itself or device object
};


// Input struct gor geometry creation
struct gkyl_gk_geometry_inp {
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


/**
 * Augment a grid with dim < 3 to 3d by adding 1 cell in the other directions
 * If dim=1, the input grid is assumed to be in z
 * If dim =2, the input grid is assumed to be in xz
 * @param grid input grid with dim <3
 * @ param geometry geometry input struct with context for augmenting grid
 */
struct gkyl_rect_grid 
gkyl_gk_geometry_augment_grid(struct gkyl_rect_grid grid, struct gkyl_gk_geometry_inp geometry);


/**
 * Augment a range with dim < 3 to 3d by adding 1 cell in the other directions
 * If dim=1, the input range is assumed to be in z
 * If dim =2, the input range is assumed to be in xz
 * @param inrange input range with dim <3
 * @param nghost number of ghost cells
 * @param ext_range output, augmented extended range
 * @param range output, augmented range
 */
void 
gkyl_gk_geometry_augment_local(const struct gkyl_range *inrange, const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);


/**
 * Evaluate and set bmag at the center of the domain
 */
void gkyl_gk_geometry_bmag_mid(struct gk_geometry* up);

/**
 * deflate geometry to lower dimensionality
 * param up_3d 3d geometry object to deflate
 * param grid deflated grid
 * param local deflated local range
 * param local_ext deflated local extended range
 * param basis deflated basis
 * param use_gpu whether or not to use gpu
 */
struct gk_geometry* gkyl_gk_geometry_deflate(const struct gk_geometry* up_3d, struct gkyl_gk_geometry_inp *geometry_inp);

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
