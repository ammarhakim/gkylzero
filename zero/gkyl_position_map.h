#pragma once

#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_comm.h>
#include <gkyl_comm_io.h>

typedef void (*mc2nu_t)(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx);

enum gkyl_position_map_id {
  GKYL_PMAP_FUNC = 0, // Function projection. User specified. Default
  GKYL_PMAP_UNIFORM_B, // Makes a uniform dB in each cell
};

struct gkyl_position_map_inp {
  enum gkyl_position_map_id id;
  mc2nu_t map_x; // Map in the first direction (x, psi, r, etc.)
  mc2nu_t map_y; // Map in the second direction (y, alpha, theta, etc.)
  mc2nu_t map_z; // Map in the third direction (z, field line length, etc.)
  void *ctx_x;  // Context for map_x.
  void *ctx_y;  // Context for map_y.
  void *ctx_z;  // Context for map_z.

  double map_strength; // Zero is uniform mapping, one is fully nonuniform mapping. In between values. Used for the mirror geometry
};

struct gkyl_position_map {
  enum gkyl_position_map_id id;
  bool is_identity; // =True if comp coords = phys coords.

  mc2nu_t map_x; // Map in the first direction (x, psi, r, etc.)
  mc2nu_t map_x_backup; // Backup map in the first direction (x, psi, r, etc.)
  mc2nu_t map_y; // Map in the second direction (y, alpha, theta, etc.)
  mc2nu_t map_y_backup; // Backup map in the second direction (y, alpha, theta, etc.)
  mc2nu_t map_z; // Map in the third direction (z, field line length, etc.)
  void *map_x_ctx;  // Context for map_x.
  void *map_x_ctx_backup;  // Context for map_x.
  void *map_y_ctx;  // Context for map_y.
  void *map_y_ctx_backup;  // Context for map_y.
  void *map_z_ctx;  // Context for map_z.

  double map_strength; // Zero is uniform mapping, one is fully nonuniform mapping. In between values. Used for the constant B mapping

  double cdim; // Number of computational dimensions.
  struct gkyl_rect_grid grid; // Position space grid.
  struct gkyl_range local, local_ext, global, global_ext; // Local & extended local position-space range.
  struct gkyl_basis basis;  // Basis for position mapping.
  struct gkyl_array *mc2nu; // Position mapping in each position direction.
  uint32_t flags;
  struct gkyl_ref_count ref_count;

  // Stuff for constant B mapping
  struct gkyl_array *bmag_global; // Global magnetic field array
  struct gkyl_comm *comm; // Communication object
  bool use_gpu; // Whether to use GPU
  struct gkyl_bmag_ctx *bmag_ctx; // Context for magnetic field calculation
  struct gkyl_position_map_const_B_ctx *constB_ctx; // Context for constant B mapping
};

struct gkyl_position_map_const_B_ctx {
  double psi, alpha;
  double theta_throat, Bmag_throat;
  double psi_min, psi_max;
  double alpha_min, alpha_max;
  double theta_min, theta_max;
  double map_strength;
  int map_order_center, map_order_expander;
};


/**
 * Create a new position map object, containing the things needed
 * to use a (univariate) position map (e.g. to use nonununiform position
 * grids).
 *
 * @param mapc2p_in Comp. to phys. mapping input object (see definition above).
 * @param grid Position space grid.
 * @param local Local position range.
 * @param local_ext Local extended position range.
 * @param global Global position range.
 * @param use_gpu Whether to create a device copy of this new object.
 * @return New position map object.
 */
struct gkyl_position_map* gkyl_position_map_new(struct gkyl_position_map_inp mapc2p_in,
  struct gkyl_rect_grid grid, struct gkyl_range local, struct gkyl_range local_ext, 
  struct gkyl_range global, struct gkyl_range global_ext, struct gkyl_basis basis);

/**
 * Set the position map object.
 * 
 * @param gpm Position map object.
 * @param mc2nu Position map array.
 * 
 * @note This function is used to set the position map array in the position map object.
 */
void gkyl_position_map_set(struct gkyl_position_map* gpm, struct gkyl_array* mc2nu);

/**
 * Evaluate the position mapping at a specific computational (position) coordinate.
 * NOTE: done on the host.
 *
 * @param gpm Gkyl position map object.
 * @param zc Computational position coordinates.
 * @param vp Resulting physical position coordinates.
 */
void
gkyl_position_map_eval_mc2nu(const struct gkyl_position_map* gpm, const double *zc, double *vp);

/**
 * Create a new pointer to the position map object.
 * Release it with gkyl_position_map_release.
 *
 * @param New pointer to the position map object.
 */
struct gkyl_position_map* gkyl_position_map_acquire(const struct gkyl_position_map* gpm);

/**
 * Release pointer to (and eventually memory associated with)
 * the position map object. 
 *
 * @param Position map object.
 */
void gkyl_position_map_release(const struct gkyl_position_map *gpm);

/**
 * Gather bmag into a global array for determining the non-uniform mapping
 * 
 * @param gpm Position map object.
 * @param bmag Local magnetic field array.
 */
void gkyl_position_map_gather_bmag_global(struct gkyl_position_map* gpm, struct gkyl_array* bmag);

/**
 * Optimize the position map object for constant B mapping.
 * 
 * @param gpm Position map object.
 */
void gkyl_position_map_optimize(struct gkyl_position_map* gpm);


void gkyl_position_map_constB_z(double t, const double *xn, double *fout, void *ctx);
