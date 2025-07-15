#pragma once

#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_comm.h>
#include <gkyl_comm_io.h>

enum gkyl_position_map_id {
  GKYL_PMAP_USER_INPUT = 0, // Function projection. User specified. Default
  GKYL_PMAP_CONSTANT_DB_POLYNOMIAL, // Makes a uniform dB in each cell. Polynomial approximation, assuming 2 local maxima in Bmag
  GKYL_PMAP_CONSTANT_DB_NUMERIC, // Makes a uniform dB in each cell, but calculates the dB numerically
};

typedef void (*mc2nu_t)(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx);

struct gkyl_position_map_inp {
  enum gkyl_position_map_id id;
  mc2nu_t maps[3]; // Position mapping in each position direction. This is defined in full 3x, 
  // not in deflated coordinates.
  void *ctxs[3]; // Context for each position mapping function.
  double map_strength; // Zero is uniform mapping, one is fully nonuniform mapping. How strong the nonuniformity is
  // Call map_strength = s, xc computational coordinate, and xnu the nonuniform coordinate
  // xnu' = xnu * s + xc * (1-s)
  double maximum_slope_at_min_B; // The maximum slope of the mapping at a magnetic field minimum. A number > 1. Hard limits on cell sizes
  double maximum_slope_at_max_B; // The maximum slope of the mapping at a magnetic field maximum. A number > 1. Hard limits on cell sizes
};

struct gkyl_position_map {
  enum gkyl_position_map_id id;
  mc2nu_t maps[3]; // Position mapping in each position direction.
  void *ctxs[3]; // Context for each position mapping function.

  double cdim; // Number of computational dimensions.
  struct gkyl_rect_grid grid; // Position space grid.
  struct gkyl_range local, local_ext, global, global_ext; // Local & extended local position-space range.
  struct gkyl_basis basis;  // Basis for position mapping.
  struct gkyl_array *mc2nu; // Position mapping in each position direction.
  uint32_t flags;
  struct gkyl_ref_count ref_count;
  bool to_optimize; // Whether to optimize the position map for constant B mapping.

  // Stuff for constant B mapping
  struct gkyl_bmag_ctx *bmag_ctx; // Context for magnetic field calculation
  struct gkyl_position_map_const_B_ctx *constB_ctx; // Context for constant B mapping
};

struct gkyl_position_map_const_B_ctx {
  mc2nu_t maps_backup[3]; // Backup of the position mapping functions.
  void *ctxs_backup[3]; // Backup of the context for each position mapping function.
  
  double psi, alpha; // The psi and alpha values for the middle flux surface to identify the 1D line we are optimizing 
  double psi_min, psi_max; // The max and min psi values for the simulation
  double alpha_min, alpha_max; // The max and min alpha values for the simulation
  double theta_min, theta_max; // The max and min theta values for the simulation
  double map_strength; // Zero is uniform mapping, one is fully nonuniform mapping. How strong the nonuniformity is
  bool enable_maximum_slope_limits_at_min_B; // Whether to enable the maximum slope limits at a magnetic field minimum
  bool enable_maximum_slope_limits_at_max_B; // Whether to enable the maximum slope limits at a magnetic field maximum
  double maximum_slope_at_min_B; // The maximum slope of the mapping at a magnetic field minimum
  double maximum_slope_at_max_B; // The maximum slope of the mapping at a magnetic field maximum

  // Polynomial-based mapping
  double theta_throat, Bmag_throat; // The theta and Bmag values at the throat of the magnetic field
  int map_order_center, map_order_expander; // The polynomial order of the center and expander maps

  // Constant B mapping
  int N_theta_boundaries; // Number of times dB/dz changes sign
  int num_extrema; // Number of extrema in the magnetic field
  double *theta_extrema; // The theta values of the extrema
  double *bmag_extrema; // The Bmag values of the extrema
  bool *min_or_max; // Whether the extrema is a minima or maxima. 1 is maxima, 0 is minima
  double dB_cell; // The change in Bmag per cell
};


/**
 * Create a new position map object. A position map is a function that maps 
 * uniform computational coordinates to non-uniform coordinates in the 
 * same coordinate space as computational coordinates. (e.g. uniform field
 * aligned -> non-uniform field aligned).
 *
 * @param pmap_info Comp. to phys. mapping input object (see definition above).
 * @param grid Position space grid.
 * @param local Local position range.
 * @param local_ext Local extended position range.
 * @param global Global position range.
 * @param use_gpu Whether to create a device copy of this new object.
 * @return New position map object.
 */
struct gkyl_position_map* gkyl_position_map_new(struct gkyl_position_map_inp pmap_info,
  struct gkyl_rect_grid grid, struct gkyl_range local, struct gkyl_range local_ext, 
  struct gkyl_range global, struct gkyl_range global_ext, struct gkyl_basis basis);

/**
 * Set the position map object. Copy the non-uniform map array to the position map object.
 * 
 * @param gpm Position map object.
 * @param mc2nu Position map array.
 * 
 * @note This function is used to set the position map array in the position map object.
 */
void gkyl_position_map_set_mc2nu(struct gkyl_position_map* gpm, struct gkyl_array* mc2nu);

/**
 * Set the magnetic field array in the position map object. This is used to set the magnetic field
 * array in the position map object.
 *
 * @param gpm Position map object.
 * @param comm Communicator object.
 * @param bmag Magnetic field array.
 */
void
gkyl_position_map_set_bmag(struct gkyl_position_map* gpm, struct gkyl_comm* comm,
  struct gkyl_array* bmag);

/**
 * Evaluate the position mapping at a specific computational (position) coordinate.
 * NOTE: done on the host.
 *
 * @param gpm Gkyl position map object.
 * @param xc Computational position coordinates.
 * @param xnu Resulting non-uniform position coordinates.
 */
void
gkyl_position_map_eval_mc2nu(const struct gkyl_position_map* gpm, const double *xc, double *xnu);

/**
 * Evaluate the slope of the position mapping at a specific computational (position) coordinate.
 * 
 * @param gpm Gkyl position map object.
 * @param ix_map Index of the map to evaluate. Calls gpm->maps[index].
 * @param x Computational position coordinates.
 * @param dx Computational position increment.
 * @param ix_comp Index in the geometry loop of which cell we are discussing
 * @param nrange Range of the computational coordinates.
 * @return Slope of the position mapping.
 */
double
gkyl_position_map_slope(const struct gkyl_position_map* gpm, int ix_map,
  double x, double dx, int ix_comp, struct gkyl_range *nrange);

/**
 * Create a new pointer to the position map object.
 * Release it with gkyl_position_map_release.
 *
 * @param gpm Position map object.
 */
struct gkyl_position_map* gkyl_position_map_acquire(const struct gkyl_position_map* gpm);

/**
 * Optimize the position map object for constant B mapping.
 * 
 * @param gpm Position map object.
 * @param grid 3D Position space grid.
 * @param global 3D Global position range.
 */
void gkyl_position_map_optimize(struct gkyl_position_map* gpm, struct gkyl_rect_grid grid,
  struct gkyl_range global);


/**
 * Release pointer to (and eventually memory associated with)
 * the position map object. 
 *
 * @param Position map object.
 */
void gkyl_position_map_release(const struct gkyl_position_map *gpm);