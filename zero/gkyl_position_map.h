#pragma once

#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_comm.h>
#include <gkyl_comm_io.h>

typedef void (*mapc2fa_t)(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx);

// Position space mappings.
struct gkyl_position_map_inp {
  mapc2fa_t mapping; // Must be 3x mapping function.
  void *ctx;  // Context for mapping.
  double numerical_mapping_fraction; // Zero is uniform mapping, one is fully nonuniform mapping. In between values. Used for the mirror geometry
};

// Object type.
typedef struct gkyl_position_map gkyl_position_map;

struct gkyl_position_map {
  bool is_identity; // =True if comp coords = phys coords.
  struct gkyl_rect_grid grid; // Position space grid.
  struct gkyl_range local, local_ext; // Local & extended local position-space range.
  struct gkyl_basis *pmap_basis;  // Basis for position mapping.
  struct gkyl_array *pmap; // Position mapping in each position direction.
  struct gkyl_position_map *on_dev; // Device copy of itself.
  // For internal/private use only:
  struct gkyl_array *pmap_ho; // Host copy of pmap.
  struct gkyl_basis pmap_basis_ho;  // Host basis for position mapping.
  uint32_t flags;
  struct gkyl_ref_count ref_count;
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
 * @param use_gpu Whether to create a device copy of this new object.
 * @return New position map object.
 */
struct gkyl_position_map* gkyl_position_map_new(struct gkyl_position_map_inp mapc2p_in,
  struct gkyl_rect_grid grid, struct gkyl_range local, 
  struct gkyl_range local_ext, struct gkyl_basis basis, bool use_gpu);

/**
 * Set the position map object.
 * 
 * @param gpm Position map object.
 * @param pmap Position map array.
 * 
 * @note This function is used to set the position map array in the position map object.
 */
void gkyl_position_map_set(struct gkyl_position_map* gpm, struct gkyl_array* pmap);


/**
 * Write the position map and its jacobian to file.
 *
 * @param gpm Gkyl position map object.
 * @param comm Communicator.
 * @param app_name Name of the app.
 */
void
gkyl_position_map_write(const struct gkyl_position_map* gpm, struct gkyl_comm* comm,
  const char* app_name);

/**
 * Evaluate the position mapping at a specific computational (position) coordinate.
 * NOTE: done on the host.
 *
 * @param gpm Gkyl position map object.
 * @param zc Computational position coordinates.
 * @param vp Resulting physical position coordinates.
 */
void
gkyl_position_map_eval_c2p(const struct gkyl_position_map* gpm, const double *zc, double *vp);

/**
 * Indicate if this position map object is allocated on the GPU.
 *
 * @param gpm Gkyl position map object.
 * @return
 */
bool gkyl_position_map_is_cu_dev(const struct gkyl_position_map* gpm);

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
