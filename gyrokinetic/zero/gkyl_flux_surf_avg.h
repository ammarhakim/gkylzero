#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_comm.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_flux_surf_avg gkyl_flux_surf_avg;

// Input parameters for the flux-surface average updater.
struct gkyl_flux_surf_avg_inp {
  const struct gkyl_rect_grid *grid; // Global configuration-space grid (2x or 3x).
  struct gkyl_basis basis; // Configuration-space basis (p=1).
  struct gkyl_basis basis_x; // 1D basis in x of the averages.
  const struct gkyl_range *local; // Local configuration-space range.
  const struct gkyl_range *local_ext; // Local extended configuration-space range.
  const struct gkyl_range *local_x; // 1D range in x of the averages.
  const struct gkyl_range *local_x_ext; // 1D extended range in x of the averages.
  const struct gkyl_array *jacobgeo; // Jacobian J weighting the average.
  struct gkyl_comm *comm; // Communicator of the ranks decomposing the averaged directions.
  bool use_gpu; // Whether to run on the GPU.
};

/**
 * Create a new updater computing the flux-surface average (FSA)
 *   <f>(x) = int J f dy dz / int J dy dz
 * of a 2x or 3x field (R operator), and the extension of a function of x to a field constant
 * along the other directions (E operator). The averaged directions may be decomposed across
 * the ranks of comm. Free using gkyl_flux_surf_avg_release.
 *
 * @param inp Input parameters.
 * @return New updater pointer.
 */
struct gkyl_flux_surf_avg *gkyl_flux_surf_avg_inew(const struct gkyl_flux_surf_avg_inp *inp);

/**
 * Compute the flux-surface average of a field.
 *
 * @param up Updater object.
 * @param fin Input field (cdim, on local).
 * @param avg Output average (1D in x, on local_x).
 */
void gkyl_flux_surf_avg_advance(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *fin, struct gkyl_array *avg
);

/**
 * Extend a function of x to a field constant along the other directions.
 *
 * @param up Updater object.
 * @param avg Input function of x (on local_x).
 * @param out Output field (cdim, on local).
 */
void gkyl_flux_surf_avg_inflate(
  struct gkyl_flux_surf_avg *up, const struct gkyl_array *avg, struct gkyl_array *out
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_flux_surf_avg_release(struct gkyl_flux_surf_avg *up);
