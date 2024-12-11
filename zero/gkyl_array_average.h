#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_array_average gkyl_array_average;

// Input of the new routine is packaged as a struct
typedef struct gkyl_array_average_inp gkyl_array_average_inp;
/**
 * The input structure contains the parameters required to perform a weighted array average:
 * 
 * @param grid       Pointer to the computational grid, used to compute the surface element.
 * @param basis      Total basis, describes the full dimensionality (ndim), polynomial type, and order.
 * @param basis_avg  Subset basis, describes the reduced dimensionality, polynomial type, and order for the output.
 * @param local      Full range of the input array, covering all dimensions of the total basis.
 * @param local_ext  Extended range of the input array, including ghost cells.
 * @param local_avg  Reduced range of the output array, covering only the non-averaged dimensions.
 * @param weights    Pointer to the array containing weights for the averaging process. (set it to NULL for integral)
 * @param avg_dim    Flag array to set which dimension is averaged
 * @param use_gpu    Boolean flag indicating whether the computation should be performed on a GPU.
 */
struct gkyl_array_average_inp {
  const struct gkyl_rect_grid *grid;    // Computational grid
  const struct gkyl_basis basis;        // Total basis for the full dimensionality
  const struct gkyl_basis basis_avg;    // Subset basis for reduced dimensionality
  const struct gkyl_range *local;       // Range for input array (total)
  const struct gkyl_range *local_ext;   // Extended range for input array (with ghosts)
  const struct gkyl_range *local_avg;   // Range for output array (reduced)
  const struct gkyl_array *weights;     // Weight array for averaging
  const int *avg_dim;                   // That indicates if a dimension is averaged or not, must be GKYL_MAX_CDIM size
  bool use_gpu;                         // Flag for GPU computation
};

/**
 * Create a new updater that computes the average of a gkyl_array, with an option to
 * set the dimensions to average
 * the input contains:
 * @param inp see gkyl_array_average_inp structure
 */
struct gkyl_array_average*
gkyl_array_average_new(const struct gkyl_array_average_inp *inp);
/**
 * Compute the array average.
 *
 * @param up array_integrate updater.
 * @param fin input gkyl_array
 * @param avgout Output gkyl_array
 */
void gkyl_array_average_advance(const struct gkyl_array_average *up, 
  const struct gkyl_array *fin, struct gkyl_array *avgout);

/**
 * Release memory associated with this updater.
 *
 * @param up array_integrate updater.
 */
void gkyl_array_average_release(struct gkyl_array_average *up);
