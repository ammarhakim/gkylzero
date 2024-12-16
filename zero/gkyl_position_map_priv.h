#include <gkyl_array.h>

// Allocate double array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

/**
 * Function that actually frees memory associated with this
 * object when the number of references has decreased to zero.
 *
 * @param ref Reference counter for this object.
 */
void gkyl_position_map_free(const struct gkyl_ref_count *ref);


/**
 * In a geometry with two symmetric magnetic field peaks, calculate the location
 * of highest mangetic field along the field line
 * 
 * @param constB_ctx Context object from the position_map
 * @param bmag_ctx Context from calc_bmag
 */
void calculate_mirror_throat_location(struct gkyl_position_map_const_B_ctx *constB_ctx, struct gkyl_bmag_ctx *bmag_ctx);

/**
 * In a geometry with two symmetric magnetic field peaks, calculate
 * the mapping that makes dB/dtheta constant along the field line
 * 
 * @param constB_ctx Context object from the position_map
 * @param bmag_ctx Context from calc_bmag
 */
void calculate_optimal_mapping(struct gkyl_position_map_const_B_ctx *constB_ctx, struct gkyl_bmag_ctx *bmag_ctx);

/**
 * Mapping for a constant magnetic field change over each cell
 * 
 * @param t Time.
 * @param xn Computational coordinates.
 * @param fout Physical coordinates.
 * @param ctx Context object.
*/
void gkyl_position_map_constB_z(double t, const double *xn, double *fout, void *ctx);