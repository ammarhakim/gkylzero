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
void
gkyl_position_map_free(const struct gkyl_ref_count *ref);