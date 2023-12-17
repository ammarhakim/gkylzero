#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_math.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_alloc_flags_priv.h>


bool
gkyl_gk_geometry_is_cu_dev(const struct gk_geometry* up)
{
  return GKYL_IS_CU_ALLOC(up->flags);
}

void
gkyl_gk_geometry_free(const struct gkyl_ref_count *ref)
{
  struct gk_geometry *up = container_of(ref, struct gk_geometry, ref_count);
  gkyl_array_release(up->bmag);
  gkyl_array_release(up->g_ij);
  gkyl_array_release(up->jacobgeo);
  gkyl_array_release(up->jacobgeo_inv);
  gkyl_array_release(up->gij);
  gkyl_array_release(up->b_i);
  gkyl_array_release(up->cmag);
  gkyl_array_release(up->jacobtot);
  gkyl_array_release(up->jacobtot_inv);
  gkyl_array_release(up->bmag_inv);
  gkyl_array_release(up->bmag_inv_sq);
  gkyl_array_release(up->gxxj);
  gkyl_array_release(up->gxyj);
  gkyl_array_release(up->gyyj);
  if (gkyl_gk_geometry_is_cu_dev(up)) 
    gkyl_cu_free(up->on_dev); 

  gkyl_free(up);
}

struct gk_geometry*
gkyl_gk_geometry_acquire(const struct gk_geometry* up)
{
  gkyl_ref_count_inc(&up->ref_count);
  return (struct gk_geometry*) up;
}

void
gkyl_gk_geometry_release(const struct gk_geometry *up)
{
  gkyl_ref_count_dec(&up->ref_count);
}



