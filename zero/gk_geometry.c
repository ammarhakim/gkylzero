#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_math.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_alloc_flags_priv.h>


bool
gkyl_gk_geometry_is_cu_dev(const struct gk_geometry* up)
{
  return GKYL_IS_CU_ALLOC(up->flags);
}

struct gk_geometry*
gkyl_gk_geometry_deflate(const struct gk_geometry* up_3d, const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_gk_geometry_deflate_cu_dev(up_3d, grid, range, range_ext, basis);
  } 
#endif 

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = *basis;
  up->range = *range;
  up->range_ext = *range_ext;
  up->grid = *grid;

  //struct gkyl_array* mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->range_ext.volume);
  // bmag, metrics and derived geo quantities
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->range_ext.volume);
  up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->range_ext.volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->range_ext.volume);
  up->b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->range_ext.volume);
  up->cmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->gxxj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->gxyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->gyyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);

  // Now fill the arrays by deflation
  int rem_dirs[3] = {0};
  if (up->grid.ndim==1) {
    rem_dirs[0] = 1;
    rem_dirs[1] = 1;
  }
  else if (up->grid.ndim==2) {
    rem_dirs[1] = 1;
  }
  struct gkyl_deflate_geo* deflator = gkyl_deflate_geo_new(&up_3d->basis, &up->basis, &up_3d->grid, &up->grid, rem_dirs, false);

  //gkyl_deflate_geo_advance(deflator, up_3d->range, up->range, mc2p, mc2p, 3); // mc2p is not currently stored. fix later so we can deflate
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->bmag, up_3d->bmag, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->g_ij, up_3d->g_ij, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->dxdz, up_3d->dxdz, 9);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->jacobgeo, up_3d->jacobgeo, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->jacobgeo_inv, up_3d->jacobgeo_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->gij, up->gij, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->b_i, up->b_i, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->cmag, up->cmag, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->jacobtot, up->jacobtot, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->jacobtot_inv, up->jacobtot_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->bmag_inv, up->bmag_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->bmag_inv_sq, up->bmag_inv_sq, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->gxxj, up->gxxj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->gxyj, up->gxyj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->range, &up->range, up_3d->gyyj, up->gyyj, 1);
  // Done deflating

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself
                   //
  //gkyl_grid_sub_array_write(&up->grid, &up->range, mc2p, "mapc2p.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->bmag, "bmag.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->g_ij, "g_ij.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->dxdz, "dxdz.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->jacobgeo, "jacobgeo.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->jacobgeo_inv, "jacogeo_inv.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->gij, "gij.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->b_i, "b_i.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->cmag, "cmag.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->jacobtot, "jacobtot.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->jacobtot_inv, "jacobtot_inv.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->bmag_inv, "bmag_inv.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->bmag_inv_sq, "bmag_inv_sq.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->gxxj, "gxxj.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->gxyj,  "gxyj.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->gyyj,  "gyyj.gkyl");

  return up;
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



