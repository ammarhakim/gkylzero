#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_fromfile.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>


struct gk_geometry*
gkyl_gk_geometry_fromfile_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_gk_geometry_fromfile_cu_dev_new(grid, range, range_ext, basis);
  } 
#endif 

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = *basis;
  up->range = *range;
  up->range_ext = *range_ext;
  up->grid = *grid;

  // bmag, metrics and derived geo quantities
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->range_ext.volume);
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->range_ext.volume);
  up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->range_ext.volume);
  up->dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->range_ext.volume);
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
  up->gxzj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->eps2= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);



  gkyl_grid_sub_array_read(&up->grid, &up->range, up->mc2p, "mapc2p.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->bmag, "bmag.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->g_ij, "g_ij.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->dxdz, "dxdz.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->dzdx, "dzdx.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->jacobgeo, "jacobgeo.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->jacobgeo_inv, "jacogeo_inv.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->gij, "gij.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->b_i, "b_i.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->cmag, "cmag.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->jacobtot, "jacobtot.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->jacobtot_inv, "jacobtot_inv.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->bmag_inv, "bmag_inv.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->bmag_inv_sq, "bmag_inv_sq.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->gxxj, "gxxj.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->gxyj,  "gxyj.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->gyyj,  "gyyj.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->gxzj,  "gxzj.gkyl");
  gkyl_grid_sub_array_read(&up->grid, &up->range, up->eps2,  "eps2.gkyl");

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself
                   
  return up;
}

