/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mapc2p.h>

}
// CPU interface to create and track a GPU object
struct gk_geometry*
gkyl_gk_geometry_mapc2p_cu_dev_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx)
{
  struct gk_geometry *up =(struct gk_geometry*) gkyl_malloc(sizeof(struct gk_geometry));

  up->basis = *basis;
  up->range = *range;
  up->range_ext = *range_ext;
  up->grid = *grid;

  struct gk_geometry *hgeo  = gkyl_gk_geometry_mapc2p_new(grid, range, range_ext, basis, mapc2p_func, mapc2p_ctx, bmag_func, bmag_ctx, false);
  struct gkyl_range nrange;

  int poly_order = basis->poly_order;
  int nodes[3] = { 1, 1, 1 };
  if (poly_order == 1){
    for (int d=0; d<grid->ndim; ++d)
      nodes[d] = grid->cells[d] + 1;
  }
  if (poly_order == 2){
    for (int d=0; d<grid->ndim; ++d)
      nodes[d] = 2*(grid->cells[d]) + 1;
  }

  gkyl_range_init_from_shape(&nrange, up->grid.ndim, nodes);

  // // Initialize the geometry object on the host side
  // mapc2p arrays, bmag, metrics and derived geo quantities
  struct gkyl_array* mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*13, nrange.volume);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange.volume);
  struct gkyl_array* mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->range_ext.volume);

  // Copy the host-side initialized geometry object to the device
  struct gkyl_array *bmag_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *g_ij_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *dxdz_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *jacobgeo_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *jacobgeo_inv_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *gij_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *b_i_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *cmag_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *jacobtot_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *jacobtot_inv_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *bmag_inv_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *bmag_inv_sq_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *gxxj_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *gxyj_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  struct gkyl_array *gyyj_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);

  gkyl_array_copy(bmag_dev, hgeo->bmag);
  gkyl_array_copy(g_ij_dev, hgeo->g_ij);
  gkyl_array_copy(dxdz_dev, hgeo->dxdz);
  gkyl_array_copy(jacobgeo_dev , hgeo->jacobgeo);
  gkyl_array_copy(jacobgeo_inv_dev, hgeo->jacobgeo_inv);
  gkyl_array_copy(gij_dev, hgeo->gij);
  gkyl_array_copy(b_i_dev, hgeo->b_i);
  gkyl_array_copy(cmag_dev, hgeo->cmag);
  gkyl_array_copy(jacobtot_dev, hgeo->jacobtot);
  gkyl_array_copy(jacobtot_inv_dev, hgeo->jacobtot_inv);
  gkyl_array_copy(bmag_inv_dev, hgeo->bmag_inv);
  gkyl_array_copy(bmag_inv_sq_dev, hgeo->bmag_inv_sq);
  gkyl_array_copy(gxxj_dev, hgeo->gxxj);
  gkyl_array_copy(gxyj_dev, hgeo->gxyj);
  gkyl_array_copy(gyyj_dev, hgeo->gyyj);

  gkyl_gk_geometry_release(hgeo);

  // this is for the memcpy below
  up->bmag  = bmag_dev->on_dev;
  up->g_ij  = g_ij_dev->on_dev;
  up->dxdz  = dxdz_dev->on_dev;
  up->jacobgeo  = jacobgeo_dev->on_dev;
  up->jacobgeo_inv = jacobgeo_inv_dev->on_dev;
  up->gij  = gij_dev->on_dev;
  up->b_i  = b_i_dev->on_dev;
  up->cmag  =  cmag_dev->on_dev;
  up->jacobtot  = jacobtot_dev->on_dev;
  up->jacobtot_inv = jacobtot_inv_dev->on_dev;
  up->bmag_inv  = bmag_inv_dev->on_dev;
  up->bmag_inv_sq = bmag_inv_sq_dev->on_dev;
  up->gxxj  = gxxj_dev->on_dev;
  up->gxyj  = gxyj_dev->on_dev;
  up->gyyj  = gyyj_dev->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);

  // Initialize the device geometry object
  struct gk_geometry *up_cu = (struct gk_geometry*) gkyl_cu_malloc(sizeof(struct gk_geometry));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gk_geometry), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;

  // geometry object should store host pointer
  up->bmag  = bmag_dev;
  up->g_ij  = g_ij_dev;
  up->dxdz  = dxdz_dev;
  up->jacobgeo  = jacobgeo_dev;
  up->jacobgeo_inv = jacobgeo_inv_dev;
  up->gij  = gij_dev;
  up->b_i  = b_i_dev;
  up->cmag  =  cmag_dev;
  up->jacobtot  = jacobtot_dev;
  up->jacobtot_inv = jacobtot_inv_dev;
  up->bmag_inv  = bmag_inv_dev;
  up->bmag_inv_sq = bmag_inv_sq_dev;
  up->gxxj  = gxxj_dev;
  up->gxyj  = gxyj_dev;
  up->gyyj  = gyyj_dev;
  
  return up;
}
