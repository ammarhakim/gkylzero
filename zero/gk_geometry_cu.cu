/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_gk_geometry.h>
}

struct gk_geom_surf*
gk_geometry_surf_cu_dev_alloc(struct gk_geom_surf* up_surf_host, struct gk_geometry* up_dev)
{
  struct gk_geom_surf *up_surf = (struct gk_geom_surf*) gkyl_malloc(sizeof(struct gk_geom_surf));
  up_surf->bmag = gkyl_array_cu_dev_new(up_surf_host->bmag->type, up_surf_host->bmag->ncomp, up_surf_host->bmag->size);
  up_surf->jacobgeo = gkyl_array_cu_dev_new(up_surf_host->jacobgeo->type, up_surf_host->jacobgeo->ncomp, up_surf_host->jacobgeo->size);
  up_surf->b_i = gkyl_array_cu_dev_new(up_surf_host->b_i->type, up_surf_host->b_i->ncomp, up_surf_host->b_i->size);
  up_surf->cmag = gkyl_array_cu_dev_new(up_surf_host->cmag->type, up_surf_host->cmag->ncomp, up_surf_host->cmag->size);
  up_surf->jacobtot_inv = gkyl_array_cu_dev_new(up_surf_host->jacobtot_inv->type, up_surf_host->jacobtot_inv->ncomp, up_surf_host->jacobtot_inv->size);
  return up_surf;
}

// CPU interface to create and track a GPU object
struct gk_geometry* 
gkyl_gk_geometry_cu_dev_new(struct gk_geometry* geo_host, struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry *up =(struct gk_geometry*) gkyl_malloc(sizeof(struct gk_geometry));

  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->global = geometry_inp->global;
  up->global_ext = geometry_inp->global_ext;
  up->grid = geometry_inp->grid;
  gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, up->basis.poly_order);

  // Copy the host-side initialized geometry object to the device
  struct gkyl_array *mc2p_dev = gkyl_array_cu_dev_new(geo_host->mc2p->type, geo_host->mc2p->ncomp, geo_host->mc2p->size);
  struct gkyl_array *mc2nu_pos_dev = gkyl_array_cu_dev_new(geo_host->mc2nu_pos->type, geo_host->mc2nu_pos->ncomp, geo_host->mc2nu_pos->size);
  struct gkyl_array *bmag_dev = gkyl_array_cu_dev_new(geo_host->bmag->type, geo_host->bmag->ncomp, geo_host->bmag->size);
  struct gkyl_array *g_ij_dev = gkyl_array_cu_dev_new(geo_host->g_ij->type, geo_host->g_ij->ncomp, geo_host->g_ij->size);
  struct gkyl_array *dxdz_dev = gkyl_array_cu_dev_new(geo_host->dxdz->type, geo_host->dxdz->ncomp, geo_host->dxdz->size);
  struct gkyl_array *dzdx_dev = gkyl_array_cu_dev_new(geo_host->dzdx->type, geo_host->dzdx->ncomp, geo_host->dzdx->size);
  struct gkyl_array *dualmag_dev = gkyl_array_cu_dev_new(geo_host->dualmag->type, geo_host->dualmag->ncomp, geo_host->dualmag->size);
  struct gkyl_array *normals_dev = gkyl_array_cu_dev_new(geo_host->normals->type, geo_host->normals->ncomp, geo_host->normals->size);
  struct gkyl_array *jacobgeo_dev = gkyl_array_cu_dev_new(geo_host->jacobgeo->type, geo_host->jacobgeo->ncomp, geo_host->jacobgeo->size);
  struct gkyl_array *jacobgeo_inv_dev = gkyl_array_cu_dev_new(geo_host->jacobgeo_inv->type, geo_host->jacobgeo_inv->ncomp, geo_host->jacobgeo_inv->size);
  struct gkyl_array *gij_dev = gkyl_array_cu_dev_new(geo_host->gij->type, geo_host->gij->ncomp, geo_host->gij->size);
  struct gkyl_array *b_i_dev = gkyl_array_cu_dev_new(geo_host->b_i->type, geo_host->b_i->ncomp, geo_host->b_i->size);
  struct gkyl_array *bcart_dev = gkyl_array_cu_dev_new(geo_host->bcart->type, geo_host->bcart->ncomp, geo_host->bcart->size);
  struct gkyl_array *cmag_dev = gkyl_array_cu_dev_new(geo_host->cmag->type, geo_host->cmag->ncomp, geo_host->cmag->size);
  struct gkyl_array *jacobtot_dev = gkyl_array_cu_dev_new(geo_host->jacobtot->type, geo_host->jacobtot->ncomp, geo_host->jacobtot->size);
  struct gkyl_array *jacobtot_inv_dev = gkyl_array_cu_dev_new(geo_host->jacobtot_inv->type, geo_host->jacobtot_inv->ncomp, geo_host->jacobtot_inv->size);
  struct gkyl_array *bmag_inv_dev = gkyl_array_cu_dev_new(geo_host->bmag_inv->type, geo_host->bmag_inv->ncomp, geo_host->bmag_inv->size);
  struct gkyl_array *bmag_inv_sq_dev = gkyl_array_cu_dev_new(geo_host->bmag_inv_sq->type, geo_host->bmag_inv_sq->ncomp, geo_host->bmag_inv_sq->size);
  struct gkyl_array *gxxj_dev = gkyl_array_cu_dev_new(geo_host->gxxj->type, geo_host->gxxj->ncomp, geo_host->gxxj->size);
  struct gkyl_array *gxyj_dev = gkyl_array_cu_dev_new(geo_host->gxyj->type, geo_host->gxyj->ncomp, geo_host->gxyj->size);
  struct gkyl_array *gyyj_dev = gkyl_array_cu_dev_new(geo_host->gyyj->type, geo_host->gyyj->ncomp, geo_host->gyyj->size);
  struct gkyl_array *gxzj_dev = gkyl_array_cu_dev_new(geo_host->gxzj->type, geo_host->gxzj->ncomp, geo_host->gxzj->size);
  struct gkyl_array *eps2_dev = gkyl_array_cu_dev_new(geo_host->eps2->type, geo_host->eps2->ncomp, geo_host->eps2->size);

  struct gk_geom_surf *geo_surf_dev[up->grid.ndim];
  for (int dir=0; dir<up->grid.ndim; ++dir)
    geo_surf_dev[dir] = gk_geometry_surf_cu_dev_alloc(geo_host->geo_surf[dir], up);

  gkyl_array_copy(mc2p_dev, geo_host->mc2p);
  gkyl_array_copy(mc2nu_pos_dev, geo_host->mc2nu_pos);
  gkyl_array_copy(bmag_dev, geo_host->bmag);
  gkyl_array_copy(g_ij_dev, geo_host->g_ij);
  gkyl_array_copy(dxdz_dev, geo_host->dxdz);
  gkyl_array_copy(dzdx_dev, geo_host->dzdx);
  gkyl_array_copy(dualmag_dev, geo_host->dualmag);
  gkyl_array_copy(normals_dev, geo_host->normals);
  gkyl_array_copy(jacobgeo_dev , geo_host->jacobgeo);
  gkyl_array_copy(jacobgeo_inv_dev, geo_host->jacobgeo_inv);
  gkyl_array_copy(gij_dev, geo_host->gij);
  gkyl_array_copy(b_i_dev, geo_host->b_i);
  gkyl_array_copy(bcart_dev, geo_host->bcart);
  gkyl_array_copy(cmag_dev, geo_host->cmag);
  gkyl_array_copy(jacobtot_dev, geo_host->jacobtot);
  gkyl_array_copy(jacobtot_inv_dev, geo_host->jacobtot_inv);
  gkyl_array_copy(bmag_inv_dev, geo_host->bmag_inv);
  gkyl_array_copy(bmag_inv_sq_dev, geo_host->bmag_inv_sq);
  gkyl_array_copy(gxxj_dev, geo_host->gxxj);
  gkyl_array_copy(gxyj_dev, geo_host->gxyj);
  gkyl_array_copy(gyyj_dev, geo_host->gyyj);
  gkyl_array_copy(gxzj_dev, geo_host->gxzj);
  gkyl_array_copy(eps2_dev, geo_host->eps2);
  for (int dir=0; dir<up->grid.ndim; ++dir) {
   gkyl_array_copy(geo_surf_dev[dir]->bmag, geo_host->geo_surf[dir]->bmag);
   gkyl_array_copy(geo_surf_dev[dir]->jacobgeo, geo_host->geo_surf[dir]->jacobgeo);
   gkyl_array_copy(geo_surf_dev[dir]->b_i, geo_host->geo_surf[dir]->b_i);
   gkyl_array_copy(geo_surf_dev[dir]->cmag, geo_host->geo_surf[dir]->cmag);
   gkyl_array_copy(geo_surf_dev[dir]->jacobtot_inv, geo_host->geo_surf[dir]->jacobtot_inv);
  }

  // this is for the memcpy below
  up->mc2p  = mc2p_dev->on_dev;
  up->mc2nu_pos  = mc2nu_pos_dev->on_dev;
  up->bmag  = bmag_dev->on_dev;
  up->g_ij  = g_ij_dev->on_dev;
  up->dxdz  = dxdz_dev->on_dev;
  up->dzdx  = dzdx_dev->on_dev;
  up->dualmag  = dualmag_dev->on_dev;
  up->normals  = normals_dev->on_dev;
  up->jacobgeo  = jacobgeo_dev->on_dev;
  up->jacobgeo_inv = jacobgeo_inv_dev->on_dev;
  up->gij  = gij_dev->on_dev;
  up->b_i  = b_i_dev->on_dev;
  up->bcart  = bcart_dev->on_dev;
  up->cmag  =  cmag_dev->on_dev;
  up->jacobtot  = jacobtot_dev->on_dev;
  up->jacobtot_inv = jacobtot_inv_dev->on_dev;
  up->bmag_inv  = bmag_inv_dev->on_dev;
  up->bmag_inv_sq = bmag_inv_sq_dev->on_dev;
  up->gxxj  = gxxj_dev->on_dev;
  up->gxyj  = gxyj_dev->on_dev;
  up->gyyj  = gyyj_dev->on_dev;
  up->gxzj  = gxzj_dev->on_dev;
  up->eps2  = eps2_dev->on_dev;
  for (int dir=0; dir<up->grid.ndim; ++dir) {
   up->geo_surf[dir]->bmag = geo_surf_dev[dir]->bmag->on_dev;
   up->geo_surf[dir]->jacobgeo = geo_surf_dev[dir]->jacobgeo->on_dev;
   up->geo_surf[dir]->b_i = geo_surf_dev[dir]->b_i->on_dev;
   up->geo_surf[dir]->cmag = geo_surf_dev[dir]->cmag->on_dev;
   up->geo_surf[dir]->jacobtot_inv = geo_surf_dev[dir]->jacobtot_inv->on_dev;
  }

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);

  // Initialize the device geometry object
  struct gk_geometry *up_cu = (struct gk_geometry*) gkyl_cu_malloc(sizeof(struct gk_geometry));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gk_geometry), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;

  // geometry object should store host pointer
  up->mc2p  = mc2p_dev;
  up->mc2nu_pos  = mc2nu_pos_dev;
  up->bmag  = bmag_dev;
  up->g_ij  = g_ij_dev;
  up->dxdz  = dxdz_dev;
  up->dzdx  = dzdx_dev;
  up->dualmag  = dualmag_dev;
  up->normals  = normals_dev;
  up->jacobgeo  = jacobgeo_dev;
  up->jacobgeo_inv = jacobgeo_inv_dev;
  up->gij  = gij_dev;
  up->b_i  = b_i_dev;
  up->bcart  = bcart_dev;
  up->cmag  =  cmag_dev;
  up->jacobtot  = jacobtot_dev;
  up->jacobtot_inv = jacobtot_inv_dev;
  up->bmag_inv  = bmag_inv_dev;
  up->bmag_inv_sq = bmag_inv_sq_dev;
  up->gxxj  = gxxj_dev;
  up->gxyj  = gxyj_dev;
  up->gyyj  = gyyj_dev;
  up->gxzj  = gxzj_dev;
  up->eps2  = eps2_dev;
  for (int dir=0; dir<up->grid.ndim; ++dir) {
   up->geo_surf[dir]->bmag = geo_surf_dev[dir]->bmag;
   up->geo_surf[dir]->jacobgeo = geo_surf_dev[dir]->jacobgeo;
   up->geo_surf[dir]->b_i = geo_surf_dev[dir]->b_i;
   up->geo_surf[dir]->cmag = geo_surf_dev[dir]->cmag;
   up->geo_surf[dir]->jacobtot_inv = geo_surf_dev[dir]->jacobtot_inv;
  }
  
  return up;
}

