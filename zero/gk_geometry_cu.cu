/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_math.h>
#include <gkyl_util.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
}

#include <cassert>

// CPU interface to create and track a GPU object
struct gk_geometry*
gkyl_gk_geometry_cu_dev_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx)
{
  struct gk_geometry *up =(struct gk_geometry*) gkyl_malloc(sizeof(struct gk_geometry));

  up->basis = basis;
  up->range = range;
  up->range_ext = range_ext;
  up->grid = grid;
  struct gkyl_range nrange;
  double dzc[3] = {0.0};

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

  gkyl_range_init_from_shape(&nrange, up->grid->ndim, nodes);

  // Initialize the geometry object on the host side
  // mapc2p arrays, bmag, metrics and derived geo quantities
  struct gkyl_array* mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim*13, nrange.volume);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim, nrange.volume);
  struct gkyl_array* mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim, up->range_ext->volume);
  struct gkyl_array* bmag = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* cmag = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* gxxj = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* gxyj = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array* gyyj = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);

  if (up->tokamak){
    const struct gkyl_tok_geo_inp *inp = mapc2p_ctx;
    const struct gkyl_tok_geo_geo_inp *ginp = bmag_ctx;
    struct gkyl_tok_geo *geo = gkyl_tok_geo_new(inp);
    gkyl_tok_geo_advance(up, &nrange, dzc, NULL, geo, NULL, bmag_ctx, 
      up->mc2p_nodal_fd, up->mc2p_nodal, up->mc2p, 
      up->bmag, up->g_ij, up->jacobgeo, up->jacobgeo_inv, up->gij, up->b_i, up->cmag, up->jacobtot, 
      up->jacobtot_inv, up->bmag_inv, up->bmag_inv_sq, up->gxxj, up->gxyj, up->gyyj);

    gkyl_calc_bmag *bcalculator = gkyl_calc_bmag_new(up->basis, geo->rzbasis, geo->fbasis, up->grid, geo->rzgrid, geo->fgrid, geo, ginp, geo->psisep, false);

    struct gkyl_array* bmagrz = gkyl_array_new(GKYL_DOUBLE, geo->rzbasis->num_basis, geo->rzlocal_ext->volume);
    struct gkyl_array* bphirz = gkyl_array_new(GKYL_DOUBLE, geo->rzbasis->num_basis, geo->rzlocal_ext->volume);
    gkyl_calc_bmag_advance(bcalculator, up->range, up->range_ext, geo->rzlocal, geo->rzlocal_ext, geo->frange, geo->frange_ext, geo->psiRZ, geo->psibyrRZ, geo->psibyr2RZ, bphirz, bmagrz, up->bmag, geo->fpoldg, mc2p);


    // now calculate the metrics
    struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(up->basis, up->grid, false);
    gkyl_calc_metric_advance(mcalc, &nrange, mc2p_nodal_fd, dzc, g_ij, up->range);

    // calculate the derived geometric quantities
    gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(up->basis, up->grid, false);
    gkyl_calc_derived_geo_advance( jcalculator, range, g_ij, bmag, 
      jacobgeo, jacobgeo_inv, gij, b_i, cmag, jacobtot, jacobtot_inv, 
      bmag_inv, bmag_inv_sq, gxxj, gxyj, gyyj);
  }
  else{
    gkyl_gk_geometry_advance(up, &nrange, dzc, mapc2p_func, mapc2p_ctx, bmag_func, bmag_ctx, 
    mc2p_nodal_fd, mc2p_nodal, mc2p, 
    bmag, g_ij, jacobgeo, jacobgeo_inv, gij, b_i, cmag, jacobtot, 
    jacobtot_inv, bmag_inv, bmag_inv_sq, gxxj, gxyj, gyyj);
  }

  // Copy the host-side initialized geometry object to the device
  struct gkyl_array *mc2p_nodal_fd_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->grid->ndim*13, nrange.volume);
  struct gkyl_array *mc2p_nodal_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->grid->ndim, nrange.volume);
  struct gkyl_array *mc2p_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->grid->ndim, up->range_ext->volume);
  struct gkyl_array *bmag_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *g_ij_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, 6*up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *jacobgeo_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *jacobgeo_inv_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *gij_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, 6*up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *b_i_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *cmag_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *jacobtot_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *jacobtot_inv_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *bmag_inv_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *bmag_inv_sq_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *gxxj_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *gxyj_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  struct gkyl_array *gyyj_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);

  gkyl_array_copy(mc2p_nodal_fd_dev, mc2p_nodal_fd);
  gkyl_array_copy(mc2p_nodal_dev, mc2p_nodal);
  gkyl_array_copy(mc2p_dev, mc2p);
  gkyl_array_copy(bmag_dev, bmag);
  gkyl_array_copy(g_ij_dev, g_ij);
  gkyl_array_copy(jacobgeo_dev , jacobgeo);
  gkyl_array_copy(jacobgeo_inv_dev, jacobgeo_inv);
  gkyl_array_copy(gij_dev, gij);
  gkyl_array_copy(b_i_dev, b_i);
  gkyl_array_copy(cmag_dev, cmag);
  gkyl_array_copy(jacobtot_dev, jacobtot);
  gkyl_array_copy(jacobtot_inv_dev, jacobtot_inv);
  gkyl_array_copy(bmag_inv_dev, bmag_inv);
  gkyl_array_copy(bmag_inv_sq_dev, bmag_inv_sq);
  gkyl_array_copy(gxxj_dev, gxxj);
  gkyl_array_copy(gxyj_dev, gxyj);
  gkyl_array_copy(gyyj_dev, gyyj);

  gkyl_array_release(mc2p_nodal_fd);
  gkyl_array_release(mc2p_nodal);
  gkyl_array_release(mc2p);
  gkyl_array_release(bmag);
  gkyl_array_release(g_ij);
  gkyl_array_release(jacobgeo);
  gkyl_array_release(jacobgeo_inv);
  gkyl_array_release(gij);
  gkyl_array_release(b_i);
  gkyl_array_release(cmag);
  gkyl_array_release(jacobtot);
  gkyl_array_release(jacobtot_inv);
  gkyl_array_release(bmag_inv);
  gkyl_array_release(bmag_inv_sq);
  gkyl_array_release(gxxj);
  gkyl_array_release(gxyj);
  gkyl_array_release(gyyj);

  // this is for the memcpy below
  up->mc2p_nodal_fd = mc2p_nodal_fd_dev->on_dev;
  up->mc2p_nodal  = mc2p_nodal_dev ->on_dev;
  up->mc2p  = mc2p_dev->on_dev;
  up->bmag  = bmag_dev->on_dev;
  up->g_ij  = g_ij_dev->on_dev;
  up->jacobgeo  = jacobgeo_dev->on_dev;
  up->jacobgeo_inv = jacobgeo_inv->on_dev;
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
  up->mc2p_nodal_fd = mc2p_nodal_fd_dev;
  up->mc2p_nodal  = mc2p_nodal_dev ;
  up->mc2p  = mc2p_dev;
  up->bmag  = bmag_dev;
  up->g_ij  = g_ij_dev;
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
       

