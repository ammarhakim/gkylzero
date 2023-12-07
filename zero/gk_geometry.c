#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_math.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
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
gkyl_gk_geometry_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx, bool tokamak, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_gk_geometry_cu_dev_new(grid, range, range_ext, basis, mapc2p_func, mapc2p_ctx, bmag_func, bmag_ctx, tokamak);
  } 
#endif 

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->tokamak = tokamak;
  up->basis = *basis;
  up->range = *range;
  up->range_ext = *range_ext;
  up->grid = *grid;
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

  gkyl_range_init_from_shape(&nrange, up->grid.ndim, nodes);
  struct gkyl_array* mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*13, nrange.volume);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange.volume);
  struct gkyl_array* mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->range_ext.volume);

  // bmag, metrics and derived geo quantities
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->range_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->range_ext.volume);
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

  if (up->tokamak){
    const struct gkyl_tok_geo_inp *inp = mapc2p_ctx;
    struct gkyl_tok_geo_geo_inp *ginp = bmag_ctx;
    ginp->cgrid = up->grid;
    ginp->cbasis = up->basis;
    struct gkyl_tok_geo *geo = gkyl_tok_geo_new(inp);
    // calculate mapc2p
    gkyl_tok_geo_advance(up, &nrange, dzc, NULL, geo, NULL, bmag_ctx, 
      mc2p_nodal_fd, mc2p_nodal, mc2p, 
      up->bmag, up->g_ij, up->jacobgeo, up->jacobgeo_inv, up->gij, up->b_i, up->cmag, up->jacobtot, 
      up->jacobtot_inv, up->bmag_inv, up->bmag_inv_sq, up->gxxj, up->gxyj, up->gyyj);
    // calculate bmag
    gkyl_calc_bmag *bcalculator = gkyl_calc_bmag_new(&up->basis, &geo->rzbasis, &geo->fbasis, &up->grid, &geo->rzgrid, &geo->fgrid, geo, ginp, geo->psisep, false);
    gkyl_calc_bmag_advance(bcalculator, &up->range, &up->range_ext, &geo->rzlocal, &geo->rzlocal_ext, &geo->frange, &geo->frange_ext, geo->psiRZ, geo->psibyrRZ, geo->psibyr2RZ, up->bmag, geo->fpoldg, mc2p);
    // now calculate the metrics
    struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid, false);
    gkyl_calc_metric_advance(mcalc, &nrange, mc2p_nodal_fd, dzc, up->g_ij, &up->range);
    // calculate the derived geometric quantities
    gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(&up->basis, &up->grid, false);
    gkyl_calc_derived_geo_advance( jcalculator, &up->range, up->g_ij, up->bmag, 
      up->jacobgeo, up->jacobgeo_inv, up->gij, up->b_i, up->cmag, up->jacobtot, up->jacobtot_inv, 
      up->bmag_inv, up->bmag_inv_sq, up->gxxj, up->gxyj, up->gyyj);
  }
  else{
    gkyl_gk_geometry_advance(up, &nrange, dzc, mapc2p_func, mapc2p_ctx, bmag_func, bmag_ctx, 
    mc2p_nodal_fd, mc2p_nodal, mc2p);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself
                   
  gkyl_grid_sub_array_write(&up->grid, &up->range, mc2p, "mapc2p.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->bmag, "bmag.gkyl");
  gkyl_grid_sub_array_write(&up->grid, &up->range, up->g_ij, "g_ij.gkyl");
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



