#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>

#include<gkyl_tok_geo.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_bmag.h>


// write out nodal coordinates 
static void
write_nodal_coordinates(const char *nm, struct gkyl_range *nrange,
  struct gkyl_array *nodes)
{
  double lower[3] = { 0.0, 0.0, 0.0 };
  double upper[3] = { 1.0, 1.0, 1.0 };
  int cells[3];
  for (int i=0; i<nrange->ndim; ++i)
    cells[i] = gkyl_range_shape(nrange, i);
  
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  gkyl_grid_sub_array_write(&grid, nrange, nodes, nm);
}

struct gk_geometry*
gkyl_gk_geometry_tok_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, void* tok_rz_ctx, void* tok_comp_ctx, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_gk_geometry_tok_cu_dev_new(grid, range, range_ext, basis, tok_rz_ctx, tok_comp_ctx);
  } 
#endif 

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
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

  const struct gkyl_tok_geo_efit_inp *inp = tok_rz_ctx;
  struct gkyl_tok_geo_grid_inp *ginp = tok_comp_ctx;
  ginp->cgrid = up->grid;
  ginp->cbasis = up->basis;
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(inp);
  // calculate mapc2p
  gkyl_tok_geo_calc(up, &nrange, dzc, NULL, geo, NULL, tok_comp_ctx, 
    mc2p_nodal_fd, mc2p_nodal, mc2p);
  // calculate bmag
  gkyl_calc_bmag *bcalculator = gkyl_calc_bmag_new(&up->basis, &geo->rzbasis, &geo->fbasis, &up->grid, &geo->rzgrid, &geo->fgrid, geo, ginp, geo->psisep, false);
  gkyl_calc_bmag_advance(bcalculator, &up->range, &up->range_ext, &geo->rzlocal, &geo->rzlocal_ext, &geo->frange, &geo->frange_ext, geo->psiRZ, geo->psibyrRZ, geo->psibyr2RZ, up->bmag, geo->fpoldg, mc2p);
  // now calculate the metrics
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid, false);
  gkyl_calc_metric_advance(mcalc, &nrange, mc2p_nodal_fd, dzc, up->g_ij, up->dxdz, &up->range);
  // calculate the derived geometric quantities
  gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(&up->basis, &up->grid, false);
  gkyl_calc_derived_geo_advance( jcalculator, &up->range, up->g_ij, up->bmag, 
    up->jacobgeo, up->jacobgeo_inv, up->gij, up->b_i, up->cmag, up->jacobtot, up->jacobtot_inv, 
    up->bmag_inv, up->bmag_inv_sq, up->gxxj, up->gxyj, up->gyyj);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself
                   
  gkyl_grid_sub_array_write(&up->grid, &up->range, mc2p, "mapc2p.gkyl");
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

