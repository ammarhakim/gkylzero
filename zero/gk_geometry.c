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
gkyl_gk_geometry_is_cu_dev(const struct gkyl_gk_geometry* up)
{
  return GKYL_IS_CU_ALLOC(up->flags);
}

void
gkyl_gk_geometry_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_gk_geometry *up = container_of(ref, struct gkyl_gk_geometry, ref_count);
  gkyl_array_release(up->mc2p_nodal_fd);
  gkyl_array_release(up->mc2p_nodal);
  gkyl_array_release(up->mc2p);

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

struct gkyl_gk_geometry*
gkyl_gk_geometry_new(const struct gkyl_rect_grid* grid, const struct gkyl_range *range, const struct gkyl_range* range_ext, 
  const struct gkyl_basis* basis, evalf_t mapc2p_func, void* mapc2p_ctx, evalf_t bmag_func, void* bmag_ctx, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_gk_geometry_new_cu_dev_new(grid, range, range_ext, basis, mapc2p_func, mapc2p_ctx, bmag_func, bmag_ctx);
  } 
#endif 

  struct gkyl_gk_geometry *up = gkyl_malloc(sizeof(*up));
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
  up->mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim*13, nrange.volume);
  up->mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim, nrange.volume);
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid->ndim, up->range_ext->volume);

  // bmag, metrics and derived geo quantities
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis->num_basis, up->range_ext->volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->jacobgeo_inv= gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis->num_basis, up->range_ext->volume);
  up->b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis->num_basis, up->range_ext->volume);
  up->cmag = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->gxxj= gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->gxyj= gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);
  up->gyyj= gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->range_ext->volume);

  gkyl_gk_geometry_advance(up, &nrange, dzc, mapc2p_func, mapc2p_ctx, bmag_func, bmag_ctx);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself
                   
  return up;
}

struct gkyl_gk_geometry*
gkyl_gk_geometry_acquire(const struct gkyl_gk_geometry* up)
{
  gkyl_ref_count_inc(&up->ref_count);
  return (struct gkyl_gk_geometry*) up;
}

void
gkyl_gk_geometry_release(const struct gkyl_gk_geometry *up)
{
  gkyl_ref_count_dec(&up->ref_count);
}



