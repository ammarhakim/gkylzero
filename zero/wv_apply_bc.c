#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_apply_bc.h>
#include <gkyl_wv_eqn.h>

struct gkyl_wv_apply_bc {
  struct gkyl_rect_grid grid;
  int dir; // direction to apply BC
  enum gkyl_edge_loc edge; // edge to apply BC
  int nghost[GKYL_MAX_DIM]; // number of ghost cells

  const struct gkyl_wv_eqn *eqn; // equation 
  const struct gkyl_wave_geom *geom; // geometry needed for BCs
  
  wv_bc_func_t bcfunc; // function pointer
  void *ctx; // context to pass to function
  
  struct gkyl_range ext_range, range; // ranges on grid
  struct gkyl_range skin, ghost; // skin and ghost ranges
};

gkyl_wv_apply_bc*
gkyl_wv_apply_bc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_wv_eqn *eqn, const struct gkyl_wave_geom *geom,
  int dir, enum gkyl_edge_loc edge, const int *nghost,
  wv_bc_func_t bcfunc, void *ctx)
{
  gkyl_wv_apply_bc *up = gkyl_malloc(sizeof(gkyl_wv_apply_bc));

  up->grid = *grid;
  up->dir = dir;
  up->edge = edge;
  gkyl_copy_int_arr(grid->ndim, nghost, up->nghost);

  up->eqn = gkyl_wv_eqn_acquire(eqn);
  up->geom = gkyl_wave_geom_acquire(geom);
  
  up->bcfunc = bcfunc;
  up->ctx = ctx;

  // compute range and extended range over grid ....
  gkyl_create_grid_ranges(grid, nghost, &up->ext_range, &up->range);
  // ... from these compute skin and ghost ranges
  gkyl_skin_ghost_ranges(&up->skin, &up->ghost, dir, edge, &up->ext_range, nghost);

  return up;
}

void
gkyl_wv_apply_bc_advance(const gkyl_wv_apply_bc *bc, double tm,
  const struct gkyl_range *update_rng, struct gkyl_array *out)
{
  enum gkyl_edge_loc edge = bc->edge;
  int dir = bc->dir, ndim = bc->grid.ndim, ncomp = out->ncomp;

  // return immediately if update region does not touch boundary
  if ( (edge == GKYL_LOWER_EDGE) && (update_rng->lower[dir] > bc->range.lower[dir]) )
    return;
  if ( (edge == GKYL_UPPER_EDGE) && (update_rng->upper[dir] < bc->range.upper[dir]) )
    return;

  // compute intersection for region to update
  struct gkyl_range up_range;
  gkyl_range_intersect(&up_range, update_rng, &bc->skin);

  int edge_idx = (edge == GKYL_LOWER_EDGE) ? bc->range.lower[dir] : bc->range.upper[dir];
  int fact = (edge == GKYL_LOWER_EDGE) ? -1 : 1;

  int gidx[GKYL_MAX_DIM]; // index into ghost cell
  
  // create iterator to walk over skin cells
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up_range);
  
  while (gkyl_range_iter_next(&iter)) {

    long sloc = gkyl_range_idx(update_rng, iter.idx);

    // compute linear index into appropriate ghost-cell
    gkyl_copy_int_arr(ndim, iter.idx, gidx);
    // this strange indexing ensures that the ghost cell index is
    // essentially "reflection" of the skin cell index; might not be
    // correct for all BCs but I am not sure how else to handle
    // multiple ghost-cell situations
    gidx[dir] = 2*edge_idx-gidx[dir]+fact;
    long gloc = gkyl_range_idx(update_rng, gidx);

    // apply boundary condition
    bc->bcfunc(tm, dir, ncomp,
      gkyl_array_fetch(out, sloc), gkyl_array_fetch(out, gloc), bc->ctx);
  }  
}

void
gkyl_wv_apply_bc_release(gkyl_wv_apply_bc* bc)
{
  gkyl_wv_eqn_release(bc->eqn);
  gkyl_wave_geom_release(bc->geom);
  gkyl_free(bc);
}
