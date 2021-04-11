#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_range.h>
#include <gkyl_rect_apply_bc.h>
#include <gkyl_rect_decomp.h>

struct gkyl_rect_apply_bc {
    struct gkyl_rect_grid grid;
    int dir; // direction to apply BC
    enum gkyl_edge_loc edge; // edge to apply BC
    int nghost[GKYL_MAX_DIM]; // number of ghost cells
    rect_bc_func_t bcfunc; // function pointer
    void *ctx; // context to pass to function

    struct gkyl_range ext_range, range; // ranges on grid
    struct gkyl_range skin, ghost; // skin and ghost ranges
};

// Increment an int vector by fact*del[d] in each direction d.
static void
incr_int_array(int ndim, int fact, const int *restrict del,
  const int *restrict inp, int *restrict out)
{
  for (int i=0; i<ndim; ++i)
    out[i] = inp[i] + fact*del[i];
}

// Create ghost and skin sub-ranges given parent (extended
// range). This code is somewhat convoluted as the skin and ghost
// ranges need to be sub-ranges of the extended range on the grid and
// not include corners. I am nor sure how to handle corners on
// physical boundaries. Also, perhaps this code could be simplified.
static void
skin_ghost_ranges_init(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *parent, const int *nghost)
{
  int ndim = parent->ndim, lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];

  if (edge == GKYL_LOWER_EDGE) {

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    up[dir] = lo[dir]+nghost[dir]-1;
    gkyl_sub_range_init(skin, parent, lo, up);

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    lo[dir] = lo[dir]-nghost[dir];
    up[dir] = lo[dir]+nghost[dir]-1;
    gkyl_sub_range_init(ghost, parent, lo, up);
  }
  else {

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    lo[dir] = up[dir]-nghost[dir]+1;
    gkyl_sub_range_init(skin, parent, lo, up);

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    up[dir] = up[dir]+nghost[dir]+1;
    lo[dir] = up[dir]-nghost[dir];
    gkyl_sub_range_init(ghost, parent, lo, up);
  }
}

gkyl_rect_apply_bc*
gkyl_rect_apply_bc_new(const struct gkyl_rect_grid *grid,
  int dir, enum gkyl_edge_loc edge, const int *nghost, rect_bc_func_t bcfunc, void *ctx)
{
  gkyl_rect_apply_bc *up = gkyl_malloc(sizeof(gkyl_rect_apply_bc));

  up->grid = *grid;
  up->dir = dir;
  up->edge = edge;
  for (int d=0; d<grid->ndim; ++d) up->nghost[d] = nghost[d];
  up->bcfunc = bcfunc;
  up->ctx = ctx;

  // compute range and extended range over grid; and skin and ghost
  // ranges for use in the BC
  gkyl_create_grid_ranges(grid, nghost, &up->ext_range, &up->range);
  skin_ghost_ranges_init(&up->skin, &up->ghost, dir, edge, &up->ext_range, nghost);

  gkyl_print_range(&up->ghost, "ghost", stdout);
  gkyl_print_range(&up->skin, "skin", stdout);

  return up;
}

void
gkyl_rect_apply_bc_advance(const gkyl_rect_apply_bc *bc, double tm,
  const struct gkyl_range *update_rng, struct gkyl_array *out)
{
}

void
gkyl_rect_apply_bc_release(gkyl_rect_apply_bc* bc)
{
  gkyl_free(bc);
}
