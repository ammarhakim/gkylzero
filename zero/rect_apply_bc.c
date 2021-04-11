#include <gkyl_rect_apply_bc.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_range.h>

struct gkyl_rect_apply_bc {
    struct gkyl_rect_grid grid;
    int dir; // direction to apply BC
    int edge; // edge to apply BC
    rect_bc_func_t bcfunc; // function pointer
    void *ctx; // context to pass to function
};

gkyl_rect_apply_bc*
gkyl_rect_apply_bc_new(const struct gkyl_rect_grid *grid,
  int dir, int edge, rect_bc_func_t bcfunc, void *ctx)
{
  gkyl_rect_apply_bc *up = gkyl_malloc(sizeof(gkyl_rect_apply_bc));

  up->grid = *grid;
  up->dir = dir;
  up->edge = edge;
  up->bcfunc = bcfunc;
  up->ctx = ctx;

  return up;
}

void
gkyl_rect_apply_bc_advance(const gkyl_rect_apply_bc *bc, double tm, struct gkyl_array *out)
{
}

void
gkyl_rect_apply_bc_release(gkyl_rect_apply_bc* bc)
{
  gkyl_free(bc);
}
