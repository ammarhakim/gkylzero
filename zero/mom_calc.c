#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc.h>

struct gkyl_mom_calc {
    struct gkyl_rect_grid grid;
    struct gkyl_mom_type *momt;
};

gkyl_mom_calc*
gkyl_mom_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt)
{
  gkyl_mom_calc *up = gkyl_malloc(sizeof(gkyl_mom_calc));
  up->grid = *grid;
  up->momt = gkyl_mom_type_aquire(momt);
  return up;
}

void
gkyl_mom_calc_advance(const gkyl_mom_calc* calc,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fin, struct gkyl_array *mout)
{
  double xc[GKYL_MAX_DIM];

  gkyl_array_clear(mout, 0.0);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_rng);
  while (gkyl_range_iter_next(&iter)) {
    
    gkyl_rect_grid_cell_center(&calc->grid, iter.idx, xc);

    long fidx = gkyl_range_idx(phase_rng, iter.idx);
    long midx = gkyl_range_idx(conf_rng, iter.idx);

    gkyl_mom_type_calc(calc->momt, xc, calc->grid.dx, iter.idx,
      gkyl_array_cfetch(fin, fidx), gkyl_array_fetch(mout, midx)
    );
  }
}

void gkyl_mom_calc_release(gkyl_mom_calc* up)
{
  gkyl_mom_type_release(up->momt);
  gkyl_free(up);
}
