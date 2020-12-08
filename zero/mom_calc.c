#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc.h>

struct gkyl_mom_calc {
    struct gkyl_rect_grid grid; // grid
    struct gkyl_mom_type *momt; // moment-type
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
gkyl_mom_calc_advance(const gkyl_mom_calc* calc, const struct gkyl_range *update_rng,
  const struct gkyl_array *fin, struct gkyl_array *mout)
{
  double xc[GKYL_MAX_DIM];
  int phase_numBasis = calc->momt->num_phase;
  int conf_numBasis = calc->momt->num_config;
  int num_mom = calc->momt->num_mom;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_rng);

  gkyl_array_clear(mout, 0.0);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&calc->grid, iter.idx, xc);

    long didx = gkyl_range_idx(update_rng, iter.idx);
    // long midx = gkyl_range_index(mom_range, iter->idx);

    // // call function to compute moment
    // gkyl_mom_type_calc(
    //   calc->momt, xc, dx, iter->idx,
    //   gkyl_array_fetch(fin, didx*phase_numBasis),
    //   gkyl_array_fetch(mout, midx*conf_numBasis*num_mom)
    // );
  }
}

void gkyl_mom_calc_release(gkyl_mom_calc* up)
{
  gkyl_mom_type_release(up->momt);
  gkyl_free(up);
}
