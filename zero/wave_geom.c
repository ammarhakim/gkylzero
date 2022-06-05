#include <assert.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>

struct gkyl_wave_geom*
gkyl_wave_geom_new(const struct gkyl_rect_grid *grid, struct gkyl_range *range,
  evalf_t mapc2p, void *ctx)
{
  struct gkyl_wave_geom *wg = gkyl_malloc(sizeof(struct gkyl_wave_geom));

  wg->range = *range;
  wg->geom = gkyl_array_new(GKYL_USER, sizeof(struct gkyl_wave_cell_geom), range->volume);

  double xc[GKYL_MAX_CDIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {

    gkyl_rect_grid_cell_center(grid, iter.idx, xc);

    struct gkyl_wave_cell_geom *geo = gkyl_array_fetch(wg->geom, gkyl_range_idx(range, iter.idx));

    // compute geometry based on grid dimensions
    switch (grid->ndim) {
      case 1:
        calc_geom_1d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;

      case 2:
        calc_geom_2d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;

      case 3:
        calc_geom_3d(grid->dx, xc, mapc2p ? mapc2p : nomapc2p, ctx, geo);
        break;
    };
  }

  wg->ref_count = gkyl_ref_count_init(wave_geom_free);

  wg->on_dev = wg;

  return wg;
}

void
wave_geom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wave_geom *wg = container_of(ref, struct gkyl_wave_geom, ref_count);
  gkyl_array_release(wg->geom);
  if (wg != wg->on_dev) { // TODO: proper gkyl way to detect device copy
    gkyl_cu_free(wg->on_dev); // wg->on_dev->geom was released with wg->geom
  }
  gkyl_free(wg);
}

struct gkyl_wave_geom *
gkyl_wave_geom_acquire(const struct gkyl_wave_geom* wg)
{
  gkyl_ref_count_inc(&wg->ref_count);
  return (struct gkyl_wave_geom*) wg;
}

void
gkyl_wave_geom_release(const struct gkyl_wave_geom *wg)
{
  gkyl_ref_count_dec(&wg->ref_count);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_wave_geom*
gkyl_wave_geom_cu_dev_new(const struct gkyl_rect_grid *grid,
                          struct gkyl_range *range,
                          evalf_t mapc2p,
                          void *ctx)
{
  assert(false);
  return 0;
}

#endif
