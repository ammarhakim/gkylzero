extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wave_geom_priv.h>
}

__global__ static void
wave_geom_set_cu_dev_ptrs(struct gkyl_wave_geom *wg_dev,
                          struct gkyl_array *geom_dev)
{
  wg_dev->geom = geom_dev;
}

// CPU interface to create and track a GPU object
struct gkyl_wave_geom*
gkyl_wave_geom_cu_dev_new(const struct gkyl_rect_grid *grid,
                          struct gkyl_range *range,
                          evalf_t mapc2p,
                          void *ctx)
{
  // STEP: CREATE HOST OBJECT
  struct gkyl_wave_geom *wg =
    (struct gkyl_wave_geom*) gkyl_malloc(sizeof(struct gkyl_wave_geom));

  // STEP: SET HOST OR COMMON HOST/DEVICE DATA IN HOST OBJECT
  wg->range = *range;

  wg->geom = gkyl_array_new(
      GKYL_USER, sizeof(struct gkyl_wave_cell_geom), range->volume);
  double xc[GKYL_MAX_CDIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(grid, iter.idx, xc);
    struct gkyl_wave_cell_geom *geo =
      (struct gkyl_wave_cell_geom*) gkyl_array_fetch(
          wg->geom, gkyl_range_idx(range, iter.idx));
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

  // STEP: COPY HOST OBJECT TO DEVICE OBJECT
  struct gkyl_wave_geom *wg_dev =
    (struct gkyl_wave_geom*) gkyl_cu_malloc(sizeof(struct gkyl_wave_geom));
  gkyl_cu_memcpy(wg_dev, wg, sizeof(struct gkyl_wave_geom), GKYL_CU_MEMCPY_H2D);

  // STEP: SET OTHER DEVICE DATA IN DEVICE OBJECT
  struct gkyl_array *geom_dev = gkyl_array_cu_dev_new(
      GKYL_USER, sizeof(struct gkyl_wave_cell_geom), range->volume);
  gkyl_array_copy(geom_dev, wg->geom);
  wave_geom_set_cu_dev_ptrs<<<1, 1>>>(wg_dev, geom_dev->on_dev);

  // STEP: KEEP POINTER TO THE DEVICE OBJECT
  wg->on_dev = wg_dev;

  return wg;
}
