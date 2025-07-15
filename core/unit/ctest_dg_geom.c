#include <acutest.h>

#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_geom.h>

#include <math.h>

static void
test_dg_1()
{
  int ndim = 1;
  double lower[] = {0.0}, upper[] = {1.0};
  int cells[] = {10};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // create range
  int nghost[3] = { 1, 1, 1 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_dg_geom *dgg
    = gkyl_dg_geom_new( &(const struct gkyl_dg_geom_inp) {
        .grid = &grid,
        .range = &ext_range,
        .nquad = 2,
      }
    );

  gkyl_dg_geom_release(dgg);
}

TEST_LIST = {
  { "test_dg_1", test_dg_1 },
  { NULL, NULL },
};
