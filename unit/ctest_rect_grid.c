#include <acutest.h>
#include <gkyl_rect_grid.h>

void test_grid_2d()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 20};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  TEST_CHECK( grid.ndim == 2 );
  for (int i=0; i<grid.ndim; ++i) {
    TEST_CHECK( grid.lower[i] == lower[i] );
    TEST_CHECK( grid.upper[i] == upper[i] );
    TEST_CHECK( grid.cells[i] == cells[i] );
    TEST_CHECK( grid.dx[i] == (upper[i]-lower[i])/cells[i] );
  }
  TEST_CHECK( grid.cellVolume == 0.075*0.2 );

  int idx[2];
  double xc[2];
  for (int i=0; i<grid.cells[0]; ++i)
    for (int j=0; j<grid.cells[1]; ++j) {
      idx[0] = i; idx[1] = j;
      gkyl_rect_grid_cell_center(&grid, idx, xc);

      TEST_CHECK( xc[0] == 1.0 + (i+0.5)*grid.dx[0] );
      TEST_CHECK( xc[1] == 1.0 + (j+0.5)*grid.dx[1] );
    }
}

TEST_LIST = {
  { "grid_2d", test_grid_2d },
  { NULL, NULL },
};
