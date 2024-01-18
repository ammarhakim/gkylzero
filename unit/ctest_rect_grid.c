#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

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

  double xc[2];  
  int x_ext[2], y_ext[2];
  
  gkyl_rect_grid_extents(&grid, 0, x_ext);
  gkyl_rect_grid_extents(&grid, 1, y_ext);  

  for (int i=x_ext[0]; i<=x_ext[1]; ++i)
    for (int j=y_ext[0]; j<=y_ext[1]; ++j) {
      gkyl_rect_grid_cell_center(&grid, (int[2]) { i, j }, xc);

      TEST_CHECK( xc[0] == 1.0 + (i-0.5)*grid.dx[0] );
      TEST_CHECK( xc[1] == 1.0 + (j-0.5)*grid.dx[1] );
    }

  double xn[2] = { 0.0, 0.0 };
  int idx[2];

  xn[0] = grid.lower[0] + 0.5*grid.dx[0];
  xn[1] = grid.lower[1] + 0.5*grid.dx[1];
  gkyl_rect_grid_coord_idx(&grid, xn, idx);
  TEST_CHECK( (idx[0] == 1) && (idx[1] == 1) );

  xn[0] = grid.lower[0] + 1.5*grid.dx[0];
  xn[1] = grid.lower[1] + 1.5*grid.dx[1];
  gkyl_rect_grid_coord_idx(&grid, xn, idx);
  TEST_CHECK( (idx[0] == 2) && (idx[1] == 2) );

  TEST_CHECK( gkyl_rect_grid_cmp(&grid, &grid) == true );

  double lower2[] = {1.0, 0.5}, upper2[] = {2.5, 5.0};
  int cells2[] = {20, 19};
  struct gkyl_rect_grid grid2;
  gkyl_rect_grid_init(&grid2, 2, lower2, upper2, cells2);

  TEST_CHECK( gkyl_rect_grid_cmp(&grid, &grid2) == false );
}

void test_grid_io()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 20};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  FILE *fp = 0;
  with_file (fp, "ctest_rect_grid.dat", "w")
    gkyl_rect_grid_write(&grid, fp);

  struct gkyl_rect_grid grid2;
  with_file (fp, "ctest_rect_grid.dat", "r")
    gkyl_rect_grid_read(&grid2, fp);

  TEST_CHECK( grid.ndim == grid2.ndim );
  for (int d=0; d<grid.ndim; ++d) {
    TEST_CHECK( grid.lower[d] == grid2.lower[d] );
    TEST_CHECK( grid.upper[d] == grid2.upper[d] );
    TEST_CHECK( grid.cells[d] == grid2.cells[d] );
    TEST_CHECK( grid.dx[d] == grid2.dx[d] );
  }
  TEST_CHECK( grid.cellVolume == grid2.cellVolume );
}

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

int cu_rect_grid_test(const struct gkyl_rect_grid grid);

void test_cu_grid_2d()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 20};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nfail = cu_rect_grid_test(grid);
  TEST_CHECK( nfail == 0 );
}
#endif

TEST_LIST = {
  { "grid_2d", test_grid_2d },
  { "grid_io", test_grid_io },
#ifdef GKYL_HAVE_CUDA
  { "cu_grid_2d", test_cu_grid_2d },
#endif  
  { NULL, NULL },
};
