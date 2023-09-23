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
}

void test_find_cell(){
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {5.0, 5.0, 10.0};
  double lower2[] = {0.0, 1.0, 3.2, 7.1}, upper2[] = {5.1, 5.0, 22.0, 8.2};
  int cells[] = {5, 5, 10};
  int cells2[] = {10, 11, 20, 3};
  double point[] = {2.5, 0.1, 4.1};
  double point2[] = {2.01, 3.1, 4.1, 8.05};
  bool pickLower = false;
  int idx = 1;
  const int *ptr = &idx;
  const int *knownIdx[3] = {NULL, ptr, NULL};
  const int *knownIdx2[4] = {NULL, NULL, NULL, NULL};
  int *cellIdx;
  int correctIdx[3]={3,idx,5};
  int correctIdx2[4]={4,6,1,3};//Got 6,9,2,6
  struct gkyl_rect_grid grid;
  struct gkyl_rect_grid grid2;
  
  printf("\nStart test_find_cell\n");
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  cellIdx = gkyl_find_cell(&grid, point, pickLower, knownIdx);
  for(int i=0;i<3;i++){
    printf("i=%i, cellIdx[i]=%i\n",i,cellIdx[i]);
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  gkyl_rect_grid_init(&grid2, 4, lower2, upper2, cells2);
  cellIdx = gkyl_find_cell(&grid2, point2, pickLower, knownIdx2);
  for(int i=0;i<4;i++){
    printf("i=%i, cellIdx[i]=%i\n",i,cellIdx[i]);
    TEST_CHECK( cellIdx[i]==correctIdx2[i]);
  }
  printf("End test_find_cell\n");
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
  { "grid_find_cell", test_find_cell },
  { "grid_io", test_grid_io },
#ifdef GKYL_HAVE_CUDA
  { "cu_grid_2d", test_cu_grid_2d },
#endif  
  { NULL, NULL },
};
