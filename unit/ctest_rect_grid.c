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

void test_find_cell_1x(){
  double lower[] = {0.0}, upper[] = {5.0};
  int cells[] = {5};
  double point[] = {2.5};
  bool pickLower = false;
  int idx = 1;
  const int *ptr = &idx;
  const int *knownIdx[1] = {NULL};
  const int *knownIdx2[1] = {ptr};
  int cellIdx[]={0};
  int correctIdx[1]={3};
  struct gkyl_rect_grid grid;
  
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx, cellIdx);
  for(int i=0;i<1;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  point[0]=2.0+1e-15;
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx, cellIdx);
  for(int i=0;i<1;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  pickLower=true;
  correctIdx[0]=2;
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx, cellIdx);
  for(int i=0;i<1;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  correctIdx[0] = idx;
  point[0] = 0.2;
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx2, cellIdx);  
  for(int i=0;i<1;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

}

void test_find_cell_2x(){
  double lower[] = {0.0, -10.0}, upper[] = {5.0, 10.0};
  int cells[] = {5, 20};
  double point[] = {2.5, 1.3};
  bool pickLower = false;
  int idx = 18;
  const int *ptr = &idx;
  const int *knownIdx[2] = {NULL,NULL};
  const int *knownIdx2[2] = {NULL,ptr};
  int cellIdx[]={0,0};
  int correctIdx[2]={3,12};
  struct gkyl_rect_grid grid;
  
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx, cellIdx);
  for(int i=0;i<2;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  point[0]=2.0;
  point[1]=1.0;
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx, cellIdx);
  for(int i=0;i<2;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  pickLower=true;
  correctIdx[0]=2;
  correctIdx[1]=11;
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx, cellIdx);
  for(int i=0;i<2;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  pickLower=false;
  correctIdx[0]=3;
  correctIdx[1]=idx;
  point[1]=7.5;
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx2, cellIdx);  
  for(int i=0;i<2;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

}

void test_find_cell_3x(){
  double lower[] = {0.0, -10.0, 1.3}, upper[] = {5.0, 10.0, 2.5};
  int cells[] = {5, 20, 100};
  double point[] = {2.5, 1.3, 1.4};
  bool pickLower = false;
  int idx = 18, idx2 = 35;
  const int *ptr = &idx, *ptr2 = &idx2;
  const int *knownIdx[] = {NULL,NULL,NULL};
  const int *knownIdx2[] = {NULL,ptr,ptr2};
  int cellIdx[]={0,0,0};
  int correctIdx[]={3,12,9};
  struct gkyl_rect_grid grid;
  
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx, cellIdx);
  for(int i=0;i<3;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  point[0]=2.0;
  point[1]=1.0+1e-15;
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx, cellIdx);
  for(int i=0;i<3;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  pickLower=true;
  correctIdx[0]=2;
  correctIdx[1]=11;
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx, cellIdx);
  for(int i=0;i<3;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

  pickLower=false;
  correctIdx[0]=3;
  correctIdx[1]=idx;
  correctIdx[2]=idx2;
  point[1]=7.5;
  point[2]=1.71;
  gkyl_rect_grid_find_cell(&grid, point, pickLower, knownIdx2, cellIdx);  
  for(int i=0;i<3;i++){
    TEST_CHECK( cellIdx[i]==correctIdx[i]);
  }

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
  { "grid_find_cell_1x", test_find_cell_1x },
  { "grid_find_cell_2x", test_find_cell_2x },
  { "grid_find_cell_3x", test_find_cell_3x },
  { "grid_io", test_grid_io },
#ifdef GKYL_HAVE_CUDA
  { "cu_grid_2d", test_cu_grid_2d },
#endif  
  { NULL, NULL },
};
