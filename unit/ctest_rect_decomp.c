#include <acutest.h>
#include <gkyl_rect_decomp.h>

void test_ranges_1d()
{
  double lower[] = { 1.0 }, upper[] = {2.5 };
  int cells[] = { 20 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  int nghost[] = { 1 };
  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  TEST_CHECK( ext_range.ndim == 1 );
  TEST_CHECK( ext_range.volume ==  22 );
  TEST_CHECK( ext_range.lower[0] == -1 );
  TEST_CHECK( ext_range.upper[0] == 20 );

  TEST_CHECK( range.ndim == 1 );
  TEST_CHECK( range.volume ==  20 );
  TEST_CHECK( range.lower[0] == 0 );
  TEST_CHECK( range.upper[0] == 19 );

  TEST_CHECK( gkyl_range_is_sub_range(&range) == 1 );
}

void test_ranges_2d()
{
  double lower[] = { 1.0, 1.0 }, upper[] = { 2.5, 5.0 };
  int cells[] = { 20, 40 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nghost[] = { 1, 0 };
  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  TEST_CHECK( ext_range.ndim == 2 );
  TEST_CHECK( ext_range.volume ==  22*40 );
  TEST_CHECK( ext_range.lower[0] == -1 );
  TEST_CHECK( ext_range.upper[0] == 20 );
  TEST_CHECK( ext_range.lower[1] == 0 );
  TEST_CHECK( ext_range.upper[1] == 39 );  

  TEST_CHECK( range.ndim == 2 );
  TEST_CHECK( range.volume ==  20*40 );
  TEST_CHECK( range.lower[0] == 0 );
  TEST_CHECK( range.upper[0] == 19 );
  TEST_CHECK( range.lower[1] == 0 );
  TEST_CHECK( range.upper[1] == 39 );  

  TEST_CHECK( gkyl_range_is_sub_range(&range) == 1 );  
}

void test_ranges_3d()
{
  double lower[] = { 1.0, 1.0, 1.0 }, upper[] = { 2.5, 5.0, 2.0 };
  int cells[] = { 20, 40, 10 };
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  int nghost[] = { 1, 0, 2 };
  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  TEST_CHECK( ext_range.ndim == 3 );
  TEST_CHECK( ext_range.volume ==  22*40*14 );
  TEST_CHECK( ext_range.lower[0] == -1 );
  TEST_CHECK( ext_range.upper[0] == 20 );
  TEST_CHECK( ext_range.lower[1] == 0 );
  TEST_CHECK( ext_range.upper[1] == 39 );
  TEST_CHECK( ext_range.lower[2] == -2 );
  TEST_CHECK( ext_range.upper[2] == 11 );  

  TEST_CHECK( range.ndim == 3 );
  TEST_CHECK( range.volume ==  20*40*10 );
  TEST_CHECK( range.lower[0] == 0 );
  TEST_CHECK( range.upper[0] == 19 );
  TEST_CHECK( range.lower[1] == 0 );
  TEST_CHECK( range.upper[1] == 39 );
  TEST_CHECK( range.lower[2] == 0 );
  TEST_CHECK( range.upper[2] == 9 );  

  TEST_CHECK( gkyl_range_is_sub_range(&range) == 1 );    
}

TEST_LIST = {
  { "ranges_1d", test_ranges_1d },
  { "ranges_2d", test_ranges_2d },
  { "ranges_3d", test_ranges_3d },
  { NULL, NULL },
};
