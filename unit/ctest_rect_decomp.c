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
  TEST_CHECK( ext_range.lower[0] == 0 );
  TEST_CHECK( ext_range.upper[0] == 21 );

  TEST_CHECK( range.ndim == 1 );
  TEST_CHECK( range.volume ==  20 );
  TEST_CHECK( range.lower[0] == 1 );
  TEST_CHECK( range.upper[0] == 20 );

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
  TEST_CHECK( ext_range.lower[0] == 0 );
  TEST_CHECK( ext_range.upper[0] == 21 );
  TEST_CHECK( ext_range.lower[1] == 1 );
  TEST_CHECK( ext_range.upper[1] == 40 );  

  TEST_CHECK( range.ndim == 2 );
  TEST_CHECK( range.volume ==  20*40 );
  TEST_CHECK( range.lower[0] == 1 );
  TEST_CHECK( range.upper[0] == 20 );
  TEST_CHECK( range.lower[1] == 1 );
  TEST_CHECK( range.upper[1] == 40 );  

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
  TEST_CHECK( ext_range.lower[0] == 0 );
  TEST_CHECK( ext_range.upper[0] == 21 );
  TEST_CHECK( ext_range.lower[1] == 1 );
  TEST_CHECK( ext_range.upper[1] == 40 );
  TEST_CHECK( ext_range.lower[2] == -1 );
  TEST_CHECK( ext_range.upper[2] == 12 ); 

  TEST_CHECK( range.ndim == 3 );
  TEST_CHECK( range.volume ==  20*40*10 );
  TEST_CHECK( range.lower[0] == 1 );
  TEST_CHECK( range.upper[0] == 20 );
  TEST_CHECK( range.lower[1] == 1 );
  TEST_CHECK( range.upper[1] == 40 );
  TEST_CHECK( range.lower[2] == 1 );
  TEST_CHECK( range.upper[2] == 10 );  

  TEST_CHECK( gkyl_range_is_sub_range(&range) == 1 );
}

void
test_ranges_from_range_2d(void)
{
  struct gkyl_range inlocal;
  gkyl_range_init(&inlocal, 2, (int[]) { 1, 2 }, (int[]) { 10, 20 });

  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&inlocal, (int[]) { 2, 1 }, &local_ext, &local);

  TEST_CHECK( local.ndim == inlocal.ndim );
  TEST_CHECK( local_ext.ndim == inlocal.ndim );
  
  for (int i=0; i<inlocal.ndim; ++i) {
    TEST_CHECK( local.lower[i] == inlocal.lower[i] );
    TEST_CHECK( local.upper[i] == inlocal.upper[i] );
  }

  TEST_CHECK( gkyl_range_is_sub_range(&local) == 1 );
}

void
test_ranges_from_range_3d(void)
{
  struct gkyl_range inlocal;
  gkyl_range_init(&inlocal, 3, (int[]) { 1, 2, 3 }, (int[]) { 10, 20, 30 });

  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&inlocal, (int[]) { 2, 1, 0 }, &local_ext, &local);

  TEST_CHECK( local.ndim == inlocal.ndim );
  TEST_CHECK( local_ext.ndim == inlocal.ndim );
  
  for (int i=0; i<inlocal.ndim; ++i) {
    TEST_CHECK( local.lower[i] == inlocal.lower[i] );
    TEST_CHECK( local.upper[i] == inlocal.upper[i] );
  }

  TEST_CHECK( gkyl_range_is_sub_range(&local) == 1 );
}

void
test_rect_decomp_2d(void)
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 100, 100 });
  
  int cuts[] = { 5, 5 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(2, cuts, &range);

  TEST_CHECK( decomp->ndim == 2 );
  TEST_CHECK( decomp->ndecomp == 25 );

  gkyl_rect_decomp_release(decomp);
}

TEST_LIST = {
  { "ranges_1d", test_ranges_1d },
  { "ranges_2d", test_ranges_2d },
  { "ranges_3d", test_ranges_3d },

  { "ranges_from_range_2d", test_ranges_from_range_2d },
  { "ranges_from_range_3d", test_ranges_from_range_3d },

  { "rect_decomp_2d", test_rect_decomp_2d },
  
  { NULL, NULL },
};
