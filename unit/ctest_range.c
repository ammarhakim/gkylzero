#include <acutest.h>
#include <gkyl_range.h>

void test_range_1()
{
  int lower[2] = {1, 1}, upper[2] = {10, 20};
  struct gkyl_range range;
  gkyl_range_init(&range, 2, lower, upper);

  TEST_CHECK( range.ndim == 2);
  TEST_CHECK( range.volume == 200);

  for (unsigned i=0; i<2; ++i) {
    TEST_CHECK( range.lower[i] == lower[i]);
    TEST_CHECK( range.upper[i] == upper[i]);
  }
}

void test_range_shape()
{
  int shape[2] = {25, 50};
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, 2, shape);

  TEST_CHECK( range.ndim == 2);
  TEST_CHECK( range.volume == 25*50);

  for (unsigned i=0; i<2; ++i) {
    TEST_CHECK( range.lower[i] == 0);
    TEST_CHECK( range.upper[i] == shape[i]-1);
  }
}

void test_shorten()
{
  int lower[3] = {1, 1, 1}, upper[3] = {20, 30, 20};
  struct gkyl_range range, shortr, shortr2;

  gkyl_range_init(&range, 3, lower, upper);
  gkyl_range_shorten(&shortr, &range, 1, 1);
  gkyl_range_shorten(&shortr2, &range, 1, 3);

  // shortr
  TEST_CHECK( shortr.lower[0] == range.lower[0] );
  TEST_CHECK( shortr.upper[0] == range.upper[0] );

  TEST_CHECK( shortr.lower[1] == range.lower[1] );
  TEST_CHECK( shortr.upper[1] == range.lower[1] ); // shortened range

  TEST_CHECK( shortr.lower[2] == range.lower[2] );
  TEST_CHECK( shortr.upper[2] == range.upper[2] );

  TEST_CHECK( shortr.volume == 20*20 );

  // shortr2
  TEST_CHECK( shortr2.lower[0] == range.lower[0] );
  TEST_CHECK( shortr2.upper[0] == range.upper[0] );

  TEST_CHECK( shortr2.lower[1] == range.lower[1] );
  TEST_CHECK( shortr2.upper[1] == range.lower[1]+3-1 ); // shortened range

  TEST_CHECK( shortr2.lower[2] == range.lower[2] );
  TEST_CHECK( shortr2.upper[2] == range.upper[2] );

  TEST_CHECK( shortr2.volume == 20*20*3 );
}

void test_skin_ghost()
{

  int lower[2] = {1, 1}, upper[3] = {10, 10};
  struct gkyl_range range;
  gkyl_range_init(&range, 2, lower, upper);

  // Test skin cell ranges
  struct gkyl_range lowerSkinRange1;
  gkyl_range_lower_skin(&lowerSkinRange1, &range, 0, 1);
  
  TEST_CHECK( lowerSkinRange1.volume == 10 );
  TEST_CHECK( lowerSkinRange1.lower[0] == 1 );
  TEST_CHECK( lowerSkinRange1.upper[0] == 1 );
  TEST_CHECK( lowerSkinRange1.lower[1] == 1 );
  TEST_CHECK( lowerSkinRange1.upper[1] == 10 );

  struct gkyl_range lowerSkinRange2;
  gkyl_range_lower_skin(&lowerSkinRange2, &range, 1, 1);
  
  TEST_CHECK( lowerSkinRange2.volume == 10 );
  TEST_CHECK( lowerSkinRange2.lower[0] == 1 );
  TEST_CHECK( lowerSkinRange2.upper[0] == 10 );
  TEST_CHECK( lowerSkinRange2.lower[1] == 1 );
  TEST_CHECK( lowerSkinRange2.upper[1] == 1 );

  struct gkyl_range upperSkinRange1;
  gkyl_range_upper_skin(&upperSkinRange1, &range, 0, 1);
  
  TEST_CHECK( upperSkinRange1.volume == 10 );
  TEST_CHECK( upperSkinRange1.lower[0] == 10 );
  TEST_CHECK( upperSkinRange1.upper[0] == 10 );
  TEST_CHECK( upperSkinRange1.lower[1] == 1 );
  TEST_CHECK( upperSkinRange1.upper[1] == 10 );

  struct gkyl_range upperSkinRange2;
  gkyl_range_upper_skin(&upperSkinRange2, &range, 1, 1);
  
  TEST_CHECK( upperSkinRange2.volume == 10 );
  TEST_CHECK( upperSkinRange2.lower[0] == 1 );
  TEST_CHECK( upperSkinRange2.upper[0] == 10 );
  TEST_CHECK( upperSkinRange2.lower[1] == 10 );
  TEST_CHECK( upperSkinRange2.upper[1] == 10 );

  // Test ghost cell ranges
  struct gkyl_range lowerGhostRange1;
  gkyl_range_lower_ghost(&lowerGhostRange1, &range, 0, 1);
  
  TEST_CHECK( lowerGhostRange1.volume == 10 );
  TEST_CHECK( lowerGhostRange1.lower[0] == 0 );
  TEST_CHECK( lowerGhostRange1.upper[0] == 0 );
  TEST_CHECK( lowerGhostRange1.lower[1] == 1 );
  TEST_CHECK( lowerGhostRange1.upper[1] == 10 );

  struct gkyl_range lowerGhostRange2;
  gkyl_range_lower_ghost(&lowerGhostRange2, &range, 1, 1);
  
  TEST_CHECK( lowerGhostRange2.volume == 10 );
  TEST_CHECK( lowerGhostRange2.lower[0] == 1 );
  TEST_CHECK( lowerGhostRange2.upper[0] == 10 );
  TEST_CHECK( lowerGhostRange2.lower[1] == 0 );
  TEST_CHECK( lowerGhostRange2.upper[1] == 0 );

  struct gkyl_range upperGhostRange1;
  gkyl_range_upper_ghost(&upperGhostRange1, &range, 0, 1);
  
  TEST_CHECK( upperGhostRange1.volume == 10 );
  TEST_CHECK( upperGhostRange1.lower[0] == 11 );
  TEST_CHECK( upperGhostRange1.upper[0] == 11 );
  TEST_CHECK( upperGhostRange1.lower[1] == 1 );
  TEST_CHECK( upperGhostRange1.upper[1] == 10 );

  struct gkyl_range upperGhostRange2;
  gkyl_range_upper_ghost(&upperGhostRange2, &range, 1, 1);
  
  TEST_CHECK( upperGhostRange2.volume == 10 );
  TEST_CHECK( upperGhostRange2.lower[0] == 1 );
  TEST_CHECK( upperGhostRange2.upper[0] == 10 );
  TEST_CHECK( upperGhostRange2.lower[1] == 11 );
  TEST_CHECK( upperGhostRange2.upper[1] == 11 );
}

TEST_LIST = {
  { "range_1", test_range_1 },
  { "range_shape",  test_range_shape },
  { "shorten", test_shorten },
  { "skin/ghost", test_skin_ghost },
  { NULL, NULL },
};
