#include <acutest.h>

#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wv_rect_apply_bc.h>

static void
bc_copy(double t, int dir, int nc, const double *skin, double *restrict ghost, void *ctx)
{
  for (int c=0; c<nc; ++c) ghost[c] = skin[c];
}

void
test_1()
{
  int ndim = 1;
  double lower[] = {-1.0}, upper[] = {1.0};
  int cells[] = {16};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[] = { 2 };

  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  gkyl_wv_rect_apply_bc *lbc = gkyl_wv_rect_apply_bc_new(&grid, 0, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_rect_apply_bc *rbc = gkyl_wv_rect_apply_bc_new(&grid, 0, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);  

  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, 1, ext_range.volume);

  gkyl_array_clear_range(distf, 1.0, range);

  // check if ghost-cells on left/right edges of domain are 0.0
  double *data = distf->data;
  TEST_CHECK( 0.0 == data[0] );
  TEST_CHECK( 0.0 == data[1] );

  TEST_CHECK( 0.0 == data[18] );
  TEST_CHECK( 0.0 == data[19] );
  
  // apply BC
  gkyl_wv_rect_apply_bc_advance(lbc, 0.0, &range, distf);
  gkyl_wv_rect_apply_bc_advance(rbc, 0.0, &range, distf);

  // check if BCs applied correctly
  TEST_CHECK( 1.0 == data[0] );
  TEST_CHECK( 1.0 == data[1] );

  TEST_CHECK( 1.0 == data[18] );
  TEST_CHECK( 1.0 == data[19] );  

  gkyl_wv_rect_apply_bc_release(lbc);
  gkyl_wv_rect_apply_bc_release(rbc);
  gkyl_array_release(distf);
}

void
test_2()
{
  int ndim = 2;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {16, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[] = { 2, 2 };

  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  gkyl_wv_rect_apply_bc *lbc = gkyl_wv_rect_apply_bc_new(&grid, 0, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_rect_apply_bc *rbc = gkyl_wv_rect_apply_bc_new(&grid, 0, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);

  gkyl_wv_rect_apply_bc *bbc = gkyl_wv_rect_apply_bc_new(&grid, 1, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_rect_apply_bc *tbc = gkyl_wv_rect_apply_bc_new(&grid, 1, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);  

  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, 1, ext_range.volume);

  // clear interior of arrat
  gkyl_array_clear_range(distf, 1.0, range);

  // check if only interior is cleared
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ext_range);

  double vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }
  TEST_CHECK( vol == range.volume );

  // apply various BCs

  gkyl_wv_rect_apply_bc_advance(lbc, 0.0, &range, distf);
  gkyl_wv_rect_apply_bc_advance(rbc, 0.0, &range, distf);

  gkyl_wv_rect_apply_bc_advance(bbc, 0.0, &range, distf);
  gkyl_wv_rect_apply_bc_advance(tbc, 0.0, &range, distf);

  gkyl_range_iter_init(&iter, &ext_range);

  vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }

  // volume should be volume of ext_range but as corners are not
  // touched by BC updater we need to subtract the volume of the 4 corners
  TEST_CHECK( vol == ext_range.volume-4*4 );

  gkyl_wv_rect_apply_bc_release(lbc);
  gkyl_wv_rect_apply_bc_release(rbc);
  gkyl_wv_rect_apply_bc_release(bbc);
  gkyl_wv_rect_apply_bc_release(tbc);
  gkyl_array_release(distf);  
}

void
test_3()
{
  int ndim = 2;
  double lower[] = {-1.0, -1.0}, upper[] = {1.0, 1.0};
  int cells[] = {16, 8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  int nghost[] = { 2, 1 };

  struct gkyl_range ext_range, range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  gkyl_wv_rect_apply_bc *lbc = gkyl_wv_rect_apply_bc_new(&grid, 0, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_rect_apply_bc *rbc = gkyl_wv_rect_apply_bc_new(&grid, 0, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);

  gkyl_wv_rect_apply_bc *bbc = gkyl_wv_rect_apply_bc_new(&grid, 1, GKYL_LOWER_EDGE, nghost, bc_copy, NULL);
  gkyl_wv_rect_apply_bc *tbc = gkyl_wv_rect_apply_bc_new(&grid, 1, GKYL_UPPER_EDGE, nghost, bc_copy, NULL);  

  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, 1, ext_range.volume);

  // clear interior of array
  gkyl_array_clear(distf, 0.0);
  gkyl_array_clear_range(distf, 1.0, range);

  // check if only interior is cleared
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ext_range);

  double vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }
  TEST_CHECK( vol == range.volume );

  /** apply BCs on restricted range: test 1 */
  struct gkyl_range sub_range;
  gkyl_sub_range_init(&sub_range, &ext_range, (int []) { 0, 2 }, (int []) { 10, 6 });

  gkyl_wv_rect_apply_bc_advance(lbc, 0.0, &sub_range, distf);
  gkyl_wv_rect_apply_bc_advance(rbc, 0.0, &sub_range, distf);

  gkyl_wv_rect_apply_bc_advance(bbc, 0.0, &sub_range, distf);
  gkyl_wv_rect_apply_bc_advance(tbc, 0.0, &sub_range, distf);

  gkyl_range_iter_init(&iter, &ext_range);

  vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }

  // volume should be volume of range + nghost*5 (as only small
  // portion of boundary is updated)
  TEST_CHECK( vol == range.volume+2*5 );

  /** apply BCs on restricted range: test 2 */

  gkyl_array_clear(distf, 0.0);
  gkyl_array_clear_range(distf, 1.0, range);
  
  gkyl_sub_range_init(&sub_range, &ext_range, (int []) { 2, 0 }, (int []) { 8, 4 });

  gkyl_wv_rect_apply_bc_advance(lbc, 0.0, &sub_range, distf);
  gkyl_wv_rect_apply_bc_advance(rbc, 0.0, &sub_range, distf);

  gkyl_wv_rect_apply_bc_advance(bbc, 0.0, &sub_range, distf);
  gkyl_wv_rect_apply_bc_advance(tbc, 0.0, &sub_range, distf);

  gkyl_range_iter_init(&iter, &ext_range);

  vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }
  
  // volume should be volume of range + nghost*7 (as only small
  // portion of boundary is updated)
  TEST_CHECK( vol == range.volume+1*7 );

  /** apply BCs on restricted range: test 3 */

  gkyl_array_clear(distf, 0.0);
  gkyl_array_clear_range(distf, 1.0, range);
  
  gkyl_sub_range_init(&sub_range, &ext_range, (int []) { 2, 1 }, (int []) { 8, 4 });

  gkyl_wv_rect_apply_bc_advance(lbc, 0.0, &sub_range, distf);
  gkyl_wv_rect_apply_bc_advance(rbc, 0.0, &sub_range, distf);

  gkyl_wv_rect_apply_bc_advance(bbc, 0.0, &sub_range, distf);
  gkyl_wv_rect_apply_bc_advance(tbc, 0.0, &sub_range, distf);

  gkyl_range_iter_init(&iter, &ext_range);

  vol = 0.0;
  while (gkyl_range_iter_next(&iter)) {
    double *f = gkyl_array_fetch(distf, gkyl_range_idx(&ext_range, iter.idx));
    vol += f[0];
  }
  
  // range does not touch boundaries
  TEST_CHECK( vol == range.volume );  

  gkyl_wv_rect_apply_bc_release(lbc);
  gkyl_wv_rect_apply_bc_release(rbc);
  gkyl_wv_rect_apply_bc_release(bbc);
  gkyl_wv_rect_apply_bc_release(tbc);
  gkyl_array_release(distf);  
}

TEST_LIST = {
  { "test_1", test_1 },
  { "test_2", test_2 },
  { "test_3", test_3 },
  { NULL, NULL },
};
