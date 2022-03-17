#include <acutest.h>

#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = x*x;
}

void
test_1()
{
  int poly_order = 1;
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  // left cell
  double *dfl = gkyl_array_fetch(distf, 0);
  TEST_CHECK( gkyl_compare(1.885618083164127, dfl[0], 1e-12) );
  TEST_CHECK( gkyl_compare(-1.632993161855453, dfl[1], 1e-12) );

  // right cell
  double *dfr = gkyl_array_fetch(distf, 1);
  TEST_CHECK( gkyl_compare(1.885618083164127, dfr[0], 1e-12) );
  TEST_CHECK( gkyl_compare(1.632993161855453, dfr[1], 1e-12) );

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
