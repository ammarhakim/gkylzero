#include <acutest.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_range.h>

void
test_mul_1()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 1, (const int[]) { 0 }, (const int[]) { 10 });

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, range.ndim, 2);
  
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);
  gkyl_array_clear(f, 0.0);
  for (size_t i=0; i<range.volume; ++i) {
    double *d = gkyl_array_fetch(f, i);
    d[0] = 1.5;
  }  

  struct gkyl_array *g = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);
  gkyl_array_clear(g, 0.0);
  for (size_t i=0; i<range.volume; ++i) {
    double *d = gkyl_array_fetch(g, i);
    d[0] = 0.5;
  }

  struct gkyl_array *fg = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);
  gkyl_array_clear(fg, 0.0);

  // fg = f*g
  gkyl_dg_mul_op(basis, 0, fg, 0, f, 0, g);

  for (size_t i=0; i<range.volume; ++i) {
    const double *d = gkyl_array_cfetch(fg, i);
    TEST_CHECK( gkyl_compare( 0.5*1.5/sqrt(2.0), d[0], 1e-15) );

    for (int k=1; k<basis.num_basis; ++k)
      TEST_CHECK( gkyl_compare( 0.0, d[k], 1e-15) );
  }

  gkyl_array_release(f);
  gkyl_array_release(g);
  gkyl_array_release(fg);
}

void
test_mul_2()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 1, (const int[]) { 0 }, (const int[]) { 10 });

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, range.ndim, 2);

  int num_basis = basis.num_basis;

  // f is a scalar field
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, range.volume);
  gkyl_array_clear(f, 0.0);
  for (size_t i=0; i<range.volume; ++i) {
    double *d = gkyl_array_fetch(f, i);
    d[0] = 1.5;
  }  

  // g is a vector field
  struct gkyl_array *g = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, range.volume);
  gkyl_array_clear(g, 0.0);
  for (size_t i=0; i<range.volume; ++i) {
    double *d = gkyl_array_fetch(g, i);
    d[0] = 0.5;
    d[num_basis] = 0.25;
  }

  // fg is a vector field
  struct gkyl_array *fg = gkyl_array_new(GKYL_DOUBLE, 2*basis.num_basis, range.volume);
  gkyl_array_clear(fg, 0.0);

  // fg = f*g. f is a scalar and g is a vector
  for (int d=0; d<2; ++d)
    gkyl_dg_mul_op(basis, d, fg, 0, f, d, g);

  for (size_t i=0; i<range.volume; ++i) {
    const double *d = gkyl_array_cfetch(fg, i);
    TEST_CHECK( gkyl_compare( 0.5*1.5/sqrt(2.0), d[0], 1e-15) );
    TEST_CHECK( gkyl_compare( 0.25*1.5/sqrt(2.0), d[num_basis], 1e-15) );

    for (int k=1; k<basis.num_basis; ++k) {
      TEST_CHECK( gkyl_compare( 0.0, d[k], 1e-15) );
      TEST_CHECK( gkyl_compare( 0.0, d[num_basis+k], 1e-15) );
    }
  }

  gkyl_array_release(f);
  gkyl_array_release(g);
  gkyl_array_release(fg);
}

TEST_LIST = {
  { "mul_1", test_mul_1 },
  { "mul_2", test_mul_2 },
  { NULL, NULL },
};
