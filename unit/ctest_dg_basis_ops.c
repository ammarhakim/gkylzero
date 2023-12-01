#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_dg_basis_ops.h>

void
test_cubic_1d(void)
{
  double val[2] = { 1.0, 2.0 };
  double grad[2] = { -1.0, -2.0 };
  double coeff[4] = { 0.0 };
  
  gkyl_dg_calc_cubic_1d(val, grad, coeff);

  struct gkyl_basis b3;
  gkyl_cart_modal_serendip(&b3, 1, 3);

  TEST_CHECK( gkyl_compare_double(val[0], b3.eval_expand((double[1]) { -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[1], b3.eval_expand((double[1]) { 1.0 }, coeff), 1.0e-15) );

  TEST_CHECK( gkyl_compare_double(grad[0], b3.eval_grad_expand(0, (double[1]) { -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grad[1], b3.eval_grad_expand(0, (double[1]) { 1.0 }, coeff), 1.0e-15) );
}

TEST_LIST = {
  { "cubic_1d", test_cubic_1d },
  { NULL, NULL },
};

