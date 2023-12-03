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

void
test_cubic_2d(void)
{
  double val[4] = { 1.0, 2.0, 3.0, 4.0 };
  double gradx[4] = { -1.0, -2.0, -3.0, -4.0 };
  double grady[4] = { 1.0, 2.0, 3.0, 4.0 };
  double gradxy[4] = { 1.5, 2.5, 3.5, 4.5 };
  double coeff[16] = { 0.0 };
  
  gkyl_dg_calc_cubic_2d(val, gradx, grady, gradxy, coeff);

  struct gkyl_basis b3;
  gkyl_cart_modal_tensor(&b3, 2, 3);

  TEST_CHECK( gkyl_compare_double(val[0], b3.eval_expand((double[2]) { -1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[1], b3.eval_expand((double[2]) { -1.0, 1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[2], b3.eval_expand((double[2]) { 1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(val[3], b3.eval_expand((double[2]) { 1.0, 1.0 }, coeff), 1.0e-15) );

  TEST_CHECK( gkyl_compare_double(gradx[0], b3.eval_grad_expand(0, (double[2]) { -1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(gradx[1], b3.eval_grad_expand(0, (double[2]) { -1.0, 1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(gradx[2], b3.eval_grad_expand(0, (double[2]) { 1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(gradx[3], b3.eval_grad_expand(0, (double[2]) { 1.0, 1.0 }, coeff), 1.0e-15) );

  TEST_CHECK( gkyl_compare_double(grady[0], b3.eval_grad_expand(1, (double[2]) { -1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grady[1], b3.eval_grad_expand(1, (double[2]) { -1.0, 1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grady[2], b3.eval_grad_expand(1, (double[2]) { 1.0, -1.0 }, coeff), 1.0e-15) );
  TEST_CHECK( gkyl_compare_double(grady[3], b3.eval_grad_expand(1, (double[2]) { 1.0, 1.0 }, coeff), 1.0e-15) );
}

TEST_LIST = {
  { "cubic_1d", test_cubic_1d },
  { "cubic_2d", test_cubic_2d },
  { NULL, NULL },
};

