#include <acutest.h>

#include <stdio.h>
#include <gkyl_dual_num.h>

static inline bool
cmp_dn(struct gkyl_dn x, struct gkyl_dn xcmp)
{
  return (x.x[0] == xcmp.x[0]) && (x.x[1] == xcmp.x[1]);
}

static void
show_dn(FILE *fp, const char *msg, struct gkyl_dn d1)
{
  fprintf(fp, "%s: [%lg, %lg]\n", msg, d1.x[0], d1.x[1]);
}

static inline
struct gkyl_dn
func_1(struct gkyl_dn x)
{
  // -2*x + x^2/(3+x^3) + 1/x
  return gdn_add(
    gdn_smul(-2.0, x),
    gdn_add(
      gdn_div(gdn_sq(x), gdn_sadd(3, gdn_cube(x))),
      gdn_inv(x)
    )
  );
}

double func_1_0(double x)
{
  return -2*x+x*x/(3+x*x*x) + 1/x;
}

double func_1_1(double x)
{
  return 2*x/(3+x*x*x) - 3*x*x*x*x/((3+x*x*x)*(3+x*x*x)) - 1/(x*x) - 2;
}

void test_basic(void)
{
  double x10 = 2.5;
  struct gkyl_dn x1 = gdn_new1(x10);

  struct gkyl_dn res = { };

  res = gdn_sq(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(x10*x10, 2*x10)) );

  res = gdn_cube(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(x10*x10*x10, 3*x10*x10)) );

  res = gdn_npow(x1, 5);
  TEST_CHECK( cmp_dn(res, gdn_new(pow(x10,5), 5*pow(x10,4))) );

  res = gdn_sqrt(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(sqrt(x10), 0.5/sqrt(x10))) );

  res = func_1(x1);
  TEST_CHECK( cmp_dn(res, gdn_new(func_1_0(x10), func_1_1(x10))) );
}

TEST_LIST = {
  { "test_basic", test_basic },
  { NULL, NULL },  
};
