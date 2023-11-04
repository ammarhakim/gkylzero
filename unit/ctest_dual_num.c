#include <acutest.h>

#include <stdio.h>

#include <gkyl_dual_num.h>
#include <gkyl_math.h>

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

// for testing root finding
struct inv_map_ctx { double s0; };

static inline struct gkyl_dn
gdn_inv_mapc2p(struct gkyl_dn x, void *ctx)
{
  struct inv_map_ctx *imctx = ctx;
  double s0 = imctx->s0;
  // sqrt(x+4) - 2 - s0
  return gdn_sadd(-2.0-s0, gdn_sqrt(gdn_sadd(4.0, x)));
}

static inline double
inv_mapc2p(double x, void *ctx)
{
  return gdn_inv_mapc2p(gdn_new0(x),ctx).x[0];
}

void test_root(void)
{
  struct inv_map_ctx imctx = { 0.75 };
  double xl = 0.0, xr = 5.0;
  double fl = inv_mapc2p(xl, &imctx);
  double fr = inv_mapc2p(xr, &imctx);

  struct gkyl_qr_res root =
    gkyl_ridders(inv_mapc2p, &imctx, xl, xr, fl, fr, 100, 1e-10);

  // compute derivative at the root
  struct gkyl_dn res = gdn_inv_mapc2p(gdn_new1(root.res), &imctx);
  // at this point 1/res.x[1] is the derivative of the mapping wrt s
  // at the root
  
  TEST_CHECK( gkyl_compare_double(root.res, 3.5625, 1e-14) );
  TEST_CHECK( gkyl_compare_double(1/res.x[1], 5.5, 1e-14) );  
}

TEST_LIST = {
  { "test_basic", test_basic },
  { "test_root", test_root },
  { NULL, NULL },  
};
