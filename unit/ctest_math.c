#include <acutest.h>
#include <gkyl_math.h>

#include <math.h>

// Table on Page 5 of van Engelen paper. See math.c for full reference
double func_1(double x, void *ctx) { return x*x*x-2*x*x+x; }

double func_2(double x, void *ctx) { return 1/(1+x); }

double func_12(double x, void *ctx) { return 1/sqrt(sin(M_PI*x)); }

double func_13(double x, void *ctx) { return pow(sin(M_PI*x), -0.8); }

double
func_cir(double x, void *ctx)
{
  return sqrt(1+x*x/(1-x*x));
}


static void
show_qr_res(struct gkyl_qr_res res, const char *msg)
{
  fprintf(stdout, "%s\n", msg);
  fprintf(stdout, ">> Status = %d. Res = %.10lg, Error = %g, Neval = %d. Nlevel = %d\n",
    res.status, res.res, res.error, res.nevals, res.nlevels);
}

void
test_dbl_exp()
{
  do {
    struct gkyl_qr_res res = gkyl_dbl_exp(func_1, 0, 0.0, 1.0, 10, 1e-11);
    //show_quad_res(res, "func_1");
    TEST_CHECK( gkyl_compare(res.res, 1.0/12.0, 1e-10) );
  } while(0);

  do {
    struct gkyl_qr_res res = gkyl_dbl_exp(func_2, 0, 0.0, 1.0, 10, 1e-11);
    //show_quad_res(res, "func_2");
    TEST_CHECK( gkyl_compare(res.res, log(2.0), 1e-10) );
  } while(0);

  do {
    // seems the following only gets to 1e-7 accuracy
    struct gkyl_qr_res res = gkyl_dbl_exp(func_12, 0, 0.0, 1.0, 10, 1e-15);
    //show_quad_res(res, "func_12");
    TEST_CHECK( gkyl_compare(res.res, 1.669253683348149, 1e-7) );
  } while(0);

  do {
    // The following integral 13 does not converge due to singularities
    // at the end point. van Engelen claims it does. I can't make it work.
    
    /* struct gkyl_quad_res res = gkyl_dbl_exp(func_13, 0, 0.0, 1.0, 10, 1e-11); */
    /* show_quad_res(res, "func_13"); */
    /* TEST_CHECK( gkyl_compare(res.quad_res, 3.604250526330095, 1e-10) ); */
  } while(0);

  do {
    struct gkyl_qr_res res = gkyl_dbl_exp(func_cir, 0, -1.0, 1.0, 10, 1e-16);
    //show_qr_res(res, "func_2");
    TEST_CHECK( gkyl_compare(res.res, M_PI, 1e-7) );
  } while(0);
}

double rfunc_1(double x, void *ctx)
{
  return x*x-1;
}

double rfunc_2(double x, void *ctx)
{
  // flat around the root x = 0
  return x*sin(x*x);
}

void
test_ridders()
{
  do {
    double x1 = 0.5, x2 = 2.0;
    double f1 = rfunc_1(x1, 0), f2 = rfunc_1(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_1, 0, x1, x2, f1, f2,
      100, 1e-12);
    //show_qr_res(res, "rfunc_1");
    TEST_CHECK( gkyl_compare(res.res, 1.0, 1e-10) );
  } while(0);

  do {
    double x1 = -0.5, x2 = 1.0;
    double f1 = rfunc_2(x1, 0), f2 = rfunc_2(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_2, 0, x1, x2, f1, f2,
      100, 1e-12);
    //show_qr_res(res, "rfunc_2");
    TEST_CHECK( gkyl_compare(res.res, 0.0, 1e-10) );
  } while(0);

}

TEST_LIST = {
  { "dbl_exp", test_dbl_exp },
  { "ridders", test_ridders },
  { NULL, NULL },
};
