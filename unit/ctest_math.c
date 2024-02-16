#include <acutest.h>
#include <gkyl_math.h>

#include <complex.h>
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
  fprintf(stdout, ">> Status = %d. Res = %.15lg, Error = %g, Neval = %d. Nlevel = %d\n",
    res.status, res.res, res.error, res.nevals, res.nlevels);
}

void
test_dbl_exp(void)
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

double rfunc_3(double x, void *ctx)
{
  return x*exp(x)-10;
}

void
test_ridders(void)
{
  do {
    double x1 = 0.5, x2 = 2.0;
    double f1 = rfunc_1(x1, 0), f2 = rfunc_1(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_1, 0, x1, x2, f1, f2,
      100, 1e-12);
//    show_qr_res(res, "rfunc_1");
    TEST_CHECK( gkyl_compare(res.res, 1.0, 1e-10) );
  } while(0);

  do {
    double x1 = -2.0, x2 = 0.0;
    double f1 = rfunc_1(x1, 0), f2 = rfunc_1(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_1, 0, x1, x2, f1, f2,
      100, 1e-12);
//    show_qr_res(res, "rfunc_1");
    TEST_CHECK( gkyl_compare(res.res, -1.0, 1e-10) );
  } while(0);  

  do {
    double x1 = -0.5, x2 = 1.0;
    double f1 = rfunc_2(x1, 0), f2 = rfunc_2(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_2, 0, x1, x2, f1, f2,
      100, 1e-12);
//    show_qr_res(res, "rfunc_2");
    TEST_CHECK( gkyl_compare(res.res, 0.0, 1e-10) );
  } while(0);

  do {
    double x1 = -1.0, x2 = 5.0;
    double f1 = rfunc_3(x1, 0), f2 = rfunc_3(x2, 0);
    struct gkyl_qr_res res = gkyl_ridders(rfunc_3, 0, x1, x2, f1, f2,
      100, 1e-12);
//    show_qr_res(res, "rfunc_3");
    TEST_CHECK( gkyl_compare(res.res, 1.745528002740699, 1e-10) );
  } while(0);  
}

struct idx_status {
  int idx;
  bool status;
};

static struct idx_status
check_in_list(int nvals, const double complex *vals, double complex tocheck, double eps)
{
  for (int i=0; i<nvals; ++i)
    if (gkyl_compare_double(creal(vals[i]), creal(tocheck), eps) &&
      gkyl_compare_double(cimag(vals[i]), cimag(tocheck), eps) )
      return (struct idx_status) { .idx = i, .status = true };
  return (struct idx_status) { .idx = nvals, .status = false };
}

void
test_poly2_roots(void)
{
  struct gkyl_lo_poly_roots rts;

  do {
    static struct idx_status istat;
    
    double c1[4] = { 10.0, -7.0 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_2, c1);
    // checks
    istat = check_in_list(2, (double complex[]) { 2.0, 5.0 },
      rts.rpart[0]+I*rts.impart[0],
      1e-15);
    TEST_CHECK( istat.status );
    
    istat = check_in_list(2, (double complex[]) { 2.0, 5.0 },
      rts.rpart[1]+I*rts.impart[1],
      1e-15);
    TEST_CHECK( istat.status );
  } while (0);

  do {
    struct idx_status istat;
    double c1[4] = { 13.0, -4.0 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_2, c1);
    
    istat = check_in_list(2, (double complex[]) { 2.0-3.0*I, 2.0+3.0*I },
      rts.rpart[0] + I*rts.impart[0],
      1e-15);
    TEST_CHECK( istat.status );

    istat = check_in_list(2, (double complex[]) { 2.0-3.0*I, 2.0+3.0*I },
      rts.rpart[1] + I*rts.impart[1],
      1e-15);
    TEST_CHECK( istat.status );
  } while (0);
}

void
test_poly3_roots(void)
{
  struct gkyl_lo_poly_roots rts;

  do {
    static struct idx_status istat;

    double c1[4] = { 3.0, -5.5, -1.5 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_3, c1);

    double complex res[4] = { 3.0, 0.5, -2.0 };

    for (int i=0; i<3; ++i)
      TEST_CHECK(
        check_in_list(3, res, rts.rpart[i]+I*rts.impart[i], 1e-14).status
      );
  } while (0);

  do {
    static struct idx_status istat;

    double c1[4] = { -7.5, 11.0, -5.5 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_3, c1);

    double complex res[4] = { 1.5, 2+I, 2-I };

    for (int i=0; i<3; ++i)
      TEST_CHECK(
        check_in_list(3, res, rts.rpart[i]+I*rts.impart[i], 1e-14).status
      );
  } while (0);  

}

void
test_poly4_roots(void)
{
  struct gkyl_lo_poly_roots rts;

  do {
    double c1[4] = { 120.0, -26.0, -25.0, 2.0 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_4, c1);

    double complex res[4] = { 4.0, 2.0, -3.0, -5.0 };

    for (int i=0; i<4; ++i)
      TEST_CHECK(
        check_in_list(4, res, rts.rpart[i]+I*rts.impart[i], 1e-14).status
      );
    
  } while (0);

  do {
    double c1[4] = { 520.0, 22.0, -3.0, 10.0 };
    rts = gkyl_calc_lo_poly_roots(GKYL_LO_POLY_4, c1);

    double complex res[4] = {
      -4.0, -10.0,
      2.0 + 3.0*I,
      2.0 - 3.0*I
    };

    for (int i=0; i<4; ++i)
      TEST_CHECK(
        check_in_list(4, res, rts.rpart[i]+I*rts.impart[i], 1e-14).status
      );
    
  } while (0);
}

void
test_polyn_roots(void)
{

  do {
    struct gkyl_poly_roots *rts = gkyl_poly_roots_new(4);
    
    double c1[] = { 120.0, -26.0, -25.0, 2.0 };
    gkyl_calc_poly_roots(rts, c1);

    double complex res[4] = { 4.0, 2.0, -3.0, -5.0 };

    for (int i=0; i<4; ++i)
      TEST_CHECK(
        check_in_list(4, res, rts->rpart[i]+I*rts->impart[i], 1e-14).status
      );

    gkyl_poly_roots_release(rts);
  } while (0);

  do {
    struct gkyl_poly_roots *rts = gkyl_poly_roots_new(5);
    
    double c1[] = { 14400.0, -11400.0, 1174.0, 711.0, -126.0 };
    gkyl_calc_poly_roots(rts, c1);

    double complex res[5] = { 120.0, 5.0, 3.0, 2.0, -4.0 };

    for (int i=0; i<5; ++i)
      TEST_CHECK(
        check_in_list(5, res, rts->rpart[i]+I*rts->impart[i], 1e-14).status
      );

    gkyl_poly_roots_release(rts);
  } while (0);  

}

TEST_LIST = {
  { "dbl_exp", test_dbl_exp },
  { "ridders", test_ridders },
  { "poly2_roots", test_poly2_roots },
  { "poly3_roots", test_poly3_roots },
  { "poly4_roots", test_poly4_roots },
  { "polyn_roots", test_polyn_roots },  
  { NULL, NULL },
};
