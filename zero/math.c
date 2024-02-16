#include <gkyl_alloc.h>
#include <gkyl_math.h>

#include <float.h>
#include <stdbool.h>
#include <complex.h>

static double ROOT_EPS = 1e-14;

// This implementation is taken from the Appendix of "Improving the
// Double Exponential Quadrature Tanh-Sinh, Sinh-Sinh and Exp-Sinh
// Formulas", Robert A. van Engelen. See
// https://www.genivia.com/qthsh.html
//
// PS: I do not think the improvements claimed in the above note are
// actually correct. Hence, at present I am using the wp34s
// implementation in Appendix B of the note.

struct gkyl_qr_res
gkyl_dbl_exp(double (*func)(double, void *), void *ctx,
  double a, double b, int n, double eps)
{
  int nev = 0;  
  double thr = 10*sqrt(eps); // too generous for larger eps, e.g. eps=1e-9
  //double thr = eps; // too generous for larger eps, e.g. eps=1e-9
  double c = (a+b)/2; // center (mean)
  double d = (b-a)/2; // half distance
  double s = func(c, ctx); nev += 1;
  double fp = 0, fm = 0;
  double p, e, v, h = 2;
  double tmax = log(2/M_PI * log((d < 1 ? 2*d : 2) / eps));
  int k = 0; // level
  do {
    double q, t;
    int j = 1;
    v = s*d*M_PI/2*h; // last sum
    p = 0;
    h /= 2;
    t = h;
    do {
      double ch = cosh(t);
      double ecs = cosh(M_PI/2 * sqrt(ch*ch - 1)); // = cosh(pi/2*sinh(t))
      double w = 1/(ecs*ecs);
      double r = sqrt(ecs*ecs - 1)/ecs;
      double x = d*r;
      if (c+x > a) {
        double y = func(c+x, ctx); nev += 1;
        if (isfinite(y))
          fp = y;
      }
      if (c-x < b) {
        double y = func(c-x, ctx); nev += 1;
        if (isfinite(y))
          fm = y;
      }
      q = ch*w*(fp+fm);
      p += q;
      j += 1+(k>0);
      t = j*h;
    } while (t <= tmax && fabs(q) > eps*fabs(p));
    s += p;
    ++k;
  } while (s && fabs(2*fabs(p) - fabs(s)) >= fabs(thr*s) && k <= n);
  s *= d*M_PI/2*h;
  e = fabs(v-s);
  if (10*e >= fabs(s)) {
    e += fabs(s);
    s = 0;
  }
  
  return (struct gkyl_qr_res) {
    .error = e,
    .res = s,
    .nevals = nev,
    .status = k>n ? 1 : 0,
    .nlevels = k
  };
}

// Helper functions for ridders
static inline double dsign(double x) { return x >= 0 ? 1 : -1; }

// See IEEE Tran. Circuit and Systems, vol CAS-26 No 11, Pg 976
// 1976. The following is almost direct implementation from the
// original paper
struct gkyl_qr_res
gkyl_ridders(double (*func)(double,void*), void *ctx,
  double xl, double xr, double fl, double fr, int max_iter, double eps)
{
  double x0 = xl, x2 = xr, f0 = fl, f2 = fr;
  double res = DBL_MAX, err = DBL_MAX;

  int nev = 0, nitr = 0, iterating = 1;
  while (iterating && nitr <= max_iter) {
    double x1 = 0.5*(x0+x2);
    double f1 = func(x1, ctx); nev += 1;
    double W = f1*f1 - f0*f2;

    double d = x2-x1;
    double x3 = x1 + dsign(f0)*f1*d/sqrt(W);
    double f3 = func(x3, ctx); nev += 1;

    if (fabs(res-x3) < eps) {
      err = fabs(res-x3);
      iterating = 0;
    }
    res = x3;

    if (f3*f0 < 0) {
      x2 = x3;
      f2 = f3;
    }
    else if (f3*f1 < 0) {
      x0 = x1 < x3 ? x1 : x3;
      f0 = x1 < x3 ? f1 : f3;

      x2 = x1 < x3 ? x3 : x1;
      f2 = x1 < x3 ? f3 : f1;
    }
    else if (f3*f2 < 0 ) {
      x0 = x3;
      f0 = f3;
    }
    nitr++;
  }
  
  return (struct gkyl_qr_res) {
    .error = err,
    .res = res,
    .nevals = nev,
    .status = nitr>max_iter ? 1 : 0,
  };
}

/////// roots of a quadratic polynomial
static struct gkyl_lo_poly_roots
quad_poly_roots(double coeff[4])
{
  double c = coeff[0], b = coeff[1], a = 1.0;

  double complex x1 = 0.0, x2 = 0.0;
  if (b>=0.0) {
    x1 = (-b-csqrt(b*b-4*a*c))/(2*a);
    x2 = 2*c/(-b-csqrt(b*b-4*a*c));
  }
  else {
    x1 = 2*c/(-b+csqrt(b*b-4*a*c));
    x2 = (-b+csqrt(b*b-4*a*c))/(2*a);
  }
  
  return (struct gkyl_lo_poly_roots) {
    .niter = 0,
    .err = { 0.0, 0.0 },
    .rpart = { creal(x1), creal(x2) },
    .impart = { cimag(x1), cimag(x2) }
  };
}

/////// roots of a cubic
static inline double complex
eval_poly3(const double coeff[4], double complex x)
{
  double complex x2 = x*x;
  double complex x3 = x2*x;
  return coeff[0] + coeff[1]*x + coeff[2]*x2 + x3;
}

static inline bool
check_converged3(double complex c1, double complex c2, double complex c3, double err[3])
{
  double eps = ROOT_EPS;
  err[0] = cabs(c1);
  err[1] = cabs(c2);
  err[2] = cabs(c3);
  
  return (err[0]+err[1]+err[2])/3.0 < eps;
}

// roots of a cubic polynomial (Durand-Kerner method)
static struct gkyl_lo_poly_roots
cubic_poly_roots(double coeff[4])
{
  double complex r1 = 0.4+0.9*I; // arbitrary complex number, not a root of unity
  double complex pn1 = r1;
  double complex qn1 = pn1*r1;
  double complex rn1 = qn1*r1;
  double complex pn = 0.0, qn = 0.0, rn = 0.0;

  double err[3] = { 0.0 };
  int max_iter = 100;
  int niter = 0;
  do {
    pn = pn1; qn = qn1; rn = rn1;
    
    pn1 = pn1 - eval_poly3(coeff, pn1)/( (pn1-qn1)*(pn1-rn1) );
    qn1 = qn1 - eval_poly3(coeff, qn1)/( (qn1-pn1)*(qn1-rn1) );
    rn1 = rn1 - eval_poly3(coeff, rn1)/( (rn1-pn1)*(rn1-qn1) );

    niter += 1;
    
  } while( !check_converged3(pn-pn1, qn-qn1, rn-rn1, err) && niter < max_iter );

  return (struct gkyl_lo_poly_roots) {
    .niter = niter,
    .err = { err[0], err[1], err[2] },
    .rpart = { creal(pn1), creal(qn1), creal(rn1) },
    .impart = { cimag(pn1), cimag(qn1), cimag(rn1) }
  };
}

/////// roots of a quartic
static inline double complex
eval_poly4(const double coeff[4], double complex x)
{
  double complex x2 = x*x;
  double complex x3 = x2*x;
  double complex x4 = x2*x2;
  return coeff[0] + coeff[1]*x + coeff[2]*x2 + coeff[3]*x3 + x4;
}

static inline bool
check_converged4(double complex c1, double complex c2, double complex c3, double complex c4,
  double err[4])
{
  double eps = ROOT_EPS;
  err[0] = cabs(c1);
  err[1] = cabs(c2);
  err[2] = cabs(c3);
  err[3] = cabs(c4);  
  
  return (err[0]+err[1]+err[2]+err[3])/4.0 < eps;
}

// roots of a quartic polynomial (Durand-Kerner method)
static struct gkyl_lo_poly_roots
quart_poly_roots(double coeff[4])
{
  double complex r1 = 0.4+0.9*I; // arbitrary complex number, not a root of unity
  double complex pn1 = r1;
  double complex qn1 = pn1*r1;
  double complex rn1 = qn1*r1;
  double complex sn1 = rn1*r1;
  double complex pn = 0.0, qn = 0.0, rn = 0.0, sn = 0.0;

  double err[4] = { 0.0 };
  int max_iter = 100;
  int niter = 0;
  do {
    pn = pn1; qn = qn1; rn = rn1; sn = sn1;
    
    pn1 = pn1 - eval_poly4(coeff, pn1)/( (pn1-qn1)*(pn1-rn1)*(pn1-sn1) );
    qn1 = qn1 - eval_poly4(coeff, qn1)/( (qn1-pn1)*(qn1-rn1)*(qn1-sn1) );
    rn1 = rn1 - eval_poly4(coeff, rn1)/( (rn1-pn1)*(rn1-qn1)*(rn1-sn1) );
    sn1 = sn1 - eval_poly4(coeff, sn1)/( (sn1-pn1)*(sn1-qn1)*(sn1-rn1) );

    niter += 1;
    
  } while( !check_converged4(pn-pn1, qn-qn1, rn-rn1, sn-sn1, err) && niter < max_iter );

  return (struct gkyl_lo_poly_roots) {
    .niter = niter,
    .err = { err[0], err[1], err[2], err[3] },
    .rpart = { creal(pn1), creal(qn1), creal(rn1), creal(sn1) },
    .impart = { cimag(pn1), cimag(qn1), cimag(rn1), cimag(sn1) }
  };
}

struct gkyl_lo_poly_roots
gkyl_calc_lo_poly_roots(enum gkyl_lo_poly_order order, double coeff[4])
{
  struct gkyl_lo_poly_roots proots = { };
  switch(order) {
    case GKYL_LO_POLY_2:
      proots = quad_poly_roots(coeff);
      break;
    case GKYL_LO_POLY_3:
      proots = cubic_poly_roots(coeff);
      break;
    case GKYL_LO_POLY_4:
      proots = quart_poly_roots(coeff);
      break;
  }
  return proots;
}

struct gkyl_poly_roots*
gkyl_poly_roots_new(int poly_order)
{
  struct gkyl_poly_roots *pr = gkyl_malloc(sizeof *pr);
  pr->impart = gkyl_malloc(sizeof(double[poly_order]));
  pr->rpart = gkyl_malloc(sizeof(double[poly_order]));
  pr->err = gkyl_malloc(sizeof(double[poly_order]));

  pr->work = gkyl_malloc(sizeof(double complex[2*poly_order]));

  pr->poly_order = poly_order;

  return pr;
}

static inline double complex
eval_poly(int poly_order, const double *coeff, double complex x)
{
  double complex xn = x;
  double complex res = coeff[0];
  for (int i=1; i<poly_order; ++i) {
    res += coeff[i]*xn;
    xn = xn*x;
  }
  return res + xn;
}

static inline bool
check_converged(int n, double complex *p1, double complex *p2, double *err)
{
  double eps = ROOT_EPS;
  double tot_err = 0.0;
  for (int i=0; i<n; ++i) {
    err[i] = cabs(p1[i]-p2[i]);
    tot_err += err[i];
  }
  return tot_err/n < eps;
}

void
gkyl_calc_poly_roots(struct gkyl_poly_roots *pr, const double *coeff)
{
  int poly_order = pr->poly_order;
  double complex *pn1 = pr->work;
  double complex *pn = pn1 + poly_order;

  double complex r1 = 0.4+0.9*I; // arbitrary complex number, not a root of unity  
  pn1[0] = r1;
  for (int i=1; i<poly_order; ++i)
    pn1[i] = pn1[i-1]*r1;

  int max_iter = 100;
  int niter = 0;
  do {
    for (int i=0; i<poly_order; ++i) {
      pn[i] = pn1[i];
      
      double complex denom = 1.0;
      for (int j=0; j<poly_order; ++j)
        if (i != j)
          denom = denom*(pn1[i]-pn1[j]);
      pn1[i] = pn1[i] - eval_poly(poly_order, coeff, pn1[i])/denom;
    }
    niter += 1;
  } while (!check_converged(poly_order, pn, pn1, pr->err) && niter < max_iter);

  pr->niter = niter;
  for (int i=0; i<poly_order; ++i) {
    pr->rpart[i] = creal(pn1[i]);
    pr->impart[i] = cimag(pn1[i]);
  }
}

void
gkyl_poly_roots_release(struct gkyl_poly_roots *pr)
{
  gkyl_free(pr->impart);
  gkyl_free(pr->rpart);
  gkyl_free(pr->err);
  gkyl_free(pr->work);
  gkyl_free(pr);
}
