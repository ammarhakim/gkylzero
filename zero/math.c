#include <gkyl_math.h>

#include <float.h>

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
// 1976. The following is almost a direct implementation from the
// original paper
struct gkyl_qr_res
gkyl_ridders(double (*func)(double,void*), void *ctx,
  double x1, double x2, double f1, double f2, int max_iter, double eps)
{
  double xl = x1, xr = x2, fl = f1, fr = f2;
  double res = DBL_MAX;

  int nev = 0, nitr = 0, iterating = 1;
  while(iterating && nitr<=max_iter) {
    double xm = 0.5*(xl+xr);
    double fm = func(xm, ctx); nev += 1;
    double W = sqrt(fm*fm - fl*fr);

    if (W == 0) { res = xm;  break; }
    
    double xnew = xm + (xm-xl)*dsign(fl-fr)*fm/W;
    if (fabs(xnew-res) < eps) break;
    res = xnew;
    double fnew = func(xnew, ctx); nev += 1;
    if (fnew == 0.0) break;

    if (fm*fnew < 0) {
      xl = xm; fl = fm;
      xr = res; fr = fnew;
    }
    else if (fl*fnew < 0) {
      xr = res; fr = fnew;      
    }
    else if (fr*fnew < fr) {
      xl = res; fl = fnew;
    }

    if (fabs(xr-xl) < eps) iterating = 0;
    nitr += 1;
  }
  
  return (struct gkyl_qr_res) {
    .error = fabs(xr-xl),
    .res = res,
    .nevals = nev,
    .status = nitr>max_iter ? 1 : 0,
  };
}
