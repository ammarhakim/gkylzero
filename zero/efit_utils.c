#include <gkyl_efit_priv.h>

static inline double sq(x) { return x*x; }

static inline double dpsidx_tensor(double x, double y, const double *cell_coeffs) {
    return 5.625*cell_coeffs[8]*(2.0*x*sq(y)-0.6666666666666666*x)+2.904737509655563*cell_coeffs[7]*(sq(y)-0.3333333333333333)+5.809475019311126*cell_coeffs[6]*x*y+1.5*cell_coeffs[3]*y+3.354101966249684*cell_coeffs[4]*x+0.8660254037844386*cell_coeffs[1];
}

static inline double dpsidy_tensor(double x, double y, const double *cell_coeffs) {
    return 5.625*cell_coeffs[8]*(2.0*sq(x)*y-0.6666666666666666*y)+5.809475019311126*cell_coeffs[7]*x*y+3.354101966249684*cell_coeffs[5]*y+2.904737509655563*cell_coeffs[6]*(sq(x)-0.3333333333333333)+1.5*cell_coeffs[3]*x+0.8660254037844386*cell_coeffs[2];
}

static inline double d2psidx2_tensor(double x, double y, const double *cell_coeffs) {
    return 5.625*cell_coeffs[8]*(2.0*sq(y)-0.6666666666666666)+5.809475019311126*cell_coeffs[6]*y+3.354101966249684*cell_coeffs[4];
}

static inline double d2psidy2_tensor(double x, double y, const double *cell_coeffs) {
    return 5.625*cell_coeffs[8]*(2.0*sq(x)-0.6666666666666666)+5.809475019311126*cell_coeffs[7]*x+3.354101966249684*cell_coeffs[5];
}

static inline double dpsidxdy_tensor(double x, double y, const double *cell_coeffs) {
    return 22.5*cell_coeffs[8]*x*y+5.809475019311126*cell_coeffs[7]*y+5.809475019311126*cell_coeffs[6]*x+1.5*cell_coeffs[3];
}

static void loss_func(double* x, double* f, const double* cell_coeffs) {
    double x0 = x[0];
    double y0 = x[1];
    f[0] = dpsidx_tensor(x0, y0, cell_coeffs);
    f[1] = dpsidy_tensor(x0, y0, cell_coeffs);
}

static void jac(double* x, double J[2][2], const double* cell_coeffs) {
    double x0 = x[0];
    double y0 = x[1];
    J[0][0] = d2psidx2_tensor(x0, y0, cell_coeffs);
    J[0][1] = dpsidxdy_tensor(x0, y0, cell_coeffs);
    J[1][0] = dpsidxdy_tensor(x0, y0, cell_coeffs);
    J[1][1] = d2psidy2_tensor(x0, y0, cell_coeffs);
}

static void print_result(int n, double x[], double dx[], double errx, double errf, int niter)
{
  double x0 = x[0];
  double y0 = x[1];
  if (x0 >= -1 && x0 <= 1 && y0 >= -1 && y0 <= 1 ) {
    printf("x = ");
    for(int i=0; i<n; i++) printf(" %g ",  x[i]);
    printf("\n");
    printf("dx = ");
    for(int i=0; i<n; i++) printf("= %g ",  dx[i]);
    printf("\n");
    printf("errx = %g\n", errx);
    printf("errf = %g\n", errf);
  }

}

bool 
newton_raphson(const double *coeffs, double *xsol)
{
  int n = 2;
  double x[2] = {0.0,0.0};
  double dx[2] = {0.0,0.0};
  double fjac[2][2];
  double fjac_inv[2][2];
  double fvec[2];
  double p[2];
  int ntrial=100;
  double errx = 0.0;
  double errf = 0.0;
  for(int niter = 0; niter < ntrial; niter++) {
    loss_func(x, fvec, coeffs);
    jac(x, fjac, coeffs);
    errf = 0.0;
    for (int i=0;i<n;i++) errf += fvec[i]*fvec[i];
    errf = sqrt(errf);
    if (errf <= 1e-19) {
      for (int i=0;i<n;i++) xsol[i] = x[i];
      return true;
    }
    for (int i=0;i<n;i++) p[i] = -fvec[i];

    double det = fjac[0][0]*fjac[1][1] - fjac[0][1]*fjac[1][0];
    fjac_inv[0][0] = fjac[1][1]/det;
    fjac_inv[1][1] = fjac[0][0]/det;
    fjac_inv[0][1] = -fjac[0][1]/det;
    fjac_inv[1][0] = -fjac[1][0]/det;

    dx[0] = fjac_inv[0][0]*p[0] + fjac_inv[0][1]*p[1];
    dx[1] = fjac_inv[1][0]*p[0] + fjac_inv[1][1]*p[1];

    errx = 0.0;
    for (int i=0;i<n;i++) {
      errx += fabs(dx[i]);
      x[i] += dx[i];
    }

    if(errx<=1e-18) {
      for (int i=0;i<n;i++) xsol[i] = x[i];
      return true;
    }

  }
  return false;

}

