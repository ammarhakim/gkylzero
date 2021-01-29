#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>

#define GUASS_QUAD_EPS 3.0e-15

// This is based on an implementation in Numerical Recipes in C book
static void
priv_gkyl_gauleg( gkyl_real x1, gkyl_real x2,  gkyl_real x[], gkyl_real w[], int n)
{
  int m, j, i;
  gkyl_real z1, z, xm, xl, pp, p3, p2, p1;

  m = (n+1)/2;
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for (i = 1; i <= m; i++) {
    z = cos(GKYL_PI*(i-0.25)/(n+0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp = n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1-p1/pp;
    } while( fabs(z-z1) > GUASS_QUAD_EPS );
    x[i] = xm-xl*z;
    x[n+1-i] = xm+xl*z;
    w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i] = w[i];
  }
}

void
gkyl_gauleg(gkyl_real x1, gkyl_real x2,  gkyl_real x[], gkyl_real w[], int n)
{
  priv_gkyl_gauleg(x1, x2, x-1, w-1, n); // actual routine assumes 1-offset arrays
}
