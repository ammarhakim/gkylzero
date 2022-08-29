// Thu Jul 28 13:04:52 2022
#include <gkyl_binop_cross_mul_ser.h>
GKYL_CU_DH
void
binop_cross_mul_1d_2d_ser_p1(const double *f, const double *g, double *fg )
{
  double tmp[4] = {0.};
  tmp[0] =  7.0710678118654757e-01*g[1]*f[1]+7.0710678118654757e-01*g[0]*f[0];
  tmp[1] =  7.0710678118654757e-01*g[1]*f[0]+7.0710678118654757e-01*g[0]*f[1];
  tmp[2] =  7.0710678118654757e-01*g[2]*f[0]+7.0710678118654757e-01*g[3]*f[1];
  tmp[3] =  7.0710678118654757e-01*g[2]*f[1]+7.0710678118654757e-01*g[3]*f[0];
 
  fg[0] = tmp[0];
  fg[1] = tmp[1];
  fg[2] = tmp[2];
  fg[3] = tmp[3];
}

