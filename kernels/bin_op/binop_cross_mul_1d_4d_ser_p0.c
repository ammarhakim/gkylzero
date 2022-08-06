// Thu Jul 28 13:04:52 2022
#include <gkyl_binop_cross_mul_ser.h>
GKYL_CU_DH
void
binop_cross_mul_1d_4d_ser_p0(const double *f, const double *g, double *fg )
{
  double tmp[1] = {0.};
  tmp[0] = 7.0710678118654757e-01*f[0]*g[0];
 
  fg[0] = tmp[0];
}

