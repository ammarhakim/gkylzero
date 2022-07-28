// Thu Jul 28 13:04:52 2022
#include <gkyl_binop_cross_mul_ser.h>
GKYL_CU_DH
void
binop_cross_mul_3d_6d_ser_p0(const double *f, const double *g, double *fg )
{
  double tmp[1] = {0.};
  tmp[0] = 3.5355339059327379e-01*g[0]*f[0];
 
  fg[0] = tmp[0];
}

