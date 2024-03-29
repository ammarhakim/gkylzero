// Thu Jul 28 13:04:18 2022
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_1d_ser_p3(const double *f, const double *g, double *fg )
{
  double tmp[4] = {0.};
  tmp[0] =  7.0710678118654757e-01*g[3]*f[3]+7.0710678118654757e-01*g[1]*f[1]+7.0710678118654757e-01*g[2]*f[2]+7.0710678118654757e-01*g[0]*f[0];
  tmp[1] =  6.2105900340811881e-01*g[3]*f[2]+6.3245553203367588e-01*g[2]*f[1]+6.3245553203367588e-01*g[1]*f[2]+7.0710678118654757e-01*g[1]*f[0]+7.0710678118654757e-01*g[0]*f[1]+6.2105900340811881e-01*g[2]*f[3];
  tmp[2] =  7.0710678118654757e-01*g[0]*f[2]+6.2105900340811881e-01*g[1]*f[3]+4.2163702135578390e-01*g[3]*f[3]+7.0710678118654757e-01*g[2]*f[0]+6.3245553203367588e-01*g[1]*f[1]+4.5175395145262565e-01*g[2]*f[2]+6.2105900340811881e-01*g[3]*f[1];
  tmp[3] =  4.2163702135578390e-01*g[3]*f[2]+6.2105900340811881e-01*g[2]*f[1]+6.2105900340811881e-01*g[1]*f[2]+7.0710678118654757e-01*g[0]*f[3]+7.0710678118654757e-01*g[3]*f[0]+4.2163702135578390e-01*g[2]*f[3];
 
  fg[0] = tmp[0];
  fg[1] = tmp[1];
  fg[2] = tmp[2];
  fg[3] = tmp[3];
  // nsum = 19, nprod = 46
}

struct gkyl_kern_op_count op_count_binop_mul_1d_ser_p3(void)
{
  return (struct gkyl_kern_op_count) { .num_sum = 19, .num_prod = 46 };
}
