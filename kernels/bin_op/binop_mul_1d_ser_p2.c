// Thu Aug 26 15:51:37 2021
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_1d_ser_p2(const double *f, const double *g, double *fg )
{
  fg[0] =  7.0710678118654757e-01*g[0]*f[0]+7.0710678118654757e-01*f[2]*g[2]+7.0710678118654757e-01*g[1]*f[1];
  fg[1] =  7.0710678118654757e-01*g[0]*f[1]+7.0710678118654757e-01*g[1]*f[0]+6.3245553203367588e-01*g[2]*f[1]+6.3245553203367588e-01*f[2]*g[1];
  fg[2] =  4.5175395145262565e-01*f[2]*g[2]+6.3245553203367588e-01*g[1]*f[1]+7.0710678118654757e-01*g[2]*f[0]+7.0710678118654757e-01*f[2]*g[0];
  // nsum = 8, nprod = 22
}

struct gkyl_kern_op_count op_count_binop_mul_1d_ser_p2(void)
{
  return (struct gkyl_kern_op_count) { .num_sum = 8, .num_prod = 22 };
}
