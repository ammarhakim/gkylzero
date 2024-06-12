// Thu Jul 28 13:04:18 2022
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_comp_par_1d_ser_p2(const double *f, const double *g, double *fg, int linc1 )
{
  switch (linc1){
    case 0:
      fg[0] =  7.0710678118654757e-01*f[0]*g[0]+7.0710678118654757e-01*f[1]*g[1]+7.0710678118654757e-01*f[2]*g[2];
    case 1:
      fg[1] =  6.3245553203367588e-01*f[2]*g[1]+6.3245553203367588e-01*f[1]*g[2]+7.0710678118654757e-01*f[0]*g[1]+7.0710678118654757e-01*f[1]*g[0];
    case 2:
      fg[2] =  7.0710678118654757e-01*f[2]*g[0]+7.0710678118654757e-01*f[0]*g[2]+6.3245553203367588e-01*f[1]*g[1]+4.5175395145262565e-01*f[2]*g[2];
  }
  // nsum = 8, nprod = 22
}

struct gkyl_kern_op_count op_count_binop_mul_comp_par_1d_ser_p2(void)
{
  return (struct gkyl_kern_op_count) { .num_sum = 8, .num_prod = 22 };
}
