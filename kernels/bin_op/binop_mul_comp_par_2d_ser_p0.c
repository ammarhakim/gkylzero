// Thu Jul 28 13:04:18 2022
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_comp_par_2d_ser_p0(const double *f, const double *g, double *fg )
{
  fg[0] = 5.0000000000000000e-01*g[0]*f[0];
   // nsum = 2, nprod = 1
}

struct gkyl_kern_op_count op_count_binop_mul_comp_par_2d_ser_p0(void)
{
  return (struct gkyl_kern_op_count) { .num_sum = 2, .num_prod = 1 };
}
