// Thu Jul 28 13:04:18 2022
#include <gkyl_binop_mul_ser.h>
GKYL_CU_DH
void
binop_mul_comp_par_1d_ser_p1(const double *f, const double *g, double *fg, int linc1 )
{
  switch (linc1){
    case 0:
      fg[0] =  7.0710678118654757e-01*f[1]*g[1]+7.0710678118654757e-01*f[0]*g[0];
    case 1:
      fg[1] =  7.0710678118654757e-01*f[1]*g[0]+7.0710678118654757e-01*f[0]*g[1];
  }
  // nsum = 2, nprod = 8
}

struct gkyl_kern_op_count op_count_binop_mul_comp_par_1d_ser_p1(void)
{
  return (struct gkyl_kern_op_count) { .num_sum = 2, .num_prod = 8 };
}
