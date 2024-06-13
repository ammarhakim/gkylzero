#include <gkyl_binop_mul_ser.h> 
 
GKYL_CU_DH
void
binop_mul_comp_par_3d_ser_p0(const double *f, const double *g, double *fg, int linc2) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  switch (linc2) { 
    case 0: 
      fg[0] = 0.3535533905932737*f[0]*g[0]; 
  } 
} 
