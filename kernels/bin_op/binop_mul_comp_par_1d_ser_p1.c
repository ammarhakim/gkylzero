#include <gkyl_binop_mul_ser.h> 
 
GKYL_CU_DH
void
binop_mul_comp_par_1d_ser_p1(const double *f, const double *g, double *fg, int linc2) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  switch (linc2) { 
    case 0: 
      fg[0] = 0.7071067811865475*f[1]*g[1]+0.7071067811865475*f[0]*g[0]; 
    case 1: 
      fg[1] = 0.7071067811865475*f[0]*g[1]+0.7071067811865475*g[0]*f[1]; 
  } 
} 
