#include <gkyl_binop_mul_ser.h> 
 
GKYL_CU_DH void binop_mul_1d_tensor_p2(const double *f, const double *g, double *fg) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  fg[0] = 0.7071067811865475*f[2]*g[2]+0.7071067811865475*f[1]*g[1]+0.7071067811865475*f[0]*g[0]; 
  fg[1] = 0.6324555320336759*f[1]*g[2]+0.6324555320336759*g[1]*f[2]+0.7071067811865475*f[0]*g[1]+0.7071067811865475*g[0]*f[1]; 
  fg[2] = 0.4517539514526256*f[2]*g[2]+0.7071067811865475*f[0]*g[2]+0.7071067811865475*g[0]*f[2]+0.6324555320336759*f[1]*g[1]; 
 
} 
