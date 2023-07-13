#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void gyrokinetic_prim_vars_1x1v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:   Input moments. 
  // prim_vars: upar = m1/m0 (first component), vtSq = 1/vdim(m2/m0 - upar^2) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1 = &moms[2]; 
  const double *m2 = &moms[4]; 
 
  double *upar = &prim_vars[0]; 
  double *vtSq = &prim_vars[2]; 
 
  double m0_inv[2] = {0.0}; 

  double uparSq[2] = {0.0}; 

  // Calculate expansions of prim_vars, which can be calculated free of aliasing errors. 
  binop_mul_1d_ser_p1(m1, m0_inv, upar); 
  binop_mul_1d_ser_p1(upar, upar, uparSq); 

  vtSq[0] = 0.7071067811865475*m0_inv[1]*m2[1]-1.0*uparSq[0]+0.7071067811865475*m0_inv[0]*m2[0]; 
  vtSq[1] = (-1.0*uparSq[1])+0.7071067811865475*m0_inv[0]*m2[1]+0.7071067811865475*m2[0]*m0_inv[1]; 

} 
 
