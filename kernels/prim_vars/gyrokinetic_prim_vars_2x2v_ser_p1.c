#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void gyrokinetic_prim_vars_2x2v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:   Input moments. 
  // prim_vars: upar = m1/m0 (first component), vtSq = 1/vdim(m2/m0 - upar^2) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1 = &moms[4]; 
  const double *m2 = &moms[8]; 
 
  double *upar = &prim_vars[0]; 
  double *vtSq = &prim_vars[4]; 
 
  double m0_inv[4] = {0.0}; 

  double uparSq[4] = {0.0}; 

  // Calculate expansions of prim_vars, which can be calculated free of aliasing errors. 
  binop_mul_2d_ser_p1(m1, m0_inv, upar); 
  binop_mul_2d_ser_p1(upar, upar, uparSq); 

  vtSq[0] = 0.1666666666666667*m0_inv[3]*m2[3]+0.1666666666666667*m0_inv[2]*m2[2]+0.1666666666666667*m0_inv[1]*m2[1]-0.3333333333333333*uparSq[0]+0.1666666666666667*m0_inv[0]*m2[0]; 
  vtSq[1] = 0.1666666666666667*m0_inv[2]*m2[3]+0.1666666666666667*m2[2]*m0_inv[3]-0.3333333333333333*uparSq[1]+0.1666666666666667*m0_inv[0]*m2[1]+0.1666666666666667*m2[0]*m0_inv[1]; 
  vtSq[2] = 0.1666666666666667*m0_inv[1]*m2[3]+0.1666666666666667*m2[1]*m0_inv[3]-0.3333333333333333*uparSq[2]+0.1666666666666667*m0_inv[0]*m2[2]+0.1666666666666667*m2[0]*m0_inv[2]; 
  vtSq[3] = (-0.3333333333333333*uparSq[3])+0.1666666666666667*m0_inv[0]*m2[3]+0.1666666666666667*m2[0]*m0_inv[3]+0.1666666666666667*m0_inv[1]*m2[2]+0.1666666666666667*m2[1]*m0_inv[2]; 

} 
 
