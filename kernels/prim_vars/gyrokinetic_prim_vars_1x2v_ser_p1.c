#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void gyrokinetic_prim_vars_1x2v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: upar = m1/m0 (first component), vtSq = 1/vdim(m2/m0 - upar^2) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1 = &moms[2]; 
  const double *m2 = &moms[4]; 
 
  double *upar = &prim_vars[0]; 
  double *vtSq = &prim_vars[2]; 
 
  double m0_inv[2] = {0.0}; 

  double uparSq[2] = {0.0}; 

  ser_1x_p1_inv(m0, m0_inv); 

  binop_mul_1d_ser_p1(m1, m0_inv, upar); 
  binop_mul_1d_ser_p1(upar, upar, uparSq); 

  vtSq[0] = 0.2357022603955158*m0_inv[1]*m2[1]-0.3333333333333333*uparSq[0]+0.2357022603955158*m0_inv[0]*m2[0]; 
  vtSq[1] = (-0.3333333333333333*uparSq[1])+0.2357022603955158*m0_inv[0]*m2[1]+0.2357022603955158*m2[0]*m0_inv[1]; 

} 
 
