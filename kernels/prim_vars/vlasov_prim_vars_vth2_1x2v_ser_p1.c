#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void vlasov_prim_vars_vth2_1x2v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: vth^2. 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[2]; 
  const double *m1y = &moms[4]; 
  const double *m2 = &moms[6]; 
 
  double *vtSq = &prim_vars[0]; 
 
  double m0Sq[2] = {0.0}; 

  double m0Sq_inv[2] = {0.0}; 

  double m1xSq[2] = {0.0}; 
  double m1ySq[2] = {0.0}; 
  double m2m0[2] = {0.0}; 

  binop_mul_1d_ser_p1(m0, m0, m0Sq); 
  binop_mul_1d_ser_p1(m0, m2, m2m0); 
  ser_1x_p1_inv(m0Sq, m0Sq_inv); 
  binop_mul_1d_ser_p1(m1x, m1x, m1xSq); 

  binop_mul_1d_ser_p1(m1y, m1y, m1ySq); 

  vtSq[0] = 0.3535533905932737*m0Sq_inv[1]*m2m0[1]-0.3535533905932737*m0Sq_inv[1]*m1ySq[1]-0.3535533905932737*m0Sq_inv[1]*m1xSq[1]+0.3535533905932737*m0Sq_inv[0]*m2m0[0]-0.3535533905932737*m0Sq_inv[0]*m1ySq[0]-0.3535533905932737*m0Sq_inv[0]*m1xSq[0]; 
  vtSq[1] = 0.3535533905932737*m0Sq_inv[0]*m2m0[1]-0.3535533905932737*m0Sq_inv[0]*m1ySq[1]-0.3535533905932737*m0Sq_inv[0]*m1xSq[1]+0.3535533905932737*m2m0[0]*m0Sq_inv[1]-0.3535533905932737*m1ySq[0]*m0Sq_inv[1]-0.3535533905932737*m1xSq[0]*m0Sq_inv[1]; 

} 
 
