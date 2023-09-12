#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void vlasov_prim_vars_1x2v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: u_i = m1i/m0 (first vdim components), vtSq = 1/vdim(m2/m0 - udrift.udrift) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[2]; 
  const double *m1y = &moms[4]; 
  const double *m2 = &moms[6]; 
 
  double *ux = &prim_vars[0]; 
  double *uy = &prim_vars[2]; 
  double *vtSq = &prim_vars[4]; 
 
  double m0_inv[2] = {0.0}; 

  double uxSq[2] = {0.0}; 
  double uySq[2] = {0.0}; 

  ser_1x_p1_inv(m0, m0_inv); 

  binop_mul_1d_ser_p1(m1x, m0_inv, ux); 
  binop_mul_1d_ser_p1(ux, ux, uxSq); 

  binop_mul_1d_ser_p1(m1y, m0_inv, uy); 
  binop_mul_1d_ser_p1(uy, uy, uySq); 

  vtSq[0] = 0.3535533905932737*m0_inv[1]*m2[1]-0.5*uySq[0]-0.5*uxSq[0]+0.3535533905932737*m0_inv[0]*m2[0]; 
  vtSq[1] = (-0.5*uySq[1])-0.5*uxSq[1]+0.3535533905932737*m0_inv[0]*m2[1]+0.3535533905932737*m2[0]*m0_inv[1]; 

} 
 
