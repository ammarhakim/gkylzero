#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void vlasov_prim_vars_1x1v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: u_i = m1i/m0 (first vdim components), vtSq = 1/vdim(m2/m0 - udrift.udrift) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[2]; 
  const double *m2 = &moms[4]; 
 
  double *ux = &prim_vars[0]; 
  double *vtSq = &prim_vars[2]; 
 
  double m0_inv[2] = {0.0}; 

  double uxSq[2] = {0.0}; 

  ser_1x_p1_inv(m0, m0_inv); 

  binop_mul_1d_ser_p1(m1x, m0_inv, ux); 
  binop_mul_1d_ser_p1(ux, ux, uxSq); 

  vtSq[0] = 0.7071067811865475*m0_inv[1]*m2[1]-1.0*uxSq[0]+0.7071067811865475*m0_inv[0]*m2[0]; 
  vtSq[1] = (-1.0*uxSq[1])+0.7071067811865475*m0_inv[0]*m2[1]+0.7071067811865475*m2[0]*m0_inv[1]; 

} 
 
