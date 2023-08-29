#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void vlasov_prim_vars_2x2v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: u_i = m1i/m0 (first vdim components), vtSq = 1/vdim(m2/m0 - udrift.udrift) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[4]; 
  const double *m1y = &moms[8]; 
  const double *m2 = &moms[12]; 
 
  double *ux = &prim_vars[0]; 
  double *uy = &prim_vars[4]; 
  double *vtSq = &prim_vars[8]; 
 
  double m0_inv[4] = {0.0}; 

  double uxSq[4] = {0.0}; 
  double uySq[4] = {0.0}; 

  ser_2x_p1_inv(m0, m0_inv); 

  binop_mul_2d_ser_p1(m1x, m0_inv, ux); 
  binop_mul_2d_ser_p1(ux, ux, uxSq); 

  binop_mul_2d_ser_p1(m1y, m0_inv, uy); 
  binop_mul_2d_ser_p1(uy, uy, uySq); 

  vtSq[0] = 0.25*m0_inv[3]*m2[3]+0.25*m0_inv[2]*m2[2]+0.25*m0_inv[1]*m2[1]-0.5*uySq[0]-0.5*uxSq[0]+0.25*m0_inv[0]*m2[0]; 
  vtSq[1] = 0.25*m0_inv[2]*m2[3]+0.25*m2[2]*m0_inv[3]-0.5*uySq[1]-0.5*uxSq[1]+0.25*m0_inv[0]*m2[1]+0.25*m2[0]*m0_inv[1]; 
  vtSq[2] = 0.25*m0_inv[1]*m2[3]+0.25*m2[1]*m0_inv[3]-0.5*uySq[2]-0.5*uxSq[2]+0.25*m0_inv[0]*m2[2]+0.25*m2[0]*m0_inv[2]; 
  vtSq[3] = (-0.5*uySq[3])-0.5*uxSq[3]+0.25*m0_inv[0]*m2[3]+0.25*m2[0]*m0_inv[3]+0.25*m0_inv[1]*m2[2]+0.25*m2[1]*m0_inv[2]; 

} 
 
