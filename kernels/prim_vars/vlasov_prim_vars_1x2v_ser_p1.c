#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void vlasov_prim_vars_1x2v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:   Input moments. 
  // prim_vars: udrift = m1/m0 (first vdim components), vtSq = 1/vdim(m2/m0 - udrift.udrift) (last component). 
 
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

  // Calculate expansions of prim_vars, which can be calculated free of aliasing errors. 
  ser_1x_p1_inv(m0, m0_inv); 

  ux[0] = 0.7071067811865475*m0_inv[1]*m1x[1]+0.7071067811865475*m0_inv[0]*m1x[0]; 
  ux[1] = 0.7071067811865475*m0_inv[0]*m1x[1]+0.7071067811865475*m1x[0]*m0_inv[1]; 
  ser_1x_p1_inv(m1x, m0_inv, ux); 
  binop_mul_1d_ser_p1(ux, ux, uxSq); 

  uy[0] = 0.7071067811865475*m0_inv[1]*m1y[1]+0.7071067811865475*m0_inv[0]*m1y[0]; 
  uy[1] = 0.7071067811865475*m0_inv[0]*m1y[1]+0.7071067811865475*m1y[0]*m0_inv[1]; 
  ser_1x_p1_inv(m1y, m0_inv, uy); 
  binop_mul_1d_ser_p1(uy, uy, uySq); 

  vtSq[0] = 0.3535533905932737*m0_inv[1]*m2[1]-0.5*uySq[0]-0.5*uxSq[0]+0.3535533905932737*m0_inv[0]*m2[0]; 
  vtSq[1] = (-0.5*uySq[1])-0.5*uxSq[1]+0.3535533905932737*m0_inv[0]*m2[1]+0.3535533905932737*m2[0]*m0_inv[1]; 

} 
 
