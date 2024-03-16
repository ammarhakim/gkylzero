#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
GKYL_CU_DH void vlasov_prim_vars_1x3v_ser_p2(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: u_i = m1i/m0 (first vdim components), vtSq = 1/vdim(m2/m0 - udrift.udrift) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[3]; 
  const double *m1y = &moms[6]; 
  const double *m1z = &moms[9]; 
  const double *m2 = &moms[12]; 
 
  double *ux = &prim_vars[0]; 
  double *uy = &prim_vars[3]; 
  double *uz = &prim_vars[6]; 
  double *vtSq = &prim_vars[9]; 
 
  double m0_inv[3] = {0.0}; 

  double uxSq[3] = {0.0}; 
  double uySq[3] = {0.0}; 
  double uzSq[3] = {0.0}; 

  ser_1x_p2_inv(m0, m0_inv); 

  binop_mul_1d_ser_p2(m1x, m0_inv, ux); 
  binop_mul_1d_ser_p2(ux, ux, uxSq); 

  binop_mul_1d_ser_p2(m1y, m0_inv, uy); 
  binop_mul_1d_ser_p2(uy, uy, uySq); 

  binop_mul_1d_ser_p2(m1z, m0_inv, uz); 
  binop_mul_1d_ser_p2(uz, uz, uzSq); 

  vtSq[0] = 0.2357022603955158*m0_inv[2]*m2[2]+0.2357022603955158*m0_inv[1]*m2[1]-0.3333333333333333*uzSq[0]-0.3333333333333333*uySq[0]-0.3333333333333333*uxSq[0]+0.2357022603955158*m0_inv[0]*m2[0]; 
  vtSq[1] = 0.21081851067789192*m0_inv[1]*m2[2]+0.21081851067789192*m2[1]*m0_inv[2]-0.3333333333333333*uzSq[1]-0.3333333333333333*uySq[1]-0.3333333333333333*uxSq[1]+0.2357022603955158*m0_inv[0]*m2[1]+0.2357022603955158*m2[0]*m0_inv[1]; 
  vtSq[2] = -(0.3333333333333333*uzSq[2])-0.3333333333333333*uySq[2]-0.3333333333333333*uxSq[2]+0.15058465048420855*m0_inv[2]*m2[2]+0.2357022603955158*m0_inv[0]*m2[2]+0.2357022603955158*m2[0]*m0_inv[2]+0.21081851067789192*m0_inv[1]*m2[1]; 

} 
 
