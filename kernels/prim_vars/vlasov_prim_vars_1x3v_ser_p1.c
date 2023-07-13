#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void vlasov_prim_vars_1x3v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:   Input moments. 
  // prim_vars: udrift = m1/m0 (first vdim components), vtSq = 1/vdim(m2/m0 - udrift.udrift) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[2]; 
  const double *m1y = &moms[4]; 
  const double *m1z = &moms[6]; 
  const double *m2 = &moms[8]; 
 
  double *ux = &prim_vars[0]; 
  double *uy = &prim_vars[2]; 
  double *uz = &prim_vars[4]; 
  double *vtSq = &prim_vars[6]; 
 
  double m0_inv[2] = {0.0}; 

  double uxSq[2] = {0.0}; 
  double uySq[2] = {0.0}; 
  double uzSq[2] = {0.0}; 

  // Calculate expansions of prim_vars, which can be calculated free of aliasing errors. 
  ser_1x_p1_inv(m0, m0_inv); 

  binop_mul_1d_ser_p1(m1x, m0_inv, ux); 
  binop_mul_1d_ser_p1(ux, ux, uxSq); 

  binop_mul_1d_ser_p1(m1y, m0_inv, uy); 
  binop_mul_1d_ser_p1(uy, uy, uySq); 

  binop_mul_1d_ser_p1(m1z, m0_inv, uz); 
  binop_mul_1d_ser_p1(uz, uz, uzSq); 

  vtSq[0] = 0.2357022603955158*m0_inv[1]*m2[1]-0.3333333333333333*uzSq[0]-0.3333333333333333*uySq[0]-0.3333333333333333*uxSq[0]+0.2357022603955158*m0_inv[0]*m2[0]; 
  vtSq[1] = (-0.3333333333333333*uzSq[1])-0.3333333333333333*uySq[1]-0.3333333333333333*uxSq[1]+0.2357022603955158*m0_inv[0]*m2[1]+0.2357022603955158*m2[0]*m0_inv[1]; 

} 
 
