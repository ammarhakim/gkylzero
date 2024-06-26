#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void vlasov_prim_vars_u_i_3x3v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: u_i = m1i/m0 (vdim components). 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[8]; 
  const double *m1y = &moms[16]; 
  const double *m1z = &moms[24]; 
 
  double *ux = &prim_vars[0]; 
  double *uy = &prim_vars[8]; 
  double *uz = &prim_vars[16]; 
  double m0_inv[8] = {0.0}; 

  ser_3x_p1_inv(m0, m0_inv); 

  binop_mul_3d_ser_p1(m1x, m0_inv, ux); 

  binop_mul_3d_ser_p1(m1y, m0_inv, uy); 

  binop_mul_3d_ser_p1(m1z, m0_inv, uz); 

} 
 
