#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void vlasov_prim_vars_u_i_1x2v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: u_i = m1i/m0 (vdim components). 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[2]; 
  const double *m1y = &moms[4]; 
 
  double *ux = &prim_vars[0]; 
  double *uy = &prim_vars[2]; 
  double m0_inv[2] = {0.0}; 

  ser_1x_p1_inv(m0, m0_inv); 

  binop_mul_1d_ser_p1(m1x, m0_inv, ux); 

  binop_mul_1d_ser_p1(m1y, m0_inv, uy); 

} 
 
