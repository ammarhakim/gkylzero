#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void gyrokinetic_prim_vars_upar_1x2v_ser_p1(const double *moms, double* prim_vars) 
{ 
  // moms:      Input moments. 
  // prim_vars: upar = m1/m0. 
 
  const double *m0 = &moms[0]; 
  const double *m1 = &moms[2]; 
 
  double *upar = &prim_vars[0]; 
 
  double m0_inv[2] = {0.0}; 

  ser_1x_p1_inv(m0, m0_inv); 

  binop_mul_1d_ser_p1(m1, m0_inv, upar); 

} 
 
