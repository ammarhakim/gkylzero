#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void transform_u_par_1x1v_ser_p1(const double *b_i, const double *moms, double* upar) 
{ 
  // moms: Input moments. 
  // b_i:  Contravariant components of field-aligned unit vector. 
  // upar: upar = udrift . bhat. 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[2]; 
 
  const double *b_x = &b_i[0]; 
 
  double m0_inv[2] = {0.0}; 

  double ux[2] = {0.0}; 
  ser_1x_p1_inv(m0, m0_inv); 
  binop_mul_1d_ser_p1(m1x, m0_inv, ux); 

  upar[0] = 0.7071067811865475*b_x[1]*ux[1]+0.7071067811865475*b_x[0]*ux[0]; 
  upar[1] = 0.7071067811865475*b_x[0]*ux[1]+0.7071067811865475*ux[0]*b_x[1]; 

} 
 
