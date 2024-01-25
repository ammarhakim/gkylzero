#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void transform_u_par_2x2v_ser_p1(const double *b_i, const double *moms, double* upar) 
{ 
  // moms: Input moments. 
  // b_i:  Contravariant components of field-aligned unit vector. 
  // upar: upar = udrift . bhat. 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[4]; 
 
  const double *b_x = &b_i[0]; 
 
  double m0_inv[4] = {0.0}; 

  double ux[4] = {0.0}; 
  ser_2x_p1_inv(m0, m0_inv); 
  binop_mul_2d_ser_p1(m1x, m0_inv, ux); 

  upar[0] = 0.5*b_x[3]*ux[3]+0.5*b_x[2]*ux[2]+0.5*b_x[1]*ux[1]+0.5*b_x[0]*ux[0]; 
  upar[1] = 0.5*b_x[2]*ux[3]+0.5*ux[2]*b_x[3]+0.5*b_x[0]*ux[1]+0.5*ux[0]*b_x[1]; 
  upar[2] = 0.5*b_x[1]*ux[3]+0.5*ux[1]*b_x[3]+0.5*b_x[0]*ux[2]+0.5*ux[0]*b_x[2]; 
  upar[3] = 0.5*b_x[0]*ux[3]+0.5*ux[0]*b_x[3]+0.5*b_x[1]*ux[2]+0.5*ux[1]*b_x[2]; 

} 
 
