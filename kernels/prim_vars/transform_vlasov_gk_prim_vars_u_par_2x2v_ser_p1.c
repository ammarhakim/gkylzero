#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void transform_vlasov_gk_prim_vars_u_par_2x2v_ser_p1(const double *b_i, const double *moms, double* upar) 
{ 
  // moms: Input moments. 
  // b_i:  Contravariant components of field-aligned unit vector. 
  // upar: upar = udrift . bhat. 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[4]; 
  const double *m1y = &moms[8]; 
 
  const double *b_x = &b_i[0]; 
  const double *b_y = &b_i[4]; 
 
  double m0_inv[4] = {0.0}; 

  double ux[4] = {0.0}; 
  double uy[4] = {0.0}; 
  ser_2x_p1_inv(m0, m0_inv); 
  binop_mul_2d_ser_p1(m1x, m0_inv, ux); 
  binop_mul_2d_ser_p1(m1y, m0_inv, uy); 

  upar[0] = 0.5*b_y[3]*uy[3]+0.5*b_x[3]*ux[3]+0.5*b_y[2]*uy[2]+0.5*b_x[2]*ux[2]+0.5*b_y[1]*uy[1]+0.5*b_x[1]*ux[1]+0.5*b_y[0]*uy[0]+0.5*b_x[0]*ux[0]; 
  upar[1] = 0.5*b_y[2]*uy[3]+0.5*b_x[2]*ux[3]+0.5*uy[2]*b_y[3]+0.5*ux[2]*b_x[3]+0.5*b_y[0]*uy[1]+0.5*b_x[0]*ux[1]+0.5*uy[0]*b_y[1]+0.5*ux[0]*b_x[1]; 
  upar[2] = 0.5*b_y[1]*uy[3]+0.5*b_x[1]*ux[3]+0.5*uy[1]*b_y[3]+0.5*ux[1]*b_x[3]+0.5*b_y[0]*uy[2]+0.5*b_x[0]*ux[2]+0.5*uy[0]*b_y[2]+0.5*ux[0]*b_x[2]; 
  upar[3] = 0.5*b_y[0]*uy[3]+0.5*b_x[0]*ux[3]+0.5*uy[0]*b_y[3]+0.5*ux[0]*b_x[3]+0.5*b_y[1]*uy[2]+0.5*b_x[1]*ux[2]+0.5*uy[1]*b_y[2]+0.5*ux[1]*b_x[2]; 

} 
 
