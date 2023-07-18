#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void transform_vlasov_gk_prim_vars_1x2v_ser_p1(const double *b_i, const double *moms, double* prim_vars) 
{ 
  // moms: Input moments. 
  // b_i:  Contravariant components of field-aligned unit vector. 
  // prim_vars: upar = udrift . bhat, vtSq = 1/vdim(m2/m0 - upar^2) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[2]; 
  const double *m1y = &moms[4]; 
  const double *m2 = &moms[6]; 
 
  const double *b_x = &b_i[0]; 
  const double *b_y = &b_i[2]; 
  const double *b_z = &b_i[4]; 
  double *upar = &prim_vars[0]; 
  double *vtSq = &prim_vars[2]; 
 
  double m0_inv[2] = {0.0}; 

  double ux[2] = {0.0}; 
  double uy[2] = {0.0}; 
  double uparSq[2] = {0.0}; 

  ser_1x_p1_inv(m0, m0_inv); 
  binop_mul_1d_ser_p1(m1x, m0_inv, ux); 
  binop_mul_1d_ser_p1(m1y, m0_inv, uy); 

  upar[0] = 0.7071067811865475*b_y[1]*uy[1]+0.7071067811865475*b_x[1]*ux[1]+0.7071067811865475*b_y[0]*uy[0]+0.7071067811865475*b_x[0]*ux[0]; 
  upar[1] = 0.7071067811865475*b_y[0]*uy[1]+0.7071067811865475*b_x[0]*ux[1]+0.7071067811865475*uy[0]*b_y[1]+0.7071067811865475*ux[0]*b_x[1]; 

  binop_mul_1d_ser_p1(upar, upar, uparSq); 
  vtSq[0] = 0.3535533905932737*m0_inv[1]*m2[1]-0.5*uparSq[0]+0.3535533905932737*m0_inv[0]*m2[0]; 
  vtSq[1] = (-0.5*uparSq[1])+0.3535533905932737*m0_inv[0]*m2[1]+0.3535533905932737*m2[0]*m0_inv[1]; 

} 
 
