#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void transform_vlasov_gk_prim_vars_2x3v_ser_p1(const double *b_i, const double *moms, double* prim_vars) 
{ 
  // moms: Input moments. 
  // b_i:  Contravariant components of field-aligned unit vector. 
  // prim_vars: upar = udrift . bhat, vtSq = 1/vdim(m2/m0 - upar^2) (last component). 
 
  const double *m0 = &moms[0]; 
  const double *m1x = &moms[4]; 
  const double *m1y = &moms[8]; 
  const double *m1z = &moms[12]; 
  const double *m2 = &moms[16]; 
 
  const double *b_x = &b_i[4]; 
  const double *b_y = &b_i[8]; 
  const double *b_z = &b_i[12]; 
  double *upar = &prim_vars[0]; 
  double *vtSq = &prim_vars[4]; 
 
  double m0_inv[4] = {0.0}; 

  double ux = {0.0}; 
  double uy = {0.0}; 
  double uz = {0.0}; 
  double uparSq[4] = {0.0}; 

  ser_2x_p1_inv(m0, m0_inv); 
  binop_mul_2d_ser_p1(m1x, m0_inv, ux); 
  binop_mul_2d_ser_p1(m1y, m0_inv, uy); 
  binop_mul_2d_ser_p1(m1z, m0_inv, uz); 

  upar[0] = 0.5*b_z[3]*uz[3]+0.5*b_y[3]*uy[3]+0.5*b_x[3]*ux[3]+0.5*b_z[2]*uz[2]+0.5*b_y[2]*uy[2]+0.5*b_x[2]*ux[2]+0.5*b_z[1]*uz[1]+0.5*b_y[1]*uy[1]+0.5*b_x[1]*ux[1]+0.5*b_z[0]*uz[0]+0.5*b_y[0]*uy[0]+0.5*b_x[0]*ux[0]; 
  upar[1] = 0.5*b_z[2]*uz[3]+0.5*b_y[2]*uy[3]+0.5*b_x[2]*ux[3]+0.5*uz[2]*b_z[3]+0.5*uy[2]*b_y[3]+0.5*ux[2]*b_x[3]+0.5*b_z[0]*uz[1]+0.5*b_y[0]*uy[1]+0.5*b_x[0]*ux[1]+0.5*uz[0]*b_z[1]+0.5*uy[0]*b_y[1]+0.5*ux[0]*b_x[1]; 
  upar[2] = 0.5*b_z[1]*uz[3]+0.5*b_y[1]*uy[3]+0.5*b_x[1]*ux[3]+0.5*uz[1]*b_z[3]+0.5*uy[1]*b_y[3]+0.5*ux[1]*b_x[3]+0.5*b_z[0]*uz[2]+0.5*b_y[0]*uy[2]+0.5*b_x[0]*ux[2]+0.5*uz[0]*b_z[2]+0.5*uy[0]*b_y[2]+0.5*ux[0]*b_x[2]; 
  upar[3] = 0.5*b_z[0]*uz[3]+0.5*b_y[0]*uy[3]+0.5*b_x[0]*ux[3]+0.5*uz[0]*b_z[3]+0.5*uy[0]*b_y[3]+0.5*ux[0]*b_x[3]+0.5*b_z[1]*uz[2]+0.5*b_y[1]*uy[2]+0.5*b_x[1]*ux[2]+0.5*uz[1]*b_z[2]+0.5*uy[1]*b_y[2]+0.5*ux[1]*b_x[2]; 

  binop_mul_2d_ser_p1(upar, upar, uparSq); 
  vtSq[0] = 0.1666666666666667*m0_inv[3]*m2[3]+0.1666666666666667*m0_inv[2]*m2[2]+0.1666666666666667*m0_inv[1]*m2[1]-0.3333333333333333*uparSq[0]+0.1666666666666667*m0_inv[0]*m2[0]; 
  vtSq[1] = 0.1666666666666667*m0_inv[2]*m2[3]+0.1666666666666667*m2[2]*m0_inv[3]-0.3333333333333333*uparSq[1]+0.1666666666666667*m0_inv[0]*m2[1]+0.1666666666666667*m2[0]*m0_inv[1]; 
  vtSq[2] = 0.1666666666666667*m0_inv[1]*m2[3]+0.1666666666666667*m2[1]*m0_inv[3]-0.3333333333333333*uparSq[2]+0.1666666666666667*m0_inv[0]*m2[2]+0.1666666666666667*m2[0]*m0_inv[2]; 
  vtSq[3] = (-0.3333333333333333*uparSq[3])+0.1666666666666667*m0_inv[0]*m2[3]+0.1666666666666667*m2[0]*m0_inv[3]+0.1666666666666667*m0_inv[1]*m2[2]+0.1666666666666667*m2[1]*m0_inv[2]; 

} 
 
