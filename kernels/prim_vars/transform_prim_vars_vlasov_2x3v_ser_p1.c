#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void transform_prim_vars_vlasov_2x3v_ser_p1(const double *b_i, const double *moms, double* prim_vars) 
{ 
  // moms: Input moments (GK). 
  // b_i:  Contravariant components of field-aligned unit vector. 
  // prim_vars: u_par_i = upar*b_i  (first vdim components), vtSq = 1/vdim(m2/m0 - upar^2) (last component)  
 
  const double *m0 = &moms[0]; 
  const double *m1 = &moms[4]; 
  const double *m2 = &moms[8]; 
  const double *b_x = &b_i[0]; 
  const double *b_y = &b_i[4]; 
  const double *b_z = &b_i[8]; 
 
  double *upar_x = &prim_vars[0]; 
  double *upar_y = &prim_vars[4]; 
  double *upar_z = &prim_vars[8]; 
  double *vtSq = &prim_vars[12]; 
 
  double m0_inv[4] = {0.0}; 

  double uparSq[4] = {0.0}; 
  ser_2x_p1_inv(m0, m0_inv); 
 
  double upar[4] = {0.0}; 
 
  binop_mul_2d_ser_p1(m1, m0_inv, upar); 
  binop_mul_2d_ser_p1(upar, upar, uparSq); 
  binop_mul_2d_ser_p1(upar, b_x, upar_x); 
  binop_mul_2d_ser_p1(upar, b_y, upar_y); 
  binop_mul_2d_ser_p1(upar, b_z, upar_z); 
 
  vtSq[0] = 0.16666666666666666*m0_inv[3]*m2[3]+0.16666666666666666*m0_inv[2]*m2[2]+0.16666666666666666*m0_inv[1]*m2[1]-0.3333333333333333*uparSq[0]+0.16666666666666666*m0_inv[0]*m2[0]; 
  vtSq[1] = 0.16666666666666666*m0_inv[2]*m2[3]+0.16666666666666666*m2[2]*m0_inv[3]-0.3333333333333333*uparSq[1]+0.16666666666666666*m0_inv[0]*m2[1]+0.16666666666666666*m2[0]*m0_inv[1]; 
  vtSq[2] = 0.16666666666666666*m0_inv[1]*m2[3]+0.16666666666666666*m2[1]*m0_inv[3]-0.3333333333333333*uparSq[2]+0.16666666666666666*m0_inv[0]*m2[2]+0.16666666666666666*m2[0]*m0_inv[2]; 
  vtSq[3] = -(0.3333333333333333*uparSq[3])+0.16666666666666666*m0_inv[0]*m2[3]+0.16666666666666666*m2[0]*m0_inv[3]+0.16666666666666666*m0_inv[1]*m2[2]+0.16666666666666666*m2[1]*m0_inv[2]; 

 
} 
 
