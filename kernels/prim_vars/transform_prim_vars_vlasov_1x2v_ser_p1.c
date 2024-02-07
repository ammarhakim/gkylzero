#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void transform_prim_vars_vlasov_1x2v_ser_p1(const double *b_i, const double *moms, double* prim_vars) 
{ 
  // moms: Input moments (GK). 
  // b_i:  Contravariant components of field-aligned unit vector. 
  // prim_vars: u_par_i = upar*b_i  (first vdim components), vtSq = 1/vdim(m2/m0 - upar^2) (last component)  
 
  const double *m0 = &moms[0]; 
  const double *m1 = &moms[2]; 
  const double *m2 = &moms[4]; 
  const double *b_x = &b_i[0]; 
  const double *b_y = &b_i[2]; 
 
  double *upar_x = &prim_vars[0]; 
  double *upar_y = &prim_vars[2]; 
  double *vtSq = &prim_vars[4]; 
 
  double m0_inv[2] = {0.0}; 

  double uparSq[2] = {0.0}; 
  ser_1x_p1_inv(m0, m0_inv); 
 
  double upar[2] = {0.0}; 
 
  binop_mul_1d_ser_p1(m1, m0_inv, upar); 
  binop_mul_1d_ser_p1(upar, upar, uparSq); 
  binop_mul_1d_ser_p1(upar, b_x, upar_x); 
  binop_mul_1d_ser_p1(upar, b_y, upar_y); 
 
  vtSq[0] = 0.35355339059327373*m0_inv[1]*m2[1]-0.5*uparSq[0]+0.35355339059327373*m0_inv[0]*m2[0]; 
  vtSq[1] = -(0.5*uparSq[1])+0.35355339059327373*m0_inv[0]*m2[1]+0.35355339059327373*m2[0]*m0_inv[1]; 

 
} 
 
