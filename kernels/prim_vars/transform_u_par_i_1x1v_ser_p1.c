#include <gkyl_dg_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void transform_u_par_i_1x1v_ser_p1(const double *b_i, const double *moms, double* u_par_i) 
{ 
  // moms: Input moments (GK). 
  // b_i:  Contravariant components of field-aligned unit vector. 
  // u_par_i: upar_i = upar*b_i. 
 
  const double *m0 = &moms[0]; 
  const double *m1 = &moms[2]; 
  const double *b_x = &b_i[0]; 
 
  double *upar_x = &u_par_i[0]; 
 
  double m0_inv[2] = {0.0}; 

  ser_1x_p1_inv(m0, m0_inv); 
 
  double upar[2] = {0.0}; 
 
  binop_mul_1d_ser_p1(m1, m0_inv, upar); 
  binop_mul_1d_ser_p1(upar, b_x, upar_x); 
 
} 
 
