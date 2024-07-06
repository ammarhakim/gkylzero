#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
GKYL_CU_DH void sr_vars_u_i_set_2x2v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i, const double *n) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // M0:    Lab frame density = GammaV*n.
  // M1i:   Lab frame flux = GammaV*n*V_drift.
  // n:     Rest-frame density.

  double n_inv[4] = {0.0}; 
  ser_2x_p1_inv(n, n_inv); 
  double temp[4] = {0.0}; 
  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs0 = gkyl_nmat_get(rhs, count); 
  gkyl_mat_clear(&rhs0, 0.0); 
  binop_mul_2d_ser_p1(n_inv, n, temp); 
  gkyl_mat_set(&rhs0,0,0,temp[0]); 
  gkyl_mat_set(&rhs0,1,0,temp[1]); 
  gkyl_mat_set(&rhs0,2,0,temp[2]); 
  gkyl_mat_set(&rhs0,3,0,temp[3]); 

  struct gkyl_mat rhs1 = gkyl_nmat_get(rhs, count+1); 
  gkyl_mat_clear(&rhs1, 0.0); 
  const double *M10 = &M1i[0]; 
  binop_mul_2d_ser_p1(n_inv, M10, temp); 
  gkyl_mat_set(&rhs1,0,0,temp[0]); 
  gkyl_mat_set(&rhs1,1,0,temp[1]); 
  gkyl_mat_set(&rhs1,2,0,temp[2]); 
  gkyl_mat_set(&rhs1,3,0,temp[3]); 

  struct gkyl_mat rhs2 = gkyl_nmat_get(rhs, count+2); 
  gkyl_mat_clear(&rhs2, 0.0); 
  const double *M11 = &M1i[4]; 
  binop_mul_2d_ser_p1(n_inv, M11, temp); 
  gkyl_mat_set(&rhs2,0,0,temp[0]); 
  gkyl_mat_set(&rhs2,1,0,temp[1]); 
  gkyl_mat_set(&rhs2,2,0,temp[2]); 
  gkyl_mat_set(&rhs2,3,0,temp[3]); 

} 
