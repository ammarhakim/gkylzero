#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void sr_vars_n_set_3x3v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // M0:    Lab frame density = Gamma*n.
  // M1i:   Lab frame flux = Gamma*n*V_drift_i.

  double M0_inv[8] = {0.0}; 
  ser_3x_p1_inv(M0, M0_inv); 
  double V_drift[8] = {0.0}; 
  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs0 = gkyl_nmat_get(rhs, count+0); 
  gkyl_mat_clear(&rhs0, 0.0); 
  const double *M10 = &M1i[0]; 
  binop_mul_3d_ser_p1(M0_inv, M10, V_drift); 
  gkyl_mat_set(&rhs0,0,0,V_drift[0]); 
  gkyl_mat_set(&rhs0,1,0,V_drift[1]); 
  gkyl_mat_set(&rhs0,2,0,V_drift[2]); 
  gkyl_mat_set(&rhs0,3,0,V_drift[3]); 
  gkyl_mat_set(&rhs0,4,0,V_drift[4]); 
  gkyl_mat_set(&rhs0,5,0,V_drift[5]); 
  gkyl_mat_set(&rhs0,6,0,V_drift[6]); 
  gkyl_mat_set(&rhs0,7,0,V_drift[7]); 

  struct gkyl_mat rhs1 = gkyl_nmat_get(rhs, count+1); 
  gkyl_mat_clear(&rhs1, 0.0); 
  const double *M11 = &M1i[8]; 
  binop_mul_3d_ser_p1(M0_inv, M11, V_drift); 
  gkyl_mat_set(&rhs1,0,0,V_drift[0]); 
  gkyl_mat_set(&rhs1,1,0,V_drift[1]); 
  gkyl_mat_set(&rhs1,2,0,V_drift[2]); 
  gkyl_mat_set(&rhs1,3,0,V_drift[3]); 
  gkyl_mat_set(&rhs1,4,0,V_drift[4]); 
  gkyl_mat_set(&rhs1,5,0,V_drift[5]); 
  gkyl_mat_set(&rhs1,6,0,V_drift[6]); 
  gkyl_mat_set(&rhs1,7,0,V_drift[7]); 

  struct gkyl_mat rhs2 = gkyl_nmat_get(rhs, count+2); 
  gkyl_mat_clear(&rhs2, 0.0); 
  const double *M12 = &M1i[16]; 
  binop_mul_3d_ser_p1(M0_inv, M12, V_drift); 
  gkyl_mat_set(&rhs2,0,0,V_drift[0]); 
  gkyl_mat_set(&rhs2,1,0,V_drift[1]); 
  gkyl_mat_set(&rhs2,2,0,V_drift[2]); 
  gkyl_mat_set(&rhs2,3,0,V_drift[3]); 
  gkyl_mat_set(&rhs2,4,0,V_drift[4]); 
  gkyl_mat_set(&rhs2,5,0,V_drift[5]); 
  gkyl_mat_set(&rhs2,6,0,V_drift[6]); 
  gkyl_mat_set(&rhs2,7,0,V_drift[7]); 

} 
