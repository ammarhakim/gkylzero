#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_surf_set_1x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *p_ij) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // p_ij:             p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_Txx_l = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_Txx_r = gkyl_nmat_get(rhs, count+1); 
  // Clear rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&rhs_Txx_l, 0.0); 
  gkyl_mat_clear(&rhs_Txx_r, 0.0); 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *Pxx = &p_ij[0]; 
  double Txx_l = 0.0; 
  double Txx_r = 0.0; 
  double rho_l = 1.58113883008419*rho[2]-1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  double rho_r = 1.58113883008419*rho[2]+1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  double Pxx_l = 1.58113883008419*Pxx[2]-1.224744871391589*Pxx[1]+0.7071067811865475*Pxx[0]; 
  double Pxx_r = 1.58113883008419*Pxx[2]+1.224744871391589*Pxx[1]+0.7071067811865475*Pxx[0]; 
  Txx_l = 3.0*Pxx_l/rho_l; 
  Txx_r = 3.0*Pxx_r/rho_r; 
 
  gkyl_mat_set(&rhs_Txx_l,0,0,Txx_l); 
  gkyl_mat_set(&rhs_Txx_r,0,0,Txx_r); 
} 
