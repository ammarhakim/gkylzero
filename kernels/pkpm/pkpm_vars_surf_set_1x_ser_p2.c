#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_surf_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       Input [rho ux, rho uy, rho uz], Fluid state vector.
  // p_ij:             p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_ux_l = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_ux_r = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uy_l = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_uy_r = gkyl_nmat_get(rhs, count+3); 
  struct gkyl_mat rhs_uz_l = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_uz_r = gkyl_nmat_get(rhs, count+5); 
  struct gkyl_mat rhs_Txx_l = gkyl_nmat_get(rhs, count+6); 
  struct gkyl_mat rhs_Txx_r = gkyl_nmat_get(rhs, count+7); 
  // Clear rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&rhs_ux_l, 0.0); 
  gkyl_mat_clear(&rhs_ux_r, 0.0); 
  gkyl_mat_clear(&rhs_uy_l, 0.0); 
  gkyl_mat_clear(&rhs_uy_r, 0.0); 
  gkyl_mat_clear(&rhs_uz_l, 0.0); 
  gkyl_mat_clear(&rhs_uz_r, 0.0); 
  gkyl_mat_clear(&rhs_Txx_l, 0.0); 
  gkyl_mat_clear(&rhs_Txx_r, 0.0); 
  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[3]; 
  const double *rhouz = &euler_pkpm[6]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *Pxx = &p_ij[0]; 
  double ux_l = 0.0; 
  double ux_r = 0.0; 
  double uy_l = 0.0; 
  double uy_r = 0.0; 
  double uz_l = 0.0; 
  double uz_r = 0.0; 
  double Txx_l = 0.0; 
  double Txx_r = 0.0; 
  double rhoux_l = 1.58113883008419*rhoux[2]-1.224744871391589*rhoux[1]+0.7071067811865475*rhoux[0]; 
  double rhoux_r = 1.58113883008419*rhoux[2]+1.224744871391589*rhoux[1]+0.7071067811865475*rhoux[0]; 
  double rhouy_l = 1.58113883008419*rhouy[2]-1.224744871391589*rhouy[1]+0.7071067811865475*rhouy[0]; 
  double rhouy_r = 1.58113883008419*rhouy[2]+1.224744871391589*rhouy[1]+0.7071067811865475*rhouy[0]; 
  double rhouz_l = 1.58113883008419*rhouz[2]-1.224744871391589*rhouz[1]+0.7071067811865475*rhouz[0]; 
  double rhouz_r = 1.58113883008419*rhouz[2]+1.224744871391589*rhouz[1]+0.7071067811865475*rhouz[0]; 
  double rho_l = 1.58113883008419*rho[2]-1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  double rho_r = 1.58113883008419*rho[2]+1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  double Pxx_l = 1.58113883008419*Pxx[2]-1.224744871391589*Pxx[1]+0.7071067811865475*Pxx[0]; 
  double Pxx_r = 1.58113883008419*Pxx[2]+1.224744871391589*Pxx[1]+0.7071067811865475*Pxx[0]; 
  ux_l = rhoux_l/rho_l; 
  ux_r = rhoux_r/rho_r; 
  uy_l = rhouy_l/rho_l; 
  uy_r = rhouy_r/rho_r; 
  uz_l = rhouz_l/rho_l; 
  uz_r = rhouz_r/rho_r; 
  Txx_l = 3.0*Pxx_l/rho_l; 
  Txx_r = 3.0*Pxx_r/rho_r; 
 
  gkyl_mat_set(&rhs_ux_l,0,0,ux_l); 
  gkyl_mat_set(&rhs_ux_r,0,0,ux_r); 
  gkyl_mat_set(&rhs_uy_l,0,0,uy_l); 
  gkyl_mat_set(&rhs_uy_r,0,0,uy_r); 
  gkyl_mat_set(&rhs_uz_l,0,0,uz_l); 
  gkyl_mat_set(&rhs_uz_r,0,0,uz_r); 
  gkyl_mat_set(&rhs_Txx_l,0,0,Txx_l); 
  gkyl_mat_set(&rhs_Txx_r,0,0,Txx_r); 
} 
