#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void pkpm_vars_surf_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // vlasov_pkpm_moms: Input [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // euler_pkpm:       Input [rho ux, rho uy, rho uz], Fluid input state vector.
  // p_ij:             p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_ux_xl = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_ux_xr = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uy_xl = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_uy_xr = gkyl_nmat_get(rhs, count+3); 
  struct gkyl_mat rhs_uz_xl = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_uz_xr = gkyl_nmat_get(rhs, count+5); 
  struct gkyl_mat rhs_Txx_xl = gkyl_nmat_get(rhs, count+6); 
  struct gkyl_mat rhs_Txx_xr = gkyl_nmat_get(rhs, count+7); 
  struct gkyl_mat rhs_ux_yl = gkyl_nmat_get(rhs, count+8); 
  struct gkyl_mat rhs_ux_yr = gkyl_nmat_get(rhs, count+9); 
  struct gkyl_mat rhs_uy_yl = gkyl_nmat_get(rhs, count+10); 
  struct gkyl_mat rhs_uy_yr = gkyl_nmat_get(rhs, count+11); 
  struct gkyl_mat rhs_uz_yl = gkyl_nmat_get(rhs, count+12); 
  struct gkyl_mat rhs_uz_yr = gkyl_nmat_get(rhs, count+13); 
  struct gkyl_mat rhs_Tyy_yl = gkyl_nmat_get(rhs, count+14); 
  struct gkyl_mat rhs_Tyy_yr = gkyl_nmat_get(rhs, count+15); 
  // Clear rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&rhs_ux_xl, 0.0); 
  gkyl_mat_clear(&rhs_ux_xr, 0.0); 
  gkyl_mat_clear(&rhs_uy_xl, 0.0); 
  gkyl_mat_clear(&rhs_uy_xr, 0.0); 
  gkyl_mat_clear(&rhs_uz_xl, 0.0); 
  gkyl_mat_clear(&rhs_uz_xr, 0.0); 
  gkyl_mat_clear(&rhs_Txx_xl, 0.0); 
  gkyl_mat_clear(&rhs_Txx_xr, 0.0); 
  gkyl_mat_clear(&rhs_ux_yl, 0.0); 
  gkyl_mat_clear(&rhs_ux_yr, 0.0); 
  gkyl_mat_clear(&rhs_uy_yl, 0.0); 
  gkyl_mat_clear(&rhs_uy_yr, 0.0); 
  gkyl_mat_clear(&rhs_uz_yl, 0.0); 
  gkyl_mat_clear(&rhs_uz_yr, 0.0); 
  gkyl_mat_clear(&rhs_Tyy_yl, 0.0); 
  gkyl_mat_clear(&rhs_Tyy_yr, 0.0); 
  const double *rhoux = &euler_pkpm[0]; 
  const double *rhouy = &euler_pkpm[4]; 
  const double *rhouz = &euler_pkpm[8]; 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *Pxx = &p_ij[0]; 
  const double *Pyy = &p_ij[12]; 
  double ux_xl[2] = {0.0}; 
  double ux_xr[2] = {0.0}; 
  double uy_xl[2] = {0.0}; 
  double uy_xr[2] = {0.0}; 
  double uz_xl[2] = {0.0}; 
  double uz_xr[2] = {0.0}; 
  double Txx_xl[2] = {0.0}; 
  double Txx_xr[2] = {0.0}; 
  double ux_yl[2] = {0.0}; 
  double ux_yr[2] = {0.0}; 
  double uy_yl[2] = {0.0}; 
  double uy_yr[2] = {0.0}; 
  double uz_yl[2] = {0.0}; 
  double uz_yr[2] = {0.0}; 
  double Tyy_yl[2] = {0.0}; 
  double Tyy_yr[2] = {0.0}; 
  double rhoux_xl[2] = {0.0}; 
  double rhoux_xr[2] = {0.0}; 
  double rhouy_xl[2] = {0.0}; 
  double rhouy_xr[2] = {0.0}; 
  double rhouz_xl[2] = {0.0}; 
  double rhouz_xr[2] = {0.0}; 
  double rho_xl[2] = {0.0}; 
  double rho_xr[2] = {0.0}; 
  double rho_inv_xl[2] = {0.0}; 
  double rho_inv_xr[2] = {0.0}; 
  double Pxx_xl[2] = {0.0}; 
  double Pxx_xr[2] = {0.0}; 
 
  rhoux_xl[0] = 0.7071067811865475*rhoux[0]-1.224744871391589*rhoux[1]; 
  rhoux_xl[1] = 0.7071067811865475*rhoux[2]-1.224744871391589*rhoux[3]; 
  rhoux_xr[0] = 1.224744871391589*rhoux[1]+0.7071067811865475*rhoux[0]; 
  rhoux_xr[1] = 1.224744871391589*rhoux[3]+0.7071067811865475*rhoux[2]; 
  rhouy_xl[0] = 0.7071067811865475*rhouy[0]-1.224744871391589*rhouy[1]; 
  rhouy_xl[1] = 0.7071067811865475*rhouy[2]-1.224744871391589*rhouy[3]; 
  rhouy_xr[0] = 1.224744871391589*rhouy[1]+0.7071067811865475*rhouy[0]; 
  rhouy_xr[1] = 1.224744871391589*rhouy[3]+0.7071067811865475*rhouy[2]; 
  rhouz_xl[0] = 0.7071067811865475*rhouz[0]-1.224744871391589*rhouz[1]; 
  rhouz_xl[1] = 0.7071067811865475*rhouz[2]-1.224744871391589*rhouz[3]; 
  rhouz_xr[0] = 1.224744871391589*rhouz[1]+0.7071067811865475*rhouz[0]; 
  rhouz_xr[1] = 1.224744871391589*rhouz[3]+0.7071067811865475*rhouz[2]; 
  rho_xl[0] = 0.7071067811865475*rho[0]-1.224744871391589*rho[1]; 
  rho_xl[1] = 0.7071067811865475*rho[2]-1.224744871391589*rho[3]; 
  rho_xr[0] = 1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  rho_xr[1] = 1.224744871391589*rho[3]+0.7071067811865475*rho[2]; 
  Pxx_xl[0] = 2.121320343559642*Pxx[0]-3.674234614174767*Pxx[1]; 
  Pxx_xl[1] = 2.121320343559642*Pxx[2]-3.674234614174767*Pxx[3]; 
  Pxx_xr[0] = 3.674234614174767*Pxx[1]+2.121320343559642*Pxx[0]; 
  Pxx_xr[1] = 3.674234614174767*Pxx[3]+2.121320343559642*Pxx[2]; 
  ser_1x_p1_inv(rho_xl, rho_inv_xl); 
  ser_1x_p1_inv(rho_xr, rho_inv_xr); 
  binop_mul_1d_ser_p1(rho_inv_xl, rhoux_xl, ux_xl); 
  binop_mul_1d_ser_p1(rho_inv_xr, rhoux_xr, ux_xr); 
  binop_mul_1d_ser_p1(rho_inv_xl, rhouy_xl, uy_xl); 
  binop_mul_1d_ser_p1(rho_inv_xr, rhouy_xr, uy_xr); 
  binop_mul_1d_ser_p1(rho_inv_xl, rhouz_xl, uz_xl); 
  binop_mul_1d_ser_p1(rho_inv_xr, rhouz_xr, uz_xr); 
  binop_mul_1d_ser_p1(rho_inv_xl, Pxx_xl, Txx_xl); 
  binop_mul_1d_ser_p1(rho_inv_xr, Pxx_xr, Txx_xr); 
 
  double rhoux_yl[2] = {0.0}; 
  double rhoux_yr[2] = {0.0}; 
  double rhouy_yl[2] = {0.0}; 
  double rhouy_yr[2] = {0.0}; 
  double rhouz_yl[2] = {0.0}; 
  double rhouz_yr[2] = {0.0}; 
  double rho_yl[2] = {0.0}; 
  double rho_yr[2] = {0.0}; 
  double rho_inv_yl[2] = {0.0}; 
  double rho_inv_yr[2] = {0.0}; 
  double Pyy_yl[2] = {0.0}; 
  double Pyy_yr[2] = {0.0}; 
 
  rhoux_yl[0] = 0.7071067811865475*rhoux[0]-1.224744871391589*rhoux[2]; 
  rhoux_yl[1] = 0.7071067811865475*rhoux[1]-1.224744871391589*rhoux[3]; 
  rhoux_yr[0] = 1.224744871391589*rhoux[2]+0.7071067811865475*rhoux[0]; 
  rhoux_yr[1] = 1.224744871391589*rhoux[3]+0.7071067811865475*rhoux[1]; 
  rhouy_yl[0] = 0.7071067811865475*rhouy[0]-1.224744871391589*rhouy[2]; 
  rhouy_yl[1] = 0.7071067811865475*rhouy[1]-1.224744871391589*rhouy[3]; 
  rhouy_yr[0] = 1.224744871391589*rhouy[2]+0.7071067811865475*rhouy[0]; 
  rhouy_yr[1] = 1.224744871391589*rhouy[3]+0.7071067811865475*rhouy[1]; 
  rhouz_yl[0] = 0.7071067811865475*rhouz[0]-1.224744871391589*rhouz[2]; 
  rhouz_yl[1] = 0.7071067811865475*rhouz[1]-1.224744871391589*rhouz[3]; 
  rhouz_yr[0] = 1.224744871391589*rhouz[2]+0.7071067811865475*rhouz[0]; 
  rhouz_yr[1] = 1.224744871391589*rhouz[3]+0.7071067811865475*rhouz[1]; 
  rho_yl[0] = 0.7071067811865475*rho[0]-1.224744871391589*rho[2]; 
  rho_yl[1] = 0.7071067811865475*rho[1]-1.224744871391589*rho[3]; 
  rho_yr[0] = 1.224744871391589*rho[2]+0.7071067811865475*rho[0]; 
  rho_yr[1] = 1.224744871391589*rho[3]+0.7071067811865475*rho[1]; 
  Pyy_yl[0] = 2.121320343559642*Pyy[0]-3.674234614174767*Pyy[2]; 
  Pyy_yl[1] = 2.121320343559642*Pyy[1]-3.674234614174767*Pyy[3]; 
  Pyy_yr[0] = 3.674234614174767*Pyy[2]+2.121320343559642*Pyy[0]; 
  Pyy_yr[1] = 3.674234614174767*Pyy[3]+2.121320343559642*Pyy[1]; 
  ser_1x_p1_inv(rho_yl, rho_inv_yl); 
  ser_1x_p1_inv(rho_yr, rho_inv_yr); 
  binop_mul_1d_ser_p1(rho_inv_yl, rhoux_yl, ux_yl); 
  binop_mul_1d_ser_p1(rho_inv_yr, rhoux_yr, ux_yr); 
  binop_mul_1d_ser_p1(rho_inv_yl, rhouy_yl, uy_yl); 
  binop_mul_1d_ser_p1(rho_inv_yr, rhouy_yr, uy_yr); 
  binop_mul_1d_ser_p1(rho_inv_yl, rhouz_yl, uz_yl); 
  binop_mul_1d_ser_p1(rho_inv_yr, rhouz_yr, uz_yr); 
  binop_mul_1d_ser_p1(rho_inv_yl, Pyy_yl, Tyy_yl); 
  binop_mul_1d_ser_p1(rho_inv_yr, Pyy_yr, Tyy_yr); 
 
 
  gkyl_mat_set(&rhs_ux_xl,0,0,ux_xl[0]); 
  gkyl_mat_set(&rhs_ux_xr,0,0,ux_xr[0]); 
  gkyl_mat_set(&rhs_uy_xl,0,0,uy_xl[0]); 
  gkyl_mat_set(&rhs_uy_xr,0,0,uy_xr[0]); 
  gkyl_mat_set(&rhs_uz_xl,0,0,uz_xl[0]); 
  gkyl_mat_set(&rhs_uz_xr,0,0,uz_xr[0]); 
  gkyl_mat_set(&rhs_Txx_xl,0,0,Txx_xl[0]); 
  gkyl_mat_set(&rhs_Txx_xr,0,0,Txx_xr[0]); 
 
  gkyl_mat_set(&rhs_ux_yl,0,0,ux_yl[0]); 
  gkyl_mat_set(&rhs_ux_yr,0,0,ux_yr[0]); 
  gkyl_mat_set(&rhs_uy_yl,0,0,uy_yl[0]); 
  gkyl_mat_set(&rhs_uy_yr,0,0,uy_yr[0]); 
  gkyl_mat_set(&rhs_uz_yl,0,0,uz_yl[0]); 
  gkyl_mat_set(&rhs_uz_yr,0,0,uz_yr[0]); 
  gkyl_mat_set(&rhs_Tyy_yl,0,0,Tyy_yl[0]); 
  gkyl_mat_set(&rhs_Tyy_yr,0,0,Tyy_yr[0]); 
 
  gkyl_mat_set(&rhs_ux_xl,1,0,ux_xl[1]); 
  gkyl_mat_set(&rhs_ux_xr,1,0,ux_xr[1]); 
  gkyl_mat_set(&rhs_uy_xl,1,0,uy_xl[1]); 
  gkyl_mat_set(&rhs_uy_xr,1,0,uy_xr[1]); 
  gkyl_mat_set(&rhs_uz_xl,1,0,uz_xl[1]); 
  gkyl_mat_set(&rhs_uz_xr,1,0,uz_xr[1]); 
  gkyl_mat_set(&rhs_Txx_xl,1,0,Txx_xl[1]); 
  gkyl_mat_set(&rhs_Txx_xr,1,0,Txx_xr[1]); 
 
  gkyl_mat_set(&rhs_ux_yl,1,0,ux_yl[1]); 
  gkyl_mat_set(&rhs_ux_yr,1,0,ux_yr[1]); 
  gkyl_mat_set(&rhs_uy_yl,1,0,uy_yl[1]); 
  gkyl_mat_set(&rhs_uy_yr,1,0,uy_yr[1]); 
  gkyl_mat_set(&rhs_uz_yl,1,0,uz_yl[1]); 
  gkyl_mat_set(&rhs_uz_yr,1,0,uz_yr[1]); 
  gkyl_mat_set(&rhs_Tyy_yl,1,0,Tyy_yl[1]); 
  gkyl_mat_set(&rhs_Tyy_yr,1,0,Tyy_yr[1]); 
 
} 
