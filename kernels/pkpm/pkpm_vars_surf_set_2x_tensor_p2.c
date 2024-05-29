#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_surf_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *p_ij) 
{ 
  // count:            integer to indicate which matrix being fetched. 
  // A:                preallocated LHS matrix. 
  // rhs:              preallocated RHS vector. 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // p_ij:             p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.

  struct gkyl_mat A_Txx_xl = gkyl_nmat_get(A, count); 
  struct gkyl_mat A_Txx_xr = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat A_Tyy_yl = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat A_Tyy_yr = gkyl_nmat_get(A, count+3); 
  struct gkyl_mat rhs_Txx_xl = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_Txx_xr = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_Tyy_yl = gkyl_nmat_get(rhs, count+2); 
  struct gkyl_mat rhs_Tyy_yr = gkyl_nmat_get(rhs, count+3); 
  // Clear A and rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&A_Txx_xl, 0.0); gkyl_mat_clear(&rhs_Txx_xl, 0.0); 
  gkyl_mat_clear(&A_Txx_xr, 0.0); gkyl_mat_clear(&rhs_Txx_xr, 0.0); 
  gkyl_mat_clear(&A_Tyy_yl, 0.0); gkyl_mat_clear(&rhs_Tyy_yl, 0.0); 
  gkyl_mat_clear(&A_Tyy_yr, 0.0); gkyl_mat_clear(&rhs_Tyy_yr, 0.0); 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *Pxx = &p_ij[0]; 
  const double *Pyy = &p_ij[27]; 
  double rho_xl[3] = {0.0}; 
  double rho_xr[3] = {0.0}; 
  rho_xl[0] = 1.58113883008419*rho[4]-1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  rho_xl[1] = 1.58113883008419*rho[6]-1.224744871391589*rho[3]+0.7071067811865475*rho[2]; 
  rho_xl[2] = 1.58113883008419*rho[8]-1.224744871391589*rho[7]+0.7071067811865475*rho[5]; 
  rho_xr[0] = 1.58113883008419*rho[4]+1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  rho_xr[1] = 1.58113883008419*rho[6]+1.224744871391589*rho[3]+0.7071067811865475*rho[2]; 
  rho_xr[2] = 1.58113883008419*rho[8]+1.224744871391589*rho[7]+0.7071067811865475*rho[5]; 
  double rho_yl[3] = {0.0}; 
  double rho_yr[3] = {0.0}; 
  rho_yl[0] = 1.58113883008419*rho[5]-1.224744871391589*rho[2]+0.7071067811865475*rho[0]; 
  rho_yl[1] = 1.58113883008419*rho[7]-1.224744871391589*rho[3]+0.7071067811865475*rho[1]; 
  rho_yl[2] = 1.58113883008419*rho[8]-1.224744871391589*rho[6]+0.7071067811865475*rho[4]; 
  rho_yr[0] = 1.58113883008419*rho[5]+1.224744871391589*rho[2]+0.7071067811865475*rho[0]; 
  rho_yr[1] = 1.58113883008419*rho[7]+1.224744871391589*rho[3]+0.7071067811865475*rho[1]; 
  rho_yr[2] = 1.58113883008419*rho[8]+1.224744871391589*rho[6]+0.7071067811865475*rho[4]; 
  double Pxx_xl[3] = {0.0}; 
  double Pxx_xr[3] = {0.0}; 
 
  Pxx_xl[0] = 4.743416490252569*Pxx[4]-3.674234614174767*Pxx[1]+2.121320343559642*Pxx[0]; 
  Pxx_xl[1] = 4.743416490252569*Pxx[6]-3.674234614174766*Pxx[3]+2.121320343559642*Pxx[2]; 
  Pxx_xl[2] = 4.743416490252571*Pxx[8]-3.674234614174768*Pxx[7]+2.121320343559643*Pxx[5]; 
  Pxx_xr[0] = 4.743416490252569*Pxx[4]+3.674234614174767*Pxx[1]+2.121320343559642*Pxx[0]; 
  Pxx_xr[1] = 4.743416490252569*Pxx[6]+3.674234614174766*Pxx[3]+2.121320343559642*Pxx[2]; 
  Pxx_xr[2] = 4.743416490252571*Pxx[8]+3.674234614174768*Pxx[7]+2.121320343559643*Pxx[5]; 
  double Pyy_yl[3] = {0.0}; 
  double Pyy_yr[3] = {0.0}; 
 
  Pyy_yl[0] = 4.743416490252569*Pyy[5]-3.674234614174767*Pyy[2]+2.121320343559642*Pyy[0]; 
  Pyy_yl[1] = 4.743416490252569*Pyy[7]-3.674234614174766*Pyy[3]+2.121320343559642*Pyy[1]; 
  Pyy_yl[2] = 4.743416490252571*Pyy[8]-3.674234614174768*Pyy[6]+2.121320343559643*Pyy[4]; 
  Pyy_yr[0] = 4.743416490252569*Pyy[5]+3.674234614174767*Pyy[2]+2.121320343559642*Pyy[0]; 
  Pyy_yr[1] = 4.743416490252569*Pyy[7]+3.674234614174766*Pyy[3]+2.121320343559642*Pyy[1]; 
  Pyy_yr[2] = 4.743416490252571*Pyy[8]+3.674234614174768*Pyy[6]+2.121320343559643*Pyy[4]; 
 
  gkyl_mat_set(&rhs_Txx_xl,0,0,Pxx_xl[0]); 
  gkyl_mat_set(&rhs_Txx_xr,0,0,Pxx_xr[0]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,0,0,Pyy_yl[0]); 
  gkyl_mat_set(&rhs_Tyy_yr,0,0,Pyy_yr[0]); 
 
  gkyl_mat_set(&rhs_Txx_xl,1,0,Pxx_xl[1]); 
  gkyl_mat_set(&rhs_Txx_xr,1,0,Pxx_xr[1]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,1,0,Pyy_yl[1]); 
  gkyl_mat_set(&rhs_Tyy_yr,1,0,Pyy_yr[1]); 
 
  gkyl_mat_set(&rhs_Txx_xl,2,0,Pxx_xl[2]); 
  gkyl_mat_set(&rhs_Txx_xr,2,0,Pxx_xr[2]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,2,0,Pyy_yl[2]); 
  gkyl_mat_set(&rhs_Tyy_yr,2,0,Pyy_yr[2]); 
 
  double temp_rho_xl = 0.0; 
  double temp_rho_xr = 0.0; 
  double temp_rho_yl = 0.0; 
  double temp_rho_yr = 0.0; 
  temp_rho_xl = 0.7071067811865475*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,0,0,temp_rho_xl); 
 
  temp_rho_xr = 0.7071067811865475*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,0,0,temp_rho_xr); 
 
  temp_rho_yl = 0.7071067811865475*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,0,0,temp_rho_yl); 
 
  temp_rho_yr = 0.7071067811865475*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,0,0,temp_rho_yr); 
 
  temp_rho_xl = 0.7071067811865475*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,0,1,temp_rho_xl); 
 
  temp_rho_xr = 0.7071067811865475*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,0,1,temp_rho_xr); 
 
  temp_rho_yl = 0.7071067811865475*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,0,1,temp_rho_yl); 
 
  temp_rho_yr = 0.7071067811865475*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,0,1,temp_rho_yr); 
 
  temp_rho_xl = 0.7071067811865475*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,0,2,temp_rho_xl); 
 
  temp_rho_xr = 0.7071067811865475*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,0,2,temp_rho_xr); 
 
  temp_rho_yl = 0.7071067811865475*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,0,2,temp_rho_yl); 
 
  temp_rho_yr = 0.7071067811865475*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,0,2,temp_rho_yr); 
 
  temp_rho_xl = 0.7071067811865475*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,1,0,temp_rho_xl); 
 
  temp_rho_xr = 0.7071067811865475*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,1,0,temp_rho_xr); 
 
  temp_rho_yl = 0.7071067811865475*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,1,0,temp_rho_yl); 
 
  temp_rho_yr = 0.7071067811865475*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,1,0,temp_rho_yr); 
 
  temp_rho_xl = 0.6324555320336759*rho_xl[2]+0.7071067811865475*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,1,1,temp_rho_xl); 
 
  temp_rho_xr = 0.6324555320336759*rho_xr[2]+0.7071067811865475*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,1,1,temp_rho_xr); 
 
  temp_rho_yl = 0.6324555320336759*rho_yl[2]+0.7071067811865475*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,1,1,temp_rho_yl); 
 
  temp_rho_yr = 0.6324555320336759*rho_yr[2]+0.7071067811865475*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,1,1,temp_rho_yr); 
 
  temp_rho_xl = 0.6324555320336759*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,1,2,temp_rho_xl); 
 
  temp_rho_xr = 0.6324555320336759*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,1,2,temp_rho_xr); 
 
  temp_rho_yl = 0.6324555320336759*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,1,2,temp_rho_yl); 
 
  temp_rho_yr = 0.6324555320336759*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,1,2,temp_rho_yr); 
 
  temp_rho_xl = 0.7071067811865475*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,2,0,temp_rho_xl); 
 
  temp_rho_xr = 0.7071067811865475*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,2,0,temp_rho_xr); 
 
  temp_rho_yl = 0.7071067811865475*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,2,0,temp_rho_yl); 
 
  temp_rho_yr = 0.7071067811865475*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,2,0,temp_rho_yr); 
 
  temp_rho_xl = 0.6324555320336759*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,2,1,temp_rho_xl); 
 
  temp_rho_xr = 0.6324555320336759*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,2,1,temp_rho_xr); 
 
  temp_rho_yl = 0.6324555320336759*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,2,1,temp_rho_yl); 
 
  temp_rho_yr = 0.6324555320336759*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,2,1,temp_rho_yr); 
 
  temp_rho_xl = 0.4517539514526256*rho_xl[2]+0.7071067811865475*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,2,2,temp_rho_xl); 
 
  temp_rho_xr = 0.4517539514526256*rho_xr[2]+0.7071067811865475*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,2,2,temp_rho_xr); 
 
  temp_rho_yl = 0.4517539514526256*rho_yl[2]+0.7071067811865475*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,2,2,temp_rho_yl); 
 
  temp_rho_yr = 0.4517539514526256*rho_yr[2]+0.7071067811865475*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,2,2,temp_rho_yr); 
 
} 
