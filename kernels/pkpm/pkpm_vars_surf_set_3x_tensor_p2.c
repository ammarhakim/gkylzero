#include <gkyl_mat.h> 
#include <gkyl_euler_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_vars_surf_set_3x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
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
  struct gkyl_mat A_Tzz_zl = gkyl_nmat_get(A, count+4); 
  struct gkyl_mat A_Tzz_zr = gkyl_nmat_get(A, count+5); 
  struct gkyl_mat rhs_Tzz_zl = gkyl_nmat_get(rhs, count+4); 
  struct gkyl_mat rhs_Tzz_zr = gkyl_nmat_get(rhs, count+5); 
  // Clear A and rhs for each component of primitive variables being solved for 
  gkyl_mat_clear(&A_Txx_xl, 0.0); gkyl_mat_clear(&rhs_Txx_xl, 0.0); 
  gkyl_mat_clear(&A_Txx_xr, 0.0); gkyl_mat_clear(&rhs_Txx_xr, 0.0); 
  gkyl_mat_clear(&A_Tyy_yl, 0.0); gkyl_mat_clear(&rhs_Tyy_yl, 0.0); 
  gkyl_mat_clear(&A_Tyy_yr, 0.0); gkyl_mat_clear(&rhs_Tyy_yr, 0.0); 
  gkyl_mat_clear(&A_Tzz_zl, 0.0); gkyl_mat_clear(&rhs_Tzz_zl, 0.0); 
  gkyl_mat_clear(&A_Tzz_zr, 0.0); gkyl_mat_clear(&rhs_Tzz_zr, 0.0); 
  const double *rho = &vlasov_pkpm_moms[0]; 
  const double *Pxx = &p_ij[0]; 
  const double *Pyy = &p_ij[81]; 
  const double *Pzz = &p_ij[135]; 
  double rho_xl[9] = {0.0}; 
  double rho_xr[9] = {0.0}; 
  rho_xl[0] = 1.58113883008419*rho[7]-1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  rho_xl[1] = 1.58113883008419*rho[11]-1.224744871391589*rho[4]+0.7071067811865475*rho[2]; 
  rho_xl[2] = 1.58113883008419*rho[13]-1.224744871391589*rho[5]+0.7071067811865475*rho[3]; 
  rho_xl[3] = 1.58113883008419*rho[17]-1.224744871391589*rho[10]+0.7071067811865475*rho[6]; 
  rho_xl[4] = 1.58113883008419*rho[20]-1.224744871391589*rho[12]+0.7071067811865475*rho[8]; 
  rho_xl[5] = 1.58113883008419*rho[21]-1.224744871391589*rho[15]+0.7071067811865475*rho[9]; 
  rho_xl[6] = 1.58113883008419*rho[23]-1.224744871391589*rho[18]+0.7071067811865475*rho[14]; 
  rho_xl[7] = 1.58113883008419*rho[24]-1.224744871391589*rho[19]+0.7071067811865475*rho[16]; 
  rho_xl[8] = 1.58113883008419*rho[26]-1.224744871391589*rho[25]+0.7071067811865475*rho[22]; 
  rho_xr[0] = 1.58113883008419*rho[7]+1.224744871391589*rho[1]+0.7071067811865475*rho[0]; 
  rho_xr[1] = 1.58113883008419*rho[11]+1.224744871391589*rho[4]+0.7071067811865475*rho[2]; 
  rho_xr[2] = 1.58113883008419*rho[13]+1.224744871391589*rho[5]+0.7071067811865475*rho[3]; 
  rho_xr[3] = 1.58113883008419*rho[17]+1.224744871391589*rho[10]+0.7071067811865475*rho[6]; 
  rho_xr[4] = 1.58113883008419*rho[20]+1.224744871391589*rho[12]+0.7071067811865475*rho[8]; 
  rho_xr[5] = 1.58113883008419*rho[21]+1.224744871391589*rho[15]+0.7071067811865475*rho[9]; 
  rho_xr[6] = 1.58113883008419*rho[23]+1.224744871391589*rho[18]+0.7071067811865475*rho[14]; 
  rho_xr[7] = 1.58113883008419*rho[24]+1.224744871391589*rho[19]+0.7071067811865475*rho[16]; 
  rho_xr[8] = 1.58113883008419*rho[26]+1.224744871391589*rho[25]+0.7071067811865475*rho[22]; 
  double rho_yl[9] = {0.0}; 
  double rho_yr[9] = {0.0}; 
  rho_yl[0] = 1.58113883008419*rho[8]-1.224744871391589*rho[2]+0.7071067811865475*rho[0]; 
  rho_yl[1] = 1.58113883008419*rho[12]-1.224744871391589*rho[4]+0.7071067811865475*rho[1]; 
  rho_yl[2] = 1.58113883008419*rho[14]-1.224744871391589*rho[6]+0.7071067811865475*rho[3]; 
  rho_yl[3] = 1.58113883008419*rho[18]-1.224744871391589*rho[10]+0.7071067811865475*rho[5]; 
  rho_yl[4] = 1.58113883008419*rho[20]-1.224744871391589*rho[11]+0.7071067811865475*rho[7]; 
  rho_yl[5] = 1.58113883008419*rho[22]-1.224744871391589*rho[16]+0.7071067811865475*rho[9]; 
  rho_yl[6] = 1.58113883008419*rho[23]-1.224744871391589*rho[17]+0.7071067811865475*rho[13]; 
  rho_yl[7] = 1.58113883008419*rho[25]-1.224744871391589*rho[19]+0.7071067811865475*rho[15]; 
  rho_yl[8] = 1.58113883008419*rho[26]-1.224744871391589*rho[24]+0.7071067811865475*rho[21]; 
  rho_yr[0] = 1.58113883008419*rho[8]+1.224744871391589*rho[2]+0.7071067811865475*rho[0]; 
  rho_yr[1] = 1.58113883008419*rho[12]+1.224744871391589*rho[4]+0.7071067811865475*rho[1]; 
  rho_yr[2] = 1.58113883008419*rho[14]+1.224744871391589*rho[6]+0.7071067811865475*rho[3]; 
  rho_yr[3] = 1.58113883008419*rho[18]+1.224744871391589*rho[10]+0.7071067811865475*rho[5]; 
  rho_yr[4] = 1.58113883008419*rho[20]+1.224744871391589*rho[11]+0.7071067811865475*rho[7]; 
  rho_yr[5] = 1.58113883008419*rho[22]+1.224744871391589*rho[16]+0.7071067811865475*rho[9]; 
  rho_yr[6] = 1.58113883008419*rho[23]+1.224744871391589*rho[17]+0.7071067811865475*rho[13]; 
  rho_yr[7] = 1.58113883008419*rho[25]+1.224744871391589*rho[19]+0.7071067811865475*rho[15]; 
  rho_yr[8] = 1.58113883008419*rho[26]+1.224744871391589*rho[24]+0.7071067811865475*rho[21]; 
  double rho_zl[9] = {0.0}; 
  double rho_zr[9] = {0.0}; 
  rho_zl[0] = 1.58113883008419*rho[9]-1.224744871391589*rho[3]+0.7071067811865475*rho[0]; 
  rho_zl[1] = 1.58113883008419*rho[15]-1.224744871391589*rho[5]+0.7071067811865475*rho[1]; 
  rho_zl[2] = 1.58113883008419*rho[16]-1.224744871391589*rho[6]+0.7071067811865475*rho[2]; 
  rho_zl[3] = 1.58113883008419*rho[19]-1.224744871391589*rho[10]+0.7071067811865475*rho[4]; 
  rho_zl[4] = 1.58113883008419*rho[21]-1.224744871391589*rho[13]+0.7071067811865475*rho[7]; 
  rho_zl[5] = 1.58113883008419*rho[22]-1.224744871391589*rho[14]+0.7071067811865475*rho[8]; 
  rho_zl[6] = 1.58113883008419*rho[24]-1.224744871391589*rho[17]+0.7071067811865475*rho[11]; 
  rho_zl[7] = 1.58113883008419*rho[25]-1.224744871391589*rho[18]+0.7071067811865475*rho[12]; 
  rho_zl[8] = 1.58113883008419*rho[26]-1.224744871391589*rho[23]+0.7071067811865475*rho[20]; 
  rho_zr[0] = 1.58113883008419*rho[9]+1.224744871391589*rho[3]+0.7071067811865475*rho[0]; 
  rho_zr[1] = 1.58113883008419*rho[15]+1.224744871391589*rho[5]+0.7071067811865475*rho[1]; 
  rho_zr[2] = 1.58113883008419*rho[16]+1.224744871391589*rho[6]+0.7071067811865475*rho[2]; 
  rho_zr[3] = 1.58113883008419*rho[19]+1.224744871391589*rho[10]+0.7071067811865475*rho[4]; 
  rho_zr[4] = 1.58113883008419*rho[21]+1.224744871391589*rho[13]+0.7071067811865475*rho[7]; 
  rho_zr[5] = 1.58113883008419*rho[22]+1.224744871391589*rho[14]+0.7071067811865475*rho[8]; 
  rho_zr[6] = 1.58113883008419*rho[24]+1.224744871391589*rho[17]+0.7071067811865475*rho[11]; 
  rho_zr[7] = 1.58113883008419*rho[25]+1.224744871391589*rho[18]+0.7071067811865475*rho[12]; 
  rho_zr[8] = 1.58113883008419*rho[26]+1.224744871391589*rho[23]+0.7071067811865475*rho[20]; 
  double Pxx_xl[9] = {0.0}; 
  double Pxx_xr[9] = {0.0}; 
 
  Pxx_xl[0] = 4.74341649025257*Pxx[7]-3.674234614174767*Pxx[1]+2.121320343559643*Pxx[0]; 
  Pxx_xl[1] = 4.74341649025257*Pxx[11]-3.674234614174767*Pxx[4]+2.121320343559643*Pxx[2]; 
  Pxx_xl[2] = 4.74341649025257*Pxx[13]-3.674234614174767*Pxx[5]+2.121320343559643*Pxx[3]; 
  Pxx_xl[3] = 4.74341649025257*Pxx[17]-3.674234614174767*Pxx[10]+2.121320343559643*Pxx[6]; 
  Pxx_xl[4] = 4.743416490252571*Pxx[20]-3.674234614174769*Pxx[12]+2.121320343559643*Pxx[8]; 
  Pxx_xl[5] = 4.743416490252571*Pxx[21]-3.674234614174769*Pxx[15]+2.121320343559643*Pxx[9]; 
  Pxx_xl[6] = 4.74341649025257*Pxx[23]-3.674234614174769*Pxx[18]+2.121320343559643*Pxx[14]; 
  Pxx_xl[7] = 4.74341649025257*Pxx[24]-3.674234614174769*Pxx[19]+2.121320343559643*Pxx[16]; 
  Pxx_xl[8] = 4.74341649025257*Pxx[26]-3.674234614174767*Pxx[25]+2.121320343559643*Pxx[22]; 
  Pxx_xr[0] = 4.74341649025257*Pxx[7]+3.674234614174767*Pxx[1]+2.121320343559643*Pxx[0]; 
  Pxx_xr[1] = 4.74341649025257*Pxx[11]+3.674234614174767*Pxx[4]+2.121320343559643*Pxx[2]; 
  Pxx_xr[2] = 4.74341649025257*Pxx[13]+3.674234614174767*Pxx[5]+2.121320343559643*Pxx[3]; 
  Pxx_xr[3] = 4.74341649025257*Pxx[17]+3.674234614174767*Pxx[10]+2.121320343559643*Pxx[6]; 
  Pxx_xr[4] = 4.743416490252571*Pxx[20]+3.674234614174769*Pxx[12]+2.121320343559643*Pxx[8]; 
  Pxx_xr[5] = 4.743416490252571*Pxx[21]+3.674234614174769*Pxx[15]+2.121320343559643*Pxx[9]; 
  Pxx_xr[6] = 4.74341649025257*Pxx[23]+3.674234614174769*Pxx[18]+2.121320343559643*Pxx[14]; 
  Pxx_xr[7] = 4.74341649025257*Pxx[24]+3.674234614174769*Pxx[19]+2.121320343559643*Pxx[16]; 
  Pxx_xr[8] = 4.74341649025257*Pxx[26]+3.674234614174767*Pxx[25]+2.121320343559643*Pxx[22]; 
  double Pyy_yl[9] = {0.0}; 
  double Pyy_yr[9] = {0.0}; 
 
  Pyy_yl[0] = 4.74341649025257*Pyy[8]-3.674234614174767*Pyy[2]+2.121320343559643*Pyy[0]; 
  Pyy_yl[1] = 4.74341649025257*Pyy[12]-3.674234614174767*Pyy[4]+2.121320343559643*Pyy[1]; 
  Pyy_yl[2] = 4.74341649025257*Pyy[14]-3.674234614174767*Pyy[6]+2.121320343559643*Pyy[3]; 
  Pyy_yl[3] = 4.74341649025257*Pyy[18]-3.674234614174767*Pyy[10]+2.121320343559643*Pyy[5]; 
  Pyy_yl[4] = 4.743416490252571*Pyy[20]-3.674234614174769*Pyy[11]+2.121320343559643*Pyy[7]; 
  Pyy_yl[5] = 4.743416490252571*Pyy[22]-3.674234614174769*Pyy[16]+2.121320343559643*Pyy[9]; 
  Pyy_yl[6] = 4.74341649025257*Pyy[23]-3.674234614174769*Pyy[17]+2.121320343559643*Pyy[13]; 
  Pyy_yl[7] = 4.74341649025257*Pyy[25]-3.674234614174769*Pyy[19]+2.121320343559643*Pyy[15]; 
  Pyy_yl[8] = 4.74341649025257*Pyy[26]-3.674234614174767*Pyy[24]+2.121320343559643*Pyy[21]; 
  Pyy_yr[0] = 4.74341649025257*Pyy[8]+3.674234614174767*Pyy[2]+2.121320343559643*Pyy[0]; 
  Pyy_yr[1] = 4.74341649025257*Pyy[12]+3.674234614174767*Pyy[4]+2.121320343559643*Pyy[1]; 
  Pyy_yr[2] = 4.74341649025257*Pyy[14]+3.674234614174767*Pyy[6]+2.121320343559643*Pyy[3]; 
  Pyy_yr[3] = 4.74341649025257*Pyy[18]+3.674234614174767*Pyy[10]+2.121320343559643*Pyy[5]; 
  Pyy_yr[4] = 4.743416490252571*Pyy[20]+3.674234614174769*Pyy[11]+2.121320343559643*Pyy[7]; 
  Pyy_yr[5] = 4.743416490252571*Pyy[22]+3.674234614174769*Pyy[16]+2.121320343559643*Pyy[9]; 
  Pyy_yr[6] = 4.74341649025257*Pyy[23]+3.674234614174769*Pyy[17]+2.121320343559643*Pyy[13]; 
  Pyy_yr[7] = 4.74341649025257*Pyy[25]+3.674234614174769*Pyy[19]+2.121320343559643*Pyy[15]; 
  Pyy_yr[8] = 4.74341649025257*Pyy[26]+3.674234614174767*Pyy[24]+2.121320343559643*Pyy[21]; 
 
  double Pzz_zl[9] = {0.0}; 
  double Pzz_zr[9] = {0.0}; 
 
  Pzz_zl[0] = 4.74341649025257*Pzz[9]-3.674234614174767*Pzz[3]+2.121320343559643*Pzz[0]; 
  Pzz_zl[1] = 4.74341649025257*Pzz[15]-3.674234614174767*Pzz[5]+2.121320343559643*Pzz[1]; 
  Pzz_zl[2] = 4.74341649025257*Pzz[16]-3.674234614174767*Pzz[6]+2.121320343559643*Pzz[2]; 
  Pzz_zl[3] = 4.74341649025257*Pzz[19]-3.674234614174767*Pzz[10]+2.121320343559643*Pzz[4]; 
  Pzz_zl[4] = 4.743416490252571*Pzz[21]-3.674234614174769*Pzz[13]+2.121320343559643*Pzz[7]; 
  Pzz_zl[5] = 4.743416490252571*Pzz[22]-3.674234614174769*Pzz[14]+2.121320343559643*Pzz[8]; 
  Pzz_zl[6] = 4.74341649025257*Pzz[24]-3.674234614174769*Pzz[17]+2.121320343559643*Pzz[11]; 
  Pzz_zl[7] = 4.74341649025257*Pzz[25]-3.674234614174769*Pzz[18]+2.121320343559643*Pzz[12]; 
  Pzz_zl[8] = 4.74341649025257*Pzz[26]-3.674234614174767*Pzz[23]+2.121320343559643*Pzz[20]; 
  Pzz_zr[0] = 4.74341649025257*Pzz[9]+3.674234614174767*Pzz[3]+2.121320343559643*Pzz[0]; 
  Pzz_zr[1] = 4.74341649025257*Pzz[15]+3.674234614174767*Pzz[5]+2.121320343559643*Pzz[1]; 
  Pzz_zr[2] = 4.74341649025257*Pzz[16]+3.674234614174767*Pzz[6]+2.121320343559643*Pzz[2]; 
  Pzz_zr[3] = 4.74341649025257*Pzz[19]+3.674234614174767*Pzz[10]+2.121320343559643*Pzz[4]; 
  Pzz_zr[4] = 4.743416490252571*Pzz[21]+3.674234614174769*Pzz[13]+2.121320343559643*Pzz[7]; 
  Pzz_zr[5] = 4.743416490252571*Pzz[22]+3.674234614174769*Pzz[14]+2.121320343559643*Pzz[8]; 
  Pzz_zr[6] = 4.74341649025257*Pzz[24]+3.674234614174769*Pzz[17]+2.121320343559643*Pzz[11]; 
  Pzz_zr[7] = 4.74341649025257*Pzz[25]+3.674234614174769*Pzz[18]+2.121320343559643*Pzz[12]; 
  Pzz_zr[8] = 4.74341649025257*Pzz[26]+3.674234614174767*Pzz[23]+2.121320343559643*Pzz[20]; 
 
  gkyl_mat_set(&rhs_Txx_xl,0,0,Pxx_xl[0]); 
  gkyl_mat_set(&rhs_Txx_xr,0,0,Pxx_xr[0]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,0,0,Pyy_yl[0]); 
  gkyl_mat_set(&rhs_Tyy_yr,0,0,Pyy_yr[0]); 
 
  gkyl_mat_set(&rhs_Tzz_zl,0,0,Pzz_zl[0]); 
  gkyl_mat_set(&rhs_Tzz_zr,0,0,Pzz_zr[0]); 
 
  gkyl_mat_set(&rhs_Txx_xl,1,0,Pxx_xl[1]); 
  gkyl_mat_set(&rhs_Txx_xr,1,0,Pxx_xr[1]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,1,0,Pyy_yl[1]); 
  gkyl_mat_set(&rhs_Tyy_yr,1,0,Pyy_yr[1]); 
 
  gkyl_mat_set(&rhs_Tzz_zl,1,0,Pzz_zl[1]); 
  gkyl_mat_set(&rhs_Tzz_zr,1,0,Pzz_zr[1]); 
 
  gkyl_mat_set(&rhs_Txx_xl,2,0,Pxx_xl[2]); 
  gkyl_mat_set(&rhs_Txx_xr,2,0,Pxx_xr[2]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,2,0,Pyy_yl[2]); 
  gkyl_mat_set(&rhs_Tyy_yr,2,0,Pyy_yr[2]); 
 
  gkyl_mat_set(&rhs_Tzz_zl,2,0,Pzz_zl[2]); 
  gkyl_mat_set(&rhs_Tzz_zr,2,0,Pzz_zr[2]); 
 
  gkyl_mat_set(&rhs_Txx_xl,3,0,Pxx_xl[3]); 
  gkyl_mat_set(&rhs_Txx_xr,3,0,Pxx_xr[3]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,3,0,Pyy_yl[3]); 
  gkyl_mat_set(&rhs_Tyy_yr,3,0,Pyy_yr[3]); 
 
  gkyl_mat_set(&rhs_Tzz_zl,3,0,Pzz_zl[3]); 
  gkyl_mat_set(&rhs_Tzz_zr,3,0,Pzz_zr[3]); 
 
  gkyl_mat_set(&rhs_Txx_xl,4,0,Pxx_xl[4]); 
  gkyl_mat_set(&rhs_Txx_xr,4,0,Pxx_xr[4]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,4,0,Pyy_yl[4]); 
  gkyl_mat_set(&rhs_Tyy_yr,4,0,Pyy_yr[4]); 
 
  gkyl_mat_set(&rhs_Tzz_zl,4,0,Pzz_zl[4]); 
  gkyl_mat_set(&rhs_Tzz_zr,4,0,Pzz_zr[4]); 
 
  gkyl_mat_set(&rhs_Txx_xl,5,0,Pxx_xl[5]); 
  gkyl_mat_set(&rhs_Txx_xr,5,0,Pxx_xr[5]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,5,0,Pyy_yl[5]); 
  gkyl_mat_set(&rhs_Tyy_yr,5,0,Pyy_yr[5]); 
 
  gkyl_mat_set(&rhs_Tzz_zl,5,0,Pzz_zl[5]); 
  gkyl_mat_set(&rhs_Tzz_zr,5,0,Pzz_zr[5]); 
 
  gkyl_mat_set(&rhs_Txx_xl,6,0,Pxx_xl[6]); 
  gkyl_mat_set(&rhs_Txx_xr,6,0,Pxx_xr[6]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,6,0,Pyy_yl[6]); 
  gkyl_mat_set(&rhs_Tyy_yr,6,0,Pyy_yr[6]); 
 
  gkyl_mat_set(&rhs_Tzz_zl,6,0,Pzz_zl[6]); 
  gkyl_mat_set(&rhs_Tzz_zr,6,0,Pzz_zr[6]); 
 
  gkyl_mat_set(&rhs_Txx_xl,7,0,Pxx_xl[7]); 
  gkyl_mat_set(&rhs_Txx_xr,7,0,Pxx_xr[7]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,7,0,Pyy_yl[7]); 
  gkyl_mat_set(&rhs_Tyy_yr,7,0,Pyy_yr[7]); 
 
  gkyl_mat_set(&rhs_Tzz_zl,7,0,Pzz_zl[7]); 
  gkyl_mat_set(&rhs_Tzz_zr,7,0,Pzz_zr[7]); 
 
  gkyl_mat_set(&rhs_Txx_xl,8,0,Pxx_xl[8]); 
  gkyl_mat_set(&rhs_Txx_xr,8,0,Pxx_xr[8]); 
 
  gkyl_mat_set(&rhs_Tyy_yl,8,0,Pyy_yl[8]); 
  gkyl_mat_set(&rhs_Tyy_yr,8,0,Pyy_yr[8]); 
 
  gkyl_mat_set(&rhs_Tzz_zl,8,0,Pzz_zl[8]); 
  gkyl_mat_set(&rhs_Tzz_zr,8,0,Pzz_zr[8]); 
 
  double temp_rho_xl = 0.0; 
  double temp_rho_xr = 0.0; 
  double temp_rho_yl = 0.0; 
  double temp_rho_yr = 0.0; 
  double temp_rho_zl = 0.0; 
  double temp_rho_zr = 0.0; 
  temp_rho_xl = 0.5*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,0,0,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,0,0,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,0,0,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,0,0,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[0]; 
  gkyl_mat_set(&A_Tzz_zl,0,0,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[0]; 
  gkyl_mat_set(&A_Tzz_zr,0,0,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,0,1,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,0,1,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,0,1,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,0,1,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,0,1,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,0,1,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,0,2,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,0,2,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,0,2,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,0,2,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,0,2,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,0,2,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,0,3,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,0,3,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,0,3,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,0,3,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,0,3,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,0,3,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[4]; 
  gkyl_mat_set(&A_Txx_xl,0,4,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[4]; 
  gkyl_mat_set(&A_Txx_xr,0,4,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[4]; 
  gkyl_mat_set(&A_Tyy_yl,0,4,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[4]; 
  gkyl_mat_set(&A_Tyy_yr,0,4,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[4]; 
  gkyl_mat_set(&A_Tzz_zl,0,4,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[4]; 
  gkyl_mat_set(&A_Tzz_zr,0,4,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[5]; 
  gkyl_mat_set(&A_Txx_xl,0,5,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[5]; 
  gkyl_mat_set(&A_Txx_xr,0,5,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[5]; 
  gkyl_mat_set(&A_Tyy_yl,0,5,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[5]; 
  gkyl_mat_set(&A_Tyy_yr,0,5,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[5]; 
  gkyl_mat_set(&A_Tzz_zl,0,5,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[5]; 
  gkyl_mat_set(&A_Tzz_zr,0,5,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[6]; 
  gkyl_mat_set(&A_Txx_xl,0,6,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[6]; 
  gkyl_mat_set(&A_Txx_xr,0,6,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[6]; 
  gkyl_mat_set(&A_Tyy_yl,0,6,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[6]; 
  gkyl_mat_set(&A_Tyy_yr,0,6,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[6]; 
  gkyl_mat_set(&A_Tzz_zl,0,6,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[6]; 
  gkyl_mat_set(&A_Tzz_zr,0,6,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[7]; 
  gkyl_mat_set(&A_Txx_xl,0,7,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[7]; 
  gkyl_mat_set(&A_Txx_xr,0,7,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[7]; 
  gkyl_mat_set(&A_Tyy_yl,0,7,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[7]; 
  gkyl_mat_set(&A_Tyy_yr,0,7,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[7]; 
  gkyl_mat_set(&A_Tzz_zl,0,7,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[7]; 
  gkyl_mat_set(&A_Tzz_zr,0,7,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[8]; 
  gkyl_mat_set(&A_Txx_xl,0,8,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[8]; 
  gkyl_mat_set(&A_Txx_xr,0,8,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[8]; 
  gkyl_mat_set(&A_Tyy_yl,0,8,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[8]; 
  gkyl_mat_set(&A_Tyy_yr,0,8,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[8]; 
  gkyl_mat_set(&A_Tzz_zl,0,8,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[8]; 
  gkyl_mat_set(&A_Tzz_zr,0,8,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,1,0,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,1,0,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,1,0,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,1,0,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,1,0,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,1,0,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[4]+0.5*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,1,1,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[4]+0.5*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,1,1,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[4]+0.5*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,1,1,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[4]+0.5*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,1,1,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[4]+0.5*rho_zl[0]; 
  gkyl_mat_set(&A_Tzz_zl,1,1,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[4]+0.5*rho_zr[0]; 
  gkyl_mat_set(&A_Tzz_zr,1,1,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,1,2,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,1,2,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,1,2,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,1,2,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,1,2,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,1,2,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[6]+0.5*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,1,3,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[6]+0.5*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,1,3,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[6]+0.5*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,1,3,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[6]+0.5*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,1,3,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[6]+0.5*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,1,3,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[6]+0.5*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,1,3,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,1,4,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,1,4,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,1,4,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,1,4,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,1,4,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,1,4,temp_rho_zr); 
 
  temp_rho_xl = 0.5000000000000001*rho_xl[7]; 
  gkyl_mat_set(&A_Txx_xl,1,5,temp_rho_xl); 
 
  temp_rho_xr = 0.5000000000000001*rho_xr[7]; 
  gkyl_mat_set(&A_Txx_xr,1,5,temp_rho_xr); 
 
  temp_rho_yl = 0.5000000000000001*rho_yl[7]; 
  gkyl_mat_set(&A_Tyy_yl,1,5,temp_rho_yl); 
 
  temp_rho_yr = 0.5000000000000001*rho_yr[7]; 
  gkyl_mat_set(&A_Tyy_yr,1,5,temp_rho_yr); 
 
  temp_rho_zl = 0.5000000000000001*rho_zl[7]; 
  gkyl_mat_set(&A_Tzz_zl,1,5,temp_rho_zl); 
 
  temp_rho_zr = 0.5000000000000001*rho_zr[7]; 
  gkyl_mat_set(&A_Tzz_zr,1,5,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,1,6,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,1,6,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,1,6,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,1,6,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,1,6,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,1,6,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[8]+0.5000000000000001*rho_xl[5]; 
  gkyl_mat_set(&A_Txx_xl,1,7,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[8]+0.5000000000000001*rho_xr[5]; 
  gkyl_mat_set(&A_Txx_xr,1,7,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[8]+0.5000000000000001*rho_yl[5]; 
  gkyl_mat_set(&A_Tyy_yl,1,7,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[8]+0.5000000000000001*rho_yr[5]; 
  gkyl_mat_set(&A_Tyy_yr,1,7,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[8]+0.5000000000000001*rho_zl[5]; 
  gkyl_mat_set(&A_Tzz_zl,1,7,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[8]+0.5000000000000001*rho_zr[5]; 
  gkyl_mat_set(&A_Tzz_zr,1,7,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[7]; 
  gkyl_mat_set(&A_Txx_xl,1,8,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[7]; 
  gkyl_mat_set(&A_Txx_xr,1,8,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[7]; 
  gkyl_mat_set(&A_Tyy_yl,1,8,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[7]; 
  gkyl_mat_set(&A_Tyy_yr,1,8,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[7]; 
  gkyl_mat_set(&A_Tzz_zl,1,8,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[7]; 
  gkyl_mat_set(&A_Tzz_zr,1,8,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,2,0,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,2,0,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,2,0,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,2,0,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,2,0,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,2,0,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,2,1,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,2,1,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,2,1,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,2,1,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,2,1,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,2,1,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[5]+0.5*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,2,2,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[5]+0.5*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,2,2,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[5]+0.5*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,2,2,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[5]+0.5*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,2,2,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[5]+0.5*rho_zl[0]; 
  gkyl_mat_set(&A_Tzz_zl,2,2,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[5]+0.5*rho_zr[0]; 
  gkyl_mat_set(&A_Tzz_zr,2,2,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[7]+0.5*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,2,3,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[7]+0.5*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,2,3,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[7]+0.5*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,2,3,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[7]+0.5*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,2,3,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[7]+0.5*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,2,3,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[7]+0.5*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,2,3,temp_rho_zr); 
 
  temp_rho_xl = 0.5000000000000001*rho_xl[6]; 
  gkyl_mat_set(&A_Txx_xl,2,4,temp_rho_xl); 
 
  temp_rho_xr = 0.5000000000000001*rho_xr[6]; 
  gkyl_mat_set(&A_Txx_xr,2,4,temp_rho_xr); 
 
  temp_rho_yl = 0.5000000000000001*rho_yl[6]; 
  gkyl_mat_set(&A_Tyy_yl,2,4,temp_rho_yl); 
 
  temp_rho_yr = 0.5000000000000001*rho_yr[6]; 
  gkyl_mat_set(&A_Tyy_yr,2,4,temp_rho_yr); 
 
  temp_rho_zl = 0.5000000000000001*rho_zl[6]; 
  gkyl_mat_set(&A_Tzz_zl,2,4,temp_rho_zl); 
 
  temp_rho_zr = 0.5000000000000001*rho_zr[6]; 
  gkyl_mat_set(&A_Tzz_zr,2,4,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,2,5,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,2,5,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,2,5,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,2,5,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,2,5,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,2,5,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[8]+0.5000000000000001*rho_xl[4]; 
  gkyl_mat_set(&A_Txx_xl,2,6,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[8]+0.5000000000000001*rho_xr[4]; 
  gkyl_mat_set(&A_Txx_xr,2,6,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[8]+0.5000000000000001*rho_yl[4]; 
  gkyl_mat_set(&A_Tyy_yl,2,6,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[8]+0.5000000000000001*rho_yr[4]; 
  gkyl_mat_set(&A_Tyy_yr,2,6,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[8]+0.5000000000000001*rho_zl[4]; 
  gkyl_mat_set(&A_Tzz_zl,2,6,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[8]+0.5000000000000001*rho_zr[4]; 
  gkyl_mat_set(&A_Tzz_zr,2,6,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,2,7,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,2,7,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,2,7,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,2,7,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,2,7,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,2,7,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[6]; 
  gkyl_mat_set(&A_Txx_xl,2,8,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[6]; 
  gkyl_mat_set(&A_Txx_xr,2,8,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[6]; 
  gkyl_mat_set(&A_Tyy_yl,2,8,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[6]; 
  gkyl_mat_set(&A_Tyy_yr,2,8,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[6]; 
  gkyl_mat_set(&A_Tzz_zl,2,8,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[6]; 
  gkyl_mat_set(&A_Tzz_zr,2,8,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,3,0,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,3,0,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,3,0,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,3,0,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,3,0,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,3,0,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[6]+0.5*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,3,1,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[6]+0.5*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,3,1,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[6]+0.5*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,3,1,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[6]+0.5*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,3,1,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[6]+0.5*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,3,1,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[6]+0.5*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,3,1,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[7]+0.5*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,3,2,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[7]+0.5*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,3,2,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[7]+0.5*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,3,2,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[7]+0.5*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,3,2,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[7]+0.5*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,3,2,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[7]+0.5*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,3,2,temp_rho_zr); 
 
  temp_rho_xl = 0.4*rho_xl[8]+0.4472135954999579*rho_xl[5]+0.4472135954999579*rho_xl[4]+0.5*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,3,3,temp_rho_xl); 
 
  temp_rho_xr = 0.4*rho_xr[8]+0.4472135954999579*rho_xr[5]+0.4472135954999579*rho_xr[4]+0.5*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,3,3,temp_rho_xr); 
 
  temp_rho_yl = 0.4*rho_yl[8]+0.4472135954999579*rho_yl[5]+0.4472135954999579*rho_yl[4]+0.5*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,3,3,temp_rho_yl); 
 
  temp_rho_yr = 0.4*rho_yr[8]+0.4472135954999579*rho_yr[5]+0.4472135954999579*rho_yr[4]+0.5*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,3,3,temp_rho_yr); 
 
  temp_rho_zl = 0.4*rho_zl[8]+0.4472135954999579*rho_zl[5]+0.4472135954999579*rho_zl[4]+0.5*rho_zl[0]; 
  gkyl_mat_set(&A_Tzz_zl,3,3,temp_rho_zl); 
 
  temp_rho_zr = 0.4*rho_zr[8]+0.4472135954999579*rho_zr[5]+0.4472135954999579*rho_zr[4]+0.5*rho_zr[0]; 
  gkyl_mat_set(&A_Tzz_zr,3,3,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,3,4,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,3,4,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,3,4,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,3,4,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,3,4,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,3,4,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,3,5,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,3,5,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,3,5,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,3,5,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,3,5,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,3,5,temp_rho_zr); 
 
  temp_rho_xl = 0.4*rho_xl[7]+0.447213595499958*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,3,6,temp_rho_xl); 
 
  temp_rho_xr = 0.4*rho_xr[7]+0.447213595499958*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,3,6,temp_rho_xr); 
 
  temp_rho_yl = 0.4*rho_yl[7]+0.447213595499958*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,3,6,temp_rho_yl); 
 
  temp_rho_yr = 0.4*rho_yr[7]+0.447213595499958*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,3,6,temp_rho_yr); 
 
  temp_rho_zl = 0.4*rho_zl[7]+0.447213595499958*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,3,6,temp_rho_zl); 
 
  temp_rho_zr = 0.4*rho_zr[7]+0.447213595499958*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,3,6,temp_rho_zr); 
 
  temp_rho_xl = 0.4*rho_xl[6]+0.447213595499958*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,3,7,temp_rho_xl); 
 
  temp_rho_xr = 0.4*rho_xr[6]+0.447213595499958*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,3,7,temp_rho_xr); 
 
  temp_rho_yl = 0.4*rho_yl[6]+0.447213595499958*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,3,7,temp_rho_yl); 
 
  temp_rho_yr = 0.4*rho_yr[6]+0.447213595499958*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,3,7,temp_rho_yr); 
 
  temp_rho_zl = 0.4*rho_zl[6]+0.447213595499958*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,3,7,temp_rho_zl); 
 
  temp_rho_zr = 0.4*rho_zr[6]+0.447213595499958*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,3,7,temp_rho_zr); 
 
  temp_rho_xl = 0.4*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,3,8,temp_rho_xl); 
 
  temp_rho_xr = 0.4*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,3,8,temp_rho_xr); 
 
  temp_rho_yl = 0.4*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,3,8,temp_rho_yl); 
 
  temp_rho_yr = 0.4*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,3,8,temp_rho_yr); 
 
  temp_rho_zl = 0.4*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,3,8,temp_rho_zl); 
 
  temp_rho_zr = 0.4*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,3,8,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[4]; 
  gkyl_mat_set(&A_Txx_xl,4,0,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[4]; 
  gkyl_mat_set(&A_Txx_xr,4,0,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[4]; 
  gkyl_mat_set(&A_Tyy_yl,4,0,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[4]; 
  gkyl_mat_set(&A_Tyy_yr,4,0,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[4]; 
  gkyl_mat_set(&A_Tzz_zl,4,0,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[4]; 
  gkyl_mat_set(&A_Tzz_zr,4,0,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,4,1,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,4,1,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,4,1,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,4,1,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,4,1,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,4,1,temp_rho_zr); 
 
  temp_rho_xl = 0.5000000000000001*rho_xl[6]; 
  gkyl_mat_set(&A_Txx_xl,4,2,temp_rho_xl); 
 
  temp_rho_xr = 0.5000000000000001*rho_xr[6]; 
  gkyl_mat_set(&A_Txx_xr,4,2,temp_rho_xr); 
 
  temp_rho_yl = 0.5000000000000001*rho_yl[6]; 
  gkyl_mat_set(&A_Tyy_yl,4,2,temp_rho_yl); 
 
  temp_rho_yr = 0.5000000000000001*rho_yr[6]; 
  gkyl_mat_set(&A_Tyy_yr,4,2,temp_rho_yr); 
 
  temp_rho_zl = 0.5000000000000001*rho_zl[6]; 
  gkyl_mat_set(&A_Tzz_zl,4,2,temp_rho_zl); 
 
  temp_rho_zr = 0.5000000000000001*rho_zr[6]; 
  gkyl_mat_set(&A_Tzz_zr,4,2,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,4,3,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,4,3,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,4,3,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,4,3,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,4,3,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,4,3,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[4]+0.5*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,4,4,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[4]+0.5*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,4,4,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[4]+0.5*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,4,4,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[4]+0.5*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,4,4,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[4]+0.5*rho_zl[0]; 
  gkyl_mat_set(&A_Tzz_zl,4,4,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[4]+0.5*rho_zr[0]; 
  gkyl_mat_set(&A_Tzz_zr,4,4,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[8]; 
  gkyl_mat_set(&A_Txx_xl,4,5,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[8]; 
  gkyl_mat_set(&A_Txx_xr,4,5,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[8]; 
  gkyl_mat_set(&A_Tyy_yl,4,5,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[8]; 
  gkyl_mat_set(&A_Tyy_yr,4,5,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[8]; 
  gkyl_mat_set(&A_Tzz_zl,4,5,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[8]; 
  gkyl_mat_set(&A_Tzz_zr,4,5,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[6]+0.5000000000000001*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,4,6,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[6]+0.5000000000000001*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,4,6,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[6]+0.5000000000000001*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,4,6,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[6]+0.5000000000000001*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,4,6,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[6]+0.5000000000000001*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,4,6,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[6]+0.5000000000000001*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,4,6,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[7]; 
  gkyl_mat_set(&A_Txx_xl,4,7,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[7]; 
  gkyl_mat_set(&A_Txx_xr,4,7,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[7]; 
  gkyl_mat_set(&A_Tyy_yl,4,7,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[7]; 
  gkyl_mat_set(&A_Tyy_yr,4,7,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[7]; 
  gkyl_mat_set(&A_Tzz_zl,4,7,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[7]; 
  gkyl_mat_set(&A_Tzz_zr,4,7,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[8]+0.5*rho_xl[5]; 
  gkyl_mat_set(&A_Txx_xl,4,8,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[8]+0.5*rho_xr[5]; 
  gkyl_mat_set(&A_Txx_xr,4,8,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[8]+0.5*rho_yl[5]; 
  gkyl_mat_set(&A_Tyy_yl,4,8,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[8]+0.5*rho_yr[5]; 
  gkyl_mat_set(&A_Tyy_yr,4,8,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[8]+0.5*rho_zl[5]; 
  gkyl_mat_set(&A_Tzz_zl,4,8,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[8]+0.5*rho_zr[5]; 
  gkyl_mat_set(&A_Tzz_zr,4,8,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[5]; 
  gkyl_mat_set(&A_Txx_xl,5,0,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[5]; 
  gkyl_mat_set(&A_Txx_xr,5,0,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[5]; 
  gkyl_mat_set(&A_Tyy_yl,5,0,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[5]; 
  gkyl_mat_set(&A_Tyy_yr,5,0,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[5]; 
  gkyl_mat_set(&A_Tzz_zl,5,0,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[5]; 
  gkyl_mat_set(&A_Tzz_zr,5,0,temp_rho_zr); 
 
  temp_rho_xl = 0.5000000000000001*rho_xl[7]; 
  gkyl_mat_set(&A_Txx_xl,5,1,temp_rho_xl); 
 
  temp_rho_xr = 0.5000000000000001*rho_xr[7]; 
  gkyl_mat_set(&A_Txx_xr,5,1,temp_rho_xr); 
 
  temp_rho_yl = 0.5000000000000001*rho_yl[7]; 
  gkyl_mat_set(&A_Tyy_yl,5,1,temp_rho_yl); 
 
  temp_rho_yr = 0.5000000000000001*rho_yr[7]; 
  gkyl_mat_set(&A_Tyy_yr,5,1,temp_rho_yr); 
 
  temp_rho_zl = 0.5000000000000001*rho_zl[7]; 
  gkyl_mat_set(&A_Tzz_zl,5,1,temp_rho_zl); 
 
  temp_rho_zr = 0.5000000000000001*rho_zr[7]; 
  gkyl_mat_set(&A_Tzz_zr,5,1,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,5,2,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,5,2,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,5,2,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,5,2,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,5,2,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,5,2,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,5,3,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,5,3,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,5,3,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,5,3,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,5,3,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,5,3,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[8]; 
  gkyl_mat_set(&A_Txx_xl,5,4,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[8]; 
  gkyl_mat_set(&A_Txx_xr,5,4,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[8]; 
  gkyl_mat_set(&A_Tyy_yl,5,4,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[8]; 
  gkyl_mat_set(&A_Tyy_yr,5,4,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[8]; 
  gkyl_mat_set(&A_Tzz_zl,5,4,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[8]; 
  gkyl_mat_set(&A_Tzz_zr,5,4,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[5]+0.5*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,5,5,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[5]+0.5*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,5,5,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[5]+0.5*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,5,5,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[5]+0.5*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,5,5,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[5]+0.5*rho_zl[0]; 
  gkyl_mat_set(&A_Tzz_zl,5,5,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[5]+0.5*rho_zr[0]; 
  gkyl_mat_set(&A_Tzz_zr,5,5,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[6]; 
  gkyl_mat_set(&A_Txx_xl,5,6,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[6]; 
  gkyl_mat_set(&A_Txx_xr,5,6,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[6]; 
  gkyl_mat_set(&A_Tyy_yl,5,6,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[6]; 
  gkyl_mat_set(&A_Tyy_yr,5,6,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[6]; 
  gkyl_mat_set(&A_Tzz_zl,5,6,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[6]; 
  gkyl_mat_set(&A_Tzz_zr,5,6,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[7]+0.5000000000000001*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,5,7,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[7]+0.5000000000000001*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,5,7,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[7]+0.5000000000000001*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,5,7,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[7]+0.5000000000000001*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,5,7,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[7]+0.5000000000000001*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,5,7,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[7]+0.5000000000000001*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,5,7,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[8]+0.5*rho_xl[4]; 
  gkyl_mat_set(&A_Txx_xl,5,8,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[8]+0.5*rho_xr[4]; 
  gkyl_mat_set(&A_Txx_xr,5,8,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[8]+0.5*rho_yl[4]; 
  gkyl_mat_set(&A_Tyy_yl,5,8,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[8]+0.5*rho_yr[4]; 
  gkyl_mat_set(&A_Tyy_yr,5,8,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[8]+0.5*rho_zl[4]; 
  gkyl_mat_set(&A_Tzz_zl,5,8,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[8]+0.5*rho_zr[4]; 
  gkyl_mat_set(&A_Tzz_zr,5,8,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[6]; 
  gkyl_mat_set(&A_Txx_xl,6,0,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[6]; 
  gkyl_mat_set(&A_Txx_xr,6,0,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[6]; 
  gkyl_mat_set(&A_Tyy_yl,6,0,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[6]; 
  gkyl_mat_set(&A_Tyy_yr,6,0,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[6]; 
  gkyl_mat_set(&A_Tzz_zl,6,0,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[6]; 
  gkyl_mat_set(&A_Tzz_zr,6,0,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,6,1,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,6,1,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,6,1,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,6,1,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,6,1,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,6,1,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[8]+0.5000000000000001*rho_xl[4]; 
  gkyl_mat_set(&A_Txx_xl,6,2,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[8]+0.5000000000000001*rho_xr[4]; 
  gkyl_mat_set(&A_Txx_xr,6,2,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[8]+0.5000000000000001*rho_yl[4]; 
  gkyl_mat_set(&A_Tyy_yl,6,2,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[8]+0.5000000000000001*rho_yr[4]; 
  gkyl_mat_set(&A_Tyy_yr,6,2,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[8]+0.5000000000000001*rho_zl[4]; 
  gkyl_mat_set(&A_Tzz_zl,6,2,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[8]+0.5000000000000001*rho_zr[4]; 
  gkyl_mat_set(&A_Tzz_zr,6,2,temp_rho_zr); 
 
  temp_rho_xl = 0.4*rho_xl[7]+0.447213595499958*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,6,3,temp_rho_xl); 
 
  temp_rho_xr = 0.4*rho_xr[7]+0.447213595499958*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,6,3,temp_rho_xr); 
 
  temp_rho_yl = 0.4*rho_yl[7]+0.447213595499958*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,6,3,temp_rho_yl); 
 
  temp_rho_yr = 0.4*rho_yr[7]+0.447213595499958*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,6,3,temp_rho_yr); 
 
  temp_rho_zl = 0.4*rho_zl[7]+0.447213595499958*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,6,3,temp_rho_zl); 
 
  temp_rho_zr = 0.4*rho_zr[7]+0.447213595499958*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,6,3,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[6]+0.5000000000000001*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,6,4,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[6]+0.5000000000000001*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,6,4,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[6]+0.5000000000000001*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,6,4,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[6]+0.5000000000000001*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,6,4,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[6]+0.5000000000000001*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,6,4,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[6]+0.5000000000000001*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,6,4,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[6]; 
  gkyl_mat_set(&A_Txx_xl,6,5,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[6]; 
  gkyl_mat_set(&A_Txx_xr,6,5,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[6]; 
  gkyl_mat_set(&A_Tyy_yl,6,5,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[6]; 
  gkyl_mat_set(&A_Tyy_yr,6,5,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[6]; 
  gkyl_mat_set(&A_Tzz_zl,6,5,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[6]; 
  gkyl_mat_set(&A_Tzz_zr,6,5,temp_rho_zr); 
 
  temp_rho_xl = 0.2857142857142857*rho_xl[8]+0.4472135954999579*rho_xl[5]+0.31943828249997*rho_xl[4]+0.5*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,6,6,temp_rho_xl); 
 
  temp_rho_xr = 0.2857142857142857*rho_xr[8]+0.4472135954999579*rho_xr[5]+0.31943828249997*rho_xr[4]+0.5*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,6,6,temp_rho_xr); 
 
  temp_rho_yl = 0.2857142857142857*rho_yl[8]+0.4472135954999579*rho_yl[5]+0.31943828249997*rho_yl[4]+0.5*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,6,6,temp_rho_yl); 
 
  temp_rho_yr = 0.2857142857142857*rho_yr[8]+0.4472135954999579*rho_yr[5]+0.31943828249997*rho_yr[4]+0.5*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,6,6,temp_rho_yr); 
 
  temp_rho_zl = 0.2857142857142857*rho_zl[8]+0.4472135954999579*rho_zl[5]+0.31943828249997*rho_zl[4]+0.5*rho_zl[0]; 
  gkyl_mat_set(&A_Tzz_zl,6,6,temp_rho_zl); 
 
  temp_rho_zr = 0.2857142857142857*rho_zr[8]+0.4472135954999579*rho_zr[5]+0.31943828249997*rho_zr[4]+0.5*rho_zr[0]; 
  gkyl_mat_set(&A_Tzz_zr,6,6,temp_rho_zr); 
 
  temp_rho_xl = 0.4*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,6,7,temp_rho_xl); 
 
  temp_rho_xr = 0.4*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,6,7,temp_rho_xr); 
 
  temp_rho_yl = 0.4*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,6,7,temp_rho_yl); 
 
  temp_rho_yr = 0.4*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,6,7,temp_rho_yr); 
 
  temp_rho_zl = 0.4*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,6,7,temp_rho_zl); 
 
  temp_rho_zr = 0.4*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,6,7,temp_rho_zr); 
 
  temp_rho_xl = 0.2857142857142857*rho_xl[6]+0.447213595499958*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,6,8,temp_rho_xl); 
 
  temp_rho_xr = 0.2857142857142857*rho_xr[6]+0.447213595499958*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,6,8,temp_rho_xr); 
 
  temp_rho_yl = 0.2857142857142857*rho_yl[6]+0.447213595499958*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,6,8,temp_rho_yl); 
 
  temp_rho_yr = 0.2857142857142857*rho_yr[6]+0.447213595499958*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,6,8,temp_rho_yr); 
 
  temp_rho_zl = 0.2857142857142857*rho_zl[6]+0.447213595499958*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,6,8,temp_rho_zl); 
 
  temp_rho_zr = 0.2857142857142857*rho_zr[6]+0.447213595499958*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,6,8,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[7]; 
  gkyl_mat_set(&A_Txx_xl,7,0,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[7]; 
  gkyl_mat_set(&A_Txx_xr,7,0,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[7]; 
  gkyl_mat_set(&A_Tyy_yl,7,0,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[7]; 
  gkyl_mat_set(&A_Tyy_yr,7,0,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[7]; 
  gkyl_mat_set(&A_Tzz_zl,7,0,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[7]; 
  gkyl_mat_set(&A_Tzz_zr,7,0,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[8]+0.5000000000000001*rho_xl[5]; 
  gkyl_mat_set(&A_Txx_xl,7,1,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[8]+0.5000000000000001*rho_xr[5]; 
  gkyl_mat_set(&A_Txx_xr,7,1,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[8]+0.5000000000000001*rho_yl[5]; 
  gkyl_mat_set(&A_Tyy_yl,7,1,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[8]+0.5000000000000001*rho_yr[5]; 
  gkyl_mat_set(&A_Tyy_yr,7,1,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[8]+0.5000000000000001*rho_zl[5]; 
  gkyl_mat_set(&A_Tzz_zl,7,1,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[8]+0.5000000000000001*rho_zr[5]; 
  gkyl_mat_set(&A_Tzz_zr,7,1,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,7,2,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,7,2,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,7,2,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,7,2,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,7,2,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,7,2,temp_rho_zr); 
 
  temp_rho_xl = 0.4*rho_xl[6]+0.447213595499958*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,7,3,temp_rho_xl); 
 
  temp_rho_xr = 0.4*rho_xr[6]+0.447213595499958*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,7,3,temp_rho_xr); 
 
  temp_rho_yl = 0.4*rho_yl[6]+0.447213595499958*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,7,3,temp_rho_yl); 
 
  temp_rho_yr = 0.4*rho_yr[6]+0.447213595499958*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,7,3,temp_rho_yr); 
 
  temp_rho_zl = 0.4*rho_zl[6]+0.447213595499958*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,7,3,temp_rho_zl); 
 
  temp_rho_zr = 0.4*rho_zr[6]+0.447213595499958*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,7,3,temp_rho_zr); 
 
  temp_rho_xl = 0.4472135954999579*rho_xl[7]; 
  gkyl_mat_set(&A_Txx_xl,7,4,temp_rho_xl); 
 
  temp_rho_xr = 0.4472135954999579*rho_xr[7]; 
  gkyl_mat_set(&A_Txx_xr,7,4,temp_rho_xr); 
 
  temp_rho_yl = 0.4472135954999579*rho_yl[7]; 
  gkyl_mat_set(&A_Tyy_yl,7,4,temp_rho_yl); 
 
  temp_rho_yr = 0.4472135954999579*rho_yr[7]; 
  gkyl_mat_set(&A_Tyy_yr,7,4,temp_rho_yr); 
 
  temp_rho_zl = 0.4472135954999579*rho_zl[7]; 
  gkyl_mat_set(&A_Tzz_zl,7,4,temp_rho_zl); 
 
  temp_rho_zr = 0.4472135954999579*rho_zr[7]; 
  gkyl_mat_set(&A_Tzz_zr,7,4,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[7]+0.5000000000000001*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,7,5,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[7]+0.5000000000000001*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,7,5,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[7]+0.5000000000000001*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,7,5,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[7]+0.5000000000000001*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,7,5,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[7]+0.5000000000000001*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,7,5,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[7]+0.5000000000000001*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,7,5,temp_rho_zr); 
 
  temp_rho_xl = 0.4*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,7,6,temp_rho_xl); 
 
  temp_rho_xr = 0.4*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,7,6,temp_rho_xr); 
 
  temp_rho_yl = 0.4*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,7,6,temp_rho_yl); 
 
  temp_rho_yr = 0.4*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,7,6,temp_rho_yr); 
 
  temp_rho_zl = 0.4*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,7,6,temp_rho_zl); 
 
  temp_rho_zr = 0.4*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,7,6,temp_rho_zr); 
 
  temp_rho_xl = 0.2857142857142857*rho_xl[8]+0.31943828249997*rho_xl[5]+0.4472135954999579*rho_xl[4]+0.5*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,7,7,temp_rho_xl); 
 
  temp_rho_xr = 0.2857142857142857*rho_xr[8]+0.31943828249997*rho_xr[5]+0.4472135954999579*rho_xr[4]+0.5*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,7,7,temp_rho_xr); 
 
  temp_rho_yl = 0.2857142857142857*rho_yl[8]+0.31943828249997*rho_yl[5]+0.4472135954999579*rho_yl[4]+0.5*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,7,7,temp_rho_yl); 
 
  temp_rho_yr = 0.2857142857142857*rho_yr[8]+0.31943828249997*rho_yr[5]+0.4472135954999579*rho_yr[4]+0.5*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,7,7,temp_rho_yr); 
 
  temp_rho_zl = 0.2857142857142857*rho_zl[8]+0.31943828249997*rho_zl[5]+0.4472135954999579*rho_zl[4]+0.5*rho_zl[0]; 
  gkyl_mat_set(&A_Tzz_zl,7,7,temp_rho_zl); 
 
  temp_rho_zr = 0.2857142857142857*rho_zr[8]+0.31943828249997*rho_zr[5]+0.4472135954999579*rho_zr[4]+0.5*rho_zr[0]; 
  gkyl_mat_set(&A_Tzz_zr,7,7,temp_rho_zr); 
 
  temp_rho_xl = 0.2857142857142857*rho_xl[7]+0.447213595499958*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,7,8,temp_rho_xl); 
 
  temp_rho_xr = 0.2857142857142857*rho_xr[7]+0.447213595499958*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,7,8,temp_rho_xr); 
 
  temp_rho_yl = 0.2857142857142857*rho_yl[7]+0.447213595499958*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,7,8,temp_rho_yl); 
 
  temp_rho_yr = 0.2857142857142857*rho_yr[7]+0.447213595499958*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,7,8,temp_rho_yr); 
 
  temp_rho_zl = 0.2857142857142857*rho_zl[7]+0.447213595499958*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,7,8,temp_rho_zl); 
 
  temp_rho_zr = 0.2857142857142857*rho_zr[7]+0.447213595499958*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,7,8,temp_rho_zr); 
 
  temp_rho_xl = 0.5*rho_xl[8]; 
  gkyl_mat_set(&A_Txx_xl,8,0,temp_rho_xl); 
 
  temp_rho_xr = 0.5*rho_xr[8]; 
  gkyl_mat_set(&A_Txx_xr,8,0,temp_rho_xr); 
 
  temp_rho_yl = 0.5*rho_yl[8]; 
  gkyl_mat_set(&A_Tyy_yl,8,0,temp_rho_yl); 
 
  temp_rho_yr = 0.5*rho_yr[8]; 
  gkyl_mat_set(&A_Tyy_yr,8,0,temp_rho_yr); 
 
  temp_rho_zl = 0.5*rho_zl[8]; 
  gkyl_mat_set(&A_Tzz_zl,8,0,temp_rho_zl); 
 
  temp_rho_zr = 0.5*rho_zr[8]; 
  gkyl_mat_set(&A_Tzz_zr,8,0,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[7]; 
  gkyl_mat_set(&A_Txx_xl,8,1,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[7]; 
  gkyl_mat_set(&A_Txx_xr,8,1,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[7]; 
  gkyl_mat_set(&A_Tyy_yl,8,1,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[7]; 
  gkyl_mat_set(&A_Tyy_yr,8,1,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[7]; 
  gkyl_mat_set(&A_Tzz_zl,8,1,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[7]; 
  gkyl_mat_set(&A_Tzz_zr,8,1,temp_rho_zr); 
 
  temp_rho_xl = 0.447213595499958*rho_xl[6]; 
  gkyl_mat_set(&A_Txx_xl,8,2,temp_rho_xl); 
 
  temp_rho_xr = 0.447213595499958*rho_xr[6]; 
  gkyl_mat_set(&A_Txx_xr,8,2,temp_rho_xr); 
 
  temp_rho_yl = 0.447213595499958*rho_yl[6]; 
  gkyl_mat_set(&A_Tyy_yl,8,2,temp_rho_yl); 
 
  temp_rho_yr = 0.447213595499958*rho_yr[6]; 
  gkyl_mat_set(&A_Tyy_yr,8,2,temp_rho_yr); 
 
  temp_rho_zl = 0.447213595499958*rho_zl[6]; 
  gkyl_mat_set(&A_Tzz_zl,8,2,temp_rho_zl); 
 
  temp_rho_zr = 0.447213595499958*rho_zr[6]; 
  gkyl_mat_set(&A_Tzz_zr,8,2,temp_rho_zr); 
 
  temp_rho_xl = 0.4*rho_xl[3]; 
  gkyl_mat_set(&A_Txx_xl,8,3,temp_rho_xl); 
 
  temp_rho_xr = 0.4*rho_xr[3]; 
  gkyl_mat_set(&A_Txx_xr,8,3,temp_rho_xr); 
 
  temp_rho_yl = 0.4*rho_yl[3]; 
  gkyl_mat_set(&A_Tyy_yl,8,3,temp_rho_yl); 
 
  temp_rho_yr = 0.4*rho_yr[3]; 
  gkyl_mat_set(&A_Tyy_yr,8,3,temp_rho_yr); 
 
  temp_rho_zl = 0.4*rho_zl[3]; 
  gkyl_mat_set(&A_Tzz_zl,8,3,temp_rho_zl); 
 
  temp_rho_zr = 0.4*rho_zr[3]; 
  gkyl_mat_set(&A_Tzz_zr,8,3,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[8]+0.5*rho_xl[5]; 
  gkyl_mat_set(&A_Txx_xl,8,4,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[8]+0.5*rho_xr[5]; 
  gkyl_mat_set(&A_Txx_xr,8,4,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[8]+0.5*rho_yl[5]; 
  gkyl_mat_set(&A_Tyy_yl,8,4,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[8]+0.5*rho_yr[5]; 
  gkyl_mat_set(&A_Tyy_yr,8,4,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[8]+0.5*rho_zl[5]; 
  gkyl_mat_set(&A_Tzz_zl,8,4,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[8]+0.5*rho_zr[5]; 
  gkyl_mat_set(&A_Tzz_zr,8,4,temp_rho_zr); 
 
  temp_rho_xl = 0.31943828249997*rho_xl[8]+0.5*rho_xl[4]; 
  gkyl_mat_set(&A_Txx_xl,8,5,temp_rho_xl); 
 
  temp_rho_xr = 0.31943828249997*rho_xr[8]+0.5*rho_xr[4]; 
  gkyl_mat_set(&A_Txx_xr,8,5,temp_rho_xr); 
 
  temp_rho_yl = 0.31943828249997*rho_yl[8]+0.5*rho_yl[4]; 
  gkyl_mat_set(&A_Tyy_yl,8,5,temp_rho_yl); 
 
  temp_rho_yr = 0.31943828249997*rho_yr[8]+0.5*rho_yr[4]; 
  gkyl_mat_set(&A_Tyy_yr,8,5,temp_rho_yr); 
 
  temp_rho_zl = 0.31943828249997*rho_zl[8]+0.5*rho_zl[4]; 
  gkyl_mat_set(&A_Tzz_zl,8,5,temp_rho_zl); 
 
  temp_rho_zr = 0.31943828249997*rho_zr[8]+0.5*rho_zr[4]; 
  gkyl_mat_set(&A_Tzz_zr,8,5,temp_rho_zr); 
 
  temp_rho_xl = 0.2857142857142857*rho_xl[6]+0.447213595499958*rho_xl[2]; 
  gkyl_mat_set(&A_Txx_xl,8,6,temp_rho_xl); 
 
  temp_rho_xr = 0.2857142857142857*rho_xr[6]+0.447213595499958*rho_xr[2]; 
  gkyl_mat_set(&A_Txx_xr,8,6,temp_rho_xr); 
 
  temp_rho_yl = 0.2857142857142857*rho_yl[6]+0.447213595499958*rho_yl[2]; 
  gkyl_mat_set(&A_Tyy_yl,8,6,temp_rho_yl); 
 
  temp_rho_yr = 0.2857142857142857*rho_yr[6]+0.447213595499958*rho_yr[2]; 
  gkyl_mat_set(&A_Tyy_yr,8,6,temp_rho_yr); 
 
  temp_rho_zl = 0.2857142857142857*rho_zl[6]+0.447213595499958*rho_zl[2]; 
  gkyl_mat_set(&A_Tzz_zl,8,6,temp_rho_zl); 
 
  temp_rho_zr = 0.2857142857142857*rho_zr[6]+0.447213595499958*rho_zr[2]; 
  gkyl_mat_set(&A_Tzz_zr,8,6,temp_rho_zr); 
 
  temp_rho_xl = 0.2857142857142857*rho_xl[7]+0.447213595499958*rho_xl[1]; 
  gkyl_mat_set(&A_Txx_xl,8,7,temp_rho_xl); 
 
  temp_rho_xr = 0.2857142857142857*rho_xr[7]+0.447213595499958*rho_xr[1]; 
  gkyl_mat_set(&A_Txx_xr,8,7,temp_rho_xr); 
 
  temp_rho_yl = 0.2857142857142857*rho_yl[7]+0.447213595499958*rho_yl[1]; 
  gkyl_mat_set(&A_Tyy_yl,8,7,temp_rho_yl); 
 
  temp_rho_yr = 0.2857142857142857*rho_yr[7]+0.447213595499958*rho_yr[1]; 
  gkyl_mat_set(&A_Tyy_yr,8,7,temp_rho_yr); 
 
  temp_rho_zl = 0.2857142857142857*rho_zl[7]+0.447213595499958*rho_zl[1]; 
  gkyl_mat_set(&A_Tzz_zl,8,7,temp_rho_zl); 
 
  temp_rho_zr = 0.2857142857142857*rho_zr[7]+0.447213595499958*rho_zr[1]; 
  gkyl_mat_set(&A_Tzz_zr,8,7,temp_rho_zr); 
 
  temp_rho_xl = 0.2040816326530612*rho_xl[8]+0.31943828249997*rho_xl[5]+0.31943828249997*rho_xl[4]+0.5*rho_xl[0]; 
  gkyl_mat_set(&A_Txx_xl,8,8,temp_rho_xl); 
 
  temp_rho_xr = 0.2040816326530612*rho_xr[8]+0.31943828249997*rho_xr[5]+0.31943828249997*rho_xr[4]+0.5*rho_xr[0]; 
  gkyl_mat_set(&A_Txx_xr,8,8,temp_rho_xr); 
 
  temp_rho_yl = 0.2040816326530612*rho_yl[8]+0.31943828249997*rho_yl[5]+0.31943828249997*rho_yl[4]+0.5*rho_yl[0]; 
  gkyl_mat_set(&A_Tyy_yl,8,8,temp_rho_yl); 
 
  temp_rho_yr = 0.2040816326530612*rho_yr[8]+0.31943828249997*rho_yr[5]+0.31943828249997*rho_yr[4]+0.5*rho_yr[0]; 
  gkyl_mat_set(&A_Tyy_yr,8,8,temp_rho_yr); 
 
  temp_rho_zl = 0.2040816326530612*rho_zl[8]+0.31943828249997*rho_zl[5]+0.31943828249997*rho_zl[4]+0.5*rho_zl[0]; 
  gkyl_mat_set(&A_Tzz_zl,8,8,temp_rho_zl); 
 
  temp_rho_zr = 0.2040816326530612*rho_zr[8]+0.31943828249997*rho_zr[5]+0.31943828249997*rho_zr[4]+0.5*rho_zr[0]; 
  gkyl_mat_set(&A_Tzz_zr,8,8,temp_rho_zr); 
 
} 
