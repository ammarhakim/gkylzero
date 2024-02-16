#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_pressure_2x_tensor_p2(const double *bvar, const double *bvar_surf, const double *vlasov_pkpm_moms, 
    double* GKYL_RESTRICT p_ij, double* GKYL_RESTRICT p_ij_surf) 
{ 
  // bvar:             Input volume expansion of magnetic field unit vector and tensor.
  //                   [bx, by, bz, bxbx, bxby, bxbz, byby, bybz, bzbz] 
  // bvar_surf:        Input surface expansion of magnetic field unit tensor and unit vector. 
  //                   [bx_xl, bx_xr, bxbx_xl, bxbx_xr, bxby_xl, bxby_xr, bxbz_xl, bxbz_xr, 
  //                    by_yl, by_yr, byby_yl, byby_yr, bxby_yl, bxby_yr, bybz_yl, bybz_yr, 
  //                    bz_zl, bz_zr, bzbz_zl, bzbz_zr, bxbz_zl, bxbz_zr, bybz_zl, bybz_zr] 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // p_ij:             Output volume expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.
  // p_ij_surf:        Output surface expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.
  //                   [Pxx_xl, Pxx_xr, Pxy_xl, Pxy_xr, Pxz_xl, Pxz_xr, 
  //                    Pxy_yl, Pxy_yr, Pyy_yl, Pyy_yr, Pyz_yl, Pyz_yr, 
  //                    Pxz_zl, Pxz_zr, Pyz_zl, Pyz_zr, Pzz_zl, Pzz_zr] 

  // Parallel/Perp pressure are first/second component of pkpm moment array and unit tensor are last six components of bvar array.
  const double *p_parallel = &vlasov_pkpm_moms[9]; 
  const double *p_perp = &vlasov_pkpm_moms[18]; 
  const double *bxbx = &bvar[27]; 
  const double *bxby = &bvar[36]; 
  const double *bxbz = &bvar[45]; 
  const double *byby = &bvar[54]; 
  const double *bybz = &bvar[63]; 
  const double *bzbz = &bvar[72]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[9]; 
  double *Pxz = &p_ij[18]; 
  double *Pyy = &p_ij[27]; 
  double *Pyz = &p_ij[36]; 
  double *Pzz = &p_ij[45]; 

  double *Pxx_xl = &p_ij_surf[0]; 
  double *Pxx_xr = &p_ij_surf[3]; 
  double *Pxy_xl = &p_ij_surf[6]; 
  double *Pxy_xr = &p_ij_surf[9]; 
  double *Pxz_xl = &p_ij_surf[12]; 
  double *Pxz_xr = &p_ij_surf[15]; 

  double *Pxy_yl = &p_ij_surf[18]; 
  double *Pxy_yr = &p_ij_surf[21]; 
  double *Pyy_yl = &p_ij_surf[24]; 
  double *Pyy_yr = &p_ij_surf[27]; 
  double *Pyz_yl = &p_ij_surf[30]; 
  double *Pyz_yr = &p_ij_surf[33]; 
 
  double DP[9] = {0.0}; 
  DP[0] = p_parallel[0] - p_perp[0]; 
  DP[1] = p_parallel[1] - p_perp[1]; 
  DP[2] = p_parallel[2] - p_perp[2]; 
  DP[3] = p_parallel[3] - p_perp[3]; 
  DP[4] = p_parallel[4] - p_perp[4]; 
  DP[5] = p_parallel[5] - p_perp[5]; 
  DP[6] = p_parallel[6] - p_perp[6]; 
  DP[7] = p_parallel[7] - p_perp[7]; 
  DP[8] = p_parallel[8] - p_perp[8]; 
  // DP b_i b_j. 
  double DP_bxbx[9] = {0.0}; 
  binop_mul_2d_tensor_p2(DP, bxbx, DP_bxbx); 
 
  double DP_bxby[9] = {0.0}; 
  binop_mul_2d_tensor_p2(DP, bxby, DP_bxby); 
 
  double DP_bxbz[9] = {0.0}; 
  binop_mul_2d_tensor_p2(DP, bxbz, DP_bxbz); 
 
  double DP_byby[9] = {0.0}; 
  binop_mul_2d_tensor_p2(DP, byby, DP_byby); 
 
  double DP_bybz[9] = {0.0}; 
  binop_mul_2d_tensor_p2(DP, bybz, DP_bybz); 
 
  double DP_bzbz[9] = {0.0}; 
  binop_mul_2d_tensor_p2(DP, bzbz, DP_bzbz); 
 
  Pxx[0] = DP_bxbx[0] + p_perp[0]; 
  Pxy[0] = DP_bxby[0]; 
  Pxz[0] = DP_bxbz[0]; 
  Pyy[0] = DP_byby[0] + p_perp[0]; 
  Pyz[0] = DP_bybz[0]; 
  Pzz[0] = DP_bzbz[0] + p_perp[0]; 
  Pxx[1] = DP_bxbx[1] + p_perp[1]; 
  Pxy[1] = DP_bxby[1]; 
  Pxz[1] = DP_bxbz[1]; 
  Pyy[1] = DP_byby[1] + p_perp[1]; 
  Pyz[1] = DP_bybz[1]; 
  Pzz[1] = DP_bzbz[1] + p_perp[1]; 
  Pxx[2] = DP_bxbx[2] + p_perp[2]; 
  Pxy[2] = DP_bxby[2]; 
  Pxz[2] = DP_bxbz[2]; 
  Pyy[2] = DP_byby[2] + p_perp[2]; 
  Pyz[2] = DP_bybz[2]; 
  Pzz[2] = DP_bzbz[2] + p_perp[2]; 
  Pxx[3] = DP_bxbx[3] + p_perp[3]; 
  Pxy[3] = DP_bxby[3]; 
  Pxz[3] = DP_bxbz[3]; 
  Pyy[3] = DP_byby[3] + p_perp[3]; 
  Pyz[3] = DP_bybz[3]; 
  Pzz[3] = DP_bzbz[3] + p_perp[3]; 
  Pxx[4] = DP_bxbx[4] + p_perp[4]; 
  Pxy[4] = DP_bxby[4]; 
  Pxz[4] = DP_bxbz[4]; 
  Pyy[4] = DP_byby[4] + p_perp[4]; 
  Pyz[4] = DP_bybz[4]; 
  Pzz[4] = DP_bzbz[4] + p_perp[4]; 
  Pxx[5] = DP_bxbx[5] + p_perp[5]; 
  Pxy[5] = DP_bxby[5]; 
  Pxz[5] = DP_bxbz[5]; 
  Pyy[5] = DP_byby[5] + p_perp[5]; 
  Pyz[5] = DP_bybz[5]; 
  Pzz[5] = DP_bzbz[5] + p_perp[5]; 
  Pxx[6] = DP_bxbx[6] + p_perp[6]; 
  Pxy[6] = DP_bxby[6]; 
  Pxz[6] = DP_bxbz[6]; 
  Pyy[6] = DP_byby[6] + p_perp[6]; 
  Pyz[6] = DP_bybz[6]; 
  Pzz[6] = DP_bzbz[6] + p_perp[6]; 
  Pxx[7] = DP_bxbx[7] + p_perp[7]; 
  Pxy[7] = DP_bxby[7]; 
  Pxz[7] = DP_bxbz[7]; 
  Pyy[7] = DP_byby[7] + p_perp[7]; 
  Pyz[7] = DP_bybz[7]; 
  Pzz[7] = DP_bzbz[7] + p_perp[7]; 
  Pxx[8] = DP_bxbx[8] + p_perp[8]; 
  Pxy[8] = DP_bxby[8]; 
  Pxz[8] = DP_bxbz[8]; 
  Pyy[8] = DP_byby[8] + p_perp[8]; 
  Pyz[8] = DP_bybz[8]; 
  Pzz[8] = DP_bzbz[8] + p_perp[8]; 
 
  Pxx_xl[0] = 1.58113883008419*Pxx[4]-1.224744871391589*Pxx[1]+0.7071067811865475*Pxx[0]; 
  Pxx_xl[1] = 1.58113883008419*Pxx[6]-1.224744871391589*Pxx[3]+0.7071067811865475*Pxx[2]; 
  Pxx_xl[2] = 1.58113883008419*Pxx[8]-1.224744871391589*Pxx[7]+0.7071067811865475*Pxx[5]; 
  Pxy_xl[0] = 1.58113883008419*Pxy[4]-1.224744871391589*Pxy[1]+0.7071067811865475*Pxy[0]; 
  Pxy_xl[1] = 1.58113883008419*Pxy[6]-1.224744871391589*Pxy[3]+0.7071067811865475*Pxy[2]; 
  Pxy_xl[2] = 1.58113883008419*Pxy[8]-1.224744871391589*Pxy[7]+0.7071067811865475*Pxy[5]; 
  Pxz_xl[0] = 1.58113883008419*Pxz[4]-1.224744871391589*Pxz[1]+0.7071067811865475*Pxz[0]; 
  Pxz_xl[1] = 1.58113883008419*Pxz[6]-1.224744871391589*Pxz[3]+0.7071067811865475*Pxz[2]; 
  Pxz_xl[2] = 1.58113883008419*Pxz[8]-1.224744871391589*Pxz[7]+0.7071067811865475*Pxz[5]; 
 
  Pxx_xr[0] = 1.58113883008419*Pxx[4]+1.224744871391589*Pxx[1]+0.7071067811865475*Pxx[0]; 
  Pxx_xr[1] = 1.58113883008419*Pxx[6]+1.224744871391589*Pxx[3]+0.7071067811865475*Pxx[2]; 
  Pxx_xr[2] = 1.58113883008419*Pxx[8]+1.224744871391589*Pxx[7]+0.7071067811865475*Pxx[5]; 
  Pxy_xr[0] = 1.58113883008419*Pxy[4]+1.224744871391589*Pxy[1]+0.7071067811865475*Pxy[0]; 
  Pxy_xr[1] = 1.58113883008419*Pxy[6]+1.224744871391589*Pxy[3]+0.7071067811865475*Pxy[2]; 
  Pxy_xr[2] = 1.58113883008419*Pxy[8]+1.224744871391589*Pxy[7]+0.7071067811865475*Pxy[5]; 
  Pxz_xr[0] = 1.58113883008419*Pxz[4]+1.224744871391589*Pxz[1]+0.7071067811865475*Pxz[0]; 
  Pxz_xr[1] = 1.58113883008419*Pxz[6]+1.224744871391589*Pxz[3]+0.7071067811865475*Pxz[2]; 
  Pxz_xr[2] = 1.58113883008419*Pxz[8]+1.224744871391589*Pxz[7]+0.7071067811865475*Pxz[5]; 
 
  Pxy_yl[0] = 1.58113883008419*Pxy[5]-1.224744871391589*Pxy[2]+0.7071067811865475*Pxy[0]; 
  Pxy_yl[1] = 1.58113883008419*Pxy[7]-1.224744871391589*Pxy[3]+0.7071067811865475*Pxy[1]; 
  Pxy_yl[2] = 1.58113883008419*Pxy[8]-1.224744871391589*Pxy[6]+0.7071067811865475*Pxy[4]; 
  Pyy_yl[0] = 1.58113883008419*Pyy[5]-1.224744871391589*Pyy[2]+0.7071067811865475*Pyy[0]; 
  Pyy_yl[1] = 1.58113883008419*Pyy[7]-1.224744871391589*Pyy[3]+0.7071067811865475*Pyy[1]; 
  Pyy_yl[2] = 1.58113883008419*Pyy[8]-1.224744871391589*Pyy[6]+0.7071067811865475*Pyy[4]; 
  Pyz_yl[0] = 1.58113883008419*Pyz[5]-1.224744871391589*Pyz[2]+0.7071067811865475*Pyz[0]; 
  Pyz_yl[1] = 1.58113883008419*Pyz[7]-1.224744871391589*Pyz[3]+0.7071067811865475*Pyz[1]; 
  Pyz_yl[2] = 1.58113883008419*Pyz[8]-1.224744871391589*Pyz[6]+0.7071067811865475*Pyz[4]; 
 
  Pxy_yr[0] = 1.58113883008419*Pxy[5]+1.224744871391589*Pxy[2]+0.7071067811865475*Pxy[0]; 
  Pxy_yr[1] = 1.58113883008419*Pxy[7]+1.224744871391589*Pxy[3]+0.7071067811865475*Pxy[1]; 
  Pxy_yr[2] = 1.58113883008419*Pxy[8]+1.224744871391589*Pxy[6]+0.7071067811865475*Pxy[4]; 
  Pyy_yr[0] = 1.58113883008419*Pyy[5]+1.224744871391589*Pyy[2]+0.7071067811865475*Pyy[0]; 
  Pyy_yr[1] = 1.58113883008419*Pyy[7]+1.224744871391589*Pyy[3]+0.7071067811865475*Pyy[1]; 
  Pyy_yr[2] = 1.58113883008419*Pyy[8]+1.224744871391589*Pyy[6]+0.7071067811865475*Pyy[4]; 
  Pyz_yr[0] = 1.58113883008419*Pyz[5]+1.224744871391589*Pyz[2]+0.7071067811865475*Pyz[0]; 
  Pyz_yr[1] = 1.58113883008419*Pyz[7]+1.224744871391589*Pyz[3]+0.7071067811865475*Pyz[1]; 
  Pyz_yr[2] = 1.58113883008419*Pyz[8]+1.224744871391589*Pyz[6]+0.7071067811865475*Pyz[4]; 
 
} 
