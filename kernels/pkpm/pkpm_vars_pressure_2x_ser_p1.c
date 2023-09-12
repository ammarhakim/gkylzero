#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_pressure_2x_ser_p1(const double *bvar, const double *bvar_surf, const double *vlasov_pkpm_moms, 
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
  const double *p_parallel = &vlasov_pkpm_moms[4]; 
  const double *p_perp = &vlasov_pkpm_moms[8]; 
  const double *bxbx = &bvar[12]; 
  const double *bxby = &bvar[16]; 
  const double *bxbz = &bvar[20]; 
  const double *byby = &bvar[24]; 
  const double *bybz = &bvar[28]; 
  const double *bzbz = &bvar[32]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[4]; 
  double *Pxz = &p_ij[8]; 
  double *Pyy = &p_ij[12]; 
  double *Pyz = &p_ij[16]; 
  double *Pzz = &p_ij[20]; 

  const double *bxbx_xl = &bvar_surf[4]; 
  const double *bxbx_xr = &bvar_surf[6]; 
  const double *bxby_xl = &bvar_surf[8]; 
  const double *bxby_xr = &bvar_surf[10]; 
  const double *bxbz_xl = &bvar_surf[12]; 
  const double *bxbz_xr = &bvar_surf[14]; 
  double *Pxx_xl = &p_ij_surf[0]; 
  double *Pxx_xr = &p_ij_surf[2]; 
  double *Pxy_xl = &p_ij_surf[4]; 
  double *Pxy_xr = &p_ij_surf[6]; 
  double *Pxz_xl = &p_ij_surf[8]; 
  double *Pxz_xr = &p_ij_surf[10]; 

  const double *byby_yl = &bvar_surf[20]; 
  const double *byby_yr = &bvar_surf[22]; 
  const double *bxby_yl = &bvar_surf[24]; 
  const double *bxby_yr = &bvar_surf[26]; 
  const double *bybz_yl = &bvar_surf[28]; 
  const double *bybz_yr = &bvar_surf[30]; 
  double *Pxy_yl = &p_ij_surf[12]; 
  double *Pxy_yr = &p_ij_surf[14]; 
  double *Pyy_yl = &p_ij_surf[16]; 
  double *Pyy_yr = &p_ij_surf[18]; 
  double *Pyz_yl = &p_ij_surf[20]; 
  double *Pyz_yr = &p_ij_surf[22]; 
 
  double DP[4] = {0.0}; 
  DP[0] = p_parallel[0] - p_perp[0]; 
  DP[1] = p_parallel[1] - p_perp[1]; 
  DP[2] = p_parallel[2] - p_perp[2]; 
  DP[3] = p_parallel[3] - p_perp[3]; 
  // DP b_i b_j. 
  double DP_bxbx[4] = {0.0}; 
  binop_mul_2d_ser_p1(DP, bxbx, DP_bxbx); 
 
  double DP_bxby[4] = {0.0}; 
  binop_mul_2d_ser_p1(DP, bxby, DP_bxby); 
 
  double DP_bxbz[4] = {0.0}; 
  binop_mul_2d_ser_p1(DP, bxbz, DP_bxbz); 
 
  double DP_byby[4] = {0.0}; 
  binop_mul_2d_ser_p1(DP, byby, DP_byby); 
 
  double DP_bybz[4] = {0.0}; 
  binop_mul_2d_ser_p1(DP, bybz, DP_bybz); 
 
  double DP_bzbz[4] = {0.0}; 
  binop_mul_2d_ser_p1(DP, bzbz, DP_bzbz); 
 
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
 
  double p_par_xl[2] = {0.0}; 
  double p_par_xr[2] = {0.0}; 
  double p_par_yl[2] = {0.0}; 
  double p_par_yr[2] = {0.0}; 
 
  double p_perp_xl[2] = {0.0}; 
  double p_perp_xr[2] = {0.0}; 
  double p_perp_yl[2] = {0.0}; 
  double p_perp_yr[2] = {0.0}; 
 
  p_par_xl[0] = 0.7071067811865475*p_parallel[0]-1.224744871391589*p_parallel[1]; 
  p_par_xl[1] = 0.7071067811865475*p_parallel[2]-1.224744871391589*p_parallel[3]; 
  p_par_xr[0] = 1.224744871391589*p_parallel[1]+0.7071067811865475*p_parallel[0]; 
  p_par_xr[1] = 1.224744871391589*p_parallel[3]+0.7071067811865475*p_parallel[2]; 
  p_perp_xl[0] = 0.7071067811865475*p_perp[0]-1.224744871391589*p_perp[1]; 
  p_perp_xl[1] = 0.7071067811865475*p_perp[2]-1.224744871391589*p_perp[3]; 
  p_perp_xr[0] = 1.224744871391589*p_perp[1]+0.7071067811865475*p_perp[0]; 
  p_perp_xr[1] = 1.224744871391589*p_perp[3]+0.7071067811865475*p_perp[2]; 
 
  p_par_yl[0] = 0.7071067811865475*p_parallel[0]-1.224744871391589*p_parallel[2]; 
  p_par_yl[1] = 0.7071067811865475*p_parallel[1]-1.224744871391589*p_parallel[3]; 
  p_par_yr[0] = 1.224744871391589*p_parallel[2]+0.7071067811865475*p_parallel[0]; 
  p_par_yr[1] = 1.224744871391589*p_parallel[3]+0.7071067811865475*p_parallel[1]; 
  p_perp_yl[0] = 0.7071067811865475*p_perp[0]-1.224744871391589*p_perp[2]; 
  p_perp_yl[1] = 0.7071067811865475*p_perp[1]-1.224744871391589*p_perp[3]; 
  p_perp_yr[0] = 1.224744871391589*p_perp[2]+0.7071067811865475*p_perp[0]; 
  p_perp_yr[1] = 1.224744871391589*p_perp[3]+0.7071067811865475*p_perp[1]; 
 
  double DP_xl[2] = {0.0}; 
  double DP_xr[2] = {0.0}; 
  double DP_yl[2] = {0.0}; 
  double DP_yr[2] = {0.0}; 
 
  DP_xl[0] = p_par_xl[0] - p_perp_xl[0]; 
  DP_xr[0] = p_par_xr[0] - p_perp_xr[0]; 
  DP_yl[0] = p_par_yl[0] - p_perp_yl[0]; 
  DP_yr[0] = p_par_yr[0] - p_perp_yr[0]; 
  DP_xl[1] = p_par_xl[1] - p_perp_xl[1]; 
  DP_xr[1] = p_par_xr[1] - p_perp_xr[1]; 
  DP_yl[1] = p_par_yl[1] - p_perp_yl[1]; 
  DP_yr[1] = p_par_yr[1] - p_perp_yr[1]; 
  // DP b_i b_j at lower x surface. 
  double DP_bxbx_xl[2] = {0.0}; 
  double DP_bxby_xl[2] = {0.0}; 
  double DP_bxbz_xl[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP_xl, bxbx_xl, DP_bxbx_xl); 
  binop_mul_1d_ser_p1(DP_xl, bxby_xl, DP_bxby_xl); 
  binop_mul_1d_ser_p1(DP_xl, bxbz_xl, DP_bxbz_xl); 
 
  // DP b_i b_j at upper x surface. 
  double DP_bxbx_xr[2] = {0.0}; 
  double DP_bxby_xr[2] = {0.0}; 
  double DP_bxbz_xr[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP_xr, bxbx_xr, DP_bxbx_xr); 
  binop_mul_1d_ser_p1(DP_xr, bxby_xr, DP_bxby_xr); 
  binop_mul_1d_ser_p1(DP_xr, bxbz_xr, DP_bxbz_xr); 
 
  // DP b_i b_j at lower y surface. 
  double DP_bxby_yl[2] = {0.0}; 
  double DP_byby_yl[2] = {0.0}; 
  double DP_bybz_yl[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP_yl, bxby_yl, DP_bxby_yl); 
  binop_mul_1d_ser_p1(DP_yl, byby_yl, DP_byby_yl); 
  binop_mul_1d_ser_p1(DP_yl, bybz_yl, DP_bybz_yl); 
 
  // DP b_i b_j at upper y surface. 
  double DP_bxby_yr[2] = {0.0}; 
  double DP_byby_yr[2] = {0.0}; 
  double DP_bybz_yr[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP_yr, bxby_yr, DP_bxby_yr); 
  binop_mul_1d_ser_p1(DP_yr, byby_yr, DP_byby_yr); 
  binop_mul_1d_ser_p1(DP_yr, bybz_yr, DP_bybz_yr); 
 
  Pxx_xl[0] = DP_bxbx_xl[0] + p_perp_xl[0]; 
  Pxx_xr[0] = DP_bxbx_xr[0] + p_perp_xr[0]; 
  Pxy_xl[0] = DP_bxby_xl[0]; 
  Pxy_xr[0] = DP_bxby_xr[0]; 
  Pxz_xl[0] = DP_bxbz_xl[0]; 
  Pxz_xr[0] = DP_bxbz_xr[0]; 
 
  Pxy_yl[0] = DP_bxby_yl[0]; 
  Pxy_yr[0] = DP_bxby_yr[0]; 
  Pyy_yl[0] = DP_byby_yl[0] + p_perp_yl[0]; 
  Pyy_yr[0] = DP_byby_yr[0] + p_perp_yr[0]; 
  Pyz_yl[0] = DP_bybz_yl[0]; 
  Pyz_yr[0] = DP_bybz_yr[0]; 
 
  Pxx_xl[1] = DP_bxbx_xl[1] + p_perp_xl[1]; 
  Pxx_xr[1] = DP_bxbx_xr[1] + p_perp_xr[1]; 
  Pxy_xl[1] = DP_bxby_xl[1]; 
  Pxy_xr[1] = DP_bxby_xr[1]; 
  Pxz_xl[1] = DP_bxbz_xl[1]; 
  Pxz_xr[1] = DP_bxbz_xr[1]; 
 
  Pxy_yl[1] = DP_bxby_yl[1]; 
  Pxy_yr[1] = DP_bxby_yr[1]; 
  Pyy_yl[1] = DP_byby_yl[1] + p_perp_yl[1]; 
  Pyy_yr[1] = DP_byby_yr[1] + p_perp_yr[1]; 
  Pyz_yl[1] = DP_bybz_yl[1]; 
  Pyz_yr[1] = DP_bybz_yr[1]; 
 
} 
