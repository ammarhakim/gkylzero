#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_pressure_3x_tensor_p2(const double *bvar, const double *vlasov_pkpm_moms, 
    double* GKYL_RESTRICT p_ij) 
{ 
  // bvar:             Input volume expansion of magnetic field unit vector and tensor.
  //                   [bx, by, bz, bxbx, bxby, bxbz, byby, bybz, bzbz] 
  // vlasov_pkpm_moms: [rho, p_parallel, p_perp], Moments computed from kinetic equation in pkpm model.
  // p_ij:             Output volume expansion of p_ij = (p_par - p_perp) b_i b_j + p_perp g_ij.

  // Parallel/Perp pressure are first/second component of pkpm moment array and unit tensor are last six components of bvar array.
  const double *p_parallel = &vlasov_pkpm_moms[27]; 
  const double *p_perp = &vlasov_pkpm_moms[54]; 
  const double *bxbx = &bvar[81]; 
  const double *bxby = &bvar[108]; 
  const double *bxbz = &bvar[135]; 
  const double *byby = &bvar[162]; 
  const double *bybz = &bvar[189]; 
  const double *bzbz = &bvar[216]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[27]; 
  double *Pxz = &p_ij[54]; 
  double *Pyy = &p_ij[81]; 
  double *Pyz = &p_ij[108]; 
  double *Pzz = &p_ij[135]; 

  double DP[27] = {0.0}; 
  DP[0] = p_parallel[0] - p_perp[0]; 
  DP[1] = p_parallel[1] - p_perp[1]; 
  DP[2] = p_parallel[2] - p_perp[2]; 
  DP[3] = p_parallel[3] - p_perp[3]; 
  DP[4] = p_parallel[4] - p_perp[4]; 
  DP[5] = p_parallel[5] - p_perp[5]; 
  DP[6] = p_parallel[6] - p_perp[6]; 
  DP[7] = p_parallel[7] - p_perp[7]; 
  DP[8] = p_parallel[8] - p_perp[8]; 
  DP[9] = p_parallel[9] - p_perp[9]; 
  DP[10] = p_parallel[10] - p_perp[10]; 
  DP[11] = p_parallel[11] - p_perp[11]; 
  DP[12] = p_parallel[12] - p_perp[12]; 
  DP[13] = p_parallel[13] - p_perp[13]; 
  DP[14] = p_parallel[14] - p_perp[14]; 
  DP[15] = p_parallel[15] - p_perp[15]; 
  DP[16] = p_parallel[16] - p_perp[16]; 
  DP[17] = p_parallel[17] - p_perp[17]; 
  DP[18] = p_parallel[18] - p_perp[18]; 
  DP[19] = p_parallel[19] - p_perp[19]; 
  DP[20] = p_parallel[20] - p_perp[20]; 
  DP[21] = p_parallel[21] - p_perp[21]; 
  DP[22] = p_parallel[22] - p_perp[22]; 
  DP[23] = p_parallel[23] - p_perp[23]; 
  DP[24] = p_parallel[24] - p_perp[24]; 
  DP[25] = p_parallel[25] - p_perp[25]; 
  DP[26] = p_parallel[26] - p_perp[26]; 
  // DP b_i b_j. 
  double DP_bxbx[27] = {0.0}; 
  binop_mul_3d_tensor_p2(DP, bxbx, DP_bxbx); 
 
  double DP_bxby[27] = {0.0}; 
  binop_mul_3d_tensor_p2(DP, bxby, DP_bxby); 
 
  double DP_bxbz[27] = {0.0}; 
  binop_mul_3d_tensor_p2(DP, bxbz, DP_bxbz); 
 
  double DP_byby[27] = {0.0}; 
  binop_mul_3d_tensor_p2(DP, byby, DP_byby); 
 
  double DP_bybz[27] = {0.0}; 
  binop_mul_3d_tensor_p2(DP, bybz, DP_bybz); 
 
  double DP_bzbz[27] = {0.0}; 
  binop_mul_3d_tensor_p2(DP, bzbz, DP_bzbz); 
 
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
  Pxx[9] = DP_bxbx[9] + p_perp[9]; 
  Pxy[9] = DP_bxby[9]; 
  Pxz[9] = DP_bxbz[9]; 
  Pyy[9] = DP_byby[9] + p_perp[9]; 
  Pyz[9] = DP_bybz[9]; 
  Pzz[9] = DP_bzbz[9] + p_perp[9]; 
  Pxx[10] = DP_bxbx[10] + p_perp[10]; 
  Pxy[10] = DP_bxby[10]; 
  Pxz[10] = DP_bxbz[10]; 
  Pyy[10] = DP_byby[10] + p_perp[10]; 
  Pyz[10] = DP_bybz[10]; 
  Pzz[10] = DP_bzbz[10] + p_perp[10]; 
  Pxx[11] = DP_bxbx[11] + p_perp[11]; 
  Pxy[11] = DP_bxby[11]; 
  Pxz[11] = DP_bxbz[11]; 
  Pyy[11] = DP_byby[11] + p_perp[11]; 
  Pyz[11] = DP_bybz[11]; 
  Pzz[11] = DP_bzbz[11] + p_perp[11]; 
  Pxx[12] = DP_bxbx[12] + p_perp[12]; 
  Pxy[12] = DP_bxby[12]; 
  Pxz[12] = DP_bxbz[12]; 
  Pyy[12] = DP_byby[12] + p_perp[12]; 
  Pyz[12] = DP_bybz[12]; 
  Pzz[12] = DP_bzbz[12] + p_perp[12]; 
  Pxx[13] = DP_bxbx[13] + p_perp[13]; 
  Pxy[13] = DP_bxby[13]; 
  Pxz[13] = DP_bxbz[13]; 
  Pyy[13] = DP_byby[13] + p_perp[13]; 
  Pyz[13] = DP_bybz[13]; 
  Pzz[13] = DP_bzbz[13] + p_perp[13]; 
  Pxx[14] = DP_bxbx[14] + p_perp[14]; 
  Pxy[14] = DP_bxby[14]; 
  Pxz[14] = DP_bxbz[14]; 
  Pyy[14] = DP_byby[14] + p_perp[14]; 
  Pyz[14] = DP_bybz[14]; 
  Pzz[14] = DP_bzbz[14] + p_perp[14]; 
  Pxx[15] = DP_bxbx[15] + p_perp[15]; 
  Pxy[15] = DP_bxby[15]; 
  Pxz[15] = DP_bxbz[15]; 
  Pyy[15] = DP_byby[15] + p_perp[15]; 
  Pyz[15] = DP_bybz[15]; 
  Pzz[15] = DP_bzbz[15] + p_perp[15]; 
  Pxx[16] = DP_bxbx[16] + p_perp[16]; 
  Pxy[16] = DP_bxby[16]; 
  Pxz[16] = DP_bxbz[16]; 
  Pyy[16] = DP_byby[16] + p_perp[16]; 
  Pyz[16] = DP_bybz[16]; 
  Pzz[16] = DP_bzbz[16] + p_perp[16]; 
  Pxx[17] = DP_bxbx[17] + p_perp[17]; 
  Pxy[17] = DP_bxby[17]; 
  Pxz[17] = DP_bxbz[17]; 
  Pyy[17] = DP_byby[17] + p_perp[17]; 
  Pyz[17] = DP_bybz[17]; 
  Pzz[17] = DP_bzbz[17] + p_perp[17]; 
  Pxx[18] = DP_bxbx[18] + p_perp[18]; 
  Pxy[18] = DP_bxby[18]; 
  Pxz[18] = DP_bxbz[18]; 
  Pyy[18] = DP_byby[18] + p_perp[18]; 
  Pyz[18] = DP_bybz[18]; 
  Pzz[18] = DP_bzbz[18] + p_perp[18]; 
  Pxx[19] = DP_bxbx[19] + p_perp[19]; 
  Pxy[19] = DP_bxby[19]; 
  Pxz[19] = DP_bxbz[19]; 
  Pyy[19] = DP_byby[19] + p_perp[19]; 
  Pyz[19] = DP_bybz[19]; 
  Pzz[19] = DP_bzbz[19] + p_perp[19]; 
  Pxx[20] = DP_bxbx[20] + p_perp[20]; 
  Pxy[20] = DP_bxby[20]; 
  Pxz[20] = DP_bxbz[20]; 
  Pyy[20] = DP_byby[20] + p_perp[20]; 
  Pyz[20] = DP_bybz[20]; 
  Pzz[20] = DP_bzbz[20] + p_perp[20]; 
  Pxx[21] = DP_bxbx[21] + p_perp[21]; 
  Pxy[21] = DP_bxby[21]; 
  Pxz[21] = DP_bxbz[21]; 
  Pyy[21] = DP_byby[21] + p_perp[21]; 
  Pyz[21] = DP_bybz[21]; 
  Pzz[21] = DP_bzbz[21] + p_perp[21]; 
  Pxx[22] = DP_bxbx[22] + p_perp[22]; 
  Pxy[22] = DP_bxby[22]; 
  Pxz[22] = DP_bxbz[22]; 
  Pyy[22] = DP_byby[22] + p_perp[22]; 
  Pyz[22] = DP_bybz[22]; 
  Pzz[22] = DP_bzbz[22] + p_perp[22]; 
  Pxx[23] = DP_bxbx[23] + p_perp[23]; 
  Pxy[23] = DP_bxby[23]; 
  Pxz[23] = DP_bxbz[23]; 
  Pyy[23] = DP_byby[23] + p_perp[23]; 
  Pyz[23] = DP_bybz[23]; 
  Pzz[23] = DP_bzbz[23] + p_perp[23]; 
  Pxx[24] = DP_bxbx[24] + p_perp[24]; 
  Pxy[24] = DP_bxby[24]; 
  Pxz[24] = DP_bxbz[24]; 
  Pyy[24] = DP_byby[24] + p_perp[24]; 
  Pyz[24] = DP_bybz[24]; 
  Pzz[24] = DP_bzbz[24] + p_perp[24]; 
  Pxx[25] = DP_bxbx[25] + p_perp[25]; 
  Pxy[25] = DP_bxby[25]; 
  Pxz[25] = DP_bxbz[25]; 
  Pyy[25] = DP_byby[25] + p_perp[25]; 
  Pyz[25] = DP_bybz[25]; 
  Pzz[25] = DP_bzbz[25] + p_perp[25]; 
  Pxx[26] = DP_bxbx[26] + p_perp[26]; 
  Pxy[26] = DP_bxby[26]; 
  Pxz[26] = DP_bxbz[26]; 
  Pyy[26] = DP_byby[26] + p_perp[26]; 
  Pyz[26] = DP_bybz[26]; 
  Pzz[26] = DP_bzbz[26] + p_perp[26]; 
 
} 