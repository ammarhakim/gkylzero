#include <gkyl_euler_pkpm_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p1(const double *bvar, const double *bvar_surf, const double *vlasov_pkpm_moms, 
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
  const double *p_parallel = &vlasov_pkpm_moms[2]; 
  const double *p_perp = &vlasov_pkpm_moms[4]; 
  const double *bxbx = &bvar[6]; 
  const double *bxby = &bvar[8]; 
  const double *bxbz = &bvar[10]; 
  const double *byby = &bvar[12]; 
  const double *bybz = &bvar[14]; 
  const double *bzbz = &bvar[16]; 

  double *Pxx = &p_ij[0]; 
  double *Pxy = &p_ij[2]; 
  double *Pxz = &p_ij[4]; 
  double *Pyy = &p_ij[6]; 
  double *Pyz = &p_ij[8]; 
  double *Pzz = &p_ij[10]; 

  const double *bxbx_xl = &bvar_surf[2]; 
  const double *bxbx_xr = &bvar_surf[3]; 
  const double *bxby_xl = &bvar_surf[4]; 
  const double *bxby_xr = &bvar_surf[5]; 
  const double *bxbz_xl = &bvar_surf[6]; 
  const double *bxbz_xr = &bvar_surf[7]; 
  double *Pxx_xl = &p_ij_surf[0]; 
  double *Pxx_xr = &p_ij_surf[1]; 
  double *Pxy_xl = &p_ij_surf[2]; 
  double *Pxy_xr = &p_ij_surf[3]; 
  double *Pxz_xl = &p_ij_surf[4]; 
  double *Pxz_xr = &p_ij_surf[5]; 

  double DP[2] = {0.0}; 
  DP[0] = p_parallel[0] - p_perp[0]; 
  DP[1] = p_parallel[1] - p_perp[1]; 
  // DP b_i b_j. 
  double DP_bxbx[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP, bxbx, DP_bxbx); 
 
  double DP_bxby[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP, bxby, DP_bxby); 
 
  double DP_bxbz[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP, bxbz, DP_bxbz); 
 
  double DP_byby[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP, byby, DP_byby); 
 
  double DP_bybz[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP, bybz, DP_bybz); 
 
  double DP_bzbz[2] = {0.0}; 
  binop_mul_1d_ser_p1(DP, bzbz, DP_bzbz); 
 
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
 
  double p_par_xl = 0.7071067811865475*p_parallel[0]-1.224744871391589*p_parallel[1]; 
  double p_par_xr = 1.224744871391589*p_parallel[1]+0.7071067811865475*p_parallel[0]; 
  double p_perp_xl = 0.7071067811865475*p_perp[0]-1.224744871391589*p_perp[1]; 
  double p_perp_xr = 1.224744871391589*p_perp[1]+0.7071067811865475*p_perp[0]; 
  Pxx_xl[0] = (p_par_xl - p_perp_xl)*bxbx_xl[0] + p_perp_xl; 
  Pxx_xr[0] = (p_par_xr - p_perp_xr)*bxbx_xr[0] + p_perp_xr; 
  Pxy_xl[0] = (p_par_xl - p_perp_xl)*bxby_xl[0]; 
  Pxy_xr[0] = (p_par_xr - p_perp_xr)*bxby_xr[0]; 
  Pxz_xl[0] = (p_par_xl - p_perp_xl)*bxbz_xl[0]; 
  Pxz_xr[0] = (p_par_xr - p_perp_xr)*bxbz_xr[0]; 
 
} 
