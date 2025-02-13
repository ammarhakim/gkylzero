#include <gkyl_mom_fpo_vlasov_kernels.h> 
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_accum_2x3v_ser_p2(const double *drag_diff_coeff_corrs, double *drag_coeff, double *drag_coeff_surf, double *diff_coeff, double *diff_coeff_surf) 
{
  // drag_diff_coeff_corrs: Corrections to be added to coeffs, function of config space only. 
  // drag_coeff: FPO drag coefficient. 
  // diff_coeff: FPO diffusion coefficient. 
 
  // Index into drag and diffusion coefficients. 
  double *ax = &drag_coeff[0]; 
  double *ay = &drag_coeff[112]; 
  double *az = &drag_coeff[224]; 
  double *Dxx = &diff_coeff[0]; 
  double *Dyy = &diff_coeff[448]; 
  double *Dzz = &diff_coeff[896]; 
 
  // Index into surface expansions. 
  double *ax_surf = &drag_coeff_surf[0]; 
  double *ay_surf = &drag_coeff_surf[48]; 
  double *az_surf = &drag_coeff_surf[96]; 
  double *Dxx_surf = &diff_coeff_surf[0]; 
  double *Dyy_surf = &diff_coeff_surf[384]; 
  double *Dzz_surf = &diff_coeff_surf[768]; 
 
  // Index into correction array. 
  const double* ax_corr = &drag_diff_coeff_corrs[0]; 
  const double* ay_corr = &drag_diff_coeff_corrs[8]; 
  const double* az_corr = &drag_diff_coeff_corrs[16]; 
  const double* D_corr = &drag_diff_coeff_corrs[24]; 
 
  ax[0] += 2.8284271247461907*ax_corr[0]; 
  ax[1] += 2.8284271247461907*ax_corr[1]; 
  ax[2] += 2.8284271247461907*ax_corr[2]; 
  ax[6] += 2.8284271247461907*ax_corr[3]; 
  ax[16] += 2.8284271247461907*ax_corr[4]; 
  ax[17] += 2.8284271247461907*ax_corr[5]; 
  ax[31] += 2.8284271247461907*ax_corr[6]; 
  ax[32] += 2.8284271247461907*ax_corr[7]; 
 
  ay[0] += 2.8284271247461907*ay_corr[0]; 
  ay[1] += 2.8284271247461907*ay_corr[1]; 
  ay[2] += 2.8284271247461907*ay_corr[2]; 
  ay[6] += 2.8284271247461907*ay_corr[3]; 
  ay[16] += 2.8284271247461907*ay_corr[4]; 
  ay[17] += 2.8284271247461907*ay_corr[5]; 
  ay[31] += 2.8284271247461907*ay_corr[6]; 
  ay[32] += 2.8284271247461907*ay_corr[7]; 
 
  az[0] += 2.8284271247461907*az_corr[0]; 
  az[1] += 2.8284271247461907*az_corr[1]; 
  az[2] += 2.8284271247461907*az_corr[2]; 
  az[6] += 2.8284271247461907*az_corr[3]; 
  az[16] += 2.8284271247461907*az_corr[4]; 
  az[17] += 2.8284271247461907*az_corr[5]; 
  az[31] += 2.8284271247461907*az_corr[6]; 
  az[32] += 2.8284271247461907*az_corr[7]; 
 
  Dxx[0] += 2.8284271247461907*D_corr[0]; 
  Dxx[1] += 2.8284271247461907*D_corr[1]; 
  Dxx[2] += 2.8284271247461907*D_corr[2]; 
  Dxx[6] += 2.8284271247461907*D_corr[3]; 
  Dxx[16] += 2.8284271247461907*D_corr[4]; 
  Dxx[17] += 2.8284271247461907*D_corr[5]; 
  Dxx[31] += 2.8284271247461907*D_corr[6]; 
  Dxx[32] += 2.8284271247461907*D_corr[7]; 
 
  Dyy[0] += 2.8284271247461907*D_corr[0]; 
  Dyy[1] += 2.8284271247461907*D_corr[1]; 
  Dyy[2] += 2.8284271247461907*D_corr[2]; 
  Dyy[6] += 2.8284271247461907*D_corr[3]; 
  Dyy[16] += 2.8284271247461907*D_corr[4]; 
  Dyy[17] += 2.8284271247461907*D_corr[5]; 
  Dyy[31] += 2.8284271247461907*D_corr[6]; 
  Dyy[32] += 2.8284271247461907*D_corr[7]; 
 
  Dzz[0] += 2.8284271247461907*D_corr[0]; 
  Dzz[1] += 2.8284271247461907*D_corr[1]; 
  Dzz[2] += 2.8284271247461907*D_corr[2]; 
  Dzz[6] += 2.8284271247461907*D_corr[3]; 
  Dzz[16] += 2.8284271247461907*D_corr[4]; 
  Dzz[17] += 2.8284271247461907*D_corr[5]; 
  Dzz[31] += 2.8284271247461907*D_corr[6]; 
  Dzz[32] += 2.8284271247461907*D_corr[7]; 
 
  ax_surf[0] += 2.0*ax_corr[0]; 
  ax_surf[1] += 2.0*ax_corr[1]; 
  ax_surf[2] += 2.0*ax_corr[2]; 
  ax_surf[5] += 2.0*ax_corr[3]; 
  ax_surf[11] += 2.0*ax_corr[4]; 
  ax_surf[12] += 2.0*ax_corr[5]; 
  ax_surf[19] += 2.0*ax_corr[6]; 
  ax_surf[20] += 2.0*ax_corr[7]; 
 
  ay_surf[0] += 2.0*ay_corr[0]; 
  ay_surf[1] += 2.0*ay_corr[1]; 
  ay_surf[2] += 2.0*ay_corr[2]; 
  ay_surf[5] += 2.0*ay_corr[3]; 
  ay_surf[11] += 2.0*ay_corr[4]; 
  ay_surf[12] += 2.0*ay_corr[5]; 
  ay_surf[19] += 2.0*ay_corr[6]; 
  ay_surf[20] += 2.0*ay_corr[7]; 
 
  az_surf[0] += 2.0*az_corr[0]; 
  az_surf[1] += 2.0*az_corr[1]; 
  az_surf[2] += 2.0*az_corr[2]; 
  az_surf[5] += 2.0*az_corr[3]; 
  az_surf[11] += 2.0*az_corr[4]; 
  az_surf[12] += 2.0*az_corr[5]; 
  az_surf[19] += 2.0*az_corr[6]; 
  az_surf[20] += 2.0*az_corr[7]; 
 
  Dxx_surf[0] += 2.0*D_corr[0]; 
  Dxx_surf[1] += 2.0*D_corr[1]; 
  Dxx_surf[2] += 2.0*D_corr[2]; 
  Dxx_surf[5] += 2.0*D_corr[3]; 
  Dxx_surf[11] += 2.0*D_corr[4]; 
  Dxx_surf[12] += 2.0*D_corr[5]; 
  Dxx_surf[19] += 2.0*D_corr[6]; 
  Dxx_surf[20] += 2.0*D_corr[7]; 
 
  Dyy_surf[0] += 2.0*D_corr[0]; 
  Dyy_surf[1] += 2.0*D_corr[1]; 
  Dyy_surf[2] += 2.0*D_corr[2]; 
  Dyy_surf[5] += 2.0*D_corr[3]; 
  Dyy_surf[11] += 2.0*D_corr[4]; 
  Dyy_surf[12] += 2.0*D_corr[5]; 
  Dyy_surf[19] += 2.0*D_corr[6]; 
  Dyy_surf[20] += 2.0*D_corr[7]; 
 
  Dzz_surf[0] += 2.0*D_corr[0]; 
  Dzz_surf[1] += 2.0*D_corr[1]; 
  Dzz_surf[2] += 2.0*D_corr[2]; 
  Dzz_surf[5] += 2.0*D_corr[3]; 
  Dzz_surf[11] += 2.0*D_corr[4]; 
  Dzz_surf[12] += 2.0*D_corr[5]; 
  Dzz_surf[19] += 2.0*D_corr[6]; 
  Dzz_surf[20] += 2.0*D_corr[7]; 
} 

