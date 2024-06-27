#include <gkyl_mom_fpo_vlasov_kernels.h> 
GKYL_CU_DH void mom_fpo_vlasov_coeff_correct_accum_3x3v_ser_p1(const double *drag_diff_coeff_corrs, double *drag_coeff, double *drag_coeff_surf, double *diff_coeff, double *diff_coeff_surf) 
{
  // drag_diff_coeff_corrs: Corrections to be added to coeffs, function of config space only. 
  // drag_coeff: FPO drag coefficient. 
  // diff_coeff: FPO diffusion coefficient. 
 
  // Index into drag and diffusion coefficients. 
  double *ax = &drag_coeff[0]; 
  double *ay = &drag_coeff[160]; 
  double *az = &drag_coeff[320]; 
  double *Dxx = &diff_coeff[0]; 
  double *Dyy = &diff_coeff[640]; 
  double *Dzz = &diff_coeff[1280]; 
 
  // Index into surface expansions. 
  double *ax_surf = &drag_coeff_surf[0]; 
  double *ay_surf = &drag_coeff_surf[32]; 
  double *az_surf = &drag_coeff_surf[64]; 
  double *Dxx_surf = &diff_coeff_surf[0]; 
  double *Dyy_surf = &diff_coeff_surf[256]; 
  double *Dzz_surf = &diff_coeff_surf[512]; 
 
  // Index into correction array. 
  const double* ax_corr = &drag_diff_coeff_corrs[0]; 
  const double* ay_corr = &drag_diff_coeff_corrs[8]; 
  const double* az_corr = &drag_diff_coeff_corrs[16]; 
  const double* D_corr = &drag_diff_coeff_corrs[24]; 
 
  ax[0] += 2.8284271247461907*ax_corr[0]; 
  ax[1] += 2.8284271247461907*ax_corr[1]; 
  ax[2] += 2.8284271247461907*ax_corr[2]; 
  ax[3] += 2.8284271247461907*ax_corr[3]; 
  ax[7] += 2.8284271247461907*ax_corr[4]; 
  ax[8] += 2.8284271247461907*ax_corr[5]; 
  ax[9] += 2.8284271247461907*ax_corr[6]; 
  ax[22] += 2.8284271247461907*ax_corr[7]; 
 
  ay[0] += 2.8284271247461907*ay_corr[0]; 
  ay[1] += 2.8284271247461907*ay_corr[1]; 
  ay[2] += 2.8284271247461907*ay_corr[2]; 
  ay[3] += 2.8284271247461907*ay_corr[3]; 
  ay[7] += 2.8284271247461907*ay_corr[4]; 
  ay[8] += 2.8284271247461907*ay_corr[5]; 
  ay[9] += 2.8284271247461907*ay_corr[6]; 
  ay[22] += 2.8284271247461907*ay_corr[7]; 
 
  az[0] += 2.8284271247461907*az_corr[0]; 
  az[1] += 2.8284271247461907*az_corr[1]; 
  az[2] += 2.8284271247461907*az_corr[2]; 
  az[3] += 2.8284271247461907*az_corr[3]; 
  az[7] += 2.8284271247461907*az_corr[4]; 
  az[8] += 2.8284271247461907*az_corr[5]; 
  az[9] += 2.8284271247461907*az_corr[6]; 
  az[22] += 2.8284271247461907*az_corr[7]; 
 
  Dxx[0] += 2.8284271247461907*D_corr[0]; 
  Dxx[1] += 2.8284271247461907*D_corr[1]; 
  Dxx[2] += 2.8284271247461907*D_corr[2]; 
  Dxx[3] += 2.8284271247461907*D_corr[3]; 
  Dxx[7] += 2.8284271247461907*D_corr[4]; 
  Dxx[8] += 2.8284271247461907*D_corr[5]; 
  Dxx[9] += 2.8284271247461907*D_corr[6]; 
  Dxx[22] += 2.8284271247461907*D_corr[7]; 
 
  Dyy[0] += 2.8284271247461907*D_corr[0]; 
  Dyy[1] += 2.8284271247461907*D_corr[1]; 
  Dyy[2] += 2.8284271247461907*D_corr[2]; 
  Dyy[3] += 2.8284271247461907*D_corr[3]; 
  Dyy[7] += 2.8284271247461907*D_corr[4]; 
  Dyy[8] += 2.8284271247461907*D_corr[5]; 
  Dyy[9] += 2.8284271247461907*D_corr[6]; 
  Dyy[22] += 2.8284271247461907*D_corr[7]; 
 
  Dzz[0] += 2.8284271247461907*D_corr[0]; 
  Dzz[1] += 2.8284271247461907*D_corr[1]; 
  Dzz[2] += 2.8284271247461907*D_corr[2]; 
  Dzz[3] += 2.8284271247461907*D_corr[3]; 
  Dzz[7] += 2.8284271247461907*D_corr[4]; 
  Dzz[8] += 2.8284271247461907*D_corr[5]; 
  Dzz[9] += 2.8284271247461907*D_corr[6]; 
  Dzz[22] += 2.8284271247461907*D_corr[7]; 
 
  ax_surf[0] += 2.0*ax_corr[0]; 
  ax_surf[1] += 2.0*ax_corr[1]; 
  ax_surf[2] += 2.0*ax_corr[2]; 
  ax_surf[3] += 2.0*ax_corr[3]; 
  ax_surf[6] += 2.0*ax_corr[4]; 
  ax_surf[7] += 2.0*ax_corr[5]; 
  ax_surf[8] += 2.0*ax_corr[6]; 
  ax_surf[16] += 2.0*ax_corr[7]; 
 
  ay_surf[0] += 2.0*ay_corr[0]; 
  ay_surf[1] += 2.0*ay_corr[1]; 
  ay_surf[2] += 2.0*ay_corr[2]; 
  ay_surf[3] += 2.0*ay_corr[3]; 
  ay_surf[6] += 2.0*ay_corr[4]; 
  ay_surf[7] += 2.0*ay_corr[5]; 
  ay_surf[8] += 2.0*ay_corr[6]; 
  ay_surf[16] += 2.0*ay_corr[7]; 
 
  az_surf[0] += 2.0*az_corr[0]; 
  az_surf[1] += 2.0*az_corr[1]; 
  az_surf[2] += 2.0*az_corr[2]; 
  az_surf[3] += 2.0*az_corr[3]; 
  az_surf[6] += 2.0*az_corr[4]; 
  az_surf[7] += 2.0*az_corr[5]; 
  az_surf[8] += 2.0*az_corr[6]; 
  az_surf[16] += 2.0*az_corr[7]; 
 
  Dxx_surf[0] += 2.0*D_corr[0]; 
  Dxx_surf[1] += 2.0*D_corr[1]; 
  Dxx_surf[2] += 2.0*D_corr[2]; 
  Dxx_surf[3] += 2.0*D_corr[3]; 
  Dxx_surf[6] += 2.0*D_corr[4]; 
  Dxx_surf[7] += 2.0*D_corr[5]; 
  Dxx_surf[8] += 2.0*D_corr[6]; 
  Dxx_surf[16] += 2.0*D_corr[7]; 
 
  Dyy_surf[0] += 2.0*D_corr[0]; 
  Dyy_surf[1] += 2.0*D_corr[1]; 
  Dyy_surf[2] += 2.0*D_corr[2]; 
  Dyy_surf[3] += 2.0*D_corr[3]; 
  Dyy_surf[6] += 2.0*D_corr[4]; 
  Dyy_surf[7] += 2.0*D_corr[5]; 
  Dyy_surf[8] += 2.0*D_corr[6]; 
  Dyy_surf[16] += 2.0*D_corr[7]; 
 
  Dzz_surf[0] += 2.0*D_corr[0]; 
  Dzz_surf[1] += 2.0*D_corr[1]; 
  Dzz_surf[2] += 2.0*D_corr[2]; 
  Dzz_surf[3] += 2.0*D_corr[3]; 
  Dzz_surf[6] += 2.0*D_corr[4]; 
  Dzz_surf[7] += 2.0*D_corr[5]; 
  Dzz_surf[8] += 2.0*D_corr[6]; 
  Dzz_surf[16] += 2.0*D_corr[7]; 
} 

