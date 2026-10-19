#include <math.h>
#include <gkyl_flash_compute_cfl_kernels.h>

double flash_compute_cfl_vlasov_pkpm_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *pkpm_max_b, const double *pkpm_lax_l, const double *pkpm_lax_r, const double *pkpm_accel_vars) 
{
  // w[NDIM]:         Cell-center coordinates. 
  // dxv[NDIM]:       Cell spacing. 
  // pkpm_accel_vars: Input pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 

  const double *pkpm_max_b_l = &pkpm_max_b[0]; 
  const double *pkpm_max_b_r = &pkpm_max_b[1]; 

  const double *pkpm_lax_dir_l = &pkpm_lax_l[0]; 
  const double *pkpm_lax_dir_r = &pkpm_lax_r[0]; 

  double max_b_l = pkpm_max_b_l[0]; 
  double max_b_r = pkpm_max_b_r[0]; 

  double max_speed_l = pkpm_lax_dir_l[0]; 
  double max_speed_r = pkpm_lax_dir_r[0]; 

  // Compute CFL frequency from streaming term 
  double max_v_par = fmax(fabs(wvpar + dvpar/2), fabs(wvpar - dvpar/2)); 
  double cfl_l = max_v_par*max_b_l + max_speed_l; 
  double cfl_r = max_v_par*max_b_r + max_speed_r; 
  return 1.5*dx1*fmax(cfl_l, cfl_r);
} 
