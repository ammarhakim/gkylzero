#include <math.h>
#include <gkyl_flash_compute_cfl_kernels.h>

double flash_compute_cfl_vlasov_pkpm_vol_1x1v_ser_p1(const double *w, const double *dxv, const double* pkpm_accel_vars) 
{
  // w[NDIM]:         Cell-center coordinates. 
  // dxv[NDIM]:       Cell spacing. 
  // pkpm_accel_vars: Input pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 

  const double dv1par = 2.0/dxv[1]; 

  const double *p_perp_div_b = &pkpm_accel_vars[0]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 5.0*fabs(0.35355339059327373*p_perp_div_b[0])*dv1par; 

  return cflFreq_mid; 
} 
