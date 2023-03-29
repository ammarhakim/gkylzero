#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_diff_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *nuVtSq, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuVtSq:    Thermal speed squared times collisionality. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvparSq4 = 4.0/(dxv[2]*dxv[2]); 

  return fabs(4.5*nuVtSq[0]*rdvparSq4); 

} 
