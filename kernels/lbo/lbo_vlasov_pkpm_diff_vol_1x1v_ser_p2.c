#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_diff_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *nuVtSq, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]:      Cell-center coordinates. 
  // dxv[2]:    Cell spacing. 
  // nuVtSq:    Thermal speed squared times collisionality. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 

  return fabs((1.272792206135785*nuVtSq[0]-1.42302494707577*nuVtSq[2])*rdvxSq4); 

} 
