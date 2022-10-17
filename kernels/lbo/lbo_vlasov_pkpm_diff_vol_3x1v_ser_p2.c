#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_diff_vol_3x1v_ser_p2(const double *w, const double *dxv, const double *nuVtSq, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]:      Cell-center coordinates. 
  // dxv[4]:    Cell spacing. 
  // nuVtSq:    Thermal speed squared times collisionality. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvparSq4 = 4.0/(dxv[3]*dxv[3]); 

  return fabs((3.181980515339463*nuVtSq[0]-3.557562367689425*(nuVtSq[9]+nuVtSq[8]+nuVtSq[7]))*rdvparSq4); 

} 
