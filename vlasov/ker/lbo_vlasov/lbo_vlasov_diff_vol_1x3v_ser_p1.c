#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuVtSqSum = &nuPrimMomsSum[6];

  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvzSq4 = 4.0/(dxv[3]*dxv[3]); 

  return fabs(6.363961030678928*nuVtSqSum[0]*rdvxSq4)+fabs(6.363961030678928*nuVtSqSum[0]*rdvySq4)+fabs(6.363961030678928*nuVtSqSum[0]*rdvzSq4); 

} 
