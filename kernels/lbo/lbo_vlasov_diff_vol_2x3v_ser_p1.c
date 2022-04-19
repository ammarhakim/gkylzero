#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[5]:      Cell-center coordinates. 
  // dxv[5]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvxSq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvySq4 = 4.0/(dxv[3]*dxv[3]); 
  const double rdvzSq4 = 4.0/(dxv[4]*dxv[4]); 

  return fabs(0.6666666666666666*nuVtSqSum[0]*rdvxSq4)+fabs(0.6666666666666666*nuVtSqSum[0]*rdvySq4)+fabs(0.6666666666666666*nuVtSqSum[0]*rdvzSq4); 

} 
