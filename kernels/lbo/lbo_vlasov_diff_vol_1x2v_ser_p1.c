#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 

  return fabs(0.9428090415820636*nuVtSqSum[0]*rdvxSq4)+fabs(0.9428090415820636*nuVtSqSum[0]*rdvySq4); 

} 
