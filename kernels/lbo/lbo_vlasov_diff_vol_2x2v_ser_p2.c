#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_vol_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]:      Cell-center coordinates. 
  // dxv[4]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvxSq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvySq4 = 4.0/(dxv[3]*dxv[3]); 

  double facDiff[8]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 
  facDiff[3] = nuVtSqSum[3]; 
  facDiff[4] = nuVtSqSum[4]; 
  facDiff[5] = nuVtSqSum[5]; 
  facDiff[6] = nuVtSqSum[6]; 
  facDiff[7] = nuVtSqSum[7]; 

  return fabs((0.9*facDiff[0]-1.006230589874905*(facDiff[5]+facDiff[4]))*rdvxSq4)+fabs((0.9*facDiff[0]-1.006230589874905*(facDiff[5]+facDiff[4]))*rdvySq4); 

} 
