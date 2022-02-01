#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_vol_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]:      Cell-center coordinates. 
  // dxv[4]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvzSq4 = 4.0/(dxv[3]*dxv[3]); 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  return fabs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4)+fabs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvySq4)+fabs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvzSq4); 

} 
