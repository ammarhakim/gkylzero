#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_vol_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // w[3]:      cell-center coordinates. 
  // dxv[3]:    cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         input distribution function.
  // out:       incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[1]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[1]   = 2.0/dxv[2]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 

  return fabs(1.333333333333333*bmag_inv[1]*nuVtSqSum[1]*rdvSq4[1]*w[2]*m_+nuVtSqSum[0]*(1.333333333333333*bmag_inv[0]*rdvSq4[1]*w[2]*m_+0.9428090415820636*rdvSq4[0])); 

} 
