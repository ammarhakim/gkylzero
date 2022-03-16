#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_vol_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // w[4]:      cell-center coordinates. 
  // dxv[4]:    cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         input distribution function.
  // out:       incremented output 
  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[2]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[1]   = 2.0/dxv[3]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 

  double facDiffVpar[8];
  // Expand nuVtSqSum in phase basis.
  facDiffVpar[0] = nuVtSqSum[0]; 
  facDiffVpar[4] = nuVtSqSum[4]; 
  facDiffVpar[5] = nuVtSqSum[5]; 

  double facDiffMu[48];
  // Expand mu diffusion coefficient in phase basis.
  facDiffMu[0] = rdv2[1]*w[3]*(bmag_inv[7]*nuVtSqSum[7]+bmag_inv[6]*nuVtSqSum[6]+bmag_inv[5]*nuVtSqSum[5]+bmag_inv[4]*nuVtSqSum[4]+bmag_inv[3]*nuVtSqSum[3]+bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0]); 
  facDiffMu[4] = 0.5773502691896258*(bmag_inv[7]*nuVtSqSum[7]+bmag_inv[6]*nuVtSqSum[6]+bmag_inv[5]*nuVtSqSum[5]+bmag_inv[4]*nuVtSqSum[4]+bmag_inv[3]*nuVtSqSum[3]+bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0]); 
  facDiffMu[5] = rdv2[1]*w[3]*((0.8*bmag_inv[6]+0.8944271909999159*bmag_inv[2])*nuVtSqSum[7]+0.8*nuVtSqSum[6]*bmag_inv[7]+0.8944271909999159*(nuVtSqSum[2]*bmag_inv[7]+bmag_inv[1]*nuVtSqSum[6]+nuVtSqSum[1]*bmag_inv[6]+bmag_inv[3]*nuVtSqSum[5]+nuVtSqSum[3]*bmag_inv[5]+bmag_inv[3]*nuVtSqSum[4])+nuVtSqSum[3]*(0.8944271909999159*bmag_inv[4]+bmag_inv[0])+nuVtSqSum[0]*bmag_inv[3]+bmag_inv[1]*nuVtSqSum[2]+nuVtSqSum[1]*bmag_inv[2]); 

  return fabs(rdv2[1]*(1.8*facDiffMu[0]-2.012461179749811*(facDiffMu[5]+facDiffMu[4]))*m_+rdvSq4[0]*(0.9*facDiffVpar[0]-1.006230589874905*(facDiffVpar[5]+facDiffVpar[4]))); 

} 
