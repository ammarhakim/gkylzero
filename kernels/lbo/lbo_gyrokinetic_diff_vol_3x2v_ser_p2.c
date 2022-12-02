#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_vol_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[5]: cell-center coordinates. 
  // dxv[5]: cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fin: input distribution function.
  // out: incremented output 

  const double *nuVtSqSum = &nuPrimMomsSum[20];

  double rdv2[2]; 
  double rdvSq4[2]; 
  rdv2[0]   = 2.0/dxv[3]; 
  rdvSq4[0] = rdv2[0]*rdv2[0]; 
  rdv2[1]   = 2.0/dxv[4]; 
  rdvSq4[1] = rdv2[1]*rdv2[1]; 

  double facDiffVpar[20];
  // Expand nuVtSqSum in phase basis.
  facDiffVpar[0] = nuVtSqSum[0]; 
  facDiffVpar[7] = nuVtSqSum[7]; 
  facDiffVpar[8] = nuVtSqSum[8]; 
  facDiffVpar[9] = nuVtSqSum[9]; 

  double facDiffMu[112];
  // Expand mu diffusion coefficient in phase basis.
  facDiffMu[0] = 0.7071067811865475*rdv2[1]*w[4]*(bmag_inv[15]*nuVtSqSum[15]+bmag_inv[13]*nuVtSqSum[13]+bmag_inv[9]*nuVtSqSum[9]+bmag_inv[7]*nuVtSqSum[7]+bmag_inv[5]*nuVtSqSum[5]+bmag_inv[3]*nuVtSqSum[3]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0]); 
  facDiffMu[7] = rdv2[1]*w[4]*((0.5656854249492381*bmag_inv[13]+0.6324555320336759*bmag_inv[3])*nuVtSqSum[15]+0.5656854249492381*nuVtSqSum[13]*bmag_inv[15]+0.6324555320336759*(nuVtSqSum[3]*bmag_inv[15]+bmag_inv[1]*nuVtSqSum[13]+nuVtSqSum[1]*bmag_inv[13]+bmag_inv[5]*nuVtSqSum[9]+nuVtSqSum[5]*bmag_inv[9]+bmag_inv[5]*nuVtSqSum[7]+nuVtSqSum[5]*bmag_inv[7])+0.7071067811865475*(bmag_inv[0]*nuVtSqSum[5]+nuVtSqSum[0]*bmag_inv[5]+bmag_inv[1]*nuVtSqSum[3]+nuVtSqSum[1]*bmag_inv[3])); 
  facDiffMu[8] = rdv2[1]*w[4]*(0.6324555320336759*bmag_inv[5]*nuVtSqSum[19]+0.7071067811865475*bmag_inv[7]*nuVtSqSum[17]+0.6324555320336759*(bmag_inv[3]*nuVtSqSum[16]+nuVtSqSum[10]*bmag_inv[15])+0.7071067811865475*(nuVtSqSum[11]*bmag_inv[13]+bmag_inv[1]*nuVtSqSum[10])+0.6324555320336759*nuVtSqSum[6]*bmag_inv[9]+0.7071067811865475*(bmag_inv[0]*nuVtSqSum[6]+nuVtSqSum[4]*bmag_inv[5]+nuVtSqSum[2]*bmag_inv[3])); 

  return fabs(rdv2[1]*(6.363961030678928*facDiffMu[0]-7.115124735378852*(facDiffMu[8]+facDiffMu[7]))*m_+rdvSq4[0]*(3.181980515339463*facDiffVpar[0]-3.557562367689425*(facDiffVpar[9]+facDiffVpar[8]+facDiffVpar[7]))); 

} 
