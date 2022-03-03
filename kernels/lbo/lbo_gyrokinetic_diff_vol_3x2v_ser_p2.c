#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_diff_vol_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // w[5]:      cell-center coordinates. 
  // dxv[5]:    cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         input distribution function.
  // out:       incremented output 
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
  facDiffMu[0] = 0.7071067811865475*rdv2[1]*w[4]*(bmag_inv[19]*nuVtSqSum[19]+bmag_inv[18]*nuVtSqSum[18]+bmag_inv[17]*nuVtSqSum[17]+bmag_inv[16]*nuVtSqSum[16]+bmag_inv[15]*nuVtSqSum[15]+bmag_inv[14]*nuVtSqSum[14]+bmag_inv[13]*nuVtSqSum[13]+bmag_inv[12]*nuVtSqSum[12]+bmag_inv[11]*nuVtSqSum[11]+bmag_inv[10]*nuVtSqSum[10]+bmag_inv[9]*nuVtSqSum[9]+bmag_inv[8]*nuVtSqSum[8]+bmag_inv[7]*nuVtSqSum[7]+bmag_inv[6]*nuVtSqSum[6]+bmag_inv[5]*nuVtSqSum[5]+bmag_inv[4]*nuVtSqSum[4]+bmag_inv[3]*nuVtSqSum[3]+bmag_inv[2]*nuVtSqSum[2]+bmag_inv[1]*nuVtSqSum[1]+bmag_inv[0]*nuVtSqSum[0]); 
  facDiffMu[7] = rdv2[1]*w[4]*((0.5656854249492381*bmag_inv[17]+0.6324555320336759*bmag_inv[6])*nuVtSqSum[19]+(0.5656854249492381*nuVtSqSum[17]+0.6324555320336759*nuVtSqSum[6])*bmag_inv[19]+0.7071067811865475*(bmag_inv[8]*nuVtSqSum[18]+nuVtSqSum[8]*bmag_inv[18])+0.6324555320336759*(bmag_inv[4]*nuVtSqSum[17]+nuVtSqSum[4]*bmag_inv[17]+bmag_inv[10]*nuVtSqSum[16]+nuVtSqSum[10]*bmag_inv[16])+(0.5656854249492381*bmag_inv[13]+0.6324555320336759*bmag_inv[3])*nuVtSqSum[15]+(0.5656854249492381*nuVtSqSum[13]+0.6324555320336759*nuVtSqSum[3])*bmag_inv[15]+0.7071067811865475*(bmag_inv[12]*nuVtSqSum[14]+nuVtSqSum[12]*bmag_inv[14])+0.6324555320336759*(bmag_inv[1]*nuVtSqSum[13]+nuVtSqSum[1]*bmag_inv[13]+bmag_inv[10]*nuVtSqSum[11]+nuVtSqSum[10]*bmag_inv[11])+0.7071067811865475*(bmag_inv[2]*nuVtSqSum[10]+nuVtSqSum[2]*bmag_inv[10])+0.6324555320336759*(bmag_inv[5]*nuVtSqSum[9]+nuVtSqSum[5]*bmag_inv[9]+bmag_inv[5]*nuVtSqSum[7]+nuVtSqSum[5]*bmag_inv[7])+0.7071067811865475*(bmag_inv[4]*nuVtSqSum[6]+nuVtSqSum[4]*bmag_inv[6]+bmag_inv[0]*nuVtSqSum[5]+nuVtSqSum[0]*bmag_inv[5]+bmag_inv[1]*nuVtSqSum[3]+nuVtSqSum[1]*bmag_inv[3])); 
  facDiffMu[8] = rdv2[1]*w[4]*((0.5656854249492381*bmag_inv[18]+0.6324555320336759*bmag_inv[5])*nuVtSqSum[19]+0.5656854249492381*nuVtSqSum[18]*bmag_inv[19]+0.6324555320336759*(nuVtSqSum[5]*bmag_inv[19]+bmag_inv[4]*nuVtSqSum[18]+nuVtSqSum[4]*bmag_inv[18])+0.7071067811865475*(bmag_inv[7]*nuVtSqSum[17]+nuVtSqSum[7]*bmag_inv[17])+(0.5656854249492381*bmag_inv[14]+0.6324555320336759*bmag_inv[3])*nuVtSqSum[16]+0.5656854249492381*nuVtSqSum[14]*bmag_inv[16]+0.6324555320336759*(nuVtSqSum[3]*bmag_inv[16]+bmag_inv[10]*nuVtSqSum[15]+nuVtSqSum[10]*bmag_inv[15]+bmag_inv[2]*nuVtSqSum[14]+nuVtSqSum[2]*bmag_inv[14])+0.7071067811865475*(bmag_inv[11]*nuVtSqSum[13]+nuVtSqSum[11]*bmag_inv[13])+0.6324555320336759*(bmag_inv[10]*nuVtSqSum[12]+nuVtSqSum[10]*bmag_inv[12])+0.7071067811865475*(bmag_inv[1]*nuVtSqSum[10]+nuVtSqSum[1]*bmag_inv[10])+0.6324555320336759*(bmag_inv[6]*nuVtSqSum[9]+nuVtSqSum[6]*bmag_inv[9]+bmag_inv[6]*nuVtSqSum[8]+nuVtSqSum[6]*bmag_inv[8])+0.7071067811865475*(bmag_inv[0]*nuVtSqSum[6]+nuVtSqSum[0]*bmag_inv[6]+bmag_inv[4]*nuVtSqSum[5]+nuVtSqSum[4]*bmag_inv[5]+bmag_inv[2]*nuVtSqSum[3]+nuVtSqSum[2]*bmag_inv[3])); 

  return fabs(rdv2[1]*(1.272792206135785*facDiffMu[0]-1.42302494707577*(facDiffMu[8]+facDiffMu[7]))*m_+rdvSq4[0]*(0.6363961030678926*facDiffVpar[0]-0.711512473537885*(facDiffVpar[9]+facDiffVpar[8]+facDiffVpar[7]))); 

} 
