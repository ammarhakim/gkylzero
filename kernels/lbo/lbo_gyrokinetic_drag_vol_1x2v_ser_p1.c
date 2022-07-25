#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_vol_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]:      cell-center coordinates. 
  // dxv[3]:    cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         input distribution function.
  // out:       incremented output 
  double rdv2[2]; 
  rdv2[0]   = 2.0/dxv[1]; 
  rdv2[1]   = 2.0/dxv[2]; 

  double alphaDrag[24]; 
  // Expand rdv2*(nu*vpar-nuUparSum) in phase basis.
  alphaDrag[0] = rdv2[0]*(2.0*nuUSum[0]-2.0*nuSum[0]*w[1]); 
  alphaDrag[1] = rdv2[0]*(2.0*nuUSum[1]-2.0*nuSum[1]*w[1]); 
  alphaDrag[2] = -0.5773502691896258*nuSum[0]*rdv2[0]*dxv[1]; 
  alphaDrag[4] = -0.5773502691896258*rdv2[0]*dxv[1]*nuSum[1]; 

  // Expand rdv2*nu*2*mu in phase basis.
  alphaDrag[12] = -4.0*nuSum[0]*rdv2[1]*w[2]; 
  alphaDrag[13] = -4.0*nuSum[1]*rdv2[1]*w[2]; 
  alphaDrag[15] = -1.154700538379252*nuSum[0]*rdv2[1]*dxv[2]; 
  alphaDrag[17] = -1.154700538379252*nuSum[1]*rdv2[1]*dxv[2]; 

  out[2] += 0.6123724356957944*(alphaDrag[4]*f[4]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[5]*alphaDrag[17]+f[3]*alphaDrag[15]+f[1]*alphaDrag[13]+f[0]*alphaDrag[12]); 
  out[4] += 0.6123724356957944*(alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.6123724356957944*(f[3]*alphaDrag[17]+f[5]*alphaDrag[15]+f[0]*alphaDrag[13]+f[1]*alphaDrag[12]); 
  out[6] += 0.6123724356957944*(f[7]*alphaDrag[17]+f[6]*alphaDrag[15]+f[4]*alphaDrag[13]+f[2]*alphaDrag[12]+alphaDrag[4]*f[7]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[7] += 0.6123724356957944*(f[6]*alphaDrag[17]+f[7]*alphaDrag[15]+f[2]*alphaDrag[13]+f[4]*alphaDrag[12]+alphaDrag[2]*f[7]+alphaDrag[4]*f[6]+alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[8] += 1.224744871391589*(alphaDrag[4]*f[9]+alphaDrag[2]*f[8])+1.369306393762915*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 1.224744871391589*(alphaDrag[2]*f[9]+alphaDrag[4]*f[8])+1.369306393762915*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[10] += 0.6123724356957944*(f[11]*alphaDrag[17]+f[10]*alphaDrag[15]+f[9]*alphaDrag[13]+f[8]*alphaDrag[12])+1.224744871391589*(alphaDrag[4]*f[11]+alphaDrag[2]*f[10])+1.369306393762915*(alphaDrag[1]*f[7]+alphaDrag[0]*f[6]+alphaDrag[4]*f[5]+alphaDrag[2]*f[3]); 
  out[11] += 0.6123724356957944*(f[10]*alphaDrag[17]+f[11]*alphaDrag[15]+f[8]*alphaDrag[13]+f[9]*alphaDrag[12])+1.224744871391589*(alphaDrag[2]*f[11]+alphaDrag[4]*f[10])+1.369306393762915*(alphaDrag[0]*f[7]+alphaDrag[1]*f[6]+alphaDrag[2]*f[5]+f[3]*alphaDrag[4]); 

  return fabs(0.1767766952966368*alphaDrag[0])+fabs(0.1767766952966368*alphaDrag[12]); 

} 
