#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double vlasov_lbo_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]:      Cell-center coordinates. 
  // dxv[2]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 

  double alphaDrag[8]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (1.414213562373095*nuUSum[0]-1.414213562373095*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (1.414213562373095*nuUSum[1]-1.414213562373095*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -0.408248290463863*nuSum[0]*dxv[1]*rdvx2; 
  alphaDrag[3] = -0.408248290463863*dxv[1]*nuSum[1]*rdvx2; 
  alphaDrag[4] = (1.414213562373095*nuUSum[2]-1.414213562373095*w[1]*nuSum[2])*rdvx2; 
  alphaDrag[6] = -0.408248290463863*dxv[1]*nuSum[2]*rdvx2; 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.8660254037844386*(alphaDrag[6]*f[6]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.7745966692414833*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 4.743416490252569*(facDiff[2]*f[4]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+1.732050807568877*alphaDrag[3]*f[7]+1.936491673103709*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6])+1.732050807568877*alphaDrag[2]*f[5]+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.5532833351724881*alphaDrag[6]*f[6]+0.8660254037844386*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6])+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*(alphaDrag[3]*f[3]+alphaDrag[1]*f[1]); 
  out[7] += (4.242640687119286*(facDiff[1]*f[4]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4+1.549193338482967*alphaDrag[6]*f[7]+1.732050807568877*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[3]*(f[5]+f[4])+f[3]*alphaDrag[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 

  return fabs(0.25*alphaDrag[0]-0.2795084971874737*alphaDrag[4])+fabs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4); 

} 
