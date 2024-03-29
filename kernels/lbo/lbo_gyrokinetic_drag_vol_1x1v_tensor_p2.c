#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]: cell-center coordinates. 
  // dxv[2]: cell spacing. 
  // m_: species mass.
  // bmag_inv: 1/(magnetic field magnitude). 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // f: input distribution function.
  // out: incremented output 

  const double *nuUSum = nuPrimMomsSum;

  double rdv2[1]; 
  rdv2[0] = 2.0/dxv[1]; 

  double alphaDrag[9]; 
  // Expand rdv2*(nu*vpar-nuUparSum) in phase basis.
  alphaDrag[0] = rdv2[0]*(1.414213562373095*nuUSum[0]-1.414213562373095*nuSum[0]*w[1]); 
  alphaDrag[1] = rdv2[0]*(1.414213562373095*nuUSum[1]-1.414213562373095*nuSum[1]*w[1]); 
  alphaDrag[2] = -0.408248290463863*nuSum[0]*rdv2[0]*dxv[1]; 
  alphaDrag[3] = -0.408248290463863*rdv2[0]*dxv[1]*nuSum[1]; 
  alphaDrag[4] = rdv2[0]*(1.414213562373095*nuUSum[2]-1.414213562373095*w[1]*nuSum[2]); 
  alphaDrag[6] = -0.408248290463863*rdv2[0]*dxv[1]*nuSum[2]; 

  out[2] += 0.8660254037844386*(alphaDrag[6]*f[6]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.7745966692414833*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 1.732050807568877*(alphaDrag[6]*f[8]+alphaDrag[3]*f[7])+1.936491673103709*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6])+1.732050807568877*alphaDrag[2]*f[5]+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.5532833351724881*alphaDrag[6]*f[6]+0.8660254037844386*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6])+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*(alphaDrag[3]*f[3]+alphaDrag[1]*f[1]); 
  out[7] += 1.549193338482967*(alphaDrag[3]*f[8]+alphaDrag[6]*f[7])+1.732050807568877*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[3]*(f[5]+f[4])+f[3]*alphaDrag[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[8] += (1.106566670344976*alphaDrag[6]+1.732050807568877*alphaDrag[2])*f[8]+1.549193338482967*alphaDrag[3]*f[7]+(1.237179148263484*alphaDrag[4]+1.936491673103709*alphaDrag[0])*f[6]+(1.732050807568877*f[5]+1.237179148263484*f[4])*alphaDrag[6]+1.936491673103709*(f[0]*alphaDrag[6]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4])+1.732050807568877*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 

  return fabs(1.25*alphaDrag[0]-1.397542485937369*alphaDrag[4]); 

} 
