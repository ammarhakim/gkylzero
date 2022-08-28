#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *nu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]:      Cell-center coordinates. 
  // dxv[2]:    Cell spacing. 
  // nu:     collisionality. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  double alphaDrag[8]; 
  // Expand dv1par*(nu*vx) in phase basis.
  alphaDrag[0] = -1.414213562373095*nu[0]*dv1par*wvpar; 
  alphaDrag[1] = -1.414213562373095*nu[1]*dv1par*wvpar; 
  alphaDrag[2] = -0.408248290463863*nu[0]*dv1par*dvpar; 
  alphaDrag[3] = -0.408248290463863*nu[1]*dv1par*dvpar; 
  alphaDrag[4] = -1.414213562373095*nu[2]*dv1par*wvpar; 
  alphaDrag[6] = -0.408248290463863*nu[2]*dv1par*dvpar; 

  out[2] += 0.8660254037844386*(alphaDrag[6]*f[6]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.7745966692414833*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 1.732050807568877*alphaDrag[3]*f[7]+1.936491673103709*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6])+1.732050807568877*alphaDrag[2]*f[5]+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.5532833351724881*alphaDrag[6]*f[6]+0.8660254037844386*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6])+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*(alphaDrag[3]*f[3]+alphaDrag[1]*f[1]); 
  out[7] += 1.549193338482967*alphaDrag[6]*f[7]+1.732050807568877*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[3]*(f[5]+f[4])+f[3]*alphaDrag[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 

  return fabs(0.25*alphaDrag[0]-0.2795084971874737*alphaDrag[4]); 

} 
