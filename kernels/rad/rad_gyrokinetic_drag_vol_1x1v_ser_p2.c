#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *vnu, const double *vsqnu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]: cell-center coordinates. 
  // dxv[2]: cell spacing. 
  // vnu: 2/pi*v*nu(v) dg field representation (v'(v||,mu) in notes) 
  // vsqnu: sqrt(mu*me/2B)*v^2*nu(v) dg field representation (v''(v||,mu) in notes) 
  // f: input distribution function.
  // out: incremented output 

  double rdv2[1]; 
  rdv2[0] = 2.0/dxv[1]; 

  double alphaDrag[8] = {0,0}; 
  alphaDrag[0] = rdv2[0]*vnu[0]; 
  alphaDrag[1] = rdv2[0]*vnu[1]; 
  alphaDrag[2] = rdv2[0]*vnu[2]; 
  alphaDrag[3] = rdv2[0]*vnu[3]; 
  alphaDrag[4] = rdv2[0]*vnu[4]; 
  alphaDrag[5] = rdv2[0]*vnu[5]; 
  alphaDrag[6] = rdv2[0]*vnu[6]; 
  alphaDrag[7] = rdv2[0]*vnu[7]; 

  out[2] += 0.8660254037844386*(alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaDrag[5]*f[7]+f[5]*alphaDrag[7])+0.7745966692414833*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 1.732050807568877*(alphaDrag[3]*f[7]+f[3]*alphaDrag[7])+1.936491673103709*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6])+1.732050807568877*(alphaDrag[2]*f[5]+f[2]*alphaDrag[5])+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.7745966692414833*alphaDrag[7]*f[7]+0.5532833351724881*alphaDrag[6]*f[6]+0.8660254037844386*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6])+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*(alphaDrag[3]*f[3]+alphaDrag[1]*f[1]); 
  out[7] += (1.549193338482967*alphaDrag[6]+1.732050807568877*alphaDrag[2])*f[7]+1.549193338482967*f[6]*alphaDrag[7]+1.732050807568877*(f[2]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 

  return fabs(1.25*alphaDrag[0]-1.397542485937369*(alphaDrag[5]+alphaDrag[4])); 

} 
