#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *nI, const double *vnu, const double *vsqnu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // nI: atomic density 
  // vnu: 2/pi*v*nu(v) dg field representation (v'(v||,mu) in notes) 
  // vsqnu: sqrt(mu*me/2B)*v^2*nu(v) dg field representation (v''(v||,mu) in notes) 
  // f: input distribution function.
  // out: incremented output 

  double rdv2[2]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdv2[1] = 2.0/dxv[2]; 

  double alphaDrag[12] = {0,0}; 
  alphaDrag[0] = 0.7071067811865475*nI[1]*vnu[1]+0.7071067811865475*nI[0]*vnu[0]; 
  alphaDrag[1] = 0.7071067811865475*nI[0]*vnu[1]+0.7071067811865475*vnu[0]*nI[1]; 
  alphaDrag[2] = 0.7071067811865475*nI[1]*vnu[4]+0.7071067811865475*nI[0]*vnu[2]; 
  alphaDrag[3] = 0.7071067811865475*nI[1]*vnu[5]+0.7071067811865475*nI[0]*vnu[3]; 
  alphaDrag[4] = 0.7071067811865475*nI[0]*vnu[4]+0.7071067811865475*nI[1]*vnu[2]; 
  alphaDrag[5] = 0.7071067811865475*nI[0]*vnu[5]+0.7071067811865475*nI[1]*vnu[3]; 
  alphaDrag[6] = 0.7071067811865475*nI[1]*vnu[7]+0.7071067811865475*nI[0]*vnu[6]; 
  alphaDrag[7] = 0.7071067811865475*nI[0]*vnu[7]+0.7071067811865475*nI[1]*vnu[6]; 
  alphaDrag[8] = 0.7071067811865475*nI[1]*vnu[9]+0.7071067811865475*nI[0]*vnu[8]; 
  alphaDrag[9] = 0.7071067811865475*nI[0]*vnu[9]+0.7071067811865475*nI[1]*vnu[8]; 
  alphaDrag[10] = 0.7071067811865475*nI[1]*vnu[11]+0.7071067811865475*nI[0]*vnu[10]; 
  alphaDrag[11] = 0.7071067811865475*nI[0]*vnu[11]+0.7071067811865475*nI[1]*vnu[10]; 

  alphaDrag[0] += 0.7071067811865475*nI[1]*vsqnu[1]+0.7071067811865475*nI[0]*vsqnu[0]; 
  alphaDrag[1] += 0.7071067811865475*nI[0]*vsqnu[1]+0.7071067811865475*vsqnu[0]*nI[1]; 
  alphaDrag[2] += 0.7071067811865475*nI[1]*vsqnu[4]+0.7071067811865475*nI[0]*vsqnu[2]; 
  alphaDrag[3] += 0.7071067811865475*nI[1]*vsqnu[5]+0.7071067811865475*nI[0]*vsqnu[3]; 
  alphaDrag[4] += 0.7071067811865475*nI[0]*vsqnu[4]+0.7071067811865475*nI[1]*vsqnu[2]; 
  alphaDrag[5] += 0.7071067811865475*nI[0]*vsqnu[5]+0.7071067811865475*nI[1]*vsqnu[3]; 
  alphaDrag[6] += 0.7071067811865475*nI[1]*vsqnu[7]+0.7071067811865475*nI[0]*vsqnu[6]; 
  alphaDrag[7] += 0.7071067811865475*nI[0]*vsqnu[7]+0.7071067811865475*nI[1]*vsqnu[6]; 
  alphaDrag[8] += 0.7071067811865475*nI[1]*vsqnu[9]+0.7071067811865475*nI[0]*vsqnu[8]; 
  alphaDrag[9] += 0.7071067811865475*nI[0]*vsqnu[9]+0.7071067811865475*nI[1]*vsqnu[8]; 
  alphaDrag[10] += 0.7071067811865475*nI[1]*vsqnu[11]+0.7071067811865475*nI[0]*vsqnu[10]; 
  alphaDrag[11] += 0.7071067811865475*nI[0]*vsqnu[11]+0.7071067811865475*nI[1]*vsqnu[10]; 

  out[2] += 1.224744871391589*(alphaDrag[11]*f[11]+alphaDrag[10]*f[10]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 1.224744871391589*(alphaDrag[10]*f[11]+f[10]*alphaDrag[11]+alphaDrag[8]*f[9]+f[8]*alphaDrag[9]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 1.224744871391589*(alphaDrag[9]*f[11]+f[9]*alphaDrag[11]+alphaDrag[8]*f[10]+f[8]*alphaDrag[10]+alphaDrag[4]*f[7]+f[4]*alphaDrag[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[7] += 1.224744871391589*(alphaDrag[8]*f[11]+f[8]*alphaDrag[11]+alphaDrag[9]*f[10]+f[9]*alphaDrag[10]+alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[4]*f[6]+f[4]*alphaDrag[6]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[8] += 2.449489742783178*(alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[6]*f[10]+f[6]*alphaDrag[10]+alphaDrag[4]*f[9]+f[4]*alphaDrag[9]+alphaDrag[2]*f[8]+f[2]*alphaDrag[8])+2.738612787525831*(alphaDrag[5]*f[7]+f[5]*alphaDrag[7]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 2.449489742783178*(alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[7]*f[10]+f[7]*alphaDrag[10]+alphaDrag[2]*f[9]+f[2]*alphaDrag[9]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+2.738612787525831*(alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[10] += 2.449489742783178*(alphaDrag[4]*f[11]+f[4]*alphaDrag[11]+alphaDrag[2]*f[10]+f[2]*alphaDrag[10]+alphaDrag[7]*f[9]+f[7]*alphaDrag[9]+alphaDrag[6]*f[8]+f[6]*alphaDrag[8])+2.738612787525831*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[4]*f[5]+f[4]*alphaDrag[5]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[11] += 2.449489742783178*(alphaDrag[2]*f[11]+f[2]*alphaDrag[11]+alphaDrag[4]*f[10]+f[4]*alphaDrag[10]+alphaDrag[6]*f[9]+f[6]*alphaDrag[9]+alphaDrag[7]*f[8]+f[7]*alphaDrag[8])+2.738612787525831*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4]); 

  return fabs(0.8838834764831842*alphaDrag[0]-0.9882117688026182*alphaDrag[8])+fabs(nI[0]*(0.375*vsqnu[0]-0.4192627457812106*vsqnu[8])); 

} 
