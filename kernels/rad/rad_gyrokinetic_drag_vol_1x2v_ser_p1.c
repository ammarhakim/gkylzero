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
  alphaDrag[0] = 0.7071067811865475*rdv2[0]*nI[1]*vnu[1]+0.7071067811865475*nI[0]*rdv2[0]*vnu[0]; 
  alphaDrag[1] = 0.7071067811865475*nI[0]*rdv2[0]*vnu[1]+0.7071067811865475*rdv2[0]*vnu[0]*nI[1]; 
  alphaDrag[2] = 0.7071067811865475*rdv2[0]*nI[1]*vnu[4]+0.7071067811865475*nI[0]*rdv2[0]*vnu[2]; 
  alphaDrag[3] = 0.7071067811865475*rdv2[0]*nI[1]*vnu[5]+0.7071067811865475*nI[0]*rdv2[0]*vnu[3]; 
  alphaDrag[4] = 0.7071067811865475*nI[0]*rdv2[0]*vnu[4]+0.7071067811865475*rdv2[0]*nI[1]*vnu[2]; 
  alphaDrag[5] = 0.7071067811865475*nI[0]*rdv2[0]*vnu[5]+0.7071067811865475*rdv2[0]*nI[1]*vnu[3]; 
  alphaDrag[6] = 0.7071067811865475*rdv2[0]*nI[1]*vnu[7]+0.7071067811865475*nI[0]*rdv2[0]*vnu[6]; 
  alphaDrag[7] = 0.7071067811865475*nI[0]*rdv2[0]*vnu[7]+0.7071067811865475*rdv2[0]*nI[1]*vnu[6]; 
  alphaDrag[8] = 0.7071067811865475*rdv2[0]*nI[1]*vnu[9]+0.7071067811865475*nI[0]*rdv2[0]*vnu[8]; 
  alphaDrag[9] = 0.7071067811865475*nI[0]*rdv2[0]*vnu[9]+0.7071067811865475*rdv2[0]*nI[1]*vnu[8]; 
  alphaDrag[10] = 0.7071067811865475*rdv2[0]*nI[1]*vnu[11]+0.7071067811865475*nI[0]*rdv2[0]*vnu[10]; 
  alphaDrag[11] = 0.7071067811865475*nI[0]*rdv2[0]*vnu[11]+0.7071067811865475*rdv2[0]*nI[1]*vnu[10]; 

  alphaDrag[0] += 0.7071067811865475*nI[1]*rdv2[1]*vsqnu[1]+0.7071067811865475*nI[0]*vsqnu[0]*rdv2[1]; 
  alphaDrag[1] += 0.7071067811865475*nI[0]*rdv2[1]*vsqnu[1]+0.7071067811865475*vsqnu[0]*nI[1]*rdv2[1]; 
  alphaDrag[2] += 0.7071067811865475*nI[1]*rdv2[1]*vsqnu[4]+0.7071067811865475*nI[0]*rdv2[1]*vsqnu[2]; 
  alphaDrag[3] += 0.7071067811865475*nI[1]*rdv2[1]*vsqnu[5]+0.7071067811865475*nI[0]*rdv2[1]*vsqnu[3]; 
  alphaDrag[4] += 0.7071067811865475*nI[0]*rdv2[1]*vsqnu[4]+0.7071067811865475*nI[1]*rdv2[1]*vsqnu[2]; 
  alphaDrag[5] += 0.7071067811865475*nI[0]*rdv2[1]*vsqnu[5]+0.7071067811865475*nI[1]*rdv2[1]*vsqnu[3]; 
  alphaDrag[6] += 0.7071067811865475*nI[1]*rdv2[1]*vsqnu[7]+0.7071067811865475*nI[0]*rdv2[1]*vsqnu[6]; 
  alphaDrag[7] += 0.7071067811865475*nI[0]*rdv2[1]*vsqnu[7]+0.7071067811865475*nI[1]*rdv2[1]*vsqnu[6]; 
  alphaDrag[8] += 0.7071067811865475*nI[1]*rdv2[1]*vsqnu[9]+0.7071067811865475*nI[0]*rdv2[1]*vsqnu[8]; 
  alphaDrag[9] += 0.7071067811865475*nI[0]*rdv2[1]*vsqnu[9]+0.7071067811865475*nI[1]*rdv2[1]*vsqnu[8]; 
  alphaDrag[10] += 0.7071067811865475*nI[1]*rdv2[1]*vsqnu[11]+0.7071067811865475*nI[0]*rdv2[1]*vsqnu[10]; 
  alphaDrag[11] += 0.7071067811865475*nI[0]*rdv2[1]*vsqnu[11]+0.7071067811865475*nI[1]*rdv2[1]*vsqnu[10]; 

  out[2] += 0.6123724356957944*(alphaDrag[11]*f[11]+alphaDrag[10]*f[10]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(alphaDrag[11]*f[11]+alphaDrag[10]*f[10]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphaDrag[10]*f[11]+f[10]*alphaDrag[11]+alphaDrag[8]*f[9]+f[8]*alphaDrag[9]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.6123724356957944*(alphaDrag[10]*f[11]+f[10]*alphaDrag[11]+alphaDrag[8]*f[9]+f[8]*alphaDrag[9]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += (0.6123724356957944*alphaDrag[9]+0.5477225575051661*alphaDrag[7])*f[11]+(0.6123724356957944*f[9]+0.5477225575051661*f[7])*alphaDrag[11]+(0.6123724356957944*alphaDrag[8]+0.5477225575051661*alphaDrag[6])*f[10]+0.6123724356957944*f[8]*alphaDrag[10]+0.5477225575051661*(f[6]*alphaDrag[10]+alphaDrag[4]*f[9]+f[4]*alphaDrag[9]+alphaDrag[2]*f[8]+f[2]*alphaDrag[8])+0.6123724356957944*((alphaDrag[5]+alphaDrag[4])*f[7]+(f[5]+f[4])*alphaDrag[7]+(alphaDrag[3]+alphaDrag[2])*f[6]+(f[3]+f[2])*alphaDrag[6]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[7] += (0.6123724356957944*alphaDrag[8]+0.5477225575051661*alphaDrag[6])*f[11]+(0.6123724356957944*f[8]+0.5477225575051661*f[6])*alphaDrag[11]+(0.6123724356957944*alphaDrag[9]+0.5477225575051661*alphaDrag[7])*f[10]+0.6123724356957944*f[9]*alphaDrag[10]+0.5477225575051661*(f[7]*alphaDrag[10]+alphaDrag[2]*f[9]+f[2]*alphaDrag[9]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+0.6123724356957944*((alphaDrag[3]+alphaDrag[2])*f[7]+(f[3]+f[2])*alphaDrag[7]+(alphaDrag[5]+alphaDrag[4])*f[6]+(f[5]+f[4])*alphaDrag[6]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[8] += 1.224744871391589*(alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[6]*f[10]+f[6]*alphaDrag[10]+alphaDrag[4]*f[9]+f[4]*alphaDrag[9]+alphaDrag[2]*f[8]+f[2]*alphaDrag[8])+1.369306393762915*(alphaDrag[5]*f[7]+f[5]*alphaDrag[7]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 1.224744871391589*(alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[7]*f[10]+f[7]*alphaDrag[10]+alphaDrag[2]*f[9]+f[2]*alphaDrag[9]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+1.369306393762915*(alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[10] += (0.3912303982179757*alphaDrag[11]+0.6123724356957944*alphaDrag[5]+1.224744871391589*alphaDrag[4])*f[11]+(0.6123724356957944*f[5]+1.224744871391589*f[4])*alphaDrag[11]+(0.3912303982179757*alphaDrag[10]+0.6123724356957944*alphaDrag[3]+1.224744871391589*alphaDrag[2])*f[10]+(0.6123724356957944*f[3]+1.224744871391589*f[2])*alphaDrag[10]+(0.3912303982179757*alphaDrag[9]+1.224744871391589*alphaDrag[7]+0.6123724356957944*alphaDrag[1])*f[9]+(1.224744871391589*f[7]+0.6123724356957944*f[1])*alphaDrag[9]+(0.3912303982179757*alphaDrag[8]+1.224744871391589*alphaDrag[6]+0.6123724356957944*alphaDrag[0])*f[8]+(1.224744871391589*f[6]+0.6123724356957944*f[0])*alphaDrag[8]+0.5477225575051661*alphaDrag[7]*f[7]+1.369306393762915*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.5477225575051661*alphaDrag[6]*f[6]+1.369306393762915*(alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[4]*f[5])+f[4]*(1.369306393762915*alphaDrag[5]+0.5477225575051661*alphaDrag[4])+1.369306393762915*alphaDrag[2]*f[3]+f[2]*(1.369306393762915*alphaDrag[3]+0.5477225575051661*alphaDrag[2]); 
  out[11] += (0.3912303982179757*alphaDrag[10]+0.6123724356957944*alphaDrag[3]+1.224744871391589*alphaDrag[2])*f[11]+(0.3912303982179757*f[10]+0.6123724356957944*f[3]+1.224744871391589*f[2])*alphaDrag[11]+(0.6123724356957944*alphaDrag[5]+1.224744871391589*alphaDrag[4])*f[10]+(0.6123724356957944*f[5]+1.224744871391589*f[4])*alphaDrag[10]+(0.3912303982179757*alphaDrag[8]+1.224744871391589*alphaDrag[6]+0.6123724356957944*alphaDrag[0])*f[9]+(0.3912303982179757*f[8]+1.224744871391589*f[6]+0.6123724356957944*f[0])*alphaDrag[9]+(1.224744871391589*alphaDrag[7]+0.6123724356957944*alphaDrag[1])*f[8]+(1.224744871391589*f[7]+0.6123724356957944*f[1])*alphaDrag[8]+(0.5477225575051661*alphaDrag[6]+1.369306393762915*alphaDrag[0])*f[7]+0.5477225575051661*f[6]*alphaDrag[7]+1.369306393762915*(f[0]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5])+(1.369306393762915*alphaDrag[3]+0.5477225575051661*alphaDrag[2])*f[4]+(1.369306393762915*f[3]+0.5477225575051661*f[2])*alphaDrag[4]; 

  return fabs(0.8838834764831842*alphaDrag[0]-0.9882117688026182*alphaDrag[8])+fabs(0.5303300858899105*alphaDrag[0]-0.592927061281571*alphaDrag[8]); 

} 
