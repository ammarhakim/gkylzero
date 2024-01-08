#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vnu, const double *vsqnu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]: cell-center coordinates. 
  // dxv[3]: cell spacing. 
  // vnu: 2/pi*v*nu(v) dg field representation (v'(v||,mu) in notes) 
  // vsqnu: sqrt(mu*me/2B)*v^2*nu(v) dg field representation (v''(v||,mu) in notes) 
  // f: input distribution function.
  // out: incremented output 

  double rdv2[2]; 
  rdv2[0] = 2.0/dxv[1]; 
  rdv2[1] = 2.0/dxv[2]; 

  double alphaDrag[12] = {0,0}; 
  alphaDrag[0] = rdv2[0]*vnu[0]; 
  alphaDrag[1] = rdv2[0]*vnu[1]; 
  alphaDrag[2] = rdv2[0]*vnu[2]; 
  alphaDrag[3] = rdv2[0]*vnu[3]; 
  alphaDrag[4] = rdv2[0]*vnu[4]; 
  alphaDrag[5] = rdv2[0]*vnu[5]; 
  alphaDrag[6] = rdv2[0]*vnu[6]; 
  alphaDrag[7] = rdv2[0]*vnu[7]; 
  alphaDrag[8] = rdv2[0]*vnu[8]; 
  alphaDrag[9] = rdv2[0]*vnu[9]; 
  alphaDrag[10] = rdv2[0]*vnu[10]; 
  alphaDrag[11] = rdv2[0]*vnu[11]; 

  betaDrag[0] += vsqnu[0]*rdv2[1]; 
  betaDrag[1] += rdv2[1]*vsqnu[1]; 
  betaDrag[2] += rdv2[1]*vsqnu[2]; 
  betaDrag[3] += rdv2[1]*vsqnu[3]; 
  betaDrag[4] += rdv2[1]*vsqnu[4]; 
  betaDrag[5] += rdv2[1]*vsqnu[5]; 
  betaDrag[6] += rdv2[1]*vsqnu[6]; 
  betaDrag[7] += rdv2[1]*vsqnu[7]; 
  betaDrag[8] += rdv2[1]*vsqnu[8]; 
  betaDrag[9] += rdv2[1]*vsqnu[9]; 
  betaDrag[10] += rdv2[1]*vsqnu[10]; 
  betaDrag[11] += rdv2[1]*vsqnu[11]; 

  out[2] += 0.6123724356957944*(alphaDrag[11]*f[11]+alphaDrag[10]*f[10]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(betaDrag[11]*f[11]+betaDrag[10]*f[10]+betaDrag[9]*f[9]+betaDrag[8]*f[8]+betaDrag[7]*f[7]+betaDrag[6]*f[6]+betaDrag[5]*f[5]+betaDrag[4]*f[4]+betaDrag[3]*f[3]+betaDrag[2]*f[2]+betaDrag[1]*f[1]+betaDrag[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphaDrag[10]*f[11]+f[10]*alphaDrag[11]+alphaDrag[8]*f[9]+f[8]*alphaDrag[9]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.6123724356957944*(betaDrag[10]*f[11]+f[10]*betaDrag[11]+betaDrag[8]*f[9]+f[8]*betaDrag[9]+betaDrag[6]*f[7]+f[6]*betaDrag[7]+betaDrag[3]*f[5]+f[3]*betaDrag[5]+betaDrag[2]*f[4]+f[2]*betaDrag[4]+betaDrag[0]*f[1]+f[0]*betaDrag[1]); 
  out[6] += 0.6123724356957944*alphaDrag[9]*f[11]+0.5477225575051661*(betaDrag[7]*f[11]+f[7]*betaDrag[11])+0.6123724356957944*(f[9]*alphaDrag[11]+alphaDrag[8]*f[10])+0.5477225575051661*(betaDrag[6]*f[10]+f[6]*betaDrag[10])+0.6123724356957944*f[8]*alphaDrag[10]+0.5477225575051661*(betaDrag[4]*f[9]+f[4]*betaDrag[9]+betaDrag[2]*f[8]+f[2]*betaDrag[8])+0.6123724356957944*((betaDrag[5]+alphaDrag[4])*f[7]+f[5]*betaDrag[7]+f[4]*alphaDrag[7]+(betaDrag[3]+alphaDrag[2])*f[6]+f[3]*betaDrag[6]+f[2]*alphaDrag[6]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+betaDrag[1]*f[4]+f[1]*betaDrag[4]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+betaDrag[0]*f[2]+f[0]*betaDrag[2]); 
  out[7] += 0.6123724356957944*alphaDrag[8]*f[11]+0.5477225575051661*(betaDrag[6]*f[11]+f[6]*betaDrag[11])+0.6123724356957944*(f[8]*alphaDrag[11]+alphaDrag[9]*f[10])+0.5477225575051661*(betaDrag[7]*f[10]+f[7]*betaDrag[10])+0.6123724356957944*f[9]*alphaDrag[10]+0.5477225575051661*(betaDrag[2]*f[9]+f[2]*betaDrag[9]+betaDrag[4]*f[8]+f[4]*betaDrag[8])+0.6123724356957944*((betaDrag[3]+alphaDrag[2])*f[7]+f[3]*betaDrag[7]+f[2]*alphaDrag[7]+(betaDrag[5]+alphaDrag[4])*f[6]+f[5]*betaDrag[6]+f[4]*alphaDrag[6]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+betaDrag[0]*f[4]+f[0]*betaDrag[4]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+betaDrag[1]*f[2]+f[1]*betaDrag[2]); 
  out[8] += 1.224744871391589*(alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[6]*f[10]+f[6]*alphaDrag[10]+alphaDrag[4]*f[9]+f[4]*alphaDrag[9]+alphaDrag[2]*f[8]+f[2]*alphaDrag[8])+1.369306393762915*(alphaDrag[5]*f[7]+f[5]*alphaDrag[7]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 1.224744871391589*(alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[7]*f[10]+f[7]*alphaDrag[10]+alphaDrag[2]*f[9]+f[2]*alphaDrag[9]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+1.369306393762915*(alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[10] += (0.3912303982179757*betaDrag[11]+0.6123724356957944*betaDrag[5]+1.224744871391589*alphaDrag[4])*f[11]+0.6123724356957944*f[5]*betaDrag[11]+1.224744871391589*f[4]*alphaDrag[11]+(0.3912303982179757*betaDrag[10]+0.6123724356957944*betaDrag[3]+1.224744871391589*alphaDrag[2])*f[10]+0.6123724356957944*f[3]*betaDrag[10]+1.224744871391589*f[2]*alphaDrag[10]+(0.3912303982179757*betaDrag[9]+1.224744871391589*alphaDrag[7])*f[9]+0.6123724356957944*(betaDrag[1]*f[9]+f[1]*betaDrag[9])+1.224744871391589*f[7]*alphaDrag[9]+(0.3912303982179757*betaDrag[8]+1.224744871391589*alphaDrag[6])*f[8]+0.6123724356957944*(betaDrag[0]*f[8]+f[0]*betaDrag[8])+1.224744871391589*f[6]*alphaDrag[8]+0.5477225575051661*betaDrag[7]*f[7]+1.369306393762915*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.5477225575051661*betaDrag[6]*f[6]+1.369306393762915*(alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[4]*f[5])+f[4]*(1.369306393762915*alphaDrag[5]+0.5477225575051661*betaDrag[4])+1.369306393762915*alphaDrag[2]*f[3]+f[2]*(1.369306393762915*alphaDrag[3]+0.5477225575051661*betaDrag[2]); 
  out[11] += (0.3912303982179757*betaDrag[10]+0.6123724356957944*betaDrag[3]+1.224744871391589*alphaDrag[2])*f[11]+(0.3912303982179757*f[10]+0.6123724356957944*f[3])*betaDrag[11]+1.224744871391589*f[2]*alphaDrag[11]+(0.6123724356957944*betaDrag[5]+1.224744871391589*alphaDrag[4])*f[10]+0.6123724356957944*f[5]*betaDrag[10]+1.224744871391589*f[4]*alphaDrag[10]+(0.3912303982179757*betaDrag[8]+1.224744871391589*alphaDrag[6]+0.6123724356957944*betaDrag[0])*f[9]+(0.3912303982179757*f[8]+0.6123724356957944*f[0])*betaDrag[9]+1.224744871391589*(f[6]*alphaDrag[9]+alphaDrag[7]*f[8])+0.6123724356957944*(betaDrag[1]*f[8]+f[1]*betaDrag[8])+f[7]*(1.224744871391589*alphaDrag[8]+0.5477225575051661*betaDrag[6]+1.369306393762915*alphaDrag[0])+0.5477225575051661*f[6]*betaDrag[7]+1.369306393762915*(f[0]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[3]*f[4])+0.5477225575051661*(betaDrag[2]*f[4]+f[2]*betaDrag[4])+1.369306393762915*f[3]*alphaDrag[4]; 

  return fabs(0.8838834764831842*alphaDrag[0]-0.9882117688026182*alphaDrag[8])+fabs(0.5303300858899105*betaDrag[0]-0.592927061281571*betaDrag[8]); 

} 
