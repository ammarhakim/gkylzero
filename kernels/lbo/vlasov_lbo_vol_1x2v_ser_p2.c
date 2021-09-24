#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double vlasov_lbo_vol_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 

  double alphaDrag[40]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.828427124746191*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = 2.0*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -0.8164965809277261*nuSum[0]*dxv[1]*rdvx2; 
  alphaDrag[7] = 2.0*nuUSum[2]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[20] = (2.0*nuUSum[3]-2.828427124746191*nuSum[0]*w[2])*rdvy2; 
  alphaDrag[21] = 2.0*nuUSum[4]*rdvy2; 
  alphaDrag[23] = -0.8164965809277261*nuSum[0]*dxv[2]*rdvy2; 
  alphaDrag[27] = 2.0*nuUSum[5]*rdvy2; 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.6123724356957944*(alphaDrag[7]*f[7]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[7]*alphaDrag[27]+f[3]*alphaDrag[23]+f[1]*alphaDrag[21]+f[0]*alphaDrag[20]); 
  out[4] += 0.5477225575051661*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7])+0.6123724356957944*(alphaDrag[2]*f[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.5477225575051661*f[1]*alphaDrag[27]+0.6123724356957944*f[5]*alphaDrag[23]+0.5477225575051661*f[7]*alphaDrag[21]+0.6123724356957944*(f[0]*alphaDrag[21]+f[1]*alphaDrag[20]); 
  out[6] += 0.6123724356957944*(f[11]*alphaDrag[27]+f[6]*alphaDrag[23]+f[4]*alphaDrag[21]+f[2]*alphaDrag[20]+alphaDrag[7]*f[13]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[8] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+1.369306393762915*alphaDrag[7]*f[11]+1.224744871391589*alphaDrag[2]*f[8]+1.369306393762915*(alphaDrag[1]*f[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 4.743416490252569*(facDiff[2]*f[7]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+1.369306393762915*f[13]*alphaDrag[27]+1.224744871391589*f[9]*alphaDrag[23]+1.369306393762915*(f[0]*alphaDrag[23]+f[5]*alphaDrag[21]+f[3]*alphaDrag[20]); 
  out[10] += 0.5477225575051661*f[4]*alphaDrag[27]+0.6123724356957944*f[10]*alphaDrag[23]+0.5477225575051661*f[11]*alphaDrag[21]+0.6123724356957944*(f[2]*alphaDrag[21]+f[4]*alphaDrag[20])+0.5477225575051661*alphaDrag[1]*f[13]+0.6123724356957944*alphaDrag[2]*f[10]+0.5477225575051661*f[5]*alphaDrag[7]+0.6123724356957944*(alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[11] += 0.6123724356957944*alphaDrag[2]*f[11]+0.3912303982179757*alphaDrag[7]*f[7]+0.6123724356957944*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7])+0.5477225575051661*alphaDrag[1]*f[1]; 
  out[12] += (4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4+1.224744871391589*(alphaDrag[2]*f[12]+alphaDrag[1]*f[11]+f[4]*alphaDrag[7])+1.369306393762915*(alphaDrag[0]*f[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[13] += 0.3912303982179757*f[7]*alphaDrag[27]+0.6123724356957944*(f[0]*alphaDrag[27]+f[13]*alphaDrag[23])+0.5477225575051661*f[1]*alphaDrag[21]+0.6123724356957944*f[7]*alphaDrag[20]; 
  out[14] += 4.743416490252569*(facDiff[2]*f[13]+facDiff[1]*f[5]+facDiff[0]*f[3])*rdvxSq4+0.6123724356957944*(f[14]*alphaDrag[23]+f[12]*alphaDrag[21]+f[8]*alphaDrag[20])+1.369306393762915*alphaDrag[7]*f[17]+1.224744871391589*alphaDrag[2]*f[14]+1.369306393762915*(alphaDrag[1]*f[10]+alphaDrag[0]*f[6]+alphaDrag[2]*f[3]); 
  out[15] += (4.242640687119286*(facDiff[1]*f[7]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvySq4+1.224744871391589*f[5]*alphaDrag[27]+(1.224744871391589*f[15]+1.369306393762915*f[1])*alphaDrag[23]+1.224744871391589*f[13]*alphaDrag[21]+1.369306393762915*(f[3]*alphaDrag[21]+f[5]*alphaDrag[20]); 
  out[16] += 4.743416490252569*(facDiff[2]*f[11]+facDiff[1]*f[4]+facDiff[0]*f[2])*rdvySq4+1.369306393762915*f[17]*alphaDrag[27]+1.224744871391589*f[16]*alphaDrag[23]+1.369306393762915*(f[2]*alphaDrag[23]+f[10]*alphaDrag[21]+f[6]*alphaDrag[20])+0.6123724356957944*(alphaDrag[2]*f[16]+alphaDrag[1]*f[15]+alphaDrag[0]*f[9]); 
  out[17] += 0.3912303982179757*f[11]*alphaDrag[27]+0.6123724356957944*(f[2]*alphaDrag[27]+f[17]*alphaDrag[23])+0.5477225575051661*f[4]*alphaDrag[21]+0.6123724356957944*(f[11]*alphaDrag[20]+alphaDrag[2]*f[17])+0.3912303982179757*alphaDrag[7]*f[13]+0.6123724356957944*(alphaDrag[0]*f[13]+f[3]*alphaDrag[7])+0.5477225575051661*alphaDrag[1]*f[5]; 
  out[18] += (4.242640687119286*facDiff[1]*f[13]+4.743416490252569*(facDiff[0]*f[5]+facDiff[1]*f[3])+4.242640687119286*facDiff[2]*f[5])*rdvxSq4+0.5477225575051661*f[12]*alphaDrag[27]+0.6123724356957944*(f[18]*alphaDrag[23]+f[8]*alphaDrag[21]+f[12]*alphaDrag[20])+1.224744871391589*(alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[7]*f[10])+1.369306393762915*(alphaDrag[0]*f[10]+alphaDrag[1]*f[6]+alphaDrag[2]*f[5]); 
  out[19] += (4.242640687119286*facDiff[1]*f[11]+4.743416490252569*(facDiff[0]*f[4]+facDiff[1]*f[2])+4.242640687119286*facDiff[2]*f[4])*rdvySq4+1.224744871391589*f[10]*alphaDrag[27]+(1.224744871391589*f[19]+1.369306393762915*f[4])*alphaDrag[23]+1.224744871391589*f[17]*alphaDrag[21]+1.369306393762915*(f[6]*alphaDrag[21]+f[10]*alphaDrag[20])+0.6123724356957944*alphaDrag[2]*f[19]+0.5477225575051661*alphaDrag[7]*f[15]+0.6123724356957944*(alphaDrag[0]*f[15]+alphaDrag[1]*f[9]); 

  return fabs(0.1767766952966368*alphaDrag[0]-0.1976423537605236*alphaDrag[7])+fabs(0.1767766952966368*alphaDrag[20]-0.1976423537605236*alphaDrag[27])+fabs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4)+fabs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvySq4); 

} 
