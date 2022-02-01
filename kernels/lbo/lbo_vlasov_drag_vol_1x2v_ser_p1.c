#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvy2 = 2.0/dxv[2]; 

  double alphaDrag[16]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.0*nuUSum[0]-2.0*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (2.0*nuUSum[1]-2.0*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -0.5773502691896258*nuSum[0]*dxv[1]*rdvx2; 
  alphaDrag[4] = -0.5773502691896258*dxv[1]*nuSum[1]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[8] = (2.0*nuUSum[2]-2.0*nuSum[0]*w[2])*rdvy2; 
  alphaDrag[9] = (2.0*nuUSum[3]-2.0*nuSum[1]*w[2])*rdvy2; 
  alphaDrag[11] = -0.5773502691896258*nuSum[0]*dxv[2]*rdvy2; 
  alphaDrag[13] = -0.5773502691896258*nuSum[1]*dxv[2]*rdvy2; 

  out[2] += 0.6123724356957944*(alphaDrag[4]*f[4]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[5]*alphaDrag[13]+f[3]*alphaDrag[11]+f[1]*alphaDrag[9]+f[0]*alphaDrag[8]); 
  out[4] += 0.6123724356957944*(alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.6123724356957944*(f[3]*alphaDrag[13]+f[5]*alphaDrag[11]+f[0]*alphaDrag[9]+f[1]*alphaDrag[8]); 
  out[6] += 0.6123724356957944*(f[7]*alphaDrag[13]+f[6]*alphaDrag[11]+f[4]*alphaDrag[9]+f[2]*alphaDrag[8]+alphaDrag[4]*f[7]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[7] += 0.6123724356957944*(f[6]*alphaDrag[13]+f[7]*alphaDrag[11]+f[2]*alphaDrag[9]+f[4]*alphaDrag[8]+alphaDrag[2]*f[7]+alphaDrag[4]*f[6]+alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 

  return fabs(0.1767766952966368*alphaDrag[0])+fabs(0.1767766952966368*alphaDrag[8]); 

} 
