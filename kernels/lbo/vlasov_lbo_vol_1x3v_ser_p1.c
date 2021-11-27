#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double vlasov_lbo_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]:      Cell-center coordinates. 
  // dxv[4]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvz2 = 2.0/dxv[3]; 
  const double rdvzSq4 = 4.0/(dxv[3]*dxv[3]); 

  double alphaDrag[48]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-2.828427124746191*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (2.828427124746191*nuUSum[1]-2.828427124746191*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -0.8164965809277261*nuSum[0]*dxv[1]*rdvx2; 
  alphaDrag[5] = -0.8164965809277261*dxv[1]*nuSum[1]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[16] = (2.828427124746191*nuUSum[2]-2.828427124746191*nuSum[0]*w[2])*rdvy2; 
  alphaDrag[17] = (2.828427124746191*nuUSum[3]-2.828427124746191*nuSum[1]*w[2])*rdvy2; 
  alphaDrag[19] = -0.8164965809277261*nuSum[0]*dxv[2]*rdvy2; 
  alphaDrag[22] = -0.8164965809277261*nuSum[1]*dxv[2]*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[32] = (2.828427124746191*nuUSum[4]-2.828427124746191*nuSum[0]*w[3])*rdvz2; 
  alphaDrag[33] = (2.828427124746191*nuUSum[5]-2.828427124746191*nuSum[1]*w[3])*rdvz2; 
  alphaDrag[36] = -0.8164965809277261*nuSum[0]*dxv[3]*rdvz2; 
  alphaDrag[40] = -0.8164965809277261*nuSum[1]*dxv[3]*rdvz2; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.4330127018922193*(alphaDrag[5]*f[5]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[6]*alphaDrag[22]+f[3]*alphaDrag[19]+f[1]*alphaDrag[17]+f[0]*alphaDrag[16]); 
  out[4] += 0.4330127018922193*(f[8]*alphaDrag[40]+f[4]*alphaDrag[36]+f[1]*alphaDrag[33]+f[0]*alphaDrag[32]); 
  out[5] += 0.4330127018922193*(alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.4330127018922193*(f[3]*alphaDrag[22]+f[6]*alphaDrag[19]+f[0]*alphaDrag[17]+f[1]*alphaDrag[16]); 
  out[7] += 0.4330127018922193*(f[11]*alphaDrag[22]+f[7]*alphaDrag[19]+f[5]*alphaDrag[17]+f[2]*alphaDrag[16]+alphaDrag[5]*f[11]+alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+alphaDrag[0]*f[3]); 
  out[8] += 0.4330127018922193*(f[4]*alphaDrag[40]+f[8]*alphaDrag[36]+f[0]*alphaDrag[33]+f[1]*alphaDrag[32]); 
  out[9] += 0.4330127018922193*(f[12]*alphaDrag[40]+f[9]*alphaDrag[36]+f[5]*alphaDrag[33]+f[2]*alphaDrag[32]+alphaDrag[5]*f[12]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[13]*alphaDrag[40]+f[10]*alphaDrag[36]+f[6]*alphaDrag[33]+f[3]*alphaDrag[32]+f[13]*alphaDrag[22]+f[10]*alphaDrag[19]+f[8]*alphaDrag[17]+f[4]*alphaDrag[16]); 
  out[11] += 0.4330127018922193*(f[7]*alphaDrag[22]+f[11]*alphaDrag[19]+f[2]*alphaDrag[17]+f[5]*alphaDrag[16]+alphaDrag[2]*f[11]+alphaDrag[5]*f[7]+alphaDrag[0]*f[6]+alphaDrag[1]*f[3]); 
  out[12] += 0.4330127018922193*(f[9]*alphaDrag[40]+f[12]*alphaDrag[36]+f[2]*alphaDrag[33]+f[5]*alphaDrag[32]+alphaDrag[2]*f[12]+alphaDrag[5]*f[9]+alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[13] += 0.4330127018922193*(f[10]*alphaDrag[40]+f[13]*alphaDrag[36]+f[3]*alphaDrag[33]+f[6]*alphaDrag[32]+f[10]*alphaDrag[22]+f[13]*alphaDrag[19]+f[4]*alphaDrag[17]+f[8]*alphaDrag[16]); 
  out[14] += 0.4330127018922193*(f[15]*alphaDrag[40]+f[14]*alphaDrag[36]+f[11]*alphaDrag[33]+f[7]*alphaDrag[32]+f[15]*alphaDrag[22]+f[14]*alphaDrag[19]+f[12]*alphaDrag[17]+f[9]*alphaDrag[16]+alphaDrag[5]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[0]*f[10]); 
  out[15] += 0.4330127018922193*(f[14]*alphaDrag[40]+f[15]*alphaDrag[36]+f[7]*alphaDrag[33]+f[11]*alphaDrag[32]+f[14]*alphaDrag[22]+f[15]*alphaDrag[19]+f[9]*alphaDrag[17]+f[12]*alphaDrag[16]+alphaDrag[2]*f[15]+alphaDrag[5]*f[14]+alphaDrag[0]*f[13]+alphaDrag[1]*f[10]); 

  return fabs(0.125*alphaDrag[0])+fabs(0.125*alphaDrag[16])+fabs(0.125*alphaDrag[32])+fabs(0.9428090415820636*nuVtSqSum[0]*rdvxSq4)+fabs(0.9428090415820636*nuVtSqSum[0]*rdvySq4)+fabs(0.9428090415820636*nuVtSqSum[0]*rdvzSq4); 

} 
