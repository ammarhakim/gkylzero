#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double vlasov_lbo_vol_2x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[5]:      Cell-center coordinates. 
  // dxv[5]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[2]; 
  const double rdvxSq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvy2 = 2.0/dxv[3]; 
  const double rdvySq4 = 4.0/(dxv[3]*dxv[3]); 
  const double rdvz2 = 2.0/dxv[4]; 
  const double rdvzSq4 = 4.0/(dxv[4]*dxv[4]); 

  double alphaDrag[96]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-5.656854249492382*nuSum[0]*w[2])*rdvx2; 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = 2.828427124746191*nuUSum[2]*rdvx2; 
  alphaDrag[3] = -1.632993161855453*nuSum[0]*dxv[2]*rdvx2; 
  alphaDrag[6] = 2.828427124746191*nuUSum[3]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[32] = (2.828427124746191*nuUSum[4]-5.656854249492382*nuSum[0]*w[3])*rdvy2; 
  alphaDrag[33] = 2.828427124746191*nuUSum[5]*rdvy2; 
  alphaDrag[34] = 2.828427124746191*nuUSum[6]*rdvy2; 
  alphaDrag[36] = -1.632993161855453*nuSum[0]*dxv[3]*rdvy2; 
  alphaDrag[38] = 2.828427124746191*nuUSum[7]*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[64] = (2.828427124746191*nuUSum[8]-5.656854249492382*nuSum[0]*w[4])*rdvz2; 
  alphaDrag[65] = 2.828427124746191*nuUSum[9]*rdvz2; 
  alphaDrag[66] = 2.828427124746191*nuUSum[10]*rdvz2; 
  alphaDrag[69] = -1.632993161855453*nuSum[0]*dxv[4]*rdvz2; 
  alphaDrag[70] = 2.828427124746191*nuUSum[11]*rdvz2; 

  // Put together updates due to drag and diffusion terms.
  out[3] += 0.3061862178478971*(alphaDrag[6]*f[6]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[4] += 0.3061862178478971*(f[6]*alphaDrag[38]+f[4]*alphaDrag[36]+f[2]*alphaDrag[34]+f[1]*alphaDrag[33]+f[0]*alphaDrag[32]); 
  out[5] += 0.3061862178478971*(f[6]*alphaDrag[70]+f[5]*alphaDrag[69]+f[2]*alphaDrag[66]+f[1]*alphaDrag[65]+f[0]*alphaDrag[64]); 
  out[7] += 0.3061862178478971*(alphaDrag[3]*f[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[8] += 0.3061862178478971*(alphaDrag[3]*f[8]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 0.3061862178478971*(f[2]*alphaDrag[38]+f[9]*alphaDrag[36]+f[6]*alphaDrag[34]+f[0]*alphaDrag[33]+f[1]*alphaDrag[32]); 
  out[10] += 0.3061862178478971*(f[1]*alphaDrag[38]+f[10]*alphaDrag[36]+f[0]*alphaDrag[34]+f[6]*alphaDrag[33]+f[2]*alphaDrag[32]); 
  out[11] += 0.3061862178478971*(f[16]*alphaDrag[38]+f[11]*alphaDrag[36]+f[8]*alphaDrag[34]+f[7]*alphaDrag[33]+f[3]*alphaDrag[32]+alphaDrag[6]*f[17]+alphaDrag[3]*f[11]+alphaDrag[2]*f[10]+alphaDrag[1]*f[9]+alphaDrag[0]*f[4]); 
  out[12] += 0.3061862178478971*(f[2]*alphaDrag[70]+f[12]*alphaDrag[69]+f[6]*alphaDrag[66]+f[0]*alphaDrag[65]+f[1]*alphaDrag[64]); 
  out[13] += 0.3061862178478971*(f[1]*alphaDrag[70]+f[13]*alphaDrag[69]+f[0]*alphaDrag[66]+f[6]*alphaDrag[65]+f[2]*alphaDrag[64]); 
  out[14] += 0.3061862178478971*(f[16]*alphaDrag[70]+f[14]*alphaDrag[69]+f[8]*alphaDrag[66]+f[7]*alphaDrag[65]+f[3]*alphaDrag[64]+alphaDrag[6]*f[20]+alphaDrag[3]*f[14]+alphaDrag[2]*f[13]+alphaDrag[1]*f[12]+alphaDrag[0]*f[5]); 
  out[15] += 0.3061862178478971*(f[17]*alphaDrag[70]+f[15]*alphaDrag[69]+f[10]*alphaDrag[66]+f[9]*alphaDrag[65]+f[4]*alphaDrag[64]+f[20]*alphaDrag[38]+f[15]*alphaDrag[36]+f[13]*alphaDrag[34]+f[12]*alphaDrag[33]+f[5]*alphaDrag[32]); 
  out[16] += 0.3061862178478971*(alphaDrag[3]*f[16]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[17] += 0.3061862178478971*(f[0]*alphaDrag[38]+f[17]*alphaDrag[36]+f[1]*alphaDrag[34]+f[2]*alphaDrag[33]+f[6]*alphaDrag[32]); 
  out[18] += 0.3061862178478971*(f[8]*alphaDrag[38]+f[18]*alphaDrag[36]+f[16]*alphaDrag[34]+f[3]*alphaDrag[33]+f[7]*alphaDrag[32]+alphaDrag[3]*f[18]+alphaDrag[2]*f[17]+alphaDrag[6]*f[10]+alphaDrag[0]*f[9]+alphaDrag[1]*f[4]); 
  out[19] += 0.3061862178478971*(f[7]*alphaDrag[38]+f[19]*alphaDrag[36]+f[3]*alphaDrag[34]+f[16]*alphaDrag[33]+f[8]*alphaDrag[32]+alphaDrag[3]*f[19]+alphaDrag[1]*f[17]+alphaDrag[0]*f[10]+alphaDrag[6]*f[9]+alphaDrag[2]*f[4]); 
  out[20] += 0.3061862178478971*(f[0]*alphaDrag[70]+f[20]*alphaDrag[69]+f[1]*alphaDrag[66]+f[2]*alphaDrag[65]+f[6]*alphaDrag[64]); 
  out[21] += 0.3061862178478971*(f[8]*alphaDrag[70]+f[21]*alphaDrag[69]+f[16]*alphaDrag[66]+f[3]*alphaDrag[65]+f[7]*alphaDrag[64]+alphaDrag[3]*f[21]+alphaDrag[2]*f[20]+alphaDrag[6]*f[13]+alphaDrag[0]*f[12]+alphaDrag[1]*f[5]); 
  out[22] += 0.3061862178478971*(f[7]*alphaDrag[70]+f[22]*alphaDrag[69]+f[3]*alphaDrag[66]+f[16]*alphaDrag[65]+f[8]*alphaDrag[64]+alphaDrag[3]*f[22]+alphaDrag[1]*f[20]+alphaDrag[0]*f[13]+alphaDrag[6]*f[12]+alphaDrag[2]*f[5]); 
  out[23] += 0.3061862178478971*(f[10]*alphaDrag[70]+f[23]*alphaDrag[69]+f[17]*alphaDrag[66]+f[4]*alphaDrag[65]+f[9]*alphaDrag[64]+f[13]*alphaDrag[38]+f[23]*alphaDrag[36]+f[20]*alphaDrag[34]+f[5]*alphaDrag[33]+f[12]*alphaDrag[32]); 
  out[24] += 0.3061862178478971*(f[9]*alphaDrag[70]+f[24]*alphaDrag[69]+f[4]*alphaDrag[66]+f[17]*alphaDrag[65]+f[10]*alphaDrag[64]+f[12]*alphaDrag[38]+f[24]*alphaDrag[36]+f[5]*alphaDrag[34]+f[20]*alphaDrag[33]+f[13]*alphaDrag[32]); 
  out[25] += 0.3061862178478971*(f[26]*alphaDrag[70]+f[25]*alphaDrag[69]+f[19]*alphaDrag[66]+f[18]*alphaDrag[65]+f[11]*alphaDrag[64]+f[27]*alphaDrag[38]+f[25]*alphaDrag[36]+f[22]*alphaDrag[34]+f[21]*alphaDrag[33]+f[14]*alphaDrag[32]+alphaDrag[6]*f[28]+alphaDrag[3]*f[25]+alphaDrag[2]*f[24]+alphaDrag[1]*f[23]+alphaDrag[0]*f[15]); 
  out[26] += 0.3061862178478971*(f[3]*alphaDrag[38]+f[26]*alphaDrag[36]+f[7]*alphaDrag[34]+f[8]*alphaDrag[33]+f[16]*alphaDrag[32]+alphaDrag[3]*f[26]+alphaDrag[0]*f[17]+alphaDrag[1]*f[10]+alphaDrag[2]*f[9]+f[4]*alphaDrag[6]); 
  out[27] += 0.3061862178478971*(f[3]*alphaDrag[70]+f[27]*alphaDrag[69]+f[7]*alphaDrag[66]+f[8]*alphaDrag[65]+f[16]*alphaDrag[64]+alphaDrag[3]*f[27]+alphaDrag[0]*f[20]+alphaDrag[1]*f[13]+alphaDrag[2]*f[12]+f[5]*alphaDrag[6]); 
  out[28] += 0.3061862178478971*(f[4]*alphaDrag[70]+f[28]*alphaDrag[69]+f[9]*alphaDrag[66]+f[10]*alphaDrag[65]+f[17]*alphaDrag[64]+f[5]*alphaDrag[38]+f[28]*alphaDrag[36]+f[12]*alphaDrag[34]+f[13]*alphaDrag[33]+f[20]*alphaDrag[32]); 
  out[29] += 0.3061862178478971*(f[19]*alphaDrag[70]+f[29]*alphaDrag[69]+f[26]*alphaDrag[66]+f[11]*alphaDrag[65]+f[18]*alphaDrag[64]+f[22]*alphaDrag[38]+f[29]*alphaDrag[36]+f[27]*alphaDrag[34]+f[14]*alphaDrag[33]+f[21]*alphaDrag[32]+alphaDrag[3]*f[29]+alphaDrag[2]*f[28]+alphaDrag[6]*f[24]+alphaDrag[0]*f[23]+alphaDrag[1]*f[15]); 
  out[30] += 0.3061862178478971*(f[18]*alphaDrag[70]+f[30]*alphaDrag[69]+f[11]*alphaDrag[66]+f[26]*alphaDrag[65]+f[19]*alphaDrag[64]+f[21]*alphaDrag[38]+f[30]*alphaDrag[36]+f[14]*alphaDrag[34]+f[27]*alphaDrag[33]+f[22]*alphaDrag[32]+alphaDrag[3]*f[30]+alphaDrag[1]*f[28]+alphaDrag[0]*f[24]+alphaDrag[6]*f[23]+alphaDrag[2]*f[15]); 
  out[31] += 0.3061862178478971*(f[11]*alphaDrag[70]+f[31]*alphaDrag[69]+f[18]*alphaDrag[66]+f[19]*alphaDrag[65]+f[26]*alphaDrag[64]+f[14]*alphaDrag[38]+f[31]*alphaDrag[36]+f[21]*alphaDrag[34]+f[22]*alphaDrag[33]+f[27]*alphaDrag[32]+alphaDrag[3]*f[31]+alphaDrag[0]*f[28]+alphaDrag[1]*f[24]+alphaDrag[2]*f[23]+alphaDrag[6]*f[15]); 

  return fabs(0.0883883476483184*alphaDrag[0])+fabs(0.0883883476483184*alphaDrag[32])+fabs(0.0883883476483184*alphaDrag[64])+fabs(0.6666666666666666*nuVtSqSum[0]*rdvxSq4)+fabs(0.6666666666666666*nuVtSqSum[0]*rdvySq4)+fabs(0.6666666666666666*nuVtSqSum[0]*rdvzSq4); 

} 
