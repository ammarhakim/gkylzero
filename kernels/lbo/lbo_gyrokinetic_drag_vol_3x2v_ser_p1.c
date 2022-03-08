#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH double lbo_gyrokinetic_drag_vol_3x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[5]:      cell-center coordinates. 
  // dxv[5]:    cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         input distribution function.
  // out:       incremented output 
  double rdv2[2]; 
  rdv2[0]   = 2.0/dxv[3]; 
  rdv2[1]   = 2.0/dxv[4]; 

  double alphaDrag[64]; 
  // Expand rdv2*(nu*vpar-nuUparSum) in phase basis.
  alphaDrag[0] = rdv2[0]*(2.0*nuUSum[0]-2.0*nuSum[0]*w[3]); 
  alphaDrag[1] = rdv2[0]*(2.0*nuUSum[1]-2.0*nuSum[1]*w[3]); 
  alphaDrag[2] = rdv2[0]*(2.0*nuUSum[2]-2.0*nuSum[2]*w[3]); 
  alphaDrag[3] = rdv2[0]*(2.0*nuUSum[3]-2.0*nuSum[3]*w[3]); 
  alphaDrag[4] = -0.5773502691896258*nuSum[0]*rdv2[0]*dxv[3]; 
  alphaDrag[6] = rdv2[0]*(2.0*nuUSum[4]-2.0*w[3]*nuSum[4]); 
  alphaDrag[7] = rdv2[0]*(2.0*nuUSum[5]-2.0*w[3]*nuSum[5]); 
  alphaDrag[8] = rdv2[0]*(2.0*nuUSum[6]-2.0*w[3]*nuSum[6]); 
  alphaDrag[9] = -0.5773502691896258*rdv2[0]*nuSum[1]*dxv[3]; 
  alphaDrag[10] = -0.5773502691896258*rdv2[0]*nuSum[2]*dxv[3]; 
  alphaDrag[11] = -0.5773502691896258*rdv2[0]*dxv[3]*nuSum[3]; 
  alphaDrag[16] = rdv2[0]*(2.0*nuUSum[7]-2.0*w[3]*nuSum[7]); 
  alphaDrag[17] = -0.5773502691896258*rdv2[0]*dxv[3]*nuSum[4]; 
  alphaDrag[18] = -0.5773502691896258*rdv2[0]*dxv[3]*nuSum[5]; 
  alphaDrag[19] = -0.5773502691896258*rdv2[0]*dxv[3]*nuSum[6]; 
  alphaDrag[26] = -0.5773502691896258*rdv2[0]*dxv[3]*nuSum[7]; 

  // Expand rdv2*nu*2*mu in phase basis.
  alphaDrag[32] = -4.0*nuSum[0]*rdv2[1]*w[4]; 
  alphaDrag[33] = -4.0*nuSum[1]*rdv2[1]*w[4]; 
  alphaDrag[34] = -4.0*rdv2[1]*nuSum[2]*w[4]; 
  alphaDrag[35] = -4.0*rdv2[1]*nuSum[3]*w[4]; 
  alphaDrag[37] = -1.154700538379252*nuSum[0]*rdv2[1]*dxv[4]; 
  alphaDrag[38] = -4.0*rdv2[1]*nuSum[4]*w[4]; 
  alphaDrag[39] = -4.0*rdv2[1]*w[4]*nuSum[5]; 
  alphaDrag[40] = -4.0*rdv2[1]*w[4]*nuSum[6]; 
  alphaDrag[44] = -1.154700538379252*nuSum[1]*rdv2[1]*dxv[4]; 
  alphaDrag[45] = -1.154700538379252*rdv2[1]*nuSum[2]*dxv[4]; 
  alphaDrag[46] = -1.154700538379252*rdv2[1]*nuSum[3]*dxv[4]; 
  alphaDrag[48] = -4.0*rdv2[1]*w[4]*nuSum[7]; 
  alphaDrag[52] = -1.154700538379252*rdv2[1]*dxv[4]*nuSum[4]; 
  alphaDrag[53] = -1.154700538379252*rdv2[1]*dxv[4]*nuSum[5]; 
  alphaDrag[54] = -1.154700538379252*rdv2[1]*dxv[4]*nuSum[6]; 
  alphaDrag[59] = -1.154700538379252*rdv2[1]*dxv[4]*nuSum[7]; 

  out[4] += 0.3061862178478971*(alphaDrag[26]*f[26]+alphaDrag[19]*f[19]+alphaDrag[18]*f[18]+alphaDrag[17]*f[17]+alphaDrag[16]*f[16]+alphaDrag[11]*f[11]+alphaDrag[10]*f[10]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[5] += 0.3061862178478971*(f[27]*alphaDrag[59]+f[22]*alphaDrag[54]+f[21]*alphaDrag[53]+f[20]*alphaDrag[52]+f[16]*alphaDrag[48]+f[14]*alphaDrag[46]+f[13]*alphaDrag[45]+f[12]*alphaDrag[44]+f[8]*alphaDrag[40]+f[7]*alphaDrag[39]+f[6]*alphaDrag[38]+f[5]*alphaDrag[37]+f[3]*alphaDrag[35]+f[2]*alphaDrag[34]+f[1]*alphaDrag[33]+f[0]*alphaDrag[32]); 
  out[9] += 0.3061862178478971*(alphaDrag[19]*f[26]+f[19]*alphaDrag[26]+alphaDrag[11]*f[18]+f[11]*alphaDrag[18]+alphaDrag[10]*f[17]+f[10]*alphaDrag[17]+alphaDrag[8]*f[16]+f[8]*alphaDrag[16]+alphaDrag[4]*f[9]+f[4]*alphaDrag[9]+alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[10] += 0.3061862178478971*(alphaDrag[18]*f[26]+f[18]*alphaDrag[26]+alphaDrag[11]*f[19]+f[11]*alphaDrag[19]+alphaDrag[9]*f[17]+f[9]*alphaDrag[17]+alphaDrag[7]*f[16]+f[7]*alphaDrag[16]+alphaDrag[4]*f[10]+f[4]*alphaDrag[10]+alphaDrag[3]*f[8]+f[3]*alphaDrag[8]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[11] += 0.3061862178478971*(alphaDrag[17]*f[26]+f[17]*alphaDrag[26]+alphaDrag[10]*f[19]+f[10]*alphaDrag[19]+alphaDrag[9]*f[18]+f[9]*alphaDrag[18]+alphaDrag[6]*f[16]+f[6]*alphaDrag[16]+alphaDrag[4]*f[11]+f[4]*alphaDrag[11]+alphaDrag[2]*f[8]+f[2]*alphaDrag[8]+alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[12] += 0.3061862178478971*(f[22]*alphaDrag[59]+f[27]*alphaDrag[54]+f[14]*alphaDrag[53]+f[13]*alphaDrag[52]+f[8]*alphaDrag[48]+f[21]*alphaDrag[46]+f[20]*alphaDrag[45]+f[5]*alphaDrag[44]+f[16]*alphaDrag[40]+f[3]*alphaDrag[39]+f[2]*alphaDrag[38]+f[12]*alphaDrag[37]+f[7]*alphaDrag[35]+f[6]*alphaDrag[34]+f[0]*alphaDrag[33]+f[1]*alphaDrag[32]); 
  out[13] += 0.3061862178478971*(f[21]*alphaDrag[59]+f[14]*alphaDrag[54]+f[27]*alphaDrag[53]+f[12]*alphaDrag[52]+f[7]*alphaDrag[48]+f[22]*alphaDrag[46]+f[5]*alphaDrag[45]+f[20]*alphaDrag[44]+f[3]*alphaDrag[40]+f[16]*alphaDrag[39]+f[1]*alphaDrag[38]+f[13]*alphaDrag[37]+f[8]*alphaDrag[35]+f[0]*alphaDrag[34]+f[6]*alphaDrag[33]+f[2]*alphaDrag[32]); 
  out[14] += 0.3061862178478971*(f[20]*alphaDrag[59]+f[13]*alphaDrag[54]+f[12]*alphaDrag[53]+f[27]*alphaDrag[52]+f[6]*alphaDrag[48]+f[5]*alphaDrag[46]+f[22]*alphaDrag[45]+f[21]*alphaDrag[44]+f[2]*alphaDrag[40]+f[1]*alphaDrag[39]+f[16]*alphaDrag[38]+f[14]*alphaDrag[37]+f[0]*alphaDrag[35]+f[8]*alphaDrag[34]+f[7]*alphaDrag[33]+f[3]*alphaDrag[32]); 
  out[15] += 0.3061862178478971*(f[31]*alphaDrag[59]+f[30]*alphaDrag[54]+f[29]*alphaDrag[53]+f[28]*alphaDrag[52]+f[26]*alphaDrag[48]+f[25]*alphaDrag[46]+f[24]*alphaDrag[45]+f[23]*alphaDrag[44]+f[19]*alphaDrag[40]+f[18]*alphaDrag[39]+f[17]*alphaDrag[38]+f[15]*alphaDrag[37]+f[11]*alphaDrag[35]+f[10]*alphaDrag[34]+f[9]*alphaDrag[33]+f[4]*alphaDrag[32]+alphaDrag[26]*f[31]+alphaDrag[19]*f[30]+alphaDrag[18]*f[29]+alphaDrag[17]*f[28]+alphaDrag[16]*f[27]+alphaDrag[11]*f[25]+alphaDrag[10]*f[24]+alphaDrag[9]*f[23]+alphaDrag[8]*f[22]+alphaDrag[7]*f[21]+alphaDrag[6]*f[20]+alphaDrag[4]*f[15]+alphaDrag[3]*f[14]+alphaDrag[2]*f[13]+alphaDrag[1]*f[12]+alphaDrag[0]*f[5]); 
  out[17] += 0.3061862178478971*(alphaDrag[11]*f[26]+f[11]*alphaDrag[26]+alphaDrag[18]*f[19]+f[18]*alphaDrag[19]+alphaDrag[4]*f[17]+f[4]*alphaDrag[17]+alphaDrag[3]*f[16]+f[3]*alphaDrag[16]+alphaDrag[9]*f[10]+f[9]*alphaDrag[10]+alphaDrag[7]*f[8]+f[7]*alphaDrag[8]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[18] += 0.3061862178478971*(alphaDrag[10]*f[26]+f[10]*alphaDrag[26]+alphaDrag[17]*f[19]+f[17]*alphaDrag[19]+alphaDrag[4]*f[18]+f[4]*alphaDrag[18]+alphaDrag[2]*f[16]+f[2]*alphaDrag[16]+alphaDrag[9]*f[11]+f[9]*alphaDrag[11]+alphaDrag[6]*f[8]+f[6]*alphaDrag[8]+alphaDrag[0]*f[7]+f[0]*alphaDrag[7]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[19] += 0.3061862178478971*(alphaDrag[9]*f[26]+f[9]*alphaDrag[26]+alphaDrag[4]*f[19]+f[4]*alphaDrag[19]+alphaDrag[17]*f[18]+f[17]*alphaDrag[18]+alphaDrag[1]*f[16]+f[1]*alphaDrag[16]+alphaDrag[10]*f[11]+f[10]*alphaDrag[11]+alphaDrag[0]*f[8]+f[0]*alphaDrag[8]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[20] += 0.3061862178478971*(f[14]*alphaDrag[59]+f[21]*alphaDrag[54]+f[22]*alphaDrag[53]+f[5]*alphaDrag[52]+f[3]*alphaDrag[48]+f[27]*alphaDrag[46]+f[12]*alphaDrag[45]+f[13]*alphaDrag[44]+f[7]*alphaDrag[40]+f[8]*alphaDrag[39]+f[0]*alphaDrag[38]+f[20]*alphaDrag[37]+f[16]*alphaDrag[35]+f[1]*alphaDrag[34]+f[2]*alphaDrag[33]+f[6]*alphaDrag[32]); 
  out[21] += 0.3061862178478971*(f[13]*alphaDrag[59]+f[20]*alphaDrag[54]+f[5]*alphaDrag[53]+f[22]*alphaDrag[52]+f[2]*alphaDrag[48]+f[12]*alphaDrag[46]+f[27]*alphaDrag[45]+f[14]*alphaDrag[44]+f[6]*alphaDrag[40]+f[0]*alphaDrag[39]+f[8]*alphaDrag[38]+f[21]*alphaDrag[37]+f[1]*alphaDrag[35]+f[16]*alphaDrag[34]+f[3]*alphaDrag[33]+f[7]*alphaDrag[32]); 
  out[22] += 0.3061862178478971*(f[12]*alphaDrag[59]+f[5]*alphaDrag[54]+f[20]*alphaDrag[53]+f[21]*alphaDrag[52]+f[1]*alphaDrag[48]+f[13]*alphaDrag[46]+f[14]*alphaDrag[45]+f[27]*alphaDrag[44]+f[0]*alphaDrag[40]+f[6]*alphaDrag[39]+f[7]*alphaDrag[38]+f[22]*alphaDrag[37]+f[2]*alphaDrag[35]+f[3]*alphaDrag[34]+f[16]*alphaDrag[33]+f[8]*alphaDrag[32]); 
  out[23] += 0.3061862178478971*(f[30]*alphaDrag[59]+f[31]*alphaDrag[54]+f[25]*alphaDrag[53]+f[24]*alphaDrag[52]+f[19]*alphaDrag[48]+f[29]*alphaDrag[46]+f[28]*alphaDrag[45]+f[15]*alphaDrag[44]+f[26]*alphaDrag[40]+f[11]*alphaDrag[39]+f[10]*alphaDrag[38]+f[23]*alphaDrag[37]+f[18]*alphaDrag[35]+f[17]*alphaDrag[34]+f[4]*alphaDrag[33]+f[9]*alphaDrag[32]+alphaDrag[19]*f[31]+alphaDrag[26]*f[30]+alphaDrag[11]*f[29]+alphaDrag[10]*f[28]+alphaDrag[8]*f[27]+alphaDrag[18]*f[25]+alphaDrag[17]*f[24]+alphaDrag[4]*f[23]+alphaDrag[16]*f[22]+alphaDrag[3]*f[21]+alphaDrag[2]*f[20]+alphaDrag[9]*f[15]+alphaDrag[7]*f[14]+alphaDrag[6]*f[13]+alphaDrag[0]*f[12]+alphaDrag[1]*f[5]); 
  out[24] += 0.3061862178478971*(f[29]*alphaDrag[59]+f[25]*alphaDrag[54]+f[31]*alphaDrag[53]+f[23]*alphaDrag[52]+f[18]*alphaDrag[48]+f[30]*alphaDrag[46]+f[15]*alphaDrag[45]+f[28]*alphaDrag[44]+f[11]*alphaDrag[40]+f[26]*alphaDrag[39]+f[9]*alphaDrag[38]+f[24]*alphaDrag[37]+f[19]*alphaDrag[35]+f[4]*alphaDrag[34]+f[17]*alphaDrag[33]+f[10]*alphaDrag[32]+alphaDrag[18]*f[31]+alphaDrag[11]*f[30]+alphaDrag[26]*f[29]+alphaDrag[9]*f[28]+alphaDrag[7]*f[27]+alphaDrag[19]*f[25]+alphaDrag[4]*f[24]+alphaDrag[17]*f[23]+alphaDrag[3]*f[22]+alphaDrag[16]*f[21]+alphaDrag[1]*f[20]+alphaDrag[10]*f[15]+alphaDrag[8]*f[14]+alphaDrag[0]*f[13]+alphaDrag[6]*f[12]+alphaDrag[2]*f[5]); 
  out[25] += 0.3061862178478971*(f[28]*alphaDrag[59]+f[24]*alphaDrag[54]+f[23]*alphaDrag[53]+f[31]*alphaDrag[52]+f[17]*alphaDrag[48]+f[15]*alphaDrag[46]+f[30]*alphaDrag[45]+f[29]*alphaDrag[44]+f[10]*alphaDrag[40]+f[9]*alphaDrag[39]+f[26]*alphaDrag[38]+f[25]*alphaDrag[37]+f[4]*alphaDrag[35]+f[19]*alphaDrag[34]+f[18]*alphaDrag[33]+f[11]*alphaDrag[32]+alphaDrag[17]*f[31]+alphaDrag[10]*f[30]+alphaDrag[9]*f[29]+alphaDrag[26]*f[28]+alphaDrag[6]*f[27]+alphaDrag[4]*f[25]+alphaDrag[19]*f[24]+alphaDrag[18]*f[23]+alphaDrag[2]*f[22]+alphaDrag[1]*f[21]+alphaDrag[16]*f[20]+alphaDrag[11]*f[15]+alphaDrag[0]*f[14]+alphaDrag[8]*f[13]+alphaDrag[7]*f[12]+alphaDrag[3]*f[5]); 
  out[26] += 0.3061862178478971*(alphaDrag[4]*f[26]+f[4]*alphaDrag[26]+alphaDrag[9]*f[19]+f[9]*alphaDrag[19]+alphaDrag[10]*f[18]+f[10]*alphaDrag[18]+alphaDrag[11]*f[17]+f[11]*alphaDrag[17]+alphaDrag[0]*f[16]+f[0]*alphaDrag[16]+alphaDrag[1]*f[8]+f[1]*alphaDrag[8]+alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]); 
  out[27] += 0.3061862178478971*(f[5]*alphaDrag[59]+f[12]*alphaDrag[54]+f[13]*alphaDrag[53]+f[14]*alphaDrag[52]+f[0]*alphaDrag[48]+f[20]*alphaDrag[46]+f[21]*alphaDrag[45]+f[22]*alphaDrag[44]+f[1]*alphaDrag[40]+f[2]*alphaDrag[39]+f[3]*alphaDrag[38]+f[27]*alphaDrag[37]+f[6]*alphaDrag[35]+f[7]*alphaDrag[34]+f[8]*alphaDrag[33]+f[16]*alphaDrag[32]); 
  out[28] += 0.3061862178478971*(f[25]*alphaDrag[59]+f[29]*alphaDrag[54]+f[30]*alphaDrag[53]+f[15]*alphaDrag[52]+f[11]*alphaDrag[48]+f[31]*alphaDrag[46]+f[23]*alphaDrag[45]+f[24]*alphaDrag[44]+f[18]*alphaDrag[40]+f[19]*alphaDrag[39]+f[4]*alphaDrag[38]+f[28]*alphaDrag[37]+f[26]*alphaDrag[35]+f[9]*alphaDrag[34]+f[10]*alphaDrag[33]+f[17]*alphaDrag[32]+alphaDrag[11]*f[31]+alphaDrag[18]*f[30]+alphaDrag[19]*f[29]+alphaDrag[4]*f[28]+alphaDrag[3]*f[27]+f[25]*alphaDrag[26]+alphaDrag[9]*f[24]+alphaDrag[10]*f[23]+alphaDrag[7]*f[22]+alphaDrag[8]*f[21]+alphaDrag[0]*f[20]+f[15]*alphaDrag[17]+f[14]*alphaDrag[16]+alphaDrag[1]*f[13]+alphaDrag[2]*f[12]+f[5]*alphaDrag[6]); 
  out[29] += 0.3061862178478971*(f[24]*alphaDrag[59]+f[28]*alphaDrag[54]+f[15]*alphaDrag[53]+f[30]*alphaDrag[52]+f[10]*alphaDrag[48]+f[23]*alphaDrag[46]+f[31]*alphaDrag[45]+f[25]*alphaDrag[44]+f[17]*alphaDrag[40]+f[4]*alphaDrag[39]+f[19]*alphaDrag[38]+f[29]*alphaDrag[37]+f[9]*alphaDrag[35]+f[26]*alphaDrag[34]+f[11]*alphaDrag[33]+f[18]*alphaDrag[32]+alphaDrag[10]*f[31]+alphaDrag[17]*f[30]+alphaDrag[4]*f[29]+alphaDrag[19]*f[28]+alphaDrag[2]*f[27]+f[24]*alphaDrag[26]+alphaDrag[9]*f[25]+alphaDrag[11]*f[23]+alphaDrag[6]*f[22]+alphaDrag[0]*f[21]+alphaDrag[8]*f[20]+f[15]*alphaDrag[18]+f[13]*alphaDrag[16]+alphaDrag[1]*f[14]+alphaDrag[3]*f[12]+f[5]*alphaDrag[7]); 
  out[30] += 0.3061862178478971*(f[23]*alphaDrag[59]+f[15]*alphaDrag[54]+f[28]*alphaDrag[53]+f[29]*alphaDrag[52]+f[9]*alphaDrag[48]+f[24]*alphaDrag[46]+f[25]*alphaDrag[45]+f[31]*alphaDrag[44]+f[4]*alphaDrag[40]+f[17]*alphaDrag[39]+f[18]*alphaDrag[38]+f[30]*alphaDrag[37]+f[10]*alphaDrag[35]+f[11]*alphaDrag[34]+f[26]*alphaDrag[33]+f[19]*alphaDrag[32]+alphaDrag[9]*f[31]+alphaDrag[4]*f[30]+alphaDrag[17]*f[29]+alphaDrag[18]*f[28]+alphaDrag[1]*f[27]+f[23]*alphaDrag[26]+alphaDrag[10]*f[25]+alphaDrag[11]*f[24]+alphaDrag[0]*f[22]+alphaDrag[6]*f[21]+alphaDrag[7]*f[20]+f[15]*alphaDrag[19]+f[12]*alphaDrag[16]+alphaDrag[2]*f[14]+alphaDrag[3]*f[13]+f[5]*alphaDrag[8]); 
  out[31] += 0.3061862178478971*(f[15]*alphaDrag[59]+f[23]*alphaDrag[54]+f[24]*alphaDrag[53]+f[25]*alphaDrag[52]+f[4]*alphaDrag[48]+f[28]*alphaDrag[46]+f[29]*alphaDrag[45]+f[30]*alphaDrag[44]+f[9]*alphaDrag[40]+f[10]*alphaDrag[39]+f[11]*alphaDrag[38]+f[31]*alphaDrag[37]+f[17]*alphaDrag[35]+f[18]*alphaDrag[34]+f[19]*alphaDrag[33]+f[26]*alphaDrag[32]+alphaDrag[4]*f[31]+alphaDrag[9]*f[30]+alphaDrag[10]*f[29]+alphaDrag[11]*f[28]+alphaDrag[0]*f[27]+f[15]*alphaDrag[26]+alphaDrag[17]*f[25]+alphaDrag[18]*f[24]+alphaDrag[19]*f[23]+alphaDrag[1]*f[22]+alphaDrag[2]*f[21]+alphaDrag[3]*f[20]+f[5]*alphaDrag[16]+alphaDrag[6]*f[14]+alphaDrag[7]*f[13]+alphaDrag[8]*f[12]); 

  return fabs(0.0883883476483184*alphaDrag[0])+fabs(0.0883883476483184*alphaDrag[32]); 

} 
