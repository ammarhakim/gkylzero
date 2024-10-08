#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_stream_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // f:         Input distribution function.
  // out:       Incremented output.
  double w3Ddx0  = w[3]/dxv[0]; 
  double dv3Ddx0 = dxv[3]/dxv[0]; 
  double w4Ddx1  = w[4]/dxv[1]; 
  double dv4Ddx1 = dxv[4]/dxv[1]; 
  double w5Ddx2  = w[5]/dxv[2]; 
  double dv5Ddx2 = dxv[5]/dxv[2]; 

  out[1] += 3.464101615137754*f[0]*w3Ddx0+f[4]*dv3Ddx0; 
  out[2] += 3.464101615137754*f[0]*w4Ddx1+f[5]*dv4Ddx1; 
  out[3] += 3.464101615137754*f[0]*w5Ddx2+f[6]*dv5Ddx2; 
  out[7] += 3.464101615137754*f[1]*w4Ddx1+3.464101615137754*f[2]*w3Ddx0+f[13]*dv4Ddx1+f[11]*dv3Ddx0; 
  out[8] += 3.464101615137754*f[1]*w5Ddx2+3.464101615137754*f[3]*w3Ddx0+f[17]*dv5Ddx2+f[12]*dv3Ddx0; 
  out[9] += 3.464101615137754*f[2]*w5Ddx2+3.464101615137754*f[3]*w4Ddx1+f[18]*dv5Ddx2+f[15]*dv4Ddx1; 
  out[10] += 3.464101615137754*f[4]*w3Ddx0+f[0]*dv3Ddx0; 
  out[11] += 3.464101615137754*f[4]*w4Ddx1+f[16]*dv4Ddx1; 
  out[12] += 3.464101615137754*f[4]*w5Ddx2+f[20]*dv5Ddx2; 
  out[13] += 3.464101615137754*f[5]*w3Ddx0+f[16]*dv3Ddx0; 
  out[14] += 3.464101615137754*f[5]*w4Ddx1+f[0]*dv4Ddx1; 
  out[15] += 3.464101615137754*f[5]*w5Ddx2+f[21]*dv5Ddx2; 
  out[17] += 3.464101615137754*f[6]*w3Ddx0+f[20]*dv3Ddx0; 
  out[18] += 3.464101615137754*f[6]*w4Ddx1+f[21]*dv4Ddx1; 
  out[19] += 3.464101615137754*f[6]*w5Ddx2+f[0]*dv5Ddx2; 
  out[22] += 3.464101615137754*f[7]*w5Ddx2+3.464101615137754*f[8]*w4Ddx1+3.464101615137754*f[9]*w3Ddx0+f[32]*dv5Ddx2+f[27]*dv4Ddx1+f[25]*dv3Ddx0; 
  out[23] += 3.464101615137754*f[10]*w4Ddx1+3.464101615137754*f[11]*w3Ddx0+f[29]*dv4Ddx1+f[2]*dv3Ddx0; 
  out[24] += 3.464101615137754*f[10]*w5Ddx2+3.464101615137754*f[12]*w3Ddx0+f[35]*dv5Ddx2+f[3]*dv3Ddx0; 
  out[25] += 3.464101615137754*f[11]*w5Ddx2+3.464101615137754*f[12]*w4Ddx1+f[36]*dv5Ddx2+f[31]*dv4Ddx1; 
  out[26] += 3.464101615137754*f[13]*w4Ddx1+3.464101615137754*f[14]*w3Ddx0+f[1]*dv4Ddx1+f[30]*dv3Ddx0; 
  out[27] += 3.464101615137754*f[13]*w5Ddx2+3.464101615137754*f[15]*w3Ddx0+f[38]*dv5Ddx2+f[31]*dv3Ddx0; 
  out[28] += 3.464101615137754*f[14]*w5Ddx2+3.464101615137754*f[15]*w4Ddx1+f[39]*dv5Ddx2+f[3]*dv4Ddx1; 
  out[29] += 3.464101615137754*f[16]*w3Ddx0+f[5]*dv3Ddx0; 
  out[30] += 3.464101615137754*f[16]*w4Ddx1+f[4]*dv4Ddx1; 
  out[31] += 3.464101615137754*f[16]*w5Ddx2+f[41]*dv5Ddx2; 
  out[32] += 3.464101615137754*f[17]*w4Ddx1+3.464101615137754*f[18]*w3Ddx0+f[38]*dv4Ddx1+f[36]*dv3Ddx0; 
  out[33] += 3.464101615137754*f[17]*w5Ddx2+3.464101615137754*f[19]*w3Ddx0+f[1]*dv5Ddx2+f[37]*dv3Ddx0; 
  out[34] += 3.464101615137754*f[18]*w5Ddx2+3.464101615137754*f[19]*w4Ddx1+f[2]*dv5Ddx2+f[40]*dv4Ddx1; 
  out[35] += 3.464101615137754*f[20]*w3Ddx0+f[6]*dv3Ddx0; 
  out[36] += 3.464101615137754*f[20]*w4Ddx1+f[41]*dv4Ddx1; 
  out[37] += 3.464101615137754*f[20]*w5Ddx2+f[4]*dv5Ddx2; 
  out[38] += 3.464101615137754*f[21]*w3Ddx0+f[41]*dv3Ddx0; 
  out[39] += 3.464101615137754*f[21]*w4Ddx1+f[6]*dv4Ddx1; 
  out[40] += 3.464101615137754*f[21]*w5Ddx2+f[5]*dv5Ddx2; 
  out[42] += 3.464101615137754*f[23]*w5Ddx2+3.464101615137754*f[24]*w4Ddx1+3.464101615137754*f[25]*w3Ddx0+f[48]*dv5Ddx2+f[45]*dv4Ddx1+f[9]*dv3Ddx0; 
  out[43] += 3.464101615137754*f[26]*w5Ddx2+3.464101615137754*f[27]*w4Ddx1+3.464101615137754*f[28]*w3Ddx0+f[51]*dv5Ddx2+f[8]*dv4Ddx1+f[46]*dv3Ddx0; 
  out[44] += 3.464101615137754*f[29]*w4Ddx1+3.464101615137754*f[30]*w3Ddx0+f[10]*dv4Ddx1+f[14]*dv3Ddx0; 
  out[45] += 3.464101615137754*f[29]*w5Ddx2+3.464101615137754*f[31]*w3Ddx0+f[54]*dv5Ddx2+f[15]*dv3Ddx0; 
  out[46] += 3.464101615137754*f[30]*w5Ddx2+3.464101615137754*f[31]*w4Ddx1+f[55]*dv5Ddx2+f[12]*dv4Ddx1; 
  out[47] += 3.464101615137754*f[32]*w5Ddx2+3.464101615137754*f[33]*w4Ddx1+3.464101615137754*f[34]*w3Ddx0+f[7]*dv5Ddx2+f[52]*dv4Ddx1+f[50]*dv3Ddx0; 
  out[48] += 3.464101615137754*f[35]*w4Ddx1+3.464101615137754*f[36]*w3Ddx0+f[54]*dv4Ddx1+f[18]*dv3Ddx0; 
  out[49] += 3.464101615137754*f[35]*w5Ddx2+3.464101615137754*f[37]*w3Ddx0+f[10]*dv5Ddx2+f[19]*dv3Ddx0; 
  out[50] += 3.464101615137754*f[36]*w5Ddx2+3.464101615137754*f[37]*w4Ddx1+f[11]*dv5Ddx2+f[56]*dv4Ddx1; 
  out[51] += 3.464101615137754*f[38]*w4Ddx1+3.464101615137754*f[39]*w3Ddx0+f[17]*dv4Ddx1+f[55]*dv3Ddx0; 
  out[52] += 3.464101615137754*f[38]*w5Ddx2+3.464101615137754*f[40]*w3Ddx0+f[13]*dv5Ddx2+f[56]*dv3Ddx0; 
  out[53] += 3.464101615137754*f[39]*w5Ddx2+3.464101615137754*f[40]*w4Ddx1+f[14]*dv5Ddx2+f[19]*dv4Ddx1; 
  out[54] += 3.464101615137754*f[41]*w3Ddx0+f[21]*dv3Ddx0; 
  out[55] += 3.464101615137754*f[41]*w4Ddx1+f[20]*dv4Ddx1; 
  out[56] += 3.464101615137754*f[41]*w5Ddx2+f[16]*dv5Ddx2; 
  out[57] += 3.464101615137754*f[44]*w5Ddx2+3.464101615137754*f[45]*w4Ddx1+3.464101615137754*f[46]*w3Ddx0+f[60]*dv5Ddx2+f[24]*dv4Ddx1+f[28]*dv3Ddx0; 
  out[58] += 3.464101615137754*f[48]*w5Ddx2+3.464101615137754*f[49]*w4Ddx1+3.464101615137754*f[50]*w3Ddx0+f[23]*dv5Ddx2+f[61]*dv4Ddx1+f[34]*dv3Ddx0; 
  out[59] += 3.464101615137754*f[51]*w5Ddx2+3.464101615137754*f[52]*w4Ddx1+3.464101615137754*f[53]*w3Ddx0+f[26]*dv5Ddx2+f[33]*dv4Ddx1+f[62]*dv3Ddx0; 
  out[60] += 3.464101615137754*f[54]*w4Ddx1+3.464101615137754*f[55]*w3Ddx0+f[35]*dv4Ddx1+f[39]*dv3Ddx0; 
  out[61] += 3.464101615137754*f[54]*w5Ddx2+3.464101615137754*f[56]*w3Ddx0+f[29]*dv5Ddx2+f[40]*dv3Ddx0; 
  out[62] += 3.464101615137754*f[55]*w5Ddx2+3.464101615137754*f[56]*w4Ddx1+f[30]*dv5Ddx2+f[37]*dv4Ddx1; 
  out[63] += 3.464101615137754*f[60]*w5Ddx2+3.464101615137754*f[61]*w4Ddx1+3.464101615137754*f[62]*w3Ddx0+f[44]*dv5Ddx2+f[49]*dv4Ddx1+f[53]*dv3Ddx0; 

  return 3.0*(fabs(w3Ddx0)+0.5*dv3Ddx0+fabs(w4Ddx1)+0.5*dv4Ddx1+fabs(w5Ddx2)+0.5*dv5Ddx2);
} 
