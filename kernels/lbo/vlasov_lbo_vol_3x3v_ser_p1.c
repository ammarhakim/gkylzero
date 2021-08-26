#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double vlasov_lbo_vol_3x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double GKYL_RESTRICT *out) 
{ 
  // w[6]:      Cell-center coordinates. 
  // dxv[6]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[3]; 
  const double rdvxSq4 = 4.0/(dxv[3]*dxv[3]); 
  const double rdvy2 = 2.0/dxv[4]; 
  const double rdvySq4 = 4.0/(dxv[4]*dxv[4]); 
  const double rdvz2 = 2.0/dxv[5]; 
  const double rdvzSq4 = 4.0/(dxv[5]*dxv[5]); 

  double alphaDrag[192]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-8.0*w[3]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = 2.828427124746191*nuUSum[2]*rdvx2; 
  alphaDrag[3] = 2.828427124746191*nuUSum[3]*rdvx2; 
  alphaDrag[4] = -2.309401076758503*dxv[3]*nuSum*rdvx2; 
  alphaDrag[7] = 2.828427124746191*nuUSum[4]*rdvx2; 
  alphaDrag[8] = 2.828427124746191*nuUSum[5]*rdvx2; 
  alphaDrag[9] = 2.828427124746191*nuUSum[6]*rdvx2; 
  alphaDrag[22] = 2.828427124746191*nuUSum[7]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[64] = (2.828427124746191*nuUSum[8]-8.0*w[4]*nuSum)*rdvy2; 
  alphaDrag[65] = 2.828427124746191*nuUSum[9]*rdvy2; 
  alphaDrag[66] = 2.828427124746191*nuUSum[10]*rdvy2; 
  alphaDrag[67] = 2.828427124746191*nuUSum[11]*rdvy2; 
  alphaDrag[69] = -2.309401076758503*dxv[4]*nuSum*rdvy2; 
  alphaDrag[71] = 2.828427124746191*nuUSum[12]*rdvy2; 
  alphaDrag[72] = 2.828427124746191*nuUSum[13]*rdvy2; 
  alphaDrag[73] = 2.828427124746191*nuUSum[14]*rdvy2; 
  alphaDrag[86] = 2.828427124746191*nuUSum[15]*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[128] = (2.828427124746191*nuUSum[16]-8.0*w[5]*nuSum)*rdvz2; 
  alphaDrag[129] = 2.828427124746191*nuUSum[17]*rdvz2; 
  alphaDrag[130] = 2.828427124746191*nuUSum[18]*rdvz2; 
  alphaDrag[131] = 2.828427124746191*nuUSum[19]*rdvz2; 
  alphaDrag[134] = -2.309401076758503*dxv[5]*nuSum*rdvz2; 
  alphaDrag[135] = 2.828427124746191*nuUSum[20]*rdvz2; 
  alphaDrag[136] = 2.828427124746191*nuUSum[21]*rdvz2; 
  alphaDrag[137] = 2.828427124746191*nuUSum[22]*rdvz2; 
  alphaDrag[150] = 2.828427124746191*nuUSum[23]*rdvz2; 

  // Put together updates due to drag and diffusion terms.
  out[4] += 0.2165063509461096*(alphaDrag[22]*f[22]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[5] += 0.2165063509461096*(f[22]*alphaDrag[86]+f[9]*alphaDrag[73]+f[8]*alphaDrag[72]+f[7]*alphaDrag[71]+f[5]*alphaDrag[69]+f[3]*alphaDrag[67]+f[2]*alphaDrag[66]+f[1]*alphaDrag[65]+f[0]*alphaDrag[64]); 
  out[6] += 0.2165063509461096*(f[22]*alphaDrag[150]+f[9]*alphaDrag[137]+f[8]*alphaDrag[136]+f[7]*alphaDrag[135]+f[6]*alphaDrag[134]+f[3]*alphaDrag[131]+f[2]*alphaDrag[130]+f[1]*alphaDrag[129]+f[0]*alphaDrag[128]); 
  out[10] += 0.2165063509461096*(alphaDrag[9]*f[22]+f[9]*alphaDrag[22]+alphaDrag[4]*f[10]+alphaDrag[3]*f[8]+f[3]*alphaDrag[8]+alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[11] += 0.2165063509461096*(alphaDrag[8]*f[22]+f[8]*alphaDrag[22]+alphaDrag[4]*f[11]+alphaDrag[3]*f[9]+f[3]*alphaDrag[9]+alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[12] += 0.2165063509461096*(alphaDrag[7]*f[22]+f[7]*alphaDrag[22]+alphaDrag[4]*f[12]+alphaDrag[2]*f[9]+f[2]*alphaDrag[9]+alphaDrag[1]*f[8]+f[1]*alphaDrag[8]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[13] += 0.2165063509461096*(f[9]*alphaDrag[86]+f[22]*alphaDrag[73]+f[3]*alphaDrag[72]+f[2]*alphaDrag[71]+f[13]*alphaDrag[69]+f[8]*alphaDrag[67]+f[7]*alphaDrag[66]+f[0]*alphaDrag[65]+f[1]*alphaDrag[64]); 
  out[14] += 0.2165063509461096*(f[8]*alphaDrag[86]+f[3]*alphaDrag[73]+f[22]*alphaDrag[72]+f[1]*alphaDrag[71]+f[14]*alphaDrag[69]+f[9]*alphaDrag[67]+f[0]*alphaDrag[66]+f[7]*alphaDrag[65]+f[2]*alphaDrag[64]); 
  out[15] += 0.2165063509461096*(f[7]*alphaDrag[86]+f[2]*alphaDrag[73]+f[1]*alphaDrag[72]+f[22]*alphaDrag[71]+f[15]*alphaDrag[69]+f[0]*alphaDrag[67]+f[9]*alphaDrag[66]+f[8]*alphaDrag[65]+f[3]*alphaDrag[64]); 
  out[16] += 0.2165063509461096*(f[42]*alphaDrag[86]+f[25]*alphaDrag[73]+f[24]*alphaDrag[72]+f[23]*alphaDrag[71]+f[16]*alphaDrag[69]+f[12]*alphaDrag[67]+f[11]*alphaDrag[66]+f[10]*alphaDrag[65]+f[4]*alphaDrag[64]+alphaDrag[22]*f[43]+alphaDrag[9]*f[28]+alphaDrag[8]*f[27]+alphaDrag[7]*f[26]+alphaDrag[4]*f[16]+alphaDrag[3]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[0]*f[5]); 
  out[17] += 0.2165063509461096*(f[9]*alphaDrag[150]+f[22]*alphaDrag[137]+f[3]*alphaDrag[136]+f[2]*alphaDrag[135]+f[17]*alphaDrag[134]+f[8]*alphaDrag[131]+f[7]*alphaDrag[130]+f[0]*alphaDrag[129]+f[1]*alphaDrag[128]); 
  out[18] += 0.2165063509461096*(f[8]*alphaDrag[150]+f[3]*alphaDrag[137]+f[22]*alphaDrag[136]+f[1]*alphaDrag[135]+f[18]*alphaDrag[134]+f[9]*alphaDrag[131]+f[0]*alphaDrag[130]+f[7]*alphaDrag[129]+f[2]*alphaDrag[128]); 
  out[19] += 0.2165063509461096*(f[7]*alphaDrag[150]+f[2]*alphaDrag[137]+f[1]*alphaDrag[136]+f[22]*alphaDrag[135]+f[19]*alphaDrag[134]+f[0]*alphaDrag[131]+f[9]*alphaDrag[130]+f[8]*alphaDrag[129]+f[3]*alphaDrag[128]); 
  out[20] += 0.2165063509461096*(f[42]*alphaDrag[150]+f[25]*alphaDrag[137]+f[24]*alphaDrag[136]+f[23]*alphaDrag[135]+f[20]*alphaDrag[134]+f[12]*alphaDrag[131]+f[11]*alphaDrag[130]+f[10]*alphaDrag[129]+f[4]*alphaDrag[128]+alphaDrag[22]*f[47]+alphaDrag[9]*f[34]+alphaDrag[8]*f[33]+alphaDrag[7]*f[32]+alphaDrag[4]*f[20]+alphaDrag[3]*f[19]+alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[0]*f[6]); 
  out[21] += 0.2165063509461096*(f[43]*alphaDrag[150]+f[28]*alphaDrag[137]+f[27]*alphaDrag[136]+f[26]*alphaDrag[135]+f[21]*alphaDrag[134]+f[15]*alphaDrag[131]+f[14]*alphaDrag[130]+f[13]*alphaDrag[129]+f[5]*alphaDrag[128]+f[47]*alphaDrag[86]+f[34]*alphaDrag[73]+f[33]*alphaDrag[72]+f[32]*alphaDrag[71]+f[21]*alphaDrag[69]+f[19]*alphaDrag[67]+f[18]*alphaDrag[66]+f[17]*alphaDrag[65]+f[6]*alphaDrag[64]); 
  out[23] += 0.2165063509461096*(alphaDrag[4]*f[23]+alphaDrag[3]*f[22]+f[3]*alphaDrag[22]+alphaDrag[8]*f[9]+f[8]*alphaDrag[9]+alphaDrag[0]*f[7]+f[0]*alphaDrag[7]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[24] += 0.2165063509461096*(alphaDrag[4]*f[24]+alphaDrag[2]*f[22]+f[2]*alphaDrag[22]+alphaDrag[7]*f[9]+f[7]*alphaDrag[9]+alphaDrag[0]*f[8]+f[0]*alphaDrag[8]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[25] += 0.2165063509461096*(alphaDrag[4]*f[25]+alphaDrag[1]*f[22]+f[1]*alphaDrag[22]+alphaDrag[0]*f[9]+f[0]*alphaDrag[9]+alphaDrag[7]*f[8]+f[7]*alphaDrag[8]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[26] += 0.2165063509461096*(f[3]*alphaDrag[86]+f[8]*alphaDrag[73]+f[9]*alphaDrag[72]+f[0]*alphaDrag[71]+f[26]*alphaDrag[69]+f[22]*alphaDrag[67]+f[1]*alphaDrag[66]+f[2]*alphaDrag[65]+f[7]*alphaDrag[64]); 
  out[27] += 0.2165063509461096*(f[2]*alphaDrag[86]+f[7]*alphaDrag[73]+f[0]*alphaDrag[72]+f[9]*alphaDrag[71]+f[27]*alphaDrag[69]+f[1]*alphaDrag[67]+f[22]*alphaDrag[66]+f[3]*alphaDrag[65]+f[8]*alphaDrag[64]); 
  out[28] += 0.2165063509461096*(f[1]*alphaDrag[86]+f[0]*alphaDrag[73]+f[7]*alphaDrag[72]+f[8]*alphaDrag[71]+f[28]*alphaDrag[69]+f[2]*alphaDrag[67]+f[3]*alphaDrag[66]+f[22]*alphaDrag[65]+f[9]*alphaDrag[64]); 
  out[29] += 0.2165063509461096*(f[25]*alphaDrag[86]+f[42]*alphaDrag[73]+f[12]*alphaDrag[72]+f[11]*alphaDrag[71]+f[29]*alphaDrag[69]+f[24]*alphaDrag[67]+f[23]*alphaDrag[66]+f[4]*alphaDrag[65]+f[10]*alphaDrag[64]+alphaDrag[9]*f[43]+alphaDrag[4]*f[29]+alphaDrag[22]*f[28]+alphaDrag[3]*f[27]+alphaDrag[2]*f[26]+alphaDrag[8]*f[15]+alphaDrag[7]*f[14]+alphaDrag[0]*f[13]+alphaDrag[1]*f[5]); 
  out[30] += 0.2165063509461096*(f[24]*alphaDrag[86]+f[12]*alphaDrag[73]+f[42]*alphaDrag[72]+f[10]*alphaDrag[71]+f[30]*alphaDrag[69]+f[25]*alphaDrag[67]+f[4]*alphaDrag[66]+f[23]*alphaDrag[65]+f[11]*alphaDrag[64]+alphaDrag[8]*f[43]+alphaDrag[4]*f[30]+alphaDrag[3]*f[28]+alphaDrag[22]*f[27]+alphaDrag[1]*f[26]+alphaDrag[9]*f[15]+alphaDrag[0]*f[14]+alphaDrag[7]*f[13]+alphaDrag[2]*f[5]); 
  out[31] += 0.2165063509461096*(f[23]*alphaDrag[86]+f[11]*alphaDrag[73]+f[10]*alphaDrag[72]+f[42]*alphaDrag[71]+f[31]*alphaDrag[69]+f[4]*alphaDrag[67]+f[25]*alphaDrag[66]+f[24]*alphaDrag[65]+f[12]*alphaDrag[64]+alphaDrag[7]*f[43]+alphaDrag[4]*f[31]+alphaDrag[2]*f[28]+alphaDrag[1]*f[27]+alphaDrag[22]*f[26]+alphaDrag[0]*f[15]+alphaDrag[9]*f[14]+alphaDrag[8]*f[13]+alphaDrag[3]*f[5]); 
  out[32] += 0.2165063509461096*(f[3]*alphaDrag[150]+f[8]*alphaDrag[137]+f[9]*alphaDrag[136]+f[0]*alphaDrag[135]+f[32]*alphaDrag[134]+f[22]*alphaDrag[131]+f[1]*alphaDrag[130]+f[2]*alphaDrag[129]+f[7]*alphaDrag[128]); 
  out[33] += 0.2165063509461096*(f[2]*alphaDrag[150]+f[7]*alphaDrag[137]+f[0]*alphaDrag[136]+f[9]*alphaDrag[135]+f[33]*alphaDrag[134]+f[1]*alphaDrag[131]+f[22]*alphaDrag[130]+f[3]*alphaDrag[129]+f[8]*alphaDrag[128]); 
  out[34] += 0.2165063509461096*(f[1]*alphaDrag[150]+f[0]*alphaDrag[137]+f[7]*alphaDrag[136]+f[8]*alphaDrag[135]+f[34]*alphaDrag[134]+f[2]*alphaDrag[131]+f[3]*alphaDrag[130]+f[22]*alphaDrag[129]+f[9]*alphaDrag[128]); 
  out[35] += 0.2165063509461096*(f[25]*alphaDrag[150]+f[42]*alphaDrag[137]+f[12]*alphaDrag[136]+f[11]*alphaDrag[135]+f[35]*alphaDrag[134]+f[24]*alphaDrag[131]+f[23]*alphaDrag[130]+f[4]*alphaDrag[129]+f[10]*alphaDrag[128]+alphaDrag[9]*f[47]+alphaDrag[4]*f[35]+alphaDrag[22]*f[34]+alphaDrag[3]*f[33]+alphaDrag[2]*f[32]+alphaDrag[8]*f[19]+alphaDrag[7]*f[18]+alphaDrag[0]*f[17]+alphaDrag[1]*f[6]); 
  out[36] += 0.2165063509461096*(f[24]*alphaDrag[150]+f[12]*alphaDrag[137]+f[42]*alphaDrag[136]+f[10]*alphaDrag[135]+f[36]*alphaDrag[134]+f[25]*alphaDrag[131]+f[4]*alphaDrag[130]+f[23]*alphaDrag[129]+f[11]*alphaDrag[128]+alphaDrag[8]*f[47]+alphaDrag[4]*f[36]+alphaDrag[3]*f[34]+alphaDrag[22]*f[33]+alphaDrag[1]*f[32]+alphaDrag[9]*f[19]+alphaDrag[0]*f[18]+alphaDrag[7]*f[17]+alphaDrag[2]*f[6]); 
  out[37] += 0.2165063509461096*(f[23]*alphaDrag[150]+f[11]*alphaDrag[137]+f[10]*alphaDrag[136]+f[42]*alphaDrag[135]+f[37]*alphaDrag[134]+f[4]*alphaDrag[131]+f[25]*alphaDrag[130]+f[24]*alphaDrag[129]+f[12]*alphaDrag[128]+alphaDrag[7]*f[47]+alphaDrag[4]*f[37]+alphaDrag[2]*f[34]+alphaDrag[1]*f[33]+alphaDrag[22]*f[32]+alphaDrag[0]*f[19]+alphaDrag[9]*f[18]+alphaDrag[8]*f[17]+alphaDrag[3]*f[6]); 
  out[38] += 0.2165063509461096*(f[28]*alphaDrag[150]+f[43]*alphaDrag[137]+f[15]*alphaDrag[136]+f[14]*alphaDrag[135]+f[38]*alphaDrag[134]+f[27]*alphaDrag[131]+f[26]*alphaDrag[130]+f[5]*alphaDrag[129]+f[13]*alphaDrag[128]+f[34]*alphaDrag[86]+f[47]*alphaDrag[73]+f[19]*alphaDrag[72]+f[18]*alphaDrag[71]+f[38]*alphaDrag[69]+f[33]*alphaDrag[67]+f[32]*alphaDrag[66]+f[6]*alphaDrag[65]+f[17]*alphaDrag[64]); 
  out[39] += 0.2165063509461096*(f[27]*alphaDrag[150]+f[15]*alphaDrag[137]+f[43]*alphaDrag[136]+f[13]*alphaDrag[135]+f[39]*alphaDrag[134]+f[28]*alphaDrag[131]+f[5]*alphaDrag[130]+f[26]*alphaDrag[129]+f[14]*alphaDrag[128]+f[33]*alphaDrag[86]+f[19]*alphaDrag[73]+f[47]*alphaDrag[72]+f[17]*alphaDrag[71]+f[39]*alphaDrag[69]+f[34]*alphaDrag[67]+f[6]*alphaDrag[66]+f[32]*alphaDrag[65]+f[18]*alphaDrag[64]); 
  out[40] += 0.2165063509461096*(f[26]*alphaDrag[150]+f[14]*alphaDrag[137]+f[13]*alphaDrag[136]+f[43]*alphaDrag[135]+f[40]*alphaDrag[134]+f[5]*alphaDrag[131]+f[28]*alphaDrag[130]+f[27]*alphaDrag[129]+f[15]*alphaDrag[128]+f[32]*alphaDrag[86]+f[18]*alphaDrag[73]+f[17]*alphaDrag[72]+f[47]*alphaDrag[71]+f[40]*alphaDrag[69]+f[6]*alphaDrag[67]+f[34]*alphaDrag[66]+f[33]*alphaDrag[65]+f[19]*alphaDrag[64]); 
  out[41] += 0.2165063509461096*(f[57]*alphaDrag[150]+f[46]*alphaDrag[137]+f[45]*alphaDrag[136]+f[44]*alphaDrag[135]+f[41]*alphaDrag[134]+f[31]*alphaDrag[131]+f[30]*alphaDrag[130]+f[29]*alphaDrag[129]+f[16]*alphaDrag[128]+f[58]*alphaDrag[86]+f[50]*alphaDrag[73]+f[49]*alphaDrag[72]+f[48]*alphaDrag[71]+f[41]*alphaDrag[69]+f[37]*alphaDrag[67]+f[36]*alphaDrag[66]+f[35]*alphaDrag[65]+f[20]*alphaDrag[64]+alphaDrag[22]*f[59]+alphaDrag[9]*f[53]+alphaDrag[8]*f[52]+alphaDrag[7]*f[51]+alphaDrag[4]*f[41]+alphaDrag[3]*f[40]+alphaDrag[2]*f[39]+alphaDrag[1]*f[38]+alphaDrag[0]*f[21]); 
  out[42] += 0.2165063509461096*(alphaDrag[4]*f[42]+alphaDrag[0]*f[22]+f[0]*alphaDrag[22]+alphaDrag[1]*f[9]+f[1]*alphaDrag[9]+alphaDrag[2]*f[8]+f[2]*alphaDrag[8]+alphaDrag[3]*f[7]+f[3]*alphaDrag[7]); 
  out[43] += 0.2165063509461096*(f[0]*alphaDrag[86]+f[1]*alphaDrag[73]+f[2]*alphaDrag[72]+f[3]*alphaDrag[71]+f[43]*alphaDrag[69]+f[7]*alphaDrag[67]+f[8]*alphaDrag[66]+f[9]*alphaDrag[65]+f[22]*alphaDrag[64]); 
  out[44] += 0.2165063509461096*(f[12]*alphaDrag[86]+f[24]*alphaDrag[73]+f[25]*alphaDrag[72]+f[4]*alphaDrag[71]+f[44]*alphaDrag[69]+f[42]*alphaDrag[67]+f[10]*alphaDrag[66]+f[11]*alphaDrag[65]+f[23]*alphaDrag[64]+alphaDrag[4]*f[44]+alphaDrag[3]*f[43]+alphaDrag[8]*f[28]+alphaDrag[9]*f[27]+alphaDrag[0]*f[26]+f[15]*alphaDrag[22]+alphaDrag[1]*f[14]+alphaDrag[2]*f[13]+f[5]*alphaDrag[7]); 
  out[45] += 0.2165063509461096*(f[11]*alphaDrag[86]+f[23]*alphaDrag[73]+f[4]*alphaDrag[72]+f[25]*alphaDrag[71]+f[45]*alphaDrag[69]+f[10]*alphaDrag[67]+f[42]*alphaDrag[66]+f[12]*alphaDrag[65]+f[24]*alphaDrag[64]+alphaDrag[4]*f[45]+alphaDrag[2]*f[43]+alphaDrag[7]*f[28]+alphaDrag[0]*f[27]+alphaDrag[9]*f[26]+f[14]*alphaDrag[22]+alphaDrag[1]*f[15]+alphaDrag[3]*f[13]+f[5]*alphaDrag[8]); 
  out[46] += 0.2165063509461096*(f[10]*alphaDrag[86]+f[4]*alphaDrag[73]+f[23]*alphaDrag[72]+f[24]*alphaDrag[71]+f[46]*alphaDrag[69]+f[11]*alphaDrag[67]+f[12]*alphaDrag[66]+f[42]*alphaDrag[65]+f[25]*alphaDrag[64]+alphaDrag[4]*f[46]+alphaDrag[1]*f[43]+alphaDrag[0]*f[28]+alphaDrag[7]*f[27]+alphaDrag[8]*f[26]+f[13]*alphaDrag[22]+alphaDrag[2]*f[15]+alphaDrag[3]*f[14]+f[5]*alphaDrag[9]); 
  out[47] += 0.2165063509461096*(f[0]*alphaDrag[150]+f[1]*alphaDrag[137]+f[2]*alphaDrag[136]+f[3]*alphaDrag[135]+f[47]*alphaDrag[134]+f[7]*alphaDrag[131]+f[8]*alphaDrag[130]+f[9]*alphaDrag[129]+f[22]*alphaDrag[128]); 
  out[48] += 0.2165063509461096*(f[12]*alphaDrag[150]+f[24]*alphaDrag[137]+f[25]*alphaDrag[136]+f[4]*alphaDrag[135]+f[48]*alphaDrag[134]+f[42]*alphaDrag[131]+f[10]*alphaDrag[130]+f[11]*alphaDrag[129]+f[23]*alphaDrag[128]+alphaDrag[4]*f[48]+alphaDrag[3]*f[47]+alphaDrag[8]*f[34]+alphaDrag[9]*f[33]+alphaDrag[0]*f[32]+f[19]*alphaDrag[22]+alphaDrag[1]*f[18]+alphaDrag[2]*f[17]+f[6]*alphaDrag[7]); 
  out[49] += 0.2165063509461096*(f[11]*alphaDrag[150]+f[23]*alphaDrag[137]+f[4]*alphaDrag[136]+f[25]*alphaDrag[135]+f[49]*alphaDrag[134]+f[10]*alphaDrag[131]+f[42]*alphaDrag[130]+f[12]*alphaDrag[129]+f[24]*alphaDrag[128]+alphaDrag[4]*f[49]+alphaDrag[2]*f[47]+alphaDrag[7]*f[34]+alphaDrag[0]*f[33]+alphaDrag[9]*f[32]+f[18]*alphaDrag[22]+alphaDrag[1]*f[19]+alphaDrag[3]*f[17]+f[6]*alphaDrag[8]); 
  out[50] += 0.2165063509461096*(f[10]*alphaDrag[150]+f[4]*alphaDrag[137]+f[23]*alphaDrag[136]+f[24]*alphaDrag[135]+f[50]*alphaDrag[134]+f[11]*alphaDrag[131]+f[12]*alphaDrag[130]+f[42]*alphaDrag[129]+f[25]*alphaDrag[128]+alphaDrag[4]*f[50]+alphaDrag[1]*f[47]+alphaDrag[0]*f[34]+alphaDrag[7]*f[33]+alphaDrag[8]*f[32]+f[17]*alphaDrag[22]+alphaDrag[2]*f[19]+alphaDrag[3]*f[18]+f[6]*alphaDrag[9]); 
  out[51] += 0.2165063509461096*(f[15]*alphaDrag[150]+f[27]*alphaDrag[137]+f[28]*alphaDrag[136]+f[5]*alphaDrag[135]+f[51]*alphaDrag[134]+f[43]*alphaDrag[131]+f[13]*alphaDrag[130]+f[14]*alphaDrag[129]+f[26]*alphaDrag[128]+f[19]*alphaDrag[86]+f[33]*alphaDrag[73]+f[34]*alphaDrag[72]+f[6]*alphaDrag[71]+f[51]*alphaDrag[69]+f[47]*alphaDrag[67]+f[17]*alphaDrag[66]+f[18]*alphaDrag[65]+f[32]*alphaDrag[64]); 
  out[52] += 0.2165063509461096*(f[14]*alphaDrag[150]+f[26]*alphaDrag[137]+f[5]*alphaDrag[136]+f[28]*alphaDrag[135]+f[52]*alphaDrag[134]+f[13]*alphaDrag[131]+f[43]*alphaDrag[130]+f[15]*alphaDrag[129]+f[27]*alphaDrag[128]+f[18]*alphaDrag[86]+f[32]*alphaDrag[73]+f[6]*alphaDrag[72]+f[34]*alphaDrag[71]+f[52]*alphaDrag[69]+f[17]*alphaDrag[67]+f[47]*alphaDrag[66]+f[19]*alphaDrag[65]+f[33]*alphaDrag[64]); 
  out[53] += 0.2165063509461096*(f[13]*alphaDrag[150]+f[5]*alphaDrag[137]+f[26]*alphaDrag[136]+f[27]*alphaDrag[135]+f[53]*alphaDrag[134]+f[14]*alphaDrag[131]+f[15]*alphaDrag[130]+f[43]*alphaDrag[129]+f[28]*alphaDrag[128]+f[17]*alphaDrag[86]+f[6]*alphaDrag[73]+f[32]*alphaDrag[72]+f[33]*alphaDrag[71]+f[53]*alphaDrag[69]+f[18]*alphaDrag[67]+f[19]*alphaDrag[66]+f[47]*alphaDrag[65]+f[34]*alphaDrag[64]); 
  out[54] += 0.2165063509461096*(f[46]*alphaDrag[150]+f[57]*alphaDrag[137]+f[31]*alphaDrag[136]+f[30]*alphaDrag[135]+f[54]*alphaDrag[134]+f[45]*alphaDrag[131]+f[44]*alphaDrag[130]+f[16]*alphaDrag[129]+f[29]*alphaDrag[128]+f[50]*alphaDrag[86]+f[58]*alphaDrag[73]+f[37]*alphaDrag[72]+f[36]*alphaDrag[71]+f[54]*alphaDrag[69]+f[49]*alphaDrag[67]+f[48]*alphaDrag[66]+f[20]*alphaDrag[65]+f[35]*alphaDrag[64]+alphaDrag[9]*f[59]+alphaDrag[4]*f[54]+alphaDrag[22]*f[53]+alphaDrag[3]*f[52]+alphaDrag[2]*f[51]+alphaDrag[8]*f[40]+alphaDrag[7]*f[39]+alphaDrag[0]*f[38]+alphaDrag[1]*f[21]); 
  out[55] += 0.2165063509461096*(f[45]*alphaDrag[150]+f[31]*alphaDrag[137]+f[57]*alphaDrag[136]+f[29]*alphaDrag[135]+f[55]*alphaDrag[134]+f[46]*alphaDrag[131]+f[16]*alphaDrag[130]+f[44]*alphaDrag[129]+f[30]*alphaDrag[128]+f[49]*alphaDrag[86]+f[37]*alphaDrag[73]+f[58]*alphaDrag[72]+f[35]*alphaDrag[71]+f[55]*alphaDrag[69]+f[50]*alphaDrag[67]+f[20]*alphaDrag[66]+f[48]*alphaDrag[65]+f[36]*alphaDrag[64]+alphaDrag[8]*f[59]+alphaDrag[4]*f[55]+alphaDrag[3]*f[53]+alphaDrag[22]*f[52]+alphaDrag[1]*f[51]+alphaDrag[9]*f[40]+alphaDrag[0]*f[39]+alphaDrag[7]*f[38]+alphaDrag[2]*f[21]); 
  out[56] += 0.2165063509461096*(f[44]*alphaDrag[150]+f[30]*alphaDrag[137]+f[29]*alphaDrag[136]+f[57]*alphaDrag[135]+f[56]*alphaDrag[134]+f[16]*alphaDrag[131]+f[46]*alphaDrag[130]+f[45]*alphaDrag[129]+f[31]*alphaDrag[128]+f[48]*alphaDrag[86]+f[36]*alphaDrag[73]+f[35]*alphaDrag[72]+f[58]*alphaDrag[71]+f[56]*alphaDrag[69]+f[20]*alphaDrag[67]+f[50]*alphaDrag[66]+f[49]*alphaDrag[65]+f[37]*alphaDrag[64]+alphaDrag[7]*f[59]+alphaDrag[4]*f[56]+alphaDrag[2]*f[53]+alphaDrag[1]*f[52]+alphaDrag[22]*f[51]+alphaDrag[0]*f[40]+alphaDrag[9]*f[39]+alphaDrag[8]*f[38]+alphaDrag[3]*f[21]); 
  out[57] += 0.2165063509461096*(f[4]*alphaDrag[86]+f[10]*alphaDrag[73]+f[11]*alphaDrag[72]+f[12]*alphaDrag[71]+f[57]*alphaDrag[69]+f[23]*alphaDrag[67]+f[24]*alphaDrag[66]+f[25]*alphaDrag[65]+f[42]*alphaDrag[64]+alphaDrag[4]*f[57]+alphaDrag[0]*f[43]+alphaDrag[1]*f[28]+alphaDrag[2]*f[27]+alphaDrag[3]*f[26]+f[5]*alphaDrag[22]+alphaDrag[7]*f[15]+alphaDrag[8]*f[14]+alphaDrag[9]*f[13]); 
  out[58] += 0.2165063509461096*(f[4]*alphaDrag[150]+f[10]*alphaDrag[137]+f[11]*alphaDrag[136]+f[12]*alphaDrag[135]+f[58]*alphaDrag[134]+f[23]*alphaDrag[131]+f[24]*alphaDrag[130]+f[25]*alphaDrag[129]+f[42]*alphaDrag[128]+alphaDrag[4]*f[58]+alphaDrag[0]*f[47]+alphaDrag[1]*f[34]+alphaDrag[2]*f[33]+alphaDrag[3]*f[32]+f[6]*alphaDrag[22]+alphaDrag[7]*f[19]+alphaDrag[8]*f[18]+alphaDrag[9]*f[17]); 
  out[59] += 0.2165063509461096*(f[5]*alphaDrag[150]+f[13]*alphaDrag[137]+f[14]*alphaDrag[136]+f[15]*alphaDrag[135]+f[59]*alphaDrag[134]+f[26]*alphaDrag[131]+f[27]*alphaDrag[130]+f[28]*alphaDrag[129]+f[43]*alphaDrag[128]+f[6]*alphaDrag[86]+f[17]*alphaDrag[73]+f[18]*alphaDrag[72]+f[19]*alphaDrag[71]+f[59]*alphaDrag[69]+f[32]*alphaDrag[67]+f[33]*alphaDrag[66]+f[34]*alphaDrag[65]+f[47]*alphaDrag[64]); 
  out[60] += 0.2165063509461096*(f[31]*alphaDrag[150]+f[45]*alphaDrag[137]+f[46]*alphaDrag[136]+f[16]*alphaDrag[135]+f[60]*alphaDrag[134]+f[57]*alphaDrag[131]+f[29]*alphaDrag[130]+f[30]*alphaDrag[129]+f[44]*alphaDrag[128]+f[37]*alphaDrag[86]+f[49]*alphaDrag[73]+f[50]*alphaDrag[72]+f[20]*alphaDrag[71]+f[60]*alphaDrag[69]+f[58]*alphaDrag[67]+f[35]*alphaDrag[66]+f[36]*alphaDrag[65]+f[48]*alphaDrag[64]+alphaDrag[4]*f[60]+alphaDrag[3]*f[59]+alphaDrag[8]*f[53]+alphaDrag[9]*f[52]+alphaDrag[0]*f[51]+alphaDrag[22]*f[40]+alphaDrag[1]*f[39]+alphaDrag[2]*f[38]+alphaDrag[7]*f[21]); 
  out[61] += 0.2165063509461096*(f[30]*alphaDrag[150]+f[44]*alphaDrag[137]+f[16]*alphaDrag[136]+f[46]*alphaDrag[135]+f[61]*alphaDrag[134]+f[29]*alphaDrag[131]+f[57]*alphaDrag[130]+f[31]*alphaDrag[129]+f[45]*alphaDrag[128]+f[36]*alphaDrag[86]+f[48]*alphaDrag[73]+f[20]*alphaDrag[72]+f[50]*alphaDrag[71]+f[61]*alphaDrag[69]+f[35]*alphaDrag[67]+f[58]*alphaDrag[66]+f[37]*alphaDrag[65]+f[49]*alphaDrag[64]+alphaDrag[4]*f[61]+alphaDrag[2]*f[59]+alphaDrag[7]*f[53]+alphaDrag[0]*f[52]+alphaDrag[9]*f[51]+alphaDrag[1]*f[40]+alphaDrag[22]*f[39]+alphaDrag[3]*f[38]+alphaDrag[8]*f[21]); 
  out[62] += 0.2165063509461096*(f[29]*alphaDrag[150]+f[16]*alphaDrag[137]+f[44]*alphaDrag[136]+f[45]*alphaDrag[135]+f[62]*alphaDrag[134]+f[30]*alphaDrag[131]+f[31]*alphaDrag[130]+f[57]*alphaDrag[129]+f[46]*alphaDrag[128]+f[35]*alphaDrag[86]+f[20]*alphaDrag[73]+f[48]*alphaDrag[72]+f[49]*alphaDrag[71]+f[62]*alphaDrag[69]+f[36]*alphaDrag[67]+f[37]*alphaDrag[66]+f[58]*alphaDrag[65]+f[50]*alphaDrag[64]+alphaDrag[4]*f[62]+alphaDrag[1]*f[59]+alphaDrag[0]*f[53]+alphaDrag[7]*f[52]+alphaDrag[8]*f[51]+alphaDrag[2]*f[40]+alphaDrag[3]*f[39]+alphaDrag[22]*f[38]+alphaDrag[9]*f[21]); 
  out[63] += 0.2165063509461096*(f[16]*alphaDrag[150]+f[29]*alphaDrag[137]+f[30]*alphaDrag[136]+f[31]*alphaDrag[135]+f[63]*alphaDrag[134]+f[44]*alphaDrag[131]+f[45]*alphaDrag[130]+f[46]*alphaDrag[129]+f[57]*alphaDrag[128]+f[20]*alphaDrag[86]+f[35]*alphaDrag[73]+f[36]*alphaDrag[72]+f[37]*alphaDrag[71]+f[63]*alphaDrag[69]+f[48]*alphaDrag[67]+f[49]*alphaDrag[66]+f[50]*alphaDrag[65]+f[58]*alphaDrag[64]+alphaDrag[4]*f[63]+alphaDrag[0]*f[59]+alphaDrag[1]*f[53]+alphaDrag[2]*f[52]+alphaDrag[3]*f[51]+alphaDrag[7]*f[40]+alphaDrag[8]*f[39]+alphaDrag[9]*f[38]+f[21]*alphaDrag[22]); 

  return fabs(0.0625*alphaDrag[0])+fabs(0.0625*alphaDrag[64])+fabs(0.0625*alphaDrag[128])+fabs(0.4714045207910317*nuVtSqSum[0]*rdvxSq4)+fabs(0.4714045207910317*nuVtSqSum[0]*rdvySq4)+fabs(0.4714045207910317*nuVtSqSum[0]*rdvzSq4); 

} 
