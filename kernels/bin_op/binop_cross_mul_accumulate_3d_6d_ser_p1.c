#include <gkyl_binop_cross_mul_ser.h> 
 
GKYL_CU_DH
void
binop_cross_mul_accumulate_3d_6d_ser_p1(double a, const double *f, const double *g, double *fg) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  fg[0] += a*(0.3535533905932737*f[7]*g[22]+0.3535533905932737*f[6]*g[9]+0.3535533905932737*f[5]*g[8]+0.3535533905932737*f[4]*g[7]+0.3535533905932737*f[3]*g[3]+0.3535533905932737*f[2]*g[2]+0.3535533905932737*f[1]*g[1]+0.3535533905932737*f[0]*g[0]); 
  fg[1] += a*(0.3535533905932737*f[6]*g[22]+0.3535533905932737*f[7]*g[9]+0.3535533905932737*f[3]*g[8]+0.3535533905932737*f[2]*g[7]+0.3535533905932737*g[3]*f[5]+0.3535533905932737*g[2]*f[4]+0.3535533905932737*f[0]*g[1]+0.3535533905932737*g[0]*f[1]); 
  fg[2] += a*(0.3535533905932737*f[5]*g[22]+0.3535533905932737*f[3]*g[9]+0.3535533905932737*f[7]*g[8]+0.3535533905932737*f[1]*g[7]+0.3535533905932737*g[3]*f[6]+0.3535533905932737*g[1]*f[4]+0.3535533905932737*f[0]*g[2]+0.3535533905932737*g[0]*f[2]); 
  fg[3] += a*(0.3535533905932737*f[4]*g[22]+0.3535533905932737*f[2]*g[9]+0.3535533905932737*f[1]*g[8]+0.3535533905932737*f[7]*g[7]+0.3535533905932737*g[2]*f[6]+0.3535533905932737*g[1]*f[5]+0.3535533905932737*f[0]*g[3]+0.3535533905932737*g[0]*f[3]); 
  fg[4] += a*(0.3535533905932737*f[7]*g[42]+0.3535533905932737*f[6]*g[25]+0.3535533905932737*f[5]*g[24]+0.3535533905932737*f[4]*g[23]+0.3535533905932737*f[3]*g[12]+0.3535533905932737*f[2]*g[11]+0.3535533905932737*f[1]*g[10]+0.3535533905932737*f[0]*g[4]); 
  fg[5] += a*(0.3535533905932737*f[7]*g[43]+0.3535533905932737*f[6]*g[28]+0.3535533905932737*f[5]*g[27]+0.3535533905932737*f[4]*g[26]+0.3535533905932737*f[3]*g[15]+0.3535533905932737*f[2]*g[14]+0.3535533905932737*f[1]*g[13]+0.3535533905932737*f[0]*g[5]); 
  fg[6] += a*(0.3535533905932737*f[7]*g[47]+0.3535533905932737*f[6]*g[34]+0.3535533905932737*f[5]*g[33]+0.3535533905932737*f[4]*g[32]+0.3535533905932737*f[3]*g[19]+0.3535533905932737*f[2]*g[18]+0.3535533905932737*f[1]*g[17]+0.3535533905932737*f[0]*g[6]); 
  fg[7] += a*(0.3535533905932737*f[3]*g[22]+0.3535533905932737*f[5]*g[9]+0.3535533905932737*f[6]*g[8]+0.3535533905932737*f[0]*g[7]+0.3535533905932737*g[3]*f[7]+0.3535533905932737*g[0]*f[4]+0.3535533905932737*f[1]*g[2]+0.3535533905932737*g[1]*f[2]); 
  fg[8] += a*(0.3535533905932737*f[2]*g[22]+0.3535533905932737*f[4]*g[9]+0.3535533905932737*f[0]*g[8]+0.3535533905932737*f[6]*g[7]+0.3535533905932737*g[2]*f[7]+0.3535533905932737*g[0]*f[5]+0.3535533905932737*f[1]*g[3]+0.3535533905932737*g[1]*f[3]); 
  fg[9] += a*(0.3535533905932737*f[1]*g[22]+0.3535533905932737*f[0]*g[9]+0.3535533905932737*f[4]*g[8]+0.3535533905932737*f[5]*g[7]+0.3535533905932737*g[1]*f[7]+0.3535533905932737*g[0]*f[6]+0.3535533905932737*f[2]*g[3]+0.3535533905932737*g[2]*f[3]); 
  fg[10] += a*(0.3535533905932737*f[6]*g[42]+0.3535533905932737*f[7]*g[25]+0.3535533905932737*f[3]*g[24]+0.3535533905932737*f[2]*g[23]+0.3535533905932737*f[5]*g[12]+0.3535533905932737*f[4]*g[11]+0.3535533905932737*f[0]*g[10]+0.3535533905932737*f[1]*g[4]); 
  fg[11] += a*(0.3535533905932737*f[5]*g[42]+0.3535533905932737*f[3]*g[25]+0.3535533905932737*f[7]*g[24]+0.3535533905932737*f[1]*g[23]+0.3535533905932737*f[6]*g[12]+0.3535533905932737*f[0]*g[11]+0.3535533905932737*f[4]*g[10]+0.3535533905932737*f[2]*g[4]); 
  fg[12] += a*(0.3535533905932737*f[4]*g[42]+0.3535533905932737*f[2]*g[25]+0.3535533905932737*f[1]*g[24]+0.3535533905932737*f[7]*g[23]+0.3535533905932737*f[0]*g[12]+0.3535533905932737*f[6]*g[11]+0.3535533905932737*f[5]*g[10]+0.3535533905932737*f[3]*g[4]); 
  fg[13] += a*(0.3535533905932737*f[6]*g[43]+0.3535533905932737*f[7]*g[28]+0.3535533905932737*f[3]*g[27]+0.3535533905932737*f[2]*g[26]+0.3535533905932737*f[5]*g[15]+0.3535533905932737*f[4]*g[14]+0.3535533905932737*f[0]*g[13]+0.3535533905932737*f[1]*g[5]); 
  fg[14] += a*(0.3535533905932737*f[5]*g[43]+0.3535533905932737*f[3]*g[28]+0.3535533905932737*f[7]*g[27]+0.3535533905932737*f[1]*g[26]+0.3535533905932737*f[6]*g[15]+0.3535533905932737*f[0]*g[14]+0.3535533905932737*f[4]*g[13]+0.3535533905932737*f[2]*g[5]); 
  fg[15] += a*(0.3535533905932737*f[4]*g[43]+0.3535533905932737*f[2]*g[28]+0.3535533905932737*f[1]*g[27]+0.3535533905932737*f[7]*g[26]+0.3535533905932737*f[0]*g[15]+0.3535533905932737*f[6]*g[14]+0.3535533905932737*f[5]*g[13]+0.3535533905932737*f[3]*g[5]); 
  fg[16] += a*(0.3535533905932737*f[7]*g[57]+0.3535533905932737*f[6]*g[46]+0.3535533905932737*f[5]*g[45]+0.3535533905932737*f[4]*g[44]+0.3535533905932737*f[3]*g[31]+0.3535533905932737*f[2]*g[30]+0.3535533905932737*f[1]*g[29]+0.3535533905932737*f[0]*g[16]); 
  fg[17] += a*(0.3535533905932737*f[6]*g[47]+0.3535533905932737*f[7]*g[34]+0.3535533905932737*f[3]*g[33]+0.3535533905932737*f[2]*g[32]+0.3535533905932737*f[5]*g[19]+0.3535533905932737*f[4]*g[18]+0.3535533905932737*f[0]*g[17]+0.3535533905932737*f[1]*g[6]); 
  fg[18] += a*(0.3535533905932737*f[5]*g[47]+0.3535533905932737*f[3]*g[34]+0.3535533905932737*f[7]*g[33]+0.3535533905932737*f[1]*g[32]+0.3535533905932737*f[6]*g[19]+0.3535533905932737*f[0]*g[18]+0.3535533905932737*f[4]*g[17]+0.3535533905932737*f[2]*g[6]); 
  fg[19] += a*(0.3535533905932737*f[4]*g[47]+0.3535533905932737*f[2]*g[34]+0.3535533905932737*f[1]*g[33]+0.3535533905932737*f[7]*g[32]+0.3535533905932737*f[0]*g[19]+0.3535533905932737*f[6]*g[18]+0.3535533905932737*f[5]*g[17]+0.3535533905932737*f[3]*g[6]); 
  fg[20] += a*(0.3535533905932737*f[7]*g[58]+0.3535533905932737*f[6]*g[50]+0.3535533905932737*f[5]*g[49]+0.3535533905932737*f[4]*g[48]+0.3535533905932737*f[3]*g[37]+0.3535533905932737*f[2]*g[36]+0.3535533905932737*f[1]*g[35]+0.3535533905932737*f[0]*g[20]); 
  fg[21] += a*(0.3535533905932737*f[7]*g[59]+0.3535533905932737*f[6]*g[53]+0.3535533905932737*f[5]*g[52]+0.3535533905932737*f[4]*g[51]+0.3535533905932737*f[3]*g[40]+0.3535533905932737*f[2]*g[39]+0.3535533905932737*f[1]*g[38]+0.3535533905932737*f[0]*g[21]); 
  fg[22] += a*(0.3535533905932737*f[0]*g[22]+0.3535533905932737*f[1]*g[9]+0.3535533905932737*f[2]*g[8]+0.3535533905932737*f[3]*g[7]+0.3535533905932737*g[0]*f[7]+0.3535533905932737*g[1]*f[6]+0.3535533905932737*g[2]*f[5]+0.3535533905932737*g[3]*f[4]); 
  fg[23] += a*(0.3535533905932737*f[3]*g[42]+0.3535533905932737*f[5]*g[25]+0.3535533905932737*f[6]*g[24]+0.3535533905932737*f[0]*g[23]+0.3535533905932737*f[7]*g[12]+0.3535533905932737*f[1]*g[11]+0.3535533905932737*f[2]*g[10]+0.3535533905932737*f[4]*g[4]); 
  fg[24] += a*(0.3535533905932737*f[2]*g[42]+0.3535533905932737*f[4]*g[25]+0.3535533905932737*f[0]*g[24]+0.3535533905932737*f[6]*g[23]+0.3535533905932737*f[1]*g[12]+0.3535533905932737*f[7]*g[11]+0.3535533905932737*f[3]*g[10]+0.3535533905932737*g[4]*f[5]); 
  fg[25] += a*(0.3535533905932737*f[1]*g[42]+0.3535533905932737*f[0]*g[25]+0.3535533905932737*f[4]*g[24]+0.3535533905932737*f[5]*g[23]+0.3535533905932737*f[2]*g[12]+0.3535533905932737*f[3]*g[11]+0.3535533905932737*f[7]*g[10]+0.3535533905932737*g[4]*f[6]); 
  fg[26] += a*(0.3535533905932737*f[3]*g[43]+0.3535533905932737*f[5]*g[28]+0.3535533905932737*f[6]*g[27]+0.3535533905932737*f[0]*g[26]+0.3535533905932737*f[7]*g[15]+0.3535533905932737*f[1]*g[14]+0.3535533905932737*f[2]*g[13]+0.3535533905932737*f[4]*g[5]); 
  fg[27] += a*(0.3535533905932737*f[2]*g[43]+0.3535533905932737*f[4]*g[28]+0.3535533905932737*f[0]*g[27]+0.3535533905932737*f[6]*g[26]+0.3535533905932737*f[1]*g[15]+0.3535533905932737*f[7]*g[14]+0.3535533905932737*f[3]*g[13]+0.3535533905932737*f[5]*g[5]); 
  fg[28] += a*(0.3535533905932737*f[1]*g[43]+0.3535533905932737*f[0]*g[28]+0.3535533905932737*f[4]*g[27]+0.3535533905932737*f[5]*g[26]+0.3535533905932737*f[2]*g[15]+0.3535533905932737*f[3]*g[14]+0.3535533905932737*f[7]*g[13]+0.3535533905932737*g[5]*f[6]); 
  fg[29] += a*(0.3535533905932737*f[6]*g[57]+0.3535533905932737*f[7]*g[46]+0.3535533905932737*f[3]*g[45]+0.3535533905932737*f[2]*g[44]+0.3535533905932737*f[5]*g[31]+0.3535533905932737*f[4]*g[30]+0.3535533905932737*f[0]*g[29]+0.3535533905932737*f[1]*g[16]); 
  fg[30] += a*(0.3535533905932737*f[5]*g[57]+0.3535533905932737*f[3]*g[46]+0.3535533905932737*f[7]*g[45]+0.3535533905932737*f[1]*g[44]+0.3535533905932737*f[6]*g[31]+0.3535533905932737*f[0]*g[30]+0.3535533905932737*f[4]*g[29]+0.3535533905932737*f[2]*g[16]); 
  fg[31] += a*(0.3535533905932737*f[4]*g[57]+0.3535533905932737*f[2]*g[46]+0.3535533905932737*f[1]*g[45]+0.3535533905932737*f[7]*g[44]+0.3535533905932737*f[0]*g[31]+0.3535533905932737*f[6]*g[30]+0.3535533905932737*f[5]*g[29]+0.3535533905932737*f[3]*g[16]); 
  fg[32] += a*(0.3535533905932737*f[3]*g[47]+0.3535533905932737*f[5]*g[34]+0.3535533905932737*f[6]*g[33]+0.3535533905932737*f[0]*g[32]+0.3535533905932737*f[7]*g[19]+0.3535533905932737*f[1]*g[18]+0.3535533905932737*f[2]*g[17]+0.3535533905932737*f[4]*g[6]); 
  fg[33] += a*(0.3535533905932737*f[2]*g[47]+0.3535533905932737*f[4]*g[34]+0.3535533905932737*f[0]*g[33]+0.3535533905932737*f[6]*g[32]+0.3535533905932737*f[1]*g[19]+0.3535533905932737*f[7]*g[18]+0.3535533905932737*f[3]*g[17]+0.3535533905932737*f[5]*g[6]); 
  fg[34] += a*(0.3535533905932737*f[1]*g[47]+0.3535533905932737*f[0]*g[34]+0.3535533905932737*f[4]*g[33]+0.3535533905932737*f[5]*g[32]+0.3535533905932737*f[2]*g[19]+0.3535533905932737*f[3]*g[18]+0.3535533905932737*f[7]*g[17]+0.3535533905932737*f[6]*g[6]); 
  fg[35] += a*(0.3535533905932737*f[6]*g[58]+0.3535533905932737*f[7]*g[50]+0.3535533905932737*f[3]*g[49]+0.3535533905932737*f[2]*g[48]+0.3535533905932737*f[5]*g[37]+0.3535533905932737*f[4]*g[36]+0.3535533905932737*f[0]*g[35]+0.3535533905932737*f[1]*g[20]); 
  fg[36] += a*(0.3535533905932737*f[5]*g[58]+0.3535533905932737*f[3]*g[50]+0.3535533905932737*f[7]*g[49]+0.3535533905932737*f[1]*g[48]+0.3535533905932737*f[6]*g[37]+0.3535533905932737*f[0]*g[36]+0.3535533905932737*f[4]*g[35]+0.3535533905932737*f[2]*g[20]); 
  fg[37] += a*(0.3535533905932737*f[4]*g[58]+0.3535533905932737*f[2]*g[50]+0.3535533905932737*f[1]*g[49]+0.3535533905932737*f[7]*g[48]+0.3535533905932737*f[0]*g[37]+0.3535533905932737*f[6]*g[36]+0.3535533905932737*f[5]*g[35]+0.3535533905932737*f[3]*g[20]); 
  fg[38] += a*(0.3535533905932737*f[6]*g[59]+0.3535533905932737*f[7]*g[53]+0.3535533905932737*f[3]*g[52]+0.3535533905932737*f[2]*g[51]+0.3535533905932737*f[5]*g[40]+0.3535533905932737*f[4]*g[39]+0.3535533905932737*f[0]*g[38]+0.3535533905932737*f[1]*g[21]); 
  fg[39] += a*(0.3535533905932737*f[5]*g[59]+0.3535533905932737*f[3]*g[53]+0.3535533905932737*f[7]*g[52]+0.3535533905932737*f[1]*g[51]+0.3535533905932737*f[6]*g[40]+0.3535533905932737*f[0]*g[39]+0.3535533905932737*f[4]*g[38]+0.3535533905932737*f[2]*g[21]); 
  fg[40] += a*(0.3535533905932737*f[4]*g[59]+0.3535533905932737*f[2]*g[53]+0.3535533905932737*f[1]*g[52]+0.3535533905932737*f[7]*g[51]+0.3535533905932737*f[0]*g[40]+0.3535533905932737*f[6]*g[39]+0.3535533905932737*f[5]*g[38]+0.3535533905932737*f[3]*g[21]); 
  fg[41] += a*(0.3535533905932737*f[7]*g[63]+0.3535533905932737*f[6]*g[62]+0.3535533905932737*f[5]*g[61]+0.3535533905932737*f[4]*g[60]+0.3535533905932737*f[3]*g[56]+0.3535533905932737*f[2]*g[55]+0.3535533905932737*f[1]*g[54]+0.3535533905932737*f[0]*g[41]); 
  fg[42] += a*(0.3535533905932737*f[0]*g[42]+0.3535533905932737*f[1]*g[25]+0.3535533905932737*f[2]*g[24]+0.3535533905932737*f[3]*g[23]+0.3535533905932737*f[4]*g[12]+0.3535533905932737*f[5]*g[11]+0.3535533905932737*f[6]*g[10]+0.3535533905932737*g[4]*f[7]); 
  fg[43] += a*(0.3535533905932737*f[0]*g[43]+0.3535533905932737*f[1]*g[28]+0.3535533905932737*f[2]*g[27]+0.3535533905932737*f[3]*g[26]+0.3535533905932737*f[4]*g[15]+0.3535533905932737*f[5]*g[14]+0.3535533905932737*f[6]*g[13]+0.3535533905932737*g[5]*f[7]); 
  fg[44] += a*(0.3535533905932737*f[3]*g[57]+0.3535533905932737*f[5]*g[46]+0.3535533905932737*f[6]*g[45]+0.3535533905932737*f[0]*g[44]+0.3535533905932737*f[7]*g[31]+0.3535533905932737*f[1]*g[30]+0.3535533905932737*f[2]*g[29]+0.3535533905932737*f[4]*g[16]); 
  fg[45] += a*(0.3535533905932737*f[2]*g[57]+0.3535533905932737*f[4]*g[46]+0.3535533905932737*f[0]*g[45]+0.3535533905932737*f[6]*g[44]+0.3535533905932737*f[1]*g[31]+0.3535533905932737*f[7]*g[30]+0.3535533905932737*f[3]*g[29]+0.3535533905932737*f[5]*g[16]); 
  fg[46] += a*(0.3535533905932737*f[1]*g[57]+0.3535533905932737*f[0]*g[46]+0.3535533905932737*f[4]*g[45]+0.3535533905932737*f[5]*g[44]+0.3535533905932737*f[2]*g[31]+0.3535533905932737*f[3]*g[30]+0.3535533905932737*f[7]*g[29]+0.3535533905932737*f[6]*g[16]); 
  fg[47] += a*(0.3535533905932737*f[0]*g[47]+0.3535533905932737*f[1]*g[34]+0.3535533905932737*f[2]*g[33]+0.3535533905932737*f[3]*g[32]+0.3535533905932737*f[4]*g[19]+0.3535533905932737*f[5]*g[18]+0.3535533905932737*f[6]*g[17]+0.3535533905932737*g[6]*f[7]); 
  fg[48] += a*(0.3535533905932737*f[3]*g[58]+0.3535533905932737*f[5]*g[50]+0.3535533905932737*f[6]*g[49]+0.3535533905932737*f[0]*g[48]+0.3535533905932737*f[7]*g[37]+0.3535533905932737*f[1]*g[36]+0.3535533905932737*f[2]*g[35]+0.3535533905932737*f[4]*g[20]); 
  fg[49] += a*(0.3535533905932737*f[2]*g[58]+0.3535533905932737*f[4]*g[50]+0.3535533905932737*f[0]*g[49]+0.3535533905932737*f[6]*g[48]+0.3535533905932737*f[1]*g[37]+0.3535533905932737*f[7]*g[36]+0.3535533905932737*f[3]*g[35]+0.3535533905932737*f[5]*g[20]); 
  fg[50] += a*(0.3535533905932737*f[1]*g[58]+0.3535533905932737*f[0]*g[50]+0.3535533905932737*f[4]*g[49]+0.3535533905932737*f[5]*g[48]+0.3535533905932737*f[2]*g[37]+0.3535533905932737*f[3]*g[36]+0.3535533905932737*f[7]*g[35]+0.3535533905932737*f[6]*g[20]); 
  fg[51] += a*(0.3535533905932737*f[3]*g[59]+0.3535533905932737*f[5]*g[53]+0.3535533905932737*f[6]*g[52]+0.3535533905932737*f[0]*g[51]+0.3535533905932737*f[7]*g[40]+0.3535533905932737*f[1]*g[39]+0.3535533905932737*f[2]*g[38]+0.3535533905932737*f[4]*g[21]); 
  fg[52] += a*(0.3535533905932737*f[2]*g[59]+0.3535533905932737*f[4]*g[53]+0.3535533905932737*f[0]*g[52]+0.3535533905932737*f[6]*g[51]+0.3535533905932737*f[1]*g[40]+0.3535533905932737*f[7]*g[39]+0.3535533905932737*f[3]*g[38]+0.3535533905932737*f[5]*g[21]); 
  fg[53] += a*(0.3535533905932737*f[1]*g[59]+0.3535533905932737*f[0]*g[53]+0.3535533905932737*f[4]*g[52]+0.3535533905932737*f[5]*g[51]+0.3535533905932737*f[2]*g[40]+0.3535533905932737*f[3]*g[39]+0.3535533905932737*f[7]*g[38]+0.3535533905932737*f[6]*g[21]); 
  fg[54] += a*(0.3535533905932737*f[6]*g[63]+0.3535533905932737*f[7]*g[62]+0.3535533905932737*f[3]*g[61]+0.3535533905932737*f[2]*g[60]+0.3535533905932737*f[5]*g[56]+0.3535533905932737*f[4]*g[55]+0.3535533905932737*f[0]*g[54]+0.3535533905932737*f[1]*g[41]); 
  fg[55] += a*(0.3535533905932737*f[5]*g[63]+0.3535533905932737*f[3]*g[62]+0.3535533905932737*f[7]*g[61]+0.3535533905932737*f[1]*g[60]+0.3535533905932737*f[6]*g[56]+0.3535533905932737*f[0]*g[55]+0.3535533905932737*f[4]*g[54]+0.3535533905932737*f[2]*g[41]); 
  fg[56] += a*(0.3535533905932737*f[4]*g[63]+0.3535533905932737*f[2]*g[62]+0.3535533905932737*f[1]*g[61]+0.3535533905932737*f[7]*g[60]+0.3535533905932737*f[0]*g[56]+0.3535533905932737*f[6]*g[55]+0.3535533905932737*f[5]*g[54]+0.3535533905932737*f[3]*g[41]); 
  fg[57] += a*(0.3535533905932737*f[0]*g[57]+0.3535533905932737*f[1]*g[46]+0.3535533905932737*f[2]*g[45]+0.3535533905932737*f[3]*g[44]+0.3535533905932737*f[4]*g[31]+0.3535533905932737*f[5]*g[30]+0.3535533905932737*f[6]*g[29]+0.3535533905932737*f[7]*g[16]); 
  fg[58] += a*(0.3535533905932737*f[0]*g[58]+0.3535533905932737*f[1]*g[50]+0.3535533905932737*f[2]*g[49]+0.3535533905932737*f[3]*g[48]+0.3535533905932737*f[4]*g[37]+0.3535533905932737*f[5]*g[36]+0.3535533905932737*f[6]*g[35]+0.3535533905932737*f[7]*g[20]); 
  fg[59] += a*(0.3535533905932737*f[0]*g[59]+0.3535533905932737*f[1]*g[53]+0.3535533905932737*f[2]*g[52]+0.3535533905932737*f[3]*g[51]+0.3535533905932737*f[4]*g[40]+0.3535533905932737*f[5]*g[39]+0.3535533905932737*f[6]*g[38]+0.3535533905932737*f[7]*g[21]); 
  fg[60] += a*(0.3535533905932737*f[3]*g[63]+0.3535533905932737*f[5]*g[62]+0.3535533905932737*f[6]*g[61]+0.3535533905932737*f[0]*g[60]+0.3535533905932737*f[7]*g[56]+0.3535533905932737*f[1]*g[55]+0.3535533905932737*f[2]*g[54]+0.3535533905932737*f[4]*g[41]); 
  fg[61] += a*(0.3535533905932737*f[2]*g[63]+0.3535533905932737*f[4]*g[62]+0.3535533905932737*f[0]*g[61]+0.3535533905932737*f[6]*g[60]+0.3535533905932737*f[1]*g[56]+0.3535533905932737*f[7]*g[55]+0.3535533905932737*f[3]*g[54]+0.3535533905932737*f[5]*g[41]); 
  fg[62] += a*(0.3535533905932737*f[1]*g[63]+0.3535533905932737*f[0]*g[62]+0.3535533905932737*f[4]*g[61]+0.3535533905932737*f[5]*g[60]+0.3535533905932737*f[2]*g[56]+0.3535533905932737*f[3]*g[55]+0.3535533905932737*f[7]*g[54]+0.3535533905932737*f[6]*g[41]); 
  fg[63] += a*(0.3535533905932737*f[0]*g[63]+0.3535533905932737*f[1]*g[62]+0.3535533905932737*f[2]*g[61]+0.3535533905932737*f[3]*g[60]+0.3535533905932737*f[4]*g[56]+0.3535533905932737*f[5]*g[55]+0.3535533905932737*f[6]*g[54]+0.3535533905932737*f[7]*g[41]); 
} 