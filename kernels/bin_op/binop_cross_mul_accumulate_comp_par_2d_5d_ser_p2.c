#include <gkyl_binop_cross_mul_ser.h> 
 
GKYL_CU_DH
void
binop_cross_mul_accumulate_comp_par_2d_5d_ser_p2(double a, const double *f, const double *g, double *fg, int linc2) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  switch (linc2) { 
    case 0: 
      fg[0] += a*(0.5*f[7]*g[32]+0.5*f[6]*g[31]+0.5*f[5]*g[17]+0.5*f[4]*g[16]+0.5*f[3]*g[6]+0.5*f[2]*g[2]+0.5*f[1]*g[1]+0.5*f[0]*g[0]); 
    case 1: 
      fg[1] += a*(0.5000000000000001*f[5]*g[32]+0.447213595499958*f[3]*g[31]+0.5000000000000001*f[7]*g[17]+0.4472135954999579*f[1]*g[16]+0.447213595499958*f[6]*g[6]+0.5*f[2]*g[6]+0.4472135954999579*g[1]*f[4]+0.5*g[2]*f[3]+0.5*f[0]*g[1]+0.5*g[0]*f[1]); 
    case 2: 
      fg[2] += a*(0.447213595499958*f[3]*g[32]+0.5000000000000001*f[4]*g[31]+0.4472135954999579*f[2]*g[17]+0.5000000000000001*f[6]*g[16]+0.447213595499958*g[6]*f[7]+0.5*f[1]*g[6]+0.4472135954999579*g[2]*f[5]+0.5*g[1]*f[3]+0.5*f[0]*g[2]+0.5*g[0]*f[2]); 
    case 3: 
      fg[3] += a*(0.5000000000000001*f[7]*g[57]+0.5000000000000001*f[6]*g[56]+0.5000000000000001*f[5]*g[34]+0.5000000000000001*f[4]*g[33]+0.5*f[3]*g[21]+0.5*f[2]*g[8]+0.5*f[1]*g[7]+0.5*f[0]*g[3]); 
    case 4: 
      fg[4] += a*(0.5000000000000001*f[7]*g[60]+0.5000000000000001*f[6]*g[59]+0.5000000000000001*f[5]*g[38]+0.5000000000000001*f[4]*g[37]+0.5*f[3]*g[22]+0.5*f[2]*g[10]+0.5*f[1]*g[9]+0.5*f[0]*g[4]); 
    case 5: 
      fg[5] += a*(0.5000000000000001*f[7]*g[69]+0.5000000000000001*f[6]*g[68]+0.5000000000000001*f[5]*g[44]+0.5000000000000001*f[4]*g[43]+0.5*f[3]*g[25]+0.5*f[2]*g[13]+0.5*f[1]*g[12]+0.5*f[0]*g[5]); 
    case 6: 
      fg[6] += a*(0.4*f[6]*g[32]+0.447213595499958*f[2]*g[32]+0.4*f[7]*g[31]+0.447213595499958*f[1]*g[31]+0.4472135954999579*f[3]*g[17]+0.4472135954999579*f[3]*g[16]+0.447213595499958*g[2]*f[7]+0.4472135954999579*f[5]*g[6]+0.4472135954999579*f[4]*g[6]+0.5*f[0]*g[6]+0.447213595499958*g[1]*f[6]+0.5*g[0]*f[3]+0.5*f[1]*g[2]+0.5*g[1]*f[2]); 
    case 7: 
      fg[7] += a*(0.5*f[5]*g[57]+0.4472135954999579*f[3]*g[56]+0.5*f[7]*g[34]+0.447213595499958*f[1]*g[33]+0.447213595499958*f[6]*g[21]+0.5*f[2]*g[21]+0.5*f[3]*g[8]+0.4472135954999579*f[4]*g[7]+0.5*f[0]*g[7]+0.5*f[1]*g[3]); 
    case 8: 
      fg[8] += a*(0.4472135954999579*f[3]*g[57]+0.5*f[4]*g[56]+0.447213595499958*f[2]*g[34]+0.5*f[6]*g[33]+0.447213595499958*f[7]*g[21]+0.5*f[1]*g[21]+0.4472135954999579*f[5]*g[8]+0.5*f[0]*g[8]+0.5*f[3]*g[7]+0.5*f[2]*g[3]); 
    case 9: 
      fg[9] += a*(0.5*f[5]*g[60]+0.4472135954999579*f[3]*g[59]+0.5*f[7]*g[38]+0.447213595499958*f[1]*g[37]+0.447213595499958*f[6]*g[22]+0.5*f[2]*g[22]+0.5*f[3]*g[10]+0.4472135954999579*f[4]*g[9]+0.5*f[0]*g[9]+0.5*f[1]*g[4]); 
    case 10: 
      fg[10] += a*(0.4472135954999579*f[3]*g[60]+0.5*f[4]*g[59]+0.447213595499958*f[2]*g[38]+0.5*f[6]*g[37]+0.447213595499958*f[7]*g[22]+0.5*f[1]*g[22]+0.4472135954999579*f[5]*g[10]+0.5*f[0]*g[10]+0.5*f[3]*g[9]+0.5*f[2]*g[4]); 
    case 11: 
      fg[11] += a*(0.5*f[7]*g[88]+0.5*f[6]*g[87]+0.5*f[5]*g[62]+0.5*f[4]*g[61]+0.5*f[3]*g[51]+0.5*f[2]*g[24]+0.5*f[1]*g[23]+0.5*f[0]*g[11]); 
    case 12: 
      fg[12] += a*(0.5*f[5]*g[69]+0.4472135954999579*f[3]*g[68]+0.5*f[7]*g[44]+0.447213595499958*f[1]*g[43]+0.447213595499958*f[6]*g[25]+0.5*f[2]*g[25]+0.5*f[3]*g[13]+0.4472135954999579*f[4]*g[12]+0.5*f[0]*g[12]+0.5*f[1]*g[5]); 
    case 13: 
      fg[13] += a*(0.4472135954999579*f[3]*g[69]+0.5*f[4]*g[68]+0.447213595499958*f[2]*g[44]+0.5*f[6]*g[43]+0.447213595499958*f[7]*g[25]+0.5*f[1]*g[25]+0.4472135954999579*f[5]*g[13]+0.5*f[0]*g[13]+0.5*f[3]*g[12]+0.5*f[2]*g[5]); 
    case 14: 
      fg[14] += a*(0.5*f[7]*g[92]+0.5*f[6]*g[91]+0.5*f[5]*g[71]+0.5*f[4]*g[70]+0.5*f[3]*g[52]+0.5*f[2]*g[27]+0.5*f[1]*g[26]+0.5*f[0]*g[14]); 
    case 15: 
      fg[15] += a*(0.5*f[7]*g[95]+0.5*f[6]*g[94]+0.5*f[5]*g[75]+0.5*f[4]*g[74]+0.5*f[3]*g[53]+0.5*f[2]*g[29]+0.5*f[1]*g[28]+0.5*f[0]*g[15]); 
    case 16: 
      fg[16] += a*(0.4472135954999579*f[7]*g[32]+0.31943828249997*f[6]*g[31]+0.5000000000000001*f[2]*g[31]+0.31943828249997*f[4]*g[16]+0.5*f[0]*g[16]+0.4472135954999579*f[3]*g[6]+0.5000000000000001*g[2]*f[6]+0.5*g[0]*f[4]+0.4472135954999579*f[1]*g[1]); 
    case 17: 
      fg[17] += a*(0.31943828249997*f[7]*g[32]+0.5000000000000001*f[1]*g[32]+0.4472135954999579*f[6]*g[31]+0.31943828249997*f[5]*g[17]+0.5*f[0]*g[17]+0.5000000000000001*g[1]*f[7]+0.4472135954999579*f[3]*g[6]+0.5*g[0]*f[5]+0.4472135954999579*f[2]*g[2]); 
    case 18: 
      fg[18] += a*(0.5*f[3]*g[58]+0.5000000000000001*f[2]*g[36]+0.5000000000000001*f[1]*g[35]+0.5*f[0]*g[18]); 
    case 19: 
      fg[19] += a*(0.5*f[3]*g[65]+0.5000000000000001*f[2]*g[41]+0.5000000000000001*f[1]*g[40]+0.5*f[0]*g[19]); 
    case 20: 
      fg[20] += a*(0.5*f[3]*g[80]+0.5000000000000001*f[2]*g[48]+0.5000000000000001*f[1]*g[47]+0.5*f[0]*g[20]); 
    case 21: 
      fg[21] += a*(0.4*f[6]*g[57]+0.4472135954999579*f[2]*g[57]+0.4*f[7]*g[56]+0.4472135954999579*f[1]*g[56]+0.447213595499958*f[3]*g[34]+0.447213595499958*f[3]*g[33]+0.4472135954999579*f[5]*g[21]+0.4472135954999579*f[4]*g[21]+0.5*f[0]*g[21]+0.447213595499958*f[7]*g[8]+0.5*f[1]*g[8]+0.447213595499958*f[6]*g[7]+0.5*f[2]*g[7]+0.5*f[3]*g[3]); 
    case 22: 
      fg[22] += a*(0.4*f[6]*g[60]+0.4472135954999579*f[2]*g[60]+0.4*f[7]*g[59]+0.4472135954999579*f[1]*g[59]+0.447213595499958*f[3]*g[38]+0.447213595499958*f[3]*g[37]+0.4472135954999579*f[5]*g[22]+0.4472135954999579*f[4]*g[22]+0.5*f[0]*g[22]+0.447213595499958*f[7]*g[10]+0.5*f[1]*g[10]+0.447213595499958*f[6]*g[9]+0.5*f[2]*g[9]+0.5*f[3]*g[4]); 
    case 23: 
      fg[23] += a*(0.5000000000000001*f[5]*g[88]+0.447213595499958*f[3]*g[87]+0.5000000000000001*f[7]*g[62]+0.4472135954999579*f[1]*g[61]+0.447213595499958*f[6]*g[51]+0.5*f[2]*g[51]+0.5*f[3]*g[24]+0.4472135954999579*f[4]*g[23]+0.5*f[0]*g[23]+0.5*f[1]*g[11]); 
    case 24: 
      fg[24] += a*(0.447213595499958*f[3]*g[88]+0.5000000000000001*f[4]*g[87]+0.4472135954999579*f[2]*g[62]+0.5000000000000001*f[6]*g[61]+0.447213595499958*f[7]*g[51]+0.5*f[1]*g[51]+0.4472135954999579*f[5]*g[24]+0.5*f[0]*g[24]+0.5*f[3]*g[23]+0.5*f[2]*g[11]); 
    case 25: 
      fg[25] += a*(0.4*f[6]*g[69]+0.4472135954999579*f[2]*g[69]+0.4*f[7]*g[68]+0.4472135954999579*f[1]*g[68]+0.447213595499958*f[3]*g[44]+0.447213595499958*f[3]*g[43]+0.4472135954999579*f[5]*g[25]+0.4472135954999579*f[4]*g[25]+0.5*f[0]*g[25]+0.447213595499958*f[7]*g[13]+0.5*f[1]*g[13]+0.447213595499958*f[6]*g[12]+0.5*f[2]*g[12]+0.5*f[3]*g[5]); 
    case 26: 
      fg[26] += a*(0.5000000000000001*f[5]*g[92]+0.447213595499958*f[3]*g[91]+0.5000000000000001*f[7]*g[71]+0.4472135954999579*f[1]*g[70]+0.447213595499958*f[6]*g[52]+0.5*f[2]*g[52]+0.5*f[3]*g[27]+0.4472135954999579*f[4]*g[26]+0.5*f[0]*g[26]+0.5*f[1]*g[14]); 
    case 27: 
      fg[27] += a*(0.447213595499958*f[3]*g[92]+0.5000000000000001*f[4]*g[91]+0.4472135954999579*f[2]*g[71]+0.5000000000000001*f[6]*g[70]+0.447213595499958*f[7]*g[52]+0.5*f[1]*g[52]+0.4472135954999579*f[5]*g[27]+0.5*f[0]*g[27]+0.5*f[3]*g[26]+0.5*f[2]*g[14]); 
    case 28: 
      fg[28] += a*(0.5000000000000001*f[5]*g[95]+0.447213595499958*f[3]*g[94]+0.5000000000000001*f[7]*g[75]+0.4472135954999579*f[1]*g[74]+0.447213595499958*f[6]*g[53]+0.5*f[2]*g[53]+0.5*f[3]*g[29]+0.4472135954999579*f[4]*g[28]+0.5*f[0]*g[28]+0.5*f[1]*g[15]); 
    case 29: 
      fg[29] += a*(0.447213595499958*f[3]*g[95]+0.5000000000000001*f[4]*g[94]+0.4472135954999579*f[2]*g[75]+0.5000000000000001*f[6]*g[74]+0.447213595499958*f[7]*g[53]+0.5*f[1]*g[53]+0.4472135954999579*f[5]*g[29]+0.5*f[0]*g[29]+0.5*f[3]*g[28]+0.5*f[2]*g[15]); 
    case 30: 
      fg[30] += a*(0.5000000000000001*f[7]*g[108]+0.5000000000000001*f[6]*g[107]+0.5000000000000001*f[5]*g[97]+0.5000000000000001*f[4]*g[96]+0.5*f[3]*g[86]+0.5*f[2]*g[55]+0.5*f[1]*g[54]+0.5*f[0]*g[30]); 
    case 31: 
      fg[31] += a*(0.4*f[3]*g[32]+0.4472135954999579*f[5]*g[31]+0.31943828249997*f[4]*g[31]+0.5*f[0]*g[31]+0.4472135954999579*f[6]*g[17]+0.31943828249997*f[6]*g[16]+0.5000000000000001*f[2]*g[16]+0.4*g[6]*f[7]+0.447213595499958*f[1]*g[6]+0.5*g[0]*f[6]+0.5000000000000001*g[2]*f[4]+0.447213595499958*g[1]*f[3]); 
    case 32: 
      fg[32] += a*(0.31943828249997*f[5]*g[32]+0.4472135954999579*f[4]*g[32]+0.5*f[0]*g[32]+0.4*f[3]*g[31]+0.31943828249997*f[7]*g[17]+0.5000000000000001*f[1]*g[17]+0.4472135954999579*f[7]*g[16]+0.5*g[0]*f[7]+0.4*f[6]*g[6]+0.447213595499958*f[2]*g[6]+0.5000000000000001*g[1]*f[5]+0.447213595499958*g[2]*f[3]); 
    case 33: 
      fg[33] += a*(0.4472135954999579*f[7]*g[57]+0.31943828249997*f[6]*g[56]+0.5000000000000001*f[2]*g[56]+0.31943828249997*f[4]*g[33]+0.5*f[0]*g[33]+0.447213595499958*f[3]*g[21]+0.5*f[6]*g[8]+0.447213595499958*f[1]*g[7]+0.5000000000000001*g[3]*f[4]); 
    case 34: 
      fg[34] += a*(0.31943828249997*f[7]*g[57]+0.5000000000000001*f[1]*g[57]+0.4472135954999579*f[6]*g[56]+0.31943828249997*f[5]*g[34]+0.5*f[0]*g[34]+0.447213595499958*f[3]*g[21]+0.447213595499958*f[2]*g[8]+0.5*f[7]*g[7]+0.5000000000000001*g[3]*f[5]); 
    case 35: 
      fg[35] += a*(0.4472135954999579*f[6]*g[58]+0.5000000000000001*f[2]*g[58]+0.5*f[3]*g[36]+0.4472135954999579*f[4]*g[35]+0.5*f[0]*g[35]+0.5000000000000001*f[1]*g[18]); 
    case 36: 
      fg[36] += a*(0.4472135954999579*f[7]*g[58]+0.5000000000000001*f[1]*g[58]+0.4472135954999579*f[5]*g[36]+0.5*f[0]*g[36]+0.5*f[3]*g[35]+0.5000000000000001*f[2]*g[18]); 
    case 37: 
      fg[37] += a*(0.4472135954999579*f[7]*g[60]+0.31943828249997*f[6]*g[59]+0.5000000000000001*f[2]*g[59]+0.31943828249997*f[4]*g[37]+0.5*f[0]*g[37]+0.447213595499958*f[3]*g[22]+0.5*f[6]*g[10]+0.447213595499958*f[1]*g[9]+0.5000000000000001*f[4]*g[4]); 
    case 38: 
      fg[38] += a*(0.31943828249997*f[7]*g[60]+0.5000000000000001*f[1]*g[60]+0.4472135954999579*f[6]*g[59]+0.31943828249997*f[5]*g[38]+0.5*f[0]*g[38]+0.447213595499958*f[3]*g[22]+0.447213595499958*f[2]*g[10]+0.5*f[7]*g[9]+0.5000000000000001*g[4]*f[5]); 
    case 39: 
      fg[39] += a*(0.5*f[3]*g[89]+0.5000000000000001*f[2]*g[64]+0.5000000000000001*f[1]*g[63]+0.5*f[0]*g[39]); 
    case 40: 
      fg[40] += a*(0.4472135954999579*f[6]*g[65]+0.5000000000000001*f[2]*g[65]+0.5*f[3]*g[41]+0.4472135954999579*f[4]*g[40]+0.5*f[0]*g[40]+0.5000000000000001*f[1]*g[19]); 
    case 41: 
      fg[41] += a*(0.4472135954999579*f[7]*g[65]+0.5000000000000001*f[1]*g[65]+0.4472135954999579*f[5]*g[41]+0.5*f[0]*g[41]+0.5*f[3]*g[40]+0.5000000000000001*f[2]*g[19]); 
    case 42: 
      fg[42] += a*(0.5*f[3]*g[90]+0.5000000000000001*f[2]*g[67]+0.5000000000000001*f[1]*g[66]+0.5*f[0]*g[42]); 
    case 43: 
      fg[43] += a*(0.4472135954999579*f[7]*g[69]+0.31943828249997*f[6]*g[68]+0.5000000000000001*f[2]*g[68]+0.31943828249997*f[4]*g[43]+0.5*f[0]*g[43]+0.447213595499958*f[3]*g[25]+0.5*f[6]*g[13]+0.447213595499958*f[1]*g[12]+0.5000000000000001*f[4]*g[5]); 
    case 44: 
      fg[44] += a*(0.31943828249997*f[7]*g[69]+0.5000000000000001*f[1]*g[69]+0.4472135954999579*f[6]*g[68]+0.31943828249997*f[5]*g[44]+0.5*f[0]*g[44]+0.447213595499958*f[3]*g[25]+0.447213595499958*f[2]*g[13]+0.5*f[7]*g[12]+0.5000000000000001*f[5]*g[5]); 
    case 45: 
      fg[45] += a*(0.5*f[3]*g[93]+0.5000000000000001*f[2]*g[73]+0.5000000000000001*f[1]*g[72]+0.5*f[0]*g[45]); 
    case 46: 
      fg[46] += a*(0.5*f[3]*g[100]+0.5000000000000001*f[2]*g[78]+0.5000000000000001*f[1]*g[77]+0.5*f[0]*g[46]); 
    case 47: 
      fg[47] += a*(0.4472135954999579*f[6]*g[80]+0.5000000000000001*f[2]*g[80]+0.5*f[3]*g[48]+0.4472135954999579*f[4]*g[47]+0.5*f[0]*g[47]+0.5000000000000001*f[1]*g[20]); 
    case 48: 
      fg[48] += a*(0.4472135954999579*f[7]*g[80]+0.5000000000000001*f[1]*g[80]+0.4472135954999579*f[5]*g[48]+0.5*f[0]*g[48]+0.5*f[3]*g[47]+0.5000000000000001*f[2]*g[20]); 
    case 49: 
      fg[49] += a*(0.5*f[3]*g[103]+0.5000000000000001*f[2]*g[82]+0.5000000000000001*f[1]*g[81]+0.5*f[0]*g[49]); 
    case 50: 
      fg[50] += a*(0.5*f[3]*g[104]+0.5000000000000001*f[2]*g[84]+0.5000000000000001*f[1]*g[83]+0.5*f[0]*g[50]); 
    case 51: 
      fg[51] += a*(0.4*f[6]*g[88]+0.447213595499958*f[2]*g[88]+0.4*f[7]*g[87]+0.447213595499958*f[1]*g[87]+0.4472135954999579*f[3]*g[62]+0.4472135954999579*f[3]*g[61]+0.4472135954999579*f[5]*g[51]+0.4472135954999579*f[4]*g[51]+0.5*f[0]*g[51]+0.447213595499958*f[7]*g[24]+0.5*f[1]*g[24]+0.447213595499958*f[6]*g[23]+0.5*f[2]*g[23]+0.5*f[3]*g[11]); 
    case 52: 
      fg[52] += a*(0.4*f[6]*g[92]+0.447213595499958*f[2]*g[92]+0.4*f[7]*g[91]+0.447213595499958*f[1]*g[91]+0.4472135954999579*f[3]*g[71]+0.4472135954999579*f[3]*g[70]+0.4472135954999579*f[5]*g[52]+0.4472135954999579*f[4]*g[52]+0.5*f[0]*g[52]+0.447213595499958*f[7]*g[27]+0.5*f[1]*g[27]+0.447213595499958*f[6]*g[26]+0.5*f[2]*g[26]+0.5*f[3]*g[14]); 
    case 53: 
      fg[53] += a*(0.4*f[6]*g[95]+0.447213595499958*f[2]*g[95]+0.4*f[7]*g[94]+0.447213595499958*f[1]*g[94]+0.4472135954999579*f[3]*g[75]+0.4472135954999579*f[3]*g[74]+0.4472135954999579*f[5]*g[53]+0.4472135954999579*f[4]*g[53]+0.5*f[0]*g[53]+0.447213595499958*f[7]*g[29]+0.5*f[1]*g[29]+0.447213595499958*f[6]*g[28]+0.5*f[2]*g[28]+0.5*f[3]*g[15]); 
    case 54: 
      fg[54] += a*(0.5*f[5]*g[108]+0.4472135954999579*f[3]*g[107]+0.5*f[7]*g[97]+0.447213595499958*f[1]*g[96]+0.447213595499958*f[6]*g[86]+0.5*f[2]*g[86]+0.5*f[3]*g[55]+0.4472135954999579*f[4]*g[54]+0.5*f[0]*g[54]+0.5*f[1]*g[30]); 
    case 55: 
      fg[55] += a*(0.4472135954999579*f[3]*g[108]+0.5*f[4]*g[107]+0.447213595499958*f[2]*g[97]+0.5*f[6]*g[96]+0.447213595499958*f[7]*g[86]+0.5*f[1]*g[86]+0.4472135954999579*f[5]*g[55]+0.5*f[0]*g[55]+0.5*f[3]*g[54]+0.5*f[2]*g[30]); 
    case 56: 
      fg[56] += a*(0.4*f[3]*g[57]+0.4472135954999579*f[5]*g[56]+0.31943828249997*f[4]*g[56]+0.5*f[0]*g[56]+0.4472135954999579*f[6]*g[34]+0.31943828249997*f[6]*g[33]+0.5000000000000001*f[2]*g[33]+0.4*f[7]*g[21]+0.4472135954999579*f[1]*g[21]+0.5*f[4]*g[8]+0.4472135954999579*f[3]*g[7]+0.5000000000000001*g[3]*f[6]); 
    case 57: 
      fg[57] += a*(0.31943828249997*f[5]*g[57]+0.4472135954999579*f[4]*g[57]+0.5*f[0]*g[57]+0.4*f[3]*g[56]+0.31943828249997*f[7]*g[34]+0.5000000000000001*f[1]*g[34]+0.4472135954999579*f[7]*g[33]+0.4*f[6]*g[21]+0.4472135954999579*f[2]*g[21]+0.4472135954999579*f[3]*g[8]+0.5*f[5]*g[7]+0.5000000000000001*g[3]*f[7]); 
    case 58: 
      fg[58] += a*(0.4472135954999579*f[5]*g[58]+0.4472135954999579*f[4]*g[58]+0.5*f[0]*g[58]+0.4472135954999579*f[7]*g[36]+0.5000000000000001*f[1]*g[36]+0.4472135954999579*f[6]*g[35]+0.5000000000000001*f[2]*g[35]+0.5*f[3]*g[18]); 
    case 59: 
      fg[59] += a*(0.4*f[3]*g[60]+0.4472135954999579*f[5]*g[59]+0.31943828249997*f[4]*g[59]+0.5*f[0]*g[59]+0.4472135954999579*f[6]*g[38]+0.31943828249997*f[6]*g[37]+0.5000000000000001*f[2]*g[37]+0.4*f[7]*g[22]+0.4472135954999579*f[1]*g[22]+0.5*f[4]*g[10]+0.4472135954999579*f[3]*g[9]+0.5000000000000001*g[4]*f[6]); 
    case 60: 
      fg[60] += a*(0.31943828249997*f[5]*g[60]+0.4472135954999579*f[4]*g[60]+0.5*f[0]*g[60]+0.4*f[3]*g[59]+0.31943828249997*f[7]*g[38]+0.5000000000000001*f[1]*g[38]+0.4472135954999579*f[7]*g[37]+0.4*f[6]*g[22]+0.4472135954999579*f[2]*g[22]+0.4472135954999579*f[3]*g[10]+0.5*f[5]*g[9]+0.5000000000000001*g[4]*f[7]); 
    case 61: 
      fg[61] += a*(0.4472135954999579*f[7]*g[88]+0.31943828249997*f[6]*g[87]+0.5000000000000001*f[2]*g[87]+0.31943828249997*f[4]*g[61]+0.5*f[0]*g[61]+0.4472135954999579*f[3]*g[51]+0.5000000000000001*f[6]*g[24]+0.4472135954999579*f[1]*g[23]+0.5*f[4]*g[11]); 
    case 62: 
      fg[62] += a*(0.31943828249997*f[7]*g[88]+0.5000000000000001*f[1]*g[88]+0.4472135954999579*f[6]*g[87]+0.31943828249997*f[5]*g[62]+0.5*f[0]*g[62]+0.4472135954999579*f[3]*g[51]+0.4472135954999579*f[2]*g[24]+0.5000000000000001*f[7]*g[23]+0.5*f[5]*g[11]); 
    case 63: 
      fg[63] += a*(0.4472135954999579*f[6]*g[89]+0.5000000000000001*f[2]*g[89]+0.5*f[3]*g[64]+0.4472135954999579*f[4]*g[63]+0.5*f[0]*g[63]+0.5000000000000001*f[1]*g[39]); 
    case 64: 
      fg[64] += a*(0.4472135954999579*f[7]*g[89]+0.5000000000000001*f[1]*g[89]+0.4472135954999579*f[5]*g[64]+0.5*f[0]*g[64]+0.5*f[3]*g[63]+0.5000000000000001*f[2]*g[39]); 
    case 65: 
      fg[65] += a*(0.4472135954999579*f[5]*g[65]+0.4472135954999579*f[4]*g[65]+0.5*f[0]*g[65]+0.4472135954999579*f[7]*g[41]+0.5000000000000001*f[1]*g[41]+0.4472135954999579*f[6]*g[40]+0.5000000000000001*f[2]*g[40]+0.5*f[3]*g[19]); 
    case 66: 
      fg[66] += a*(0.4472135954999579*f[6]*g[90]+0.5000000000000001*f[2]*g[90]+0.5*f[3]*g[67]+0.4472135954999579*f[4]*g[66]+0.5*f[0]*g[66]+0.5000000000000001*f[1]*g[42]); 
    case 67: 
      fg[67] += a*(0.4472135954999579*f[7]*g[90]+0.5000000000000001*f[1]*g[90]+0.4472135954999579*f[5]*g[67]+0.5*f[0]*g[67]+0.5*f[3]*g[66]+0.5000000000000001*f[2]*g[42]); 
    case 68: 
      fg[68] += a*(0.4*f[3]*g[69]+0.4472135954999579*f[5]*g[68]+0.31943828249997*f[4]*g[68]+0.5*f[0]*g[68]+0.4472135954999579*f[6]*g[44]+0.31943828249997*f[6]*g[43]+0.5000000000000001*f[2]*g[43]+0.4*f[7]*g[25]+0.4472135954999579*f[1]*g[25]+0.5*f[4]*g[13]+0.4472135954999579*f[3]*g[12]+0.5000000000000001*g[5]*f[6]); 
    case 69: 
      fg[69] += a*(0.31943828249997*f[5]*g[69]+0.4472135954999579*f[4]*g[69]+0.5*f[0]*g[69]+0.4*f[3]*g[68]+0.31943828249997*f[7]*g[44]+0.5000000000000001*f[1]*g[44]+0.4472135954999579*f[7]*g[43]+0.4*f[6]*g[25]+0.4472135954999579*f[2]*g[25]+0.4472135954999579*f[3]*g[13]+0.5*f[5]*g[12]+0.5000000000000001*g[5]*f[7]); 
    case 70: 
      fg[70] += a*(0.4472135954999579*f[7]*g[92]+0.31943828249997*f[6]*g[91]+0.5000000000000001*f[2]*g[91]+0.31943828249997*f[4]*g[70]+0.5*f[0]*g[70]+0.4472135954999579*f[3]*g[52]+0.5000000000000001*f[6]*g[27]+0.4472135954999579*f[1]*g[26]+0.5*f[4]*g[14]); 
    case 71: 
      fg[71] += a*(0.31943828249997*f[7]*g[92]+0.5000000000000001*f[1]*g[92]+0.4472135954999579*f[6]*g[91]+0.31943828249997*f[5]*g[71]+0.5*f[0]*g[71]+0.4472135954999579*f[3]*g[52]+0.4472135954999579*f[2]*g[27]+0.5000000000000001*f[7]*g[26]+0.5*f[5]*g[14]); 
    case 72: 
      fg[72] += a*(0.4472135954999579*f[6]*g[93]+0.5000000000000001*f[2]*g[93]+0.5*f[3]*g[73]+0.4472135954999579*f[4]*g[72]+0.5*f[0]*g[72]+0.5000000000000001*f[1]*g[45]); 
    case 73: 
      fg[73] += a*(0.4472135954999579*f[7]*g[93]+0.5000000000000001*f[1]*g[93]+0.4472135954999579*f[5]*g[73]+0.5*f[0]*g[73]+0.5*f[3]*g[72]+0.5000000000000001*f[2]*g[45]); 
    case 74: 
      fg[74] += a*(0.4472135954999579*f[7]*g[95]+0.31943828249997*f[6]*g[94]+0.5000000000000001*f[2]*g[94]+0.31943828249997*f[4]*g[74]+0.5*f[0]*g[74]+0.4472135954999579*f[3]*g[53]+0.5000000000000001*f[6]*g[29]+0.4472135954999579*f[1]*g[28]+0.5*f[4]*g[15]); 
    case 75: 
      fg[75] += a*(0.31943828249997*f[7]*g[95]+0.5000000000000001*f[1]*g[95]+0.4472135954999579*f[6]*g[94]+0.31943828249997*f[5]*g[75]+0.5*f[0]*g[75]+0.4472135954999579*f[3]*g[53]+0.4472135954999579*f[2]*g[29]+0.5000000000000001*f[7]*g[28]+0.5*f[5]*g[15]); 
    case 76: 
      fg[76] += a*(0.5*f[3]*g[109]+0.5000000000000001*f[2]*g[99]+0.5000000000000001*f[1]*g[98]+0.5*f[0]*g[76]); 
    case 77: 
      fg[77] += a*(0.4472135954999579*f[6]*g[100]+0.5000000000000001*f[2]*g[100]+0.5*f[3]*g[78]+0.4472135954999579*f[4]*g[77]+0.5*f[0]*g[77]+0.5000000000000001*f[1]*g[46]); 
    case 78: 
      fg[78] += a*(0.4472135954999579*f[7]*g[100]+0.5000000000000001*f[1]*g[100]+0.4472135954999579*f[5]*g[78]+0.5*f[0]*g[78]+0.5*f[3]*g[77]+0.5000000000000001*f[2]*g[46]); 
    case 79: 
      fg[79] += a*(0.5*f[3]*g[110]+0.5000000000000001*f[2]*g[102]+0.5000000000000001*f[1]*g[101]+0.5*f[0]*g[79]); 
    case 80: 
      fg[80] += a*(0.4472135954999579*f[5]*g[80]+0.4472135954999579*f[4]*g[80]+0.5*f[0]*g[80]+0.4472135954999579*f[7]*g[48]+0.5000000000000001*f[1]*g[48]+0.4472135954999579*f[6]*g[47]+0.5000000000000001*f[2]*g[47]+0.5*f[3]*g[20]); 
    case 81: 
      fg[81] += a*(0.4472135954999579*f[6]*g[103]+0.5000000000000001*f[2]*g[103]+0.5*f[3]*g[82]+0.4472135954999579*f[4]*g[81]+0.5*f[0]*g[81]+0.5000000000000001*f[1]*g[49]); 
    case 82: 
      fg[82] += a*(0.4472135954999579*f[7]*g[103]+0.5000000000000001*f[1]*g[103]+0.4472135954999579*f[5]*g[82]+0.5*f[0]*g[82]+0.5*f[3]*g[81]+0.5000000000000001*f[2]*g[49]); 
    case 83: 
      fg[83] += a*(0.4472135954999579*f[6]*g[104]+0.5000000000000001*f[2]*g[104]+0.5*f[3]*g[84]+0.4472135954999579*f[4]*g[83]+0.5*f[0]*g[83]+0.5000000000000001*f[1]*g[50]); 
    case 84: 
      fg[84] += a*(0.4472135954999579*f[7]*g[104]+0.5000000000000001*f[1]*g[104]+0.4472135954999579*f[5]*g[84]+0.5*f[0]*g[84]+0.5*f[3]*g[83]+0.5000000000000001*f[2]*g[50]); 
    case 85: 
      fg[85] += a*(0.5*f[3]*g[111]+0.5000000000000001*f[2]*g[106]+0.5000000000000001*f[1]*g[105]+0.5*f[0]*g[85]); 
    case 86: 
      fg[86] += a*(0.4*f[6]*g[108]+0.4472135954999579*f[2]*g[108]+0.4*f[7]*g[107]+0.4472135954999579*f[1]*g[107]+0.447213595499958*f[3]*g[97]+0.447213595499958*f[3]*g[96]+0.4472135954999579*f[5]*g[86]+0.4472135954999579*f[4]*g[86]+0.5*f[0]*g[86]+0.447213595499958*f[7]*g[55]+0.5*f[1]*g[55]+0.447213595499958*f[6]*g[54]+0.5*f[2]*g[54]+0.5*f[3]*g[30]); 
    case 87: 
      fg[87] += a*(0.4*f[3]*g[88]+0.4472135954999579*f[5]*g[87]+0.31943828249997*f[4]*g[87]+0.5*f[0]*g[87]+0.4472135954999579*f[6]*g[62]+0.31943828249997*f[6]*g[61]+0.5000000000000001*f[2]*g[61]+0.4*f[7]*g[51]+0.447213595499958*f[1]*g[51]+0.5000000000000001*f[4]*g[24]+0.447213595499958*f[3]*g[23]+0.5*f[6]*g[11]); 
    case 88: 
      fg[88] += a*(0.31943828249997*f[5]*g[88]+0.4472135954999579*f[4]*g[88]+0.5*f[0]*g[88]+0.4*f[3]*g[87]+0.31943828249997*f[7]*g[62]+0.5000000000000001*f[1]*g[62]+0.4472135954999579*f[7]*g[61]+0.4*f[6]*g[51]+0.447213595499958*f[2]*g[51]+0.447213595499958*f[3]*g[24]+0.5000000000000001*f[5]*g[23]+0.5*f[7]*g[11]); 
    case 89: 
      fg[89] += a*(0.4472135954999579*f[5]*g[89]+0.4472135954999579*f[4]*g[89]+0.5*f[0]*g[89]+0.4472135954999579*f[7]*g[64]+0.5000000000000001*f[1]*g[64]+0.4472135954999579*f[6]*g[63]+0.5000000000000001*f[2]*g[63]+0.5*f[3]*g[39]); 
    case 90: 
      fg[90] += a*(0.4472135954999579*f[5]*g[90]+0.4472135954999579*f[4]*g[90]+0.5*f[0]*g[90]+0.4472135954999579*f[7]*g[67]+0.5000000000000001*f[1]*g[67]+0.4472135954999579*f[6]*g[66]+0.5000000000000001*f[2]*g[66]+0.5*f[3]*g[42]); 
    case 91: 
      fg[91] += a*(0.4*f[3]*g[92]+0.4472135954999579*f[5]*g[91]+0.31943828249997*f[4]*g[91]+0.5*f[0]*g[91]+0.4472135954999579*f[6]*g[71]+0.31943828249997*f[6]*g[70]+0.5000000000000001*f[2]*g[70]+0.4*f[7]*g[52]+0.447213595499958*f[1]*g[52]+0.5000000000000001*f[4]*g[27]+0.447213595499958*f[3]*g[26]+0.5*f[6]*g[14]); 
    case 92: 
      fg[92] += a*(0.31943828249997*f[5]*g[92]+0.4472135954999579*f[4]*g[92]+0.5*f[0]*g[92]+0.4*f[3]*g[91]+0.31943828249997*f[7]*g[71]+0.5000000000000001*f[1]*g[71]+0.4472135954999579*f[7]*g[70]+0.4*f[6]*g[52]+0.447213595499958*f[2]*g[52]+0.447213595499958*f[3]*g[27]+0.5000000000000001*f[5]*g[26]+0.5*f[7]*g[14]); 
    case 93: 
      fg[93] += a*(0.4472135954999579*f[5]*g[93]+0.4472135954999579*f[4]*g[93]+0.5*f[0]*g[93]+0.4472135954999579*f[7]*g[73]+0.5000000000000001*f[1]*g[73]+0.4472135954999579*f[6]*g[72]+0.5000000000000001*f[2]*g[72]+0.5*f[3]*g[45]); 
    case 94: 
      fg[94] += a*(0.4*f[3]*g[95]+0.4472135954999579*f[5]*g[94]+0.31943828249997*f[4]*g[94]+0.5*f[0]*g[94]+0.4472135954999579*f[6]*g[75]+0.31943828249997*f[6]*g[74]+0.5000000000000001*f[2]*g[74]+0.4*f[7]*g[53]+0.447213595499958*f[1]*g[53]+0.5000000000000001*f[4]*g[29]+0.447213595499958*f[3]*g[28]+0.5*f[6]*g[15]); 
    case 95: 
      fg[95] += a*(0.31943828249997*f[5]*g[95]+0.4472135954999579*f[4]*g[95]+0.5*f[0]*g[95]+0.4*f[3]*g[94]+0.31943828249997*f[7]*g[75]+0.5000000000000001*f[1]*g[75]+0.4472135954999579*f[7]*g[74]+0.4*f[6]*g[53]+0.447213595499958*f[2]*g[53]+0.447213595499958*f[3]*g[29]+0.5000000000000001*f[5]*g[28]+0.5*f[7]*g[15]); 
    case 96: 
      fg[96] += a*(0.4472135954999579*f[7]*g[108]+0.31943828249997*f[6]*g[107]+0.5000000000000001*f[2]*g[107]+0.31943828249997*f[4]*g[96]+0.5*f[0]*g[96]+0.447213595499958*f[3]*g[86]+0.5*f[6]*g[55]+0.447213595499958*f[1]*g[54]+0.5000000000000001*f[4]*g[30]); 
    case 97: 
      fg[97] += a*(0.31943828249997*f[7]*g[108]+0.5000000000000001*f[1]*g[108]+0.4472135954999579*f[6]*g[107]+0.31943828249997*f[5]*g[97]+0.5*f[0]*g[97]+0.447213595499958*f[3]*g[86]+0.447213595499958*f[2]*g[55]+0.5*f[7]*g[54]+0.5000000000000001*f[5]*g[30]); 
    case 98: 
      fg[98] += a*(0.4472135954999579*f[6]*g[109]+0.5000000000000001*f[2]*g[109]+0.5*f[3]*g[99]+0.4472135954999579*f[4]*g[98]+0.5*f[0]*g[98]+0.5000000000000001*f[1]*g[76]); 
    case 99: 
      fg[99] += a*(0.4472135954999579*f[7]*g[109]+0.5000000000000001*f[1]*g[109]+0.4472135954999579*f[5]*g[99]+0.5*f[0]*g[99]+0.5*f[3]*g[98]+0.5000000000000001*f[2]*g[76]); 
    case 100: 
      fg[100] += a*(0.4472135954999579*f[5]*g[100]+0.4472135954999579*f[4]*g[100]+0.5*f[0]*g[100]+0.4472135954999579*f[7]*g[78]+0.5000000000000001*f[1]*g[78]+0.4472135954999579*f[6]*g[77]+0.5000000000000001*f[2]*g[77]+0.5*f[3]*g[46]); 
    case 101: 
      fg[101] += a*(0.4472135954999579*f[6]*g[110]+0.5000000000000001*f[2]*g[110]+0.5*f[3]*g[102]+0.4472135954999579*f[4]*g[101]+0.5*f[0]*g[101]+0.5000000000000001*f[1]*g[79]); 
    case 102: 
      fg[102] += a*(0.4472135954999579*f[7]*g[110]+0.5000000000000001*f[1]*g[110]+0.4472135954999579*f[5]*g[102]+0.5*f[0]*g[102]+0.5*f[3]*g[101]+0.5000000000000001*f[2]*g[79]); 
    case 103: 
      fg[103] += a*(0.4472135954999579*f[5]*g[103]+0.4472135954999579*f[4]*g[103]+0.5*f[0]*g[103]+0.4472135954999579*f[7]*g[82]+0.5000000000000001*f[1]*g[82]+0.4472135954999579*f[6]*g[81]+0.5000000000000001*f[2]*g[81]+0.5*f[3]*g[49]); 
    case 104: 
      fg[104] += a*(0.4472135954999579*f[5]*g[104]+0.4472135954999579*f[4]*g[104]+0.5*f[0]*g[104]+0.4472135954999579*f[7]*g[84]+0.5000000000000001*f[1]*g[84]+0.4472135954999579*f[6]*g[83]+0.5000000000000001*f[2]*g[83]+0.5*f[3]*g[50]); 
    case 105: 
      fg[105] += a*(0.4472135954999579*f[6]*g[111]+0.5000000000000001*f[2]*g[111]+0.5*f[3]*g[106]+0.4472135954999579*f[4]*g[105]+0.5*f[0]*g[105]+0.5000000000000001*f[1]*g[85]); 
    case 106: 
      fg[106] += a*(0.4472135954999579*f[7]*g[111]+0.5000000000000001*f[1]*g[111]+0.4472135954999579*f[5]*g[106]+0.5*f[0]*g[106]+0.5*f[3]*g[105]+0.5000000000000001*f[2]*g[85]); 
    case 107: 
      fg[107] += a*(0.4*f[3]*g[108]+0.4472135954999579*f[5]*g[107]+0.31943828249997*f[4]*g[107]+0.5*f[0]*g[107]+0.4472135954999579*f[6]*g[97]+0.31943828249997*f[6]*g[96]+0.5000000000000001*f[2]*g[96]+0.4*f[7]*g[86]+0.4472135954999579*f[1]*g[86]+0.5*f[4]*g[55]+0.4472135954999579*f[3]*g[54]+0.5000000000000001*f[6]*g[30]); 
    case 108: 
      fg[108] += a*(0.31943828249997*f[5]*g[108]+0.4472135954999579*f[4]*g[108]+0.5*f[0]*g[108]+0.4*f[3]*g[107]+0.31943828249997*f[7]*g[97]+0.5000000000000001*f[1]*g[97]+0.4472135954999579*f[7]*g[96]+0.4*f[6]*g[86]+0.4472135954999579*f[2]*g[86]+0.4472135954999579*f[3]*g[55]+0.5*f[5]*g[54]+0.5000000000000001*f[7]*g[30]); 
    case 109: 
      fg[109] += a*(0.4472135954999579*f[5]*g[109]+0.4472135954999579*f[4]*g[109]+0.5*f[0]*g[109]+0.4472135954999579*f[7]*g[99]+0.5000000000000001*f[1]*g[99]+0.4472135954999579*f[6]*g[98]+0.5000000000000001*f[2]*g[98]+0.5*f[3]*g[76]); 
    case 110: 
      fg[110] += a*(0.4472135954999579*f[5]*g[110]+0.4472135954999579*f[4]*g[110]+0.5*f[0]*g[110]+0.4472135954999579*f[7]*g[102]+0.5000000000000001*f[1]*g[102]+0.4472135954999579*f[6]*g[101]+0.5000000000000001*f[2]*g[101]+0.5*f[3]*g[79]); 
    case 111: 
      fg[111] += a*(0.4472135954999579*f[5]*g[111]+0.4472135954999579*f[4]*g[111]+0.5*f[0]*g[111]+0.4472135954999579*f[7]*g[106]+0.5000000000000001*f[1]*g[106]+0.4472135954999579*f[6]*g[105]+0.5000000000000001*f[2]*g[105]+0.5*f[3]*g[85]); 
  } 
} 