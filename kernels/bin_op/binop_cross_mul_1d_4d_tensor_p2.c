#include <gkyl_binop_cross_mul_tensor.h> 
 
GKYL_CU_DH
void
binop_cross_mul_1d_4d_tensor_p2(const double *f, const double *g, double *fg) 
{ 
  // f:  First input DG field. 
  // g:  Second input DG field. 
  // fg: Output DG field f*g using weak multiplication. 
 
  double tmp[81] = {0.0}; 
  tmp[0] = 0.7071067811865475*f[2]*g[11]+0.7071067811865475*f[1]*g[1]+0.7071067811865475*f[0]*g[0]; 
  tmp[1] = 0.6324555320336759*f[1]*g[11]+0.6324555320336759*g[1]*f[2]+0.7071067811865475*f[0]*g[1]+0.7071067811865475*g[0]*f[1]; 
  tmp[2] = 0.7071067811865475*f[2]*g[19]+0.7071067811865475*f[1]*g[5]+0.7071067811865475*f[0]*g[2]; 
  tmp[3] = 0.7071067811865475*f[2]*g[21]+0.7071067811865475*f[1]*g[6]+0.7071067811865475*f[0]*g[3]; 
  tmp[4] = 0.7071067811865475*f[2]*g[25]+0.7071067811865475*f[1]*g[8]+0.7071067811865475*f[0]*g[4]; 
  tmp[5] = 0.632455532033676*f[1]*g[19]+0.6324555320336759*f[2]*g[5]+0.7071067811865475*f[0]*g[5]+0.7071067811865475*f[1]*g[2]; 
  tmp[6] = 0.632455532033676*f[1]*g[21]+0.6324555320336759*f[2]*g[6]+0.7071067811865475*f[0]*g[6]+0.7071067811865475*f[1]*g[3]; 
  tmp[7] = 0.7071067811865475*f[2]*g[32]+0.7071067811865475*f[1]*g[15]+0.7071067811865475*f[0]*g[7]; 
  tmp[8] = 0.632455532033676*f[1]*g[25]+0.6324555320336759*f[2]*g[8]+0.7071067811865475*f[0]*g[8]+0.7071067811865475*f[1]*g[4]; 
  tmp[9] = 0.7071067811865475*f[2]*g[35]+0.7071067811865475*f[1]*g[16]+0.7071067811865475*f[0]*g[9]; 
  tmp[10] = 0.7071067811865475*f[2]*g[37]+0.7071067811865475*f[1]*g[17]+0.7071067811865475*f[0]*g[10]; 
  tmp[11] = 0.4517539514526256*f[2]*g[11]+0.7071067811865475*f[0]*g[11]+0.7071067811865475*g[0]*f[2]+0.6324555320336759*f[1]*g[1]; 
  tmp[12] = 0.7071067811865475*f[2]*g[44]+0.7071067811865475*f[1]*g[20]+0.7071067811865475*f[0]*g[12]; 
  tmp[13] = 0.7071067811865475*f[2]*g[45]+0.7071067811865475*f[1]*g[23]+0.7071067811865475*f[0]*g[13]; 
  tmp[14] = 0.7071067811865475*f[2]*g[47]+0.7071067811865475*f[1]*g[28]+0.7071067811865475*f[0]*g[14]; 
  tmp[15] = 0.6324555320336759*f[1]*g[32]+0.6324555320336759*f[2]*g[15]+0.7071067811865475*f[0]*g[15]+0.7071067811865475*f[1]*g[7]; 
  tmp[16] = 0.6324555320336759*f[1]*g[35]+0.6324555320336759*f[2]*g[16]+0.7071067811865475*f[0]*g[16]+0.7071067811865475*f[1]*g[9]; 
  tmp[17] = 0.6324555320336759*f[1]*g[37]+0.6324555320336759*f[2]*g[17]+0.7071067811865475*f[0]*g[17]+0.7071067811865475*f[1]*g[10]; 
  tmp[18] = 0.7071067811865475*f[2]*g[50]+0.7071067811865475*f[1]*g[31]+0.7071067811865475*f[0]*g[18]; 
  tmp[19] = 0.4517539514526256*f[2]*g[19]+0.7071067811865475*f[0]*g[19]+0.632455532033676*f[1]*g[5]+0.7071067811865475*f[2]*g[2]; 
  tmp[20] = 0.632455532033676*f[1]*g[44]+0.6324555320336759*f[2]*g[20]+0.7071067811865475*f[0]*g[20]+0.7071067811865475*f[1]*g[12]; 
  tmp[21] = 0.4517539514526256*f[2]*g[21]+0.7071067811865475*f[0]*g[21]+0.632455532033676*f[1]*g[6]+0.7071067811865475*f[2]*g[3]; 
  tmp[22] = 0.7071067811865475*f[2]*g[54]+0.7071067811865475*f[1]*g[33]+0.7071067811865475*f[0]*g[22]; 
  tmp[23] = 0.632455532033676*f[1]*g[45]+0.6324555320336759*f[2]*g[23]+0.7071067811865475*f[0]*g[23]+0.7071067811865475*f[1]*g[13]; 
  tmp[24] = 0.7071067811865475*f[2]*g[55]+0.7071067811865475*f[1]*g[34]+0.7071067811865475*f[0]*g[24]; 
  tmp[25] = 0.4517539514526256*f[2]*g[25]+0.7071067811865475*f[0]*g[25]+0.632455532033676*f[1]*g[8]+0.7071067811865475*f[2]*g[4]; 
  tmp[26] = 0.7071067811865475*f[2]*g[57]+0.7071067811865475*f[1]*g[36]+0.7071067811865475*f[0]*g[26]; 
  tmp[27] = 0.7071067811865475*f[2]*g[58]+0.7071067811865475*f[1]*g[39]+0.7071067811865475*f[0]*g[27]; 
  tmp[28] = 0.632455532033676*f[1]*g[47]+0.6324555320336759*f[2]*g[28]+0.7071067811865475*f[0]*g[28]+0.7071067811865475*f[1]*g[14]; 
  tmp[29] = 0.7071067811865475*f[2]*g[60]+0.7071067811865475*f[1]*g[41]+0.7071067811865475*f[0]*g[29]; 
  tmp[30] = 0.7071067811865475*f[2]*g[62]+0.7071067811865475*f[1]*g[42]+0.7071067811865475*f[0]*g[30]; 
  tmp[31] = 0.632455532033676*f[1]*g[50]+0.6324555320336759*f[2]*g[31]+0.7071067811865475*f[0]*g[31]+0.7071067811865475*f[1]*g[18]; 
  tmp[32] = 0.4517539514526256*f[2]*g[32]+0.7071067811865475*f[0]*g[32]+0.6324555320336759*f[1]*g[15]+0.7071067811865475*f[2]*g[7]; 
  tmp[33] = 0.6324555320336759*f[1]*g[54]+0.6324555320336759*f[2]*g[33]+0.7071067811865475*f[0]*g[33]+0.7071067811865475*f[1]*g[22]; 
  tmp[34] = 0.6324555320336759*f[1]*g[55]+0.6324555320336759*f[2]*g[34]+0.7071067811865475*f[0]*g[34]+0.7071067811865475*f[1]*g[24]; 
  tmp[35] = 0.4517539514526256*f[2]*g[35]+0.7071067811865475*f[0]*g[35]+0.6324555320336759*f[1]*g[16]+0.7071067811865475*f[2]*g[9]; 
  tmp[36] = 0.6324555320336759*f[1]*g[57]+0.6324555320336759*f[2]*g[36]+0.7071067811865475*f[0]*g[36]+0.7071067811865475*f[1]*g[26]; 
  tmp[37] = 0.4517539514526256*f[2]*g[37]+0.7071067811865475*f[0]*g[37]+0.6324555320336759*f[1]*g[17]+0.7071067811865475*f[2]*g[10]; 
  tmp[38] = 0.7071067811865475*f[2]*g[66]+0.7071067811865475*f[1]*g[51]+0.7071067811865475*f[0]*g[38]; 
  tmp[39] = 0.6324555320336759*f[1]*g[58]+0.6324555320336759*f[2]*g[39]+0.7071067811865475*f[0]*g[39]+0.7071067811865475*f[1]*g[27]; 
  tmp[40] = 0.7071067811865475*f[2]*g[67]+0.7071067811865475*f[1]*g[52]+0.7071067811865475*f[0]*g[40]; 
  tmp[41] = 0.6324555320336759*f[1]*g[60]+0.6324555320336759*f[2]*g[41]+0.7071067811865475*f[0]*g[41]+0.7071067811865475*f[1]*g[29]; 
  tmp[42] = 0.6324555320336759*f[1]*g[62]+0.6324555320336759*f[2]*g[42]+0.7071067811865475*f[0]*g[42]+0.7071067811865475*f[1]*g[30]; 
  tmp[43] = 0.7071067811865475*f[2]*g[69]+0.7071067811865475*f[1]*g[53]+0.7071067811865475*f[0]*g[43]; 
  tmp[44] = 0.4517539514526256*f[2]*g[44]+0.7071067811865475*f[0]*g[44]+0.632455532033676*f[1]*g[20]+0.7071067811865475*f[2]*g[12]; 
  tmp[45] = 0.4517539514526256*f[2]*g[45]+0.7071067811865475*f[0]*g[45]+0.632455532033676*f[1]*g[23]+0.7071067811865475*f[2]*g[13]; 
  tmp[46] = 0.7071067811865475*f[2]*g[72]+0.7071067811865475*f[1]*g[56]+0.7071067811865475*f[0]*g[46]; 
  tmp[47] = 0.4517539514526256*f[2]*g[47]+0.7071067811865475*f[0]*g[47]+0.632455532033676*f[1]*g[28]+0.7071067811865475*f[2]*g[14]; 
  tmp[48] = 0.7071067811865475*f[2]*g[73]+0.7071067811865475*f[1]*g[61]+0.7071067811865475*f[0]*g[48]; 
  tmp[49] = 0.7071067811865475*f[2]*g[74]+0.7071067811865475*f[1]*g[64]+0.7071067811865475*f[0]*g[49]; 
  tmp[50] = 0.4517539514526256*f[2]*g[50]+0.7071067811865475*f[0]*g[50]+0.632455532033676*f[1]*g[31]+0.7071067811865475*f[2]*g[18]; 
  tmp[51] = 0.632455532033676*f[1]*g[66]+0.6324555320336759*f[2]*g[51]+0.7071067811865475*f[0]*g[51]+0.7071067811865475*f[1]*g[38]; 
  tmp[52] = 0.632455532033676*f[1]*g[67]+0.6324555320336759*f[2]*g[52]+0.7071067811865475*f[0]*g[52]+0.7071067811865475*f[1]*g[40]; 
  tmp[53] = 0.632455532033676*f[1]*g[69]+0.6324555320336759*f[2]*g[53]+0.7071067811865475*f[0]*g[53]+0.7071067811865475*f[1]*g[43]; 
  tmp[54] = 0.4517539514526256*f[2]*g[54]+0.7071067811865475*f[0]*g[54]+0.6324555320336759*f[1]*g[33]+0.7071067811865475*f[2]*g[22]; 
  tmp[55] = 0.4517539514526256*f[2]*g[55]+0.7071067811865475*f[0]*g[55]+0.6324555320336759*f[1]*g[34]+0.7071067811865475*f[2]*g[24]; 
  tmp[56] = 0.6324555320336759*f[1]*g[72]+0.6324555320336759*f[2]*g[56]+0.7071067811865475*f[0]*g[56]+0.7071067811865475*f[1]*g[46]; 
  tmp[57] = 0.4517539514526256*f[2]*g[57]+0.7071067811865475*f[0]*g[57]+0.6324555320336759*f[1]*g[36]+0.7071067811865475*f[2]*g[26]; 
  tmp[58] = 0.4517539514526256*f[2]*g[58]+0.7071067811865475*f[0]*g[58]+0.6324555320336759*f[1]*g[39]+0.7071067811865475*f[2]*g[27]; 
  tmp[59] = 0.7071067811865475*f[2]*g[76]+0.7071067811865475*f[1]*g[68]+0.7071067811865475*f[0]*g[59]; 
  tmp[60] = 0.4517539514526256*f[2]*g[60]+0.7071067811865475*f[0]*g[60]+0.6324555320336759*f[1]*g[41]+0.7071067811865475*f[2]*g[29]; 
  tmp[61] = 0.6324555320336759*f[1]*g[73]+0.6324555320336759*f[2]*g[61]+0.7071067811865475*f[0]*g[61]+0.7071067811865475*f[1]*g[48]; 
  tmp[62] = 0.4517539514526256*f[2]*g[62]+0.7071067811865475*f[0]*g[62]+0.6324555320336759*f[1]*g[42]+0.7071067811865475*f[2]*g[30]; 
  tmp[63] = 0.7071067811865475*f[2]*g[77]+0.7071067811865475*f[1]*g[70]+0.7071067811865475*f[0]*g[63]; 
  tmp[64] = 0.6324555320336759*f[1]*g[74]+0.6324555320336759*f[2]*g[64]+0.7071067811865475*f[0]*g[64]+0.7071067811865475*f[1]*g[49]; 
  tmp[65] = 0.7071067811865475*f[2]*g[78]+0.7071067811865475*f[1]*g[71]+0.7071067811865475*f[0]*g[65]; 
  tmp[66] = 0.4517539514526256*f[2]*g[66]+0.7071067811865475*f[0]*g[66]+0.632455532033676*f[1]*g[51]+0.7071067811865475*f[2]*g[38]; 
  tmp[67] = 0.4517539514526256*f[2]*g[67]+0.7071067811865475*f[0]*g[67]+0.632455532033676*f[1]*g[52]+0.7071067811865475*f[2]*g[40]; 
  tmp[68] = 0.632455532033676*f[1]*g[76]+0.6324555320336759*f[2]*g[68]+0.7071067811865475*f[0]*g[68]+0.7071067811865475*f[1]*g[59]; 
  tmp[69] = 0.4517539514526256*f[2]*g[69]+0.7071067811865475*f[0]*g[69]+0.632455532033676*f[1]*g[53]+0.7071067811865475*f[2]*g[43]; 
  tmp[70] = 0.632455532033676*f[1]*g[77]+0.6324555320336759*f[2]*g[70]+0.7071067811865475*f[0]*g[70]+0.7071067811865475*f[1]*g[63]; 
  tmp[71] = 0.632455532033676*f[1]*g[78]+0.6324555320336759*f[2]*g[71]+0.7071067811865475*f[0]*g[71]+0.7071067811865475*f[1]*g[65]; 
  tmp[72] = 0.4517539514526256*f[2]*g[72]+0.7071067811865475*f[0]*g[72]+0.6324555320336759*f[1]*g[56]+0.7071067811865475*f[2]*g[46]; 
  tmp[73] = 0.4517539514526256*f[2]*g[73]+0.7071067811865475*f[0]*g[73]+0.6324555320336759*f[1]*g[61]+0.7071067811865475*f[2]*g[48]; 
  tmp[74] = 0.4517539514526256*f[2]*g[74]+0.7071067811865475*f[0]*g[74]+0.6324555320336759*f[1]*g[64]+0.7071067811865475*f[2]*g[49]; 
  tmp[75] = 0.7071067811865475*f[2]*g[80]+0.7071067811865475*f[1]*g[79]+0.7071067811865475*f[0]*g[75]; 
  tmp[76] = 0.4517539514526256*f[2]*g[76]+0.7071067811865475*f[0]*g[76]+0.632455532033676*f[1]*g[68]+0.7071067811865475*f[2]*g[59]; 
  tmp[77] = 0.4517539514526256*f[2]*g[77]+0.7071067811865475*f[0]*g[77]+0.632455532033676*f[1]*g[70]+0.7071067811865475*f[2]*g[63]; 
  tmp[78] = 0.4517539514526256*f[2]*g[78]+0.7071067811865475*f[0]*g[78]+0.632455532033676*f[1]*g[71]+0.7071067811865475*f[2]*g[65]; 
  tmp[79] = 0.632455532033676*f[1]*g[80]+0.6324555320336759*f[2]*g[79]+0.7071067811865475*f[0]*g[79]+0.7071067811865475*f[1]*g[75]; 
  tmp[80] = 0.4517539514526256*f[2]*g[80]+0.7071067811865475*f[0]*g[80]+0.632455532033676*f[1]*g[79]+0.7071067811865475*f[2]*g[75]; 
 
  fg[0] = tmp[0]; 
  fg[1] = tmp[1]; 
  fg[2] = tmp[2]; 
  fg[3] = tmp[3]; 
  fg[4] = tmp[4]; 
  fg[5] = tmp[5]; 
  fg[6] = tmp[6]; 
  fg[7] = tmp[7]; 
  fg[8] = tmp[8]; 
  fg[9] = tmp[9]; 
  fg[10] = tmp[10]; 
  fg[11] = tmp[11]; 
  fg[12] = tmp[12]; 
  fg[13] = tmp[13]; 
  fg[14] = tmp[14]; 
  fg[15] = tmp[15]; 
  fg[16] = tmp[16]; 
  fg[17] = tmp[17]; 
  fg[18] = tmp[18]; 
  fg[19] = tmp[19]; 
  fg[20] = tmp[20]; 
  fg[21] = tmp[21]; 
  fg[22] = tmp[22]; 
  fg[23] = tmp[23]; 
  fg[24] = tmp[24]; 
  fg[25] = tmp[25]; 
  fg[26] = tmp[26]; 
  fg[27] = tmp[27]; 
  fg[28] = tmp[28]; 
  fg[29] = tmp[29]; 
  fg[30] = tmp[30]; 
  fg[31] = tmp[31]; 
  fg[32] = tmp[32]; 
  fg[33] = tmp[33]; 
  fg[34] = tmp[34]; 
  fg[35] = tmp[35]; 
  fg[36] = tmp[36]; 
  fg[37] = tmp[37]; 
  fg[38] = tmp[38]; 
  fg[39] = tmp[39]; 
  fg[40] = tmp[40]; 
  fg[41] = tmp[41]; 
  fg[42] = tmp[42]; 
  fg[43] = tmp[43]; 
  fg[44] = tmp[44]; 
  fg[45] = tmp[45]; 
  fg[46] = tmp[46]; 
  fg[47] = tmp[47]; 
  fg[48] = tmp[48]; 
  fg[49] = tmp[49]; 
  fg[50] = tmp[50]; 
  fg[51] = tmp[51]; 
  fg[52] = tmp[52]; 
  fg[53] = tmp[53]; 
  fg[54] = tmp[54]; 
  fg[55] = tmp[55]; 
  fg[56] = tmp[56]; 
  fg[57] = tmp[57]; 
  fg[58] = tmp[58]; 
  fg[59] = tmp[59]; 
  fg[60] = tmp[60]; 
  fg[61] = tmp[61]; 
  fg[62] = tmp[62]; 
  fg[63] = tmp[63]; 
  fg[64] = tmp[64]; 
  fg[65] = tmp[65]; 
  fg[66] = tmp[66]; 
  fg[67] = tmp[67]; 
  fg[68] = tmp[68]; 
  fg[69] = tmp[69]; 
  fg[70] = tmp[70]; 
  fg[71] = tmp[71]; 
  fg[72] = tmp[72]; 
  fg[73] = tmp[73]; 
  fg[74] = tmp[74]; 
  fg[75] = tmp[75]; 
  fg[76] = tmp[76]; 
  fg[77] = tmp[77]; 
  fg[78] = tmp[78]; 
  fg[79] = tmp[79]; 
  fg[80] = tmp[80]; 
} 