#include <gkyl_inflate_surf_kernels.h>

GKYL_CU_DH void inflate_surfx_3x2v_ser_p2(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[2] = 1.4142135623730951*fld_deflated[1]; 
  fld[3] = 1.4142135623730951*fld_deflated[2]; 
  fld[4] = 1.4142135623730951*fld_deflated[3]; 
  fld[5] = 1.4142135623730951*fld_deflated[4]; 
  fld[8] = 1.4142135623730951*fld_deflated[5]; 
  fld[10] = 1.4142135623730951*fld_deflated[6]; 
  fld[11] = 1.4142135623730951*fld_deflated[7]; 
  fld[13] = 1.4142135623730951*fld_deflated[8]; 
  fld[14] = 1.4142135623730951*fld_deflated[9]; 
  fld[15] = 1.4142135623730951*fld_deflated[10]; 
  fld[17] = 1.4142135623730951*fld_deflated[11]; 
  fld[18] = 1.4142135623730951*fld_deflated[12]; 
  fld[19] = 1.4142135623730951*fld_deflated[13]; 
  fld[20] = 1.4142135623730951*fld_deflated[14]; 
  fld[24] = 1.4142135623730951*fld_deflated[15]; 
  fld[27] = 1.4142135623730951*fld_deflated[16]; 
  fld[29] = 1.4142135623730951*fld_deflated[17]; 
  fld[30] = 1.4142135623730951*fld_deflated[18]; 
  fld[34] = 1.4142135623730951*fld_deflated[19]; 
  fld[36] = 1.4142135623730951*fld_deflated[20]; 
  fld[38] = 1.4142135623730951*fld_deflated[21]; 
  fld[39] = 1.4142135623730951*fld_deflated[22]; 
  fld[41] = 1.4142135623730951*fld_deflated[23]; 
  fld[42] = 1.4142135623730951*fld_deflated[24]; 
  fld[44] = 1.4142135623730951*fld_deflated[25]; 
  fld[45] = 1.4142135623730951*fld_deflated[26]; 
  fld[46] = 1.4142135623730951*fld_deflated[27]; 
  fld[48] = 1.4142135623730951*fld_deflated[28]; 
  fld[49] = 1.4142135623730951*fld_deflated[29]; 
  fld[50] = 1.4142135623730951*fld_deflated[30]; 
  fld[55] = 1.4142135623730951*fld_deflated[31]; 
  fld[62] = 1.4142135623730951*fld_deflated[32]; 
  fld[64] = 1.4142135623730951*fld_deflated[33]; 
  fld[67] = 1.4142135623730951*fld_deflated[34]; 
  fld[71] = 1.4142135623730951*fld_deflated[35]; 
  fld[73] = 1.4142135623730951*fld_deflated[36]; 
  fld[75] = 1.4142135623730951*fld_deflated[37]; 
  fld[76] = 1.4142135623730951*fld_deflated[38]; 
  fld[78] = 1.4142135623730951*fld_deflated[39]; 
  fld[79] = 1.4142135623730951*fld_deflated[40]; 
  fld[82] = 1.4142135623730951*fld_deflated[41]; 
  fld[84] = 1.4142135623730951*fld_deflated[42]; 
  fld[85] = 1.4142135623730951*fld_deflated[43]; 
  fld[97] = 1.4142135623730951*fld_deflated[44]; 
  fld[99] = 1.4142135623730951*fld_deflated[45]; 
  fld[102] = 1.4142135623730951*fld_deflated[46]; 
  fld[106] = 1.4142135623730951*fld_deflated[47]; 
}
 
GKYL_CU_DH void inflate_surfy_3x2v_ser_p2(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[1] = 1.4142135623730951*fld_deflated[1]; 
  fld[3] = 1.4142135623730951*fld_deflated[2]; 
  fld[4] = 1.4142135623730951*fld_deflated[3]; 
  fld[5] = 1.4142135623730951*fld_deflated[4]; 
  fld[7] = 1.4142135623730951*fld_deflated[5]; 
  fld[9] = 1.4142135623730951*fld_deflated[6]; 
  fld[11] = 1.4142135623730951*fld_deflated[7]; 
  fld[12] = 1.4142135623730951*fld_deflated[8]; 
  fld[14] = 1.4142135623730951*fld_deflated[9]; 
  fld[15] = 1.4142135623730951*fld_deflated[10]; 
  fld[16] = 1.4142135623730951*fld_deflated[11]; 
  fld[18] = 1.4142135623730951*fld_deflated[12]; 
  fld[19] = 1.4142135623730951*fld_deflated[13]; 
  fld[20] = 1.4142135623730951*fld_deflated[14]; 
  fld[23] = 1.4142135623730951*fld_deflated[15]; 
  fld[26] = 1.4142135623730951*fld_deflated[16]; 
  fld[28] = 1.4142135623730951*fld_deflated[17]; 
  fld[30] = 1.4142135623730951*fld_deflated[18]; 
  fld[33] = 1.4142135623730951*fld_deflated[19]; 
  fld[35] = 1.4142135623730951*fld_deflated[20]; 
  fld[37] = 1.4142135623730951*fld_deflated[21]; 
  fld[39] = 1.4142135623730951*fld_deflated[22]; 
  fld[40] = 1.4142135623730951*fld_deflated[23]; 
  fld[42] = 1.4142135623730951*fld_deflated[24]; 
  fld[43] = 1.4142135623730951*fld_deflated[25]; 
  fld[45] = 1.4142135623730951*fld_deflated[26]; 
  fld[46] = 1.4142135623730951*fld_deflated[27]; 
  fld[47] = 1.4142135623730951*fld_deflated[28]; 
  fld[49] = 1.4142135623730951*fld_deflated[29]; 
  fld[50] = 1.4142135623730951*fld_deflated[30]; 
  fld[54] = 1.4142135623730951*fld_deflated[31]; 
  fld[61] = 1.4142135623730951*fld_deflated[32]; 
  fld[63] = 1.4142135623730951*fld_deflated[33]; 
  fld[66] = 1.4142135623730951*fld_deflated[34]; 
  fld[70] = 1.4142135623730951*fld_deflated[35]; 
  fld[72] = 1.4142135623730951*fld_deflated[36]; 
  fld[74] = 1.4142135623730951*fld_deflated[37]; 
  fld[76] = 1.4142135623730951*fld_deflated[38]; 
  fld[77] = 1.4142135623730951*fld_deflated[39]; 
  fld[79] = 1.4142135623730951*fld_deflated[40]; 
  fld[81] = 1.4142135623730951*fld_deflated[41]; 
  fld[83] = 1.4142135623730951*fld_deflated[42]; 
  fld[85] = 1.4142135623730951*fld_deflated[43]; 
  fld[96] = 1.4142135623730951*fld_deflated[44]; 
  fld[98] = 1.4142135623730951*fld_deflated[45]; 
  fld[101] = 1.4142135623730951*fld_deflated[46]; 
  fld[105] = 1.4142135623730951*fld_deflated[47]; 
}
 
GKYL_CU_DH void inflate_surfz_3x2v_ser_p2(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[1] = 1.4142135623730951*fld_deflated[1]; 
  fld[2] = 1.4142135623730951*fld_deflated[2]; 
  fld[4] = 1.4142135623730951*fld_deflated[3]; 
  fld[5] = 1.4142135623730951*fld_deflated[4]; 
  fld[6] = 1.4142135623730951*fld_deflated[5]; 
  fld[9] = 1.4142135623730951*fld_deflated[6]; 
  fld[10] = 1.4142135623730951*fld_deflated[7]; 
  fld[12] = 1.4142135623730951*fld_deflated[8]; 
  fld[13] = 1.4142135623730951*fld_deflated[9]; 
  fld[15] = 1.4142135623730951*fld_deflated[10]; 
  fld[16] = 1.4142135623730951*fld_deflated[11]; 
  fld[17] = 1.4142135623730951*fld_deflated[12]; 
  fld[19] = 1.4142135623730951*fld_deflated[13]; 
  fld[20] = 1.4142135623730951*fld_deflated[14]; 
  fld[22] = 1.4142135623730951*fld_deflated[15]; 
  fld[25] = 1.4142135623730951*fld_deflated[16]; 
  fld[28] = 1.4142135623730951*fld_deflated[17]; 
  fld[29] = 1.4142135623730951*fld_deflated[18]; 
  fld[31] = 1.4142135623730951*fld_deflated[19]; 
  fld[32] = 1.4142135623730951*fld_deflated[20]; 
  fld[37] = 1.4142135623730951*fld_deflated[21]; 
  fld[38] = 1.4142135623730951*fld_deflated[22]; 
  fld[40] = 1.4142135623730951*fld_deflated[23]; 
  fld[41] = 1.4142135623730951*fld_deflated[24]; 
  fld[43] = 1.4142135623730951*fld_deflated[25]; 
  fld[44] = 1.4142135623730951*fld_deflated[26]; 
  fld[46] = 1.4142135623730951*fld_deflated[27]; 
  fld[47] = 1.4142135623730951*fld_deflated[28]; 
  fld[48] = 1.4142135623730951*fld_deflated[29]; 
  fld[50] = 1.4142135623730951*fld_deflated[30]; 
  fld[53] = 1.4142135623730951*fld_deflated[31]; 
  fld[59] = 1.4142135623730951*fld_deflated[32]; 
  fld[60] = 1.4142135623730951*fld_deflated[33]; 
  fld[65] = 1.4142135623730951*fld_deflated[34]; 
  fld[68] = 1.4142135623730951*fld_deflated[35]; 
  fld[69] = 1.4142135623730951*fld_deflated[36]; 
  fld[74] = 1.4142135623730951*fld_deflated[37]; 
  fld[75] = 1.4142135623730951*fld_deflated[38]; 
  fld[77] = 1.4142135623730951*fld_deflated[39]; 
  fld[78] = 1.4142135623730951*fld_deflated[40]; 
  fld[80] = 1.4142135623730951*fld_deflated[41]; 
  fld[83] = 1.4142135623730951*fld_deflated[42]; 
  fld[84] = 1.4142135623730951*fld_deflated[43]; 
  fld[94] = 1.4142135623730951*fld_deflated[44]; 
  fld[95] = 1.4142135623730951*fld_deflated[45]; 
  fld[100] = 1.4142135623730951*fld_deflated[46]; 
  fld[104] = 1.4142135623730951*fld_deflated[47]; 
}
 
