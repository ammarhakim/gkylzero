#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_vmap_M1i_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *jacob_vel_inv, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *jacob_vel_inv1 = &jacob_vel_inv[3]; 
  const double dv12 = 2.0/dxv[3]; 
  const double *jacob_vel_inv2 = &jacob_vel_inv[6]; 
  double p0_over_gamma[27] = {0.0}; 
  p0_over_gamma[0] = 2.738612787525831*jacob_vel_inv0[1]*gamma[7]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[1]*dv10; 
  p0_over_gamma[1] = 2.449489742783178*jacob_vel_inv0[2]*gamma[7]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[7]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[1]*dv10; 
  p0_over_gamma[2] = 2.738612787525831*jacob_vel_inv0[1]*gamma[11]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[4]*dv10; 
  p0_over_gamma[3] = 2.738612787525831*jacob_vel_inv0[1]*gamma[13]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[5]*dv10; 
  p0_over_gamma[4] = 2.449489742783178*jacob_vel_inv0[2]*gamma[11]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[11]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[4]*dv10; 
  p0_over_gamma[5] = 2.449489742783178*jacob_vel_inv0[2]*gamma[13]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[13]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[5]*dv10; 
  p0_over_gamma[6] = 2.738612787525831*jacob_vel_inv0[1]*gamma[17]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[10]*dv10; 
  p0_over_gamma[7] = 2.449489742783178*jacob_vel_inv0[1]*gamma[7]*dv10+1.224744871391589*gamma[1]*jacob_vel_inv0[2]*dv10; 
  p0_over_gamma[8] = 2.738612787525831*jacob_vel_inv0[1]*gamma[20]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[12]*dv10; 
  p0_over_gamma[9] = 2.738612787525831*jacob_vel_inv0[1]*gamma[21]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[15]*dv10; 
  p0_over_gamma[10] = 2.449489742783178*jacob_vel_inv0[2]*gamma[17]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[17]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[10]*dv10; 
  p0_over_gamma[11] = 2.449489742783178*jacob_vel_inv0[1]*gamma[11]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[4]*dv10; 
  p0_over_gamma[12] = 2.449489742783178*jacob_vel_inv0[2]*gamma[20]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[20]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[12]*dv10; 
  p0_over_gamma[13] = 2.449489742783178*jacob_vel_inv0[1]*gamma[13]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[5]*dv10; 
  p0_over_gamma[14] = 2.738612787525831*jacob_vel_inv0[1]*gamma[23]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[18]*dv10; 
  p0_over_gamma[15] = 2.449489742783178*jacob_vel_inv0[2]*gamma[21]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[21]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[15]*dv10; 
  p0_over_gamma[16] = 2.738612787525831*jacob_vel_inv0[1]*gamma[24]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[19]*dv10; 
  p0_over_gamma[17] = 2.449489742783178*jacob_vel_inv0[1]*gamma[17]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[10]*dv10; 
  p0_over_gamma[18] = 2.449489742783178*jacob_vel_inv0[2]*gamma[23]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[23]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[18]*dv10; 
  p0_over_gamma[19] = 2.449489742783178*jacob_vel_inv0[2]*gamma[24]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[24]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[19]*dv10; 
  p0_over_gamma[20] = 2.449489742783178*jacob_vel_inv0[1]*gamma[20]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[12]*dv10; 
  p0_over_gamma[21] = 2.449489742783178*jacob_vel_inv0[1]*gamma[21]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[15]*dv10; 
  p0_over_gamma[22] = 2.738612787525831*jacob_vel_inv0[1]*gamma[26]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[25]*dv10; 
  p0_over_gamma[23] = 2.449489742783178*jacob_vel_inv0[1]*gamma[23]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[18]*dv10; 
  p0_over_gamma[24] = 2.449489742783178*jacob_vel_inv0[1]*gamma[24]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[19]*dv10; 
  p0_over_gamma[25] = 2.449489742783178*jacob_vel_inv0[2]*gamma[26]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[26]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[25]*dv10; 
  p0_over_gamma[26] = 2.449489742783178*jacob_vel_inv0[1]*gamma[26]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[25]*dv10; 

  double p1_over_gamma[27] = {0.0}; 
  p1_over_gamma[0] = 2.738612787525831*jacob_vel_inv1[1]*gamma[8]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[2]*dv11; 
  p1_over_gamma[1] = 2.738612787525831*jacob_vel_inv1[1]*gamma[12]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[4]*dv11; 
  p1_over_gamma[2] = 2.449489742783178*jacob_vel_inv1[2]*gamma[8]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[8]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[2]*dv11; 
  p1_over_gamma[3] = 2.738612787525831*jacob_vel_inv1[1]*gamma[14]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[6]*dv11; 
  p1_over_gamma[4] = 2.449489742783178*jacob_vel_inv1[2]*gamma[12]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[12]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[4]*dv11; 
  p1_over_gamma[5] = 2.738612787525831*jacob_vel_inv1[1]*gamma[18]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[10]*dv11; 
  p1_over_gamma[6] = 2.449489742783178*jacob_vel_inv1[2]*gamma[14]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[14]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[6]*dv11; 
  p1_over_gamma[7] = 2.738612787525831*jacob_vel_inv1[1]*gamma[20]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[11]*dv11; 
  p1_over_gamma[8] = 2.449489742783178*jacob_vel_inv1[1]*gamma[8]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[2]*dv11; 
  p1_over_gamma[9] = 2.738612787525831*jacob_vel_inv1[1]*gamma[22]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[16]*dv11; 
  p1_over_gamma[10] = 2.449489742783178*jacob_vel_inv1[2]*gamma[18]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[18]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[10]*dv11; 
  p1_over_gamma[11] = 2.449489742783178*jacob_vel_inv1[2]*gamma[20]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[20]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[11]*dv11; 
  p1_over_gamma[12] = 2.449489742783178*jacob_vel_inv1[1]*gamma[12]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[4]*dv11; 
  p1_over_gamma[13] = 2.738612787525831*jacob_vel_inv1[1]*gamma[23]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[17]*dv11; 
  p1_over_gamma[14] = 2.449489742783178*jacob_vel_inv1[1]*gamma[14]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[6]*dv11; 
  p1_over_gamma[15] = 2.738612787525831*jacob_vel_inv1[1]*gamma[25]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[19]*dv11; 
  p1_over_gamma[16] = 2.449489742783178*jacob_vel_inv1[2]*gamma[22]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[22]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[16]*dv11; 
  p1_over_gamma[17] = 2.449489742783178*jacob_vel_inv1[2]*gamma[23]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[23]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[17]*dv11; 
  p1_over_gamma[18] = 2.449489742783178*jacob_vel_inv1[1]*gamma[18]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[10]*dv11; 
  p1_over_gamma[19] = 2.449489742783178*jacob_vel_inv1[2]*gamma[25]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[25]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[19]*dv11; 
  p1_over_gamma[20] = 2.449489742783178*jacob_vel_inv1[1]*gamma[20]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[11]*dv11; 
  p1_over_gamma[21] = 2.738612787525831*jacob_vel_inv1[1]*gamma[26]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[24]*dv11; 
  p1_over_gamma[22] = 2.449489742783178*jacob_vel_inv1[1]*gamma[22]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[16]*dv11; 
  p1_over_gamma[23] = 2.449489742783178*jacob_vel_inv1[1]*gamma[23]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[17]*dv11; 
  p1_over_gamma[24] = 2.449489742783178*jacob_vel_inv1[2]*gamma[26]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[26]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[24]*dv11; 
  p1_over_gamma[25] = 2.449489742783178*jacob_vel_inv1[1]*gamma[25]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[19]*dv11; 
  p1_over_gamma[26] = 2.449489742783178*jacob_vel_inv1[1]*gamma[26]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[24]*dv11; 

  double p2_over_gamma[27] = {0.0}; 
  p2_over_gamma[0] = 2.738612787525831*jacob_vel_inv2[1]*gamma[9]*dv12+1.224744871391589*jacob_vel_inv2[0]*gamma[3]*dv12; 
  p2_over_gamma[1] = 2.738612787525831*jacob_vel_inv2[1]*gamma[15]*dv12+1.224744871391589*jacob_vel_inv2[0]*gamma[5]*dv12; 
  p2_over_gamma[2] = 2.738612787525831*jacob_vel_inv2[1]*gamma[16]*dv12+1.224744871391589*jacob_vel_inv2[0]*gamma[6]*dv12; 
  p2_over_gamma[3] = 2.449489742783178*jacob_vel_inv2[2]*gamma[9]*dv12+2.738612787525831*jacob_vel_inv2[0]*gamma[9]*dv12+1.224744871391589*jacob_vel_inv2[1]*gamma[3]*dv12; 
  p2_over_gamma[4] = 2.738612787525831*jacob_vel_inv2[1]*gamma[19]*dv12+1.224744871391589*jacob_vel_inv2[0]*gamma[10]*dv12; 
  p2_over_gamma[5] = 2.449489742783178*jacob_vel_inv2[2]*gamma[15]*dv12+2.738612787525831*jacob_vel_inv2[0]*gamma[15]*dv12+1.224744871391589*jacob_vel_inv2[1]*gamma[5]*dv12; 
  p2_over_gamma[6] = 2.449489742783178*jacob_vel_inv2[2]*gamma[16]*dv12+2.738612787525831*jacob_vel_inv2[0]*gamma[16]*dv12+1.224744871391589*jacob_vel_inv2[1]*gamma[6]*dv12; 
  p2_over_gamma[7] = 2.738612787525831*jacob_vel_inv2[1]*gamma[21]*dv12+1.224744871391589*jacob_vel_inv2[0]*gamma[13]*dv12; 
  p2_over_gamma[8] = 2.738612787525831*jacob_vel_inv2[1]*gamma[22]*dv12+1.224744871391589*jacob_vel_inv2[0]*gamma[14]*dv12; 
  p2_over_gamma[9] = 2.449489742783178*jacob_vel_inv2[1]*gamma[9]*dv12+1.224744871391589*jacob_vel_inv2[2]*gamma[3]*dv12; 
  p2_over_gamma[10] = 2.449489742783178*jacob_vel_inv2[2]*gamma[19]*dv12+2.738612787525831*jacob_vel_inv2[0]*gamma[19]*dv12+1.224744871391589*jacob_vel_inv2[1]*gamma[10]*dv12; 
  p2_over_gamma[11] = 2.738612787525831*jacob_vel_inv2[1]*gamma[24]*dv12+1.224744871391589*jacob_vel_inv2[0]*gamma[17]*dv12; 
  p2_over_gamma[12] = 2.738612787525831*jacob_vel_inv2[1]*gamma[25]*dv12+1.224744871391589*jacob_vel_inv2[0]*gamma[18]*dv12; 
  p2_over_gamma[13] = 2.449489742783178*jacob_vel_inv2[2]*gamma[21]*dv12+2.738612787525831*jacob_vel_inv2[0]*gamma[21]*dv12+1.224744871391589*jacob_vel_inv2[1]*gamma[13]*dv12; 
  p2_over_gamma[14] = 2.449489742783178*jacob_vel_inv2[2]*gamma[22]*dv12+2.738612787525831*jacob_vel_inv2[0]*gamma[22]*dv12+1.224744871391589*jacob_vel_inv2[1]*gamma[14]*dv12; 
  p2_over_gamma[15] = 2.449489742783178*jacob_vel_inv2[1]*gamma[15]*dv12+1.224744871391589*jacob_vel_inv2[2]*gamma[5]*dv12; 
  p2_over_gamma[16] = 2.449489742783178*jacob_vel_inv2[1]*gamma[16]*dv12+1.224744871391589*jacob_vel_inv2[2]*gamma[6]*dv12; 
  p2_over_gamma[17] = 2.449489742783178*jacob_vel_inv2[2]*gamma[24]*dv12+2.738612787525831*jacob_vel_inv2[0]*gamma[24]*dv12+1.224744871391589*jacob_vel_inv2[1]*gamma[17]*dv12; 
  p2_over_gamma[18] = 2.449489742783178*jacob_vel_inv2[2]*gamma[25]*dv12+2.738612787525831*jacob_vel_inv2[0]*gamma[25]*dv12+1.224744871391589*jacob_vel_inv2[1]*gamma[18]*dv12; 
  p2_over_gamma[19] = 2.449489742783178*jacob_vel_inv2[1]*gamma[19]*dv12+1.224744871391589*jacob_vel_inv2[2]*gamma[10]*dv12; 
  p2_over_gamma[20] = 2.738612787525831*jacob_vel_inv2[1]*gamma[26]*dv12+1.224744871391589*jacob_vel_inv2[0]*gamma[23]*dv12; 
  p2_over_gamma[21] = 2.449489742783178*jacob_vel_inv2[1]*gamma[21]*dv12+1.224744871391589*jacob_vel_inv2[2]*gamma[13]*dv12; 
  p2_over_gamma[22] = 2.449489742783178*jacob_vel_inv2[1]*gamma[22]*dv12+1.224744871391589*jacob_vel_inv2[2]*gamma[14]*dv12; 
  p2_over_gamma[23] = 2.449489742783178*jacob_vel_inv2[2]*gamma[26]*dv12+2.738612787525831*jacob_vel_inv2[0]*gamma[26]*dv12+1.224744871391589*jacob_vel_inv2[1]*gamma[23]*dv12; 
  p2_over_gamma[24] = 2.449489742783178*jacob_vel_inv2[1]*gamma[24]*dv12+1.224744871391589*jacob_vel_inv2[2]*gamma[17]*dv12; 
  p2_over_gamma[25] = 2.449489742783178*jacob_vel_inv2[1]*gamma[25]*dv12+1.224744871391589*jacob_vel_inv2[2]*gamma[18]*dv12; 
  p2_over_gamma[26] = 2.449489742783178*jacob_vel_inv2[1]*gamma[26]*dv12+1.224744871391589*jacob_vel_inv2[2]*gamma[23]*dv12; 

  out[0] += (p0_over_gamma[26]*f[75]+p0_over_gamma[25]*f[65]+p0_over_gamma[24]*f[63]+p0_over_gamma[23]*f[59]+p0_over_gamma[22]*f[49]+p0_over_gamma[21]*f[48]+p0_over_gamma[20]*f[46]+p0_over_gamma[19]*f[43]+p0_over_gamma[18]*f[40]+p0_over_gamma[17]*f[38]+p0_over_gamma[16]*f[30]+p0_over_gamma[15]*f[29]+p0_over_gamma[14]*f[27]+p0_over_gamma[13]*f[26]+p0_over_gamma[12]*f[24]+p0_over_gamma[11]*f[22]+p0_over_gamma[10]*f[18]+p0_over_gamma[9]*f[14]+p0_over_gamma[8]*f[13]+p0_over_gamma[7]*f[12]+p0_over_gamma[6]*f[10]+p0_over_gamma[5]*f[9]+p0_over_gamma[4]*f[7]+p0_over_gamma[3]*f[4]+p0_over_gamma[2]*f[3]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (1.0*p0_over_gamma[26]*f[79]+p0_over_gamma[25]*f[71]+p0_over_gamma[24]*f[70]+p0_over_gamma[23]*f[68]+p0_over_gamma[22]*f[64]+p0_over_gamma[21]*f[61]+p0_over_gamma[20]*f[56]+1.0*p0_over_gamma[19]*f[53]+1.0*p0_over_gamma[18]*f[52]+1.0*p0_over_gamma[17]*f[51]+1.0*p0_over_gamma[16]*f[42]+1.0*p0_over_gamma[15]*f[41]+1.0*p0_over_gamma[14]*f[39]+1.0*p0_over_gamma[13]*f[36]+1.0*p0_over_gamma[12]*f[34]+1.0*p0_over_gamma[11]*f[33]+p0_over_gamma[10]*f[31]+1.0*p0_over_gamma[9]*f[28]+1.0*p0_over_gamma[8]*f[23]+1.0*p0_over_gamma[7]*f[20]+p0_over_gamma[6]*f[17]+p0_over_gamma[5]*f[16]+p0_over_gamma[4]*f[15]+p0_over_gamma[3]*f[8]+p0_over_gamma[2]*f[6]+p0_over_gamma[1]*f[5]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p0_over_gamma[26]*f[80]+1.0*p0_over_gamma[25]*f[78]+1.0*p0_over_gamma[24]*f[77]+1.0*p0_over_gamma[23]*f[76]+p0_over_gamma[22]*f[74]+p0_over_gamma[21]*f[73]+p0_over_gamma[20]*f[72]+p0_over_gamma[19]*f[69]+p0_over_gamma[18]*f[67]+p0_over_gamma[17]*f[66]+1.0*p0_over_gamma[16]*f[62]+1.0*p0_over_gamma[15]*f[60]+1.0*p0_over_gamma[14]*f[58]+1.0*p0_over_gamma[13]*f[57]+1.0*p0_over_gamma[12]*f[55]+1.0*p0_over_gamma[11]*f[54]+1.0*p0_over_gamma[10]*f[50]+p0_over_gamma[9]*f[47]+p0_over_gamma[8]*f[45]+p0_over_gamma[7]*f[44]+p0_over_gamma[6]*f[37]+p0_over_gamma[5]*f[35]+p0_over_gamma[4]*f[32]+1.0*p0_over_gamma[3]*f[25]+1.0*p0_over_gamma[2]*f[21]+1.0*p0_over_gamma[1]*f[19]+p0_over_gamma[0]*f[11])*volFact; 
  out[3] += (p1_over_gamma[26]*f[75]+p1_over_gamma[25]*f[65]+p1_over_gamma[24]*f[63]+p1_over_gamma[23]*f[59]+p1_over_gamma[22]*f[49]+p1_over_gamma[21]*f[48]+p1_over_gamma[20]*f[46]+p1_over_gamma[19]*f[43]+p1_over_gamma[18]*f[40]+p1_over_gamma[17]*f[38]+p1_over_gamma[16]*f[30]+p1_over_gamma[15]*f[29]+p1_over_gamma[14]*f[27]+p1_over_gamma[13]*f[26]+p1_over_gamma[12]*f[24]+p1_over_gamma[11]*f[22]+p1_over_gamma[10]*f[18]+p1_over_gamma[9]*f[14]+p1_over_gamma[8]*f[13]+p1_over_gamma[7]*f[12]+p1_over_gamma[6]*f[10]+p1_over_gamma[5]*f[9]+p1_over_gamma[4]*f[7]+p1_over_gamma[3]*f[4]+p1_over_gamma[2]*f[3]+p1_over_gamma[1]*f[2]+f[0]*p1_over_gamma[0])*volFact; 
  out[4] += (1.0*p1_over_gamma[26]*f[79]+p1_over_gamma[25]*f[71]+p1_over_gamma[24]*f[70]+p1_over_gamma[23]*f[68]+p1_over_gamma[22]*f[64]+p1_over_gamma[21]*f[61]+p1_over_gamma[20]*f[56]+1.0*p1_over_gamma[19]*f[53]+1.0*p1_over_gamma[18]*f[52]+1.0*p1_over_gamma[17]*f[51]+1.0*p1_over_gamma[16]*f[42]+1.0*p1_over_gamma[15]*f[41]+1.0*p1_over_gamma[14]*f[39]+1.0*p1_over_gamma[13]*f[36]+1.0*p1_over_gamma[12]*f[34]+1.0*p1_over_gamma[11]*f[33]+p1_over_gamma[10]*f[31]+1.0*p1_over_gamma[9]*f[28]+1.0*p1_over_gamma[8]*f[23]+1.0*p1_over_gamma[7]*f[20]+p1_over_gamma[6]*f[17]+p1_over_gamma[5]*f[16]+p1_over_gamma[4]*f[15]+p1_over_gamma[3]*f[8]+p1_over_gamma[2]*f[6]+p1_over_gamma[1]*f[5]+p1_over_gamma[0]*f[1])*volFact; 
  out[5] += (p1_over_gamma[26]*f[80]+1.0*p1_over_gamma[25]*f[78]+1.0*p1_over_gamma[24]*f[77]+1.0*p1_over_gamma[23]*f[76]+p1_over_gamma[22]*f[74]+p1_over_gamma[21]*f[73]+p1_over_gamma[20]*f[72]+p1_over_gamma[19]*f[69]+p1_over_gamma[18]*f[67]+p1_over_gamma[17]*f[66]+1.0*p1_over_gamma[16]*f[62]+1.0*p1_over_gamma[15]*f[60]+1.0*p1_over_gamma[14]*f[58]+1.0*p1_over_gamma[13]*f[57]+1.0*p1_over_gamma[12]*f[55]+1.0*p1_over_gamma[11]*f[54]+1.0*p1_over_gamma[10]*f[50]+p1_over_gamma[9]*f[47]+p1_over_gamma[8]*f[45]+p1_over_gamma[7]*f[44]+p1_over_gamma[6]*f[37]+p1_over_gamma[5]*f[35]+p1_over_gamma[4]*f[32]+1.0*p1_over_gamma[3]*f[25]+1.0*p1_over_gamma[2]*f[21]+1.0*p1_over_gamma[1]*f[19]+p1_over_gamma[0]*f[11])*volFact; 
  out[6] += (p2_over_gamma[26]*f[75]+p2_over_gamma[25]*f[65]+p2_over_gamma[24]*f[63]+p2_over_gamma[23]*f[59]+p2_over_gamma[22]*f[49]+p2_over_gamma[21]*f[48]+p2_over_gamma[20]*f[46]+p2_over_gamma[19]*f[43]+p2_over_gamma[18]*f[40]+p2_over_gamma[17]*f[38]+p2_over_gamma[16]*f[30]+p2_over_gamma[15]*f[29]+p2_over_gamma[14]*f[27]+p2_over_gamma[13]*f[26]+p2_over_gamma[12]*f[24]+p2_over_gamma[11]*f[22]+p2_over_gamma[10]*f[18]+p2_over_gamma[9]*f[14]+p2_over_gamma[8]*f[13]+p2_over_gamma[7]*f[12]+p2_over_gamma[6]*f[10]+p2_over_gamma[5]*f[9]+p2_over_gamma[4]*f[7]+p2_over_gamma[3]*f[4]+p2_over_gamma[2]*f[3]+p2_over_gamma[1]*f[2]+f[0]*p2_over_gamma[0])*volFact; 
  out[7] += (1.0*p2_over_gamma[26]*f[79]+p2_over_gamma[25]*f[71]+p2_over_gamma[24]*f[70]+p2_over_gamma[23]*f[68]+p2_over_gamma[22]*f[64]+p2_over_gamma[21]*f[61]+p2_over_gamma[20]*f[56]+1.0*p2_over_gamma[19]*f[53]+1.0*p2_over_gamma[18]*f[52]+1.0*p2_over_gamma[17]*f[51]+1.0*p2_over_gamma[16]*f[42]+1.0*p2_over_gamma[15]*f[41]+1.0*p2_over_gamma[14]*f[39]+1.0*p2_over_gamma[13]*f[36]+1.0*p2_over_gamma[12]*f[34]+1.0*p2_over_gamma[11]*f[33]+p2_over_gamma[10]*f[31]+1.0*p2_over_gamma[9]*f[28]+1.0*p2_over_gamma[8]*f[23]+1.0*p2_over_gamma[7]*f[20]+p2_over_gamma[6]*f[17]+p2_over_gamma[5]*f[16]+p2_over_gamma[4]*f[15]+p2_over_gamma[3]*f[8]+p2_over_gamma[2]*f[6]+p2_over_gamma[1]*f[5]+p2_over_gamma[0]*f[1])*volFact; 
  out[8] += (p2_over_gamma[26]*f[80]+1.0*p2_over_gamma[25]*f[78]+1.0*p2_over_gamma[24]*f[77]+1.0*p2_over_gamma[23]*f[76]+p2_over_gamma[22]*f[74]+p2_over_gamma[21]*f[73]+p2_over_gamma[20]*f[72]+p2_over_gamma[19]*f[69]+p2_over_gamma[18]*f[67]+p2_over_gamma[17]*f[66]+1.0*p2_over_gamma[16]*f[62]+1.0*p2_over_gamma[15]*f[60]+1.0*p2_over_gamma[14]*f[58]+1.0*p2_over_gamma[13]*f[57]+1.0*p2_over_gamma[12]*f[55]+1.0*p2_over_gamma[11]*f[54]+1.0*p2_over_gamma[10]*f[50]+p2_over_gamma[9]*f[47]+p2_over_gamma[8]*f[45]+p2_over_gamma[7]*f[44]+p2_over_gamma[6]*f[37]+p2_over_gamma[5]*f[35]+p2_over_gamma[4]*f[32]+1.0*p2_over_gamma[3]*f[25]+1.0*p2_over_gamma[2]*f[21]+1.0*p2_over_gamma[1]*f[19]+p2_over_gamma[0]*f[11])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_M3i_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double *p0 = &vmap[0]; 
  const double *p1 = &vmap[4]; 
  const double *p2 = &vmap[8]; 
  out[0] += (2.0*p0[2]*f[12]+2.0*p0[1]*f[2]+2.0*f[0]*p0[0])*volFact; 
  out[1] += (2.0*p0[2]*f[20]+2.0*p0[1]*f[5]+2.0*p0[0]*f[1])*volFact; 
  out[2] += (2.0*p0[2]*f[44]+2.0*p0[1]*f[19]+2.0*p0[0]*f[11])*volFact; 
  out[3] += (2.0*p1[2]*f[13]+2.0*p1[1]*f[3]+2.0*f[0]*p1[0])*volFact; 
  out[4] += (2.0*p1[2]*f[23]+2.0*p1[1]*f[6]+2.0*p1[0]*f[1])*volFact; 
  out[5] += (2.0*p1[2]*f[45]+2.0*p1[1]*f[21]+2.0*p1[0]*f[11])*volFact; 
  out[6] += (2.0*p2[2]*f[14]+2.0*p2[1]*f[4]+2.0*f[0]*p2[0])*volFact; 
  out[7] += (2.0*p2[2]*f[28]+2.0*p2[1]*f[8]+2.0*p2[0]*f[1])*volFact; 
  out[8] += (2.0*p2[2]*f[47]+2.0*p2[1]*f[25]+2.0*p2[0]*f[11])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_int_mom_1x3v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *vmap, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double *p0 = &vmap[0]; 
  const double *p1 = &vmap[4]; 
  const double *p2 = &vmap[8]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (1.414213562373095*gamma[26]*f[75]+1.414213562373095*gamma[25]*f[65]+1.414213562373095*gamma[24]*f[63]+1.414213562373095*gamma[23]*f[59]+1.414213562373095*gamma[22]*f[49]+1.414213562373095*gamma[21]*f[48]+1.414213562373095*gamma[20]*f[46]+1.414213562373095*gamma[19]*f[43]+1.414213562373095*gamma[18]*f[40]+1.414213562373095*gamma[17]*f[38]+1.414213562373095*gamma[16]*f[30]+1.414213562373095*gamma[15]*f[29]+1.414213562373095*gamma[14]*f[27]+1.414213562373095*gamma[13]*f[26]+1.414213562373095*gamma[12]*f[24]+1.414213562373095*gamma[11]*f[22]+1.414213562373095*gamma[10]*f[18]+1.414213562373095*gamma[9]*f[14]+1.414213562373095*gamma[8]*f[13]+1.414213562373095*gamma[7]*f[12]+1.414213562373095*gamma[6]*f[10]+1.414213562373095*gamma[5]*f[9]+1.414213562373095*gamma[4]*f[7]+1.414213562373095*gamma[3]*f[4]+1.414213562373095*gamma[2]*f[3]+1.414213562373095*gamma[1]*f[2]+1.414213562373095*f[0]*gamma[0])*volFact; 
  out[2] += (2.828427124746191*p0[2]*f[12]+2.828427124746191*p0[1]*f[2]+2.828427124746191*f[0]*p0[0])*volFact; 
  out[3] += (2.828427124746191*p1[2]*f[13]+2.828427124746191*p1[1]*f[3]+2.828427124746191*f[0]*p1[0])*volFact; 
  out[4] += (2.828427124746191*p2[2]*f[14]+2.828427124746191*p2[1]*f[4]+2.828427124746191*f[0]*p2[0])*volFact; 
} 
