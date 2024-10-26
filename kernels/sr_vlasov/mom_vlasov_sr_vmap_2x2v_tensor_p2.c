#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_vmap_M1i_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *jacob_vel_inv, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double dv10 = 2.0/dxv[2]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  const double dv11 = 2.0/dxv[3]; 
  const double *jacob_vel_inv1 = &jacob_vel_inv[3]; 
  double p0_over_gamma[9] = {0.0}; 
  p0_over_gamma[0] = 2.738612787525831*jacob_vel_inv0[1]*gamma[4]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[1]*dv10; 
  p0_over_gamma[1] = 2.449489742783178*jacob_vel_inv0[2]*gamma[4]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[4]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[1]*dv10; 
  p0_over_gamma[2] = 2.738612787525831*jacob_vel_inv0[1]*gamma[6]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[3]*dv10; 
  p0_over_gamma[3] = 2.449489742783178*jacob_vel_inv0[2]*gamma[6]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[6]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[3]*dv10; 
  p0_over_gamma[4] = 2.449489742783178*jacob_vel_inv0[1]*gamma[4]*dv10+1.224744871391589*gamma[1]*jacob_vel_inv0[2]*dv10; 
  p0_over_gamma[5] = 2.738612787525831*jacob_vel_inv0[1]*gamma[8]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[7]*dv10; 
  p0_over_gamma[6] = 2.449489742783178*jacob_vel_inv0[1]*gamma[6]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[3]*dv10; 
  p0_over_gamma[7] = 2.449489742783178*jacob_vel_inv0[2]*gamma[8]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[8]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[7]*dv10; 
  p0_over_gamma[8] = 2.449489742783178*jacob_vel_inv0[1]*gamma[8]*dv10+1.224744871391589*jacob_vel_inv0[2]*gamma[7]*dv10; 

  double p1_over_gamma[9] = {0.0}; 
  p1_over_gamma[0] = 2.738612787525831*jacob_vel_inv1[1]*gamma[5]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[2]*dv11; 
  p1_over_gamma[1] = 2.738612787525831*jacob_vel_inv1[1]*gamma[7]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[3]*dv11; 
  p1_over_gamma[2] = 2.449489742783178*jacob_vel_inv1[2]*gamma[5]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[5]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[2]*dv11; 
  p1_over_gamma[3] = 2.449489742783178*jacob_vel_inv1[2]*gamma[7]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[7]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[3]*dv11; 
  p1_over_gamma[4] = 2.738612787525831*jacob_vel_inv1[1]*gamma[8]*dv11+1.224744871391589*jacob_vel_inv1[0]*gamma[6]*dv11; 
  p1_over_gamma[5] = 2.449489742783178*jacob_vel_inv1[1]*gamma[5]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[2]*dv11; 
  p1_over_gamma[6] = 2.449489742783178*jacob_vel_inv1[2]*gamma[8]*dv11+2.738612787525831*jacob_vel_inv1[0]*gamma[8]*dv11+1.224744871391589*jacob_vel_inv1[1]*gamma[6]*dv11; 
  p1_over_gamma[7] = 2.449489742783178*jacob_vel_inv1[1]*gamma[7]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[3]*dv11; 
  p1_over_gamma[8] = 2.449489742783178*jacob_vel_inv1[1]*gamma[8]*dv11+1.224744871391589*jacob_vel_inv1[2]*gamma[6]*dv11; 

  out[0] += (p0_over_gamma[8]*f[49]+p0_over_gamma[7]*f[30]+p0_over_gamma[6]*f[27]+p0_over_gamma[5]*f[14]+p0_over_gamma[4]*f[13]+p0_over_gamma[3]*f[10]+p0_over_gamma[2]*f[4]+p0_over_gamma[1]*f[3]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (p0_over_gamma[8]*f[64]+1.0*p0_over_gamma[7]*f[42]+1.0*p0_over_gamma[6]*f[39]+1.0*p0_over_gamma[5]*f[28]+1.0*p0_over_gamma[4]*f[23]+p0_over_gamma[3]*f[17]+p0_over_gamma[2]*f[8]+p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (p0_over_gamma[8]*f[65]+1.0*p0_over_gamma[7]*f[43]+1.0*p0_over_gamma[6]*f[40]+1.0*p0_over_gamma[5]*f[29]+1.0*p0_over_gamma[4]*f[24]+p0_over_gamma[3]*f[18]+p0_over_gamma[2]*f[9]+p0_over_gamma[1]*f[7]+p0_over_gamma[0]*f[2])*volFact; 
  out[3] += (p0_over_gamma[8]*f[71]+p0_over_gamma[7]*f[53]+p0_over_gamma[6]*f[52]+p0_over_gamma[5]*f[41]+p0_over_gamma[4]*f[34]+p0_over_gamma[3]*f[31]+p0_over_gamma[2]*f[16]+p0_over_gamma[1]*f[15]+p0_over_gamma[0]*f[5])*volFact; 
  out[4] += (p0_over_gamma[8]*f[74]+1.0*p0_over_gamma[7]*f[62]+1.0*p0_over_gamma[6]*f[58]+p0_over_gamma[5]*f[47]+p0_over_gamma[4]*f[45]+p0_over_gamma[3]*f[37]+1.0*p0_over_gamma[2]*f[25]+1.0*p0_over_gamma[1]*f[21]+p0_over_gamma[0]*f[11])*volFact; 
  out[5] += (p0_over_gamma[8]*f[75]+1.0*p0_over_gamma[7]*f[63]+1.0*p0_over_gamma[6]*f[59]+p0_over_gamma[5]*f[48]+p0_over_gamma[4]*f[46]+p0_over_gamma[3]*f[38]+1.0*p0_over_gamma[2]*f[26]+1.0*p0_over_gamma[1]*f[22]+p0_over_gamma[0]*f[12])*volFact; 
  out[6] += (p0_over_gamma[8]*f[78]+p0_over_gamma[7]*f[69]+p0_over_gamma[6]*f[67]+1.0*p0_over_gamma[5]*f[60]+1.0*p0_over_gamma[4]*f[55]+p0_over_gamma[3]*f[50]+1.0*p0_over_gamma[2]*f[35]+1.0*p0_over_gamma[1]*f[32]+p0_over_gamma[0]*f[19])*volFact; 
  out[7] += (p0_over_gamma[8]*f[79]+p0_over_gamma[7]*f[70]+p0_over_gamma[6]*f[68]+1.0*p0_over_gamma[5]*f[61]+1.0*p0_over_gamma[4]*f[56]+p0_over_gamma[3]*f[51]+1.0*p0_over_gamma[2]*f[36]+1.0*p0_over_gamma[1]*f[33]+p0_over_gamma[0]*f[20])*volFact; 
  out[8] += (p0_over_gamma[8]*f[80]+p0_over_gamma[7]*f[77]+p0_over_gamma[6]*f[76]+p0_over_gamma[5]*f[73]+p0_over_gamma[4]*f[72]+p0_over_gamma[3]*f[66]+p0_over_gamma[2]*f[57]+p0_over_gamma[1]*f[54]+p0_over_gamma[0]*f[44])*volFact; 
  out[9] += (p1_over_gamma[8]*f[49]+p1_over_gamma[7]*f[30]+p1_over_gamma[6]*f[27]+p1_over_gamma[5]*f[14]+p1_over_gamma[4]*f[13]+p1_over_gamma[3]*f[10]+p1_over_gamma[2]*f[4]+p1_over_gamma[1]*f[3]+f[0]*p1_over_gamma[0])*volFact; 
  out[10] += (p1_over_gamma[8]*f[64]+1.0*p1_over_gamma[7]*f[42]+1.0*p1_over_gamma[6]*f[39]+1.0*p1_over_gamma[5]*f[28]+1.0*p1_over_gamma[4]*f[23]+p1_over_gamma[3]*f[17]+p1_over_gamma[2]*f[8]+p1_over_gamma[1]*f[6]+p1_over_gamma[0]*f[1])*volFact; 
  out[11] += (p1_over_gamma[8]*f[65]+1.0*p1_over_gamma[7]*f[43]+1.0*p1_over_gamma[6]*f[40]+1.0*p1_over_gamma[5]*f[29]+1.0*p1_over_gamma[4]*f[24]+p1_over_gamma[3]*f[18]+p1_over_gamma[2]*f[9]+p1_over_gamma[1]*f[7]+p1_over_gamma[0]*f[2])*volFact; 
  out[12] += (p1_over_gamma[8]*f[71]+p1_over_gamma[7]*f[53]+p1_over_gamma[6]*f[52]+p1_over_gamma[5]*f[41]+p1_over_gamma[4]*f[34]+p1_over_gamma[3]*f[31]+p1_over_gamma[2]*f[16]+p1_over_gamma[1]*f[15]+p1_over_gamma[0]*f[5])*volFact; 
  out[13] += (p1_over_gamma[8]*f[74]+1.0*p1_over_gamma[7]*f[62]+1.0*p1_over_gamma[6]*f[58]+p1_over_gamma[5]*f[47]+p1_over_gamma[4]*f[45]+p1_over_gamma[3]*f[37]+1.0*p1_over_gamma[2]*f[25]+1.0*p1_over_gamma[1]*f[21]+p1_over_gamma[0]*f[11])*volFact; 
  out[14] += (p1_over_gamma[8]*f[75]+1.0*p1_over_gamma[7]*f[63]+1.0*p1_over_gamma[6]*f[59]+p1_over_gamma[5]*f[48]+p1_over_gamma[4]*f[46]+p1_over_gamma[3]*f[38]+1.0*p1_over_gamma[2]*f[26]+1.0*p1_over_gamma[1]*f[22]+p1_over_gamma[0]*f[12])*volFact; 
  out[15] += (p1_over_gamma[8]*f[78]+p1_over_gamma[7]*f[69]+p1_over_gamma[6]*f[67]+1.0*p1_over_gamma[5]*f[60]+1.0*p1_over_gamma[4]*f[55]+p1_over_gamma[3]*f[50]+1.0*p1_over_gamma[2]*f[35]+1.0*p1_over_gamma[1]*f[32]+p1_over_gamma[0]*f[19])*volFact; 
  out[16] += (p1_over_gamma[8]*f[79]+p1_over_gamma[7]*f[70]+p1_over_gamma[6]*f[68]+1.0*p1_over_gamma[5]*f[61]+1.0*p1_over_gamma[4]*f[56]+p1_over_gamma[3]*f[51]+1.0*p1_over_gamma[2]*f[36]+1.0*p1_over_gamma[1]*f[33]+p1_over_gamma[0]*f[20])*volFact; 
  out[17] += (p1_over_gamma[8]*f[80]+p1_over_gamma[7]*f[77]+p1_over_gamma[6]*f[76]+p1_over_gamma[5]*f[73]+p1_over_gamma[4]*f[72]+p1_over_gamma[3]*f[66]+p1_over_gamma[2]*f[57]+p1_over_gamma[1]*f[54]+p1_over_gamma[0]*f[44])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_M3i_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]*dxv[3]/4; 
  const double *p0 = &vmap[0]; 
  const double *p1 = &vmap[4]; 
  out[0] += (1.414213562373095*p0[2]*f[13]+1.414213562373095*p0[1]*f[3]+1.414213562373095*f[0]*p0[0])*volFact; 
  out[1] += (1.414213562373095*p0[2]*f[23]+1.414213562373095*p0[1]*f[6]+1.414213562373095*p0[0]*f[1])*volFact; 
  out[2] += (1.414213562373095*p0[2]*f[24]+1.414213562373095*p0[1]*f[7]+1.414213562373095*p0[0]*f[2])*volFact; 
  out[3] += (1.414213562373095*p0[2]*f[34]+1.414213562373095*p0[1]*f[15]+1.414213562373095*p0[0]*f[5])*volFact; 
  out[4] += (1.414213562373095*p0[2]*f[45]+1.414213562373095*p0[1]*f[21]+1.414213562373095*p0[0]*f[11])*volFact; 
  out[5] += (1.414213562373095*p0[2]*f[46]+1.414213562373095*p0[1]*f[22]+1.414213562373095*p0[0]*f[12])*volFact; 
  out[6] += (1.414213562373095*p0[2]*f[55]+1.414213562373095*p0[1]*f[32]+1.414213562373095*p0[0]*f[19])*volFact; 
  out[7] += (1.414213562373095*p0[2]*f[56]+1.414213562373095*p0[1]*f[33]+1.414213562373095*p0[0]*f[20])*volFact; 
  out[8] += (1.414213562373095*p0[2]*f[72]+1.414213562373095*p0[1]*f[54]+1.414213562373095*p0[0]*f[44])*volFact; 
  out[9] += (1.414213562373095*p1[2]*f[14]+1.414213562373095*p1[1]*f[4]+1.414213562373095*f[0]*p1[0])*volFact; 
  out[10] += (1.414213562373095*p1[2]*f[28]+1.414213562373095*p1[1]*f[8]+1.414213562373095*p1[0]*f[1])*volFact; 
  out[11] += (1.414213562373095*p1[2]*f[29]+1.414213562373095*p1[1]*f[9]+1.414213562373095*p1[0]*f[2])*volFact; 
  out[12] += (1.414213562373095*p1[2]*f[41]+1.414213562373095*p1[1]*f[16]+1.414213562373095*p1[0]*f[5])*volFact; 
  out[13] += (1.414213562373095*p1[2]*f[47]+1.414213562373095*p1[1]*f[25]+1.414213562373095*p1[0]*f[11])*volFact; 
  out[14] += (1.414213562373095*p1[2]*f[48]+1.414213562373095*p1[1]*f[26]+1.414213562373095*p1[0]*f[12])*volFact; 
  out[15] += (1.414213562373095*p1[2]*f[60]+1.414213562373095*p1[1]*f[35]+1.414213562373095*p1[0]*f[19])*volFact; 
  out[16] += (1.414213562373095*p1[2]*f[61]+1.414213562373095*p1[1]*f[36]+1.414213562373095*p1[0]*f[20])*volFact; 
  out[17] += (1.414213562373095*p1[2]*f[73]+1.414213562373095*p1[1]*f[57]+1.414213562373095*p1[0]*f[44])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_int_mom_2x2v_tensor_p2(const double *w, const double *dxv, const int *idx, const double *vmap, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double *p0 = &vmap[0]; 
  const double *p1 = &vmap[4]; 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.0*gamma[8]*f[49]+2.0*gamma[7]*f[30]+2.0*gamma[6]*f[27]+2.0*gamma[5]*f[14]+2.0*gamma[4]*f[13]+2.0*gamma[3]*f[10]+2.0*gamma[2]*f[4]+2.0*gamma[1]*f[3]+2.0*f[0]*gamma[0])*volFact; 
  out[2] += (2.828427124746191*p0[2]*f[13]+2.828427124746191*p0[1]*f[3]+2.828427124746191*f[0]*p0[0])*volFact; 
  out[3] += (2.828427124746191*p1[2]*f[14]+2.828427124746191*p1[1]*f[4]+2.828427124746191*f[0]*p1[0])*volFact; 
} 
