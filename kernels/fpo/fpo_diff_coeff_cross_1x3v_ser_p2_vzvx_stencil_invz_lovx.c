#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_invz_lovx(const double *dxv, double gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // diff_coeff: Output array for diffusion tensor. 

  double dv1_pv1 = 2.0/dxv[3]; 
  double dv1_pv2 = 2.0/dxv[1]; 
  double dv1_sq = 4.0/dxv[3]/dxv[1]; 
 
  const double* GCL = fpo_g_stencil[0]; 
  const double* GCC = fpo_g_stencil[1]; 
  const double* GCR = fpo_g_stencil[2]; 
  const double* GTL = fpo_g_stencil[3]; 
  const double* GTC = fpo_g_stencil[4]; 
  const double* GTR = fpo_g_stencil[5]; 

  const double* g_surf_CL = fpo_g_surf_stencil[0]; 
  const double* g_surf_CL_pv2 = &g_surf_CL[0]; 
  const double* g_surf_CC = fpo_g_surf_stencil[1]; 
  const double* g_surf_CC_pv2 = &g_surf_CC[0]; 
  const double* g_surf_CR = fpo_g_surf_stencil[2]; 
  const double* g_surf_CR_pv2 = &g_surf_CR[0]; 
  
  const double* g_surf_CC_pv1 = &g_surf_CC[40]; 
  const double* dgdpv1_surf_CC_pv2 = &fpo_dgdv_surf[40]; 
  const double* dgdpv2_surf_CC_pv1 = &fpo_dgdv_surf[120]; 
  const double* dgdpv1_surf_CC_pv1 = &fpo_dgdv_surf[160]; 
  
  double surft1_upper[20], surft1_lower[20]; 
  double surft2_upper[20], surft2_lower[20]; 
  
  double *diff_coeff_vxvy = &diff_coeff[48]; 
  double *diff_coeff_vxvz = &diff_coeff[96]; 
  double *diff_coeff_vyvx = &diff_coeff[144]; 
  double *diff_coeff_vyvz = &diff_coeff[240]; 
  double *diff_coeff_vzvx = &diff_coeff[288]; 
  double *diff_coeff_vzvy = &diff_coeff[336]; 
  
  double *out = diff_coeff_vzvx; 
  
  surft1_upper[0] = -(0.12168640803947765*GTR[29])+0.12168640803947765*GTL[29]+0.12168640803947765*GCR[29]-0.12168640803947765*GCL[29]-0.12168640803947765*GTR[26]-0.12168640803947765*GTL[26]+0.2433728160789553*GTC[26]-0.12168640803947765*GCR[26]-0.12168640803947765*GCL[26]+0.2433728160789553*GCC[26]+0.08646852977022904*GTR[14]-0.08646852977022904*GTL[14]+0.08646852977022904*GCR[14]-0.08646852977022904*GCL[14]+0.08646852977022904*GTR[12]-0.08646852977022904*GTL[12]+0.08646852977022904*GCR[12]-0.08646852977022904*GCL[12]+0.1750503603816304*GTR[9]+0.1750503603816304*GTL[9]-0.3501007207632608*GTC[9]-0.1750503603816304*GCR[9]-0.1750503603816304*GCL[9]+0.3501007207632608*GCC[9]-0.12438815100070813*GTR[4]-0.12438815100070813*GTL[4]+0.24877630200141632*GTC[4]-0.12438815100070813*GCR[4]-0.12438815100070813*GCL[4]+0.24877630200141632*GCC[4]-0.12438815100070813*GTR[2]+0.12438815100070813*GTL[2]+0.12438815100070813*GCR[2]-0.12438815100070813*GCL[2]+0.0883883476483184*GTR[0]-0.0883883476483184*GTL[0]+0.0883883476483184*GCR[0]-0.0883883476483184*GCL[0]; 
  surft1_upper[1] = -(0.12168640803947765*GTR[41])+0.12168640803947765*GTL[41]+0.12168640803947765*GCR[41]-0.12168640803947765*GCL[41]-0.12168640803947765*GTR[36]-0.12168640803947765*GTL[36]+0.2433728160789553*GTC[36]-0.12168640803947765*GCR[36]-0.12168640803947765*GCL[36]+0.2433728160789553*GCC[36]+0.08646852977022904*GTR[28]-0.08646852977022904*GTL[28]+0.08646852977022904*GCR[28]-0.08646852977022904*GCL[28]+0.08646852977022904*GTR[20]-0.08646852977022904*GTL[20]+0.08646852977022904*GCR[20]-0.08646852977022904*GCL[20]+0.1750503603816304*GTR[16]+0.1750503603816304*GTL[16]-0.3501007207632608*GTC[16]-0.1750503603816304*GCR[16]-0.1750503603816304*GCL[16]+0.3501007207632608*GCC[16]-0.12438815100070813*GTR[8]-0.12438815100070813*GTL[8]+0.24877630200141632*GTC[8]-0.12438815100070813*GCR[8]-0.12438815100070813*GCL[8]+0.24877630200141632*GCC[8]-0.12438815100070813*GTR[5]+0.12438815100070813*GTL[5]+0.12438815100070813*GCR[5]-0.12438815100070813*GCL[5]+0.0883883476483184*GTR[1]-0.0883883476483184*GTL[1]+0.0883883476483184*GCR[1]-0.0883883476483184*GCL[1]; 
  surft1_upper[2] = -(0.12168640803947765*GTR[43])+0.12168640803947765*GTL[43]+0.12168640803947765*GCR[43]-0.12168640803947765*GCL[43]-0.12168640803947765*GTR[38]-0.12168640803947765*GTL[38]+0.2433728160789553*GTC[38]-0.12168640803947765*GCR[38]-0.12168640803947765*GCL[38]+0.2433728160789553*GCC[38]+0.08646852977022904*GTR[30]-0.08646852977022904*GTL[30]+0.08646852977022904*GCR[30]-0.08646852977022904*GCL[30]+0.08646852977022904*GTR[22]-0.08646852977022904*GTL[22]+0.08646852977022904*GCR[22]-0.08646852977022904*GCL[22]+0.1750503603816304*GTR[18]+0.1750503603816304*GTL[18]-0.3501007207632608*GTC[18]-0.1750503603816304*GCR[18]-0.1750503603816304*GCL[18]+0.3501007207632608*GCC[18]-0.12438815100070813*GTR[10]-0.12438815100070813*GTL[10]+0.24877630200141632*GTC[10]-0.12438815100070813*GCR[10]-0.12438815100070813*GCL[10]+0.24877630200141632*GCC[10]-0.12438815100070813*GTR[7]+0.12438815100070813*GTL[7]+0.12438815100070813*GCR[7]-0.12438815100070813*GCL[7]+0.0883883476483184*GTR[3]-0.0883883476483184*GTL[3]+0.0883883476483184*GCR[3]-0.0883883476483184*GCL[3]; 
  surft1_upper[3] = -(0.2107670413149332*GTR[29])-0.2107670413149332*GTL[29]-0.4215340826298664*GTC[29]+0.2107670413149332*GCR[29]+0.2107670413149332*GCL[29]+0.4215340826298664*GCC[29]-0.2107670413149332*GTR[26]+0.2107670413149332*GTL[26]-0.2107670413149332*GCR[26]+0.2107670413149332*GCL[26]+0.1497678868178187*GTR[14]+0.1497678868178187*GTL[14]+0.29953577363563744*GTC[14]+0.1497678868178187*GCR[14]+0.1497678868178187*GCL[14]+0.29953577363563744*GCC[14]+0.1497678868178187*GTR[12]+0.1497678868178187*GTL[12]-0.29953577363563744*GTC[12]+0.1497678868178187*GCR[12]+0.1497678868178187*GCL[12]-0.29953577363563744*GCC[12]+0.30319611806422586*GTR[9]-0.30319611806422586*GTL[9]-0.30319611806422586*GCR[9]+0.30319611806422586*GCL[9]-0.21544659739277597*GTR[4]+0.21544659739277597*GTL[4]-0.21544659739277597*GCR[4]+0.21544659739277597*GCL[4]-0.21544659739277597*GTR[2]-0.21544659739277597*GTL[2]+0.43089319478555205*GTC[2]+0.21544659739277597*GCR[2]+0.21544659739277597*GCL[2]-0.43089319478555205*GCC[2]+0.15309310892394856*GTR[0]+0.15309310892394856*GTL[0]-0.3061862178478971*GTC[0]+0.15309310892394856*GCR[0]+0.15309310892394856*GCL[0]-0.3061862178478971*GCC[0]; 
  surft1_upper[4] = -(0.12168640803947765*GTR[47])+0.12168640803947765*GTL[47]+0.12168640803947765*GCR[47]-0.12168640803947765*GCL[47]-0.12168640803947765*GTR[45]-0.12168640803947765*GTL[45]+0.2433728160789553*GTC[45]-0.12168640803947765*GCR[45]-0.12168640803947765*GCL[45]+0.2433728160789553*GCC[45]+0.08646852977022904*GTR[42]-0.08646852977022904*GTL[42]+0.08646852977022904*GCR[42]-0.08646852977022904*GCL[42]+0.08646852977022904*GTR[33]-0.08646852977022904*GTL[33]+0.08646852977022904*GCR[33]-0.08646852977022904*GCL[33]+0.1750503603816304*GTR[31]+0.1750503603816304*GTL[31]-0.3501007207632608*GTC[31]-0.1750503603816304*GCR[31]-0.1750503603816304*GCL[31]+0.3501007207632608*GCC[31]-0.12438815100070813*GTR[17]-0.12438815100070813*GTL[17]+0.24877630200141632*GTC[17]-0.12438815100070813*GCR[17]-0.12438815100070813*GCL[17]+0.24877630200141632*GCC[17]-0.12438815100070813*GTR[15]+0.12438815100070813*GTL[15]+0.12438815100070813*GCR[15]-0.12438815100070813*GCL[15]+0.0883883476483184*GTR[6]-0.0883883476483184*GTL[6]+0.0883883476483184*GCR[6]-0.0883883476483184*GCL[6]; 
  surft1_upper[5] = -(0.21076704131493318*GTR[41])-0.21076704131493318*GTL[41]-0.42153408262986636*GTC[41]+0.21076704131493318*GCR[41]+0.21076704131493318*GCL[41]+0.42153408262986636*GCC[41]-0.21076704131493318*GTR[36]+0.21076704131493318*GTL[36]-0.21076704131493318*GCR[36]+0.21076704131493318*GCL[36]+0.1497678868178187*GTR[28]+0.1497678868178187*GTL[28]+0.29953577363563744*GTC[28]+0.1497678868178187*GCR[28]+0.1497678868178187*GCL[28]+0.29953577363563744*GCC[28]+0.1497678868178187*GTR[20]+0.1497678868178187*GTL[20]-0.29953577363563744*GTC[20]+0.1497678868178187*GCR[20]+0.1497678868178187*GCL[20]-0.29953577363563744*GCC[20]+0.30319611806422586*GTR[16]-0.30319611806422586*GTL[16]-0.30319611806422586*GCR[16]+0.30319611806422586*GCL[16]-0.21544659739277597*GTR[8]+0.21544659739277597*GTL[8]-0.21544659739277597*GCR[8]+0.21544659739277597*GCL[8]-0.21544659739277597*GTR[5]-0.21544659739277597*GTL[5]+0.43089319478555205*GTC[5]+0.21544659739277597*GCR[5]+0.21544659739277597*GCL[5]-0.43089319478555205*GCC[5]+0.15309310892394856*GTR[1]+0.15309310892394856*GTL[1]-0.3061862178478971*GTC[1]+0.15309310892394856*GCR[1]+0.15309310892394856*GCL[1]-0.3061862178478971*GCC[1]; 
  surft1_upper[6] = -(0.21076704131493318*GTR[43])-0.21076704131493318*GTL[43]-0.42153408262986636*GTC[43]+0.21076704131493318*GCR[43]+0.21076704131493318*GCL[43]+0.42153408262986636*GCC[43]-0.21076704131493318*GTR[38]+0.21076704131493318*GTL[38]-0.21076704131493318*GCR[38]+0.21076704131493318*GCL[38]+0.1497678868178187*GTR[30]+0.1497678868178187*GTL[30]+0.29953577363563744*GTC[30]+0.1497678868178187*GCR[30]+0.1497678868178187*GCL[30]+0.29953577363563744*GCC[30]+0.1497678868178187*GTR[22]+0.1497678868178187*GTL[22]-0.29953577363563744*GTC[22]+0.1497678868178187*GCR[22]+0.1497678868178187*GCL[22]-0.29953577363563744*GCC[22]+0.30319611806422586*GTR[18]-0.30319611806422586*GTL[18]-0.30319611806422586*GCR[18]+0.30319611806422586*GCL[18]-0.21544659739277597*GTR[10]+0.21544659739277597*GTL[10]-0.21544659739277597*GCR[10]+0.21544659739277597*GCL[10]-0.21544659739277597*GTR[7]-0.21544659739277597*GTL[7]+0.43089319478555205*GTC[7]+0.21544659739277597*GCR[7]+0.21544659739277597*GCL[7]-0.43089319478555205*GCC[7]+0.15309310892394856*GTR[3]+0.15309310892394856*GTL[3]-0.3061862178478971*GTC[3]+0.15309310892394856*GCR[3]+0.15309310892394856*GCL[3]-0.3061862178478971*GCC[3]; 
  surft1_upper[7] = 0.1750503603816304*GTR[35]+0.1750503603816304*GTL[35]-0.3501007207632608*GTC[35]-0.1750503603816304*GCR[35]-0.1750503603816304*GCL[35]+0.3501007207632608*GCC[35]-0.12438815100070813*GTR[25]-0.12438815100070813*GTL[25]+0.24877630200141632*GTC[25]-0.12438815100070813*GCR[25]-0.12438815100070813*GCL[25]+0.24877630200141632*GCC[25]-0.12438815100070813*GTR[19]+0.12438815100070813*GTL[19]+0.12438815100070813*GCR[19]-0.12438815100070813*GCL[19]+0.0883883476483184*GTR[11]-0.0883883476483184*GTL[11]+0.0883883476483184*GCR[11]-0.0883883476483184*GCL[11]; 
  surft1_upper[8] = 0.1750503603816304*GTR[40]+0.1750503603816304*GTL[40]-0.3501007207632608*GTC[40]-0.1750503603816304*GCR[40]-0.1750503603816304*GCL[40]+0.3501007207632608*GCC[40]-0.12438815100070813*GTR[27]-0.12438815100070813*GTL[27]+0.24877630200141632*GTC[27]-0.12438815100070813*GCR[27]-0.12438815100070813*GCL[27]+0.24877630200141632*GCC[27]-0.12438815100070813*GTR[24]+0.12438815100070813*GTL[24]+0.12438815100070813*GCR[24]-0.12438815100070813*GCL[24]+0.0883883476483184*GTR[13]-0.0883883476483184*GTL[13]+0.0883883476483184*GCR[13]-0.0883883476483184*GCL[13]; 
  surft1_upper[9] = -(0.27209908031404895*GTR[29])+0.27209908031404895*GTL[29]+0.27209908031404895*GCR[29]-0.27209908031404895*GCL[29]-0.27209908031404895*GTR[26]-0.27209908031404895*GTL[26]-0.7953665424564508*GTC[26]-0.27209908031404895*GCR[26]-0.27209908031404895*GCL[26]-0.7953665424564508*GCC[26]+0.1933495104806964*GTR[14]-0.1933495104806964*GTL[14]+0.1933495104806964*GCR[14]-0.1933495104806964*GCL[14]+0.1933495104806964*GTR[12]-0.1933495104806964*GTL[12]+0.1933495104806964*GCR[12]-0.1933495104806964*GCL[12]+0.3914245052991616*GTR[9]+0.3914245052991616*GTL[9]+1.1441639385667801*GTC[9]-0.3914245052991616*GCR[9]-0.3914245052991616*GCL[9]-1.1441639385667801*GCC[9]-0.2781403612330919*GTR[4]-0.2781403612330919*GTL[4]-0.8130256712967302*GTC[4]-0.2781403612330919*GCR[4]-0.2781403612330919*GCL[4]-0.8130256712967302*GCC[4]-0.2781403612330919*GTR[2]+0.2781403612330919*GTL[2]+0.2781403612330919*GCR[2]-0.2781403612330919*GCL[2]+0.19764235376052364*GTR[0]-0.19764235376052364*GTL[0]+0.19764235376052364*GCR[0]-0.19764235376052364*GCL[0]; 
  surft1_upper[10] = -(0.2107670413149332*GTR[47])-0.2107670413149332*GTL[47]-0.4215340826298664*GTC[47]+0.2107670413149332*GCR[47]+0.2107670413149332*GCL[47]+0.4215340826298664*GCC[47]-0.2107670413149332*GTR[45]+0.2107670413149332*GTL[45]-0.2107670413149332*GCR[45]+0.2107670413149332*GCL[45]+0.1497678868178187*GTR[42]+0.1497678868178187*GTL[42]+0.29953577363563744*GTC[42]+0.1497678868178187*GCR[42]+0.1497678868178187*GCL[42]+0.29953577363563744*GCC[42]+0.1497678868178187*GTR[33]+0.1497678868178187*GTL[33]-0.29953577363563744*GTC[33]+0.1497678868178187*GCR[33]+0.1497678868178187*GCL[33]-0.29953577363563744*GCC[33]+0.30319611806422586*GTR[31]-0.30319611806422586*GTL[31]-0.30319611806422586*GCR[31]+0.30319611806422586*GCL[31]-0.21544659739277597*GTR[17]+0.21544659739277597*GTL[17]-0.21544659739277597*GCR[17]+0.21544659739277597*GCL[17]-0.21544659739277597*GTR[15]-0.21544659739277597*GTL[15]+0.43089319478555205*GTC[15]+0.21544659739277597*GCR[15]+0.21544659739277597*GCL[15]-0.43089319478555205*GCC[15]+0.15309310892394856*GTR[6]+0.15309310892394856*GTL[6]-0.3061862178478971*GTC[6]+0.15309310892394856*GCR[6]+0.15309310892394856*GCL[6]-0.3061862178478971*GCC[6]; 
  surft1_upper[11] = 0.1750503603816304*GTR[44]+0.1750503603816304*GTL[44]-0.3501007207632608*GTC[44]-0.1750503603816304*GCR[44]-0.1750503603816304*GCL[44]+0.3501007207632608*GCC[44]-0.12438815100070813*GTR[37]-0.12438815100070813*GTL[37]+0.24877630200141632*GTC[37]-0.12438815100070813*GCR[37]-0.12438815100070813*GCL[37]+0.24877630200141632*GCC[37]-0.12438815100070813*GTR[32]+0.12438815100070813*GTL[32]+0.12438815100070813*GCR[32]-0.12438815100070813*GCL[32]+0.0883883476483184*GTR[21]-0.0883883476483184*GTL[21]+0.0883883476483184*GCR[21]-0.0883883476483184*GCL[21]; 
  surft1_upper[12] = 0.1750503603816304*GTR[46]+0.1750503603816304*GTL[46]-0.3501007207632608*GTC[46]-0.1750503603816304*GCR[46]-0.1750503603816304*GCL[46]+0.3501007207632608*GCC[46]-0.12438815100070813*GTR[39]-0.12438815100070813*GTL[39]+0.24877630200141632*GTC[39]-0.12438815100070813*GCR[39]-0.12438815100070813*GCL[39]+0.24877630200141632*GCC[39]-0.12438815100070813*GTR[34]+0.12438815100070813*GTL[34]+0.12438815100070813*GCR[34]-0.12438815100070813*GCL[34]+0.0883883476483184*GTR[23]-0.0883883476483184*GTL[23]+0.0883883476483184*GCR[23]-0.0883883476483184*GCL[23]; 
  surft1_upper[13] = 0.303196118064226*GTR[35]-0.303196118064226*GTL[35]-0.303196118064226*GCR[35]+0.303196118064226*GCL[35]-0.21544659739277597*GTR[25]+0.21544659739277597*GTL[25]-0.21544659739277597*GCR[25]+0.21544659739277597*GCL[25]-0.21544659739277597*GTR[19]-0.21544659739277597*GTL[19]+0.43089319478555205*GTC[19]+0.21544659739277597*GCR[19]+0.21544659739277597*GCL[19]-0.43089319478555205*GCC[19]+0.15309310892394856*GTR[11]+0.15309310892394856*GTL[11]-0.3061862178478971*GTC[11]+0.15309310892394856*GCR[11]+0.15309310892394856*GCL[11]-0.3061862178478971*GCC[11]; 
  surft1_upper[14] = 0.303196118064226*GTR[40]-0.303196118064226*GTL[40]-0.303196118064226*GCR[40]+0.303196118064226*GCL[40]-0.21544659739277597*GTR[27]+0.21544659739277597*GTL[27]-0.21544659739277597*GCR[27]+0.21544659739277597*GCL[27]-0.21544659739277597*GTR[24]-0.21544659739277597*GTL[24]+0.43089319478555205*GTC[24]+0.21544659739277597*GCR[24]+0.21544659739277597*GCL[24]-0.43089319478555205*GCC[24]+0.15309310892394856*GTR[13]+0.15309310892394856*GTL[13]-0.3061862178478971*GTC[13]+0.15309310892394856*GCR[13]+0.15309310892394856*GCL[13]-0.3061862178478971*GCC[13]; 
  surft1_upper[15] = -(0.27209908031404895*GTR[41])+0.27209908031404895*GTL[41]+0.27209908031404895*GCR[41]-0.27209908031404895*GCL[41]-0.27209908031404895*GTR[36]-0.27209908031404895*GTL[36]-0.7953665424564508*GTC[36]-0.27209908031404895*GCR[36]-0.27209908031404895*GCL[36]-0.7953665424564508*GCC[36]+0.1933495104806964*GTR[28]-0.1933495104806964*GTL[28]+0.1933495104806964*GCR[28]-0.1933495104806964*GCL[28]+0.1933495104806964*GTR[20]-0.1933495104806964*GTL[20]+0.1933495104806964*GCR[20]-0.1933495104806964*GCL[20]+0.39142450529916156*GTR[16]+0.39142450529916156*GTL[16]+1.14416393856678*GTC[16]-0.39142450529916156*GCR[16]-0.39142450529916156*GCL[16]-1.14416393856678*GCC[16]-0.2781403612330919*GTR[8]-0.2781403612330919*GTL[8]-0.8130256712967302*GTC[8]-0.2781403612330919*GCR[8]-0.2781403612330919*GCL[8]-0.8130256712967302*GCC[8]-0.2781403612330919*GTR[5]+0.2781403612330919*GTL[5]+0.2781403612330919*GCR[5]-0.2781403612330919*GCL[5]+0.19764235376052366*GTR[1]-0.19764235376052366*GTL[1]+0.19764235376052366*GCR[1]-0.19764235376052366*GCL[1]; 
  surft1_upper[16] = -(0.27209908031404895*GTR[43])+0.27209908031404895*GTL[43]+0.27209908031404895*GCR[43]-0.27209908031404895*GCL[43]-0.27209908031404895*GTR[38]-0.27209908031404895*GTL[38]-0.7953665424564508*GTC[38]-0.27209908031404895*GCR[38]-0.27209908031404895*GCL[38]-0.7953665424564508*GCC[38]+0.1933495104806964*GTR[30]-0.1933495104806964*GTL[30]+0.1933495104806964*GCR[30]-0.1933495104806964*GCL[30]+0.1933495104806964*GTR[22]-0.1933495104806964*GTL[22]+0.1933495104806964*GCR[22]-0.1933495104806964*GCL[22]+0.39142450529916156*GTR[18]+0.39142450529916156*GTL[18]+1.14416393856678*GTC[18]-0.39142450529916156*GCR[18]-0.39142450529916156*GCL[18]-1.14416393856678*GCC[18]-0.2781403612330919*GTR[10]-0.2781403612330919*GTL[10]-0.8130256712967302*GTC[10]-0.2781403612330919*GCR[10]-0.2781403612330919*GCL[10]-0.8130256712967302*GCC[10]-0.2781403612330919*GTR[7]+0.2781403612330919*GTL[7]+0.2781403612330919*GCR[7]-0.2781403612330919*GCL[7]+0.19764235376052366*GTR[3]-0.19764235376052366*GTL[3]+0.19764235376052366*GCR[3]-0.19764235376052366*GCL[3]; 
  surft1_upper[17] = 0.303196118064226*GTR[44]-0.303196118064226*GTL[44]-0.303196118064226*GCR[44]+0.303196118064226*GCL[44]-0.21544659739277597*GTR[37]+0.21544659739277597*GTL[37]-0.21544659739277597*GCR[37]+0.21544659739277597*GCL[37]-0.21544659739277597*GTR[32]-0.21544659739277597*GTL[32]+0.43089319478555205*GTC[32]+0.21544659739277597*GCR[32]+0.21544659739277597*GCL[32]-0.43089319478555205*GCC[32]+0.15309310892394856*GTR[21]+0.15309310892394856*GTL[21]-0.3061862178478971*GTC[21]+0.15309310892394856*GCR[21]+0.15309310892394856*GCL[21]-0.3061862178478971*GCC[21]; 
  surft1_upper[18] = 0.303196118064226*GTR[46]-0.303196118064226*GTL[46]-0.303196118064226*GCR[46]+0.303196118064226*GCL[46]-0.21544659739277597*GTR[39]+0.21544659739277597*GTL[39]-0.21544659739277597*GCR[39]+0.21544659739277597*GCL[39]-0.21544659739277597*GTR[34]-0.21544659739277597*GTL[34]+0.43089319478555205*GTC[34]+0.21544659739277597*GCR[34]+0.21544659739277597*GCL[34]-0.43089319478555205*GCC[34]+0.15309310892394856*GTR[23]+0.15309310892394856*GTL[23]-0.3061862178478971*GTC[23]+0.15309310892394856*GCR[23]+0.15309310892394856*GCL[23]-0.3061862178478971*GCC[23]; 
  surft1_upper[19] = -(0.27209908031404895*GTR[47])+0.27209908031404895*GTL[47]+0.27209908031404895*GCR[47]-0.27209908031404895*GCL[47]-0.27209908031404895*GTR[45]-0.27209908031404895*GTL[45]-0.7953665424564508*GTC[45]-0.27209908031404895*GCR[45]-0.27209908031404895*GCL[45]-0.7953665424564508*GCC[45]+0.1933495104806964*GTR[42]-0.1933495104806964*GTL[42]+0.1933495104806964*GCR[42]-0.1933495104806964*GCL[42]+0.1933495104806964*GTR[33]-0.1933495104806964*GTL[33]+0.1933495104806964*GCR[33]-0.1933495104806964*GCL[33]+0.3914245052991616*GTR[31]+0.3914245052991616*GTL[31]+1.1441639385667801*GTC[31]-0.3914245052991616*GCR[31]-0.3914245052991616*GCL[31]-1.1441639385667801*GCC[31]-0.2781403612330919*GTR[17]-0.2781403612330919*GTL[17]-0.8130256712967302*GTC[17]-0.2781403612330919*GCR[17]-0.2781403612330919*GCL[17]-0.8130256712967302*GCC[17]-0.2781403612330919*GTR[15]+0.2781403612330919*GTL[15]+0.2781403612330919*GCR[15]-0.2781403612330919*GCL[15]+0.19764235376052364*GTR[6]-0.19764235376052364*GTL[6]+0.19764235376052364*GCR[6]-0.19764235376052364*GCL[6]; 
  surft1_lower[0] = dgdpv1_surf_CC_pv2[0]/dv1_pv1; 
  surft1_lower[1] = dgdpv1_surf_CC_pv2[1]/dv1_pv1; 
  surft1_lower[2] = dgdpv1_surf_CC_pv2[2]/dv1_pv1; 
  surft1_lower[3] = dgdpv1_surf_CC_pv2[3]/dv1_pv1; 
  surft1_lower[4] = dgdpv1_surf_CC_pv2[4]/dv1_pv1; 
  surft1_lower[5] = dgdpv1_surf_CC_pv2[5]/dv1_pv1; 
  surft1_lower[6] = dgdpv1_surf_CC_pv2[6]/dv1_pv1; 
  surft1_lower[7] = dgdpv1_surf_CC_pv2[7]/dv1_pv1; 
  surft1_lower[8] = dgdpv1_surf_CC_pv2[8]/dv1_pv1; 
  surft1_lower[9] = dgdpv1_surf_CC_pv2[9]/dv1_pv1; 
  surft1_lower[10] = dgdpv1_surf_CC_pv2[10]/dv1_pv1; 
  surft1_lower[11] = dgdpv1_surf_CC_pv2[11]/dv1_pv1; 
  surft1_lower[12] = dgdpv1_surf_CC_pv2[12]/dv1_pv1; 
  surft1_lower[13] = dgdpv1_surf_CC_pv2[13]/dv1_pv1; 
  surft1_lower[14] = dgdpv1_surf_CC_pv2[14]/dv1_pv1; 
  surft1_lower[15] = dgdpv1_surf_CC_pv2[15]/dv1_pv1; 
  surft1_lower[16] = dgdpv1_surf_CC_pv2[16]/dv1_pv1; 
  surft1_lower[17] = dgdpv1_surf_CC_pv2[17]/dv1_pv1; 
  surft1_lower[18] = dgdpv1_surf_CC_pv2[18]/dv1_pv1; 
  surft1_lower[19] = dgdpv1_surf_CC_pv2[19]/dv1_pv1; 

  surft2_upper[0] = 0.34587411908091625*GCR[14]+0.34587411908091625*GCC[14]-0.49755260400283263*GCR[4]+0.49755260400283263*GCC[4]+0.3535533905932737*GCR[0]+0.3535533905932737*GCC[0]; 
  surft2_upper[1] = 0.34587411908091625*GCR[28]+0.34587411908091625*GCC[28]-0.49755260400283263*GCR[8]+0.49755260400283263*GCC[8]+0.3535533905932737*GCR[1]+0.3535533905932737*GCC[1]; 
  surft2_upper[2] = 0.34587411908091625*GCR[29]+0.34587411908091625*GCC[29]-0.49755260400283263*GCR[9]+0.49755260400283263*GCC[9]+0.3535533905932737*GCR[2]+0.3535533905932737*GCC[2]; 
  surft2_upper[3] = 0.34587411908091625*GCR[30]+0.34587411908091625*GCC[30]-0.49755260400283263*GCR[10]+0.49755260400283263*GCC[10]+0.3535533905932737*GCR[3]+0.3535533905932737*GCC[3]; 
  surft2_upper[4] = 0.34587411908091625*GCR[41]+0.34587411908091625*GCC[41]-0.49755260400283263*GCR[16]+0.49755260400283263*GCC[16]+0.3535533905932737*GCR[5]+0.3535533905932737*GCC[5]; 
  surft2_upper[5] = 0.34587411908091625*GCR[42]+0.34587411908091625*GCC[42]-0.49755260400283263*GCR[17]+0.49755260400283263*GCC[17]+0.3535533905932737*GCR[6]+0.3535533905932737*GCC[6]; 
  surft2_upper[6] = 0.34587411908091625*GCR[43]+0.34587411908091625*GCC[43]-0.49755260400283263*GCR[18]+0.49755260400283263*GCC[18]+0.3535533905932737*GCR[7]+0.3535533905932737*GCC[7]; 
  surft2_upper[7] = -(0.49755260400283263*GCR[25])+0.49755260400283263*GCC[25]+0.3535533905932737*GCR[11]+0.3535533905932737*GCC[11]; 
  surft2_upper[8] = -(0.49755260400283263*GCR[26])+0.49755260400283263*GCC[26]+0.3535533905932737*GCR[12]+0.3535533905932737*GCC[12]; 
  surft2_upper[9] = -(0.49755260400283263*GCR[27])+0.49755260400283263*GCC[27]+0.3535533905932737*GCR[13]+0.3535533905932737*GCC[13]; 
  surft2_upper[10] = 0.34587411908091625*GCR[47]+0.34587411908091625*GCC[47]-0.49755260400283263*GCR[31]+0.49755260400283263*GCC[31]+0.3535533905932737*GCR[15]+0.3535533905932737*GCC[15]; 
  surft2_upper[11] = -(0.49755260400283263*GCR[35])+0.49755260400283263*GCC[35]+0.3535533905932737*GCR[19]+0.3535533905932737*GCC[19]; 
  surft2_upper[12] = -(0.49755260400283263*GCR[36])+0.49755260400283263*GCC[36]+0.3535533905932737*GCR[20]+0.3535533905932737*GCC[20]; 
  surft2_upper[13] = -(0.49755260400283263*GCR[37])+0.49755260400283263*GCC[37]+0.3535533905932737*GCR[21]+0.3535533905932737*GCC[21]; 
  surft2_upper[14] = -(0.49755260400283263*GCR[38])+0.49755260400283263*GCC[38]+0.3535533905932737*GCR[22]+0.3535533905932737*GCC[22]; 
  surft2_upper[15] = -(0.49755260400283263*GCR[39])+0.49755260400283263*GCC[39]+0.3535533905932737*GCR[23]+0.3535533905932737*GCC[23]; 
  surft2_upper[16] = -(0.49755260400283263*GCR[40])+0.49755260400283263*GCC[40]+0.3535533905932737*GCR[24]+0.3535533905932737*GCC[24]; 
  surft2_upper[17] = -(0.49755260400283263*GCR[44])+0.49755260400283263*GCC[44]+0.3535533905932737*GCR[32]+0.3535533905932737*GCC[32]; 
  surft2_upper[18] = -(0.49755260400283263*GCR[45])+0.49755260400283263*GCC[45]+0.3535533905932737*GCR[33]+0.3535533905932737*GCC[33]; 
  surft2_upper[19] = -(0.49755260400283263*GCR[46])+0.49755260400283263*GCC[46]+0.3535533905932737*GCR[34]+0.3535533905932737*GCC[34]; 
  surft2_lower[0] = 0.34587411908091625*GCL[14]+0.34587411908091625*GCC[14]+0.49755260400283263*GCL[4]-0.49755260400283263*GCC[4]+0.3535533905932737*GCL[0]+0.3535533905932737*GCC[0]; 
  surft2_lower[1] = 0.34587411908091625*GCL[28]+0.34587411908091625*GCC[28]+0.49755260400283263*GCL[8]-0.49755260400283263*GCC[8]+0.3535533905932737*GCL[1]+0.3535533905932737*GCC[1]; 
  surft2_lower[2] = 0.34587411908091625*GCL[29]+0.34587411908091625*GCC[29]+0.49755260400283263*GCL[9]-0.49755260400283263*GCC[9]+0.3535533905932737*GCL[2]+0.3535533905932737*GCC[2]; 
  surft2_lower[3] = 0.34587411908091625*GCL[30]+0.34587411908091625*GCC[30]+0.49755260400283263*GCL[10]-0.49755260400283263*GCC[10]+0.3535533905932737*GCL[3]+0.3535533905932737*GCC[3]; 
  surft2_lower[4] = 0.34587411908091625*GCL[41]+0.34587411908091625*GCC[41]+0.49755260400283263*GCL[16]-0.49755260400283263*GCC[16]+0.3535533905932737*GCL[5]+0.3535533905932737*GCC[5]; 
  surft2_lower[5] = 0.34587411908091625*GCL[42]+0.34587411908091625*GCC[42]+0.49755260400283263*GCL[17]-0.49755260400283263*GCC[17]+0.3535533905932737*GCL[6]+0.3535533905932737*GCC[6]; 
  surft2_lower[6] = 0.34587411908091625*GCL[43]+0.34587411908091625*GCC[43]+0.49755260400283263*GCL[18]-0.49755260400283263*GCC[18]+0.3535533905932737*GCL[7]+0.3535533905932737*GCC[7]; 
  surft2_lower[7] = 0.49755260400283263*GCL[25]-0.49755260400283263*GCC[25]+0.3535533905932737*GCL[11]+0.3535533905932737*GCC[11]; 
  surft2_lower[8] = 0.49755260400283263*GCL[26]-0.49755260400283263*GCC[26]+0.3535533905932737*GCL[12]+0.3535533905932737*GCC[12]; 
  surft2_lower[9] = 0.49755260400283263*GCL[27]-0.49755260400283263*GCC[27]+0.3535533905932737*GCL[13]+0.3535533905932737*GCC[13]; 
  surft2_lower[10] = 0.34587411908091625*GCL[47]+0.34587411908091625*GCC[47]+0.49755260400283263*GCL[31]-0.49755260400283263*GCC[31]+0.3535533905932737*GCL[15]+0.3535533905932737*GCC[15]; 
  surft2_lower[11] = 0.49755260400283263*GCL[35]-0.49755260400283263*GCC[35]+0.3535533905932737*GCL[19]+0.3535533905932737*GCC[19]; 
  surft2_lower[12] = 0.49755260400283263*GCL[36]-0.49755260400283263*GCC[36]+0.3535533905932737*GCL[20]+0.3535533905932737*GCC[20]; 
  surft2_lower[13] = 0.49755260400283263*GCL[37]-0.49755260400283263*GCC[37]+0.3535533905932737*GCL[21]+0.3535533905932737*GCC[21]; 
  surft2_lower[14] = 0.49755260400283263*GCL[38]-0.49755260400283263*GCC[38]+0.3535533905932737*GCL[22]+0.3535533905932737*GCC[22]; 
  surft2_lower[15] = 0.49755260400283263*GCL[39]-0.49755260400283263*GCC[39]+0.3535533905932737*GCL[23]+0.3535533905932737*GCC[23]; 
  surft2_lower[16] = 0.49755260400283263*GCL[40]-0.49755260400283263*GCC[40]+0.3535533905932737*GCL[24]+0.3535533905932737*GCC[24]; 
  surft2_lower[17] = 0.49755260400283263*GCL[44]-0.49755260400283263*GCC[44]+0.3535533905932737*GCL[32]+0.3535533905932737*GCC[32]; 
  surft2_lower[18] = 0.49755260400283263*GCL[45]-0.49755260400283263*GCC[45]+0.3535533905932737*GCL[33]+0.3535533905932737*GCC[33]; 
  surft2_lower[19] = 0.49755260400283263*GCL[46]-0.49755260400283263*GCC[46]+0.3535533905932737*GCL[34]+0.3535533905932737*GCC[34]; 

  out[0] = 0.7071067811865475*surft1_upper[0]*dv1_sq*gamma-0.7071067811865475*surft1_lower[0]*dv1_sq*gamma; 
  out[1] = 0.7071067811865475*surft1_upper[1]*dv1_sq*gamma-0.7071067811865475*surft1_lower[1]*dv1_sq*gamma; 
  out[2] = -(1.224744871391589*surft2_upper[0]*dv1_sq*gamma)+1.224744871391589*surft2_lower[0]*dv1_sq*gamma+1.224744871391589*surft1_upper[0]*dv1_sq*gamma+1.224744871391589*surft1_lower[0]*dv1_sq*gamma; 
  out[3] = 0.7071067811865475*surft1_upper[2]*dv1_sq*gamma-0.7071067811865475*surft1_lower[2]*dv1_sq*gamma; 
  out[4] = 0.7071067811865475*surft1_upper[3]*dv1_sq*gamma-0.7071067811865475*surft1_lower[3]*dv1_sq*gamma; 
  out[5] = -(1.224744871391589*surft2_upper[1]*dv1_sq*gamma)+1.224744871391589*surft2_lower[1]*dv1_sq*gamma+1.224744871391589*surft1_upper[1]*dv1_sq*gamma+1.224744871391589*surft1_lower[1]*dv1_sq*gamma; 
  out[6] = 0.7071067811865475*surft1_upper[4]*dv1_sq*gamma-0.7071067811865475*surft1_lower[4]*dv1_sq*gamma; 
  out[7] = -(1.224744871391589*surft2_upper[3]*dv1_sq*gamma)+1.224744871391589*surft2_lower[3]*dv1_sq*gamma+1.224744871391589*surft1_upper[2]*dv1_sq*gamma+1.224744871391589*surft1_lower[2]*dv1_sq*gamma; 
  out[8] = 0.7071067811865475*surft1_upper[5]*dv1_sq*gamma-0.7071067811865475*surft1_lower[5]*dv1_sq*gamma; 
  out[9] = 1.224744871391589*surft1_upper[3]*dv1_sq*gamma+1.224744871391589*surft1_lower[3]*dv1_sq*gamma-2.1213203435596424*surft2_upper[0]*dv1_sq*gamma-2.1213203435596424*surft2_lower[0]*dv1_sq*gamma+3.0*GCC[0]*dv1_sq*gamma; 
  out[10] = 0.7071067811865475*surft1_upper[6]*dv1_sq*gamma-0.7071067811865475*surft1_lower[6]*dv1_sq*gamma; 
  out[11] = 0.7071067811865475*surft1_upper[7]*dv1_sq*gamma-0.7071067811865475*surft1_lower[7]*dv1_sq*gamma; 
  out[12] = -(2.7386127875258306*surft2_upper[2]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[2]*dv1_sq*gamma+1.5811388300841895*surft1_upper[0]*dv1_sq*gamma-1.5811388300841895*surft1_lower[0]*dv1_sq*gamma; 
  out[13] = 0.7071067811865475*surft1_upper[8]*dv1_sq*gamma-0.7071067811865475*surft1_lower[8]*dv1_sq*gamma; 
  out[14] = 0.7071067811865475*surft1_upper[9]*dv1_sq*gamma-0.7071067811865475*surft1_lower[9]*dv1_sq*gamma; 
  out[15] = -(1.224744871391589*surft2_upper[5]*dv1_sq*gamma)+1.224744871391589*surft2_lower[5]*dv1_sq*gamma+1.224744871391589*surft1_upper[4]*dv1_sq*gamma+1.224744871391589*surft1_lower[4]*dv1_sq*gamma; 
  out[16] = 1.224744871391589*surft1_upper[5]*dv1_sq*gamma+1.224744871391589*surft1_lower[5]*dv1_sq*gamma-2.1213203435596424*surft2_upper[1]*dv1_sq*gamma-2.1213203435596424*surft2_lower[1]*dv1_sq*gamma+3.0*GCC[1]*dv1_sq*gamma; 
  out[17] = 0.7071067811865475*surft1_upper[10]*dv1_sq*gamma-0.7071067811865475*surft1_lower[10]*dv1_sq*gamma; 
  out[18] = 1.224744871391589*surft1_upper[6]*dv1_sq*gamma+1.224744871391589*surft1_lower[6]*dv1_sq*gamma-2.1213203435596424*surft2_upper[3]*dv1_sq*gamma-2.1213203435596424*surft2_lower[3]*dv1_sq*gamma+3.0*GCC[3]*dv1_sq*gamma; 
  out[19] = -(1.224744871391589*surft2_upper[7]*dv1_sq*gamma)+1.224744871391589*surft2_lower[7]*dv1_sq*gamma+1.224744871391589*surft1_upper[7]*dv1_sq*gamma+1.224744871391589*surft1_lower[7]*dv1_sq*gamma; 
  out[20] = -(2.7386127875258306*surft2_upper[4]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[4]*dv1_sq*gamma+1.5811388300841898*surft1_upper[1]*dv1_sq*gamma-1.5811388300841898*surft1_lower[1]*dv1_sq*gamma; 
  out[21] = 0.7071067811865475*surft1_upper[11]*dv1_sq*gamma-0.7071067811865475*surft1_lower[11]*dv1_sq*gamma; 
  out[22] = -(2.7386127875258306*surft2_upper[6]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[6]*dv1_sq*gamma+1.5811388300841898*surft1_upper[2]*dv1_sq*gamma-1.5811388300841898*surft1_lower[2]*dv1_sq*gamma; 
  out[23] = 0.7071067811865475*surft1_upper[12]*dv1_sq*gamma-0.7071067811865475*surft1_lower[12]*dv1_sq*gamma; 
  out[24] = -(1.224744871391589*surft2_upper[9]*dv1_sq*gamma)+1.224744871391589*surft2_lower[9]*dv1_sq*gamma+1.224744871391589*surft1_upper[8]*dv1_sq*gamma+1.224744871391589*surft1_lower[8]*dv1_sq*gamma; 
  out[25] = 0.7071067811865475*surft1_upper[13]*dv1_sq*gamma-0.7071067811865475*surft1_lower[13]*dv1_sq*gamma; 
  out[26] = 1.5811388300841898*surft1_upper[3]*dv1_sq*gamma-1.5811388300841898*surft1_lower[3]*dv1_sq*gamma-4.743416490252569*surft2_upper[2]*dv1_sq*gamma-4.743416490252569*surft2_lower[2]*dv1_sq*gamma+6.7082039324993685*GCC[2]*dv1_sq*gamma; 
  out[27] = 0.7071067811865475*surft1_upper[14]*dv1_sq*gamma-0.7071067811865475*surft1_lower[14]*dv1_sq*gamma; 
  out[28] = 0.7071067811865475*surft1_upper[15]*dv1_sq*gamma-0.7071067811865475*surft1_lower[15]*dv1_sq*gamma; 
  out[29] = 1.224744871391589*surft1_upper[9]*dv1_sq*gamma+1.224744871391589*surft1_lower[9]*dv1_sq*gamma+6.7082039324993685*GCC[4]*dv1_sq*gamma-2.7386127875258306*surft2_upper[0]*dv1_sq*gamma+2.7386127875258306*surft2_lower[0]*dv1_sq*gamma; 
  out[30] = 0.7071067811865475*surft1_upper[16]*dv1_sq*gamma-0.7071067811865475*surft1_lower[16]*dv1_sq*gamma; 
  out[31] = 1.224744871391589*surft1_upper[10]*dv1_sq*gamma+1.224744871391589*surft1_lower[10]*dv1_sq*gamma+3.0*GCC[6]*dv1_sq*gamma-2.1213203435596424*surft2_upper[5]*dv1_sq*gamma-2.1213203435596424*surft2_lower[5]*dv1_sq*gamma; 
  out[32] = -(1.224744871391589*surft2_upper[13]*dv1_sq*gamma)+1.224744871391589*surft2_lower[13]*dv1_sq*gamma+1.224744871391589*surft1_upper[11]*dv1_sq*gamma+1.224744871391589*surft1_lower[11]*dv1_sq*gamma; 
  out[33] = -(2.7386127875258306*surft2_upper[10]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[10]*dv1_sq*gamma+1.5811388300841895*surft1_upper[4]*dv1_sq*gamma-1.5811388300841895*surft1_lower[4]*dv1_sq*gamma; 
  out[34] = -(1.224744871391589*surft2_upper[15]*dv1_sq*gamma)+1.224744871391589*surft2_lower[15]*dv1_sq*gamma+1.224744871391589*surft1_upper[12]*dv1_sq*gamma+1.224744871391589*surft1_lower[12]*dv1_sq*gamma; 
  out[35] = 1.224744871391589*surft1_upper[13]*dv1_sq*gamma+1.224744871391589*surft1_lower[13]*dv1_sq*gamma+3.0*GCC[11]*dv1_sq*gamma-2.1213203435596424*surft2_upper[7]*dv1_sq*gamma-2.1213203435596424*surft2_lower[7]*dv1_sq*gamma; 
  out[36] = 1.5811388300841895*surft1_upper[5]*dv1_sq*gamma-1.5811388300841895*surft1_lower[5]*dv1_sq*gamma+6.708203932499369*GCC[5]*dv1_sq*gamma-4.743416490252569*surft2_upper[4]*dv1_sq*gamma-4.743416490252569*surft2_lower[4]*dv1_sq*gamma; 
  out[37] = 0.7071067811865475*surft1_upper[17]*dv1_sq*gamma-0.7071067811865475*surft1_lower[17]*dv1_sq*gamma; 
  out[38] = 6.708203932499369*GCC[7]*dv1_sq*gamma-4.743416490252569*surft2_upper[6]*dv1_sq*gamma-4.743416490252569*surft2_lower[6]*dv1_sq*gamma+1.5811388300841895*surft1_upper[6]*dv1_sq*gamma-1.5811388300841895*surft1_lower[6]*dv1_sq*gamma; 
  out[39] = 0.7071067811865475*surft1_upper[18]*dv1_sq*gamma-0.7071067811865475*surft1_lower[18]*dv1_sq*gamma; 
  out[40] = 1.224744871391589*surft1_upper[14]*dv1_sq*gamma+1.224744871391589*surft1_lower[14]*dv1_sq*gamma+3.0*GCC[13]*dv1_sq*gamma-2.1213203435596424*surft2_upper[9]*dv1_sq*gamma-2.1213203435596424*surft2_lower[9]*dv1_sq*gamma; 
  out[41] = 1.224744871391589*surft1_upper[15]*dv1_sq*gamma+1.224744871391589*surft1_lower[15]*dv1_sq*gamma+6.708203932499369*GCC[8]*dv1_sq*gamma-2.7386127875258306*surft2_upper[1]*dv1_sq*gamma+2.7386127875258306*surft2_lower[1]*dv1_sq*gamma; 
  out[42] = 0.7071067811865475*surft1_upper[19]*dv1_sq*gamma-0.7071067811865475*surft1_lower[19]*dv1_sq*gamma; 
  out[43] = 1.224744871391589*surft1_upper[16]*dv1_sq*gamma+1.224744871391589*surft1_lower[16]*dv1_sq*gamma+6.708203932499369*GCC[10]*dv1_sq*gamma-2.7386127875258306*surft2_upper[3]*dv1_sq*gamma+2.7386127875258306*surft2_lower[3]*dv1_sq*gamma; 
  out[44] = 3.0*GCC[21]*dv1_sq*gamma+1.224744871391589*surft1_upper[17]*dv1_sq*gamma+1.224744871391589*surft1_lower[17]*dv1_sq*gamma-2.1213203435596424*surft2_upper[13]*dv1_sq*gamma-2.1213203435596424*surft2_lower[13]*dv1_sq*gamma; 
  out[45] = 6.7082039324993685*GCC[15]*dv1_sq*gamma-4.743416490252569*surft2_upper[10]*dv1_sq*gamma-4.743416490252569*surft2_lower[10]*dv1_sq*gamma+1.5811388300841898*surft1_upper[10]*dv1_sq*gamma-1.5811388300841898*surft1_lower[10]*dv1_sq*gamma; 
  out[46] = 3.0*GCC[23]*dv1_sq*gamma+1.224744871391589*surft1_upper[18]*dv1_sq*gamma+1.224744871391589*surft1_lower[18]*dv1_sq*gamma-2.1213203435596424*surft2_upper[15]*dv1_sq*gamma-2.1213203435596424*surft2_lower[15]*dv1_sq*gamma; 
  out[47] = 1.224744871391589*surft1_upper[19]*dv1_sq*gamma+1.224744871391589*surft1_lower[19]*dv1_sq*gamma+6.7082039324993685*GCC[17]*dv1_sq*gamma-2.7386127875258306*surft2_upper[5]*dv1_sq*gamma+2.7386127875258306*surft2_lower[5]*dv1_sq*gamma; 
} 
