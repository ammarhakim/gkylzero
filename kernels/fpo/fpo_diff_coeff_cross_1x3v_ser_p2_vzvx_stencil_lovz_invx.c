#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_lovz_invx(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[3]; 
  double dv2 = 2.0/dxv[1]; 
  double dv1_sq = 4.0/dxv[3]/dxv[1]; 
 
  double surft1_lo[20] = {0.0}; 
  double surft1_up[20] = {0.0}; 
  double surft2_lo[20] = {0.0}; 
  double surft2_up[20] = {0.0}; 
  double vol[48] = {0.0}; 
  double *out = &diff_coeff[288]; 

  const double* GCL = fpo_g_stencil[0]; 
  const double* GTL = fpo_g_stencil[1]; 
  const double* GCC = fpo_g_stencil[2]; 
  const double* G_surf_CC_vx = &fpo_g_surf_stencil[2][0]; 
  const double* GTC = fpo_g_stencil[3]; 
  const double* GCR = fpo_g_stencil[4]; 
  const double* GTR = fpo_g_stencil[5]; 
  const double* dGdvx_surf_CC_vz = &fpo_dgdv_surf[80]; 

  surft1_lo[0] = dGdvx_surf_CC_vz[0]/dv2; 
  surft1_lo[1] = dGdvx_surf_CC_vz[1]/dv2; 
  surft1_lo[2] = dGdvx_surf_CC_vz[2]/dv2; 
  surft1_lo[3] = dGdvx_surf_CC_vz[3]/dv2; 
  surft1_lo[4] = dGdvx_surf_CC_vz[4]/dv2; 
  surft1_lo[5] = dGdvx_surf_CC_vz[5]/dv2; 
  surft1_lo[6] = dGdvx_surf_CC_vz[6]/dv2; 
  surft1_lo[7] = dGdvx_surf_CC_vz[7]/dv2; 
  surft1_lo[8] = dGdvx_surf_CC_vz[8]/dv2; 
  surft1_lo[9] = dGdvx_surf_CC_vz[9]/dv2; 
  surft1_lo[10] = dGdvx_surf_CC_vz[10]/dv2; 
  surft1_lo[11] = dGdvx_surf_CC_vz[11]/dv2; 
  surft1_lo[12] = dGdvx_surf_CC_vz[12]/dv2; 
  surft1_lo[13] = dGdvx_surf_CC_vz[13]/dv2; 
  surft1_lo[14] = dGdvx_surf_CC_vz[14]/dv2; 
  surft1_lo[15] = dGdvx_surf_CC_vz[15]/dv2; 
  surft1_lo[16] = dGdvx_surf_CC_vz[16]/dv2; 
  surft1_lo[17] = dGdvx_surf_CC_vz[17]/dv2; 
  surft1_lo[18] = dGdvx_surf_CC_vz[18]/dv2; 
  surft1_lo[19] = dGdvx_surf_CC_vz[19]/dv2; 
  surft1_up[0] = -(0.12168640803947765*GTR[29])-0.12168640803947765*GTL[29]+0.2433728160789553*GTC[29]-0.12168640803947765*GCR[29]-0.12168640803947765*GCL[29]+0.2433728160789553*GCC[29]-0.12168640803947765*GTR[26]+0.12168640803947765*GTL[26]+0.12168640803947765*GCR[26]-0.12168640803947765*GCL[26]+0.08646852977022904*GTR[14]-0.08646852977022904*GTL[14]+0.08646852977022904*GCR[14]-0.08646852977022904*GCL[14]+0.08646852977022904*GTR[12]-0.08646852977022904*GTL[12]+0.08646852977022904*GCR[12]-0.08646852977022904*GCL[12]+0.1750503603816304*GTR[9]+0.1750503603816304*GTL[9]-0.3501007207632608*GTC[9]-0.1750503603816304*GCR[9]-0.1750503603816304*GCL[9]+0.3501007207632608*GCC[9]-0.12438815100070813*GTR[4]+0.12438815100070813*GTL[4]+0.12438815100070813*GCR[4]-0.12438815100070813*GCL[4]-0.12438815100070813*GTR[2]-0.12438815100070813*GTL[2]+0.24877630200141632*GTC[2]-0.12438815100070813*GCR[2]-0.12438815100070813*GCL[2]+0.24877630200141632*GCC[2]+0.0883883476483184*GTR[0]-0.0883883476483184*GTL[0]+0.0883883476483184*GCR[0]-0.0883883476483184*GCL[0]; 
  surft1_up[1] = -(0.12168640803947765*GTR[41])-0.12168640803947765*GTL[41]+0.2433728160789553*GTC[41]-0.12168640803947765*GCR[41]-0.12168640803947765*GCL[41]+0.2433728160789553*GCC[41]-0.12168640803947765*GTR[36]+0.12168640803947765*GTL[36]+0.12168640803947765*GCR[36]-0.12168640803947765*GCL[36]+0.08646852977022904*GTR[28]-0.08646852977022904*GTL[28]+0.08646852977022904*GCR[28]-0.08646852977022904*GCL[28]+0.08646852977022904*GTR[20]-0.08646852977022904*GTL[20]+0.08646852977022904*GCR[20]-0.08646852977022904*GCL[20]+0.1750503603816304*GTR[16]+0.1750503603816304*GTL[16]-0.3501007207632608*GTC[16]-0.1750503603816304*GCR[16]-0.1750503603816304*GCL[16]+0.3501007207632608*GCC[16]-0.12438815100070813*GTR[8]+0.12438815100070813*GTL[8]+0.12438815100070813*GCR[8]-0.12438815100070813*GCL[8]-0.12438815100070813*GTR[5]-0.12438815100070813*GTL[5]+0.24877630200141632*GTC[5]-0.12438815100070813*GCR[5]-0.12438815100070813*GCL[5]+0.24877630200141632*GCC[5]+0.0883883476483184*GTR[1]-0.0883883476483184*GTL[1]+0.0883883476483184*GCR[1]-0.0883883476483184*GCL[1]; 
  surft1_up[2] = -(0.2107670413149332*GTR[29])+0.2107670413149332*GTL[29]-0.2107670413149332*GCR[29]+0.2107670413149332*GCL[29]-0.2107670413149332*GTR[26]-0.2107670413149332*GTL[26]-0.4215340826298664*GTC[26]+0.2107670413149332*GCR[26]+0.2107670413149332*GCL[26]+0.4215340826298664*GCC[26]+0.1497678868178187*GTR[14]+0.1497678868178187*GTL[14]-0.29953577363563744*GTC[14]+0.1497678868178187*GCR[14]+0.1497678868178187*GCL[14]-0.29953577363563744*GCC[14]+0.1497678868178187*GTR[12]+0.1497678868178187*GTL[12]+0.29953577363563744*GTC[12]+0.1497678868178187*GCR[12]+0.1497678868178187*GCL[12]+0.29953577363563744*GCC[12]+0.30319611806422586*GTR[9]-0.30319611806422586*GTL[9]-0.30319611806422586*GCR[9]+0.30319611806422586*GCL[9]-0.21544659739277597*GTR[4]-0.21544659739277597*GTL[4]+0.43089319478555205*GTC[4]+0.21544659739277597*GCR[4]+0.21544659739277597*GCL[4]-0.43089319478555205*GCC[4]-0.21544659739277597*GTR[2]+0.21544659739277597*GTL[2]-0.21544659739277597*GCR[2]+0.21544659739277597*GCL[2]+0.15309310892394856*GTR[0]+0.15309310892394856*GTL[0]-0.3061862178478971*GTC[0]+0.15309310892394856*GCR[0]+0.15309310892394856*GCL[0]-0.3061862178478971*GCC[0]; 
  surft1_up[3] = -(0.12168640803947765*GTR[43])-0.12168640803947765*GTL[43]+0.2433728160789553*GTC[43]-0.12168640803947765*GCR[43]-0.12168640803947765*GCL[43]+0.2433728160789553*GCC[43]-0.12168640803947765*GTR[38]+0.12168640803947765*GTL[38]+0.12168640803947765*GCR[38]-0.12168640803947765*GCL[38]+0.08646852977022904*GTR[30]-0.08646852977022904*GTL[30]+0.08646852977022904*GCR[30]-0.08646852977022904*GCL[30]+0.08646852977022904*GTR[22]-0.08646852977022904*GTL[22]+0.08646852977022904*GCR[22]-0.08646852977022904*GCL[22]+0.1750503603816304*GTR[18]+0.1750503603816304*GTL[18]-0.3501007207632608*GTC[18]-0.1750503603816304*GCR[18]-0.1750503603816304*GCL[18]+0.3501007207632608*GCC[18]-0.12438815100070813*GTR[10]+0.12438815100070813*GTL[10]+0.12438815100070813*GCR[10]-0.12438815100070813*GCL[10]-0.12438815100070813*GTR[7]-0.12438815100070813*GTL[7]+0.24877630200141632*GTC[7]-0.12438815100070813*GCR[7]-0.12438815100070813*GCL[7]+0.24877630200141632*GCC[7]+0.0883883476483184*GTR[3]-0.0883883476483184*GTL[3]+0.0883883476483184*GCR[3]-0.0883883476483184*GCL[3]; 
  surft1_up[4] = -(0.21076704131493318*GTR[41])+0.21076704131493318*GTL[41]-0.21076704131493318*GCR[41]+0.21076704131493318*GCL[41]-0.21076704131493318*GTR[36]-0.21076704131493318*GTL[36]-0.42153408262986636*GTC[36]+0.21076704131493318*GCR[36]+0.21076704131493318*GCL[36]+0.42153408262986636*GCC[36]+0.1497678868178187*GTR[28]+0.1497678868178187*GTL[28]-0.29953577363563744*GTC[28]+0.1497678868178187*GCR[28]+0.1497678868178187*GCL[28]-0.29953577363563744*GCC[28]+0.1497678868178187*GTR[20]+0.1497678868178187*GTL[20]+0.29953577363563744*GTC[20]+0.1497678868178187*GCR[20]+0.1497678868178187*GCL[20]+0.29953577363563744*GCC[20]+0.30319611806422586*GTR[16]-0.30319611806422586*GTL[16]-0.30319611806422586*GCR[16]+0.30319611806422586*GCL[16]-0.21544659739277597*GTR[8]-0.21544659739277597*GTL[8]+0.43089319478555205*GTC[8]+0.21544659739277597*GCR[8]+0.21544659739277597*GCL[8]-0.43089319478555205*GCC[8]-0.21544659739277597*GTR[5]+0.21544659739277597*GTL[5]-0.21544659739277597*GCR[5]+0.21544659739277597*GCL[5]+0.15309310892394856*GTR[1]+0.15309310892394856*GTL[1]-0.3061862178478971*GTC[1]+0.15309310892394856*GCR[1]+0.15309310892394856*GCL[1]-0.3061862178478971*GCC[1]; 
  surft1_up[5] = -(0.12168640803947765*GTR[47])-0.12168640803947765*GTL[47]+0.2433728160789553*GTC[47]-0.12168640803947765*GCR[47]-0.12168640803947765*GCL[47]+0.2433728160789553*GCC[47]-0.12168640803947765*GTR[45]+0.12168640803947765*GTL[45]+0.12168640803947765*GCR[45]-0.12168640803947765*GCL[45]+0.08646852977022904*GTR[42]-0.08646852977022904*GTL[42]+0.08646852977022904*GCR[42]-0.08646852977022904*GCL[42]+0.08646852977022904*GTR[33]-0.08646852977022904*GTL[33]+0.08646852977022904*GCR[33]-0.08646852977022904*GCL[33]+0.1750503603816304*GTR[31]+0.1750503603816304*GTL[31]-0.3501007207632608*GTC[31]-0.1750503603816304*GCR[31]-0.1750503603816304*GCL[31]+0.3501007207632608*GCC[31]-0.12438815100070813*GTR[17]+0.12438815100070813*GTL[17]+0.12438815100070813*GCR[17]-0.12438815100070813*GCL[17]-0.12438815100070813*GTR[15]-0.12438815100070813*GTL[15]+0.24877630200141632*GTC[15]-0.12438815100070813*GCR[15]-0.12438815100070813*GCL[15]+0.24877630200141632*GCC[15]+0.0883883476483184*GTR[6]-0.0883883476483184*GTL[6]+0.0883883476483184*GCR[6]-0.0883883476483184*GCL[6]; 
  surft1_up[6] = -(0.21076704131493318*GTR[43])+0.21076704131493318*GTL[43]-0.21076704131493318*GCR[43]+0.21076704131493318*GCL[43]-0.21076704131493318*GTR[38]-0.21076704131493318*GTL[38]-0.42153408262986636*GTC[38]+0.21076704131493318*GCR[38]+0.21076704131493318*GCL[38]+0.42153408262986636*GCC[38]+0.1497678868178187*GTR[30]+0.1497678868178187*GTL[30]-0.29953577363563744*GTC[30]+0.1497678868178187*GCR[30]+0.1497678868178187*GCL[30]-0.29953577363563744*GCC[30]+0.1497678868178187*GTR[22]+0.1497678868178187*GTL[22]+0.29953577363563744*GTC[22]+0.1497678868178187*GCR[22]+0.1497678868178187*GCL[22]+0.29953577363563744*GCC[22]+0.30319611806422586*GTR[18]-0.30319611806422586*GTL[18]-0.30319611806422586*GCR[18]+0.30319611806422586*GCL[18]-0.21544659739277597*GTR[10]-0.21544659739277597*GTL[10]+0.43089319478555205*GTC[10]+0.21544659739277597*GCR[10]+0.21544659739277597*GCL[10]-0.43089319478555205*GCC[10]-0.21544659739277597*GTR[7]+0.21544659739277597*GTL[7]-0.21544659739277597*GCR[7]+0.21544659739277597*GCL[7]+0.15309310892394856*GTR[3]+0.15309310892394856*GTL[3]-0.3061862178478971*GTC[3]+0.15309310892394856*GCR[3]+0.15309310892394856*GCL[3]-0.3061862178478971*GCC[3]; 
  surft1_up[7] = 0.1750503603816304*GTR[35]+0.1750503603816304*GTL[35]-0.3501007207632608*GTC[35]-0.1750503603816304*GCR[35]-0.1750503603816304*GCL[35]+0.3501007207632608*GCC[35]-0.12438815100070813*GTR[25]+0.12438815100070813*GTL[25]+0.12438815100070813*GCR[25]-0.12438815100070813*GCL[25]-0.12438815100070813*GTR[19]-0.12438815100070813*GTL[19]+0.24877630200141632*GTC[19]-0.12438815100070813*GCR[19]-0.12438815100070813*GCL[19]+0.24877630200141632*GCC[19]+0.0883883476483184*GTR[11]-0.0883883476483184*GTL[11]+0.0883883476483184*GCR[11]-0.0883883476483184*GCL[11]; 
  surft1_up[8] = -(0.27209908031404895*GTR[29])-0.27209908031404895*GTL[29]-0.7953665424564508*GTC[29]-0.27209908031404895*GCR[29]-0.27209908031404895*GCL[29]-0.7953665424564508*GCC[29]-0.27209908031404895*GTR[26]+0.27209908031404895*GTL[26]+0.27209908031404895*GCR[26]-0.27209908031404895*GCL[26]+0.1933495104806964*GTR[14]-0.1933495104806964*GTL[14]+0.1933495104806964*GCR[14]-0.1933495104806964*GCL[14]+0.1933495104806964*GTR[12]-0.1933495104806964*GTL[12]+0.1933495104806964*GCR[12]-0.1933495104806964*GCL[12]+0.3914245052991616*GTR[9]+0.3914245052991616*GTL[9]+1.1441639385667801*GTC[9]-0.3914245052991616*GCR[9]-0.3914245052991616*GCL[9]-1.1441639385667801*GCC[9]-0.2781403612330919*GTR[4]+0.2781403612330919*GTL[4]+0.2781403612330919*GCR[4]-0.2781403612330919*GCL[4]-0.2781403612330919*GTR[2]-0.2781403612330919*GTL[2]-0.8130256712967302*GTC[2]-0.2781403612330919*GCR[2]-0.2781403612330919*GCL[2]-0.8130256712967302*GCC[2]+0.19764235376052364*GTR[0]-0.19764235376052364*GTL[0]+0.19764235376052364*GCR[0]-0.19764235376052364*GCL[0]; 
  surft1_up[9] = 0.1750503603816304*GTR[40]+0.1750503603816304*GTL[40]-0.3501007207632608*GTC[40]-0.1750503603816304*GCR[40]-0.1750503603816304*GCL[40]+0.3501007207632608*GCC[40]-0.12438815100070813*GTR[27]+0.12438815100070813*GTL[27]+0.12438815100070813*GCR[27]-0.12438815100070813*GCL[27]-0.12438815100070813*GTR[24]-0.12438815100070813*GTL[24]+0.24877630200141632*GTC[24]-0.12438815100070813*GCR[24]-0.12438815100070813*GCL[24]+0.24877630200141632*GCC[24]+0.0883883476483184*GTR[13]-0.0883883476483184*GTL[13]+0.0883883476483184*GCR[13]-0.0883883476483184*GCL[13]; 
  surft1_up[10] = -(0.2107670413149332*GTR[47])+0.2107670413149332*GTL[47]-0.2107670413149332*GCR[47]+0.2107670413149332*GCL[47]-0.2107670413149332*GTR[45]-0.2107670413149332*GTL[45]-0.4215340826298664*GTC[45]+0.2107670413149332*GCR[45]+0.2107670413149332*GCL[45]+0.4215340826298664*GCC[45]+0.1497678868178187*GTR[42]+0.1497678868178187*GTL[42]-0.29953577363563744*GTC[42]+0.1497678868178187*GCR[42]+0.1497678868178187*GCL[42]-0.29953577363563744*GCC[42]+0.1497678868178187*GTR[33]+0.1497678868178187*GTL[33]+0.29953577363563744*GTC[33]+0.1497678868178187*GCR[33]+0.1497678868178187*GCL[33]+0.29953577363563744*GCC[33]+0.30319611806422586*GTR[31]-0.30319611806422586*GTL[31]-0.30319611806422586*GCR[31]+0.30319611806422586*GCL[31]-0.21544659739277597*GTR[17]-0.21544659739277597*GTL[17]+0.43089319478555205*GTC[17]+0.21544659739277597*GCR[17]+0.21544659739277597*GCL[17]-0.43089319478555205*GCC[17]-0.21544659739277597*GTR[15]+0.21544659739277597*GTL[15]-0.21544659739277597*GCR[15]+0.21544659739277597*GCL[15]+0.15309310892394856*GTR[6]+0.15309310892394856*GTL[6]-0.3061862178478971*GTC[6]+0.15309310892394856*GCR[6]+0.15309310892394856*GCL[6]-0.3061862178478971*GCC[6]; 
  surft1_up[11] = 0.303196118064226*GTR[35]-0.303196118064226*GTL[35]-0.303196118064226*GCR[35]+0.303196118064226*GCL[35]-0.21544659739277597*GTR[25]-0.21544659739277597*GTL[25]+0.43089319478555205*GTC[25]+0.21544659739277597*GCR[25]+0.21544659739277597*GCL[25]-0.43089319478555205*GCC[25]-0.21544659739277597*GTR[19]+0.21544659739277597*GTL[19]-0.21544659739277597*GCR[19]+0.21544659739277597*GCL[19]+0.15309310892394856*GTR[11]+0.15309310892394856*GTL[11]-0.3061862178478971*GTC[11]+0.15309310892394856*GCR[11]+0.15309310892394856*GCL[11]-0.3061862178478971*GCC[11]; 
  surft1_up[12] = -(0.27209908031404895*GTR[41])-0.27209908031404895*GTL[41]-0.7953665424564508*GTC[41]-0.27209908031404895*GCR[41]-0.27209908031404895*GCL[41]-0.7953665424564508*GCC[41]-0.27209908031404895*GTR[36]+0.27209908031404895*GTL[36]+0.27209908031404895*GCR[36]-0.27209908031404895*GCL[36]+0.1933495104806964*GTR[28]-0.1933495104806964*GTL[28]+0.1933495104806964*GCR[28]-0.1933495104806964*GCL[28]+0.1933495104806964*GTR[20]-0.1933495104806964*GTL[20]+0.1933495104806964*GCR[20]-0.1933495104806964*GCL[20]+0.39142450529916156*GTR[16]+0.39142450529916156*GTL[16]+1.14416393856678*GTC[16]-0.39142450529916156*GCR[16]-0.39142450529916156*GCL[16]-1.14416393856678*GCC[16]-0.2781403612330919*GTR[8]+0.2781403612330919*GTL[8]+0.2781403612330919*GCR[8]-0.2781403612330919*GCL[8]-0.2781403612330919*GTR[5]-0.2781403612330919*GTL[5]-0.8130256712967302*GTC[5]-0.2781403612330919*GCR[5]-0.2781403612330919*GCL[5]-0.8130256712967302*GCC[5]+0.19764235376052366*GTR[1]-0.19764235376052366*GTL[1]+0.19764235376052366*GCR[1]-0.19764235376052366*GCL[1]; 
  surft1_up[13] = 0.1750503603816304*GTR[44]+0.1750503603816304*GTL[44]-0.3501007207632608*GTC[44]-0.1750503603816304*GCR[44]-0.1750503603816304*GCL[44]+0.3501007207632608*GCC[44]-0.12438815100070813*GTR[37]+0.12438815100070813*GTL[37]+0.12438815100070813*GCR[37]-0.12438815100070813*GCL[37]-0.12438815100070813*GTR[32]-0.12438815100070813*GTL[32]+0.24877630200141632*GTC[32]-0.12438815100070813*GCR[32]-0.12438815100070813*GCL[32]+0.24877630200141632*GCC[32]+0.0883883476483184*GTR[21]-0.0883883476483184*GTL[21]+0.0883883476483184*GCR[21]-0.0883883476483184*GCL[21]; 
  surft1_up[14] = -(0.27209908031404895*GTR[43])-0.27209908031404895*GTL[43]-0.7953665424564508*GTC[43]-0.27209908031404895*GCR[43]-0.27209908031404895*GCL[43]-0.7953665424564508*GCC[43]-0.27209908031404895*GTR[38]+0.27209908031404895*GTL[38]+0.27209908031404895*GCR[38]-0.27209908031404895*GCL[38]+0.1933495104806964*GTR[30]-0.1933495104806964*GTL[30]+0.1933495104806964*GCR[30]-0.1933495104806964*GCL[30]+0.1933495104806964*GTR[22]-0.1933495104806964*GTL[22]+0.1933495104806964*GCR[22]-0.1933495104806964*GCL[22]+0.39142450529916156*GTR[18]+0.39142450529916156*GTL[18]+1.14416393856678*GTC[18]-0.39142450529916156*GCR[18]-0.39142450529916156*GCL[18]-1.14416393856678*GCC[18]-0.2781403612330919*GTR[10]+0.2781403612330919*GTL[10]+0.2781403612330919*GCR[10]-0.2781403612330919*GCL[10]-0.2781403612330919*GTR[7]-0.2781403612330919*GTL[7]-0.8130256712967302*GTC[7]-0.2781403612330919*GCR[7]-0.2781403612330919*GCL[7]-0.8130256712967302*GCC[7]+0.19764235376052366*GTR[3]-0.19764235376052366*GTL[3]+0.19764235376052366*GCR[3]-0.19764235376052366*GCL[3]; 
  surft1_up[15] = 0.1750503603816304*GTR[46]+0.1750503603816304*GTL[46]-0.3501007207632608*GTC[46]-0.1750503603816304*GCR[46]-0.1750503603816304*GCL[46]+0.3501007207632608*GCC[46]-0.12438815100070813*GTR[39]+0.12438815100070813*GTL[39]+0.12438815100070813*GCR[39]-0.12438815100070813*GCL[39]-0.12438815100070813*GTR[34]-0.12438815100070813*GTL[34]+0.24877630200141632*GTC[34]-0.12438815100070813*GCR[34]-0.12438815100070813*GCL[34]+0.24877630200141632*GCC[34]+0.0883883476483184*GTR[23]-0.0883883476483184*GTL[23]+0.0883883476483184*GCR[23]-0.0883883476483184*GCL[23]; 
  surft1_up[16] = 0.303196118064226*GTR[40]-0.303196118064226*GTL[40]-0.303196118064226*GCR[40]+0.303196118064226*GCL[40]-0.21544659739277597*GTR[27]-0.21544659739277597*GTL[27]+0.43089319478555205*GTC[27]+0.21544659739277597*GCR[27]+0.21544659739277597*GCL[27]-0.43089319478555205*GCC[27]-0.21544659739277597*GTR[24]+0.21544659739277597*GTL[24]-0.21544659739277597*GCR[24]+0.21544659739277597*GCL[24]+0.15309310892394856*GTR[13]+0.15309310892394856*GTL[13]-0.3061862178478971*GTC[13]+0.15309310892394856*GCR[13]+0.15309310892394856*GCL[13]-0.3061862178478971*GCC[13]; 
  surft1_up[17] = 0.303196118064226*GTR[44]-0.303196118064226*GTL[44]-0.303196118064226*GCR[44]+0.303196118064226*GCL[44]-0.21544659739277597*GTR[37]-0.21544659739277597*GTL[37]+0.43089319478555205*GTC[37]+0.21544659739277597*GCR[37]+0.21544659739277597*GCL[37]-0.43089319478555205*GCC[37]-0.21544659739277597*GTR[32]+0.21544659739277597*GTL[32]-0.21544659739277597*GCR[32]+0.21544659739277597*GCL[32]+0.15309310892394856*GTR[21]+0.15309310892394856*GTL[21]-0.3061862178478971*GTC[21]+0.15309310892394856*GCR[21]+0.15309310892394856*GCL[21]-0.3061862178478971*GCC[21]; 
  surft1_up[18] = -(0.27209908031404895*GTR[47])-0.27209908031404895*GTL[47]-0.7953665424564508*GTC[47]-0.27209908031404895*GCR[47]-0.27209908031404895*GCL[47]-0.7953665424564508*GCC[47]-0.27209908031404895*GTR[45]+0.27209908031404895*GTL[45]+0.27209908031404895*GCR[45]-0.27209908031404895*GCL[45]+0.1933495104806964*GTR[42]-0.1933495104806964*GTL[42]+0.1933495104806964*GCR[42]-0.1933495104806964*GCL[42]+0.1933495104806964*GTR[33]-0.1933495104806964*GTL[33]+0.1933495104806964*GCR[33]-0.1933495104806964*GCL[33]+0.3914245052991616*GTR[31]+0.3914245052991616*GTL[31]+1.1441639385667801*GTC[31]-0.3914245052991616*GCR[31]-0.3914245052991616*GCL[31]-1.1441639385667801*GCC[31]-0.2781403612330919*GTR[17]+0.2781403612330919*GTL[17]+0.2781403612330919*GCR[17]-0.2781403612330919*GCL[17]-0.2781403612330919*GTR[15]-0.2781403612330919*GTL[15]-0.8130256712967302*GTC[15]-0.2781403612330919*GCR[15]-0.2781403612330919*GCL[15]-0.8130256712967302*GCC[15]+0.19764235376052364*GTR[6]-0.19764235376052364*GTL[6]+0.19764235376052364*GCR[6]-0.19764235376052364*GCL[6]; 
  surft1_up[19] = 0.303196118064226*GTR[46]-0.303196118064226*GTL[46]-0.303196118064226*GCR[46]+0.303196118064226*GCL[46]-0.21544659739277597*GTR[39]-0.21544659739277597*GTL[39]+0.43089319478555205*GTC[39]+0.21544659739277597*GCR[39]+0.21544659739277597*GCL[39]-0.43089319478555205*GCC[39]-0.21544659739277597*GTR[34]+0.21544659739277597*GTL[34]-0.21544659739277597*GCR[34]+0.21544659739277597*GCL[34]+0.15309310892394856*GTR[23]+0.15309310892394856*GTL[23]-0.3061862178478971*GTC[23]+0.15309310892394856*GCR[23]+0.15309310892394856*GCL[23]-0.3061862178478971*GCC[23]; 

  surft2_lo[0] = 0.34587411908091625*GCL[12]+0.34587411908091625*GCC[12]+0.49755260400283263*GCL[2]-0.49755260400283263*GCC[2]+0.3535533905932737*GCL[0]+0.3535533905932737*GCC[0]; 
  surft2_lo[1] = 0.34587411908091625*GCL[20]+0.34587411908091625*GCC[20]+0.49755260400283263*GCL[5]-0.49755260400283263*GCC[5]+0.3535533905932737*GCL[1]+0.3535533905932737*GCC[1]; 
  surft2_lo[2] = 0.34587411908091625*GCL[22]+0.34587411908091625*GCC[22]+0.49755260400283263*GCL[7]-0.49755260400283263*GCC[7]+0.3535533905932737*GCL[3]+0.3535533905932737*GCC[3]; 
  surft2_lo[3] = 0.34587411908091625*GCL[26]+0.34587411908091625*GCC[26]+0.49755260400283263*GCL[9]-0.49755260400283263*GCC[9]+0.3535533905932737*GCL[4]+0.3535533905932737*GCC[4]; 
  surft2_lo[4] = 0.34587411908091625*GCL[33]+0.34587411908091625*GCC[33]+0.49755260400283263*GCL[15]-0.49755260400283263*GCC[15]+0.3535533905932737*GCL[6]+0.3535533905932737*GCC[6]; 
  surft2_lo[5] = 0.34587411908091625*GCL[36]+0.34587411908091625*GCC[36]+0.49755260400283263*GCL[16]-0.49755260400283263*GCC[16]+0.3535533905932737*GCL[8]+0.3535533905932737*GCC[8]; 
  surft2_lo[6] = 0.34587411908091625*GCL[38]+0.34587411908091625*GCC[38]+0.49755260400283263*GCL[18]-0.49755260400283263*GCC[18]+0.3535533905932737*GCL[10]+0.3535533905932737*GCC[10]; 
  surft2_lo[7] = 0.49755260400283263*GCL[19]-0.49755260400283263*GCC[19]+0.3535533905932737*GCL[11]+0.3535533905932737*GCC[11]; 
  surft2_lo[8] = 0.49755260400283263*GCL[24]-0.49755260400283263*GCC[24]+0.3535533905932737*GCL[13]+0.3535533905932737*GCC[13]; 
  surft2_lo[9] = 0.49755260400283263*GCL[29]-0.49755260400283263*GCC[29]+0.3535533905932737*GCL[14]+0.3535533905932737*GCC[14]; 
  surft2_lo[10] = 0.34587411908091625*GCL[45]+0.34587411908091625*GCC[45]+0.49755260400283263*GCL[31]-0.49755260400283263*GCC[31]+0.3535533905932737*GCL[17]+0.3535533905932737*GCC[17]; 
  surft2_lo[11] = 0.49755260400283263*GCL[32]-0.49755260400283263*GCC[32]+0.3535533905932737*GCL[21]+0.3535533905932737*GCC[21]; 
  surft2_lo[12] = 0.49755260400283263*GCL[34]-0.49755260400283263*GCC[34]+0.3535533905932737*GCL[23]+0.3535533905932737*GCC[23]; 
  surft2_lo[13] = 0.49755260400283263*GCL[35]-0.49755260400283263*GCC[35]+0.3535533905932737*GCL[25]+0.3535533905932737*GCC[25]; 
  surft2_lo[14] = 0.49755260400283263*GCL[40]-0.49755260400283263*GCC[40]+0.3535533905932737*GCL[27]+0.3535533905932737*GCC[27]; 
  surft2_lo[15] = 0.49755260400283263*GCL[41]-0.49755260400283263*GCC[41]+0.3535533905932737*GCL[28]+0.3535533905932737*GCC[28]; 
  surft2_lo[16] = 0.49755260400283263*GCL[43]-0.49755260400283263*GCC[43]+0.3535533905932737*GCL[30]+0.3535533905932737*GCC[30]; 
  surft2_lo[17] = 0.49755260400283263*GCL[44]-0.49755260400283263*GCC[44]+0.3535533905932737*GCL[37]+0.3535533905932737*GCC[37]; 
  surft2_lo[18] = 0.49755260400283263*GCL[46]-0.49755260400283263*GCC[46]+0.3535533905932737*GCL[39]+0.3535533905932737*GCC[39]; 
  surft2_lo[19] = 0.49755260400283263*GCL[47]-0.49755260400283263*GCC[47]+0.3535533905932737*GCL[42]+0.3535533905932737*GCC[42]; 
  surft2_up[0] = 0.34587411908091625*GCR[12]+0.34587411908091625*GCC[12]-0.49755260400283263*GCR[2]+0.49755260400283263*GCC[2]+0.3535533905932737*GCR[0]+0.3535533905932737*GCC[0]; 
  surft2_up[1] = 0.34587411908091625*GCR[20]+0.34587411908091625*GCC[20]-0.49755260400283263*GCR[5]+0.49755260400283263*GCC[5]+0.3535533905932737*GCR[1]+0.3535533905932737*GCC[1]; 
  surft2_up[2] = 0.34587411908091625*GCR[22]+0.34587411908091625*GCC[22]-0.49755260400283263*GCR[7]+0.49755260400283263*GCC[7]+0.3535533905932737*GCR[3]+0.3535533905932737*GCC[3]; 
  surft2_up[3] = 0.34587411908091625*GCR[26]+0.34587411908091625*GCC[26]-0.49755260400283263*GCR[9]+0.49755260400283263*GCC[9]+0.3535533905932737*GCR[4]+0.3535533905932737*GCC[4]; 
  surft2_up[4] = 0.34587411908091625*GCR[33]+0.34587411908091625*GCC[33]-0.49755260400283263*GCR[15]+0.49755260400283263*GCC[15]+0.3535533905932737*GCR[6]+0.3535533905932737*GCC[6]; 
  surft2_up[5] = 0.34587411908091625*GCR[36]+0.34587411908091625*GCC[36]-0.49755260400283263*GCR[16]+0.49755260400283263*GCC[16]+0.3535533905932737*GCR[8]+0.3535533905932737*GCC[8]; 
  surft2_up[6] = 0.34587411908091625*GCR[38]+0.34587411908091625*GCC[38]-0.49755260400283263*GCR[18]+0.49755260400283263*GCC[18]+0.3535533905932737*GCR[10]+0.3535533905932737*GCC[10]; 
  surft2_up[7] = -(0.49755260400283263*GCR[19])+0.49755260400283263*GCC[19]+0.3535533905932737*GCR[11]+0.3535533905932737*GCC[11]; 
  surft2_up[8] = -(0.49755260400283263*GCR[24])+0.49755260400283263*GCC[24]+0.3535533905932737*GCR[13]+0.3535533905932737*GCC[13]; 
  surft2_up[9] = -(0.49755260400283263*GCR[29])+0.49755260400283263*GCC[29]+0.3535533905932737*GCR[14]+0.3535533905932737*GCC[14]; 
  surft2_up[10] = 0.34587411908091625*GCR[45]+0.34587411908091625*GCC[45]-0.49755260400283263*GCR[31]+0.49755260400283263*GCC[31]+0.3535533905932737*GCR[17]+0.3535533905932737*GCC[17]; 
  surft2_up[11] = -(0.49755260400283263*GCR[32])+0.49755260400283263*GCC[32]+0.3535533905932737*GCR[21]+0.3535533905932737*GCC[21]; 
  surft2_up[12] = -(0.49755260400283263*GCR[34])+0.49755260400283263*GCC[34]+0.3535533905932737*GCR[23]+0.3535533905932737*GCC[23]; 
  surft2_up[13] = -(0.49755260400283263*GCR[35])+0.49755260400283263*GCC[35]+0.3535533905932737*GCR[25]+0.3535533905932737*GCC[25]; 
  surft2_up[14] = -(0.49755260400283263*GCR[40])+0.49755260400283263*GCC[40]+0.3535533905932737*GCR[27]+0.3535533905932737*GCC[27]; 
  surft2_up[15] = -(0.49755260400283263*GCR[41])+0.49755260400283263*GCC[41]+0.3535533905932737*GCR[28]+0.3535533905932737*GCC[28]; 
  surft2_up[16] = -(0.49755260400283263*GCR[43])+0.49755260400283263*GCC[43]+0.3535533905932737*GCR[30]+0.3535533905932737*GCC[30]; 
  surft2_up[17] = -(0.49755260400283263*GCR[44])+0.49755260400283263*GCC[44]+0.3535533905932737*GCR[37]+0.3535533905932737*GCC[37]; 
  surft2_up[18] = -(0.49755260400283263*GCR[46])+0.49755260400283263*GCC[46]+0.3535533905932737*GCR[39]+0.3535533905932737*GCC[39]; 
  surft2_up[19] = -(0.49755260400283263*GCR[47])+0.49755260400283263*GCC[47]+0.3535533905932737*GCR[42]+0.3535533905932737*GCC[42]; 

  vol[9] = 3.0*GCC[0]; 
  vol[16] = 3.0*GCC[1]; 
  vol[18] = 3.0*GCC[3]; 
  vol[26] = 6.7082039324993685*GCC[2]; 
  vol[29] = 6.7082039324993685*GCC[4]; 
  vol[31] = 3.0*GCC[6]; 
  vol[35] = 3.0*GCC[11]; 
  vol[36] = 6.708203932499369*GCC[5]; 
  vol[38] = 6.708203932499369*GCC[7]; 
  vol[40] = 3.0*GCC[13]; 
  vol[41] = 6.708203932499369*GCC[8]; 
  vol[43] = 6.708203932499369*GCC[10]; 
  vol[44] = 3.0*GCC[21]; 
  vol[45] = 6.7082039324993685*GCC[15]; 
  vol[46] = 3.0*GCC[23]; 
  vol[47] = 6.7082039324993685*GCC[17]; 

  out[0] = (vol[0]+0.7071067811865475*surft1_up[0]-0.7071067811865475*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[1] = (vol[1]+0.7071067811865475*surft1_up[1]-0.7071067811865475*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[2] = (vol[2]+0.7071067811865475*surft1_up[2]-0.7071067811865475*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[3] = (vol[3]+0.7071067811865475*surft1_up[3]-0.7071067811865475*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[4] = (vol[4]-1.224744871391589*surft2_up[0]+1.224744871391589*surft2_lo[0]+1.224744871391589*surft1_up[0]+1.224744871391589*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[5] = (vol[5]+0.7071067811865475*surft1_up[4]-0.7071067811865475*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[6] = (vol[6]+0.7071067811865475*surft1_up[5]-0.7071067811865475*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[7] = (vol[7]+0.7071067811865475*surft1_up[6]-0.7071067811865475*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[8] = (vol[8]-1.224744871391589*surft2_up[1]+1.224744871391589*surft2_lo[1]+1.224744871391589*surft1_up[1]+1.224744871391589*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[9] = (vol[9]+1.224744871391589*surft1_up[2]+1.224744871391589*surft1_lo[2]-2.1213203435596424*surft2_up[0]-2.1213203435596424*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[10] = (vol[10]+1.224744871391589*surft1_up[3]+1.224744871391589*surft1_lo[3]-1.224744871391589*surft2_up[2]+1.224744871391589*surft2_lo[2])*dv1_sq*gamma_avg; 
  out[11] = (vol[11]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[12] = (vol[12]+0.7071067811865475*surft1_up[8]-0.7071067811865475*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[13] = (vol[13]+0.7071067811865475*surft1_up[9]-0.7071067811865475*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]-2.7386127875258306*surft2_up[3]+2.7386127875258306*surft2_lo[3]+1.5811388300841895*surft1_up[0]-1.5811388300841895*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+0.7071067811865475*surft1_up[10]-0.7071067811865475*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[16] = (vol[16]+1.224744871391589*surft1_up[4]+1.224744871391589*surft1_lo[4]-2.1213203435596424*surft2_up[1]-2.1213203435596424*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[17] = (vol[17]+1.224744871391589*surft1_up[5]+1.224744871391589*surft1_lo[5]-1.224744871391589*surft2_up[4]+1.224744871391589*surft2_lo[4])*dv1_sq*gamma_avg; 
  out[18] = (vol[18]+1.224744871391589*surft1_up[6]+1.224744871391589*surft1_lo[6]-2.1213203435596424*surft2_up[2]-2.1213203435596424*surft2_lo[2])*dv1_sq*gamma_avg; 
  out[19] = (vol[19]+0.7071067811865475*surft1_up[11]-0.7071067811865475*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[20] = (vol[20]+0.7071067811865475*surft1_up[12]-0.7071067811865475*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[21] = (vol[21]+0.7071067811865475*surft1_up[13]-0.7071067811865475*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[22] = (vol[22]+0.7071067811865475*surft1_up[14]-0.7071067811865475*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[23] = (vol[23]+0.7071067811865475*surft1_up[15]-0.7071067811865475*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[24] = (vol[24]+0.7071067811865475*surft1_up[16]-0.7071067811865475*surft1_lo[16])*dv1_sq*gamma_avg; 
  out[25] = (vol[25]-1.224744871391589*surft2_up[7]+1.224744871391589*surft2_lo[7]+1.224744871391589*surft1_up[7]+1.224744871391589*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[26] = (vol[26]+1.224744871391589*surft1_up[8]+1.224744871391589*surft1_lo[8]-2.7386127875258306*surft2_up[0]+2.7386127875258306*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[27] = (vol[27]+1.224744871391589*surft1_up[9]+1.224744871391589*surft1_lo[9]-1.224744871391589*surft2_up[8]+1.224744871391589*surft2_lo[8])*dv1_sq*gamma_avg; 
  out[28] = (vol[28]-2.7386127875258306*surft2_up[5]+2.7386127875258306*surft2_lo[5]+1.5811388300841898*surft1_up[1]-1.5811388300841898*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[29] = (vol[29]-4.743416490252569*surft2_up[3]-4.743416490252569*surft2_lo[3]+1.5811388300841898*surft1_up[2]-1.5811388300841898*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[30] = (vol[30]-2.7386127875258306*surft2_up[6]+2.7386127875258306*surft2_lo[6]+1.5811388300841898*surft1_up[3]-1.5811388300841898*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[31] = (vol[31]+1.224744871391589*surft1_up[10]+1.224744871391589*surft1_lo[10]-2.1213203435596424*surft2_up[4]-2.1213203435596424*surft2_lo[4])*dv1_sq*gamma_avg; 
  out[32] = (vol[32]+0.7071067811865475*surft1_up[17]-0.7071067811865475*surft1_lo[17])*dv1_sq*gamma_avg; 
  out[33] = (vol[33]+0.7071067811865475*surft1_up[18]-0.7071067811865475*surft1_lo[18])*dv1_sq*gamma_avg; 
  out[34] = (vol[34]+0.7071067811865475*surft1_up[19]-0.7071067811865475*surft1_lo[19])*dv1_sq*gamma_avg; 
  out[35] = (vol[35]+1.224744871391589*surft1_up[11]+1.224744871391589*surft1_lo[11]-2.1213203435596424*surft2_up[7]-2.1213203435596424*surft2_lo[7])*dv1_sq*gamma_avg; 
  out[36] = (vol[36]+1.224744871391589*surft1_up[12]+1.224744871391589*surft1_lo[12]-2.7386127875258306*surft2_up[1]+2.7386127875258306*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[37] = (vol[37]+1.224744871391589*surft1_up[13]+1.224744871391589*surft1_lo[13]-1.224744871391589*surft2_up[11]+1.224744871391589*surft2_lo[11])*dv1_sq*gamma_avg; 
  out[38] = (vol[38]+1.224744871391589*surft1_up[14]+1.224744871391589*surft1_lo[14]-2.7386127875258306*surft2_up[2]+2.7386127875258306*surft2_lo[2])*dv1_sq*gamma_avg; 
  out[39] = (vol[39]+1.224744871391589*surft1_up[15]+1.224744871391589*surft1_lo[15]-1.224744871391589*surft2_up[12]+1.224744871391589*surft2_lo[12])*dv1_sq*gamma_avg; 
  out[40] = (vol[40]+1.224744871391589*surft1_up[16]+1.224744871391589*surft1_lo[16]-2.1213203435596424*surft2_up[8]-2.1213203435596424*surft2_lo[8])*dv1_sq*gamma_avg; 
  out[41] = (vol[41]-4.743416490252569*surft2_up[5]-4.743416490252569*surft2_lo[5]+1.5811388300841895*surft1_up[4]-1.5811388300841895*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[42] = (vol[42]-2.7386127875258306*surft2_up[10]+2.7386127875258306*surft2_lo[10]+1.5811388300841895*surft1_up[5]-1.5811388300841895*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[43] = (vol[43]-4.743416490252569*surft2_up[6]-4.743416490252569*surft2_lo[6]+1.5811388300841895*surft1_up[6]-1.5811388300841895*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[44] = (vol[44]+1.224744871391589*surft1_up[17]+1.224744871391589*surft1_lo[17]-2.1213203435596424*surft2_up[11]-2.1213203435596424*surft2_lo[11])*dv1_sq*gamma_avg; 
  out[45] = (vol[45]+1.224744871391589*surft1_up[18]+1.224744871391589*surft1_lo[18]-2.7386127875258306*surft2_up[4]+2.7386127875258306*surft2_lo[4])*dv1_sq*gamma_avg; 
  out[46] = (vol[46]+1.224744871391589*surft1_up[19]+1.224744871391589*surft1_lo[19]-2.1213203435596424*surft2_up[12]-2.1213203435596424*surft2_lo[12])*dv1_sq*gamma_avg; 
  out[47] = (vol[47]-4.743416490252569*surft2_up[10]-4.743416490252569*surft2_lo[10]+1.5811388300841898*surft1_up[10]-1.5811388300841898*surft1_lo[10])*dv1_sq*gamma_avg; 
} 