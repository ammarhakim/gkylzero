#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p2_upvx_invy(const double *dxv, double gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // diff_coeff: Output array for diffusion tensor. 

  double dv1_pv1 = 2.0/dxv[1]; 
  double dv1_pv2 = 2.0/dxv[2]; 
  double dv1_sq = 4.0/dxv[1]/dxv[2]; 
 
  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GTL = fpo_g_stencil[2]; 
  const double* GBC = fpo_g_stencil[3]; 
  const double* GCC = fpo_g_stencil[4]; 
  const double* GTC = fpo_g_stencil[5]; 

  const double* g_surf_CL = fpo_g_surf_stencil[1]; 
  const double* g_surf_CL_pv2 = &g_surf_CL[20]; 
  const double* g_surf_CC = fpo_g_surf_stencil[4]; 
  const double* g_surf_CC_pv2 = &g_surf_CC[20]; 
  
  const double* g_surf_CC_pv1 = &g_surf_CC[0]; 
  const double* dgdpv1_surf_CC_pv2 = &fpo_dgdv_surf[60]; 
  const double* dgdpv2_surf_CC_pv1 = &fpo_dgdv_surf[20]; 
  const double* dgdpv1_surf_CC_pv1 = &fpo_dgdv_surf[0]; 
  
  double surft1_upper[20], surft1_lower[20]; 
  double surft2_upper[20], surft2_lower[20]; 
  
  double *diff_coeff_vxvy = &diff_coeff[48]; 
  double *diff_coeff_vxvz = &diff_coeff[96]; 
  double *diff_coeff_vyvx = &diff_coeff[144]; 
  double *diff_coeff_vyvz = &diff_coeff[240]; 
  double *diff_coeff_vzvx = &diff_coeff[288]; 
  double *diff_coeff_vzvy = &diff_coeff[336]; 
  
  double *out = diff_coeff_vxvy; 
  
  surft1_upper[0] = dgdpv2_surf_CC_pv1[0]/dv1_pv2; 
  surft1_upper[1] = dgdpv2_surf_CC_pv1[1]/dv1_pv2; 
  surft1_upper[2] = dgdpv2_surf_CC_pv1[2]/dv1_pv2; 
  surft1_upper[3] = dgdpv2_surf_CC_pv1[3]/dv1_pv2; 
  surft1_upper[4] = dgdpv2_surf_CC_pv1[4]/dv1_pv2; 
  surft1_upper[5] = dgdpv2_surf_CC_pv1[5]/dv1_pv2; 
  surft1_upper[6] = dgdpv2_surf_CC_pv1[6]/dv1_pv2; 
  surft1_upper[7] = dgdpv2_surf_CC_pv1[7]/dv1_pv2; 
  surft1_upper[8] = dgdpv2_surf_CC_pv1[8]/dv1_pv2; 
  surft1_upper[9] = dgdpv2_surf_CC_pv1[9]/dv1_pv2; 
  surft1_upper[10] = dgdpv2_surf_CC_pv1[10]/dv1_pv2; 
  surft1_upper[11] = dgdpv2_surf_CC_pv1[11]/dv1_pv2; 
  surft1_upper[12] = dgdpv2_surf_CC_pv1[12]/dv1_pv2; 
  surft1_upper[13] = dgdpv2_surf_CC_pv1[13]/dv1_pv2; 
  surft1_upper[14] = dgdpv2_surf_CC_pv1[14]/dv1_pv2; 
  surft1_upper[15] = dgdpv2_surf_CC_pv1[15]/dv1_pv2; 
  surft1_upper[16] = dgdpv2_surf_CC_pv1[16]/dv1_pv2; 
  surft1_upper[17] = dgdpv2_surf_CC_pv1[17]/dv1_pv2; 
  surft1_upper[18] = dgdpv2_surf_CC_pv1[18]/dv1_pv2; 
  surft1_upper[19] = dgdpv2_surf_CC_pv1[19]/dv1_pv2; 
  surft1_lower[0] = 0.12168640803947765*GTL[24]-0.12168640803947765*GTC[24]-0.12168640803947765*GBL[24]+0.12168640803947765*GBC[24]-0.12168640803947765*GTL[22]-0.12168640803947765*GTC[22]+0.2433728160789553*GCL[22]+0.2433728160789553*GCC[22]-0.12168640803947765*GBL[22]-0.12168640803947765*GBC[22]+0.08646852977022904*GTL[13]+0.08646852977022904*GTC[13]-0.08646852977022904*GBL[13]-0.08646852977022904*GBC[13]+0.08646852977022904*GTL[12]+0.08646852977022904*GTC[12]-0.08646852977022904*GBL[12]-0.08646852977022904*GBC[12]-0.1750503603816304*GTL[7]+0.1750503603816304*GTC[7]+0.3501007207632608*GCL[7]-0.3501007207632608*GCC[7]-0.1750503603816304*GBL[7]+0.1750503603816304*GBC[7]-0.12438815100070813*GTL[3]-0.12438815100070813*GTC[3]+0.24877630200141632*GCL[3]+0.24877630200141632*GCC[3]-0.12438815100070813*GBL[3]-0.12438815100070813*GBC[3]+0.12438815100070813*GTL[2]-0.12438815100070813*GTC[2]-0.12438815100070813*GBL[2]+0.12438815100070813*GBC[2]+0.0883883476483184*GTL[0]+0.0883883476483184*GTC[0]-0.0883883476483184*GBL[0]-0.0883883476483184*GBC[0]; 
  surft1_lower[1] = 0.12168640803947765*GTL[34]-0.12168640803947765*GTC[34]-0.12168640803947765*GBL[34]+0.12168640803947765*GBC[34]-0.12168640803947765*GTL[33]-0.12168640803947765*GTC[33]+0.2433728160789553*GCL[33]+0.2433728160789553*GCC[33]-0.12168640803947765*GBL[33]-0.12168640803947765*GBC[33]+0.08646852977022904*GTL[23]+0.08646852977022904*GTC[23]-0.08646852977022904*GBL[23]-0.08646852977022904*GBC[23]+0.08646852977022904*GTL[20]+0.08646852977022904*GTC[20]-0.08646852977022904*GBL[20]-0.08646852977022904*GBC[20]-0.1750503603816304*GTL[15]+0.1750503603816304*GTC[15]+0.3501007207632608*GCL[15]-0.3501007207632608*GCC[15]-0.1750503603816304*GBL[15]+0.1750503603816304*GBC[15]-0.12438815100070813*GTL[6]-0.12438815100070813*GTC[6]+0.24877630200141632*GCL[6]+0.24877630200141632*GCC[6]-0.12438815100070813*GBL[6]-0.12438815100070813*GBC[6]+0.12438815100070813*GTL[5]-0.12438815100070813*GTC[5]-0.12438815100070813*GBL[5]+0.12438815100070813*GBC[5]+0.0883883476483184*GTL[1]+0.0883883476483184*GTC[1]-0.0883883476483184*GBL[1]-0.0883883476483184*GBC[1]; 
  surft1_lower[2] = 0.2107670413149332*GTL[24]-0.2107670413149332*GTC[24]+0.4215340826298664*GCL[24]-0.4215340826298664*GCC[24]+0.2107670413149332*GBL[24]-0.2107670413149332*GBC[24]-0.2107670413149332*GTL[22]-0.2107670413149332*GTC[22]+0.2107670413149332*GBL[22]+0.2107670413149332*GBC[22]+0.1497678868178187*GTL[13]+0.1497678868178187*GTC[13]+0.29953577363563744*GCL[13]+0.29953577363563744*GCC[13]+0.1497678868178187*GBL[13]+0.1497678868178187*GBC[13]+0.1497678868178187*GTL[12]+0.1497678868178187*GTC[12]-0.29953577363563744*GCL[12]-0.29953577363563744*GCC[12]+0.1497678868178187*GBL[12]+0.1497678868178187*GBC[12]-0.30319611806422586*GTL[7]+0.30319611806422586*GTC[7]+0.30319611806422586*GBL[7]-0.30319611806422586*GBC[7]-0.21544659739277597*GTL[3]-0.21544659739277597*GTC[3]+0.21544659739277597*GBL[3]+0.21544659739277597*GBC[3]+0.21544659739277597*GTL[2]-0.21544659739277597*GTC[2]-0.43089319478555205*GCL[2]+0.43089319478555205*GCC[2]+0.21544659739277597*GBL[2]-0.21544659739277597*GBC[2]+0.15309310892394856*GTL[0]+0.15309310892394856*GTC[0]-0.3061862178478971*GCL[0]-0.3061862178478971*GCC[0]+0.15309310892394856*GBL[0]+0.15309310892394856*GBC[0]; 
  surft1_lower[3] = 0.12168640803947765*GTL[40]-0.12168640803947765*GTC[40]-0.12168640803947765*GBL[40]+0.12168640803947765*GBC[40]-0.12168640803947765*GTL[38]-0.12168640803947765*GTC[38]+0.2433728160789553*GCL[38]+0.2433728160789553*GCC[38]-0.12168640803947765*GBL[38]-0.12168640803947765*GBC[38]+0.08646852977022904*GTL[27]+0.08646852977022904*GTC[27]-0.08646852977022904*GBL[27]-0.08646852977022904*GBC[27]+0.08646852977022904*GTL[26]+0.08646852977022904*GTC[26]-0.08646852977022904*GBL[26]-0.08646852977022904*GBC[26]-0.1750503603816304*GTL[18]+0.1750503603816304*GTC[18]+0.3501007207632608*GCL[18]-0.3501007207632608*GCC[18]-0.1750503603816304*GBL[18]+0.1750503603816304*GBC[18]-0.12438815100070813*GTL[10]-0.12438815100070813*GTC[10]+0.24877630200141632*GCL[10]+0.24877630200141632*GCC[10]-0.12438815100070813*GBL[10]-0.12438815100070813*GBC[10]+0.12438815100070813*GTL[9]-0.12438815100070813*GTC[9]-0.12438815100070813*GBL[9]+0.12438815100070813*GBC[9]+0.0883883476483184*GTL[4]+0.0883883476483184*GTC[4]-0.0883883476483184*GBL[4]-0.0883883476483184*GBC[4]; 
  surft1_lower[4] = 0.21076704131493318*GTL[34]-0.21076704131493318*GTC[34]+0.42153408262986636*GCL[34]-0.42153408262986636*GCC[34]+0.21076704131493318*GBL[34]-0.21076704131493318*GBC[34]-0.21076704131493318*GTL[33]-0.21076704131493318*GTC[33]+0.21076704131493318*GBL[33]+0.21076704131493318*GBC[33]+0.1497678868178187*GTL[23]+0.1497678868178187*GTC[23]+0.29953577363563744*GCL[23]+0.29953577363563744*GCC[23]+0.1497678868178187*GBL[23]+0.1497678868178187*GBC[23]+0.1497678868178187*GTL[20]+0.1497678868178187*GTC[20]-0.29953577363563744*GCL[20]-0.29953577363563744*GCC[20]+0.1497678868178187*GBL[20]+0.1497678868178187*GBC[20]-0.30319611806422586*GTL[15]+0.30319611806422586*GTC[15]+0.30319611806422586*GBL[15]-0.30319611806422586*GBC[15]-0.21544659739277597*GTL[6]-0.21544659739277597*GTC[6]+0.21544659739277597*GBL[6]+0.21544659739277597*GBC[6]+0.21544659739277597*GTL[5]-0.21544659739277597*GTC[5]-0.43089319478555205*GCL[5]+0.43089319478555205*GCC[5]+0.21544659739277597*GBL[5]-0.21544659739277597*GBC[5]+0.15309310892394856*GTL[1]+0.15309310892394856*GTC[1]-0.3061862178478971*GCL[1]-0.3061862178478971*GCC[1]+0.15309310892394856*GBL[1]+0.15309310892394856*GBC[1]; 
  surft1_lower[5] = 0.12168640803947765*GTL[46]-0.12168640803947765*GTC[46]-0.12168640803947765*GBL[46]+0.12168640803947765*GBC[46]-0.12168640803947765*GTL[45]-0.12168640803947765*GTC[45]+0.2433728160789553*GCL[45]+0.2433728160789553*GCC[45]-0.12168640803947765*GBL[45]-0.12168640803947765*GBC[45]+0.08646852977022904*GTL[39]+0.08646852977022904*GTC[39]-0.08646852977022904*GBL[39]-0.08646852977022904*GBC[39]+0.08646852977022904*GTL[36]+0.08646852977022904*GTC[36]-0.08646852977022904*GBL[36]-0.08646852977022904*GBC[36]-0.1750503603816304*GTL[31]+0.1750503603816304*GTC[31]+0.3501007207632608*GCL[31]-0.3501007207632608*GCC[31]-0.1750503603816304*GBL[31]+0.1750503603816304*GBC[31]-0.12438815100070813*GTL[17]-0.12438815100070813*GTC[17]+0.24877630200141632*GCL[17]+0.24877630200141632*GCC[17]-0.12438815100070813*GBL[17]-0.12438815100070813*GBC[17]+0.12438815100070813*GTL[16]-0.12438815100070813*GTC[16]-0.12438815100070813*GBL[16]+0.12438815100070813*GBC[16]+0.0883883476483184*GTL[8]+0.0883883476483184*GTC[8]-0.0883883476483184*GBL[8]-0.0883883476483184*GBC[8]; 
  surft1_lower[6] = 0.21076704131493318*GTL[40]-0.21076704131493318*GTC[40]+0.42153408262986636*GCL[40]-0.42153408262986636*GCC[40]+0.21076704131493318*GBL[40]-0.21076704131493318*GBC[40]-0.21076704131493318*GTL[38]-0.21076704131493318*GTC[38]+0.21076704131493318*GBL[38]+0.21076704131493318*GBC[38]+0.1497678868178187*GTL[27]+0.1497678868178187*GTC[27]+0.29953577363563744*GCL[27]+0.29953577363563744*GCC[27]+0.1497678868178187*GBL[27]+0.1497678868178187*GBC[27]+0.1497678868178187*GTL[26]+0.1497678868178187*GTC[26]-0.29953577363563744*GCL[26]-0.29953577363563744*GCC[26]+0.1497678868178187*GBL[26]+0.1497678868178187*GBC[26]-0.30319611806422586*GTL[18]+0.30319611806422586*GTC[18]+0.30319611806422586*GBL[18]-0.30319611806422586*GBC[18]-0.21544659739277597*GTL[10]-0.21544659739277597*GTC[10]+0.21544659739277597*GBL[10]+0.21544659739277597*GBC[10]+0.21544659739277597*GTL[9]-0.21544659739277597*GTC[9]-0.43089319478555205*GCL[9]+0.43089319478555205*GCC[9]+0.21544659739277597*GBL[9]-0.21544659739277597*GBC[9]+0.15309310892394856*GTL[4]+0.15309310892394856*GTC[4]-0.3061862178478971*GCL[4]-0.3061862178478971*GCC[4]+0.15309310892394856*GBL[4]+0.15309310892394856*GBC[4]; 
  surft1_lower[7] = -(0.1750503603816304*GTL[32])+0.1750503603816304*GTC[32]+0.3501007207632608*GCL[32]-0.3501007207632608*GCC[32]-0.1750503603816304*GBL[32]+0.1750503603816304*GBC[32]-0.12438815100070813*GTL[21]-0.12438815100070813*GTC[21]+0.24877630200141632*GCL[21]+0.24877630200141632*GCC[21]-0.12438815100070813*GBL[21]-0.12438815100070813*GBC[21]+0.12438815100070813*GTL[19]-0.12438815100070813*GTC[19]-0.12438815100070813*GBL[19]+0.12438815100070813*GBC[19]+0.0883883476483184*GTL[11]+0.0883883476483184*GTC[11]-0.0883883476483184*GBL[11]-0.0883883476483184*GBC[11]; 
  surft1_lower[8] = 0.27209908031404895*GTL[24]-0.27209908031404895*GTC[24]-0.27209908031404895*GBL[24]+0.27209908031404895*GBC[24]-0.27209908031404895*GTL[22]-0.27209908031404895*GTC[22]-0.7953665424564508*GCL[22]-0.7953665424564508*GCC[22]-0.27209908031404895*GBL[22]-0.27209908031404895*GBC[22]+0.1933495104806964*GTL[13]+0.1933495104806964*GTC[13]-0.1933495104806964*GBL[13]-0.1933495104806964*GBC[13]+0.1933495104806964*GTL[12]+0.1933495104806964*GTC[12]-0.1933495104806964*GBL[12]-0.1933495104806964*GBC[12]-0.3914245052991616*GTL[7]+0.3914245052991616*GTC[7]-1.1441639385667801*GCL[7]+1.1441639385667801*GCC[7]-0.3914245052991616*GBL[7]+0.3914245052991616*GBC[7]-0.2781403612330919*GTL[3]-0.2781403612330919*GTC[3]-0.8130256712967302*GCL[3]-0.8130256712967302*GCC[3]-0.2781403612330919*GBL[3]-0.2781403612330919*GBC[3]+0.2781403612330919*GTL[2]-0.2781403612330919*GTC[2]-0.2781403612330919*GBL[2]+0.2781403612330919*GBC[2]+0.19764235376052364*GTL[0]+0.19764235376052364*GTC[0]-0.19764235376052364*GBL[0]-0.19764235376052364*GBC[0]; 
  surft1_lower[9] = -(0.1750503603816304*GTL[43])+0.1750503603816304*GTC[43]+0.3501007207632608*GCL[43]-0.3501007207632608*GCC[43]-0.1750503603816304*GBL[43]+0.1750503603816304*GBC[43]-0.12438815100070813*GTL[30]-0.12438815100070813*GTC[30]+0.24877630200141632*GCL[30]+0.24877630200141632*GCC[30]-0.12438815100070813*GBL[30]-0.12438815100070813*GBC[30]+0.12438815100070813*GTL[29]-0.12438815100070813*GTC[29]-0.12438815100070813*GBL[29]+0.12438815100070813*GBC[29]+0.0883883476483184*GTL[14]+0.0883883476483184*GTC[14]-0.0883883476483184*GBL[14]-0.0883883476483184*GBC[14]; 
  surft1_lower[10] = 0.2107670413149332*GTL[46]-0.2107670413149332*GTC[46]+0.4215340826298664*GCL[46]-0.4215340826298664*GCC[46]+0.2107670413149332*GBL[46]-0.2107670413149332*GBC[46]-0.2107670413149332*GTL[45]-0.2107670413149332*GTC[45]+0.2107670413149332*GBL[45]+0.2107670413149332*GBC[45]+0.1497678868178187*GTL[39]+0.1497678868178187*GTC[39]+0.29953577363563744*GCL[39]+0.29953577363563744*GCC[39]+0.1497678868178187*GBL[39]+0.1497678868178187*GBC[39]+0.1497678868178187*GTL[36]+0.1497678868178187*GTC[36]-0.29953577363563744*GCL[36]-0.29953577363563744*GCC[36]+0.1497678868178187*GBL[36]+0.1497678868178187*GBC[36]-0.30319611806422586*GTL[31]+0.30319611806422586*GTC[31]+0.30319611806422586*GBL[31]-0.30319611806422586*GBC[31]-0.21544659739277597*GTL[17]-0.21544659739277597*GTC[17]+0.21544659739277597*GBL[17]+0.21544659739277597*GBC[17]+0.21544659739277597*GTL[16]-0.21544659739277597*GTC[16]-0.43089319478555205*GCL[16]+0.43089319478555205*GCC[16]+0.21544659739277597*GBL[16]-0.21544659739277597*GBC[16]+0.15309310892394856*GTL[8]+0.15309310892394856*GTC[8]-0.3061862178478971*GCL[8]-0.3061862178478971*GCC[8]+0.15309310892394856*GBL[8]+0.15309310892394856*GBC[8]; 
  surft1_lower[11] = -(0.303196118064226*GTL[32])+0.303196118064226*GTC[32]+0.303196118064226*GBL[32]-0.303196118064226*GBC[32]-0.21544659739277597*GTL[21]-0.21544659739277597*GTC[21]+0.21544659739277597*GBL[21]+0.21544659739277597*GBC[21]+0.21544659739277597*GTL[19]-0.21544659739277597*GTC[19]-0.43089319478555205*GCL[19]+0.43089319478555205*GCC[19]+0.21544659739277597*GBL[19]-0.21544659739277597*GBC[19]+0.15309310892394856*GTL[11]+0.15309310892394856*GTC[11]-0.3061862178478971*GCL[11]-0.3061862178478971*GCC[11]+0.15309310892394856*GBL[11]+0.15309310892394856*GBC[11]; 
  surft1_lower[12] = 0.27209908031404895*GTL[34]-0.27209908031404895*GTC[34]-0.27209908031404895*GBL[34]+0.27209908031404895*GBC[34]-0.27209908031404895*GTL[33]-0.27209908031404895*GTC[33]-0.7953665424564508*GCL[33]-0.7953665424564508*GCC[33]-0.27209908031404895*GBL[33]-0.27209908031404895*GBC[33]+0.1933495104806964*GTL[23]+0.1933495104806964*GTC[23]-0.1933495104806964*GBL[23]-0.1933495104806964*GBC[23]+0.1933495104806964*GTL[20]+0.1933495104806964*GTC[20]-0.1933495104806964*GBL[20]-0.1933495104806964*GBC[20]-0.39142450529916156*GTL[15]+0.39142450529916156*GTC[15]-1.14416393856678*GCL[15]+1.14416393856678*GCC[15]-0.39142450529916156*GBL[15]+0.39142450529916156*GBC[15]-0.2781403612330919*GTL[6]-0.2781403612330919*GTC[6]-0.8130256712967302*GCL[6]-0.8130256712967302*GCC[6]-0.2781403612330919*GBL[6]-0.2781403612330919*GBC[6]+0.2781403612330919*GTL[5]-0.2781403612330919*GTC[5]-0.2781403612330919*GBL[5]+0.2781403612330919*GBC[5]+0.19764235376052366*GTL[1]+0.19764235376052366*GTC[1]-0.19764235376052366*GBL[1]-0.19764235376052366*GBC[1]; 
  surft1_lower[13] = -(0.1750503603816304*GTL[44])+0.1750503603816304*GTC[44]+0.3501007207632608*GCL[44]-0.3501007207632608*GCC[44]-0.1750503603816304*GBL[44]+0.1750503603816304*GBC[44]-0.12438815100070813*GTL[37]-0.12438815100070813*GTC[37]+0.24877630200141632*GCL[37]+0.24877630200141632*GCC[37]-0.12438815100070813*GBL[37]-0.12438815100070813*GBC[37]+0.12438815100070813*GTL[35]-0.12438815100070813*GTC[35]-0.12438815100070813*GBL[35]+0.12438815100070813*GBC[35]+0.0883883476483184*GTL[25]+0.0883883476483184*GTC[25]-0.0883883476483184*GBL[25]-0.0883883476483184*GBC[25]; 
  surft1_lower[14] = 0.27209908031404895*GTL[40]-0.27209908031404895*GTC[40]-0.27209908031404895*GBL[40]+0.27209908031404895*GBC[40]-0.27209908031404895*GTL[38]-0.27209908031404895*GTC[38]-0.7953665424564508*GCL[38]-0.7953665424564508*GCC[38]-0.27209908031404895*GBL[38]-0.27209908031404895*GBC[38]+0.1933495104806964*GTL[27]+0.1933495104806964*GTC[27]-0.1933495104806964*GBL[27]-0.1933495104806964*GBC[27]+0.1933495104806964*GTL[26]+0.1933495104806964*GTC[26]-0.1933495104806964*GBL[26]-0.1933495104806964*GBC[26]-0.39142450529916156*GTL[18]+0.39142450529916156*GTC[18]-1.14416393856678*GCL[18]+1.14416393856678*GCC[18]-0.39142450529916156*GBL[18]+0.39142450529916156*GBC[18]-0.2781403612330919*GTL[10]-0.2781403612330919*GTC[10]-0.8130256712967302*GCL[10]-0.8130256712967302*GCC[10]-0.2781403612330919*GBL[10]-0.2781403612330919*GBC[10]+0.2781403612330919*GTL[9]-0.2781403612330919*GTC[9]-0.2781403612330919*GBL[9]+0.2781403612330919*GBC[9]+0.19764235376052366*GTL[4]+0.19764235376052366*GTC[4]-0.19764235376052366*GBL[4]-0.19764235376052366*GBC[4]; 
  surft1_lower[15] = -(0.1750503603816304*GTL[47])+0.1750503603816304*GTC[47]+0.3501007207632608*GCL[47]-0.3501007207632608*GCC[47]-0.1750503603816304*GBL[47]+0.1750503603816304*GBC[47]-0.12438815100070813*GTL[42]-0.12438815100070813*GTC[42]+0.24877630200141632*GCL[42]+0.24877630200141632*GCC[42]-0.12438815100070813*GBL[42]-0.12438815100070813*GBC[42]+0.12438815100070813*GTL[41]-0.12438815100070813*GTC[41]-0.12438815100070813*GBL[41]+0.12438815100070813*GBC[41]+0.0883883476483184*GTL[28]+0.0883883476483184*GTC[28]-0.0883883476483184*GBL[28]-0.0883883476483184*GBC[28]; 
  surft1_lower[16] = -(0.303196118064226*GTL[43])+0.303196118064226*GTC[43]+0.303196118064226*GBL[43]-0.303196118064226*GBC[43]-0.21544659739277597*GTL[30]-0.21544659739277597*GTC[30]+0.21544659739277597*GBL[30]+0.21544659739277597*GBC[30]+0.21544659739277597*GTL[29]-0.21544659739277597*GTC[29]-0.43089319478555205*GCL[29]+0.43089319478555205*GCC[29]+0.21544659739277597*GBL[29]-0.21544659739277597*GBC[29]+0.15309310892394856*GTL[14]+0.15309310892394856*GTC[14]-0.3061862178478971*GCL[14]-0.3061862178478971*GCC[14]+0.15309310892394856*GBL[14]+0.15309310892394856*GBC[14]; 
  surft1_lower[17] = -(0.303196118064226*GTL[44])+0.303196118064226*GTC[44]+0.303196118064226*GBL[44]-0.303196118064226*GBC[44]-0.21544659739277597*GTL[37]-0.21544659739277597*GTC[37]+0.21544659739277597*GBL[37]+0.21544659739277597*GBC[37]+0.21544659739277597*GTL[35]-0.21544659739277597*GTC[35]-0.43089319478555205*GCL[35]+0.43089319478555205*GCC[35]+0.21544659739277597*GBL[35]-0.21544659739277597*GBC[35]+0.15309310892394856*GTL[25]+0.15309310892394856*GTC[25]-0.3061862178478971*GCL[25]-0.3061862178478971*GCC[25]+0.15309310892394856*GBL[25]+0.15309310892394856*GBC[25]; 
  surft1_lower[18] = 0.27209908031404895*GTL[46]-0.27209908031404895*GTC[46]-0.27209908031404895*GBL[46]+0.27209908031404895*GBC[46]-0.27209908031404895*GTL[45]-0.27209908031404895*GTC[45]-0.7953665424564508*GCL[45]-0.7953665424564508*GCC[45]-0.27209908031404895*GBL[45]-0.27209908031404895*GBC[45]+0.1933495104806964*GTL[39]+0.1933495104806964*GTC[39]-0.1933495104806964*GBL[39]-0.1933495104806964*GBC[39]+0.1933495104806964*GTL[36]+0.1933495104806964*GTC[36]-0.1933495104806964*GBL[36]-0.1933495104806964*GBC[36]-0.3914245052991616*GTL[31]+0.3914245052991616*GTC[31]-1.1441639385667801*GCL[31]+1.1441639385667801*GCC[31]-0.3914245052991616*GBL[31]+0.3914245052991616*GBC[31]-0.2781403612330919*GTL[17]-0.2781403612330919*GTC[17]-0.8130256712967302*GCL[17]-0.8130256712967302*GCC[17]-0.2781403612330919*GBL[17]-0.2781403612330919*GBC[17]+0.2781403612330919*GTL[16]-0.2781403612330919*GTC[16]-0.2781403612330919*GBL[16]+0.2781403612330919*GBC[16]+0.19764235376052364*GTL[8]+0.19764235376052364*GTC[8]-0.19764235376052364*GBL[8]-0.19764235376052364*GBC[8]; 
  surft1_lower[19] = -(0.303196118064226*GTL[47])+0.303196118064226*GTC[47]+0.303196118064226*GBL[47]-0.303196118064226*GBC[47]-0.21544659739277597*GTL[42]-0.21544659739277597*GTC[42]+0.21544659739277597*GBL[42]+0.21544659739277597*GBC[42]+0.21544659739277597*GTL[41]-0.21544659739277597*GTC[41]-0.43089319478555205*GCL[41]+0.43089319478555205*GCC[41]+0.21544659739277597*GBL[41]-0.21544659739277597*GBC[41]+0.15309310892394856*GTL[28]+0.15309310892394856*GTC[28]-0.3061862178478971*GCL[28]-0.3061862178478971*GCC[28]+0.15309310892394856*GBL[28]+0.15309310892394856*GBC[28]; 

  surft2_upper[0] = 0.34587411908091625*GTC[13]+0.34587411908091625*GCC[13]-0.49755260400283263*GTC[3]+0.49755260400283263*GCC[3]+0.3535533905932737*GTC[0]+0.3535533905932737*GCC[0]; 
  surft2_upper[1] = 0.34587411908091625*GTC[23]+0.34587411908091625*GCC[23]-0.49755260400283263*GTC[6]+0.49755260400283263*GCC[6]+0.3535533905932737*GTC[1]+0.3535533905932737*GCC[1]; 
  surft2_upper[2] = 0.34587411908091625*GTC[24]+0.34587411908091625*GCC[24]-0.49755260400283263*GTC[7]+0.49755260400283263*GCC[7]+0.3535533905932737*GTC[2]+0.3535533905932737*GCC[2]; 
  surft2_upper[3] = 0.34587411908091625*GTC[27]+0.34587411908091625*GCC[27]-0.49755260400283263*GTC[10]+0.49755260400283263*GCC[10]+0.3535533905932737*GTC[4]+0.3535533905932737*GCC[4]; 
  surft2_upper[4] = 0.34587411908091625*GTC[34]+0.34587411908091625*GCC[34]-0.49755260400283263*GTC[15]+0.49755260400283263*GCC[15]+0.3535533905932737*GTC[5]+0.3535533905932737*GCC[5]; 
  surft2_upper[5] = 0.34587411908091625*GTC[39]+0.34587411908091625*GCC[39]-0.49755260400283263*GTC[17]+0.49755260400283263*GCC[17]+0.3535533905932737*GTC[8]+0.3535533905932737*GCC[8]; 
  surft2_upper[6] = 0.34587411908091625*GTC[40]+0.34587411908091625*GCC[40]-0.49755260400283263*GTC[18]+0.49755260400283263*GCC[18]+0.3535533905932737*GTC[9]+0.3535533905932737*GCC[9]; 
  surft2_upper[7] = -(0.49755260400283263*GTC[21])+0.49755260400283263*GCC[21]+0.3535533905932737*GTC[11]+0.3535533905932737*GCC[11]; 
  surft2_upper[8] = -(0.49755260400283263*GTC[22])+0.49755260400283263*GCC[22]+0.3535533905932737*GTC[12]+0.3535533905932737*GCC[12]; 
  surft2_upper[9] = -(0.49755260400283263*GTC[30])+0.49755260400283263*GCC[30]+0.3535533905932737*GTC[14]+0.3535533905932737*GCC[14]; 
  surft2_upper[10] = 0.34587411908091625*GTC[46]+0.34587411908091625*GCC[46]-0.49755260400283263*GTC[31]+0.49755260400283263*GCC[31]+0.3535533905932737*GTC[16]+0.3535533905932737*GCC[16]; 
  surft2_upper[11] = -(0.49755260400283263*GTC[32])+0.49755260400283263*GCC[32]+0.3535533905932737*GTC[19]+0.3535533905932737*GCC[19]; 
  surft2_upper[12] = -(0.49755260400283263*GTC[33])+0.49755260400283263*GCC[33]+0.3535533905932737*GTC[20]+0.3535533905932737*GCC[20]; 
  surft2_upper[13] = -(0.49755260400283263*GTC[37])+0.49755260400283263*GCC[37]+0.3535533905932737*GTC[25]+0.3535533905932737*GCC[25]; 
  surft2_upper[14] = -(0.49755260400283263*GTC[38])+0.49755260400283263*GCC[38]+0.3535533905932737*GTC[26]+0.3535533905932737*GCC[26]; 
  surft2_upper[15] = -(0.49755260400283263*GTC[42])+0.49755260400283263*GCC[42]+0.3535533905932737*GTC[28]+0.3535533905932737*GCC[28]; 
  surft2_upper[16] = -(0.49755260400283263*GTC[43])+0.49755260400283263*GCC[43]+0.3535533905932737*GTC[29]+0.3535533905932737*GCC[29]; 
  surft2_upper[17] = -(0.49755260400283263*GTC[44])+0.49755260400283263*GCC[44]+0.3535533905932737*GTC[35]+0.3535533905932737*GCC[35]; 
  surft2_upper[18] = -(0.49755260400283263*GTC[45])+0.49755260400283263*GCC[45]+0.3535533905932737*GTC[36]+0.3535533905932737*GCC[36]; 
  surft2_upper[19] = -(0.49755260400283263*GTC[47])+0.49755260400283263*GCC[47]+0.3535533905932737*GTC[41]+0.3535533905932737*GCC[41]; 
  surft2_lower[0] = 0.34587411908091625*GCC[13]+0.34587411908091625*GBC[13]-0.49755260400283263*GCC[3]+0.49755260400283263*GBC[3]+0.3535533905932737*GCC[0]+0.3535533905932737*GBC[0]; 
  surft2_lower[1] = 0.34587411908091625*GCC[23]+0.34587411908091625*GBC[23]-0.49755260400283263*GCC[6]+0.49755260400283263*GBC[6]+0.3535533905932737*GCC[1]+0.3535533905932737*GBC[1]; 
  surft2_lower[2] = 0.34587411908091625*GCC[24]+0.34587411908091625*GBC[24]-0.49755260400283263*GCC[7]+0.49755260400283263*GBC[7]+0.3535533905932737*GCC[2]+0.3535533905932737*GBC[2]; 
  surft2_lower[3] = 0.34587411908091625*GCC[27]+0.34587411908091625*GBC[27]-0.49755260400283263*GCC[10]+0.49755260400283263*GBC[10]+0.3535533905932737*GCC[4]+0.3535533905932737*GBC[4]; 
  surft2_lower[4] = 0.34587411908091625*GCC[34]+0.34587411908091625*GBC[34]-0.49755260400283263*GCC[15]+0.49755260400283263*GBC[15]+0.3535533905932737*GCC[5]+0.3535533905932737*GBC[5]; 
  surft2_lower[5] = 0.34587411908091625*GCC[39]+0.34587411908091625*GBC[39]-0.49755260400283263*GCC[17]+0.49755260400283263*GBC[17]+0.3535533905932737*GCC[8]+0.3535533905932737*GBC[8]; 
  surft2_lower[6] = 0.34587411908091625*GCC[40]+0.34587411908091625*GBC[40]-0.49755260400283263*GCC[18]+0.49755260400283263*GBC[18]+0.3535533905932737*GCC[9]+0.3535533905932737*GBC[9]; 
  surft2_lower[7] = -(0.49755260400283263*GCC[21])+0.49755260400283263*GBC[21]+0.3535533905932737*GCC[11]+0.3535533905932737*GBC[11]; 
  surft2_lower[8] = -(0.49755260400283263*GCC[22])+0.49755260400283263*GBC[22]+0.3535533905932737*GCC[12]+0.3535533905932737*GBC[12]; 
  surft2_lower[9] = -(0.49755260400283263*GCC[30])+0.49755260400283263*GBC[30]+0.3535533905932737*GCC[14]+0.3535533905932737*GBC[14]; 
  surft2_lower[10] = 0.34587411908091625*GCC[46]+0.34587411908091625*GBC[46]-0.49755260400283263*GCC[31]+0.49755260400283263*GBC[31]+0.3535533905932737*GCC[16]+0.3535533905932737*GBC[16]; 
  surft2_lower[11] = -(0.49755260400283263*GCC[32])+0.49755260400283263*GBC[32]+0.3535533905932737*GCC[19]+0.3535533905932737*GBC[19]; 
  surft2_lower[12] = -(0.49755260400283263*GCC[33])+0.49755260400283263*GBC[33]+0.3535533905932737*GCC[20]+0.3535533905932737*GBC[20]; 
  surft2_lower[13] = -(0.49755260400283263*GCC[37])+0.49755260400283263*GBC[37]+0.3535533905932737*GCC[25]+0.3535533905932737*GBC[25]; 
  surft2_lower[14] = -(0.49755260400283263*GCC[38])+0.49755260400283263*GBC[38]+0.3535533905932737*GCC[26]+0.3535533905932737*GBC[26]; 
  surft2_lower[15] = -(0.49755260400283263*GCC[42])+0.49755260400283263*GBC[42]+0.3535533905932737*GCC[28]+0.3535533905932737*GBC[28]; 
  surft2_lower[16] = -(0.49755260400283263*GCC[43])+0.49755260400283263*GBC[43]+0.3535533905932737*GCC[29]+0.3535533905932737*GBC[29]; 
  surft2_lower[17] = -(0.49755260400283263*GCC[44])+0.49755260400283263*GBC[44]+0.3535533905932737*GCC[35]+0.3535533905932737*GBC[35]; 
  surft2_lower[18] = -(0.49755260400283263*GCC[45])+0.49755260400283263*GBC[45]+0.3535533905932737*GCC[36]+0.3535533905932737*GBC[36]; 
  surft2_lower[19] = -(0.49755260400283263*GCC[47])+0.49755260400283263*GBC[47]+0.3535533905932737*GCC[41]+0.3535533905932737*GBC[41]; 

  out[0] = 0.7071067811865475*surft1_upper[0]*dv1_sq*gamma-0.7071067811865475*surft1_lower[0]*dv1_sq*gamma; 
  out[1] = 0.7071067811865475*surft1_upper[1]*dv1_sq*gamma-0.7071067811865475*surft1_lower[1]*dv1_sq*gamma; 
  out[2] = -(1.224744871391589*surft2_upper[0]*dv1_sq*gamma)+1.224744871391589*surft2_lower[0]*dv1_sq*gamma+1.224744871391589*surft1_upper[0]*dv1_sq*gamma+1.224744871391589*surft1_lower[0]*dv1_sq*gamma; 
  out[3] = 0.7071067811865475*surft1_upper[2]*dv1_sq*gamma-0.7071067811865475*surft1_lower[2]*dv1_sq*gamma; 
  out[4] = 0.7071067811865475*surft1_upper[3]*dv1_sq*gamma-0.7071067811865475*surft1_lower[3]*dv1_sq*gamma; 
  out[5] = -(1.224744871391589*surft2_upper[1]*dv1_sq*gamma)+1.224744871391589*surft2_lower[1]*dv1_sq*gamma+1.224744871391589*surft1_upper[1]*dv1_sq*gamma+1.224744871391589*surft1_lower[1]*dv1_sq*gamma; 
  out[6] = 0.7071067811865475*surft1_upper[4]*dv1_sq*gamma-0.7071067811865475*surft1_lower[4]*dv1_sq*gamma; 
  out[7] = 1.224744871391589*surft1_upper[2]*dv1_sq*gamma+1.224744871391589*surft1_lower[2]*dv1_sq*gamma-2.1213203435596424*surft2_upper[0]*dv1_sq*gamma-2.1213203435596424*surft2_lower[0]*dv1_sq*gamma+3.0*GCC[0]*dv1_sq*gamma; 
  out[8] = 0.7071067811865475*surft1_upper[5]*dv1_sq*gamma-0.7071067811865475*surft1_lower[5]*dv1_sq*gamma; 
  out[9] = -(1.224744871391589*surft2_upper[3]*dv1_sq*gamma)+1.224744871391589*surft2_lower[3]*dv1_sq*gamma+1.224744871391589*surft1_upper[3]*dv1_sq*gamma+1.224744871391589*surft1_lower[3]*dv1_sq*gamma; 
  out[10] = 0.7071067811865475*surft1_upper[6]*dv1_sq*gamma-0.7071067811865475*surft1_lower[6]*dv1_sq*gamma; 
  out[11] = 0.7071067811865475*surft1_upper[7]*dv1_sq*gamma-0.7071067811865475*surft1_lower[7]*dv1_sq*gamma; 
  out[12] = -(2.7386127875258306*surft2_upper[2]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[2]*dv1_sq*gamma+1.5811388300841895*surft1_upper[0]*dv1_sq*gamma-1.5811388300841895*surft1_lower[0]*dv1_sq*gamma; 
  out[13] = 0.7071067811865475*surft1_upper[8]*dv1_sq*gamma-0.7071067811865475*surft1_lower[8]*dv1_sq*gamma; 
  out[14] = 0.7071067811865475*surft1_upper[9]*dv1_sq*gamma-0.7071067811865475*surft1_lower[9]*dv1_sq*gamma; 
  out[15] = 1.224744871391589*surft1_upper[4]*dv1_sq*gamma+1.224744871391589*surft1_lower[4]*dv1_sq*gamma-2.1213203435596424*surft2_upper[1]*dv1_sq*gamma-2.1213203435596424*surft2_lower[1]*dv1_sq*gamma+3.0*GCC[1]*dv1_sq*gamma; 
  out[16] = -(1.224744871391589*surft2_upper[5]*dv1_sq*gamma)+1.224744871391589*surft2_lower[5]*dv1_sq*gamma+1.224744871391589*surft1_upper[5]*dv1_sq*gamma+1.224744871391589*surft1_lower[5]*dv1_sq*gamma; 
  out[17] = 0.7071067811865475*surft1_upper[10]*dv1_sq*gamma-0.7071067811865475*surft1_lower[10]*dv1_sq*gamma; 
  out[18] = 1.224744871391589*surft1_upper[6]*dv1_sq*gamma+1.224744871391589*surft1_lower[6]*dv1_sq*gamma+3.0*GCC[4]*dv1_sq*gamma-2.1213203435596424*surft2_upper[3]*dv1_sq*gamma-2.1213203435596424*surft2_lower[3]*dv1_sq*gamma; 
  out[19] = -(1.224744871391589*surft2_upper[7]*dv1_sq*gamma)+1.224744871391589*surft2_lower[7]*dv1_sq*gamma+1.224744871391589*surft1_upper[7]*dv1_sq*gamma+1.224744871391589*surft1_lower[7]*dv1_sq*gamma; 
  out[20] = -(2.7386127875258306*surft2_upper[4]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[4]*dv1_sq*gamma+1.5811388300841898*surft1_upper[1]*dv1_sq*gamma-1.5811388300841898*surft1_lower[1]*dv1_sq*gamma; 
  out[21] = 0.7071067811865475*surft1_upper[11]*dv1_sq*gamma-0.7071067811865475*surft1_lower[11]*dv1_sq*gamma; 
  out[22] = -(4.743416490252569*surft2_upper[2]*dv1_sq*gamma)-4.743416490252569*surft2_lower[2]*dv1_sq*gamma+1.5811388300841898*surft1_upper[2]*dv1_sq*gamma-1.5811388300841898*surft1_lower[2]*dv1_sq*gamma+6.7082039324993685*GCC[2]*dv1_sq*gamma; 
  out[23] = 0.7071067811865475*surft1_upper[12]*dv1_sq*gamma-0.7071067811865475*surft1_lower[12]*dv1_sq*gamma; 
  out[24] = 1.224744871391589*surft1_upper[8]*dv1_sq*gamma+1.224744871391589*surft1_lower[8]*dv1_sq*gamma+6.7082039324993685*GCC[3]*dv1_sq*gamma-2.7386127875258306*surft2_upper[0]*dv1_sq*gamma+2.7386127875258306*surft2_lower[0]*dv1_sq*gamma; 
  out[25] = 0.7071067811865475*surft1_upper[13]*dv1_sq*gamma-0.7071067811865475*surft1_lower[13]*dv1_sq*gamma; 
  out[26] = -(2.7386127875258306*surft2_upper[6]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[6]*dv1_sq*gamma+1.5811388300841898*surft1_upper[3]*dv1_sq*gamma-1.5811388300841898*surft1_lower[3]*dv1_sq*gamma; 
  out[27] = 0.7071067811865475*surft1_upper[14]*dv1_sq*gamma-0.7071067811865475*surft1_lower[14]*dv1_sq*gamma; 
  out[28] = 0.7071067811865475*surft1_upper[15]*dv1_sq*gamma-0.7071067811865475*surft1_lower[15]*dv1_sq*gamma; 
  out[29] = -(1.224744871391589*surft2_upper[9]*dv1_sq*gamma)+1.224744871391589*surft2_lower[9]*dv1_sq*gamma+1.224744871391589*surft1_upper[9]*dv1_sq*gamma+1.224744871391589*surft1_lower[9]*dv1_sq*gamma; 
  out[30] = 0.7071067811865475*surft1_upper[16]*dv1_sq*gamma-0.7071067811865475*surft1_lower[16]*dv1_sq*gamma; 
  out[31] = 1.224744871391589*surft1_upper[10]*dv1_sq*gamma+1.224744871391589*surft1_lower[10]*dv1_sq*gamma+3.0*GCC[8]*dv1_sq*gamma-2.1213203435596424*surft2_upper[5]*dv1_sq*gamma-2.1213203435596424*surft2_lower[5]*dv1_sq*gamma; 
  out[32] = 1.224744871391589*surft1_upper[11]*dv1_sq*gamma+1.224744871391589*surft1_lower[11]*dv1_sq*gamma+3.0*GCC[11]*dv1_sq*gamma-2.1213203435596424*surft2_upper[7]*dv1_sq*gamma-2.1213203435596424*surft2_lower[7]*dv1_sq*gamma; 
  out[33] = 6.708203932499369*GCC[5]*dv1_sq*gamma-4.743416490252569*surft2_upper[4]*dv1_sq*gamma-4.743416490252569*surft2_lower[4]*dv1_sq*gamma+1.5811388300841895*surft1_upper[4]*dv1_sq*gamma-1.5811388300841895*surft1_lower[4]*dv1_sq*gamma; 
  out[34] = 1.224744871391589*surft1_upper[12]*dv1_sq*gamma+1.224744871391589*surft1_lower[12]*dv1_sq*gamma+6.708203932499369*GCC[6]*dv1_sq*gamma-2.7386127875258306*surft2_upper[1]*dv1_sq*gamma+2.7386127875258306*surft2_lower[1]*dv1_sq*gamma; 
  out[35] = -(1.224744871391589*surft2_upper[13]*dv1_sq*gamma)+1.224744871391589*surft2_lower[13]*dv1_sq*gamma+1.224744871391589*surft1_upper[13]*dv1_sq*gamma+1.224744871391589*surft1_lower[13]*dv1_sq*gamma; 
  out[36] = -(2.7386127875258306*surft2_upper[10]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[10]*dv1_sq*gamma+1.5811388300841895*surft1_upper[5]*dv1_sq*gamma-1.5811388300841895*surft1_lower[5]*dv1_sq*gamma; 
  out[37] = 0.7071067811865475*surft1_upper[17]*dv1_sq*gamma-0.7071067811865475*surft1_lower[17]*dv1_sq*gamma; 
  out[38] = 6.708203932499369*GCC[9]*dv1_sq*gamma-4.743416490252569*surft2_upper[6]*dv1_sq*gamma-4.743416490252569*surft2_lower[6]*dv1_sq*gamma+1.5811388300841895*surft1_upper[6]*dv1_sq*gamma-1.5811388300841895*surft1_lower[6]*dv1_sq*gamma; 
  out[39] = 0.7071067811865475*surft1_upper[18]*dv1_sq*gamma-0.7071067811865475*surft1_lower[18]*dv1_sq*gamma; 
  out[40] = 1.224744871391589*surft1_upper[14]*dv1_sq*gamma+1.224744871391589*surft1_lower[14]*dv1_sq*gamma+6.708203932499369*GCC[10]*dv1_sq*gamma-2.7386127875258306*surft2_upper[3]*dv1_sq*gamma+2.7386127875258306*surft2_lower[3]*dv1_sq*gamma; 
  out[41] = -(1.224744871391589*surft2_upper[15]*dv1_sq*gamma)+1.224744871391589*surft2_lower[15]*dv1_sq*gamma+1.224744871391589*surft1_upper[15]*dv1_sq*gamma+1.224744871391589*surft1_lower[15]*dv1_sq*gamma; 
  out[42] = 0.7071067811865475*surft1_upper[19]*dv1_sq*gamma-0.7071067811865475*surft1_lower[19]*dv1_sq*gamma; 
  out[43] = 1.224744871391589*surft1_upper[16]*dv1_sq*gamma+1.224744871391589*surft1_lower[16]*dv1_sq*gamma+3.0*GCC[14]*dv1_sq*gamma-2.1213203435596424*surft2_upper[9]*dv1_sq*gamma-2.1213203435596424*surft2_lower[9]*dv1_sq*gamma; 
  out[44] = 3.0*GCC[25]*dv1_sq*gamma+1.224744871391589*surft1_upper[17]*dv1_sq*gamma+1.224744871391589*surft1_lower[17]*dv1_sq*gamma-2.1213203435596424*surft2_upper[13]*dv1_sq*gamma-2.1213203435596424*surft2_lower[13]*dv1_sq*gamma; 
  out[45] = 6.7082039324993685*GCC[16]*dv1_sq*gamma-4.743416490252569*surft2_upper[10]*dv1_sq*gamma-4.743416490252569*surft2_lower[10]*dv1_sq*gamma+1.5811388300841898*surft1_upper[10]*dv1_sq*gamma-1.5811388300841898*surft1_lower[10]*dv1_sq*gamma; 
  out[46] = 1.224744871391589*surft1_upper[18]*dv1_sq*gamma+1.224744871391589*surft1_lower[18]*dv1_sq*gamma+6.7082039324993685*GCC[17]*dv1_sq*gamma-2.7386127875258306*surft2_upper[5]*dv1_sq*gamma+2.7386127875258306*surft2_lower[5]*dv1_sq*gamma; 
  out[47] = 3.0*GCC[28]*dv1_sq*gamma+1.224744871391589*surft1_upper[19]*dv1_sq*gamma+1.224744871391589*surft1_lower[19]*dv1_sq*gamma-2.1213203435596424*surft2_upper[15]*dv1_sq*gamma-2.1213203435596424*surft2_lower[15]*dv1_sq*gamma; 
} 
