#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p2_upvy_upvx(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[2]; 
  double dv2 = 2.0/dxv[1]; 
  double dv1_sq = 4.0/dxv[2]/dxv[1]; 
 
  double surft1_lo[20] = {0.0}; 
  double surft1_up[20] = {0.0}; 
  double surft2_lo[20] = {0.0}; 
  double surft2_up[20] = {0.0}; 
  double vol[48] = {0.0}; 
  double *out = &diff_coeff[144]; 

  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GBC = fpo_g_stencil[2]; 
  const double* GCC = fpo_g_stencil[3]; 
  const double* G_surf_CC_vx = &fpo_g_surf_stencil[3][0]; 
  const double* dGdvx_surf_CC_vy = &fpo_dgdv_surf[40]; 

  surft1_lo[0] = -(0.12168640803947765*GCL[24])+1.490992801802391*GCC[24]-0.12168640803947765*GBL[24]+0.12168640803947765*GBC[24]+0.12168640803947765*GCL[22]-1.2476199857234356*GCC[22]-0.12168640803947765*GBL[22]-0.12168640803947765*GBC[22]-0.08646852977022904*GCL[13]+0.704100885271865*GCC[13]-0.08646852977022904*GBL[13]-0.08646852977022904*GBC[13]-0.08646852977022904*GCL[12]+0.704100885271865*GCC[12]-0.08646852977022904*GBL[12]-0.08646852977022904*GBC[12]+0.1750503603816304*GCL[7]-1.23571053216145*GCC[7]-0.1750503603816304*GBL[7]+0.1750503603816304*GBC[7]+0.12438815100070813*GCL[3]-0.48798428469508565*GCC[3]-0.12438815100070813*GBL[3]-0.12438815100070813*GBC[3]-0.12438815100070813*GCL[2]+0.736760586696502*GCC[2]-0.12438815100070813*GBL[2]+0.12438815100070813*GBC[2]-0.0883883476483184*GCL[0]+0.26516504294495524*GCC[0]-0.0883883476483184*GBL[0]-0.0883883476483184*GBC[0]; 
  surft1_lo[1] = -(0.12168640803947765*GCL[34])+1.490992801802391*GCC[34]-0.12168640803947765*GBL[34]+0.12168640803947765*GBC[34]+0.12168640803947765*GCL[33]-1.2476199857234354*GCC[33]-0.12168640803947765*GBL[33]-0.12168640803947765*GBC[33]-0.08646852977022904*GCL[23]+0.704100885271865*GCC[23]-0.08646852977022904*GBL[23]-0.08646852977022904*GBC[23]-0.08646852977022904*GCL[20]+0.704100885271865*GCC[20]-0.08646852977022904*GBL[20]-0.08646852977022904*GBC[20]+0.1750503603816304*GCL[15]-1.23571053216145*GCC[15]-0.1750503603816304*GBL[15]+0.1750503603816304*GBC[15]+0.12438815100070813*GCL[6]-0.48798428469508565*GCC[6]-0.12438815100070813*GBL[6]-0.12438815100070813*GBC[6]-0.12438815100070813*GCL[5]+0.736760586696502*GCC[5]-0.12438815100070813*GBL[5]+0.12438815100070813*GBC[5]-0.0883883476483184*GCL[1]+0.26516504294495524*GCC[1]-0.0883883476483184*GBL[1]-0.0883883476483184*GBC[1]; 
  surft1_lo[2] = 0.2107670413149332*GCL[24]+2.1609412038113476*GCC[24]+0.2107670413149332*GBL[24]-0.2107670413149332*GBC[24]-0.2107670413149332*GCL[22]-2.582475286441214*GCC[22]+0.2107670413149332*GBL[22]+0.2107670413149332*GBC[22]+0.1497678868178187*GCL[13]+0.9200027333094577*GCC[13]+0.1497678868178187*GBL[13]-0.44930366045345604*GBC[13]+0.1497678868178187*GCL[12]+1.5190742805807327*GCC[12]+0.1497678868178187*GBL[12]+0.1497678868178187*GBC[12]-0.30319611806422586*GCL[7]-1.5339211890231543*GCC[7]+0.30319611806422586*GBL[7]-0.30319611806422586*GBC[7]-0.21544659739277597*GCL[3]-0.41432037960149226*GCC[3]+0.21544659739277597*GBL[3]-0.6463397921783279*GBC[3]+0.21544659739277597*GCL[2]+0.8452135743870443*GCC[2]+0.21544659739277597*GBL[2]-0.21544659739277597*GBC[2]+0.15309310892394856*GCL[0]+0.15309310892394856*GCC[0]+0.15309310892394856*GBL[0]-0.45927932677184563*GBC[0]; 
  surft1_lo[3] = -(0.12168640803947765*GCL[40])+1.490992801802391*GCC[40]-0.12168640803947765*GBL[40]+0.12168640803947765*GBC[40]+0.12168640803947765*GCL[38]-1.2476199857234354*GCC[38]-0.12168640803947765*GBL[38]-0.12168640803947765*GBC[38]-0.08646852977022904*GCL[27]+0.704100885271865*GCC[27]-0.08646852977022904*GBL[27]-0.08646852977022904*GBC[27]-0.08646852977022904*GCL[26]+0.704100885271865*GCC[26]-0.08646852977022904*GBL[26]-0.08646852977022904*GBC[26]+0.1750503603816304*GCL[18]-1.23571053216145*GCC[18]-0.1750503603816304*GBL[18]+0.1750503603816304*GBC[18]+0.12438815100070813*GCL[10]-0.48798428469508565*GCC[10]-0.12438815100070813*GBL[10]-0.12438815100070813*GBC[10]-0.12438815100070813*GCL[9]+0.736760586696502*GCC[9]-0.12438815100070813*GBL[9]+0.12438815100070813*GBC[9]-0.0883883476483184*GCL[4]+0.26516504294495524*GCC[4]-0.0883883476483184*GBL[4]-0.0883883476483184*GBC[4]; 
  surft1_lo[4] = 0.21076704131493318*GCL[34]+2.160941203811348*GCC[34]+0.21076704131493318*GBL[34]-0.21076704131493318*GBC[34]-0.21076704131493318*GCL[33]-2.5824752864412144*GCC[33]+0.21076704131493318*GBL[33]+0.21076704131493318*GBC[33]+0.1497678868178187*GCL[23]+0.9200027333094578*GCC[23]+0.1497678868178187*GBL[23]-0.44930366045345616*GBC[23]+0.1497678868178187*GCL[20]+1.5190742805807327*GCC[20]+0.1497678868178187*GBL[20]+0.1497678868178187*GBC[20]-0.30319611806422586*GCL[15]-1.5339211890231543*GCC[15]+0.30319611806422586*GBL[15]-0.30319611806422586*GBC[15]-0.21544659739277597*GCL[6]-0.41432037960149226*GCC[6]+0.21544659739277597*GBL[6]-0.6463397921783279*GBC[6]+0.21544659739277597*GCL[5]+0.8452135743870443*GCC[5]+0.21544659739277597*GBL[5]-0.21544659739277597*GBC[5]+0.15309310892394856*GCL[1]+0.15309310892394856*GCC[1]+0.15309310892394856*GBL[1]-0.45927932677184563*GBC[1]; 
  surft1_lo[5] = -(0.12168640803947765*GCL[46])+1.490992801802391*GCC[46]-0.12168640803947765*GBL[46]+0.12168640803947765*GBC[46]+0.12168640803947765*GCL[45]-1.2476199857234356*GCC[45]-0.12168640803947765*GBL[45]-0.12168640803947765*GBC[45]-0.08646852977022904*GCL[39]+0.704100885271865*GCC[39]-0.08646852977022904*GBL[39]-0.08646852977022904*GBC[39]-0.08646852977022904*GCL[36]+0.704100885271865*GCC[36]-0.08646852977022904*GBL[36]-0.08646852977022904*GBC[36]+0.1750503603816304*GCL[31]-1.23571053216145*GCC[31]-0.1750503603816304*GBL[31]+0.1750503603816304*GBC[31]+0.12438815100070813*GCL[17]-0.48798428469508565*GCC[17]-0.12438815100070813*GBL[17]-0.12438815100070813*GBC[17]-0.12438815100070813*GCL[16]+0.736760586696502*GCC[16]-0.12438815100070813*GBL[16]+0.12438815100070813*GBC[16]-0.0883883476483184*GCL[8]+0.26516504294495524*GCC[8]-0.0883883476483184*GBL[8]-0.0883883476483184*GBC[8]; 
  surft1_lo[6] = 0.21076704131493318*GCL[40]+2.160941203811348*GCC[40]+0.21076704131493318*GBL[40]-0.21076704131493318*GBC[40]-0.21076704131493318*GCL[38]-2.5824752864412144*GCC[38]+0.21076704131493318*GBL[38]+0.21076704131493318*GBC[38]+0.1497678868178187*GCL[27]+0.9200027333094578*GCC[27]+0.1497678868178187*GBL[27]-0.44930366045345616*GBC[27]+0.1497678868178187*GCL[26]+1.5190742805807327*GCC[26]+0.1497678868178187*GBL[26]+0.1497678868178187*GBC[26]-0.30319611806422586*GCL[18]-1.5339211890231543*GCC[18]+0.30319611806422586*GBL[18]-0.30319611806422586*GBC[18]-0.21544659739277597*GCL[10]-0.41432037960149226*GCC[10]+0.21544659739277597*GBL[10]-0.6463397921783279*GBC[10]+0.21544659739277597*GCL[9]+0.8452135743870443*GCC[9]+0.21544659739277597*GBL[9]-0.21544659739277597*GBC[9]+0.15309310892394856*GCL[4]+0.15309310892394856*GCC[4]+0.15309310892394856*GBL[4]-0.45927932677184563*GBC[4]; 
  surft1_lo[7] = 0.1750503603816304*GCL[32]-1.23571053216145*GCC[32]-0.1750503603816304*GBL[32]+0.1750503603816304*GBC[32]+0.12438815100070813*GCL[21]-0.48798428469508576*GCC[21]-0.12438815100070813*GBL[21]-0.12438815100070813*GBC[21]-0.12438815100070813*GCL[19]+0.736760586696502*GCC[19]-0.12438815100070813*GBL[19]+0.12438815100070813*GBC[19]-0.0883883476483184*GCL[11]+0.26516504294495524*GCC[11]-0.0883883476483184*GBL[11]-0.0883883476483184*GBC[11]; 
  surft1_lo[8] = -(0.27209908031404895*GCL[24])+1.9943965557084689*GCC[24]-0.27209908031404895*GBL[24]-1.0674656227704997*GBC[24]+0.27209908031404895*GCL[22]-2.789763098164919*GCC[22]-0.27209908031404895*GBL[22]-0.27209908031404895*GBC[22]-0.1933495104806964*GCL[13]+1.5744174424856705*GCC[13]-0.1933495104806964*GBL[13]-0.1933495104806964*GBC[13]-0.1933495104806964*GCL[12]+1.5744174424856705*GCC[12]-0.1933495104806964*GBL[12]-0.1933495104806964*GBC[12]+0.3914245052991616*GCL[7]-0.8361198012603392*GCC[7]-0.3914245052991616*GBL[7]-1.5355884438659417*GBC[7]+0.2781403612330919*GCL[3]-1.0911660325298218*GCC[3]-0.2781403612330919*GBL[3]-0.2781403612330919*GBC[3]-0.2781403612330919*GCL[2]+0.2781403612330919*GCC[2]-0.2781403612330919*GBL[2]-1.0911660325298218*GBC[2]-0.19764235376052364*GCL[0]+0.592927061281571*GCC[0]-0.19764235376052364*GBL[0]-0.19764235376052364*GBC[0]; 
  surft1_lo[9] = 0.1750503603816304*GCL[43]-1.23571053216145*GCC[43]-0.1750503603816304*GBL[43]+0.1750503603816304*GBC[43]+0.12438815100070813*GCL[30]-0.48798428469508576*GCC[30]-0.12438815100070813*GBL[30]-0.12438815100070813*GBC[30]-0.12438815100070813*GCL[29]+0.736760586696502*GCC[29]-0.12438815100070813*GBL[29]+0.12438815100070813*GBC[29]-0.0883883476483184*GCL[14]+0.26516504294495524*GCC[14]-0.0883883476483184*GBL[14]-0.0883883476483184*GBC[14]; 
  surft1_lo[10] = 0.2107670413149332*GCL[46]+2.1609412038113476*GCC[46]+0.2107670413149332*GBL[46]-0.2107670413149332*GBC[46]-0.2107670413149332*GCL[45]-2.582475286441214*GCC[45]+0.2107670413149332*GBL[45]+0.2107670413149332*GBC[45]+0.1497678868178187*GCL[39]+0.9200027333094577*GCC[39]+0.1497678868178187*GBL[39]-0.44930366045345604*GBC[39]+0.1497678868178187*GCL[36]+1.5190742805807327*GCC[36]+0.1497678868178187*GBL[36]+0.1497678868178187*GBC[36]-0.30319611806422586*GCL[31]-1.5339211890231543*GCC[31]+0.30319611806422586*GBL[31]-0.30319611806422586*GBC[31]-0.21544659739277597*GCL[17]-0.41432037960149226*GCC[17]+0.21544659739277597*GBL[17]-0.6463397921783279*GBC[17]+0.21544659739277597*GCL[16]+0.8452135743870443*GCC[16]+0.21544659739277597*GBL[16]-0.21544659739277597*GBC[16]+0.15309310892394856*GCL[8]+0.15309310892394856*GCC[8]+0.15309310892394856*GBL[8]-0.45927932677184563*GBC[8]; 
  surft1_lo[11] = -(0.303196118064226*GCL[32])-1.533921189023155*GCC[32]+0.303196118064226*GBL[32]-0.303196118064226*GBC[32]-0.21544659739277597*GCL[21]-0.41432037960149226*GCC[21]+0.21544659739277597*GBL[21]-0.6463397921783279*GBC[21]+0.21544659739277597*GCL[19]+0.8452135743870443*GCC[19]+0.21544659739277597*GBL[19]-0.21544659739277597*GBC[19]+0.15309310892394856*GCL[11]+0.15309310892394856*GCC[11]+0.15309310892394856*GBL[11]-0.4592793267718458*GBC[11]; 
  surft1_lo[12] = -(0.27209908031404895*GCL[34])+1.9943965557084689*GCC[34]-0.27209908031404895*GBL[34]-1.0674656227704997*GBC[34]+0.27209908031404895*GCL[33]-2.789763098164919*GCC[33]-0.27209908031404895*GBL[33]-0.27209908031404895*GBC[33]-0.1933495104806964*GCL[23]+1.5744174424856705*GCC[23]-0.1933495104806964*GBL[23]-0.1933495104806964*GBC[23]-0.1933495104806964*GCL[20]+1.5744174424856705*GCC[20]-0.1933495104806964*GBL[20]-0.1933495104806964*GBC[20]+0.39142450529916156*GCL[15]-0.8361198012603392*GCC[15]-0.39142450529916156*GBL[15]-1.5355884438659415*GBC[15]+0.2781403612330919*GCL[6]-1.091166032529822*GCC[6]-0.2781403612330919*GBL[6]-0.2781403612330919*GBC[6]-0.2781403612330919*GCL[5]+0.2781403612330919*GCC[5]-0.2781403612330919*GBL[5]-1.091166032529822*GBC[5]-0.19764235376052366*GCL[1]+0.5929270612815709*GCC[1]-0.19764235376052366*GBL[1]-0.19764235376052366*GBC[1]; 
  surft1_lo[13] = 0.1750503603816304*GCL[44]-1.23571053216145*GCC[44]-0.1750503603816304*GBL[44]+0.1750503603816304*GBC[44]+0.12438815100070813*GCL[37]-0.48798428469508576*GCC[37]-0.12438815100070813*GBL[37]-0.12438815100070813*GBC[37]-0.12438815100070813*GCL[35]+0.736760586696502*GCC[35]-0.12438815100070813*GBL[35]+0.12438815100070813*GBC[35]-0.0883883476483184*GCL[25]+0.26516504294495524*GCC[25]-0.0883883476483184*GBL[25]-0.0883883476483184*GBC[25]; 
  surft1_lo[14] = -(0.27209908031404895*GCL[40])+1.9943965557084689*GCC[40]-0.27209908031404895*GBL[40]-1.0674656227704997*GBC[40]+0.27209908031404895*GCL[38]-2.789763098164919*GCC[38]-0.27209908031404895*GBL[38]-0.27209908031404895*GBC[38]-0.1933495104806964*GCL[27]+1.5744174424856705*GCC[27]-0.1933495104806964*GBL[27]-0.1933495104806964*GBC[27]-0.1933495104806964*GCL[26]+1.5744174424856705*GCC[26]-0.1933495104806964*GBL[26]-0.1933495104806964*GBC[26]+0.39142450529916156*GCL[18]-0.8361198012603392*GCC[18]-0.39142450529916156*GBL[18]-1.5355884438659415*GBC[18]+0.2781403612330919*GCL[10]-1.091166032529822*GCC[10]-0.2781403612330919*GBL[10]-0.2781403612330919*GBC[10]-0.2781403612330919*GCL[9]+0.2781403612330919*GCC[9]-0.2781403612330919*GBL[9]-1.091166032529822*GBC[9]-0.19764235376052366*GCL[4]+0.5929270612815709*GCC[4]-0.19764235376052366*GBL[4]-0.19764235376052366*GBC[4]; 
  surft1_lo[15] = 0.1750503603816304*GCL[47]-1.23571053216145*GCC[47]-0.1750503603816304*GBL[47]+0.1750503603816304*GBC[47]+0.12438815100070813*GCL[42]-0.48798428469508576*GCC[42]-0.12438815100070813*GBL[42]-0.12438815100070813*GBC[42]-0.12438815100070813*GCL[41]+0.736760586696502*GCC[41]-0.12438815100070813*GBL[41]+0.12438815100070813*GBC[41]-0.0883883476483184*GCL[28]+0.26516504294495524*GCC[28]-0.0883883476483184*GBL[28]-0.0883883476483184*GBC[28]; 
  surft1_lo[16] = -(0.303196118064226*GCL[43])-1.533921189023155*GCC[43]+0.303196118064226*GBL[43]-0.303196118064226*GBC[43]-0.21544659739277597*GCL[30]-0.41432037960149226*GCC[30]+0.21544659739277597*GBL[30]-0.6463397921783279*GBC[30]+0.21544659739277597*GCL[29]+0.8452135743870443*GCC[29]+0.21544659739277597*GBL[29]-0.21544659739277597*GBC[29]+0.15309310892394856*GCL[14]+0.15309310892394856*GCC[14]+0.15309310892394856*GBL[14]-0.4592793267718458*GBC[14]; 
  surft1_lo[17] = -(0.303196118064226*GCL[44])-1.533921189023155*GCC[44]+0.303196118064226*GBL[44]-0.303196118064226*GBC[44]-0.21544659739277597*GCL[37]-0.41432037960149226*GCC[37]+0.21544659739277597*GBL[37]-0.6463397921783279*GBC[37]+0.21544659739277597*GCL[35]+0.8452135743870443*GCC[35]+0.21544659739277597*GBL[35]-0.21544659739277597*GBC[35]+0.15309310892394856*GCL[25]+0.15309310892394856*GCC[25]+0.15309310892394856*GBL[25]-0.4592793267718458*GBC[25]; 
  surft1_lo[18] = -(0.27209908031404895*GCL[46])+1.9943965557084689*GCC[46]-0.27209908031404895*GBL[46]-1.0674656227704997*GBC[46]+0.27209908031404895*GCL[45]-2.789763098164919*GCC[45]-0.27209908031404895*GBL[45]-0.27209908031404895*GBC[45]-0.1933495104806964*GCL[39]+1.5744174424856705*GCC[39]-0.1933495104806964*GBL[39]-0.1933495104806964*GBC[39]-0.1933495104806964*GCL[36]+1.5744174424856705*GCC[36]-0.1933495104806964*GBL[36]-0.1933495104806964*GBC[36]+0.3914245052991616*GCL[31]-0.8361198012603392*GCC[31]-0.3914245052991616*GBL[31]-1.5355884438659417*GBC[31]+0.2781403612330919*GCL[17]-1.0911660325298218*GCC[17]-0.2781403612330919*GBL[17]-0.2781403612330919*GBC[17]-0.2781403612330919*GCL[16]+0.2781403612330919*GCC[16]-0.2781403612330919*GBL[16]-1.0911660325298218*GBC[16]-0.19764235376052364*GCL[8]+0.592927061281571*GCC[8]-0.19764235376052364*GBL[8]-0.19764235376052364*GBC[8]; 
  surft1_lo[19] = -(0.303196118064226*GCL[47])-1.533921189023155*GCC[47]+0.303196118064226*GBL[47]-0.303196118064226*GBC[47]-0.21544659739277597*GCL[42]-0.41432037960149226*GCC[42]+0.21544659739277597*GBL[42]-0.6463397921783279*GBC[42]+0.21544659739277597*GCL[41]+0.8452135743870443*GCC[41]+0.21544659739277597*GBL[41]-0.21544659739277597*GBC[41]+0.15309310892394856*GCL[28]+0.15309310892394856*GCC[28]+0.15309310892394856*GBL[28]-0.4592793267718458*GBC[28]; 
  surft1_up[0] = dGdvx_surf_CC_vy[0]/dv2; 
  surft1_up[1] = dGdvx_surf_CC_vy[1]/dv2; 
  surft1_up[2] = dGdvx_surf_CC_vy[2]/dv2; 
  surft1_up[3] = dGdvx_surf_CC_vy[3]/dv2; 
  surft1_up[4] = dGdvx_surf_CC_vy[4]/dv2; 
  surft1_up[5] = dGdvx_surf_CC_vy[5]/dv2; 
  surft1_up[6] = dGdvx_surf_CC_vy[6]/dv2; 
  surft1_up[7] = dGdvx_surf_CC_vy[7]/dv2; 
  surft1_up[8] = dGdvx_surf_CC_vy[8]/dv2; 
  surft1_up[9] = dGdvx_surf_CC_vy[9]/dv2; 
  surft1_up[10] = dGdvx_surf_CC_vy[10]/dv2; 
  surft1_up[11] = dGdvx_surf_CC_vy[11]/dv2; 
  surft1_up[12] = dGdvx_surf_CC_vy[12]/dv2; 
  surft1_up[13] = dGdvx_surf_CC_vy[13]/dv2; 
  surft1_up[14] = dGdvx_surf_CC_vy[14]/dv2; 
  surft1_up[15] = dGdvx_surf_CC_vy[15]/dv2; 
  surft1_up[16] = dGdvx_surf_CC_vy[16]/dv2; 
  surft1_up[17] = dGdvx_surf_CC_vy[17]/dv2; 
  surft1_up[18] = dGdvx_surf_CC_vy[18]/dv2; 
  surft1_up[19] = dGdvx_surf_CC_vy[19]/dv2; 

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
  surft2_up[0] = G_surf_CC_vx[0]; 
  surft2_up[1] = G_surf_CC_vx[1]; 
  surft2_up[2] = G_surf_CC_vx[2]; 
  surft2_up[3] = G_surf_CC_vx[3]; 
  surft2_up[4] = G_surf_CC_vx[4]; 
  surft2_up[5] = G_surf_CC_vx[5]; 
  surft2_up[6] = G_surf_CC_vx[6]; 
  surft2_up[7] = G_surf_CC_vx[7]; 
  surft2_up[8] = G_surf_CC_vx[8]; 
  surft2_up[9] = G_surf_CC_vx[9]; 
  surft2_up[10] = G_surf_CC_vx[10]; 
  surft2_up[11] = G_surf_CC_vx[11]; 
  surft2_up[12] = G_surf_CC_vx[12]; 
  surft2_up[13] = G_surf_CC_vx[13]; 
  surft2_up[14] = G_surf_CC_vx[14]; 
  surft2_up[15] = G_surf_CC_vx[15]; 
  surft2_up[16] = G_surf_CC_vx[16]; 
  surft2_up[17] = G_surf_CC_vx[17]; 
  surft2_up[18] = G_surf_CC_vx[18]; 
  surft2_up[19] = G_surf_CC_vx[19]; 

  vol[7] = 3.0*GCC[0]; 
  vol[15] = 3.0*GCC[1]; 
  vol[18] = 3.0*GCC[4]; 
  vol[22] = 6.7082039324993685*GCC[2]; 
  vol[24] = 6.7082039324993685*GCC[3]; 
  vol[31] = 3.0*GCC[8]; 
  vol[32] = 3.0*GCC[11]; 
  vol[33] = 6.708203932499369*GCC[5]; 
  vol[34] = 6.708203932499369*GCC[6]; 
  vol[38] = 6.708203932499369*GCC[9]; 
  vol[40] = 6.708203932499369*GCC[10]; 
  vol[43] = 3.0*GCC[14]; 
  vol[44] = 3.0*GCC[25]; 
  vol[45] = 6.7082039324993685*GCC[16]; 
  vol[46] = 6.7082039324993685*GCC[17]; 
  vol[47] = 3.0*GCC[28]; 

  out[0] = (vol[0]+0.7071067811865475*surft1_up[0]-0.7071067811865475*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[1] = (vol[1]+0.7071067811865475*surft1_up[1]-0.7071067811865475*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[2] = (vol[2]+0.7071067811865475*surft1_up[2]-0.7071067811865475*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[3] = (vol[3]-1.224744871391589*surft2_up[0]+1.224744871391589*surft2_lo[0]+1.224744871391589*surft1_up[0]+1.224744871391589*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[4] = (vol[4]+0.7071067811865475*surft1_up[3]-0.7071067811865475*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[5] = (vol[5]+0.7071067811865475*surft1_up[4]-0.7071067811865475*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[6] = (vol[6]-1.224744871391589*surft2_up[1]+1.224744871391589*surft2_lo[1]+1.224744871391589*surft1_up[1]+1.224744871391589*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[7] = (vol[7]+1.224744871391589*surft1_up[2]+1.224744871391589*surft1_lo[2]-2.1213203435596424*surft2_up[0]-2.1213203435596424*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[8] = (vol[8]+0.7071067811865475*surft1_up[5]-0.7071067811865475*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[9] = (vol[9]+0.7071067811865475*surft1_up[6]-0.7071067811865475*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[10] = (vol[10]-1.224744871391589*surft2_up[3]+1.224744871391589*surft2_lo[3]+1.224744871391589*surft1_up[3]+1.224744871391589*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[11] = (vol[11]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[12] = (vol[12]+0.7071067811865475*surft1_up[8]-0.7071067811865475*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[13] = (vol[13]-2.7386127875258306*surft2_up[2]+2.7386127875258306*surft2_lo[2]+1.5811388300841895*surft1_up[0]-1.5811388300841895*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]+0.7071067811865475*surft1_up[9]-0.7071067811865475*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+1.224744871391589*surft1_up[4]+1.224744871391589*surft1_lo[4]-2.1213203435596424*surft2_up[1]-2.1213203435596424*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[16] = (vol[16]+0.7071067811865475*surft1_up[10]-0.7071067811865475*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[17] = (vol[17]-1.224744871391589*surft2_up[5]+1.224744871391589*surft2_lo[5]+1.224744871391589*surft1_up[5]+1.224744871391589*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[18] = (vol[18]+1.224744871391589*surft1_up[6]+1.224744871391589*surft1_lo[6]-2.1213203435596424*surft2_up[3]-2.1213203435596424*surft2_lo[3])*dv1_sq*gamma_avg; 
  out[19] = (vol[19]+0.7071067811865475*surft1_up[11]-0.7071067811865475*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[20] = (vol[20]+0.7071067811865475*surft1_up[12]-0.7071067811865475*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[21] = (vol[21]-1.224744871391589*surft2_up[7]+1.224744871391589*surft2_lo[7]+1.224744871391589*surft1_up[7]+1.224744871391589*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[22] = (vol[22]+1.224744871391589*surft1_up[8]+1.224744871391589*surft1_lo[8]-2.7386127875258306*surft2_up[0]+2.7386127875258306*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[23] = (vol[23]-2.7386127875258306*surft2_up[4]+2.7386127875258306*surft2_lo[4]+1.5811388300841898*surft1_up[1]-1.5811388300841898*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[24] = (vol[24]-4.743416490252569*surft2_up[2]-4.743416490252569*surft2_lo[2]+1.5811388300841898*surft1_up[2]-1.5811388300841898*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[25] = (vol[25]+0.7071067811865475*surft1_up[13]-0.7071067811865475*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[26] = (vol[26]+0.7071067811865475*surft1_up[14]-0.7071067811865475*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[27] = (vol[27]-2.7386127875258306*surft2_up[6]+2.7386127875258306*surft2_lo[6]+1.5811388300841898*surft1_up[3]-1.5811388300841898*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[28] = (vol[28]+0.7071067811865475*surft1_up[15]-0.7071067811865475*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[29] = (vol[29]+0.7071067811865475*surft1_up[16]-0.7071067811865475*surft1_lo[16])*dv1_sq*gamma_avg; 
  out[30] = (vol[30]-1.224744871391589*surft2_up[9]+1.224744871391589*surft2_lo[9]+1.224744871391589*surft1_up[9]+1.224744871391589*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[31] = (vol[31]+1.224744871391589*surft1_up[10]+1.224744871391589*surft1_lo[10]-2.1213203435596424*surft2_up[5]-2.1213203435596424*surft2_lo[5])*dv1_sq*gamma_avg; 
  out[32] = (vol[32]+1.224744871391589*surft1_up[11]+1.224744871391589*surft1_lo[11]-2.1213203435596424*surft2_up[7]-2.1213203435596424*surft2_lo[7])*dv1_sq*gamma_avg; 
  out[33] = (vol[33]+1.224744871391589*surft1_up[12]+1.224744871391589*surft1_lo[12]-2.7386127875258306*surft2_up[1]+2.7386127875258306*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[34] = (vol[34]-4.743416490252569*surft2_up[4]-4.743416490252569*surft2_lo[4]+1.5811388300841895*surft1_up[4]-1.5811388300841895*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[35] = (vol[35]+0.7071067811865475*surft1_up[17]-0.7071067811865475*surft1_lo[17])*dv1_sq*gamma_avg; 
  out[36] = (vol[36]+0.7071067811865475*surft1_up[18]-0.7071067811865475*surft1_lo[18])*dv1_sq*gamma_avg; 
  out[37] = (vol[37]-1.224744871391589*surft2_up[13]+1.224744871391589*surft2_lo[13]+1.224744871391589*surft1_up[13]+1.224744871391589*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[38] = (vol[38]+1.224744871391589*surft1_up[14]+1.224744871391589*surft1_lo[14]-2.7386127875258306*surft2_up[3]+2.7386127875258306*surft2_lo[3])*dv1_sq*gamma_avg; 
  out[39] = (vol[39]-2.7386127875258306*surft2_up[10]+2.7386127875258306*surft2_lo[10]+1.5811388300841895*surft1_up[5]-1.5811388300841895*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[40] = (vol[40]-4.743416490252569*surft2_up[6]-4.743416490252569*surft2_lo[6]+1.5811388300841895*surft1_up[6]-1.5811388300841895*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[41] = (vol[41]+0.7071067811865475*surft1_up[19]-0.7071067811865475*surft1_lo[19])*dv1_sq*gamma_avg; 
  out[42] = (vol[42]-1.224744871391589*surft2_up[15]+1.224744871391589*surft2_lo[15]+1.224744871391589*surft1_up[15]+1.224744871391589*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[43] = (vol[43]+1.224744871391589*surft1_up[16]+1.224744871391589*surft1_lo[16]-2.1213203435596424*surft2_up[9]-2.1213203435596424*surft2_lo[9])*dv1_sq*gamma_avg; 
  out[44] = (vol[44]+1.224744871391589*surft1_up[17]+1.224744871391589*surft1_lo[17]-2.1213203435596424*surft2_up[13]-2.1213203435596424*surft2_lo[13])*dv1_sq*gamma_avg; 
  out[45] = (vol[45]+1.224744871391589*surft1_up[18]+1.224744871391589*surft1_lo[18]-2.7386127875258306*surft2_up[5]+2.7386127875258306*surft2_lo[5])*dv1_sq*gamma_avg; 
  out[46] = (vol[46]-4.743416490252569*surft2_up[10]-4.743416490252569*surft2_lo[10]+1.5811388300841898*surft1_up[10]-1.5811388300841898*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[47] = (vol[47]+1.224744871391589*surft1_up[19]+1.224744871391589*surft1_lo[19]-2.1213203435596424*surft2_up[15]-2.1213203435596424*surft2_lo[15])*dv1_sq*gamma_avg; 
} 