#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p2_lovz_lovx(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
 double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1_pv1 = 2.0/dxv[3]; 
  double dv1_pv2 = 2.0/dxv[1]; 
  double dv1_sq = 4.0/dxv[3]/dxv[1]; 
 
  const double* GCC = fpo_g_stencil[0]; 
  const double* GCR = fpo_g_stencil[1]; 
  const double* GTC = fpo_g_stencil[2]; 
  const double* GTR = fpo_g_stencil[3]; 

  const double* g_surf_CC = fpo_g_surf_stencil[0]; 
  const double* g_surf_CC_pv2 = &g_surf_CC[0]; 
  const double* g_surf_CR = fpo_g_surf_stencil[1]; 
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
  
  surft1_upper[0] = (0.11180339887498948*dgdpv1_surf_CC_pv1[8])/dv1_pv1+(0.08660254037844387*dgdpv1_surf_CC_pv1[2])/dv1_pv1+(0.05*dgdpv1_surf_CC_pv1[0])/dv1_pv1-0.20999597273098428*GTR[29]+0.5020433520257306*GTC[29]+0.20999597273098428*GCR[29]-0.5020433520257306*GCC[29]-0.20218664720405521*GTR[26]+0.2845589849538555*GTC[26]-0.20218664720405521*GCR[26]+0.2845589849538555*GCC[26]+0.14921997708919524*GTR[14]-0.35674444853774495*GTC[14]+0.14921997708919524*GCR[14]-0.35674444853774495*GCC[14]+0.14051136087662214*GTR[12]-0.14051136087662214*GTC[12]+0.14051136087662214*GCR[12]-0.14051136087662214*GCC[12]+0.2908529064802475*GTR[9]-0.4093485350462743*GTC[9]-0.2908529064802475*GCR[9]+0.4093485350462743*GCC[9]-0.20667569704733044*GTR[4]+0.29087690695550217*GTC[4]-0.20667569704733044*GCR[4]+0.29087690695550217*GCC[4]-0.2021307453761506*GTR[2]+0.2021307453761506*GTC[2]+0.2021307453761506*GCR[2]-0.2021307453761506*GCC[2]+0.14363106492851735*GTR[0]-0.14363106492851735*GTC[0]+0.14363106492851735*GCR[0]-0.14363106492851735*GCC[0]; 
  surft1_upper[1] = (0.11180339887498951*dgdpv1_surf_CC_pv1[12])/dv1_pv1+(0.08660254037844387*dgdpv1_surf_CC_pv1[4])/dv1_pv1+(0.05*dgdpv1_surf_CC_pv1[1])/dv1_pv1-0.20999597273098425*GTR[41]+0.5020433520257305*GTC[41]+0.20999597273098425*GCR[41]-0.5020433520257305*GCC[41]-0.2021866472040551*GTR[36]+0.28455898495385545*GTC[36]-0.2021866472040551*GCR[36]+0.28455898495385545*GCC[36]+0.1492199770891953*GTR[28]-0.35674444853774506*GTC[28]+0.1492199770891953*GCR[28]-0.35674444853774506*GCC[28]+0.14051136087662214*GTR[20]-0.14051136087662214*GTC[20]+0.14051136087662214*GCR[20]-0.14051136087662214*GCC[20]+0.2908529064802475*GTR[16]-0.4093485350462743*GTC[16]-0.2908529064802475*GCR[16]+0.4093485350462743*GCC[16]-0.20667569704733044*GTR[8]+0.29087690695550217*GTC[8]-0.20667569704733044*GCR[8]+0.29087690695550217*GCC[8]-0.2021307453761506*GTR[5]+0.2021307453761506*GTC[5]+0.2021307453761506*GCR[5]-0.2021307453761506*GCC[5]+0.14363106492851735*GTR[1]-0.14363106492851735*GTC[1]+0.14363106492851735*GCR[1]-0.14363106492851735*GCC[1]; 
  surft1_upper[2] = (0.11180339887498951*dgdpv1_surf_CC_pv1[14])/dv1_pv1+(0.08660254037844387*dgdpv1_surf_CC_pv1[6])/dv1_pv1+(0.05*dgdpv1_surf_CC_pv1[3])/dv1_pv1-0.20999597273098425*GTR[43]+0.5020433520257305*GTC[43]+0.20999597273098425*GCR[43]-0.5020433520257305*GCC[43]-0.2021866472040551*GTR[38]+0.28455898495385545*GTC[38]-0.2021866472040551*GCR[38]+0.28455898495385545*GCC[38]+0.1492199770891953*GTR[30]-0.35674444853774506*GTC[30]+0.1492199770891953*GCR[30]-0.35674444853774506*GCC[30]+0.14051136087662214*GTR[22]-0.14051136087662214*GTC[22]+0.14051136087662214*GCR[22]-0.14051136087662214*GCC[22]+0.2908529064802475*GTR[18]-0.4093485350462743*GTC[18]-0.2908529064802475*GCR[18]+0.4093485350462743*GCC[18]-0.20667569704733044*GTR[10]+0.29087690695550217*GTC[10]-0.20667569704733044*GCR[10]+0.29087690695550217*GCC[10]-0.2021307453761506*GTR[7]+0.2021307453761506*GTC[7]+0.2021307453761506*GCR[7]-0.2021307453761506*GCC[7]+0.14363106492851735*GTR[3]-0.14363106492851735*GTC[3]+0.14363106492851735*GCR[3]-0.14363106492851735*GCC[3]; 
  surft1_upper[3] = -((0.19364916731037082*dgdpv1_surf_CC_pv1[8])/dv1_pv1)-(0.15*dgdpv1_surf_CC_pv1[2])/dv1_pv1-(0.08660254037844387*dgdpv1_surf_CC_pv1[0])/dv1_pv1-0.05781038847495312*GTR[29]-1.2910986759406198*GTC[29]+0.05781038847495312*GCR[29]+1.2910986759406198*GCC[29]-0.07133653706043892*GTR[26]-0.07133653706043892*GTC[26]-0.07133653706043892*GCR[26]-0.07133653706043892*GCC[26]+0.04107919181288743*GTR[14]+0.9174352838211527*GTC[14]+0.04107919181288743*GCR[14]+0.9174352838211527*GCC[14]+0.05616295755668199*GTR[12]-0.05616295755668199*GTC[12]+0.05616295755668199*GCR[12]-0.05616295755668199*GCC[12]+0.10262022457558416*GTR[9]+0.10262022457558416*GTC[9]-0.10262022457558416*GCR[9]-0.10262022457558416*GCC[9]-0.07292038680986264*GTR[4]-0.07292038680986264*GTC[4]-0.07292038680986264*GCR[4]-0.07292038680986264*GCC[4]-0.08079247402229095*GTR[2]+0.08079247402229095*GTC[2]+0.08079247402229095*GCR[2]-0.08079247402229095*GCC[2]+0.057409915846480676*GTR[0]-0.057409915846480676*GTC[0]+0.057409915846480676*GCR[0]-0.057409915846480676*GCC[0]; 
  surft1_upper[4] = (0.11180339887498948*dgdpv1_surf_CC_pv1[18])/dv1_pv1+(0.08660254037844387*dgdpv1_surf_CC_pv1[10])/dv1_pv1+(0.05*dgdpv1_surf_CC_pv1[5])/dv1_pv1-0.20999597273098428*GTR[47]+0.5020433520257306*GTC[47]+0.20999597273098428*GCR[47]-0.5020433520257306*GCC[47]-0.20218664720405521*GTR[45]+0.2845589849538555*GTC[45]-0.20218664720405521*GCR[45]+0.2845589849538555*GCC[45]+0.14921997708919524*GTR[42]-0.35674444853774495*GTC[42]+0.14921997708919524*GCR[42]-0.35674444853774495*GCC[42]+0.14051136087662214*GTR[33]-0.14051136087662214*GTC[33]+0.14051136087662214*GCR[33]-0.14051136087662214*GCC[33]+0.2908529064802475*GTR[31]-0.4093485350462743*GTC[31]-0.2908529064802475*GCR[31]+0.4093485350462743*GCC[31]-0.20667569704733044*GTR[17]+0.29087690695550217*GTC[17]-0.20667569704733044*GCR[17]+0.29087690695550217*GCC[17]-0.2021307453761506*GTR[15]+0.2021307453761506*GTC[15]+0.2021307453761506*GCR[15]-0.2021307453761506*GCC[15]+0.14363106492851735*GTR[6]-0.14363106492851735*GTC[6]+0.14363106492851735*GCR[6]-0.14363106492851735*GCC[6]; 
  surft1_upper[5] = -((0.19364916731037085*dgdpv1_surf_CC_pv1[12])/dv1_pv1)-(0.15*dgdpv1_surf_CC_pv1[4])/dv1_pv1-(0.08660254037844387*dgdpv1_surf_CC_pv1[1])/dv1_pv1-0.057810388474953116*GTR[41]-1.2910986759406196*GTC[41]+0.057810388474953116*GCR[41]+1.2910986759406196*GCC[41]-0.07133653706043892*GTR[36]-0.07133653706043892*GTC[36]-0.07133653706043892*GCR[36]-0.07133653706043892*GCC[36]+0.04107919181288744*GTR[28]+0.9174352838211529*GTC[28]+0.04107919181288744*GCR[28]+0.9174352838211529*GCC[28]+0.05616295755668199*GTR[20]-0.05616295755668199*GTC[20]+0.05616295755668199*GCR[20]-0.05616295755668199*GCC[20]+0.10262022457558416*GTR[16]+0.10262022457558416*GTC[16]-0.10262022457558416*GCR[16]-0.10262022457558416*GCC[16]-0.07292038680986264*GTR[8]-0.07292038680986264*GTC[8]-0.07292038680986264*GCR[8]-0.07292038680986264*GCC[8]-0.08079247402229095*GTR[5]+0.08079247402229095*GTC[5]+0.08079247402229095*GCR[5]-0.08079247402229095*GCC[5]+0.057409915846480676*GTR[1]-0.057409915846480676*GTC[1]+0.057409915846480676*GCR[1]-0.057409915846480676*GCC[1]; 
  surft1_upper[6] = -((0.19364916731037085*dgdpv1_surf_CC_pv1[14])/dv1_pv1)-(0.15*dgdpv1_surf_CC_pv1[6])/dv1_pv1-(0.08660254037844387*dgdpv1_surf_CC_pv1[3])/dv1_pv1-0.057810388474953116*GTR[43]-1.2910986759406196*GTC[43]+0.057810388474953116*GCR[43]+1.2910986759406196*GCC[43]-0.07133653706043892*GTR[38]-0.07133653706043892*GTC[38]-0.07133653706043892*GCR[38]-0.07133653706043892*GCC[38]+0.04107919181288744*GTR[30]+0.9174352838211529*GTC[30]+0.04107919181288744*GCR[30]+0.9174352838211529*GCC[30]+0.05616295755668199*GTR[22]-0.05616295755668199*GTC[22]+0.05616295755668199*GCR[22]-0.05616295755668199*GCC[22]+0.10262022457558416*GTR[18]+0.10262022457558416*GTC[18]-0.10262022457558416*GCR[18]-0.10262022457558416*GCC[18]-0.07292038680986264*GTR[10]-0.07292038680986264*GTC[10]-0.07292038680986264*GCR[10]-0.07292038680986264*GCC[10]-0.08079247402229095*GTR[7]+0.08079247402229095*GTC[7]+0.08079247402229095*GCR[7]-0.08079247402229095*GCC[7]+0.057409915846480676*GTR[3]-0.057409915846480676*GTC[3]+0.057409915846480676*GCR[3]-0.057409915846480676*GCC[3]; 
  surft1_upper[7] = (0.08660254037844385*dgdpv1_surf_CC_pv1[11])/dv1_pv1+(0.05*dgdpv1_surf_CC_pv1[7])/dv1_pv1+0.2908529064802475*GTR[35]-0.4093485350462743*GTC[35]-0.2908529064802475*GCR[35]+0.4093485350462743*GCC[35]-0.2066756970473305*GTR[25]+0.2908769069555022*GTC[25]-0.2066756970473305*GCR[25]+0.2908769069555022*GCC[25]-0.2021307453761506*GTR[19]+0.2021307453761506*GTC[19]+0.2021307453761506*GCR[19]-0.2021307453761506*GCC[19]+0.14363106492851735*GTR[11]-0.14363106492851735*GTC[11]+0.14363106492851735*GCR[11]-0.14363106492851735*GCC[11]; 
  surft1_upper[8] = (0.08660254037844385*dgdpv1_surf_CC_pv1[16])/dv1_pv1+(0.05*dgdpv1_surf_CC_pv1[9])/dv1_pv1+0.2908529064802475*GTR[40]-0.4093485350462743*GTC[40]-0.2908529064802475*GCR[40]+0.4093485350462743*GCC[40]-0.2066756970473305*GTR[27]+0.2908769069555022*GTC[27]-0.2066756970473305*GCR[27]+0.2908769069555022*GCC[27]-0.2021307453761506*GTR[24]+0.2021307453761506*GTC[24]+0.2021307453761506*GCR[24]-0.2021307453761506*GCC[24]+0.14363106492851735*GTR[13]-0.14363106492851735*GTC[13]+0.14363106492851735*GCR[13]-0.14363106492851735*GCC[13]; 
  surft1_upper[9] = (0.25*dgdpv1_surf_CC_pv1[8])/dv1_pv1+(0.19364916731037082*dgdpv1_surf_CC_pv1[2])/dv1_pv1+(0.11180339887498948*dgdpv1_surf_CC_pv1[0])/dv1_pv1-0.4695652700276729*GTR[29]+1.1226030627813903*GTC[29]+0.4695652700276729*GCR[29]-1.1226030627813903*GCC[29]-0.4521030872910352*GTR[26]-0.7032714691193882*GTC[26]-0.4521030872910352*GCR[26]-0.7032714691193882*GCC[26]+0.3336660123724018*GTR[14]-0.7977048375260732*GTC[14]+0.3336660123724018*GCR[14]-0.7977048375260732*GCC[14]+0.3141929545311315*GTR[12]-0.3141929545311315*GTC[12]+0.3141929545311315*GCR[12]-0.3141929545311315*GCC[12]+0.6503668703432224*GTR[9]+1.0116817983116795*GTC[9]-0.6503668703432224*GCR[9]-1.0116817983116795*GCC[9]-0.46214090789498335*GTR[4]-0.71888585672553*GTC[4]-0.46214090789498335*GCR[4]-0.71888585672553*GCC[4]-0.45197808700377406*GTR[2]+0.45197808700377406*GTC[2]+0.45197808700377406*GCR[2]-0.45197808700377406*GCC[2]+0.3211688248608508*GTR[0]-0.3211688248608508*GTC[0]+0.3211688248608508*GCR[0]-0.3211688248608508*GCC[0]; 
  surft1_upper[10] = -((0.19364916731037082*dgdpv1_surf_CC_pv1[18])/dv1_pv1)-(0.15*dgdpv1_surf_CC_pv1[10])/dv1_pv1-(0.08660254037844387*dgdpv1_surf_CC_pv1[5])/dv1_pv1-0.05781038847495312*GTR[47]-1.2910986759406198*GTC[47]+0.05781038847495312*GCR[47]+1.2910986759406198*GCC[47]-0.07133653706043892*GTR[45]-0.07133653706043892*GTC[45]-0.07133653706043892*GCR[45]-0.07133653706043892*GCC[45]+0.04107919181288743*GTR[42]+0.9174352838211527*GTC[42]+0.04107919181288743*GCR[42]+0.9174352838211527*GCC[42]+0.05616295755668199*GTR[33]-0.05616295755668199*GTC[33]+0.05616295755668199*GCR[33]-0.05616295755668199*GCC[33]+0.10262022457558416*GTR[31]+0.10262022457558416*GTC[31]-0.10262022457558416*GCR[31]-0.10262022457558416*GCC[31]-0.07292038680986264*GTR[17]-0.07292038680986264*GTC[17]-0.07292038680986264*GCR[17]-0.07292038680986264*GCC[17]-0.08079247402229095*GTR[15]+0.08079247402229095*GTC[15]+0.08079247402229095*GCR[15]-0.08079247402229095*GCC[15]+0.057409915846480676*GTR[6]-0.057409915846480676*GTC[6]+0.057409915846480676*GCR[6]-0.057409915846480676*GCC[6]; 
  surft1_upper[11] = (0.08660254037844385*dgdpv1_surf_CC_pv1[17])/dv1_pv1+(0.05*dgdpv1_surf_CC_pv1[13])/dv1_pv1+0.2908529064802475*GTR[44]-0.4093485350462743*GTC[44]-0.2908529064802475*GCR[44]+0.4093485350462743*GCC[44]-0.2066756970473305*GTR[37]+0.2908769069555022*GTC[37]-0.2066756970473305*GCR[37]+0.2908769069555022*GCC[37]-0.2021307453761506*GTR[32]+0.2021307453761506*GTC[32]+0.2021307453761506*GCR[32]-0.2021307453761506*GCC[32]+0.14363106492851735*GTR[21]-0.14363106492851735*GTC[21]+0.14363106492851735*GCR[21]-0.14363106492851735*GCC[21]; 
  surft1_upper[12] = (0.08660254037844385*dgdpv1_surf_CC_pv1[19])/dv1_pv1+(0.05*dgdpv1_surf_CC_pv1[15])/dv1_pv1+0.2908529064802475*GTR[46]-0.4093485350462743*GTC[46]-0.2908529064802475*GCR[46]+0.4093485350462743*GCC[46]-0.2066756970473305*GTR[39]+0.2908769069555022*GTC[39]-0.2066756970473305*GCR[39]+0.2908769069555022*GCC[39]-0.2021307453761506*GTR[34]+0.2021307453761506*GTC[34]+0.2021307453761506*GCR[34]-0.2021307453761506*GCC[34]+0.14363106492851735*GTR[23]-0.14363106492851735*GTC[23]+0.14363106492851735*GCR[23]-0.14363106492851735*GCC[23]; 
  surft1_upper[13] = -((0.15*dgdpv1_surf_CC_pv1[11])/dv1_pv1)-(0.08660254037844385*dgdpv1_surf_CC_pv1[7])/dv1_pv1+0.10262022457558415*GTR[35]+0.10262022457558415*GTC[35]-0.10262022457558415*GCR[35]-0.10262022457558415*GCC[35]-0.07292038680986264*GTR[25]-0.07292038680986264*GTC[25]-0.07292038680986264*GCR[25]-0.07292038680986264*GCC[25]-0.08079247402229095*GTR[19]+0.08079247402229095*GTC[19]+0.08079247402229095*GCR[19]-0.08079247402229095*GCC[19]+0.05740991584648068*GTR[11]-0.05740991584648068*GTC[11]+0.05740991584648068*GCR[11]-0.05740991584648068*GCC[11]; 
  surft1_upper[14] = -((0.15*dgdpv1_surf_CC_pv1[16])/dv1_pv1)-(0.08660254037844385*dgdpv1_surf_CC_pv1[9])/dv1_pv1+0.10262022457558415*GTR[40]+0.10262022457558415*GTC[40]-0.10262022457558415*GCR[40]-0.10262022457558415*GCC[40]-0.07292038680986264*GTR[27]-0.07292038680986264*GTC[27]-0.07292038680986264*GCR[27]-0.07292038680986264*GCC[27]-0.08079247402229095*GTR[24]+0.08079247402229095*GTC[24]+0.08079247402229095*GCR[24]-0.08079247402229095*GCC[24]+0.05740991584648068*GTR[13]-0.05740991584648068*GTC[13]+0.05740991584648068*GCR[13]-0.05740991584648068*GCC[13]; 
  surft1_upper[15] = (0.25*dgdpv1_surf_CC_pv1[12])/dv1_pv1+(0.19364916731037085*dgdpv1_surf_CC_pv1[4])/dv1_pv1+(0.11180339887498951*dgdpv1_surf_CC_pv1[1])/dv1_pv1-0.4695652700276729*GTR[41]+1.1226030627813903*GTC[41]+0.4695652700276729*GCR[41]-1.1226030627813903*GCC[41]-0.4521030872910352*GTR[36]-0.7032714691193882*GTC[36]-0.4521030872910352*GCR[36]-0.7032714691193882*GCC[36]+0.3336660123724018*GTR[28]-0.7977048375260732*GTC[28]+0.3336660123724018*GCR[28]-0.7977048375260732*GCC[28]+0.3141929545311315*GTR[20]-0.3141929545311315*GTC[20]+0.3141929545311315*GCR[20]-0.3141929545311315*GCC[20]+0.6503668703432222*GTR[16]+1.0116817983116793*GTC[16]-0.6503668703432222*GCR[16]-1.0116817983116793*GCC[16]-0.46214090789498363*GTR[8]-0.7188858567255303*GTC[8]-0.46214090789498363*GCR[8]-0.7188858567255303*GCC[8]-0.45197808700377406*GTR[5]+0.45197808700377406*GTC[5]+0.45197808700377406*GCR[5]-0.45197808700377406*GCC[5]+0.3211688248608508*GTR[1]-0.3211688248608508*GTC[1]+0.3211688248608508*GCR[1]-0.3211688248608508*GCC[1]; 
  surft1_upper[16] = (0.25*dgdpv1_surf_CC_pv1[14])/dv1_pv1+(0.19364916731037085*dgdpv1_surf_CC_pv1[6])/dv1_pv1+(0.11180339887498951*dgdpv1_surf_CC_pv1[3])/dv1_pv1-0.4695652700276729*GTR[43]+1.1226030627813903*GTC[43]+0.4695652700276729*GCR[43]-1.1226030627813903*GCC[43]-0.4521030872910352*GTR[38]-0.7032714691193882*GTC[38]-0.4521030872910352*GCR[38]-0.7032714691193882*GCC[38]+0.3336660123724018*GTR[30]-0.7977048375260732*GTC[30]+0.3336660123724018*GCR[30]-0.7977048375260732*GCC[30]+0.3141929545311315*GTR[22]-0.3141929545311315*GTC[22]+0.3141929545311315*GCR[22]-0.3141929545311315*GCC[22]+0.6503668703432222*GTR[18]+1.0116817983116793*GTC[18]-0.6503668703432222*GCR[18]-1.0116817983116793*GCC[18]-0.46214090789498363*GTR[10]-0.7188858567255303*GTC[10]-0.46214090789498363*GCR[10]-0.7188858567255303*GCC[10]-0.45197808700377406*GTR[7]+0.45197808700377406*GTC[7]+0.45197808700377406*GCR[7]-0.45197808700377406*GCC[7]+0.3211688248608508*GTR[3]-0.3211688248608508*GTC[3]+0.3211688248608508*GCR[3]-0.3211688248608508*GCC[3]; 
  surft1_upper[17] = -((0.15*dgdpv1_surf_CC_pv1[17])/dv1_pv1)-(0.08660254037844385*dgdpv1_surf_CC_pv1[13])/dv1_pv1+0.10262022457558415*GTR[44]+0.10262022457558415*GTC[44]-0.10262022457558415*GCR[44]-0.10262022457558415*GCC[44]-0.07292038680986264*GTR[37]-0.07292038680986264*GTC[37]-0.07292038680986264*GCR[37]-0.07292038680986264*GCC[37]-0.08079247402229095*GTR[32]+0.08079247402229095*GTC[32]+0.08079247402229095*GCR[32]-0.08079247402229095*GCC[32]+0.05740991584648068*GTR[21]-0.05740991584648068*GTC[21]+0.05740991584648068*GCR[21]-0.05740991584648068*GCC[21]; 
  surft1_upper[18] = -((0.15*dgdpv1_surf_CC_pv1[19])/dv1_pv1)-(0.08660254037844385*dgdpv1_surf_CC_pv1[15])/dv1_pv1+0.10262022457558415*GTR[46]+0.10262022457558415*GTC[46]-0.10262022457558415*GCR[46]-0.10262022457558415*GCC[46]-0.07292038680986264*GTR[39]-0.07292038680986264*GTC[39]-0.07292038680986264*GCR[39]-0.07292038680986264*GCC[39]-0.08079247402229095*GTR[34]+0.08079247402229095*GTC[34]+0.08079247402229095*GCR[34]-0.08079247402229095*GCC[34]+0.05740991584648068*GTR[23]-0.05740991584648068*GTC[23]+0.05740991584648068*GCR[23]-0.05740991584648068*GCC[23]; 
  surft1_upper[19] = (0.25*dgdpv1_surf_CC_pv1[18])/dv1_pv1+(0.19364916731037082*dgdpv1_surf_CC_pv1[10])/dv1_pv1+(0.11180339887498948*dgdpv1_surf_CC_pv1[5])/dv1_pv1-0.4695652700276729*GTR[47]+1.1226030627813903*GTC[47]+0.4695652700276729*GCR[47]-1.1226030627813903*GCC[47]-0.4521030872910352*GTR[45]-0.7032714691193882*GTC[45]-0.4521030872910352*GCR[45]-0.7032714691193882*GCC[45]+0.3336660123724018*GTR[42]-0.7977048375260732*GTC[42]+0.3336660123724018*GCR[42]-0.7977048375260732*GCC[42]+0.3141929545311315*GTR[33]-0.3141929545311315*GTC[33]+0.3141929545311315*GCR[33]-0.3141929545311315*GCC[33]+0.6503668703432224*GTR[31]+1.0116817983116795*GTC[31]-0.6503668703432224*GCR[31]-1.0116817983116795*GCC[31]-0.46214090789498335*GTR[17]-0.71888585672553*GTC[17]-0.46214090789498335*GCR[17]-0.71888585672553*GCC[17]-0.45197808700377406*GTR[15]+0.45197808700377406*GTC[15]+0.45197808700377406*GCR[15]-0.45197808700377406*GCC[15]+0.3211688248608508*GTR[6]-0.3211688248608508*GTC[6]+0.3211688248608508*GCR[6]-0.3211688248608508*GCC[6]; 
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
  surft2_lower[0] = g_surf_CC_pv1[0]; 
  surft2_lower[1] = g_surf_CC_pv1[1]; 
  surft2_lower[2] = g_surf_CC_pv1[2]; 
  surft2_lower[3] = g_surf_CC_pv1[3]; 
  surft2_lower[4] = g_surf_CC_pv1[4]; 
  surft2_lower[5] = g_surf_CC_pv1[5]; 
  surft2_lower[6] = g_surf_CC_pv1[6]; 
  surft2_lower[7] = g_surf_CC_pv1[7]; 
  surft2_lower[8] = g_surf_CC_pv1[8]; 
  surft2_lower[9] = g_surf_CC_pv1[9]; 
  surft2_lower[10] = g_surf_CC_pv1[10]; 
  surft2_lower[11] = g_surf_CC_pv1[11]; 
  surft2_lower[12] = g_surf_CC_pv1[12]; 
  surft2_lower[13] = g_surf_CC_pv1[13]; 
  surft2_lower[14] = g_surf_CC_pv1[14]; 
  surft2_lower[15] = g_surf_CC_pv1[15]; 
  surft2_lower[16] = g_surf_CC_pv1[16]; 
  surft2_lower[17] = g_surf_CC_pv1[17]; 
  surft2_lower[18] = g_surf_CC_pv1[18]; 
  surft2_lower[19] = g_surf_CC_pv1[19]; 

  out[0] = 0.7071067811865475*surft1_upper[0]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[0]*dv1_sq*gamma_avg; 
  out[1] = 0.7071067811865475*surft1_upper[1]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[1]*dv1_sq*gamma_avg; 
  out[2] = -(1.224744871391589*surft2_upper[0]*dv1_sq*gamma_avg)+1.224744871391589*surft2_lower[0]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[0]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[0]*dv1_sq*gamma_avg; 
  out[3] = 0.7071067811865475*surft1_upper[2]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[2]*dv1_sq*gamma_avg; 
  out[4] = 0.7071067811865475*surft1_upper[3]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[3]*dv1_sq*gamma_avg; 
  out[5] = -(1.224744871391589*surft2_upper[1]*dv1_sq*gamma_avg)+1.224744871391589*surft2_lower[1]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[1]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[1]*dv1_sq*gamma_avg; 
  out[6] = 0.7071067811865475*surft1_upper[4]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[4]*dv1_sq*gamma_avg; 
  out[7] = -(1.224744871391589*surft2_upper[3]*dv1_sq*gamma_avg)+1.224744871391589*surft2_lower[3]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[2]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[2]*dv1_sq*gamma_avg; 
  out[8] = 0.7071067811865475*surft1_upper[5]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[5]*dv1_sq*gamma_avg; 
  out[9] = 1.224744871391589*surft1_upper[3]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[3]*dv1_sq*gamma_avg-2.1213203435596424*surft2_upper[0]*dv1_sq*gamma_avg-2.1213203435596424*surft2_lower[0]*dv1_sq*gamma_avg+3.0*GCC[0]*dv1_sq*gamma_avg; 
  out[10] = 0.7071067811865475*surft1_upper[6]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[6]*dv1_sq*gamma_avg; 
  out[11] = 0.7071067811865475*surft1_upper[7]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[7]*dv1_sq*gamma_avg; 
  out[12] = -(2.7386127875258306*surft2_upper[2]*dv1_sq*gamma_avg)+2.7386127875258306*surft2_lower[2]*dv1_sq*gamma_avg+1.5811388300841895*surft1_upper[0]*dv1_sq*gamma_avg-1.5811388300841895*surft1_lower[0]*dv1_sq*gamma_avg; 
  out[13] = 0.7071067811865475*surft1_upper[8]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[8]*dv1_sq*gamma_avg; 
  out[14] = 0.7071067811865475*surft1_upper[9]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[9]*dv1_sq*gamma_avg; 
  out[15] = -(1.224744871391589*surft2_upper[5]*dv1_sq*gamma_avg)+1.224744871391589*surft2_lower[5]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[4]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[4]*dv1_sq*gamma_avg; 
  out[16] = 1.224744871391589*surft1_upper[5]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[5]*dv1_sq*gamma_avg-2.1213203435596424*surft2_upper[1]*dv1_sq*gamma_avg-2.1213203435596424*surft2_lower[1]*dv1_sq*gamma_avg+3.0*GCC[1]*dv1_sq*gamma_avg; 
  out[17] = 0.7071067811865475*surft1_upper[10]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[10]*dv1_sq*gamma_avg; 
  out[18] = 1.224744871391589*surft1_upper[6]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[6]*dv1_sq*gamma_avg-2.1213203435596424*surft2_upper[3]*dv1_sq*gamma_avg-2.1213203435596424*surft2_lower[3]*dv1_sq*gamma_avg+3.0*GCC[3]*dv1_sq*gamma_avg; 
  out[19] = -(1.224744871391589*surft2_upper[7]*dv1_sq*gamma_avg)+1.224744871391589*surft2_lower[7]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[7]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[7]*dv1_sq*gamma_avg; 
  out[20] = -(2.7386127875258306*surft2_upper[4]*dv1_sq*gamma_avg)+2.7386127875258306*surft2_lower[4]*dv1_sq*gamma_avg+1.5811388300841898*surft1_upper[1]*dv1_sq*gamma_avg-1.5811388300841898*surft1_lower[1]*dv1_sq*gamma_avg; 
  out[21] = 0.7071067811865475*surft1_upper[11]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[11]*dv1_sq*gamma_avg; 
  out[22] = -(2.7386127875258306*surft2_upper[6]*dv1_sq*gamma_avg)+2.7386127875258306*surft2_lower[6]*dv1_sq*gamma_avg+1.5811388300841898*surft1_upper[2]*dv1_sq*gamma_avg-1.5811388300841898*surft1_lower[2]*dv1_sq*gamma_avg; 
  out[23] = 0.7071067811865475*surft1_upper[12]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[12]*dv1_sq*gamma_avg; 
  out[24] = -(1.224744871391589*surft2_upper[9]*dv1_sq*gamma_avg)+1.224744871391589*surft2_lower[9]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[8]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[8]*dv1_sq*gamma_avg; 
  out[25] = 0.7071067811865475*surft1_upper[13]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[13]*dv1_sq*gamma_avg; 
  out[26] = 1.5811388300841898*surft1_upper[3]*dv1_sq*gamma_avg-1.5811388300841898*surft1_lower[3]*dv1_sq*gamma_avg-4.743416490252569*surft2_upper[2]*dv1_sq*gamma_avg-4.743416490252569*surft2_lower[2]*dv1_sq*gamma_avg+6.7082039324993685*GCC[2]*dv1_sq*gamma_avg; 
  out[27] = 0.7071067811865475*surft1_upper[14]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[14]*dv1_sq*gamma_avg; 
  out[28] = 0.7071067811865475*surft1_upper[15]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[15]*dv1_sq*gamma_avg; 
  out[29] = 1.224744871391589*surft1_upper[9]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[9]*dv1_sq*gamma_avg+6.7082039324993685*GCC[4]*dv1_sq*gamma_avg-2.7386127875258306*surft2_upper[0]*dv1_sq*gamma_avg+2.7386127875258306*surft2_lower[0]*dv1_sq*gamma_avg; 
  out[30] = 0.7071067811865475*surft1_upper[16]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[16]*dv1_sq*gamma_avg; 
  out[31] = 1.224744871391589*surft1_upper[10]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[10]*dv1_sq*gamma_avg+3.0*GCC[6]*dv1_sq*gamma_avg-2.1213203435596424*surft2_upper[5]*dv1_sq*gamma_avg-2.1213203435596424*surft2_lower[5]*dv1_sq*gamma_avg; 
  out[32] = -(1.224744871391589*surft2_upper[13]*dv1_sq*gamma_avg)+1.224744871391589*surft2_lower[13]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[11]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[11]*dv1_sq*gamma_avg; 
  out[33] = -(2.7386127875258306*surft2_upper[10]*dv1_sq*gamma_avg)+2.7386127875258306*surft2_lower[10]*dv1_sq*gamma_avg+1.5811388300841895*surft1_upper[4]*dv1_sq*gamma_avg-1.5811388300841895*surft1_lower[4]*dv1_sq*gamma_avg; 
  out[34] = -(1.224744871391589*surft2_upper[15]*dv1_sq*gamma_avg)+1.224744871391589*surft2_lower[15]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[12]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[12]*dv1_sq*gamma_avg; 
  out[35] = 1.224744871391589*surft1_upper[13]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[13]*dv1_sq*gamma_avg+3.0*GCC[11]*dv1_sq*gamma_avg-2.1213203435596424*surft2_upper[7]*dv1_sq*gamma_avg-2.1213203435596424*surft2_lower[7]*dv1_sq*gamma_avg; 
  out[36] = 1.5811388300841895*surft1_upper[5]*dv1_sq*gamma_avg-1.5811388300841895*surft1_lower[5]*dv1_sq*gamma_avg+6.708203932499369*GCC[5]*dv1_sq*gamma_avg-4.743416490252569*surft2_upper[4]*dv1_sq*gamma_avg-4.743416490252569*surft2_lower[4]*dv1_sq*gamma_avg; 
  out[37] = 0.7071067811865475*surft1_upper[17]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[17]*dv1_sq*gamma_avg; 
  out[38] = 6.708203932499369*GCC[7]*dv1_sq*gamma_avg-4.743416490252569*surft2_upper[6]*dv1_sq*gamma_avg-4.743416490252569*surft2_lower[6]*dv1_sq*gamma_avg+1.5811388300841895*surft1_upper[6]*dv1_sq*gamma_avg-1.5811388300841895*surft1_lower[6]*dv1_sq*gamma_avg; 
  out[39] = 0.7071067811865475*surft1_upper[18]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[18]*dv1_sq*gamma_avg; 
  out[40] = 1.224744871391589*surft1_upper[14]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[14]*dv1_sq*gamma_avg+3.0*GCC[13]*dv1_sq*gamma_avg-2.1213203435596424*surft2_upper[9]*dv1_sq*gamma_avg-2.1213203435596424*surft2_lower[9]*dv1_sq*gamma_avg; 
  out[41] = 1.224744871391589*surft1_upper[15]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[15]*dv1_sq*gamma_avg+6.708203932499369*GCC[8]*dv1_sq*gamma_avg-2.7386127875258306*surft2_upper[1]*dv1_sq*gamma_avg+2.7386127875258306*surft2_lower[1]*dv1_sq*gamma_avg; 
  out[42] = 0.7071067811865475*surft1_upper[19]*dv1_sq*gamma_avg-0.7071067811865475*surft1_lower[19]*dv1_sq*gamma_avg; 
  out[43] = 1.224744871391589*surft1_upper[16]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[16]*dv1_sq*gamma_avg+6.708203932499369*GCC[10]*dv1_sq*gamma_avg-2.7386127875258306*surft2_upper[3]*dv1_sq*gamma_avg+2.7386127875258306*surft2_lower[3]*dv1_sq*gamma_avg; 
  out[44] = 3.0*GCC[21]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[17]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[17]*dv1_sq*gamma_avg-2.1213203435596424*surft2_upper[13]*dv1_sq*gamma_avg-2.1213203435596424*surft2_lower[13]*dv1_sq*gamma_avg; 
  out[45] = 6.7082039324993685*GCC[15]*dv1_sq*gamma_avg-4.743416490252569*surft2_upper[10]*dv1_sq*gamma_avg-4.743416490252569*surft2_lower[10]*dv1_sq*gamma_avg+1.5811388300841898*surft1_upper[10]*dv1_sq*gamma_avg-1.5811388300841898*surft1_lower[10]*dv1_sq*gamma_avg; 
  out[46] = 3.0*GCC[23]*dv1_sq*gamma_avg+1.224744871391589*surft1_upper[18]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[18]*dv1_sq*gamma_avg-2.1213203435596424*surft2_upper[15]*dv1_sq*gamma_avg-2.1213203435596424*surft2_lower[15]*dv1_sq*gamma_avg; 
  out[47] = 1.224744871391589*surft1_upper[19]*dv1_sq*gamma_avg+1.224744871391589*surft1_lower[19]*dv1_sq*gamma_avg+6.7082039324993685*GCC[17]*dv1_sq*gamma_avg-2.7386127875258306*surft2_upper[5]*dv1_sq*gamma_avg+2.7386127875258306*surft2_lower[5]*dv1_sq*gamma_avg; 
} 
