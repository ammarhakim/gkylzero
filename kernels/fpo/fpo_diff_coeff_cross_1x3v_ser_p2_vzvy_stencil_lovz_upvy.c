#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvy_ser_p2_lovz_upvy(const double *dxv, double gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // diff_coeff: Output array for diffusion tensor. 

  double dv1_pv1 = 2.0/dxv[3]; 
  double dv1_pv2 = 2.0/dxv[2]; 
  double dv1_sq = 4.0/dxv[3]/dxv[2]; 
 
  const double* GBC = fpo_g_stencil[0]; 
  const double* GBR = fpo_g_stencil[1]; 
  const double* GCC = fpo_g_stencil[2]; 
  const double* GCR = fpo_g_stencil[3]; 

  const double* g_surf_CC = fpo_g_surf_stencil[2]; 
  const double* g_surf_CC_pv2 = &g_surf_CC[20]; 
  const double* g_surf_CR = fpo_g_surf_stencil[3]; 
  const double* g_surf_CR_pv2 = &g_surf_CR[20]; 
  
  const double* g_surf_CC_pv1 = &g_surf_CC[40]; 
  const double* dgdpv1_surf_CC_pv2 = &fpo_dgdv_surf[100]; 
  const double* dgdpv2_surf_CC_pv1 = &fpo_dgdv_surf[140]; 
  const double* dgdpv1_surf_CC_pv1 = &fpo_dgdv_surf[160]; 
  
  double surft1_upper[20], surft1_lower[20]; 
  double surft2_upper[20], surft2_lower[20]; 
  
  double *diff_coeff_vxvy = &diff_coeff[48]; 
  double *diff_coeff_vxvz = &diff_coeff[96]; 
  double *diff_coeff_vyvx = &diff_coeff[144]; 
  double *diff_coeff_vyvz = &diff_coeff[240]; 
  double *diff_coeff_vzvx = &diff_coeff[288]; 
  double *diff_coeff_vzvy = &diff_coeff[336]; 
  
  double *out = diff_coeff_vzvy; 
  
  surft1_upper[0] = dgdpv1_surf_CC_pv2[0]/dv1_pv1; 
  surft1_upper[1] = dgdpv1_surf_CC_pv2[1]/dv1_pv1; 
  surft1_upper[2] = dgdpv1_surf_CC_pv2[2]/dv1_pv1; 
  surft1_upper[3] = dgdpv1_surf_CC_pv2[3]/dv1_pv1; 
  surft1_upper[4] = dgdpv1_surf_CC_pv2[4]/dv1_pv1; 
  surft1_upper[5] = dgdpv1_surf_CC_pv2[5]/dv1_pv1; 
  surft1_upper[6] = dgdpv1_surf_CC_pv2[6]/dv1_pv1; 
  surft1_upper[7] = dgdpv1_surf_CC_pv2[7]/dv1_pv1; 
  surft1_upper[8] = dgdpv1_surf_CC_pv2[8]/dv1_pv1; 
  surft1_upper[9] = dgdpv1_surf_CC_pv2[9]/dv1_pv1; 
  surft1_upper[10] = dgdpv1_surf_CC_pv2[10]/dv1_pv1; 
  surft1_upper[11] = dgdpv1_surf_CC_pv2[11]/dv1_pv1; 
  surft1_upper[12] = dgdpv1_surf_CC_pv2[12]/dv1_pv1; 
  surft1_upper[13] = dgdpv1_surf_CC_pv2[13]/dv1_pv1; 
  surft1_upper[14] = dgdpv1_surf_CC_pv2[14]/dv1_pv1; 
  surft1_upper[15] = dgdpv1_surf_CC_pv2[15]/dv1_pv1; 
  surft1_upper[16] = dgdpv1_surf_CC_pv2[16]/dv1_pv1; 
  surft1_upper[17] = dgdpv1_surf_CC_pv2[17]/dv1_pv1; 
  surft1_upper[18] = dgdpv1_surf_CC_pv2[18]/dv1_pv1; 
  surft1_upper[19] = dgdpv1_surf_CC_pv2[19]/dv1_pv1; 
  surft1_lower[0] = (0.09316949906249124*dgdpv1_surf_CC_pv1[9])/dv1_pv1-(0.07216878364870323*dgdpv1_surf_CC_pv1[3])/dv1_pv1+(0.041666666666666664*dgdpv1_surf_CC_pv1[0])/dv1_pv1-0.26249496591373034*GCR[30]+0.6275541900321632*GCC[30]+0.26249496591373034*GBR[30]-0.6275541900321632*GBC[30]-0.24493289823330752*GCR[27]+0.23869256961589846*GCC[27]-0.24493289823330752*GBR[27]+0.23869256961589846*GBC[27]+0.186524971361494*GCR[14]-0.44593056067218106*GCC[14]+0.186524971361494*GBR[14]-0.44593056067218106*GBC[14]+0.16753277642981868*GCR[13]-0.16753277642981868*GCC[13]+0.16753277642981868*GBR[13]-0.16753277642981868*GBC[13]+0.35234495615276884*GCR[10]-0.3433680145947365*GCC[10]-0.35234495615276884*GBR[10]+0.3433680145947365*GBC[10]-0.2503710218860407*GCR[4]+0.2439921423475428*GCC[4]-0.2503710218860407*GBR[4]+0.2439921423475428*GBC[4]-0.2410020425638719*GCR[3]+0.2410020425638719*GCC[3]+0.2410020425638719*GBR[3]-0.2410020425638719*GBC[3]+0.1712524235686168*GCR[0]-0.1712524235686168*GCC[0]+0.1712524235686168*GBR[0]-0.1712524235686168*GBC[0]; 
  surft1_lower[1] = (0.09316949906249125*dgdpv1_surf_CC_pv1[15])/dv1_pv1-(0.07216878364870323*dgdpv1_surf_CC_pv1[5])/dv1_pv1+(0.041666666666666664*dgdpv1_surf_CC_pv1[1])/dv1_pv1-0.2624949659137303*GCR[42]+0.6275541900321632*GCC[42]+0.2624949659137303*GBR[42]-0.6275541900321632*GBC[42]-0.24493289823330755*GCR[39]+0.23869256961589838*GCC[39]-0.24493289823330755*GBR[39]+0.23869256961589838*GBC[39]+0.18652497136149407*GCR[28]-0.44593056067218123*GCC[28]+0.18652497136149407*GBR[28]-0.44593056067218123*GBC[28]+0.1675327764298187*GCR[23]-0.1675327764298187*GCC[23]+0.1675327764298187*GBR[23]-0.1675327764298187*GBC[23]+0.35234495615276884*GCR[17]-0.3433680145947365*GCC[17]-0.35234495615276884*GBR[17]+0.3433680145947365*GBC[17]-0.2503710218860407*GCR[8]+0.2439921423475428*GCC[8]-0.2503710218860407*GBR[8]+0.2439921423475428*GBC[8]-0.2410020425638719*GCR[6]+0.2410020425638719*GCC[6]+0.2410020425638719*GBR[6]-0.2410020425638719*GBC[6]+0.1712524235686168*GCR[1]-0.1712524235686168*GCC[1]+0.1712524235686168*GBR[1]-0.1712524235686168*GBC[1]; 
  surft1_lower[2] = (0.09316949906249125*dgdpv1_surf_CC_pv1[16])/dv1_pv1-(0.07216878364870323*dgdpv1_surf_CC_pv1[6])/dv1_pv1+(0.041666666666666664*dgdpv1_surf_CC_pv1[2])/dv1_pv1-0.2624949659137303*GCR[43]+0.6275541900321632*GCC[43]+0.2624949659137303*GBR[43]-0.6275541900321632*GBC[43]-0.24493289823330755*GCR[40]+0.23869256961589838*GCC[40]-0.24493289823330755*GBR[40]+0.23869256961589838*GBC[40]+0.18652497136149407*GCR[29]-0.44593056067218123*GCC[29]+0.18652497136149407*GBR[29]-0.44593056067218123*GBC[29]+0.1675327764298187*GCR[24]-0.1675327764298187*GCC[24]+0.1675327764298187*GBR[24]-0.1675327764298187*GBC[24]+0.35234495615276884*GCR[18]-0.3433680145947365*GCC[18]-0.35234495615276884*GBR[18]+0.3433680145947365*GBC[18]-0.2503710218860407*GCR[9]+0.2439921423475428*GCC[9]-0.2503710218860407*GBR[9]+0.2439921423475428*GBC[9]-0.2410020425638719*GCR[7]+0.2410020425638719*GCC[7]+0.2410020425638719*GBR[7]-0.2410020425638719*GBC[7]+0.1712524235686168*GCR[2]-0.1712524235686168*GCC[2]+0.1712524235686168*GBR[2]-0.1712524235686168*GBC[2]; 
  surft1_lower[3] = -((0.16137430609197573*dgdpv1_surf_CC_pv1[9])/dv1_pv1)+(0.125*dgdpv1_surf_CC_pv1[3])/dv1_pv1-(0.07216878364870323*dgdpv1_surf_CC_pv1[0])/dv1_pv1+0.03312053506377521*GCR[30]-1.5084898242683071*GCC[30]-0.03312053506377521*GBR[30]+1.5084898242683071*GBC[30]+0.0027021415553196565*GCR[27]+0.008106424665958968*GCC[27]+0.0027021415553196565*GBR[27]+0.008106424665958968*GBC[27]-0.023534953642800078*GCR[14]+1.0719101613675308*GCC[14]-0.023534953642800078*GBR[14]+1.0719101613675308*GBC[14]+0.009360492926113665*GCR[13]-0.009360492926113665*GCC[13]+0.009360492926113665*GBR[13]-0.009360492926113665*GBC[13]-0.003887129718772127*GCR[10]-0.01166138915631638*GCC[10]+0.003887129718772127*GBR[10]+0.01166138915631638*GBC[10]+0.002762135864009948*GCR[4]+0.008286407592029844*GCC[4]+0.002762135864009948*GBR[4]+0.008286407592029844*GBC[4]-0.013465412337048493*GCR[3]+0.013465412337048493*GCC[3]+0.013465412337048493*GBR[3]-0.013465412337048493*GBC[3]+0.009568319307746778*GCR[0]-0.009568319307746778*GCC[0]+0.009568319307746778*GBR[0]-0.009568319307746778*GBC[0]; 
  surft1_lower[4] = (0.09316949906249124*dgdpv1_surf_CC_pv1[19])/dv1_pv1-(0.07216878364870323*dgdpv1_surf_CC_pv1[10])/dv1_pv1+(0.041666666666666664*dgdpv1_surf_CC_pv1[4])/dv1_pv1-0.26249496591373034*GCR[47]+0.6275541900321632*GCC[47]+0.26249496591373034*GBR[47]-0.6275541900321632*GBC[47]-0.24493289823330752*GCR[46]+0.23869256961589846*GCC[46]-0.24493289823330752*GBR[46]+0.23869256961589846*GBC[46]+0.186524971361494*GCR[41]-0.44593056067218106*GCC[41]+0.186524971361494*GBR[41]-0.44593056067218106*GBC[41]+0.16753277642981868*GCR[34]-0.16753277642981868*GCC[34]+0.16753277642981868*GBR[34]-0.16753277642981868*GBC[34]+0.35234495615276884*GCR[31]-0.3433680145947365*GCC[31]-0.35234495615276884*GBR[31]+0.3433680145947365*GBC[31]-0.2503710218860407*GCR[16]+0.2439921423475428*GCC[16]-0.2503710218860407*GBR[16]+0.2439921423475428*GBC[16]-0.2410020425638719*GCR[15]+0.2410020425638719*GCC[15]+0.2410020425638719*GBR[15]-0.2410020425638719*GBC[15]+0.1712524235686168*GCR[5]-0.1712524235686168*GCC[5]+0.1712524235686168*GBR[5]-0.1712524235686168*GBC[5]; 
  surft1_lower[5] = -((0.1613743060919757*dgdpv1_surf_CC_pv1[15])/dv1_pv1)+(0.125*dgdpv1_surf_CC_pv1[5])/dv1_pv1-(0.07216878364870323*dgdpv1_surf_CC_pv1[1])/dv1_pv1+0.03312053506377521*GCR[42]-1.5084898242683071*GCC[42]-0.03312053506377521*GBR[42]+1.5084898242683071*GBC[42]+0.0027021415553196556*GCR[39]+0.008106424665958968*GCC[39]+0.0027021415553196556*GBR[39]+0.008106424665958968*GBC[39]-0.02353495364280008*GCR[28]+1.071910161367531*GCC[28]-0.02353495364280008*GBR[28]+1.071910161367531*GBC[28]+0.009360492926113666*GCR[23]-0.009360492926113666*GCC[23]+0.009360492926113666*GBR[23]-0.009360492926113666*GBC[23]-0.003887129718772127*GCR[17]-0.01166138915631638*GCC[17]+0.003887129718772127*GBR[17]+0.01166138915631638*GBC[17]+0.002762135864009948*GCR[8]+0.008286407592029844*GCC[8]+0.002762135864009948*GBR[8]+0.008286407592029844*GBC[8]-0.013465412337048493*GCR[6]+0.013465412337048493*GCC[6]+0.013465412337048493*GBR[6]-0.013465412337048493*GBC[6]+0.009568319307746778*GCR[1]-0.009568319307746778*GCC[1]+0.009568319307746778*GBR[1]-0.009568319307746778*GBC[1]; 
  surft1_lower[6] = -((0.1613743060919757*dgdpv1_surf_CC_pv1[16])/dv1_pv1)+(0.125*dgdpv1_surf_CC_pv1[6])/dv1_pv1-(0.07216878364870323*dgdpv1_surf_CC_pv1[2])/dv1_pv1+0.03312053506377521*GCR[43]-1.5084898242683071*GCC[43]-0.03312053506377521*GBR[43]+1.5084898242683071*GBC[43]+0.0027021415553196556*GCR[40]+0.008106424665958968*GCC[40]+0.0027021415553196556*GBR[40]+0.008106424665958968*GBC[40]-0.02353495364280008*GCR[29]+1.071910161367531*GCC[29]-0.02353495364280008*GBR[29]+1.071910161367531*GBC[29]+0.009360492926113666*GCR[24]-0.009360492926113666*GCC[24]+0.009360492926113666*GBR[24]-0.009360492926113666*GBC[24]-0.003887129718772127*GCR[18]-0.01166138915631638*GCC[18]+0.003887129718772127*GBR[18]+0.01166138915631638*GBC[18]+0.002762135864009948*GCR[9]+0.008286407592029844*GCC[9]+0.002762135864009948*GBR[9]+0.008286407592029844*GBC[9]-0.013465412337048493*GCR[7]+0.013465412337048493*GCC[7]+0.013465412337048493*GBR[7]-0.013465412337048493*GBC[7]+0.009568319307746778*GCR[2]-0.009568319307746778*GCC[2]+0.009568319307746778*GBR[2]-0.009568319307746778*GBC[2]; 
  surft1_lower[7] = -((0.07216878364870322*dgdpv1_surf_CC_pv1[13])/dv1_pv1)+(0.041666666666666664*dgdpv1_surf_CC_pv1[7])/dv1_pv1+0.35234495615276884*GCR[37]-0.3433680145947365*GCC[37]-0.35234495615276884*GBR[37]+0.3433680145947365*GBC[37]-0.2503710218860407*GCR[25]+0.24399214234754282*GCC[25]-0.2503710218860407*GBR[25]+0.24399214234754282*GBC[25]-0.24100204256387192*GCR[21]+0.24100204256387192*GCC[21]+0.24100204256387192*GBR[21]-0.24100204256387192*GBC[21]+0.1712524235686168*GCR[11]-0.1712524235686168*GCC[11]+0.1712524235686168*GBR[11]-0.1712524235686168*GBC[11]; 
  surft1_lower[8] = -((0.07216878364870322*dgdpv1_surf_CC_pv1[14])/dv1_pv1)+(0.041666666666666664*dgdpv1_surf_CC_pv1[8])/dv1_pv1+0.35234495615276884*GCR[38]-0.3433680145947365*GCC[38]-0.35234495615276884*GBR[38]+0.3433680145947365*GBC[38]-0.2503710218860407*GCR[26]+0.24399214234754282*GCC[26]-0.2503710218860407*GBR[26]+0.24399214234754282*GBC[26]-0.24100204256387192*GCR[22]+0.24100204256387192*GCC[22]+0.24100204256387192*GBR[22]-0.24100204256387192*GBC[22]+0.1712524235686168*GCR[12]-0.1712524235686168*GCC[12]+0.1712524235686168*GBR[12]-0.1712524235686168*GBC[12]; 
  surft1_lower[9] = (0.20833333333333334*dgdpv1_surf_CC_pv1[9])/dv1_pv1-(0.16137430609197573*dgdpv1_surf_CC_pv1[3])/dv1_pv1+(0.09316949906249124*dgdpv1_surf_CC_pv1[0])/dv1_pv1-0.5869565875345911*GCR[30]+1.4032538284767377*GCC[30]+0.5869565875345911*GBR[30]-1.4032538284767377*GBC[30]-0.5476866103757138*GCR[27]-0.8058318916992987*GCC[27]-0.5476866103757138*GBR[27]-0.8058318916992987*GBC[27]+0.4170825154655021*GCR[14]-0.9971310469075912*GCC[14]+0.4170825154655021*GBR[14]-0.9971310469075912*GBC[14]+0.37461467655634906*GCR[13]-0.37461467655634906*GCC[13]+0.37461467655634906*GBR[13]-0.37461467655634906*GBC[13]+0.787867273486774*GCR[10]+1.1592187272321324*GCC[10]-0.787867273486774*GBR[10]-1.1592187272321324*GBC[10]-0.5598466245332747*GCR[4]-0.8237233774980027*GCC[4]-0.5598466245332747*GBR[4]-0.8237233774980027*GBC[4]-0.5388969498891153*GCR[3]+0.5388969498891153*GCC[3]+0.5388969498891153*GBR[3]-0.5388969498891153*GBC[3]+0.3829320604110143*GCR[0]-0.3829320604110143*GCC[0]+0.3829320604110143*GBR[0]-0.3829320604110143*GBC[0]; 
  surft1_lower[10] = -((0.16137430609197573*dgdpv1_surf_CC_pv1[19])/dv1_pv1)+(0.125*dgdpv1_surf_CC_pv1[10])/dv1_pv1-(0.07216878364870323*dgdpv1_surf_CC_pv1[4])/dv1_pv1+0.03312053506377521*GCR[47]-1.5084898242683071*GCC[47]-0.03312053506377521*GBR[47]+1.5084898242683071*GBC[47]+0.0027021415553196565*GCR[46]+0.008106424665958968*GCC[46]+0.0027021415553196565*GBR[46]+0.008106424665958968*GBC[46]-0.023534953642800078*GCR[41]+1.0719101613675308*GCC[41]-0.023534953642800078*GBR[41]+1.0719101613675308*GBC[41]+0.009360492926113665*GCR[34]-0.009360492926113665*GCC[34]+0.009360492926113665*GBR[34]-0.009360492926113665*GBC[34]-0.003887129718772127*GCR[31]-0.01166138915631638*GCC[31]+0.003887129718772127*GBR[31]+0.01166138915631638*GBC[31]+0.002762135864009948*GCR[16]+0.008286407592029844*GCC[16]+0.002762135864009948*GBR[16]+0.008286407592029844*GBC[16]-0.013465412337048493*GCR[15]+0.013465412337048493*GCC[15]+0.013465412337048493*GBR[15]-0.013465412337048493*GBC[15]+0.009568319307746778*GCR[5]-0.009568319307746778*GCC[5]+0.009568319307746778*GBR[5]-0.009568319307746778*GBC[5]; 
  surft1_lower[11] = -((0.07216878364870322*dgdpv1_surf_CC_pv1[17])/dv1_pv1)+(0.041666666666666664*dgdpv1_surf_CC_pv1[11])/dv1_pv1+0.35234495615276884*GCR[44]-0.3433680145947365*GCC[44]-0.35234495615276884*GBR[44]+0.3433680145947365*GBC[44]-0.2503710218860407*GCR[35]+0.24399214234754282*GCC[35]-0.2503710218860407*GBR[35]+0.24399214234754282*GBC[35]-0.24100204256387192*GCR[32]+0.24100204256387192*GCC[32]+0.24100204256387192*GBR[32]-0.24100204256387192*GBC[32]+0.1712524235686168*GCR[19]-0.1712524235686168*GCC[19]+0.1712524235686168*GBR[19]-0.1712524235686168*GBC[19]; 
  surft1_lower[12] = -((0.07216878364870322*dgdpv1_surf_CC_pv1[18])/dv1_pv1)+(0.041666666666666664*dgdpv1_surf_CC_pv1[12])/dv1_pv1+0.35234495615276884*GCR[45]-0.3433680145947365*GCC[45]-0.35234495615276884*GBR[45]+0.3433680145947365*GBC[45]-0.2503710218860407*GCR[36]+0.24399214234754282*GCC[36]-0.2503710218860407*GBR[36]+0.24399214234754282*GBC[36]-0.24100204256387192*GCR[33]+0.24100204256387192*GCC[33]+0.24100204256387192*GBR[33]-0.24100204256387192*GBC[33]+0.1712524235686168*GCR[20]-0.1712524235686168*GCC[20]+0.1712524235686168*GBR[20]-0.1712524235686168*GBC[20]; 
  surft1_lower[13] = (0.125*dgdpv1_surf_CC_pv1[13])/dv1_pv1-(0.07216878364870322*dgdpv1_surf_CC_pv1[7])/dv1_pv1-0.0038871297187721273*GCR[37]-0.01166138915631638*GCC[37]+0.0038871297187721273*GBR[37]+0.01166138915631638*GBC[37]+0.002762135864009948*GCR[25]+0.008286407592029844*GCC[25]+0.002762135864009948*GBR[25]+0.008286407592029844*GBC[25]-0.013465412337048493*GCR[21]+0.013465412337048493*GCC[21]+0.013465412337048493*GBR[21]-0.013465412337048493*GBC[21]+0.00956831930774678*GCR[11]-0.00956831930774678*GCC[11]+0.00956831930774678*GBR[11]-0.00956831930774678*GBC[11]; 
  surft1_lower[14] = (0.125*dgdpv1_surf_CC_pv1[14])/dv1_pv1-(0.07216878364870322*dgdpv1_surf_CC_pv1[8])/dv1_pv1-0.0038871297187721273*GCR[38]-0.01166138915631638*GCC[38]+0.0038871297187721273*GBR[38]+0.01166138915631638*GBC[38]+0.002762135864009948*GCR[26]+0.008286407592029844*GCC[26]+0.002762135864009948*GBR[26]+0.008286407592029844*GBC[26]-0.013465412337048493*GCR[22]+0.013465412337048493*GCC[22]+0.013465412337048493*GBR[22]-0.013465412337048493*GBC[22]+0.00956831930774678*GCR[12]-0.00956831930774678*GCC[12]+0.00956831930774678*GBR[12]-0.00956831930774678*GBC[12]; 
  surft1_lower[15] = (0.20833333333333334*dgdpv1_surf_CC_pv1[15])/dv1_pv1-(0.1613743060919757*dgdpv1_surf_CC_pv1[5])/dv1_pv1+(0.09316949906249125*dgdpv1_surf_CC_pv1[1])/dv1_pv1-0.5869565875345911*GCR[42]+1.4032538284767377*GCC[42]+0.5869565875345911*GBR[42]-1.4032538284767377*GBC[42]-0.5476866103757138*GCR[39]-0.8058318916992987*GCC[39]-0.5476866103757138*GBR[39]-0.8058318916992987*GBC[39]+0.4170825154655021*GCR[28]-0.9971310469075912*GCC[28]+0.4170825154655021*GBR[28]-0.9971310469075912*GBC[28]+0.37461467655634906*GCR[23]-0.37461467655634906*GCC[23]+0.37461467655634906*GBR[23]-0.37461467655634906*GBC[23]+0.787867273486774*GCR[17]+1.1592187272321322*GCC[17]-0.787867273486774*GBR[17]-1.1592187272321322*GBC[17]-0.5598466245332747*GCR[8]-0.8237233774980027*GCC[8]-0.5598466245332747*GBR[8]-0.8237233774980027*GBC[8]-0.5388969498891153*GCR[6]+0.5388969498891153*GCC[6]+0.5388969498891153*GBR[6]-0.5388969498891153*GBC[6]+0.38293206041101435*GCR[1]-0.38293206041101435*GCC[1]+0.38293206041101435*GBR[1]-0.38293206041101435*GBC[1]; 
  surft1_lower[16] = (0.20833333333333334*dgdpv1_surf_CC_pv1[16])/dv1_pv1-(0.1613743060919757*dgdpv1_surf_CC_pv1[6])/dv1_pv1+(0.09316949906249125*dgdpv1_surf_CC_pv1[2])/dv1_pv1-0.5869565875345911*GCR[43]+1.4032538284767377*GCC[43]+0.5869565875345911*GBR[43]-1.4032538284767377*GBC[43]-0.5476866103757138*GCR[40]-0.8058318916992987*GCC[40]-0.5476866103757138*GBR[40]-0.8058318916992987*GBC[40]+0.4170825154655021*GCR[29]-0.9971310469075912*GCC[29]+0.4170825154655021*GBR[29]-0.9971310469075912*GBC[29]+0.37461467655634906*GCR[24]-0.37461467655634906*GCC[24]+0.37461467655634906*GBR[24]-0.37461467655634906*GBC[24]+0.787867273486774*GCR[18]+1.1592187272321322*GCC[18]-0.787867273486774*GBR[18]-1.1592187272321322*GBC[18]-0.5598466245332747*GCR[9]-0.8237233774980027*GCC[9]-0.5598466245332747*GBR[9]-0.8237233774980027*GBC[9]-0.5388969498891153*GCR[7]+0.5388969498891153*GCC[7]+0.5388969498891153*GBR[7]-0.5388969498891153*GBC[7]+0.38293206041101435*GCR[2]-0.38293206041101435*GCC[2]+0.38293206041101435*GBR[2]-0.38293206041101435*GBC[2]; 
  surft1_lower[17] = (0.125*dgdpv1_surf_CC_pv1[17])/dv1_pv1-(0.07216878364870322*dgdpv1_surf_CC_pv1[11])/dv1_pv1-0.0038871297187721273*GCR[44]-0.01166138915631638*GCC[44]+0.0038871297187721273*GBR[44]+0.01166138915631638*GBC[44]+0.002762135864009948*GCR[35]+0.008286407592029844*GCC[35]+0.002762135864009948*GBR[35]+0.008286407592029844*GBC[35]-0.013465412337048493*GCR[32]+0.013465412337048493*GCC[32]+0.013465412337048493*GBR[32]-0.013465412337048493*GBC[32]+0.00956831930774678*GCR[19]-0.00956831930774678*GCC[19]+0.00956831930774678*GBR[19]-0.00956831930774678*GBC[19]; 
  surft1_lower[18] = (0.125*dgdpv1_surf_CC_pv1[18])/dv1_pv1-(0.07216878364870322*dgdpv1_surf_CC_pv1[12])/dv1_pv1-0.0038871297187721273*GCR[45]-0.01166138915631638*GCC[45]+0.0038871297187721273*GBR[45]+0.01166138915631638*GBC[45]+0.002762135864009948*GCR[36]+0.008286407592029844*GCC[36]+0.002762135864009948*GBR[36]+0.008286407592029844*GBC[36]-0.013465412337048493*GCR[33]+0.013465412337048493*GCC[33]+0.013465412337048493*GBR[33]-0.013465412337048493*GBC[33]+0.00956831930774678*GCR[20]-0.00956831930774678*GCC[20]+0.00956831930774678*GBR[20]-0.00956831930774678*GBC[20]; 
  surft1_lower[19] = (0.20833333333333334*dgdpv1_surf_CC_pv1[19])/dv1_pv1-(0.16137430609197573*dgdpv1_surf_CC_pv1[10])/dv1_pv1+(0.09316949906249124*dgdpv1_surf_CC_pv1[4])/dv1_pv1-0.5869565875345911*GCR[47]+1.4032538284767377*GCC[47]+0.5869565875345911*GBR[47]-1.4032538284767377*GBC[47]-0.5476866103757138*GCR[46]-0.8058318916992987*GCC[46]-0.5476866103757138*GBR[46]-0.8058318916992987*GBC[46]+0.4170825154655021*GCR[41]-0.9971310469075912*GCC[41]+0.4170825154655021*GBR[41]-0.9971310469075912*GBC[41]+0.37461467655634906*GCR[34]-0.37461467655634906*GCC[34]+0.37461467655634906*GBR[34]-0.37461467655634906*GBC[34]+0.787867273486774*GCR[31]+1.1592187272321324*GCC[31]-0.787867273486774*GBR[31]-1.1592187272321324*GBC[31]-0.5598466245332747*GCR[16]-0.8237233774980027*GCC[16]-0.5598466245332747*GBR[16]-0.8237233774980027*GBC[16]-0.5388969498891153*GCR[15]+0.5388969498891153*GCC[15]+0.5388969498891153*GBR[15]-0.5388969498891153*GBC[15]+0.3829320604110143*GCR[5]-0.3829320604110143*GCC[5]+0.3829320604110143*GBR[5]-0.3829320604110143*GBC[5]; 

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

  out[0] = 0.7071067811865475*surft1_upper[0]*dv1_sq*gamma-0.7071067811865475*surft1_lower[0]*dv1_sq*gamma; 
  out[1] = 0.7071067811865475*surft1_upper[1]*dv1_sq*gamma-0.7071067811865475*surft1_lower[1]*dv1_sq*gamma; 
  out[2] = 0.7071067811865475*surft1_upper[2]*dv1_sq*gamma-0.7071067811865475*surft1_lower[2]*dv1_sq*gamma; 
  out[3] = -(1.224744871391589*surft2_upper[0]*dv1_sq*gamma)+1.224744871391589*surft2_lower[0]*dv1_sq*gamma+1.224744871391589*surft1_upper[0]*dv1_sq*gamma+1.224744871391589*surft1_lower[0]*dv1_sq*gamma; 
  out[4] = 0.7071067811865475*surft1_upper[3]*dv1_sq*gamma-0.7071067811865475*surft1_lower[3]*dv1_sq*gamma; 
  out[5] = 0.7071067811865475*surft1_upper[4]*dv1_sq*gamma-0.7071067811865475*surft1_lower[4]*dv1_sq*gamma; 
  out[6] = -(1.224744871391589*surft2_upper[1]*dv1_sq*gamma)+1.224744871391589*surft2_lower[1]*dv1_sq*gamma+1.224744871391589*surft1_upper[1]*dv1_sq*gamma+1.224744871391589*surft1_lower[1]*dv1_sq*gamma; 
  out[7] = -(1.224744871391589*surft2_upper[2]*dv1_sq*gamma)+1.224744871391589*surft2_lower[2]*dv1_sq*gamma+1.224744871391589*surft1_upper[2]*dv1_sq*gamma+1.224744871391589*surft1_lower[2]*dv1_sq*gamma; 
  out[8] = 0.7071067811865475*surft1_upper[5]*dv1_sq*gamma-0.7071067811865475*surft1_lower[5]*dv1_sq*gamma; 
  out[9] = 0.7071067811865475*surft1_upper[6]*dv1_sq*gamma-0.7071067811865475*surft1_lower[6]*dv1_sq*gamma; 
  out[10] = 1.224744871391589*surft1_upper[3]*dv1_sq*gamma+1.224744871391589*surft1_lower[3]*dv1_sq*gamma-2.1213203435596424*surft2_upper[0]*dv1_sq*gamma-2.1213203435596424*surft2_lower[0]*dv1_sq*gamma+3.0*GCC[0]*dv1_sq*gamma; 
  out[11] = 0.7071067811865475*surft1_upper[7]*dv1_sq*gamma-0.7071067811865475*surft1_lower[7]*dv1_sq*gamma; 
  out[12] = 0.7071067811865475*surft1_upper[8]*dv1_sq*gamma-0.7071067811865475*surft1_lower[8]*dv1_sq*gamma; 
  out[13] = -(2.7386127875258306*surft2_upper[3]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[3]*dv1_sq*gamma+1.5811388300841895*surft1_upper[0]*dv1_sq*gamma-1.5811388300841895*surft1_lower[0]*dv1_sq*gamma; 
  out[14] = 0.7071067811865475*surft1_upper[9]*dv1_sq*gamma-0.7071067811865475*surft1_lower[9]*dv1_sq*gamma; 
  out[15] = -(1.224744871391589*surft2_upper[4]*dv1_sq*gamma)+1.224744871391589*surft2_lower[4]*dv1_sq*gamma+1.224744871391589*surft1_upper[4]*dv1_sq*gamma+1.224744871391589*surft1_lower[4]*dv1_sq*gamma; 
  out[16] = 0.7071067811865475*surft1_upper[10]*dv1_sq*gamma-0.7071067811865475*surft1_lower[10]*dv1_sq*gamma; 
  out[17] = 1.224744871391589*surft1_upper[5]*dv1_sq*gamma+1.224744871391589*surft1_lower[5]*dv1_sq*gamma-2.1213203435596424*surft2_upper[1]*dv1_sq*gamma-2.1213203435596424*surft2_lower[1]*dv1_sq*gamma+3.0*GCC[1]*dv1_sq*gamma; 
  out[18] = 1.224744871391589*surft1_upper[6]*dv1_sq*gamma+1.224744871391589*surft1_lower[6]*dv1_sq*gamma-2.1213203435596424*surft2_upper[2]*dv1_sq*gamma-2.1213203435596424*surft2_lower[2]*dv1_sq*gamma+3.0*GCC[2]*dv1_sq*gamma; 
  out[19] = 0.7071067811865475*surft1_upper[11]*dv1_sq*gamma-0.7071067811865475*surft1_lower[11]*dv1_sq*gamma; 
  out[20] = 0.7071067811865475*surft1_upper[12]*dv1_sq*gamma-0.7071067811865475*surft1_lower[12]*dv1_sq*gamma; 
  out[21] = -(1.224744871391589*surft2_upper[7]*dv1_sq*gamma)+1.224744871391589*surft2_lower[7]*dv1_sq*gamma+1.224744871391589*surft1_upper[7]*dv1_sq*gamma+1.224744871391589*surft1_lower[7]*dv1_sq*gamma; 
  out[22] = -(1.224744871391589*surft2_upper[8]*dv1_sq*gamma)+1.224744871391589*surft2_lower[8]*dv1_sq*gamma+1.224744871391589*surft1_upper[8]*dv1_sq*gamma+1.224744871391589*surft1_lower[8]*dv1_sq*gamma; 
  out[23] = -(2.7386127875258306*surft2_upper[5]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[5]*dv1_sq*gamma+1.5811388300841898*surft1_upper[1]*dv1_sq*gamma-1.5811388300841898*surft1_lower[1]*dv1_sq*gamma; 
  out[24] = -(2.7386127875258306*surft2_upper[6]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[6]*dv1_sq*gamma+1.5811388300841898*surft1_upper[2]*dv1_sq*gamma-1.5811388300841898*surft1_lower[2]*dv1_sq*gamma; 
  out[25] = 0.7071067811865475*surft1_upper[13]*dv1_sq*gamma-0.7071067811865475*surft1_lower[13]*dv1_sq*gamma; 
  out[26] = 0.7071067811865475*surft1_upper[14]*dv1_sq*gamma-0.7071067811865475*surft1_lower[14]*dv1_sq*gamma; 
  out[27] = -(4.743416490252569*surft2_upper[3]*dv1_sq*gamma)-4.743416490252569*surft2_lower[3]*dv1_sq*gamma+1.5811388300841898*surft1_upper[3]*dv1_sq*gamma-1.5811388300841898*surft1_lower[3]*dv1_sq*gamma+6.7082039324993685*GCC[3]*dv1_sq*gamma; 
  out[28] = 0.7071067811865475*surft1_upper[15]*dv1_sq*gamma-0.7071067811865475*surft1_lower[15]*dv1_sq*gamma; 
  out[29] = 0.7071067811865475*surft1_upper[16]*dv1_sq*gamma-0.7071067811865475*surft1_lower[16]*dv1_sq*gamma; 
  out[30] = 1.224744871391589*surft1_upper[9]*dv1_sq*gamma+1.224744871391589*surft1_lower[9]*dv1_sq*gamma+6.7082039324993685*GCC[4]*dv1_sq*gamma-2.7386127875258306*surft2_upper[0]*dv1_sq*gamma+2.7386127875258306*surft2_lower[0]*dv1_sq*gamma; 
  out[31] = 1.224744871391589*surft1_upper[10]*dv1_sq*gamma+1.224744871391589*surft1_lower[10]*dv1_sq*gamma+3.0*GCC[5]*dv1_sq*gamma-2.1213203435596424*surft2_upper[4]*dv1_sq*gamma-2.1213203435596424*surft2_lower[4]*dv1_sq*gamma; 
  out[32] = -(1.224744871391589*surft2_upper[11]*dv1_sq*gamma)+1.224744871391589*surft2_lower[11]*dv1_sq*gamma+1.224744871391589*surft1_upper[11]*dv1_sq*gamma+1.224744871391589*surft1_lower[11]*dv1_sq*gamma; 
  out[33] = -(1.224744871391589*surft2_upper[12]*dv1_sq*gamma)+1.224744871391589*surft2_lower[12]*dv1_sq*gamma+1.224744871391589*surft1_upper[12]*dv1_sq*gamma+1.224744871391589*surft1_lower[12]*dv1_sq*gamma; 
  out[34] = -(2.7386127875258306*surft2_upper[10]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[10]*dv1_sq*gamma+1.5811388300841895*surft1_upper[4]*dv1_sq*gamma-1.5811388300841895*surft1_lower[4]*dv1_sq*gamma; 
  out[35] = 0.7071067811865475*surft1_upper[17]*dv1_sq*gamma-0.7071067811865475*surft1_lower[17]*dv1_sq*gamma; 
  out[36] = 0.7071067811865475*surft1_upper[18]*dv1_sq*gamma-0.7071067811865475*surft1_lower[18]*dv1_sq*gamma; 
  out[37] = 1.224744871391589*surft1_upper[13]*dv1_sq*gamma+1.224744871391589*surft1_lower[13]*dv1_sq*gamma+3.0*GCC[11]*dv1_sq*gamma-2.1213203435596424*surft2_upper[7]*dv1_sq*gamma-2.1213203435596424*surft2_lower[7]*dv1_sq*gamma; 
  out[38] = 1.224744871391589*surft1_upper[14]*dv1_sq*gamma+1.224744871391589*surft1_lower[14]*dv1_sq*gamma+3.0*GCC[12]*dv1_sq*gamma-2.1213203435596424*surft2_upper[8]*dv1_sq*gamma-2.1213203435596424*surft2_lower[8]*dv1_sq*gamma; 
  out[39] = 6.708203932499369*GCC[6]*dv1_sq*gamma-4.743416490252569*surft2_upper[5]*dv1_sq*gamma-4.743416490252569*surft2_lower[5]*dv1_sq*gamma+1.5811388300841895*surft1_upper[5]*dv1_sq*gamma-1.5811388300841895*surft1_lower[5]*dv1_sq*gamma; 
  out[40] = 6.708203932499369*GCC[7]*dv1_sq*gamma-4.743416490252569*surft2_upper[6]*dv1_sq*gamma-4.743416490252569*surft2_lower[6]*dv1_sq*gamma+1.5811388300841895*surft1_upper[6]*dv1_sq*gamma-1.5811388300841895*surft1_lower[6]*dv1_sq*gamma; 
  out[41] = 0.7071067811865475*surft1_upper[19]*dv1_sq*gamma-0.7071067811865475*surft1_lower[19]*dv1_sq*gamma; 
  out[42] = 1.224744871391589*surft1_upper[15]*dv1_sq*gamma+1.224744871391589*surft1_lower[15]*dv1_sq*gamma+6.708203932499369*GCC[8]*dv1_sq*gamma-2.7386127875258306*surft2_upper[1]*dv1_sq*gamma+2.7386127875258306*surft2_lower[1]*dv1_sq*gamma; 
  out[43] = 1.224744871391589*surft1_upper[16]*dv1_sq*gamma+1.224744871391589*surft1_lower[16]*dv1_sq*gamma+6.708203932499369*GCC[9]*dv1_sq*gamma-2.7386127875258306*surft2_upper[2]*dv1_sq*gamma+2.7386127875258306*surft2_lower[2]*dv1_sq*gamma; 
  out[44] = 3.0*GCC[19]*dv1_sq*gamma+1.224744871391589*surft1_upper[17]*dv1_sq*gamma+1.224744871391589*surft1_lower[17]*dv1_sq*gamma-2.1213203435596424*surft2_upper[11]*dv1_sq*gamma-2.1213203435596424*surft2_lower[11]*dv1_sq*gamma; 
  out[45] = 3.0*GCC[20]*dv1_sq*gamma+1.224744871391589*surft1_upper[18]*dv1_sq*gamma+1.224744871391589*surft1_lower[18]*dv1_sq*gamma-2.1213203435596424*surft2_upper[12]*dv1_sq*gamma-2.1213203435596424*surft2_lower[12]*dv1_sq*gamma; 
  out[46] = 6.7082039324993685*GCC[15]*dv1_sq*gamma-4.743416490252569*surft2_upper[10]*dv1_sq*gamma-4.743416490252569*surft2_lower[10]*dv1_sq*gamma+1.5811388300841898*surft1_upper[10]*dv1_sq*gamma-1.5811388300841898*surft1_lower[10]*dv1_sq*gamma; 
  out[47] = 1.224744871391589*surft1_upper[19]*dv1_sq*gamma+1.224744871391589*surft1_lower[19]*dv1_sq*gamma+6.7082039324993685*GCC[16]*dv1_sq*gamma-2.7386127875258306*surft2_upper[4]*dv1_sq*gamma+2.7386127875258306*surft2_lower[4]*dv1_sq*gamma; 
} 
