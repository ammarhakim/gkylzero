#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_lovx_lovy(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1_pv1 = 2.0/dxv[1]; 
  double dv1_pv2 = 2.0/dxv[2]; 
  double dv1_sq = 4.0/dxv[1]/dxv[2]; 
 
  const double* GCC = fpo_g_stencil[0]; 
  const double* GTC = fpo_g_stencil[1]; 
  const double* GCR = fpo_g_stencil[2]; 
  const double* GTR = fpo_g_stencil[3]; 

  const double* g_surf_CC = fpo_g_surf_stencil[0]; 
  const double* g_surf_CC_pv2 = &g_surf_CC[8]; 
  const double* g_surf_CR = fpo_g_surf_stencil[2]; 
  const double* g_surf_CR_pv2 = &g_surf_CR[8]; 
  
  const double* g_surf_CC_pv1 = &g_surf_CC[0]; 
  const double* dgdpv1_surf_CC_pv2 = &fpo_dgdv_surf[24]; 
  const double* dgdpv2_surf_CC_pv1 = &fpo_dgdv_surf[8]; 
  const double* dgdpv1_surf_CC_pv1 = &fpo_dgdv_surf[0]; 
  
  double surft1_upper[8], surft1_lower[8]; 
  double surft2_upper[8], surft2_lower[8]; 
  
  double *diff_coeff_vxvy = &diff_coeff[40]; 
  double *diff_coeff_vxvz = &diff_coeff[80]; 
  double *diff_coeff_vyvx = &diff_coeff[120]; 
  double *diff_coeff_vyvz = &diff_coeff[200]; 
  double *diff_coeff_vzvx = &diff_coeff[240]; 
  double *diff_coeff_vzvy = &diff_coeff[280]; 
  
  double *out = diff_coeff_vxvy; 
  
  surft1_upper[0] = (0.14433756729740646*dgdpv1_surf_CC_pv1[2])/dv1_pv1+(0.08333333333333333*dgdpv1_surf_CC_pv1[0])/dv1_pv1+0.036828478186799324*GTR[7]-0.5524271728019898*GTC[7]-0.036828478186799324*GCR[7]+0.5524271728019898*GCC[7]-0.03827327723098713*GTR[3]+0.03827327723098713*GTC[3]+0.03827327723098713*GCR[3]-0.03827327723098713*GCC[3]-0.03189439769248927*GTR[2]+0.478415965387339*GTC[2]-0.03189439769248927*GCR[2]+0.478415965387339*GCC[2]+0.033145630368119385*GTR[0]-0.033145630368119385*GTC[0]+0.033145630368119385*GCR[0]-0.033145630368119385*GCC[0]; 
  surft1_upper[1] = (0.14433756729740646*dgdpv1_surf_CC_pv1[4])/dv1_pv1+(0.08333333333333333*dgdpv1_surf_CC_pv1[1])/dv1_pv1+0.036828478186799324*GTR[11]-0.5524271728019898*GTC[11]-0.036828478186799324*GCR[11]+0.5524271728019898*GCC[11]-0.03827327723098713*GTR[6]+0.03827327723098713*GTC[6]+0.03827327723098713*GCR[6]-0.03827327723098713*GCC[6]-0.03189439769248927*GTR[5]+0.478415965387339*GTC[5]-0.03189439769248927*GCR[5]+0.478415965387339*GCC[5]+0.033145630368119385*GTR[1]-0.033145630368119385*GTC[1]+0.033145630368119385*GCR[1]-0.033145630368119385*GCC[1]; 
  surft1_upper[2] = -((0.25*dgdpv1_surf_CC_pv1[2])/dv1_pv1)-(0.14433756729740646*dgdpv1_surf_CC_pv1[0])/dv1_pv1+0.34445949507888407*GTR[7]+0.5485836403108156*GTC[7]-0.34445949507888407*GCR[7]-0.5485836403108156*GCC[7]-0.2872621298570347*GTR[3]+0.2872621298570347*GTC[3]+0.2872621298570347*GCR[3]-0.2872621298570347*GCC[3]-0.2983106733130745*GTR[2]-0.4750873686097112*GTC[2]-0.2983106733130745*GCR[2]-0.4750873686097112*GCC[2]+0.24877630200141632*GTR[0]-0.24877630200141632*GTC[0]+0.24877630200141632*GCR[0]-0.24877630200141632*GCC[0]; 
  surft1_upper[3] = (0.14433756729740646*dgdpv1_surf_CC_pv1[6])/dv1_pv1+(0.08333333333333333*dgdpv1_surf_CC_pv1[3])/dv1_pv1+0.036828478186799324*GTR[14]-0.5524271728019898*GTC[14]-0.036828478186799324*GCR[14]+0.5524271728019898*GCC[14]-0.03827327723098713*GTR[10]+0.03827327723098713*GTC[10]+0.03827327723098713*GCR[10]-0.03827327723098713*GCC[10]-0.03189439769248927*GTR[9]+0.478415965387339*GTC[9]-0.03189439769248927*GCR[9]+0.478415965387339*GCC[9]+0.033145630368119385*GTR[4]-0.033145630368119385*GTC[4]+0.033145630368119385*GCR[4]-0.033145630368119385*GCC[4]; 
  surft1_upper[4] = -((0.25*dgdpv1_surf_CC_pv1[4])/dv1_pv1)-(0.14433756729740646*dgdpv1_surf_CC_pv1[1])/dv1_pv1+0.34445949507888407*GTR[11]+0.5485836403108156*GTC[11]-0.34445949507888407*GCR[11]-0.5485836403108156*GCC[11]-0.2872621298570347*GTR[6]+0.2872621298570347*GTC[6]+0.2872621298570347*GCR[6]-0.2872621298570347*GCC[6]-0.2983106733130745*GTR[5]-0.4750873686097112*GTC[5]-0.2983106733130745*GCR[5]-0.4750873686097112*GCC[5]+0.24877630200141632*GTR[1]-0.24877630200141632*GTC[1]+0.24877630200141632*GCR[1]-0.24877630200141632*GCC[1]; 
  surft1_upper[5] = (0.14433756729740646*dgdpv1_surf_CC_pv1[7])/dv1_pv1+(0.08333333333333333*dgdpv1_surf_CC_pv1[5])/dv1_pv1+0.036828478186799324*GTR[15]-0.5524271728019898*GTC[15]-0.036828478186799324*GCR[15]+0.5524271728019898*GCC[15]-0.03827327723098713*GTR[13]+0.03827327723098713*GTC[13]+0.03827327723098713*GCR[13]-0.03827327723098713*GCC[13]-0.03189439769248927*GTR[12]+0.478415965387339*GTC[12]-0.03189439769248927*GCR[12]+0.478415965387339*GCC[12]+0.033145630368119385*GTR[8]-0.033145630368119385*GTC[8]+0.033145630368119385*GCR[8]-0.033145630368119385*GCC[8]; 
  surft1_upper[6] = -((0.25*dgdpv1_surf_CC_pv1[6])/dv1_pv1)-(0.14433756729740646*dgdpv1_surf_CC_pv1[3])/dv1_pv1+0.34445949507888407*GTR[14]+0.5485836403108156*GTC[14]-0.34445949507888407*GCR[14]-0.5485836403108156*GCC[14]-0.2872621298570347*GTR[10]+0.2872621298570347*GTC[10]+0.2872621298570347*GCR[10]-0.2872621298570347*GCC[10]-0.2983106733130745*GTR[9]-0.4750873686097112*GTC[9]-0.2983106733130745*GCR[9]-0.4750873686097112*GCC[9]+0.24877630200141632*GTR[4]-0.24877630200141632*GTC[4]+0.24877630200141632*GCR[4]-0.24877630200141632*GCC[4]; 
  surft1_upper[7] = -((0.25*dgdpv1_surf_CC_pv1[7])/dv1_pv1)-(0.14433756729740646*dgdpv1_surf_CC_pv1[5])/dv1_pv1+0.34445949507888407*GTR[15]+0.5485836403108156*GTC[15]-0.34445949507888407*GCR[15]-0.5485836403108156*GCC[15]-0.2872621298570347*GTR[13]+0.2872621298570347*GTC[13]+0.2872621298570347*GCR[13]-0.2872621298570347*GCC[13]-0.2983106733130745*GTR[12]-0.4750873686097112*GTC[12]-0.2983106733130745*GCR[12]-0.4750873686097112*GCC[12]+0.24877630200141632*GTR[8]-0.24877630200141632*GTC[8]+0.24877630200141632*GCR[8]-0.24877630200141632*GCC[8]; 
  surft1_lower[0] = dgdpv1_surf_CC_pv2[0]/dv1_pv1; 
  surft1_lower[1] = dgdpv1_surf_CC_pv2[1]/dv1_pv1; 
  surft1_lower[2] = dgdpv1_surf_CC_pv2[2]/dv1_pv1; 
  surft1_lower[3] = dgdpv1_surf_CC_pv2[3]/dv1_pv1; 
  surft1_lower[4] = dgdpv1_surf_CC_pv2[4]/dv1_pv1; 
  surft1_lower[5] = dgdpv1_surf_CC_pv2[5]/dv1_pv1; 
  surft1_lower[6] = dgdpv1_surf_CC_pv2[6]/dv1_pv1; 
  surft1_lower[7] = dgdpv1_surf_CC_pv2[7]/dv1_pv1; 

  surft2_upper[0] = -(0.408248290463863*GCR[2])+0.408248290463863*GCC[2]+0.3535533905932737*GCR[0]+0.3535533905932737*GCC[0]; 
  surft2_upper[1] = -(0.408248290463863*GCR[5])+0.408248290463863*GCC[5]+0.3535533905932737*GCR[1]+0.3535533905932737*GCC[1]; 
  surft2_upper[2] = -(0.408248290463863*GCR[7])+0.408248290463863*GCC[7]+0.3535533905932737*GCR[3]+0.3535533905932737*GCC[3]; 
  surft2_upper[3] = -(0.408248290463863*GCR[9])+0.408248290463863*GCC[9]+0.3535533905932737*GCR[4]+0.3535533905932737*GCC[4]; 
  surft2_upper[4] = -(0.408248290463863*GCR[11])+0.408248290463863*GCC[11]+0.3535533905932737*GCR[6]+0.3535533905932737*GCC[6]; 
  surft2_upper[5] = -(0.408248290463863*GCR[12])+0.408248290463863*GCC[12]+0.3535533905932737*GCR[8]+0.3535533905932737*GCC[8]; 
  surft2_upper[6] = -(0.408248290463863*GCR[14])+0.408248290463863*GCC[14]+0.3535533905932737*GCR[10]+0.3535533905932737*GCC[10]; 
  surft2_upper[7] = -(0.408248290463863*GCR[15])+0.408248290463863*GCC[15]+0.3535533905932737*GCR[13]+0.3535533905932737*GCC[13]; 
  surft2_lower[0] = g_surf_CC_pv1[0]; 
  surft2_lower[1] = g_surf_CC_pv1[1]; 
  surft2_lower[2] = g_surf_CC_pv1[2]; 
  surft2_lower[3] = g_surf_CC_pv1[3]; 
  surft2_lower[4] = g_surf_CC_pv1[4]; 
  surft2_lower[5] = g_surf_CC_pv1[5]; 
  surft2_lower[6] = g_surf_CC_pv1[6]; 
  surft2_lower[7] = g_surf_CC_pv1[7]; 

  out[0] = (0.7071067811865475*surft1_upper[0]-0.7071067811865475*surft1_lower[0])*dv1_sq*gamma_avg; 
  out[1] = (0.7071067811865475*surft1_upper[1]-0.7071067811865475*surft1_lower[1])*dv1_sq*gamma_avg; 
  out[2] = (0.7071067811865475*surft1_upper[2]-0.7071067811865475*surft1_lower[2])*dv1_sq*gamma_avg; 
  out[3] = (-(1.224744871391589*surft2_upper[0])+1.224744871391589*surft2_lower[0]+1.224744871391589*surft1_upper[0]+1.224744871391589*surft1_lower[0])*dv1_sq*gamma_avg; 
  out[4] = (0.7071067811865475*surft1_upper[3]-0.7071067811865475*surft1_lower[3])*dv1_sq*gamma_avg; 
  out[5] = (0.7071067811865475*surft1_upper[4]-0.7071067811865475*surft1_lower[4])*dv1_sq*gamma_avg; 
  out[6] = (-(1.224744871391589*surft2_upper[1])+1.224744871391589*surft2_lower[1]+1.224744871391589*surft1_upper[1]+1.224744871391589*surft1_lower[1])*dv1_sq*gamma_avg; 
  out[7] = (1.224744871391589*surft1_upper[2]+1.224744871391589*surft1_lower[2]-2.1213203435596424*surft2_upper[0]-2.1213203435596424*surft2_lower[0]+3.0*GCC[0])*dv1_sq*gamma_avg; 
  out[8] = (0.7071067811865475*surft1_upper[5]-0.7071067811865475*surft1_lower[5])*dv1_sq*gamma_avg; 
  out[9] = (0.7071067811865475*surft1_upper[6]-0.7071067811865475*surft1_lower[6])*dv1_sq*gamma_avg; 
  out[10] = (-(1.224744871391589*surft2_upper[3])+1.224744871391589*surft2_lower[3]+1.224744871391589*surft1_upper[3]+1.224744871391589*surft1_lower[3])*dv1_sq*gamma_avg; 
  out[11] = (1.224744871391589*surft1_upper[4]+1.224744871391589*surft1_lower[4]-2.1213203435596424*surft2_upper[1]-2.1213203435596424*surft2_lower[1]+3.0*GCC[1])*dv1_sq*gamma_avg; 
  out[12] = (0.7071067811865475*surft1_upper[7]-0.7071067811865475*surft1_lower[7])*dv1_sq*gamma_avg; 
  out[13] = (-(1.224744871391589*surft2_upper[5])+1.224744871391589*surft2_lower[5]+1.224744871391589*surft1_upper[5]+1.224744871391589*surft1_lower[5])*dv1_sq*gamma_avg; 
  out[14] = (1.224744871391589*surft1_upper[6]+1.224744871391589*surft1_lower[6]+3.0*GCC[4]-2.1213203435596424*surft2_upper[3]-2.1213203435596424*surft2_lower[3])*dv1_sq*gamma_avg; 
  out[15] = (3.0*GCC[8]+1.224744871391589*surft1_upper[7]+1.224744871391589*surft1_lower[7]-2.1213203435596424*surft2_upper[5]-2.1213203435596424*surft2_lower[5])*dv1_sq*gamma_avg; 
  out[18] = (6.7082039324993685*GCC[2]-2.7386127875258306*surft2_upper[0]+2.7386127875258306*surft2_lower[0])*dv1_sq*gamma_avg; 
  out[20] = (6.708203932499369*GCC[5]-2.7386127875258306*surft2_upper[1]+2.7386127875258306*surft2_lower[1])*dv1_sq*gamma_avg; 
  out[22] = (6.708203932499369*GCC[9]-2.7386127875258306*surft2_upper[3]+2.7386127875258306*surft2_lower[3])*dv1_sq*gamma_avg; 
  out[23] = (6.7082039324993685*GCC[12]-2.7386127875258306*surft2_upper[5]+2.7386127875258306*surft2_lower[5])*dv1_sq*gamma_avg; 
  out[24] = (-(2.7386127875258306*surft2_upper[2])+2.7386127875258306*surft2_lower[2]+1.5811388300841895*surft1_upper[0]-1.5811388300841895*surft1_lower[0])*dv1_sq*gamma_avg; 
  out[25] = (-(2.7386127875258306*surft2_upper[4])+2.7386127875258306*surft2_lower[4]+1.5811388300841898*surft1_upper[1]-1.5811388300841898*surft1_lower[1])*dv1_sq*gamma_avg; 
  out[26] = (6.7082039324993685*GCC[3]-4.743416490252569*surft2_upper[2]-4.743416490252569*surft2_lower[2]+1.5811388300841898*surft1_upper[2]-1.5811388300841898*surft1_lower[2])*dv1_sq*gamma_avg; 
  out[27] = (-(2.7386127875258306*surft2_upper[6])+2.7386127875258306*surft2_lower[6]+1.5811388300841898*surft1_upper[3]-1.5811388300841898*surft1_lower[3])*dv1_sq*gamma_avg; 
  out[28] = (6.708203932499369*GCC[6]-4.743416490252569*surft2_upper[4]-4.743416490252569*surft2_lower[4]+1.5811388300841895*surft1_upper[4]-1.5811388300841895*surft1_lower[4])*dv1_sq*gamma_avg; 
  out[29] = (-(2.7386127875258306*surft2_upper[7])+2.7386127875258306*surft2_lower[7]+1.5811388300841895*surft1_upper[5]-1.5811388300841895*surft1_lower[5])*dv1_sq*gamma_avg; 
  out[30] = (6.708203932499369*GCC[10]-4.743416490252569*surft2_upper[6]-4.743416490252569*surft2_lower[6]+1.5811388300841895*surft1_upper[6]-1.5811388300841895*surft1_lower[6])*dv1_sq*gamma_avg; 
  out[31] = (6.7082039324993685*GCC[13]-4.743416490252569*surft2_upper[7]-4.743416490252569*surft2_lower[7]+1.5811388300841898*surft1_upper[7]-1.5811388300841898*surft1_lower[7])*dv1_sq*gamma_avg; 
  out[38] = 3.0*GCC[32]*dv1_sq*gamma_avg; 
  out[39] = 3.0*GCC[33]*dv1_sq*gamma_avg; 
} 

