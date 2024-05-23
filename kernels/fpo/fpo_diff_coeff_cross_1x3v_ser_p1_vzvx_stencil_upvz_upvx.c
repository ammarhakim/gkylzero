#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vzvx_ser_p1_upvz_upvx(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
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
 
  const double* GBL = fpo_g_stencil[0]; 
  const double* GBC = fpo_g_stencil[1]; 
  const double* GCL = fpo_g_stencil[2]; 
  const double* GCC = fpo_g_stencil[3]; 

  const double* g_surf_CL = fpo_g_surf_stencil[2]; 
  const double* g_surf_CL_pv2 = &g_surf_CL[0]; 
  const double* g_surf_CC = fpo_g_surf_stencil[3]; 
  const double* g_surf_CC_pv2 = &g_surf_CC[0]; 
  
  const double* g_surf_CC_pv1 = &g_surf_CC[16]; 
  const double* dgdpv1_surf_CC_pv2 = &fpo_dgdv_surf[16]; 
  const double* dgdpv2_surf_CC_pv1 = &fpo_dgdv_surf[48]; 
  const double* dgdpv1_surf_CC_pv1 = &fpo_dgdv_surf[64]; 
  
  double surft1_upper[8], surft1_lower[8]; 
  double surft2_upper[8], surft2_lower[8]; 
  
  double *diff_coeff_vxvy = &diff_coeff[40]; 
  double *diff_coeff_vxvz = &diff_coeff[80]; 
  double *diff_coeff_vyvx = &diff_coeff[120]; 
  double *diff_coeff_vyvz = &diff_coeff[200]; 
  double *diff_coeff_vzvx = &diff_coeff[240]; 
  double *diff_coeff_vzvy = &diff_coeff[280]; 
  
  double *out = diff_coeff_vzvx; 
  
  surft1_upper[0] = dgdpv1_surf_CC_pv2[0]/dv1_pv1; 
  surft1_upper[1] = dgdpv1_surf_CC_pv2[1]/dv1_pv1; 
  surft1_upper[2] = dgdpv1_surf_CC_pv2[2]/dv1_pv1; 
  surft1_upper[3] = dgdpv1_surf_CC_pv2[3]/dv1_pv1; 
  surft1_upper[4] = dgdpv1_surf_CC_pv2[4]/dv1_pv1; 
  surft1_upper[5] = dgdpv1_surf_CC_pv2[5]/dv1_pv1; 
  surft1_upper[6] = dgdpv1_surf_CC_pv2[6]/dv1_pv1; 
  surft1_upper[7] = dgdpv1_surf_CC_pv2[7]/dv1_pv1; 
  surft1_lower[0] = -((0.14433756729740646*dgdpv1_surf_CC_pv1[2])/dv1_pv1)+(0.08333333333333333*dgdpv1_surf_CC_pv1[0])/dv1_pv1+0.036828478186799324*GCL[9]-0.5524271728019898*GCC[9]-0.036828478186799324*GBL[9]+0.5524271728019898*GBC[9]-0.03189439769248927*GCL[4]+0.478415965387339*GCC[4]-0.03189439769248927*GBL[4]+0.478415965387339*GBC[4]+0.03827327723098713*GCL[2]-0.03827327723098713*GCC[2]-0.03827327723098713*GBL[2]+0.03827327723098713*GBC[2]-0.033145630368119385*GCL[0]+0.033145630368119385*GCC[0]-0.033145630368119385*GBL[0]+0.033145630368119385*GBC[0]; 
  surft1_lower[1] = -((0.14433756729740646*dgdpv1_surf_CC_pv1[4])/dv1_pv1)+(0.08333333333333333*dgdpv1_surf_CC_pv1[1])/dv1_pv1+0.036828478186799324*GCL[12]-0.5524271728019898*GCC[12]-0.036828478186799324*GBL[12]+0.5524271728019898*GBC[12]-0.03189439769248927*GCL[8]+0.478415965387339*GCC[8]-0.03189439769248927*GBL[8]+0.478415965387339*GBC[8]+0.03827327723098713*GCL[5]-0.03827327723098713*GCC[5]-0.03827327723098713*GBL[5]+0.03827327723098713*GBC[5]-0.033145630368119385*GCL[1]+0.033145630368119385*GCC[1]-0.033145630368119385*GBL[1]+0.033145630368119385*GBC[1]; 
  surft1_lower[2] = -((0.14433756729740646*dgdpv1_surf_CC_pv1[6])/dv1_pv1)+(0.08333333333333333*dgdpv1_surf_CC_pv1[3])/dv1_pv1+0.036828478186799324*GCL[14]-0.5524271728019898*GCC[14]-0.036828478186799324*GBL[14]+0.5524271728019898*GBC[14]-0.03189439769248927*GCL[10]+0.478415965387339*GCC[10]-0.03189439769248927*GBL[10]+0.478415965387339*GBC[10]+0.03827327723098713*GCL[7]-0.03827327723098713*GCC[7]-0.03827327723098713*GBL[7]+0.03827327723098713*GBC[7]-0.033145630368119385*GCL[3]+0.033145630368119385*GCC[3]-0.033145630368119385*GBL[3]+0.033145630368119385*GBC[3]; 
  surft1_lower[3] = -((0.25*dgdpv1_surf_CC_pv1[2])/dv1_pv1)+(0.14433756729740646*dgdpv1_surf_CC_pv1[0])/dv1_pv1-0.34445949507888407*GCL[9]-0.5485836403108156*GCC[9]+0.34445949507888407*GBL[9]+0.5485836403108156*GBC[9]+0.2983106733130745*GCL[4]+0.4750873686097112*GCC[4]+0.2983106733130745*GBL[4]+0.4750873686097112*GBC[4]-0.2872621298570347*GCL[2]+0.2872621298570347*GCC[2]+0.2872621298570347*GBL[2]-0.2872621298570347*GBC[2]+0.24877630200141632*GCL[0]-0.24877630200141632*GCC[0]+0.24877630200141632*GBL[0]-0.24877630200141632*GBC[0]; 
  surft1_lower[4] = -((0.14433756729740646*dgdpv1_surf_CC_pv1[7])/dv1_pv1)+(0.08333333333333333*dgdpv1_surf_CC_pv1[5])/dv1_pv1+0.036828478186799324*GCL[15]-0.5524271728019898*GCC[15]-0.036828478186799324*GBL[15]+0.5524271728019898*GBC[15]-0.03189439769248927*GCL[13]+0.478415965387339*GCC[13]-0.03189439769248927*GBL[13]+0.478415965387339*GBC[13]+0.03827327723098713*GCL[11]-0.03827327723098713*GCC[11]-0.03827327723098713*GBL[11]+0.03827327723098713*GBC[11]-0.033145630368119385*GCL[6]+0.033145630368119385*GCC[6]-0.033145630368119385*GBL[6]+0.033145630368119385*GBC[6]; 
  surft1_lower[5] = -((0.25*dgdpv1_surf_CC_pv1[4])/dv1_pv1)+(0.14433756729740646*dgdpv1_surf_CC_pv1[1])/dv1_pv1-0.34445949507888407*GCL[12]-0.5485836403108156*GCC[12]+0.34445949507888407*GBL[12]+0.5485836403108156*GBC[12]+0.2983106733130745*GCL[8]+0.4750873686097112*GCC[8]+0.2983106733130745*GBL[8]+0.4750873686097112*GBC[8]-0.2872621298570347*GCL[5]+0.2872621298570347*GCC[5]+0.2872621298570347*GBL[5]-0.2872621298570347*GBC[5]+0.24877630200141632*GCL[1]-0.24877630200141632*GCC[1]+0.24877630200141632*GBL[1]-0.24877630200141632*GBC[1]; 
  surft1_lower[6] = -((0.25*dgdpv1_surf_CC_pv1[6])/dv1_pv1)+(0.14433756729740646*dgdpv1_surf_CC_pv1[3])/dv1_pv1-0.34445949507888407*GCL[14]-0.5485836403108156*GCC[14]+0.34445949507888407*GBL[14]+0.5485836403108156*GBC[14]+0.2983106733130745*GCL[10]+0.4750873686097112*GCC[10]+0.2983106733130745*GBL[10]+0.4750873686097112*GBC[10]-0.2872621298570347*GCL[7]+0.2872621298570347*GCC[7]+0.2872621298570347*GBL[7]-0.2872621298570347*GBC[7]+0.24877630200141632*GCL[3]-0.24877630200141632*GCC[3]+0.24877630200141632*GBL[3]-0.24877630200141632*GBC[3]; 
  surft1_lower[7] = -((0.25*dgdpv1_surf_CC_pv1[7])/dv1_pv1)+(0.14433756729740646*dgdpv1_surf_CC_pv1[5])/dv1_pv1-0.34445949507888407*GCL[15]-0.5485836403108156*GCC[15]+0.34445949507888407*GBL[15]+0.5485836403108156*GBC[15]+0.2983106733130745*GCL[13]+0.4750873686097112*GCC[13]+0.2983106733130745*GBL[13]+0.4750873686097112*GBC[13]-0.2872621298570347*GCL[11]+0.2872621298570347*GCC[11]+0.2872621298570347*GBL[11]-0.2872621298570347*GBC[11]+0.24877630200141632*GCL[6]-0.24877630200141632*GCC[6]+0.24877630200141632*GBL[6]-0.24877630200141632*GBC[6]; 

  surft2_upper[0] = g_surf_CC_pv1[0]; 
  surft2_upper[1] = g_surf_CC_pv1[1]; 
  surft2_upper[2] = g_surf_CC_pv1[2]; 
  surft2_upper[3] = g_surf_CC_pv1[3]; 
  surft2_upper[4] = g_surf_CC_pv1[4]; 
  surft2_upper[5] = g_surf_CC_pv1[5]; 
  surft2_upper[6] = g_surf_CC_pv1[6]; 
  surft2_upper[7] = g_surf_CC_pv1[7]; 
  surft2_lower[0] = 0.408248290463863*GCL[4]-0.408248290463863*GCC[4]+0.3535533905932737*GCL[0]+0.3535533905932737*GCC[0]; 
  surft2_lower[1] = 0.408248290463863*GCL[8]-0.408248290463863*GCC[8]+0.3535533905932737*GCL[1]+0.3535533905932737*GCC[1]; 
  surft2_lower[2] = 0.408248290463863*GCL[9]-0.408248290463863*GCC[9]+0.3535533905932737*GCL[2]+0.3535533905932737*GCC[2]; 
  surft2_lower[3] = 0.408248290463863*GCL[10]-0.408248290463863*GCC[10]+0.3535533905932737*GCL[3]+0.3535533905932737*GCC[3]; 
  surft2_lower[4] = 0.408248290463863*GCL[12]-0.408248290463863*GCC[12]+0.3535533905932737*GCL[5]+0.3535533905932737*GCC[5]; 
  surft2_lower[5] = 0.408248290463863*GCL[13]-0.408248290463863*GCC[13]+0.3535533905932737*GCL[6]+0.3535533905932737*GCC[6]; 
  surft2_lower[6] = 0.408248290463863*GCL[14]-0.408248290463863*GCC[14]+0.3535533905932737*GCL[7]+0.3535533905932737*GCC[7]; 
  surft2_lower[7] = 0.408248290463863*GCL[15]-0.408248290463863*GCC[15]+0.3535533905932737*GCL[11]+0.3535533905932737*GCC[11]; 

  out[0] = (0.7071067811865475*surft1_upper[0]-0.7071067811865475*surft1_lower[0])*dv1_sq*gamma_avg; 
  out[1] = (0.7071067811865475*surft1_upper[1]-0.7071067811865475*surft1_lower[1])*dv1_sq*gamma_avg; 
  out[2] = (-(1.224744871391589*surft2_upper[0])+1.224744871391589*surft2_lower[0]+1.224744871391589*surft1_upper[0]+1.224744871391589*surft1_lower[0])*dv1_sq*gamma_avg; 
  out[3] = (0.7071067811865475*surft1_upper[2]-0.7071067811865475*surft1_lower[2])*dv1_sq*gamma_avg; 
  out[4] = (0.7071067811865475*surft1_upper[3]-0.7071067811865475*surft1_lower[3])*dv1_sq*gamma_avg; 
  out[5] = (-(1.224744871391589*surft2_upper[1])+1.224744871391589*surft2_lower[1]+1.224744871391589*surft1_upper[1]+1.224744871391589*surft1_lower[1])*dv1_sq*gamma_avg; 
  out[6] = (0.7071067811865475*surft1_upper[4]-0.7071067811865475*surft1_lower[4])*dv1_sq*gamma_avg; 
  out[7] = (-(1.224744871391589*surft2_upper[3])+1.224744871391589*surft2_lower[3]+1.224744871391589*surft1_upper[2]+1.224744871391589*surft1_lower[2])*dv1_sq*gamma_avg; 
  out[8] = (0.7071067811865475*surft1_upper[5]-0.7071067811865475*surft1_lower[5])*dv1_sq*gamma_avg; 
  out[9] = (1.224744871391589*surft1_upper[3]+1.224744871391589*surft1_lower[3]-2.1213203435596424*surft2_upper[0]-2.1213203435596424*surft2_lower[0]+3.0*GCC[0])*dv1_sq*gamma_avg; 
  out[10] = (0.7071067811865475*surft1_upper[6]-0.7071067811865475*surft1_lower[6])*dv1_sq*gamma_avg; 
  out[11] = (-(1.224744871391589*surft2_upper[5])+1.224744871391589*surft2_lower[5]+1.224744871391589*surft1_upper[4]+1.224744871391589*surft1_lower[4])*dv1_sq*gamma_avg; 
  out[12] = (1.224744871391589*surft1_upper[5]+1.224744871391589*surft1_lower[5]-2.1213203435596424*surft2_upper[1]-2.1213203435596424*surft2_lower[1]+3.0*GCC[1])*dv1_sq*gamma_avg; 
  out[13] = (0.7071067811865475*surft1_upper[7]-0.7071067811865475*surft1_lower[7])*dv1_sq*gamma_avg; 
  out[14] = (1.224744871391589*surft1_upper[6]+1.224744871391589*surft1_lower[6]-2.1213203435596424*surft2_upper[3]-2.1213203435596424*surft2_lower[3]+3.0*GCC[3])*dv1_sq*gamma_avg; 
  out[15] = (1.224744871391589*surft1_upper[7]+1.224744871391589*surft1_lower[7]+3.0*GCC[6]-2.1213203435596424*surft2_upper[5]-2.1213203435596424*surft2_lower[5])*dv1_sq*gamma_avg; 
  out[16] = (-(2.7386127875258306*surft2_upper[2])+2.7386127875258306*surft2_lower[2]+1.5811388300841895*surft1_upper[0]-1.5811388300841895*surft1_lower[0])*dv1_sq*gamma_avg; 
  out[17] = (-(2.7386127875258306*surft2_upper[4])+2.7386127875258306*surft2_lower[4]+1.5811388300841898*surft1_upper[1]-1.5811388300841898*surft1_lower[1])*dv1_sq*gamma_avg; 
  out[18] = (-(2.7386127875258306*surft2_upper[6])+2.7386127875258306*surft2_lower[6]+1.5811388300841898*surft1_upper[2]-1.5811388300841898*surft1_lower[2])*dv1_sq*gamma_avg; 
  out[19] = (1.5811388300841898*surft1_upper[3]-1.5811388300841898*surft1_lower[3]-4.743416490252569*surft2_upper[2]-4.743416490252569*surft2_lower[2]+6.7082039324993685*GCC[2])*dv1_sq*gamma_avg; 
  out[20] = (-(2.7386127875258306*surft2_upper[7])+2.7386127875258306*surft2_lower[7]+1.5811388300841895*surft1_upper[4]-1.5811388300841895*surft1_lower[4])*dv1_sq*gamma_avg; 
  out[21] = (1.5811388300841895*surft1_upper[5]-1.5811388300841895*surft1_lower[5]+6.708203932499369*GCC[5]-4.743416490252569*surft2_upper[4]-4.743416490252569*surft2_lower[4])*dv1_sq*gamma_avg; 
  out[22] = (6.708203932499369*GCC[7]-4.743416490252569*surft2_upper[6]-4.743416490252569*surft2_lower[6]+1.5811388300841895*surft1_upper[6]-1.5811388300841895*surft1_lower[6])*dv1_sq*gamma_avg; 
  out[23] = (6.7082039324993685*GCC[11]-4.743416490252569*surft2_upper[7]-4.743416490252569*surft2_lower[7]+1.5811388300841898*surft1_upper[7]-1.5811388300841898*surft1_lower[7])*dv1_sq*gamma_avg; 
  out[30] = 3.0*GCC[24]*dv1_sq*gamma_avg; 
  out[31] = 3.0*GCC[25]*dv1_sq*gamma_avg; 
  out[34] = (6.7082039324993685*GCC[4]-2.7386127875258306*surft2_upper[0]+2.7386127875258306*surft2_lower[0])*dv1_sq*gamma_avg; 
  out[36] = (6.708203932499369*GCC[8]-2.7386127875258306*surft2_upper[1]+2.7386127875258306*surft2_lower[1])*dv1_sq*gamma_avg; 
  out[38] = (6.708203932499369*GCC[10]-2.7386127875258306*surft2_upper[3]+2.7386127875258306*surft2_lower[3])*dv1_sq*gamma_avg; 
  out[39] = (6.7082039324993685*GCC[13]-2.7386127875258306*surft2_upper[5]+2.7386127875258306*surft2_lower[5])*dv1_sq*gamma_avg; 
} 

