#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_lovy_upvx(const double *dxv, double gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // diff_coeff: Output array for diffusion tensor. 

  double dv1_pv1 = 2.0/dxv[2]; 
  double dv1_pv2 = 2.0/dxv[1]; 
  double dv1_sq = 4.0/dxv[2]/dxv[1]; 
 
  const double* GBC = fpo_g_stencil[0]; 
  const double* GBR = fpo_g_stencil[1]; 
  const double* GCC = fpo_g_stencil[2]; 
  const double* GCR = fpo_g_stencil[3]; 

  const double* g_surf_CC = fpo_g_surf_stencil[2]; 
  const double* g_surf_CC_pv2 = &g_surf_CC[0]; 
  const double* g_surf_CR = fpo_g_surf_stencil[3]; 
  const double* g_surf_CR_pv2 = &g_surf_CR[0]; 
  
  const double* g_surf_CC_pv1 = &g_surf_CC[8]; 
  const double* dgdpv1_surf_CC_pv2 = &fpo_dgdv_surf[8]; 
  const double* dgdpv2_surf_CC_pv1 = &fpo_dgdv_surf[24]; 
  const double* dgdpv1_surf_CC_pv1 = &fpo_dgdv_surf[32]; 
  
  double surft1_upper[8], surft1_lower[8]; 
  double surft2_upper[8], surft2_lower[8]; 
  
  double *diff_coeff_vxvy = &diff_coeff[40]; 
  double *diff_coeff_vxvz = &diff_coeff[80]; 
  double *diff_coeff_vyvx = &diff_coeff[120]; 
  double *diff_coeff_vyvz = &diff_coeff[200]; 
  double *diff_coeff_vzvx = &diff_coeff[240]; 
  double *diff_coeff_vzvy = &diff_coeff[280]; 
  
  double *out = diff_coeff_vyvx; 
  
  surft1_upper[0] = dgdpv1_surf_CC_pv2[0]/dv1_pv1; 
  surft1_upper[1] = dgdpv1_surf_CC_pv2[1]/dv1_pv1; 
  surft1_upper[2] = dgdpv1_surf_CC_pv2[2]/dv1_pv1; 
  surft1_upper[3] = dgdpv1_surf_CC_pv2[3]/dv1_pv1; 
  surft1_upper[4] = dgdpv1_surf_CC_pv2[4]/dv1_pv1; 
  surft1_upper[5] = dgdpv1_surf_CC_pv2[5]/dv1_pv1; 
  surft1_upper[6] = dgdpv1_surf_CC_pv2[6]/dv1_pv1; 
  surft1_upper[7] = dgdpv1_surf_CC_pv2[7]/dv1_pv1; 
  surft1_lower[0] = -((0.11547005383792518*dgdpv1_surf_CC_pv1[2])/dv1_pv1)+(0.06666666666666667*dgdpv1_surf_CC_pv1[0])/dv1_pv1+0.011785113019775787*GCR[7]-0.6010407640085651*GCC[7]-0.011785113019775787*GBR[7]+0.6010407640085651*GBC[7]-0.010206207261596573*GCR[3]+0.5205165703414251*GCC[3]-0.010206207261596573*GBR[3]+0.5205165703414251*GBC[3]-0.020412414523193145*GCR[2]+0.020412414523193145*GCC[2]+0.020412414523193145*GBR[2]-0.020412414523193145*GBC[2]+0.01767766952966368*GCR[0]-0.01767766952966368*GCC[0]+0.01767766952966368*GBR[0]-0.01767766952966368*GBC[0]; 
  surft1_lower[1] = -((0.11547005383792518*dgdpv1_surf_CC_pv1[4])/dv1_pv1)+(0.06666666666666667*dgdpv1_surf_CC_pv1[1])/dv1_pv1+0.011785113019775787*GCR[11]-0.6010407640085651*GCC[11]-0.011785113019775787*GBR[11]+0.6010407640085651*GBC[11]-0.010206207261596573*GCR[6]+0.5205165703414251*GCC[6]-0.010206207261596573*GBR[6]+0.5205165703414251*GBC[6]-0.020412414523193145*GCR[5]+0.020412414523193145*GCC[5]+0.020412414523193145*GBR[5]-0.020412414523193145*GBC[5]+0.01767766952966368*GCR[1]-0.01767766952966368*GCC[1]+0.01767766952966368*GBR[1]-0.01767766952966368*GBC[1]; 
  surft1_lower[2] = (0.2*dgdpv1_surf_CC_pv1[2])/dv1_pv1-(0.11547005383792518*dgdpv1_surf_CC_pv1[0])/dv1_pv1+0.3878358759406697*GCR[7]+0.6327848502189876*GCC[7]-0.3878358759406697*GBR[7]-0.6327848502189876*GBC[7]-0.3358757210636099*GCR[3]-0.5480077554195741*GCC[3]-0.3358757210636099*GBR[3]-0.5480077554195741*GBC[3]-0.31819805153394626*GCR[2]+0.31819805153394626*GCC[2]+0.31819805153394626*GBR[2]-0.31819805153394626*GBC[2]+0.2755675960631073*GCR[0]-0.2755675960631073*GCC[0]+0.2755675960631073*GBR[0]-0.2755675960631073*GBC[0]; 
  surft1_lower[3] = -((0.11547005383792518*dgdpv1_surf_CC_pv1[6])/dv1_pv1)+(0.06666666666666667*dgdpv1_surf_CC_pv1[3])/dv1_pv1+0.011785113019775787*GCR[14]-0.6010407640085651*GCC[14]-0.011785113019775787*GBR[14]+0.6010407640085651*GBC[14]-0.010206207261596573*GCR[10]+0.5205165703414251*GCC[10]-0.010206207261596573*GBR[10]+0.5205165703414251*GBC[10]-0.020412414523193145*GCR[9]+0.020412414523193145*GCC[9]+0.020412414523193145*GBR[9]-0.020412414523193145*GBC[9]+0.01767766952966368*GCR[4]-0.01767766952966368*GCC[4]+0.01767766952966368*GBR[4]-0.01767766952966368*GBC[4]; 
  surft1_lower[4] = (0.2*dgdpv1_surf_CC_pv1[4])/dv1_pv1-(0.11547005383792518*dgdpv1_surf_CC_pv1[1])/dv1_pv1+0.3878358759406697*GCR[11]+0.6327848502189876*GCC[11]-0.3878358759406697*GBR[11]-0.6327848502189876*GBC[11]-0.3358757210636099*GCR[6]-0.5480077554195741*GCC[6]-0.3358757210636099*GBR[6]-0.5480077554195741*GBC[6]-0.31819805153394626*GCR[5]+0.31819805153394626*GCC[5]+0.31819805153394626*GBR[5]-0.31819805153394626*GBC[5]+0.2755675960631073*GCR[1]-0.2755675960631073*GCC[1]+0.2755675960631073*GBR[1]-0.2755675960631073*GBC[1]; 
  surft1_lower[5] = -((0.11547005383792518*dgdpv1_surf_CC_pv1[7])/dv1_pv1)+(0.06666666666666667*dgdpv1_surf_CC_pv1[5])/dv1_pv1+0.011785113019775787*GCR[15]-0.6010407640085651*GCC[15]-0.011785113019775787*GBR[15]+0.6010407640085651*GBC[15]-0.010206207261596573*GCR[13]+0.5205165703414251*GCC[13]-0.010206207261596573*GBR[13]+0.5205165703414251*GBC[13]-0.020412414523193145*GCR[12]+0.020412414523193145*GCC[12]+0.020412414523193145*GBR[12]-0.020412414523193145*GBC[12]+0.01767766952966368*GCR[8]-0.01767766952966368*GCC[8]+0.01767766952966368*GBR[8]-0.01767766952966368*GBC[8]; 
  surft1_lower[6] = (0.2*dgdpv1_surf_CC_pv1[6])/dv1_pv1-(0.11547005383792518*dgdpv1_surf_CC_pv1[3])/dv1_pv1+0.3878358759406697*GCR[14]+0.6327848502189876*GCC[14]-0.3878358759406697*GBR[14]-0.6327848502189876*GBC[14]-0.3358757210636099*GCR[10]-0.5480077554195741*GCC[10]-0.3358757210636099*GBR[10]-0.5480077554195741*GBC[10]-0.31819805153394626*GCR[9]+0.31819805153394626*GCC[9]+0.31819805153394626*GBR[9]-0.31819805153394626*GBC[9]+0.2755675960631073*GCR[4]-0.2755675960631073*GCC[4]+0.2755675960631073*GBR[4]-0.2755675960631073*GBC[4]; 
  surft1_lower[7] = (0.2*dgdpv1_surf_CC_pv1[7])/dv1_pv1-(0.11547005383792518*dgdpv1_surf_CC_pv1[5])/dv1_pv1+0.3878358759406697*GCR[15]+0.6327848502189876*GCC[15]-0.3878358759406697*GBR[15]-0.6327848502189876*GBC[15]-0.3358757210636099*GCR[13]-0.5480077554195741*GCC[13]-0.3358757210636099*GBR[13]-0.5480077554195741*GBC[13]-0.31819805153394626*GCR[12]+0.31819805153394626*GCC[12]+0.31819805153394626*GBR[12]-0.31819805153394626*GBC[12]+0.2755675960631073*GCR[8]-0.2755675960631073*GCC[8]+0.2755675960631073*GBR[8]-0.2755675960631073*GBC[8]; 

  surft2_upper[0] = -(0.408248290463863*GCR[3])+0.408248290463863*GCC[3]+0.3535533905932737*GCR[0]+0.3535533905932737*GCC[0]; 
  surft2_upper[1] = -(0.408248290463863*GCR[6])+0.408248290463863*GCC[6]+0.3535533905932737*GCR[1]+0.3535533905932737*GCC[1]; 
  surft2_upper[2] = -(0.408248290463863*GCR[7])+0.408248290463863*GCC[7]+0.3535533905932737*GCR[2]+0.3535533905932737*GCC[2]; 
  surft2_upper[3] = -(0.408248290463863*GCR[10])+0.408248290463863*GCC[10]+0.3535533905932737*GCR[4]+0.3535533905932737*GCC[4]; 
  surft2_upper[4] = -(0.408248290463863*GCR[11])+0.408248290463863*GCC[11]+0.3535533905932737*GCR[5]+0.3535533905932737*GCC[5]; 
  surft2_upper[5] = -(0.408248290463863*GCR[13])+0.408248290463863*GCC[13]+0.3535533905932737*GCR[8]+0.3535533905932737*GCC[8]; 
  surft2_upper[6] = -(0.408248290463863*GCR[14])+0.408248290463863*GCC[14]+0.3535533905932737*GCR[9]+0.3535533905932737*GCC[9]; 
  surft2_upper[7] = -(0.408248290463863*GCR[15])+0.408248290463863*GCC[15]+0.3535533905932737*GCR[12]+0.3535533905932737*GCC[12]; 
  surft2_lower[0] = g_surf_CC_pv1[0]; 
  surft2_lower[1] = g_surf_CC_pv1[1]; 
  surft2_lower[2] = g_surf_CC_pv1[2]; 
  surft2_lower[3] = g_surf_CC_pv1[3]; 
  surft2_lower[4] = g_surf_CC_pv1[4]; 
  surft2_lower[5] = g_surf_CC_pv1[5]; 
  surft2_lower[6] = g_surf_CC_pv1[6]; 
  surft2_lower[7] = g_surf_CC_pv1[7]; 

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
  out[11] = 1.224744871391589*surft1_upper[4]*dv1_sq*gamma+1.224744871391589*surft1_lower[4]*dv1_sq*gamma-2.1213203435596424*surft2_upper[1]*dv1_sq*gamma-2.1213203435596424*surft2_lower[1]*dv1_sq*gamma+3.0*GCC[1]*dv1_sq*gamma; 
  out[12] = -(1.224744871391589*surft2_upper[5]*dv1_sq*gamma)+1.224744871391589*surft2_lower[5]*dv1_sq*gamma+1.224744871391589*surft1_upper[5]*dv1_sq*gamma+1.224744871391589*surft1_lower[5]*dv1_sq*gamma; 
  out[13] = 0.7071067811865475*surft1_upper[7]*dv1_sq*gamma-0.7071067811865475*surft1_lower[7]*dv1_sq*gamma; 
  out[14] = 1.224744871391589*surft1_upper[6]*dv1_sq*gamma+1.224744871391589*surft1_lower[6]*dv1_sq*gamma+3.0*GCC[4]*dv1_sq*gamma-2.1213203435596424*surft2_upper[3]*dv1_sq*gamma-2.1213203435596424*surft2_lower[3]*dv1_sq*gamma; 
  out[15] = 3.0*GCC[8]*dv1_sq*gamma+1.224744871391589*surft1_upper[7]*dv1_sq*gamma+1.224744871391589*surft1_lower[7]*dv1_sq*gamma-2.1213203435596424*surft2_upper[5]*dv1_sq*gamma-2.1213203435596424*surft2_lower[5]*dv1_sq*gamma; 
  out[16] = -(2.7386127875258306*surft2_upper[2]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[2]*dv1_sq*gamma+1.5811388300841895*surft1_upper[0]*dv1_sq*gamma-1.5811388300841895*surft1_lower[0]*dv1_sq*gamma; 
  out[17] = -(2.7386127875258306*surft2_upper[4]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[4]*dv1_sq*gamma+1.5811388300841898*surft1_upper[1]*dv1_sq*gamma-1.5811388300841898*surft1_lower[1]*dv1_sq*gamma; 
  out[18] = -(4.743416490252569*surft2_upper[2]*dv1_sq*gamma)-4.743416490252569*surft2_lower[2]*dv1_sq*gamma+1.5811388300841898*surft1_upper[2]*dv1_sq*gamma-1.5811388300841898*surft1_lower[2]*dv1_sq*gamma+6.7082039324993685*GCC[2]*dv1_sq*gamma; 
  out[19] = -(2.7386127875258306*surft2_upper[6]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[6]*dv1_sq*gamma+1.5811388300841898*surft1_upper[3]*dv1_sq*gamma-1.5811388300841898*surft1_lower[3]*dv1_sq*gamma; 
  out[20] = 6.708203932499369*GCC[5]*dv1_sq*gamma-4.743416490252569*surft2_upper[4]*dv1_sq*gamma-4.743416490252569*surft2_lower[4]*dv1_sq*gamma+1.5811388300841895*surft1_upper[4]*dv1_sq*gamma-1.5811388300841895*surft1_lower[4]*dv1_sq*gamma; 
  out[21] = -(2.7386127875258306*surft2_upper[7]*dv1_sq*gamma)+2.7386127875258306*surft2_lower[7]*dv1_sq*gamma+1.5811388300841895*surft1_upper[5]*dv1_sq*gamma-1.5811388300841895*surft1_lower[5]*dv1_sq*gamma; 
  out[22] = 6.708203932499369*GCC[9]*dv1_sq*gamma-4.743416490252569*surft2_upper[6]*dv1_sq*gamma-4.743416490252569*surft2_lower[6]*dv1_sq*gamma+1.5811388300841895*surft1_upper[6]*dv1_sq*gamma-1.5811388300841895*surft1_lower[6]*dv1_sq*gamma; 
  out[23] = 6.7082039324993685*GCC[12]*dv1_sq*gamma-4.743416490252569*surft2_upper[7]*dv1_sq*gamma-4.743416490252569*surft2_lower[7]*dv1_sq*gamma+1.5811388300841898*surft1_upper[7]*dv1_sq*gamma-1.5811388300841898*surft1_lower[7]*dv1_sq*gamma; 
  out[26] = 6.7082039324993685*GCC[3]*dv1_sq*gamma-2.7386127875258306*surft2_upper[0]*dv1_sq*gamma+2.7386127875258306*surft2_lower[0]*dv1_sq*gamma; 
  out[28] = 6.708203932499369*GCC[6]*dv1_sq*gamma-2.7386127875258306*surft2_upper[1]*dv1_sq*gamma+2.7386127875258306*surft2_lower[1]*dv1_sq*gamma; 
  out[30] = 6.708203932499369*GCC[10]*dv1_sq*gamma-2.7386127875258306*surft2_upper[3]*dv1_sq*gamma+2.7386127875258306*surft2_lower[3]*dv1_sq*gamma; 
  out[31] = 6.7082039324993685*GCC[13]*dv1_sq*gamma-2.7386127875258306*surft2_upper[5]*dv1_sq*gamma+2.7386127875258306*surft2_lower[5]*dv1_sq*gamma; 
  out[38] = 3.0*GCC[32]*dv1_sq*gamma; 
  out[39] = 3.0*GCC[33]*dv1_sq*gamma; 
} 
