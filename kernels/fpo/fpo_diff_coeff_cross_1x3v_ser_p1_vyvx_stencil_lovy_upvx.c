#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vyvx_ser_p1_lovy_upvx(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
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
 
  double surft1_lo[8] = {0.0}; 
  double surft1_up[8] = {0.0}; 
  double surft2_lo[8] = {0.0}; 
  double surft2_up[8] = {0.0}; 
  double vol[40] = {0.0}; 
  double *out = &diff_coeff[120]; 

  const double* GCL = fpo_g_stencil[0]; 
  const double* GTL = fpo_g_stencil[1]; 
  const double* GCC = fpo_g_stencil[2]; 
  const double* G_surf_CC_vx = &fpo_g_surf_stencil[2][0]; 
  const double* GTC = fpo_g_stencil[3]; 
  const double* dGdvx_surf_CC_vy = &fpo_dgdv_surf[24]; 

  surft1_lo[0] = dGdvx_surf_CC_vy[0]/dv2; 
  surft1_lo[1] = dGdvx_surf_CC_vy[1]/dv2; 
  surft1_lo[2] = dGdvx_surf_CC_vy[2]/dv2; 
  surft1_lo[3] = dGdvx_surf_CC_vy[3]/dv2; 
  surft1_lo[4] = dGdvx_surf_CC_vy[4]/dv2; 
  surft1_lo[5] = dGdvx_surf_CC_vy[5]/dv2; 
  surft1_lo[6] = dGdvx_surf_CC_vy[6]/dv2; 
  surft1_lo[7] = dGdvx_surf_CC_vy[7]/dv2; 
  surft1_up[0] = -(0.7071067811865475*GTC[7])+0.7071067811865475*GCC[7]+0.6123724356957944*GTC[2]+0.6123724356957944*GCC[2]; 
  surft1_up[1] = -(0.7071067811865475*GTC[11])+0.7071067811865475*GCC[11]+0.6123724356957944*GTC[5]+0.6123724356957944*GCC[5]; 
  surft1_up[3] = -(0.7071067811865475*GTC[14])+0.7071067811865475*GCC[14]+0.6123724356957944*GTC[9]+0.6123724356957944*GCC[9]; 
  surft1_up[5] = -(0.7071067811865475*GTC[15])+0.7071067811865475*GCC[15]+0.6123724356957944*GTC[12]+0.6123724356957944*GCC[12]; 

  surft2_lo[0] = 0.408248290463863*GCL[2]-0.408248290463863*GCC[2]+0.3535533905932737*GCL[0]+0.3535533905932737*GCC[0]; 
  surft2_lo[1] = 0.408248290463863*GCL[5]-0.408248290463863*GCC[5]+0.3535533905932737*GCL[1]+0.3535533905932737*GCC[1]; 
  surft2_lo[2] = 0.408248290463863*GCL[7]-0.408248290463863*GCC[7]+0.3535533905932737*GCL[3]+0.3535533905932737*GCC[3]; 
  surft2_lo[3] = 0.408248290463863*GCL[9]-0.408248290463863*GCC[9]+0.3535533905932737*GCL[4]+0.3535533905932737*GCC[4]; 
  surft2_lo[4] = 0.408248290463863*GCL[11]-0.408248290463863*GCC[11]+0.3535533905932737*GCL[6]+0.3535533905932737*GCC[6]; 
  surft2_lo[5] = 0.408248290463863*GCL[12]-0.408248290463863*GCC[12]+0.3535533905932737*GCL[8]+0.3535533905932737*GCC[8]; 
  surft2_lo[6] = 0.408248290463863*GCL[14]-0.408248290463863*GCC[14]+0.3535533905932737*GCL[10]+0.3535533905932737*GCC[10]; 
  surft2_lo[7] = 0.408248290463863*GCL[15]-0.408248290463863*GCC[15]+0.3535533905932737*GCL[13]+0.3535533905932737*GCC[13]; 
  surft2_up[0] = 1.5811388300841895*GCC[16]+1.224744871391589*GCC[2]+0.7071067811865475*GCC[0]; 
  surft2_up[1] = 1.5811388300841898*GCC[17]+1.224744871391589*GCC[5]+0.7071067811865475*GCC[1]; 
  surft2_up[2] = 1.5811388300841898*GCC[18]+1.224744871391589*GCC[7]+0.7071067811865475*GCC[3]; 
  surft2_up[3] = 1.5811388300841898*GCC[19]+1.224744871391589*GCC[9]+0.7071067811865475*GCC[4]; 
  surft2_up[4] = 1.5811388300841895*GCC[20]+1.224744871391589*GCC[11]+0.7071067811865475*GCC[6]; 
  surft2_up[5] = 1.5811388300841895*GCC[21]+1.224744871391589*GCC[12]+0.7071067811865475*GCC[8]; 
  surft2_up[6] = 1.5811388300841895*GCC[22]+1.224744871391589*GCC[14]+0.7071067811865475*GCC[10]; 
  surft2_up[7] = 1.5811388300841898*GCC[23]+1.224744871391589*GCC[15]+0.7071067811865475*GCC[13]; 

  vol[7] = 3.0*GCC[0]; 
  vol[11] = 3.0*GCC[1]; 
  vol[14] = 3.0*GCC[4]; 
  vol[15] = 3.0*GCC[8]; 
  vol[18] = 6.7082039324993685*GCC[2]; 
  vol[20] = 6.708203932499369*GCC[5]; 
  vol[22] = 6.708203932499369*GCC[9]; 
  vol[23] = 6.7082039324993685*GCC[12]; 
  vol[26] = 6.7082039324993685*GCC[3]; 
  vol[28] = 6.708203932499369*GCC[6]; 
  vol[30] = 6.708203932499369*GCC[10]; 
  vol[31] = 6.7082039324993685*GCC[13]; 
  vol[38] = 3.0*GCC[32]; 
  vol[39] = 3.0*GCC[33]; 

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
  out[11] = (vol[11]+1.224744871391589*surft1_up[4]+1.224744871391589*surft1_lo[4]-2.1213203435596424*surft2_up[1]-2.1213203435596424*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[12] = (vol[12]+0.7071067811865475*surft1_up[7]-0.7071067811865475*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[13] = (vol[13]-1.224744871391589*surft2_up[5]+1.224744871391589*surft2_lo[5]+1.224744871391589*surft1_up[5]+1.224744871391589*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[14] = (vol[14]+1.224744871391589*surft1_up[6]+1.224744871391589*surft1_lo[6]-2.1213203435596424*surft2_up[3]-2.1213203435596424*surft2_lo[3])*dv1_sq*gamma_avg; 
  out[15] = (vol[15]+1.224744871391589*surft1_up[7]+1.224744871391589*surft1_lo[7]-2.1213203435596424*surft2_up[5]-2.1213203435596424*surft2_lo[5])*dv1_sq*gamma_avg; 
  out[16] = vol[16]*dv1_sq*gamma_avg; 
  out[17] = vol[17]*dv1_sq*gamma_avg; 
  out[18] = (vol[18]-2.7386127875258306*surft2_up[0]+2.7386127875258306*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[19] = vol[19]*dv1_sq*gamma_avg; 
  out[20] = (vol[20]-2.7386127875258306*surft2_up[1]+2.7386127875258306*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[21] = vol[21]*dv1_sq*gamma_avg; 
  out[22] = (vol[22]-2.7386127875258306*surft2_up[3]+2.7386127875258306*surft2_lo[3])*dv1_sq*gamma_avg; 
  out[23] = (vol[23]-2.7386127875258306*surft2_up[5]+2.7386127875258306*surft2_lo[5])*dv1_sq*gamma_avg; 
  out[24] = (vol[24]-2.7386127875258306*surft2_up[2]+2.7386127875258306*surft2_lo[2]+1.5811388300841895*surft1_up[0]-1.5811388300841895*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[25] = (vol[25]-2.7386127875258306*surft2_up[4]+2.7386127875258306*surft2_lo[4]+1.5811388300841898*surft1_up[1]-1.5811388300841898*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[26] = (vol[26]-4.743416490252569*surft2_up[2]-4.743416490252569*surft2_lo[2]+1.5811388300841898*surft1_up[2]-1.5811388300841898*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[27] = (vol[27]-2.7386127875258306*surft2_up[6]+2.7386127875258306*surft2_lo[6]+1.5811388300841898*surft1_up[3]-1.5811388300841898*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[28] = (vol[28]-4.743416490252569*surft2_up[4]-4.743416490252569*surft2_lo[4]+1.5811388300841895*surft1_up[4]-1.5811388300841895*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[29] = (vol[29]-2.7386127875258306*surft2_up[7]+2.7386127875258306*surft2_lo[7]+1.5811388300841895*surft1_up[5]-1.5811388300841895*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[30] = (vol[30]-4.743416490252569*surft2_up[6]-4.743416490252569*surft2_lo[6]+1.5811388300841895*surft1_up[6]-1.5811388300841895*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[31] = (vol[31]-4.743416490252569*surft2_up[7]-4.743416490252569*surft2_lo[7]+1.5811388300841898*surft1_up[7]-1.5811388300841898*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[32] = vol[32]*dv1_sq*gamma_avg; 
  out[33] = vol[33]*dv1_sq*gamma_avg; 
  out[34] = vol[34]*dv1_sq*gamma_avg; 
  out[35] = vol[35]*dv1_sq*gamma_avg; 
  out[36] = vol[36]*dv1_sq*gamma_avg; 
  out[37] = vol[37]*dv1_sq*gamma_avg; 
  out[38] = vol[38]*dv1_sq*gamma_avg; 
  out[39] = vol[39]*dv1_sq*gamma_avg; 
} 