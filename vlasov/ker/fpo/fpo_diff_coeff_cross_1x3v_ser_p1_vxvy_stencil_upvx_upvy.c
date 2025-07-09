#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_diff_coeff_cross_1x3v_vxvy_ser_p1_upvx_upvy(const double *dxv, const double *gamma, const double* fpo_g_stencil[9], const double* fpo_g_surf_stencil[9], const double* fpo_dgdv_surf, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[9]: 9 cell stencil of Rosenbluth potential G. 
  // fpo_g_surf_stencil[9]: 9 cell stencil of surface projection of G. 
  // fpo_dgdv_surf: Surface expansion of dG/dv in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[1]; 
  double dv2 = 2.0/dxv[2]; 
  double dv1_sq = 4.0/dxv[1]/dxv[2]; 
 
  double surft1_lo[16] = {0.0}; 
  double surft1_up[16] = {0.0}; 
  double surft2_lo[16] = {0.0}; 
  double surft2_up[16] = {0.0}; 
  double vol[40] = {0.0}; 
  double *out = &diff_coeff[40]; 

  const double* GBL = fpo_g_stencil[0]; 
  const double* GCL = fpo_g_stencil[1]; 
  const double* GTL = fpo_g_stencil[2]; 
  const double* GBC = fpo_g_stencil[3]; 
  const double* GCC = fpo_g_stencil[4]; 
  const double* G_surf_CC_vx = &fpo_g_surf_stencil[4][0]; 
  const double* G_surf_CC_vy = &fpo_g_surf_stencil[4][16]; 
  const double* GTC = fpo_g_stencil[5]; 
  const double* GBR = fpo_g_stencil[6]; 
  const double* GCR = fpo_g_stencil[7]; 
  const double* GTR = fpo_g_stencil[8]; 
  const double* dGdvx_surf_CC_vy = &fpo_dgdv_surf[32]; 

  surft1_lo[0] = 0.599071547271275*GCC[26]+0.599071547271275*GBC[26]-0.8617863895711042*GCC[7]+0.8617863895711042*GBC[7]+0.6123724356957944*GCC[2]+0.6123724356957944*GBC[2]; 
  surft1_lo[1] = 0.599071547271275*GCC[28]+0.599071547271275*GBC[28]-0.8617863895711042*GCC[11]+0.8617863895711042*GBC[11]+0.6123724356957944*GCC[5]+0.6123724356957944*GBC[5]; 
  surft1_lo[2] = -(1.9270129491651047*GCC[18])+1.9270129491651047*GBC[18]+1.369306393762915*GCC[16]+1.369306393762915*GBC[16]; 
  surft1_lo[3] = 0.599071547271275*GCC[30]+0.599071547271275*GBC[30]-0.8617863895711042*GCC[14]+0.8617863895711042*GBC[14]+0.6123724356957944*GCC[9]+0.6123724356957944*GBC[9]; 
  surft1_lo[4] = -(1.927012949165105*GCC[20])+1.927012949165105*GBC[20]+1.369306393762915*GCC[17]+1.369306393762915*GBC[17]; 
  surft1_lo[5] = 0.599071547271275*GCC[31]+0.599071547271275*GBC[31]-0.8617863895711042*GCC[15]+0.8617863895711042*GBC[15]+0.6123724356957944*GCC[12]+0.6123724356957944*GBC[12]; 
  surft1_lo[6] = -(1.927012949165105*GCC[22])+1.927012949165105*GBC[22]+1.369306393762915*GCC[19]+1.369306393762915*GBC[19]; 
  surft1_lo[7] = -(1.9270129491651047*GCC[23])+1.9270129491651047*GBC[23]+1.369306393762915*GCC[21]+1.369306393762915*GBC[21]; 
  surft1_lo[12] = -(0.8617863895711042*GCC[38])+0.8617863895711042*GBC[38]+0.6123724356957944*GCC[34]+0.6123724356957944*GBC[34]; 
  surft1_lo[13] = -(0.8617863895711042*GCC[39])+0.8617863895711042*GBC[39]+0.6123724356957944*GCC[36]+0.6123724356957944*GBC[36]; 
  surft1_up[0] = dGdvx_surf_CC_vy[0]/dv1; 
  surft1_up[1] = dGdvx_surf_CC_vy[1]/dv1; 
  surft1_up[2] = dGdvx_surf_CC_vy[2]/dv1; 
  surft1_up[3] = dGdvx_surf_CC_vy[3]/dv1; 
  surft1_up[4] = dGdvx_surf_CC_vy[4]/dv1; 
  surft1_up[5] = dGdvx_surf_CC_vy[5]/dv1; 
  surft1_up[6] = dGdvx_surf_CC_vy[6]/dv1; 
  surft1_up[7] = dGdvx_surf_CC_vy[7]/dv1; 
  surft1_up[8] = dGdvx_surf_CC_vy[8]/dv1; 
  surft1_up[9] = dGdvx_surf_CC_vy[9]/dv1; 
  surft1_up[10] = dGdvx_surf_CC_vy[10]/dv1; 
  surft1_up[11] = dGdvx_surf_CC_vy[11]/dv1; 
  surft1_up[12] = dGdvx_surf_CC_vy[12]/dv1; 
  surft1_up[13] = dGdvx_surf_CC_vy[13]/dv1; 
  surft1_up[14] = dGdvx_surf_CC_vy[14]/dv1; 
  surft1_up[15] = dGdvx_surf_CC_vy[15]/dv1; 

  surft2_lo[0] = 0.34587411908091625*GCL[16]+0.34587411908091625*GCC[16]+0.49755260400283263*GCL[2]-0.49755260400283263*GCC[2]+0.3535533905932737*GCL[0]+0.3535533905932737*GCC[0]; 
  surft2_lo[1] = 0.34587411908091625*GCL[17]+0.34587411908091625*GCC[17]+0.49755260400283263*GCL[5]-0.49755260400283263*GCC[5]+0.3535533905932737*GCL[1]+0.3535533905932737*GCC[1]; 
  surft2_lo[2] = 0.34587411908091625*GCL[18]+0.34587411908091625*GCC[18]+0.49755260400283263*GCL[7]-0.49755260400283263*GCC[7]+0.3535533905932737*GCL[3]+0.3535533905932737*GCC[3]; 
  surft2_lo[3] = 0.34587411908091625*GCL[19]+0.34587411908091625*GCC[19]+0.49755260400283263*GCL[9]-0.49755260400283263*GCC[9]+0.3535533905932737*GCL[4]+0.3535533905932737*GCC[4]; 
  surft2_lo[4] = 0.34587411908091625*GCL[20]+0.34587411908091625*GCC[20]+0.49755260400283263*GCL[11]-0.49755260400283263*GCC[11]+0.3535533905932737*GCL[6]+0.3535533905932737*GCC[6]; 
  surft2_lo[5] = 0.34587411908091625*GCL[21]+0.34587411908091625*GCC[21]+0.49755260400283263*GCL[12]-0.49755260400283263*GCC[12]+0.3535533905932737*GCL[8]+0.3535533905932737*GCC[8]; 
  surft2_lo[6] = 0.34587411908091625*GCL[22]+0.34587411908091625*GCC[22]+0.49755260400283263*GCL[14]-0.49755260400283263*GCC[14]+0.3535533905932737*GCL[10]+0.3535533905932737*GCC[10]; 
  surft2_lo[7] = 0.34587411908091625*GCL[23]+0.34587411908091625*GCC[23]+0.49755260400283263*GCL[15]-0.49755260400283263*GCC[15]+0.3535533905932737*GCL[13]+0.3535533905932737*GCC[13]; 
  surft2_lo[8] = 0.49755260400283263*GCL[26]-0.49755260400283263*GCC[26]+0.3535533905932737*GCL[24]+0.3535533905932737*GCC[24]; 
  surft2_lo[9] = 0.49755260400283263*GCL[28]-0.49755260400283263*GCC[28]+0.3535533905932737*GCL[25]+0.3535533905932737*GCC[25]; 
  surft2_lo[10] = 0.49755260400283263*GCL[30]-0.49755260400283263*GCC[30]+0.3535533905932737*GCL[27]+0.3535533905932737*GCC[27]; 
  surft2_lo[11] = 0.49755260400283263*GCL[31]-0.49755260400283263*GCC[31]+0.3535533905932737*GCL[29]+0.3535533905932737*GCC[29]; 
  surft2_lo[12] = 0.49755260400283263*GCL[34]-0.49755260400283263*GCC[34]+0.3535533905932737*GCL[32]+0.3535533905932737*GCC[32]; 
  surft2_lo[13] = 0.49755260400283263*GCL[36]-0.49755260400283263*GCC[36]+0.3535533905932737*GCL[33]+0.3535533905932737*GCC[33]; 
  surft2_lo[14] = 0.49755260400283263*GCL[38]-0.49755260400283263*GCC[38]+0.3535533905932737*GCL[35]+0.3535533905932737*GCC[35]; 
  surft2_lo[15] = 0.49755260400283263*GCL[39]-0.49755260400283263*GCC[39]+0.3535533905932737*GCL[37]+0.3535533905932737*GCC[37]; 
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
  out[16] = (vol[16]+0.7071067811865475*surft1_up[8]-0.7071067811865475*surft1_lo[8])*dv1_sq*gamma_avg; 
  out[17] = (vol[17]+0.7071067811865475*surft1_up[9]-0.7071067811865475*surft1_lo[9])*dv1_sq*gamma_avg; 
  out[18] = (vol[18]+1.224744871391589*surft1_up[8]+1.224744871391589*surft1_lo[8]-2.7386127875258306*surft2_up[0]+2.7386127875258306*surft2_lo[0])*dv1_sq*gamma_avg; 
  out[19] = (vol[19]+0.7071067811865475*surft1_up[10]-0.7071067811865475*surft1_lo[10])*dv1_sq*gamma_avg; 
  out[20] = (vol[20]+1.224744871391589*surft1_up[9]+1.224744871391589*surft1_lo[9]-2.7386127875258306*surft2_up[1]+2.7386127875258306*surft2_lo[1])*dv1_sq*gamma_avg; 
  out[21] = (vol[21]+0.7071067811865475*surft1_up[11]-0.7071067811865475*surft1_lo[11])*dv1_sq*gamma_avg; 
  out[22] = (vol[22]+1.224744871391589*surft1_up[10]+1.224744871391589*surft1_lo[10]-2.7386127875258306*surft2_up[3]+2.7386127875258306*surft2_lo[3])*dv1_sq*gamma_avg; 
  out[23] = (vol[23]+1.224744871391589*surft1_up[11]+1.224744871391589*surft1_lo[11]-2.7386127875258306*surft2_up[5]+2.7386127875258306*surft2_lo[5])*dv1_sq*gamma_avg; 
  out[24] = (vol[24]-2.7386127875258306*surft2_up[2]+2.7386127875258306*surft2_lo[2]+1.5811388300841895*surft1_up[0]-1.5811388300841895*surft1_lo[0])*dv1_sq*gamma_avg; 
  out[25] = (vol[25]-2.7386127875258306*surft2_up[4]+2.7386127875258306*surft2_lo[4]+1.5811388300841898*surft1_up[1]-1.5811388300841898*surft1_lo[1])*dv1_sq*gamma_avg; 
  out[26] = (vol[26]-4.743416490252569*surft2_up[2]-4.743416490252569*surft2_lo[2]+1.5811388300841898*surft1_up[2]-1.5811388300841898*surft1_lo[2])*dv1_sq*gamma_avg; 
  out[27] = (vol[27]-2.7386127875258306*surft2_up[6]+2.7386127875258306*surft2_lo[6]+1.5811388300841898*surft1_up[3]-1.5811388300841898*surft1_lo[3])*dv1_sq*gamma_avg; 
  out[28] = (vol[28]-4.743416490252569*surft2_up[4]-4.743416490252569*surft2_lo[4]+1.5811388300841895*surft1_up[4]-1.5811388300841895*surft1_lo[4])*dv1_sq*gamma_avg; 
  out[29] = (vol[29]-2.7386127875258306*surft2_up[7]+2.7386127875258306*surft2_lo[7]+1.5811388300841895*surft1_up[5]-1.5811388300841895*surft1_lo[5])*dv1_sq*gamma_avg; 
  out[30] = (vol[30]-4.743416490252569*surft2_up[6]-4.743416490252569*surft2_lo[6]+1.5811388300841895*surft1_up[6]-1.5811388300841895*surft1_lo[6])*dv1_sq*gamma_avg; 
  out[31] = (vol[31]-4.743416490252569*surft2_up[7]-4.743416490252569*surft2_lo[7]+1.5811388300841898*surft1_up[7]-1.5811388300841898*surft1_lo[7])*dv1_sq*gamma_avg; 
  out[32] = (vol[32]+0.7071067811865475*surft1_up[12]-0.7071067811865475*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[33] = (vol[33]+0.7071067811865475*surft1_up[13]-0.7071067811865475*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[34] = (vol[34]+0.7071067811865475*surft1_up[14]-0.7071067811865475*surft1_lo[14])*dv1_sq*gamma_avg; 
  out[35] = (vol[35]-1.224744871391589*surft2_up[12]+1.224744871391589*surft2_lo[12]+1.224744871391589*surft1_up[12]+1.224744871391589*surft1_lo[12])*dv1_sq*gamma_avg; 
  out[36] = (vol[36]+0.7071067811865475*surft1_up[15]-0.7071067811865475*surft1_lo[15])*dv1_sq*gamma_avg; 
  out[37] = (vol[37]-1.224744871391589*surft2_up[13]+1.224744871391589*surft2_lo[13]+1.224744871391589*surft1_up[13]+1.224744871391589*surft1_lo[13])*dv1_sq*gamma_avg; 
  out[38] = (vol[38]+1.224744871391589*surft1_up[14]+1.224744871391589*surft1_lo[14]-2.1213203435596424*surft2_up[12]-2.1213203435596424*surft2_lo[12])*dv1_sq*gamma_avg; 
  out[39] = (vol[39]+1.224744871391589*surft1_up[15]+1.224744871391589*surft1_lo[15]-2.1213203435596424*surft2_up[13]-2.1213203435596424*surft2_lo[13])*dv1_sq*gamma_avg; 
} 
