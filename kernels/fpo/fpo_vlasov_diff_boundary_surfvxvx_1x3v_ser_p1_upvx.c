#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_vlasov_diff_boundary_surfvxvx_1x3v_ser_p1_upvx(const double *dxv, const double* diff_coeff_stencil[9], const double* f_stencil[9], double* out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_stencil[9]: 9-cell stencil of diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 


  double dv1_sq = 4.0/dxv[1]/dxv[1]; 
 
  double D_rec_lo[8], D_rec_up[8]; 
  double f_rec_lo[8], f_rec_up[8]; 
  double df_rec_lo[8], df_rec_up[8]; 
  double surft1_lo[8], surft1_up[8]; 
  double surft2_lo[8], surft2_up[8]; 

  const double* DL = &diff_coeff_stencil[0][0]; 
  const double* fL = f_stencil[0]; 
  const double* DC = &diff_coeff_stencil[1][0]; 
  const double* fC = f_stencil[1]; 

  D_rec_lo[0] = -(0.408248290463863*DL[2])+0.408248290463863*DC[2]+0.3535533905932737*DL[0]+0.3535533905932737*DC[0]; 
  D_rec_lo[1] = -(0.408248290463863*DL[5])+0.408248290463863*DC[5]+0.3535533905932737*DL[1]+0.3535533905932737*DC[1]; 
  D_rec_lo[2] = -(0.408248290463863*DL[7])+0.408248290463863*DC[7]+0.3535533905932737*DL[3]+0.3535533905932737*DC[3]; 
  D_rec_lo[3] = -(0.408248290463863*DL[9])+0.408248290463863*DC[9]+0.3535533905932737*DL[4]+0.3535533905932737*DC[4]; 
  D_rec_lo[4] = -(0.408248290463863*DL[11])+0.408248290463863*DC[11]+0.3535533905932737*DL[6]+0.3535533905932737*DC[6]; 
  D_rec_lo[5] = -(0.408248290463863*DL[12])+0.408248290463863*DC[12]+0.3535533905932737*DL[8]+0.3535533905932737*DC[8]; 
  D_rec_lo[6] = -(0.408248290463863*DL[14])+0.408248290463863*DC[14]+0.3535533905932737*DL[10]+0.3535533905932737*DC[10]; 
  D_rec_lo[7] = -(0.408248290463863*DL[15])+0.408248290463863*DC[15]+0.3535533905932737*DL[13]+0.3535533905932737*DC[13]; 
  D_rec_up[0] = 1.5811388300841895*DC[16]+1.224744871391589*DC[2]+0.7071067811865475*DC[0]; 
  D_rec_up[1] = 1.5811388300841898*DC[17]+1.224744871391589*DC[5]+0.7071067811865475*DC[1]; 
  D_rec_up[2] = 1.5811388300841898*DC[18]+1.224744871391589*DC[7]+0.7071067811865475*DC[3]; 
  D_rec_up[3] = 1.5811388300841898*DC[19]+1.224744871391589*DC[9]+0.7071067811865475*DC[4]; 
  D_rec_up[4] = 1.5811388300841895*DC[20]+1.224744871391589*DC[11]+0.7071067811865475*DC[6]; 
  D_rec_up[5] = 1.5811388300841895*DC[21]+1.224744871391589*DC[12]+0.7071067811865475*DC[8]; 
  D_rec_up[6] = 1.5811388300841895*DC[22]+1.224744871391589*DC[14]+0.7071067811865475*DC[10]; 
  D_rec_up[7] = 1.5811388300841898*DC[23]+1.224744871391589*DC[15]+0.7071067811865475*DC[13]; 

  f_rec_lo[0] = -(0.408248290463863*fL[2])+0.408248290463863*fC[2]+0.3535533905932737*fL[0]+0.3535533905932737*fC[0]; 
  f_rec_lo[1] = -(0.408248290463863*fL[5])+0.408248290463863*fC[5]+0.3535533905932737*fL[1]+0.3535533905932737*fC[1]; 
  f_rec_lo[2] = -(0.408248290463863*fL[7])+0.408248290463863*fC[7]+0.3535533905932737*fL[3]+0.3535533905932737*fC[3]; 
  f_rec_lo[3] = -(0.408248290463863*fL[9])+0.408248290463863*fC[9]+0.3535533905932737*fL[4]+0.3535533905932737*fC[4]; 
  f_rec_lo[4] = -(0.408248290463863*fL[11])+0.408248290463863*fC[11]+0.3535533905932737*fL[6]+0.3535533905932737*fC[6]; 
  f_rec_lo[5] = -(0.408248290463863*fL[12])+0.408248290463863*fC[12]+0.3535533905932737*fL[8]+0.3535533905932737*fC[8]; 
  f_rec_lo[6] = -(0.408248290463863*fL[14])+0.408248290463863*fC[14]+0.3535533905932737*fL[10]+0.3535533905932737*fC[10]; 
  f_rec_lo[7] = -(0.408248290463863*fL[15])+0.408248290463863*fC[15]+0.3535533905932737*fL[13]+0.3535533905932737*fC[13]; 
  f_rec_up[0] = 1.5811388300841895*fC[16]+1.224744871391589*fC[2]+0.7071067811865475*fC[0]; 
  f_rec_up[1] = 1.5811388300841898*fC[17]+1.224744871391589*fC[5]+0.7071067811865475*fC[1]; 
  f_rec_up[2] = 1.5811388300841898*fC[18]+1.224744871391589*fC[7]+0.7071067811865475*fC[3]; 
  f_rec_up[3] = 1.5811388300841898*fC[19]+1.224744871391589*fC[9]+0.7071067811865475*fC[4]; 
  f_rec_up[4] = 1.5811388300841895*fC[20]+1.224744871391589*fC[11]+0.7071067811865475*fC[6]; 
  f_rec_up[5] = 1.5811388300841895*fC[21]+1.224744871391589*fC[12]+0.7071067811865475*fC[8]; 
  f_rec_up[6] = 1.5811388300841895*fC[22]+1.224744871391589*fC[14]+0.7071067811865475*fC[10]; 
  f_rec_up[7] = 1.5811388300841898*fC[23]+1.224744871391589*fC[15]+0.7071067811865475*fC[13]; 

  df_rec_lo[0] = -(0.7654655446197428*fL[2])-0.7654655446197428*fC[2]+0.7954951288348656*fL[0]-0.7954951288348656*fC[0]; 
  df_rec_lo[1] = -(0.7654655446197428*fL[5])-0.7654655446197428*fC[5]+0.7954951288348656*fL[1]-0.7954951288348656*fC[1]; 
  df_rec_lo[2] = -(0.7654655446197428*fL[7])-0.7654655446197428*fC[7]+0.7954951288348656*fL[3]-0.7954951288348656*fC[3]; 
  df_rec_lo[3] = -(0.7654655446197428*fL[9])-0.7654655446197428*fC[9]+0.7954951288348656*fL[4]-0.7954951288348656*fC[4]; 
  df_rec_lo[4] = -(0.7654655446197428*fL[11])-0.7654655446197428*fC[11]+0.7954951288348656*fL[6]-0.7954951288348656*fC[6]; 
  df_rec_lo[5] = -(0.7654655446197428*fL[12])-0.7654655446197428*fC[12]+0.7954951288348656*fL[8]-0.7954951288348656*fC[8]; 
  df_rec_lo[6] = -(0.7654655446197428*fL[14])-0.7654655446197428*fC[14]+0.7954951288348656*fL[10]-0.7954951288348656*fC[10]; 
  df_rec_lo[7] = -(0.7654655446197428*fL[15])-0.7654655446197428*fC[15]+0.7954951288348656*fL[13]-0.7954951288348656*fC[13]; 

  surft1_lo[0] = 0.3535533905932737*D_rec_lo[7]*df_rec_lo[7]+0.3535533905932737*D_rec_lo[6]*df_rec_lo[6]+0.3535533905932737*D_rec_lo[5]*df_rec_lo[5]+0.3535533905932737*D_rec_lo[4]*df_rec_lo[4]+0.3535533905932737*D_rec_lo[3]*df_rec_lo[3]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[2]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[1]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[0]; 
  surft1_lo[1] = 0.3535533905932737*D_rec_lo[6]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[6]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[3]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[3]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[4]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[4]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[1]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[1]; 
  surft1_lo[2] = 0.3535533905932737*D_rec_lo[5]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[5]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[3]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[3]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[4]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[4]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[2]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[2]; 
  surft1_lo[3] = 0.3535533905932737*D_rec_lo[4]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[4]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[3]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[3]; 
  surft1_lo[4] = 0.3535533905932737*D_rec_lo[3]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[3]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[5]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[5]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[4]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[4]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[2]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[2]; 
  surft1_lo[5] = 0.3535533905932737*D_rec_lo[2]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[4]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[4]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[3]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[3]; 
  surft1_lo[6] = 0.3535533905932737*D_rec_lo[1]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[4]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[4]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[3]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[3]; 
  surft1_lo[7] = 0.3535533905932737*D_rec_lo[0]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[3]*df_rec_lo[4]+0.3535533905932737*df_rec_lo[3]*D_rec_lo[4]; 

  surft2_lo[0] = 0.3535533905932737*D_rec_lo[7]*f_rec_lo[7]+0.3535533905932737*D_rec_lo[6]*f_rec_lo[6]+0.3535533905932737*D_rec_lo[5]*f_rec_lo[5]+0.3535533905932737*D_rec_lo[4]*f_rec_lo[4]+0.3535533905932737*D_rec_lo[3]*f_rec_lo[3]+0.3535533905932737*D_rec_lo[2]*f_rec_lo[2]+0.3535533905932737*D_rec_lo[1]*f_rec_lo[1]+0.3535533905932737*D_rec_lo[0]*f_rec_lo[0]; 
  surft2_lo[1] = 0.3535533905932737*D_rec_lo[6]*f_rec_lo[7]+0.3535533905932737*f_rec_lo[6]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[3]*f_rec_lo[5]+0.3535533905932737*f_rec_lo[3]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[2]*f_rec_lo[4]+0.3535533905932737*f_rec_lo[2]*D_rec_lo[4]+0.3535533905932737*D_rec_lo[0]*f_rec_lo[1]+0.3535533905932737*f_rec_lo[0]*D_rec_lo[1]; 
  surft2_lo[2] = 0.3535533905932737*D_rec_lo[5]*f_rec_lo[7]+0.3535533905932737*f_rec_lo[5]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[3]*f_rec_lo[6]+0.3535533905932737*f_rec_lo[3]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[1]*f_rec_lo[4]+0.3535533905932737*f_rec_lo[1]*D_rec_lo[4]+0.3535533905932737*D_rec_lo[0]*f_rec_lo[2]+0.3535533905932737*f_rec_lo[0]*D_rec_lo[2]; 
  surft2_lo[3] = 0.3535533905932737*D_rec_lo[4]*f_rec_lo[7]+0.3535533905932737*f_rec_lo[4]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[2]*f_rec_lo[6]+0.3535533905932737*f_rec_lo[2]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[1]*f_rec_lo[5]+0.3535533905932737*f_rec_lo[1]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[0]*f_rec_lo[3]+0.3535533905932737*f_rec_lo[0]*D_rec_lo[3]; 
  surft2_lo[4] = 0.3535533905932737*D_rec_lo[3]*f_rec_lo[7]+0.3535533905932737*f_rec_lo[3]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[5]*f_rec_lo[6]+0.3535533905932737*f_rec_lo[5]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[0]*f_rec_lo[4]+0.3535533905932737*f_rec_lo[0]*D_rec_lo[4]+0.3535533905932737*D_rec_lo[1]*f_rec_lo[2]+0.3535533905932737*f_rec_lo[1]*D_rec_lo[2]; 
  surft2_lo[5] = 0.3535533905932737*D_rec_lo[2]*f_rec_lo[7]+0.3535533905932737*f_rec_lo[2]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[4]*f_rec_lo[6]+0.3535533905932737*f_rec_lo[4]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[0]*f_rec_lo[5]+0.3535533905932737*f_rec_lo[0]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[1]*f_rec_lo[3]+0.3535533905932737*f_rec_lo[1]*D_rec_lo[3]; 
  surft2_lo[6] = 0.3535533905932737*D_rec_lo[1]*f_rec_lo[7]+0.3535533905932737*f_rec_lo[1]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[0]*f_rec_lo[6]+0.3535533905932737*f_rec_lo[0]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[4]*f_rec_lo[5]+0.3535533905932737*f_rec_lo[4]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[2]*f_rec_lo[3]+0.3535533905932737*f_rec_lo[2]*D_rec_lo[3]; 
  surft2_lo[7] = 0.3535533905932737*D_rec_lo[0]*f_rec_lo[7]+0.3535533905932737*f_rec_lo[0]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[1]*f_rec_lo[6]+0.3535533905932737*f_rec_lo[1]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[2]*f_rec_lo[5]+0.3535533905932737*f_rec_lo[2]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[3]*f_rec_lo[4]+0.3535533905932737*f_rec_lo[3]*D_rec_lo[4]; 
  surft2_up[0] = 0.3535533905932737*D_rec_up[7]*f_rec_up[7]+0.3535533905932737*D_rec_up[6]*f_rec_up[6]+0.3535533905932737*D_rec_up[5]*f_rec_up[5]+0.3535533905932737*D_rec_up[4]*f_rec_up[4]+0.3535533905932737*D_rec_up[3]*f_rec_up[3]+0.3535533905932737*D_rec_up[2]*f_rec_up[2]+0.3535533905932737*D_rec_up[1]*f_rec_up[1]+0.3535533905932737*D_rec_up[0]*f_rec_up[0]; 
  surft2_up[1] = 0.3535533905932737*D_rec_up[6]*f_rec_up[7]+0.3535533905932737*f_rec_up[6]*D_rec_up[7]+0.3535533905932737*D_rec_up[3]*f_rec_up[5]+0.3535533905932737*f_rec_up[3]*D_rec_up[5]+0.3535533905932737*D_rec_up[2]*f_rec_up[4]+0.3535533905932737*f_rec_up[2]*D_rec_up[4]+0.3535533905932737*D_rec_up[0]*f_rec_up[1]+0.3535533905932737*f_rec_up[0]*D_rec_up[1]; 
  surft2_up[2] = 0.3535533905932737*D_rec_up[5]*f_rec_up[7]+0.3535533905932737*f_rec_up[5]*D_rec_up[7]+0.3535533905932737*D_rec_up[3]*f_rec_up[6]+0.3535533905932737*f_rec_up[3]*D_rec_up[6]+0.3535533905932737*D_rec_up[1]*f_rec_up[4]+0.3535533905932737*f_rec_up[1]*D_rec_up[4]+0.3535533905932737*D_rec_up[0]*f_rec_up[2]+0.3535533905932737*f_rec_up[0]*D_rec_up[2]; 
  surft2_up[3] = 0.3535533905932737*D_rec_up[4]*f_rec_up[7]+0.3535533905932737*f_rec_up[4]*D_rec_up[7]+0.3535533905932737*D_rec_up[2]*f_rec_up[6]+0.3535533905932737*f_rec_up[2]*D_rec_up[6]+0.3535533905932737*D_rec_up[1]*f_rec_up[5]+0.3535533905932737*f_rec_up[1]*D_rec_up[5]+0.3535533905932737*D_rec_up[0]*f_rec_up[3]+0.3535533905932737*f_rec_up[0]*D_rec_up[3]; 
  surft2_up[4] = 0.3535533905932737*D_rec_up[3]*f_rec_up[7]+0.3535533905932737*f_rec_up[3]*D_rec_up[7]+0.3535533905932737*D_rec_up[5]*f_rec_up[6]+0.3535533905932737*f_rec_up[5]*D_rec_up[6]+0.3535533905932737*D_rec_up[0]*f_rec_up[4]+0.3535533905932737*f_rec_up[0]*D_rec_up[4]+0.3535533905932737*D_rec_up[1]*f_rec_up[2]+0.3535533905932737*f_rec_up[1]*D_rec_up[2]; 
  surft2_up[5] = 0.3535533905932737*D_rec_up[2]*f_rec_up[7]+0.3535533905932737*f_rec_up[2]*D_rec_up[7]+0.3535533905932737*D_rec_up[4]*f_rec_up[6]+0.3535533905932737*f_rec_up[4]*D_rec_up[6]+0.3535533905932737*D_rec_up[0]*f_rec_up[5]+0.3535533905932737*f_rec_up[0]*D_rec_up[5]+0.3535533905932737*D_rec_up[1]*f_rec_up[3]+0.3535533905932737*f_rec_up[1]*D_rec_up[3]; 
  surft2_up[6] = 0.3535533905932737*D_rec_up[1]*f_rec_up[7]+0.3535533905932737*f_rec_up[1]*D_rec_up[7]+0.3535533905932737*D_rec_up[0]*f_rec_up[6]+0.3535533905932737*f_rec_up[0]*D_rec_up[6]+0.3535533905932737*D_rec_up[4]*f_rec_up[5]+0.3535533905932737*f_rec_up[4]*D_rec_up[5]+0.3535533905932737*D_rec_up[2]*f_rec_up[3]+0.3535533905932737*f_rec_up[2]*D_rec_up[3]; 
  surft2_up[7] = 0.3535533905932737*D_rec_up[0]*f_rec_up[7]+0.3535533905932737*f_rec_up[0]*D_rec_up[7]+0.3535533905932737*D_rec_up[1]*f_rec_up[6]+0.3535533905932737*f_rec_up[1]*D_rec_up[6]+0.3535533905932737*D_rec_up[2]*f_rec_up[5]+0.3535533905932737*f_rec_up[2]*D_rec_up[5]+0.3535533905932737*D_rec_up[3]*f_rec_up[4]+0.3535533905932737*f_rec_up[3]*D_rec_up[4]; 

  out[0] += 0.35355339059327373*surft1_up[0]*dv1_sq-0.35355339059327373*surft1_lo[0]*dv1_sq; 
  out[1] += 0.35355339059327373*surft1_up[1]*dv1_sq-0.35355339059327373*surft1_lo[1]*dv1_sq; 
  out[2] += -(0.6123724356957945*surft2_up[0]*dv1_sq)+0.6123724356957945*surft2_lo[0]*dv1_sq+0.6123724356957945*surft1_up[0]*dv1_sq+0.6123724356957945*surft1_lo[0]*dv1_sq; 
  out[3] += 0.35355339059327373*surft1_up[2]*dv1_sq-0.35355339059327373*surft1_lo[2]*dv1_sq; 
  out[4] += 0.35355339059327373*surft1_up[3]*dv1_sq-0.35355339059327373*surft1_lo[3]*dv1_sq; 
  out[5] += -(0.6123724356957945*surft2_up[1]*dv1_sq)+0.6123724356957945*surft2_lo[1]*dv1_sq+0.6123724356957945*surft1_up[1]*dv1_sq+0.6123724356957945*surft1_lo[1]*dv1_sq; 
  out[6] += 0.35355339059327373*surft1_up[4]*dv1_sq-0.35355339059327373*surft1_lo[4]*dv1_sq; 
  out[7] += -(0.6123724356957945*surft2_up[2]*dv1_sq)+0.6123724356957945*surft2_lo[2]*dv1_sq+0.6123724356957945*surft1_up[2]*dv1_sq+0.6123724356957945*surft1_lo[2]*dv1_sq; 
  out[8] += 0.35355339059327373*surft1_up[5]*dv1_sq-0.35355339059327373*surft1_lo[5]*dv1_sq; 
  out[9] += -(0.6123724356957945*surft2_up[3]*dv1_sq)+0.6123724356957945*surft2_lo[3]*dv1_sq+0.6123724356957945*surft1_up[3]*dv1_sq+0.6123724356957945*surft1_lo[3]*dv1_sq; 
  out[10] += 0.35355339059327373*surft1_up[6]*dv1_sq-0.35355339059327373*surft1_lo[6]*dv1_sq; 
  out[11] += -(0.6123724356957945*surft2_up[4]*dv1_sq)+0.6123724356957945*surft2_lo[4]*dv1_sq+0.6123724356957945*surft1_up[4]*dv1_sq+0.6123724356957945*surft1_lo[4]*dv1_sq; 
  out[12] += -(0.6123724356957945*surft2_up[5]*dv1_sq)+0.6123724356957945*surft2_lo[5]*dv1_sq+0.6123724356957945*surft1_up[5]*dv1_sq+0.6123724356957945*surft1_lo[5]*dv1_sq; 
  out[13] += 0.35355339059327373*surft1_up[7]*dv1_sq-0.35355339059327373*surft1_lo[7]*dv1_sq; 
  out[14] += -(0.6123724356957945*surft2_up[6]*dv1_sq)+0.6123724356957945*surft2_lo[6]*dv1_sq+0.6123724356957945*surft1_up[6]*dv1_sq+0.6123724356957945*surft1_lo[6]*dv1_sq; 
  out[15] += -(0.6123724356957945*surft2_up[7]*dv1_sq)+0.6123724356957945*surft2_lo[7]*dv1_sq+0.6123724356957945*surft1_up[7]*dv1_sq+0.6123724356957945*surft1_lo[7]*dv1_sq; 
  out[16] += -(2.3717082451262845*surft2_up[0]*dv1_sq)-2.3717082451262845*surft2_lo[0]*dv1_sq+0.7905694150420948*surft1_up[0]*dv1_sq-0.7905694150420948*surft1_lo[0]*dv1_sq; 
  out[17] += -(2.3717082451262845*surft2_up[1]*dv1_sq)-2.3717082451262845*surft2_lo[1]*dv1_sq+0.7905694150420949*surft1_up[1]*dv1_sq-0.7905694150420949*surft1_lo[1]*dv1_sq; 
  out[18] += -(2.3717082451262845*surft2_up[2]*dv1_sq)-2.3717082451262845*surft2_lo[2]*dv1_sq+0.7905694150420949*surft1_up[2]*dv1_sq-0.7905694150420949*surft1_lo[2]*dv1_sq; 
  out[19] += -(2.3717082451262845*surft2_up[3]*dv1_sq)-2.3717082451262845*surft2_lo[3]*dv1_sq+0.7905694150420949*surft1_up[3]*dv1_sq-0.7905694150420949*surft1_lo[3]*dv1_sq; 
  out[20] += -(2.3717082451262845*surft2_up[4]*dv1_sq)-2.3717082451262845*surft2_lo[4]*dv1_sq+0.7905694150420948*surft1_up[4]*dv1_sq-0.7905694150420948*surft1_lo[4]*dv1_sq; 
  out[21] += -(2.3717082451262845*surft2_up[5]*dv1_sq)-2.3717082451262845*surft2_lo[5]*dv1_sq+0.7905694150420948*surft1_up[5]*dv1_sq-0.7905694150420948*surft1_lo[5]*dv1_sq; 
  out[22] += -(2.3717082451262845*surft2_up[6]*dv1_sq)-2.3717082451262845*surft2_lo[6]*dv1_sq+0.7905694150420948*surft1_up[6]*dv1_sq-0.7905694150420948*surft1_lo[6]*dv1_sq; 
  out[23] += -(2.3717082451262845*surft2_up[7]*dv1_sq)-2.3717082451262845*surft2_lo[7]*dv1_sq+0.7905694150420949*surft1_up[7]*dv1_sq-0.7905694150420949*surft1_lo[7]*dv1_sq; 
} 