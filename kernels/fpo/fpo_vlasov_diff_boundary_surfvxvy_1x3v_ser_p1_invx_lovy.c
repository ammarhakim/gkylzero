#include <gkyl_fpo_vlasov_kernels.h> 
 
GKYL_CU_DH void fpo_vlasov_diff_boundary_surfvxvy_1x3v_ser_p1_invx_lovy(const double *dxv, const double* diff_coeff_stencil[9], const double* f_stencil[9], double* out) { 
  // dxv[NDIM]: Cell spacing in each direction. 
  // diff_coeff_stencil[3]: 3-cell stencil of diffusion tensor. 
  // f_stencil[9]: 9-cell stencil of distribution function. 
  // out: Incremented output. 


  double dv1_sq = 4.0/dxv[1]/dxv[2]; 
 
  double D_rec_lo[8], D_rec_up[8]; 
  double f_rec_lo[8], f_rec_up[8]; 
  double df_rec_lo[8], df_rec_up[8]; 
  double surft1_lo[8], surft1_up[8]; 
  double surft2_lo[8], surft2_up[8]; 

  const double* DCL = &diff_coeff_stencil[0][120]; 
  const double* fCL = f_stencil[0]; 
  const double* DTL = &diff_coeff_stencil[1][120]; 
  const double* fTL = f_stencil[1]; 
  const double* DCC = &diff_coeff_stencil[2][120]; 
  const double* fCC = f_stencil[2]; 
  const double* DTC = &diff_coeff_stencil[3][120]; 
  const double* fTC = f_stencil[3]; 
  const double* DCR = &diff_coeff_stencil[4][120]; 
  const double* fCR = f_stencil[4]; 
  const double* DTR = &diff_coeff_stencil[5][120]; 
  const double* fTR = f_stencil[5]; 

  D_rec_lo[0] = 0.408248290463863*DCL[2]-0.408248290463863*DCC[2]+0.3535533905932737*DCL[0]+0.3535533905932737*DCC[0]; 
  D_rec_lo[1] = 0.408248290463863*DCL[5]-0.408248290463863*DCC[5]+0.3535533905932737*DCL[1]+0.3535533905932737*DCC[1]; 
  D_rec_lo[2] = 0.408248290463863*DCL[7]-0.408248290463863*DCC[7]+0.3535533905932737*DCL[3]+0.3535533905932737*DCC[3]; 
  D_rec_lo[3] = 0.408248290463863*DCL[9]-0.408248290463863*DCC[9]+0.3535533905932737*DCL[4]+0.3535533905932737*DCC[4]; 
  D_rec_lo[4] = 0.408248290463863*DCL[11]-0.408248290463863*DCC[11]+0.3535533905932737*DCL[6]+0.3535533905932737*DCC[6]; 
  D_rec_lo[5] = 0.408248290463863*DCL[12]-0.408248290463863*DCC[12]+0.3535533905932737*DCL[8]+0.3535533905932737*DCC[8]; 
  D_rec_lo[6] = 0.408248290463863*DCL[14]-0.408248290463863*DCC[14]+0.3535533905932737*DCL[10]+0.3535533905932737*DCC[10]; 
  D_rec_lo[7] = 0.408248290463863*DCL[15]-0.408248290463863*DCC[15]+0.3535533905932737*DCL[13]+0.3535533905932737*DCC[13]; 
  D_rec_up[0] = -(0.408248290463863*DCR[2])+0.408248290463863*DCC[2]+0.3535533905932737*DCR[0]+0.3535533905932737*DCC[0]; 
  D_rec_up[1] = -(0.408248290463863*DCR[5])+0.408248290463863*DCC[5]+0.3535533905932737*DCR[1]+0.3535533905932737*DCC[1]; 
  D_rec_up[2] = -(0.408248290463863*DCR[7])+0.408248290463863*DCC[7]+0.3535533905932737*DCR[3]+0.3535533905932737*DCC[3]; 
  D_rec_up[3] = -(0.408248290463863*DCR[9])+0.408248290463863*DCC[9]+0.3535533905932737*DCR[4]+0.3535533905932737*DCC[4]; 
  D_rec_up[4] = -(0.408248290463863*DCR[11])+0.408248290463863*DCC[11]+0.3535533905932737*DCR[6]+0.3535533905932737*DCC[6]; 
  D_rec_up[5] = -(0.408248290463863*DCR[12])+0.408248290463863*DCC[12]+0.3535533905932737*DCR[8]+0.3535533905932737*DCC[8]; 
  D_rec_up[6] = -(0.408248290463863*DCR[14])+0.408248290463863*DCC[14]+0.3535533905932737*DCR[10]+0.3535533905932737*DCC[10]; 
  D_rec_up[7] = -(0.408248290463863*DCR[15])+0.408248290463863*DCC[15]+0.3535533905932737*DCR[13]+0.3535533905932737*DCC[13]; 

  f_rec_lo[0] = 0.408248290463863*fCL[2]-0.408248290463863*fCC[2]+0.3535533905932737*fCL[0]+0.3535533905932737*fCC[0]; 
  f_rec_lo[1] = 0.408248290463863*fCL[5]-0.408248290463863*fCC[5]+0.3535533905932737*fCL[1]+0.3535533905932737*fCC[1]; 
  f_rec_lo[2] = 0.408248290463863*fCL[7]-0.408248290463863*fCC[7]+0.3535533905932737*fCL[3]+0.3535533905932737*fCC[3]; 
  f_rec_lo[3] = 0.408248290463863*fCL[9]-0.408248290463863*fCC[9]+0.3535533905932737*fCL[4]+0.3535533905932737*fCC[4]; 
  f_rec_lo[4] = 0.408248290463863*fCL[11]-0.408248290463863*fCC[11]+0.3535533905932737*fCL[6]+0.3535533905932737*fCC[6]; 
  f_rec_lo[5] = 0.408248290463863*fCL[12]-0.408248290463863*fCC[12]+0.3535533905932737*fCL[8]+0.3535533905932737*fCC[8]; 
  f_rec_lo[6] = 0.408248290463863*fCL[14]-0.408248290463863*fCC[14]+0.3535533905932737*fCL[10]+0.3535533905932737*fCC[10]; 
  f_rec_lo[7] = 0.408248290463863*fCL[15]-0.408248290463863*fCC[15]+0.3535533905932737*fCL[13]+0.3535533905932737*fCC[13]; 
  f_rec_up[0] = -(0.408248290463863*fCR[2])+0.408248290463863*fCC[2]+0.3535533905932737*fCR[0]+0.3535533905932737*fCC[0]; 
  f_rec_up[1] = -(0.408248290463863*fCR[5])+0.408248290463863*fCC[5]+0.3535533905932737*fCR[1]+0.3535533905932737*fCC[1]; 
  f_rec_up[2] = -(0.408248290463863*fCR[7])+0.408248290463863*fCC[7]+0.3535533905932737*fCR[3]+0.3535533905932737*fCC[3]; 
  f_rec_up[3] = -(0.408248290463863*fCR[9])+0.408248290463863*fCC[9]+0.3535533905932737*fCR[4]+0.3535533905932737*fCC[4]; 
  f_rec_up[4] = -(0.408248290463863*fCR[11])+0.408248290463863*fCC[11]+0.3535533905932737*fCR[6]+0.3535533905932737*fCC[6]; 
  f_rec_up[5] = -(0.408248290463863*fCR[12])+0.408248290463863*fCC[12]+0.3535533905932737*fCR[8]+0.3535533905932737*fCC[8]; 
  f_rec_up[6] = -(0.408248290463863*fCR[14])+0.408248290463863*fCC[14]+0.3535533905932737*fCR[10]+0.3535533905932737*fCC[10]; 
  f_rec_up[7] = -(0.408248290463863*fCR[15])+0.408248290463863*fCC[15]+0.3535533905932737*fCR[13]+0.3535533905932737*fCC[13]; 

  df_rec_lo[0] = -(0.11785113019775789*fTL[7])+0.11785113019775789*fTC[7]+0.11785113019775789*fCL[7]-0.11785113019775789*fCC[7]-0.10206207261596573*fTL[3]-0.10206207261596573*fTC[3]+0.10206207261596573*fCL[3]+0.10206207261596573*fCC[3]-0.8660254037844386*f_rec_lo[2]+0.10206207261596573*fTL[2]-0.10206207261596573*fTC[2]+0.10206207261596573*fCL[2]-0.10206207261596573*fCC[2]-0.5*f_rec_lo[0]+0.0883883476483184*fTL[0]+0.0883883476483184*fTC[0]+0.0883883476483184*fCL[0]+0.0883883476483184*fCC[0]; 
  df_rec_lo[1] = -(0.11785113019775789*fTL[11])+0.11785113019775789*fTC[11]+0.11785113019775789*fCL[11]-0.11785113019775789*fCC[11]-0.10206207261596573*fTL[6]-0.10206207261596573*fTC[6]+0.10206207261596573*fCL[6]+0.10206207261596573*fCC[6]+0.10206207261596573*fTL[5]-0.10206207261596573*fTC[5]+0.10206207261596573*fCL[5]-0.10206207261596573*fCC[5]-0.8660254037844386*f_rec_lo[4]-0.5*f_rec_lo[1]+0.0883883476483184*fTL[1]+0.0883883476483184*fTC[1]+0.0883883476483184*fCL[1]+0.0883883476483184*fCC[1]; 
  df_rec_lo[2] = -(0.20412414523193148*fTL[7])+0.20412414523193148*fTC[7]+0.20412414523193148*fCL[7]-0.20412414523193148*fCC[7]-0.1767766952966368*fTL[3]-0.1767766952966368*fTC[3]+0.1767766952966368*fCL[3]+0.1767766952966368*fCC[3]+1.5*f_rec_lo[2]+0.1767766952966368*fTL[2]-0.1767766952966368*fTC[2]-0.5303300858899105*fCL[2]+0.5303300858899105*fCC[2]+0.8660254037844386*f_rec_lo[0]+0.15309310892394856*fTL[0]+0.15309310892394856*fTC[0]-0.45927932677184563*fCL[0]-0.45927932677184563*fCC[0]; 
  df_rec_lo[3] = -(0.11785113019775789*fTL[14])+0.11785113019775789*fTC[14]+0.11785113019775789*fCL[14]-0.11785113019775789*fCC[14]-0.10206207261596573*fTL[10]-0.10206207261596573*fTC[10]+0.10206207261596573*fCL[10]+0.10206207261596573*fCC[10]+0.10206207261596573*fTL[9]-0.10206207261596573*fTC[9]+0.10206207261596573*fCL[9]-0.10206207261596573*fCC[9]-0.8660254037844386*f_rec_lo[6]+0.0883883476483184*fTL[4]+0.0883883476483184*fTC[4]+0.0883883476483184*fCL[4]+0.0883883476483184*fCC[4]-0.5*f_rec_lo[3]; 
  df_rec_lo[4] = -(0.20412414523193148*fTL[11])+0.20412414523193148*fTC[11]+0.20412414523193148*fCL[11]-0.20412414523193148*fCC[11]-0.1767766952966368*fTL[6]-0.1767766952966368*fTC[6]+0.1767766952966368*fCL[6]+0.1767766952966368*fCC[6]+0.1767766952966368*fTL[5]-0.1767766952966368*fTC[5]-0.5303300858899105*fCL[5]+0.5303300858899105*fCC[5]+1.5*f_rec_lo[4]+0.8660254037844386*f_rec_lo[1]+0.15309310892394856*fTL[1]+0.15309310892394856*fTC[1]-0.45927932677184563*fCL[1]-0.45927932677184563*fCC[1]; 
  df_rec_lo[5] = -(0.11785113019775789*fTL[15])+0.11785113019775789*fTC[15]+0.11785113019775789*fCL[15]-0.11785113019775789*fCC[15]-0.10206207261596573*fTL[13]-0.10206207261596573*fTC[13]+0.10206207261596573*fCL[13]+0.10206207261596573*fCC[13]+0.10206207261596573*fTL[12]-0.10206207261596573*fTC[12]+0.10206207261596573*fCL[12]-0.10206207261596573*fCC[12]+0.0883883476483184*fTL[8]+0.0883883476483184*fTC[8]+0.0883883476483184*fCL[8]+0.0883883476483184*fCC[8]-0.8660254037844386*f_rec_lo[7]-0.5*f_rec_lo[5]; 
  df_rec_lo[6] = -(0.20412414523193148*fTL[14])+0.20412414523193148*fTC[14]+0.20412414523193148*fCL[14]-0.20412414523193148*fCC[14]-0.1767766952966368*fTL[10]-0.1767766952966368*fTC[10]+0.1767766952966368*fCL[10]+0.1767766952966368*fCC[10]+0.1767766952966368*fTL[9]-0.1767766952966368*fTC[9]-0.5303300858899105*fCL[9]+0.5303300858899105*fCC[9]+1.5*f_rec_lo[6]+0.15309310892394856*fTL[4]+0.15309310892394856*fTC[4]-0.45927932677184563*fCL[4]-0.45927932677184563*fCC[4]+0.8660254037844386*f_rec_lo[3]; 
  df_rec_lo[7] = -(0.20412414523193148*fTL[15])+0.20412414523193148*fTC[15]+0.20412414523193148*fCL[15]-0.20412414523193148*fCC[15]-0.1767766952966368*fTL[13]-0.1767766952966368*fTC[13]+0.1767766952966368*fCL[13]+0.1767766952966368*fCC[13]+0.1767766952966368*fTL[12]-0.1767766952966368*fTC[12]-0.5303300858899105*fCL[12]+0.5303300858899105*fCC[12]+0.15309310892394856*fTL[8]+0.15309310892394856*fTC[8]-0.45927932677184563*fCL[8]-0.45927932677184563*fCC[8]+1.5*f_rec_lo[7]+0.8660254037844386*f_rec_lo[5]; 
  df_rec_up[0] = 0.11785113019775789*fTR[7]-0.11785113019775789*fTC[7]-0.11785113019775789*fCR[7]+0.11785113019775789*fCC[7]-0.10206207261596573*fTR[3]-0.10206207261596573*fTC[3]+0.10206207261596573*fCR[3]+0.10206207261596573*fCC[3]-0.8660254037844386*f_rec_up[2]-0.10206207261596573*fTR[2]+0.10206207261596573*fTC[2]-0.10206207261596573*fCR[2]+0.10206207261596573*fCC[2]-0.5*f_rec_up[0]+0.0883883476483184*fTR[0]+0.0883883476483184*fTC[0]+0.0883883476483184*fCR[0]+0.0883883476483184*fCC[0]; 
  df_rec_up[1] = 0.11785113019775789*fTR[11]-0.11785113019775789*fTC[11]-0.11785113019775789*fCR[11]+0.11785113019775789*fCC[11]-0.10206207261596573*fTR[6]-0.10206207261596573*fTC[6]+0.10206207261596573*fCR[6]+0.10206207261596573*fCC[6]-0.10206207261596573*fTR[5]+0.10206207261596573*fTC[5]-0.10206207261596573*fCR[5]+0.10206207261596573*fCC[5]-0.8660254037844386*f_rec_up[4]-0.5*f_rec_up[1]+0.0883883476483184*fTR[1]+0.0883883476483184*fTC[1]+0.0883883476483184*fCR[1]+0.0883883476483184*fCC[1]; 
  df_rec_up[2] = 0.20412414523193148*fTR[7]-0.20412414523193148*fTC[7]-0.20412414523193148*fCR[7]+0.20412414523193148*fCC[7]-0.1767766952966368*fTR[3]-0.1767766952966368*fTC[3]+0.1767766952966368*fCR[3]+0.1767766952966368*fCC[3]+1.5*f_rec_up[2]-0.1767766952966368*fTR[2]+0.1767766952966368*fTC[2]+0.5303300858899105*fCR[2]-0.5303300858899105*fCC[2]+0.8660254037844386*f_rec_up[0]+0.15309310892394856*fTR[0]+0.15309310892394856*fTC[0]-0.45927932677184563*fCR[0]-0.45927932677184563*fCC[0]; 
  df_rec_up[3] = 0.11785113019775789*fTR[14]-0.11785113019775789*fTC[14]-0.11785113019775789*fCR[14]+0.11785113019775789*fCC[14]-0.10206207261596573*fTR[10]-0.10206207261596573*fTC[10]+0.10206207261596573*fCR[10]+0.10206207261596573*fCC[10]-0.10206207261596573*fTR[9]+0.10206207261596573*fTC[9]-0.10206207261596573*fCR[9]+0.10206207261596573*fCC[9]-0.8660254037844386*f_rec_up[6]+0.0883883476483184*fTR[4]+0.0883883476483184*fTC[4]+0.0883883476483184*fCR[4]+0.0883883476483184*fCC[4]-0.5*f_rec_up[3]; 
  df_rec_up[4] = 0.20412414523193148*fTR[11]-0.20412414523193148*fTC[11]-0.20412414523193148*fCR[11]+0.20412414523193148*fCC[11]-0.1767766952966368*fTR[6]-0.1767766952966368*fTC[6]+0.1767766952966368*fCR[6]+0.1767766952966368*fCC[6]-0.1767766952966368*fTR[5]+0.1767766952966368*fTC[5]+0.5303300858899105*fCR[5]-0.5303300858899105*fCC[5]+1.5*f_rec_up[4]+0.8660254037844386*f_rec_up[1]+0.15309310892394856*fTR[1]+0.15309310892394856*fTC[1]-0.45927932677184563*fCR[1]-0.45927932677184563*fCC[1]; 
  df_rec_up[5] = 0.11785113019775789*fTR[15]-0.11785113019775789*fTC[15]-0.11785113019775789*fCR[15]+0.11785113019775789*fCC[15]-0.10206207261596573*fTR[13]-0.10206207261596573*fTC[13]+0.10206207261596573*fCR[13]+0.10206207261596573*fCC[13]-0.10206207261596573*fTR[12]+0.10206207261596573*fTC[12]-0.10206207261596573*fCR[12]+0.10206207261596573*fCC[12]+0.0883883476483184*fTR[8]+0.0883883476483184*fTC[8]+0.0883883476483184*fCR[8]+0.0883883476483184*fCC[8]-0.8660254037844386*f_rec_up[7]-0.5*f_rec_up[5]; 
  df_rec_up[6] = 0.20412414523193148*fTR[14]-0.20412414523193148*fTC[14]-0.20412414523193148*fCR[14]+0.20412414523193148*fCC[14]-0.1767766952966368*fTR[10]-0.1767766952966368*fTC[10]+0.1767766952966368*fCR[10]+0.1767766952966368*fCC[10]-0.1767766952966368*fTR[9]+0.1767766952966368*fTC[9]+0.5303300858899105*fCR[9]-0.5303300858899105*fCC[9]+1.5*f_rec_up[6]+0.15309310892394856*fTR[4]+0.15309310892394856*fTC[4]-0.45927932677184563*fCR[4]-0.45927932677184563*fCC[4]+0.8660254037844386*f_rec_up[3]; 
  df_rec_up[7] = 0.20412414523193148*fTR[15]-0.20412414523193148*fTC[15]-0.20412414523193148*fCR[15]+0.20412414523193148*fCC[15]-0.1767766952966368*fTR[13]-0.1767766952966368*fTC[13]+0.1767766952966368*fCR[13]+0.1767766952966368*fCC[13]-0.1767766952966368*fTR[12]+0.1767766952966368*fTC[12]+0.5303300858899105*fCR[12]-0.5303300858899105*fCC[12]+0.15309310892394856*fTR[8]+0.15309310892394856*fTC[8]-0.45927932677184563*fCR[8]-0.45927932677184563*fCC[8]+1.5*f_rec_up[7]+0.8660254037844386*f_rec_up[5]; 

  surft1_lo[0] = 0.3535533905932737*D_rec_lo[7]*df_rec_lo[7]+0.3535533905932737*D_rec_lo[6]*df_rec_lo[6]+0.3535533905932737*D_rec_lo[5]*df_rec_lo[5]+0.3535533905932737*D_rec_lo[4]*df_rec_lo[4]+0.3535533905932737*D_rec_lo[3]*df_rec_lo[3]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[2]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[1]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[0]; 
  surft1_lo[1] = 0.3535533905932737*D_rec_lo[6]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[6]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[3]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[3]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[4]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[4]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[1]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[1]; 
  surft1_lo[2] = 0.3535533905932737*D_rec_lo[5]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[5]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[3]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[3]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[4]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[4]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[2]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[2]; 
  surft1_lo[3] = 0.3535533905932737*D_rec_lo[4]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[4]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[3]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[3]; 
  surft1_lo[4] = 0.3535533905932737*D_rec_lo[3]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[3]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[5]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[5]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[4]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[4]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[2]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[2]; 
  surft1_lo[5] = 0.3535533905932737*D_rec_lo[2]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[4]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[4]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[3]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[3]; 
  surft1_lo[6] = 0.3535533905932737*D_rec_lo[1]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[0]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[4]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[4]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[3]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[3]; 
  surft1_lo[7] = 0.3535533905932737*D_rec_lo[0]*df_rec_lo[7]+0.3535533905932737*df_rec_lo[0]*D_rec_lo[7]+0.3535533905932737*D_rec_lo[1]*df_rec_lo[6]+0.3535533905932737*df_rec_lo[1]*D_rec_lo[6]+0.3535533905932737*D_rec_lo[2]*df_rec_lo[5]+0.3535533905932737*df_rec_lo[2]*D_rec_lo[5]+0.3535533905932737*D_rec_lo[3]*df_rec_lo[4]+0.3535533905932737*df_rec_lo[3]*D_rec_lo[4]; 
  surft1_up[0] = 0.3535533905932737*D_rec_up[7]*df_rec_up[7]+0.3535533905932737*D_rec_up[6]*df_rec_up[6]+0.3535533905932737*D_rec_up[5]*df_rec_up[5]+0.3535533905932737*D_rec_up[4]*df_rec_up[4]+0.3535533905932737*D_rec_up[3]*df_rec_up[3]+0.3535533905932737*D_rec_up[2]*df_rec_up[2]+0.3535533905932737*D_rec_up[1]*df_rec_up[1]+0.3535533905932737*D_rec_up[0]*df_rec_up[0]; 
  surft1_up[1] = 0.3535533905932737*D_rec_up[6]*df_rec_up[7]+0.3535533905932737*df_rec_up[6]*D_rec_up[7]+0.3535533905932737*D_rec_up[3]*df_rec_up[5]+0.3535533905932737*df_rec_up[3]*D_rec_up[5]+0.3535533905932737*D_rec_up[2]*df_rec_up[4]+0.3535533905932737*df_rec_up[2]*D_rec_up[4]+0.3535533905932737*D_rec_up[0]*df_rec_up[1]+0.3535533905932737*df_rec_up[0]*D_rec_up[1]; 
  surft1_up[2] = 0.3535533905932737*D_rec_up[5]*df_rec_up[7]+0.3535533905932737*df_rec_up[5]*D_rec_up[7]+0.3535533905932737*D_rec_up[3]*df_rec_up[6]+0.3535533905932737*df_rec_up[3]*D_rec_up[6]+0.3535533905932737*D_rec_up[1]*df_rec_up[4]+0.3535533905932737*df_rec_up[1]*D_rec_up[4]+0.3535533905932737*D_rec_up[0]*df_rec_up[2]+0.3535533905932737*df_rec_up[0]*D_rec_up[2]; 
  surft1_up[3] = 0.3535533905932737*D_rec_up[4]*df_rec_up[7]+0.3535533905932737*df_rec_up[4]*D_rec_up[7]+0.3535533905932737*D_rec_up[2]*df_rec_up[6]+0.3535533905932737*df_rec_up[2]*D_rec_up[6]+0.3535533905932737*D_rec_up[1]*df_rec_up[5]+0.3535533905932737*df_rec_up[1]*D_rec_up[5]+0.3535533905932737*D_rec_up[0]*df_rec_up[3]+0.3535533905932737*df_rec_up[0]*D_rec_up[3]; 
  surft1_up[4] = 0.3535533905932737*D_rec_up[3]*df_rec_up[7]+0.3535533905932737*df_rec_up[3]*D_rec_up[7]+0.3535533905932737*D_rec_up[5]*df_rec_up[6]+0.3535533905932737*df_rec_up[5]*D_rec_up[6]+0.3535533905932737*D_rec_up[0]*df_rec_up[4]+0.3535533905932737*df_rec_up[0]*D_rec_up[4]+0.3535533905932737*D_rec_up[1]*df_rec_up[2]+0.3535533905932737*df_rec_up[1]*D_rec_up[2]; 
  surft1_up[5] = 0.3535533905932737*D_rec_up[2]*df_rec_up[7]+0.3535533905932737*df_rec_up[2]*D_rec_up[7]+0.3535533905932737*D_rec_up[4]*df_rec_up[6]+0.3535533905932737*df_rec_up[4]*D_rec_up[6]+0.3535533905932737*D_rec_up[0]*df_rec_up[5]+0.3535533905932737*df_rec_up[0]*D_rec_up[5]+0.3535533905932737*D_rec_up[1]*df_rec_up[3]+0.3535533905932737*df_rec_up[1]*D_rec_up[3]; 
  surft1_up[6] = 0.3535533905932737*D_rec_up[1]*df_rec_up[7]+0.3535533905932737*df_rec_up[1]*D_rec_up[7]+0.3535533905932737*D_rec_up[0]*df_rec_up[6]+0.3535533905932737*df_rec_up[0]*D_rec_up[6]+0.3535533905932737*D_rec_up[4]*df_rec_up[5]+0.3535533905932737*df_rec_up[4]*D_rec_up[5]+0.3535533905932737*D_rec_up[2]*df_rec_up[3]+0.3535533905932737*df_rec_up[2]*D_rec_up[3]; 
  surft1_up[7] = 0.3535533905932737*D_rec_up[0]*df_rec_up[7]+0.3535533905932737*df_rec_up[0]*D_rec_up[7]+0.3535533905932737*D_rec_up[1]*df_rec_up[6]+0.3535533905932737*df_rec_up[1]*D_rec_up[6]+0.3535533905932737*D_rec_up[2]*df_rec_up[5]+0.3535533905932737*df_rec_up[2]*D_rec_up[5]+0.3535533905932737*D_rec_up[3]*df_rec_up[4]+0.3535533905932737*df_rec_up[3]*D_rec_up[4]; 

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