#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_no_by_surfx_3x2v_ser_p1(const double *w, const double *dxv, 
  const double *alpha_surf_l, const double *alpha_surf_r, 
  const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
  const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // alpha_surf_l: Surface expansion of phase space flux on the left.
  // alpha_surf_r: Surface expansion of phase space flux on the right.
  // sgn_alpha_surf_l: sign(alpha_surf_l) at quadrature points.
  // sgn_alpha_surf_r: sign(alpha_surf_r) at quadrature points.
  // const_sgn_alpha_l: Boolean array true if sign(alpha_surf_l) is only one sign, either +1 or -1.
  // const_sgn_alpha_r: Boolean array true if sign(alpha_surf_r) is only one sign, either +1 or -1.
  // fl,fc,fr: distribution function in left, center and right cells.
  // out: output increment in center cell.

  double rdx2 = 2.0/dxv[0];

  const double *alphaL = &alpha_surf_l[0];
  const double *alphaR = &alpha_surf_r[0];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_l[0];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_r[0];
  const int *const_sgn_alphaL = &const_sgn_alpha_l[0];
  const int *const_sgn_alphaR = &const_sgn_alpha_r[0];

  double fUpL[24] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = 1.224744871391589*fl[1]+0.7071067811865475*fl[0]; 
  fUpL[1] = 1.224744871391589*fl[6]+0.7071067811865475*fl[2]; 
  fUpL[2] = 1.224744871391589*fl[7]+0.7071067811865475*fl[3]; 
  fUpL[3] = 1.224744871391589*fl[9]+0.7071067811865475*fl[4]; 
  fUpL[4] = 1.224744871391589*fl[12]+0.7071067811865475*fl[5]; 
  fUpL[5] = 1.224744871391589*fl[16]+0.7071067811865475*fl[8]; 
  fUpL[6] = 1.224744871391589*fl[17]+0.7071067811865475*fl[10]; 
  fUpL[7] = 1.224744871391589*fl[18]+0.7071067811865475*fl[11]; 
  fUpL[8] = 1.224744871391589*fl[20]+0.7071067811865475*fl[13]; 
  fUpL[9] = 1.224744871391589*fl[21]+0.7071067811865475*fl[14]; 
  fUpL[10] = 1.224744871391589*fl[23]+0.7071067811865475*fl[15]; 
  fUpL[11] = 1.224744871391589*fl[26]+0.7071067811865475*fl[19]; 
  fUpL[12] = 1.224744871391589*fl[27]+0.7071067811865475*fl[22]; 
  fUpL[13] = 1.224744871391589*fl[28]+0.7071067811865475*fl[24]; 
  fUpL[14] = 1.224744871391589*fl[29]+0.7071067811865475*fl[25]; 
  fUpL[15] = 1.224744871391589*fl[31]+0.7071067811865475*fl[30]; 
  fUpL[16] = 1.224744871391589*fl[33]+0.7071067811865475*fl[32]; 
  fUpL[17] = 1.224744871391589*fl[37]+0.7071067811865475*fl[34]; 
  fUpL[18] = 1.224744871391589*fl[38]+0.7071067811865475*fl[35]; 
  fUpL[19] = 1.224744871391589*fl[40]+0.7071067811865475*fl[36]; 
  fUpL[20] = 1.224744871391589*fl[43]+0.7071067811865475*fl[39]; 
  fUpL[21] = 1.224744871391589*fl[44]+0.7071067811865475*fl[41]; 
  fUpL[22] = 1.224744871391589*fl[45]+0.7071067811865475*fl[42]; 
  fUpL[23] = 1.224744871391589*fl[47]+0.7071067811865475*fl[46]; 
    } else { 
  fUpL[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[1]; 
  fUpL[1] = 0.7071067811865475*fc[2]-1.224744871391589*fc[6]; 
  fUpL[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[7]; 
  fUpL[3] = 0.7071067811865475*fc[4]-1.224744871391589*fc[9]; 
  fUpL[4] = 0.7071067811865475*fc[5]-1.224744871391589*fc[12]; 
  fUpL[5] = 0.7071067811865475*fc[8]-1.224744871391589*fc[16]; 
  fUpL[6] = 0.7071067811865475*fc[10]-1.224744871391589*fc[17]; 
  fUpL[7] = 0.7071067811865475*fc[11]-1.224744871391589*fc[18]; 
  fUpL[8] = 0.7071067811865475*fc[13]-1.224744871391589*fc[20]; 
  fUpL[9] = 0.7071067811865475*fc[14]-1.224744871391589*fc[21]; 
  fUpL[10] = 0.7071067811865475*fc[15]-1.224744871391589*fc[23]; 
  fUpL[11] = 0.7071067811865475*fc[19]-1.224744871391589*fc[26]; 
  fUpL[12] = 0.7071067811865475*fc[22]-1.224744871391589*fc[27]; 
  fUpL[13] = 0.7071067811865475*fc[24]-1.224744871391589*fc[28]; 
  fUpL[14] = 0.7071067811865475*fc[25]-1.224744871391589*fc[29]; 
  fUpL[15] = 0.7071067811865475*fc[30]-1.224744871391589*fc[31]; 
  fUpL[16] = 0.7071067811865475*fc[32]-1.224744871391589*fc[33]; 
  fUpL[17] = 0.7071067811865475*fc[34]-1.224744871391589*fc[37]; 
  fUpL[18] = 0.7071067811865475*fc[35]-1.224744871391589*fc[38]; 
  fUpL[19] = 0.7071067811865475*fc[36]-1.224744871391589*fc[40]; 
  fUpL[20] = 0.7071067811865475*fc[39]-1.224744871391589*fc[43]; 
  fUpL[21] = 0.7071067811865475*fc[41]-1.224744871391589*fc[44]; 
  fUpL[22] = 0.7071067811865475*fc[42]-1.224744871391589*fc[45]; 
  fUpL[23] = 0.7071067811865475*fc[46]-1.224744871391589*fc[47]; 
    } 
  } else { 
  double f_lr[24] = {0.};
  double f_cl[24] = {0.};
  double sgn_alphaUpL[2] = {0.};
  sgn_alphaUpL[0] = 0.7071067811865475*sgn_alpha_surfL[1]+0.7071067811865475*sgn_alpha_surfL[0]; 
  sgn_alphaUpL[1] = 0.7071067811865475*sgn_alpha_surfL[1]-0.7071067811865475*sgn_alpha_surfL[0]; 

  f_lr[0] = 1.224744871391589*fl[1]+0.7071067811865475*fl[0]; 
  f_lr[1] = 1.224744871391589*fl[6]+0.7071067811865475*fl[2]; 
  f_lr[2] = 1.224744871391589*fl[7]+0.7071067811865475*fl[3]; 
  f_lr[3] = 1.224744871391589*fl[9]+0.7071067811865475*fl[4]; 
  f_lr[4] = 1.224744871391589*fl[12]+0.7071067811865475*fl[5]; 
  f_lr[5] = 1.224744871391589*fl[16]+0.7071067811865475*fl[8]; 
  f_lr[6] = 1.224744871391589*fl[17]+0.7071067811865475*fl[10]; 
  f_lr[7] = 1.224744871391589*fl[18]+0.7071067811865475*fl[11]; 
  f_lr[8] = 1.224744871391589*fl[20]+0.7071067811865475*fl[13]; 
  f_lr[9] = 1.224744871391589*fl[21]+0.7071067811865475*fl[14]; 
  f_lr[10] = 1.224744871391589*fl[23]+0.7071067811865475*fl[15]; 
  f_lr[11] = 1.224744871391589*fl[26]+0.7071067811865475*fl[19]; 
  f_lr[12] = 1.224744871391589*fl[27]+0.7071067811865475*fl[22]; 
  f_lr[13] = 1.224744871391589*fl[28]+0.7071067811865475*fl[24]; 
  f_lr[14] = 1.224744871391589*fl[29]+0.7071067811865475*fl[25]; 
  f_lr[15] = 1.224744871391589*fl[31]+0.7071067811865475*fl[30]; 
  f_lr[16] = 1.224744871391589*fl[33]+0.7071067811865475*fl[32]; 
  f_lr[17] = 1.224744871391589*fl[37]+0.7071067811865475*fl[34]; 
  f_lr[18] = 1.224744871391589*fl[38]+0.7071067811865475*fl[35]; 
  f_lr[19] = 1.224744871391589*fl[40]+0.7071067811865475*fl[36]; 
  f_lr[20] = 1.224744871391589*fl[43]+0.7071067811865475*fl[39]; 
  f_lr[21] = 1.224744871391589*fl[44]+0.7071067811865475*fl[41]; 
  f_lr[22] = 1.224744871391589*fl[45]+0.7071067811865475*fl[42]; 
  f_lr[23] = 1.224744871391589*fl[47]+0.7071067811865475*fl[46]; 

  f_cl[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[1]; 
  f_cl[1] = 0.7071067811865475*fc[2]-1.224744871391589*fc[6]; 
  f_cl[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[7]; 
  f_cl[3] = 0.7071067811865475*fc[4]-1.224744871391589*fc[9]; 
  f_cl[4] = 0.7071067811865475*fc[5]-1.224744871391589*fc[12]; 
  f_cl[5] = 0.7071067811865475*fc[8]-1.224744871391589*fc[16]; 
  f_cl[6] = 0.7071067811865475*fc[10]-1.224744871391589*fc[17]; 
  f_cl[7] = 0.7071067811865475*fc[11]-1.224744871391589*fc[18]; 
  f_cl[8] = 0.7071067811865475*fc[13]-1.224744871391589*fc[20]; 
  f_cl[9] = 0.7071067811865475*fc[14]-1.224744871391589*fc[21]; 
  f_cl[10] = 0.7071067811865475*fc[15]-1.224744871391589*fc[23]; 
  f_cl[11] = 0.7071067811865475*fc[19]-1.224744871391589*fc[26]; 
  f_cl[12] = 0.7071067811865475*fc[22]-1.224744871391589*fc[27]; 
  f_cl[13] = 0.7071067811865475*fc[24]-1.224744871391589*fc[28]; 
  f_cl[14] = 0.7071067811865475*fc[25]-1.224744871391589*fc[29]; 
  f_cl[15] = 0.7071067811865475*fc[30]-1.224744871391589*fc[31]; 
  f_cl[16] = 0.7071067811865475*fc[32]-1.224744871391589*fc[33]; 
  f_cl[17] = 0.7071067811865475*fc[34]-1.224744871391589*fc[37]; 
  f_cl[18] = 0.7071067811865475*fc[35]-1.224744871391589*fc[38]; 
  f_cl[19] = 0.7071067811865475*fc[36]-1.224744871391589*fc[40]; 
  f_cl[20] = 0.7071067811865475*fc[39]-1.224744871391589*fc[43]; 
  f_cl[21] = 0.7071067811865475*fc[41]-1.224744871391589*fc[44]; 
  f_cl[22] = 0.7071067811865475*fc[42]-1.224744871391589*fc[45]; 
  f_cl[23] = 0.7071067811865475*fc[46]-1.224744871391589*fc[47]; 

  fUpL[0] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[2]-0.3535533905932737*f_cl[2])+(0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[5]-0.3535533905932737*f_cl[5])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[2]+(0.3535533905932737*f_lr[0]-0.3535533905932737*f_cl[0])*sgn_alphaUpL[1]; 
  fUpL[3] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[7]-0.3535533905932737*f_cl[7])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[3]; 
  fUpL[4] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[9]-0.3535533905932737*f_cl[9])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[4]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[4]; 
  fUpL[5] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[5]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[5]+(0.3535533905932737*f_lr[1]-0.3535533905932737*f_cl[1])*sgn_alphaUpL[1]; 
  fUpL[6] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[11]-0.3535533905932737*f_cl[11])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[6]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[6]; 
  fUpL[7] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[7]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[7]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[3]-0.3535533905932737*f_cl[3]); 
  fUpL[8] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[12]-0.3535533905932737*f_cl[12])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[8]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[8]; 
  fUpL[9] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[9]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[9]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[4]-0.3535533905932737*f_cl[4]); 
  fUpL[10] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[14]-0.3535533905932737*f_cl[14])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[10]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[10]; 
  fUpL[11] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[11]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[11]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[6]-0.3535533905932737*f_cl[6]); 
  fUpL[12] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[12]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[12]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[8]-0.3535533905932737*f_cl[8]); 
  fUpL[13] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[15]-0.3535533905932737*f_cl[15])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[13]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[13]; 
  fUpL[14] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[14]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[14]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[10]-0.3535533905932737*f_cl[10]); 
  fUpL[15] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[15]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[15]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[13]-0.3535533905932737*f_cl[13]); 
  fUpL[16] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[18]-0.3535533905932737*f_cl[18])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[16]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[16]; 
  fUpL[17] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[20]-0.3535533905932737*f_cl[20])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[17]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[17]; 
  fUpL[18] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[18]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[18]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[16]-0.3535533905932737*f_cl[16]); 
  fUpL[19] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[22]-0.3535533905932737*f_cl[22])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[19]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[19]; 
  fUpL[20] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[20]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[20]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[17]-0.3535533905932737*f_cl[17]); 
  fUpL[21] = sgn_alphaUpL[1]*(0.3535533905932737*f_lr[23]-0.3535533905932737*f_cl[23])+(0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[21]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[21]; 
  fUpL[22] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[22]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[22]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[19]-0.3535533905932737*f_cl[19]); 
  fUpL[23] = (0.3535533905932737*sgn_alphaUpL[0]+0.5)*f_lr[23]+(0.5-0.3535533905932737*sgn_alphaUpL[0])*f_cl[23]+sgn_alphaUpL[1]*(0.3535533905932737*f_lr[21]-0.3535533905932737*f_cl[21]); 

  } 
  double fUpR[24] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = 1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  fUpR[1] = 1.224744871391589*fc[6]+0.7071067811865475*fc[2]; 
  fUpR[2] = 1.224744871391589*fc[7]+0.7071067811865475*fc[3]; 
  fUpR[3] = 1.224744871391589*fc[9]+0.7071067811865475*fc[4]; 
  fUpR[4] = 1.224744871391589*fc[12]+0.7071067811865475*fc[5]; 
  fUpR[5] = 1.224744871391589*fc[16]+0.7071067811865475*fc[8]; 
  fUpR[6] = 1.224744871391589*fc[17]+0.7071067811865475*fc[10]; 
  fUpR[7] = 1.224744871391589*fc[18]+0.7071067811865475*fc[11]; 
  fUpR[8] = 1.224744871391589*fc[20]+0.7071067811865475*fc[13]; 
  fUpR[9] = 1.224744871391589*fc[21]+0.7071067811865475*fc[14]; 
  fUpR[10] = 1.224744871391589*fc[23]+0.7071067811865475*fc[15]; 
  fUpR[11] = 1.224744871391589*fc[26]+0.7071067811865475*fc[19]; 
  fUpR[12] = 1.224744871391589*fc[27]+0.7071067811865475*fc[22]; 
  fUpR[13] = 1.224744871391589*fc[28]+0.7071067811865475*fc[24]; 
  fUpR[14] = 1.224744871391589*fc[29]+0.7071067811865475*fc[25]; 
  fUpR[15] = 1.224744871391589*fc[31]+0.7071067811865475*fc[30]; 
  fUpR[16] = 1.224744871391589*fc[33]+0.7071067811865475*fc[32]; 
  fUpR[17] = 1.224744871391589*fc[37]+0.7071067811865475*fc[34]; 
  fUpR[18] = 1.224744871391589*fc[38]+0.7071067811865475*fc[35]; 
  fUpR[19] = 1.224744871391589*fc[40]+0.7071067811865475*fc[36]; 
  fUpR[20] = 1.224744871391589*fc[43]+0.7071067811865475*fc[39]; 
  fUpR[21] = 1.224744871391589*fc[44]+0.7071067811865475*fc[41]; 
  fUpR[22] = 1.224744871391589*fc[45]+0.7071067811865475*fc[42]; 
  fUpR[23] = 1.224744871391589*fc[47]+0.7071067811865475*fc[46]; 
    } else { 
  fUpR[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[1]; 
  fUpR[1] = 0.7071067811865475*fr[2]-1.224744871391589*fr[6]; 
  fUpR[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[7]; 
  fUpR[3] = 0.7071067811865475*fr[4]-1.224744871391589*fr[9]; 
  fUpR[4] = 0.7071067811865475*fr[5]-1.224744871391589*fr[12]; 
  fUpR[5] = 0.7071067811865475*fr[8]-1.224744871391589*fr[16]; 
  fUpR[6] = 0.7071067811865475*fr[10]-1.224744871391589*fr[17]; 
  fUpR[7] = 0.7071067811865475*fr[11]-1.224744871391589*fr[18]; 
  fUpR[8] = 0.7071067811865475*fr[13]-1.224744871391589*fr[20]; 
  fUpR[9] = 0.7071067811865475*fr[14]-1.224744871391589*fr[21]; 
  fUpR[10] = 0.7071067811865475*fr[15]-1.224744871391589*fr[23]; 
  fUpR[11] = 0.7071067811865475*fr[19]-1.224744871391589*fr[26]; 
  fUpR[12] = 0.7071067811865475*fr[22]-1.224744871391589*fr[27]; 
  fUpR[13] = 0.7071067811865475*fr[24]-1.224744871391589*fr[28]; 
  fUpR[14] = 0.7071067811865475*fr[25]-1.224744871391589*fr[29]; 
  fUpR[15] = 0.7071067811865475*fr[30]-1.224744871391589*fr[31]; 
  fUpR[16] = 0.7071067811865475*fr[32]-1.224744871391589*fr[33]; 
  fUpR[17] = 0.7071067811865475*fr[34]-1.224744871391589*fr[37]; 
  fUpR[18] = 0.7071067811865475*fr[35]-1.224744871391589*fr[38]; 
  fUpR[19] = 0.7071067811865475*fr[36]-1.224744871391589*fr[40]; 
  fUpR[20] = 0.7071067811865475*fr[39]-1.224744871391589*fr[43]; 
  fUpR[21] = 0.7071067811865475*fr[41]-1.224744871391589*fr[44]; 
  fUpR[22] = 0.7071067811865475*fr[42]-1.224744871391589*fr[45]; 
  fUpR[23] = 0.7071067811865475*fr[46]-1.224744871391589*fr[47]; 
    } 
  } else { 
  double f_cr[24] = {0.};
  double f_rl[24] = {0.};
  double sgn_alphaUpR[2] = {0.};
  sgn_alphaUpR[0] = 0.7071067811865475*sgn_alpha_surfR[1]+0.7071067811865475*sgn_alpha_surfR[0]; 
  sgn_alphaUpR[1] = 0.7071067811865475*sgn_alpha_surfR[1]-0.7071067811865475*sgn_alpha_surfR[0]; 

  f_cr[0] = 1.224744871391589*fc[1]+0.7071067811865475*fc[0]; 
  f_cr[1] = 1.224744871391589*fc[6]+0.7071067811865475*fc[2]; 
  f_cr[2] = 1.224744871391589*fc[7]+0.7071067811865475*fc[3]; 
  f_cr[3] = 1.224744871391589*fc[9]+0.7071067811865475*fc[4]; 
  f_cr[4] = 1.224744871391589*fc[12]+0.7071067811865475*fc[5]; 
  f_cr[5] = 1.224744871391589*fc[16]+0.7071067811865475*fc[8]; 
  f_cr[6] = 1.224744871391589*fc[17]+0.7071067811865475*fc[10]; 
  f_cr[7] = 1.224744871391589*fc[18]+0.7071067811865475*fc[11]; 
  f_cr[8] = 1.224744871391589*fc[20]+0.7071067811865475*fc[13]; 
  f_cr[9] = 1.224744871391589*fc[21]+0.7071067811865475*fc[14]; 
  f_cr[10] = 1.224744871391589*fc[23]+0.7071067811865475*fc[15]; 
  f_cr[11] = 1.224744871391589*fc[26]+0.7071067811865475*fc[19]; 
  f_cr[12] = 1.224744871391589*fc[27]+0.7071067811865475*fc[22]; 
  f_cr[13] = 1.224744871391589*fc[28]+0.7071067811865475*fc[24]; 
  f_cr[14] = 1.224744871391589*fc[29]+0.7071067811865475*fc[25]; 
  f_cr[15] = 1.224744871391589*fc[31]+0.7071067811865475*fc[30]; 
  f_cr[16] = 1.224744871391589*fc[33]+0.7071067811865475*fc[32]; 
  f_cr[17] = 1.224744871391589*fc[37]+0.7071067811865475*fc[34]; 
  f_cr[18] = 1.224744871391589*fc[38]+0.7071067811865475*fc[35]; 
  f_cr[19] = 1.224744871391589*fc[40]+0.7071067811865475*fc[36]; 
  f_cr[20] = 1.224744871391589*fc[43]+0.7071067811865475*fc[39]; 
  f_cr[21] = 1.224744871391589*fc[44]+0.7071067811865475*fc[41]; 
  f_cr[22] = 1.224744871391589*fc[45]+0.7071067811865475*fc[42]; 
  f_cr[23] = 1.224744871391589*fc[47]+0.7071067811865475*fc[46]; 

  f_rl[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[1]; 
  f_rl[1] = 0.7071067811865475*fr[2]-1.224744871391589*fr[6]; 
  f_rl[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[7]; 
  f_rl[3] = 0.7071067811865475*fr[4]-1.224744871391589*fr[9]; 
  f_rl[4] = 0.7071067811865475*fr[5]-1.224744871391589*fr[12]; 
  f_rl[5] = 0.7071067811865475*fr[8]-1.224744871391589*fr[16]; 
  f_rl[6] = 0.7071067811865475*fr[10]-1.224744871391589*fr[17]; 
  f_rl[7] = 0.7071067811865475*fr[11]-1.224744871391589*fr[18]; 
  f_rl[8] = 0.7071067811865475*fr[13]-1.224744871391589*fr[20]; 
  f_rl[9] = 0.7071067811865475*fr[14]-1.224744871391589*fr[21]; 
  f_rl[10] = 0.7071067811865475*fr[15]-1.224744871391589*fr[23]; 
  f_rl[11] = 0.7071067811865475*fr[19]-1.224744871391589*fr[26]; 
  f_rl[12] = 0.7071067811865475*fr[22]-1.224744871391589*fr[27]; 
  f_rl[13] = 0.7071067811865475*fr[24]-1.224744871391589*fr[28]; 
  f_rl[14] = 0.7071067811865475*fr[25]-1.224744871391589*fr[29]; 
  f_rl[15] = 0.7071067811865475*fr[30]-1.224744871391589*fr[31]; 
  f_rl[16] = 0.7071067811865475*fr[32]-1.224744871391589*fr[33]; 
  f_rl[17] = 0.7071067811865475*fr[34]-1.224744871391589*fr[37]; 
  f_rl[18] = 0.7071067811865475*fr[35]-1.224744871391589*fr[38]; 
  f_rl[19] = 0.7071067811865475*fr[36]-1.224744871391589*fr[40]; 
  f_rl[20] = 0.7071067811865475*fr[39]-1.224744871391589*fr[43]; 
  f_rl[21] = 0.7071067811865475*fr[41]-1.224744871391589*fr[44]; 
  f_rl[22] = 0.7071067811865475*fr[42]-1.224744871391589*fr[45]; 
  f_rl[23] = 0.7071067811865475*fr[46]-1.224744871391589*fr[47]; 

  fUpR[0] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[2]-0.3535533905932737*f_rl[2])+(0.3535533905932737*f_cr[0]-0.3535533905932737*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[5]-0.3535533905932737*f_rl[5])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[1]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[2]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[2]+(0.3535533905932737*f_cr[0]-0.3535533905932737*f_rl[0])*sgn_alphaUpR[1]; 
  fUpR[3] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[7]-0.3535533905932737*f_rl[7])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[3]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[3]; 
  fUpR[4] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[9]-0.3535533905932737*f_rl[9])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[4]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[4]; 
  fUpR[5] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[5]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[5]+(0.3535533905932737*f_cr[1]-0.3535533905932737*f_rl[1])*sgn_alphaUpR[1]; 
  fUpR[6] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[11]-0.3535533905932737*f_rl[11])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[6]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[6]; 
  fUpR[7] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[7]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[7]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[3]-0.3535533905932737*f_rl[3]); 
  fUpR[8] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[12]-0.3535533905932737*f_rl[12])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[8]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[8]; 
  fUpR[9] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[9]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[9]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[4]-0.3535533905932737*f_rl[4]); 
  fUpR[10] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[14]-0.3535533905932737*f_rl[14])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[10]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[10]; 
  fUpR[11] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[11]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[11]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[6]-0.3535533905932737*f_rl[6]); 
  fUpR[12] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[12]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[12]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[8]-0.3535533905932737*f_rl[8]); 
  fUpR[13] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[15]-0.3535533905932737*f_rl[15])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[13]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[13]; 
  fUpR[14] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[14]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[14]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[10]-0.3535533905932737*f_rl[10]); 
  fUpR[15] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[15]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[15]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[13]-0.3535533905932737*f_rl[13]); 
  fUpR[16] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[18]-0.3535533905932737*f_rl[18])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[16]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[16]; 
  fUpR[17] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[20]-0.3535533905932737*f_rl[20])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[17]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[17]; 
  fUpR[18] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[18]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[18]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[16]-0.3535533905932737*f_rl[16]); 
  fUpR[19] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[22]-0.3535533905932737*f_rl[22])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[19]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[19]; 
  fUpR[20] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[20]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[20]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[17]-0.3535533905932737*f_rl[17]); 
  fUpR[21] = sgn_alphaUpR[1]*(0.3535533905932737*f_cr[23]-0.3535533905932737*f_rl[23])+(0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[21]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[21]; 
  fUpR[22] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[22]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[22]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[19]-0.3535533905932737*f_rl[19]); 
  fUpR[23] = (0.5-0.3535533905932737*sgn_alphaUpR[0])*f_rl[23]+(0.3535533905932737*sgn_alphaUpR[0]+0.5)*f_cr[23]+sgn_alphaUpR[1]*(0.3535533905932737*f_cr[21]-0.3535533905932737*f_rl[21]); 

  } 
  double GhatL[24] = {0.};
  double GhatR[24] = {0.};
  GhatL[0] = 0.25*(alphaL[2]*fUpL[2]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.25*(alphaL[2]*fUpL[5]+alphaL[0]*fUpL[1]); 
  GhatL[2] = 0.25*(alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.25*(alphaL[2]*fUpL[7]+alphaL[0]*fUpL[3]); 
  GhatL[4] = 0.25*(alphaL[2]*fUpL[9]+alphaL[0]*fUpL[4]); 
  GhatL[5] = 0.25*(alphaL[0]*fUpL[5]+fUpL[1]*alphaL[2]); 
  GhatL[6] = 0.25*(alphaL[2]*fUpL[11]+alphaL[0]*fUpL[6]); 
  GhatL[7] = 0.25*(alphaL[0]*fUpL[7]+alphaL[2]*fUpL[3]); 
  GhatL[8] = 0.25*(alphaL[2]*fUpL[12]+alphaL[0]*fUpL[8]); 
  GhatL[9] = 0.25*(alphaL[0]*fUpL[9]+alphaL[2]*fUpL[4]); 
  GhatL[10] = 0.25*(alphaL[2]*fUpL[14]+alphaL[0]*fUpL[10]); 
  GhatL[11] = 0.25*(alphaL[0]*fUpL[11]+alphaL[2]*fUpL[6]); 
  GhatL[12] = 0.25*(alphaL[0]*fUpL[12]+alphaL[2]*fUpL[8]); 
  GhatL[13] = 0.25*(alphaL[2]*fUpL[15]+alphaL[0]*fUpL[13]); 
  GhatL[14] = 0.25*(alphaL[0]*fUpL[14]+alphaL[2]*fUpL[10]); 
  GhatL[15] = 0.25*(alphaL[0]*fUpL[15]+alphaL[2]*fUpL[13]); 
  GhatL[16] = 0.2500000000000001*alphaL[2]*fUpL[18]+0.25*alphaL[0]*fUpL[16]; 
  GhatL[17] = 0.2500000000000001*alphaL[2]*fUpL[20]+0.25*alphaL[0]*fUpL[17]; 
  GhatL[18] = 0.25*alphaL[0]*fUpL[18]+0.2500000000000001*alphaL[2]*fUpL[16]; 
  GhatL[19] = 0.2500000000000001*alphaL[2]*fUpL[22]+0.25*alphaL[0]*fUpL[19]; 
  GhatL[20] = 0.25*alphaL[0]*fUpL[20]+0.2500000000000001*alphaL[2]*fUpL[17]; 
  GhatL[21] = 0.2500000000000001*alphaL[2]*fUpL[23]+0.25*alphaL[0]*fUpL[21]; 
  GhatL[22] = 0.25*alphaL[0]*fUpL[22]+0.2500000000000001*alphaL[2]*fUpL[19]; 
  GhatL[23] = 0.25*alphaL[0]*fUpL[23]+0.2500000000000001*alphaL[2]*fUpL[21]; 

  GhatR[0] = 0.25*(alphaR[2]*fUpR[2]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.25*(alphaR[2]*fUpR[5]+alphaR[0]*fUpR[1]); 
  GhatR[2] = 0.25*(alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.25*(alphaR[2]*fUpR[7]+alphaR[0]*fUpR[3]); 
  GhatR[4] = 0.25*(alphaR[2]*fUpR[9]+alphaR[0]*fUpR[4]); 
  GhatR[5] = 0.25*(alphaR[0]*fUpR[5]+fUpR[1]*alphaR[2]); 
  GhatR[6] = 0.25*(alphaR[2]*fUpR[11]+alphaR[0]*fUpR[6]); 
  GhatR[7] = 0.25*(alphaR[0]*fUpR[7]+alphaR[2]*fUpR[3]); 
  GhatR[8] = 0.25*(alphaR[2]*fUpR[12]+alphaR[0]*fUpR[8]); 
  GhatR[9] = 0.25*(alphaR[0]*fUpR[9]+alphaR[2]*fUpR[4]); 
  GhatR[10] = 0.25*(alphaR[2]*fUpR[14]+alphaR[0]*fUpR[10]); 
  GhatR[11] = 0.25*(alphaR[0]*fUpR[11]+alphaR[2]*fUpR[6]); 
  GhatR[12] = 0.25*(alphaR[0]*fUpR[12]+alphaR[2]*fUpR[8]); 
  GhatR[13] = 0.25*(alphaR[2]*fUpR[15]+alphaR[0]*fUpR[13]); 
  GhatR[14] = 0.25*(alphaR[0]*fUpR[14]+alphaR[2]*fUpR[10]); 
  GhatR[15] = 0.25*(alphaR[0]*fUpR[15]+alphaR[2]*fUpR[13]); 
  GhatR[16] = 0.2500000000000001*alphaR[2]*fUpR[18]+0.25*alphaR[0]*fUpR[16]; 
  GhatR[17] = 0.2500000000000001*alphaR[2]*fUpR[20]+0.25*alphaR[0]*fUpR[17]; 
  GhatR[18] = 0.25*alphaR[0]*fUpR[18]+0.2500000000000001*alphaR[2]*fUpR[16]; 
  GhatR[19] = 0.2500000000000001*alphaR[2]*fUpR[22]+0.25*alphaR[0]*fUpR[19]; 
  GhatR[20] = 0.25*alphaR[0]*fUpR[20]+0.2500000000000001*alphaR[2]*fUpR[17]; 
  GhatR[21] = 0.2500000000000001*alphaR[2]*fUpR[23]+0.25*alphaR[0]*fUpR[21]; 
  GhatR[22] = 0.25*alphaR[0]*fUpR[22]+0.2500000000000001*alphaR[2]*fUpR[19]; 
  GhatR[23] = 0.25*alphaR[0]*fUpR[23]+0.2500000000000001*alphaR[2]*fUpR[21]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdx2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdx2; 
  out[5] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdx2; 
  out[6] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 
  out[7] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdx2; 
  out[8] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdx2; 
  out[9] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdx2; 
  out[10] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdx2; 
  out[11] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdx2; 
  out[12] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdx2; 
  out[13] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdx2; 
  out[14] += (0.7071067811865475*GhatL[9]-0.7071067811865475*GhatR[9])*rdx2; 
  out[15] += (0.7071067811865475*GhatL[10]-0.7071067811865475*GhatR[10])*rdx2; 
  out[16] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdx2; 
  out[17] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdx2; 
  out[18] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdx2; 
  out[19] += (0.7071067811865475*GhatL[11]-0.7071067811865475*GhatR[11])*rdx2; 
  out[20] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdx2; 
  out[21] += ((-1.224744871391589*GhatR[9])-1.224744871391589*GhatL[9])*rdx2; 
  out[22] += (0.7071067811865475*GhatL[12]-0.7071067811865475*GhatR[12])*rdx2; 
  out[23] += ((-1.224744871391589*GhatR[10])-1.224744871391589*GhatL[10])*rdx2; 
  out[24] += (0.7071067811865475*GhatL[13]-0.7071067811865475*GhatR[13])*rdx2; 
  out[25] += (0.7071067811865475*GhatL[14]-0.7071067811865475*GhatR[14])*rdx2; 
  out[26] += ((-1.224744871391589*GhatR[11])-1.224744871391589*GhatL[11])*rdx2; 
  out[27] += ((-1.224744871391589*GhatR[12])-1.224744871391589*GhatL[12])*rdx2; 
  out[28] += ((-1.224744871391589*GhatR[13])-1.224744871391589*GhatL[13])*rdx2; 
  out[29] += ((-1.224744871391589*GhatR[14])-1.224744871391589*GhatL[14])*rdx2; 
  out[30] += (0.7071067811865475*GhatL[15]-0.7071067811865475*GhatR[15])*rdx2; 
  out[31] += ((-1.224744871391589*GhatR[15])-1.224744871391589*GhatL[15])*rdx2; 
  out[32] += (0.7071067811865475*GhatL[16]-0.7071067811865475*GhatR[16])*rdx2; 
  out[33] += ((-1.224744871391589*GhatR[16])-1.224744871391589*GhatL[16])*rdx2; 
  out[34] += (0.7071067811865475*GhatL[17]-0.7071067811865475*GhatR[17])*rdx2; 
  out[35] += (0.7071067811865475*GhatL[18]-0.7071067811865475*GhatR[18])*rdx2; 
  out[36] += (0.7071067811865475*GhatL[19]-0.7071067811865475*GhatR[19])*rdx2; 
  out[37] += ((-1.224744871391589*GhatR[17])-1.224744871391589*GhatL[17])*rdx2; 
  out[38] += ((-1.224744871391589*GhatR[18])-1.224744871391589*GhatL[18])*rdx2; 
  out[39] += (0.7071067811865475*GhatL[20]-0.7071067811865475*GhatR[20])*rdx2; 
  out[40] += ((-1.224744871391589*GhatR[19])-1.224744871391589*GhatL[19])*rdx2; 
  out[41] += (0.7071067811865475*GhatL[21]-0.7071067811865475*GhatR[21])*rdx2; 
  out[42] += (0.7071067811865475*GhatL[22]-0.7071067811865475*GhatR[22])*rdx2; 
  out[43] += ((-1.224744871391589*GhatR[20])-1.224744871391589*GhatL[20])*rdx2; 
  out[44] += ((-1.224744871391589*GhatR[21])-1.224744871391589*GhatL[21])*rdx2; 
  out[45] += ((-1.224744871391589*GhatR[22])-1.224744871391589*GhatL[22])*rdx2; 
  out[46] += (0.7071067811865475*GhatL[23]-0.7071067811865475*GhatR[23])*rdx2; 
  out[47] += ((-1.224744871391589*GhatR[23])-1.224744871391589*GhatL[23])*rdx2; 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 0.375*rdx2*cflFreq; 

} 