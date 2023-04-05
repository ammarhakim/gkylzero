#include <gkyl_lbo_vlasov_pkpm_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_drag_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nu, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:    cell-center coordinates. 
  // dxv[3]:  cell spacing. 
  // nu:       collisionalities added (self and cross species collisionalities). 
  // fl/fc/fr: Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells 
  // out:      incremented distribution function in cell 

  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[12]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[12]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[12]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 

  double alphaDrSurf_l[4] = {0.0}; 
  alphaDrSurf_l[0] = nu[0]*wvpar-0.5*nu[0]*dvpar; 
  alphaDrSurf_l[1] = nu[1]*wvpar-0.5*nu[1]*dvpar; 
  alphaDrSurf_l[2] = nu[2]*wvpar-0.5*nu[2]*dvpar; 
  alphaDrSurf_l[3] = nu[3]*wvpar-0.5*nu[3]*dvpar; 

  double alphaDrSurf_r[4] = {0.0}; 
  alphaDrSurf_r[0] = nu[0]*wvpar+0.5*nu[0]*dvpar; 
  alphaDrSurf_r[1] = nu[1]*wvpar+0.5*nu[1]*dvpar; 
  alphaDrSurf_r[2] = nu[2]*wvpar+0.5*nu[2]*dvpar; 
  alphaDrSurf_r[3] = nu[3]*wvpar+0.5*nu[3]*dvpar; 

  double Ghat_F_0_r[4]; 
  double Ghat_F_0_l[4]; 
  double Ghat_G_1_r[4]; 
  double Ghat_G_1_l[4]; 
  if (wvpar>0) { 

  Ghat_F_0_r[0] = 0.7905694150420947*(alphaDrSurf_r[3]*F_0r[11]+alphaDrSurf_r[2]*F_0r[10]+alphaDrSurf_r[1]*F_0r[9]+alphaDrSurf_r[0]*F_0r[8])-0.6123724356957944*(alphaDrSurf_r[3]*F_0r[7]+alphaDrSurf_r[2]*F_0r[6]+alphaDrSurf_r[1]*F_0r[5])+0.3535533905932737*alphaDrSurf_r[3]*F_0r[4]-0.6123724356957944*alphaDrSurf_r[0]*F_0r[3]+0.3535533905932737*(F_0r[2]*alphaDrSurf_r[2]+F_0r[1]*alphaDrSurf_r[1]+F_0r[0]*alphaDrSurf_r[0]); 
  Ghat_F_0_r[1] = 0.7905694150420947*(alphaDrSurf_r[2]*F_0r[11]+alphaDrSurf_r[3]*F_0r[10]+alphaDrSurf_r[0]*F_0r[9]+alphaDrSurf_r[1]*F_0r[8])-0.6123724356957944*(alphaDrSurf_r[2]*F_0r[7]+alphaDrSurf_r[3]*F_0r[6]+alphaDrSurf_r[0]*F_0r[5])+0.3535533905932737*(alphaDrSurf_r[2]*F_0r[4]+F_0r[2]*alphaDrSurf_r[3])-0.6123724356957944*alphaDrSurf_r[1]*F_0r[3]+0.3535533905932737*(F_0r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*F_0r[1]); 
  Ghat_F_0_r[2] = 0.7905694150420947*(alphaDrSurf_r[1]*F_0r[11]+alphaDrSurf_r[0]*F_0r[10]+alphaDrSurf_r[3]*F_0r[9]+alphaDrSurf_r[2]*F_0r[8])-0.6123724356957944*(alphaDrSurf_r[1]*F_0r[7]+alphaDrSurf_r[0]*F_0r[6]+alphaDrSurf_r[3]*F_0r[5])+0.3535533905932737*(alphaDrSurf_r[1]*F_0r[4]+F_0r[1]*alphaDrSurf_r[3])-0.6123724356957944*alphaDrSurf_r[2]*F_0r[3]+0.3535533905932737*(F_0r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*F_0r[2]); 
  Ghat_F_0_r[3] = 0.7905694150420947*(alphaDrSurf_r[0]*F_0r[11]+alphaDrSurf_r[1]*F_0r[10]+alphaDrSurf_r[2]*F_0r[9]+alphaDrSurf_r[3]*F_0r[8])-0.6123724356957944*(alphaDrSurf_r[0]*F_0r[7]+alphaDrSurf_r[1]*F_0r[6]+alphaDrSurf_r[2]*F_0r[5])+0.3535533905932737*alphaDrSurf_r[0]*F_0r[4]-0.6123724356957944*F_0r[3]*alphaDrSurf_r[3]+0.3535533905932737*(F_0r[0]*alphaDrSurf_r[3]+F_0r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*F_0r[2]); 
  Ghat_G_1_r[0] = 0.7905694150420947*(alphaDrSurf_r[3]*G_1r[11]+alphaDrSurf_r[2]*G_1r[10]+alphaDrSurf_r[1]*G_1r[9]+alphaDrSurf_r[0]*G_1r[8])-0.6123724356957944*(alphaDrSurf_r[3]*G_1r[7]+alphaDrSurf_r[2]*G_1r[6]+alphaDrSurf_r[1]*G_1r[5])+0.3535533905932737*alphaDrSurf_r[3]*G_1r[4]-0.6123724356957944*alphaDrSurf_r[0]*G_1r[3]+0.3535533905932737*(G_1r[2]*alphaDrSurf_r[2]+G_1r[1]*alphaDrSurf_r[1]+G_1r[0]*alphaDrSurf_r[0]); 
  Ghat_G_1_r[1] = 0.7905694150420947*(alphaDrSurf_r[2]*G_1r[11]+alphaDrSurf_r[3]*G_1r[10]+alphaDrSurf_r[0]*G_1r[9]+alphaDrSurf_r[1]*G_1r[8])-0.6123724356957944*(alphaDrSurf_r[2]*G_1r[7]+alphaDrSurf_r[3]*G_1r[6]+alphaDrSurf_r[0]*G_1r[5])+0.3535533905932737*(alphaDrSurf_r[2]*G_1r[4]+G_1r[2]*alphaDrSurf_r[3])-0.6123724356957944*alphaDrSurf_r[1]*G_1r[3]+0.3535533905932737*(G_1r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*G_1r[1]); 
  Ghat_G_1_r[2] = 0.7905694150420947*(alphaDrSurf_r[1]*G_1r[11]+alphaDrSurf_r[0]*G_1r[10]+alphaDrSurf_r[3]*G_1r[9]+alphaDrSurf_r[2]*G_1r[8])-0.6123724356957944*(alphaDrSurf_r[1]*G_1r[7]+alphaDrSurf_r[0]*G_1r[6]+alphaDrSurf_r[3]*G_1r[5])+0.3535533905932737*(alphaDrSurf_r[1]*G_1r[4]+G_1r[1]*alphaDrSurf_r[3])-0.6123724356957944*alphaDrSurf_r[2]*G_1r[3]+0.3535533905932737*(G_1r[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*G_1r[2]); 
  Ghat_G_1_r[3] = 0.7905694150420947*(alphaDrSurf_r[0]*G_1r[11]+alphaDrSurf_r[1]*G_1r[10]+alphaDrSurf_r[2]*G_1r[9]+alphaDrSurf_r[3]*G_1r[8])-0.6123724356957944*(alphaDrSurf_r[0]*G_1r[7]+alphaDrSurf_r[1]*G_1r[6]+alphaDrSurf_r[2]*G_1r[5])+0.3535533905932737*alphaDrSurf_r[0]*G_1r[4]-0.6123724356957944*G_1r[3]*alphaDrSurf_r[3]+0.3535533905932737*(G_1r[0]*alphaDrSurf_r[3]+G_1r[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*G_1r[2]); 

  Ghat_F_0_l[0] = 0.7905694150420947*(alphaDrSurf_l[3]*F_0c[11]+alphaDrSurf_l[2]*F_0c[10]+alphaDrSurf_l[1]*F_0c[9]+alphaDrSurf_l[0]*F_0c[8])-0.6123724356957944*(alphaDrSurf_l[3]*F_0c[7]+alphaDrSurf_l[2]*F_0c[6]+alphaDrSurf_l[1]*F_0c[5])+0.3535533905932737*alphaDrSurf_l[3]*F_0c[4]-0.6123724356957944*alphaDrSurf_l[0]*F_0c[3]+0.3535533905932737*(F_0c[2]*alphaDrSurf_l[2]+F_0c[1]*alphaDrSurf_l[1]+F_0c[0]*alphaDrSurf_l[0]); 
  Ghat_F_0_l[1] = 0.7905694150420947*(alphaDrSurf_l[2]*F_0c[11]+alphaDrSurf_l[3]*F_0c[10]+alphaDrSurf_l[0]*F_0c[9]+alphaDrSurf_l[1]*F_0c[8])-0.6123724356957944*(alphaDrSurf_l[2]*F_0c[7]+alphaDrSurf_l[3]*F_0c[6]+alphaDrSurf_l[0]*F_0c[5])+0.3535533905932737*(alphaDrSurf_l[2]*F_0c[4]+F_0c[2]*alphaDrSurf_l[3])-0.6123724356957944*alphaDrSurf_l[1]*F_0c[3]+0.3535533905932737*(F_0c[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*F_0c[1]); 
  Ghat_F_0_l[2] = 0.7905694150420947*(alphaDrSurf_l[1]*F_0c[11]+alphaDrSurf_l[0]*F_0c[10]+alphaDrSurf_l[3]*F_0c[9]+alphaDrSurf_l[2]*F_0c[8])-0.6123724356957944*(alphaDrSurf_l[1]*F_0c[7]+alphaDrSurf_l[0]*F_0c[6]+alphaDrSurf_l[3]*F_0c[5])+0.3535533905932737*(alphaDrSurf_l[1]*F_0c[4]+F_0c[1]*alphaDrSurf_l[3])-0.6123724356957944*alphaDrSurf_l[2]*F_0c[3]+0.3535533905932737*(F_0c[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*F_0c[2]); 
  Ghat_F_0_l[3] = 0.7905694150420947*(alphaDrSurf_l[0]*F_0c[11]+alphaDrSurf_l[1]*F_0c[10]+alphaDrSurf_l[2]*F_0c[9]+alphaDrSurf_l[3]*F_0c[8])-0.6123724356957944*(alphaDrSurf_l[0]*F_0c[7]+alphaDrSurf_l[1]*F_0c[6]+alphaDrSurf_l[2]*F_0c[5])+0.3535533905932737*alphaDrSurf_l[0]*F_0c[4]-0.6123724356957944*F_0c[3]*alphaDrSurf_l[3]+0.3535533905932737*(F_0c[0]*alphaDrSurf_l[3]+F_0c[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*F_0c[2]); 
  Ghat_G_1_l[0] = 0.7905694150420947*(alphaDrSurf_l[3]*G_1c[11]+alphaDrSurf_l[2]*G_1c[10]+alphaDrSurf_l[1]*G_1c[9]+alphaDrSurf_l[0]*G_1c[8])-0.6123724356957944*(alphaDrSurf_l[3]*G_1c[7]+alphaDrSurf_l[2]*G_1c[6]+alphaDrSurf_l[1]*G_1c[5])+0.3535533905932737*alphaDrSurf_l[3]*G_1c[4]-0.6123724356957944*alphaDrSurf_l[0]*G_1c[3]+0.3535533905932737*(G_1c[2]*alphaDrSurf_l[2]+G_1c[1]*alphaDrSurf_l[1]+G_1c[0]*alphaDrSurf_l[0]); 
  Ghat_G_1_l[1] = 0.7905694150420947*(alphaDrSurf_l[2]*G_1c[11]+alphaDrSurf_l[3]*G_1c[10]+alphaDrSurf_l[0]*G_1c[9]+alphaDrSurf_l[1]*G_1c[8])-0.6123724356957944*(alphaDrSurf_l[2]*G_1c[7]+alphaDrSurf_l[3]*G_1c[6]+alphaDrSurf_l[0]*G_1c[5])+0.3535533905932737*(alphaDrSurf_l[2]*G_1c[4]+G_1c[2]*alphaDrSurf_l[3])-0.6123724356957944*alphaDrSurf_l[1]*G_1c[3]+0.3535533905932737*(G_1c[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*G_1c[1]); 
  Ghat_G_1_l[2] = 0.7905694150420947*(alphaDrSurf_l[1]*G_1c[11]+alphaDrSurf_l[0]*G_1c[10]+alphaDrSurf_l[3]*G_1c[9]+alphaDrSurf_l[2]*G_1c[8])-0.6123724356957944*(alphaDrSurf_l[1]*G_1c[7]+alphaDrSurf_l[0]*G_1c[6]+alphaDrSurf_l[3]*G_1c[5])+0.3535533905932737*(alphaDrSurf_l[1]*G_1c[4]+G_1c[1]*alphaDrSurf_l[3])-0.6123724356957944*alphaDrSurf_l[2]*G_1c[3]+0.3535533905932737*(G_1c[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*G_1c[2]); 
  Ghat_G_1_l[3] = 0.7905694150420947*(alphaDrSurf_l[0]*G_1c[11]+alphaDrSurf_l[1]*G_1c[10]+alphaDrSurf_l[2]*G_1c[9]+alphaDrSurf_l[3]*G_1c[8])-0.6123724356957944*(alphaDrSurf_l[0]*G_1c[7]+alphaDrSurf_l[1]*G_1c[6]+alphaDrSurf_l[2]*G_1c[5])+0.3535533905932737*alphaDrSurf_l[0]*G_1c[4]-0.6123724356957944*G_1c[3]*alphaDrSurf_l[3]+0.3535533905932737*(G_1c[0]*alphaDrSurf_l[3]+G_1c[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*G_1c[2]); 

  } else { 

  Ghat_F_0_r[0] = 0.7905694150420947*(alphaDrSurf_r[3]*F_0c[11]+alphaDrSurf_r[2]*F_0c[10]+alphaDrSurf_r[1]*F_0c[9]+alphaDrSurf_r[0]*F_0c[8])+0.6123724356957944*(alphaDrSurf_r[3]*F_0c[7]+alphaDrSurf_r[2]*F_0c[6]+alphaDrSurf_r[1]*F_0c[5])+0.3535533905932737*alphaDrSurf_r[3]*F_0c[4]+0.6123724356957944*alphaDrSurf_r[0]*F_0c[3]+0.3535533905932737*(F_0c[2]*alphaDrSurf_r[2]+F_0c[1]*alphaDrSurf_r[1]+F_0c[0]*alphaDrSurf_r[0]); 
  Ghat_F_0_r[1] = 0.7905694150420947*(alphaDrSurf_r[2]*F_0c[11]+alphaDrSurf_r[3]*F_0c[10]+alphaDrSurf_r[0]*F_0c[9]+alphaDrSurf_r[1]*F_0c[8])+0.6123724356957944*(alphaDrSurf_r[2]*F_0c[7]+alphaDrSurf_r[3]*F_0c[6]+alphaDrSurf_r[0]*F_0c[5])+0.3535533905932737*(alphaDrSurf_r[2]*F_0c[4]+F_0c[2]*alphaDrSurf_r[3])+0.6123724356957944*alphaDrSurf_r[1]*F_0c[3]+0.3535533905932737*(F_0c[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*F_0c[1]); 
  Ghat_F_0_r[2] = 0.7905694150420947*(alphaDrSurf_r[1]*F_0c[11]+alphaDrSurf_r[0]*F_0c[10]+alphaDrSurf_r[3]*F_0c[9]+alphaDrSurf_r[2]*F_0c[8])+0.6123724356957944*(alphaDrSurf_r[1]*F_0c[7]+alphaDrSurf_r[0]*F_0c[6]+alphaDrSurf_r[3]*F_0c[5])+0.3535533905932737*(alphaDrSurf_r[1]*F_0c[4]+F_0c[1]*alphaDrSurf_r[3])+0.6123724356957944*alphaDrSurf_r[2]*F_0c[3]+0.3535533905932737*(F_0c[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*F_0c[2]); 
  Ghat_F_0_r[3] = 0.7905694150420947*(alphaDrSurf_r[0]*F_0c[11]+alphaDrSurf_r[1]*F_0c[10]+alphaDrSurf_r[2]*F_0c[9]+alphaDrSurf_r[3]*F_0c[8])+0.6123724356957944*(alphaDrSurf_r[0]*F_0c[7]+alphaDrSurf_r[1]*F_0c[6]+alphaDrSurf_r[2]*F_0c[5])+0.3535533905932737*alphaDrSurf_r[0]*F_0c[4]+0.6123724356957944*F_0c[3]*alphaDrSurf_r[3]+0.3535533905932737*(F_0c[0]*alphaDrSurf_r[3]+F_0c[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*F_0c[2]); 
  Ghat_G_1_r[0] = 0.7905694150420947*(alphaDrSurf_r[3]*G_1c[11]+alphaDrSurf_r[2]*G_1c[10]+alphaDrSurf_r[1]*G_1c[9]+alphaDrSurf_r[0]*G_1c[8])+0.6123724356957944*(alphaDrSurf_r[3]*G_1c[7]+alphaDrSurf_r[2]*G_1c[6]+alphaDrSurf_r[1]*G_1c[5])+0.3535533905932737*alphaDrSurf_r[3]*G_1c[4]+0.6123724356957944*alphaDrSurf_r[0]*G_1c[3]+0.3535533905932737*(G_1c[2]*alphaDrSurf_r[2]+G_1c[1]*alphaDrSurf_r[1]+G_1c[0]*alphaDrSurf_r[0]); 
  Ghat_G_1_r[1] = 0.7905694150420947*(alphaDrSurf_r[2]*G_1c[11]+alphaDrSurf_r[3]*G_1c[10]+alphaDrSurf_r[0]*G_1c[9]+alphaDrSurf_r[1]*G_1c[8])+0.6123724356957944*(alphaDrSurf_r[2]*G_1c[7]+alphaDrSurf_r[3]*G_1c[6]+alphaDrSurf_r[0]*G_1c[5])+0.3535533905932737*(alphaDrSurf_r[2]*G_1c[4]+G_1c[2]*alphaDrSurf_r[3])+0.6123724356957944*alphaDrSurf_r[1]*G_1c[3]+0.3535533905932737*(G_1c[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*G_1c[1]); 
  Ghat_G_1_r[2] = 0.7905694150420947*(alphaDrSurf_r[1]*G_1c[11]+alphaDrSurf_r[0]*G_1c[10]+alphaDrSurf_r[3]*G_1c[9]+alphaDrSurf_r[2]*G_1c[8])+0.6123724356957944*(alphaDrSurf_r[1]*G_1c[7]+alphaDrSurf_r[0]*G_1c[6]+alphaDrSurf_r[3]*G_1c[5])+0.3535533905932737*(alphaDrSurf_r[1]*G_1c[4]+G_1c[1]*alphaDrSurf_r[3])+0.6123724356957944*alphaDrSurf_r[2]*G_1c[3]+0.3535533905932737*(G_1c[0]*alphaDrSurf_r[2]+alphaDrSurf_r[0]*G_1c[2]); 
  Ghat_G_1_r[3] = 0.7905694150420947*(alphaDrSurf_r[0]*G_1c[11]+alphaDrSurf_r[1]*G_1c[10]+alphaDrSurf_r[2]*G_1c[9]+alphaDrSurf_r[3]*G_1c[8])+0.6123724356957944*(alphaDrSurf_r[0]*G_1c[7]+alphaDrSurf_r[1]*G_1c[6]+alphaDrSurf_r[2]*G_1c[5])+0.3535533905932737*alphaDrSurf_r[0]*G_1c[4]+0.6123724356957944*G_1c[3]*alphaDrSurf_r[3]+0.3535533905932737*(G_1c[0]*alphaDrSurf_r[3]+G_1c[1]*alphaDrSurf_r[2]+alphaDrSurf_r[1]*G_1c[2]); 

  Ghat_F_0_l[0] = 0.7905694150420947*(alphaDrSurf_l[3]*F_0l[11]+alphaDrSurf_l[2]*F_0l[10]+alphaDrSurf_l[1]*F_0l[9]+alphaDrSurf_l[0]*F_0l[8])+0.6123724356957944*(alphaDrSurf_l[3]*F_0l[7]+alphaDrSurf_l[2]*F_0l[6]+alphaDrSurf_l[1]*F_0l[5])+0.3535533905932737*alphaDrSurf_l[3]*F_0l[4]+0.6123724356957944*alphaDrSurf_l[0]*F_0l[3]+0.3535533905932737*(F_0l[2]*alphaDrSurf_l[2]+F_0l[1]*alphaDrSurf_l[1]+F_0l[0]*alphaDrSurf_l[0]); 
  Ghat_F_0_l[1] = 0.7905694150420947*(alphaDrSurf_l[2]*F_0l[11]+alphaDrSurf_l[3]*F_0l[10]+alphaDrSurf_l[0]*F_0l[9]+alphaDrSurf_l[1]*F_0l[8])+0.6123724356957944*(alphaDrSurf_l[2]*F_0l[7]+alphaDrSurf_l[3]*F_0l[6]+alphaDrSurf_l[0]*F_0l[5])+0.3535533905932737*(alphaDrSurf_l[2]*F_0l[4]+F_0l[2]*alphaDrSurf_l[3])+0.6123724356957944*alphaDrSurf_l[1]*F_0l[3]+0.3535533905932737*(F_0l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*F_0l[1]); 
  Ghat_F_0_l[2] = 0.7905694150420947*(alphaDrSurf_l[1]*F_0l[11]+alphaDrSurf_l[0]*F_0l[10]+alphaDrSurf_l[3]*F_0l[9]+alphaDrSurf_l[2]*F_0l[8])+0.6123724356957944*(alphaDrSurf_l[1]*F_0l[7]+alphaDrSurf_l[0]*F_0l[6]+alphaDrSurf_l[3]*F_0l[5])+0.3535533905932737*(alphaDrSurf_l[1]*F_0l[4]+F_0l[1]*alphaDrSurf_l[3])+0.6123724356957944*alphaDrSurf_l[2]*F_0l[3]+0.3535533905932737*(F_0l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*F_0l[2]); 
  Ghat_F_0_l[3] = 0.7905694150420947*(alphaDrSurf_l[0]*F_0l[11]+alphaDrSurf_l[1]*F_0l[10]+alphaDrSurf_l[2]*F_0l[9]+alphaDrSurf_l[3]*F_0l[8])+0.6123724356957944*(alphaDrSurf_l[0]*F_0l[7]+alphaDrSurf_l[1]*F_0l[6]+alphaDrSurf_l[2]*F_0l[5])+0.3535533905932737*alphaDrSurf_l[0]*F_0l[4]+0.6123724356957944*F_0l[3]*alphaDrSurf_l[3]+0.3535533905932737*(F_0l[0]*alphaDrSurf_l[3]+F_0l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*F_0l[2]); 
  Ghat_G_1_l[0] = 0.7905694150420947*(alphaDrSurf_l[3]*G_1l[11]+alphaDrSurf_l[2]*G_1l[10]+alphaDrSurf_l[1]*G_1l[9]+alphaDrSurf_l[0]*G_1l[8])+0.6123724356957944*(alphaDrSurf_l[3]*G_1l[7]+alphaDrSurf_l[2]*G_1l[6]+alphaDrSurf_l[1]*G_1l[5])+0.3535533905932737*alphaDrSurf_l[3]*G_1l[4]+0.6123724356957944*alphaDrSurf_l[0]*G_1l[3]+0.3535533905932737*(G_1l[2]*alphaDrSurf_l[2]+G_1l[1]*alphaDrSurf_l[1]+G_1l[0]*alphaDrSurf_l[0]); 
  Ghat_G_1_l[1] = 0.7905694150420947*(alphaDrSurf_l[2]*G_1l[11]+alphaDrSurf_l[3]*G_1l[10]+alphaDrSurf_l[0]*G_1l[9]+alphaDrSurf_l[1]*G_1l[8])+0.6123724356957944*(alphaDrSurf_l[2]*G_1l[7]+alphaDrSurf_l[3]*G_1l[6]+alphaDrSurf_l[0]*G_1l[5])+0.3535533905932737*(alphaDrSurf_l[2]*G_1l[4]+G_1l[2]*alphaDrSurf_l[3])+0.6123724356957944*alphaDrSurf_l[1]*G_1l[3]+0.3535533905932737*(G_1l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*G_1l[1]); 
  Ghat_G_1_l[2] = 0.7905694150420947*(alphaDrSurf_l[1]*G_1l[11]+alphaDrSurf_l[0]*G_1l[10]+alphaDrSurf_l[3]*G_1l[9]+alphaDrSurf_l[2]*G_1l[8])+0.6123724356957944*(alphaDrSurf_l[1]*G_1l[7]+alphaDrSurf_l[0]*G_1l[6]+alphaDrSurf_l[3]*G_1l[5])+0.3535533905932737*(alphaDrSurf_l[1]*G_1l[4]+G_1l[1]*alphaDrSurf_l[3])+0.6123724356957944*alphaDrSurf_l[2]*G_1l[3]+0.3535533905932737*(G_1l[0]*alphaDrSurf_l[2]+alphaDrSurf_l[0]*G_1l[2]); 
  Ghat_G_1_l[3] = 0.7905694150420947*(alphaDrSurf_l[0]*G_1l[11]+alphaDrSurf_l[1]*G_1l[10]+alphaDrSurf_l[2]*G_1l[9]+alphaDrSurf_l[3]*G_1l[8])+0.6123724356957944*(alphaDrSurf_l[0]*G_1l[7]+alphaDrSurf_l[1]*G_1l[6]+alphaDrSurf_l[2]*G_1l[5])+0.3535533905932737*alphaDrSurf_l[0]*G_1l[4]+0.6123724356957944*G_1l[3]*alphaDrSurf_l[3]+0.3535533905932737*(G_1l[0]*alphaDrSurf_l[3]+G_1l[1]*alphaDrSurf_l[2]+alphaDrSurf_l[1]*G_1l[2]); 

  } 
  out_F_0[0] += (0.7071067811865475*Ghat_F_0_r[0]-0.7071067811865475*Ghat_F_0_l[0])*dv1par; 
  out_F_0[1] += (0.7071067811865475*Ghat_F_0_r[1]-0.7071067811865475*Ghat_F_0_l[1])*dv1par; 
  out_F_0[2] += (0.7071067811865475*Ghat_F_0_r[2]-0.7071067811865475*Ghat_F_0_l[2])*dv1par; 
  out_F_0[3] += 1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dv1par; 
  out_F_0[4] += (0.7071067811865475*Ghat_F_0_r[3]-0.7071067811865475*Ghat_F_0_l[3])*dv1par; 
  out_F_0[5] += 1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dv1par; 
  out_F_0[6] += 1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2])*dv1par; 
  out_F_0[7] += 1.224744871391589*(Ghat_F_0_r[3]+Ghat_F_0_l[3])*dv1par; 
  out_F_0[8] += (1.58113883008419*Ghat_F_0_r[0]-1.58113883008419*Ghat_F_0_l[0])*dv1par; 
  out_F_0[9] += (1.58113883008419*Ghat_F_0_r[1]-1.58113883008419*Ghat_F_0_l[1])*dv1par; 
  out_F_0[10] += (1.58113883008419*Ghat_F_0_r[2]-1.58113883008419*Ghat_F_0_l[2])*dv1par; 
  out_F_0[11] += (1.58113883008419*Ghat_F_0_r[3]-1.58113883008419*Ghat_F_0_l[3])*dv1par; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_r[0]-0.7071067811865475*Ghat_G_1_l[0])*dv1par; 
  out_G_1[1] += (0.7071067811865475*Ghat_G_1_r[1]-0.7071067811865475*Ghat_G_1_l[1])*dv1par; 
  out_G_1[2] += (0.7071067811865475*Ghat_G_1_r[2]-0.7071067811865475*Ghat_G_1_l[2])*dv1par; 
  out_G_1[3] += 1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dv1par; 
  out_G_1[4] += (0.7071067811865475*Ghat_G_1_r[3]-0.7071067811865475*Ghat_G_1_l[3])*dv1par; 
  out_G_1[5] += 1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dv1par; 
  out_G_1[6] += 1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dv1par; 
  out_G_1[7] += 1.224744871391589*(Ghat_G_1_r[3]+Ghat_G_1_l[3])*dv1par; 
  out_G_1[8] += (1.58113883008419*Ghat_G_1_r[0]-1.58113883008419*Ghat_G_1_l[0])*dv1par; 
  out_G_1[9] += (1.58113883008419*Ghat_G_1_r[1]-1.58113883008419*Ghat_G_1_l[1])*dv1par; 
  out_G_1[10] += (1.58113883008419*Ghat_G_1_r[2]-1.58113883008419*Ghat_G_1_l[2])*dv1par; 
  out_G_1[11] += (1.58113883008419*Ghat_G_1_r[3]-1.58113883008419*Ghat_G_1_l[3])*dv1par; 
} 
