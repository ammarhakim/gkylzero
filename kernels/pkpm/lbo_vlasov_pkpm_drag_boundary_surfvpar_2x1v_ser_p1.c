#include <gkyl_lbo_vlasov_pkpm_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nu, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:       Cell-center coordinates. 
  // dxv[3]:     Cell spacing. 
  // nu:          Collisionality. 
  // fSkin/fEdge: Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in skin cell/last edge cell 
  // out:         Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0Skin = &fSkin[0]; 
  const double *G_1Skin = &fSkin[12]; 
  const double *F_0Edge = &fEdge[0]; 
  const double *G_1Edge = &fEdge[12]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 

  double alphaDrSurf[4] = {0.0}; 
  double Ghat_F_0[4] = {0.0}; 
  double Ghat_G_1[4] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nu[0]*wvpar+0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar+0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar+0.5*nu[2]*dvpar; 
  alphaDrSurf[3] = nu[3]*wvpar+0.5*nu[3]*dvpar; 

  Ghat_F_0[0] = 0.7905694150420947*alphaDrSurf[3]*F_0Skin[11]+0.7905694150420948*alphaDrSurf[2]*F_0Skin[10]+0.7905694150420948*alphaDrSurf[1]*F_0Skin[9]+0.7905694150420947*alphaDrSurf[0]*F_0Skin[8]+0.6123724356957944*alphaDrSurf[3]*F_0Skin[7]+0.6123724356957944*alphaDrSurf[2]*F_0Skin[6]+0.6123724356957944*alphaDrSurf[1]*F_0Skin[5]+0.3535533905932737*alphaDrSurf[3]*F_0Skin[4]+0.6123724356957944*alphaDrSurf[0]*F_0Skin[3]+0.3535533905932737*F_0Skin[2]*alphaDrSurf[2]+0.3535533905932737*F_0Skin[1]*alphaDrSurf[1]+0.3535533905932737*F_0Skin[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.7905694150420947*alphaDrSurf[2]*F_0Skin[11]+0.7905694150420948*alphaDrSurf[3]*F_0Skin[10]+0.7905694150420948*alphaDrSurf[0]*F_0Skin[9]+0.7905694150420947*alphaDrSurf[1]*F_0Skin[8]+0.6123724356957944*alphaDrSurf[2]*F_0Skin[7]+0.6123724356957944*alphaDrSurf[3]*F_0Skin[6]+0.6123724356957944*alphaDrSurf[0]*F_0Skin[5]+0.3535533905932737*alphaDrSurf[2]*F_0Skin[4]+0.3535533905932737*F_0Skin[2]*alphaDrSurf[3]+0.6123724356957944*alphaDrSurf[1]*F_0Skin[3]+0.3535533905932737*F_0Skin[0]*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0]*F_0Skin[1]; 
  Ghat_F_0[2] = 0.7905694150420947*alphaDrSurf[1]*F_0Skin[11]+0.7905694150420948*alphaDrSurf[0]*F_0Skin[10]+0.7905694150420948*alphaDrSurf[3]*F_0Skin[9]+0.7905694150420947*alphaDrSurf[2]*F_0Skin[8]+0.6123724356957944*alphaDrSurf[1]*F_0Skin[7]+0.6123724356957944*alphaDrSurf[0]*F_0Skin[6]+0.6123724356957944*alphaDrSurf[3]*F_0Skin[5]+0.3535533905932737*alphaDrSurf[1]*F_0Skin[4]+0.3535533905932737*F_0Skin[1]*alphaDrSurf[3]+0.6123724356957944*alphaDrSurf[2]*F_0Skin[3]+0.3535533905932737*F_0Skin[0]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0]*F_0Skin[2]; 
  Ghat_F_0[3] = 0.7905694150420947*alphaDrSurf[0]*F_0Skin[11]+0.7905694150420948*alphaDrSurf[1]*F_0Skin[10]+0.7905694150420948*alphaDrSurf[2]*F_0Skin[9]+0.7905694150420947*alphaDrSurf[3]*F_0Skin[8]+0.6123724356957944*alphaDrSurf[0]*F_0Skin[7]+0.6123724356957944*alphaDrSurf[1]*F_0Skin[6]+0.6123724356957944*alphaDrSurf[2]*F_0Skin[5]+0.3535533905932737*alphaDrSurf[0]*F_0Skin[4]+0.6123724356957944*F_0Skin[3]*alphaDrSurf[3]+0.3535533905932737*F_0Skin[0]*alphaDrSurf[3]+0.3535533905932737*F_0Skin[1]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]*F_0Skin[2]; 
  Ghat_G_1[0] = 0.7905694150420947*alphaDrSurf[3]*G_1Skin[11]+0.7905694150420948*alphaDrSurf[2]*G_1Skin[10]+0.7905694150420948*alphaDrSurf[1]*G_1Skin[9]+0.7905694150420947*alphaDrSurf[0]*G_1Skin[8]+0.6123724356957944*alphaDrSurf[3]*G_1Skin[7]+0.6123724356957944*alphaDrSurf[2]*G_1Skin[6]+0.6123724356957944*alphaDrSurf[1]*G_1Skin[5]+0.3535533905932737*alphaDrSurf[3]*G_1Skin[4]+0.6123724356957944*alphaDrSurf[0]*G_1Skin[3]+0.3535533905932737*G_1Skin[2]*alphaDrSurf[2]+0.3535533905932737*G_1Skin[1]*alphaDrSurf[1]+0.3535533905932737*G_1Skin[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.7905694150420947*alphaDrSurf[2]*G_1Skin[11]+0.7905694150420948*alphaDrSurf[3]*G_1Skin[10]+0.7905694150420948*alphaDrSurf[0]*G_1Skin[9]+0.7905694150420947*alphaDrSurf[1]*G_1Skin[8]+0.6123724356957944*alphaDrSurf[2]*G_1Skin[7]+0.6123724356957944*alphaDrSurf[3]*G_1Skin[6]+0.6123724356957944*alphaDrSurf[0]*G_1Skin[5]+0.3535533905932737*alphaDrSurf[2]*G_1Skin[4]+0.3535533905932737*G_1Skin[2]*alphaDrSurf[3]+0.6123724356957944*alphaDrSurf[1]*G_1Skin[3]+0.3535533905932737*G_1Skin[0]*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0]*G_1Skin[1]; 
  Ghat_G_1[2] = 0.7905694150420947*alphaDrSurf[1]*G_1Skin[11]+0.7905694150420948*alphaDrSurf[0]*G_1Skin[10]+0.7905694150420948*alphaDrSurf[3]*G_1Skin[9]+0.7905694150420947*alphaDrSurf[2]*G_1Skin[8]+0.6123724356957944*alphaDrSurf[1]*G_1Skin[7]+0.6123724356957944*alphaDrSurf[0]*G_1Skin[6]+0.6123724356957944*alphaDrSurf[3]*G_1Skin[5]+0.3535533905932737*alphaDrSurf[1]*G_1Skin[4]+0.3535533905932737*G_1Skin[1]*alphaDrSurf[3]+0.6123724356957944*alphaDrSurf[2]*G_1Skin[3]+0.3535533905932737*G_1Skin[0]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0]*G_1Skin[2]; 
  Ghat_G_1[3] = 0.7905694150420947*alphaDrSurf[0]*G_1Skin[11]+0.7905694150420948*alphaDrSurf[1]*G_1Skin[10]+0.7905694150420948*alphaDrSurf[2]*G_1Skin[9]+0.7905694150420947*alphaDrSurf[3]*G_1Skin[8]+0.6123724356957944*alphaDrSurf[0]*G_1Skin[7]+0.6123724356957944*alphaDrSurf[1]*G_1Skin[6]+0.6123724356957944*alphaDrSurf[2]*G_1Skin[5]+0.3535533905932737*alphaDrSurf[0]*G_1Skin[4]+0.6123724356957944*G_1Skin[3]*alphaDrSurf[3]+0.3535533905932737*G_1Skin[0]*alphaDrSurf[3]+0.3535533905932737*G_1Skin[1]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]*G_1Skin[2]; 

  out_F_0[0] += 0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += 0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += 0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += 0.7071067811865475*Ghat_F_0[3]*dv1par; 
  out_F_0[5] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += 1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[8] += 1.58113883008419*Ghat_F_0[0]*dv1par; 
  out_F_0[9] += 1.58113883008419*Ghat_F_0[1]*dv1par; 
  out_F_0[10] += 1.58113883008419*Ghat_F_0[2]*dv1par; 
  out_F_0[11] += 1.58113883008419*Ghat_F_0[3]*dv1par; 
  out_G_1[0] += 0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += 0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += 0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += 0.7071067811865475*Ghat_G_1[3]*dv1par; 
  out_G_1[5] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += 1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[8] += 1.58113883008419*Ghat_G_1[0]*dv1par; 
  out_G_1[9] += 1.58113883008419*Ghat_G_1[1]*dv1par; 
  out_G_1[10] += 1.58113883008419*Ghat_G_1[2]*dv1par; 
  out_G_1[11] += 1.58113883008419*Ghat_G_1[3]*dv1par; 

  } else { 

  alphaDrSurf[0] = nu[0]*wvpar-0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar-0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar-0.5*nu[2]*dvpar; 
  alphaDrSurf[3] = nu[3]*wvpar-0.5*nu[3]*dvpar; 

  Ghat_F_0[0] = 0.7905694150420947*alphaDrSurf[3]*F_0Skin[11]+0.7905694150420948*alphaDrSurf[2]*F_0Skin[10]+0.7905694150420948*alphaDrSurf[1]*F_0Skin[9]+0.7905694150420947*alphaDrSurf[0]*F_0Skin[8]-0.6123724356957944*alphaDrSurf[3]*F_0Skin[7]-0.6123724356957944*alphaDrSurf[2]*F_0Skin[6]-0.6123724356957944*alphaDrSurf[1]*F_0Skin[5]+0.3535533905932737*alphaDrSurf[3]*F_0Skin[4]-0.6123724356957944*alphaDrSurf[0]*F_0Skin[3]+0.3535533905932737*F_0Skin[2]*alphaDrSurf[2]+0.3535533905932737*F_0Skin[1]*alphaDrSurf[1]+0.3535533905932737*F_0Skin[0]*alphaDrSurf[0]; 
  Ghat_F_0[1] = 0.7905694150420947*alphaDrSurf[2]*F_0Skin[11]+0.7905694150420948*alphaDrSurf[3]*F_0Skin[10]+0.7905694150420948*alphaDrSurf[0]*F_0Skin[9]+0.7905694150420947*alphaDrSurf[1]*F_0Skin[8]-0.6123724356957944*alphaDrSurf[2]*F_0Skin[7]-0.6123724356957944*alphaDrSurf[3]*F_0Skin[6]-0.6123724356957944*alphaDrSurf[0]*F_0Skin[5]+0.3535533905932737*alphaDrSurf[2]*F_0Skin[4]+0.3535533905932737*F_0Skin[2]*alphaDrSurf[3]-0.6123724356957944*alphaDrSurf[1]*F_0Skin[3]+0.3535533905932737*F_0Skin[0]*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0]*F_0Skin[1]; 
  Ghat_F_0[2] = 0.7905694150420947*alphaDrSurf[1]*F_0Skin[11]+0.7905694150420948*alphaDrSurf[0]*F_0Skin[10]+0.7905694150420948*alphaDrSurf[3]*F_0Skin[9]+0.7905694150420947*alphaDrSurf[2]*F_0Skin[8]-0.6123724356957944*alphaDrSurf[1]*F_0Skin[7]-0.6123724356957944*alphaDrSurf[0]*F_0Skin[6]-0.6123724356957944*alphaDrSurf[3]*F_0Skin[5]+0.3535533905932737*alphaDrSurf[1]*F_0Skin[4]+0.3535533905932737*F_0Skin[1]*alphaDrSurf[3]-0.6123724356957944*alphaDrSurf[2]*F_0Skin[3]+0.3535533905932737*F_0Skin[0]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0]*F_0Skin[2]; 
  Ghat_F_0[3] = 0.7905694150420947*alphaDrSurf[0]*F_0Skin[11]+0.7905694150420948*alphaDrSurf[1]*F_0Skin[10]+0.7905694150420948*alphaDrSurf[2]*F_0Skin[9]+0.7905694150420947*alphaDrSurf[3]*F_0Skin[8]-0.6123724356957944*alphaDrSurf[0]*F_0Skin[7]-0.6123724356957944*alphaDrSurf[1]*F_0Skin[6]-0.6123724356957944*alphaDrSurf[2]*F_0Skin[5]+0.3535533905932737*alphaDrSurf[0]*F_0Skin[4]-0.6123724356957944*F_0Skin[3]*alphaDrSurf[3]+0.3535533905932737*F_0Skin[0]*alphaDrSurf[3]+0.3535533905932737*F_0Skin[1]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]*F_0Skin[2]; 
  Ghat_G_1[0] = 0.7905694150420947*alphaDrSurf[3]*G_1Skin[11]+0.7905694150420948*alphaDrSurf[2]*G_1Skin[10]+0.7905694150420948*alphaDrSurf[1]*G_1Skin[9]+0.7905694150420947*alphaDrSurf[0]*G_1Skin[8]-0.6123724356957944*alphaDrSurf[3]*G_1Skin[7]-0.6123724356957944*alphaDrSurf[2]*G_1Skin[6]-0.6123724356957944*alphaDrSurf[1]*G_1Skin[5]+0.3535533905932737*alphaDrSurf[3]*G_1Skin[4]-0.6123724356957944*alphaDrSurf[0]*G_1Skin[3]+0.3535533905932737*G_1Skin[2]*alphaDrSurf[2]+0.3535533905932737*G_1Skin[1]*alphaDrSurf[1]+0.3535533905932737*G_1Skin[0]*alphaDrSurf[0]; 
  Ghat_G_1[1] = 0.7905694150420947*alphaDrSurf[2]*G_1Skin[11]+0.7905694150420948*alphaDrSurf[3]*G_1Skin[10]+0.7905694150420948*alphaDrSurf[0]*G_1Skin[9]+0.7905694150420947*alphaDrSurf[1]*G_1Skin[8]-0.6123724356957944*alphaDrSurf[2]*G_1Skin[7]-0.6123724356957944*alphaDrSurf[3]*G_1Skin[6]-0.6123724356957944*alphaDrSurf[0]*G_1Skin[5]+0.3535533905932737*alphaDrSurf[2]*G_1Skin[4]+0.3535533905932737*G_1Skin[2]*alphaDrSurf[3]-0.6123724356957944*alphaDrSurf[1]*G_1Skin[3]+0.3535533905932737*G_1Skin[0]*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0]*G_1Skin[1]; 
  Ghat_G_1[2] = 0.7905694150420947*alphaDrSurf[1]*G_1Skin[11]+0.7905694150420948*alphaDrSurf[0]*G_1Skin[10]+0.7905694150420948*alphaDrSurf[3]*G_1Skin[9]+0.7905694150420947*alphaDrSurf[2]*G_1Skin[8]-0.6123724356957944*alphaDrSurf[1]*G_1Skin[7]-0.6123724356957944*alphaDrSurf[0]*G_1Skin[6]-0.6123724356957944*alphaDrSurf[3]*G_1Skin[5]+0.3535533905932737*alphaDrSurf[1]*G_1Skin[4]+0.3535533905932737*G_1Skin[1]*alphaDrSurf[3]-0.6123724356957944*alphaDrSurf[2]*G_1Skin[3]+0.3535533905932737*G_1Skin[0]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0]*G_1Skin[2]; 
  Ghat_G_1[3] = 0.7905694150420947*alphaDrSurf[0]*G_1Skin[11]+0.7905694150420948*alphaDrSurf[1]*G_1Skin[10]+0.7905694150420948*alphaDrSurf[2]*G_1Skin[9]+0.7905694150420947*alphaDrSurf[3]*G_1Skin[8]-0.6123724356957944*alphaDrSurf[0]*G_1Skin[7]-0.6123724356957944*alphaDrSurf[1]*G_1Skin[6]-0.6123724356957944*alphaDrSurf[2]*G_1Skin[5]+0.3535533905932737*alphaDrSurf[0]*G_1Skin[4]-0.6123724356957944*G_1Skin[3]*alphaDrSurf[3]+0.3535533905932737*G_1Skin[0]*alphaDrSurf[3]+0.3535533905932737*G_1Skin[1]*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]*G_1Skin[2]; 

  out_F_0[0] += -0.7071067811865475*Ghat_F_0[0]*dv1par; 
  out_F_0[1] += -0.7071067811865475*Ghat_F_0[1]*dv1par; 
  out_F_0[2] += -0.7071067811865475*Ghat_F_0[2]*dv1par; 
  out_F_0[3] += 1.224744871391589*Ghat_F_0[0]*dv1par; 
  out_F_0[4] += -0.7071067811865475*Ghat_F_0[3]*dv1par; 
  out_F_0[5] += 1.224744871391589*Ghat_F_0[1]*dv1par; 
  out_F_0[6] += 1.224744871391589*Ghat_F_0[2]*dv1par; 
  out_F_0[7] += 1.224744871391589*Ghat_F_0[3]*dv1par; 
  out_F_0[8] += -1.58113883008419*Ghat_F_0[0]*dv1par; 
  out_F_0[9] += -1.58113883008419*Ghat_F_0[1]*dv1par; 
  out_F_0[10] += -1.58113883008419*Ghat_F_0[2]*dv1par; 
  out_F_0[11] += -1.58113883008419*Ghat_F_0[3]*dv1par; 
  out_G_1[0] += -0.7071067811865475*Ghat_G_1[0]*dv1par; 
  out_G_1[1] += -0.7071067811865475*Ghat_G_1[1]*dv1par; 
  out_G_1[2] += -0.7071067811865475*Ghat_G_1[2]*dv1par; 
  out_G_1[3] += 1.224744871391589*Ghat_G_1[0]*dv1par; 
  out_G_1[4] += -0.7071067811865475*Ghat_G_1[3]*dv1par; 
  out_G_1[5] += 1.224744871391589*Ghat_G_1[1]*dv1par; 
  out_G_1[6] += 1.224744871391589*Ghat_G_1[2]*dv1par; 
  out_G_1[7] += 1.224744871391589*Ghat_G_1[3]*dv1par; 
  out_G_1[8] += -1.58113883008419*Ghat_G_1[0]*dv1par; 
  out_G_1[9] += -1.58113883008419*Ghat_G_1[1]*dv1par; 
  out_G_1[10] += -1.58113883008419*Ghat_G_1[2]*dv1par; 
  out_G_1[11] += -1.58113883008419*Ghat_G_1[3]*dv1par; 

  }

  return 0.;

} 
