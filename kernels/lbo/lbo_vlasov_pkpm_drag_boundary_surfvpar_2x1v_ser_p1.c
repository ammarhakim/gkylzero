#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_drag_boundary_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nu, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nu:         collisionality. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 

  double alphaDrSurf[4] = {0.0}; 
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nu[0]*wvpar+0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar+0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar+0.5*nu[2]*dvpar; 
  alphaDrSurf[3] = nu[3]*wvpar+0.5*nu[3]*dvpar; 

  Ghat[0] = 0.7905694150420947*alphaDrSurf[3]*fSkin[11]+0.7905694150420948*alphaDrSurf[2]*fSkin[10]+0.7905694150420948*alphaDrSurf[1]*fSkin[9]+0.7905694150420947*alphaDrSurf[0]*fSkin[8]+0.6123724356957944*alphaDrSurf[3]*fSkin[7]+0.6123724356957944*alphaDrSurf[2]*fSkin[6]+0.6123724356957944*alphaDrSurf[1]*fSkin[5]+0.3535533905932737*alphaDrSurf[3]*fSkin[4]+0.6123724356957944*alphaDrSurf[0]*fSkin[3]+0.3535533905932737*alphaDrSurf[2]*fSkin[2]+0.3535533905932737*alphaDrSurf[1]*fSkin[1]+0.3535533905932737*alphaDrSurf[0]*fSkin[0]; 
  Ghat[1] = 0.7905694150420947*alphaDrSurf[2]*fSkin[11]+0.7905694150420948*alphaDrSurf[3]*fSkin[10]+0.7905694150420948*alphaDrSurf[0]*fSkin[9]+0.7905694150420947*alphaDrSurf[1]*fSkin[8]+0.6123724356957944*alphaDrSurf[2]*fSkin[7]+0.6123724356957944*alphaDrSurf[3]*fSkin[6]+0.6123724356957944*alphaDrSurf[0]*fSkin[5]+0.3535533905932737*alphaDrSurf[2]*fSkin[4]+0.6123724356957944*alphaDrSurf[1]*fSkin[3]+0.3535533905932737*fSkin[2]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[0]*fSkin[1]+0.3535533905932737*fSkin[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.7905694150420947*alphaDrSurf[1]*fSkin[11]+0.7905694150420948*alphaDrSurf[0]*fSkin[10]+0.7905694150420948*alphaDrSurf[3]*fSkin[9]+0.7905694150420947*alphaDrSurf[2]*fSkin[8]+0.6123724356957944*alphaDrSurf[1]*fSkin[7]+0.6123724356957944*alphaDrSurf[0]*fSkin[6]+0.6123724356957944*alphaDrSurf[3]*fSkin[5]+0.3535533905932737*alphaDrSurf[1]*fSkin[4]+0.6123724356957944*alphaDrSurf[2]*fSkin[3]+0.3535533905932737*fSkin[1]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[0]*fSkin[2]+0.3535533905932737*fSkin[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.7905694150420947*alphaDrSurf[0]*fSkin[11]+0.7905694150420948*alphaDrSurf[1]*fSkin[10]+0.7905694150420948*alphaDrSurf[2]*fSkin[9]+0.7905694150420947*alphaDrSurf[3]*fSkin[8]+0.6123724356957944*alphaDrSurf[0]*fSkin[7]+0.6123724356957944*alphaDrSurf[1]*fSkin[6]+0.6123724356957944*alphaDrSurf[2]*fSkin[5]+0.3535533905932737*alphaDrSurf[0]*fSkin[4]+0.6123724356957944*alphaDrSurf[3]*fSkin[3]+0.3535533905932737*fSkin[0]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[1]*fSkin[2]+0.3535533905932737*fSkin[1]*alphaDrSurf[2]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += 0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += 1.224744871391589*Ghat[0]*dv1par; 
  out[4] += 0.7071067811865475*Ghat[3]*dv1par; 
  out[5] += 1.224744871391589*Ghat[1]*dv1par; 
  out[6] += 1.224744871391589*Ghat[2]*dv1par; 
  out[7] += 1.224744871391589*Ghat[3]*dv1par; 
  out[8] += 1.58113883008419*Ghat[0]*dv1par; 
  out[9] += 1.58113883008419*Ghat[1]*dv1par; 
  out[10] += 1.58113883008419*Ghat[2]*dv1par; 
  out[11] += 1.58113883008419*Ghat[3]*dv1par; 

  } else { 

  alphaDrSurf[0] = nu[0]*wvpar-0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar-0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar-0.5*nu[2]*dvpar; 
  alphaDrSurf[3] = nu[3]*wvpar-0.5*nu[3]*dvpar; 

  Ghat[0] = 0.7905694150420947*alphaDrSurf[3]*fSkin[11]+0.7905694150420948*alphaDrSurf[2]*fSkin[10]+0.7905694150420948*alphaDrSurf[1]*fSkin[9]+0.7905694150420947*alphaDrSurf[0]*fSkin[8]-0.6123724356957944*alphaDrSurf[3]*fSkin[7]-0.6123724356957944*alphaDrSurf[2]*fSkin[6]-0.6123724356957944*alphaDrSurf[1]*fSkin[5]+0.3535533905932737*alphaDrSurf[3]*fSkin[4]-0.6123724356957944*alphaDrSurf[0]*fSkin[3]+0.3535533905932737*alphaDrSurf[2]*fSkin[2]+0.3535533905932737*alphaDrSurf[1]*fSkin[1]+0.3535533905932737*alphaDrSurf[0]*fSkin[0]; 
  Ghat[1] = 0.7905694150420947*alphaDrSurf[2]*fSkin[11]+0.7905694150420948*alphaDrSurf[3]*fSkin[10]+0.7905694150420948*alphaDrSurf[0]*fSkin[9]+0.7905694150420947*alphaDrSurf[1]*fSkin[8]-0.6123724356957944*alphaDrSurf[2]*fSkin[7]-0.6123724356957944*alphaDrSurf[3]*fSkin[6]-0.6123724356957944*alphaDrSurf[0]*fSkin[5]+0.3535533905932737*alphaDrSurf[2]*fSkin[4]-0.6123724356957944*alphaDrSurf[1]*fSkin[3]+0.3535533905932737*fSkin[2]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[0]*fSkin[1]+0.3535533905932737*fSkin[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.7905694150420947*alphaDrSurf[1]*fSkin[11]+0.7905694150420948*alphaDrSurf[0]*fSkin[10]+0.7905694150420948*alphaDrSurf[3]*fSkin[9]+0.7905694150420947*alphaDrSurf[2]*fSkin[8]-0.6123724356957944*alphaDrSurf[1]*fSkin[7]-0.6123724356957944*alphaDrSurf[0]*fSkin[6]-0.6123724356957944*alphaDrSurf[3]*fSkin[5]+0.3535533905932737*alphaDrSurf[1]*fSkin[4]-0.6123724356957944*alphaDrSurf[2]*fSkin[3]+0.3535533905932737*fSkin[1]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[0]*fSkin[2]+0.3535533905932737*fSkin[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.7905694150420947*alphaDrSurf[0]*fSkin[11]+0.7905694150420948*alphaDrSurf[1]*fSkin[10]+0.7905694150420948*alphaDrSurf[2]*fSkin[9]+0.7905694150420947*alphaDrSurf[3]*fSkin[8]-0.6123724356957944*alphaDrSurf[0]*fSkin[7]-0.6123724356957944*alphaDrSurf[1]*fSkin[6]-0.6123724356957944*alphaDrSurf[2]*fSkin[5]+0.3535533905932737*alphaDrSurf[0]*fSkin[4]-0.6123724356957944*alphaDrSurf[3]*fSkin[3]+0.3535533905932737*fSkin[0]*alphaDrSurf[3]+0.3535533905932737*alphaDrSurf[1]*fSkin[2]+0.3535533905932737*fSkin[1]*alphaDrSurf[2]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += 1.224744871391589*Ghat[0]*dv1par; 
  out[4] += -0.7071067811865475*Ghat[3]*dv1par; 
  out[5] += 1.224744871391589*Ghat[1]*dv1par; 
  out[6] += 1.224744871391589*Ghat[2]*dv1par; 
  out[7] += 1.224744871391589*Ghat[3]*dv1par; 
  out[8] += -1.58113883008419*Ghat[0]*dv1par; 
  out[9] += -1.58113883008419*Ghat[1]*dv1par; 
  out[10] += -1.58113883008419*Ghat[2]*dv1par; 
  out[11] += -1.58113883008419*Ghat[3]*dv1par; 

  } 
} 
