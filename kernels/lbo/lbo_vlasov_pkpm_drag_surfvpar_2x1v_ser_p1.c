#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_drag_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, const double *nu, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nu:         collisionalities added (self and cross species collisionalities). 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 

  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 

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

  double Ghat_r[4]; 
  double Ghat_l[4]; 
  if (wvpar>0) { 

  Ghat_r[0] = 0.7905694150420947*(alphaDrSurf_r[3]*fr[11]+alphaDrSurf_r[2]*fr[10]+alphaDrSurf_r[1]*fr[9]+alphaDrSurf_r[0]*fr[8])-0.6123724356957944*(alphaDrSurf_r[3]*fr[7]+alphaDrSurf_r[2]*fr[6]+alphaDrSurf_r[1]*fr[5])+0.3535533905932737*alphaDrSurf_r[3]*fr[4]-0.6123724356957944*alphaDrSurf_r[0]*fr[3]+0.3535533905932737*(alphaDrSurf_r[2]*fr[2]+alphaDrSurf_r[1]*fr[1]+alphaDrSurf_r[0]*fr[0]); 
  Ghat_r[1] = 0.7905694150420947*(alphaDrSurf_r[2]*fr[11]+alphaDrSurf_r[3]*fr[10]+alphaDrSurf_r[0]*fr[9]+alphaDrSurf_r[1]*fr[8])-0.6123724356957944*(alphaDrSurf_r[2]*fr[7]+alphaDrSurf_r[3]*fr[6]+alphaDrSurf_r[0]*fr[5])+0.3535533905932737*alphaDrSurf_r[2]*fr[4]-0.6123724356957944*alphaDrSurf_r[1]*fr[3]+0.3535533905932737*(fr[2]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*fr[1]+fr[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.7905694150420947*(alphaDrSurf_r[1]*fr[11]+alphaDrSurf_r[0]*fr[10]+alphaDrSurf_r[3]*fr[9]+alphaDrSurf_r[2]*fr[8])-0.6123724356957944*(alphaDrSurf_r[1]*fr[7]+alphaDrSurf_r[0]*fr[6]+alphaDrSurf_r[3]*fr[5])+0.3535533905932737*alphaDrSurf_r[1]*fr[4]-0.6123724356957944*alphaDrSurf_r[2]*fr[3]+0.3535533905932737*(fr[1]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*fr[2]+fr[0]*alphaDrSurf_r[2]); 
  Ghat_r[3] = 0.7905694150420947*(alphaDrSurf_r[0]*fr[11]+alphaDrSurf_r[1]*fr[10]+alphaDrSurf_r[2]*fr[9]+alphaDrSurf_r[3]*fr[8])-0.6123724356957944*(alphaDrSurf_r[0]*fr[7]+alphaDrSurf_r[1]*fr[6]+alphaDrSurf_r[2]*fr[5])+0.3535533905932737*alphaDrSurf_r[0]*fr[4]-0.6123724356957944*alphaDrSurf_r[3]*fr[3]+0.3535533905932737*(fr[0]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*fr[2]+fr[1]*alphaDrSurf_r[2]); 

  Ghat_l[0] = 0.7905694150420947*(alphaDrSurf_l[3]*fc[11]+alphaDrSurf_l[2]*fc[10]+alphaDrSurf_l[1]*fc[9]+alphaDrSurf_l[0]*fc[8])-0.6123724356957944*(alphaDrSurf_l[3]*fc[7]+alphaDrSurf_l[2]*fc[6]+alphaDrSurf_l[1]*fc[5])+0.3535533905932737*alphaDrSurf_l[3]*fc[4]-0.6123724356957944*alphaDrSurf_l[0]*fc[3]+0.3535533905932737*(alphaDrSurf_l[2]*fc[2]+alphaDrSurf_l[1]*fc[1]+alphaDrSurf_l[0]*fc[0]); 
  Ghat_l[1] = 0.7905694150420947*(alphaDrSurf_l[2]*fc[11]+alphaDrSurf_l[3]*fc[10]+alphaDrSurf_l[0]*fc[9]+alphaDrSurf_l[1]*fc[8])-0.6123724356957944*(alphaDrSurf_l[2]*fc[7]+alphaDrSurf_l[3]*fc[6]+alphaDrSurf_l[0]*fc[5])+0.3535533905932737*alphaDrSurf_l[2]*fc[4]-0.6123724356957944*alphaDrSurf_l[1]*fc[3]+0.3535533905932737*(fc[2]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*fc[1]+fc[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.7905694150420947*(alphaDrSurf_l[1]*fc[11]+alphaDrSurf_l[0]*fc[10]+alphaDrSurf_l[3]*fc[9]+alphaDrSurf_l[2]*fc[8])-0.6123724356957944*(alphaDrSurf_l[1]*fc[7]+alphaDrSurf_l[0]*fc[6]+alphaDrSurf_l[3]*fc[5])+0.3535533905932737*alphaDrSurf_l[1]*fc[4]-0.6123724356957944*alphaDrSurf_l[2]*fc[3]+0.3535533905932737*(fc[1]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*fc[2]+fc[0]*alphaDrSurf_l[2]); 
  Ghat_l[3] = 0.7905694150420947*(alphaDrSurf_l[0]*fc[11]+alphaDrSurf_l[1]*fc[10]+alphaDrSurf_l[2]*fc[9]+alphaDrSurf_l[3]*fc[8])-0.6123724356957944*(alphaDrSurf_l[0]*fc[7]+alphaDrSurf_l[1]*fc[6]+alphaDrSurf_l[2]*fc[5])+0.3535533905932737*alphaDrSurf_l[0]*fc[4]-0.6123724356957944*alphaDrSurf_l[3]*fc[3]+0.3535533905932737*(fc[0]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*fc[2]+fc[1]*alphaDrSurf_l[2]); 

  } else { 

  Ghat_r[0] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_r[3]*fc[11]+alphaDrSurf_r[2]*fc[10]+alphaDrSurf_r[1]*fc[9]+alphaDrSurf_r[0]*fc[8])+7.348469228349534*(alphaDrSurf_r[3]*fc[7]+alphaDrSurf_r[2]*fc[6]+alphaDrSurf_r[1]*fc[5])+4.242640687119286*alphaDrSurf_r[3]*fc[4]+7.348469228349534*alphaDrSurf_r[0]*fc[3]+4.242640687119286*(alphaDrSurf_r[2]*fc[2]+alphaDrSurf_r[1]*fc[1]+alphaDrSurf_r[0]*fc[0])); 
  Ghat_r[1] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_r[2]*fc[11]+alphaDrSurf_r[3]*fc[10]+alphaDrSurf_r[0]*fc[9]+alphaDrSurf_r[1]*fc[8])+7.348469228349534*(alphaDrSurf_r[2]*fc[7]+alphaDrSurf_r[3]*fc[6]+alphaDrSurf_r[0]*fc[5])+4.242640687119286*alphaDrSurf_r[2]*fc[4]+7.348469228349534*alphaDrSurf_r[1]*fc[3]+4.242640687119286*(fc[2]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*fc[1]+fc[0]*alphaDrSurf_r[1])); 
  Ghat_r[2] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_r[1]*fc[11]+alphaDrSurf_r[0]*fc[10]+alphaDrSurf_r[3]*fc[9]+alphaDrSurf_r[2]*fc[8])+7.348469228349534*(alphaDrSurf_r[1]*fc[7]+alphaDrSurf_r[0]*fc[6]+alphaDrSurf_r[3]*fc[5])+4.242640687119286*alphaDrSurf_r[1]*fc[4]+7.348469228349534*alphaDrSurf_r[2]*fc[3]+4.242640687119286*(fc[1]*alphaDrSurf_r[3]+alphaDrSurf_r[0]*fc[2]+fc[0]*alphaDrSurf_r[2])); 
  Ghat_r[3] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_r[0]*fc[11]+alphaDrSurf_r[1]*fc[10]+alphaDrSurf_r[2]*fc[9]+alphaDrSurf_r[3]*fc[8])+7.348469228349534*(alphaDrSurf_r[0]*fc[7]+alphaDrSurf_r[1]*fc[6]+alphaDrSurf_r[2]*fc[5])+4.242640687119286*alphaDrSurf_r[0]*fc[4]+7.348469228349534*alphaDrSurf_r[3]*fc[3]+4.242640687119286*(fc[0]*alphaDrSurf_r[3]+alphaDrSurf_r[1]*fc[2]+fc[1]*alphaDrSurf_r[2])); 

  Ghat_l[0] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_l[3]*fl[11]+alphaDrSurf_l[2]*fl[10]+alphaDrSurf_l[1]*fl[9]+alphaDrSurf_l[0]*fl[8])+7.348469228349534*(alphaDrSurf_l[3]*fl[7]+alphaDrSurf_l[2]*fl[6]+alphaDrSurf_l[1]*fl[5])+4.242640687119286*alphaDrSurf_l[3]*fl[4]+7.348469228349534*alphaDrSurf_l[0]*fl[3]+4.242640687119286*(alphaDrSurf_l[2]*fl[2]+alphaDrSurf_l[1]*fl[1]+alphaDrSurf_l[0]*fl[0])); 
  Ghat_l[1] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_l[2]*fl[11]+alphaDrSurf_l[3]*fl[10]+alphaDrSurf_l[0]*fl[9]+alphaDrSurf_l[1]*fl[8])+7.348469228349534*(alphaDrSurf_l[2]*fl[7]+alphaDrSurf_l[3]*fl[6]+alphaDrSurf_l[0]*fl[5])+4.242640687119286*alphaDrSurf_l[2]*fl[4]+7.348469228349534*alphaDrSurf_l[1]*fl[3]+4.242640687119286*(fl[2]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*fl[1]+fl[0]*alphaDrSurf_l[1])); 
  Ghat_l[2] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_l[1]*fl[11]+alphaDrSurf_l[0]*fl[10]+alphaDrSurf_l[3]*fl[9]+alphaDrSurf_l[2]*fl[8])+7.348469228349534*(alphaDrSurf_l[1]*fl[7]+alphaDrSurf_l[0]*fl[6]+alphaDrSurf_l[3]*fl[5])+4.242640687119286*alphaDrSurf_l[1]*fl[4]+7.348469228349534*alphaDrSurf_l[2]*fl[3]+4.242640687119286*(fl[1]*alphaDrSurf_l[3]+alphaDrSurf_l[0]*fl[2]+fl[0]*alphaDrSurf_l[2])); 
  Ghat_l[3] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_l[0]*fl[11]+alphaDrSurf_l[1]*fl[10]+alphaDrSurf_l[2]*fl[9]+alphaDrSurf_l[3]*fl[8])+7.348469228349534*(alphaDrSurf_l[0]*fl[7]+alphaDrSurf_l[1]*fl[6]+alphaDrSurf_l[2]*fl[5])+4.242640687119286*alphaDrSurf_l[0]*fl[4]+7.348469228349534*alphaDrSurf_l[3]*fl[3]+4.242640687119286*(fl[0]*alphaDrSurf_l[3]+alphaDrSurf_l[1]*fl[2]+fl[1]*alphaDrSurf_l[2])); 

  } 
  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*dv1par; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*dv1par; 
  out[2] += (0.7071067811865475*Ghat_r[2]-0.7071067811865475*Ghat_l[2])*dv1par; 
  out[3] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[4] += (0.7071067811865475*Ghat_r[3]-0.7071067811865475*Ghat_l[3])*dv1par; 
  out[5] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv1par; 
  out[6] += 1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv1par; 
  out[7] += 1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv1par; 
  out[8] += (1.58113883008419*Ghat_r[0]-1.58113883008419*Ghat_l[0])*dv1par; 
  out[9] += (1.58113883008419*Ghat_r[1]-1.58113883008419*Ghat_l[1])*dv1par; 
  out[10] += (1.58113883008419*Ghat_r[2]-1.58113883008419*Ghat_l[2])*dv1par; 
  out[11] += (1.58113883008419*Ghat_r[3]-1.58113883008419*Ghat_l[3])*dv1par; 
} 
