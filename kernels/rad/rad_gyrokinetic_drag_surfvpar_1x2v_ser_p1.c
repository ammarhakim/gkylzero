#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double *nuField, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // nuField:   2/pi*v*nu(v) field dg representation (v'(v||,mu) in notes)
  // fl/fc/fr:  distribution function in cells 
  // out:       incremented distribution function in cell 

  double rdv2 = 2.0/dxv[1]; 

  double Ghat_r[4] = {0.0}; 
  double Ghat_l[4] = {0.0}; 
  double alphaDrSurf_l[4] = {0.0}; 
  alphaDrSurf_l[0] = 1.58113883008419*nuField[8]-1.224744871391589*nuField[2]+0.7071067811865475*nuField[0]; 
  alphaDrSurf_l[1] = 1.58113883008419*nuField[9]-1.224744871391589*nuField[4]+0.7071067811865475*nuField[1]; 
  alphaDrSurf_l[2] = 1.58113883008419*nuField[10]-1.224744871391589*nuField[6]+0.7071067811865475*nuField[3]; 
  alphaDrSurf_l[3] = 1.58113883008419*nuField[11]-1.224744871391589*nuField[7]+0.7071067811865475*nuField[5]; 

  double alphaDrSurf_r[4] = {0.0}; 
  alphaDrSurf_r[0] = 1.58113883008419*nuField[8]+1.224744871391589*nuField[2]+0.7071067811865475*nuField[0]; 
  alphaDrSurf_r[1] = 1.58113883008419*nuField[9]+1.224744871391589*nuField[4]+0.7071067811865475*nuField[1]; 
  alphaDrSurf_r[2] = 1.58113883008419*nuField[10]+1.224744871391589*nuField[6]+0.7071067811865475*nuField[3]; 
  alphaDrSurf_r[3] = 1.58113883008419*nuField[11]+1.224744871391589*nuField[7]+0.7071067811865475*nuField[5]; 

  if (w[1]>0) {

  Ghat_l[0] = 0.7905694150420947*alphaDrSurf_l[3]*fc[11]+0.7905694150420948*alphaDrSurf_l[2]*fc[10]+0.7905694150420948*alphaDrSurf_l[1]*fc[9]+0.7905694150420947*alphaDrSurf_l[0]*fc[8]-0.6123724356957944*alphaDrSurf_l[3]*fc[7]-0.6123724356957944*alphaDrSurf_l[2]*fc[6]+0.3535533905932737*alphaDrSurf_l[3]*fc[5]-0.6123724356957944*alphaDrSurf_l[1]*fc[4]+0.3535533905932737*alphaDrSurf_l[2]*fc[3]-0.6123724356957944*alphaDrSurf_l[0]*fc[2]+0.3535533905932737*alphaDrSurf_l[1]*fc[1]+0.3535533905932737*alphaDrSurf_l[0]*fc[0]; 
  Ghat_l[1] = 0.7905694150420947*alphaDrSurf_l[2]*fc[11]+0.7905694150420948*alphaDrSurf_l[3]*fc[10]+0.7905694150420948*alphaDrSurf_l[0]*fc[9]+0.7905694150420947*alphaDrSurf_l[1]*fc[8]-0.6123724356957944*alphaDrSurf_l[2]*fc[7]-0.6123724356957944*alphaDrSurf_l[3]*fc[6]+0.3535533905932737*alphaDrSurf_l[2]*fc[5]-0.6123724356957944*alphaDrSurf_l[0]*fc[4]+0.3535533905932737*alphaDrSurf_l[3]*fc[3]-0.6123724356957944*alphaDrSurf_l[1]*fc[2]+0.3535533905932737*alphaDrSurf_l[0]*fc[1]+0.3535533905932737*fc[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = 0.7905694150420947*alphaDrSurf_l[1]*fc[11]+0.7905694150420948*alphaDrSurf_l[0]*fc[10]+0.7905694150420948*alphaDrSurf_l[3]*fc[9]+0.7905694150420947*alphaDrSurf_l[2]*fc[8]-0.6123724356957944*alphaDrSurf_l[1]*fc[7]-0.6123724356957944*alphaDrSurf_l[0]*fc[6]+0.3535533905932737*alphaDrSurf_l[1]*fc[5]-0.6123724356957944*alphaDrSurf_l[3]*fc[4]+0.3535533905932737*alphaDrSurf_l[0]*fc[3]+0.3535533905932737*fc[1]*alphaDrSurf_l[3]-0.6123724356957944*alphaDrSurf_l[2]*fc[2]+0.3535533905932737*fc[0]*alphaDrSurf_l[2]; 
  Ghat_l[3] = 0.7905694150420947*alphaDrSurf_l[0]*fc[11]+0.7905694150420948*alphaDrSurf_l[1]*fc[10]+0.7905694150420948*alphaDrSurf_l[2]*fc[9]+0.7905694150420947*alphaDrSurf_l[3]*fc[8]-0.6123724356957944*alphaDrSurf_l[0]*fc[7]-0.6123724356957944*alphaDrSurf_l[1]*fc[6]+0.3535533905932737*alphaDrSurf_l[0]*fc[5]-0.6123724356957944*alphaDrSurf_l[2]*fc[4]+0.3535533905932737*alphaDrSurf_l[1]*fc[3]-0.6123724356957944*fc[2]*alphaDrSurf_l[3]+0.3535533905932737*fc[0]*alphaDrSurf_l[3]+0.3535533905932737*fc[1]*alphaDrSurf_l[2]; 

  Ghat_r[0] = 0.7905694150420947*alphaDrSurf_r[3]*fr[11]+0.7905694150420948*alphaDrSurf_r[2]*fr[10]+0.7905694150420948*alphaDrSurf_r[1]*fr[9]+0.7905694150420947*alphaDrSurf_r[0]*fr[8]-0.6123724356957944*alphaDrSurf_r[3]*fr[7]-0.6123724356957944*alphaDrSurf_r[2]*fr[6]+0.3535533905932737*alphaDrSurf_r[3]*fr[5]-0.6123724356957944*alphaDrSurf_r[1]*fr[4]+0.3535533905932737*alphaDrSurf_r[2]*fr[3]-0.6123724356957944*alphaDrSurf_r[0]*fr[2]+0.3535533905932737*alphaDrSurf_r[1]*fr[1]+0.3535533905932737*alphaDrSurf_r[0]*fr[0]; 
  Ghat_r[1] = 0.7905694150420947*alphaDrSurf_r[2]*fr[11]+0.7905694150420948*alphaDrSurf_r[3]*fr[10]+0.7905694150420948*alphaDrSurf_r[0]*fr[9]+0.7905694150420947*alphaDrSurf_r[1]*fr[8]-0.6123724356957944*alphaDrSurf_r[2]*fr[7]-0.6123724356957944*alphaDrSurf_r[3]*fr[6]+0.3535533905932737*alphaDrSurf_r[2]*fr[5]-0.6123724356957944*alphaDrSurf_r[0]*fr[4]+0.3535533905932737*alphaDrSurf_r[3]*fr[3]-0.6123724356957944*alphaDrSurf_r[1]*fr[2]+0.3535533905932737*alphaDrSurf_r[0]*fr[1]+0.3535533905932737*fr[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = 0.7905694150420947*alphaDrSurf_r[1]*fr[11]+0.7905694150420948*alphaDrSurf_r[0]*fr[10]+0.7905694150420948*alphaDrSurf_r[3]*fr[9]+0.7905694150420947*alphaDrSurf_r[2]*fr[8]-0.6123724356957944*alphaDrSurf_r[1]*fr[7]-0.6123724356957944*alphaDrSurf_r[0]*fr[6]+0.3535533905932737*alphaDrSurf_r[1]*fr[5]-0.6123724356957944*alphaDrSurf_r[3]*fr[4]+0.3535533905932737*alphaDrSurf_r[0]*fr[3]+0.3535533905932737*fr[1]*alphaDrSurf_r[3]-0.6123724356957944*alphaDrSurf_r[2]*fr[2]+0.3535533905932737*fr[0]*alphaDrSurf_r[2]; 
  Ghat_r[3] = 0.7905694150420947*alphaDrSurf_r[0]*fr[11]+0.7905694150420948*alphaDrSurf_r[1]*fr[10]+0.7905694150420948*alphaDrSurf_r[2]*fr[9]+0.7905694150420947*alphaDrSurf_r[3]*fr[8]-0.6123724356957944*alphaDrSurf_r[0]*fr[7]-0.6123724356957944*alphaDrSurf_r[1]*fr[6]+0.3535533905932737*alphaDrSurf_r[0]*fr[5]-0.6123724356957944*alphaDrSurf_r[2]*fr[4]+0.3535533905932737*alphaDrSurf_r[1]*fr[3]-0.6123724356957944*fr[2]*alphaDrSurf_r[3]+0.3535533905932737*fr[0]*alphaDrSurf_r[3]+0.3535533905932737*fr[1]*alphaDrSurf_r[2]; 

 } else { 

  Ghat_l[0] = 0.7905694150420947*alphaDrSurf_l[3]*fl[11]+0.7905694150420947*alphaDrSurf_l[2]*fl[10]+0.7905694150420947*alphaDrSurf_l[1]*fl[9]+0.7905694150420947*alphaDrSurf_l[0]*fl[8]+0.6123724356957944*alphaDrSurf_l[3]*fl[7]+0.6123724356957944*alphaDrSurf_l[2]*fl[6]+0.3535533905932737*alphaDrSurf_l[3]*fl[5]+0.6123724356957944*alphaDrSurf_l[1]*fl[4]+0.3535533905932737*alphaDrSurf_l[2]*fl[3]+0.6123724356957944*alphaDrSurf_l[0]*fl[2]+0.3535533905932737*alphaDrSurf_l[1]*fl[1]+0.3535533905932737*alphaDrSurf_l[0]*fl[0]; 
  Ghat_l[1] = 0.7905694150420947*alphaDrSurf_l[2]*fl[11]+0.7905694150420947*alphaDrSurf_l[3]*fl[10]+0.7905694150420947*alphaDrSurf_l[0]*fl[9]+0.7905694150420947*alphaDrSurf_l[1]*fl[8]+0.6123724356957944*alphaDrSurf_l[2]*fl[7]+0.6123724356957944*alphaDrSurf_l[3]*fl[6]+0.3535533905932737*alphaDrSurf_l[2]*fl[5]+0.6123724356957944*alphaDrSurf_l[0]*fl[4]+0.3535533905932737*alphaDrSurf_l[3]*fl[3]+0.6123724356957944*alphaDrSurf_l[1]*fl[2]+0.3535533905932737*alphaDrSurf_l[0]*fl[1]+0.3535533905932737*fl[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = 0.7905694150420947*alphaDrSurf_l[1]*fl[11]+0.7905694150420947*alphaDrSurf_l[0]*fl[10]+0.7905694150420947*alphaDrSurf_l[3]*fl[9]+0.7905694150420947*alphaDrSurf_l[2]*fl[8]+0.6123724356957944*alphaDrSurf_l[1]*fl[7]+0.6123724356957944*alphaDrSurf_l[0]*fl[6]+0.3535533905932737*alphaDrSurf_l[1]*fl[5]+0.6123724356957944*alphaDrSurf_l[3]*fl[4]+0.3535533905932737*alphaDrSurf_l[0]*fl[3]+0.3535533905932737*fl[1]*alphaDrSurf_l[3]+0.6123724356957944*alphaDrSurf_l[2]*fl[2]+0.3535533905932737*fl[0]*alphaDrSurf_l[2]; 
  Ghat_l[3] = 0.7905694150420947*alphaDrSurf_l[0]*fl[11]+0.7905694150420947*alphaDrSurf_l[1]*fl[10]+0.7905694150420947*alphaDrSurf_l[2]*fl[9]+0.7905694150420947*alphaDrSurf_l[3]*fl[8]+0.6123724356957944*alphaDrSurf_l[0]*fl[7]+0.6123724356957944*alphaDrSurf_l[1]*fl[6]+0.3535533905932737*alphaDrSurf_l[0]*fl[5]+0.6123724356957944*alphaDrSurf_l[2]*fl[4]+0.3535533905932737*alphaDrSurf_l[1]*fl[3]+0.6123724356957944*fl[2]*alphaDrSurf_l[3]+0.3535533905932737*fl[0]*alphaDrSurf_l[3]+0.3535533905932737*fl[1]*alphaDrSurf_l[2]; 

  Ghat_r[0] = 0.7905694150420947*alphaDrSurf_r[3]*fc[11]+0.7905694150420947*alphaDrSurf_r[2]*fc[10]+0.7905694150420947*alphaDrSurf_r[1]*fc[9]+0.7905694150420947*alphaDrSurf_r[0]*fc[8]+0.6123724356957944*alphaDrSurf_r[3]*fc[7]+0.6123724356957944*alphaDrSurf_r[2]*fc[6]+0.3535533905932737*alphaDrSurf_r[3]*fc[5]+0.6123724356957944*alphaDrSurf_r[1]*fc[4]+0.3535533905932737*alphaDrSurf_r[2]*fc[3]+0.6123724356957944*alphaDrSurf_r[0]*fc[2]+0.3535533905932737*alphaDrSurf_r[1]*fc[1]+0.3535533905932737*alphaDrSurf_r[0]*fc[0]; 
  Ghat_r[1] = 0.7905694150420947*alphaDrSurf_r[2]*fc[11]+0.7905694150420947*alphaDrSurf_r[3]*fc[10]+0.7905694150420947*alphaDrSurf_r[0]*fc[9]+0.7905694150420947*alphaDrSurf_r[1]*fc[8]+0.6123724356957944*alphaDrSurf_r[2]*fc[7]+0.6123724356957944*alphaDrSurf_r[3]*fc[6]+0.3535533905932737*alphaDrSurf_r[2]*fc[5]+0.6123724356957944*alphaDrSurf_r[0]*fc[4]+0.3535533905932737*alphaDrSurf_r[3]*fc[3]+0.6123724356957944*alphaDrSurf_r[1]*fc[2]+0.3535533905932737*alphaDrSurf_r[0]*fc[1]+0.3535533905932737*fc[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = 0.7905694150420947*alphaDrSurf_r[1]*fc[11]+0.7905694150420947*alphaDrSurf_r[0]*fc[10]+0.7905694150420947*alphaDrSurf_r[3]*fc[9]+0.7905694150420947*alphaDrSurf_r[2]*fc[8]+0.6123724356957944*alphaDrSurf_r[1]*fc[7]+0.6123724356957944*alphaDrSurf_r[0]*fc[6]+0.3535533905932737*alphaDrSurf_r[1]*fc[5]+0.6123724356957944*alphaDrSurf_r[3]*fc[4]+0.3535533905932737*alphaDrSurf_r[0]*fc[3]+0.3535533905932737*fc[1]*alphaDrSurf_r[3]+0.6123724356957944*alphaDrSurf_r[2]*fc[2]+0.3535533905932737*fc[0]*alphaDrSurf_r[2]; 
  Ghat_r[3] = 0.7905694150420947*alphaDrSurf_r[0]*fc[11]+0.7905694150420947*alphaDrSurf_r[1]*fc[10]+0.7905694150420947*alphaDrSurf_r[2]*fc[9]+0.7905694150420947*alphaDrSurf_r[3]*fc[8]+0.6123724356957944*alphaDrSurf_r[0]*fc[7]+0.6123724356957944*alphaDrSurf_r[1]*fc[6]+0.3535533905932737*alphaDrSurf_r[0]*fc[5]+0.6123724356957944*alphaDrSurf_r[2]*fc[4]+0.3535533905932737*alphaDrSurf_r[1]*fc[3]+0.6123724356957944*fc[2]*alphaDrSurf_r[3]+0.3535533905932737*fc[0]*alphaDrSurf_r[3]+0.3535533905932737*fc[1]*alphaDrSurf_r[2]; 

 } 
  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*rdv2; 
  out[2] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[3] += (0.7071067811865475*Ghat_r[2]-0.7071067811865475*Ghat_l[2])*rdv2; 
  out[4] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[5] += (0.7071067811865475*Ghat_r[3]-0.7071067811865475*Ghat_l[3])*rdv2; 
  out[6] += 1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdv2; 
  out[7] += 1.224744871391589*(Ghat_r[3]+Ghat_l[3])*rdv2; 
  out[8] += (1.58113883008419*Ghat_r[0]-1.58113883008419*Ghat_l[0])*rdv2; 
  out[9] += (1.58113883008419*Ghat_r[1]-1.58113883008419*Ghat_l[1])*rdv2; 
  out[10] += (1.58113883008419*Ghat_r[2]-1.58113883008419*Ghat_l[2])*rdv2; 
  out[11] += (1.58113883008419*Ghat_r[3]-1.58113883008419*Ghat_l[3])*rdv2; 

  return 0.;

} 
