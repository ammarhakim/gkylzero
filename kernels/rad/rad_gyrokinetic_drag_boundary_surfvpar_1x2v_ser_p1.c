#include <gkyl_rad_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_1x2v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double *nI, const double *nuField, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // nI:        ion density. 
  // nuField:       2/pi*v*nu(v) field dg representation (v'(v||,mu) in notes
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 

  double rdv2 = 2.0/dxv[1]; 
  double alphaDrSurf[4] = {0.0}; 
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  if (w[1]>0) { 

  Ghat[0] = 0.7905694150420947*(alphaDrSurf[3]*fEdge[11]+alphaDrSurf[2]*fEdge[10]+alphaDrSurf[1]*fEdge[9]+alphaDrSurf[0]*fEdge[8])-0.6123724356957944*(alphaDrSurf[3]*fEdge[7]+alphaDrSurf[2]*fEdge[6])+0.3535533905932737*alphaDrSurf[3]*fEdge[5]-0.6123724356957944*alphaDrSurf[1]*fEdge[4]+0.3535533905932737*alphaDrSurf[2]*fEdge[3]-0.6123724356957944*alphaDrSurf[0]*fEdge[2]+0.3535533905932737*(alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0]); 
  Ghat[1] = 0.7905694150420947*(alphaDrSurf[2]*fEdge[11]+alphaDrSurf[3]*fEdge[10]+alphaDrSurf[0]*fEdge[9]+alphaDrSurf[1]*fEdge[8])-0.6123724356957944*(alphaDrSurf[2]*fEdge[7]+alphaDrSurf[3]*fEdge[6])+0.3535533905932737*alphaDrSurf[2]*fEdge[5]-0.6123724356957944*alphaDrSurf[0]*fEdge[4]+0.3535533905932737*alphaDrSurf[3]*fEdge[3]-0.6123724356957944*alphaDrSurf[1]*fEdge[2]+0.3535533905932737*(alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.7905694150420947*(alphaDrSurf[1]*fEdge[11]+alphaDrSurf[0]*fEdge[10]+alphaDrSurf[3]*fEdge[9]+alphaDrSurf[2]*fEdge[8])-0.6123724356957944*(alphaDrSurf[1]*fEdge[7]+alphaDrSurf[0]*fEdge[6])+0.3535533905932737*alphaDrSurf[1]*fEdge[5]-0.6123724356957944*alphaDrSurf[3]*fEdge[4]+0.3535533905932737*(alphaDrSurf[0]*fEdge[3]+fEdge[1]*alphaDrSurf[3])+alphaDrSurf[2]*(0.3535533905932737*fEdge[0]-0.6123724356957944*fEdge[2]); 
  Ghat[3] = 0.7905694150420947*(alphaDrSurf[0]*fEdge[11]+alphaDrSurf[1]*fEdge[10]+alphaDrSurf[2]*fEdge[9]+alphaDrSurf[3]*fEdge[8])-0.6123724356957944*(alphaDrSurf[0]*fEdge[7]+alphaDrSurf[1]*fEdge[6])+0.3535533905932737*alphaDrSurf[0]*fEdge[5]-0.6123724356957944*alphaDrSurf[2]*fEdge[4]+0.3535533905932737*alphaDrSurf[1]*fEdge[3]-0.6123724356957944*fEdge[2]*alphaDrSurf[3]+0.3535533905932737*(fEdge[0]*alphaDrSurf[3]+fEdge[1]*alphaDrSurf[2]); 

  } else { 

  Ghat[0] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf[3]*fSkin[11]+alphaDrSurf[2]*fSkin[10]+alphaDrSurf[1]*fSkin[9]+alphaDrSurf[0]*fSkin[8])+7.348469228349534*(alphaDrSurf[3]*fSkin[7]+alphaDrSurf[2]*fSkin[6])+4.242640687119286*alphaDrSurf[3]*fSkin[5]+7.348469228349534*alphaDrSurf[1]*fSkin[4]+4.242640687119286*alphaDrSurf[2]*fSkin[3]+7.348469228349534*alphaDrSurf[0]*fSkin[2]+4.242640687119286*(alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])); 
  Ghat[1] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf[2]*fSkin[11]+alphaDrSurf[3]*fSkin[10]+alphaDrSurf[0]*fSkin[9]+alphaDrSurf[1]*fSkin[8])+7.348469228349534*(alphaDrSurf[2]*fSkin[7]+alphaDrSurf[3]*fSkin[6])+4.242640687119286*alphaDrSurf[2]*fSkin[5]+7.348469228349534*alphaDrSurf[0]*fSkin[4]+4.242640687119286*alphaDrSurf[3]*fSkin[3]+7.348469228349534*alphaDrSurf[1]*fSkin[2]+4.242640687119286*(alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])); 
  Ghat[2] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf[1]*fSkin[11]+alphaDrSurf[0]*fSkin[10]+alphaDrSurf[3]*fSkin[9]+alphaDrSurf[2]*fSkin[8])+7.348469228349534*(alphaDrSurf[1]*fSkin[7]+alphaDrSurf[0]*fSkin[6])+4.242640687119286*alphaDrSurf[1]*fSkin[5]+7.348469228349534*alphaDrSurf[3]*fSkin[4]+4.242640687119286*(alphaDrSurf[0]*fSkin[3]+fSkin[1]*alphaDrSurf[3])+alphaDrSurf[2]*(7.348469228349534*fSkin[2]+4.242640687119286*fSkin[0])); 
  Ghat[3] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf[0]*fSkin[11]+alphaDrSurf[1]*fSkin[10]+alphaDrSurf[2]*fSkin[9]+alphaDrSurf[3]*fSkin[8])+7.348469228349534*(alphaDrSurf[0]*fSkin[7]+alphaDrSurf[1]*fSkin[6])+4.242640687119286*alphaDrSurf[0]*fSkin[5]+7.348469228349534*alphaDrSurf[2]*fSkin[4]+4.242640687119286*alphaDrSurf[1]*fSkin[3]+(7.348469228349534*fSkin[2]+4.242640687119286*fSkin[0])*alphaDrSurf[3]+4.242640687119286*fSkin[1]*alphaDrSurf[2]); 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -1.224744871391589*Ghat[0]*rdv2; 
  out[3] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[4] += -1.224744871391589*Ghat[1]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[6] += -1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -1.224744871391589*Ghat[3]*rdv2; 
  out[8] += -1.58113883008419*Ghat[0]*rdv2; 
  out[9] += -1.58113883008419*Ghat[1]*rdv2; 
  out[10] += -1.58113883008419*Ghat[2]*rdv2; 
  out[11] += -1.58113883008419*Ghat[3]*rdv2; 

  } else { 

  if (w[1]>0) { 

  Ghat[0] = 0.7905694150420947*(alphaDrSurf[3]*fSkin[11]+alphaDrSurf[2]*fSkin[10]+alphaDrSurf[1]*fSkin[9]+alphaDrSurf[0]*fSkin[8])-0.6123724356957944*(alphaDrSurf[3]*fSkin[7]+alphaDrSurf[2]*fSkin[6])+0.3535533905932737*alphaDrSurf[3]*fSkin[5]-0.6123724356957944*alphaDrSurf[1]*fSkin[4]+0.3535533905932737*alphaDrSurf[2]*fSkin[3]-0.6123724356957944*alphaDrSurf[0]*fSkin[2]+0.3535533905932737*(alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0]); 
  Ghat[1] = 0.7905694150420947*(alphaDrSurf[2]*fSkin[11]+alphaDrSurf[3]*fSkin[10]+alphaDrSurf[0]*fSkin[9]+alphaDrSurf[1]*fSkin[8])-0.6123724356957944*(alphaDrSurf[2]*fSkin[7]+alphaDrSurf[3]*fSkin[6])+0.3535533905932737*alphaDrSurf[2]*fSkin[5]-0.6123724356957944*alphaDrSurf[0]*fSkin[4]+0.3535533905932737*alphaDrSurf[3]*fSkin[3]-0.6123724356957944*alphaDrSurf[1]*fSkin[2]+0.3535533905932737*(alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.7905694150420947*(alphaDrSurf[1]*fSkin[11]+alphaDrSurf[0]*fSkin[10]+alphaDrSurf[3]*fSkin[9]+alphaDrSurf[2]*fSkin[8])-0.6123724356957944*(alphaDrSurf[1]*fSkin[7]+alphaDrSurf[0]*fSkin[6])+0.3535533905932737*alphaDrSurf[1]*fSkin[5]-0.6123724356957944*alphaDrSurf[3]*fSkin[4]+0.3535533905932737*(alphaDrSurf[0]*fSkin[3]+fSkin[1]*alphaDrSurf[3])+alphaDrSurf[2]*(0.3535533905932737*fSkin[0]-0.6123724356957944*fSkin[2]); 
  Ghat[3] = 0.7905694150420947*(alphaDrSurf[0]*fSkin[11]+alphaDrSurf[1]*fSkin[10]+alphaDrSurf[2]*fSkin[9]+alphaDrSurf[3]*fSkin[8])-0.6123724356957944*(alphaDrSurf[0]*fSkin[7]+alphaDrSurf[1]*fSkin[6])+0.3535533905932737*alphaDrSurf[0]*fSkin[5]-0.6123724356957944*alphaDrSurf[2]*fSkin[4]+0.3535533905932737*alphaDrSurf[1]*fSkin[3]-0.6123724356957944*fSkin[2]*alphaDrSurf[3]+0.3535533905932737*(fSkin[0]*alphaDrSurf[3]+fSkin[1]*alphaDrSurf[2]); 

  } else { 

  Ghat[0] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf[3]*fEdge[11]+alphaDrSurf[2]*fEdge[10]+alphaDrSurf[1]*fEdge[9]+alphaDrSurf[0]*fEdge[8])+7.348469228349534*(alphaDrSurf[3]*fEdge[7]+alphaDrSurf[2]*fEdge[6])+4.242640687119286*alphaDrSurf[3]*fEdge[5]+7.348469228349534*alphaDrSurf[1]*fEdge[4]+4.242640687119286*alphaDrSurf[2]*fEdge[3]+7.348469228349534*alphaDrSurf[0]*fEdge[2]+4.242640687119286*(alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])); 
  Ghat[1] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf[2]*fEdge[11]+alphaDrSurf[3]*fEdge[10]+alphaDrSurf[0]*fEdge[9]+alphaDrSurf[1]*fEdge[8])+7.348469228349534*(alphaDrSurf[2]*fEdge[7]+alphaDrSurf[3]*fEdge[6])+4.242640687119286*alphaDrSurf[2]*fEdge[5]+7.348469228349534*alphaDrSurf[0]*fEdge[4]+4.242640687119286*alphaDrSurf[3]*fEdge[3]+7.348469228349534*alphaDrSurf[1]*fEdge[2]+4.242640687119286*(alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])); 
  Ghat[2] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf[1]*fEdge[11]+alphaDrSurf[0]*fEdge[10]+alphaDrSurf[3]*fEdge[9]+alphaDrSurf[2]*fEdge[8])+7.348469228349534*(alphaDrSurf[1]*fEdge[7]+alphaDrSurf[0]*fEdge[6])+4.242640687119286*alphaDrSurf[1]*fEdge[5]+7.348469228349534*alphaDrSurf[3]*fEdge[4]+4.242640687119286*(alphaDrSurf[0]*fEdge[3]+fEdge[1]*alphaDrSurf[3])+alphaDrSurf[2]*(7.348469228349534*fEdge[2]+4.242640687119286*fEdge[0])); 
  Ghat[3] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf[0]*fEdge[11]+alphaDrSurf[1]*fEdge[10]+alphaDrSurf[2]*fEdge[9]+alphaDrSurf[3]*fEdge[8])+7.348469228349534*(alphaDrSurf[0]*fEdge[7]+alphaDrSurf[1]*fEdge[6])+4.242640687119286*alphaDrSurf[0]*fEdge[5]+7.348469228349534*alphaDrSurf[2]*fEdge[4]+4.242640687119286*alphaDrSurf[1]*fEdge[3]+(7.348469228349534*fEdge[2]+4.242640687119286*fEdge[0])*alphaDrSurf[3]+4.242640687119286*fEdge[1]*alphaDrSurf[2]); 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[4] += -1.224744871391589*Ghat[1]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[6] += -1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -1.224744871391589*Ghat[3]*rdv2; 
  out[8] += 1.58113883008419*Ghat[0]*rdv2; 
  out[9] += 1.58113883008419*Ghat[1]*rdv2; 
  out[10] += 1.58113883008419*Ghat[2]*rdv2; 
  out[11] += 1.58113883008419*Ghat[3]*rdv2; 

  } 
  return 0.;

} 
