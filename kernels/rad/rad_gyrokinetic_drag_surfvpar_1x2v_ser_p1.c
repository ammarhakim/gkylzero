#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double *nI, const double *nuField, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // nI:        ion density. 
  // nuField:   2/pi*v*nu(v) field dg representation (v'(v||,mu) in notes
  // fl/fc/fr:  distribution function in cells 
  // out:       incremented distribution function in cell 

  double rdv2 = 2.0/dxv[1]; 

  double Ghat_r[4] = {0.0}; 
  double Ghat_l[4] = {0.0}; 
  double alphaDrSurf_l[4] = {0.0}; 
  alphaDrSurf_l[0] = 1.118033988749895*nI[1]*nuField[9]+1.118033988749895*nI[0]*nuField[8]-0.8660254037844386*nI[1]*nuField[4]-0.8660254037844386*nI[0]*nuField[2]+0.5*nI[1]*nuField[1]+0.5*nI[0]*nuField[0]; 
  alphaDrSurf_l[1] = 1.118033988749895*nI[0]*nuField[9]+1.118033988749895*nI[1]*nuField[8]-0.8660254037844386*nI[0]*nuField[4]-0.8660254037844386*nI[1]*nuField[2]+0.5*nI[0]*nuField[1]+0.5*nuField[0]*nI[1]; 
  alphaDrSurf_l[2] = 1.118033988749895*nI[1]*nuField[11]+1.118033988749895*nI[0]*nuField[10]-0.8660254037844386*nI[1]*nuField[7]-0.8660254037844386*nI[0]*nuField[6]+0.5*nI[1]*nuField[5]+0.5*nI[0]*nuField[3]; 
  alphaDrSurf_l[3] = 1.118033988749895*nI[0]*nuField[11]+1.118033988749895*nI[1]*nuField[10]-0.8660254037844386*nI[0]*nuField[7]-0.8660254037844386*nI[1]*nuField[6]+0.5*nI[0]*nuField[5]+0.5*nI[1]*nuField[3]; 

  double alphaDrSurf_r[4] = {0.0}; 
  alphaDrSurf_r[0] = 1.118033988749895*nI[1]*nuField[9]+1.118033988749895*nI[0]*nuField[8]+0.8660254037844386*nI[1]*nuField[4]+0.8660254037844386*nI[0]*nuField[2]+0.5*nI[1]*nuField[1]+0.5*nI[0]*nuField[0]; 
  alphaDrSurf_r[1] = 1.118033988749895*nI[0]*nuField[9]+1.118033988749895*nI[1]*nuField[8]+0.8660254037844386*nI[0]*nuField[4]+0.8660254037844386*nI[1]*nuField[2]+0.5*nI[0]*nuField[1]+0.5*nuField[0]*nI[1]; 
  alphaDrSurf_r[2] = 1.118033988749895*nI[1]*nuField[11]+1.118033988749895*nI[0]*nuField[10]+0.8660254037844386*nI[1]*nuField[7]+0.8660254037844386*nI[0]*nuField[6]+0.5*nI[1]*nuField[5]+0.5*nI[0]*nuField[3]; 
  alphaDrSurf_r[3] = 1.118033988749895*nI[0]*nuField[11]+1.118033988749895*nI[1]*nuField[10]+0.8660254037844386*nI[0]*nuField[7]+0.8660254037844386*nI[1]*nuField[6]+0.5*nI[0]*nuField[5]+0.5*nI[1]*nuField[3]; 

  Ghat_l[0] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_l[3]*fc[11]+alphaDrSurf_l[2]*fc[10]+alphaDrSurf_l[1]*fc[9]+alphaDrSurf_l[0]*fc[8])-7.348469228349534*(alphaDrSurf_l[3]*fc[7]+alphaDrSurf_l[2]*fc[6])+4.242640687119286*alphaDrSurf_l[3]*fc[5]-7.348469228349534*alphaDrSurf_l[1]*fc[4]+4.242640687119286*alphaDrSurf_l[2]*fc[3]-7.348469228349534*alphaDrSurf_l[0]*fc[2]+4.242640687119286*(alphaDrSurf_l[1]*fc[1]+alphaDrSurf_l[0]*fc[0])); 
  Ghat_l[1] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_l[2]*fc[11]+alphaDrSurf_l[3]*fc[10]+alphaDrSurf_l[0]*fc[9]+alphaDrSurf_l[1]*fc[8])-7.348469228349534*(alphaDrSurf_l[2]*fc[7]+alphaDrSurf_l[3]*fc[6])+4.242640687119286*alphaDrSurf_l[2]*fc[5]-7.348469228349534*alphaDrSurf_l[0]*fc[4]+4.242640687119286*alphaDrSurf_l[3]*fc[3]-7.348469228349534*alphaDrSurf_l[1]*fc[2]+4.242640687119286*(alphaDrSurf_l[0]*fc[1]+fc[0]*alphaDrSurf_l[1])); 
  Ghat_l[2] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_l[1]*fc[11]+alphaDrSurf_l[0]*fc[10]+alphaDrSurf_l[3]*fc[9]+alphaDrSurf_l[2]*fc[8])-7.348469228349534*(alphaDrSurf_l[1]*fc[7]+alphaDrSurf_l[0]*fc[6])+4.242640687119286*alphaDrSurf_l[1]*fc[5]-7.348469228349534*alphaDrSurf_l[3]*fc[4]+4.242640687119286*(alphaDrSurf_l[0]*fc[3]+fc[1]*alphaDrSurf_l[3])+alphaDrSurf_l[2]*(4.242640687119286*fc[0]-7.348469228349534*fc[2])); 
  Ghat_l[3] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_l[0]*fc[11]+alphaDrSurf_l[1]*fc[10]+alphaDrSurf_l[2]*fc[9]+alphaDrSurf_l[3]*fc[8])-7.348469228349534*(alphaDrSurf_l[0]*fc[7]+alphaDrSurf_l[1]*fc[6])+4.242640687119286*alphaDrSurf_l[0]*fc[5]-7.348469228349534*alphaDrSurf_l[2]*fc[4]+4.242640687119286*alphaDrSurf_l[1]*fc[3]+(4.242640687119286*fc[0]-7.348469228349534*fc[2])*alphaDrSurf_l[3]+4.242640687119286*fc[1]*alphaDrSurf_l[2]); 

  Ghat_r[0] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_r[3]*fr[11]+alphaDrSurf_r[2]*fr[10]+alphaDrSurf_r[1]*fr[9]+alphaDrSurf_r[0]*fr[8])-7.348469228349534*(alphaDrSurf_r[3]*fr[7]+alphaDrSurf_r[2]*fr[6])+4.242640687119286*alphaDrSurf_r[3]*fr[5]-7.348469228349534*alphaDrSurf_r[1]*fr[4]+4.242640687119286*alphaDrSurf_r[2]*fr[3]-7.348469228349534*alphaDrSurf_r[0]*fr[2]+4.242640687119286*(alphaDrSurf_r[1]*fr[1]+alphaDrSurf_r[0]*fr[0])); 
  Ghat_r[1] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_r[2]*fr[11]+alphaDrSurf_r[3]*fr[10]+alphaDrSurf_r[0]*fr[9]+alphaDrSurf_r[1]*fr[8])-7.348469228349534*(alphaDrSurf_r[2]*fr[7]+alphaDrSurf_r[3]*fr[6])+4.242640687119286*alphaDrSurf_r[2]*fr[5]-7.348469228349534*alphaDrSurf_r[0]*fr[4]+4.242640687119286*alphaDrSurf_r[3]*fr[3]-7.348469228349534*alphaDrSurf_r[1]*fr[2]+4.242640687119286*(alphaDrSurf_r[0]*fr[1]+fr[0]*alphaDrSurf_r[1])); 
  Ghat_r[2] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_r[1]*fr[11]+alphaDrSurf_r[0]*fr[10]+alphaDrSurf_r[3]*fr[9]+alphaDrSurf_r[2]*fr[8])-7.348469228349534*(alphaDrSurf_r[1]*fr[7]+alphaDrSurf_r[0]*fr[6])+4.242640687119286*alphaDrSurf_r[1]*fr[5]-7.348469228349534*alphaDrSurf_r[3]*fr[4]+4.242640687119286*(alphaDrSurf_r[0]*fr[3]+fr[1]*alphaDrSurf_r[3])+alphaDrSurf_r[2]*(4.242640687119286*fr[0]-7.348469228349534*fr[2])); 
  Ghat_r[3] = 0.08333333333333333*(9.48683298050514*(alphaDrSurf_r[0]*fr[11]+alphaDrSurf_r[1]*fr[10]+alphaDrSurf_r[2]*fr[9]+alphaDrSurf_r[3]*fr[8])-7.348469228349534*(alphaDrSurf_r[0]*fr[7]+alphaDrSurf_r[1]*fr[6])+4.242640687119286*alphaDrSurf_r[0]*fr[5]-7.348469228349534*alphaDrSurf_r[2]*fr[4]+4.242640687119286*alphaDrSurf_r[1]*fr[3]+(4.242640687119286*fr[0]-7.348469228349534*fr[2])*alphaDrSurf_r[3]+4.242640687119286*fr[1]*alphaDrSurf_r[2]); 

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
