#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nI, const double *nuField, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:     cell-center coordinates. 
  // dxv[2]:   cell spacing. 
  // nI:        ion density. 
  // nuField:   2/pi*v*nu(v) field dg representation (v'(v||,mu) in notes
  // fl/fc/fr:  distribution function in cells 
  // out:       incremented distribution function in cell 

  double rdv2 = 2.0/dxv[1]; 

  double Ghat_r[2] = {0.0}; 
  double Ghat_l[2] = {0.0}; 
  double alphaDrSurf_l[2] = {0.0}; 
  alphaDrSurf_l[0] = 1.1180339887498951*nI[1]*nuField[5]+1.118033988749895*nI[0]*nuField[4]-0.8660254037844386*nI[1]*nuField[3]-0.8660254037844386*nI[0]*nuField[2]+0.5*nI[1]*nuField[1]+0.5*nI[0]*nuField[0]; 
  alphaDrSurf_l[1] = 1.1180339887498951*nI[0]*nuField[5]+1.118033988749895*nI[1]*nuField[4]-0.8660254037844386*nI[0]*nuField[3]-0.8660254037844386*nI[1]*nuField[2]+0.5*nI[0]*nuField[1]+0.5*nuField[0]*nI[1]; 

  double alphaDrSurf_r[2] = {0.0}; 
  alphaDrSurf_r[0] = 1.1180339887498951*nI[1]*nuField[5]+1.118033988749895*nI[0]*nuField[4]+0.8660254037844386*nI[1]*nuField[3]+0.8660254037844386*nI[0]*nuField[2]+0.5*nI[1]*nuField[1]+0.5*nI[0]*nuField[0]; 
  alphaDrSurf_r[1] = 1.1180339887498951*nI[0]*nuField[5]+1.118033988749895*nI[1]*nuField[4]+0.8660254037844386*nI[0]*nuField[3]+0.8660254037844386*nI[1]*nuField[2]+0.5*nI[0]*nuField[1]+0.5*nuField[0]*nI[1]; 

  Ghat_l[0] = 0.16666666666666666*(6.7082039324993685*alphaDrSurf_l[1]*fc[5]+6.708203932499369*alphaDrSurf_l[0]*fc[4]-5.196152422706631*(alphaDrSurf_l[1]*fc[3]+alphaDrSurf_l[0]*fc[2])+3.0*(alphaDrSurf_l[1]*fc[1]+alphaDrSurf_l[0]*fc[0])); 
  Ghat_l[1] = 0.16666666666666666*(6.7082039324993685*alphaDrSurf_l[0]*fc[5]+6.708203932499369*alphaDrSurf_l[1]*fc[4]-5.196152422706631*(alphaDrSurf_l[0]*fc[3]+alphaDrSurf_l[1]*fc[2])+3.0*(alphaDrSurf_l[0]*fc[1]+fc[0]*alphaDrSurf_l[1])); 

  Ghat_r[0] = 0.16666666666666666*(6.7082039324993685*alphaDrSurf_r[1]*fr[5]+6.708203932499369*alphaDrSurf_r[0]*fr[4]-5.196152422706631*(alphaDrSurf_r[1]*fr[3]+alphaDrSurf_r[0]*fr[2])+3.0*(alphaDrSurf_r[1]*fr[1]+alphaDrSurf_r[0]*fr[0])); 
  Ghat_r[1] = 0.16666666666666666*(6.7082039324993685*alphaDrSurf_r[0]*fr[5]+6.708203932499369*alphaDrSurf_r[1]*fr[4]-5.196152422706631*(alphaDrSurf_r[0]*fr[3]+alphaDrSurf_r[1]*fr[2])+3.0*(alphaDrSurf_r[0]*fr[1]+fr[0]*alphaDrSurf_r[1])); 

  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*rdv2; 
  out[2] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[3] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[4] += (1.5811388300841895*Ghat_r[0]-1.5811388300841895*Ghat_l[0])*rdv2; 
  out[5] += (1.5811388300841898*Ghat_r[1]-1.5811388300841898*Ghat_l[1])*rdv2; 

  return 0.;

} 
