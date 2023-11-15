#include <gkyl_rad_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nI, const double *nuField, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:     cell-center coordinates. 
  // dxv[2]:   cell spacing. 
  // nI:        ion density. 
  // nuField:       2/pi*v*nu(v) field dg representation (v'(v||,mu) in notes
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 

  double rdv2 = 2.0/dxv[1]; 
  double alphaDrSurf[2] = {0.0}; 
  double Ghat[2] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.16666666666666666*(6.7082039324993685*nI[1]*nuField[5]+6.708203932499369*nI[0]*nuField[4]+5.196152422706631*(nI[1]*nuField[3]+nI[0]*nuField[2])+3.0*(nI[1]*nuField[1]+nI[0]*nuField[0])); 
  alphaDrSurf[1] = 0.16666666666666666*(6.7082039324993685*nI[0]*nuField[5]+6.708203932499369*nI[1]*nuField[4]+5.196152422706631*(nI[0]*nuField[3]+nI[1]*nuField[2])+3.0*(nI[0]*nuField[1]+nuField[0]*nI[1])); 

  Ghat[0] = 0.16666666666666666*(6.7082039324993685*alphaDrSurf[1]*fEdge[5]+6.708203932499369*alphaDrSurf[0]*fEdge[4]-5.196152422706631*(alphaDrSurf[1]*fEdge[3]+alphaDrSurf[0]*fEdge[2])+3.0*(alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0])); 
  Ghat[1] = 0.16666666666666666*(6.7082039324993685*alphaDrSurf[0]*fEdge[5]+6.708203932499369*alphaDrSurf[1]*fEdge[4]-5.196152422706631*(alphaDrSurf[0]*fEdge[3]+alphaDrSurf[1]*fEdge[2])+3.0*(alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1])); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += 1.5811388300841895*Ghat[0]*rdv2; 
  out[5] += 1.5811388300841898*Ghat[1]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.16666666666666666*(6.7082039324993685*nI[1]*nuField[5]+6.708203932499369*nI[0]*nuField[4]-5.196152422706631*(nI[1]*nuField[3]+nI[0]*nuField[2])+3.0*(nI[1]*nuField[1]+nI[0]*nuField[0])); 
  alphaDrSurf[1] = 0.16666666666666666*(6.7082039324993685*nI[0]*nuField[5]+6.708203932499369*nI[1]*nuField[4]-5.196152422706631*(nI[0]*nuField[3]+nI[1]*nuField[2])+3.0*(nI[0]*nuField[1]+nuField[0]*nI[1])); 

  Ghat[0] = 0.16666666666666666*(6.7082039324993685*alphaDrSurf[1]*fSkin[5]+6.708203932499369*alphaDrSurf[0]*fSkin[4]-5.196152422706631*(alphaDrSurf[1]*fSkin[3]+alphaDrSurf[0]*fSkin[2])+3.0*(alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0])); 
  Ghat[1] = 0.16666666666666666*(6.7082039324993685*alphaDrSurf[0]*fSkin[5]+6.708203932499369*alphaDrSurf[1]*fSkin[4]-5.196152422706631*(alphaDrSurf[0]*fSkin[3]+alphaDrSurf[1]*fSkin[2])+3.0*(alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1])); 

  out[0] += -(0.7071067811865475*Ghat[0]*rdv2); 
  out[1] += -(0.7071067811865475*Ghat[1]*rdv2); 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += -(1.5811388300841895*Ghat[0]*rdv2); 
  out[5] += -(1.5811388300841898*Ghat[1]*rdv2); 

  } 

  return 0.;

} 
