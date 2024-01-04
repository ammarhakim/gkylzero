#include <gkyl_rad_gyrokinetic_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double *nI, const double *nuField, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:     cell-center coordinates. 
  // dxv[2]:   cell spacing. 
  // nI:        ion density. 
  // nuField:       2/pi*v*nu(v) field dg representation (v'(v||,mu) in notes
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 

  double rdv2 = 2.0/dxv[1]; 
  double alphaDrSurf[3] = {0.0}; 
  double Ghat[3] = {0.0}; 

  if (edge == -1) { 

  if (w[1]>0) { 

  Ghat[0] = 1.118033988749895*alphaDrSurf[1]*fEdge[7]-0.8660254037844386*alphaDrSurf[2]*fEdge[6]+1.118033988749895*alphaDrSurf[0]*fEdge[5]+0.5*alphaDrSurf[2]*fEdge[4]-0.8660254037844386*(alphaDrSurf[1]*fEdge[3]+alphaDrSurf[0]*fEdge[2])+0.5*(alphaDrSurf[1]*fEdge[1]+alphaDrSurf[0]*fEdge[0]); 
  Ghat[1] = (alphaDrSurf[2]+1.118033988749895*alphaDrSurf[0])*fEdge[7]+alphaDrSurf[1]*((-0.7745966692414833*fEdge[6])+1.118033988749895*fEdge[5]+0.4472135954999579*fEdge[4])-0.7745966692414833*alphaDrSurf[2]*fEdge[3]-0.8660254037844386*(alphaDrSurf[0]*fEdge[3]+alphaDrSurf[1]*fEdge[2])+0.4472135954999579*fEdge[1]*alphaDrSurf[2]+0.5*(alphaDrSurf[0]*fEdge[1]+fEdge[0]*alphaDrSurf[1]); 
  Ghat[2] = alphaDrSurf[1]*fEdge[7]+((-0.5532833351724881*alphaDrSurf[2])-0.8660254037844386*alphaDrSurf[0])*fEdge[6]+1.118033988749895*alphaDrSurf[2]*fEdge[5]+(0.31943828249997*alphaDrSurf[2]+0.5*alphaDrSurf[0])*fEdge[4]-0.7745966692414833*alphaDrSurf[1]*fEdge[3]+alphaDrSurf[2]*(0.5*fEdge[0]-0.8660254037844386*fEdge[2])+0.4472135954999579*alphaDrSurf[1]*fEdge[1]; 

  } else { 

  Ghat[0] = 1.118033988749895*alphaDrSurf[1]*fSkin[7]+0.8660254037844387*alphaDrSurf[2]*fSkin[6]+1.118033988749895*alphaDrSurf[0]*fSkin[5]+0.5*alphaDrSurf[2]*fSkin[4]+0.8660254037844386*alphaDrSurf[1]*fSkin[3]+0.8660254037844386*alphaDrSurf[0]*fSkin[2]+0.5*alphaDrSurf[1]*fSkin[1]+0.5*alphaDrSurf[0]*fSkin[0]; 
  Ghat[1] = 1.0*alphaDrSurf[2]*fSkin[7]+1.118033988749895*alphaDrSurf[0]*fSkin[7]+0.7745966692414834*alphaDrSurf[1]*fSkin[6]+1.118033988749895*alphaDrSurf[1]*fSkin[5]+0.4472135954999579*alphaDrSurf[1]*fSkin[4]+0.7745966692414833*alphaDrSurf[2]*fSkin[3]+0.8660254037844386*alphaDrSurf[0]*fSkin[3]+0.8660254037844386*alphaDrSurf[1]*fSkin[2]+0.4472135954999579*fSkin[1]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*fSkin[1]+0.5*fSkin[0]*alphaDrSurf[1]; 
  Ghat[2] = 1.0*alphaDrSurf[1]*fSkin[7]+0.5532833351724881*alphaDrSurf[2]*fSkin[6]+0.8660254037844387*alphaDrSurf[0]*fSkin[6]+1.118033988749895*alphaDrSurf[2]*fSkin[5]+0.31943828249997*alphaDrSurf[2]*fSkin[4]+0.5*alphaDrSurf[0]*fSkin[4]+0.7745966692414833*alphaDrSurf[1]*fSkin[3]+0.8660254037844386*alphaDrSurf[2]*fSkin[2]+0.5*fSkin[0]*alphaDrSurf[2]+0.4472135954999579*alphaDrSurf[1]*fSkin[1]; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -1.224744871391589*Ghat[0]*rdv2; 
  out[3] += -1.224744871391589*Ghat[1]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[5] += -1.58113883008419*Ghat[0]*rdv2; 
  out[6] += -1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -1.58113883008419*Ghat[1]*rdv2; 

  } else { 

  if (w[1]>0) { 

  Ghat[0] = 1.118033988749895*alphaDrSurf[1]*fSkin[7]-0.8660254037844386*alphaDrSurf[2]*fSkin[6]+1.118033988749895*alphaDrSurf[0]*fSkin[5]+0.5*alphaDrSurf[2]*fSkin[4]-0.8660254037844386*(alphaDrSurf[1]*fSkin[3]+alphaDrSurf[0]*fSkin[2])+0.5*(alphaDrSurf[1]*fSkin[1]+alphaDrSurf[0]*fSkin[0]); 
  Ghat[1] = (alphaDrSurf[2]+1.118033988749895*alphaDrSurf[0])*fSkin[7]+alphaDrSurf[1]*((-0.7745966692414833*fSkin[6])+1.118033988749895*fSkin[5]+0.4472135954999579*fSkin[4])-0.7745966692414833*alphaDrSurf[2]*fSkin[3]-0.8660254037844386*(alphaDrSurf[0]*fSkin[3]+alphaDrSurf[1]*fSkin[2])+0.4472135954999579*fSkin[1]*alphaDrSurf[2]+0.5*(alphaDrSurf[0]*fSkin[1]+fSkin[0]*alphaDrSurf[1]); 
  Ghat[2] = alphaDrSurf[1]*fSkin[7]+((-0.5532833351724881*alphaDrSurf[2])-0.8660254037844386*alphaDrSurf[0])*fSkin[6]+1.118033988749895*alphaDrSurf[2]*fSkin[5]+(0.31943828249997*alphaDrSurf[2]+0.5*alphaDrSurf[0])*fSkin[4]-0.7745966692414833*alphaDrSurf[1]*fSkin[3]+alphaDrSurf[2]*(0.5*fSkin[0]-0.8660254037844386*fSkin[2])+0.4472135954999579*alphaDrSurf[1]*fSkin[1]; 

  } else { 

  Ghat[0] = 1.118033988749895*alphaDrSurf[1]*fEdge[7]+0.8660254037844387*alphaDrSurf[2]*fEdge[6]+1.118033988749895*alphaDrSurf[0]*fEdge[5]+0.5*alphaDrSurf[2]*fEdge[4]+0.8660254037844386*alphaDrSurf[1]*fEdge[3]+0.8660254037844386*alphaDrSurf[0]*fEdge[2]+0.5*alphaDrSurf[1]*fEdge[1]+0.5*alphaDrSurf[0]*fEdge[0]; 
  Ghat[1] = 1.0*alphaDrSurf[2]*fEdge[7]+1.118033988749895*alphaDrSurf[0]*fEdge[7]+0.7745966692414834*alphaDrSurf[1]*fEdge[6]+1.118033988749895*alphaDrSurf[1]*fEdge[5]+0.4472135954999579*alphaDrSurf[1]*fEdge[4]+0.7745966692414833*alphaDrSurf[2]*fEdge[3]+0.8660254037844386*alphaDrSurf[0]*fEdge[3]+0.8660254037844386*alphaDrSurf[1]*fEdge[2]+0.4472135954999579*fEdge[1]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*fEdge[1]+0.5*fEdge[0]*alphaDrSurf[1]; 
  Ghat[2] = 1.0*alphaDrSurf[1]*fEdge[7]+0.5532833351724881*alphaDrSurf[2]*fEdge[6]+0.8660254037844387*alphaDrSurf[0]*fEdge[6]+1.118033988749895*alphaDrSurf[2]*fEdge[5]+0.31943828249997*alphaDrSurf[2]*fEdge[4]+0.5*alphaDrSurf[0]*fEdge[4]+0.7745966692414833*alphaDrSurf[1]*fEdge[3]+0.8660254037844386*alphaDrSurf[2]*fEdge[2]+0.5*fEdge[0]*alphaDrSurf[2]+0.4472135954999579*alphaDrSurf[1]*fEdge[1]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -1.224744871391589*Ghat[0]*rdv2; 
  out[3] += -1.224744871391589*Ghat[1]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[5] += 1.58113883008419*Ghat[0]*rdv2; 
  out[6] += -1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.58113883008419*Ghat[1]*rdv2; 

  } 
  return 0.;

} 
