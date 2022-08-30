#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH void lbo_vlasov_pkpm_drag_boundary_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, const double *nu, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nu:         collisionality. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 

  double alphaDrSurf[3] = {0.0}; 
  double Ghat[3] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nu[0]*wvpar+0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar+0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar+0.5*nu[2]*dvpar; 

  Ghat[0] = 1.118033988749895*alphaDrSurf[1]*fSkin[7]+0.8660254037844387*alphaDrSurf[2]*fSkin[6]+1.118033988749895*alphaDrSurf[0]*fSkin[5]+0.5*alphaDrSurf[2]*fSkin[4]+0.8660254037844386*alphaDrSurf[1]*fSkin[3]+0.8660254037844386*alphaDrSurf[0]*fSkin[2]+0.5*alphaDrSurf[1]*fSkin[1]+0.5*alphaDrSurf[0]*fSkin[0]; 
  Ghat[1] = 1.0*alphaDrSurf[2]*fSkin[7]+1.118033988749895*alphaDrSurf[0]*fSkin[7]+0.7745966692414834*alphaDrSurf[1]*fSkin[6]+1.118033988749895*alphaDrSurf[1]*fSkin[5]+0.4472135954999579*alphaDrSurf[1]*fSkin[4]+0.7745966692414833*alphaDrSurf[2]*fSkin[3]+0.8660254037844386*alphaDrSurf[0]*fSkin[3]+0.8660254037844386*alphaDrSurf[1]*fSkin[2]+0.4472135954999579*fSkin[1]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*fSkin[1]+0.5*fSkin[0]*alphaDrSurf[1]; 
  Ghat[2] = 1.0*alphaDrSurf[1]*fSkin[7]+0.5532833351724881*alphaDrSurf[2]*fSkin[6]+0.8660254037844387*alphaDrSurf[0]*fSkin[6]+1.118033988749895*alphaDrSurf[2]*fSkin[5]+0.31943828249997*alphaDrSurf[2]*fSkin[4]+0.5*alphaDrSurf[0]*fSkin[4]+0.7745966692414833*alphaDrSurf[1]*fSkin[3]+0.8660254037844386*alphaDrSurf[2]*fSkin[2]+0.5*fSkin[0]*alphaDrSurf[2]+0.4472135954999579*alphaDrSurf[1]*fSkin[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += 1.224744871391589*Ghat[0]*dv1par; 
  out[3] += 1.224744871391589*Ghat[1]*dv1par; 
  out[4] += 0.7071067811865475*Ghat[2]*dv1par; 
  out[5] += 1.58113883008419*Ghat[0]*dv1par; 
  out[6] += 1.224744871391589*Ghat[2]*dv1par; 
  out[7] += 1.58113883008419*Ghat[1]*dv1par; 

  } else { 

  alphaDrSurf[0] = nu[0]*wvpar-0.5*nu[0]*dvpar; 
  alphaDrSurf[1] = nu[1]*wvpar-0.5*nu[1]*dvpar; 
  alphaDrSurf[2] = nu[2]*wvpar-0.5*nu[2]*dvpar; 

  Ghat[0] = 1.118033988749895*alphaDrSurf[1]*fSkin[7]-0.8660254037844387*alphaDrSurf[2]*fSkin[6]+1.118033988749895*alphaDrSurf[0]*fSkin[5]+0.5*alphaDrSurf[2]*fSkin[4]-0.8660254037844386*alphaDrSurf[1]*fSkin[3]-0.8660254037844386*alphaDrSurf[0]*fSkin[2]+0.5*alphaDrSurf[1]*fSkin[1]+0.5*alphaDrSurf[0]*fSkin[0]; 
  Ghat[1] = 1.0*alphaDrSurf[2]*fSkin[7]+1.118033988749895*alphaDrSurf[0]*fSkin[7]-0.7745966692414834*alphaDrSurf[1]*fSkin[6]+1.118033988749895*alphaDrSurf[1]*fSkin[5]+0.4472135954999579*alphaDrSurf[1]*fSkin[4]-0.7745966692414833*alphaDrSurf[2]*fSkin[3]-0.8660254037844386*alphaDrSurf[0]*fSkin[3]-0.8660254037844386*alphaDrSurf[1]*fSkin[2]+0.4472135954999579*fSkin[1]*alphaDrSurf[2]+0.5*alphaDrSurf[0]*fSkin[1]+0.5*fSkin[0]*alphaDrSurf[1]; 
  Ghat[2] = 1.0*alphaDrSurf[1]*fSkin[7]-0.5532833351724881*alphaDrSurf[2]*fSkin[6]-0.8660254037844387*alphaDrSurf[0]*fSkin[6]+1.118033988749895*alphaDrSurf[2]*fSkin[5]+0.31943828249997*alphaDrSurf[2]*fSkin[4]+0.5*alphaDrSurf[0]*fSkin[4]-0.7745966692414833*alphaDrSurf[1]*fSkin[3]-0.8660254037844386*alphaDrSurf[2]*fSkin[2]+0.5*fSkin[0]*alphaDrSurf[2]+0.4472135954999579*alphaDrSurf[1]*fSkin[1]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += 1.224744871391589*Ghat[0]*dv1par; 
  out[3] += 1.224744871391589*Ghat[1]*dv1par; 
  out[4] += -0.7071067811865475*Ghat[2]*dv1par; 
  out[5] += -1.58113883008419*Ghat[0]*dv1par; 
  out[6] += 1.224744871391589*Ghat[2]*dv1par; 
  out[7] += -1.58113883008419*Ghat[1]*dv1par; 

  } 
} 
