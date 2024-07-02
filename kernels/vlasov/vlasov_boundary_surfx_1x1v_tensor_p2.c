#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_tensor_p2(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // alpha_geo:   Fields used only for general geometry.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[3]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fSkin[4]+1.224744871391589*fSkin[1]+0.7071067811865475*fSkin[0])*wv+(0.4564354645876384*fSkin[6]+0.3535533905932737*fSkin[3]+0.2041241452319315*fSkin[2])*dv; 
  Ghat[1] = (1.58113883008419*fSkin[6]+1.224744871391589*fSkin[3]+0.7071067811865475*fSkin[2])*wv+(0.408248290463863*fSkin[8]+0.3162277660168379*fSkin[7]+0.1825741858350554*fSkin[5]+0.4564354645876384*fSkin[4]+0.3535533905932737*fSkin[1]+0.2041241452319315*fSkin[0])*dv; 
  Ghat[2] = (1.58113883008419*fSkin[8]+1.224744871391589*fSkin[7]+0.7071067811865475*fSkin[5])*wv+(0.408248290463863*fSkin[6]+0.3162277660168379*fSkin[3]+0.1825741858350554*fSkin[2])*dv; 

  } else { 

  Ghat[0] = 1.58113883008419*fEdge[4]*wv-1.224744871391589*fEdge[1]*wv+0.7071067811865475*fEdge[0]*wv+0.4564354645876383*fEdge[6]*dv-0.3535533905932737*fEdge[3]*dv+0.2041241452319315*fEdge[2]*dv; 
  Ghat[1] = 1.58113883008419*fEdge[6]*wv-1.224744871391589*fEdge[3]*wv+0.7071067811865475*fEdge[2]*wv+0.408248290463863*fEdge[8]*dv-0.3162277660168379*fEdge[7]*dv+0.1825741858350554*fEdge[5]*dv+0.4564354645876384*fEdge[4]*dv-0.3535533905932737*fEdge[1]*dv+0.2041241452319315*fEdge[0]*dv; 
  Ghat[2] = 1.58113883008419*fEdge[8]*wv-1.224744871391589*fEdge[7]*wv+0.7071067811865475*fEdge[5]*wv+0.408248290463863*fEdge[6]*dv-0.3162277660168379*fEdge[3]*dv+0.1825741858350554*fEdge[2]*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += -1.58113883008419*Ghat[0]*dx10; 
  out[5] += -0.7071067811865475*Ghat[2]*dx10; 
  out[6] += -1.58113883008419*Ghat[1]*dx10; 
  out[7] += -1.224744871391589*Ghat[2]*dx10; 
  out[8] += -1.58113883008419*Ghat[2]*dx10; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fEdge[4]+1.224744871391589*fEdge[1]+0.7071067811865475*fEdge[0])*wv+(0.4564354645876384*fEdge[6]+0.3535533905932737*fEdge[3]+0.2041241452319315*fEdge[2])*dv; 
  Ghat[1] = (1.58113883008419*fEdge[6]+1.224744871391589*fEdge[3]+0.7071067811865475*fEdge[2])*wv+(0.408248290463863*fEdge[8]+0.3162277660168379*fEdge[7]+0.1825741858350554*fEdge[5]+0.4564354645876384*fEdge[4]+0.3535533905932737*fEdge[1]+0.2041241452319315*fEdge[0])*dv; 
  Ghat[2] = (1.58113883008419*fEdge[8]+1.224744871391589*fEdge[7]+0.7071067811865475*fEdge[5])*wv+(0.408248290463863*fEdge[6]+0.3162277660168379*fEdge[3]+0.1825741858350554*fEdge[2])*dv; 

  } else { 

  Ghat[0] = 1.58113883008419*fSkin[4]*wv-1.224744871391589*fSkin[1]*wv+0.7071067811865475*fSkin[0]*wv+0.4564354645876383*fSkin[6]*dv-0.3535533905932737*fSkin[3]*dv+0.2041241452319315*fSkin[2]*dv; 
  Ghat[1] = 1.58113883008419*fSkin[6]*wv-1.224744871391589*fSkin[3]*wv+0.7071067811865475*fSkin[2]*wv+0.408248290463863*fSkin[8]*dv-0.3162277660168379*fSkin[7]*dv+0.1825741858350554*fSkin[5]*dv+0.4564354645876384*fSkin[4]*dv-0.3535533905932737*fSkin[1]*dv+0.2041241452319315*fSkin[0]*dv; 
  Ghat[2] = 1.58113883008419*fSkin[8]*wv-1.224744871391589*fSkin[7]*wv+0.7071067811865475*fSkin[5]*wv+0.408248290463863*fSkin[6]*dv-0.3162277660168379*fSkin[3]*dv+0.1825741858350554*fSkin[2]*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += 1.58113883008419*Ghat[0]*dx10; 
  out[5] += 0.7071067811865475*Ghat[2]*dx10; 
  out[6] += 1.58113883008419*Ghat[1]*dx10; 
  out[7] += -1.224744871391589*Ghat[2]*dx10; 
  out[8] += 1.58113883008419*Ghat[2]*dx10; 

  } 
  return 0.;

} 
