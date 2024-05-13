#include <gkyl_vlasov_poisson_kernels.h> 
GKYL_CU_DH double vlasov_poisson_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
    const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat[3]; 

  if (edge == -1) { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fskin[4]+1.224744871391589*fskin[1]+0.7071067811865475*fskin[0])*wv+(0.4564354645876384*fskin[6]+0.3535533905932737*fskin[3]+0.2041241452319315*fskin[2])*dv; 
  Ghat[1] = (1.58113883008419*fskin[6]+1.224744871391589*fskin[3]+0.7071067811865475*fskin[2])*wv+(0.3162277660168379*fskin[7]+0.1825741858350554*fskin[5]+0.4564354645876384*fskin[4]+0.3535533905932737*fskin[1]+0.2041241452319315*fskin[0])*dv; 
  Ghat[2] = (1.224744871391589*fskin[7]+0.7071067811865475*fskin[5])*wv+(0.408248290463863*fskin[6]+0.3162277660168379*fskin[3]+0.1825741858350554*fskin[2])*dv; 

  } else { 

  Ghat[0] = 1.58113883008419*fedge[4]*wv-1.224744871391589*fedge[1]*wv+0.7071067811865475*fedge[0]*wv+0.4564354645876383*fedge[6]*dv-0.3535533905932737*fedge[3]*dv+0.2041241452319315*fedge[2]*dv; 
  Ghat[1] = 1.58113883008419*fedge[6]*wv-1.224744871391589*fedge[3]*wv+0.7071067811865475*fedge[2]*wv-0.3162277660168379*fedge[7]*dv+0.1825741858350554*fedge[5]*dv+0.4564354645876384*fedge[4]*dv-0.3535533905932737*fedge[1]*dv+0.2041241452319315*fedge[0]*dv; 
  Ghat[2] = (-1.224744871391589*fedge[7]*wv)+0.7071067811865475*fedge[5]*wv+0.408248290463863*fedge[6]*dv-0.3162277660168379*fedge[3]*dv+0.1825741858350554*fedge[2]*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += -1.58113883008419*Ghat[0]*dx10; 
  out[5] += -0.7071067811865475*Ghat[2]*dx10; 
  out[6] += -1.58113883008419*Ghat[1]*dx10; 
  out[7] += -1.224744871391589*Ghat[2]*dx10; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.58113883008419*fedge[4]+1.224744871391589*fedge[1]+0.7071067811865475*fedge[0])*wv+(0.4564354645876384*fedge[6]+0.3535533905932737*fedge[3]+0.2041241452319315*fedge[2])*dv; 
  Ghat[1] = (1.58113883008419*fedge[6]+1.224744871391589*fedge[3]+0.7071067811865475*fedge[2])*wv+(0.3162277660168379*fedge[7]+0.1825741858350554*fedge[5]+0.4564354645876384*fedge[4]+0.3535533905932737*fedge[1]+0.2041241452319315*fedge[0])*dv; 
  Ghat[2] = (1.224744871391589*fedge[7]+0.7071067811865475*fedge[5])*wv+(0.408248290463863*fedge[6]+0.3162277660168379*fedge[3]+0.1825741858350554*fedge[2])*dv; 

  } else { 

  Ghat[0] = 1.58113883008419*fskin[4]*wv-1.224744871391589*fskin[1]*wv+0.7071067811865475*fskin[0]*wv+0.4564354645876383*fskin[6]*dv-0.3535533905932737*fskin[3]*dv+0.2041241452319315*fskin[2]*dv; 
  Ghat[1] = 1.58113883008419*fskin[6]*wv-1.224744871391589*fskin[3]*wv+0.7071067811865475*fskin[2]*wv-0.3162277660168379*fskin[7]*dv+0.1825741858350554*fskin[5]*dv+0.4564354645876384*fskin[4]*dv-0.3535533905932737*fskin[1]*dv+0.2041241452319315*fskin[0]*dv; 
  Ghat[2] = (-1.224744871391589*fskin[7]*wv)+0.7071067811865475*fskin[5]*wv+0.408248290463863*fskin[6]*dv-0.3162277660168379*fskin[3]*dv+0.1825741858350554*fskin[2]*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 
  out[4] += 1.58113883008419*Ghat[0]*dx10; 
  out[5] += 0.7071067811865475*Ghat[2]*dx10; 
  out[6] += 1.58113883008419*Ghat[1]*dx10; 
  out[7] += -1.224744871391589*Ghat[2]*dx10; 

  } 
  return 0.;

} 
