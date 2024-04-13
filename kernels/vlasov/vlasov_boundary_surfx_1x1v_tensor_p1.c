#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_tensor_p1(const double *w, const double *dxv, const double *alpha_geo, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
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

  Ghat[0] = (1.224744871391589*fSkin[1]+0.7071067811865475*fSkin[0])*wv+(0.3535533905932737*fSkin[3]+0.2041241452319315*fSkin[2])*dv; 
  Ghat[1] = (1.224744871391589*fSkin[3]+0.7071067811865475*fSkin[2])*wv+(0.3535533905932737*fSkin[1]+0.2041241452319315*fSkin[0])*dv; 
  Ghat[2] = (0.3162277660168379*fSkin[3]+0.1825741858350554*fSkin[2])*dv; 

  } else { 

  Ghat[0] = -0.1178511301977579*((10.39230484541326*fEdge[1]-6.0*fEdge[0])*wv+(3.0*fEdge[3]-1.732050807568877*fEdge[2])*dv); 
  Ghat[1] = -0.1178511301977579*((10.39230484541326*fEdge[3]-6.0*fEdge[2])*wv+(3.0*fEdge[1]-1.732050807568877*fEdge[0])*dv); 
  Ghat[2] = -0.04714045207910316*(6.708203932499369*fEdge[3]-3.872983346207417*fEdge[2])*dv; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += -0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 

  } else { 

  if (wv>0) { 

  Ghat[0] = (1.224744871391589*fEdge[1]+0.7071067811865475*fEdge[0])*wv+(0.3535533905932737*fEdge[3]+0.2041241452319315*fEdge[2])*dv; 
  Ghat[1] = (1.224744871391589*fEdge[3]+0.7071067811865475*fEdge[2])*wv+(0.3535533905932737*fEdge[1]+0.2041241452319315*fEdge[0])*dv; 
  Ghat[2] = (0.3162277660168379*fEdge[3]+0.1825741858350554*fEdge[2])*dv; 

  } else { 

  Ghat[0] = -0.1178511301977579*((10.39230484541326*fSkin[1]-6.0*fSkin[0])*wv+(3.0*fSkin[3]-1.732050807568877*fSkin[2])*dv); 
  Ghat[1] = -0.1178511301977579*((10.39230484541326*fSkin[3]-6.0*fSkin[2])*wv+(3.0*fSkin[1]-1.732050807568877*fSkin[0])*dv); 
  Ghat[2] = -0.04714045207910316*(6.708203932499369*fSkin[3]-3.872983346207417*fSkin[2])*dv; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*dx10; 
  out[1] += -1.224744871391589*Ghat[0]*dx10; 
  out[2] += 0.7071067811865475*Ghat[1]*dx10; 
  out[3] += -1.224744871391589*Ghat[1]*dx10; 

  } 
  return 0.;

} 
