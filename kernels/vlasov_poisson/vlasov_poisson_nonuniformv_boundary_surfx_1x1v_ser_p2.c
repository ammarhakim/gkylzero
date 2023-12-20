#include <gkyl_vlasov_poisson_kernels.h> 
GKYL_CU_DH double vlasov_poisson_nonuniformv_boundary_surfx_1x1v_ser_p2(const double *w, const double *dxv, const double *vcoord, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // vcoord:      Discrete (DG) velocity coordinate.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double rdx2 = 2/dxv[0]; 
  double Ghat[3]; 

  if (edge == -1) { 

  if (0.7071067811865475*vcoord[0]-0.7905694150420947*vcoord[2]>0) { 

  Ghat[0] = 0.8660254037844386*vcoord[2]*fSkin[7]+1.118033988749895*vcoord[1]*fSkin[6]+0.5*vcoord[2]*fSkin[5]+1.118033988749895*vcoord[0]*fSkin[4]+vcoord[1]*(0.8660254037844386*fSkin[3]+0.5*fSkin[2])+vcoord[0]*(0.8660254037844386*fSkin[1]+0.5*fSkin[0]); 
  Ghat[1] = 0.7745966692414833*vcoord[1]*fSkin[7]+(vcoord[2]+1.118033988749895*vcoord[0])*fSkin[6]+vcoord[1]*(0.4472135954999579*fSkin[5]+1.118033988749895*fSkin[4])+(0.7745966692414833*vcoord[2]+0.8660254037844386*vcoord[0])*fSkin[3]+fSkin[2]*(0.4472135954999579*vcoord[2]+0.5*vcoord[0])+(0.8660254037844386*fSkin[1]+0.5*fSkin[0])*vcoord[1]; 
  Ghat[2] = (0.5532833351724881*vcoord[2]+0.8660254037844386*vcoord[0])*fSkin[7]+vcoord[1]*fSkin[6]+(0.31943828249997*vcoord[2]+0.5*vcoord[0])*fSkin[5]+1.118033988749895*vcoord[2]*fSkin[4]+0.7745966692414833*vcoord[1]*fSkin[3]+(0.8660254037844386*fSkin[1]+0.5*fSkin[0])*vcoord[2]+0.4472135954999579*vcoord[1]*fSkin[2]; 

  } else { 

  Ghat[0] = (-0.8660254037844387*vcoord[2]*fEdge[7])+1.118033988749895*vcoord[1]*fEdge[6]+0.5*vcoord[2]*fEdge[5]+1.118033988749895*vcoord[0]*fEdge[4]-0.8660254037844386*vcoord[1]*fEdge[3]+0.5*vcoord[1]*fEdge[2]-0.8660254037844386*vcoord[0]*fEdge[1]+0.5*fEdge[0]*vcoord[0]; 
  Ghat[1] = (-0.7745966692414834*vcoord[1]*fEdge[7])+1.0*vcoord[2]*fEdge[6]+1.118033988749895*vcoord[0]*fEdge[6]+0.4472135954999579*vcoord[1]*fEdge[5]+1.118033988749895*vcoord[1]*fEdge[4]-0.7745966692414833*vcoord[2]*fEdge[3]-0.8660254037844386*vcoord[0]*fEdge[3]+0.4472135954999579*fEdge[2]*vcoord[2]+0.5*vcoord[0]*fEdge[2]-0.8660254037844386*fEdge[1]*vcoord[1]+0.5*fEdge[0]*vcoord[1]; 
  Ghat[2] = (-0.5532833351724881*vcoord[2]*fEdge[7])-0.8660254037844387*vcoord[0]*fEdge[7]+1.0*vcoord[1]*fEdge[6]+0.31943828249997*vcoord[2]*fEdge[5]+0.5*vcoord[0]*fEdge[5]+1.118033988749895*vcoord[2]*fEdge[4]-0.7745966692414833*vcoord[1]*fEdge[3]-0.8660254037844386*fEdge[1]*vcoord[2]+0.5*fEdge[0]*vcoord[2]+0.4472135954999579*vcoord[1]*fEdge[2]; 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*rdx2; 
  out[1] += -1.224744871391589*Ghat[0]*rdx2; 
  out[2] += -0.7071067811865475*Ghat[1]*rdx2; 
  out[3] += -1.224744871391589*Ghat[1]*rdx2; 
  out[4] += -1.58113883008419*Ghat[0]*rdx2; 
  out[5] += -0.7071067811865475*Ghat[2]*rdx2; 
  out[6] += -1.58113883008419*Ghat[1]*rdx2; 
  out[7] += -1.224744871391589*Ghat[2]*rdx2; 

  } else { 

  if (0.7071067811865475*vcoord[0]-0.7905694150420947*vcoord[2]>0) { 

  Ghat[0] = 0.8660254037844386*vcoord[2]*fEdge[7]+1.118033988749895*vcoord[1]*fEdge[6]+0.5*vcoord[2]*fEdge[5]+1.118033988749895*vcoord[0]*fEdge[4]+vcoord[1]*(0.8660254037844386*fEdge[3]+0.5*fEdge[2])+vcoord[0]*(0.8660254037844386*fEdge[1]+0.5*fEdge[0]); 
  Ghat[1] = 0.7745966692414833*vcoord[1]*fEdge[7]+(vcoord[2]+1.118033988749895*vcoord[0])*fEdge[6]+vcoord[1]*(0.4472135954999579*fEdge[5]+1.118033988749895*fEdge[4])+(0.7745966692414833*vcoord[2]+0.8660254037844386*vcoord[0])*fEdge[3]+fEdge[2]*(0.4472135954999579*vcoord[2]+0.5*vcoord[0])+(0.8660254037844386*fEdge[1]+0.5*fEdge[0])*vcoord[1]; 
  Ghat[2] = (0.5532833351724881*vcoord[2]+0.8660254037844386*vcoord[0])*fEdge[7]+vcoord[1]*fEdge[6]+(0.31943828249997*vcoord[2]+0.5*vcoord[0])*fEdge[5]+1.118033988749895*vcoord[2]*fEdge[4]+0.7745966692414833*vcoord[1]*fEdge[3]+(0.8660254037844386*fEdge[1]+0.5*fEdge[0])*vcoord[2]+0.4472135954999579*vcoord[1]*fEdge[2]; 

  } else { 

  Ghat[0] = (-0.8660254037844387*vcoord[2]*fSkin[7])+1.118033988749895*vcoord[1]*fSkin[6]+0.5*vcoord[2]*fSkin[5]+1.118033988749895*vcoord[0]*fSkin[4]-0.8660254037844386*vcoord[1]*fSkin[3]+0.5*vcoord[1]*fSkin[2]-0.8660254037844386*vcoord[0]*fSkin[1]+0.5*fSkin[0]*vcoord[0]; 
  Ghat[1] = (-0.7745966692414834*vcoord[1]*fSkin[7])+1.0*vcoord[2]*fSkin[6]+1.118033988749895*vcoord[0]*fSkin[6]+0.4472135954999579*vcoord[1]*fSkin[5]+1.118033988749895*vcoord[1]*fSkin[4]-0.7745966692414833*vcoord[2]*fSkin[3]-0.8660254037844386*vcoord[0]*fSkin[3]+0.4472135954999579*fSkin[2]*vcoord[2]+0.5*vcoord[0]*fSkin[2]-0.8660254037844386*fSkin[1]*vcoord[1]+0.5*fSkin[0]*vcoord[1]; 
  Ghat[2] = (-0.5532833351724881*vcoord[2]*fSkin[7])-0.8660254037844387*vcoord[0]*fSkin[7]+1.0*vcoord[1]*fSkin[6]+0.31943828249997*vcoord[2]*fSkin[5]+0.5*vcoord[0]*fSkin[5]+1.118033988749895*vcoord[2]*fSkin[4]-0.7745966692414833*vcoord[1]*fSkin[3]-0.8660254037844386*fSkin[1]*vcoord[2]+0.5*fSkin[0]*vcoord[2]+0.4472135954999579*vcoord[1]*fSkin[2]; 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*rdx2; 
  out[1] += -1.224744871391589*Ghat[0]*rdx2; 
  out[2] += 0.7071067811865475*Ghat[1]*rdx2; 
  out[3] += -1.224744871391589*Ghat[1]*rdx2; 
  out[4] += 1.58113883008419*Ghat[0]*rdx2; 
  out[5] += 0.7071067811865475*Ghat[2]*rdx2; 
  out[6] += 1.58113883008419*Ghat[1]*rdx2; 
  out[7] += -1.224744871391589*Ghat[2]*rdx2; 

  } 
  return 0.;

} 
