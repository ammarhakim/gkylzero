#include <gkyl_vlasov_poisson_kernels.h> 
GKYL_CU_DH double vlasov_poisson_nonuniformv_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *vcoord, int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
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

  Ghat[0] = vcoord[2]*(0.8660254037844386*fSkin[5]+0.5*fSkin[4])+vcoord[1]*(0.8660254037844386*fSkin[3]+0.5*fSkin[2])+vcoord[0]*(0.8660254037844386*fSkin[1]+0.5*fSkin[0]); 
  Ghat[1] = vcoord[1]*(0.7745966692414833*fSkin[5]+0.4472135954999579*fSkin[4])+(0.7745966692414833*vcoord[2]+0.8660254037844386*vcoord[0])*fSkin[3]+fSkin[2]*(0.4472135954999579*vcoord[2]+0.5*vcoord[0])+(0.8660254037844386*fSkin[1]+0.5*fSkin[0])*vcoord[1]; 
  Ghat[2] = (0.5532833351724881*vcoord[2]+0.8660254037844386*vcoord[0])*fSkin[5]+(0.31943828249997*vcoord[2]+0.5*vcoord[0])*fSkin[4]+0.7745966692414833*vcoord[1]*fSkin[3]+(0.8660254037844386*fSkin[1]+0.5*fSkin[0])*vcoord[2]+0.4472135954999579*vcoord[1]*fSkin[2]; 

  } else { 

  Ghat[0] = -0.1*(vcoord[2]*(8.660254037844387*fEdge[5]-5.0*fEdge[4])+vcoord[1]*(8.660254037844386*fEdge[3]-5.0*fEdge[2])+vcoord[0]*(8.660254037844386*fEdge[1]-5.0*fEdge[0])); 
  Ghat[1] = -0.1*(vcoord[1]*(7.745966692414834*fEdge[5]-4.47213595499958*fEdge[4])+(7.745966692414834*vcoord[2]+8.660254037844386*vcoord[0])*fEdge[3]+fEdge[2]*((-4.47213595499958*vcoord[2])-5.0*vcoord[0])+(8.660254037844386*fEdge[1]-5.0*fEdge[0])*vcoord[1]); 
  Ghat[2] = -0.01428571428571429*((38.72983346207417*vcoord[2]+60.62177826491071*vcoord[0])*fEdge[5]+((-22.3606797749979*vcoord[2])-35.0*vcoord[0])*fEdge[4]+54.22176684690384*vcoord[1]*fEdge[3]+(60.6217782649107*fEdge[1]-35.0*fEdge[0])*vcoord[2]-31.30495168499706*vcoord[1]*fEdge[2]); 

  } 

  out[0] += -0.7071067811865475*Ghat[0]*rdx2; 
  out[1] += -1.224744871391589*Ghat[0]*rdx2; 
  out[2] += -0.7071067811865475*Ghat[1]*rdx2; 
  out[3] += -1.224744871391589*Ghat[1]*rdx2; 
  out[4] += -0.7071067811865475*Ghat[2]*rdx2; 
  out[5] += -1.224744871391589*Ghat[2]*rdx2; 

  } else { 

  if (0.7071067811865475*vcoord[0]-0.7905694150420947*vcoord[2]>0) { 

  Ghat[0] = vcoord[2]*(0.8660254037844386*fEdge[5]+0.5*fEdge[4])+vcoord[1]*(0.8660254037844386*fEdge[3]+0.5*fEdge[2])+vcoord[0]*(0.8660254037844386*fEdge[1]+0.5*fEdge[0]); 
  Ghat[1] = vcoord[1]*(0.7745966692414833*fEdge[5]+0.4472135954999579*fEdge[4])+(0.7745966692414833*vcoord[2]+0.8660254037844386*vcoord[0])*fEdge[3]+fEdge[2]*(0.4472135954999579*vcoord[2]+0.5*vcoord[0])+(0.8660254037844386*fEdge[1]+0.5*fEdge[0])*vcoord[1]; 
  Ghat[2] = (0.5532833351724881*vcoord[2]+0.8660254037844386*vcoord[0])*fEdge[5]+(0.31943828249997*vcoord[2]+0.5*vcoord[0])*fEdge[4]+0.7745966692414833*vcoord[1]*fEdge[3]+(0.8660254037844386*fEdge[1]+0.5*fEdge[0])*vcoord[2]+0.4472135954999579*vcoord[1]*fEdge[2]; 

  } else { 

  Ghat[0] = -0.1*(vcoord[2]*(8.660254037844387*fSkin[5]-5.0*fSkin[4])+vcoord[1]*(8.660254037844386*fSkin[3]-5.0*fSkin[2])+vcoord[0]*(8.660254037844386*fSkin[1]-5.0*fSkin[0])); 
  Ghat[1] = -0.1*(vcoord[1]*(7.745966692414834*fSkin[5]-4.47213595499958*fSkin[4])+(7.745966692414834*vcoord[2]+8.660254037844386*vcoord[0])*fSkin[3]+fSkin[2]*((-4.47213595499958*vcoord[2])-5.0*vcoord[0])+(8.660254037844386*fSkin[1]-5.0*fSkin[0])*vcoord[1]); 
  Ghat[2] = -0.01428571428571429*((38.72983346207417*vcoord[2]+60.62177826491071*vcoord[0])*fSkin[5]+((-22.3606797749979*vcoord[2])-35.0*vcoord[0])*fSkin[4]+54.22176684690384*vcoord[1]*fSkin[3]+(60.6217782649107*fSkin[1]-35.0*fSkin[0])*vcoord[2]-31.30495168499706*vcoord[1]*fSkin[2]); 

  } 

  out[0] += 0.7071067811865475*Ghat[0]*rdx2; 
  out[1] += -1.224744871391589*Ghat[0]*rdx2; 
  out[2] += 0.7071067811865475*Ghat[1]*rdx2; 
  out[3] += -1.224744871391589*Ghat[1]*rdx2; 
  out[4] += 0.7071067811865475*Ghat[2]*rdx2; 
  out[5] += -1.224744871391589*Ghat[2]*rdx2; 

  } 
  return 0.;

} 
