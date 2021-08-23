#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double *E0 = &qmem[0]; 

  double Ghat[3]; 
  double alpha[3]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 
  alpha[2] = E0[2]; 

  double fUpwindQuad[3];
  double fUpwind[3];

  if (edge == -1) { 

  fUpwindQuad[0] = (-1.5*fSkin[7])+copysignf(1.0,0.6324555320336759*alpha[2]-0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0])*((-1.5*fSkin[7])+1.5*fEdge[7]+0.7745966692414833*(fSkin[6]+fEdge[6])+1.118033988749895*fSkin[5]-1.118033988749895*fEdge[5]+0.4472135954999579*fSkin[4]-0.4472135954999579*fEdge[4]-1.161895003862225*(fSkin[3]+fEdge[3])+0.8660254037844386*(fSkin[2]+fEdge[2])-0.6708203932499369*fSkin[1]+0.6708203932499369*fEdge[1]+0.5*fSkin[0]-0.5*fEdge[0])-1.5*fEdge[7]+0.7745966692414833*fSkin[6]-0.7745966692414833*fEdge[6]+1.118033988749895*(fSkin[5]+fEdge[5])+0.4472135954999579*(fSkin[4]+fEdge[4])-1.161895003862225*fSkin[3]+1.161895003862225*fEdge[3]+0.8660254037844386*fSkin[2]-0.8660254037844386*fEdge[2]-0.6708203932499369*(fSkin[1]+fEdge[1])+0.5*(fSkin[0]+fEdge[0]); 
  fUpwindQuad[1] = copysignf(1.0,0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2])*((-0.9682458365518543*(fSkin[6]+fEdge[6]))+1.118033988749895*fSkin[5]-1.118033988749895*fEdge[5]-0.5590169943749475*fSkin[4]+0.5590169943749475*fEdge[4]+0.8660254037844386*(fSkin[2]+fEdge[2])+0.5*fSkin[0]-0.5*fEdge[0])-0.9682458365518543*fSkin[6]+0.9682458365518543*fEdge[6]+1.118033988749895*(fSkin[5]+fEdge[5])-0.5590169943749475*(fSkin[4]+fEdge[4])+0.8660254037844386*fSkin[2]-0.8660254037844386*fEdge[2]+0.5*(fSkin[0]+fEdge[0]); 
  fUpwindQuad[2] = copysignf(1.0,0.6324555320336759*alpha[2]+0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0])*(1.5*fSkin[7]-1.5*fEdge[7]+0.7745966692414833*(fSkin[6]+fEdge[6])+1.118033988749895*fSkin[5]-1.118033988749895*fEdge[5]+0.4472135954999579*fSkin[4]-0.4472135954999579*fEdge[4]+1.161895003862225*(fSkin[3]+fEdge[3])+0.8660254037844386*(fSkin[2]+fEdge[2])+0.6708203932499369*fSkin[1]-0.6708203932499369*fEdge[1]+0.5*fSkin[0]-0.5*fEdge[0])+1.5*(fSkin[7]+fEdge[7])+0.7745966692414833*fSkin[6]-0.7745966692414833*fEdge[6]+1.118033988749895*(fSkin[5]+fEdge[5])+0.4472135954999579*(fSkin[4]+fEdge[4])+1.161895003862225*fSkin[3]-1.161895003862225*fEdge[3]+0.8660254037844386*fSkin[2]-0.8660254037844386*fEdge[2]+0.6708203932499369*(fSkin[1]+fEdge[1])+0.5*(fSkin[0]+fEdge[0]); 

  fUpwind[0] = 0.0392837100659193*(5.0*fUpwindQuad[2]+8.0*fUpwindQuad[1]+5.0*fUpwindQuad[0]); 
  fUpwind[1] = 0.2635231383473649*(fUpwindQuad[2]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.1756820922315766*(fUpwindQuad[2]-2.0*fUpwindQuad[1]+fUpwindQuad[0]); 

  Ghat[0] = 0.7071067811865475*(alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1414213562373095*(4.47213595499958*(alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2])+5.0*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1])); 
  Ghat[2] = 0.02020305089104421*((22.3606797749979*alpha[2]+35.0*alpha[0])*fUpwind[2]+35.0*fUpwind[0]*alpha[2]+31.30495168499706*alpha[1]*fUpwind[1]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -1.224744871391589*Ghat[1]*dv10; 
  out[4] += -0.7071067811865475*Ghat[2]*dv10; 
  out[5] += -1.58113883008419*Ghat[0]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.58113883008419*Ghat[1]*dv10; 

  } else { 

  fUpwindQuad[0] = copysignf(1.0,0.6324555320336759*alpha[2]-0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0])*(1.5*fSkin[7]-1.5*fEdge[7]+0.7745966692414833*(fSkin[6]+fEdge[6])-1.118033988749895*fSkin[5]+1.118033988749895*fEdge[5]-0.4472135954999579*fSkin[4]+0.4472135954999579*fEdge[4]-1.161895003862225*(fSkin[3]+fEdge[3])+0.8660254037844386*(fSkin[2]+fEdge[2])+0.6708203932499369*fSkin[1]-0.6708203932499369*fEdge[1]-0.5*fSkin[0]+0.5*fEdge[0])-1.5*(fSkin[7]+fEdge[7])-0.7745966692414833*fSkin[6]+0.7745966692414833*fEdge[6]+1.118033988749895*(fSkin[5]+fEdge[5])+0.4472135954999579*(fSkin[4]+fEdge[4])+1.161895003862225*fSkin[3]-1.161895003862225*fEdge[3]-0.8660254037844386*fSkin[2]+0.8660254037844386*fEdge[2]-0.6708203932499369*(fSkin[1]+fEdge[1])+0.5*(fSkin[0]+fEdge[0]); 
  fUpwindQuad[1] = copysignf(1.0,0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2])*((-0.9682458365518543*(fSkin[6]+fEdge[6]))-1.118033988749895*fSkin[5]+1.118033988749895*fEdge[5]+0.5590169943749475*fSkin[4]-0.5590169943749475*fEdge[4]+0.8660254037844386*(fSkin[2]+fEdge[2])-0.5*fSkin[0]+0.5*fEdge[0])+0.9682458365518543*fSkin[6]-0.9682458365518543*fEdge[6]+1.118033988749895*(fSkin[5]+fEdge[5])-0.5590169943749475*(fSkin[4]+fEdge[4])-0.8660254037844386*fSkin[2]+0.8660254037844386*fEdge[2]+0.5*(fSkin[0]+fEdge[0]); 
  fUpwindQuad[2] = 1.5*fSkin[7]+copysignf(1.0,0.6324555320336759*alpha[2]+0.9486832980505137*alpha[1]+0.7071067811865475*alpha[0])*((-1.5*fSkin[7])+1.5*fEdge[7]+0.7745966692414833*(fSkin[6]+fEdge[6])-1.118033988749895*fSkin[5]+1.118033988749895*fEdge[5]-0.4472135954999579*fSkin[4]+0.4472135954999579*fEdge[4]+1.161895003862225*(fSkin[3]+fEdge[3])+0.8660254037844386*(fSkin[2]+fEdge[2])-0.6708203932499369*fSkin[1]+0.6708203932499369*fEdge[1]-0.5*fSkin[0]+0.5*fEdge[0])+1.5*fEdge[7]-0.7745966692414833*fSkin[6]+0.7745966692414833*fEdge[6]+1.118033988749895*(fSkin[5]+fEdge[5])+0.4472135954999579*(fSkin[4]+fEdge[4])-1.161895003862225*fSkin[3]+1.161895003862225*fEdge[3]-0.8660254037844386*fSkin[2]+0.8660254037844386*fEdge[2]+0.6708203932499369*(fSkin[1]+fEdge[1])+0.5*(fSkin[0]+fEdge[0]); 

  fUpwind[0] = 0.0392837100659193*(5.0*fUpwindQuad[2]+8.0*fUpwindQuad[1]+5.0*fUpwindQuad[0]); 
  fUpwind[1] = 0.2635231383473649*(fUpwindQuad[2]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.1756820922315766*(fUpwindQuad[2]-2.0*fUpwindQuad[1]+fUpwindQuad[0]); 

  Ghat[0] = 0.7071067811865475*(alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1414213562373095*(4.47213595499958*(alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2])+5.0*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1])); 
  Ghat[2] = 0.02020305089104421*((22.3606797749979*alpha[2]+35.0*alpha[0])*fUpwind[2]+35.0*fUpwind[0]*alpha[2]+31.30495168499706*alpha[1]*fUpwind[1]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -1.224744871391589*Ghat[1]*dv10; 
  out[4] += 0.7071067811865475*Ghat[2]*dv10; 
  out[5] += 1.58113883008419*Ghat[0]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += 1.58113883008419*Ghat[1]*dv10; 

  } 
} 
