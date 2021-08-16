#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[10]; 

  double Ghat[4]; 
  double alpha[4]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 

  double fUpwindQuad[4];
  double fUpwind[4];

  if (edge == -1) { 

  if (alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = 0.6123724356957944*fSkin[7]-0.6123724356957944*fSkin[6]+0.3535533905932737*fSkin[5]-0.6123724356957944*fSkin[4]-0.3535533905932737*fSkin[3]+0.6123724356957944*fSkin[2]-0.3535533905932737*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[0] = (-0.6123724356957944*fEdge[7])+0.6123724356957944*fEdge[6]+0.3535533905932737*fEdge[5]+0.6123724356957944*fEdge[4]-0.3535533905932737*fEdge[3]-0.6123724356957944*fEdge[2]-0.3535533905932737*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 
  if ((-alpha[3])-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = (-0.6123724356957944*(fSkin[7]+fSkin[6]))-0.3535533905932737*fSkin[5]+0.6123724356957944*fSkin[4]-0.3535533905932737*fSkin[3]+0.6123724356957944*fSkin[2]+0.3535533905932737*(fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[1] = 0.6123724356957944*(fEdge[7]+fEdge[6])-0.3535533905932737*fEdge[5]-0.6123724356957944*fEdge[4]-0.3535533905932737*fEdge[3]-0.6123724356957944*fEdge[2]+0.3535533905932737*(fEdge[1]+fEdge[0]); 
  } 
  if ((-alpha[3])+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = (-0.6123724356957944*fSkin[7])+0.6123724356957944*fSkin[6]-0.3535533905932737*fSkin[5]-0.6123724356957944*fSkin[4]+0.3535533905932737*fSkin[3]+0.6123724356957944*fSkin[2]-0.3535533905932737*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 

    fUpwindQuad[2] = 0.6123724356957944*fEdge[7]-0.6123724356957944*fEdge[6]-0.3535533905932737*fEdge[5]+0.6123724356957944*fEdge[4]+0.3535533905932737*fEdge[3]-0.6123724356957944*fEdge[2]-0.3535533905932737*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 
  if (alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = 0.6123724356957944*(fSkin[7]+fSkin[6])+0.3535533905932737*fSkin[5]+0.6123724356957944*fSkin[4]+0.3535533905932737*fSkin[3]+0.6123724356957944*fSkin[2]+0.3535533905932737*(fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[3] = (-0.6123724356957944*(fEdge[7]+fEdge[6]))+0.3535533905932737*fEdge[5]-0.6123724356957944*fEdge[4]+0.3535533905932737*fEdge[3]-0.6123724356957944*fEdge[2]+0.3535533905932737*(fEdge[1]+fEdge[0]); 
  } 

  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] = 0.5*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += -0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10; 

  } else { 

  if (alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = 0.6123724356957944*fEdge[7]-0.6123724356957944*fEdge[6]+0.3535533905932737*fEdge[5]-0.6123724356957944*fEdge[4]-0.3535533905932737*fEdge[3]+0.6123724356957944*fEdge[2]-0.3535533905932737*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[0] = (-0.6123724356957944*fSkin[7])+0.6123724356957944*fSkin[6]+0.3535533905932737*fSkin[5]+0.6123724356957944*fSkin[4]-0.3535533905932737*fSkin[3]-0.6123724356957944*fSkin[2]-0.3535533905932737*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 
  if ((-alpha[3])-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = (-0.6123724356957944*(fEdge[7]+fEdge[6]))-0.3535533905932737*fEdge[5]+0.6123724356957944*fEdge[4]-0.3535533905932737*fEdge[3]+0.6123724356957944*fEdge[2]+0.3535533905932737*(fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[1] = 0.6123724356957944*(fSkin[7]+fSkin[6])-0.3535533905932737*fSkin[5]-0.6123724356957944*fSkin[4]-0.3535533905932737*fSkin[3]-0.6123724356957944*fSkin[2]+0.3535533905932737*(fSkin[1]+fSkin[0]); 
  } 
  if ((-alpha[3])+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = (-0.6123724356957944*fEdge[7])+0.6123724356957944*fEdge[6]-0.3535533905932737*fEdge[5]-0.6123724356957944*fEdge[4]+0.3535533905932737*fEdge[3]+0.6123724356957944*fEdge[2]-0.3535533905932737*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 

    fUpwindQuad[2] = 0.6123724356957944*fSkin[7]-0.6123724356957944*fSkin[6]-0.3535533905932737*fSkin[5]+0.6123724356957944*fSkin[4]+0.3535533905932737*fSkin[3]-0.6123724356957944*fSkin[2]-0.3535533905932737*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 
  if (alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = 0.6123724356957944*(fEdge[7]+fEdge[6])+0.3535533905932737*fEdge[5]+0.6123724356957944*fEdge[4]+0.3535533905932737*fEdge[3]+0.6123724356957944*fEdge[2]+0.3535533905932737*(fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[3] = (-0.6123724356957944*(fSkin[7]+fSkin[6]))+0.3535533905932737*fSkin[5]-0.6123724356957944*fSkin[4]+0.3535533905932737*fSkin[3]-0.6123724356957944*fSkin[2]+0.3535533905932737*(fSkin[1]+fSkin[0]); 
  } 

  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] = 0.5*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += 0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10; 

  } 
} 
