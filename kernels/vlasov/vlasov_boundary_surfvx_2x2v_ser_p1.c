#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[20]; 

  double Ghat[8]; 
  double alpha[8]; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 1.414213562373095*(B2[2]*wv2+E0[2]); 
  alpha[3] = 0.408248290463863*B2[0]*dv2; 
  alpha[4] = 1.414213562373095*(B2[3]*wv2+E0[3]); 
  alpha[5] = 0.408248290463863*B2[1]*dv2; 
  alpha[6] = 0.408248290463863*B2[2]*dv2; 
  alpha[7] = 0.408248290463863*B2[3]*dv2; 

  double fUpwindQuad[8];
  double fUpwind[8];

  if (edge == -1) { 

  if ((-alpha[7])+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*(fSkin[14]+fSkin[13])-0.25*fSkin[12]+0.4330127018922193*fSkin[11]-0.4330127018922193*fSkin[10]+0.25*(fSkin[9]+fSkin[8])-0.4330127018922193*(fSkin[7]+fSkin[6])+0.25*fSkin[5]-0.25*fSkin[4]+0.4330127018922193*fSkin[3]-0.25*(fSkin[2]+fSkin[1])+0.25*fSkin[0]; 
  } else { 

    fUpwindQuad[0] = 0.4330127018922193*fEdge[15]-0.4330127018922193*(fEdge[14]+fEdge[13])-0.25*fEdge[12]-0.4330127018922193*fEdge[11]+0.4330127018922193*fEdge[10]+0.25*(fEdge[9]+fEdge[8])+0.4330127018922193*(fEdge[7]+fEdge[6])+0.25*fEdge[5]-0.25*fEdge[4]-0.4330127018922193*fEdge[3]-0.25*(fEdge[2]+fEdge[1])+0.25*fEdge[0]; 
  } 
  if (alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = 0.4330127018922193*(fSkin[15]+fSkin[14])-0.4330127018922193*fSkin[13]+0.25*fSkin[12]-0.4330127018922193*(fSkin[11]+fSkin[10])+0.25*fSkin[9]-0.25*fSkin[8]-0.4330127018922193*fSkin[7]+0.4330127018922193*fSkin[6]-0.25*(fSkin[5]+fSkin[4])+0.4330127018922193*fSkin[3]-0.25*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[1] = (-0.4330127018922193*(fEdge[15]+fEdge[14]))+0.4330127018922193*fEdge[13]+0.25*fEdge[12]+0.4330127018922193*(fEdge[11]+fEdge[10])+0.25*fEdge[9]-0.25*fEdge[8]+0.4330127018922193*fEdge[7]-0.4330127018922193*fEdge[6]-0.25*(fEdge[5]+fEdge[4])-0.4330127018922193*fEdge[3]-0.25*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } 
  if (alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = 0.4330127018922193*fSkin[15]-0.4330127018922193*fSkin[14]+0.4330127018922193*fSkin[13]+0.25*fSkin[12]-0.4330127018922193*(fSkin[11]+fSkin[10])-0.25*fSkin[9]+0.25*fSkin[8]+0.4330127018922193*fSkin[7]-0.4330127018922193*fSkin[6]-0.25*(fSkin[5]+fSkin[4])+0.4330127018922193*fSkin[3]+0.25*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } else { 

    fUpwindQuad[2] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*fEdge[14]-0.4330127018922193*fEdge[13]+0.25*fEdge[12]+0.4330127018922193*(fEdge[11]+fEdge[10])-0.25*fEdge[9]+0.25*fEdge[8]-0.4330127018922193*fEdge[7]+0.4330127018922193*fEdge[6]-0.25*(fEdge[5]+fEdge[4])-0.4330127018922193*fEdge[3]+0.25*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } 
  if ((-alpha[7])-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = (-0.4330127018922193*(fSkin[15]+fSkin[14]+fSkin[13]))-0.25*fSkin[12]+0.4330127018922193*fSkin[11]-0.4330127018922193*fSkin[10]-0.25*(fSkin[9]+fSkin[8])+0.4330127018922193*(fSkin[7]+fSkin[6])+0.25*fSkin[5]-0.25*fSkin[4]+0.4330127018922193*fSkin[3]+0.25*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[3] = 0.4330127018922193*(fEdge[15]+fEdge[14]+fEdge[13])-0.25*fEdge[12]-0.4330127018922193*fEdge[11]+0.4330127018922193*fEdge[10]-0.25*(fEdge[9]+fEdge[8])-0.4330127018922193*(fEdge[7]+fEdge[6])+0.25*fEdge[5]-0.25*fEdge[4]-0.4330127018922193*fEdge[3]+0.25*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } 
  if (alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = 0.4330127018922193*fSkin[15]-0.4330127018922193*(fSkin[14]+fSkin[13])+0.25*fSkin[12]+0.4330127018922193*(fSkin[11]+fSkin[10])-0.25*(fSkin[9]+fSkin[8])-0.4330127018922193*(fSkin[7]+fSkin[6])+0.25*(fSkin[5]+fSkin[4])+0.4330127018922193*fSkin[3]-0.25*(fSkin[2]+fSkin[1])+0.25*fSkin[0]; 
  } else { 

    fUpwindQuad[4] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*(fEdge[14]+fEdge[13])+0.25*fEdge[12]-0.4330127018922193*(fEdge[11]+fEdge[10])-0.25*(fEdge[9]+fEdge[8])+0.4330127018922193*(fEdge[7]+fEdge[6])+0.25*(fEdge[5]+fEdge[4])-0.4330127018922193*fEdge[3]-0.25*(fEdge[2]+fEdge[1])+0.25*fEdge[0]; 
  } 
  if ((-alpha[7])-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = (-0.4330127018922193*(fSkin[15]+fSkin[14]))+0.4330127018922193*fSkin[13]-0.25*fSkin[12]-0.4330127018922193*fSkin[11]+0.4330127018922193*fSkin[10]-0.25*fSkin[9]+0.25*fSkin[8]-0.4330127018922193*fSkin[7]+0.4330127018922193*fSkin[6]-0.25*fSkin[5]+0.25*fSkin[4]+0.4330127018922193*fSkin[3]-0.25*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[5] = 0.4330127018922193*(fEdge[15]+fEdge[14])-0.4330127018922193*fEdge[13]-0.25*fEdge[12]+0.4330127018922193*fEdge[11]-0.4330127018922193*fEdge[10]-0.25*fEdge[9]+0.25*fEdge[8]+0.4330127018922193*fEdge[7]-0.4330127018922193*fEdge[6]-0.25*fEdge[5]+0.25*fEdge[4]-0.4330127018922193*fEdge[3]-0.25*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } 
  if ((-alpha[7])+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*fSkin[14]-0.4330127018922193*fSkin[13]-0.25*fSkin[12]-0.4330127018922193*fSkin[11]+0.4330127018922193*fSkin[10]+0.25*fSkin[9]-0.25*fSkin[8]+0.4330127018922193*fSkin[7]-0.4330127018922193*fSkin[6]-0.25*fSkin[5]+0.25*fSkin[4]+0.4330127018922193*fSkin[3]+0.25*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } else { 

    fUpwindQuad[6] = 0.4330127018922193*fEdge[15]-0.4330127018922193*fEdge[14]+0.4330127018922193*fEdge[13]-0.25*fEdge[12]+0.4330127018922193*fEdge[11]-0.4330127018922193*fEdge[10]+0.25*fEdge[9]-0.25*fEdge[8]-0.4330127018922193*fEdge[7]+0.4330127018922193*fEdge[6]-0.25*fEdge[5]+0.25*fEdge[4]-0.4330127018922193*fEdge[3]+0.25*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } 
  if (alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[7] = 0.4330127018922193*(fSkin[15]+fSkin[14]+fSkin[13])+0.25*fSkin[12]+0.4330127018922193*(fSkin[11]+fSkin[10])+0.25*(fSkin[9]+fSkin[8])+0.4330127018922193*(fSkin[7]+fSkin[6])+0.25*(fSkin[5]+fSkin[4])+0.4330127018922193*fSkin[3]+0.25*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } else { 

    fUpwindQuad[7] = (-0.4330127018922193*(fEdge[15]+fEdge[14]+fEdge[13]))+0.25*fEdge[12]-0.4330127018922193*(fEdge[11]+fEdge[10])+0.25*(fEdge[9]+fEdge[8])-0.4330127018922193*(fEdge[7]+fEdge[6])+0.25*(fEdge[5]+fEdge[4])-0.4330127018922193*fEdge[3]+0.25*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } 

  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[5] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[6] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Ghat[0] = 0.3535533905932737*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] = 0.3535533905932737*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] = 0.3535533905932737*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -1.224744871391589*Ghat[1]*dv10; 
  out[7] += -1.224744871391589*Ghat[2]*dv10; 
  out[8] += -0.7071067811865475*Ghat[5]*dv10; 
  out[9] += -0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += -0.7071067811865475*Ghat[7]*dv10; 
  out[13] += -1.224744871391589*Ghat[5]*dv10; 
  out[14] += -1.224744871391589*Ghat[6]*dv10; 
  out[15] += -1.224744871391589*Ghat[7]*dv10; 

  } else { 

  if ((-alpha[7])+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*(fEdge[14]+fEdge[13])-0.25*fEdge[12]+0.4330127018922193*fEdge[11]-0.4330127018922193*fEdge[10]+0.25*(fEdge[9]+fEdge[8])-0.4330127018922193*(fEdge[7]+fEdge[6])+0.25*fEdge[5]-0.25*fEdge[4]+0.4330127018922193*fEdge[3]-0.25*(fEdge[2]+fEdge[1])+0.25*fEdge[0]; 
  } else { 

    fUpwindQuad[0] = 0.4330127018922193*fSkin[15]-0.4330127018922193*(fSkin[14]+fSkin[13])-0.25*fSkin[12]-0.4330127018922193*fSkin[11]+0.4330127018922193*fSkin[10]+0.25*(fSkin[9]+fSkin[8])+0.4330127018922193*(fSkin[7]+fSkin[6])+0.25*fSkin[5]-0.25*fSkin[4]-0.4330127018922193*fSkin[3]-0.25*(fSkin[2]+fSkin[1])+0.25*fSkin[0]; 
  } 
  if (alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = 0.4330127018922193*(fEdge[15]+fEdge[14])-0.4330127018922193*fEdge[13]+0.25*fEdge[12]-0.4330127018922193*(fEdge[11]+fEdge[10])+0.25*fEdge[9]-0.25*fEdge[8]-0.4330127018922193*fEdge[7]+0.4330127018922193*fEdge[6]-0.25*(fEdge[5]+fEdge[4])+0.4330127018922193*fEdge[3]-0.25*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[1] = (-0.4330127018922193*(fSkin[15]+fSkin[14]))+0.4330127018922193*fSkin[13]+0.25*fSkin[12]+0.4330127018922193*(fSkin[11]+fSkin[10])+0.25*fSkin[9]-0.25*fSkin[8]+0.4330127018922193*fSkin[7]-0.4330127018922193*fSkin[6]-0.25*(fSkin[5]+fSkin[4])-0.4330127018922193*fSkin[3]-0.25*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } 
  if (alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = 0.4330127018922193*fEdge[15]-0.4330127018922193*fEdge[14]+0.4330127018922193*fEdge[13]+0.25*fEdge[12]-0.4330127018922193*(fEdge[11]+fEdge[10])-0.25*fEdge[9]+0.25*fEdge[8]+0.4330127018922193*fEdge[7]-0.4330127018922193*fEdge[6]-0.25*(fEdge[5]+fEdge[4])+0.4330127018922193*fEdge[3]+0.25*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } else { 

    fUpwindQuad[2] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*fSkin[14]-0.4330127018922193*fSkin[13]+0.25*fSkin[12]+0.4330127018922193*(fSkin[11]+fSkin[10])-0.25*fSkin[9]+0.25*fSkin[8]-0.4330127018922193*fSkin[7]+0.4330127018922193*fSkin[6]-0.25*(fSkin[5]+fSkin[4])-0.4330127018922193*fSkin[3]+0.25*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } 
  if ((-alpha[7])-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = (-0.4330127018922193*(fEdge[15]+fEdge[14]+fEdge[13]))-0.25*fEdge[12]+0.4330127018922193*fEdge[11]-0.4330127018922193*fEdge[10]-0.25*(fEdge[9]+fEdge[8])+0.4330127018922193*(fEdge[7]+fEdge[6])+0.25*fEdge[5]-0.25*fEdge[4]+0.4330127018922193*fEdge[3]+0.25*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[3] = 0.4330127018922193*(fSkin[15]+fSkin[14]+fSkin[13])-0.25*fSkin[12]-0.4330127018922193*fSkin[11]+0.4330127018922193*fSkin[10]-0.25*(fSkin[9]+fSkin[8])-0.4330127018922193*(fSkin[7]+fSkin[6])+0.25*fSkin[5]-0.25*fSkin[4]-0.4330127018922193*fSkin[3]+0.25*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } 
  if (alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[4] = 0.4330127018922193*fEdge[15]-0.4330127018922193*(fEdge[14]+fEdge[13])+0.25*fEdge[12]+0.4330127018922193*(fEdge[11]+fEdge[10])-0.25*(fEdge[9]+fEdge[8])-0.4330127018922193*(fEdge[7]+fEdge[6])+0.25*(fEdge[5]+fEdge[4])+0.4330127018922193*fEdge[3]-0.25*(fEdge[2]+fEdge[1])+0.25*fEdge[0]; 
  } else { 

    fUpwindQuad[4] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*(fSkin[14]+fSkin[13])+0.25*fSkin[12]-0.4330127018922193*(fSkin[11]+fSkin[10])-0.25*(fSkin[9]+fSkin[8])+0.4330127018922193*(fSkin[7]+fSkin[6])+0.25*(fSkin[5]+fSkin[4])-0.4330127018922193*fSkin[3]-0.25*(fSkin[2]+fSkin[1])+0.25*fSkin[0]; 
  } 
  if ((-alpha[7])-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[5] = (-0.4330127018922193*(fEdge[15]+fEdge[14]))+0.4330127018922193*fEdge[13]-0.25*fEdge[12]-0.4330127018922193*fEdge[11]+0.4330127018922193*fEdge[10]-0.25*fEdge[9]+0.25*fEdge[8]-0.4330127018922193*fEdge[7]+0.4330127018922193*fEdge[6]-0.25*fEdge[5]+0.25*fEdge[4]+0.4330127018922193*fEdge[3]-0.25*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[5] = 0.4330127018922193*(fSkin[15]+fSkin[14])-0.4330127018922193*fSkin[13]-0.25*fSkin[12]+0.4330127018922193*fSkin[11]-0.4330127018922193*fSkin[10]-0.25*fSkin[9]+0.25*fSkin[8]+0.4330127018922193*fSkin[7]-0.4330127018922193*fSkin[6]-0.25*fSkin[5]+0.25*fSkin[4]-0.4330127018922193*fSkin[3]-0.25*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } 
  if ((-alpha[7])+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[6] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*fEdge[14]-0.4330127018922193*fEdge[13]-0.25*fEdge[12]-0.4330127018922193*fEdge[11]+0.4330127018922193*fEdge[10]+0.25*fEdge[9]-0.25*fEdge[8]+0.4330127018922193*fEdge[7]-0.4330127018922193*fEdge[6]-0.25*fEdge[5]+0.25*fEdge[4]+0.4330127018922193*fEdge[3]+0.25*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } else { 

    fUpwindQuad[6] = 0.4330127018922193*fSkin[15]-0.4330127018922193*fSkin[14]+0.4330127018922193*fSkin[13]-0.25*fSkin[12]+0.4330127018922193*fSkin[11]-0.4330127018922193*fSkin[10]+0.25*fSkin[9]-0.25*fSkin[8]-0.4330127018922193*fSkin[7]+0.4330127018922193*fSkin[6]-0.25*fSkin[5]+0.25*fSkin[4]-0.4330127018922193*fSkin[3]+0.25*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } 
  if (alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[7] = 0.4330127018922193*(fEdge[15]+fEdge[14]+fEdge[13])+0.25*fEdge[12]+0.4330127018922193*(fEdge[11]+fEdge[10])+0.25*(fEdge[9]+fEdge[8])+0.4330127018922193*(fEdge[7]+fEdge[6])+0.25*(fEdge[5]+fEdge[4])+0.4330127018922193*fEdge[3]+0.25*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } else { 

    fUpwindQuad[7] = (-0.4330127018922193*(fSkin[15]+fSkin[14]+fSkin[13]))+0.25*fSkin[12]-0.4330127018922193*(fSkin[11]+fSkin[10])+0.25*(fSkin[9]+fSkin[8])-0.4330127018922193*(fSkin[7]+fSkin[6])+0.25*(fSkin[5]+fSkin[4])-0.4330127018922193*fSkin[3]+0.25*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } 

  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[5] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[6] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Ghat[0] = 0.3535533905932737*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] = 0.3535533905932737*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] = 0.3535533905932737*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -1.224744871391589*Ghat[1]*dv10; 
  out[7] += -1.224744871391589*Ghat[2]*dv10; 
  out[8] += 0.7071067811865475*Ghat[5]*dv10; 
  out[9] += 0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += 0.7071067811865475*Ghat[7]*dv10; 
  out[13] += -1.224744871391589*Ghat[5]*dv10; 
  out[14] += -1.224744871391589*Ghat[6]*dv10; 
  out[15] += -1.224744871391589*Ghat[7]*dv10; 

  } 
} 
