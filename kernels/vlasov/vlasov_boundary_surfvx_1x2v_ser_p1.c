#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_1x2v_p1_surfvx_quad.h> 
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

  double alpha[4] = {0.0}; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 

  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  if (alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[0] = ser_1x2v_p1_surfvx_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_1x2v_p1_surfvx_quad_0(-1, fEdge); 
  } 
  if ((-alpha[3])-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = ser_1x2v_p1_surfvx_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_1x2v_p1_surfvx_quad_1(-1, fEdge); 
  } 
  if ((-alpha[3])+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x2v_p1_surfvx_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_1x2v_p1_surfvx_quad_2(-1, fEdge); 
  } 
  if (alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x2v_p1_surfvx_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_1x2v_p1_surfvx_quad_3(-1, fEdge); 
  } 

  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] += 0.5*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 

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

    fUpwindQuad[0] = ser_1x2v_p1_surfvx_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = ser_1x2v_p1_surfvx_quad_0(-1, fSkin); 
  } 
  if ((-alpha[3])-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[1] = ser_1x2v_p1_surfvx_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = ser_1x2v_p1_surfvx_quad_1(-1, fSkin); 
  } 
  if ((-alpha[3])+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x2v_p1_surfvx_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = ser_1x2v_p1_surfvx_quad_2(-1, fSkin); 
  } 
  if (alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x2v_p1_surfvx_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = ser_1x2v_p1_surfvx_quad_3(-1, fSkin); 
  } 

  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] += 0.5*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 

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
