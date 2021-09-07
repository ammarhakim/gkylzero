#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_1x2v_p1_surfvy_quad.h> 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[2]; 
  double alpha[4] = {0.0}; 

  alpha[0] = -2.449489742783178*A1[1]*dx10*wv1; 
  alpha[2] = -0.7071067811865475*A1[1]*dv1*dx10; 

  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[0] = ser_1x2v_p1_surfvy_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_1x2v_p1_surfvy_quad_0(-1, fEdge); 
  } 
  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[1] = ser_1x2v_p1_surfvy_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_1x2v_p1_surfvy_quad_1(-1, fEdge); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x2v_p1_surfvy_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_1x2v_p1_surfvy_quad_2(-1, fEdge); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x2v_p1_surfvy_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_1x2v_p1_surfvy_quad_3(-1, fEdge); 
  } 

  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] += 0.5*(alpha[2]*fUpwind[2]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.5*(alpha[2]*fUpwind[3]+alpha[0]*fUpwind[1]); 
  Ghat[2] += 0.5*(alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.5*(alpha[0]*fUpwind[3]+fUpwind[1]*alpha[2]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += -0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -1.224744871391589*Ghat[1]*dv11; 
  out[6] += -1.224744871391589*Ghat[2]*dv11; 
  out[7] += -1.224744871391589*Ghat[3]*dv11; 

  } else { 

  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[0] = ser_1x2v_p1_surfvy_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = ser_1x2v_p1_surfvy_quad_0(-1, fSkin); 
  } 
  if (alpha[0]-alpha[2] > 0) { 

    fUpwindQuad[1] = ser_1x2v_p1_surfvy_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = ser_1x2v_p1_surfvy_quad_1(-1, fSkin); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[2] = ser_1x2v_p1_surfvy_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = ser_1x2v_p1_surfvy_quad_2(-1, fSkin); 
  } 
  if (alpha[2]+alpha[0] > 0) { 

    fUpwindQuad[3] = ser_1x2v_p1_surfvy_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = ser_1x2v_p1_surfvy_quad_3(-1, fSkin); 
  } 

  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Ghat[0] += 0.5*(alpha[2]*fUpwind[2]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.5*(alpha[2]*fUpwind[3]+alpha[0]*fUpwind[1]); 
  Ghat[2] += 0.5*(alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.5*(alpha[0]*fUpwind[3]+fUpwind[1]*alpha[2]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += 0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -1.224744871391589*Ghat[1]*dv11; 
  out[6] += -1.224744871391589*Ghat[2]*dv11; 
  out[7] += -1.224744871391589*Ghat[3]*dv11; 

  } 
} 
