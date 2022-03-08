#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx3_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind.h> 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  double alpha[8] = {0.0}; 

  alpha[0] = -2.449489742783178*phi[1]*dx10; 
  alpha[2] = -2.449489742783178*phi[3]*dx10; 

  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_quad_0_r(fSkin); 
    fUpwindQuad[4] = ser_4x_p1_surfx3_quad_4_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_quad_0_l(fEdge); 
    fUpwindQuad[4] = ser_4x_p1_surfx3_quad_4_l(fEdge); 
  } 
  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_quad_1_r(fSkin); 
    fUpwindQuad[5] = ser_4x_p1_surfx3_quad_5_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_quad_1_l(fEdge); 
    fUpwindQuad[5] = ser_4x_p1_surfx3_quad_5_l(fEdge); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_quad_2_r(fSkin); 
    fUpwindQuad[6] = ser_4x_p1_surfx3_quad_6_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_quad_2_l(fEdge); 
    fUpwindQuad[6] = ser_4x_p1_surfx3_quad_6_l(fEdge); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_quad_3_r(fSkin); 
    fUpwindQuad[7] = ser_4x_p1_surfx3_quad_7_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_quad_3_l(fEdge); 
    fUpwindQuad[7] = ser_4x_p1_surfx3_quad_7_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.3535533905932737*(alpha[2]*fUpwind[2]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*(alpha[2]*fUpwind[4]+alpha[0]*fUpwind[1]); 
  Ghat[2] += 0.3535533905932737*(alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[2]*fUpwind[6]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.3535533905932737*(alpha[0]*fUpwind[4]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3535533905932737*(alpha[2]*fUpwind[7]+alpha[0]*fUpwind[5]); 
  Ghat[6] += 0.3535533905932737*(alpha[0]*fUpwind[6]+alpha[2]*fUpwind[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[2]*fUpwind[5]); 

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

  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_quad_0_r(fEdge); 
    fUpwindQuad[4] = ser_4x_p1_surfx3_quad_4_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_quad_0_l(fSkin); 
    fUpwindQuad[4] = ser_4x_p1_surfx3_quad_4_l(fSkin); 
  } 
  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_quad_1_r(fEdge); 
    fUpwindQuad[5] = ser_4x_p1_surfx3_quad_5_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_quad_1_l(fSkin); 
    fUpwindQuad[5] = ser_4x_p1_surfx3_quad_5_l(fSkin); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_quad_2_r(fEdge); 
    fUpwindQuad[6] = ser_4x_p1_surfx3_quad_6_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_quad_2_l(fSkin); 
    fUpwindQuad[6] = ser_4x_p1_surfx3_quad_6_l(fSkin); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_quad_3_r(fEdge); 
    fUpwindQuad[7] = ser_4x_p1_surfx3_quad_7_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_quad_3_l(fSkin); 
    fUpwindQuad[7] = ser_4x_p1_surfx3_quad_7_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.3535533905932737*(alpha[2]*fUpwind[2]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*(alpha[2]*fUpwind[4]+alpha[0]*fUpwind[1]); 
  Ghat[2] += 0.3535533905932737*(alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3535533905932737*(alpha[2]*fUpwind[6]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.3535533905932737*(alpha[0]*fUpwind[4]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3535533905932737*(alpha[2]*fUpwind[7]+alpha[0]*fUpwind[5]); 
  Ghat[6] += 0.3535533905932737*(alpha[0]*fUpwind[6]+alpha[2]*fUpwind[3]); 
  Ghat[7] += 0.3535533905932737*(alpha[0]*fUpwind[7]+alpha[2]*fUpwind[5]); 

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
