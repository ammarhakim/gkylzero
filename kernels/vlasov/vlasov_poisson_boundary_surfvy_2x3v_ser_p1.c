#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p1_surfx4_quad.h> 
#include <gkyl_basis_ser_5x_p1_upwind.h> 
GKYL_CU_DH void vlasov_poisson_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  double alpha[16] = {0.0}; 

  alpha[0] = -3.464101615137754*phi[2]*dx11; 
  alpha[1] = -3.464101615137754*phi[3]*dx11; 

  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_quad_0_r(fSkin); 
    fUpwindQuad[4] = ser_5x_p1_surfx4_quad_4_r(fSkin); 
    fUpwindQuad[8] = ser_5x_p1_surfx4_quad_8_r(fSkin); 
    fUpwindQuad[12] = ser_5x_p1_surfx4_quad_12_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_quad_0_l(fEdge); 
    fUpwindQuad[4] = ser_5x_p1_surfx4_quad_4_l(fEdge); 
    fUpwindQuad[8] = ser_5x_p1_surfx4_quad_8_l(fEdge); 
    fUpwindQuad[12] = ser_5x_p1_surfx4_quad_12_l(fEdge); 
  } 
  if (alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_quad_1_r(fSkin); 
    fUpwindQuad[5] = ser_5x_p1_surfx4_quad_5_r(fSkin); 
    fUpwindQuad[9] = ser_5x_p1_surfx4_quad_9_r(fSkin); 
    fUpwindQuad[13] = ser_5x_p1_surfx4_quad_13_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_quad_1_l(fEdge); 
    fUpwindQuad[5] = ser_5x_p1_surfx4_quad_5_l(fEdge); 
    fUpwindQuad[9] = ser_5x_p1_surfx4_quad_9_l(fEdge); 
    fUpwindQuad[13] = ser_5x_p1_surfx4_quad_13_l(fEdge); 
  } 
  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_quad_2_r(fSkin); 
    fUpwindQuad[6] = ser_5x_p1_surfx4_quad_6_r(fSkin); 
    fUpwindQuad[10] = ser_5x_p1_surfx4_quad_10_r(fSkin); 
    fUpwindQuad[14] = ser_5x_p1_surfx4_quad_14_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_quad_2_l(fEdge); 
    fUpwindQuad[6] = ser_5x_p1_surfx4_quad_6_l(fEdge); 
    fUpwindQuad[10] = ser_5x_p1_surfx4_quad_10_l(fEdge); 
    fUpwindQuad[14] = ser_5x_p1_surfx4_quad_14_l(fEdge); 
  } 
  if (alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_quad_3_r(fSkin); 
    fUpwindQuad[7] = ser_5x_p1_surfx4_quad_7_r(fSkin); 
    fUpwindQuad[11] = ser_5x_p1_surfx4_quad_11_r(fSkin); 
    fUpwindQuad[15] = ser_5x_p1_surfx4_quad_15_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_quad_3_l(fEdge); 
    fUpwindQuad[7] = ser_5x_p1_surfx4_quad_7_l(fEdge); 
    fUpwindQuad[11] = ser_5x_p1_surfx4_quad_11_l(fEdge); 
    fUpwindQuad[15] = ser_5x_p1_surfx4_quad_15_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.25*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.25*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.25*(alpha[1]*fUpwind[5]+alpha[0]*fUpwind[2]); 
  Ghat[3] += 0.25*(alpha[1]*fUpwind[6]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.25*(alpha[1]*fUpwind[8]+alpha[0]*fUpwind[4]); 
  Ghat[5] += 0.25*(alpha[0]*fUpwind[5]+alpha[1]*fUpwind[2]); 
  Ghat[6] += 0.25*(alpha[0]*fUpwind[6]+alpha[1]*fUpwind[3]); 
  Ghat[7] += 0.25*(alpha[1]*fUpwind[11]+alpha[0]*fUpwind[7]); 
  Ghat[8] += 0.25*(alpha[0]*fUpwind[8]+alpha[1]*fUpwind[4]); 
  Ghat[9] += 0.25*(alpha[1]*fUpwind[12]+alpha[0]*fUpwind[9]); 
  Ghat[10] += 0.25*(alpha[1]*fUpwind[13]+alpha[0]*fUpwind[10]); 
  Ghat[11] += 0.25*(alpha[0]*fUpwind[11]+alpha[1]*fUpwind[7]); 
  Ghat[12] += 0.25*(alpha[0]*fUpwind[12]+alpha[1]*fUpwind[9]); 
  Ghat[13] += 0.25*(alpha[0]*fUpwind[13]+alpha[1]*fUpwind[10]); 
  Ghat[14] += 0.25*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]); 
  Ghat[15] += 0.25*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += -0.7071067811865475*Ghat[4]*dv11; 
  out[6] += -0.7071067811865475*Ghat[5]*dv11; 
  out[7] += -0.7071067811865475*Ghat[6]*dv11; 
  out[8] += -0.7071067811865475*Ghat[7]*dv11; 
  out[9] += -1.224744871391589*Ghat[1]*dv11; 
  out[10] += -1.224744871391589*Ghat[2]*dv11; 
  out[11] += -1.224744871391589*Ghat[3]*dv11; 
  out[12] += -0.7071067811865475*Ghat[8]*dv11; 
  out[13] += -0.7071067811865475*Ghat[9]*dv11; 
  out[14] += -0.7071067811865475*Ghat[10]*dv11; 
  out[15] += -1.224744871391589*Ghat[4]*dv11; 
  out[16] += -0.7071067811865475*Ghat[11]*dv11; 
  out[17] += -1.224744871391589*Ghat[5]*dv11; 
  out[18] += -1.224744871391589*Ghat[6]*dv11; 
  out[19] += -1.224744871391589*Ghat[7]*dv11; 
  out[20] += -0.7071067811865475*Ghat[12]*dv11; 
  out[21] += -0.7071067811865475*Ghat[13]*dv11; 
  out[22] += -0.7071067811865475*Ghat[14]*dv11; 
  out[23] += -1.224744871391589*Ghat[8]*dv11; 
  out[24] += -1.224744871391589*Ghat[9]*dv11; 
  out[25] += -1.224744871391589*Ghat[10]*dv11; 
  out[26] += -1.224744871391589*Ghat[11]*dv11; 
  out[27] += -0.7071067811865475*Ghat[15]*dv11; 
  out[28] += -1.224744871391589*Ghat[12]*dv11; 
  out[29] += -1.224744871391589*Ghat[13]*dv11; 
  out[30] += -1.224744871391589*Ghat[14]*dv11; 
  out[31] += -1.224744871391589*Ghat[15]*dv11; 

  } else { 

  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_quad_0_r(fEdge); 
    fUpwindQuad[4] = ser_5x_p1_surfx4_quad_4_r(fEdge); 
    fUpwindQuad[8] = ser_5x_p1_surfx4_quad_8_r(fEdge); 
    fUpwindQuad[12] = ser_5x_p1_surfx4_quad_12_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_quad_0_l(fSkin); 
    fUpwindQuad[4] = ser_5x_p1_surfx4_quad_4_l(fSkin); 
    fUpwindQuad[8] = ser_5x_p1_surfx4_quad_8_l(fSkin); 
    fUpwindQuad[12] = ser_5x_p1_surfx4_quad_12_l(fSkin); 
  } 
  if (alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_quad_1_r(fEdge); 
    fUpwindQuad[5] = ser_5x_p1_surfx4_quad_5_r(fEdge); 
    fUpwindQuad[9] = ser_5x_p1_surfx4_quad_9_r(fEdge); 
    fUpwindQuad[13] = ser_5x_p1_surfx4_quad_13_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_quad_1_l(fSkin); 
    fUpwindQuad[5] = ser_5x_p1_surfx4_quad_5_l(fSkin); 
    fUpwindQuad[9] = ser_5x_p1_surfx4_quad_9_l(fSkin); 
    fUpwindQuad[13] = ser_5x_p1_surfx4_quad_13_l(fSkin); 
  } 
  if (alpha[0]-alpha[1] > 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_quad_2_r(fEdge); 
    fUpwindQuad[6] = ser_5x_p1_surfx4_quad_6_r(fEdge); 
    fUpwindQuad[10] = ser_5x_p1_surfx4_quad_10_r(fEdge); 
    fUpwindQuad[14] = ser_5x_p1_surfx4_quad_14_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_quad_2_l(fSkin); 
    fUpwindQuad[6] = ser_5x_p1_surfx4_quad_6_l(fSkin); 
    fUpwindQuad[10] = ser_5x_p1_surfx4_quad_10_l(fSkin); 
    fUpwindQuad[14] = ser_5x_p1_surfx4_quad_14_l(fSkin); 
  } 
  if (alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_quad_3_r(fEdge); 
    fUpwindQuad[7] = ser_5x_p1_surfx4_quad_7_r(fEdge); 
    fUpwindQuad[11] = ser_5x_p1_surfx4_quad_11_r(fEdge); 
    fUpwindQuad[15] = ser_5x_p1_surfx4_quad_15_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_quad_3_l(fSkin); 
    fUpwindQuad[7] = ser_5x_p1_surfx4_quad_7_l(fSkin); 
    fUpwindQuad[11] = ser_5x_p1_surfx4_quad_11_l(fSkin); 
    fUpwindQuad[15] = ser_5x_p1_surfx4_quad_15_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.25*(alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.25*(alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.25*(alpha[1]*fUpwind[5]+alpha[0]*fUpwind[2]); 
  Ghat[3] += 0.25*(alpha[1]*fUpwind[6]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.25*(alpha[1]*fUpwind[8]+alpha[0]*fUpwind[4]); 
  Ghat[5] += 0.25*(alpha[0]*fUpwind[5]+alpha[1]*fUpwind[2]); 
  Ghat[6] += 0.25*(alpha[0]*fUpwind[6]+alpha[1]*fUpwind[3]); 
  Ghat[7] += 0.25*(alpha[1]*fUpwind[11]+alpha[0]*fUpwind[7]); 
  Ghat[8] += 0.25*(alpha[0]*fUpwind[8]+alpha[1]*fUpwind[4]); 
  Ghat[9] += 0.25*(alpha[1]*fUpwind[12]+alpha[0]*fUpwind[9]); 
  Ghat[10] += 0.25*(alpha[1]*fUpwind[13]+alpha[0]*fUpwind[10]); 
  Ghat[11] += 0.25*(alpha[0]*fUpwind[11]+alpha[1]*fUpwind[7]); 
  Ghat[12] += 0.25*(alpha[0]*fUpwind[12]+alpha[1]*fUpwind[9]); 
  Ghat[13] += 0.25*(alpha[0]*fUpwind[13]+alpha[1]*fUpwind[10]); 
  Ghat[14] += 0.25*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]); 
  Ghat[15] += 0.25*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += 0.7071067811865475*Ghat[3]*dv11; 
  out[4] += -1.224744871391589*Ghat[0]*dv11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11; 
  out[6] += 0.7071067811865475*Ghat[5]*dv11; 
  out[7] += 0.7071067811865475*Ghat[6]*dv11; 
  out[8] += 0.7071067811865475*Ghat[7]*dv11; 
  out[9] += -1.224744871391589*Ghat[1]*dv11; 
  out[10] += -1.224744871391589*Ghat[2]*dv11; 
  out[11] += -1.224744871391589*Ghat[3]*dv11; 
  out[12] += 0.7071067811865475*Ghat[8]*dv11; 
  out[13] += 0.7071067811865475*Ghat[9]*dv11; 
  out[14] += 0.7071067811865475*Ghat[10]*dv11; 
  out[15] += -1.224744871391589*Ghat[4]*dv11; 
  out[16] += 0.7071067811865475*Ghat[11]*dv11; 
  out[17] += -1.224744871391589*Ghat[5]*dv11; 
  out[18] += -1.224744871391589*Ghat[6]*dv11; 
  out[19] += -1.224744871391589*Ghat[7]*dv11; 
  out[20] += 0.7071067811865475*Ghat[12]*dv11; 
  out[21] += 0.7071067811865475*Ghat[13]*dv11; 
  out[22] += 0.7071067811865475*Ghat[14]*dv11; 
  out[23] += -1.224744871391589*Ghat[8]*dv11; 
  out[24] += -1.224744871391589*Ghat[9]*dv11; 
  out[25] += -1.224744871391589*Ghat[10]*dv11; 
  out[26] += -1.224744871391589*Ghat[11]*dv11; 
  out[27] += 0.7071067811865475*Ghat[15]*dv11; 
  out[28] += -1.224744871391589*Ghat[12]*dv11; 
  out[29] += -1.224744871391589*Ghat[13]*dv11; 
  out[30] += -1.224744871391589*Ghat[14]*dv11; 
  out[31] += -1.224744871391589*Ghat[15]*dv11; 

  } 
} 
