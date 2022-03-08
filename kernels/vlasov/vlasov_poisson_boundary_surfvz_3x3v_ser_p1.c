#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_6x_p1_surfx6_quad.h> 
#include <gkyl_basis_ser_6x_p1_upwind.h> 
GKYL_CU_DH void vlasov_poisson_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv12 = 2/dxv[5]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dx12 = 2/dxv[2]; 
  double alpha[32] = {0.0}; 

  alpha[0] = -3.464101615137754*phi[3]*dx12; 
  alpha[1] = -3.464101615137754*phi[5]*dx12; 
  alpha[2] = -3.464101615137754*phi[6]*dx12; 
  alpha[6] = -3.464101615137754*phi[7]*dx12; 

  double fUpwindQuad[32] = {0.0};
  double fUpwind[32] = {0.0};
  double Ghat[32] = {0.0}; 

  if (edge == -1) { 

  if (alpha[6]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx6_quad_0_r(fSkin); 
    fUpwindQuad[8] = ser_6x_p1_surfx6_quad_8_r(fSkin); 
    fUpwindQuad[16] = ser_6x_p1_surfx6_quad_16_r(fSkin); 
    fUpwindQuad[24] = ser_6x_p1_surfx6_quad_24_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx6_quad_0_l(fEdge); 
    fUpwindQuad[8] = ser_6x_p1_surfx6_quad_8_l(fEdge); 
    fUpwindQuad[16] = ser_6x_p1_surfx6_quad_16_l(fEdge); 
    fUpwindQuad[24] = ser_6x_p1_surfx6_quad_24_l(fEdge); 
  } 
  if ((-alpha[6])-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_6x_p1_surfx6_quad_1_r(fSkin); 
    fUpwindQuad[9] = ser_6x_p1_surfx6_quad_9_r(fSkin); 
    fUpwindQuad[17] = ser_6x_p1_surfx6_quad_17_r(fSkin); 
    fUpwindQuad[25] = ser_6x_p1_surfx6_quad_25_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_6x_p1_surfx6_quad_1_l(fEdge); 
    fUpwindQuad[9] = ser_6x_p1_surfx6_quad_9_l(fEdge); 
    fUpwindQuad[17] = ser_6x_p1_surfx6_quad_17_l(fEdge); 
    fUpwindQuad[25] = ser_6x_p1_surfx6_quad_25_l(fEdge); 
  } 
  if ((-alpha[6])+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_6x_p1_surfx6_quad_2_r(fSkin); 
    fUpwindQuad[10] = ser_6x_p1_surfx6_quad_10_r(fSkin); 
    fUpwindQuad[18] = ser_6x_p1_surfx6_quad_18_r(fSkin); 
    fUpwindQuad[26] = ser_6x_p1_surfx6_quad_26_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_6x_p1_surfx6_quad_2_l(fEdge); 
    fUpwindQuad[10] = ser_6x_p1_surfx6_quad_10_l(fEdge); 
    fUpwindQuad[18] = ser_6x_p1_surfx6_quad_18_l(fEdge); 
    fUpwindQuad[26] = ser_6x_p1_surfx6_quad_26_l(fEdge); 
  } 
  if (alpha[6]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_6x_p1_surfx6_quad_3_r(fSkin); 
    fUpwindQuad[11] = ser_6x_p1_surfx6_quad_11_r(fSkin); 
    fUpwindQuad[19] = ser_6x_p1_surfx6_quad_19_r(fSkin); 
    fUpwindQuad[27] = ser_6x_p1_surfx6_quad_27_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_6x_p1_surfx6_quad_3_l(fEdge); 
    fUpwindQuad[11] = ser_6x_p1_surfx6_quad_11_l(fEdge); 
    fUpwindQuad[19] = ser_6x_p1_surfx6_quad_19_l(fEdge); 
    fUpwindQuad[27] = ser_6x_p1_surfx6_quad_27_l(fEdge); 
  } 
  if (alpha[6]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx6_quad_4_r(fSkin); 
    fUpwindQuad[12] = ser_6x_p1_surfx6_quad_12_r(fSkin); 
    fUpwindQuad[20] = ser_6x_p1_surfx6_quad_20_r(fSkin); 
    fUpwindQuad[28] = ser_6x_p1_surfx6_quad_28_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx6_quad_4_l(fEdge); 
    fUpwindQuad[12] = ser_6x_p1_surfx6_quad_12_l(fEdge); 
    fUpwindQuad[20] = ser_6x_p1_surfx6_quad_20_l(fEdge); 
    fUpwindQuad[28] = ser_6x_p1_surfx6_quad_28_l(fEdge); 
  } 
  if ((-alpha[6])-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_6x_p1_surfx6_quad_5_r(fSkin); 
    fUpwindQuad[13] = ser_6x_p1_surfx6_quad_13_r(fSkin); 
    fUpwindQuad[21] = ser_6x_p1_surfx6_quad_21_r(fSkin); 
    fUpwindQuad[29] = ser_6x_p1_surfx6_quad_29_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_6x_p1_surfx6_quad_5_l(fEdge); 
    fUpwindQuad[13] = ser_6x_p1_surfx6_quad_13_l(fEdge); 
    fUpwindQuad[21] = ser_6x_p1_surfx6_quad_21_l(fEdge); 
    fUpwindQuad[29] = ser_6x_p1_surfx6_quad_29_l(fEdge); 
  } 
  if ((-alpha[6])+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_6x_p1_surfx6_quad_6_r(fSkin); 
    fUpwindQuad[14] = ser_6x_p1_surfx6_quad_14_r(fSkin); 
    fUpwindQuad[22] = ser_6x_p1_surfx6_quad_22_r(fSkin); 
    fUpwindQuad[30] = ser_6x_p1_surfx6_quad_30_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_6x_p1_surfx6_quad_6_l(fEdge); 
    fUpwindQuad[14] = ser_6x_p1_surfx6_quad_14_l(fEdge); 
    fUpwindQuad[22] = ser_6x_p1_surfx6_quad_22_l(fEdge); 
    fUpwindQuad[30] = ser_6x_p1_surfx6_quad_30_l(fEdge); 
  } 
  if (alpha[6]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_6x_p1_surfx6_quad_7_r(fSkin); 
    fUpwindQuad[15] = ser_6x_p1_surfx6_quad_15_r(fSkin); 
    fUpwindQuad[23] = ser_6x_p1_surfx6_quad_23_r(fSkin); 
    fUpwindQuad[31] = ser_6x_p1_surfx6_quad_31_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_6x_p1_surfx6_quad_7_l(fEdge); 
    fUpwindQuad[15] = ser_6x_p1_surfx6_quad_15_l(fEdge); 
    fUpwindQuad[23] = ser_6x_p1_surfx6_quad_23_l(fEdge); 
    fUpwindQuad[31] = ser_6x_p1_surfx6_quad_31_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_6x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.1767766952966368*(alpha[6]*fUpwind[6]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.1767766952966368*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.1767766952966368*(alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.1767766952966368*(alpha[6]*fUpwind[16]+alpha[2]*fUpwind[8]+alpha[1]*fUpwind[7]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.1767766952966368*(alpha[6]*fUpwind[17]+alpha[2]*fUpwind[10]+alpha[1]*fUpwind[9]+alpha[0]*fUpwind[4]); 
  Ghat[5] += 0.1767766952966368*(alpha[6]*fUpwind[20]+alpha[2]*fUpwind[13]+alpha[1]*fUpwind[12]+alpha[0]*fUpwind[5]); 
  Ghat[6] += 0.1767766952966368*(alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[7] += 0.1767766952966368*(alpha[2]*fUpwind[16]+alpha[6]*fUpwind[8]+alpha[0]*fUpwind[7]+alpha[1]*fUpwind[3]); 
  Ghat[8] += 0.1767766952966368*(alpha[1]*fUpwind[16]+alpha[0]*fUpwind[8]+alpha[6]*fUpwind[7]+alpha[2]*fUpwind[3]); 
  Ghat[9] += 0.1767766952966368*(alpha[2]*fUpwind[17]+alpha[6]*fUpwind[10]+alpha[0]*fUpwind[9]+alpha[1]*fUpwind[4]); 
  Ghat[10] += 0.1767766952966368*(alpha[1]*fUpwind[17]+alpha[0]*fUpwind[10]+alpha[6]*fUpwind[9]+alpha[2]*fUpwind[4]); 
  Ghat[11] += 0.1767766952966368*(alpha[6]*fUpwind[26]+alpha[2]*fUpwind[19]+alpha[1]*fUpwind[18]+alpha[0]*fUpwind[11]); 
  Ghat[12] += 0.1767766952966368*(alpha[2]*fUpwind[20]+alpha[6]*fUpwind[13]+alpha[0]*fUpwind[12]+alpha[1]*fUpwind[5]); 
  Ghat[13] += 0.1767766952966368*(alpha[1]*fUpwind[20]+alpha[0]*fUpwind[13]+alpha[6]*fUpwind[12]+alpha[2]*fUpwind[5]); 
  Ghat[14] += 0.1767766952966368*(alpha[6]*fUpwind[27]+alpha[2]*fUpwind[22]+alpha[1]*fUpwind[21]+alpha[0]*fUpwind[14]); 
  Ghat[15] += 0.1767766952966368*(alpha[6]*fUpwind[28]+alpha[2]*fUpwind[24]+alpha[1]*fUpwind[23]+alpha[0]*fUpwind[15]); 
  Ghat[16] += 0.1767766952966368*(alpha[0]*fUpwind[16]+alpha[1]*fUpwind[8]+alpha[2]*fUpwind[7]+fUpwind[3]*alpha[6]); 
  Ghat[17] += 0.1767766952966368*(alpha[0]*fUpwind[17]+alpha[1]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[4]*alpha[6]); 
  Ghat[18] += 0.1767766952966368*(alpha[2]*fUpwind[26]+alpha[6]*fUpwind[19]+alpha[0]*fUpwind[18]+alpha[1]*fUpwind[11]); 
  Ghat[19] += 0.1767766952966368*(alpha[1]*fUpwind[26]+alpha[0]*fUpwind[19]+alpha[6]*fUpwind[18]+alpha[2]*fUpwind[11]); 
  Ghat[20] += 0.1767766952966368*(alpha[0]*fUpwind[20]+alpha[1]*fUpwind[13]+alpha[2]*fUpwind[12]+fUpwind[5]*alpha[6]); 
  Ghat[21] += 0.1767766952966368*(alpha[2]*fUpwind[27]+alpha[6]*fUpwind[22]+alpha[0]*fUpwind[21]+alpha[1]*fUpwind[14]); 
  Ghat[22] += 0.1767766952966368*(alpha[1]*fUpwind[27]+alpha[0]*fUpwind[22]+alpha[6]*fUpwind[21]+alpha[2]*fUpwind[14]); 
  Ghat[23] += 0.1767766952966368*(alpha[2]*fUpwind[28]+alpha[6]*fUpwind[24]+alpha[0]*fUpwind[23]+alpha[1]*fUpwind[15]); 
  Ghat[24] += 0.1767766952966368*(alpha[1]*fUpwind[28]+alpha[0]*fUpwind[24]+alpha[6]*fUpwind[23]+alpha[2]*fUpwind[15]); 
  Ghat[25] += 0.1767766952966368*(alpha[6]*fUpwind[31]+alpha[2]*fUpwind[30]+alpha[1]*fUpwind[29]+alpha[0]*fUpwind[25]); 
  Ghat[26] += 0.1767766952966368*(alpha[0]*fUpwind[26]+alpha[1]*fUpwind[19]+alpha[2]*fUpwind[18]+alpha[6]*fUpwind[11]); 
  Ghat[27] += 0.1767766952966368*(alpha[0]*fUpwind[27]+alpha[1]*fUpwind[22]+alpha[2]*fUpwind[21]+alpha[6]*fUpwind[14]); 
  Ghat[28] += 0.1767766952966368*(alpha[0]*fUpwind[28]+alpha[1]*fUpwind[24]+alpha[2]*fUpwind[23]+alpha[6]*fUpwind[15]); 
  Ghat[29] += 0.1767766952966368*(alpha[2]*fUpwind[31]+alpha[6]*fUpwind[30]+alpha[0]*fUpwind[29]+alpha[1]*fUpwind[25]); 
  Ghat[30] += 0.1767766952966368*(alpha[1]*fUpwind[31]+alpha[0]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[2]*fUpwind[25]); 
  Ghat[31] += 0.1767766952966368*(alpha[0]*fUpwind[31]+alpha[1]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[6]*fUpwind[25]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv12; 
  out[1] += -0.7071067811865475*Ghat[1]*dv12; 
  out[2] += -0.7071067811865475*Ghat[2]*dv12; 
  out[3] += -0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -0.7071067811865475*Ghat[4]*dv12; 
  out[5] += -0.7071067811865475*Ghat[5]*dv12; 
  out[6] += -1.224744871391589*Ghat[0]*dv12; 
  out[7] += -0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -0.7071067811865475*Ghat[7]*dv12; 
  out[9] += -0.7071067811865475*Ghat[8]*dv12; 
  out[10] += -0.7071067811865475*Ghat[9]*dv12; 
  out[11] += -0.7071067811865475*Ghat[10]*dv12; 
  out[12] += -0.7071067811865475*Ghat[11]*dv12; 
  out[13] += -0.7071067811865475*Ghat[12]*dv12; 
  out[14] += -0.7071067811865475*Ghat[13]*dv12; 
  out[15] += -0.7071067811865475*Ghat[14]*dv12; 
  out[16] += -0.7071067811865475*Ghat[15]*dv12; 
  out[17] += -1.224744871391589*Ghat[1]*dv12; 
  out[18] += -1.224744871391589*Ghat[2]*dv12; 
  out[19] += -1.224744871391589*Ghat[3]*dv12; 
  out[20] += -1.224744871391589*Ghat[4]*dv12; 
  out[21] += -1.224744871391589*Ghat[5]*dv12; 
  out[22] += -0.7071067811865475*Ghat[16]*dv12; 
  out[23] += -0.7071067811865475*Ghat[17]*dv12; 
  out[24] += -0.7071067811865475*Ghat[18]*dv12; 
  out[25] += -0.7071067811865475*Ghat[19]*dv12; 
  out[26] += -0.7071067811865475*Ghat[20]*dv12; 
  out[27] += -0.7071067811865475*Ghat[21]*dv12; 
  out[28] += -0.7071067811865475*Ghat[22]*dv12; 
  out[29] += -0.7071067811865475*Ghat[23]*dv12; 
  out[30] += -0.7071067811865475*Ghat[24]*dv12; 
  out[31] += -0.7071067811865475*Ghat[25]*dv12; 
  out[32] += -1.224744871391589*Ghat[6]*dv12; 
  out[33] += -1.224744871391589*Ghat[7]*dv12; 
  out[34] += -1.224744871391589*Ghat[8]*dv12; 
  out[35] += -1.224744871391589*Ghat[9]*dv12; 
  out[36] += -1.224744871391589*Ghat[10]*dv12; 
  out[37] += -1.224744871391589*Ghat[11]*dv12; 
  out[38] += -1.224744871391589*Ghat[12]*dv12; 
  out[39] += -1.224744871391589*Ghat[13]*dv12; 
  out[40] += -1.224744871391589*Ghat[14]*dv12; 
  out[41] += -1.224744871391589*Ghat[15]*dv12; 
  out[42] += -0.7071067811865475*Ghat[26]*dv12; 
  out[43] += -0.7071067811865475*Ghat[27]*dv12; 
  out[44] += -0.7071067811865475*Ghat[28]*dv12; 
  out[45] += -0.7071067811865475*Ghat[29]*dv12; 
  out[46] += -0.7071067811865475*Ghat[30]*dv12; 
  out[47] += -1.224744871391589*Ghat[16]*dv12; 
  out[48] += -1.224744871391589*Ghat[17]*dv12; 
  out[49] += -1.224744871391589*Ghat[18]*dv12; 
  out[50] += -1.224744871391589*Ghat[19]*dv12; 
  out[51] += -1.224744871391589*Ghat[20]*dv12; 
  out[52] += -1.224744871391589*Ghat[21]*dv12; 
  out[53] += -1.224744871391589*Ghat[22]*dv12; 
  out[54] += -1.224744871391589*Ghat[23]*dv12; 
  out[55] += -1.224744871391589*Ghat[24]*dv12; 
  out[56] += -1.224744871391589*Ghat[25]*dv12; 
  out[57] += -0.7071067811865475*Ghat[31]*dv12; 
  out[58] += -1.224744871391589*Ghat[26]*dv12; 
  out[59] += -1.224744871391589*Ghat[27]*dv12; 
  out[60] += -1.224744871391589*Ghat[28]*dv12; 
  out[61] += -1.224744871391589*Ghat[29]*dv12; 
  out[62] += -1.224744871391589*Ghat[30]*dv12; 
  out[63] += -1.224744871391589*Ghat[31]*dv12; 

  } else { 

  if (alpha[6]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx6_quad_0_r(fEdge); 
    fUpwindQuad[8] = ser_6x_p1_surfx6_quad_8_r(fEdge); 
    fUpwindQuad[16] = ser_6x_p1_surfx6_quad_16_r(fEdge); 
    fUpwindQuad[24] = ser_6x_p1_surfx6_quad_24_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx6_quad_0_l(fSkin); 
    fUpwindQuad[8] = ser_6x_p1_surfx6_quad_8_l(fSkin); 
    fUpwindQuad[16] = ser_6x_p1_surfx6_quad_16_l(fSkin); 
    fUpwindQuad[24] = ser_6x_p1_surfx6_quad_24_l(fSkin); 
  } 
  if ((-alpha[6])-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_6x_p1_surfx6_quad_1_r(fEdge); 
    fUpwindQuad[9] = ser_6x_p1_surfx6_quad_9_r(fEdge); 
    fUpwindQuad[17] = ser_6x_p1_surfx6_quad_17_r(fEdge); 
    fUpwindQuad[25] = ser_6x_p1_surfx6_quad_25_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_6x_p1_surfx6_quad_1_l(fSkin); 
    fUpwindQuad[9] = ser_6x_p1_surfx6_quad_9_l(fSkin); 
    fUpwindQuad[17] = ser_6x_p1_surfx6_quad_17_l(fSkin); 
    fUpwindQuad[25] = ser_6x_p1_surfx6_quad_25_l(fSkin); 
  } 
  if ((-alpha[6])+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_6x_p1_surfx6_quad_2_r(fEdge); 
    fUpwindQuad[10] = ser_6x_p1_surfx6_quad_10_r(fEdge); 
    fUpwindQuad[18] = ser_6x_p1_surfx6_quad_18_r(fEdge); 
    fUpwindQuad[26] = ser_6x_p1_surfx6_quad_26_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_6x_p1_surfx6_quad_2_l(fSkin); 
    fUpwindQuad[10] = ser_6x_p1_surfx6_quad_10_l(fSkin); 
    fUpwindQuad[18] = ser_6x_p1_surfx6_quad_18_l(fSkin); 
    fUpwindQuad[26] = ser_6x_p1_surfx6_quad_26_l(fSkin); 
  } 
  if (alpha[6]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_6x_p1_surfx6_quad_3_r(fEdge); 
    fUpwindQuad[11] = ser_6x_p1_surfx6_quad_11_r(fEdge); 
    fUpwindQuad[19] = ser_6x_p1_surfx6_quad_19_r(fEdge); 
    fUpwindQuad[27] = ser_6x_p1_surfx6_quad_27_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_6x_p1_surfx6_quad_3_l(fSkin); 
    fUpwindQuad[11] = ser_6x_p1_surfx6_quad_11_l(fSkin); 
    fUpwindQuad[19] = ser_6x_p1_surfx6_quad_19_l(fSkin); 
    fUpwindQuad[27] = ser_6x_p1_surfx6_quad_27_l(fSkin); 
  } 
  if (alpha[6]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx6_quad_4_r(fEdge); 
    fUpwindQuad[12] = ser_6x_p1_surfx6_quad_12_r(fEdge); 
    fUpwindQuad[20] = ser_6x_p1_surfx6_quad_20_r(fEdge); 
    fUpwindQuad[28] = ser_6x_p1_surfx6_quad_28_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx6_quad_4_l(fSkin); 
    fUpwindQuad[12] = ser_6x_p1_surfx6_quad_12_l(fSkin); 
    fUpwindQuad[20] = ser_6x_p1_surfx6_quad_20_l(fSkin); 
    fUpwindQuad[28] = ser_6x_p1_surfx6_quad_28_l(fSkin); 
  } 
  if ((-alpha[6])-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_6x_p1_surfx6_quad_5_r(fEdge); 
    fUpwindQuad[13] = ser_6x_p1_surfx6_quad_13_r(fEdge); 
    fUpwindQuad[21] = ser_6x_p1_surfx6_quad_21_r(fEdge); 
    fUpwindQuad[29] = ser_6x_p1_surfx6_quad_29_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_6x_p1_surfx6_quad_5_l(fSkin); 
    fUpwindQuad[13] = ser_6x_p1_surfx6_quad_13_l(fSkin); 
    fUpwindQuad[21] = ser_6x_p1_surfx6_quad_21_l(fSkin); 
    fUpwindQuad[29] = ser_6x_p1_surfx6_quad_29_l(fSkin); 
  } 
  if ((-alpha[6])+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_6x_p1_surfx6_quad_6_r(fEdge); 
    fUpwindQuad[14] = ser_6x_p1_surfx6_quad_14_r(fEdge); 
    fUpwindQuad[22] = ser_6x_p1_surfx6_quad_22_r(fEdge); 
    fUpwindQuad[30] = ser_6x_p1_surfx6_quad_30_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_6x_p1_surfx6_quad_6_l(fSkin); 
    fUpwindQuad[14] = ser_6x_p1_surfx6_quad_14_l(fSkin); 
    fUpwindQuad[22] = ser_6x_p1_surfx6_quad_22_l(fSkin); 
    fUpwindQuad[30] = ser_6x_p1_surfx6_quad_30_l(fSkin); 
  } 
  if (alpha[6]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_6x_p1_surfx6_quad_7_r(fEdge); 
    fUpwindQuad[15] = ser_6x_p1_surfx6_quad_15_r(fEdge); 
    fUpwindQuad[23] = ser_6x_p1_surfx6_quad_23_r(fEdge); 
    fUpwindQuad[31] = ser_6x_p1_surfx6_quad_31_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_6x_p1_surfx6_quad_7_l(fSkin); 
    fUpwindQuad[15] = ser_6x_p1_surfx6_quad_15_l(fSkin); 
    fUpwindQuad[23] = ser_6x_p1_surfx6_quad_23_l(fSkin); 
    fUpwindQuad[31] = ser_6x_p1_surfx6_quad_31_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_6x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.1767766952966368*(alpha[6]*fUpwind[6]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.1767766952966368*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.1767766952966368*(alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.1767766952966368*(alpha[6]*fUpwind[16]+alpha[2]*fUpwind[8]+alpha[1]*fUpwind[7]+alpha[0]*fUpwind[3]); 
  Ghat[4] += 0.1767766952966368*(alpha[6]*fUpwind[17]+alpha[2]*fUpwind[10]+alpha[1]*fUpwind[9]+alpha[0]*fUpwind[4]); 
  Ghat[5] += 0.1767766952966368*(alpha[6]*fUpwind[20]+alpha[2]*fUpwind[13]+alpha[1]*fUpwind[12]+alpha[0]*fUpwind[5]); 
  Ghat[6] += 0.1767766952966368*(alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[7] += 0.1767766952966368*(alpha[2]*fUpwind[16]+alpha[6]*fUpwind[8]+alpha[0]*fUpwind[7]+alpha[1]*fUpwind[3]); 
  Ghat[8] += 0.1767766952966368*(alpha[1]*fUpwind[16]+alpha[0]*fUpwind[8]+alpha[6]*fUpwind[7]+alpha[2]*fUpwind[3]); 
  Ghat[9] += 0.1767766952966368*(alpha[2]*fUpwind[17]+alpha[6]*fUpwind[10]+alpha[0]*fUpwind[9]+alpha[1]*fUpwind[4]); 
  Ghat[10] += 0.1767766952966368*(alpha[1]*fUpwind[17]+alpha[0]*fUpwind[10]+alpha[6]*fUpwind[9]+alpha[2]*fUpwind[4]); 
  Ghat[11] += 0.1767766952966368*(alpha[6]*fUpwind[26]+alpha[2]*fUpwind[19]+alpha[1]*fUpwind[18]+alpha[0]*fUpwind[11]); 
  Ghat[12] += 0.1767766952966368*(alpha[2]*fUpwind[20]+alpha[6]*fUpwind[13]+alpha[0]*fUpwind[12]+alpha[1]*fUpwind[5]); 
  Ghat[13] += 0.1767766952966368*(alpha[1]*fUpwind[20]+alpha[0]*fUpwind[13]+alpha[6]*fUpwind[12]+alpha[2]*fUpwind[5]); 
  Ghat[14] += 0.1767766952966368*(alpha[6]*fUpwind[27]+alpha[2]*fUpwind[22]+alpha[1]*fUpwind[21]+alpha[0]*fUpwind[14]); 
  Ghat[15] += 0.1767766952966368*(alpha[6]*fUpwind[28]+alpha[2]*fUpwind[24]+alpha[1]*fUpwind[23]+alpha[0]*fUpwind[15]); 
  Ghat[16] += 0.1767766952966368*(alpha[0]*fUpwind[16]+alpha[1]*fUpwind[8]+alpha[2]*fUpwind[7]+fUpwind[3]*alpha[6]); 
  Ghat[17] += 0.1767766952966368*(alpha[0]*fUpwind[17]+alpha[1]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[4]*alpha[6]); 
  Ghat[18] += 0.1767766952966368*(alpha[2]*fUpwind[26]+alpha[6]*fUpwind[19]+alpha[0]*fUpwind[18]+alpha[1]*fUpwind[11]); 
  Ghat[19] += 0.1767766952966368*(alpha[1]*fUpwind[26]+alpha[0]*fUpwind[19]+alpha[6]*fUpwind[18]+alpha[2]*fUpwind[11]); 
  Ghat[20] += 0.1767766952966368*(alpha[0]*fUpwind[20]+alpha[1]*fUpwind[13]+alpha[2]*fUpwind[12]+fUpwind[5]*alpha[6]); 
  Ghat[21] += 0.1767766952966368*(alpha[2]*fUpwind[27]+alpha[6]*fUpwind[22]+alpha[0]*fUpwind[21]+alpha[1]*fUpwind[14]); 
  Ghat[22] += 0.1767766952966368*(alpha[1]*fUpwind[27]+alpha[0]*fUpwind[22]+alpha[6]*fUpwind[21]+alpha[2]*fUpwind[14]); 
  Ghat[23] += 0.1767766952966368*(alpha[2]*fUpwind[28]+alpha[6]*fUpwind[24]+alpha[0]*fUpwind[23]+alpha[1]*fUpwind[15]); 
  Ghat[24] += 0.1767766952966368*(alpha[1]*fUpwind[28]+alpha[0]*fUpwind[24]+alpha[6]*fUpwind[23]+alpha[2]*fUpwind[15]); 
  Ghat[25] += 0.1767766952966368*(alpha[6]*fUpwind[31]+alpha[2]*fUpwind[30]+alpha[1]*fUpwind[29]+alpha[0]*fUpwind[25]); 
  Ghat[26] += 0.1767766952966368*(alpha[0]*fUpwind[26]+alpha[1]*fUpwind[19]+alpha[2]*fUpwind[18]+alpha[6]*fUpwind[11]); 
  Ghat[27] += 0.1767766952966368*(alpha[0]*fUpwind[27]+alpha[1]*fUpwind[22]+alpha[2]*fUpwind[21]+alpha[6]*fUpwind[14]); 
  Ghat[28] += 0.1767766952966368*(alpha[0]*fUpwind[28]+alpha[1]*fUpwind[24]+alpha[2]*fUpwind[23]+alpha[6]*fUpwind[15]); 
  Ghat[29] += 0.1767766952966368*(alpha[2]*fUpwind[31]+alpha[6]*fUpwind[30]+alpha[0]*fUpwind[29]+alpha[1]*fUpwind[25]); 
  Ghat[30] += 0.1767766952966368*(alpha[1]*fUpwind[31]+alpha[0]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[2]*fUpwind[25]); 
  Ghat[31] += 0.1767766952966368*(alpha[0]*fUpwind[31]+alpha[1]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[6]*fUpwind[25]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv12; 
  out[1] += 0.7071067811865475*Ghat[1]*dv12; 
  out[2] += 0.7071067811865475*Ghat[2]*dv12; 
  out[3] += 0.7071067811865475*Ghat[3]*dv12; 
  out[4] += 0.7071067811865475*Ghat[4]*dv12; 
  out[5] += 0.7071067811865475*Ghat[5]*dv12; 
  out[6] += -1.224744871391589*Ghat[0]*dv12; 
  out[7] += 0.7071067811865475*Ghat[6]*dv12; 
  out[8] += 0.7071067811865475*Ghat[7]*dv12; 
  out[9] += 0.7071067811865475*Ghat[8]*dv12; 
  out[10] += 0.7071067811865475*Ghat[9]*dv12; 
  out[11] += 0.7071067811865475*Ghat[10]*dv12; 
  out[12] += 0.7071067811865475*Ghat[11]*dv12; 
  out[13] += 0.7071067811865475*Ghat[12]*dv12; 
  out[14] += 0.7071067811865475*Ghat[13]*dv12; 
  out[15] += 0.7071067811865475*Ghat[14]*dv12; 
  out[16] += 0.7071067811865475*Ghat[15]*dv12; 
  out[17] += -1.224744871391589*Ghat[1]*dv12; 
  out[18] += -1.224744871391589*Ghat[2]*dv12; 
  out[19] += -1.224744871391589*Ghat[3]*dv12; 
  out[20] += -1.224744871391589*Ghat[4]*dv12; 
  out[21] += -1.224744871391589*Ghat[5]*dv12; 
  out[22] += 0.7071067811865475*Ghat[16]*dv12; 
  out[23] += 0.7071067811865475*Ghat[17]*dv12; 
  out[24] += 0.7071067811865475*Ghat[18]*dv12; 
  out[25] += 0.7071067811865475*Ghat[19]*dv12; 
  out[26] += 0.7071067811865475*Ghat[20]*dv12; 
  out[27] += 0.7071067811865475*Ghat[21]*dv12; 
  out[28] += 0.7071067811865475*Ghat[22]*dv12; 
  out[29] += 0.7071067811865475*Ghat[23]*dv12; 
  out[30] += 0.7071067811865475*Ghat[24]*dv12; 
  out[31] += 0.7071067811865475*Ghat[25]*dv12; 
  out[32] += -1.224744871391589*Ghat[6]*dv12; 
  out[33] += -1.224744871391589*Ghat[7]*dv12; 
  out[34] += -1.224744871391589*Ghat[8]*dv12; 
  out[35] += -1.224744871391589*Ghat[9]*dv12; 
  out[36] += -1.224744871391589*Ghat[10]*dv12; 
  out[37] += -1.224744871391589*Ghat[11]*dv12; 
  out[38] += -1.224744871391589*Ghat[12]*dv12; 
  out[39] += -1.224744871391589*Ghat[13]*dv12; 
  out[40] += -1.224744871391589*Ghat[14]*dv12; 
  out[41] += -1.224744871391589*Ghat[15]*dv12; 
  out[42] += 0.7071067811865475*Ghat[26]*dv12; 
  out[43] += 0.7071067811865475*Ghat[27]*dv12; 
  out[44] += 0.7071067811865475*Ghat[28]*dv12; 
  out[45] += 0.7071067811865475*Ghat[29]*dv12; 
  out[46] += 0.7071067811865475*Ghat[30]*dv12; 
  out[47] += -1.224744871391589*Ghat[16]*dv12; 
  out[48] += -1.224744871391589*Ghat[17]*dv12; 
  out[49] += -1.224744871391589*Ghat[18]*dv12; 
  out[50] += -1.224744871391589*Ghat[19]*dv12; 
  out[51] += -1.224744871391589*Ghat[20]*dv12; 
  out[52] += -1.224744871391589*Ghat[21]*dv12; 
  out[53] += -1.224744871391589*Ghat[22]*dv12; 
  out[54] += -1.224744871391589*Ghat[23]*dv12; 
  out[55] += -1.224744871391589*Ghat[24]*dv12; 
  out[56] += -1.224744871391589*Ghat[25]*dv12; 
  out[57] += 0.7071067811865475*Ghat[31]*dv12; 
  out[58] += -1.224744871391589*Ghat[26]*dv12; 
  out[59] += -1.224744871391589*Ghat[27]*dv12; 
  out[60] += -1.224744871391589*Ghat[28]*dv12; 
  out[61] += -1.224744871391589*Ghat[29]*dv12; 
  out[62] += -1.224744871391589*Ghat[30]*dv12; 
  out[63] += -1.224744871391589*Ghat[31]*dv12; 

  } 
} 
