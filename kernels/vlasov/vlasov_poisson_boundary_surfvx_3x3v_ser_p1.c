#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_6x_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_6x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_3x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[3]; 
  const double dv1 = dxv[3], wv1 = w[3]; 
  const double dv2 = dxv[4], wv2 = w[4]; 
  const double dv3 = dxv[5], wv3 = w[5]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double dx12 = 2/dxv[2]; 
  double alpha[32] = {0.0}; 

  alpha[0] = -3.464101615137754*phi[1]*dx10; 
  alpha[2] = -3.464101615137754*phi[4]*dx10; 
  alpha[3] = -3.464101615137754*phi[5]*dx10; 
  alpha[8] = -3.464101615137754*phi[7]*dx10; 

  double fUpwindQuad[32] = {0.0};
  double fUpwind[32] = {0.0};
  double Ghat[32] = {0.0}; 

  if (edge == -1) { 

  if (alpha[8]-alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx4_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_6x_p1_surfx4_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = ser_6x_p1_surfx4_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = ser_6x_p1_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx4_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_6x_p1_surfx4_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = ser_6x_p1_surfx4_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = ser_6x_p1_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if (alpha[8]-alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx4_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = ser_6x_p1_surfx4_eval_quad_node_5_r(fSkin); 
    fUpwindQuad[6] = ser_6x_p1_surfx4_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = ser_6x_p1_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx4_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = ser_6x_p1_surfx4_eval_quad_node_5_l(fEdge); 
    fUpwindQuad[6] = ser_6x_p1_surfx4_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = ser_6x_p1_surfx4_eval_quad_node_7_l(fEdge); 
  } 
  if (alpha[8]-alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_6x_p1_surfx4_eval_quad_node_8_r(fSkin); 
    fUpwindQuad[9] = ser_6x_p1_surfx4_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = ser_6x_p1_surfx4_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = ser_6x_p1_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_6x_p1_surfx4_eval_quad_node_8_l(fEdge); 
    fUpwindQuad[9] = ser_6x_p1_surfx4_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = ser_6x_p1_surfx4_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = ser_6x_p1_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if (alpha[8]-alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_6x_p1_surfx4_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = ser_6x_p1_surfx4_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = ser_6x_p1_surfx4_eval_quad_node_14_r(fSkin); 
    fUpwindQuad[15] = ser_6x_p1_surfx4_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_6x_p1_surfx4_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = ser_6x_p1_surfx4_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = ser_6x_p1_surfx4_eval_quad_node_14_l(fEdge); 
    fUpwindQuad[15] = ser_6x_p1_surfx4_eval_quad_node_15_l(fEdge); 
  } 
  if ((-alpha[8])+alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[16] = ser_6x_p1_surfx4_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = ser_6x_p1_surfx4_eval_quad_node_17_r(fSkin); 
    fUpwindQuad[18] = ser_6x_p1_surfx4_eval_quad_node_18_r(fSkin); 
    fUpwindQuad[19] = ser_6x_p1_surfx4_eval_quad_node_19_r(fSkin); 
  } else { 
    fUpwindQuad[16] = ser_6x_p1_surfx4_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = ser_6x_p1_surfx4_eval_quad_node_17_l(fEdge); 
    fUpwindQuad[18] = ser_6x_p1_surfx4_eval_quad_node_18_l(fEdge); 
    fUpwindQuad[19] = ser_6x_p1_surfx4_eval_quad_node_19_l(fEdge); 
  } 
  if ((-alpha[8])+alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[20] = ser_6x_p1_surfx4_eval_quad_node_20_r(fSkin); 
    fUpwindQuad[21] = ser_6x_p1_surfx4_eval_quad_node_21_r(fSkin); 
    fUpwindQuad[22] = ser_6x_p1_surfx4_eval_quad_node_22_r(fSkin); 
    fUpwindQuad[23] = ser_6x_p1_surfx4_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[20] = ser_6x_p1_surfx4_eval_quad_node_20_l(fEdge); 
    fUpwindQuad[21] = ser_6x_p1_surfx4_eval_quad_node_21_l(fEdge); 
    fUpwindQuad[22] = ser_6x_p1_surfx4_eval_quad_node_22_l(fEdge); 
    fUpwindQuad[23] = ser_6x_p1_surfx4_eval_quad_node_23_l(fEdge); 
  } 
  if ((-alpha[8])+alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[24] = ser_6x_p1_surfx4_eval_quad_node_24_r(fSkin); 
    fUpwindQuad[25] = ser_6x_p1_surfx4_eval_quad_node_25_r(fSkin); 
    fUpwindQuad[26] = ser_6x_p1_surfx4_eval_quad_node_26_r(fSkin); 
    fUpwindQuad[27] = ser_6x_p1_surfx4_eval_quad_node_27_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_6x_p1_surfx4_eval_quad_node_24_l(fEdge); 
    fUpwindQuad[25] = ser_6x_p1_surfx4_eval_quad_node_25_l(fEdge); 
    fUpwindQuad[26] = ser_6x_p1_surfx4_eval_quad_node_26_l(fEdge); 
    fUpwindQuad[27] = ser_6x_p1_surfx4_eval_quad_node_27_l(fEdge); 
  } 
  if ((-alpha[8])+alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[28] = ser_6x_p1_surfx4_eval_quad_node_28_r(fSkin); 
    fUpwindQuad[29] = ser_6x_p1_surfx4_eval_quad_node_29_r(fSkin); 
    fUpwindQuad[30] = ser_6x_p1_surfx4_eval_quad_node_30_r(fSkin); 
    fUpwindQuad[31] = ser_6x_p1_surfx4_eval_quad_node_31_r(fSkin); 
  } else { 
    fUpwindQuad[28] = ser_6x_p1_surfx4_eval_quad_node_28_l(fEdge); 
    fUpwindQuad[29] = ser_6x_p1_surfx4_eval_quad_node_29_l(fEdge); 
    fUpwindQuad[30] = ser_6x_p1_surfx4_eval_quad_node_30_l(fEdge); 
    fUpwindQuad[31] = ser_6x_p1_surfx4_eval_quad_node_31_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*(alpha[8]*fUpwind[8]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1767766952966368*(alpha[8]*fUpwind[16]+alpha[3]*fUpwind[7]+alpha[2]*fUpwind[6]+alpha[0]*fUpwind[1]); 
  Ghat[2] = 0.1767766952966368*(alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.1767766952966368*(alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.1767766952966368*(alpha[8]*fUpwind[19]+alpha[3]*fUpwind[11]+alpha[2]*fUpwind[10]+alpha[0]*fUpwind[4]); 
  Ghat[5] = 0.1767766952966368*(alpha[8]*fUpwind[22]+alpha[3]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[0]*fUpwind[5]); 
  Ghat[6] = 0.1767766952966368*(alpha[3]*fUpwind[16]+fUpwind[7]*alpha[8]+alpha[0]*fUpwind[6]+fUpwind[1]*alpha[2]); 
  Ghat[7] = 0.1767766952966368*(alpha[2]*fUpwind[16]+fUpwind[6]*alpha[8]+alpha[0]*fUpwind[7]+fUpwind[1]*alpha[3]); 
  Ghat[8] = 0.1767766952966368*(alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[9] = 0.1767766952966368*(alpha[8]*fUpwind[26]+alpha[3]*fUpwind[18]+alpha[2]*fUpwind[17]+alpha[0]*fUpwind[9]); 
  Ghat[10] = 0.1767766952966368*(alpha[3]*fUpwind[19]+alpha[8]*fUpwind[11]+alpha[0]*fUpwind[10]+alpha[2]*fUpwind[4]); 
  Ghat[11] = 0.1767766952966368*(alpha[2]*fUpwind[19]+alpha[0]*fUpwind[11]+alpha[8]*fUpwind[10]+alpha[3]*fUpwind[4]); 
  Ghat[12] = 0.1767766952966368*(alpha[8]*fUpwind[27]+alpha[3]*fUpwind[21]+alpha[2]*fUpwind[20]+alpha[0]*fUpwind[12]); 
  Ghat[13] = 0.1767766952966368*(alpha[3]*fUpwind[22]+alpha[8]*fUpwind[14]+alpha[0]*fUpwind[13]+alpha[2]*fUpwind[5]); 
  Ghat[14] = 0.1767766952966368*(alpha[2]*fUpwind[22]+alpha[0]*fUpwind[14]+alpha[8]*fUpwind[13]+alpha[3]*fUpwind[5]); 
  Ghat[15] = 0.1767766952966368*(alpha[8]*fUpwind[30]+alpha[3]*fUpwind[25]+alpha[2]*fUpwind[24]+alpha[0]*fUpwind[15]); 
  Ghat[16] = 0.1767766952966368*(alpha[0]*fUpwind[16]+fUpwind[1]*alpha[8]+alpha[2]*fUpwind[7]+alpha[3]*fUpwind[6]); 
  Ghat[17] = 0.1767766952966368*(alpha[3]*fUpwind[26]+alpha[8]*fUpwind[18]+alpha[0]*fUpwind[17]+alpha[2]*fUpwind[9]); 
  Ghat[18] = 0.1767766952966368*(alpha[2]*fUpwind[26]+alpha[0]*fUpwind[18]+alpha[8]*fUpwind[17]+alpha[3]*fUpwind[9]); 
  Ghat[19] = 0.1767766952966368*(alpha[0]*fUpwind[19]+alpha[2]*fUpwind[11]+alpha[3]*fUpwind[10]+fUpwind[4]*alpha[8]); 
  Ghat[20] = 0.1767766952966368*(alpha[3]*fUpwind[27]+alpha[8]*fUpwind[21]+alpha[0]*fUpwind[20]+alpha[2]*fUpwind[12]); 
  Ghat[21] = 0.1767766952966368*(alpha[2]*fUpwind[27]+alpha[0]*fUpwind[21]+alpha[8]*fUpwind[20]+alpha[3]*fUpwind[12]); 
  Ghat[22] = 0.1767766952966368*(alpha[0]*fUpwind[22]+alpha[2]*fUpwind[14]+alpha[3]*fUpwind[13]+fUpwind[5]*alpha[8]); 
  Ghat[23] = 0.1767766952966368*(alpha[8]*fUpwind[31]+alpha[3]*fUpwind[29]+alpha[2]*fUpwind[28]+alpha[0]*fUpwind[23]); 
  Ghat[24] = 0.1767766952966368*(alpha[3]*fUpwind[30]+alpha[8]*fUpwind[25]+alpha[0]*fUpwind[24]+alpha[2]*fUpwind[15]); 
  Ghat[25] = 0.1767766952966368*(alpha[2]*fUpwind[30]+alpha[0]*fUpwind[25]+alpha[8]*fUpwind[24]+alpha[3]*fUpwind[15]); 
  Ghat[26] = 0.1767766952966368*(alpha[0]*fUpwind[26]+alpha[2]*fUpwind[18]+alpha[3]*fUpwind[17]+alpha[8]*fUpwind[9]); 
  Ghat[27] = 0.1767766952966368*(alpha[0]*fUpwind[27]+alpha[2]*fUpwind[21]+alpha[3]*fUpwind[20]+alpha[8]*fUpwind[12]); 
  Ghat[28] = 0.1767766952966368*(alpha[3]*fUpwind[31]+alpha[8]*fUpwind[29]+alpha[0]*fUpwind[28]+alpha[2]*fUpwind[23]); 
  Ghat[29] = 0.1767766952966368*(alpha[2]*fUpwind[31]+alpha[0]*fUpwind[29]+alpha[8]*fUpwind[28]+alpha[3]*fUpwind[23]); 
  Ghat[30] = 0.1767766952966368*(alpha[0]*fUpwind[30]+alpha[2]*fUpwind[25]+alpha[3]*fUpwind[24]+alpha[8]*fUpwind[15]); 
  Ghat[31] = 0.1767766952966368*(alpha[0]*fUpwind[31]+alpha[2]*fUpwind[29]+alpha[3]*fUpwind[28]+alpha[8]*fUpwind[23]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -0.7071067811865475*Ghat[3]*dv10; 
  out[4] += -1.224744871391589*Ghat[0]*dv10; 
  out[5] += -0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -0.7071067811865475*Ghat[6]*dv10; 
  out[8] += -0.7071067811865475*Ghat[7]*dv10; 
  out[9] += -0.7071067811865475*Ghat[8]*dv10; 
  out[10] += -1.224744871391589*Ghat[1]*dv10; 
  out[11] += -1.224744871391589*Ghat[2]*dv10; 
  out[12] += -1.224744871391589*Ghat[3]*dv10; 
  out[13] += -0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -0.7071067811865475*Ghat[10]*dv10; 
  out[15] += -0.7071067811865475*Ghat[11]*dv10; 
  out[16] += -1.224744871391589*Ghat[4]*dv10; 
  out[17] += -0.7071067811865475*Ghat[12]*dv10; 
  out[18] += -0.7071067811865475*Ghat[13]*dv10; 
  out[19] += -0.7071067811865475*Ghat[14]*dv10; 
  out[20] += -1.224744871391589*Ghat[5]*dv10; 
  out[21] += -0.7071067811865475*Ghat[15]*dv10; 
  out[22] += -0.7071067811865475*Ghat[16]*dv10; 
  out[23] += -1.224744871391589*Ghat[6]*dv10; 
  out[24] += -1.224744871391589*Ghat[7]*dv10; 
  out[25] += -1.224744871391589*Ghat[8]*dv10; 
  out[26] += -0.7071067811865475*Ghat[17]*dv10; 
  out[27] += -0.7071067811865475*Ghat[18]*dv10; 
  out[28] += -0.7071067811865475*Ghat[19]*dv10; 
  out[29] += -1.224744871391589*Ghat[9]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -1.224744871391589*Ghat[11]*dv10; 
  out[32] += -0.7071067811865475*Ghat[20]*dv10; 
  out[33] += -0.7071067811865475*Ghat[21]*dv10; 
  out[34] += -0.7071067811865475*Ghat[22]*dv10; 
  out[35] += -1.224744871391589*Ghat[12]*dv10; 
  out[36] += -1.224744871391589*Ghat[13]*dv10; 
  out[37] += -1.224744871391589*Ghat[14]*dv10; 
  out[38] += -0.7071067811865475*Ghat[23]*dv10; 
  out[39] += -0.7071067811865475*Ghat[24]*dv10; 
  out[40] += -0.7071067811865475*Ghat[25]*dv10; 
  out[41] += -1.224744871391589*Ghat[15]*dv10; 
  out[42] += -1.224744871391589*Ghat[16]*dv10; 
  out[43] += -0.7071067811865475*Ghat[26]*dv10; 
  out[44] += -1.224744871391589*Ghat[17]*dv10; 
  out[45] += -1.224744871391589*Ghat[18]*dv10; 
  out[46] += -1.224744871391589*Ghat[19]*dv10; 
  out[47] += -0.7071067811865475*Ghat[27]*dv10; 
  out[48] += -1.224744871391589*Ghat[20]*dv10; 
  out[49] += -1.224744871391589*Ghat[21]*dv10; 
  out[50] += -1.224744871391589*Ghat[22]*dv10; 
  out[51] += -0.7071067811865475*Ghat[28]*dv10; 
  out[52] += -0.7071067811865475*Ghat[29]*dv10; 
  out[53] += -0.7071067811865475*Ghat[30]*dv10; 
  out[54] += -1.224744871391589*Ghat[23]*dv10; 
  out[55] += -1.224744871391589*Ghat[24]*dv10; 
  out[56] += -1.224744871391589*Ghat[25]*dv10; 
  out[57] += -1.224744871391589*Ghat[26]*dv10; 
  out[58] += -1.224744871391589*Ghat[27]*dv10; 
  out[59] += -0.7071067811865475*Ghat[31]*dv10; 
  out[60] += -1.224744871391589*Ghat[28]*dv10; 
  out[61] += -1.224744871391589*Ghat[29]*dv10; 
  out[62] += -1.224744871391589*Ghat[30]*dv10; 
  out[63] += -1.224744871391589*Ghat[31]*dv10; 

  } else { 

  if (alpha[8]-alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx4_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_6x_p1_surfx4_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = ser_6x_p1_surfx4_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = ser_6x_p1_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx4_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_6x_p1_surfx4_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = ser_6x_p1_surfx4_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = ser_6x_p1_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if (alpha[8]-alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx4_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = ser_6x_p1_surfx4_eval_quad_node_5_r(fEdge); 
    fUpwindQuad[6] = ser_6x_p1_surfx4_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = ser_6x_p1_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx4_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = ser_6x_p1_surfx4_eval_quad_node_5_l(fSkin); 
    fUpwindQuad[6] = ser_6x_p1_surfx4_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = ser_6x_p1_surfx4_eval_quad_node_7_l(fSkin); 
  } 
  if (alpha[8]-alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_6x_p1_surfx4_eval_quad_node_8_r(fEdge); 
    fUpwindQuad[9] = ser_6x_p1_surfx4_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = ser_6x_p1_surfx4_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = ser_6x_p1_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_6x_p1_surfx4_eval_quad_node_8_l(fSkin); 
    fUpwindQuad[9] = ser_6x_p1_surfx4_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = ser_6x_p1_surfx4_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = ser_6x_p1_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if (alpha[8]-alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_6x_p1_surfx4_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = ser_6x_p1_surfx4_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = ser_6x_p1_surfx4_eval_quad_node_14_r(fEdge); 
    fUpwindQuad[15] = ser_6x_p1_surfx4_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_6x_p1_surfx4_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = ser_6x_p1_surfx4_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = ser_6x_p1_surfx4_eval_quad_node_14_l(fSkin); 
    fUpwindQuad[15] = ser_6x_p1_surfx4_eval_quad_node_15_l(fSkin); 
  } 
  if ((-alpha[8])+alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[16] = ser_6x_p1_surfx4_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = ser_6x_p1_surfx4_eval_quad_node_17_r(fEdge); 
    fUpwindQuad[18] = ser_6x_p1_surfx4_eval_quad_node_18_r(fEdge); 
    fUpwindQuad[19] = ser_6x_p1_surfx4_eval_quad_node_19_r(fEdge); 
  } else { 
    fUpwindQuad[16] = ser_6x_p1_surfx4_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = ser_6x_p1_surfx4_eval_quad_node_17_l(fSkin); 
    fUpwindQuad[18] = ser_6x_p1_surfx4_eval_quad_node_18_l(fSkin); 
    fUpwindQuad[19] = ser_6x_p1_surfx4_eval_quad_node_19_l(fSkin); 
  } 
  if ((-alpha[8])+alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[20] = ser_6x_p1_surfx4_eval_quad_node_20_r(fEdge); 
    fUpwindQuad[21] = ser_6x_p1_surfx4_eval_quad_node_21_r(fEdge); 
    fUpwindQuad[22] = ser_6x_p1_surfx4_eval_quad_node_22_r(fEdge); 
    fUpwindQuad[23] = ser_6x_p1_surfx4_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[20] = ser_6x_p1_surfx4_eval_quad_node_20_l(fSkin); 
    fUpwindQuad[21] = ser_6x_p1_surfx4_eval_quad_node_21_l(fSkin); 
    fUpwindQuad[22] = ser_6x_p1_surfx4_eval_quad_node_22_l(fSkin); 
    fUpwindQuad[23] = ser_6x_p1_surfx4_eval_quad_node_23_l(fSkin); 
  } 
  if ((-alpha[8])+alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[24] = ser_6x_p1_surfx4_eval_quad_node_24_r(fEdge); 
    fUpwindQuad[25] = ser_6x_p1_surfx4_eval_quad_node_25_r(fEdge); 
    fUpwindQuad[26] = ser_6x_p1_surfx4_eval_quad_node_26_r(fEdge); 
    fUpwindQuad[27] = ser_6x_p1_surfx4_eval_quad_node_27_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_6x_p1_surfx4_eval_quad_node_24_l(fSkin); 
    fUpwindQuad[25] = ser_6x_p1_surfx4_eval_quad_node_25_l(fSkin); 
    fUpwindQuad[26] = ser_6x_p1_surfx4_eval_quad_node_26_l(fSkin); 
    fUpwindQuad[27] = ser_6x_p1_surfx4_eval_quad_node_27_l(fSkin); 
  } 
  if ((-alpha[8])+alpha[3]-alpha[2]+alpha[0] > 0) { 
    fUpwindQuad[28] = ser_6x_p1_surfx4_eval_quad_node_28_r(fEdge); 
    fUpwindQuad[29] = ser_6x_p1_surfx4_eval_quad_node_29_r(fEdge); 
    fUpwindQuad[30] = ser_6x_p1_surfx4_eval_quad_node_30_r(fEdge); 
    fUpwindQuad[31] = ser_6x_p1_surfx4_eval_quad_node_31_r(fEdge); 
  } else { 
    fUpwindQuad[28] = ser_6x_p1_surfx4_eval_quad_node_28_l(fSkin); 
    fUpwindQuad[29] = ser_6x_p1_surfx4_eval_quad_node_29_l(fSkin); 
    fUpwindQuad[30] = ser_6x_p1_surfx4_eval_quad_node_30_l(fSkin); 
    fUpwindQuad[31] = ser_6x_p1_surfx4_eval_quad_node_31_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*(alpha[8]*fUpwind[8]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.1767766952966368*(alpha[8]*fUpwind[16]+alpha[3]*fUpwind[7]+alpha[2]*fUpwind[6]+alpha[0]*fUpwind[1]); 
  Ghat[2] = 0.1767766952966368*(alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.1767766952966368*(alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.1767766952966368*(alpha[8]*fUpwind[19]+alpha[3]*fUpwind[11]+alpha[2]*fUpwind[10]+alpha[0]*fUpwind[4]); 
  Ghat[5] = 0.1767766952966368*(alpha[8]*fUpwind[22]+alpha[3]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[0]*fUpwind[5]); 
  Ghat[6] = 0.1767766952966368*(alpha[3]*fUpwind[16]+fUpwind[7]*alpha[8]+alpha[0]*fUpwind[6]+fUpwind[1]*alpha[2]); 
  Ghat[7] = 0.1767766952966368*(alpha[2]*fUpwind[16]+fUpwind[6]*alpha[8]+alpha[0]*fUpwind[7]+fUpwind[1]*alpha[3]); 
  Ghat[8] = 0.1767766952966368*(alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[9] = 0.1767766952966368*(alpha[8]*fUpwind[26]+alpha[3]*fUpwind[18]+alpha[2]*fUpwind[17]+alpha[0]*fUpwind[9]); 
  Ghat[10] = 0.1767766952966368*(alpha[3]*fUpwind[19]+alpha[8]*fUpwind[11]+alpha[0]*fUpwind[10]+alpha[2]*fUpwind[4]); 
  Ghat[11] = 0.1767766952966368*(alpha[2]*fUpwind[19]+alpha[0]*fUpwind[11]+alpha[8]*fUpwind[10]+alpha[3]*fUpwind[4]); 
  Ghat[12] = 0.1767766952966368*(alpha[8]*fUpwind[27]+alpha[3]*fUpwind[21]+alpha[2]*fUpwind[20]+alpha[0]*fUpwind[12]); 
  Ghat[13] = 0.1767766952966368*(alpha[3]*fUpwind[22]+alpha[8]*fUpwind[14]+alpha[0]*fUpwind[13]+alpha[2]*fUpwind[5]); 
  Ghat[14] = 0.1767766952966368*(alpha[2]*fUpwind[22]+alpha[0]*fUpwind[14]+alpha[8]*fUpwind[13]+alpha[3]*fUpwind[5]); 
  Ghat[15] = 0.1767766952966368*(alpha[8]*fUpwind[30]+alpha[3]*fUpwind[25]+alpha[2]*fUpwind[24]+alpha[0]*fUpwind[15]); 
  Ghat[16] = 0.1767766952966368*(alpha[0]*fUpwind[16]+fUpwind[1]*alpha[8]+alpha[2]*fUpwind[7]+alpha[3]*fUpwind[6]); 
  Ghat[17] = 0.1767766952966368*(alpha[3]*fUpwind[26]+alpha[8]*fUpwind[18]+alpha[0]*fUpwind[17]+alpha[2]*fUpwind[9]); 
  Ghat[18] = 0.1767766952966368*(alpha[2]*fUpwind[26]+alpha[0]*fUpwind[18]+alpha[8]*fUpwind[17]+alpha[3]*fUpwind[9]); 
  Ghat[19] = 0.1767766952966368*(alpha[0]*fUpwind[19]+alpha[2]*fUpwind[11]+alpha[3]*fUpwind[10]+fUpwind[4]*alpha[8]); 
  Ghat[20] = 0.1767766952966368*(alpha[3]*fUpwind[27]+alpha[8]*fUpwind[21]+alpha[0]*fUpwind[20]+alpha[2]*fUpwind[12]); 
  Ghat[21] = 0.1767766952966368*(alpha[2]*fUpwind[27]+alpha[0]*fUpwind[21]+alpha[8]*fUpwind[20]+alpha[3]*fUpwind[12]); 
  Ghat[22] = 0.1767766952966368*(alpha[0]*fUpwind[22]+alpha[2]*fUpwind[14]+alpha[3]*fUpwind[13]+fUpwind[5]*alpha[8]); 
  Ghat[23] = 0.1767766952966368*(alpha[8]*fUpwind[31]+alpha[3]*fUpwind[29]+alpha[2]*fUpwind[28]+alpha[0]*fUpwind[23]); 
  Ghat[24] = 0.1767766952966368*(alpha[3]*fUpwind[30]+alpha[8]*fUpwind[25]+alpha[0]*fUpwind[24]+alpha[2]*fUpwind[15]); 
  Ghat[25] = 0.1767766952966368*(alpha[2]*fUpwind[30]+alpha[0]*fUpwind[25]+alpha[8]*fUpwind[24]+alpha[3]*fUpwind[15]); 
  Ghat[26] = 0.1767766952966368*(alpha[0]*fUpwind[26]+alpha[2]*fUpwind[18]+alpha[3]*fUpwind[17]+alpha[8]*fUpwind[9]); 
  Ghat[27] = 0.1767766952966368*(alpha[0]*fUpwind[27]+alpha[2]*fUpwind[21]+alpha[3]*fUpwind[20]+alpha[8]*fUpwind[12]); 
  Ghat[28] = 0.1767766952966368*(alpha[3]*fUpwind[31]+alpha[8]*fUpwind[29]+alpha[0]*fUpwind[28]+alpha[2]*fUpwind[23]); 
  Ghat[29] = 0.1767766952966368*(alpha[2]*fUpwind[31]+alpha[0]*fUpwind[29]+alpha[8]*fUpwind[28]+alpha[3]*fUpwind[23]); 
  Ghat[30] = 0.1767766952966368*(alpha[0]*fUpwind[30]+alpha[2]*fUpwind[25]+alpha[3]*fUpwind[24]+alpha[8]*fUpwind[15]); 
  Ghat[31] = 0.1767766952966368*(alpha[0]*fUpwind[31]+alpha[2]*fUpwind[29]+alpha[3]*fUpwind[28]+alpha[8]*fUpwind[23]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += 0.7071067811865475*Ghat[3]*dv10; 
  out[4] += -1.224744871391589*Ghat[0]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += 0.7071067811865475*Ghat[5]*dv10; 
  out[7] += 0.7071067811865475*Ghat[6]*dv10; 
  out[8] += 0.7071067811865475*Ghat[7]*dv10; 
  out[9] += 0.7071067811865475*Ghat[8]*dv10; 
  out[10] += -1.224744871391589*Ghat[1]*dv10; 
  out[11] += -1.224744871391589*Ghat[2]*dv10; 
  out[12] += -1.224744871391589*Ghat[3]*dv10; 
  out[13] += 0.7071067811865475*Ghat[9]*dv10; 
  out[14] += 0.7071067811865475*Ghat[10]*dv10; 
  out[15] += 0.7071067811865475*Ghat[11]*dv10; 
  out[16] += -1.224744871391589*Ghat[4]*dv10; 
  out[17] += 0.7071067811865475*Ghat[12]*dv10; 
  out[18] += 0.7071067811865475*Ghat[13]*dv10; 
  out[19] += 0.7071067811865475*Ghat[14]*dv10; 
  out[20] += -1.224744871391589*Ghat[5]*dv10; 
  out[21] += 0.7071067811865475*Ghat[15]*dv10; 
  out[22] += 0.7071067811865475*Ghat[16]*dv10; 
  out[23] += -1.224744871391589*Ghat[6]*dv10; 
  out[24] += -1.224744871391589*Ghat[7]*dv10; 
  out[25] += -1.224744871391589*Ghat[8]*dv10; 
  out[26] += 0.7071067811865475*Ghat[17]*dv10; 
  out[27] += 0.7071067811865475*Ghat[18]*dv10; 
  out[28] += 0.7071067811865475*Ghat[19]*dv10; 
  out[29] += -1.224744871391589*Ghat[9]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -1.224744871391589*Ghat[11]*dv10; 
  out[32] += 0.7071067811865475*Ghat[20]*dv10; 
  out[33] += 0.7071067811865475*Ghat[21]*dv10; 
  out[34] += 0.7071067811865475*Ghat[22]*dv10; 
  out[35] += -1.224744871391589*Ghat[12]*dv10; 
  out[36] += -1.224744871391589*Ghat[13]*dv10; 
  out[37] += -1.224744871391589*Ghat[14]*dv10; 
  out[38] += 0.7071067811865475*Ghat[23]*dv10; 
  out[39] += 0.7071067811865475*Ghat[24]*dv10; 
  out[40] += 0.7071067811865475*Ghat[25]*dv10; 
  out[41] += -1.224744871391589*Ghat[15]*dv10; 
  out[42] += -1.224744871391589*Ghat[16]*dv10; 
  out[43] += 0.7071067811865475*Ghat[26]*dv10; 
  out[44] += -1.224744871391589*Ghat[17]*dv10; 
  out[45] += -1.224744871391589*Ghat[18]*dv10; 
  out[46] += -1.224744871391589*Ghat[19]*dv10; 
  out[47] += 0.7071067811865475*Ghat[27]*dv10; 
  out[48] += -1.224744871391589*Ghat[20]*dv10; 
  out[49] += -1.224744871391589*Ghat[21]*dv10; 
  out[50] += -1.224744871391589*Ghat[22]*dv10; 
  out[51] += 0.7071067811865475*Ghat[28]*dv10; 
  out[52] += 0.7071067811865475*Ghat[29]*dv10; 
  out[53] += 0.7071067811865475*Ghat[30]*dv10; 
  out[54] += -1.224744871391589*Ghat[23]*dv10; 
  out[55] += -1.224744871391589*Ghat[24]*dv10; 
  out[56] += -1.224744871391589*Ghat[25]*dv10; 
  out[57] += -1.224744871391589*Ghat[26]*dv10; 
  out[58] += -1.224744871391589*Ghat[27]*dv10; 
  out[59] += 0.7071067811865475*Ghat[31]*dv10; 
  out[60] += -1.224744871391589*Ghat[28]*dv10; 
  out[61] += -1.224744871391589*Ghat[29]*dv10; 
  out[62] += -1.224744871391589*Ghat[30]*dv10; 
  out[63] += -1.224744871391589*Ghat[31]*dv10; 

  } 
} 
