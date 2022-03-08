#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p1_surfx3_quad.h> 
#include <gkyl_basis_ser_5x_p1_upwind.h> 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
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
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[4]; 
  const double *A2 = &vecA[8]; 
  double alpha[16] = {0.0}; 

  alpha[0] = 3.464101615137754*A2[1]*dx10*wv3-3.464101615137754*A0[2]*dx11*wv2+3.464101615137754*A1[1]*dx10*wv2-3.464101615137754*phi[1]*dx10; 
  alpha[1] = -3.464101615137754*A0[3]*dx11*wv2; 
  alpha[2] = 3.464101615137754*A2[3]*dx10*wv3+3.464101615137754*A1[3]*dx10*wv2-3.464101615137754*phi[3]*dx10; 
  alpha[3] = A1[1]*dv2*dx10-1.0*A0[2]*dv2*dx11; 
  alpha[4] = A2[1]*dv3*dx10; 
  alpha[6] = -1.0*A0[3]*dv2*dx11; 
  alpha[7] = A1[3]*dv2*dx10; 
  alpha[9] = A2[3]*dv3*dx10; 

  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  if (alpha[9]+alpha[7]+alpha[6]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx3_quad_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx3_quad_0_l(fEdge); 
  } 
  if (alpha[9]+alpha[7]-alpha[6]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx3_quad_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx3_quad_1_l(fEdge); 
  } 
  if ((-alpha[9])-alpha[7]+alpha[6]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx3_quad_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx3_quad_2_l(fEdge); 
  } 
  if ((-alpha[9])-alpha[7]-alpha[6]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx3_quad_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx3_quad_3_l(fEdge); 
  } 
  if (alpha[9]-alpha[7]-alpha[6]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p1_surfx3_quad_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_5x_p1_surfx3_quad_4_l(fEdge); 
  } 
  if (alpha[9]-alpha[7]+alpha[6]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p1_surfx3_quad_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_5x_p1_surfx3_quad_5_l(fEdge); 
  } 
  if ((-alpha[9])+alpha[7]-alpha[6]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p1_surfx3_quad_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_5x_p1_surfx3_quad_6_l(fEdge); 
  } 
  if ((-alpha[9])+alpha[7]+alpha[6]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p1_surfx3_quad_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_5x_p1_surfx3_quad_7_l(fEdge); 
  } 
  if ((-alpha[9])+alpha[7]+alpha[6]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p1_surfx3_quad_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_5x_p1_surfx3_quad_8_l(fEdge); 
  } 
  if ((-alpha[9])+alpha[7]-alpha[6]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p1_surfx3_quad_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_5x_p1_surfx3_quad_9_l(fEdge); 
  } 
  if (alpha[9]-alpha[7]+alpha[6]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p1_surfx3_quad_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_5x_p1_surfx3_quad_10_l(fEdge); 
  } 
  if (alpha[9]-alpha[7]-alpha[6]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p1_surfx3_quad_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_5x_p1_surfx3_quad_11_l(fEdge); 
  } 
  if ((-alpha[9])-alpha[7]-alpha[6]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p1_surfx3_quad_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_5x_p1_surfx3_quad_12_l(fEdge); 
  } 
  if ((-alpha[9])-alpha[7]+alpha[6]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p1_surfx3_quad_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_5x_p1_surfx3_quad_13_l(fEdge); 
  } 
  if (alpha[9]+alpha[7]-alpha[6]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p1_surfx3_quad_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_5x_p1_surfx3_quad_14_l(fEdge); 
  } 
  if (alpha[9]+alpha[7]+alpha[6]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p1_surfx3_quad_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_5x_p1_surfx3_quad_15_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.25*(alpha[9]*fUpwind[9]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.25*(alpha[9]*fUpwind[12]+alpha[7]*fUpwind[11]+alpha[4]*fUpwind[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.25*(alpha[6]*fUpwind[11]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.25*(alpha[9]*fUpwind[14]+alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.25*(alpha[7]*fUpwind[14]+alpha[6]*fUpwind[13]+alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] += 0.25*(alpha[4]*fUpwind[12]+alpha[3]*fUpwind[11]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] += 0.25*(alpha[9]*fUpwind[15]+alpha[4]*fUpwind[13]+alpha[2]*fUpwind[11]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += 0.25*(alpha[4]*fUpwind[14]+alpha[1]*fUpwind[11]+alpha[9]*fUpwind[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] += 0.25*(alpha[7]*fUpwind[15]+alpha[3]*fUpwind[13]+alpha[2]*fUpwind[12]+alpha[6]*fUpwind[10]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] += 0.25*(alpha[6]*fUpwind[15]+alpha[3]*fUpwind[14]+alpha[1]*fUpwind[12]+alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] += 0.25*(alpha[2]*fUpwind[14]+alpha[1]*fUpwind[13]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.25*(alpha[4]*fUpwind[15]+alpha[9]*fUpwind[13]+alpha[0]*fUpwind[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]); 
  Ghat[12] += 0.25*(alpha[3]*fUpwind[15]+alpha[6]*fUpwind[14]+alpha[7]*fUpwind[13]+alpha[0]*fUpwind[12]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+alpha[4]*fUpwind[5]); 
  Ghat[13] += 0.25*(alpha[2]*fUpwind[15]+alpha[0]*fUpwind[13]+alpha[7]*fUpwind[12]+alpha[9]*fUpwind[11]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[14] += 0.25*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]+alpha[6]*fUpwind[12]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[15] += 0.25*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[3]*fUpwind[12]+alpha[4]*fUpwind[11]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += -0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += -0.7071067811865475*Ghat[8]*dv10; 
  out[13] += -0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += -0.7071067811865475*Ghat[10]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += -0.7071067811865475*Ghat[11]*dv10; 
  out[18] += -1.224744871391589*Ghat[6]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 
  out[20] += -0.7071067811865475*Ghat[12]*dv10; 
  out[21] += -1.224744871391589*Ghat[8]*dv10; 
  out[22] += -1.224744871391589*Ghat[9]*dv10; 
  out[23] += -0.7071067811865475*Ghat[13]*dv10; 
  out[24] += -0.7071067811865475*Ghat[14]*dv10; 
  out[25] += -1.224744871391589*Ghat[10]*dv10; 
  out[26] += -1.224744871391589*Ghat[11]*dv10; 
  out[27] += -1.224744871391589*Ghat[12]*dv10; 
  out[28] += -0.7071067811865475*Ghat[15]*dv10; 
  out[29] += -1.224744871391589*Ghat[13]*dv10; 
  out[30] += -1.224744871391589*Ghat[14]*dv10; 
  out[31] += -1.224744871391589*Ghat[15]*dv10; 

  } else { 

  if (alpha[9]+alpha[7]+alpha[6]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx3_quad_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx3_quad_0_l(fSkin); 
  } 
  if (alpha[9]+alpha[7]-alpha[6]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx3_quad_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx3_quad_1_l(fSkin); 
  } 
  if ((-alpha[9])-alpha[7]+alpha[6]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx3_quad_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx3_quad_2_l(fSkin); 
  } 
  if ((-alpha[9])-alpha[7]-alpha[6]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx3_quad_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx3_quad_3_l(fSkin); 
  } 
  if (alpha[9]-alpha[7]-alpha[6]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p1_surfx3_quad_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_5x_p1_surfx3_quad_4_l(fSkin); 
  } 
  if (alpha[9]-alpha[7]+alpha[6]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p1_surfx3_quad_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_5x_p1_surfx3_quad_5_l(fSkin); 
  } 
  if ((-alpha[9])+alpha[7]-alpha[6]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p1_surfx3_quad_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_5x_p1_surfx3_quad_6_l(fSkin); 
  } 
  if ((-alpha[9])+alpha[7]+alpha[6]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p1_surfx3_quad_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_5x_p1_surfx3_quad_7_l(fSkin); 
  } 
  if ((-alpha[9])+alpha[7]+alpha[6]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p1_surfx3_quad_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_5x_p1_surfx3_quad_8_l(fSkin); 
  } 
  if ((-alpha[9])+alpha[7]-alpha[6]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p1_surfx3_quad_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_5x_p1_surfx3_quad_9_l(fSkin); 
  } 
  if (alpha[9]-alpha[7]+alpha[6]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p1_surfx3_quad_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_5x_p1_surfx3_quad_10_l(fSkin); 
  } 
  if (alpha[9]-alpha[7]-alpha[6]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p1_surfx3_quad_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_5x_p1_surfx3_quad_11_l(fSkin); 
  } 
  if ((-alpha[9])-alpha[7]-alpha[6]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p1_surfx3_quad_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_5x_p1_surfx3_quad_12_l(fSkin); 
  } 
  if ((-alpha[9])-alpha[7]+alpha[6]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p1_surfx3_quad_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_5x_p1_surfx3_quad_13_l(fSkin); 
  } 
  if (alpha[9]+alpha[7]-alpha[6]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p1_surfx3_quad_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_5x_p1_surfx3_quad_14_l(fSkin); 
  } 
  if (alpha[9]+alpha[7]+alpha[6]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p1_surfx3_quad_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_5x_p1_surfx3_quad_15_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p1_upwind(fUpwindQuad, fUpwind); 

  Ghat[0] += 0.25*(alpha[9]*fUpwind[9]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.25*(alpha[9]*fUpwind[12]+alpha[7]*fUpwind[11]+alpha[4]*fUpwind[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.25*(alpha[6]*fUpwind[11]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.25*(alpha[9]*fUpwind[14]+alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.25*(alpha[7]*fUpwind[14]+alpha[6]*fUpwind[13]+alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] += 0.25*(alpha[4]*fUpwind[12]+alpha[3]*fUpwind[11]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] += 0.25*(alpha[9]*fUpwind[15]+alpha[4]*fUpwind[13]+alpha[2]*fUpwind[11]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] += 0.25*(alpha[4]*fUpwind[14]+alpha[1]*fUpwind[11]+alpha[9]*fUpwind[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] += 0.25*(alpha[7]*fUpwind[15]+alpha[3]*fUpwind[13]+alpha[2]*fUpwind[12]+alpha[6]*fUpwind[10]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] += 0.25*(alpha[6]*fUpwind[15]+alpha[3]*fUpwind[14]+alpha[1]*fUpwind[12]+alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] += 0.25*(alpha[2]*fUpwind[14]+alpha[1]*fUpwind[13]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.25*(alpha[4]*fUpwind[15]+alpha[9]*fUpwind[13]+alpha[0]*fUpwind[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]); 
  Ghat[12] += 0.25*(alpha[3]*fUpwind[15]+alpha[6]*fUpwind[14]+alpha[7]*fUpwind[13]+alpha[0]*fUpwind[12]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+alpha[4]*fUpwind[5]); 
  Ghat[13] += 0.25*(alpha[2]*fUpwind[15]+alpha[0]*fUpwind[13]+alpha[7]*fUpwind[12]+alpha[9]*fUpwind[11]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[14] += 0.25*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]+alpha[6]*fUpwind[12]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[15] += 0.25*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[3]*fUpwind[12]+alpha[4]*fUpwind[11]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += 0.7071067811865475*Ghat[5]*dv10; 
  out[7] += -1.224744871391589*Ghat[1]*dv10; 
  out[8] += -1.224744871391589*Ghat[2]*dv10; 
  out[9] += 0.7071067811865475*Ghat[6]*dv10; 
  out[10] += 0.7071067811865475*Ghat[7]*dv10; 
  out[11] += -1.224744871391589*Ghat[3]*dv10; 
  out[12] += 0.7071067811865475*Ghat[8]*dv10; 
  out[13] += 0.7071067811865475*Ghat[9]*dv10; 
  out[14] += -1.224744871391589*Ghat[4]*dv10; 
  out[15] += 0.7071067811865475*Ghat[10]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += 0.7071067811865475*Ghat[11]*dv10; 
  out[18] += -1.224744871391589*Ghat[6]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 
  out[20] += 0.7071067811865475*Ghat[12]*dv10; 
  out[21] += -1.224744871391589*Ghat[8]*dv10; 
  out[22] += -1.224744871391589*Ghat[9]*dv10; 
  out[23] += 0.7071067811865475*Ghat[13]*dv10; 
  out[24] += 0.7071067811865475*Ghat[14]*dv10; 
  out[25] += -1.224744871391589*Ghat[10]*dv10; 
  out[26] += -1.224744871391589*Ghat[11]*dv10; 
  out[27] += -1.224744871391589*Ghat[12]*dv10; 
  out[28] += 0.7071067811865475*Ghat[15]*dv10; 
  out[29] += -1.224744871391589*Ghat[13]*dv10; 
  out[30] += -1.224744871391589*Ghat[14]*dv10; 
  out[31] += -1.224744871391589*Ghat[15]*dv10; 

  } 
} 
