#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_1x3v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_1x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_ext_phiA_boundary_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double *pots, const double *EBext, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // pots:        potentials phi_tot=phi+phi_ext and A_ext (scaled by q/m).
  // EBext:       external E and B fields (scaled by q/m).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 

  const double dv12 = 2/dxv[3]; 
  const double dx10 = 2/dxv[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 

  double alpha[16] = {0.0}; 

  const double *phi = &pots[0]; 

  const double *Ax = &pots[2]; 
  const double *Ay = &pots[4]; 
  const double *Az = &pots[6]; 

  alpha[0] = -(3.4641016151377544*Az[1]*dx10*wv1); 
  alpha[2] = -(1.0*Az[1]*dv1*dx10); 

  double fUpwindQuad[18] = {0.0};
  double fUpwind[16] = {0.0};
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx4_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx4_eval_quad_node_0_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx4_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx4_eval_quad_node_1_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx4_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx4_eval_quad_node_4_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx4_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx4_eval_quad_node_6_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx4_eval_quad_node_7_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx4_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx4_eval_quad_node_8_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx4_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx4_eval_quad_node_9_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx4_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx4_eval_quad_node_10_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx4_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx4_eval_quad_node_12_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx4_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx4_eval_quad_node_13_l(fEdge); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx4_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx4_eval_quad_node_14_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx4_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx4_eval_quad_node_15_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx4_eval_quad_node_16_r(fSkin); 
  } else { 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx4_eval_quad_node_16_l(fEdge); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx4_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx4_eval_quad_node_17_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[2]*fUpwind[2]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[2]*fUpwind[4]+0.3535533905932737*alpha[0]*fUpwind[1]; 
  Ghat[2] = 0.3162277660168379*alpha[2]*fUpwind[8]+0.3535533905932737*alpha[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.3535533905932737*alpha[2]*fUpwind[6]+0.3535533905932737*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.31622776601683794*alpha[2]*fUpwind[9]+0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alpha[2]; 
  Ghat[5] = 0.3535533905932737*alpha[2]*fUpwind[7]+0.3535533905932737*alpha[0]*fUpwind[5]; 
  Ghat[6] = 0.31622776601683794*alpha[2]*fUpwind[10]+0.3535533905932737*alpha[0]*fUpwind[6]+0.3535533905932737*alpha[2]*fUpwind[3]; 
  Ghat[7] = 0.3162277660168379*alpha[2]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[7]+0.3535533905932737*alpha[2]*fUpwind[5]; 
  Ghat[8] = 0.3535533905932737*alpha[0]*fUpwind[8]+0.3162277660168379*alpha[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alpha[0]*fUpwind[9]+0.31622776601683794*alpha[2]*fUpwind[4]; 
  Ghat[10] = 0.3535533905932737*alpha[0]*fUpwind[10]+0.31622776601683794*alpha[2]*fUpwind[6]; 
  Ghat[11] = 0.3535533905932737*alpha[0]*fUpwind[11]+0.3162277660168379*alpha[2]*fUpwind[7]; 
  Ghat[12] = 0.3535533905932737*alpha[2]*fUpwind[14]+0.3535533905932737*alpha[0]*fUpwind[12]; 
  Ghat[13] = 0.3535533905932737*alpha[2]*fUpwind[15]+0.3535533905932737*alpha[0]*fUpwind[13]; 
  Ghat[14] = 0.3535533905932737*alpha[0]*fUpwind[14]+0.3535533905932737*alpha[2]*fUpwind[12]; 
  Ghat[15] = 0.3535533905932737*alpha[0]*fUpwind[15]+0.3535533905932737*alpha[2]*fUpwind[13]; 

  out[0] += -(0.7071067811865475*Ghat[0]*dv12); 
  out[1] += -(0.7071067811865475*Ghat[1]*dv12); 
  out[2] += -(0.7071067811865475*Ghat[2]*dv12); 
  out[3] += -(0.7071067811865475*Ghat[3]*dv12); 
  out[4] += -(1.224744871391589*Ghat[0]*dv12); 
  out[5] += -(0.7071067811865475*Ghat[4]*dv12); 
  out[6] += -(0.7071067811865475*Ghat[5]*dv12); 
  out[7] += -(0.7071067811865475*Ghat[6]*dv12); 
  out[8] += -(1.224744871391589*Ghat[1]*dv12); 
  out[9] += -(1.224744871391589*Ghat[2]*dv12); 
  out[10] += -(1.224744871391589*Ghat[3]*dv12); 
  out[11] += -(0.7071067811865475*Ghat[7]*dv12); 
  out[12] += -(1.224744871391589*Ghat[4]*dv12); 
  out[13] += -(1.224744871391589*Ghat[5]*dv12); 
  out[14] += -(1.224744871391589*Ghat[6]*dv12); 
  out[15] += -(1.224744871391589*Ghat[7]*dv12); 
  out[16] += -(0.7071067811865475*Ghat[8]*dv12); 
  out[17] += -(0.7071067811865475*Ghat[9]*dv12); 
  out[18] += -(0.7071067811865475*Ghat[10]*dv12); 
  out[19] += -(1.224744871391589*Ghat[8]*dv12); 
  out[20] += -(0.7071067811865475*Ghat[11]*dv12); 
  out[21] += -(1.224744871391589*Ghat[9]*dv12); 
  out[22] += -(1.224744871391589*Ghat[10]*dv12); 
  out[23] += -(1.224744871391589*Ghat[11]*dv12); 
  out[24] += -(0.7071067811865475*Ghat[12]*dv12); 
  out[25] += -(0.7071067811865475*Ghat[13]*dv12); 
  out[26] += -(0.7071067811865475*Ghat[14]*dv12); 
  out[27] += -(1.224744871391589*Ghat[12]*dv12); 
  out[28] += -(0.7071067811865475*Ghat[15]*dv12); 
  out[29] += -(1.224744871391589*Ghat[13]*dv12); 
  out[30] += -(1.224744871391589*Ghat[14]*dv12); 
  out[31] += -(1.224744871391589*Ghat[15]*dv12); 
  out[32] += -(1.5811388300841895*Ghat[0]*dv12); 
  out[33] += -(1.5811388300841898*Ghat[1]*dv12); 
  out[34] += -(1.5811388300841898*Ghat[2]*dv12); 
  out[35] += -(1.5811388300841898*Ghat[3]*dv12); 
  out[36] += -(1.5811388300841895*Ghat[4]*dv12); 
  out[37] += -(1.5811388300841895*Ghat[5]*dv12); 
  out[38] += -(1.5811388300841895*Ghat[6]*dv12); 
  out[39] += -(1.5811388300841898*Ghat[7]*dv12); 

  } else { 

  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx4_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x3v_p1_surfx4_eval_quad_node_0_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx4_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_1x3v_p1_surfx4_eval_quad_node_1_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_1x3v_p1_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_1x3v_p1_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx4_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = hyb_1x3v_p1_surfx4_eval_quad_node_4_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = hyb_1x3v_p1_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx4_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = hyb_1x3v_p1_surfx4_eval_quad_node_6_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = hyb_1x3v_p1_surfx4_eval_quad_node_7_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx4_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = hyb_1x3v_p1_surfx4_eval_quad_node_8_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx4_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_1x3v_p1_surfx4_eval_quad_node_9_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx4_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = hyb_1x3v_p1_surfx4_eval_quad_node_10_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.4743416490252568*alpha[2] > 0) { 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = hyb_1x3v_p1_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx4_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = hyb_1x3v_p1_surfx4_eval_quad_node_12_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx4_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = hyb_1x3v_p1_surfx4_eval_quad_node_13_l(fSkin); 
  } 
  if (0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx4_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = hyb_1x3v_p1_surfx4_eval_quad_node_14_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx4_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = hyb_1x3v_p1_surfx4_eval_quad_node_15_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx4_eval_quad_node_16_r(fEdge); 
  } else { 
    fUpwindQuad[16] = hyb_1x3v_p1_surfx4_eval_quad_node_16_l(fSkin); 
  } 
  if (0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx4_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[17] = hyb_1x3v_p1_surfx4_eval_quad_node_17_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alpha[2]*fUpwind[2]+0.3535533905932737*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alpha[2]*fUpwind[4]+0.3535533905932737*alpha[0]*fUpwind[1]; 
  Ghat[2] = 0.3162277660168379*alpha[2]*fUpwind[8]+0.3535533905932737*alpha[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.3535533905932737*alpha[2]*fUpwind[6]+0.3535533905932737*alpha[0]*fUpwind[3]; 
  Ghat[4] = 0.31622776601683794*alpha[2]*fUpwind[9]+0.3535533905932737*alpha[0]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alpha[2]; 
  Ghat[5] = 0.3535533905932737*alpha[2]*fUpwind[7]+0.3535533905932737*alpha[0]*fUpwind[5]; 
  Ghat[6] = 0.31622776601683794*alpha[2]*fUpwind[10]+0.3535533905932737*alpha[0]*fUpwind[6]+0.3535533905932737*alpha[2]*fUpwind[3]; 
  Ghat[7] = 0.3162277660168379*alpha[2]*fUpwind[11]+0.3535533905932737*alpha[0]*fUpwind[7]+0.3535533905932737*alpha[2]*fUpwind[5]; 
  Ghat[8] = 0.3535533905932737*alpha[0]*fUpwind[8]+0.3162277660168379*alpha[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alpha[0]*fUpwind[9]+0.31622776601683794*alpha[2]*fUpwind[4]; 
  Ghat[10] = 0.3535533905932737*alpha[0]*fUpwind[10]+0.31622776601683794*alpha[2]*fUpwind[6]; 
  Ghat[11] = 0.3535533905932737*alpha[0]*fUpwind[11]+0.3162277660168379*alpha[2]*fUpwind[7]; 
  Ghat[12] = 0.3535533905932737*alpha[2]*fUpwind[14]+0.3535533905932737*alpha[0]*fUpwind[12]; 
  Ghat[13] = 0.3535533905932737*alpha[2]*fUpwind[15]+0.3535533905932737*alpha[0]*fUpwind[13]; 
  Ghat[14] = 0.3535533905932737*alpha[0]*fUpwind[14]+0.3535533905932737*alpha[2]*fUpwind[12]; 
  Ghat[15] = 0.3535533905932737*alpha[0]*fUpwind[15]+0.3535533905932737*alpha[2]*fUpwind[13]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv12; 
  out[1] += 0.7071067811865475*Ghat[1]*dv12; 
  out[2] += 0.7071067811865475*Ghat[2]*dv12; 
  out[3] += 0.7071067811865475*Ghat[3]*dv12; 
  out[4] += -(1.224744871391589*Ghat[0]*dv12); 
  out[5] += 0.7071067811865475*Ghat[4]*dv12; 
  out[6] += 0.7071067811865475*Ghat[5]*dv12; 
  out[7] += 0.7071067811865475*Ghat[6]*dv12; 
  out[8] += -(1.224744871391589*Ghat[1]*dv12); 
  out[9] += -(1.224744871391589*Ghat[2]*dv12); 
  out[10] += -(1.224744871391589*Ghat[3]*dv12); 
  out[11] += 0.7071067811865475*Ghat[7]*dv12; 
  out[12] += -(1.224744871391589*Ghat[4]*dv12); 
  out[13] += -(1.224744871391589*Ghat[5]*dv12); 
  out[14] += -(1.224744871391589*Ghat[6]*dv12); 
  out[15] += -(1.224744871391589*Ghat[7]*dv12); 
  out[16] += 0.7071067811865475*Ghat[8]*dv12; 
  out[17] += 0.7071067811865475*Ghat[9]*dv12; 
  out[18] += 0.7071067811865475*Ghat[10]*dv12; 
  out[19] += -(1.224744871391589*Ghat[8]*dv12); 
  out[20] += 0.7071067811865475*Ghat[11]*dv12; 
  out[21] += -(1.224744871391589*Ghat[9]*dv12); 
  out[22] += -(1.224744871391589*Ghat[10]*dv12); 
  out[23] += -(1.224744871391589*Ghat[11]*dv12); 
  out[24] += 0.7071067811865475*Ghat[12]*dv12; 
  out[25] += 0.7071067811865475*Ghat[13]*dv12; 
  out[26] += 0.7071067811865475*Ghat[14]*dv12; 
  out[27] += -(1.224744871391589*Ghat[12]*dv12); 
  out[28] += 0.7071067811865475*Ghat[15]*dv12; 
  out[29] += -(1.224744871391589*Ghat[13]*dv12); 
  out[30] += -(1.224744871391589*Ghat[14]*dv12); 
  out[31] += -(1.224744871391589*Ghat[15]*dv12); 
  out[32] += 1.5811388300841895*Ghat[0]*dv12; 
  out[33] += 1.5811388300841898*Ghat[1]*dv12; 
  out[34] += 1.5811388300841898*Ghat[2]*dv12; 
  out[35] += 1.5811388300841898*Ghat[3]*dv12; 
  out[36] += 1.5811388300841895*Ghat[4]*dv12; 
  out[37] += 1.5811388300841895*Ghat[5]*dv12; 
  out[38] += 1.5811388300841895*Ghat[6]*dv12; 
  out[39] += 1.5811388300841898*Ghat[7]*dv12; 

  } 
  return 0.;

} 
