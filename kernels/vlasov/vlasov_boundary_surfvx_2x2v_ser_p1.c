#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x2v_p1_upwind_quad_to_modal.h> 
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

  double alpha[12] = {0.0}; 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 1.414213562373095*(B2[2]*wv2+E0[2]); 
  alpha[3] = 0.408248290463863*B2[0]*dv2; 
  alpha[4] = 1.414213562373095*(B2[3]*wv2+E0[3]); 
  alpha[5] = 0.408248290463863*B2[1]*dv2; 
  alpha[6] = 0.408248290463863*B2[2]*dv2; 
  alpha[7] = 0.408248290463863*B2[3]*dv2; 

  double fUpwindQuad[12] = {0.0};
  double fUpwind[12] = {0.0};
  double Ghat[12] = {0.0}; 

  if (edge == -1) { 

  if ((-141200537061855*alpha[7])+141200537061855*alpha[6]+141200537061855*alpha[5]+105244666443262*alpha[4]-141200537061855*alpha[3]-105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if (105244666443262*alpha[4]-105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (141200537061855*alpha[7]-141200537061855*alpha[6]-141200537061855*alpha[5]+105244666443262*alpha[4]+141200537061855*alpha[3]-105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (141200537061855*alpha[7]-141200537061855*alpha[6]+141200537061855*alpha[5]-105244666443262*alpha[4]-141200537061855*alpha[3]+105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if ((-105244666443262*alpha[4])+105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_l(fEdge); 
  } 
  if ((-141200537061855*alpha[7])+141200537061855*alpha[6]-141200537061855*alpha[5]-105244666443262*alpha[4]+141200537061855*alpha[3]+105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if (141200537061855*alpha[7]+141200537061855*alpha[6]-141200537061855*alpha[5]-105244666443262*alpha[4]-141200537061855*alpha[3]-105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_l(fEdge); 
  } 
  if ((-105244666443262*alpha[4])-105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_l(fEdge); 
  } 
  if ((-141200537061855*alpha[7])-141200537061855*alpha[6]+141200537061855*alpha[5]-105244666443262*alpha[4]+141200537061855*alpha[3]-105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if ((-141200537061855*alpha[7])-141200537061855*alpha[6]-141200537061855*alpha[5]+105244666443262*alpha[4]-141200537061855*alpha[3]+105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_l(fEdge); 
  } 
  if (105244666443262*alpha[4]+105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_l(fEdge); 
  } 
  if (141200537061855*alpha[7]+141200537061855*alpha[6]+141200537061855*alpha[5]+105244666443262*alpha[4]+141200537061855*alpha[3]+105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.02357022603955158*(13.41640786499874*alpha[7]*fUpwind[11]+13.41640786499874*(alpha[6]*fUpwind[10]+alpha[5]*fUpwind[9])+13.41640786499874*alpha[3]*fUpwind[8]+15.0*(alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3])); 
  Ghat[4] = 0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] = 0.02357022603955158*(13.41640786499874*alpha[6]*fUpwind[11]+13.41640786499874*(alpha[7]*fUpwind[10]+alpha[3]*fUpwind[9])+13.41640786499874*alpha[5]*fUpwind[8]+15.0*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3])); 
  Ghat[6] = 0.02357022603955158*(13.41640786499874*alpha[5]*fUpwind[11]+13.41640786499874*(alpha[3]*fUpwind[10]+alpha[7]*fUpwind[9])+13.41640786499874*alpha[6]*fUpwind[8]+15.0*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3])); 
  Ghat[7] = 0.02357022603955158*(13.41640786499874*alpha[3]*fUpwind[11]+13.41640786499874*(alpha[5]*fUpwind[10]+alpha[6]*fUpwind[9])+13.41640786499874*alpha[7]*fUpwind[8]+15.0*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4])); 
  Ghat[8] = 0.02357022603955158*(15.0*alpha[4]*fUpwind[11]+15.0*(alpha[2]*fUpwind[10]+alpha[1]*fUpwind[9])+15.0*alpha[0]*fUpwind[8]+13.41640786499874*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[3]*fUpwind[3])); 
  Ghat[9] = 0.02357022603955158*(15.0*alpha[2]*fUpwind[11]+15.0*(alpha[4]*fUpwind[10]+alpha[0]*fUpwind[9])+15.0*alpha[1]*fUpwind[8]+13.41640786499874*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5])); 
  Ghat[10] = 0.02357022603955158*(15.0*alpha[1]*fUpwind[11]+15.0*(alpha[0]*fUpwind[10]+alpha[4]*fUpwind[9])+15.0*alpha[2]*fUpwind[8]+13.41640786499874*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])); 
  Ghat[11] = 0.02357022603955158*(15.0*alpha[0]*fUpwind[11]+15.0*(alpha[1]*fUpwind[10]+alpha[2]*fUpwind[9])+15.0*alpha[4]*fUpwind[8]+13.41640786499874*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6])); 

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
  out[16] += -1.58113883008419*Ghat[0]*dv10; 
  out[17] += -1.58113883008419*Ghat[1]*dv10; 
  out[18] += -1.58113883008419*Ghat[2]*dv10; 
  out[19] += -1.58113883008419*Ghat[3]*dv10; 
  out[20] += -1.58113883008419*Ghat[4]*dv10; 
  out[21] += -1.58113883008419*Ghat[5]*dv10; 
  out[22] += -1.58113883008419*Ghat[6]*dv10; 
  out[23] += -1.58113883008419*Ghat[7]*dv10; 
  out[24] += -0.7071067811865475*Ghat[8]*dv10; 
  out[25] += -0.7071067811865475*Ghat[9]*dv10; 
  out[26] += -0.7071067811865475*Ghat[10]*dv10; 
  out[27] += -1.224744871391589*Ghat[8]*dv10; 
  out[28] += -0.7071067811865475*Ghat[11]*dv10; 
  out[29] += -1.224744871391589*Ghat[9]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -1.224744871391589*Ghat[11]*dv10; 

  } else { 

  if ((-141200537061855*alpha[7])+141200537061855*alpha[6]+141200537061855*alpha[5]+105244666443262*alpha[4]-141200537061855*alpha[3]-105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if (105244666443262*alpha[4]-105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (141200537061855*alpha[7]-141200537061855*alpha[6]-141200537061855*alpha[5]+105244666443262*alpha[4]+141200537061855*alpha[3]-105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (141200537061855*alpha[7]-141200537061855*alpha[6]+141200537061855*alpha[5]-105244666443262*alpha[4]-141200537061855*alpha[3]+105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if ((-105244666443262*alpha[4])+105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_l(fSkin); 
  } 
  if ((-141200537061855*alpha[7])+141200537061855*alpha[6]-141200537061855*alpha[5]-105244666443262*alpha[4]+141200537061855*alpha[3]+105244666443262*alpha[2]-105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if (141200537061855*alpha[7]+141200537061855*alpha[6]-141200537061855*alpha[5]-105244666443262*alpha[4]-141200537061855*alpha[3]-105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_l(fSkin); 
  } 
  if ((-105244666443262*alpha[4])-105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_l(fSkin); 
  } 
  if ((-141200537061855*alpha[7])-141200537061855*alpha[6]+141200537061855*alpha[5]-105244666443262*alpha[4]+141200537061855*alpha[3]-105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if ((-141200537061855*alpha[7])-141200537061855*alpha[6]-141200537061855*alpha[5]+105244666443262*alpha[4]-141200537061855*alpha[3]+105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_l(fSkin); 
  } 
  if (105244666443262*alpha[4]+105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_l(fSkin); 
  } 
  if (141200537061855*alpha[7]+141200537061855*alpha[6]+141200537061855*alpha[5]+105244666443262*alpha[4]+141200537061855*alpha[3]+105244666443262*alpha[2]+105244666443262*alpha[1]+105244666443262*alpha[0] > 0) { 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.02357022603955158*(13.41640786499874*alpha[7]*fUpwind[11]+13.41640786499874*(alpha[6]*fUpwind[10]+alpha[5]*fUpwind[9])+13.41640786499874*alpha[3]*fUpwind[8]+15.0*(alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3])); 
  Ghat[4] = 0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] = 0.02357022603955158*(13.41640786499874*alpha[6]*fUpwind[11]+13.41640786499874*(alpha[7]*fUpwind[10]+alpha[3]*fUpwind[9])+13.41640786499874*alpha[5]*fUpwind[8]+15.0*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3])); 
  Ghat[6] = 0.02357022603955158*(13.41640786499874*alpha[5]*fUpwind[11]+13.41640786499874*(alpha[3]*fUpwind[10]+alpha[7]*fUpwind[9])+13.41640786499874*alpha[6]*fUpwind[8]+15.0*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3])); 
  Ghat[7] = 0.02357022603955158*(13.41640786499874*alpha[3]*fUpwind[11]+13.41640786499874*(alpha[5]*fUpwind[10]+alpha[6]*fUpwind[9])+13.41640786499874*alpha[7]*fUpwind[8]+15.0*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4])); 
  Ghat[8] = 0.02357022603955158*(15.0*alpha[4]*fUpwind[11]+15.0*(alpha[2]*fUpwind[10]+alpha[1]*fUpwind[9])+15.0*alpha[0]*fUpwind[8]+13.41640786499874*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[3]*fUpwind[3])); 
  Ghat[9] = 0.02357022603955158*(15.0*alpha[2]*fUpwind[11]+15.0*(alpha[4]*fUpwind[10]+alpha[0]*fUpwind[9])+15.0*alpha[1]*fUpwind[8]+13.41640786499874*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5])); 
  Ghat[10] = 0.02357022603955158*(15.0*alpha[1]*fUpwind[11]+15.0*(alpha[0]*fUpwind[10]+alpha[4]*fUpwind[9])+15.0*alpha[2]*fUpwind[8]+13.41640786499874*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])); 
  Ghat[11] = 0.02357022603955158*(15.0*alpha[0]*fUpwind[11]+15.0*(alpha[1]*fUpwind[10]+alpha[2]*fUpwind[9])+15.0*alpha[4]*fUpwind[8]+13.41640786499874*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6])); 

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
  out[16] += 1.58113883008419*Ghat[0]*dv10; 
  out[17] += 1.58113883008419*Ghat[1]*dv10; 
  out[18] += 1.58113883008419*Ghat[2]*dv10; 
  out[19] += 1.58113883008419*Ghat[3]*dv10; 
  out[20] += 1.58113883008419*Ghat[4]*dv10; 
  out[21] += 1.58113883008419*Ghat[5]*dv10; 
  out[22] += 1.58113883008419*Ghat[6]*dv10; 
  out[23] += 1.58113883008419*Ghat[7]*dv10; 
  out[24] += 0.7071067811865475*Ghat[8]*dv10; 
  out[25] += 0.7071067811865475*Ghat[9]*dv10; 
  out[26] += 0.7071067811865475*Ghat[10]*dv10; 
  out[27] += -1.224744871391589*Ghat[8]*dv10; 
  out[28] += 0.7071067811865475*Ghat[11]*dv10; 
  out[29] += -1.224744871391589*Ghat[9]*dv10; 
  out[30] += -1.224744871391589*Ghat[10]*dv10; 
  out[31] += -1.224744871391589*Ghat[11]*dv10; 

  } 
} 
