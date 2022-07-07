#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_5x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
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
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E0 = &qmem[0]; 
  const double *B0 = &qmem[12]; 
  const double *B1 = &qmem[16]; 
  const double *B2 = &qmem[20]; 

  double alpha[32] = {0.0}; 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 2.0*(B2[2]*wv2+E0[2])-2.0*B1[2]*wv3; 
  alpha[3] = 0.5773502691896258*B2[0]*dv2; 
  alpha[4] = -0.5773502691896258*B1[0]*dv3; 
  alpha[5] = 2.0*(B2[3]*wv2+E0[3])-2.0*B1[3]*wv3; 
  alpha[6] = 0.5773502691896258*B2[1]*dv2; 
  alpha[7] = 0.5773502691896258*B2[2]*dv2; 
  alpha[8] = -0.5773502691896258*B1[1]*dv3; 
  alpha[9] = -0.5773502691896258*B1[2]*dv3; 
  alpha[11] = 0.5773502691896258*B2[3]*dv2; 
  alpha[12] = -0.5773502691896258*B1[3]*dv3; 

  double fUpwindQuad[16] = {0.0};
  double fUpwind[32] = {0.0};
  double Ghat[32] = {0.0}; 

  if (edge == -1) { 

  if ((-alpha[12])-alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if (alpha[12]-alpha[11]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if ((-alpha[12])+alpha[11]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (alpha[12]+alpha[11]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if (alpha[12]+alpha[11]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p1_surfx3_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_5x_p1_surfx3_eval_quad_node_4_l(fEdge); 
  } 
  if ((-alpha[12])+alpha[11]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p1_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_5x_p1_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if (alpha[12]-alpha[11]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p1_surfx3_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_5x_p1_surfx3_eval_quad_node_6_l(fEdge); 
  } 
  if ((-alpha[12])-alpha[11]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p1_surfx3_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_5x_p1_surfx3_eval_quad_node_7_l(fEdge); 
  } 
  if (alpha[12]+alpha[11]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p1_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_5x_p1_surfx3_eval_quad_node_8_l(fEdge); 
  } 
  if ((-alpha[12])+alpha[11]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p1_surfx3_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_5x_p1_surfx3_eval_quad_node_9_l(fEdge); 
  } 
  if (alpha[12]-alpha[11]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p1_surfx3_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_5x_p1_surfx3_eval_quad_node_10_l(fEdge); 
  } 
  if ((-alpha[12])-alpha[11]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p1_surfx3_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_5x_p1_surfx3_eval_quad_node_11_l(fEdge); 
  } 
  if ((-alpha[12])-alpha[11]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p1_surfx3_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_5x_p1_surfx3_eval_quad_node_12_l(fEdge); 
  } 
  if (alpha[12]-alpha[11]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p1_surfx3_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_5x_p1_surfx3_eval_quad_node_13_l(fEdge); 
  } 
  if ((-alpha[12])+alpha[11]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p1_surfx3_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_5x_p1_surfx3_eval_quad_node_14_l(fEdge); 
  } 
  if (alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p1_surfx3_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_5x_p1_surfx3_eval_quad_node_15_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*(alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.25*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.01666666666666667*(13.41640786499874*alpha[11]*fUpwind[20]+13.41640786499874*(alpha[7]*fUpwind[18]+alpha[6]*fUpwind[17])+13.41640786499874*alpha[3]*fUpwind[16]+15.0*(alpha[12]*fUpwind[15]+alpha[9]*fUpwind[14]+alpha[8]*fUpwind[13]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]+alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3])); 
  Ghat[4] = 0.01666666666666667*(13.41640786499874*alpha[12]*fUpwind[28]+13.41640786499874*(alpha[9]*fUpwind[26]+alpha[8]*fUpwind[25])+13.41640786499874*alpha[4]*fUpwind[24]+15.0*(alpha[11]*fUpwind[15]+alpha[7]*fUpwind[14]+alpha[6]*fUpwind[13]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4])); 
  Ghat[5] = 0.25*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] = 0.01666666666666667*(13.41640786499874*alpha[7]*fUpwind[20]+13.41640786499874*(alpha[11]*fUpwind[18]+alpha[3]*fUpwind[17])+13.41640786499874*alpha[6]*fUpwind[16]+15.0*(alpha[9]*fUpwind[15]+alpha[12]*fUpwind[14]+alpha[4]*fUpwind[13]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[8]*fUpwind[10]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3])); 
  Ghat[7] = 0.01666666666666667*(13.41640786499874*alpha[6]*fUpwind[20]+13.41640786499874*(alpha[3]*fUpwind[18]+alpha[11]*fUpwind[17])+13.41640786499874*alpha[7]*fUpwind[16]+15.0*(alpha[8]*fUpwind[15]+alpha[4]*fUpwind[14]+alpha[12]*fUpwind[13]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[9]*fUpwind[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3])); 
  Ghat[8] = 0.01666666666666667*(13.41640786499874*alpha[9]*fUpwind[28]+13.41640786499874*(alpha[12]*fUpwind[26]+alpha[4]*fUpwind[25])+13.41640786499874*alpha[8]*fUpwind[24]+15.0*(alpha[7]*fUpwind[15]+alpha[11]*fUpwind[14]+alpha[3]*fUpwind[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[6]*fUpwind[10]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4])); 
  Ghat[9] = 0.01666666666666667*(13.41640786499874*alpha[8]*fUpwind[28]+13.41640786499874*(alpha[4]*fUpwind[26]+alpha[12]*fUpwind[25])+13.41640786499874*alpha[9]*fUpwind[24]+15.0*(alpha[6]*fUpwind[15]+alpha[3]*fUpwind[14]+alpha[11]*fUpwind[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4])); 
  Ghat[10] = 0.01666666666666667*(13.41640786499874*alpha[12]*fUpwind[31]+13.41640786499874*(alpha[9]*fUpwind[30]+alpha[8]*fUpwind[29])+13.41640786499874*(alpha[4]*fUpwind[27]+alpha[11]*fUpwind[23])+13.41640786499874*(alpha[7]*fUpwind[22]+alpha[6]*fUpwind[21])+13.41640786499874*alpha[3]*fUpwind[19]+15.0*(alpha[5]*fUpwind[15]+alpha[2]*fUpwind[14]+alpha[1]*fUpwind[13]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4])); 
  Ghat[11] = 0.01666666666666667*(13.41640786499874*alpha[3]*fUpwind[20]+13.41640786499874*(alpha[6]*fUpwind[18]+alpha[7]*fUpwind[17])+13.41640786499874*alpha[11]*fUpwind[16]+15.0*(alpha[4]*fUpwind[15]+alpha[8]*fUpwind[14]+alpha[9]*fUpwind[13]+fUpwind[10]*alpha[12]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5])); 
  Ghat[12] = 0.01666666666666667*(13.41640786499874*alpha[4]*fUpwind[28]+13.41640786499874*(alpha[8]*fUpwind[26]+alpha[9]*fUpwind[25])+13.41640786499874*alpha[12]*fUpwind[24]+15.0*(alpha[3]*fUpwind[15]+alpha[6]*fUpwind[14]+alpha[7]*fUpwind[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+fUpwind[10]*alpha[11]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5])); 
  Ghat[13] = 0.01666666666666667*(13.41640786499874*alpha[9]*fUpwind[31]+13.41640786499874*(alpha[12]*fUpwind[30]+alpha[4]*fUpwind[29])+13.41640786499874*(alpha[8]*fUpwind[27]+alpha[7]*fUpwind[23])+13.41640786499874*(alpha[11]*fUpwind[22]+alpha[3]*fUpwind[21])+13.41640786499874*alpha[6]*fUpwind[19]+15.0*(alpha[2]*fUpwind[15]+alpha[5]*fUpwind[14]+alpha[0]*fUpwind[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6])); 
  Ghat[14] = 0.01666666666666667*(13.41640786499874*alpha[8]*fUpwind[31]+13.41640786499874*(alpha[4]*fUpwind[30]+alpha[12]*fUpwind[29])+13.41640786499874*(alpha[9]*fUpwind[27]+alpha[6]*fUpwind[23])+13.41640786499874*(alpha[3]*fUpwind[22]+alpha[11]*fUpwind[21])+13.41640786499874*alpha[7]*fUpwind[19]+15.0*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]+alpha[5]*fUpwind[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7])); 
  Ghat[15] = 0.01666666666666667*(13.41640786499874*alpha[4]*fUpwind[31]+13.41640786499874*(alpha[8]*fUpwind[30]+alpha[9]*fUpwind[29])+13.41640786499874*(alpha[12]*fUpwind[27]+alpha[3]*fUpwind[23])+13.41640786499874*(alpha[6]*fUpwind[22]+alpha[7]*fUpwind[21])+13.41640786499874*alpha[11]*fUpwind[19]+15.0*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[5]*fUpwind[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8])); 
  Ghat[16] = 0.01666666666666667*(15.0*alpha[12]*fUpwind[23]+15.0*(alpha[9]*fUpwind[22]+alpha[8]*fUpwind[21]+alpha[5]*fUpwind[20])+15.0*(alpha[4]*fUpwind[19]+alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17])+15.0*alpha[0]*fUpwind[16]+13.41640786499874*(alpha[11]*fUpwind[11]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[3]*fUpwind[3])); 
  Ghat[17] = 0.01666666666666667*(15.0*alpha[9]*fUpwind[23]+15.0*(alpha[12]*fUpwind[22]+alpha[4]*fUpwind[21]+alpha[2]*fUpwind[20])+15.0*(alpha[8]*fUpwind[19]+alpha[5]*fUpwind[18]+alpha[0]*fUpwind[17])+15.0*alpha[1]*fUpwind[16]+13.41640786499874*(alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])); 
  Ghat[18] = 0.01666666666666667*(15.0*alpha[8]*fUpwind[23]+15.0*(alpha[4]*fUpwind[22]+alpha[12]*fUpwind[21]+alpha[1]*fUpwind[20])+15.0*(alpha[9]*fUpwind[19]+alpha[0]*fUpwind[18]+alpha[5]*fUpwind[17])+15.0*alpha[2]*fUpwind[16]+13.41640786499874*(alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7])); 
  Ghat[19] = 0.01666666666666667*(15.0*alpha[5]*fUpwind[23]+15.0*(alpha[2]*fUpwind[22]+alpha[1]*fUpwind[21]+alpha[12]*fUpwind[20])+15.0*(alpha[0]*fUpwind[19]+alpha[9]*fUpwind[18]+alpha[8]*fUpwind[17])+15.0*alpha[4]*fUpwind[16]+13.41640786499874*(alpha[11]*fUpwind[15]+alpha[7]*fUpwind[14]+alpha[6]*fUpwind[13]+alpha[3]*fUpwind[10])); 
  Ghat[20] = 0.01666666666666667*(15.0*alpha[4]*fUpwind[23]+15.0*(alpha[8]*fUpwind[22]+alpha[9]*fUpwind[21]+alpha[0]*fUpwind[20])+15.0*(alpha[12]*fUpwind[19]+alpha[1]*fUpwind[18]+alpha[2]*fUpwind[17])+15.0*alpha[5]*fUpwind[16]+13.41640786499874*(alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7])); 
  Ghat[21] = 0.01666666666666667*(15.0*alpha[2]*fUpwind[23]+15.0*(alpha[5]*fUpwind[22]+alpha[0]*fUpwind[21]+alpha[9]*fUpwind[20])+15.0*(alpha[1]*fUpwind[19]+alpha[12]*fUpwind[18]+alpha[4]*fUpwind[17])+15.0*alpha[8]*fUpwind[16]+13.41640786499874*(alpha[7]*fUpwind[15]+alpha[11]*fUpwind[14]+alpha[3]*fUpwind[13]+alpha[6]*fUpwind[10])); 
  Ghat[22] = 0.01666666666666667*(15.0*alpha[1]*fUpwind[23]+15.0*(alpha[0]*fUpwind[22]+alpha[5]*fUpwind[21]+alpha[8]*fUpwind[20])+15.0*(alpha[2]*fUpwind[19]+alpha[4]*fUpwind[18]+alpha[12]*fUpwind[17])+15.0*alpha[9]*fUpwind[16]+13.41640786499874*(alpha[6]*fUpwind[15]+alpha[3]*fUpwind[14]+alpha[11]*fUpwind[13]+alpha[7]*fUpwind[10])); 
  Ghat[23] = 0.01666666666666667*(15.0*alpha[0]*fUpwind[23]+15.0*(alpha[1]*fUpwind[22]+alpha[2]*fUpwind[21]+alpha[4]*fUpwind[20])+15.0*(alpha[5]*fUpwind[19]+alpha[8]*fUpwind[18]+alpha[9]*fUpwind[17])+15.0*alpha[12]*fUpwind[16]+13.41640786499874*(alpha[3]*fUpwind[15]+alpha[6]*fUpwind[14]+alpha[7]*fUpwind[13]+fUpwind[10]*alpha[11])); 
  Ghat[24] = 0.01666666666666667*(15.0*alpha[11]*fUpwind[31]+15.0*(alpha[7]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[5]*fUpwind[28])+15.0*(alpha[3]*fUpwind[27]+alpha[2]*fUpwind[26]+alpha[1]*fUpwind[25])+15.0*alpha[0]*fUpwind[24]+13.41640786499874*(alpha[12]*fUpwind[12]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[4]*fUpwind[4])); 
  Ghat[25] = 0.01666666666666667*(15.0*alpha[7]*fUpwind[31]+15.0*(alpha[11]*fUpwind[30]+alpha[3]*fUpwind[29]+alpha[2]*fUpwind[28])+15.0*(alpha[6]*fUpwind[27]+alpha[5]*fUpwind[26]+alpha[0]*fUpwind[25])+15.0*alpha[1]*fUpwind[24]+13.41640786499874*(alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8])); 
  Ghat[26] = 0.01666666666666667*(15.0*alpha[6]*fUpwind[31]+15.0*(alpha[3]*fUpwind[30]+alpha[11]*fUpwind[29]+alpha[1]*fUpwind[28])+15.0*(alpha[7]*fUpwind[27]+alpha[0]*fUpwind[26]+alpha[5]*fUpwind[25])+15.0*alpha[2]*fUpwind[24]+13.41640786499874*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9])); 
  Ghat[27] = 0.01666666666666667*(15.0*alpha[5]*fUpwind[31]+15.0*(alpha[2]*fUpwind[30]+alpha[1]*fUpwind[29]+alpha[11]*fUpwind[28])+15.0*(alpha[0]*fUpwind[27]+alpha[7]*fUpwind[26]+alpha[6]*fUpwind[25])+15.0*alpha[3]*fUpwind[24]+13.41640786499874*(alpha[12]*fUpwind[15]+alpha[9]*fUpwind[14]+alpha[8]*fUpwind[13]+alpha[4]*fUpwind[10])); 
  Ghat[28] = 0.01666666666666667*(15.0*alpha[3]*fUpwind[31]+15.0*(alpha[6]*fUpwind[30]+alpha[7]*fUpwind[29]+alpha[0]*fUpwind[28])+15.0*(alpha[11]*fUpwind[27]+alpha[1]*fUpwind[26]+alpha[2]*fUpwind[25])+15.0*alpha[5]*fUpwind[24]+13.41640786499874*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9])); 
  Ghat[29] = 0.01666666666666667*(15.0*alpha[2]*fUpwind[31]+15.0*(alpha[5]*fUpwind[30]+alpha[0]*fUpwind[29]+alpha[7]*fUpwind[28])+15.0*(alpha[1]*fUpwind[27]+alpha[11]*fUpwind[26]+alpha[3]*fUpwind[25])+15.0*alpha[6]*fUpwind[24]+13.41640786499874*(alpha[9]*fUpwind[15]+alpha[12]*fUpwind[14]+alpha[4]*fUpwind[13]+alpha[8]*fUpwind[10])); 
  Ghat[30] = 0.01666666666666667*(15.0*alpha[1]*fUpwind[31]+15.0*(alpha[0]*fUpwind[30]+alpha[5]*fUpwind[29]+alpha[6]*fUpwind[28])+15.0*(alpha[2]*fUpwind[27]+alpha[3]*fUpwind[26]+alpha[11]*fUpwind[25])+15.0*alpha[7]*fUpwind[24]+13.41640786499874*(alpha[8]*fUpwind[15]+alpha[4]*fUpwind[14]+alpha[12]*fUpwind[13]+alpha[9]*fUpwind[10])); 
  Ghat[31] = 0.01666666666666667*(15.0*alpha[0]*fUpwind[31]+15.0*(alpha[1]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[3]*fUpwind[28])+15.0*(alpha[5]*fUpwind[27]+alpha[6]*fUpwind[26]+alpha[7]*fUpwind[25])+15.0*alpha[11]*fUpwind[24]+13.41640786499874*(alpha[4]*fUpwind[15]+alpha[8]*fUpwind[14]+alpha[9]*fUpwind[13]+fUpwind[10]*alpha[12])); 

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
  out[32] += -1.58113883008419*Ghat[0]*dv10; 
  out[33] += -1.58113883008419*Ghat[1]*dv10; 
  out[34] += -1.58113883008419*Ghat[2]*dv10; 
  out[35] += -1.58113883008419*Ghat[3]*dv10; 
  out[36] += -1.58113883008419*Ghat[4]*dv10; 
  out[37] += -1.58113883008419*Ghat[5]*dv10; 
  out[38] += -1.58113883008419*Ghat[6]*dv10; 
  out[39] += -1.58113883008419*Ghat[7]*dv10; 
  out[40] += -1.58113883008419*Ghat[8]*dv10; 
  out[41] += -1.58113883008419*Ghat[9]*dv10; 
  out[42] += -1.58113883008419*Ghat[10]*dv10; 
  out[43] += -1.58113883008419*Ghat[11]*dv10; 
  out[44] += -1.58113883008419*Ghat[12]*dv10; 
  out[45] += -1.58113883008419*Ghat[13]*dv10; 
  out[46] += -1.58113883008419*Ghat[14]*dv10; 
  out[47] += -1.58113883008419*Ghat[15]*dv10; 
  out[48] += -0.7071067811865475*Ghat[16]*dv10; 
  out[49] += -0.7071067811865475*Ghat[17]*dv10; 
  out[50] += -0.7071067811865475*Ghat[18]*dv10; 
  out[51] += -1.224744871391589*Ghat[16]*dv10; 
  out[52] += -0.7071067811865475*Ghat[19]*dv10; 
  out[53] += -0.7071067811865475*Ghat[20]*dv10; 
  out[54] += -1.224744871391589*Ghat[17]*dv10; 
  out[55] += -1.224744871391589*Ghat[18]*dv10; 
  out[56] += -0.7071067811865475*Ghat[21]*dv10; 
  out[57] += -0.7071067811865475*Ghat[22]*dv10; 
  out[58] += -1.224744871391589*Ghat[19]*dv10; 
  out[59] += -1.224744871391589*Ghat[20]*dv10; 
  out[60] += -0.7071067811865475*Ghat[23]*dv10; 
  out[61] += -1.224744871391589*Ghat[21]*dv10; 
  out[62] += -1.224744871391589*Ghat[22]*dv10; 
  out[63] += -1.224744871391589*Ghat[23]*dv10; 
  out[64] += -0.7071067811865475*Ghat[24]*dv10; 
  out[65] += -0.7071067811865475*Ghat[25]*dv10; 
  out[66] += -0.7071067811865475*Ghat[26]*dv10; 
  out[67] += -1.224744871391589*Ghat[24]*dv10; 
  out[68] += -0.7071067811865475*Ghat[27]*dv10; 
  out[69] += -0.7071067811865475*Ghat[28]*dv10; 
  out[70] += -1.224744871391589*Ghat[25]*dv10; 
  out[71] += -1.224744871391589*Ghat[26]*dv10; 
  out[72] += -0.7071067811865475*Ghat[29]*dv10; 
  out[73] += -0.7071067811865475*Ghat[30]*dv10; 
  out[74] += -1.224744871391589*Ghat[27]*dv10; 
  out[75] += -1.224744871391589*Ghat[28]*dv10; 
  out[76] += -0.7071067811865475*Ghat[31]*dv10; 
  out[77] += -1.224744871391589*Ghat[29]*dv10; 
  out[78] += -1.224744871391589*Ghat[30]*dv10; 
  out[79] += -1.224744871391589*Ghat[31]*dv10; 

  } else { 

  if ((-alpha[12])-alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if (alpha[12]-alpha[11]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if ((-alpha[12])+alpha[11]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (alpha[12]+alpha[11]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if (alpha[12]+alpha[11]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p1_surfx3_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_5x_p1_surfx3_eval_quad_node_4_l(fSkin); 
  } 
  if ((-alpha[12])+alpha[11]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p1_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_5x_p1_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if (alpha[12]-alpha[11]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p1_surfx3_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_5x_p1_surfx3_eval_quad_node_6_l(fSkin); 
  } 
  if ((-alpha[12])-alpha[11]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p1_surfx3_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_5x_p1_surfx3_eval_quad_node_7_l(fSkin); 
  } 
  if (alpha[12]+alpha[11]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p1_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_5x_p1_surfx3_eval_quad_node_8_l(fSkin); 
  } 
  if ((-alpha[12])+alpha[11]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p1_surfx3_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_5x_p1_surfx3_eval_quad_node_9_l(fSkin); 
  } 
  if (alpha[12]-alpha[11]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p1_surfx3_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_5x_p1_surfx3_eval_quad_node_10_l(fSkin); 
  } 
  if ((-alpha[12])-alpha[11]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p1_surfx3_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_5x_p1_surfx3_eval_quad_node_11_l(fSkin); 
  } 
  if ((-alpha[12])-alpha[11]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p1_surfx3_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_5x_p1_surfx3_eval_quad_node_12_l(fSkin); 
  } 
  if (alpha[12]-alpha[11]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p1_surfx3_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_5x_p1_surfx3_eval_quad_node_13_l(fSkin); 
  } 
  if ((-alpha[12])+alpha[11]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p1_surfx3_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_5x_p1_surfx3_eval_quad_node_14_l(fSkin); 
  } 
  if (alpha[12]+alpha[11]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p1_surfx3_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_5x_p1_surfx3_eval_quad_node_15_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*(alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.25*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.01666666666666667*(13.41640786499874*alpha[11]*fUpwind[20]+13.41640786499874*(alpha[7]*fUpwind[18]+alpha[6]*fUpwind[17])+13.41640786499874*alpha[3]*fUpwind[16]+15.0*(alpha[12]*fUpwind[15]+alpha[9]*fUpwind[14]+alpha[8]*fUpwind[13]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]+alpha[4]*fUpwind[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3])); 
  Ghat[4] = 0.01666666666666667*(13.41640786499874*alpha[12]*fUpwind[28]+13.41640786499874*(alpha[9]*fUpwind[26]+alpha[8]*fUpwind[25])+13.41640786499874*alpha[4]*fUpwind[24]+15.0*(alpha[11]*fUpwind[15]+alpha[7]*fUpwind[14]+alpha[6]*fUpwind[13]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[3]*fUpwind[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4])); 
  Ghat[5] = 0.25*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] = 0.01666666666666667*(13.41640786499874*alpha[7]*fUpwind[20]+13.41640786499874*(alpha[11]*fUpwind[18]+alpha[3]*fUpwind[17])+13.41640786499874*alpha[6]*fUpwind[16]+15.0*(alpha[9]*fUpwind[15]+alpha[12]*fUpwind[14]+alpha[4]*fUpwind[13]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[8]*fUpwind[10]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3])); 
  Ghat[7] = 0.01666666666666667*(13.41640786499874*alpha[6]*fUpwind[20]+13.41640786499874*(alpha[3]*fUpwind[18]+alpha[11]*fUpwind[17])+13.41640786499874*alpha[7]*fUpwind[16]+15.0*(alpha[8]*fUpwind[15]+alpha[4]*fUpwind[14]+alpha[12]*fUpwind[13]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[9]*fUpwind[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3])); 
  Ghat[8] = 0.01666666666666667*(13.41640786499874*alpha[9]*fUpwind[28]+13.41640786499874*(alpha[12]*fUpwind[26]+alpha[4]*fUpwind[25])+13.41640786499874*alpha[8]*fUpwind[24]+15.0*(alpha[7]*fUpwind[15]+alpha[11]*fUpwind[14]+alpha[3]*fUpwind[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[6]*fUpwind[10]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4])); 
  Ghat[9] = 0.01666666666666667*(13.41640786499874*alpha[8]*fUpwind[28]+13.41640786499874*(alpha[4]*fUpwind[26]+alpha[12]*fUpwind[25])+13.41640786499874*alpha[9]*fUpwind[24]+15.0*(alpha[6]*fUpwind[15]+alpha[3]*fUpwind[14]+alpha[11]*fUpwind[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[7]*fUpwind[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4])); 
  Ghat[10] = 0.01666666666666667*(13.41640786499874*alpha[12]*fUpwind[31]+13.41640786499874*(alpha[9]*fUpwind[30]+alpha[8]*fUpwind[29])+13.41640786499874*(alpha[4]*fUpwind[27]+alpha[11]*fUpwind[23])+13.41640786499874*(alpha[7]*fUpwind[22]+alpha[6]*fUpwind[21])+13.41640786499874*alpha[3]*fUpwind[19]+15.0*(alpha[5]*fUpwind[15]+alpha[2]*fUpwind[14]+alpha[1]*fUpwind[13]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]+alpha[0]*fUpwind[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4])); 
  Ghat[11] = 0.01666666666666667*(13.41640786499874*alpha[3]*fUpwind[20]+13.41640786499874*(alpha[6]*fUpwind[18]+alpha[7]*fUpwind[17])+13.41640786499874*alpha[11]*fUpwind[16]+15.0*(alpha[4]*fUpwind[15]+alpha[8]*fUpwind[14]+alpha[9]*fUpwind[13]+fUpwind[10]*alpha[12]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5])); 
  Ghat[12] = 0.01666666666666667*(13.41640786499874*alpha[4]*fUpwind[28]+13.41640786499874*(alpha[8]*fUpwind[26]+alpha[9]*fUpwind[25])+13.41640786499874*alpha[12]*fUpwind[24]+15.0*(alpha[3]*fUpwind[15]+alpha[6]*fUpwind[14]+alpha[7]*fUpwind[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+fUpwind[10]*alpha[11]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5])); 
  Ghat[13] = 0.01666666666666667*(13.41640786499874*alpha[9]*fUpwind[31]+13.41640786499874*(alpha[12]*fUpwind[30]+alpha[4]*fUpwind[29])+13.41640786499874*(alpha[8]*fUpwind[27]+alpha[7]*fUpwind[23])+13.41640786499874*(alpha[11]*fUpwind[22]+alpha[3]*fUpwind[21])+13.41640786499874*alpha[6]*fUpwind[19]+15.0*(alpha[2]*fUpwind[15]+alpha[5]*fUpwind[14]+alpha[0]*fUpwind[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[1]*fUpwind[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6])); 
  Ghat[14] = 0.01666666666666667*(13.41640786499874*alpha[8]*fUpwind[31]+13.41640786499874*(alpha[4]*fUpwind[30]+alpha[12]*fUpwind[29])+13.41640786499874*(alpha[9]*fUpwind[27]+alpha[6]*fUpwind[23])+13.41640786499874*(alpha[3]*fUpwind[22]+alpha[11]*fUpwind[21])+13.41640786499874*alpha[7]*fUpwind[19]+15.0*(alpha[1]*fUpwind[15]+alpha[0]*fUpwind[14]+alpha[5]*fUpwind[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[2]*fUpwind[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7])); 
  Ghat[15] = 0.01666666666666667*(13.41640786499874*alpha[4]*fUpwind[31]+13.41640786499874*(alpha[8]*fUpwind[30]+alpha[9]*fUpwind[29])+13.41640786499874*(alpha[12]*fUpwind[27]+alpha[3]*fUpwind[23])+13.41640786499874*(alpha[6]*fUpwind[22]+alpha[7]*fUpwind[21])+13.41640786499874*alpha[11]*fUpwind[19]+15.0*(alpha[0]*fUpwind[15]+alpha[1]*fUpwind[14]+alpha[2]*fUpwind[13]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[5]*fUpwind[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8])); 
  Ghat[16] = 0.01666666666666667*(15.0*alpha[12]*fUpwind[23]+15.0*(alpha[9]*fUpwind[22]+alpha[8]*fUpwind[21]+alpha[5]*fUpwind[20])+15.0*(alpha[4]*fUpwind[19]+alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17])+15.0*alpha[0]*fUpwind[16]+13.41640786499874*(alpha[11]*fUpwind[11]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[3]*fUpwind[3])); 
  Ghat[17] = 0.01666666666666667*(15.0*alpha[9]*fUpwind[23]+15.0*(alpha[12]*fUpwind[22]+alpha[4]*fUpwind[21]+alpha[2]*fUpwind[20])+15.0*(alpha[8]*fUpwind[19]+alpha[5]*fUpwind[18]+alpha[0]*fUpwind[17])+15.0*alpha[1]*fUpwind[16]+13.41640786499874*(alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6])); 
  Ghat[18] = 0.01666666666666667*(15.0*alpha[8]*fUpwind[23]+15.0*(alpha[4]*fUpwind[22]+alpha[12]*fUpwind[21]+alpha[1]*fUpwind[20])+15.0*(alpha[9]*fUpwind[19]+alpha[0]*fUpwind[18]+alpha[5]*fUpwind[17])+15.0*alpha[2]*fUpwind[16]+13.41640786499874*(alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7])); 
  Ghat[19] = 0.01666666666666667*(15.0*alpha[5]*fUpwind[23]+15.0*(alpha[2]*fUpwind[22]+alpha[1]*fUpwind[21]+alpha[12]*fUpwind[20])+15.0*(alpha[0]*fUpwind[19]+alpha[9]*fUpwind[18]+alpha[8]*fUpwind[17])+15.0*alpha[4]*fUpwind[16]+13.41640786499874*(alpha[11]*fUpwind[15]+alpha[7]*fUpwind[14]+alpha[6]*fUpwind[13]+alpha[3]*fUpwind[10])); 
  Ghat[20] = 0.01666666666666667*(15.0*alpha[4]*fUpwind[23]+15.0*(alpha[8]*fUpwind[22]+alpha[9]*fUpwind[21]+alpha[0]*fUpwind[20])+15.0*(alpha[12]*fUpwind[19]+alpha[1]*fUpwind[18]+alpha[2]*fUpwind[17])+15.0*alpha[5]*fUpwind[16]+13.41640786499874*(alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7])); 
  Ghat[21] = 0.01666666666666667*(15.0*alpha[2]*fUpwind[23]+15.0*(alpha[5]*fUpwind[22]+alpha[0]*fUpwind[21]+alpha[9]*fUpwind[20])+15.0*(alpha[1]*fUpwind[19]+alpha[12]*fUpwind[18]+alpha[4]*fUpwind[17])+15.0*alpha[8]*fUpwind[16]+13.41640786499874*(alpha[7]*fUpwind[15]+alpha[11]*fUpwind[14]+alpha[3]*fUpwind[13]+alpha[6]*fUpwind[10])); 
  Ghat[22] = 0.01666666666666667*(15.0*alpha[1]*fUpwind[23]+15.0*(alpha[0]*fUpwind[22]+alpha[5]*fUpwind[21]+alpha[8]*fUpwind[20])+15.0*(alpha[2]*fUpwind[19]+alpha[4]*fUpwind[18]+alpha[12]*fUpwind[17])+15.0*alpha[9]*fUpwind[16]+13.41640786499874*(alpha[6]*fUpwind[15]+alpha[3]*fUpwind[14]+alpha[11]*fUpwind[13]+alpha[7]*fUpwind[10])); 
  Ghat[23] = 0.01666666666666667*(15.0*alpha[0]*fUpwind[23]+15.0*(alpha[1]*fUpwind[22]+alpha[2]*fUpwind[21]+alpha[4]*fUpwind[20])+15.0*(alpha[5]*fUpwind[19]+alpha[8]*fUpwind[18]+alpha[9]*fUpwind[17])+15.0*alpha[12]*fUpwind[16]+13.41640786499874*(alpha[3]*fUpwind[15]+alpha[6]*fUpwind[14]+alpha[7]*fUpwind[13]+fUpwind[10]*alpha[11])); 
  Ghat[24] = 0.01666666666666667*(15.0*alpha[11]*fUpwind[31]+15.0*(alpha[7]*fUpwind[30]+alpha[6]*fUpwind[29]+alpha[5]*fUpwind[28])+15.0*(alpha[3]*fUpwind[27]+alpha[2]*fUpwind[26]+alpha[1]*fUpwind[25])+15.0*alpha[0]*fUpwind[24]+13.41640786499874*(alpha[12]*fUpwind[12]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[4]*fUpwind[4])); 
  Ghat[25] = 0.01666666666666667*(15.0*alpha[7]*fUpwind[31]+15.0*(alpha[11]*fUpwind[30]+alpha[3]*fUpwind[29]+alpha[2]*fUpwind[28])+15.0*(alpha[6]*fUpwind[27]+alpha[5]*fUpwind[26]+alpha[0]*fUpwind[25])+15.0*alpha[1]*fUpwind[24]+13.41640786499874*(alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8])); 
  Ghat[26] = 0.01666666666666667*(15.0*alpha[6]*fUpwind[31]+15.0*(alpha[3]*fUpwind[30]+alpha[11]*fUpwind[29]+alpha[1]*fUpwind[28])+15.0*(alpha[7]*fUpwind[27]+alpha[0]*fUpwind[26]+alpha[5]*fUpwind[25])+15.0*alpha[2]*fUpwind[24]+13.41640786499874*(alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9])); 
  Ghat[27] = 0.01666666666666667*(15.0*alpha[5]*fUpwind[31]+15.0*(alpha[2]*fUpwind[30]+alpha[1]*fUpwind[29]+alpha[11]*fUpwind[28])+15.0*(alpha[0]*fUpwind[27]+alpha[7]*fUpwind[26]+alpha[6]*fUpwind[25])+15.0*alpha[3]*fUpwind[24]+13.41640786499874*(alpha[12]*fUpwind[15]+alpha[9]*fUpwind[14]+alpha[8]*fUpwind[13]+alpha[4]*fUpwind[10])); 
  Ghat[28] = 0.01666666666666667*(15.0*alpha[3]*fUpwind[31]+15.0*(alpha[6]*fUpwind[30]+alpha[7]*fUpwind[29]+alpha[0]*fUpwind[28])+15.0*(alpha[11]*fUpwind[27]+alpha[1]*fUpwind[26]+alpha[2]*fUpwind[25])+15.0*alpha[5]*fUpwind[24]+13.41640786499874*(alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9])); 
  Ghat[29] = 0.01666666666666667*(15.0*alpha[2]*fUpwind[31]+15.0*(alpha[5]*fUpwind[30]+alpha[0]*fUpwind[29]+alpha[7]*fUpwind[28])+15.0*(alpha[1]*fUpwind[27]+alpha[11]*fUpwind[26]+alpha[3]*fUpwind[25])+15.0*alpha[6]*fUpwind[24]+13.41640786499874*(alpha[9]*fUpwind[15]+alpha[12]*fUpwind[14]+alpha[4]*fUpwind[13]+alpha[8]*fUpwind[10])); 
  Ghat[30] = 0.01666666666666667*(15.0*alpha[1]*fUpwind[31]+15.0*(alpha[0]*fUpwind[30]+alpha[5]*fUpwind[29]+alpha[6]*fUpwind[28])+15.0*(alpha[2]*fUpwind[27]+alpha[3]*fUpwind[26]+alpha[11]*fUpwind[25])+15.0*alpha[7]*fUpwind[24]+13.41640786499874*(alpha[8]*fUpwind[15]+alpha[4]*fUpwind[14]+alpha[12]*fUpwind[13]+alpha[9]*fUpwind[10])); 
  Ghat[31] = 0.01666666666666667*(15.0*alpha[0]*fUpwind[31]+15.0*(alpha[1]*fUpwind[30]+alpha[2]*fUpwind[29]+alpha[3]*fUpwind[28])+15.0*(alpha[5]*fUpwind[27]+alpha[6]*fUpwind[26]+alpha[7]*fUpwind[25])+15.0*alpha[11]*fUpwind[24]+13.41640786499874*(alpha[4]*fUpwind[15]+alpha[8]*fUpwind[14]+alpha[9]*fUpwind[13]+fUpwind[10]*alpha[12])); 

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
  out[32] += 1.58113883008419*Ghat[0]*dv10; 
  out[33] += 1.58113883008419*Ghat[1]*dv10; 
  out[34] += 1.58113883008419*Ghat[2]*dv10; 
  out[35] += 1.58113883008419*Ghat[3]*dv10; 
  out[36] += 1.58113883008419*Ghat[4]*dv10; 
  out[37] += 1.58113883008419*Ghat[5]*dv10; 
  out[38] += 1.58113883008419*Ghat[6]*dv10; 
  out[39] += 1.58113883008419*Ghat[7]*dv10; 
  out[40] += 1.58113883008419*Ghat[8]*dv10; 
  out[41] += 1.58113883008419*Ghat[9]*dv10; 
  out[42] += 1.58113883008419*Ghat[10]*dv10; 
  out[43] += 1.58113883008419*Ghat[11]*dv10; 
  out[44] += 1.58113883008419*Ghat[12]*dv10; 
  out[45] += 1.58113883008419*Ghat[13]*dv10; 
  out[46] += 1.58113883008419*Ghat[14]*dv10; 
  out[47] += 1.58113883008419*Ghat[15]*dv10; 
  out[48] += 0.7071067811865475*Ghat[16]*dv10; 
  out[49] += 0.7071067811865475*Ghat[17]*dv10; 
  out[50] += 0.7071067811865475*Ghat[18]*dv10; 
  out[51] += -1.224744871391589*Ghat[16]*dv10; 
  out[52] += 0.7071067811865475*Ghat[19]*dv10; 
  out[53] += 0.7071067811865475*Ghat[20]*dv10; 
  out[54] += -1.224744871391589*Ghat[17]*dv10; 
  out[55] += -1.224744871391589*Ghat[18]*dv10; 
  out[56] += 0.7071067811865475*Ghat[21]*dv10; 
  out[57] += 0.7071067811865475*Ghat[22]*dv10; 
  out[58] += -1.224744871391589*Ghat[19]*dv10; 
  out[59] += -1.224744871391589*Ghat[20]*dv10; 
  out[60] += 0.7071067811865475*Ghat[23]*dv10; 
  out[61] += -1.224744871391589*Ghat[21]*dv10; 
  out[62] += -1.224744871391589*Ghat[22]*dv10; 
  out[63] += -1.224744871391589*Ghat[23]*dv10; 
  out[64] += 0.7071067811865475*Ghat[24]*dv10; 
  out[65] += 0.7071067811865475*Ghat[25]*dv10; 
  out[66] += 0.7071067811865475*Ghat[26]*dv10; 
  out[67] += -1.224744871391589*Ghat[24]*dv10; 
  out[68] += 0.7071067811865475*Ghat[27]*dv10; 
  out[69] += 0.7071067811865475*Ghat[28]*dv10; 
  out[70] += -1.224744871391589*Ghat[25]*dv10; 
  out[71] += -1.224744871391589*Ghat[26]*dv10; 
  out[72] += 0.7071067811865475*Ghat[29]*dv10; 
  out[73] += 0.7071067811865475*Ghat[30]*dv10; 
  out[74] += -1.224744871391589*Ghat[27]*dv10; 
  out[75] += -1.224744871391589*Ghat[28]*dv10; 
  out[76] += 0.7071067811865475*Ghat[31]*dv10; 
  out[77] += -1.224744871391589*Ghat[29]*dv10; 
  out[78] += -1.224744871391589*Ghat[30]*dv10; 
  out[79] += -1.224744871391589*Ghat[31]*dv10; 

  } 
} 
