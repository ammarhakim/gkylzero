#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_ser_5x_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_5x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_boundary_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *E1 = &qmem[4]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *B0 = &qmem[12]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *B1 = &qmem[16]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
  const double *B2 = &qmem[20]; 

  double alpha[16] = {0.0}; 

  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};
  double Ghat[16] = {0.0}; 

  if (edge == -1) { 

  alpha[0] = 1.224744871391589*B0[0]*p2_over_gamma[2]-1.224744871391589*B2[0]*p0_over_gamma[2]+0.7071067811865475*B0[0]*p2_over_gamma[0]-0.7071067811865475*B2[0]*p0_over_gamma[0]+2.0*E1[0]; 
  alpha[1] = 1.224744871391589*B0[1]*p2_over_gamma[2]-1.224744871391589*B2[1]*p0_over_gamma[2]+2.0*E1[1]-0.7071067811865475*p0_over_gamma[0]*B2[1]+0.7071067811865475*p2_over_gamma[0]*B0[1]; 
  alpha[2] = 1.224744871391589*B0[2]*p2_over_gamma[2]-1.224744871391589*B2[2]*p0_over_gamma[2]+2.0*E1[2]-0.7071067811865475*p0_over_gamma[0]*B2[2]+0.7071067811865475*p2_over_gamma[0]*B0[2]; 
  alpha[3] = 1.224744871391589*B0[0]*p2_over_gamma[4]-1.224744871391589*B2[0]*p0_over_gamma[4]+0.7071067811865475*B0[0]*p2_over_gamma[1]-0.7071067811865475*B2[0]*p0_over_gamma[1]; 
  alpha[4] = 1.224744871391589*B0[0]*p2_over_gamma[6]-1.224744871391589*B2[0]*p0_over_gamma[6]+0.7071067811865475*B0[0]*p2_over_gamma[3]-0.7071067811865475*B2[0]*p0_over_gamma[3]; 
  alpha[5] = 2.0*E1[3]-1.224744871391589*p0_over_gamma[2]*B2[3]-0.7071067811865475*p0_over_gamma[0]*B2[3]+1.224744871391589*p2_over_gamma[2]*B0[3]+0.7071067811865475*p2_over_gamma[0]*B0[3]; 
  alpha[6] = 1.224744871391589*B0[1]*p2_over_gamma[4]-1.224744871391589*B2[1]*p0_over_gamma[4]+0.7071067811865475*B0[1]*p2_over_gamma[1]-0.7071067811865475*B2[1]*p0_over_gamma[1]; 
  alpha[7] = 1.224744871391589*B0[2]*p2_over_gamma[4]-1.224744871391589*B2[2]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[2]+0.7071067811865475*p2_over_gamma[1]*B0[2]; 
  alpha[8] = 1.224744871391589*B0[1]*p2_over_gamma[6]-1.224744871391589*B2[1]*p0_over_gamma[6]+0.7071067811865475*B0[1]*p2_over_gamma[3]-0.7071067811865475*B2[1]*p0_over_gamma[3]; 
  alpha[9] = 1.224744871391589*B0[2]*p2_over_gamma[6]-1.224744871391589*B2[2]*p0_over_gamma[6]+0.7071067811865475*B0[2]*p2_over_gamma[3]-0.7071067811865475*B2[2]*p0_over_gamma[3]; 
  alpha[10] = 1.224744871391589*B0[0]*p2_over_gamma[7]-1.224744871391589*B2[0]*p0_over_gamma[7]+0.7071067811865475*B0[0]*p2_over_gamma[5]-0.7071067811865475*B2[0]*p0_over_gamma[5]; 
  alpha[11] = 1.224744871391589*B0[3]*p2_over_gamma[4]-1.224744871391589*B2[3]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[3]+0.7071067811865475*p2_over_gamma[1]*B0[3]; 
  alpha[12] = 1.224744871391589*B0[3]*p2_over_gamma[6]-1.224744871391589*B2[3]*p0_over_gamma[6]+0.7071067811865475*B0[3]*p2_over_gamma[3]-0.7071067811865475*B2[3]*p0_over_gamma[3]; 
  alpha[13] = 1.224744871391589*B0[1]*p2_over_gamma[7]-1.224744871391589*B2[1]*p0_over_gamma[7]+0.7071067811865475*B0[1]*p2_over_gamma[5]-0.7071067811865475*B2[1]*p0_over_gamma[5]; 
  alpha[14] = 1.224744871391589*B0[2]*p2_over_gamma[7]-1.224744871391589*B2[2]*p0_over_gamma[7]+0.7071067811865475*B0[2]*p2_over_gamma[5]-0.7071067811865475*B2[2]*p0_over_gamma[5]; 
  alpha[15] = 1.224744871391589*B0[3]*p2_over_gamma[7]-1.224744871391589*B2[3]*p0_over_gamma[7]+0.7071067811865475*B0[3]*p2_over_gamma[5]-0.7071067811865475*B2[3]*p0_over_gamma[5]; 

  if (alpha[15]-alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_eval_quad_node_0_l(fEdge); 
  } 
  if ((-alpha[15])+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_eval_quad_node_1_l(fEdge); 
  } 
  if ((-alpha[15])+alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if (alpha[15]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if ((-alpha[15])+alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p1_surfx4_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_5x_p1_surfx4_eval_quad_node_4_l(fEdge); 
  } 
  if (alpha[15]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p1_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_5x_p1_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if (alpha[15]-alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p1_surfx4_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_5x_p1_surfx4_eval_quad_node_6_l(fEdge); 
  } 
  if ((-alpha[15])+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p1_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_5x_p1_surfx4_eval_quad_node_7_l(fEdge); 
  } 
  if ((-alpha[15])-alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p1_surfx4_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_5x_p1_surfx4_eval_quad_node_8_l(fEdge); 
  } 
  if (alpha[15]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p1_surfx4_eval_quad_node_9_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_5x_p1_surfx4_eval_quad_node_9_l(fEdge); 
  } 
  if (alpha[15]+alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p1_surfx4_eval_quad_node_10_r(fSkin); 
  } else { 
    fUpwindQuad[10] = ser_5x_p1_surfx4_eval_quad_node_10_l(fEdge); 
  } 
  if ((-alpha[15])-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p1_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[11] = ser_5x_p1_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if (alpha[15]+alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p1_surfx4_eval_quad_node_12_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_5x_p1_surfx4_eval_quad_node_12_l(fEdge); 
  } 
  if ((-alpha[15])-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p1_surfx4_eval_quad_node_13_r(fSkin); 
  } else { 
    fUpwindQuad[13] = ser_5x_p1_surfx4_eval_quad_node_13_l(fEdge); 
  } 
  if ((-alpha[15])-alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p1_surfx4_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[14] = ser_5x_p1_surfx4_eval_quad_node_14_l(fEdge); 
  } 
  if (alpha[15]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p1_surfx4_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_5x_p1_surfx4_eval_quad_node_15_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*(alpha[15]*fUpwind[15]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alpha[14]*fUpwind[15]+fUpwind[14]*alpha[15]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.25*(alpha[13]*fUpwind[15]+fUpwind[13]*alpha[15]+alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.25*(alpha[12]*fUpwind[15]+fUpwind[12]*alpha[15]+alpha[9]*fUpwind[14]+fUpwind[9]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.25*(alpha[11]*fUpwind[15]+fUpwind[11]*alpha[15]+alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] = 0.25*(alpha[10]*fUpwind[15]+fUpwind[10]*alpha[15]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] = 0.25*(alpha[9]*fUpwind[15]+fUpwind[9]*alpha[15]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[8]*fUpwind[10]+fUpwind[8]*alpha[10]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] = 0.25*(alpha[8]*fUpwind[15]+fUpwind[8]*alpha[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] = 0.25*(alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15]+alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14]+alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] = 0.25*(alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15]+alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[7]*fUpwind[10]+fUpwind[7]*alpha[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] = 0.25*(alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14]+alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]+alpha[0]*fUpwind[10]+fUpwind[0]*alpha[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] = 0.25*(alpha[4]*fUpwind[15]+fUpwind[4]*alpha[15]+alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[9]*fUpwind[13]+fUpwind[9]*alpha[13]+alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[12] = 0.25*(alpha[3]*fUpwind[15]+fUpwind[3]*alpha[15]+alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[13] = 0.25*(alpha[2]*fUpwind[15]+fUpwind[2]*alpha[15]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14]+alpha[0]*fUpwind[13]+fUpwind[0]*alpha[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[14] = 0.25*(alpha[1]*fUpwind[15]+fUpwind[1]*alpha[15]+alpha[0]*fUpwind[14]+fUpwind[0]*alpha[14]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[15] = 0.25*(alpha[0]*fUpwind[15]+fUpwind[0]*alpha[15]+alpha[1]*fUpwind[14]+fUpwind[1]*alpha[14]+alpha[2]*fUpwind[13]+fUpwind[2]*alpha[13]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]); 

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

  alpha[0] = (-1.224744871391589*B0[0]*p2_over_gamma[2])+1.224744871391589*B2[0]*p0_over_gamma[2]+0.7071067811865475*B0[0]*p2_over_gamma[0]-0.7071067811865475*B2[0]*p0_over_gamma[0]+2.0*E1[0]; 
  alpha[1] = (-1.224744871391589*B0[1]*p2_over_gamma[2])+1.224744871391589*B2[1]*p0_over_gamma[2]+2.0*E1[1]-0.7071067811865475*p0_over_gamma[0]*B2[1]+0.7071067811865475*p2_over_gamma[0]*B0[1]; 
  alpha[2] = (-1.224744871391589*B0[2]*p2_over_gamma[2])+1.224744871391589*B2[2]*p0_over_gamma[2]+2.0*E1[2]-0.7071067811865475*p0_over_gamma[0]*B2[2]+0.7071067811865475*p2_over_gamma[0]*B0[2]; 
  alpha[3] = (-1.224744871391589*B0[0]*p2_over_gamma[4])+1.224744871391589*B2[0]*p0_over_gamma[4]+0.7071067811865475*B0[0]*p2_over_gamma[1]-0.7071067811865475*B2[0]*p0_over_gamma[1]; 
  alpha[4] = (-1.224744871391589*B0[0]*p2_over_gamma[6])+1.224744871391589*B2[0]*p0_over_gamma[6]+0.7071067811865475*B0[0]*p2_over_gamma[3]-0.7071067811865475*B2[0]*p0_over_gamma[3]; 
  alpha[5] = 2.0*E1[3]+1.224744871391589*p0_over_gamma[2]*B2[3]-0.7071067811865475*p0_over_gamma[0]*B2[3]-1.224744871391589*p2_over_gamma[2]*B0[3]+0.7071067811865475*p2_over_gamma[0]*B0[3]; 
  alpha[6] = (-1.224744871391589*B0[1]*p2_over_gamma[4])+1.224744871391589*B2[1]*p0_over_gamma[4]+0.7071067811865475*B0[1]*p2_over_gamma[1]-0.7071067811865475*B2[1]*p0_over_gamma[1]; 
  alpha[7] = (-1.224744871391589*B0[2]*p2_over_gamma[4])+1.224744871391589*B2[2]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[2]+0.7071067811865475*p2_over_gamma[1]*B0[2]; 
  alpha[8] = (-1.224744871391589*B0[1]*p2_over_gamma[6])+1.224744871391589*B2[1]*p0_over_gamma[6]+0.7071067811865475*B0[1]*p2_over_gamma[3]-0.7071067811865475*B2[1]*p0_over_gamma[3]; 
  alpha[9] = (-1.224744871391589*B0[2]*p2_over_gamma[6])+1.224744871391589*B2[2]*p0_over_gamma[6]+0.7071067811865475*B0[2]*p2_over_gamma[3]-0.7071067811865475*B2[2]*p0_over_gamma[3]; 
  alpha[10] = (-1.224744871391589*B0[0]*p2_over_gamma[7])+1.224744871391589*B2[0]*p0_over_gamma[7]+0.7071067811865475*B0[0]*p2_over_gamma[5]-0.7071067811865475*B2[0]*p0_over_gamma[5]; 
  alpha[11] = (-1.224744871391589*B0[3]*p2_over_gamma[4])+1.224744871391589*B2[3]*p0_over_gamma[4]-0.7071067811865475*p0_over_gamma[1]*B2[3]+0.7071067811865475*p2_over_gamma[1]*B0[3]; 
  alpha[12] = (-1.224744871391589*B0[3]*p2_over_gamma[6])+1.224744871391589*B2[3]*p0_over_gamma[6]+0.7071067811865475*B0[3]*p2_over_gamma[3]-0.7071067811865475*B2[3]*p0_over_gamma[3]; 
  alpha[13] = (-1.224744871391589*B0[1]*p2_over_gamma[7])+1.224744871391589*B2[1]*p0_over_gamma[7]+0.7071067811865475*B0[1]*p2_over_gamma[5]-0.7071067811865475*B2[1]*p0_over_gamma[5]; 
  alpha[14] = (-1.224744871391589*B0[2]*p2_over_gamma[7])+1.224744871391589*B2[2]*p0_over_gamma[7]+0.7071067811865475*B0[2]*p2_over_gamma[5]-0.7071067811865475*B2[2]*p0_over_gamma[5]; 
  alpha[15] = (-1.224744871391589*B0[3]*p2_over_gamma[7])+1.224744871391589*B2[3]*p0_over_gamma[7]+0.7071067811865475*B0[3]*p2_over_gamma[5]-0.7071067811865475*B2[3]*p0_over_gamma[5]; 

  if (alpha[15]-alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p1_surfx4_eval_quad_node_0_l(fSkin); 
  } 
  if ((-alpha[15])+alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_5x_p1_surfx4_eval_quad_node_1_l(fSkin); 
  } 
  if ((-alpha[15])+alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_5x_p1_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if (alpha[15]-alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_5x_p1_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if ((-alpha[15])+alpha[14]-alpha[13]+alpha[12]+alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_5x_p1_surfx4_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_5x_p1_surfx4_eval_quad_node_4_l(fSkin); 
  } 
  if (alpha[15]-alpha[14]+alpha[13]-alpha[12]+alpha[11]-alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_5x_p1_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_5x_p1_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if (alpha[15]-alpha[14]+alpha[13]+alpha[12]-alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_5x_p1_surfx4_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_5x_p1_surfx4_eval_quad_node_6_l(fSkin); 
  } 
  if ((-alpha[15])+alpha[14]-alpha[13]-alpha[12]-alpha[11]+alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_5x_p1_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_5x_p1_surfx4_eval_quad_node_7_l(fSkin); 
  } 
  if ((-alpha[15])-alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]-alpha[8]+alpha[7]-alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[8] = ser_5x_p1_surfx4_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_5x_p1_surfx4_eval_quad_node_8_l(fSkin); 
  } 
  if (alpha[15]+alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]+alpha[8]+alpha[7]-alpha[6]-alpha[5]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[9] = ser_5x_p1_surfx4_eval_quad_node_9_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_5x_p1_surfx4_eval_quad_node_9_l(fSkin); 
  } 
  if (alpha[15]+alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]-alpha[8]-alpha[7]+alpha[6]-alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[10] = ser_5x_p1_surfx4_eval_quad_node_10_r(fEdge); 
  } else { 
    fUpwindQuad[10] = ser_5x_p1_surfx4_eval_quad_node_10_l(fSkin); 
  } 
  if ((-alpha[15])-alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]+alpha[8]-alpha[7]+alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[11] = ser_5x_p1_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[11] = ser_5x_p1_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if (alpha[15]+alpha[14]+alpha[13]-alpha[12]-alpha[11]+alpha[10]-alpha[9]-alpha[8]-alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[12] = ser_5x_p1_surfx4_eval_quad_node_12_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_5x_p1_surfx4_eval_quad_node_12_l(fSkin); 
  } 
  if ((-alpha[15])-alpha[14]-alpha[13]+alpha[12]-alpha[11]-alpha[10]+alpha[9]+alpha[8]-alpha[7]-alpha[6]+alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[13] = ser_5x_p1_surfx4_eval_quad_node_13_r(fEdge); 
  } else { 
    fUpwindQuad[13] = ser_5x_p1_surfx4_eval_quad_node_13_l(fSkin); 
  } 
  if ((-alpha[15])-alpha[14]-alpha[13]-alpha[12]+alpha[11]-alpha[10]-alpha[9]-alpha[8]+alpha[7]+alpha[6]+alpha[5]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[14] = ser_5x_p1_surfx4_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[14] = ser_5x_p1_surfx4_eval_quad_node_14_l(fSkin); 
  } 
  if (alpha[15]+alpha[14]+alpha[13]+alpha[12]+alpha[11]+alpha[10]+alpha[9]+alpha[8]+alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[15] = ser_5x_p1_surfx4_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_5x_p1_surfx4_eval_quad_node_15_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*(alpha[15]*fUpwind[15]+alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[12]*fUpwind[12]+alpha[11]*fUpwind[11]+alpha[10]*fUpwind[10]+alpha[9]*fUpwind[9]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.25*(alpha[14]*fUpwind[15]+fUpwind[14]*alpha[15]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[9]*fUpwind[12]+fUpwind[9]*alpha[12]+alpha[7]*fUpwind[11]+fUpwind[7]*alpha[11]+alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.25*(alpha[13]*fUpwind[15]+fUpwind[13]*alpha[15]+alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[8]*fUpwind[12]+fUpwind[8]*alpha[12]+alpha[6]*fUpwind[11]+fUpwind[6]*alpha[11]+alpha[4]*fUpwind[9]+fUpwind[4]*alpha[9]+alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.25*(alpha[12]*fUpwind[15]+fUpwind[12]*alpha[15]+alpha[9]*fUpwind[14]+fUpwind[9]*alpha[14]+alpha[8]*fUpwind[13]+fUpwind[8]*alpha[13]+alpha[5]*fUpwind[11]+fUpwind[5]*alpha[11]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10]+alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.25*(alpha[11]*fUpwind[15]+fUpwind[11]*alpha[15]+alpha[7]*fUpwind[14]+fUpwind[7]*alpha[14]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13]+alpha[5]*fUpwind[12]+fUpwind[5]*alpha[12]+alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10]+alpha[2]*fUpwind[9]+fUpwind[2]*alpha[9]+alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]); 
  Ghat[5] = 0.25*(alpha[10]*fUpwind[15]+fUpwind[10]*alpha[15]+alpha[13]*fUpwind[14]+fUpwind[13]*alpha[14]+alpha[4]*fUpwind[12]+fUpwind[4]*alpha[12]+alpha[3]*fUpwind[11]+fUpwind[3]*alpha[11]+alpha[8]*fUpwind[9]+fUpwind[8]*alpha[9]+alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[6] = 0.25*(alpha[9]*fUpwind[15]+fUpwind[9]*alpha[15]+alpha[12]*fUpwind[14]+fUpwind[12]*alpha[14]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[2]*fUpwind[11]+fUpwind[2]*alpha[11]+alpha[8]*fUpwind[10]+fUpwind[8]*alpha[10]+alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[7] = 0.25*(alpha[8]*fUpwind[15]+fUpwind[8]*alpha[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[12]*fUpwind[13]+fUpwind[12]*alpha[13]+alpha[1]*fUpwind[11]+fUpwind[1]*alpha[11]+alpha[9]*fUpwind[10]+fUpwind[9]*alpha[10]+alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[8] = 0.25*(alpha[7]*fUpwind[15]+fUpwind[7]*alpha[15]+alpha[11]*fUpwind[14]+fUpwind[11]*alpha[14]+alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[2]*fUpwind[12]+fUpwind[2]*alpha[12]+alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10]+alpha[5]*fUpwind[9]+fUpwind[5]*alpha[9]+alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[9] = 0.25*(alpha[6]*fUpwind[15]+fUpwind[6]*alpha[15]+alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[11]*fUpwind[13]+fUpwind[11]*alpha[13]+alpha[1]*fUpwind[12]+fUpwind[1]*alpha[12]+alpha[7]*fUpwind[10]+fUpwind[7]*alpha[10]+alpha[0]*fUpwind[9]+fUpwind[0]*alpha[9]+alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[10] = 0.25*(alpha[5]*fUpwind[15]+fUpwind[5]*alpha[15]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14]+alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[11]*fUpwind[12]+fUpwind[11]*alpha[12]+alpha[0]*fUpwind[10]+fUpwind[0]*alpha[10]+alpha[7]*fUpwind[9]+fUpwind[7]*alpha[9]+alpha[6]*fUpwind[8]+fUpwind[6]*alpha[8]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] = 0.25*(alpha[4]*fUpwind[15]+fUpwind[4]*alpha[15]+alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[9]*fUpwind[13]+fUpwind[9]*alpha[13]+alpha[10]*fUpwind[12]+fUpwind[10]*alpha[12]+alpha[0]*fUpwind[11]+fUpwind[0]*alpha[11]+alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[12] = 0.25*(alpha[3]*fUpwind[15]+fUpwind[3]*alpha[15]+alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[0]*fUpwind[12]+fUpwind[0]*alpha[12]+alpha[10]*fUpwind[11]+fUpwind[10]*alpha[11]+alpha[1]*fUpwind[9]+fUpwind[1]*alpha[9]+alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[13] = 0.25*(alpha[2]*fUpwind[15]+fUpwind[2]*alpha[15]+alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14]+alpha[0]*fUpwind[13]+fUpwind[0]*alpha[13]+alpha[7]*fUpwind[12]+fUpwind[7]*alpha[12]+alpha[9]*fUpwind[11]+fUpwind[9]*alpha[11]+alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10]+alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[14] = 0.25*(alpha[1]*fUpwind[15]+fUpwind[1]*alpha[15]+alpha[0]*fUpwind[14]+fUpwind[0]*alpha[14]+alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13]+alpha[6]*fUpwind[12]+fUpwind[6]*alpha[12]+alpha[8]*fUpwind[11]+fUpwind[8]*alpha[11]+alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10]+alpha[3]*fUpwind[9]+fUpwind[3]*alpha[9]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]); 
  Ghat[15] = 0.25*(alpha[0]*fUpwind[15]+fUpwind[0]*alpha[15]+alpha[1]*fUpwind[14]+fUpwind[1]*alpha[14]+alpha[2]*fUpwind[13]+fUpwind[2]*alpha[13]+alpha[3]*fUpwind[12]+fUpwind[3]*alpha[12]+alpha[4]*fUpwind[11]+fUpwind[4]*alpha[11]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10]+alpha[6]*fUpwind[9]+fUpwind[6]*alpha[9]+alpha[7]*fUpwind[8]+fUpwind[7]*alpha[8]); 

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
  return 0.;

} 
