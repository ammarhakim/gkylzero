#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_boundary_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *E1 = &qmem[2]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *B0 = &qmem[6]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double *B1 = &qmem[8]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
  const double *B2 = &qmem[10]; 

  double alpha[8] = {0.0}; 

  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  alpha[0] = 1.224744871391589*B0[0]*p2_over_gamma[2]-1.224744871391589*B2[0]*p0_over_gamma[2]+0.7071067811865475*B0[0]*p2_over_gamma[0]-0.7071067811865475*B2[0]*p0_over_gamma[0]+2.0*E1[0]; 
  alpha[1] = 1.224744871391589*B0[1]*p2_over_gamma[2]-1.224744871391589*B2[1]*p0_over_gamma[2]+2.0*E1[1]-0.7071067811865475*p0_over_gamma[0]*B2[1]+0.7071067811865475*p2_over_gamma[0]*B0[1]; 
  alpha[2] = 1.224744871391589*B0[0]*p2_over_gamma[4]-1.224744871391589*B2[0]*p0_over_gamma[4]+0.7071067811865475*B0[0]*p2_over_gamma[1]-0.7071067811865475*B2[0]*p0_over_gamma[1]; 
  alpha[3] = 1.224744871391589*B0[0]*p2_over_gamma[6]-1.224744871391589*B2[0]*p0_over_gamma[6]+0.7071067811865475*B0[0]*p2_over_gamma[3]-0.7071067811865475*B2[0]*p0_over_gamma[3]; 
  alpha[4] = 1.224744871391589*B0[1]*p2_over_gamma[4]-1.224744871391589*B2[1]*p0_over_gamma[4]+0.7071067811865475*B0[1]*p2_over_gamma[1]-0.7071067811865475*B2[1]*p0_over_gamma[1]; 
  alpha[5] = 1.224744871391589*B0[1]*p2_over_gamma[6]-1.224744871391589*B2[1]*p0_over_gamma[6]+0.7071067811865475*B0[1]*p2_over_gamma[3]-0.7071067811865475*B2[1]*p0_over_gamma[3]; 
  alpha[6] = 1.224744871391589*B0[0]*p2_over_gamma[7]-1.224744871391589*B2[0]*p0_over_gamma[7]+0.7071067811865475*B0[0]*p2_over_gamma[5]-0.7071067811865475*B2[0]*p0_over_gamma[5]; 
  alpha[7] = 1.224744871391589*B0[1]*p2_over_gamma[7]-1.224744871391589*B2[1]*p0_over_gamma[7]+0.7071067811865475*B0[1]*p2_over_gamma[5]-0.7071067811865475*B2[1]*p0_over_gamma[5]; 

  if ((-alpha[7])+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if (alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if ((-alpha[7])+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if (alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p1_surfx3_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_4x_p1_surfx3_eval_quad_node_4_l(fEdge); 
  } 
  if ((-alpha[7])-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p1_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_4x_p1_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if ((-alpha[7])-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p1_surfx3_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_4x_p1_surfx3_eval_quad_node_6_l(fEdge); 
  } 
  if (alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p1_surfx3_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_4x_p1_surfx3_eval_quad_node_7_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] = 0.3535533905932737*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] = 0.3535533905932737*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += -0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -0.7071067811865475*Ghat[4]*dv11; 
  out[6] += -1.224744871391589*Ghat[1]*dv11; 
  out[7] += -1.224744871391589*Ghat[2]*dv11; 
  out[8] += -0.7071067811865475*Ghat[5]*dv11; 
  out[9] += -0.7071067811865475*Ghat[6]*dv11; 
  out[10] += -1.224744871391589*Ghat[3]*dv11; 
  out[11] += -1.224744871391589*Ghat[4]*dv11; 
  out[12] += -0.7071067811865475*Ghat[7]*dv11; 
  out[13] += -1.224744871391589*Ghat[5]*dv11; 
  out[14] += -1.224744871391589*Ghat[6]*dv11; 
  out[15] += -1.224744871391589*Ghat[7]*dv11; 

  } else { 

  alpha[0] = (-1.224744871391589*B0[0]*p2_over_gamma[2])+1.224744871391589*B2[0]*p0_over_gamma[2]+0.7071067811865475*B0[0]*p2_over_gamma[0]-0.7071067811865475*B2[0]*p0_over_gamma[0]+2.0*E1[0]; 
  alpha[1] = (-1.224744871391589*B0[1]*p2_over_gamma[2])+1.224744871391589*B2[1]*p0_over_gamma[2]+2.0*E1[1]-0.7071067811865475*p0_over_gamma[0]*B2[1]+0.7071067811865475*p2_over_gamma[0]*B0[1]; 
  alpha[2] = (-1.224744871391589*B0[0]*p2_over_gamma[4])+1.224744871391589*B2[0]*p0_over_gamma[4]+0.7071067811865475*B0[0]*p2_over_gamma[1]-0.7071067811865475*B2[0]*p0_over_gamma[1]; 
  alpha[3] = (-1.224744871391589*B0[0]*p2_over_gamma[6])+1.224744871391589*B2[0]*p0_over_gamma[6]+0.7071067811865475*B0[0]*p2_over_gamma[3]-0.7071067811865475*B2[0]*p0_over_gamma[3]; 
  alpha[4] = (-1.224744871391589*B0[1]*p2_over_gamma[4])+1.224744871391589*B2[1]*p0_over_gamma[4]+0.7071067811865475*B0[1]*p2_over_gamma[1]-0.7071067811865475*B2[1]*p0_over_gamma[1]; 
  alpha[5] = (-1.224744871391589*B0[1]*p2_over_gamma[6])+1.224744871391589*B2[1]*p0_over_gamma[6]+0.7071067811865475*B0[1]*p2_over_gamma[3]-0.7071067811865475*B2[1]*p0_over_gamma[3]; 
  alpha[6] = (-1.224744871391589*B0[0]*p2_over_gamma[7])+1.224744871391589*B2[0]*p0_over_gamma[7]+0.7071067811865475*B0[0]*p2_over_gamma[5]-0.7071067811865475*B2[0]*p0_over_gamma[5]; 
  alpha[7] = (-1.224744871391589*B0[1]*p2_over_gamma[7])+1.224744871391589*B2[1]*p0_over_gamma[7]+0.7071067811865475*B0[1]*p2_over_gamma[5]-0.7071067811865475*B2[1]*p0_over_gamma[5]; 

  if ((-alpha[7])+alpha[6]+alpha[5]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p1_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if (alpha[7]-alpha[6]-alpha[5]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_4x_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (alpha[7]-alpha[6]+alpha[5]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_4x_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if ((-alpha[7])+alpha[6]-alpha[5]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_4x_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if (alpha[7]+alpha[6]-alpha[5]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[4] = ser_4x_p1_surfx3_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_4x_p1_surfx3_eval_quad_node_4_l(fSkin); 
  } 
  if ((-alpha[7])-alpha[6]+alpha[5]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[5] = ser_4x_p1_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_4x_p1_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if ((-alpha[7])-alpha[6]-alpha[5]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[6] = ser_4x_p1_surfx3_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_4x_p1_surfx3_eval_quad_node_6_l(fSkin); 
  } 
  if (alpha[7]+alpha[6]+alpha[5]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[7] = ser_4x_p1_surfx3_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_4x_p1_surfx3_eval_quad_node_7_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7]+alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.3535533905932737*(alpha[5]*fUpwind[7]+fUpwind[5]*alpha[7]+alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.3535533905932737*(alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7]+alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] = 0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7]+alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] = 0.3535533905932737*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7]+alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] = 0.3535533905932737*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7]+alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] = 0.3535533905932737*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += 0.7071067811865475*Ghat[3]*dv11; 
  out[5] += 0.7071067811865475*Ghat[4]*dv11; 
  out[6] += -1.224744871391589*Ghat[1]*dv11; 
  out[7] += -1.224744871391589*Ghat[2]*dv11; 
  out[8] += 0.7071067811865475*Ghat[5]*dv11; 
  out[9] += 0.7071067811865475*Ghat[6]*dv11; 
  out[10] += -1.224744871391589*Ghat[3]*dv11; 
  out[11] += -1.224744871391589*Ghat[4]*dv11; 
  out[12] += 0.7071067811865475*Ghat[7]*dv11; 
  out[13] += -1.224744871391589*Ghat[5]*dv11; 
  out[14] += -1.224744871391589*Ghat[6]*dv11; 
  out[15] += -1.224744871391589*Ghat[7]*dv11; 

  } 
  return 0.;

} 
