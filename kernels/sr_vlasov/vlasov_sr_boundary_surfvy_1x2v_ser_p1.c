#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_sr_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
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
  const double *E1 = &qmem[2]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double *p1_over_gamma = &p_over_gamma[4]; 
  const double *B2 = &qmem[10]; 

  double alpha[4] = {0.0}; 

  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  alpha[0] = (-1.224744871391589*B2[0]*p0_over_gamma[2])-0.7071067811865475*B2[0]*p0_over_gamma[0]+1.414213562373095*E1[0]; 
  alpha[1] = (-1.224744871391589*B2[1]*p0_over_gamma[2])+1.414213562373095*E1[1]-0.7071067811865475*p0_over_gamma[0]*B2[1]; 
  alpha[2] = (-1.224744871391589*B2[0]*p0_over_gamma[3])-0.7071067811865475*B2[0]*p0_over_gamma[1]; 
  alpha[3] = (-1.224744871391589*B2[1]*p0_over_gamma[3])-0.7071067811865475*B2[1]*p0_over_gamma[1]; 

  if (alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_3x_p1_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_3x_p1_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if ((-alpha[3])+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_3x_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_3x_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if ((-alpha[3])-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_3x_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_3x_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_3x_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_3x_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv11; 
  out[1] += -0.7071067811865475*Ghat[1]*dv11; 
  out[2] += -0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += -0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -1.224744871391589*Ghat[1]*dv11; 
  out[6] += -1.224744871391589*Ghat[2]*dv11; 
  out[7] += -1.224744871391589*Ghat[3]*dv11; 

  } else { 

  alpha[0] = 1.224744871391589*B2[0]*p0_over_gamma[2]-0.7071067811865475*B2[0]*p0_over_gamma[0]+1.414213562373095*E1[0]; 
  alpha[1] = 1.224744871391589*B2[1]*p0_over_gamma[2]+1.414213562373095*E1[1]-0.7071067811865475*p0_over_gamma[0]*B2[1]; 
  alpha[2] = 1.224744871391589*B2[0]*p0_over_gamma[3]-0.7071067811865475*B2[0]*p0_over_gamma[1]; 
  alpha[3] = 1.224744871391589*B2[1]*p0_over_gamma[3]-0.7071067811865475*B2[1]*p0_over_gamma[1]; 

  if (alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[0] = ser_3x_p1_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_3x_p1_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if ((-alpha[3])+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[1] = ser_3x_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_3x_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if ((-alpha[3])-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[2] = ser_3x_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_3x_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad[3] = ser_3x_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_3x_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] = 0.5*(alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] = 0.5*(alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv11; 
  out[1] += 0.7071067811865475*Ghat[1]*dv11; 
  out[2] += 0.7071067811865475*Ghat[2]*dv11; 
  out[3] += -1.224744871391589*Ghat[0]*dv11; 
  out[4] += 0.7071067811865475*Ghat[3]*dv11; 
  out[5] += -1.224744871391589*Ghat[1]*dv11; 
  out[6] += -1.224744871391589*Ghat[2]*dv11; 
  out[7] += -1.224744871391589*Ghat[3]*dv11; 

  } 
} 
