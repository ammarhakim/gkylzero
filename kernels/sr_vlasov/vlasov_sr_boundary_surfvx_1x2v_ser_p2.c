#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_boundary_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // gamma:       Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double *E0 = &qmem[0]; 
  double p1_over_gamma_l[3] = {0.0}; 
  double p1_over_gamma_r[3] = {0.0}; 
  p1_over_gamma_l[0] = 2.738612787525831*gamma[6]*dv11-2.121320343559642*gamma[3]*dv11+1.224744871391589*gamma[2]*dv11; 
  p1_over_gamma_l[1] = 2.738612787525831*gamma[5]*dv11-4.743416490252569*gamma[7]*dv11; 
  p1_over_gamma_r[0] = 2.738612787525831*gamma[6]*dv11+2.121320343559642*gamma[3]*dv11+1.224744871391589*gamma[2]*dv11; 
  p1_over_gamma_r[1] = 4.743416490252569*gamma[7]*dv11+2.738612787525831*gamma[5]*dv11; 

  const double *B2 = &qmem[15]; 

  double alpha[8] = {0.0}; 

  double fUpwindQuad[9] = {0.0};
  double fUpwind[8] = {0.0};
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  alpha[0] = B2[0]*p1_over_gamma_r[0]+1.414213562373095*E0[0]; 
  alpha[1] = 1.414213562373095*E0[1]+p1_over_gamma_r[0]*B2[1]; 
  alpha[2] = B2[0]*p1_over_gamma_r[1]; 
  alpha[3] = B2[1]*p1_over_gamma_r[1]; 
  alpha[4] = 1.414213562373095*E0[2]+p1_over_gamma_r[0]*B2[2]; 
  alpha[6] = 1.0*p1_over_gamma_r[1]*B2[2]; 

  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]+0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]-0.9*alpha[3]+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(fEdge); 
  } 
  if (0.75*alpha[6]-0.5590169943749475*alpha[4]-0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(fEdge); 
  } 
  if (0.5*alpha[0]-0.5590169943749475*alpha[4] > 0) { 
    fUpwindQuad[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(fEdge); 
  } 
  if ((-0.75*alpha[6])-0.5590169943749475*alpha[4]+0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(fEdge); 
  } 
  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]-0.9*alpha[3]-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(fEdge); 
  } 
  if (0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(fEdge); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]+0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[6]*fUpwind[6]+0.5*alpha[4]*fUpwind[4]+0.5*alpha[3]*fUpwind[3]+0.5*alpha[2]*fUpwind[2]+0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.447213595499958*alpha[3]*fUpwind[6]+0.447213595499958*fUpwind[3]*alpha[6]+0.4472135954999579*alpha[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alpha[4]+0.5*alpha[2]*fUpwind[3]+0.5*fUpwind[2]*alpha[3]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.447213595499958*alpha[3]*fUpwind[7]+0.5000000000000001*alpha[4]*fUpwind[6]+0.5000000000000001*fUpwind[4]*alpha[6]+0.4472135954999579*alpha[2]*fUpwind[5]+0.5*alpha[1]*fUpwind[3]+0.5*fUpwind[1]*alpha[3]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.4*alpha[6]*fUpwind[7]+0.447213595499958*alpha[2]*fUpwind[7]+0.447213595499958*alpha[1]*fUpwind[6]+0.447213595499958*fUpwind[1]*alpha[6]+0.4472135954999579*alpha[3]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[4]+0.4472135954999579*fUpwind[3]*alpha[4]+0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[0]*alpha[3]+0.5*alpha[1]*fUpwind[2]+0.5*fUpwind[1]*alpha[2]; 
  Ghat[4] = 0.31943828249997*alpha[6]*fUpwind[6]+0.5000000000000001*alpha[2]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alpha[6]+0.31943828249997*alpha[4]*fUpwind[4]+0.5*alpha[0]*fUpwind[4]+0.5*fUpwind[0]*alpha[4]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[1]*fUpwind[1]; 
  Ghat[5] = 0.5000000000000001*alpha[1]*fUpwind[7]+0.4472135954999579*alpha[6]*fUpwind[6]+0.5*alpha[0]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[2]*fUpwind[2]; 
  Ghat[6] = 0.4*alpha[3]*fUpwind[7]+0.31943828249997*alpha[4]*fUpwind[6]+0.5*alpha[0]*fUpwind[6]+0.4472135954999579*fUpwind[5]*alpha[6]+0.31943828249997*fUpwind[4]*alpha[6]+0.5*fUpwind[0]*alpha[6]+0.5000000000000001*alpha[2]*fUpwind[4]+0.5000000000000001*fUpwind[2]*alpha[4]+0.447213595499958*alpha[1]*fUpwind[3]+0.447213595499958*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.4472135954999579*alpha[4]*fUpwind[7]+0.5*alpha[0]*fUpwind[7]+0.4*alpha[3]*fUpwind[6]+0.4*fUpwind[3]*alpha[6]+0.5000000000000001*alpha[1]*fUpwind[5]+0.447213595499958*alpha[2]*fUpwind[3]+0.447213595499958*fUpwind[2]*alpha[3]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += -0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -0.7071067811865475*Ghat[4]*dv10; 
  out[8] += -1.58113883008419*Ghat[0]*dv10; 
  out[9] += -0.7071067811865475*Ghat[5]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += -1.58113883008419*Ghat[1]*dv10; 
  out[13] += -0.7071067811865475*Ghat[6]*dv10; 
  out[14] += -1.58113883008419*Ghat[2]*dv10; 
  out[15] += -0.7071067811865475*Ghat[7]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += -1.224744871391589*Ghat[6]*dv10; 
  out[18] += -1.58113883008419*Ghat[3]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 

  } else { 

  alpha[0] = B2[0]*p1_over_gamma_l[0]+1.414213562373095*E0[0]; 
  alpha[1] = 1.414213562373095*E0[1]+p1_over_gamma_l[0]*B2[1]; 
  alpha[2] = B2[0]*p1_over_gamma_l[1]; 
  alpha[3] = B2[1]*p1_over_gamma_l[1]; 
  alpha[4] = 1.414213562373095*E0[2]+p1_over_gamma_l[0]*B2[2]; 
  alpha[6] = 1.0*p1_over_gamma_l[1]*B2[2]; 

  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]+0.9*alpha[3]-0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.4472135954999579*alpha[4]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]-0.9*alpha[3]+0.6708203932499369*alpha[2]-0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(fSkin); 
  } 
  if (0.75*alpha[6]-0.5590169943749475*alpha[4]-0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(fSkin); 
  } 
  if (0.5*alpha[0]-0.5590169943749475*alpha[4] > 0) { 
    fUpwindQuad[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(fSkin); 
  } 
  if ((-0.75*alpha[6])-0.5590169943749475*alpha[4]+0.6708203932499369*alpha[2]+0.5*alpha[0] > 0) { 
    fUpwindQuad[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(fSkin); 
  } 
  if ((-0.5999999999999999*alpha[6])+0.4472135954999579*alpha[4]-0.9*alpha[3]-0.6708203932499369*alpha[2]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(fSkin); 
  } 
  if (0.4472135954999579*alpha[4]+0.6708203932499369*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(fSkin); 
  } 
  if (0.5999999999999999*alpha[6]+0.4472135954999579*alpha[4]+0.9*alpha[3]+0.6708203932499369*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[6]*fUpwind[6]+0.5*alpha[4]*fUpwind[4]+0.5*alpha[3]*fUpwind[3]+0.5*alpha[2]*fUpwind[2]+0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.447213595499958*alpha[3]*fUpwind[6]+0.447213595499958*fUpwind[3]*alpha[6]+0.4472135954999579*alpha[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alpha[4]+0.5*alpha[2]*fUpwind[3]+0.5*fUpwind[2]*alpha[3]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.447213595499958*alpha[3]*fUpwind[7]+0.5000000000000001*alpha[4]*fUpwind[6]+0.5000000000000001*fUpwind[4]*alpha[6]+0.4472135954999579*alpha[2]*fUpwind[5]+0.5*alpha[1]*fUpwind[3]+0.5*fUpwind[1]*alpha[3]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.4*alpha[6]*fUpwind[7]+0.447213595499958*alpha[2]*fUpwind[7]+0.447213595499958*alpha[1]*fUpwind[6]+0.447213595499958*fUpwind[1]*alpha[6]+0.4472135954999579*alpha[3]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[4]+0.4472135954999579*fUpwind[3]*alpha[4]+0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[0]*alpha[3]+0.5*alpha[1]*fUpwind[2]+0.5*fUpwind[1]*alpha[2]; 
  Ghat[4] = 0.31943828249997*alpha[6]*fUpwind[6]+0.5000000000000001*alpha[2]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alpha[6]+0.31943828249997*alpha[4]*fUpwind[4]+0.5*alpha[0]*fUpwind[4]+0.5*fUpwind[0]*alpha[4]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[1]*fUpwind[1]; 
  Ghat[5] = 0.5000000000000001*alpha[1]*fUpwind[7]+0.4472135954999579*alpha[6]*fUpwind[6]+0.5*alpha[0]*fUpwind[5]+0.4472135954999579*alpha[3]*fUpwind[3]+0.4472135954999579*alpha[2]*fUpwind[2]; 
  Ghat[6] = 0.4*alpha[3]*fUpwind[7]+0.31943828249997*alpha[4]*fUpwind[6]+0.5*alpha[0]*fUpwind[6]+0.4472135954999579*fUpwind[5]*alpha[6]+0.31943828249997*fUpwind[4]*alpha[6]+0.5*fUpwind[0]*alpha[6]+0.5000000000000001*alpha[2]*fUpwind[4]+0.5000000000000001*fUpwind[2]*alpha[4]+0.447213595499958*alpha[1]*fUpwind[3]+0.447213595499958*fUpwind[1]*alpha[3]; 
  Ghat[7] = 0.4472135954999579*alpha[4]*fUpwind[7]+0.5*alpha[0]*fUpwind[7]+0.4*alpha[3]*fUpwind[6]+0.4*fUpwind[3]*alpha[6]+0.5000000000000001*alpha[1]*fUpwind[5]+0.447213595499958*alpha[2]*fUpwind[3]+0.447213595499958*fUpwind[2]*alpha[3]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += 0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += 0.7071067811865475*Ghat[4]*dv10; 
  out[8] += 1.58113883008419*Ghat[0]*dv10; 
  out[9] += 0.7071067811865475*Ghat[5]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += -1.224744871391589*Ghat[4]*dv10; 
  out[12] += 1.58113883008419*Ghat[1]*dv10; 
  out[13] += 0.7071067811865475*Ghat[6]*dv10; 
  out[14] += 1.58113883008419*Ghat[2]*dv10; 
  out[15] += 0.7071067811865475*Ghat[7]*dv10; 
  out[16] += -1.224744871391589*Ghat[5]*dv10; 
  out[17] += -1.224744871391589*Ghat[6]*dv10; 
  out[18] += 1.58113883008419*Ghat[3]*dv10; 
  out[19] += -1.224744871391589*Ghat[7]*dv10; 

  } 
  return 0.;

} 
