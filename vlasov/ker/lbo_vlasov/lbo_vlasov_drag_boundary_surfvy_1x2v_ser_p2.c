#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[9]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUy = &nuPrimMomsSum[3]; 

  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[9] = {0.0};
  double fUpwind[8] = {0.0};;
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[2]+1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[2]+1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[1]); 
  alphaDrSurf[4] = 0.5*(2.828427124746191*nuSum[2]*w[2]-2.828427124746191*sumNuUy[2]+1.414213562373095*dxv[2]*nuSum[2]); 

  if (0.4472135954999572*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_3x_p2_surfx3_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_3x_p2_surfx3_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = ser_3x_p2_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_3x_p2_surfx3_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_3x_p2_surfx3_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = ser_3x_p2_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (0.4472135954999572*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_3x_p2_surfx3_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = ser_3x_p2_surfx3_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = ser_3x_p2_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_3x_p2_surfx3_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = ser_3x_p2_surfx3_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = ser_3x_p2_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if (0.4472135954999572*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_3x_p2_surfx3_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = ser_3x_p2_surfx3_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = ser_3x_p2_surfx3_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_3x_p2_surfx3_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = ser_3x_p2_surfx3_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = ser_3x_p2_surfx3_eval_quad_node_8_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[3] = 0.447213595499958*alphaDrSurf[1]*fUpwind[6]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*alphaDrSurf[1]*fUpwind[2]; 
  Ghat[4] = 0.31943828249997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[5] = 0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[5]; 
  Ghat[6] = 0.31943828249997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.447213595499958*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[7] = 0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += 1.58113883008419*Ghat[0]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[4]*rdv2; 
  out[14] += 1.224744871391589*Ghat[5]*rdv2; 
  out[15] += 1.58113883008419*Ghat[1]*rdv2; 
  out[16] += 1.58113883008419*Ghat[2]*rdv2; 
  out[17] += 1.224744871391589*Ghat[6]*rdv2; 
  out[18] += 1.224744871391589*Ghat[7]*rdv2; 
  out[19] += 1.58113883008419*Ghat[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[2]-1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[2]-1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[1]); 
  alphaDrSurf[4] = 0.5*(2.828427124746191*nuSum[2]*w[2]-2.828427124746191*sumNuUy[2]-1.414213562373095*dxv[2]*nuSum[2]); 

  if (0.4472135954999572*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_3x_p2_surfx3_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_3x_p2_surfx3_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = ser_3x_p2_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_3x_p2_surfx3_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_3x_p2_surfx3_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = ser_3x_p2_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (0.4472135954999572*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_3x_p2_surfx3_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = ser_3x_p2_surfx3_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = ser_3x_p2_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_3x_p2_surfx3_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = ser_3x_p2_surfx3_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = ser_3x_p2_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if (0.4472135954999572*alphaDrSurf[4]-0.6708203932499357*alphaDrSurf[1]+0.5*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_3x_p2_surfx3_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = ser_3x_p2_surfx3_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = ser_3x_p2_surfx3_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_3x_p2_surfx3_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = ser_3x_p2_surfx3_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = ser_3x_p2_surfx3_eval_quad_node_8_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.4472135954999579*alphaDrSurf[1]*fUpwind[4]+0.4472135954999579*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.5000000000000001*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[1]*fUpwind[3]+0.5*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[3] = 0.447213595499958*alphaDrSurf[1]*fUpwind[6]+0.4472135954999579*fUpwind[3]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[3]+0.5*alphaDrSurf[1]*fUpwind[2]; 
  Ghat[4] = 0.31943828249997*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.4472135954999579*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[5] = 0.5000000000000001*alphaDrSurf[1]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[5]; 
  Ghat[6] = 0.31943828249997*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[6]+0.5000000000000001*fUpwind[2]*alphaDrSurf[4]+0.447213595499958*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[7] = 0.4472135954999579*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[7]+0.5000000000000001*alphaDrSurf[1]*fUpwind[5]; 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += -1.58113883008419*Ghat[0]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[4]*rdv2; 
  out[14] += 1.224744871391589*Ghat[5]*rdv2; 
  out[15] += -1.58113883008419*Ghat[1]*rdv2; 
  out[16] += -1.58113883008419*Ghat[2]*rdv2; 
  out[17] += 1.224744871391589*Ghat[6]*rdv2; 
  out[18] += 1.224744871391589*Ghat[7]*rdv2; 
  out[19] += -1.58113883008419*Ghat[3]*rdv2; 

  } 

  return 0.;

} 
