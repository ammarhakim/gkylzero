#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void lbo_vlasov_drag_boundary_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  const double *sumNuUy = &nuUSum[2]; 

  double alphaDrSurf[4] = {0.0}; 
  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};;
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[2]+1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[2]+1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = ser_3x_p1_surfx3_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_3x_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_3x_p1_surfx3_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_3x_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[2] = ser_3x_p1_surfx3_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = ser_3x_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_3x_p1_surfx3_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = ser_3x_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.5*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.5*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[2]-1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[2]-1.414213562373095*dxv[2])-2.828427124746191*sumNuUy[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = ser_3x_p1_surfx3_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_3x_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_3x_p1_surfx3_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_3x_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[2] = ser_3x_p1_surfx3_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = ser_3x_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_3x_p1_surfx3_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = ser_3x_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.5*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.5*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 1.224744871391589*Ghat[1]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 

  } 
} 
