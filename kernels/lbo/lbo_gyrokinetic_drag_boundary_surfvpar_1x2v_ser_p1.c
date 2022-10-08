#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_1x2v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[4]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 

  const double *nuUSum = nuPrimMomsSum;

  double rdv2 = 2.0/dxv[1]; 


  double alphaDrSurf[4] = {0.0}; 
  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};;
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[1]+1.414213562373095*dxv[1])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.828427124746191*nuSum[1]*w[1]-2.828427124746191*nuUSum[1]+1.414213562373095*dxv[1]*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_1x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.5*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.5*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[4] += 1.224744871391589*Ghat[1]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += 1.58113883008419*Ghat[0]*rdv2; 
  out[9] += 1.58113883008419*Ghat[1]*rdv2; 
  out[10] += 1.58113883008419*Ghat[2]*rdv2; 
  out[11] += 1.58113883008419*Ghat[3]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[1]-1.414213562373095*dxv[1])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.828427124746191*nuSum[1]*w[1]-2.828427124746191*nuUSum[1]-1.414213562373095*dxv[1]*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_1x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.5*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.5*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[3] = 0.5*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[4] += 1.224744871391589*Ghat[1]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[6] += 1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 1.224744871391589*Ghat[3]*rdv2; 
  out[8] += -1.58113883008419*Ghat[0]*rdv2; 
  out[9] += -1.58113883008419*Ghat[1]*rdv2; 
  out[10] += -1.58113883008419*Ghat[2]*rdv2; 
  out[11] += -1.58113883008419*Ghat[3]*rdv2; 

  } 
} 
