#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_2x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_gkhyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfvpar_2x2v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:     cell-center coordinates. 
  // dxv[4]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[8]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 


  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};;
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[2]+dxv[2])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[2]+dxv[2])-2.0*nuUSum[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(2.0*nuSum[2]*w[2]-2.0*nuUSum[2]+dxv[2]*nuSum[2]); 
  alphaDrSurf[4] = -0.7071067811865475*(2.0*nuUSum[3]+((-2.0*w[2])-1.0*dxv[2])*nuSum[3]); 

  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = gkhyb_2x2v_p1_surfx3_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = gkhyb_2x2v_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[0] = gkhyb_2x2v_p1_surfx3_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = gkhyb_2x2v_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = gkhyb_2x2v_p1_surfx3_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = gkhyb_2x2v_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[2] = gkhyb_2x2v_p1_surfx3_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = gkhyb_2x2v_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = gkhyb_2x2v_p1_surfx3_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = gkhyb_2x2v_p1_surfx3_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[4] = gkhyb_2x2v_p1_surfx3_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = gkhyb_2x2v_p1_surfx3_eval_quad_node_5_l(fEdge); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = gkhyb_2x2v_p1_surfx3_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = gkhyb_2x2v_p1_surfx3_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[6] = gkhyb_2x2v_p1_surfx3_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = gkhyb_2x2v_p1_surfx3_eval_quad_node_7_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_2x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alphaDrSurf[4]*fUpwind[4]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alphaDrSurf[2]*fUpwind[4]+fUpwind[2]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[4]+fUpwind[1]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.3535533905932737*(alphaDrSurf[4]*fUpwind[7]+alphaDrSurf[2]*fUpwind[6]+alphaDrSurf[1]*fUpwind[5]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[4] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[4]+fUpwind[0]*alphaDrSurf[4]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[5] = 0.3535533905932737*(alphaDrSurf[2]*fUpwind[7]+alphaDrSurf[4]*fUpwind[6]+alphaDrSurf[0]*fUpwind[5]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[6] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[7]+alphaDrSurf[0]*fUpwind[6]+alphaDrSurf[4]*fUpwind[5]+alphaDrSurf[2]*fUpwind[3]); 
  Ghat[7] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[7]+alphaDrSurf[1]*fUpwind[6]+alphaDrSurf[2]*fUpwind[5]+fUpwind[3]*alphaDrSurf[4]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 1.224744871391589*Ghat[1]*rdv2; 
  out[7] += 1.224744871391589*Ghat[2]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 1.224744871391589*Ghat[4]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat[7]*rdv2; 
  out[16] += 1.58113883008419*Ghat[0]*rdv2; 
  out[17] += 1.58113883008419*Ghat[1]*rdv2; 
  out[18] += 1.58113883008419*Ghat[2]*rdv2; 
  out[19] += 1.58113883008419*Ghat[3]*rdv2; 
  out[20] += 1.58113883008419*Ghat[4]*rdv2; 
  out[21] += 1.58113883008419*Ghat[5]*rdv2; 
  out[22] += 1.58113883008419*Ghat[6]*rdv2; 
  out[23] += 1.58113883008419*Ghat[7]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[2]-1.0*dxv[2])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[2]-1.0*dxv[2])-2.0*nuUSum[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(2.0*nuSum[2]*w[2]-2.0*nuUSum[2]-1.0*dxv[2]*nuSum[2]); 
  alphaDrSurf[4] = -0.7071067811865475*(2.0*nuUSum[3]+(dxv[2]-2.0*w[2])*nuSum[3]); 

  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = gkhyb_2x2v_p1_surfx3_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = gkhyb_2x2v_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[0] = gkhyb_2x2v_p1_surfx3_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = gkhyb_2x2v_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = gkhyb_2x2v_p1_surfx3_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = gkhyb_2x2v_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[2] = gkhyb_2x2v_p1_surfx3_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = gkhyb_2x2v_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = gkhyb_2x2v_p1_surfx3_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = gkhyb_2x2v_p1_surfx3_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[4] = gkhyb_2x2v_p1_surfx3_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = gkhyb_2x2v_p1_surfx3_eval_quad_node_5_l(fSkin); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = gkhyb_2x2v_p1_surfx3_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = gkhyb_2x2v_p1_surfx3_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[6] = gkhyb_2x2v_p1_surfx3_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = gkhyb_2x2v_p1_surfx3_eval_quad_node_7_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_2x2v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*(alphaDrSurf[4]*fUpwind[4]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.3535533905932737*(alphaDrSurf[2]*fUpwind[4]+fUpwind[2]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[4]+fUpwind[1]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.3535533905932737*(alphaDrSurf[4]*fUpwind[7]+alphaDrSurf[2]*fUpwind[6]+alphaDrSurf[1]*fUpwind[5]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[4] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[4]+fUpwind[0]*alphaDrSurf[4]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[5] = 0.3535533905932737*(alphaDrSurf[2]*fUpwind[7]+alphaDrSurf[4]*fUpwind[6]+alphaDrSurf[0]*fUpwind[5]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[6] = 0.3535533905932737*(alphaDrSurf[1]*fUpwind[7]+alphaDrSurf[0]*fUpwind[6]+alphaDrSurf[4]*fUpwind[5]+alphaDrSurf[2]*fUpwind[3]); 
  Ghat[7] = 0.3535533905932737*(alphaDrSurf[0]*fUpwind[7]+alphaDrSurf[1]*fUpwind[6]+alphaDrSurf[2]*fUpwind[5]+fUpwind[3]*alphaDrSurf[4]); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 1.224744871391589*Ghat[1]*rdv2; 
  out[7] += 1.224744871391589*Ghat[2]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[9] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 1.224744871391589*Ghat[4]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[13] += 1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 1.224744871391589*Ghat[7]*rdv2; 
  out[16] += -1.58113883008419*Ghat[0]*rdv2; 
  out[17] += -1.58113883008419*Ghat[1]*rdv2; 
  out[18] += -1.58113883008419*Ghat[2]*rdv2; 
  out[19] += -1.58113883008419*Ghat[3]*rdv2; 
  out[20] += -1.58113883008419*Ghat[4]*rdv2; 
  out[21] += -1.58113883008419*Ghat[5]*rdv2; 
  out[22] += -1.58113883008419*Ghat[6]*rdv2; 
  out[23] += -1.58113883008419*Ghat[7]*rdv2; 

  } 
} 
