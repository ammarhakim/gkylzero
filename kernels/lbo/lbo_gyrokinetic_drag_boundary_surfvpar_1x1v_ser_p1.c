#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_gkhyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:     cell-center coordinates. 
  // dxv[2]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[2]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 


  double alphaDrSurf[2] = {0.0}; 
  double fUpwindQuad[2] = {0.0};
  double fUpwind[2] = {0.0};;
  double Ghat[2] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]+dxv[1])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*nuUSum[1]+dxv[1]*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_1x1v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.7071067811865475*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += 1.58113883008419*Ghat[0]*rdv2; 
  out[5] += 1.58113883008419*Ghat[1]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]-1.0*dxv[1])-2.0*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*nuUSum[1]-1.0*dxv[1]*nuSum[1]); 

  if (alphaDrSurf[0]-alphaDrSurf[1] < 0) { 
    fUpwindQuad[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  gkhyb_1x1v_p1_vpardir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.7071067811865475*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*Ghat[1]*rdv2; 
  out[4] += -1.58113883008419*Ghat[0]*rdv2; 
  out[5] += -1.58113883008419*Ghat[1]*rdv2; 

  } 
} 
