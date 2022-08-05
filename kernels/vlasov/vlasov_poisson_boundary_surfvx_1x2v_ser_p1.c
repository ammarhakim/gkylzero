#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_poisson_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  double alpha[4] = {0.0}; 

  alpha[0] = -2.449489742783178*phi[1]*dx10; 

  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  if (alpha[0] > 0) { 
    fUpwindQuad[0] = ser_3x_p1_surfx2_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_3x_p1_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_3x_p1_surfx2_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_3x_p1_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if (alpha[0] > 0) { 
    fUpwindQuad[2] = ser_3x_p1_surfx2_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = ser_3x_p1_surfx2_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_3x_p1_surfx2_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = ser_3x_p1_surfx2_eval_quad_node_3_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alpha[0]*fUpwind[1]; 
  Ghat[2] = 0.5*alpha[0]*fUpwind[2]; 
  Ghat[3] = 0.5*alpha[0]*fUpwind[3]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += -0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += -0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10; 

  } else { 

  if (alpha[0] > 0) { 
    fUpwindQuad[0] = ser_3x_p1_surfx2_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_3x_p1_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_3x_p1_surfx2_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_3x_p1_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if (alpha[0] > 0) { 
    fUpwindQuad[2] = ser_3x_p1_surfx2_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = ser_3x_p1_surfx2_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_3x_p1_surfx2_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = ser_3x_p1_surfx2_eval_quad_node_3_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alpha[0]*fUpwind[1]; 
  Ghat[2] = 0.5*alpha[0]*fUpwind[2]; 
  Ghat[3] = 0.5*alpha[0]*fUpwind[3]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -1.224744871391589*Ghat[0]*dv10; 
  out[3] += 0.7071067811865475*Ghat[2]*dv10; 
  out[4] += -1.224744871391589*Ghat[1]*dv10; 
  out[5] += 0.7071067811865475*Ghat[3]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.224744871391589*Ghat[3]*dv10; 

  } 
} 
