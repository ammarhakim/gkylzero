#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_hyb_2x1v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_vmap_boundary_surfvx_2x1v_ser_p1(const double *w, const double *dxv, const double *jacob_vel_inv, const double *gamma, const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // jacob_vel_inv[VDIM]: Inverse velocity space Jacobian in each direction (unused in uniform grid simulations).
  // gamma:       Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:        q/m*EM fields.
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2.0/dxv[2]; 
  const double *E0 = &qmem[0]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  double alpha[4] = {0.0}; 

  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  alpha[0] = 1.58113883008419*E0[0]*jacob_vel_inv0[2]+1.224744871391589*E0[0]*jacob_vel_inv0[1]+0.7071067811865475*E0[0]*jacob_vel_inv0[0]; 
  alpha[1] = 1.58113883008419*E0[1]*jacob_vel_inv0[2]+1.224744871391589*E0[1]*jacob_vel_inv0[1]+0.7071067811865475*jacob_vel_inv0[0]*E0[1]; 
  alpha[2] = 1.58113883008419*E0[2]*jacob_vel_inv0[2]+1.224744871391589*jacob_vel_inv0[1]*E0[2]+0.7071067811865475*jacob_vel_inv0[0]*E0[2]; 
  alpha[3] = 1.58113883008419*jacob_vel_inv0[2]*E0[3]+1.224744871391589*jacob_vel_inv0[1]*E0[3]+0.7071067811865475*jacob_vel_inv0[0]*E0[3]; 

  if (0.5*alpha[3]-0.5*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(fEdge); 
  } 
  if ((-0.5*alpha[3])+0.5*alpha[2]-0.5*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(fEdge); 
  } 
  if (0.5*(alpha[1]+alpha[0])-0.5*(alpha[3]+alpha[2]) > 0) { 
    fUpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(fEdge); 
  } 
  if (0.5*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[3]*fUpwind[3]+0.5*alpha[2]*fUpwind[2]+0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alpha[2]*fUpwind[3]+0.5*fUpwind[2]*alpha[3]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.5*alpha[1]*fUpwind[3]+0.5*fUpwind[1]*alpha[3]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[0]*alpha[3]+0.5*alpha[1]*fUpwind[2]+0.5*fUpwind[1]*alpha[2]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -1.224744871391589*Ghat[1]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += -1.58113883008419*Ghat[0]*dv10; 
  out[8] += -1.224744871391589*Ghat[3]*dv10; 
  out[9] += -1.58113883008419*Ghat[1]*dv10; 
  out[10] += -1.58113883008419*Ghat[2]*dv10; 
  out[11] += -1.58113883008419*Ghat[3]*dv10; 

  } else { 

  alpha[0] = 1.58113883008419*E0[0]*jacob_vel_inv0[2]-1.224744871391589*E0[0]*jacob_vel_inv0[1]+0.7071067811865475*E0[0]*jacob_vel_inv0[0]; 
  alpha[1] = 1.58113883008419*E0[1]*jacob_vel_inv0[2]-1.224744871391589*E0[1]*jacob_vel_inv0[1]+0.7071067811865475*jacob_vel_inv0[0]*E0[1]; 
  alpha[2] = 1.58113883008419*E0[2]*jacob_vel_inv0[2]-1.224744871391589*jacob_vel_inv0[1]*E0[2]+0.7071067811865475*jacob_vel_inv0[0]*E0[2]; 
  alpha[3] = 1.58113883008419*jacob_vel_inv0[2]*E0[3]-1.224744871391589*jacob_vel_inv0[1]*E0[3]+0.7071067811865475*jacob_vel_inv0[0]*E0[3]; 

  if (0.5*alpha[3]-0.5*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) { 
    fUpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(fSkin); 
  } 
  if ((-0.5*alpha[3])+0.5*alpha[2]-0.5*alpha[1]+0.5*alpha[0] > 0) { 
    fUpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(fSkin); 
  } 
  if (0.5*(alpha[1]+alpha[0])-0.5*(alpha[3]+alpha[2]) > 0) { 
    fUpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(fSkin); 
  } 
  if (0.5*(alpha[3]+alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.5*alpha[3]*fUpwind[3]+0.5*alpha[2]*fUpwind[2]+0.5*alpha[1]*fUpwind[1]+0.5*alpha[0]*fUpwind[0]; 
  Ghat[1] = 0.5*alpha[2]*fUpwind[3]+0.5*fUpwind[2]*alpha[3]+0.5*alpha[0]*fUpwind[1]+0.5*fUpwind[0]*alpha[1]; 
  Ghat[2] = 0.5*alpha[1]*fUpwind[3]+0.5*fUpwind[1]*alpha[3]+0.5*alpha[0]*fUpwind[2]+0.5*fUpwind[0]*alpha[2]; 
  Ghat[3] = 0.5*alpha[0]*fUpwind[3]+0.5*fUpwind[0]*alpha[3]+0.5*alpha[1]*fUpwind[2]+0.5*fUpwind[1]*alpha[2]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -1.224744871391589*Ghat[1]*dv10; 
  out[6] += -1.224744871391589*Ghat[2]*dv10; 
  out[7] += 1.58113883008419*Ghat[0]*dv10; 
  out[8] += -1.224744871391589*Ghat[3]*dv10; 
  out[9] += 1.58113883008419*Ghat[1]*dv10; 
  out[10] += 1.58113883008419*Ghat[2]*dv10; 
  out[11] += 1.58113883008419*Ghat[3]*dv10; 

  } 
  return 0.;

} 
