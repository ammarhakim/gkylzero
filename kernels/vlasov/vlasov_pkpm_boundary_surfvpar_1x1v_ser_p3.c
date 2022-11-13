#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p3_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p3_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_1x1v_ser_p3(const double *w, const double *dxv, 
     const double *bb_grad_u, const double *p_force, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_force:   total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.
  // bb_grad_u: bb : grad(u).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:       Incremented distribution function in center cell.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  double alphaSurf[4] = {0.0}; 
  double fUpwindQuad[4] = {0.0};
  double fUpwind[4] = {0.0};;
  double Ghat[4] = {0.0}; 

  if (edge == -1) { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force[3]; 

  if ((-0.5701294036773671*alphaSurf[3])+0.9681844646844029*alphaSurf[2]-1.054672281193885*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = ser_2x_p3_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_2x_p3_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.7702725556588816*alphaSurf[3]-0.5164305132317774*alphaSurf[2]-0.4163900395009129*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[1] = ser_2x_p3_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_2x_p3_surfx2_eval_quad_node_1_l(fEdge); 
  } 
  if ((-0.7702725556588816*alphaSurf[3])-0.5164305132317774*alphaSurf[2]+0.4163900395009129*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = ser_2x_p3_surfx2_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_2x_p3_surfx2_eval_quad_node_2_l(fEdge); 
  } 
  if (0.5701294036773671*alphaSurf[3]+0.9681844646844029*alphaSurf[2]+1.054672281193885*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[3] = ser_2x_p3_surfx2_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_2x_p3_surfx2_eval_quad_node_3_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[3]*fUpwind[3]+0.7071067811865475*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6210590034081186*alphaSurf[2]*fUpwind[3]+0.6210590034081186*fUpwind[2]*alphaSurf[3]+0.6324555320336759*alphaSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaSurf[2]+0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.421637021355784*alphaSurf[3]*fUpwind[3]+0.6210590034081186*alphaSurf[1]*fUpwind[3]+0.6210590034081186*fUpwind[1]*alphaSurf[3]+0.4517539514526256*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaSurf[2]+0.6324555320336759*alphaSurf[1]*fUpwind[1]; 
  Ghat[3] = 0.421637021355784*alphaSurf[2]*fUpwind[3]+0.7071067811865475*alphaSurf[0]*fUpwind[3]+0.421637021355784*fUpwind[2]*alphaSurf[3]+0.7071067811865475*fUpwind[0]*alphaSurf[3]+0.6210590034081186*alphaSurf[1]*fUpwind[2]+0.6210590034081186*fUpwind[1]*alphaSurf[2]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += -0.7071067811865475*Ghat[2]*dv1par; 
  out[5] += -1.58113883008419*Ghat[0]*dv1par; 
  out[6] += -1.224744871391589*Ghat[2]*dv1par; 
  out[7] += -1.58113883008419*Ghat[1]*dv1par; 
  out[8] += -0.7071067811865475*Ghat[3]*dv1par; 
  out[9] += -1.870828693386971*Ghat[0]*dv1par; 
  out[10] += -1.224744871391589*Ghat[3]*dv1par; 
  out[11] += -1.870828693386971*Ghat[1]*dv1par; 

  } else { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force[3]; 

  if ((-0.5701294036773671*alphaSurf[3])+0.9681844646844029*alphaSurf[2]-1.054672281193885*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = ser_2x_p3_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x_p3_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.7702725556588816*alphaSurf[3]-0.5164305132317774*alphaSurf[2]-0.4163900395009129*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[1] = ser_2x_p3_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x_p3_surfx2_eval_quad_node_1_l(fSkin); 
  } 
  if ((-0.7702725556588816*alphaSurf[3])-0.5164305132317774*alphaSurf[2]+0.4163900395009129*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = ser_2x_p3_surfx2_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_2x_p3_surfx2_eval_quad_node_2_l(fSkin); 
  } 
  if (0.5701294036773671*alphaSurf[3]+0.9681844646844029*alphaSurf[2]+1.054672281193885*alphaSurf[1]+0.7071067811865475*alphaSurf[0] > 0) { 
    fUpwindQuad[3] = ser_2x_p3_surfx2_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_2x_p3_surfx2_eval_quad_node_3_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p3_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[3]*fUpwind[3]+0.7071067811865475*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.6210590034081186*alphaSurf[2]*fUpwind[3]+0.6210590034081186*fUpwind[2]*alphaSurf[3]+0.6324555320336759*alphaSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaSurf[2]+0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.421637021355784*alphaSurf[3]*fUpwind[3]+0.6210590034081186*alphaSurf[1]*fUpwind[3]+0.6210590034081186*fUpwind[1]*alphaSurf[3]+0.4517539514526256*alphaSurf[2]*fUpwind[2]+0.7071067811865475*alphaSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaSurf[2]+0.6324555320336759*alphaSurf[1]*fUpwind[1]; 
  Ghat[3] = 0.421637021355784*alphaSurf[2]*fUpwind[3]+0.7071067811865475*alphaSurf[0]*fUpwind[3]+0.421637021355784*fUpwind[2]*alphaSurf[3]+0.7071067811865475*fUpwind[0]*alphaSurf[3]+0.6210590034081186*alphaSurf[1]*fUpwind[2]+0.6210590034081186*fUpwind[1]*alphaSurf[2]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += 0.7071067811865475*Ghat[2]*dv1par; 
  out[5] += 1.58113883008419*Ghat[0]*dv1par; 
  out[6] += -1.224744871391589*Ghat[2]*dv1par; 
  out[7] += 1.58113883008419*Ghat[1]*dv1par; 
  out[8] += 0.7071067811865475*Ghat[3]*dv1par; 
  out[9] += -1.870828693386971*Ghat[0]*dv1par; 
  out[10] += -1.224744871391589*Ghat[3]*dv1par; 
  out[11] += -1.870828693386971*Ghat[1]*dv1par; 

  } 
} 
