#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *bb_grad_u, const double *p_force, 
     const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:     Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // p_force:     total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.
  // bb_grad_u:   bb : grad(u).
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Incremented distribution function in center cell.
  const double dx0 = 2.0/dxv[0]; 
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  double alphaSurf[2] = {0.0}; 
  double fUpwindQuad[2] = {0.0};
  double fUpwind[2] = {0.0};;
  double Ghat[2] = {0.0}; 

  if (edge == -1) { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 

  if (0.7071067811865475*alphaSurf[0]-0.7071067811865475*alphaSurf[1] > 0) { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fEdge); 
  } 
  if (0.7071067811865475*(alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += -1.58113883008419*Ghat[0]*dv1par; 
  out[5] += -1.58113883008419*Ghat[1]*dv1par; 

  } else { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 

  if (0.7071067811865475*alphaSurf[0]-0.7071067811865475*alphaSurf[1] > 0) { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_1x1v_p1_surfx2_eval_quad_node_0_l(fSkin); 
  } 
  if (0.7071067811865475*(alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_1x1v_p1_surfx2_eval_quad_node_1_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.7071067811865475*alphaSurf[1]*fUpwind[1]+0.7071067811865475*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.7071067811865475*alphaSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaSurf[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -1.224744871391589*Ghat[0]*dv1par; 
  out[3] += -1.224744871391589*Ghat[1]*dv1par; 
  out[4] += 1.58113883008419*Ghat[0]*dv1par; 
  out[5] += 1.58113883008419*Ghat[1]*dv1par; 

  } 
} 
