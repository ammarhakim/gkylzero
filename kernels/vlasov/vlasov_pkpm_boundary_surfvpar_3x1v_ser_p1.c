#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_3x1v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_3x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_boundary_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, 
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
  const double dx1 = 2.0/dxv[1]; 
  const double dx2 = 2.0/dxv[2]; 
  const double dv1par = 2.0/dxv[3]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  double alphaSurf[8] = {0.0}; 
  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};;
  double Ghat[8] = {0.0}; 

  if (edge == -1) { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf[4] = (-1.0*bb_grad_u[4]*wvpar)-0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf[5] = (-1.0*bb_grad_u[5]*wvpar)-0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf[6] = (-1.0*bb_grad_u[6]*wvpar)-0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf[7] = (-1.0*bb_grad_u[7]*wvpar)-0.5*bb_grad_u[7]*dvpar+p_force[7]; 

  if ((-0.3535533905932737*alphaSurf[7])+0.3535533905932737*(alphaSurf[6]+alphaSurf[5]+alphaSurf[4])-0.3535533905932737*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(fSkin); 
  } else { 
    fUpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(fEdge); 
  } 
  if (0.3535533905932737*alphaSurf[7]-0.3535533905932737*(alphaSurf[6]+alphaSurf[5])+0.3535533905932737*(alphaSurf[4]+alphaSurf[3])-0.3535533905932737*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(fSkin); 
  } else { 
    fUpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(fEdge); 
  } 
  if (0.3535533905932737*alphaSurf[7]-0.3535533905932737*alphaSurf[6]+0.3535533905932737*alphaSurf[5]-0.3535533905932737*(alphaSurf[4]+alphaSurf[3])+0.3535533905932737*alphaSurf[2]-0.3535533905932737*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if ((-0.3535533905932737*alphaSurf[7])+0.3535533905932737*alphaSurf[6]-0.3535533905932737*(alphaSurf[5]+alphaSurf[4])+0.3535533905932737*(alphaSurf[3]+alphaSurf[2])-0.3535533905932737*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(fEdge); 
  } 
  if (0.3535533905932737*(alphaSurf[7]+alphaSurf[6])-0.3535533905932737*(alphaSurf[5]+alphaSurf[4]+alphaSurf[3]+alphaSurf[2])+0.3535533905932737*(alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(fSkin); 
  } else { 
    fUpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(fEdge); 
  } 
  if ((-0.3535533905932737*(alphaSurf[7]+alphaSurf[6]))+0.3535533905932737*alphaSurf[5]-0.3535533905932737*alphaSurf[4]+0.3535533905932737*alphaSurf[3]-0.3535533905932737*alphaSurf[2]+0.3535533905932737*(alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if ((-0.3535533905932737*(alphaSurf[7]+alphaSurf[6]+alphaSurf[5]))+0.3535533905932737*alphaSurf[4]-0.3535533905932737*alphaSurf[3]+0.3535533905932737*(alphaSurf[2]+alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(fSkin); 
  } else { 
    fUpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(fEdge); 
  } 
  if (0.3535533905932737*(alphaSurf[7]+alphaSurf[6]+alphaSurf[5]+alphaSurf[4]+alphaSurf[3]+alphaSurf[2]+alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaSurf[7]*fUpwind[7]+0.3535533905932737*alphaSurf[6]*fUpwind[6]+0.3535533905932737*alphaSurf[5]*fUpwind[5]+0.3535533905932737*alphaSurf[4]*fUpwind[4]+0.3535533905932737*alphaSurf[3]*fUpwind[3]+0.3535533905932737*alphaSurf[2]*fUpwind[2]+0.3535533905932737*alphaSurf[1]*fUpwind[1]+0.3535533905932737*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alphaSurf[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alphaSurf[7]+0.3535533905932737*alphaSurf[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaSurf[5]+0.3535533905932737*alphaSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaSurf[4]+0.3535533905932737*alphaSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.3535533905932737*alphaSurf[5]*fUpwind[7]+0.3535533905932737*fUpwind[5]*alphaSurf[7]+0.3535533905932737*alphaSurf[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alphaSurf[6]+0.3535533905932737*alphaSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaSurf[4]+0.3535533905932737*alphaSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaSurf[2]; 
  Ghat[3] = 0.3535533905932737*alphaSurf[4]*fUpwind[7]+0.3535533905932737*fUpwind[4]*alphaSurf[7]+0.3535533905932737*alphaSurf[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alphaSurf[6]+0.3535533905932737*alphaSurf[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alphaSurf[5]+0.3535533905932737*alphaSurf[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alphaSurf[3]; 
  Ghat[4] = 0.3535533905932737*alphaSurf[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alphaSurf[7]+0.3535533905932737*alphaSurf[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alphaSurf[6]+0.3535533905932737*alphaSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaSurf[4]+0.3535533905932737*alphaSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaSurf[2]; 
  Ghat[5] = 0.3535533905932737*alphaSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaSurf[7]+0.3535533905932737*alphaSurf[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alphaSurf[6]+0.3535533905932737*alphaSurf[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alphaSurf[5]+0.3535533905932737*alphaSurf[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alphaSurf[3]; 
  Ghat[6] = 0.3535533905932737*alphaSurf[1]*fUpwind[7]+0.3535533905932737*fUpwind[1]*alphaSurf[7]+0.3535533905932737*alphaSurf[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alphaSurf[6]+0.3535533905932737*alphaSurf[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alphaSurf[5]+0.3535533905932737*alphaSurf[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alphaSurf[3]; 
  Ghat[7] = 0.3535533905932737*alphaSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaSurf[7]+0.3535533905932737*alphaSurf[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alphaSurf[6]+0.3535533905932737*alphaSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alphaSurf[5]+0.3535533905932737*alphaSurf[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alphaSurf[4]; 

  out[0] += -0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += -0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += -0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += -0.7071067811865475*Ghat[3]*dv1par; 
  out[4] += -1.224744871391589*Ghat[0]*dv1par; 
  out[5] += -0.7071067811865475*Ghat[4]*dv1par; 
  out[6] += -0.7071067811865475*Ghat[5]*dv1par; 
  out[7] += -0.7071067811865475*Ghat[6]*dv1par; 
  out[8] += -1.224744871391589*Ghat[1]*dv1par; 
  out[9] += -1.224744871391589*Ghat[2]*dv1par; 
  out[10] += -1.224744871391589*Ghat[3]*dv1par; 
  out[11] += -0.7071067811865475*Ghat[7]*dv1par; 
  out[12] += -1.224744871391589*Ghat[4]*dv1par; 
  out[13] += -1.224744871391589*Ghat[5]*dv1par; 
  out[14] += -1.224744871391589*Ghat[6]*dv1par; 
  out[15] += -1.224744871391589*Ghat[7]*dv1par; 
  out[16] += -1.58113883008419*Ghat[0]*dv1par; 
  out[17] += -1.58113883008419*Ghat[1]*dv1par; 
  out[18] += -1.58113883008419*Ghat[2]*dv1par; 
  out[19] += -1.58113883008419*Ghat[3]*dv1par; 
  out[20] += -1.58113883008419*Ghat[4]*dv1par; 
  out[21] += -1.58113883008419*Ghat[5]*dv1par; 
  out[22] += -1.58113883008419*Ghat[6]*dv1par; 
  out[23] += -1.58113883008419*Ghat[7]*dv1par; 

  } else { 

  alphaSurf[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf[4] = (-1.0*bb_grad_u[4]*wvpar)+0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf[5] = (-1.0*bb_grad_u[5]*wvpar)+0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf[6] = (-1.0*bb_grad_u[6]*wvpar)+0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf[7] = (-1.0*bb_grad_u[7]*wvpar)+0.5*bb_grad_u[7]*dvpar+p_force[7]; 

  if ((-0.3535533905932737*alphaSurf[7])+0.3535533905932737*(alphaSurf[6]+alphaSurf[5]+alphaSurf[4])-0.3535533905932737*(alphaSurf[3]+alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(fEdge); 
  } else { 
    fUpwindQuad[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(fSkin); 
  } 
  if (0.3535533905932737*alphaSurf[7]-0.3535533905932737*(alphaSurf[6]+alphaSurf[5])+0.3535533905932737*(alphaSurf[4]+alphaSurf[3])-0.3535533905932737*(alphaSurf[2]+alphaSurf[1])+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(fEdge); 
  } else { 
    fUpwindQuad[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(fSkin); 
  } 
  if (0.3535533905932737*alphaSurf[7]-0.3535533905932737*alphaSurf[6]+0.3535533905932737*alphaSurf[5]-0.3535533905932737*(alphaSurf[4]+alphaSurf[3])+0.3535533905932737*alphaSurf[2]-0.3535533905932737*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if ((-0.3535533905932737*alphaSurf[7])+0.3535533905932737*alphaSurf[6]-0.3535533905932737*(alphaSurf[5]+alphaSurf[4])+0.3535533905932737*(alphaSurf[3]+alphaSurf[2])-0.3535533905932737*alphaSurf[1]+0.3535533905932737*alphaSurf[0] > 0) { 
    fUpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(fSkin); 
  } 
  if (0.3535533905932737*(alphaSurf[7]+alphaSurf[6])-0.3535533905932737*(alphaSurf[5]+alphaSurf[4]+alphaSurf[3]+alphaSurf[2])+0.3535533905932737*(alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(fEdge); 
  } else { 
    fUpwindQuad[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(fSkin); 
  } 
  if ((-0.3535533905932737*(alphaSurf[7]+alphaSurf[6]))+0.3535533905932737*alphaSurf[5]-0.3535533905932737*alphaSurf[4]+0.3535533905932737*alphaSurf[3]-0.3535533905932737*alphaSurf[2]+0.3535533905932737*(alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if ((-0.3535533905932737*(alphaSurf[7]+alphaSurf[6]+alphaSurf[5]))+0.3535533905932737*alphaSurf[4]-0.3535533905932737*alphaSurf[3]+0.3535533905932737*(alphaSurf[2]+alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(fEdge); 
  } else { 
    fUpwindQuad[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(fSkin); 
  } 
  if (0.3535533905932737*(alphaSurf[7]+alphaSurf[6]+alphaSurf[5]+alphaSurf[4]+alphaSurf[3]+alphaSurf[2]+alphaSurf[1]+alphaSurf[0]) > 0) { 
    fUpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaSurf[7]*fUpwind[7]+0.3535533905932737*alphaSurf[6]*fUpwind[6]+0.3535533905932737*alphaSurf[5]*fUpwind[5]+0.3535533905932737*alphaSurf[4]*fUpwind[4]+0.3535533905932737*alphaSurf[3]*fUpwind[3]+0.3535533905932737*alphaSurf[2]*fUpwind[2]+0.3535533905932737*alphaSurf[1]*fUpwind[1]+0.3535533905932737*alphaSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alphaSurf[6]*fUpwind[7]+0.3535533905932737*fUpwind[6]*alphaSurf[7]+0.3535533905932737*alphaSurf[3]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaSurf[5]+0.3535533905932737*alphaSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaSurf[4]+0.3535533905932737*alphaSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaSurf[1]; 
  Ghat[2] = 0.3535533905932737*alphaSurf[5]*fUpwind[7]+0.3535533905932737*fUpwind[5]*alphaSurf[7]+0.3535533905932737*alphaSurf[3]*fUpwind[6]+0.3535533905932737*fUpwind[3]*alphaSurf[6]+0.3535533905932737*alphaSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaSurf[4]+0.3535533905932737*alphaSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaSurf[2]; 
  Ghat[3] = 0.3535533905932737*alphaSurf[4]*fUpwind[7]+0.3535533905932737*fUpwind[4]*alphaSurf[7]+0.3535533905932737*alphaSurf[2]*fUpwind[6]+0.3535533905932737*fUpwind[2]*alphaSurf[6]+0.3535533905932737*alphaSurf[1]*fUpwind[5]+0.3535533905932737*fUpwind[1]*alphaSurf[5]+0.3535533905932737*alphaSurf[0]*fUpwind[3]+0.3535533905932737*fUpwind[0]*alphaSurf[3]; 
  Ghat[4] = 0.3535533905932737*alphaSurf[3]*fUpwind[7]+0.3535533905932737*fUpwind[3]*alphaSurf[7]+0.3535533905932737*alphaSurf[5]*fUpwind[6]+0.3535533905932737*fUpwind[5]*alphaSurf[6]+0.3535533905932737*alphaSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaSurf[4]+0.3535533905932737*alphaSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaSurf[2]; 
  Ghat[5] = 0.3535533905932737*alphaSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaSurf[7]+0.3535533905932737*alphaSurf[4]*fUpwind[6]+0.3535533905932737*fUpwind[4]*alphaSurf[6]+0.3535533905932737*alphaSurf[0]*fUpwind[5]+0.3535533905932737*fUpwind[0]*alphaSurf[5]+0.3535533905932737*alphaSurf[1]*fUpwind[3]+0.3535533905932737*fUpwind[1]*alphaSurf[3]; 
  Ghat[6] = 0.3535533905932737*alphaSurf[1]*fUpwind[7]+0.3535533905932737*fUpwind[1]*alphaSurf[7]+0.3535533905932737*alphaSurf[0]*fUpwind[6]+0.3535533905932737*fUpwind[0]*alphaSurf[6]+0.3535533905932737*alphaSurf[4]*fUpwind[5]+0.3535533905932737*fUpwind[4]*alphaSurf[5]+0.3535533905932737*alphaSurf[2]*fUpwind[3]+0.3535533905932737*fUpwind[2]*alphaSurf[3]; 
  Ghat[7] = 0.3535533905932737*alphaSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaSurf[7]+0.3535533905932737*alphaSurf[1]*fUpwind[6]+0.3535533905932737*fUpwind[1]*alphaSurf[6]+0.3535533905932737*alphaSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[2]*alphaSurf[5]+0.3535533905932737*alphaSurf[3]*fUpwind[4]+0.3535533905932737*fUpwind[3]*alphaSurf[4]; 

  out[0] += 0.7071067811865475*Ghat[0]*dv1par; 
  out[1] += 0.7071067811865475*Ghat[1]*dv1par; 
  out[2] += 0.7071067811865475*Ghat[2]*dv1par; 
  out[3] += 0.7071067811865475*Ghat[3]*dv1par; 
  out[4] += -1.224744871391589*Ghat[0]*dv1par; 
  out[5] += 0.7071067811865475*Ghat[4]*dv1par; 
  out[6] += 0.7071067811865475*Ghat[5]*dv1par; 
  out[7] += 0.7071067811865475*Ghat[6]*dv1par; 
  out[8] += -1.224744871391589*Ghat[1]*dv1par; 
  out[9] += -1.224744871391589*Ghat[2]*dv1par; 
  out[10] += -1.224744871391589*Ghat[3]*dv1par; 
  out[11] += 0.7071067811865475*Ghat[7]*dv1par; 
  out[12] += -1.224744871391589*Ghat[4]*dv1par; 
  out[13] += -1.224744871391589*Ghat[5]*dv1par; 
  out[14] += -1.224744871391589*Ghat[6]*dv1par; 
  out[15] += -1.224744871391589*Ghat[7]*dv1par; 
  out[16] += 1.58113883008419*Ghat[0]*dv1par; 
  out[17] += 1.58113883008419*Ghat[1]*dv1par; 
  out[18] += 1.58113883008419*Ghat[2]*dv1par; 
  out[19] += 1.58113883008419*Ghat[3]*dv1par; 
  out[20] += 1.58113883008419*Ghat[4]*dv1par; 
  out[21] += 1.58113883008419*Ghat[5]*dv1par; 
  out[22] += 1.58113883008419*Ghat[6]*dv1par; 
  out[23] += 1.58113883008419*Ghat[7]*dv1par; 

  } 
} 
