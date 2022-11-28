#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *bb_grad_u, const double *p_force, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_force:   total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.
  // bb_grad_u: bb : grad(u).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  double alphaSurf_l[3] = {0.0}; 
  alphaSurf_l[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf_l[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf_l[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force[2]; 

  double alphaSurf_r[3] = {0.0}; 
  alphaSurf_r[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf_r[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf_r[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force[2]; 

  double fUpwindQuad_l[3] = {0.0};
  double fUpwindQuad_r[3] = {0.0};
  double fUpwind_l[3] = {0.0};
  double fUpwind_r[3] = {0.0};
  double Ghat_l[3] = {0.0}; 
  double Ghat_r[3] = {0.0}; 

  if (0.6324555320336759*alphaSurf_l[2]-0.9486832980505137*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fc); 
  } 
  if (0.6324555320336759*alphaSurf_r[2]-0.9486832980505137*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(fr); 
  } 
  if (0.7071067811865475*alphaSurf_l[0]-0.7905694150420947*alphaSurf_l[2] > 0) { 
    fUpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fc); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-0.7905694150420947*alphaSurf_r[2] > 0) { 
    fUpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(fr); 
  } 
  if (0.6324555320336759*alphaSurf_l[2]+0.9486832980505137*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fc); 
  } 
  if (0.6324555320336759*alphaSurf_r[2]+0.9486832980505137*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_2x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.7071067811865475*alphaSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.6324555320336759*alphaSurf_l[1]*fUpwind_l[2]+0.6324555320336759*fUpwind_l[1]*alphaSurf_l[2]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.4517539514526256*alphaSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaSurf_l[0]*fUpwind_l[2]+0.7071067811865475*fUpwind_l[0]*alphaSurf_l[2]+0.6324555320336759*alphaSurf_l[1]*fUpwind_l[1]; 

  Ghat_r[0] = 0.7071067811865475*alphaSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.6324555320336759*alphaSurf_r[1]*fUpwind_r[2]+0.6324555320336759*fUpwind_r[1]*alphaSurf_r[2]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.4517539514526256*alphaSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaSurf_r[0]*fUpwind_r[2]+0.7071067811865475*fUpwind_r[0]*alphaSurf_r[2]+0.6324555320336759*alphaSurf_r[1]*fUpwind_r[1]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv1par; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv1par; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv1par; 
  out[4] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv1par; 
  out[5] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv1par; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv1par; 
  out[7] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv1par; 

} 
