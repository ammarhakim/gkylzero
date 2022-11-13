#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfvpar_2x1v_ser_p2(const double *w, const double *dxv, 
     const double *bb_grad_u, const double *p_force, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_force:   total pressure force = 1/rho (b . div(P) + p_perp div(b)) for Euler PKPM.
  // bb_grad_u: bb : grad(u).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  double alphaSurf_l[8] = {0.0}; 
  alphaSurf_l[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf_l[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf_l[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf_l[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf_l[4] = (-1.0*bb_grad_u[4]*wvpar)+0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf_l[5] = (-1.0*bb_grad_u[5]*wvpar)+0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf_l[6] = (-1.0*bb_grad_u[6]*wvpar)+0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf_l[7] = (-1.0*bb_grad_u[7]*wvpar)+0.5*bb_grad_u[7]*dvpar+p_force[7]; 

  double alphaSurf_r[8] = {0.0}; 
  alphaSurf_r[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force[0]; 
  alphaSurf_r[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force[1]; 
  alphaSurf_r[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force[2]; 
  alphaSurf_r[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force[3]; 
  alphaSurf_r[4] = (-1.0*bb_grad_u[4]*wvpar)-0.5*bb_grad_u[4]*dvpar+p_force[4]; 
  alphaSurf_r[5] = (-1.0*bb_grad_u[5]*wvpar)-0.5*bb_grad_u[5]*dvpar+p_force[5]; 
  alphaSurf_r[6] = (-1.0*bb_grad_u[6]*wvpar)-0.5*bb_grad_u[6]*dvpar+p_force[6]; 
  alphaSurf_r[7] = (-1.0*bb_grad_u[7]*wvpar)-0.5*bb_grad_u[7]*dvpar+p_force[7]; 

  double fUpwindQuad_l[9] = {0.0};
  double fUpwindQuad_r[9] = {0.0};
  double fUpwind_l[8] = {0.0};
  double fUpwind_r[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 

  if ((-0.5999999999999995*alphaSurf_l[7])-0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*(alphaSurf_l[5]+alphaSurf_l[4])+0.9*alphaSurf_l[3]-0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx3_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p2_surfx3_eval_quad_node_0_l(fc); 
  } 
  if ((-0.5999999999999995*alphaSurf_r[7])-0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*(alphaSurf_r[5]+alphaSurf_r[4])+0.9*alphaSurf_r[3]-0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_3x_p2_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_3x_p2_surfx3_eval_quad_node_0_l(fr); 
  } 
  if (0.75*alphaSurf_l[7]-0.5590169943749475*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx3_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p2_surfx3_eval_quad_node_1_l(fc); 
  } 
  if (0.75*alphaSurf_r[7]-0.5590169943749475*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_3x_p2_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_3x_p2_surfx3_eval_quad_node_1_l(fr); 
  } 
  if ((-0.5999999999999995*alphaSurf_l[7])+0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*(alphaSurf_l[5]+alphaSurf_l[4])-0.9*alphaSurf_l[3]+0.6708203932499369*alphaSurf_l[2]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx3_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_3x_p2_surfx3_eval_quad_node_2_l(fc); 
  } 
  if ((-0.5999999999999995*alphaSurf_r[7])+0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*(alphaSurf_r[5]+alphaSurf_r[4])-0.9*alphaSurf_r[3]+0.6708203932499369*alphaSurf_r[2]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_3x_p2_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_3x_p2_surfx3_eval_quad_node_2_l(fr); 
  } 
  if (0.75*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[5]-0.5590169943749475*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_3x_p2_surfx3_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_3x_p2_surfx3_eval_quad_node_3_l(fc); 
  } 
  if (0.75*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[5]-0.5590169943749475*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_3x_p2_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_3x_p2_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (0.5*alphaSurf_l[0]-0.5590169943749475*(alphaSurf_l[5]+alphaSurf_l[4]) > 0) { 
    fUpwindQuad_l[4] = ser_3x_p2_surfx3_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_3x_p2_surfx3_eval_quad_node_4_l(fc); 
  } 
  if (0.5*alphaSurf_r[0]-0.5590169943749475*(alphaSurf_r[5]+alphaSurf_r[4]) > 0) { 
    fUpwindQuad_r[4] = ser_3x_p2_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_3x_p2_surfx3_eval_quad_node_4_l(fr); 
  } 
  if ((-0.75*alphaSurf_l[6])+0.4472135954999579*alphaSurf_l[5]-0.5590169943749475*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = ser_3x_p2_surfx3_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_3x_p2_surfx3_eval_quad_node_5_l(fc); 
  } 
  if ((-0.75*alphaSurf_r[6])+0.4472135954999579*alphaSurf_r[5]-0.5590169943749475*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = ser_3x_p2_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_3x_p2_surfx3_eval_quad_node_5_l(fr); 
  } 
  if (0.5999999999999995*alphaSurf_l[7]-0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*(alphaSurf_l[5]+alphaSurf_l[4])-0.9*alphaSurf_l[3]-0.6708203932499369*alphaSurf_l[2]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[6] = ser_3x_p2_surfx3_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_3x_p2_surfx3_eval_quad_node_6_l(fc); 
  } 
  if (0.5999999999999995*alphaSurf_r[7]-0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*(alphaSurf_r[5]+alphaSurf_r[4])-0.9*alphaSurf_r[3]-0.6708203932499369*alphaSurf_r[2]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[6] = ser_3x_p2_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_3x_p2_surfx3_eval_quad_node_6_l(fr); 
  } 
  if ((-0.75*alphaSurf_l[7])-0.5590169943749475*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = ser_3x_p2_surfx3_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_3x_p2_surfx3_eval_quad_node_7_l(fc); 
  } 
  if ((-0.75*alphaSurf_r[7])-0.5590169943749475*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = ser_3x_p2_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_3x_p2_surfx3_eval_quad_node_7_l(fr); 
  } 
  if (0.5999999999999995*alphaSurf_l[7]+0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*(alphaSurf_l[5]+alphaSurf_l[4])+0.9*alphaSurf_l[3]+0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    fUpwindQuad_l[8] = ser_3x_p2_surfx3_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_3x_p2_surfx3_eval_quad_node_8_l(fc); 
  } 
  if (0.5999999999999995*alphaSurf_r[7]+0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*(alphaSurf_r[5]+alphaSurf_r[4])+0.9*alphaSurf_r[3]+0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    fUpwindQuad_r[8] = ser_3x_p2_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_3x_p2_surfx3_eval_quad_node_8_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.5*alphaSurf_l[7]*fUpwind_l[7]+0.5*alphaSurf_l[6]*fUpwind_l[6]+0.5*alphaSurf_l[5]*fUpwind_l[5]+0.5*alphaSurf_l[4]*fUpwind_l[4]+0.5*alphaSurf_l[3]*fUpwind_l[3]+0.5*alphaSurf_l[2]*fUpwind_l[2]+0.5*alphaSurf_l[1]*fUpwind_l[1]+0.5*alphaSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.5000000000000001*alphaSurf_l[5]*fUpwind_l[7]+0.5000000000000001*fUpwind_l[5]*alphaSurf_l[7]+0.447213595499958*alphaSurf_l[3]*fUpwind_l[6]+0.447213595499958*fUpwind_l[3]*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alphaSurf_l[4]+0.5*alphaSurf_l[2]*fUpwind_l[3]+0.5*fUpwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaSurf_l[1]; 
  Ghat_l[2] = 0.447213595499958*alphaSurf_l[3]*fUpwind_l[7]+0.447213595499958*fUpwind_l[3]*alphaSurf_l[7]+0.5000000000000001*alphaSurf_l[4]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[4]*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[2]*fUpwind_l[5]+0.4472135954999579*fUpwind_l[2]*alphaSurf_l[5]+0.5*alphaSurf_l[1]*fUpwind_l[3]+0.5*fUpwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*fUpwind_l[2]+0.5*fUpwind_l[0]*alphaSurf_l[2]; 
  Ghat_l[3] = 0.4*alphaSurf_l[6]*fUpwind_l[7]+0.447213595499958*alphaSurf_l[2]*fUpwind_l[7]+0.4*fUpwind_l[6]*alphaSurf_l[7]+0.447213595499958*fUpwind_l[2]*alphaSurf_l[7]+0.447213595499958*alphaSurf_l[1]*fUpwind_l[6]+0.447213595499958*fUpwind_l[1]*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[5]+0.4472135954999579*fUpwind_l[3]*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[3]*alphaSurf_l[4]+0.5*alphaSurf_l[0]*fUpwind_l[3]+0.5*fUpwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*fUpwind_l[2]+0.5*fUpwind_l[1]*alphaSurf_l[2]; 
  Ghat_l[4] = 0.4472135954999579*alphaSurf_l[7]*fUpwind_l[7]+0.31943828249997*alphaSurf_l[6]*fUpwind_l[6]+0.5000000000000001*alphaSurf_l[2]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alphaSurf_l[6]+0.31943828249997*alphaSurf_l[4]*fUpwind_l[4]+0.5*alphaSurf_l[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[3]+0.4472135954999579*alphaSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[5] = 0.31943828249997*alphaSurf_l[7]*fUpwind_l[7]+0.5000000000000001*alphaSurf_l[1]*fUpwind_l[7]+0.5000000000000001*fUpwind_l[1]*alphaSurf_l[7]+0.4472135954999579*alphaSurf_l[6]*fUpwind_l[6]+0.31943828249997*alphaSurf_l[5]*fUpwind_l[5]+0.5*alphaSurf_l[0]*fUpwind_l[5]+0.5*fUpwind_l[0]*alphaSurf_l[5]+0.4472135954999579*alphaSurf_l[3]*fUpwind_l[3]+0.4472135954999579*alphaSurf_l[2]*fUpwind_l[2]; 
  Ghat_l[6] = 0.4*alphaSurf_l[3]*fUpwind_l[7]+0.4*fUpwind_l[3]*alphaSurf_l[7]+0.4472135954999579*alphaSurf_l[5]*fUpwind_l[6]+0.31943828249997*alphaSurf_l[4]*fUpwind_l[6]+0.5*alphaSurf_l[0]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[5]*alphaSurf_l[6]+0.31943828249997*fUpwind_l[4]*alphaSurf_l[6]+0.5*fUpwind_l[0]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[2]*fUpwind_l[4]+0.5000000000000001*fUpwind_l[2]*alphaSurf_l[4]+0.447213595499958*alphaSurf_l[1]*fUpwind_l[3]+0.447213595499958*fUpwind_l[1]*alphaSurf_l[3]; 
  Ghat_l[7] = 0.31943828249997*alphaSurf_l[5]*fUpwind_l[7]+0.4472135954999579*alphaSurf_l[4]*fUpwind_l[7]+0.5*alphaSurf_l[0]*fUpwind_l[7]+0.31943828249997*fUpwind_l[5]*alphaSurf_l[7]+0.4472135954999579*fUpwind_l[4]*alphaSurf_l[7]+0.5*fUpwind_l[0]*alphaSurf_l[7]+0.4*alphaSurf_l[3]*fUpwind_l[6]+0.4*fUpwind_l[3]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[1]*fUpwind_l[5]+0.5000000000000001*fUpwind_l[1]*alphaSurf_l[5]+0.447213595499958*alphaSurf_l[2]*fUpwind_l[3]+0.447213595499958*fUpwind_l[2]*alphaSurf_l[3]; 

  Ghat_r[0] = 0.5*alphaSurf_r[7]*fUpwind_r[7]+0.5*alphaSurf_r[6]*fUpwind_r[6]+0.5*alphaSurf_r[5]*fUpwind_r[5]+0.5*alphaSurf_r[4]*fUpwind_r[4]+0.5*alphaSurf_r[3]*fUpwind_r[3]+0.5*alphaSurf_r[2]*fUpwind_r[2]+0.5*alphaSurf_r[1]*fUpwind_r[1]+0.5*alphaSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.5000000000000001*alphaSurf_r[5]*fUpwind_r[7]+0.5000000000000001*fUpwind_r[5]*alphaSurf_r[7]+0.447213595499958*alphaSurf_r[3]*fUpwind_r[6]+0.447213595499958*fUpwind_r[3]*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alphaSurf_r[4]+0.5*alphaSurf_r[2]*fUpwind_r[3]+0.5*fUpwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaSurf_r[1]; 
  Ghat_r[2] = 0.447213595499958*alphaSurf_r[3]*fUpwind_r[7]+0.447213595499958*fUpwind_r[3]*alphaSurf_r[7]+0.5000000000000001*alphaSurf_r[4]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[4]*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[2]*fUpwind_r[5]+0.4472135954999579*fUpwind_r[2]*alphaSurf_r[5]+0.5*alphaSurf_r[1]*fUpwind_r[3]+0.5*fUpwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*fUpwind_r[2]+0.5*fUpwind_r[0]*alphaSurf_r[2]; 
  Ghat_r[3] = 0.4*alphaSurf_r[6]*fUpwind_r[7]+0.447213595499958*alphaSurf_r[2]*fUpwind_r[7]+0.4*fUpwind_r[6]*alphaSurf_r[7]+0.447213595499958*fUpwind_r[2]*alphaSurf_r[7]+0.447213595499958*alphaSurf_r[1]*fUpwind_r[6]+0.447213595499958*fUpwind_r[1]*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[5]+0.4472135954999579*fUpwind_r[3]*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[3]*alphaSurf_r[4]+0.5*alphaSurf_r[0]*fUpwind_r[3]+0.5*fUpwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*fUpwind_r[2]+0.5*fUpwind_r[1]*alphaSurf_r[2]; 
  Ghat_r[4] = 0.4472135954999579*alphaSurf_r[7]*fUpwind_r[7]+0.31943828249997*alphaSurf_r[6]*fUpwind_r[6]+0.5000000000000001*alphaSurf_r[2]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alphaSurf_r[6]+0.31943828249997*alphaSurf_r[4]*fUpwind_r[4]+0.5*alphaSurf_r[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[3]+0.4472135954999579*alphaSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[5] = 0.31943828249997*alphaSurf_r[7]*fUpwind_r[7]+0.5000000000000001*alphaSurf_r[1]*fUpwind_r[7]+0.5000000000000001*fUpwind_r[1]*alphaSurf_r[7]+0.4472135954999579*alphaSurf_r[6]*fUpwind_r[6]+0.31943828249997*alphaSurf_r[5]*fUpwind_r[5]+0.5*alphaSurf_r[0]*fUpwind_r[5]+0.5*fUpwind_r[0]*alphaSurf_r[5]+0.4472135954999579*alphaSurf_r[3]*fUpwind_r[3]+0.4472135954999579*alphaSurf_r[2]*fUpwind_r[2]; 
  Ghat_r[6] = 0.4*alphaSurf_r[3]*fUpwind_r[7]+0.4*fUpwind_r[3]*alphaSurf_r[7]+0.4472135954999579*alphaSurf_r[5]*fUpwind_r[6]+0.31943828249997*alphaSurf_r[4]*fUpwind_r[6]+0.5*alphaSurf_r[0]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[5]*alphaSurf_r[6]+0.31943828249997*fUpwind_r[4]*alphaSurf_r[6]+0.5*fUpwind_r[0]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[2]*fUpwind_r[4]+0.5000000000000001*fUpwind_r[2]*alphaSurf_r[4]+0.447213595499958*alphaSurf_r[1]*fUpwind_r[3]+0.447213595499958*fUpwind_r[1]*alphaSurf_r[3]; 
  Ghat_r[7] = 0.31943828249997*alphaSurf_r[5]*fUpwind_r[7]+0.4472135954999579*alphaSurf_r[4]*fUpwind_r[7]+0.5*alphaSurf_r[0]*fUpwind_r[7]+0.31943828249997*fUpwind_r[5]*alphaSurf_r[7]+0.4472135954999579*fUpwind_r[4]*alphaSurf_r[7]+0.5*fUpwind_r[0]*alphaSurf_r[7]+0.4*alphaSurf_r[3]*fUpwind_r[6]+0.4*fUpwind_r[3]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[1]*fUpwind_r[5]+0.5000000000000001*fUpwind_r[1]*alphaSurf_r[5]+0.447213595499958*alphaSurf_r[2]*fUpwind_r[3]+0.447213595499958*fUpwind_r[2]*alphaSurf_r[3]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv1par; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv1par; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv1par; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv1par; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv1par; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv1par; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv1par; 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv1par; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv1par; 
  out[9] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv1par; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv1par; 
  out[11] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv1par; 
  out[12] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv1par; 
  out[13] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv1par; 
  out[14] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv1par; 
  out[15] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv1par; 
  out[16] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv1par; 
  out[17] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv1par; 
  out[18] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv1par; 
  out[19] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv1par; 

} 
