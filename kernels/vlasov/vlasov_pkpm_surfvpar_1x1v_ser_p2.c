#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfvpar_1x1v_ser_p2(const double *w, const double *dxv, 
     const double *bb_grad_u, const double *p_force, const double *div_b, const double *p_perp_div_b, 
     const double *g_dist_sourcel, const double *g_dist_sourcec, const double *g_dist_sourcer, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:            Cell-center coordinates.
  // dxv[NDIM]:          Cell spacing.
  // p_force:            total pressure force = 1/rho (div(p_parallel b_hat) - p_perp div(b)).
  // bb_grad_u:          bb : grad(u).
  // div_b:              divergence of the magnetic field unit vector.
  // p_perp_div_b:       2*p_perp/rho*div(b).
  // g_dist_sourcel/c/r: 2*T_perp/m*G - T_perp/m*F_0 in left/center/right cells.
  // fl/fc/fr:           Input Distribution function [F_0, T_perp/m G = T_perp/m (F_0 - F_1)] in left/center/right cells.
  // out:                Incremented distribution function in center cell.
  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[8]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[8]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[8]; 

  const double *F_0_sourcel = &g_dist_sourcel[0]; 
  const double *G_1_sourcel = &g_dist_sourcel[8]; 
  const double *F_0_sourcec = &g_dist_sourcec[0]; 
  const double *G_1_sourcec = &g_dist_sourcec[8]; 
  const double *F_0_sourcer = &g_dist_sourcer[0]; 
  const double *G_1_sourcer = &g_dist_sourcer[8]; 

  const double *p_force_F_0 = &p_force[0]; 
  const double *p_force_G_1 = &p_force[3]; 

  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[8]; 

  double alphaSurf_F_0_l[3] = {0.0}; 
  alphaSurf_F_0_l[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force_F_0[0]; 
  alphaSurf_F_0_l[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force_F_0[1]; 
  alphaSurf_F_0_l[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force_F_0[2]; 

  double alphaSurf_F_0_r[3] = {0.0}; 
  alphaSurf_F_0_r[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force_F_0[0]; 
  alphaSurf_F_0_r[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force_F_0[1]; 
  alphaSurf_F_0_r[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force_F_0[2]; 

  double alphaSurf_G_1_l[3] = {0.0}; 
  alphaSurf_G_1_l[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force_G_1[0]; 
  alphaSurf_G_1_l[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force_G_1[1]; 
  alphaSurf_G_1_l[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force_G_1[2]; 

  double alphaSurf_G_1_r[3] = {0.0}; 
  alphaSurf_G_1_r[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force_G_1[0]; 
  alphaSurf_G_1_r[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force_G_1[1]; 
  alphaSurf_G_1_r[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force_G_1[2]; 

  double div_b_Surf[3] = {0.0}; 
  div_b_Surf[0] = div_b[0]; 
  div_b_Surf[1] = div_b[1]; 
  div_b_Surf[2] = div_b[2]; 

  double F_0_UpwindQuad_l[3] = {0.0};
  double F_0_UpwindQuad_r[3] = {0.0};
  double F_0_Upwind_l[3] = {0.0};
  double F_0_Upwind_r[3] = {0.0};
  double Ghat_F_0_l[3] = {0.0}; 
  double Ghat_F_0_r[3] = {0.0}; 
  double G_1_UpwindQuad_l[3] = {0.0};
  double G_1_UpwindQuad_r[3] = {0.0};
  double G_1_Upwind_l[3] = {0.0};
  double G_1_Upwind_r[3] = {0.0};
  double Ghat_G_1_l[3] = {0.0}; 
  double Ghat_G_1_r[3] = {0.0}; 

  if (0.6324555320336759*alphaSurf_F_0_l[2]-0.9486832980505137*alphaSurf_F_0_l[1]+0.7071067811865475*alphaSurf_F_0_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(F_0c); 
  } 
  if (0.6324555320336759*alphaSurf_F_0_r[2]-0.9486832980505137*alphaSurf_F_0_r[1]+0.7071067811865475*alphaSurf_F_0_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(F_0r); 
  } 
  if (0.7071067811865475*alphaSurf_F_0_l[0]-0.7905694150420947*alphaSurf_F_0_l[2] > 0) { 
    F_0_UpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(F_0c); 
  } 
  if (0.7071067811865475*alphaSurf_F_0_r[0]-0.7905694150420947*alphaSurf_F_0_r[2] > 0) { 
    F_0_UpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(F_0r); 
  } 
  if (0.6324555320336759*alphaSurf_F_0_l[2]+0.9486832980505137*alphaSurf_F_0_l[1]+0.7071067811865475*alphaSurf_F_0_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(F_0c); 
  } 
  if (0.6324555320336759*alphaSurf_F_0_r[2]+0.9486832980505137*alphaSurf_F_0_r[1]+0.7071067811865475*alphaSurf_F_0_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(F_0r); 
  } 
  if (0.6324555320336759*alphaSurf_G_1_l[2]-0.9486832980505137*alphaSurf_G_1_l[1]+0.7071067811865475*alphaSurf_G_1_l[0] > 0) { 
    G_1_UpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(G_1c); 
  } 
  if (0.6324555320336759*alphaSurf_G_1_r[2]-0.9486832980505137*alphaSurf_G_1_r[1]+0.7071067811865475*alphaSurf_G_1_r[0] > 0) { 
    G_1_UpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(G_1r); 
  } 
  if (0.7071067811865475*alphaSurf_G_1_l[0]-0.7905694150420947*alphaSurf_G_1_l[2] > 0) { 
    G_1_UpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(G_1c); 
  } 
  if (0.7071067811865475*alphaSurf_G_1_r[0]-0.7905694150420947*alphaSurf_G_1_r[2] > 0) { 
    G_1_UpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(G_1r); 
  } 
  if (0.6324555320336759*alphaSurf_G_1_l[2]+0.9486832980505137*alphaSurf_G_1_l[1]+0.7071067811865475*alphaSurf_G_1_l[0] > 0) { 
    G_1_UpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(G_1c); 
  } 
  if (0.6324555320336759*alphaSurf_G_1_r[2]+0.9486832980505137*alphaSurf_G_1_r[1]+0.7071067811865475*alphaSurf_G_1_r[0] > 0) { 
    G_1_UpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  ser_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  ser_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  ser_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  double F_0_div_b_UpwindQuad_l[3] = {0.0};
  double F_0_div_b_UpwindQuad_r[3] = {0.0};
  double F_0_div_b_Upwind_l[3] = {0.0};
  double F_0_div_b_Upwind_r[3] = {0.0};
  double Ghat_F_0_div_b_l[3] = {0.0}; 
  double Ghat_F_0_div_b_r[3] = {0.0}; 

  double G_1_div_b_UpwindQuad_l[3] = {0.0};
  double G_1_div_b_UpwindQuad_r[3] = {0.0};
  double G_1_div_b_Upwind_l[3] = {0.0};
  double G_1_div_b_Upwind_r[3] = {0.0};
  double Ghat_G_1_div_b_l[3] = {0.0}; 
  double Ghat_G_1_div_b_r[3] = {0.0}; 

  if (0.6324555320336759*div_b_Surf[2]-0.9486832980505137*div_b_Surf[1]+0.7071067811865475*div_b_Surf[0] < 0) { 
    F_0_div_b_UpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[0] = ser_2x_p2_surfx2_eval_quad_node_0_l(G_1_sourcer); 
  } 
  if (0.7071067811865475*div_b_Surf[0]-0.7905694150420947*div_b_Surf[2] < 0) { 
    F_0_div_b_UpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[1] = ser_2x_p2_surfx2_eval_quad_node_1_l(G_1_sourcer); 
  } 
  if (0.6324555320336759*div_b_Surf[2]+0.9486832980505137*div_b_Surf[1]+0.7071067811865475*div_b_Surf[0] < 0) { 
    F_0_div_b_UpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[2] = ser_2x_p2_surfx2_eval_quad_node_2_l(G_1_sourcer); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(F_0_div_b_UpwindQuad_l, F_0_div_b_Upwind_l); 
  ser_2x_p2_upwind_quad_to_modal(F_0_div_b_UpwindQuad_r, F_0_div_b_Upwind_r); 
  ser_2x_p2_upwind_quad_to_modal(G_1_div_b_UpwindQuad_l, G_1_div_b_Upwind_l); 
  ser_2x_p2_upwind_quad_to_modal(G_1_div_b_UpwindQuad_r, G_1_div_b_Upwind_r); 

  Ghat_F_0_l[0] = 0.7071067811865475*F_0_Upwind_l[2]*alphaSurf_F_0_l[2]+0.7071067811865475*F_0_Upwind_l[1]*alphaSurf_F_0_l[1]+0.7071067811865475*F_0_Upwind_l[0]*alphaSurf_F_0_l[0]; 
  Ghat_F_0_l[1] = 0.6324555320336759*F_0_Upwind_l[1]*alphaSurf_F_0_l[2]+0.6324555320336759*alphaSurf_F_0_l[1]*F_0_Upwind_l[2]+0.7071067811865475*F_0_Upwind_l[0]*alphaSurf_F_0_l[1]+0.7071067811865475*alphaSurf_F_0_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.4517539514526256*F_0_Upwind_l[2]*alphaSurf_F_0_l[2]+0.7071067811865475*F_0_Upwind_l[0]*alphaSurf_F_0_l[2]+0.7071067811865475*alphaSurf_F_0_l[0]*F_0_Upwind_l[2]+0.6324555320336759*F_0_Upwind_l[1]*alphaSurf_F_0_l[1]; 
  Ghat_G_1_l[0] = 0.7071067811865475*G_1_Upwind_l[2]*alphaSurf_G_1_l[2]+0.7071067811865475*G_1_Upwind_l[1]*alphaSurf_G_1_l[1]+0.7071067811865475*G_1_Upwind_l[0]*alphaSurf_G_1_l[0]; 
  Ghat_G_1_l[1] = 0.6324555320336759*G_1_Upwind_l[1]*alphaSurf_G_1_l[2]+0.6324555320336759*alphaSurf_G_1_l[1]*G_1_Upwind_l[2]+0.7071067811865475*G_1_Upwind_l[0]*alphaSurf_G_1_l[1]+0.7071067811865475*alphaSurf_G_1_l[0]*G_1_Upwind_l[1]; 
  Ghat_G_1_l[2] = 0.4517539514526256*G_1_Upwind_l[2]*alphaSurf_G_1_l[2]+0.7071067811865475*G_1_Upwind_l[0]*alphaSurf_G_1_l[2]+0.7071067811865475*alphaSurf_G_1_l[0]*G_1_Upwind_l[2]+0.6324555320336759*G_1_Upwind_l[1]*alphaSurf_G_1_l[1]; 
  Ghat_F_0_div_b_l[0] = 0.7071067811865475*F_0_div_b_Upwind_l[2]*div_b_Surf[2]+0.7071067811865475*F_0_div_b_Upwind_l[1]*div_b_Surf[1]+0.7071067811865475*F_0_div_b_Upwind_l[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b_l[1] = 0.6324555320336759*F_0_div_b_Upwind_l[1]*div_b_Surf[2]+0.6324555320336759*div_b_Surf[1]*F_0_div_b_Upwind_l[2]+0.7071067811865475*F_0_div_b_Upwind_l[0]*div_b_Surf[1]+0.7071067811865475*div_b_Surf[0]*F_0_div_b_Upwind_l[1]; 
  Ghat_F_0_div_b_l[2] = 0.4517539514526256*F_0_div_b_Upwind_l[2]*div_b_Surf[2]+0.7071067811865475*F_0_div_b_Upwind_l[0]*div_b_Surf[2]+0.7071067811865475*div_b_Surf[0]*F_0_div_b_Upwind_l[2]+0.6324555320336759*F_0_div_b_Upwind_l[1]*div_b_Surf[1]; 
  Ghat_G_1_div_b_l[0] = 0.7071067811865475*G_1_div_b_Upwind_l[2]*div_b_Surf[2]+0.7071067811865475*G_1_div_b_Upwind_l[1]*div_b_Surf[1]+0.7071067811865475*G_1_div_b_Upwind_l[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b_l[1] = 0.6324555320336759*G_1_div_b_Upwind_l[1]*div_b_Surf[2]+0.6324555320336759*div_b_Surf[1]*G_1_div_b_Upwind_l[2]+0.7071067811865475*G_1_div_b_Upwind_l[0]*div_b_Surf[1]+0.7071067811865475*div_b_Surf[0]*G_1_div_b_Upwind_l[1]; 
  Ghat_G_1_div_b_l[2] = 0.4517539514526256*G_1_div_b_Upwind_l[2]*div_b_Surf[2]+0.7071067811865475*G_1_div_b_Upwind_l[0]*div_b_Surf[2]+0.7071067811865475*div_b_Surf[0]*G_1_div_b_Upwind_l[2]+0.6324555320336759*G_1_div_b_Upwind_l[1]*div_b_Surf[1]; 

  Ghat_F_0_r[0] = 0.7071067811865475*F_0_Upwind_r[2]*alphaSurf_F_0_r[2]+0.7071067811865475*F_0_Upwind_r[1]*alphaSurf_F_0_r[1]+0.7071067811865475*F_0_Upwind_r[0]*alphaSurf_F_0_r[0]; 
  Ghat_F_0_r[1] = 0.6324555320336759*F_0_Upwind_r[1]*alphaSurf_F_0_r[2]+0.6324555320336759*alphaSurf_F_0_r[1]*F_0_Upwind_r[2]+0.7071067811865475*F_0_Upwind_r[0]*alphaSurf_F_0_r[1]+0.7071067811865475*alphaSurf_F_0_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.4517539514526256*F_0_Upwind_r[2]*alphaSurf_F_0_r[2]+0.7071067811865475*F_0_Upwind_r[0]*alphaSurf_F_0_r[2]+0.7071067811865475*alphaSurf_F_0_r[0]*F_0_Upwind_r[2]+0.6324555320336759*F_0_Upwind_r[1]*alphaSurf_F_0_r[1]; 
  Ghat_G_1_r[0] = 0.7071067811865475*G_1_Upwind_r[2]*alphaSurf_G_1_r[2]+0.7071067811865475*G_1_Upwind_r[1]*alphaSurf_G_1_r[1]+0.7071067811865475*G_1_Upwind_r[0]*alphaSurf_G_1_r[0]; 
  Ghat_G_1_r[1] = 0.6324555320336759*G_1_Upwind_r[1]*alphaSurf_G_1_r[2]+0.6324555320336759*alphaSurf_G_1_r[1]*G_1_Upwind_r[2]+0.7071067811865475*G_1_Upwind_r[0]*alphaSurf_G_1_r[1]+0.7071067811865475*alphaSurf_G_1_r[0]*G_1_Upwind_r[1]; 
  Ghat_G_1_r[2] = 0.4517539514526256*G_1_Upwind_r[2]*alphaSurf_G_1_r[2]+0.7071067811865475*G_1_Upwind_r[0]*alphaSurf_G_1_r[2]+0.7071067811865475*alphaSurf_G_1_r[0]*G_1_Upwind_r[2]+0.6324555320336759*G_1_Upwind_r[1]*alphaSurf_G_1_r[1]; 
  Ghat_F_0_div_b_r[0] = 0.7071067811865475*F_0_div_b_Upwind_r[2]*div_b_Surf[2]+0.7071067811865475*F_0_div_b_Upwind_r[1]*div_b_Surf[1]+0.7071067811865475*F_0_div_b_Upwind_r[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b_r[1] = 0.6324555320336759*F_0_div_b_Upwind_r[1]*div_b_Surf[2]+0.6324555320336759*div_b_Surf[1]*F_0_div_b_Upwind_r[2]+0.7071067811865475*F_0_div_b_Upwind_r[0]*div_b_Surf[1]+0.7071067811865475*div_b_Surf[0]*F_0_div_b_Upwind_r[1]; 
  Ghat_F_0_div_b_r[2] = 0.4517539514526256*F_0_div_b_Upwind_r[2]*div_b_Surf[2]+0.7071067811865475*F_0_div_b_Upwind_r[0]*div_b_Surf[2]+0.7071067811865475*div_b_Surf[0]*F_0_div_b_Upwind_r[2]+0.6324555320336759*F_0_div_b_Upwind_r[1]*div_b_Surf[1]; 
  Ghat_G_1_div_b_r[0] = 0.7071067811865475*G_1_div_b_Upwind_r[2]*div_b_Surf[2]+0.7071067811865475*G_1_div_b_Upwind_r[1]*div_b_Surf[1]+0.7071067811865475*G_1_div_b_Upwind_r[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b_r[1] = 0.6324555320336759*G_1_div_b_Upwind_r[1]*div_b_Surf[2]+0.6324555320336759*div_b_Surf[1]*G_1_div_b_Upwind_r[2]+0.7071067811865475*G_1_div_b_Upwind_r[0]*div_b_Surf[1]+0.7071067811865475*div_b_Surf[0]*G_1_div_b_Upwind_r[1]; 
  Ghat_G_1_div_b_r[2] = 0.4517539514526256*G_1_div_b_Upwind_r[2]*div_b_Surf[2]+0.7071067811865475*G_1_div_b_Upwind_r[0]*div_b_Surf[2]+0.7071067811865475*div_b_Surf[0]*G_1_div_b_Upwind_r[2]+0.6324555320336759*G_1_div_b_Upwind_r[1]*div_b_Surf[1]; 

  out_F_0[0] += ((-0.7071067811865475*Ghat_F_0_r[0])+0.7071067811865475*(Ghat_F_0_l[0]+Ghat_F_0_div_b_r[0])-0.7071067811865475*Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[1] += ((-0.7071067811865475*Ghat_F_0_r[1])+0.7071067811865475*(Ghat_F_0_l[1]+Ghat_F_0_div_b_r[1])-0.7071067811865475*Ghat_F_0_div_b_l[1])*dv1par; 
  out_F_0[2] += (1.224744871391589*(Ghat_F_0_div_b_r[0]+Ghat_F_0_div_b_l[0])-1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0]))*dv1par; 
  out_F_0[3] += (1.224744871391589*(Ghat_F_0_div_b_r[1]+Ghat_F_0_div_b_l[1])-1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1]))*dv1par; 
  out_F_0[4] += ((-0.7071067811865475*Ghat_F_0_r[2])+0.7071067811865475*(Ghat_F_0_l[2]+Ghat_F_0_div_b_r[2])-0.7071067811865475*Ghat_F_0_div_b_l[2])*dv1par; 
  out_F_0[5] += ((-1.58113883008419*Ghat_F_0_r[0])+1.58113883008419*(Ghat_F_0_l[0]+Ghat_F_0_div_b_r[0])-1.58113883008419*Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[6] += (1.224744871391589*(Ghat_F_0_div_b_r[2]+Ghat_F_0_div_b_l[2])-1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2]))*dv1par; 
  out_F_0[7] += ((-1.58113883008419*Ghat_F_0_r[1])+1.58113883008419*(Ghat_F_0_l[1]+Ghat_F_0_div_b_r[1])-1.58113883008419*Ghat_F_0_div_b_l[1])*dv1par; 
  out_G_1[0] += ((-0.7071067811865475*Ghat_G_1_r[0])+0.7071067811865475*(Ghat_G_1_l[0]+Ghat_G_1_div_b_r[0])-0.7071067811865475*Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[1] += ((-0.7071067811865475*Ghat_G_1_r[1])+0.7071067811865475*(Ghat_G_1_l[1]+Ghat_G_1_div_b_r[1])-0.7071067811865475*Ghat_G_1_div_b_l[1])*dv1par; 
  out_G_1[2] += (1.224744871391589*(Ghat_G_1_div_b_r[0]+Ghat_G_1_div_b_l[0])-1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0]))*dv1par; 
  out_G_1[3] += (1.224744871391589*(Ghat_G_1_div_b_r[1]+Ghat_G_1_div_b_l[1])-1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1]))*dv1par; 
  out_G_1[4] += ((-0.7071067811865475*Ghat_G_1_r[2])+0.7071067811865475*(Ghat_G_1_l[2]+Ghat_G_1_div_b_r[2])-0.7071067811865475*Ghat_G_1_div_b_l[2])*dv1par; 
  out_G_1[5] += ((-1.58113883008419*Ghat_G_1_r[0])+1.58113883008419*(Ghat_G_1_l[0]+Ghat_G_1_div_b_r[0])-1.58113883008419*Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[6] += (1.224744871391589*(Ghat_G_1_div_b_r[2]+Ghat_G_1_div_b_l[2])-1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2]))*dv1par; 
  out_G_1[7] += ((-1.58113883008419*Ghat_G_1_r[1])+1.58113883008419*(Ghat_G_1_l[1]+Ghat_G_1_div_b_r[1])-1.58113883008419*Ghat_G_1_div_b_l[1])*dv1par; 

} 
