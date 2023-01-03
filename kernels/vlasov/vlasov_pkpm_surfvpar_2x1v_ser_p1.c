#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_2x1v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfvpar_2x1v_ser_p1(const double *w, const double *dxv, 
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
  const double dv1par = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[12]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[12]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[12]; 

  const double *F_0_sourcel = &g_dist_sourcel[0]; 
  const double *G_1_sourcel = &g_dist_sourcel[12]; 
  const double *F_0_sourcec = &g_dist_sourcec[0]; 
  const double *G_1_sourcec = &g_dist_sourcec[12]; 
  const double *F_0_sourcer = &g_dist_sourcer[0]; 
  const double *G_1_sourcer = &g_dist_sourcer[12]; 

  const double *p_force_F_0 = &p_force[0]; 
  const double *p_force_G_1 = &p_force[4]; 

  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 

  double alphaSurf_F_0_l[4] = {0.0}; 
  alphaSurf_F_0_l[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force_F_0[0]; 
  alphaSurf_F_0_l[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force_F_0[1]; 
  alphaSurf_F_0_l[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force_F_0[2]; 
  alphaSurf_F_0_l[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force_F_0[3]; 

  double alphaSurf_F_0_r[4] = {0.0}; 
  alphaSurf_F_0_r[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force_F_0[0]; 
  alphaSurf_F_0_r[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force_F_0[1]; 
  alphaSurf_F_0_r[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force_F_0[2]; 
  alphaSurf_F_0_r[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force_F_0[3]; 

  double alphaSurf_G_1_l[4] = {0.0}; 
  alphaSurf_G_1_l[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force_G_1[0]; 
  alphaSurf_G_1_l[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force_G_1[1]; 
  alphaSurf_G_1_l[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force_G_1[2]; 
  alphaSurf_G_1_l[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force_G_1[3]; 

  double alphaSurf_G_1_r[4] = {0.0}; 
  alphaSurf_G_1_r[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force_G_1[0]; 
  alphaSurf_G_1_r[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force_G_1[1]; 
  alphaSurf_G_1_r[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force_G_1[2]; 
  alphaSurf_G_1_r[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force_G_1[3]; 

  double div_b_Surf[4] = {0.0}; 
  div_b_Surf[0] = div_b[0]; 
  div_b_Surf[1] = div_b[1]; 
  div_b_Surf[2] = div_b[2]; 
  div_b_Surf[3] = div_b[3]; 

  double F_0_UpwindQuad_l[4] = {0.0};
  double F_0_UpwindQuad_r[4] = {0.0};
  double F_0_Upwind_l[4] = {0.0};
  double F_0_Upwind_r[4] = {0.0};
  double Ghat_F_0_l[4] = {0.0}; 
  double Ghat_F_0_r[4] = {0.0}; 
  double G_1_UpwindQuad_l[4] = {0.0};
  double G_1_UpwindQuad_r[4] = {0.0};
  double G_1_Upwind_l[4] = {0.0};
  double G_1_Upwind_r[4] = {0.0};
  double Ghat_G_1_l[4] = {0.0}; 
  double Ghat_G_1_r[4] = {0.0}; 

  if (0.5*alphaSurf_F_0_l[3]-0.5*(alphaSurf_F_0_l[2]+alphaSurf_F_0_l[1])+0.5*alphaSurf_F_0_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0c); 
  } 
  if (0.5*alphaSurf_F_0_r[3]-0.5*(alphaSurf_F_0_r[2]+alphaSurf_F_0_r[1])+0.5*alphaSurf_F_0_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0r); 
  } 
  if ((-0.5*alphaSurf_F_0_l[3])+0.5*alphaSurf_F_0_l[2]-0.5*alphaSurf_F_0_l[1]+0.5*alphaSurf_F_0_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0c); 
  } 
  if ((-0.5*alphaSurf_F_0_r[3])+0.5*alphaSurf_F_0_r[2]-0.5*alphaSurf_F_0_r[1]+0.5*alphaSurf_F_0_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0r); 
  } 
  if (0.5*(alphaSurf_F_0_l[1]+alphaSurf_F_0_l[0])-0.5*(alphaSurf_F_0_l[3]+alphaSurf_F_0_l[2]) > 0) { 
    F_0_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0c); 
  } 
  if (0.5*(alphaSurf_F_0_r[1]+alphaSurf_F_0_r[0])-0.5*(alphaSurf_F_0_r[3]+alphaSurf_F_0_r[2]) > 0) { 
    F_0_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0r); 
  } 
  if (0.5*(alphaSurf_F_0_l[3]+alphaSurf_F_0_l[2]+alphaSurf_F_0_l[1]+alphaSurf_F_0_l[0]) > 0) { 
    F_0_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0c); 
  } 
  if (0.5*(alphaSurf_F_0_r[3]+alphaSurf_F_0_r[2]+alphaSurf_F_0_r[1]+alphaSurf_F_0_r[0]) > 0) { 
    F_0_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0r); 
  } 
  if (0.5*alphaSurf_G_1_l[3]-0.5*(alphaSurf_G_1_l[2]+alphaSurf_G_1_l[1])+0.5*alphaSurf_G_1_l[0] > 0) { 
    G_1_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1c); 
  } 
  if (0.5*alphaSurf_G_1_r[3]-0.5*(alphaSurf_G_1_r[2]+alphaSurf_G_1_r[1])+0.5*alphaSurf_G_1_r[0] > 0) { 
    G_1_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1r); 
  } 
  if ((-0.5*alphaSurf_G_1_l[3])+0.5*alphaSurf_G_1_l[2]-0.5*alphaSurf_G_1_l[1]+0.5*alphaSurf_G_1_l[0] > 0) { 
    G_1_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1c); 
  } 
  if ((-0.5*alphaSurf_G_1_r[3])+0.5*alphaSurf_G_1_r[2]-0.5*alphaSurf_G_1_r[1]+0.5*alphaSurf_G_1_r[0] > 0) { 
    G_1_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1r); 
  } 
  if (0.5*(alphaSurf_G_1_l[1]+alphaSurf_G_1_l[0])-0.5*(alphaSurf_G_1_l[3]+alphaSurf_G_1_l[2]) > 0) { 
    G_1_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1c); 
  } 
  if (0.5*(alphaSurf_G_1_r[1]+alphaSurf_G_1_r[0])-0.5*(alphaSurf_G_1_r[3]+alphaSurf_G_1_r[2]) > 0) { 
    G_1_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1r); 
  } 
  if (0.5*(alphaSurf_G_1_l[3]+alphaSurf_G_1_l[2]+alphaSurf_G_1_l[1]+alphaSurf_G_1_l[0]) > 0) { 
    G_1_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1c); 
  } 
  if (0.5*(alphaSurf_G_1_r[3]+alphaSurf_G_1_r[2]+alphaSurf_G_1_r[1]+alphaSurf_G_1_r[0]) > 0) { 
    G_1_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  double F_0_div_b_UpwindQuad_l[4] = {0.0};
  double F_0_div_b_UpwindQuad_r[4] = {0.0};
  double F_0_div_b_Upwind_l[4] = {0.0};
  double F_0_div_b_Upwind_r[4] = {0.0};
  double Ghat_F_0_div_b_l[4] = {0.0}; 
  double Ghat_F_0_div_b_r[4] = {0.0}; 

  double G_1_div_b_UpwindQuad_l[4] = {0.0};
  double G_1_div_b_UpwindQuad_r[4] = {0.0};
  double G_1_div_b_Upwind_l[4] = {0.0};
  double G_1_div_b_Upwind_r[4] = {0.0};
  double Ghat_G_1_div_b_l[4] = {0.0}; 
  double Ghat_G_1_div_b_r[4] = {0.0}; 

  if (0.5*div_b_Surf[3]-0.5*(div_b_Surf[2]+div_b_Surf[1])+0.5*div_b_Surf[0] < 0) { 
    F_0_div_b_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[0] = hyb_2x1v_p1_surfx3_eval_quad_node_0_l(G_1_sourcer); 
  } 
  if ((-0.5*div_b_Surf[3])+0.5*div_b_Surf[2]-0.5*div_b_Surf[1]+0.5*div_b_Surf[0] < 0) { 
    F_0_div_b_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[1] = hyb_2x1v_p1_surfx3_eval_quad_node_1_l(G_1_sourcer); 
  } 
  if (0.5*(div_b_Surf[1]+div_b_Surf[0])-0.5*(div_b_Surf[3]+div_b_Surf[2]) < 0) { 
    F_0_div_b_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[2] = hyb_2x1v_p1_surfx3_eval_quad_node_2_l(G_1_sourcer); 
  } 
  if (0.5*(div_b_Surf[3]+div_b_Surf[2]+div_b_Surf[1]+div_b_Surf[0]) < 0) { 
    F_0_div_b_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[3] = hyb_2x1v_p1_surfx3_eval_quad_node_3_l(G_1_sourcer); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_div_b_UpwindQuad_l, F_0_div_b_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(F_0_div_b_UpwindQuad_r, F_0_div_b_Upwind_r); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_div_b_UpwindQuad_l, G_1_div_b_Upwind_l); 
  hyb_2x1v_p1_vdir_upwind_quad_to_modal(G_1_div_b_UpwindQuad_r, G_1_div_b_Upwind_r); 

  Ghat_F_0_l[0] = 0.5*F_0_Upwind_l[3]*alphaSurf_F_0_l[3]+0.5*F_0_Upwind_l[2]*alphaSurf_F_0_l[2]+0.5*F_0_Upwind_l[1]*alphaSurf_F_0_l[1]+0.5*F_0_Upwind_l[0]*alphaSurf_F_0_l[0]; 
  Ghat_F_0_l[1] = 0.5*F_0_Upwind_l[2]*alphaSurf_F_0_l[3]+0.5*alphaSurf_F_0_l[2]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_F_0_l[1]+0.5*alphaSurf_F_0_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.5*F_0_Upwind_l[1]*alphaSurf_F_0_l[3]+0.5*alphaSurf_F_0_l[1]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_F_0_l[2]+0.5*alphaSurf_F_0_l[0]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[3] = 0.5*F_0_Upwind_l[0]*alphaSurf_F_0_l[3]+0.5*alphaSurf_F_0_l[0]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[1]*alphaSurf_F_0_l[2]+0.5*alphaSurf_F_0_l[1]*F_0_Upwind_l[2]; 
  Ghat_G_1_l[0] = 0.5*G_1_Upwind_l[3]*alphaSurf_G_1_l[3]+0.5*G_1_Upwind_l[2]*alphaSurf_G_1_l[2]+0.5*G_1_Upwind_l[1]*alphaSurf_G_1_l[1]+0.5*G_1_Upwind_l[0]*alphaSurf_G_1_l[0]; 
  Ghat_G_1_l[1] = 0.5*G_1_Upwind_l[2]*alphaSurf_G_1_l[3]+0.5*alphaSurf_G_1_l[2]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_G_1_l[1]+0.5*alphaSurf_G_1_l[0]*G_1_Upwind_l[1]; 
  Ghat_G_1_l[2] = 0.5*G_1_Upwind_l[1]*alphaSurf_G_1_l[3]+0.5*alphaSurf_G_1_l[1]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_G_1_l[2]+0.5*alphaSurf_G_1_l[0]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[3] = 0.5*G_1_Upwind_l[0]*alphaSurf_G_1_l[3]+0.5*alphaSurf_G_1_l[0]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[1]*alphaSurf_G_1_l[2]+0.5*alphaSurf_G_1_l[1]*G_1_Upwind_l[2]; 
  Ghat_F_0_div_b_l[0] = 0.5*F_0_div_b_Upwind_l[3]*div_b_Surf[3]+0.5*F_0_div_b_Upwind_l[2]*div_b_Surf[2]+0.5*F_0_div_b_Upwind_l[1]*div_b_Surf[1]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b_l[1] = 0.5*F_0_div_b_Upwind_l[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*F_0_div_b_Upwind_l[3]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[1]; 
  Ghat_F_0_div_b_l[2] = 0.5*F_0_div_b_Upwind_l[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*F_0_div_b_Upwind_l[3]+0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[2]; 
  Ghat_F_0_div_b_l[3] = 0.5*F_0_div_b_Upwind_l[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_l[3]+0.5*F_0_div_b_Upwind_l[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*F_0_div_b_Upwind_l[2]; 
  Ghat_G_1_div_b_l[0] = 0.5*G_1_div_b_Upwind_l[3]*div_b_Surf[3]+0.5*G_1_div_b_Upwind_l[2]*div_b_Surf[2]+0.5*G_1_div_b_Upwind_l[1]*div_b_Surf[1]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b_l[1] = 0.5*G_1_div_b_Upwind_l[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*G_1_div_b_Upwind_l[3]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[1]; 
  Ghat_G_1_div_b_l[2] = 0.5*G_1_div_b_Upwind_l[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*G_1_div_b_Upwind_l[3]+0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[2]; 
  Ghat_G_1_div_b_l[3] = 0.5*G_1_div_b_Upwind_l[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_l[3]+0.5*G_1_div_b_Upwind_l[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*G_1_div_b_Upwind_l[2]; 

  Ghat_F_0_r[0] = 0.5*F_0_Upwind_r[3]*alphaSurf_F_0_r[3]+0.5*F_0_Upwind_r[2]*alphaSurf_F_0_r[2]+0.5*F_0_Upwind_r[1]*alphaSurf_F_0_r[1]+0.5*F_0_Upwind_r[0]*alphaSurf_F_0_r[0]; 
  Ghat_F_0_r[1] = 0.5*F_0_Upwind_r[2]*alphaSurf_F_0_r[3]+0.5*alphaSurf_F_0_r[2]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_F_0_r[1]+0.5*alphaSurf_F_0_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.5*F_0_Upwind_r[1]*alphaSurf_F_0_r[3]+0.5*alphaSurf_F_0_r[1]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_F_0_r[2]+0.5*alphaSurf_F_0_r[0]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[3] = 0.5*F_0_Upwind_r[0]*alphaSurf_F_0_r[3]+0.5*alphaSurf_F_0_r[0]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[1]*alphaSurf_F_0_r[2]+0.5*alphaSurf_F_0_r[1]*F_0_Upwind_r[2]; 
  Ghat_G_1_r[0] = 0.5*G_1_Upwind_r[3]*alphaSurf_G_1_r[3]+0.5*G_1_Upwind_r[2]*alphaSurf_G_1_r[2]+0.5*G_1_Upwind_r[1]*alphaSurf_G_1_r[1]+0.5*G_1_Upwind_r[0]*alphaSurf_G_1_r[0]; 
  Ghat_G_1_r[1] = 0.5*G_1_Upwind_r[2]*alphaSurf_G_1_r[3]+0.5*alphaSurf_G_1_r[2]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_G_1_r[1]+0.5*alphaSurf_G_1_r[0]*G_1_Upwind_r[1]; 
  Ghat_G_1_r[2] = 0.5*G_1_Upwind_r[1]*alphaSurf_G_1_r[3]+0.5*alphaSurf_G_1_r[1]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_G_1_r[2]+0.5*alphaSurf_G_1_r[0]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[3] = 0.5*G_1_Upwind_r[0]*alphaSurf_G_1_r[3]+0.5*alphaSurf_G_1_r[0]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[1]*alphaSurf_G_1_r[2]+0.5*alphaSurf_G_1_r[1]*G_1_Upwind_r[2]; 
  Ghat_F_0_div_b_r[0] = 0.5*F_0_div_b_Upwind_r[3]*div_b_Surf[3]+0.5*F_0_div_b_Upwind_r[2]*div_b_Surf[2]+0.5*F_0_div_b_Upwind_r[1]*div_b_Surf[1]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b_r[1] = 0.5*F_0_div_b_Upwind_r[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*F_0_div_b_Upwind_r[3]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[1]; 
  Ghat_F_0_div_b_r[2] = 0.5*F_0_div_b_Upwind_r[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*F_0_div_b_Upwind_r[3]+0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[2]; 
  Ghat_F_0_div_b_r[3] = 0.5*F_0_div_b_Upwind_r[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*F_0_div_b_Upwind_r[3]+0.5*F_0_div_b_Upwind_r[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*F_0_div_b_Upwind_r[2]; 
  Ghat_G_1_div_b_r[0] = 0.5*G_1_div_b_Upwind_r[3]*div_b_Surf[3]+0.5*G_1_div_b_Upwind_r[2]*div_b_Surf[2]+0.5*G_1_div_b_Upwind_r[1]*div_b_Surf[1]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b_r[1] = 0.5*G_1_div_b_Upwind_r[2]*div_b_Surf[3]+0.5*div_b_Surf[2]*G_1_div_b_Upwind_r[3]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[1]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[1]; 
  Ghat_G_1_div_b_r[2] = 0.5*G_1_div_b_Upwind_r[1]*div_b_Surf[3]+0.5*div_b_Surf[1]*G_1_div_b_Upwind_r[3]+0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[2]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[2]; 
  Ghat_G_1_div_b_r[3] = 0.5*G_1_div_b_Upwind_r[0]*div_b_Surf[3]+0.5*div_b_Surf[0]*G_1_div_b_Upwind_r[3]+0.5*G_1_div_b_Upwind_r[1]*div_b_Surf[2]+0.5*div_b_Surf[1]*G_1_div_b_Upwind_r[2]; 

  out_F_0[0] += ((-0.7071067811865475*Ghat_F_0_r[0])+0.7071067811865475*(Ghat_F_0_l[0]+Ghat_F_0_div_b_r[0])-0.7071067811865475*Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[1] += ((-0.7071067811865475*Ghat_F_0_r[1])+0.7071067811865475*(Ghat_F_0_l[1]+Ghat_F_0_div_b_r[1])-0.7071067811865475*Ghat_F_0_div_b_l[1])*dv1par; 
  out_F_0[2] += ((-0.7071067811865475*Ghat_F_0_r[2])+0.7071067811865475*(Ghat_F_0_l[2]+Ghat_F_0_div_b_r[2])-0.7071067811865475*Ghat_F_0_div_b_l[2])*dv1par; 
  out_F_0[3] += (1.224744871391589*(Ghat_F_0_div_b_r[0]+Ghat_F_0_div_b_l[0])-1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0]))*dv1par; 
  out_F_0[4] += ((-0.7071067811865475*Ghat_F_0_r[3])+0.7071067811865475*(Ghat_F_0_l[3]+Ghat_F_0_div_b_r[3])-0.7071067811865475*Ghat_F_0_div_b_l[3])*dv1par; 
  out_F_0[5] += (1.224744871391589*(Ghat_F_0_div_b_r[1]+Ghat_F_0_div_b_l[1])-1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1]))*dv1par; 
  out_F_0[6] += (1.224744871391589*(Ghat_F_0_div_b_r[2]+Ghat_F_0_div_b_l[2])-1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2]))*dv1par; 
  out_F_0[7] += (1.224744871391589*(Ghat_F_0_div_b_r[3]+Ghat_F_0_div_b_l[3])-1.224744871391589*(Ghat_F_0_r[3]+Ghat_F_0_l[3]))*dv1par; 
  out_F_0[8] += ((-1.58113883008419*Ghat_F_0_r[0])+1.58113883008419*(Ghat_F_0_l[0]+Ghat_F_0_div_b_r[0])-1.58113883008419*Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[9] += ((-1.58113883008419*Ghat_F_0_r[1])+1.58113883008419*(Ghat_F_0_l[1]+Ghat_F_0_div_b_r[1])-1.58113883008419*Ghat_F_0_div_b_l[1])*dv1par; 
  out_F_0[10] += ((-1.58113883008419*Ghat_F_0_r[2])+1.58113883008419*(Ghat_F_0_l[2]+Ghat_F_0_div_b_r[2])-1.58113883008419*Ghat_F_0_div_b_l[2])*dv1par; 
  out_F_0[11] += ((-1.58113883008419*Ghat_F_0_r[3])+1.58113883008419*(Ghat_F_0_l[3]+Ghat_F_0_div_b_r[3])-1.58113883008419*Ghat_F_0_div_b_l[3])*dv1par; 
  out_G_1[0] += ((-0.7071067811865475*Ghat_G_1_r[0])+0.7071067811865475*(Ghat_G_1_l[0]+Ghat_G_1_div_b_r[0])-0.7071067811865475*Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[1] += ((-0.7071067811865475*Ghat_G_1_r[1])+0.7071067811865475*(Ghat_G_1_l[1]+Ghat_G_1_div_b_r[1])-0.7071067811865475*Ghat_G_1_div_b_l[1])*dv1par; 
  out_G_1[2] += ((-0.7071067811865475*Ghat_G_1_r[2])+0.7071067811865475*(Ghat_G_1_l[2]+Ghat_G_1_div_b_r[2])-0.7071067811865475*Ghat_G_1_div_b_l[2])*dv1par; 
  out_G_1[3] += (1.224744871391589*(Ghat_G_1_div_b_r[0]+Ghat_G_1_div_b_l[0])-1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0]))*dv1par; 
  out_G_1[4] += ((-0.7071067811865475*Ghat_G_1_r[3])+0.7071067811865475*(Ghat_G_1_l[3]+Ghat_G_1_div_b_r[3])-0.7071067811865475*Ghat_G_1_div_b_l[3])*dv1par; 
  out_G_1[5] += (1.224744871391589*(Ghat_G_1_div_b_r[1]+Ghat_G_1_div_b_l[1])-1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1]))*dv1par; 
  out_G_1[6] += (1.224744871391589*(Ghat_G_1_div_b_r[2]+Ghat_G_1_div_b_l[2])-1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2]))*dv1par; 
  out_G_1[7] += (1.224744871391589*(Ghat_G_1_div_b_r[3]+Ghat_G_1_div_b_l[3])-1.224744871391589*(Ghat_G_1_r[3]+Ghat_G_1_l[3]))*dv1par; 
  out_G_1[8] += ((-1.58113883008419*Ghat_G_1_r[0])+1.58113883008419*(Ghat_G_1_l[0]+Ghat_G_1_div_b_r[0])-1.58113883008419*Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[9] += ((-1.58113883008419*Ghat_G_1_r[1])+1.58113883008419*(Ghat_G_1_l[1]+Ghat_G_1_div_b_r[1])-1.58113883008419*Ghat_G_1_div_b_l[1])*dv1par; 
  out_G_1[10] += ((-1.58113883008419*Ghat_G_1_r[2])+1.58113883008419*(Ghat_G_1_l[2]+Ghat_G_1_div_b_r[2])-1.58113883008419*Ghat_G_1_div_b_l[2])*dv1par; 
  out_G_1[11] += ((-1.58113883008419*Ghat_G_1_r[3])+1.58113883008419*(Ghat_G_1_l[3]+Ghat_G_1_div_b_r[3])-1.58113883008419*Ghat_G_1_div_b_l[3])*dv1par; 

} 
