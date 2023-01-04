#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_3x1v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_3x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfvpar_3x1v_ser_p1(const double *w, const double *dxv, 
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
  // g_dist_sourcel/c/r: [T_perp/m G, 2 (T_perp/m)^2 (F_2 - F_1)] in left/center/right cells.
  // fl/fc/fr:           Input Distribution function [F_0, T_perp/m G = T_perp/m (F_0 - F_1)] in left/center/right cells.
  // out:                Incremented distribution function in center cell.
  const double dv1par = 2.0/dxv[3]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[24]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[24]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[24]; 

  const double *F_0_sourcel = &g_dist_sourcel[0]; 
  const double *G_1_sourcel = &g_dist_sourcel[24]; 
  const double *F_0_sourcec = &g_dist_sourcec[0]; 
  const double *G_1_sourcec = &g_dist_sourcec[24]; 
  const double *F_0_sourcer = &g_dist_sourcer[0]; 
  const double *G_1_sourcer = &g_dist_sourcer[24]; 

  const double *p_force_F_0 = &p_force[0]; 
  const double *p_force_G_1 = &p_force[8]; 

  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[24]; 

  double alphaSurf_F_0_l[8] = {0.0}; 
  alphaSurf_F_0_l[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force_F_0[0]; 
  alphaSurf_F_0_l[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force_F_0[1]; 
  alphaSurf_F_0_l[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force_F_0[2]; 
  alphaSurf_F_0_l[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force_F_0[3]; 
  alphaSurf_F_0_l[4] = (-1.0*bb_grad_u[4]*wvpar)+0.5*bb_grad_u[4]*dvpar+p_force_F_0[4]; 
  alphaSurf_F_0_l[5] = (-1.0*bb_grad_u[5]*wvpar)+0.5*bb_grad_u[5]*dvpar+p_force_F_0[5]; 
  alphaSurf_F_0_l[6] = (-1.0*bb_grad_u[6]*wvpar)+0.5*bb_grad_u[6]*dvpar+p_force_F_0[6]; 
  alphaSurf_F_0_l[7] = (-1.0*bb_grad_u[7]*wvpar)+0.5*bb_grad_u[7]*dvpar+p_force_F_0[7]; 

  double alphaSurf_F_0_r[8] = {0.0}; 
  alphaSurf_F_0_r[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force_F_0[0]; 
  alphaSurf_F_0_r[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force_F_0[1]; 
  alphaSurf_F_0_r[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force_F_0[2]; 
  alphaSurf_F_0_r[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force_F_0[3]; 
  alphaSurf_F_0_r[4] = (-1.0*bb_grad_u[4]*wvpar)-0.5*bb_grad_u[4]*dvpar+p_force_F_0[4]; 
  alphaSurf_F_0_r[5] = (-1.0*bb_grad_u[5]*wvpar)-0.5*bb_grad_u[5]*dvpar+p_force_F_0[5]; 
  alphaSurf_F_0_r[6] = (-1.0*bb_grad_u[6]*wvpar)-0.5*bb_grad_u[6]*dvpar+p_force_F_0[6]; 
  alphaSurf_F_0_r[7] = (-1.0*bb_grad_u[7]*wvpar)-0.5*bb_grad_u[7]*dvpar+p_force_F_0[7]; 

  double alphaSurf_G_1_l[8] = {0.0}; 
  alphaSurf_G_1_l[0] = (-1.0*bb_grad_u[0]*wvpar)+0.5*bb_grad_u[0]*dvpar+p_force_G_1[0]; 
  alphaSurf_G_1_l[1] = (-1.0*bb_grad_u[1]*wvpar)+0.5*bb_grad_u[1]*dvpar+p_force_G_1[1]; 
  alphaSurf_G_1_l[2] = (-1.0*bb_grad_u[2]*wvpar)+0.5*bb_grad_u[2]*dvpar+p_force_G_1[2]; 
  alphaSurf_G_1_l[3] = (-1.0*bb_grad_u[3]*wvpar)+0.5*bb_grad_u[3]*dvpar+p_force_G_1[3]; 
  alphaSurf_G_1_l[4] = (-1.0*bb_grad_u[4]*wvpar)+0.5*bb_grad_u[4]*dvpar+p_force_G_1[4]; 
  alphaSurf_G_1_l[5] = (-1.0*bb_grad_u[5]*wvpar)+0.5*bb_grad_u[5]*dvpar+p_force_G_1[5]; 
  alphaSurf_G_1_l[6] = (-1.0*bb_grad_u[6]*wvpar)+0.5*bb_grad_u[6]*dvpar+p_force_G_1[6]; 
  alphaSurf_G_1_l[7] = (-1.0*bb_grad_u[7]*wvpar)+0.5*bb_grad_u[7]*dvpar+p_force_G_1[7]; 

  double alphaSurf_G_1_r[8] = {0.0}; 
  alphaSurf_G_1_r[0] = (-1.0*bb_grad_u[0]*wvpar)-0.5*bb_grad_u[0]*dvpar+p_force_G_1[0]; 
  alphaSurf_G_1_r[1] = (-1.0*bb_grad_u[1]*wvpar)-0.5*bb_grad_u[1]*dvpar+p_force_G_1[1]; 
  alphaSurf_G_1_r[2] = (-1.0*bb_grad_u[2]*wvpar)-0.5*bb_grad_u[2]*dvpar+p_force_G_1[2]; 
  alphaSurf_G_1_r[3] = (-1.0*bb_grad_u[3]*wvpar)-0.5*bb_grad_u[3]*dvpar+p_force_G_1[3]; 
  alphaSurf_G_1_r[4] = (-1.0*bb_grad_u[4]*wvpar)-0.5*bb_grad_u[4]*dvpar+p_force_G_1[4]; 
  alphaSurf_G_1_r[5] = (-1.0*bb_grad_u[5]*wvpar)-0.5*bb_grad_u[5]*dvpar+p_force_G_1[5]; 
  alphaSurf_G_1_r[6] = (-1.0*bb_grad_u[6]*wvpar)-0.5*bb_grad_u[6]*dvpar+p_force_G_1[6]; 
  alphaSurf_G_1_r[7] = (-1.0*bb_grad_u[7]*wvpar)-0.5*bb_grad_u[7]*dvpar+p_force_G_1[7]; 

  double div_b_Surf[8] = {0.0}; 
  div_b_Surf[0] = div_b[0]; 
  div_b_Surf[1] = div_b[1]; 
  div_b_Surf[2] = div_b[2]; 
  div_b_Surf[3] = div_b[3]; 
  div_b_Surf[4] = div_b[4]; 
  div_b_Surf[5] = div_b[5]; 
  div_b_Surf[6] = div_b[6]; 
  div_b_Surf[7] = div_b[7]; 

  double F_0_UpwindQuad_l[8] = {0.0};
  double F_0_UpwindQuad_r[8] = {0.0};
  double F_0_Upwind_l[8] = {0.0};
  double F_0_Upwind_r[8] = {0.0};
  double Ghat_F_0_l[8] = {0.0}; 
  double Ghat_F_0_r[8] = {0.0}; 
  double G_1_UpwindQuad_l[8] = {0.0};
  double G_1_UpwindQuad_r[8] = {0.0};
  double G_1_Upwind_l[8] = {0.0};
  double G_1_Upwind_r[8] = {0.0};
  double Ghat_G_1_l[8] = {0.0}; 
  double Ghat_G_1_r[8] = {0.0}; 

  if ((-0.3535533905932737*alphaSurf_F_0_l[7])+0.3535533905932737*(alphaSurf_F_0_l[6]+alphaSurf_F_0_l[5]+alphaSurf_F_0_l[4])-0.3535533905932737*(alphaSurf_F_0_l[3]+alphaSurf_F_0_l[2]+alphaSurf_F_0_l[1])+0.3535533905932737*alphaSurf_F_0_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(F_0c); 
  } 
  if ((-0.3535533905932737*alphaSurf_F_0_r[7])+0.3535533905932737*(alphaSurf_F_0_r[6]+alphaSurf_F_0_r[5]+alphaSurf_F_0_r[4])-0.3535533905932737*(alphaSurf_F_0_r[3]+alphaSurf_F_0_r[2]+alphaSurf_F_0_r[1])+0.3535533905932737*alphaSurf_F_0_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(F_0r); 
  } 
  if (0.3535533905932737*alphaSurf_F_0_l[7]-0.3535533905932737*(alphaSurf_F_0_l[6]+alphaSurf_F_0_l[5])+0.3535533905932737*(alphaSurf_F_0_l[4]+alphaSurf_F_0_l[3])-0.3535533905932737*(alphaSurf_F_0_l[2]+alphaSurf_F_0_l[1])+0.3535533905932737*alphaSurf_F_0_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(F_0c); 
  } 
  if (0.3535533905932737*alphaSurf_F_0_r[7]-0.3535533905932737*(alphaSurf_F_0_r[6]+alphaSurf_F_0_r[5])+0.3535533905932737*(alphaSurf_F_0_r[4]+alphaSurf_F_0_r[3])-0.3535533905932737*(alphaSurf_F_0_r[2]+alphaSurf_F_0_r[1])+0.3535533905932737*alphaSurf_F_0_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(F_0r); 
  } 
  if (0.3535533905932737*alphaSurf_F_0_l[7]-0.3535533905932737*alphaSurf_F_0_l[6]+0.3535533905932737*alphaSurf_F_0_l[5]-0.3535533905932737*(alphaSurf_F_0_l[4]+alphaSurf_F_0_l[3])+0.3535533905932737*alphaSurf_F_0_l[2]-0.3535533905932737*alphaSurf_F_0_l[1]+0.3535533905932737*alphaSurf_F_0_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(F_0c); 
  } 
  if (0.3535533905932737*alphaSurf_F_0_r[7]-0.3535533905932737*alphaSurf_F_0_r[6]+0.3535533905932737*alphaSurf_F_0_r[5]-0.3535533905932737*(alphaSurf_F_0_r[4]+alphaSurf_F_0_r[3])+0.3535533905932737*alphaSurf_F_0_r[2]-0.3535533905932737*alphaSurf_F_0_r[1]+0.3535533905932737*alphaSurf_F_0_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(F_0r); 
  } 
  if ((-0.3535533905932737*alphaSurf_F_0_l[7])+0.3535533905932737*alphaSurf_F_0_l[6]-0.3535533905932737*(alphaSurf_F_0_l[5]+alphaSurf_F_0_l[4])+0.3535533905932737*(alphaSurf_F_0_l[3]+alphaSurf_F_0_l[2])-0.3535533905932737*alphaSurf_F_0_l[1]+0.3535533905932737*alphaSurf_F_0_l[0] > 0) { 
    F_0_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(F_0c); 
  } 
  if ((-0.3535533905932737*alphaSurf_F_0_r[7])+0.3535533905932737*alphaSurf_F_0_r[6]-0.3535533905932737*(alphaSurf_F_0_r[5]+alphaSurf_F_0_r[4])+0.3535533905932737*(alphaSurf_F_0_r[3]+alphaSurf_F_0_r[2])-0.3535533905932737*alphaSurf_F_0_r[1]+0.3535533905932737*alphaSurf_F_0_r[0] > 0) { 
    F_0_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(F_0r); 
  } 
  if (0.3535533905932737*(alphaSurf_F_0_l[7]+alphaSurf_F_0_l[6])-0.3535533905932737*(alphaSurf_F_0_l[5]+alphaSurf_F_0_l[4]+alphaSurf_F_0_l[3]+alphaSurf_F_0_l[2])+0.3535533905932737*(alphaSurf_F_0_l[1]+alphaSurf_F_0_l[0]) > 0) { 
    F_0_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(F_0c); 
  } 
  if (0.3535533905932737*(alphaSurf_F_0_r[7]+alphaSurf_F_0_r[6])-0.3535533905932737*(alphaSurf_F_0_r[5]+alphaSurf_F_0_r[4]+alphaSurf_F_0_r[3]+alphaSurf_F_0_r[2])+0.3535533905932737*(alphaSurf_F_0_r[1]+alphaSurf_F_0_r[0]) > 0) { 
    F_0_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(F_0r); 
  } 
  if ((-0.3535533905932737*(alphaSurf_F_0_l[7]+alphaSurf_F_0_l[6]))+0.3535533905932737*alphaSurf_F_0_l[5]-0.3535533905932737*alphaSurf_F_0_l[4]+0.3535533905932737*alphaSurf_F_0_l[3]-0.3535533905932737*alphaSurf_F_0_l[2]+0.3535533905932737*(alphaSurf_F_0_l[1]+alphaSurf_F_0_l[0]) > 0) { 
    F_0_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(F_0c); 
  } 
  if ((-0.3535533905932737*(alphaSurf_F_0_r[7]+alphaSurf_F_0_r[6]))+0.3535533905932737*alphaSurf_F_0_r[5]-0.3535533905932737*alphaSurf_F_0_r[4]+0.3535533905932737*alphaSurf_F_0_r[3]-0.3535533905932737*alphaSurf_F_0_r[2]+0.3535533905932737*(alphaSurf_F_0_r[1]+alphaSurf_F_0_r[0]) > 0) { 
    F_0_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(F_0r); 
  } 
  if ((-0.3535533905932737*(alphaSurf_F_0_l[7]+alphaSurf_F_0_l[6]+alphaSurf_F_0_l[5]))+0.3535533905932737*alphaSurf_F_0_l[4]-0.3535533905932737*alphaSurf_F_0_l[3]+0.3535533905932737*(alphaSurf_F_0_l[2]+alphaSurf_F_0_l[1]+alphaSurf_F_0_l[0]) > 0) { 
    F_0_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(F_0c); 
  } 
  if ((-0.3535533905932737*(alphaSurf_F_0_r[7]+alphaSurf_F_0_r[6]+alphaSurf_F_0_r[5]))+0.3535533905932737*alphaSurf_F_0_r[4]-0.3535533905932737*alphaSurf_F_0_r[3]+0.3535533905932737*(alphaSurf_F_0_r[2]+alphaSurf_F_0_r[1]+alphaSurf_F_0_r[0]) > 0) { 
    F_0_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(F_0r); 
  } 
  if (0.3535533905932737*(alphaSurf_F_0_l[7]+alphaSurf_F_0_l[6]+alphaSurf_F_0_l[5]+alphaSurf_F_0_l[4]+alphaSurf_F_0_l[3]+alphaSurf_F_0_l[2]+alphaSurf_F_0_l[1]+alphaSurf_F_0_l[0]) > 0) { 
    F_0_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(F_0c); 
  } 
  if (0.3535533905932737*(alphaSurf_F_0_r[7]+alphaSurf_F_0_r[6]+alphaSurf_F_0_r[5]+alphaSurf_F_0_r[4]+alphaSurf_F_0_r[3]+alphaSurf_F_0_r[2]+alphaSurf_F_0_r[1]+alphaSurf_F_0_r[0]) > 0) { 
    F_0_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(F_0r); 
  } 
  if ((-0.3535533905932737*alphaSurf_G_1_l[7])+0.3535533905932737*(alphaSurf_G_1_l[6]+alphaSurf_G_1_l[5]+alphaSurf_G_1_l[4])-0.3535533905932737*(alphaSurf_G_1_l[3]+alphaSurf_G_1_l[2]+alphaSurf_G_1_l[1])+0.3535533905932737*alphaSurf_G_1_l[0] > 0) { 
    G_1_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(G_1c); 
  } 
  if ((-0.3535533905932737*alphaSurf_G_1_r[7])+0.3535533905932737*(alphaSurf_G_1_r[6]+alphaSurf_G_1_r[5]+alphaSurf_G_1_r[4])-0.3535533905932737*(alphaSurf_G_1_r[3]+alphaSurf_G_1_r[2]+alphaSurf_G_1_r[1])+0.3535533905932737*alphaSurf_G_1_r[0] > 0) { 
    G_1_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(G_1r); 
  } 
  if (0.3535533905932737*alphaSurf_G_1_l[7]-0.3535533905932737*(alphaSurf_G_1_l[6]+alphaSurf_G_1_l[5])+0.3535533905932737*(alphaSurf_G_1_l[4]+alphaSurf_G_1_l[3])-0.3535533905932737*(alphaSurf_G_1_l[2]+alphaSurf_G_1_l[1])+0.3535533905932737*alphaSurf_G_1_l[0] > 0) { 
    G_1_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(G_1c); 
  } 
  if (0.3535533905932737*alphaSurf_G_1_r[7]-0.3535533905932737*(alphaSurf_G_1_r[6]+alphaSurf_G_1_r[5])+0.3535533905932737*(alphaSurf_G_1_r[4]+alphaSurf_G_1_r[3])-0.3535533905932737*(alphaSurf_G_1_r[2]+alphaSurf_G_1_r[1])+0.3535533905932737*alphaSurf_G_1_r[0] > 0) { 
    G_1_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(G_1r); 
  } 
  if (0.3535533905932737*alphaSurf_G_1_l[7]-0.3535533905932737*alphaSurf_G_1_l[6]+0.3535533905932737*alphaSurf_G_1_l[5]-0.3535533905932737*(alphaSurf_G_1_l[4]+alphaSurf_G_1_l[3])+0.3535533905932737*alphaSurf_G_1_l[2]-0.3535533905932737*alphaSurf_G_1_l[1]+0.3535533905932737*alphaSurf_G_1_l[0] > 0) { 
    G_1_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(G_1c); 
  } 
  if (0.3535533905932737*alphaSurf_G_1_r[7]-0.3535533905932737*alphaSurf_G_1_r[6]+0.3535533905932737*alphaSurf_G_1_r[5]-0.3535533905932737*(alphaSurf_G_1_r[4]+alphaSurf_G_1_r[3])+0.3535533905932737*alphaSurf_G_1_r[2]-0.3535533905932737*alphaSurf_G_1_r[1]+0.3535533905932737*alphaSurf_G_1_r[0] > 0) { 
    G_1_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(G_1r); 
  } 
  if ((-0.3535533905932737*alphaSurf_G_1_l[7])+0.3535533905932737*alphaSurf_G_1_l[6]-0.3535533905932737*(alphaSurf_G_1_l[5]+alphaSurf_G_1_l[4])+0.3535533905932737*(alphaSurf_G_1_l[3]+alphaSurf_G_1_l[2])-0.3535533905932737*alphaSurf_G_1_l[1]+0.3535533905932737*alphaSurf_G_1_l[0] > 0) { 
    G_1_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(G_1c); 
  } 
  if ((-0.3535533905932737*alphaSurf_G_1_r[7])+0.3535533905932737*alphaSurf_G_1_r[6]-0.3535533905932737*(alphaSurf_G_1_r[5]+alphaSurf_G_1_r[4])+0.3535533905932737*(alphaSurf_G_1_r[3]+alphaSurf_G_1_r[2])-0.3535533905932737*alphaSurf_G_1_r[1]+0.3535533905932737*alphaSurf_G_1_r[0] > 0) { 
    G_1_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(G_1r); 
  } 
  if (0.3535533905932737*(alphaSurf_G_1_l[7]+alphaSurf_G_1_l[6])-0.3535533905932737*(alphaSurf_G_1_l[5]+alphaSurf_G_1_l[4]+alphaSurf_G_1_l[3]+alphaSurf_G_1_l[2])+0.3535533905932737*(alphaSurf_G_1_l[1]+alphaSurf_G_1_l[0]) > 0) { 
    G_1_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(G_1c); 
  } 
  if (0.3535533905932737*(alphaSurf_G_1_r[7]+alphaSurf_G_1_r[6])-0.3535533905932737*(alphaSurf_G_1_r[5]+alphaSurf_G_1_r[4]+alphaSurf_G_1_r[3]+alphaSurf_G_1_r[2])+0.3535533905932737*(alphaSurf_G_1_r[1]+alphaSurf_G_1_r[0]) > 0) { 
    G_1_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(G_1r); 
  } 
  if ((-0.3535533905932737*(alphaSurf_G_1_l[7]+alphaSurf_G_1_l[6]))+0.3535533905932737*alphaSurf_G_1_l[5]-0.3535533905932737*alphaSurf_G_1_l[4]+0.3535533905932737*alphaSurf_G_1_l[3]-0.3535533905932737*alphaSurf_G_1_l[2]+0.3535533905932737*(alphaSurf_G_1_l[1]+alphaSurf_G_1_l[0]) > 0) { 
    G_1_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(G_1c); 
  } 
  if ((-0.3535533905932737*(alphaSurf_G_1_r[7]+alphaSurf_G_1_r[6]))+0.3535533905932737*alphaSurf_G_1_r[5]-0.3535533905932737*alphaSurf_G_1_r[4]+0.3535533905932737*alphaSurf_G_1_r[3]-0.3535533905932737*alphaSurf_G_1_r[2]+0.3535533905932737*(alphaSurf_G_1_r[1]+alphaSurf_G_1_r[0]) > 0) { 
    G_1_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(G_1r); 
  } 
  if ((-0.3535533905932737*(alphaSurf_G_1_l[7]+alphaSurf_G_1_l[6]+alphaSurf_G_1_l[5]))+0.3535533905932737*alphaSurf_G_1_l[4]-0.3535533905932737*alphaSurf_G_1_l[3]+0.3535533905932737*(alphaSurf_G_1_l[2]+alphaSurf_G_1_l[1]+alphaSurf_G_1_l[0]) > 0) { 
    G_1_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(G_1c); 
  } 
  if ((-0.3535533905932737*(alphaSurf_G_1_r[7]+alphaSurf_G_1_r[6]+alphaSurf_G_1_r[5]))+0.3535533905932737*alphaSurf_G_1_r[4]-0.3535533905932737*alphaSurf_G_1_r[3]+0.3535533905932737*(alphaSurf_G_1_r[2]+alphaSurf_G_1_r[1]+alphaSurf_G_1_r[0]) > 0) { 
    G_1_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(G_1r); 
  } 
  if (0.3535533905932737*(alphaSurf_G_1_l[7]+alphaSurf_G_1_l[6]+alphaSurf_G_1_l[5]+alphaSurf_G_1_l[4]+alphaSurf_G_1_l[3]+alphaSurf_G_1_l[2]+alphaSurf_G_1_l[1]+alphaSurf_G_1_l[0]) > 0) { 
    G_1_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(G_1l); 
  } else { 
    G_1_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(G_1c); 
  } 
  if (0.3535533905932737*(alphaSurf_G_1_r[7]+alphaSurf_G_1_r[6]+alphaSurf_G_1_r[5]+alphaSurf_G_1_r[4]+alphaSurf_G_1_r[3]+alphaSurf_G_1_r[2]+alphaSurf_G_1_r[1]+alphaSurf_G_1_r[0]) > 0) { 
    G_1_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(G_1c); 
  } else { 
    G_1_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  double F_0_div_b_UpwindQuad_l[8] = {0.0};
  double F_0_div_b_UpwindQuad_r[8] = {0.0};
  double F_0_div_b_Upwind_l[8] = {0.0};
  double F_0_div_b_Upwind_r[8] = {0.0};
  double Ghat_F_0_div_b_l[8] = {0.0}; 
  double Ghat_F_0_div_b_r[8] = {0.0}; 

  double G_1_div_b_UpwindQuad_l[8] = {0.0};
  double G_1_div_b_UpwindQuad_r[8] = {0.0};
  double G_1_div_b_Upwind_l[8] = {0.0};
  double G_1_div_b_Upwind_r[8] = {0.0};
  double Ghat_G_1_div_b_l[8] = {0.0}; 
  double Ghat_G_1_div_b_r[8] = {0.0}; 

  if ((-0.3535533905932737*div_b_Surf[7])+0.3535533905932737*(div_b_Surf[6]+div_b_Surf[5]+div_b_Surf[4])-0.3535533905932737*(div_b_Surf[3]+div_b_Surf[2]+div_b_Surf[1])+0.3535533905932737*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[0] = hyb_3x1v_p1_surfx4_eval_quad_node_0_l(G_1_sourcer); 
  } 
  if (0.3535533905932737*div_b_Surf[7]-0.3535533905932737*(div_b_Surf[6]+div_b_Surf[5])+0.3535533905932737*(div_b_Surf[4]+div_b_Surf[3])-0.3535533905932737*(div_b_Surf[2]+div_b_Surf[1])+0.3535533905932737*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[1] = hyb_3x1v_p1_surfx4_eval_quad_node_1_l(G_1_sourcer); 
  } 
  if (0.3535533905932737*div_b_Surf[7]-0.3535533905932737*div_b_Surf[6]+0.3535533905932737*div_b_Surf[5]-0.3535533905932737*(div_b_Surf[4]+div_b_Surf[3])+0.3535533905932737*div_b_Surf[2]-0.3535533905932737*div_b_Surf[1]+0.3535533905932737*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[2] = hyb_3x1v_p1_surfx4_eval_quad_node_2_l(G_1_sourcer); 
  } 
  if ((-0.3535533905932737*div_b_Surf[7])+0.3535533905932737*div_b_Surf[6]-0.3535533905932737*(div_b_Surf[5]+div_b_Surf[4])+0.3535533905932737*(div_b_Surf[3]+div_b_Surf[2])-0.3535533905932737*div_b_Surf[1]+0.3535533905932737*div_b_Surf[0] > 0) { 
    F_0_div_b_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[3] = hyb_3x1v_p1_surfx4_eval_quad_node_3_l(G_1_sourcer); 
  } 
  if (0.3535533905932737*(div_b_Surf[7]+div_b_Surf[6])-0.3535533905932737*(div_b_Surf[5]+div_b_Surf[4]+div_b_Surf[3]+div_b_Surf[2])+0.3535533905932737*(div_b_Surf[1]+div_b_Surf[0]) > 0) { 
    F_0_div_b_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[4] = hyb_3x1v_p1_surfx4_eval_quad_node_4_l(G_1_sourcer); 
  } 
  if ((-0.3535533905932737*(div_b_Surf[7]+div_b_Surf[6]))+0.3535533905932737*div_b_Surf[5]-0.3535533905932737*div_b_Surf[4]+0.3535533905932737*div_b_Surf[3]-0.3535533905932737*div_b_Surf[2]+0.3535533905932737*(div_b_Surf[1]+div_b_Surf[0]) > 0) { 
    F_0_div_b_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[5] = hyb_3x1v_p1_surfx4_eval_quad_node_5_l(G_1_sourcer); 
  } 
  if ((-0.3535533905932737*(div_b_Surf[7]+div_b_Surf[6]+div_b_Surf[5]))+0.3535533905932737*div_b_Surf[4]-0.3535533905932737*div_b_Surf[3]+0.3535533905932737*(div_b_Surf[2]+div_b_Surf[1]+div_b_Surf[0]) > 0) { 
    F_0_div_b_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[6] = hyb_3x1v_p1_surfx4_eval_quad_node_6_l(G_1_sourcer); 
  } 
  if (0.3535533905932737*(div_b_Surf[7]+div_b_Surf[6]+div_b_Surf[5]+div_b_Surf[4]+div_b_Surf[3]+div_b_Surf[2]+div_b_Surf[1]+div_b_Surf[0]) > 0) { 
    F_0_div_b_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(F_0_sourcel); 
    F_0_div_b_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(F_0_sourcec); 
    G_1_div_b_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(G_1_sourcel); 
    G_1_div_b_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_r(G_1_sourcec); 
  } else { 
    F_0_div_b_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(F_0_sourcec); 
    F_0_div_b_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(F_0_sourcer); 
    G_1_div_b_UpwindQuad_l[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(G_1_sourcec); 
    G_1_div_b_UpwindQuad_r[7] = hyb_3x1v_p1_surfx4_eval_quad_node_7_l(G_1_sourcer); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(F_0_div_b_UpwindQuad_l, F_0_div_b_Upwind_l); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(F_0_div_b_UpwindQuad_r, F_0_div_b_Upwind_r); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(G_1_div_b_UpwindQuad_l, G_1_div_b_Upwind_l); 
  hyb_3x1v_p1_vdir_upwind_quad_to_modal(G_1_div_b_UpwindQuad_r, G_1_div_b_Upwind_r); 

  Ghat_F_0_l[0] = 0.3535533905932737*F_0_Upwind_l[7]*alphaSurf_F_0_l[7]+0.3535533905932737*F_0_Upwind_l[6]*alphaSurf_F_0_l[6]+0.3535533905932737*F_0_Upwind_l[5]*alphaSurf_F_0_l[5]+0.3535533905932737*F_0_Upwind_l[4]*alphaSurf_F_0_l[4]+0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_F_0_l[3]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_F_0_l[2]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_F_0_l[1]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_F_0_l[0]; 
  Ghat_F_0_l[1] = 0.3535533905932737*F_0_Upwind_l[6]*alphaSurf_F_0_l[7]+0.3535533905932737*alphaSurf_F_0_l[6]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_F_0_l[5]+0.3535533905932737*alphaSurf_F_0_l[3]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_F_0_l[4]+0.3535533905932737*alphaSurf_F_0_l[2]*F_0_Upwind_l[4]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_F_0_l[1]+0.3535533905932737*alphaSurf_F_0_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.3535533905932737*F_0_Upwind_l[5]*alphaSurf_F_0_l[7]+0.3535533905932737*alphaSurf_F_0_l[5]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_F_0_l[6]+0.3535533905932737*alphaSurf_F_0_l[3]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_F_0_l[4]+0.3535533905932737*alphaSurf_F_0_l[1]*F_0_Upwind_l[4]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_F_0_l[2]+0.3535533905932737*alphaSurf_F_0_l[0]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[3] = 0.3535533905932737*F_0_Upwind_l[4]*alphaSurf_F_0_l[7]+0.3535533905932737*alphaSurf_F_0_l[4]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_F_0_l[6]+0.3535533905932737*alphaSurf_F_0_l[2]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_F_0_l[5]+0.3535533905932737*alphaSurf_F_0_l[1]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_F_0_l[3]+0.3535533905932737*alphaSurf_F_0_l[0]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[4] = 0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_F_0_l[7]+0.3535533905932737*alphaSurf_F_0_l[3]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[5]*alphaSurf_F_0_l[6]+0.3535533905932737*alphaSurf_F_0_l[5]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_F_0_l[4]+0.3535533905932737*alphaSurf_F_0_l[0]*F_0_Upwind_l[4]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_F_0_l[2]+0.3535533905932737*alphaSurf_F_0_l[1]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[5] = 0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_F_0_l[7]+0.3535533905932737*alphaSurf_F_0_l[2]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[4]*alphaSurf_F_0_l[6]+0.3535533905932737*alphaSurf_F_0_l[4]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_F_0_l[5]+0.3535533905932737*alphaSurf_F_0_l[0]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_F_0_l[3]+0.3535533905932737*alphaSurf_F_0_l[1]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[6] = 0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_F_0_l[7]+0.3535533905932737*alphaSurf_F_0_l[1]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_F_0_l[6]+0.3535533905932737*alphaSurf_F_0_l[0]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[4]*alphaSurf_F_0_l[5]+0.3535533905932737*alphaSurf_F_0_l[4]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_F_0_l[3]+0.3535533905932737*alphaSurf_F_0_l[2]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[7] = 0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_F_0_l[7]+0.3535533905932737*alphaSurf_F_0_l[0]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_F_0_l[6]+0.3535533905932737*alphaSurf_F_0_l[1]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_F_0_l[5]+0.3535533905932737*alphaSurf_F_0_l[2]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_F_0_l[4]+0.3535533905932737*alphaSurf_F_0_l[3]*F_0_Upwind_l[4]; 
  Ghat_G_1_l[0] = 0.3535533905932737*G_1_Upwind_l[7]*alphaSurf_G_1_l[7]+0.3535533905932737*G_1_Upwind_l[6]*alphaSurf_G_1_l[6]+0.3535533905932737*G_1_Upwind_l[5]*alphaSurf_G_1_l[5]+0.3535533905932737*G_1_Upwind_l[4]*alphaSurf_G_1_l[4]+0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_G_1_l[3]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_G_1_l[2]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_G_1_l[1]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_G_1_l[0]; 
  Ghat_G_1_l[1] = 0.3535533905932737*G_1_Upwind_l[6]*alphaSurf_G_1_l[7]+0.3535533905932737*alphaSurf_G_1_l[6]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_G_1_l[5]+0.3535533905932737*alphaSurf_G_1_l[3]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_G_1_l[4]+0.3535533905932737*alphaSurf_G_1_l[2]*G_1_Upwind_l[4]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_G_1_l[1]+0.3535533905932737*alphaSurf_G_1_l[0]*G_1_Upwind_l[1]; 
  Ghat_G_1_l[2] = 0.3535533905932737*G_1_Upwind_l[5]*alphaSurf_G_1_l[7]+0.3535533905932737*alphaSurf_G_1_l[5]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_G_1_l[6]+0.3535533905932737*alphaSurf_G_1_l[3]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_G_1_l[4]+0.3535533905932737*alphaSurf_G_1_l[1]*G_1_Upwind_l[4]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_G_1_l[2]+0.3535533905932737*alphaSurf_G_1_l[0]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[3] = 0.3535533905932737*G_1_Upwind_l[4]*alphaSurf_G_1_l[7]+0.3535533905932737*alphaSurf_G_1_l[4]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_G_1_l[6]+0.3535533905932737*alphaSurf_G_1_l[2]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_G_1_l[5]+0.3535533905932737*alphaSurf_G_1_l[1]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_G_1_l[3]+0.3535533905932737*alphaSurf_G_1_l[0]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[4] = 0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_G_1_l[7]+0.3535533905932737*alphaSurf_G_1_l[3]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[5]*alphaSurf_G_1_l[6]+0.3535533905932737*alphaSurf_G_1_l[5]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_G_1_l[4]+0.3535533905932737*alphaSurf_G_1_l[0]*G_1_Upwind_l[4]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_G_1_l[2]+0.3535533905932737*alphaSurf_G_1_l[1]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[5] = 0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_G_1_l[7]+0.3535533905932737*alphaSurf_G_1_l[2]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[4]*alphaSurf_G_1_l[6]+0.3535533905932737*alphaSurf_G_1_l[4]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_G_1_l[5]+0.3535533905932737*alphaSurf_G_1_l[0]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_G_1_l[3]+0.3535533905932737*alphaSurf_G_1_l[1]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[6] = 0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_G_1_l[7]+0.3535533905932737*alphaSurf_G_1_l[1]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_G_1_l[6]+0.3535533905932737*alphaSurf_G_1_l[0]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[4]*alphaSurf_G_1_l[5]+0.3535533905932737*alphaSurf_G_1_l[4]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_G_1_l[3]+0.3535533905932737*alphaSurf_G_1_l[2]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[7] = 0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_G_1_l[7]+0.3535533905932737*alphaSurf_G_1_l[0]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_G_1_l[6]+0.3535533905932737*alphaSurf_G_1_l[1]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_G_1_l[5]+0.3535533905932737*alphaSurf_G_1_l[2]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_G_1_l[4]+0.3535533905932737*alphaSurf_G_1_l[3]*G_1_Upwind_l[4]; 
  Ghat_F_0_div_b_l[0] = 0.3535533905932737*F_0_div_b_Upwind_l[7]*div_b_Surf[7]+0.3535533905932737*F_0_div_b_Upwind_l[6]*div_b_Surf[6]+0.3535533905932737*F_0_div_b_Upwind_l[5]*div_b_Surf[5]+0.3535533905932737*F_0_div_b_Upwind_l[4]*div_b_Surf[4]+0.3535533905932737*F_0_div_b_Upwind_l[3]*div_b_Surf[3]+0.3535533905932737*F_0_div_b_Upwind_l[2]*div_b_Surf[2]+0.3535533905932737*F_0_div_b_Upwind_l[1]*div_b_Surf[1]+0.3535533905932737*F_0_div_b_Upwind_l[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b_l[1] = 0.3535533905932737*F_0_div_b_Upwind_l[6]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[6]*F_0_div_b_Upwind_l[7]+0.3535533905932737*F_0_div_b_Upwind_l[3]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[3]*F_0_div_b_Upwind_l[5]+0.3535533905932737*F_0_div_b_Upwind_l[2]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_l[4]+0.3535533905932737*F_0_div_b_Upwind_l[0]*div_b_Surf[1]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_l[1]; 
  Ghat_F_0_div_b_l[2] = 0.3535533905932737*F_0_div_b_Upwind_l[5]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[5]*F_0_div_b_Upwind_l[7]+0.3535533905932737*F_0_div_b_Upwind_l[3]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[3]*F_0_div_b_Upwind_l[6]+0.3535533905932737*F_0_div_b_Upwind_l[1]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_l[4]+0.3535533905932737*F_0_div_b_Upwind_l[0]*div_b_Surf[2]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_l[2]; 
  Ghat_F_0_div_b_l[3] = 0.3535533905932737*F_0_div_b_Upwind_l[4]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[4]*F_0_div_b_Upwind_l[7]+0.3535533905932737*F_0_div_b_Upwind_l[2]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_l[6]+0.3535533905932737*F_0_div_b_Upwind_l[1]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_l[5]+0.3535533905932737*F_0_div_b_Upwind_l[0]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_l[3]; 
  Ghat_F_0_div_b_l[4] = 0.3535533905932737*F_0_div_b_Upwind_l[3]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[3]*F_0_div_b_Upwind_l[7]+0.3535533905932737*F_0_div_b_Upwind_l[5]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[5]*F_0_div_b_Upwind_l[6]+0.3535533905932737*F_0_div_b_Upwind_l[0]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_l[4]+0.3535533905932737*F_0_div_b_Upwind_l[1]*div_b_Surf[2]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_l[2]; 
  Ghat_F_0_div_b_l[5] = 0.3535533905932737*F_0_div_b_Upwind_l[2]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_l[7]+0.3535533905932737*F_0_div_b_Upwind_l[4]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[4]*F_0_div_b_Upwind_l[6]+0.3535533905932737*F_0_div_b_Upwind_l[0]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_l[5]+0.3535533905932737*F_0_div_b_Upwind_l[1]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_l[3]; 
  Ghat_F_0_div_b_l[6] = 0.3535533905932737*F_0_div_b_Upwind_l[1]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_l[7]+0.3535533905932737*F_0_div_b_Upwind_l[0]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_l[6]+0.3535533905932737*F_0_div_b_Upwind_l[4]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[4]*F_0_div_b_Upwind_l[5]+0.3535533905932737*F_0_div_b_Upwind_l[2]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_l[3]; 
  Ghat_F_0_div_b_l[7] = 0.3535533905932737*F_0_div_b_Upwind_l[0]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_l[7]+0.3535533905932737*F_0_div_b_Upwind_l[1]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_l[6]+0.3535533905932737*F_0_div_b_Upwind_l[2]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_l[5]+0.3535533905932737*F_0_div_b_Upwind_l[3]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[3]*F_0_div_b_Upwind_l[4]; 
  Ghat_G_1_div_b_l[0] = 0.3535533905932737*G_1_div_b_Upwind_l[7]*div_b_Surf[7]+0.3535533905932737*G_1_div_b_Upwind_l[6]*div_b_Surf[6]+0.3535533905932737*G_1_div_b_Upwind_l[5]*div_b_Surf[5]+0.3535533905932737*G_1_div_b_Upwind_l[4]*div_b_Surf[4]+0.3535533905932737*G_1_div_b_Upwind_l[3]*div_b_Surf[3]+0.3535533905932737*G_1_div_b_Upwind_l[2]*div_b_Surf[2]+0.3535533905932737*G_1_div_b_Upwind_l[1]*div_b_Surf[1]+0.3535533905932737*G_1_div_b_Upwind_l[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b_l[1] = 0.3535533905932737*G_1_div_b_Upwind_l[6]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[6]*G_1_div_b_Upwind_l[7]+0.3535533905932737*G_1_div_b_Upwind_l[3]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[3]*G_1_div_b_Upwind_l[5]+0.3535533905932737*G_1_div_b_Upwind_l[2]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_l[4]+0.3535533905932737*G_1_div_b_Upwind_l[0]*div_b_Surf[1]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_l[1]; 
  Ghat_G_1_div_b_l[2] = 0.3535533905932737*G_1_div_b_Upwind_l[5]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[5]*G_1_div_b_Upwind_l[7]+0.3535533905932737*G_1_div_b_Upwind_l[3]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[3]*G_1_div_b_Upwind_l[6]+0.3535533905932737*G_1_div_b_Upwind_l[1]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_l[4]+0.3535533905932737*G_1_div_b_Upwind_l[0]*div_b_Surf[2]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_l[2]; 
  Ghat_G_1_div_b_l[3] = 0.3535533905932737*G_1_div_b_Upwind_l[4]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[4]*G_1_div_b_Upwind_l[7]+0.3535533905932737*G_1_div_b_Upwind_l[2]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_l[6]+0.3535533905932737*G_1_div_b_Upwind_l[1]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_l[5]+0.3535533905932737*G_1_div_b_Upwind_l[0]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_l[3]; 
  Ghat_G_1_div_b_l[4] = 0.3535533905932737*G_1_div_b_Upwind_l[3]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[3]*G_1_div_b_Upwind_l[7]+0.3535533905932737*G_1_div_b_Upwind_l[5]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[5]*G_1_div_b_Upwind_l[6]+0.3535533905932737*G_1_div_b_Upwind_l[0]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_l[4]+0.3535533905932737*G_1_div_b_Upwind_l[1]*div_b_Surf[2]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_l[2]; 
  Ghat_G_1_div_b_l[5] = 0.3535533905932737*G_1_div_b_Upwind_l[2]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_l[7]+0.3535533905932737*G_1_div_b_Upwind_l[4]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[4]*G_1_div_b_Upwind_l[6]+0.3535533905932737*G_1_div_b_Upwind_l[0]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_l[5]+0.3535533905932737*G_1_div_b_Upwind_l[1]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_l[3]; 
  Ghat_G_1_div_b_l[6] = 0.3535533905932737*G_1_div_b_Upwind_l[1]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_l[7]+0.3535533905932737*G_1_div_b_Upwind_l[0]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_l[6]+0.3535533905932737*G_1_div_b_Upwind_l[4]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[4]*G_1_div_b_Upwind_l[5]+0.3535533905932737*G_1_div_b_Upwind_l[2]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_l[3]; 
  Ghat_G_1_div_b_l[7] = 0.3535533905932737*G_1_div_b_Upwind_l[0]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_l[7]+0.3535533905932737*G_1_div_b_Upwind_l[1]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_l[6]+0.3535533905932737*G_1_div_b_Upwind_l[2]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_l[5]+0.3535533905932737*G_1_div_b_Upwind_l[3]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[3]*G_1_div_b_Upwind_l[4]; 

  Ghat_F_0_r[0] = 0.3535533905932737*F_0_Upwind_r[7]*alphaSurf_F_0_r[7]+0.3535533905932737*F_0_Upwind_r[6]*alphaSurf_F_0_r[6]+0.3535533905932737*F_0_Upwind_r[5]*alphaSurf_F_0_r[5]+0.3535533905932737*F_0_Upwind_r[4]*alphaSurf_F_0_r[4]+0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_F_0_r[3]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_F_0_r[2]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_F_0_r[1]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_F_0_r[0]; 
  Ghat_F_0_r[1] = 0.3535533905932737*F_0_Upwind_r[6]*alphaSurf_F_0_r[7]+0.3535533905932737*alphaSurf_F_0_r[6]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_F_0_r[5]+0.3535533905932737*alphaSurf_F_0_r[3]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_F_0_r[4]+0.3535533905932737*alphaSurf_F_0_r[2]*F_0_Upwind_r[4]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_F_0_r[1]+0.3535533905932737*alphaSurf_F_0_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.3535533905932737*F_0_Upwind_r[5]*alphaSurf_F_0_r[7]+0.3535533905932737*alphaSurf_F_0_r[5]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_F_0_r[6]+0.3535533905932737*alphaSurf_F_0_r[3]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_F_0_r[4]+0.3535533905932737*alphaSurf_F_0_r[1]*F_0_Upwind_r[4]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_F_0_r[2]+0.3535533905932737*alphaSurf_F_0_r[0]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[3] = 0.3535533905932737*F_0_Upwind_r[4]*alphaSurf_F_0_r[7]+0.3535533905932737*alphaSurf_F_0_r[4]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_F_0_r[6]+0.3535533905932737*alphaSurf_F_0_r[2]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_F_0_r[5]+0.3535533905932737*alphaSurf_F_0_r[1]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_F_0_r[3]+0.3535533905932737*alphaSurf_F_0_r[0]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[4] = 0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_F_0_r[7]+0.3535533905932737*alphaSurf_F_0_r[3]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[5]*alphaSurf_F_0_r[6]+0.3535533905932737*alphaSurf_F_0_r[5]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_F_0_r[4]+0.3535533905932737*alphaSurf_F_0_r[0]*F_0_Upwind_r[4]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_F_0_r[2]+0.3535533905932737*alphaSurf_F_0_r[1]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[5] = 0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_F_0_r[7]+0.3535533905932737*alphaSurf_F_0_r[2]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[4]*alphaSurf_F_0_r[6]+0.3535533905932737*alphaSurf_F_0_r[4]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_F_0_r[5]+0.3535533905932737*alphaSurf_F_0_r[0]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_F_0_r[3]+0.3535533905932737*alphaSurf_F_0_r[1]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[6] = 0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_F_0_r[7]+0.3535533905932737*alphaSurf_F_0_r[1]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_F_0_r[6]+0.3535533905932737*alphaSurf_F_0_r[0]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[4]*alphaSurf_F_0_r[5]+0.3535533905932737*alphaSurf_F_0_r[4]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_F_0_r[3]+0.3535533905932737*alphaSurf_F_0_r[2]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[7] = 0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_F_0_r[7]+0.3535533905932737*alphaSurf_F_0_r[0]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_F_0_r[6]+0.3535533905932737*alphaSurf_F_0_r[1]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_F_0_r[5]+0.3535533905932737*alphaSurf_F_0_r[2]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_F_0_r[4]+0.3535533905932737*alphaSurf_F_0_r[3]*F_0_Upwind_r[4]; 
  Ghat_G_1_r[0] = 0.3535533905932737*G_1_Upwind_r[7]*alphaSurf_G_1_r[7]+0.3535533905932737*G_1_Upwind_r[6]*alphaSurf_G_1_r[6]+0.3535533905932737*G_1_Upwind_r[5]*alphaSurf_G_1_r[5]+0.3535533905932737*G_1_Upwind_r[4]*alphaSurf_G_1_r[4]+0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_G_1_r[3]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_G_1_r[2]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_G_1_r[1]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_G_1_r[0]; 
  Ghat_G_1_r[1] = 0.3535533905932737*G_1_Upwind_r[6]*alphaSurf_G_1_r[7]+0.3535533905932737*alphaSurf_G_1_r[6]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_G_1_r[5]+0.3535533905932737*alphaSurf_G_1_r[3]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_G_1_r[4]+0.3535533905932737*alphaSurf_G_1_r[2]*G_1_Upwind_r[4]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_G_1_r[1]+0.3535533905932737*alphaSurf_G_1_r[0]*G_1_Upwind_r[1]; 
  Ghat_G_1_r[2] = 0.3535533905932737*G_1_Upwind_r[5]*alphaSurf_G_1_r[7]+0.3535533905932737*alphaSurf_G_1_r[5]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_G_1_r[6]+0.3535533905932737*alphaSurf_G_1_r[3]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_G_1_r[4]+0.3535533905932737*alphaSurf_G_1_r[1]*G_1_Upwind_r[4]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_G_1_r[2]+0.3535533905932737*alphaSurf_G_1_r[0]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[3] = 0.3535533905932737*G_1_Upwind_r[4]*alphaSurf_G_1_r[7]+0.3535533905932737*alphaSurf_G_1_r[4]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_G_1_r[6]+0.3535533905932737*alphaSurf_G_1_r[2]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_G_1_r[5]+0.3535533905932737*alphaSurf_G_1_r[1]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_G_1_r[3]+0.3535533905932737*alphaSurf_G_1_r[0]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[4] = 0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_G_1_r[7]+0.3535533905932737*alphaSurf_G_1_r[3]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[5]*alphaSurf_G_1_r[6]+0.3535533905932737*alphaSurf_G_1_r[5]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_G_1_r[4]+0.3535533905932737*alphaSurf_G_1_r[0]*G_1_Upwind_r[4]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_G_1_r[2]+0.3535533905932737*alphaSurf_G_1_r[1]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[5] = 0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_G_1_r[7]+0.3535533905932737*alphaSurf_G_1_r[2]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[4]*alphaSurf_G_1_r[6]+0.3535533905932737*alphaSurf_G_1_r[4]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_G_1_r[5]+0.3535533905932737*alphaSurf_G_1_r[0]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_G_1_r[3]+0.3535533905932737*alphaSurf_G_1_r[1]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[6] = 0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_G_1_r[7]+0.3535533905932737*alphaSurf_G_1_r[1]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_G_1_r[6]+0.3535533905932737*alphaSurf_G_1_r[0]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[4]*alphaSurf_G_1_r[5]+0.3535533905932737*alphaSurf_G_1_r[4]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_G_1_r[3]+0.3535533905932737*alphaSurf_G_1_r[2]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[7] = 0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_G_1_r[7]+0.3535533905932737*alphaSurf_G_1_r[0]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_G_1_r[6]+0.3535533905932737*alphaSurf_G_1_r[1]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_G_1_r[5]+0.3535533905932737*alphaSurf_G_1_r[2]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_G_1_r[4]+0.3535533905932737*alphaSurf_G_1_r[3]*G_1_Upwind_r[4]; 
  Ghat_F_0_div_b_r[0] = 0.3535533905932737*F_0_div_b_Upwind_r[7]*div_b_Surf[7]+0.3535533905932737*F_0_div_b_Upwind_r[6]*div_b_Surf[6]+0.3535533905932737*F_0_div_b_Upwind_r[5]*div_b_Surf[5]+0.3535533905932737*F_0_div_b_Upwind_r[4]*div_b_Surf[4]+0.3535533905932737*F_0_div_b_Upwind_r[3]*div_b_Surf[3]+0.3535533905932737*F_0_div_b_Upwind_r[2]*div_b_Surf[2]+0.3535533905932737*F_0_div_b_Upwind_r[1]*div_b_Surf[1]+0.3535533905932737*F_0_div_b_Upwind_r[0]*div_b_Surf[0]; 
  Ghat_F_0_div_b_r[1] = 0.3535533905932737*F_0_div_b_Upwind_r[6]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[6]*F_0_div_b_Upwind_r[7]+0.3535533905932737*F_0_div_b_Upwind_r[3]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[3]*F_0_div_b_Upwind_r[5]+0.3535533905932737*F_0_div_b_Upwind_r[2]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_r[4]+0.3535533905932737*F_0_div_b_Upwind_r[0]*div_b_Surf[1]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_r[1]; 
  Ghat_F_0_div_b_r[2] = 0.3535533905932737*F_0_div_b_Upwind_r[5]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[5]*F_0_div_b_Upwind_r[7]+0.3535533905932737*F_0_div_b_Upwind_r[3]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[3]*F_0_div_b_Upwind_r[6]+0.3535533905932737*F_0_div_b_Upwind_r[1]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_r[4]+0.3535533905932737*F_0_div_b_Upwind_r[0]*div_b_Surf[2]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_r[2]; 
  Ghat_F_0_div_b_r[3] = 0.3535533905932737*F_0_div_b_Upwind_r[4]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[4]*F_0_div_b_Upwind_r[7]+0.3535533905932737*F_0_div_b_Upwind_r[2]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_r[6]+0.3535533905932737*F_0_div_b_Upwind_r[1]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_r[5]+0.3535533905932737*F_0_div_b_Upwind_r[0]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_r[3]; 
  Ghat_F_0_div_b_r[4] = 0.3535533905932737*F_0_div_b_Upwind_r[3]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[3]*F_0_div_b_Upwind_r[7]+0.3535533905932737*F_0_div_b_Upwind_r[5]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[5]*F_0_div_b_Upwind_r[6]+0.3535533905932737*F_0_div_b_Upwind_r[0]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_r[4]+0.3535533905932737*F_0_div_b_Upwind_r[1]*div_b_Surf[2]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_r[2]; 
  Ghat_F_0_div_b_r[5] = 0.3535533905932737*F_0_div_b_Upwind_r[2]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_r[7]+0.3535533905932737*F_0_div_b_Upwind_r[4]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[4]*F_0_div_b_Upwind_r[6]+0.3535533905932737*F_0_div_b_Upwind_r[0]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_r[5]+0.3535533905932737*F_0_div_b_Upwind_r[1]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_r[3]; 
  Ghat_F_0_div_b_r[6] = 0.3535533905932737*F_0_div_b_Upwind_r[1]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_r[7]+0.3535533905932737*F_0_div_b_Upwind_r[0]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_r[6]+0.3535533905932737*F_0_div_b_Upwind_r[4]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[4]*F_0_div_b_Upwind_r[5]+0.3535533905932737*F_0_div_b_Upwind_r[2]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_r[3]; 
  Ghat_F_0_div_b_r[7] = 0.3535533905932737*F_0_div_b_Upwind_r[0]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[0]*F_0_div_b_Upwind_r[7]+0.3535533905932737*F_0_div_b_Upwind_r[1]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[1]*F_0_div_b_Upwind_r[6]+0.3535533905932737*F_0_div_b_Upwind_r[2]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[2]*F_0_div_b_Upwind_r[5]+0.3535533905932737*F_0_div_b_Upwind_r[3]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[3]*F_0_div_b_Upwind_r[4]; 
  Ghat_G_1_div_b_r[0] = 0.3535533905932737*G_1_div_b_Upwind_r[7]*div_b_Surf[7]+0.3535533905932737*G_1_div_b_Upwind_r[6]*div_b_Surf[6]+0.3535533905932737*G_1_div_b_Upwind_r[5]*div_b_Surf[5]+0.3535533905932737*G_1_div_b_Upwind_r[4]*div_b_Surf[4]+0.3535533905932737*G_1_div_b_Upwind_r[3]*div_b_Surf[3]+0.3535533905932737*G_1_div_b_Upwind_r[2]*div_b_Surf[2]+0.3535533905932737*G_1_div_b_Upwind_r[1]*div_b_Surf[1]+0.3535533905932737*G_1_div_b_Upwind_r[0]*div_b_Surf[0]; 
  Ghat_G_1_div_b_r[1] = 0.3535533905932737*G_1_div_b_Upwind_r[6]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[6]*G_1_div_b_Upwind_r[7]+0.3535533905932737*G_1_div_b_Upwind_r[3]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[3]*G_1_div_b_Upwind_r[5]+0.3535533905932737*G_1_div_b_Upwind_r[2]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_r[4]+0.3535533905932737*G_1_div_b_Upwind_r[0]*div_b_Surf[1]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_r[1]; 
  Ghat_G_1_div_b_r[2] = 0.3535533905932737*G_1_div_b_Upwind_r[5]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[5]*G_1_div_b_Upwind_r[7]+0.3535533905932737*G_1_div_b_Upwind_r[3]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[3]*G_1_div_b_Upwind_r[6]+0.3535533905932737*G_1_div_b_Upwind_r[1]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_r[4]+0.3535533905932737*G_1_div_b_Upwind_r[0]*div_b_Surf[2]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_r[2]; 
  Ghat_G_1_div_b_r[3] = 0.3535533905932737*G_1_div_b_Upwind_r[4]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[4]*G_1_div_b_Upwind_r[7]+0.3535533905932737*G_1_div_b_Upwind_r[2]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_r[6]+0.3535533905932737*G_1_div_b_Upwind_r[1]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_r[5]+0.3535533905932737*G_1_div_b_Upwind_r[0]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_r[3]; 
  Ghat_G_1_div_b_r[4] = 0.3535533905932737*G_1_div_b_Upwind_r[3]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[3]*G_1_div_b_Upwind_r[7]+0.3535533905932737*G_1_div_b_Upwind_r[5]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[5]*G_1_div_b_Upwind_r[6]+0.3535533905932737*G_1_div_b_Upwind_r[0]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_r[4]+0.3535533905932737*G_1_div_b_Upwind_r[1]*div_b_Surf[2]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_r[2]; 
  Ghat_G_1_div_b_r[5] = 0.3535533905932737*G_1_div_b_Upwind_r[2]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_r[7]+0.3535533905932737*G_1_div_b_Upwind_r[4]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[4]*G_1_div_b_Upwind_r[6]+0.3535533905932737*G_1_div_b_Upwind_r[0]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_r[5]+0.3535533905932737*G_1_div_b_Upwind_r[1]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_r[3]; 
  Ghat_G_1_div_b_r[6] = 0.3535533905932737*G_1_div_b_Upwind_r[1]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_r[7]+0.3535533905932737*G_1_div_b_Upwind_r[0]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_r[6]+0.3535533905932737*G_1_div_b_Upwind_r[4]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[4]*G_1_div_b_Upwind_r[5]+0.3535533905932737*G_1_div_b_Upwind_r[2]*div_b_Surf[3]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_r[3]; 
  Ghat_G_1_div_b_r[7] = 0.3535533905932737*G_1_div_b_Upwind_r[0]*div_b_Surf[7]+0.3535533905932737*div_b_Surf[0]*G_1_div_b_Upwind_r[7]+0.3535533905932737*G_1_div_b_Upwind_r[1]*div_b_Surf[6]+0.3535533905932737*div_b_Surf[1]*G_1_div_b_Upwind_r[6]+0.3535533905932737*G_1_div_b_Upwind_r[2]*div_b_Surf[5]+0.3535533905932737*div_b_Surf[2]*G_1_div_b_Upwind_r[5]+0.3535533905932737*G_1_div_b_Upwind_r[3]*div_b_Surf[4]+0.3535533905932737*div_b_Surf[3]*G_1_div_b_Upwind_r[4]; 

  out_F_0[0] += ((-0.7071067811865475*Ghat_F_0_r[0])+0.7071067811865475*Ghat_F_0_l[0]-0.7071067811865475*Ghat_F_0_div_b_r[0]+0.7071067811865475*Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[1] += ((-0.7071067811865475*Ghat_F_0_r[1])+0.7071067811865475*Ghat_F_0_l[1]-0.7071067811865475*Ghat_F_0_div_b_r[1]+0.7071067811865475*Ghat_F_0_div_b_l[1])*dv1par; 
  out_F_0[2] += ((-0.7071067811865475*Ghat_F_0_r[2])+0.7071067811865475*Ghat_F_0_l[2]-0.7071067811865475*Ghat_F_0_div_b_r[2]+0.7071067811865475*Ghat_F_0_div_b_l[2])*dv1par; 
  out_F_0[3] += ((-0.7071067811865475*Ghat_F_0_r[3])+0.7071067811865475*Ghat_F_0_l[3]-0.7071067811865475*Ghat_F_0_div_b_r[3]+0.7071067811865475*Ghat_F_0_div_b_l[3])*dv1par; 
  out_F_0[4] += -1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0]+Ghat_F_0_div_b_r[0]+Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[5] += ((-0.7071067811865475*Ghat_F_0_r[4])+0.7071067811865475*Ghat_F_0_l[4]-0.7071067811865475*Ghat_F_0_div_b_r[4]+0.7071067811865475*Ghat_F_0_div_b_l[4])*dv1par; 
  out_F_0[6] += ((-0.7071067811865475*Ghat_F_0_r[5])+0.7071067811865475*Ghat_F_0_l[5]-0.7071067811865475*Ghat_F_0_div_b_r[5]+0.7071067811865475*Ghat_F_0_div_b_l[5])*dv1par; 
  out_F_0[7] += ((-0.7071067811865475*Ghat_F_0_r[6])+0.7071067811865475*Ghat_F_0_l[6]-0.7071067811865475*Ghat_F_0_div_b_r[6]+0.7071067811865475*Ghat_F_0_div_b_l[6])*dv1par; 
  out_F_0[8] += -1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1]+Ghat_F_0_div_b_r[1]+Ghat_F_0_div_b_l[1])*dv1par; 
  out_F_0[9] += -1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2]+Ghat_F_0_div_b_r[2]+Ghat_F_0_div_b_l[2])*dv1par; 
  out_F_0[10] += -1.224744871391589*(Ghat_F_0_r[3]+Ghat_F_0_l[3]+Ghat_F_0_div_b_r[3]+Ghat_F_0_div_b_l[3])*dv1par; 
  out_F_0[11] += ((-0.7071067811865475*Ghat_F_0_r[7])+0.7071067811865475*Ghat_F_0_l[7]-0.7071067811865475*Ghat_F_0_div_b_r[7]+0.7071067811865475*Ghat_F_0_div_b_l[7])*dv1par; 
  out_F_0[12] += -1.224744871391589*(Ghat_F_0_r[4]+Ghat_F_0_l[4]+Ghat_F_0_div_b_r[4]+Ghat_F_0_div_b_l[4])*dv1par; 
  out_F_0[13] += -1.224744871391589*(Ghat_F_0_r[5]+Ghat_F_0_l[5]+Ghat_F_0_div_b_r[5]+Ghat_F_0_div_b_l[5])*dv1par; 
  out_F_0[14] += -1.224744871391589*(Ghat_F_0_r[6]+Ghat_F_0_l[6]+Ghat_F_0_div_b_r[6]+Ghat_F_0_div_b_l[6])*dv1par; 
  out_F_0[15] += -1.224744871391589*(Ghat_F_0_r[7]+Ghat_F_0_l[7]+Ghat_F_0_div_b_r[7]+Ghat_F_0_div_b_l[7])*dv1par; 
  out_F_0[16] += ((-1.58113883008419*Ghat_F_0_r[0])+1.58113883008419*Ghat_F_0_l[0]-1.58113883008419*Ghat_F_0_div_b_r[0]+1.58113883008419*Ghat_F_0_div_b_l[0])*dv1par; 
  out_F_0[17] += ((-1.58113883008419*Ghat_F_0_r[1])+1.58113883008419*Ghat_F_0_l[1]-1.58113883008419*Ghat_F_0_div_b_r[1]+1.58113883008419*Ghat_F_0_div_b_l[1])*dv1par; 
  out_F_0[18] += ((-1.58113883008419*Ghat_F_0_r[2])+1.58113883008419*Ghat_F_0_l[2]-1.58113883008419*Ghat_F_0_div_b_r[2]+1.58113883008419*Ghat_F_0_div_b_l[2])*dv1par; 
  out_F_0[19] += ((-1.58113883008419*Ghat_F_0_r[3])+1.58113883008419*Ghat_F_0_l[3]-1.58113883008419*Ghat_F_0_div_b_r[3]+1.58113883008419*Ghat_F_0_div_b_l[3])*dv1par; 
  out_F_0[20] += ((-1.58113883008419*Ghat_F_0_r[4])+1.58113883008419*Ghat_F_0_l[4]-1.58113883008419*Ghat_F_0_div_b_r[4]+1.58113883008419*Ghat_F_0_div_b_l[4])*dv1par; 
  out_F_0[21] += ((-1.58113883008419*Ghat_F_0_r[5])+1.58113883008419*Ghat_F_0_l[5]-1.58113883008419*Ghat_F_0_div_b_r[5]+1.58113883008419*Ghat_F_0_div_b_l[5])*dv1par; 
  out_F_0[22] += ((-1.58113883008419*Ghat_F_0_r[6])+1.58113883008419*Ghat_F_0_l[6]-1.58113883008419*Ghat_F_0_div_b_r[6]+1.58113883008419*Ghat_F_0_div_b_l[6])*dv1par; 
  out_F_0[23] += ((-1.58113883008419*Ghat_F_0_r[7])+1.58113883008419*Ghat_F_0_l[7]-1.58113883008419*Ghat_F_0_div_b_r[7]+1.58113883008419*Ghat_F_0_div_b_l[7])*dv1par; 
  out_G_1[0] += ((-0.7071067811865475*Ghat_G_1_r[0])+0.7071067811865475*Ghat_G_1_l[0]-0.7071067811865475*Ghat_G_1_div_b_r[0]+0.7071067811865475*Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[1] += ((-0.7071067811865475*Ghat_G_1_r[1])+0.7071067811865475*Ghat_G_1_l[1]-0.7071067811865475*Ghat_G_1_div_b_r[1]+0.7071067811865475*Ghat_G_1_div_b_l[1])*dv1par; 
  out_G_1[2] += ((-0.7071067811865475*Ghat_G_1_r[2])+0.7071067811865475*Ghat_G_1_l[2]-0.7071067811865475*Ghat_G_1_div_b_r[2]+0.7071067811865475*Ghat_G_1_div_b_l[2])*dv1par; 
  out_G_1[3] += ((-0.7071067811865475*Ghat_G_1_r[3])+0.7071067811865475*Ghat_G_1_l[3]-0.7071067811865475*Ghat_G_1_div_b_r[3]+0.7071067811865475*Ghat_G_1_div_b_l[3])*dv1par; 
  out_G_1[4] += -1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0]+Ghat_G_1_div_b_r[0]+Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[5] += ((-0.7071067811865475*Ghat_G_1_r[4])+0.7071067811865475*Ghat_G_1_l[4]-0.7071067811865475*Ghat_G_1_div_b_r[4]+0.7071067811865475*Ghat_G_1_div_b_l[4])*dv1par; 
  out_G_1[6] += ((-0.7071067811865475*Ghat_G_1_r[5])+0.7071067811865475*Ghat_G_1_l[5]-0.7071067811865475*Ghat_G_1_div_b_r[5]+0.7071067811865475*Ghat_G_1_div_b_l[5])*dv1par; 
  out_G_1[7] += ((-0.7071067811865475*Ghat_G_1_r[6])+0.7071067811865475*Ghat_G_1_l[6]-0.7071067811865475*Ghat_G_1_div_b_r[6]+0.7071067811865475*Ghat_G_1_div_b_l[6])*dv1par; 
  out_G_1[8] += -1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1]+Ghat_G_1_div_b_r[1]+Ghat_G_1_div_b_l[1])*dv1par; 
  out_G_1[9] += -1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2]+Ghat_G_1_div_b_r[2]+Ghat_G_1_div_b_l[2])*dv1par; 
  out_G_1[10] += -1.224744871391589*(Ghat_G_1_r[3]+Ghat_G_1_l[3]+Ghat_G_1_div_b_r[3]+Ghat_G_1_div_b_l[3])*dv1par; 
  out_G_1[11] += ((-0.7071067811865475*Ghat_G_1_r[7])+0.7071067811865475*Ghat_G_1_l[7]-0.7071067811865475*Ghat_G_1_div_b_r[7]+0.7071067811865475*Ghat_G_1_div_b_l[7])*dv1par; 
  out_G_1[12] += -1.224744871391589*(Ghat_G_1_r[4]+Ghat_G_1_l[4]+Ghat_G_1_div_b_r[4]+Ghat_G_1_div_b_l[4])*dv1par; 
  out_G_1[13] += -1.224744871391589*(Ghat_G_1_r[5]+Ghat_G_1_l[5]+Ghat_G_1_div_b_r[5]+Ghat_G_1_div_b_l[5])*dv1par; 
  out_G_1[14] += -1.224744871391589*(Ghat_G_1_r[6]+Ghat_G_1_l[6]+Ghat_G_1_div_b_r[6]+Ghat_G_1_div_b_l[6])*dv1par; 
  out_G_1[15] += -1.224744871391589*(Ghat_G_1_r[7]+Ghat_G_1_l[7]+Ghat_G_1_div_b_r[7]+Ghat_G_1_div_b_l[7])*dv1par; 
  out_G_1[16] += ((-1.58113883008419*Ghat_G_1_r[0])+1.58113883008419*Ghat_G_1_l[0]-1.58113883008419*Ghat_G_1_div_b_r[0]+1.58113883008419*Ghat_G_1_div_b_l[0])*dv1par; 
  out_G_1[17] += ((-1.58113883008419*Ghat_G_1_r[1])+1.58113883008419*Ghat_G_1_l[1]-1.58113883008419*Ghat_G_1_div_b_r[1]+1.58113883008419*Ghat_G_1_div_b_l[1])*dv1par; 
  out_G_1[18] += ((-1.58113883008419*Ghat_G_1_r[2])+1.58113883008419*Ghat_G_1_l[2]-1.58113883008419*Ghat_G_1_div_b_r[2]+1.58113883008419*Ghat_G_1_div_b_l[2])*dv1par; 
  out_G_1[19] += ((-1.58113883008419*Ghat_G_1_r[3])+1.58113883008419*Ghat_G_1_l[3]-1.58113883008419*Ghat_G_1_div_b_r[3]+1.58113883008419*Ghat_G_1_div_b_l[3])*dv1par; 
  out_G_1[20] += ((-1.58113883008419*Ghat_G_1_r[4])+1.58113883008419*Ghat_G_1_l[4]-1.58113883008419*Ghat_G_1_div_b_r[4]+1.58113883008419*Ghat_G_1_div_b_l[4])*dv1par; 
  out_G_1[21] += ((-1.58113883008419*Ghat_G_1_r[5])+1.58113883008419*Ghat_G_1_l[5]-1.58113883008419*Ghat_G_1_div_b_r[5]+1.58113883008419*Ghat_G_1_div_b_l[5])*dv1par; 
  out_G_1[22] += ((-1.58113883008419*Ghat_G_1_r[6])+1.58113883008419*Ghat_G_1_l[6]-1.58113883008419*Ghat_G_1_div_b_r[6]+1.58113883008419*Ghat_G_1_div_b_l[6])*dv1par; 
  out_G_1[23] += ((-1.58113883008419*Ghat_G_1_r[7])+1.58113883008419*Ghat_G_1_l[7]-1.58113883008419*Ghat_G_1_div_b_r[7]+1.58113883008419*Ghat_G_1_div_b_l[7])*dv1par; 

} 
