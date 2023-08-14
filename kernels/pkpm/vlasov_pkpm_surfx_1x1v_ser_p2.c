#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_ser_2x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_pkpm_surfx_1x1v_ser_p2(const double *w, const double *dxv, 
    const double *bvar_l, const double *bvar_c, const double *bvar_r, 
    const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r, 
    const double *fl, const double *fc, const double *fr, 
    const double *pkpm_lax, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:              Cell-center coordinates.
  // dxv[NDIM]:            Cell spacing.
  // bvar_l/c/r:           Input magnetic field unit vector in left/center/right cells.
  // pkpm_prim_surf_l/c/r: Input surface primitive variables [u_i, 3*T_ii/m] in left/center/right cells in each direction.
  // fl/fc/fr:             Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells.
  // pkpm_lax:             Surface expansion of pkpm Lax penalization: lambda_i = |u_i| + sqrt(3.0*T_ii/m).
  // out:                  Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *bl = &bvar_l[0]; 
  const double *bc = &bvar_c[0]; 
  const double *br = &bvar_r[0]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[8]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[8]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[8]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[8]; 
  double alpha_l[8] = {0.0}; 
  double alpha_c[8] = {0.0}; 
  double alpha_r[8] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar; 
  alpha_l[2] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[3] = 0.408248290463863*bl[1]*dvpar; 
  alpha_l[4] = 1.414213562373095*bl[2]*wvpar; 
  alpha_l[6] = 0.408248290463863*bl[2]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar; 
  alpha_c[2] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[3] = 0.408248290463863*bc[1]*dvpar; 
  alpha_c[4] = 1.414213562373095*bc[2]*wvpar; 
  alpha_c[6] = 0.408248290463863*bc[2]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar; 
  alpha_r[2] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[3] = 0.408248290463863*br[1]*dvpar; 
  alpha_r[4] = 1.414213562373095*br[2]*wvpar; 
  alpha_r[6] = 0.408248290463863*br[2]*dvpar; 

  double alphaSurf_l[3] = {0.0}; 
  alphaSurf_l[0] = 0.3458741190809163*alpha_l[4]+0.3458741190809163*alpha_c[4]+0.4975526040028326*alpha_l[1]-0.4975526040028326*alpha_c[1]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.3458741190809163*alpha_l[6]+0.3458741190809163*alpha_c[6]+0.4975526040028326*alpha_l[3]-0.4975526040028326*alpha_c[3]+0.3535533905932737*alpha_l[2]+0.3535533905932737*alpha_c[2]; 

  double alphaSurf_r[3] = {0.0}; 
  alphaSurf_r[0] = 0.3458741190809163*alpha_r[4]+0.3458741190809163*alpha_c[4]-0.4975526040028326*alpha_r[1]+0.4975526040028326*alpha_c[1]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = 0.3458741190809163*alpha_r[6]+0.3458741190809163*alpha_c[6]-0.4975526040028326*alpha_r[3]+0.4975526040028326*alpha_c[3]+0.3535533905932737*alpha_r[2]+0.3535533905932737*alpha_c[2]; 

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

  if (0.7071067811865475*alphaSurf_l[0]-0.9486832980505137*alphaSurf_l[1] > 0) { 
    F_0_UpwindQuad_l[0] = ser_2x_p2_surfx1_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = ser_2x_p2_surfx1_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = ser_2x_p2_surfx1_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = ser_2x_p2_surfx1_eval_quad_node_0_l(G_1c); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-0.9486832980505137*alphaSurf_r[1] > 0) { 
    F_0_UpwindQuad_r[0] = ser_2x_p2_surfx1_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = ser_2x_p2_surfx1_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = ser_2x_p2_surfx1_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = ser_2x_p2_surfx1_eval_quad_node_0_l(G_1r); 
  } 
  if (0.7071067811865475*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = ser_2x_p2_surfx1_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = ser_2x_p2_surfx1_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = ser_2x_p2_surfx1_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = ser_2x_p2_surfx1_eval_quad_node_1_l(G_1c); 
  } 
  if (0.7071067811865475*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = ser_2x_p2_surfx1_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = ser_2x_p2_surfx1_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = ser_2x_p2_surfx1_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = ser_2x_p2_surfx1_eval_quad_node_1_l(G_1r); 
  } 
  if (0.9486832980505137*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = ser_2x_p2_surfx1_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = ser_2x_p2_surfx1_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = ser_2x_p2_surfx1_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = ser_2x_p2_surfx1_eval_quad_node_2_l(G_1c); 
  } 
  if (0.9486832980505137*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = ser_2x_p2_surfx1_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = ser_2x_p2_surfx1_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = ser_2x_p2_surfx1_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = ser_2x_p2_surfx1_eval_quad_node_2_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  ser_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  ser_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  ser_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  Ghat_F_0_l[0] = 0.7071067811865475*F_0_Upwind_l[1]*alphaSurf_l[1]+0.7071067811865475*F_0_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_F_0_l[1] = 0.6324555320336759*alphaSurf_l[1]*F_0_Upwind_l[2]+0.7071067811865475*F_0_Upwind_l[0]*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.7071067811865475*alphaSurf_l[0]*F_0_Upwind_l[2]+0.6324555320336759*F_0_Upwind_l[1]*alphaSurf_l[1]; 
  Ghat_G_1_l[0] = 0.7071067811865475*G_1_Upwind_l[1]*alphaSurf_l[1]+0.7071067811865475*G_1_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_G_1_l[1] = 0.6324555320336759*alphaSurf_l[1]*G_1_Upwind_l[2]+0.7071067811865475*G_1_Upwind_l[0]*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0]*G_1_Upwind_l[1]; 
  Ghat_G_1_l[2] = 0.7071067811865475*alphaSurf_l[0]*G_1_Upwind_l[2]+0.6324555320336759*G_1_Upwind_l[1]*alphaSurf_l[1]; 

  Ghat_F_0_r[0] = 0.7071067811865475*F_0_Upwind_r[1]*alphaSurf_r[1]+0.7071067811865475*F_0_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_F_0_r[1] = 0.6324555320336759*alphaSurf_r[1]*F_0_Upwind_r[2]+0.7071067811865475*F_0_Upwind_r[0]*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.7071067811865475*alphaSurf_r[0]*F_0_Upwind_r[2]+0.6324555320336759*F_0_Upwind_r[1]*alphaSurf_r[1]; 
  Ghat_G_1_r[0] = 0.7071067811865475*G_1_Upwind_r[1]*alphaSurf_r[1]+0.7071067811865475*G_1_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_G_1_r[1] = 0.6324555320336759*alphaSurf_r[1]*G_1_Upwind_r[2]+0.7071067811865475*G_1_Upwind_r[0]*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0]*G_1_Upwind_r[1]; 
  Ghat_G_1_r[2] = 0.7071067811865475*alphaSurf_r[0]*G_1_Upwind_r[2]+0.6324555320336759*G_1_Upwind_r[1]*alphaSurf_r[1]; 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_l[0]-0.7071067811865475*Ghat_F_0_r[0])*dx1; 
  out_F_0[1] += -1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dx1; 
  out_F_0[2] += (0.7071067811865475*Ghat_F_0_l[1]-0.7071067811865475*Ghat_F_0_r[1])*dx1; 
  out_F_0[3] += -1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dx1; 
  out_F_0[4] += (1.58113883008419*Ghat_F_0_l[0]-1.58113883008419*Ghat_F_0_r[0])*dx1; 
  out_F_0[5] += (0.7071067811865475*Ghat_F_0_l[2]-0.7071067811865475*Ghat_F_0_r[2])*dx1; 
  out_F_0[6] += (1.58113883008419*Ghat_F_0_l[1]-1.58113883008419*Ghat_F_0_r[1])*dx1; 
  out_F_0[7] += -1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2])*dx1; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_l[0]-0.7071067811865475*Ghat_G_1_r[0])*dx1; 
  out_G_1[1] += -1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dx1; 
  out_G_1[2] += (0.7071067811865475*Ghat_G_1_l[1]-0.7071067811865475*Ghat_G_1_r[1])*dx1; 
  out_G_1[3] += -1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dx1; 
  out_G_1[4] += (1.58113883008419*Ghat_G_1_l[0]-1.58113883008419*Ghat_G_1_r[0])*dx1; 
  out_G_1[5] += (0.7071067811865475*Ghat_G_1_l[2]-0.7071067811865475*Ghat_G_1_r[2])*dx1; 
  out_G_1[6] += (1.58113883008419*Ghat_G_1_l[1]-1.58113883008419*Ghat_G_1_r[1])*dx1; 
  out_G_1[7] += -1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dx1; 

  const double *u_surf_lr = &pkpm_prim_surf_l[1]; 
  const double *u_surf_cl = &pkpm_prim_surf_c[0]; 
  const double *u_surf_cr = &pkpm_prim_surf_c[1]; 
  const double *u_surf_rl = &pkpm_prim_surf_r[0]; 

  const double *pkpm_lax_l = &pkpm_lax[0]; 
  const double *pkpm_lax_r = &pkpm_lax[1]; 

  double F_0_lr[3] = {0.0}; 
  double F_0_cl[3] = {0.0}; 
  double F_0_cr[3] = {0.0}; 
  double F_0_rl[3] = {0.0}; 
  double Ghat_F_0_u_l[3] = {0.0}; 
  double Ghat_F_0_u_r[3] = {0.0}; 

  double G_1_lr[3] = {0.0}; 
  double G_1_cl[3] = {0.0}; 
  double G_1_cr[3] = {0.0}; 
  double G_1_rl[3] = {0.0}; 
  double Ghat_G_1_u_l[3] = {0.0}; 
  double Ghat_G_1_u_r[3] = {0.0}; 

  F_0_lr[0] = 1.58113883008419*F_0l[4]+1.224744871391589*F_0l[1]+0.7071067811865475*F_0l[0]; 
  F_0_lr[1] = 1.58113883008419*F_0l[6]+1.224744871391589*F_0l[3]+0.7071067811865475*F_0l[2]; 
  F_0_lr[2] = 1.224744871391589*F_0l[7]+0.7071067811865475*F_0l[5]; 

  F_0_cl[0] = 1.58113883008419*F_0c[4]-1.224744871391589*F_0c[1]+0.7071067811865475*F_0c[0]; 
  F_0_cl[1] = 1.58113883008419*F_0c[6]-1.224744871391589*F_0c[3]+0.7071067811865475*F_0c[2]; 
  F_0_cl[2] = 0.7071067811865475*F_0c[5]-1.224744871391589*F_0c[7]; 

  F_0_cr[0] = 1.58113883008419*F_0c[4]+1.224744871391589*F_0c[1]+0.7071067811865475*F_0c[0]; 
  F_0_cr[1] = 1.58113883008419*F_0c[6]+1.224744871391589*F_0c[3]+0.7071067811865475*F_0c[2]; 
  F_0_cr[2] = 1.224744871391589*F_0c[7]+0.7071067811865475*F_0c[5]; 

  F_0_rl[0] = 1.58113883008419*F_0r[4]-1.224744871391589*F_0r[1]+0.7071067811865475*F_0r[0]; 
  F_0_rl[1] = 1.58113883008419*F_0r[6]-1.224744871391589*F_0r[3]+0.7071067811865475*F_0r[2]; 
  F_0_rl[2] = 0.7071067811865475*F_0r[5]-1.224744871391589*F_0r[7]; 

  G_1_lr[0] = 1.58113883008419*G_1l[4]+1.224744871391589*G_1l[1]+0.7071067811865475*G_1l[0]; 
  G_1_lr[1] = 1.58113883008419*G_1l[6]+1.224744871391589*G_1l[3]+0.7071067811865475*G_1l[2]; 
  G_1_lr[2] = 1.224744871391589*G_1l[7]+0.7071067811865475*G_1l[5]; 

  G_1_cl[0] = 1.58113883008419*G_1c[4]-1.224744871391589*G_1c[1]+0.7071067811865475*G_1c[0]; 
  G_1_cl[1] = 1.58113883008419*G_1c[6]-1.224744871391589*G_1c[3]+0.7071067811865475*G_1c[2]; 
  G_1_cl[2] = 0.7071067811865475*G_1c[5]-1.224744871391589*G_1c[7]; 

  G_1_cr[0] = 1.58113883008419*G_1c[4]+1.224744871391589*G_1c[1]+0.7071067811865475*G_1c[0]; 
  G_1_cr[1] = 1.58113883008419*G_1c[6]+1.224744871391589*G_1c[3]+0.7071067811865475*G_1c[2]; 
  G_1_cr[2] = 1.224744871391589*G_1c[7]+0.7071067811865475*G_1c[5]; 

  G_1_rl[0] = 1.58113883008419*G_1r[4]-1.224744871391589*G_1r[1]+0.7071067811865475*G_1r[0]; 
  G_1_rl[1] = 1.58113883008419*G_1r[6]-1.224744871391589*G_1r[3]+0.7071067811865475*G_1r[2]; 
  G_1_rl[2] = 0.7071067811865475*G_1r[5]-1.224744871391589*G_1r[7]; 

  double ul_r = u_surf_lr[0]; 
  double uc_l = u_surf_cl[0]; 
  double uc_r = u_surf_cr[0]; 
  double ur_l = u_surf_rl[0]; 
  double avg_u_l = 0.5*(ul_r + uc_l); 
  double avg_u_r = 0.5*(uc_r + ur_l); 

  double max_speed_l = pkpm_lax_l[0]; 
  double max_speed_r = pkpm_lax_r[0]; 

  Ghat_F_0_u_l[0] = 0.5*F_0_lr[0]*max_speed_l-0.5*F_0_cl[0]*max_speed_l+0.5*F_0_lr[0]*avg_u_l+0.5*F_0_cl[0]*avg_u_l; 
  Ghat_F_0_u_l[1] = 0.5*F_0_lr[1]*max_speed_l-0.5*F_0_cl[1]*max_speed_l+0.5*F_0_lr[1]*avg_u_l+0.5*F_0_cl[1]*avg_u_l; 
  Ghat_F_0_u_l[2] = 0.5*F_0_lr[2]*max_speed_l-0.5*F_0_cl[2]*max_speed_l+0.5*F_0_lr[2]*avg_u_l+0.5*F_0_cl[2]*avg_u_l; 
  Ghat_G_1_u_l[0] = 0.5*G_1_lr[0]*max_speed_l-0.5*G_1_cl[0]*max_speed_l+0.5*G_1_lr[0]*avg_u_l+0.5*G_1_cl[0]*avg_u_l; 
  Ghat_G_1_u_l[1] = 0.5*G_1_lr[1]*max_speed_l-0.5*G_1_cl[1]*max_speed_l+0.5*G_1_lr[1]*avg_u_l+0.5*G_1_cl[1]*avg_u_l; 
  Ghat_G_1_u_l[2] = 0.5*G_1_lr[2]*max_speed_l-0.5*G_1_cl[2]*max_speed_l+0.5*G_1_lr[2]*avg_u_l+0.5*G_1_cl[2]*avg_u_l; 

  Ghat_F_0_u_r[0] = (-0.5*F_0_rl[0]*max_speed_r)+0.5*F_0_cr[0]*max_speed_r+0.5*F_0_rl[0]*avg_u_r+0.5*F_0_cr[0]*avg_u_r; 
  Ghat_F_0_u_r[1] = (-0.5*F_0_rl[1]*max_speed_r)+0.5*F_0_cr[1]*max_speed_r+0.5*F_0_rl[1]*avg_u_r+0.5*F_0_cr[1]*avg_u_r; 
  Ghat_F_0_u_r[2] = (-0.5*F_0_rl[2]*max_speed_r)+0.5*F_0_cr[2]*max_speed_r+0.5*F_0_rl[2]*avg_u_r+0.5*F_0_cr[2]*avg_u_r; 
  Ghat_G_1_u_r[0] = (-0.5*G_1_rl[0]*max_speed_r)+0.5*G_1_cr[0]*max_speed_r+0.5*G_1_rl[0]*avg_u_r+0.5*G_1_cr[0]*avg_u_r; 
  Ghat_G_1_u_r[1] = (-0.5*G_1_rl[1]*max_speed_r)+0.5*G_1_cr[1]*max_speed_r+0.5*G_1_rl[1]*avg_u_r+0.5*G_1_cr[1]*avg_u_r; 
  Ghat_G_1_u_r[2] = (-0.5*G_1_rl[2]*max_speed_r)+0.5*G_1_cr[2]*max_speed_r+0.5*G_1_rl[2]*avg_u_r+0.5*G_1_cr[2]*avg_u_r; 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_u_l[0]-0.7071067811865475*Ghat_F_0_u_r[0])*dx1; 
  out_F_0[1] += -1.224744871391589*(Ghat_F_0_u_r[0]+Ghat_F_0_u_l[0])*dx1; 
  out_F_0[2] += (0.7071067811865475*Ghat_F_0_u_l[1]-0.7071067811865475*Ghat_F_0_u_r[1])*dx1; 
  out_F_0[3] += -1.224744871391589*(Ghat_F_0_u_r[1]+Ghat_F_0_u_l[1])*dx1; 
  out_F_0[4] += (1.58113883008419*Ghat_F_0_u_l[0]-1.58113883008419*Ghat_F_0_u_r[0])*dx1; 
  out_F_0[5] += (0.7071067811865475*Ghat_F_0_u_l[2]-0.7071067811865475*Ghat_F_0_u_r[2])*dx1; 
  out_F_0[6] += (1.58113883008419*Ghat_F_0_u_l[1]-1.58113883008419*Ghat_F_0_u_r[1])*dx1; 
  out_F_0[7] += -1.224744871391589*(Ghat_F_0_u_r[2]+Ghat_F_0_u_l[2])*dx1; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_u_l[0]-0.7071067811865475*Ghat_G_1_u_r[0])*dx1; 
  out_G_1[1] += -1.224744871391589*(Ghat_G_1_u_r[0]+Ghat_G_1_u_l[0])*dx1; 
  out_G_1[2] += (0.7071067811865475*Ghat_G_1_u_l[1]-0.7071067811865475*Ghat_G_1_u_r[1])*dx1; 
  out_G_1[3] += -1.224744871391589*(Ghat_G_1_u_r[1]+Ghat_G_1_u_l[1])*dx1; 
  out_G_1[4] += (1.58113883008419*Ghat_G_1_u_l[0]-1.58113883008419*Ghat_G_1_u_r[0])*dx1; 
  out_G_1[5] += (0.7071067811865475*Ghat_G_1_u_l[2]-0.7071067811865475*Ghat_G_1_u_r[2])*dx1; 
  out_G_1[6] += (1.58113883008419*Ghat_G_1_u_l[1]-1.58113883008419*Ghat_G_1_u_r[1])*dx1; 
  out_G_1[7] += -1.224744871391589*(Ghat_G_1_u_r[2]+Ghat_G_1_u_l[2])*dx1; 

  return 0.;

} 
