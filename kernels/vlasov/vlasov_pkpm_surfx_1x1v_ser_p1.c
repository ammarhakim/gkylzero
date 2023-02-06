#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_hyb_1x1v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_hyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_1x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *u_il, const double *u_ic, const double *u_ir, 
     const double *T_ijl, const double *T_ijc, const double *T_ijr,
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                 Cell-center coordinates.
  // dxv[NDIM]:               Cell spacing.
  // bvarl/bvarc/bvarr:       Input magnetic field unit vector in left/center/right cells.
  // u_il/u_ic/u_ir:          Input bulk velocity (ux,uy,uz) in left/center/right cells.
  // T_ijl/T_ijc/T_ijr:       Input Temperature tensor/mass (for penalization) in left/center/right cells.
  // fl/fc/fr:                Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells.
  // out:                     Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *ul = &u_il[0]; 
  const double *uc = &u_ic[0]; 
  const double *ur = &u_ir[0]; 
  const double *bl = &bvarl[0]; 
  const double *bc = &bvarc[0]; 
  const double *br = &bvarr[0]; 
  // Get thermal velocity in direction of update for penalization vth^2 = 3.0*T_ii/m. 
  const double *vth_sql = &T_ijl[0]; 
  const double *vth_sqc = &T_ijc[0]; 
  const double *vth_sqr = &T_ijr[0]; 

  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[6]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[6]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[6]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 
  double alpha_l[6] = {0.0}; 
  double alpha_c[6] = {0.0}; 
  double alpha_r[6] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar; 
  alpha_l[2] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[3] = 0.408248290463863*bl[1]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar; 
  alpha_c[2] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[3] = 0.408248290463863*bc[1]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar; 
  alpha_r[2] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[3] = 0.408248290463863*br[1]*dvpar; 

  double alphaSurf_l[3] = {0.0}; 
  alphaSurf_l[0] = 0.408248290463863*alpha_l[1]-0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.408248290463863*alpha_l[3]-0.408248290463863*alpha_c[3]+0.3535533905932737*alpha_l[2]+0.3535533905932737*alpha_c[2]; 

  double alphaSurf_r[3] = {0.0}; 
  alphaSurf_r[0] = (-0.408248290463863*alpha_r[1])+0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = (-0.408248290463863*alpha_r[3])+0.408248290463863*alpha_c[3]+0.3535533905932737*alpha_r[2]+0.3535533905932737*alpha_c[2]; 

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
    F_0_UpwindQuad_l[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(G_1c); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-0.9486832980505137*alphaSurf_r[1] > 0) { 
    F_0_UpwindQuad_r[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(G_1r); 
  } 
  if (0.7071067811865475*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(G_1c); 
  } 
  if (0.7071067811865475*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(G_1r); 
  } 
  if (0.9486832980505137*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(G_1c); 
  } 
  if (0.9486832980505137*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

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
  out_F_0[4] += (0.7071067811865475*Ghat_F_0_l[2]-0.7071067811865475*Ghat_F_0_r[2])*dx1; 
  out_F_0[5] += -1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2])*dx1; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_l[0]-0.7071067811865475*Ghat_G_1_r[0])*dx1; 
  out_G_1[1] += -1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dx1; 
  out_G_1[2] += (0.7071067811865475*Ghat_G_1_l[1]-0.7071067811865475*Ghat_G_1_r[1])*dx1; 
  out_G_1[3] += -1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dx1; 
  out_G_1[4] += (0.7071067811865475*Ghat_G_1_l[2]-0.7071067811865475*Ghat_G_1_r[2])*dx1; 
  out_G_1[5] += -1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dx1; 

  double alpha_u_l[6] = {0.0}; 
  double alpha_u_c[6] = {0.0}; 
  double alpha_u_r[6] = {0.0}; 
  double vth_sq_l[6] = {0.0}; 
  double vth_sq_c[6] = {0.0}; 
  double vth_sq_r[6] = {0.0}; 
  alpha_u_l[0] = 1.414213562373095*ul[0]; 
  alpha_u_l[1] = 1.414213562373095*ul[1]; 

  alpha_u_c[0] = 1.414213562373095*uc[0]; 
  alpha_u_c[1] = 1.414213562373095*uc[1]; 

  alpha_u_r[0] = 1.414213562373095*ur[0]; 
  alpha_u_r[1] = 1.414213562373095*ur[1]; 

  vth_sq_l[0] = 1.414213562373095*vth_sql[0]; 
  vth_sq_l[1] = 1.414213562373095*vth_sql[1]; 

  vth_sq_c[0] = 1.414213562373095*vth_sqc[0]; 
  vth_sq_c[1] = 1.414213562373095*vth_sqc[1]; 

  vth_sq_r[0] = 1.414213562373095*vth_sqr[0]; 
  vth_sq_r[1] = 1.414213562373095*vth_sqr[1]; 

  double lax_F_0_quad_l[3] = {0.0};
  double lax_F_0_quad_r[3] = {0.0};
  double lax_F_0_modal_l[3] = {0.0};
  double lax_F_0_modal_r[3] = {0.0};
  double lax_G_1_quad_l[3] = {0.0};
  double lax_G_1_quad_r[3] = {0.0};
  double lax_G_1_modal_l[3] = {0.0};
  double lax_G_1_modal_r[3] = {0.0};

  double alpha_l_r = 0.0; 
  double alpha_c_l = 0.0; 
  double alpha_c_r = 0.0; 
  double alpha_r_l = 0.0; 
  double alphaQuad_l = 0.0; 
  double alphaQuad_r = 0.0; 

  double vth_sq_l_r = 0.0; 
  double vth_sq_c_l = 0.0; 
  double vth_sq_c_r = 0.0; 
  double vth_sq_r_l = 0.0; 
  double vthQuad_l = 0.0; 
  double vthQuad_r = 0.0; 

  double max_speed_l = 0.0; 
  double max_speed_r = 0.0; 

  double F_0_l_r = 0.0; 
  double F_0_c_l = 0.0; 
  double F_0_c_r = 0.0; 
  double F_0_r_l = 0.0; 
  double G_1_l_r = 0.0; 
  double G_1_c_l = 0.0; 
  double G_1_c_r = 0.0; 
  double G_1_r_l = 0.0; 

  alpha_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(alpha_u_l); 
  alpha_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(alpha_u_c); 
  alpha_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(alpha_u_c); 
  alpha_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(vth_sq_l); 
  vth_sq_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(vth_sq_c); 
  vth_sq_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(vth_sq_c); 
  vth_sq_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  F_0_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(F_0l); 
  F_0_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(F_0c); 
  F_0_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(F_0c); 
  F_0_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(F_0r); 
  lax_F_0_quad_l[0] = 0.5*(alpha_l_r*F_0_l_r + alpha_c_l*F_0_c_l) - 0.5*max_speed_l*(F_0_c_l - F_0_l_r); 
  lax_F_0_quad_r[0] = 0.5*(alpha_c_r*F_0_c_r + alpha_r_l*F_0_r_l) - 0.5*max_speed_r*(F_0_r_l - F_0_c_r); 
  G_1_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(G_1l); 
  G_1_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(G_1c); 
  G_1_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_0_r(G_1c); 
  G_1_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_0_l(G_1r); 
  lax_G_1_quad_l[0] = 0.5*(alpha_l_r*G_1_l_r + alpha_c_l*G_1_c_l) - 0.5*max_speed_l*(G_1_c_l - G_1_l_r); 
  lax_G_1_quad_r[0] = 0.5*(alpha_c_r*G_1_c_r + alpha_r_l*G_1_r_l) - 0.5*max_speed_r*(G_1_r_l - G_1_c_r); 

  alpha_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(alpha_u_l); 
  alpha_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(alpha_u_c); 
  alpha_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(alpha_u_c); 
  alpha_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(vth_sq_l); 
  vth_sq_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(vth_sq_c); 
  vth_sq_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(vth_sq_c); 
  vth_sq_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  F_0_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(F_0l); 
  F_0_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(F_0c); 
  F_0_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(F_0c); 
  F_0_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(F_0r); 
  lax_F_0_quad_l[1] = 0.5*(alpha_l_r*F_0_l_r + alpha_c_l*F_0_c_l) - 0.5*max_speed_l*(F_0_c_l - F_0_l_r); 
  lax_F_0_quad_r[1] = 0.5*(alpha_c_r*F_0_c_r + alpha_r_l*F_0_r_l) - 0.5*max_speed_r*(F_0_r_l - F_0_c_r); 
  G_1_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(G_1l); 
  G_1_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(G_1c); 
  G_1_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_1_r(G_1c); 
  G_1_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_1_l(G_1r); 
  lax_G_1_quad_l[1] = 0.5*(alpha_l_r*G_1_l_r + alpha_c_l*G_1_c_l) - 0.5*max_speed_l*(G_1_c_l - G_1_l_r); 
  lax_G_1_quad_r[1] = 0.5*(alpha_c_r*G_1_c_r + alpha_r_l*G_1_r_l) - 0.5*max_speed_r*(G_1_r_l - G_1_c_r); 

  alpha_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(alpha_u_l); 
  alpha_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(alpha_u_c); 
  alpha_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(alpha_u_c); 
  alpha_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(alpha_u_r); 
  alphaQuad_l = fmax(fabs(alpha_l_r), fabs(alpha_c_l)); 
  alphaQuad_r = fmax(fabs(alpha_c_r), fabs(alpha_r_l)); 
  vth_sq_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(vth_sq_l); 
  vth_sq_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(vth_sq_c); 
  vth_sq_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(vth_sq_c); 
  vth_sq_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(vth_sq_r); 
  vthQuad_l = fmax(sqrt(fabs(vth_sq_l_r)), sqrt(fabs(vth_sq_c_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sq_c_r)), sqrt(fabs(vth_sq_r_l))); 
  max_speed_l = alphaQuad_l + vthQuad_l; 
  max_speed_r = alphaQuad_r + vthQuad_r; 
  F_0_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(F_0l); 
  F_0_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(F_0c); 
  F_0_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(F_0c); 
  F_0_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(F_0r); 
  lax_F_0_quad_l[2] = 0.5*(alpha_l_r*F_0_l_r + alpha_c_l*F_0_c_l) - 0.5*max_speed_l*(F_0_c_l - F_0_l_r); 
  lax_F_0_quad_r[2] = 0.5*(alpha_c_r*F_0_c_r + alpha_r_l*F_0_r_l) - 0.5*max_speed_r*(F_0_r_l - F_0_c_r); 
  G_1_l_r = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(G_1l); 
  G_1_c_l = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(G_1c); 
  G_1_c_r = hyb_1x1v_p1_surfx1_eval_quad_node_2_r(G_1c); 
  G_1_r_l = hyb_1x1v_p1_surfx1_eval_quad_node_2_l(G_1r); 
  lax_G_1_quad_l[2] = 0.5*(alpha_l_r*G_1_l_r + alpha_c_l*G_1_c_l) - 0.5*max_speed_l*(G_1_c_l - G_1_l_r); 
  lax_G_1_quad_r[2] = 0.5*(alpha_c_r*G_1_c_r + alpha_r_l*G_1_r_l) - 0.5*max_speed_r*(G_1_r_l - G_1_c_r); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(lax_F_0_quad_l, lax_F_0_modal_l); 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(lax_F_0_quad_r, lax_F_0_modal_r); 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(lax_G_1_quad_l, lax_G_1_modal_l); 
  hyb_1x1v_p1_xdir_upwind_quad_to_modal(lax_G_1_quad_r, lax_G_1_modal_r); 

  out_F_0[0] += (0.7071067811865475*lax_F_0_modal_l[0]-0.7071067811865475*lax_F_0_modal_r[0])*dx1; 
  out_F_0[1] += -1.224744871391589*(lax_F_0_modal_r[0]+lax_F_0_modal_l[0])*dx1; 
  out_F_0[2] += (0.7071067811865475*lax_F_0_modal_l[1]-0.7071067811865475*lax_F_0_modal_r[1])*dx1; 
  out_F_0[3] += -1.224744871391589*(lax_F_0_modal_r[1]+lax_F_0_modal_l[1])*dx1; 
  out_F_0[4] += (0.7071067811865475*lax_F_0_modal_l[2]-0.7071067811865475*lax_F_0_modal_r[2])*dx1; 
  out_F_0[5] += -1.224744871391589*(lax_F_0_modal_r[2]+lax_F_0_modal_l[2])*dx1; 
  out_G_1[0] += (0.7071067811865475*lax_G_1_modal_l[0]-0.7071067811865475*lax_G_1_modal_r[0])*dx1; 
  out_G_1[1] += -1.224744871391589*(lax_G_1_modal_r[0]+lax_G_1_modal_l[0])*dx1; 
  out_G_1[2] += (0.7071067811865475*lax_G_1_modal_l[1]-0.7071067811865475*lax_G_1_modal_r[1])*dx1; 
  out_G_1[3] += -1.224744871391589*(lax_G_1_modal_r[1]+lax_G_1_modal_l[1])*dx1; 
  out_G_1[4] += (0.7071067811865475*lax_G_1_modal_l[2]-0.7071067811865475*lax_G_1_modal_r[2])*dx1; 
  out_G_1[5] += -1.224744871391589*(lax_G_1_modal_r[2]+lax_G_1_modal_l[2])*dx1; 


} 
