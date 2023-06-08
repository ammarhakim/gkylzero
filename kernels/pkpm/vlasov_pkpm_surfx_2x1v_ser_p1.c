#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_hyb_2x1v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_hyb_2x1v_p1_upwind_quad_to_modal.h> 
#include <gkyl_basis_ser_2x_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_2x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_2x1v_ser_p1(const double *w, const double *dxv, 
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
  const double dvpar = dxv[2], wvpar = w[2]; 
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
  const double *G_1l = &fl[12]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[12]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[12]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[12]; 
  double alpha_l[8] = {0.0}; 
  double alpha_c[8] = {0.0}; 
  double alpha_r[8] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar; 
  alpha_l[2] = 1.414213562373095*bl[2]*wvpar; 
  alpha_l[3] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[4] = 1.414213562373095*bl[3]*wvpar; 
  alpha_l[5] = 0.408248290463863*bl[1]*dvpar; 
  alpha_l[6] = 0.408248290463863*bl[2]*dvpar; 
  alpha_l[7] = 0.408248290463863*bl[3]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar; 
  alpha_c[2] = 1.414213562373095*bc[2]*wvpar; 
  alpha_c[3] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[4] = 1.414213562373095*bc[3]*wvpar; 
  alpha_c[5] = 0.408248290463863*bc[1]*dvpar; 
  alpha_c[6] = 0.408248290463863*bc[2]*dvpar; 
  alpha_c[7] = 0.408248290463863*bc[3]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar; 
  alpha_r[2] = 1.414213562373095*br[2]*wvpar; 
  alpha_r[3] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[4] = 1.414213562373095*br[3]*wvpar; 
  alpha_r[5] = 0.408248290463863*br[1]*dvpar; 
  alpha_r[6] = 0.408248290463863*br[2]*dvpar; 
  alpha_r[7] = 0.408248290463863*br[3]*dvpar; 

  double alphaSurf_l[4] = {0.0}; 
  alphaSurf_l[0] = 0.408248290463863*alpha_l[1]-0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.408248290463863*alpha_l[4]-0.408248290463863*alpha_c[4]+0.3535533905932737*alpha_l[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_l[2] = 0.408248290463863*alpha_l[5]-0.408248290463863*alpha_c[5]+0.3535533905932737*alpha_l[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_l[3] = 0.408248290463863*alpha_l[7]-0.408248290463863*alpha_c[7]+0.3535533905932737*alpha_l[6]+0.3535533905932737*alpha_c[6]; 

  double alphaSurf_r[4] = {0.0}; 
  alphaSurf_r[0] = (-0.408248290463863*alpha_r[1])+0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = (-0.408248290463863*alpha_r[4])+0.408248290463863*alpha_c[4]+0.3535533905932737*alpha_r[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_r[2] = (-0.408248290463863*alpha_r[5])+0.408248290463863*alpha_c[5]+0.3535533905932737*alpha_r[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_r[3] = (-0.408248290463863*alpha_r[7])+0.408248290463863*alpha_c[7]+0.3535533905932737*alpha_r[6]+0.3535533905932737*alpha_c[6]; 

  double F_0_UpwindQuad_l[6] = {0.0};
  double F_0_UpwindQuad_r[6] = {0.0};
  double F_0_Upwind_l[6] = {0.0};
  double F_0_Upwind_r[6] = {0.0};
  double Ghat_F_0_l[6] = {0.0}; 
  double Ghat_F_0_r[6] = {0.0}; 
  double G_1_UpwindQuad_l[6] = {0.0};
  double G_1_UpwindQuad_r[6] = {0.0};
  double G_1_Upwind_l[6] = {0.0};
  double G_1_Upwind_r[6] = {0.0};
  double Ghat_G_1_l[6] = {0.0}; 
  double Ghat_G_1_r[6] = {0.0}; 

  if (0.6708203932499369*alphaSurf_l[3]-0.6708203932499369*alphaSurf_l[2]-0.5*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(G_1c); 
  } 
  if (0.6708203932499369*alphaSurf_r[3]-0.6708203932499369*alphaSurf_r[2]-0.5*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = hyb_2x1v_p1_surfx1_eval_quad_node_0_l(G_1r); 
  } 
  if (0.5*alphaSurf_l[0]-0.5*alphaSurf_l[1] > 0) { 
    F_0_UpwindQuad_l[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(G_1c); 
  } 
  if (0.5*alphaSurf_r[0]-0.5*alphaSurf_r[1] > 0) { 
    F_0_UpwindQuad_r[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = hyb_2x1v_p1_surfx1_eval_quad_node_1_l(G_1r); 
  } 
  if ((-0.6708203932499369*alphaSurf_l[3])+0.6708203932499369*alphaSurf_l[2]-0.5*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(G_1c); 
  } 
  if ((-0.6708203932499369*alphaSurf_r[3])+0.6708203932499369*alphaSurf_r[2]-0.5*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = hyb_2x1v_p1_surfx1_eval_quad_node_2_l(G_1r); 
  } 
  if (0.5*(alphaSurf_l[1]+alphaSurf_l[0])-0.6708203932499369*(alphaSurf_l[3]+alphaSurf_l[2]) > 0) { 
    F_0_UpwindQuad_l[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(F_0l); 
    G_1_UpwindQuad_l[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(F_0c); 
    G_1_UpwindQuad_l[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(G_1c); 
  } 
  if (0.5*(alphaSurf_r[1]+alphaSurf_r[0])-0.6708203932499369*(alphaSurf_r[3]+alphaSurf_r[2]) > 0) { 
    F_0_UpwindQuad_r[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(F_0c); 
    G_1_UpwindQuad_r[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(F_0r); 
    G_1_UpwindQuad_r[3] = hyb_2x1v_p1_surfx1_eval_quad_node_3_l(G_1r); 
  } 
  if (0.5*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(F_0l); 
    G_1_UpwindQuad_l[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(F_0c); 
    G_1_UpwindQuad_l[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(G_1c); 
  } 
  if (0.5*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(F_0c); 
    G_1_UpwindQuad_r[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(F_0r); 
    G_1_UpwindQuad_r[4] = hyb_2x1v_p1_surfx1_eval_quad_node_4_l(G_1r); 
  } 
  if (0.6708203932499369*(alphaSurf_l[3]+alphaSurf_l[2])+0.5*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(F_0l); 
    G_1_UpwindQuad_l[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(F_0c); 
    G_1_UpwindQuad_l[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(G_1c); 
  } 
  if (0.6708203932499369*(alphaSurf_r[3]+alphaSurf_r[2])+0.5*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(F_0c); 
    G_1_UpwindQuad_r[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(F_0r); 
    G_1_UpwindQuad_r[5] = hyb_2x1v_p1_surfx1_eval_quad_node_5_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  hyb_2x1v_p1_xdir_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  Ghat_F_0_l[0] = 0.5*F_0_Upwind_l[3]*alphaSurf_l[3]+0.5*F_0_Upwind_l[2]*alphaSurf_l[2]+0.5*F_0_Upwind_l[1]*alphaSurf_l[1]+0.5*F_0_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_F_0_l[1] = 0.5*F_0_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[5]+0.4472135954999579*alphaSurf_l[2]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[3] = 0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[5]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[4] = 0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[5]+0.5*alphaSurf_l[0]*F_0_Upwind_l[4]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_F_0_l[5] = 0.5*alphaSurf_l[0]*F_0_Upwind_l[5]+0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[4]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[3]; 
  Ghat_G_1_l[0] = 0.5*G_1_Upwind_l[3]*alphaSurf_l[3]+0.5*G_1_Upwind_l[2]*alphaSurf_l[2]+0.5*G_1_Upwind_l[1]*alphaSurf_l[1]+0.5*G_1_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_G_1_l[1] = 0.5*G_1_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*G_1_Upwind_l[1]; 
  Ghat_G_1_l[2] = 0.447213595499958*alphaSurf_l[3]*G_1_Upwind_l[5]+0.4472135954999579*alphaSurf_l[2]*G_1_Upwind_l[4]+0.5*G_1_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[3] = 0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[5]+0.4472135954999579*alphaSurf_l[3]*G_1_Upwind_l[4]+0.5*G_1_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[4] = 0.5000000000000001*alphaSurf_l[1]*G_1_Upwind_l[5]+0.5*alphaSurf_l[0]*G_1_Upwind_l[4]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*G_1_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_G_1_l[5] = 0.5*alphaSurf_l[0]*G_1_Upwind_l[5]+0.5000000000000001*alphaSurf_l[1]*G_1_Upwind_l[4]+0.447213595499958*G_1_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[3]; 

  Ghat_F_0_r[0] = 0.5*F_0_Upwind_r[3]*alphaSurf_r[3]+0.5*F_0_Upwind_r[2]*alphaSurf_r[2]+0.5*F_0_Upwind_r[1]*alphaSurf_r[1]+0.5*F_0_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_F_0_r[1] = 0.5*F_0_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[5]+0.4472135954999579*alphaSurf_r[2]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[3] = 0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[5]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[4] = 0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[5]+0.5*alphaSurf_r[0]*F_0_Upwind_r[4]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_F_0_r[5] = 0.5*alphaSurf_r[0]*F_0_Upwind_r[5]+0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[4]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[3]; 
  Ghat_G_1_r[0] = 0.5*G_1_Upwind_r[3]*alphaSurf_r[3]+0.5*G_1_Upwind_r[2]*alphaSurf_r[2]+0.5*G_1_Upwind_r[1]*alphaSurf_r[1]+0.5*G_1_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_G_1_r[1] = 0.5*G_1_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*G_1_Upwind_r[1]; 
  Ghat_G_1_r[2] = 0.447213595499958*alphaSurf_r[3]*G_1_Upwind_r[5]+0.4472135954999579*alphaSurf_r[2]*G_1_Upwind_r[4]+0.5*G_1_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[3] = 0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[5]+0.4472135954999579*alphaSurf_r[3]*G_1_Upwind_r[4]+0.5*G_1_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[4] = 0.5000000000000001*alphaSurf_r[1]*G_1_Upwind_r[5]+0.5*alphaSurf_r[0]*G_1_Upwind_r[4]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*G_1_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_G_1_r[5] = 0.5*alphaSurf_r[0]*G_1_Upwind_r[5]+0.5000000000000001*alphaSurf_r[1]*G_1_Upwind_r[4]+0.447213595499958*G_1_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[3]; 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_l[0]-0.7071067811865475*Ghat_F_0_r[0])*dx1; 
  out_F_0[1] += -1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dx1; 
  out_F_0[2] += (0.7071067811865475*Ghat_F_0_l[1]-0.7071067811865475*Ghat_F_0_r[1])*dx1; 
  out_F_0[3] += (0.7071067811865475*Ghat_F_0_l[2]-0.7071067811865475*Ghat_F_0_r[2])*dx1; 
  out_F_0[4] += -1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dx1; 
  out_F_0[5] += -1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2])*dx1; 
  out_F_0[6] += (0.7071067811865475*Ghat_F_0_l[3]-0.7071067811865475*Ghat_F_0_r[3])*dx1; 
  out_F_0[7] += -1.224744871391589*(Ghat_F_0_r[3]+Ghat_F_0_l[3])*dx1; 
  out_F_0[8] += (0.7071067811865475*Ghat_F_0_l[4]-0.7071067811865475*Ghat_F_0_r[4])*dx1; 
  out_F_0[9] += -1.224744871391589*(Ghat_F_0_r[4]+Ghat_F_0_l[4])*dx1; 
  out_F_0[10] += (0.7071067811865475*Ghat_F_0_l[5]-0.7071067811865475*Ghat_F_0_r[5])*dx1; 
  out_F_0[11] += -1.224744871391589*(Ghat_F_0_r[5]+Ghat_F_0_l[5])*dx1; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_l[0]-0.7071067811865475*Ghat_G_1_r[0])*dx1; 
  out_G_1[1] += -1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dx1; 
  out_G_1[2] += (0.7071067811865475*Ghat_G_1_l[1]-0.7071067811865475*Ghat_G_1_r[1])*dx1; 
  out_G_1[3] += (0.7071067811865475*Ghat_G_1_l[2]-0.7071067811865475*Ghat_G_1_r[2])*dx1; 
  out_G_1[4] += -1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dx1; 
  out_G_1[5] += -1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dx1; 
  out_G_1[6] += (0.7071067811865475*Ghat_G_1_l[3]-0.7071067811865475*Ghat_G_1_r[3])*dx1; 
  out_G_1[7] += -1.224744871391589*(Ghat_G_1_r[3]+Ghat_G_1_l[3])*dx1; 
  out_G_1[8] += (0.7071067811865475*Ghat_G_1_l[4]-0.7071067811865475*Ghat_G_1_r[4])*dx1; 
  out_G_1[9] += -1.224744871391589*(Ghat_G_1_r[4]+Ghat_G_1_l[4])*dx1; 
  out_G_1[10] += (0.7071067811865475*Ghat_G_1_l[5]-0.7071067811865475*Ghat_G_1_r[5])*dx1; 
  out_G_1[11] += -1.224744871391589*(Ghat_G_1_r[5]+Ghat_G_1_l[5])*dx1; 

  double F_0l_rx[2] = {0.0}; 
  double F_0c_lx[2] = {0.0}; 
  double F_0c_rx[2] = {0.0}; 
  double F_0r_lx[2] = {0.0}; 

  double G_1l_rx[2] = {0.0}; 
  double G_1c_lx[2] = {0.0}; 
  double G_1c_rx[2] = {0.0}; 
  double G_1r_lx[2] = {0.0}; 

  double avg_F_0_l[2] = {0.0}; 
  double avg_F_0_r[2] = {0.0}; 
  double avg_G_1_l[2] = {0.0}; 
  double avg_G_1_r[2] = {0.0}; 
  double avg_u_l[2] = {0.0}; 
  double avg_u_r[2] = {0.0}; 
  double ul_r = 0.0; 
  double uc_l = 0.0; 
  double uc_r = 0.0; 
  double ur_l = 0.0; 
  double uQuad_l = 0.0; 
  double uQuad_r = 0.0; 

  double vth_sql_r = 0.0; 
  double vth_sqc_l = 0.0; 
  double vth_sqc_r = 0.0; 
  double vth_sqr_l = 0.0; 
  double vthQuad_l = 0.0; 
  double vthQuad_r = 0.0; 

  double max_speed_quad_l[2] = {0.0}; 
  double max_speed_quad_r[2] = {0.0}; 
  double max_speed_modal_l[2] = {0.0}; 
  double max_speed_modal_r[2] = {0.0}; 
  avg_u_l[0] = 0.6123724356957944*ul[1]-0.6123724356957944*uc[1]+0.3535533905932737*ul[0]+0.3535533905932737*uc[0]; 
  avg_u_l[1] = 0.6123724356957944*ul[3]-0.6123724356957944*uc[3]+0.3535533905932737*ul[2]+0.3535533905932737*uc[2]; 

  avg_u_r[0] = (-0.6123724356957944*ur[1])+0.6123724356957944*uc[1]+0.3535533905932737*ur[0]+0.3535533905932737*uc[0]; 
  avg_u_r[1] = (-0.6123724356957944*ur[3])+0.6123724356957944*uc[3]+0.3535533905932737*ur[2]+0.3535533905932737*uc[2]; 

  ul_r = ser_2x_p1_surfx1_eval_quad_node_0_r(ul); 
  uc_l = ser_2x_p1_surfx1_eval_quad_node_0_l(uc); 
  uc_r = ser_2x_p1_surfx1_eval_quad_node_0_r(uc); 
  ur_l = ser_2x_p1_surfx1_eval_quad_node_0_l(ur); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_2x_p1_surfx1_eval_quad_node_0_r(vth_sql); 
  vth_sqc_l = ser_2x_p1_surfx1_eval_quad_node_0_l(vth_sqc); 
  vth_sqc_r = ser_2x_p1_surfx1_eval_quad_node_0_r(vth_sqc); 
  vth_sqr_l = ser_2x_p1_surfx1_eval_quad_node_0_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[0] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[0] = uQuad_r + vthQuad_r; 
  ul_r = ser_2x_p1_surfx1_eval_quad_node_1_r(ul); 
  uc_l = ser_2x_p1_surfx1_eval_quad_node_1_l(uc); 
  uc_r = ser_2x_p1_surfx1_eval_quad_node_1_r(uc); 
  ur_l = ser_2x_p1_surfx1_eval_quad_node_1_l(ur); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_2x_p1_surfx1_eval_quad_node_1_r(vth_sql); 
  vth_sqc_l = ser_2x_p1_surfx1_eval_quad_node_1_l(vth_sqc); 
  vth_sqc_r = ser_2x_p1_surfx1_eval_quad_node_1_r(vth_sqc); 
  vth_sqr_l = ser_2x_p1_surfx1_eval_quad_node_1_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[1] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[1] = uQuad_r + vthQuad_r; 
  ser_2x_p1_upwind_quad_to_modal(max_speed_quad_l, max_speed_modal_l); 
  ser_2x_p1_upwind_quad_to_modal(max_speed_quad_r, max_speed_modal_r); 
  F_0l_rx[0] = 0.6123724356957944*F_0l[1]+0.3535533905932737*F_0l[0]; 
  F_0l_rx[1] = 0.6123724356957944*F_0l[4]+0.3535533905932737*F_0l[2]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[0]-0.6123724356957944*F_0c[1]; 
  F_0c_lx[1] = 0.3535533905932737*F_0c[2]-0.6123724356957944*F_0c[4]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[1]+0.3535533905932737*F_0c[0]; 
  F_0c_rx[1] = 0.6123724356957944*F_0c[4]+0.3535533905932737*F_0c[2]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[0]-0.6123724356957944*F_0r[1]; 
  F_0r_lx[1] = 0.3535533905932737*F_0r[2]-0.6123724356957944*F_0r[4]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[1]+0.3535533905932737*G_1l[0]; 
  G_1l_rx[1] = 0.6123724356957944*G_1l[4]+0.3535533905932737*G_1l[2]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[0]-0.6123724356957944*G_1c[1]; 
  G_1c_lx[1] = 0.3535533905932737*G_1c[2]-0.6123724356957944*G_1c[4]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[1]+0.3535533905932737*G_1c[0]; 
  G_1c_rx[1] = 0.6123724356957944*G_1c[4]+0.3535533905932737*G_1c[2]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[0]-0.6123724356957944*G_1r[1]; 
  G_1r_lx[1] = 0.3535533905932737*G_1r[2]-0.6123724356957944*G_1r[4]; 

  out_F_0[0] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[0] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[1])-0.6123724356957944*F_0l[0]; 
  F_0l_rx[1] = (-1.060660171779821*F_0l[4])-0.6123724356957944*F_0l[2]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[1]-0.6123724356957944*F_0c[0]; 
  F_0c_lx[1] = 1.060660171779821*F_0c[4]-0.6123724356957944*F_0c[2]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[1]+0.6123724356957944*F_0c[0]; 
  F_0c_rx[1] = 1.060660171779821*F_0c[4]+0.6123724356957944*F_0c[2]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[0]-1.060660171779821*F_0r[1]; 
  F_0r_lx[1] = 0.6123724356957944*F_0r[2]-1.060660171779821*F_0r[4]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[1])-0.6123724356957944*G_1l[0]; 
  G_1l_rx[1] = (-1.060660171779821*G_1l[4])-0.6123724356957944*G_1l[2]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[1]-0.6123724356957944*G_1c[0]; 
  G_1c_lx[1] = 1.060660171779821*G_1c[4]-0.6123724356957944*G_1c[2]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[1]+0.6123724356957944*G_1c[0]; 
  G_1c_rx[1] = 1.060660171779821*G_1c[4]+0.6123724356957944*G_1c[2]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[0]-1.060660171779821*G_1r[1]; 
  G_1r_lx[1] = 0.6123724356957944*G_1r[2]-1.060660171779821*G_1r[4]; 

  out_F_0[1] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[1] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.6123724356957944*F_0l[4]+0.3535533905932737*F_0l[2]; 
  F_0l_rx[1] = 0.6123724356957944*F_0l[1]+0.3535533905932737*F_0l[0]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[2]-0.6123724356957944*F_0c[4]; 
  F_0c_lx[1] = 0.3535533905932737*F_0c[0]-0.6123724356957944*F_0c[1]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[4]+0.3535533905932737*F_0c[2]; 
  F_0c_rx[1] = 0.6123724356957944*F_0c[1]+0.3535533905932737*F_0c[0]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[2]-0.6123724356957944*F_0r[4]; 
  F_0r_lx[1] = 0.3535533905932737*F_0r[0]-0.6123724356957944*F_0r[1]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[4]+0.3535533905932737*G_1l[2]; 
  G_1l_rx[1] = 0.6123724356957944*G_1l[1]+0.3535533905932737*G_1l[0]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[2]-0.6123724356957944*G_1c[4]; 
  G_1c_lx[1] = 0.3535533905932737*G_1c[0]-0.6123724356957944*G_1c[1]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[4]+0.3535533905932737*G_1c[2]; 
  G_1c_rx[1] = 0.6123724356957944*G_1c[1]+0.3535533905932737*G_1c[0]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[2]-0.6123724356957944*G_1r[4]; 
  G_1r_lx[1] = 0.3535533905932737*G_1r[0]-0.6123724356957944*G_1r[1]; 

  out_F_0[2] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[2] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.6123724356957944*F_0l[5]+0.3535533905932737*F_0l[3]; 
  F_0l_rx[1] = 0.6123724356957944*F_0l[7]+0.3535533905932737*F_0l[6]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[3]-0.6123724356957944*F_0c[5]; 
  F_0c_lx[1] = 0.3535533905932737*F_0c[6]-0.6123724356957944*F_0c[7]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[5]+0.3535533905932737*F_0c[3]; 
  F_0c_rx[1] = 0.6123724356957944*F_0c[7]+0.3535533905932737*F_0c[6]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[3]-0.6123724356957944*F_0r[5]; 
  F_0r_lx[1] = 0.3535533905932737*F_0r[6]-0.6123724356957944*F_0r[7]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[5]+0.3535533905932737*G_1l[3]; 
  G_1l_rx[1] = 0.6123724356957944*G_1l[7]+0.3535533905932737*G_1l[6]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[3]-0.6123724356957944*G_1c[5]; 
  G_1c_lx[1] = 0.3535533905932737*G_1c[6]-0.6123724356957944*G_1c[7]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[5]+0.3535533905932737*G_1c[3]; 
  G_1c_rx[1] = 0.6123724356957944*G_1c[7]+0.3535533905932737*G_1c[6]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[3]-0.6123724356957944*G_1r[5]; 
  G_1r_lx[1] = 0.3535533905932737*G_1r[6]-0.6123724356957944*G_1r[7]; 

  out_F_0[3] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[3] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[4])-0.6123724356957944*F_0l[2]; 
  F_0l_rx[1] = (-1.060660171779821*F_0l[1])-0.6123724356957944*F_0l[0]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[4]-0.6123724356957944*F_0c[2]; 
  F_0c_lx[1] = 1.060660171779821*F_0c[1]-0.6123724356957944*F_0c[0]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[4]+0.6123724356957944*F_0c[2]; 
  F_0c_rx[1] = 1.060660171779821*F_0c[1]+0.6123724356957944*F_0c[0]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[2]-1.060660171779821*F_0r[4]; 
  F_0r_lx[1] = 0.6123724356957944*F_0r[0]-1.060660171779821*F_0r[1]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[4])-0.6123724356957944*G_1l[2]; 
  G_1l_rx[1] = (-1.060660171779821*G_1l[1])-0.6123724356957944*G_1l[0]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[4]-0.6123724356957944*G_1c[2]; 
  G_1c_lx[1] = 1.060660171779821*G_1c[1]-0.6123724356957944*G_1c[0]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[4]+0.6123724356957944*G_1c[2]; 
  G_1c_rx[1] = 1.060660171779821*G_1c[1]+0.6123724356957944*G_1c[0]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[2]-1.060660171779821*G_1r[4]; 
  G_1r_lx[1] = 0.6123724356957944*G_1r[0]-1.060660171779821*G_1r[1]; 

  out_F_0[4] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[4] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[5])-0.6123724356957944*F_0l[3]; 
  F_0l_rx[1] = (-1.060660171779821*F_0l[7])-0.6123724356957944*F_0l[6]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[5]-0.6123724356957944*F_0c[3]; 
  F_0c_lx[1] = 1.060660171779821*F_0c[7]-0.6123724356957944*F_0c[6]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[5]+0.6123724356957944*F_0c[3]; 
  F_0c_rx[1] = 1.060660171779821*F_0c[7]+0.6123724356957944*F_0c[6]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[3]-1.060660171779821*F_0r[5]; 
  F_0r_lx[1] = 0.6123724356957944*F_0r[6]-1.060660171779821*F_0r[7]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[5])-0.6123724356957944*G_1l[3]; 
  G_1l_rx[1] = (-1.060660171779821*G_1l[7])-0.6123724356957944*G_1l[6]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[5]-0.6123724356957944*G_1c[3]; 
  G_1c_lx[1] = 1.060660171779821*G_1c[7]-0.6123724356957944*G_1c[6]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[5]+0.6123724356957944*G_1c[3]; 
  G_1c_rx[1] = 1.060660171779821*G_1c[7]+0.6123724356957944*G_1c[6]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[3]-1.060660171779821*G_1r[5]; 
  G_1r_lx[1] = 0.6123724356957944*G_1r[6]-1.060660171779821*G_1r[7]; 

  out_F_0[5] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[5] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.6123724356957944*F_0l[7]+0.3535533905932737*F_0l[6]; 
  F_0l_rx[1] = 0.6123724356957944*F_0l[5]+0.3535533905932737*F_0l[3]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[6]-0.6123724356957944*F_0c[7]; 
  F_0c_lx[1] = 0.3535533905932737*F_0c[3]-0.6123724356957944*F_0c[5]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[7]+0.3535533905932737*F_0c[6]; 
  F_0c_rx[1] = 0.6123724356957944*F_0c[5]+0.3535533905932737*F_0c[3]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[6]-0.6123724356957944*F_0r[7]; 
  F_0r_lx[1] = 0.3535533905932737*F_0r[3]-0.6123724356957944*F_0r[5]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[7]+0.3535533905932737*G_1l[6]; 
  G_1l_rx[1] = 0.6123724356957944*G_1l[5]+0.3535533905932737*G_1l[3]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[6]-0.6123724356957944*G_1c[7]; 
  G_1c_lx[1] = 0.3535533905932737*G_1c[3]-0.6123724356957944*G_1c[5]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[7]+0.3535533905932737*G_1c[6]; 
  G_1c_rx[1] = 0.6123724356957944*G_1c[5]+0.3535533905932737*G_1c[3]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[6]-0.6123724356957944*G_1r[7]; 
  G_1r_lx[1] = 0.3535533905932737*G_1r[3]-0.6123724356957944*G_1r[5]; 

  out_F_0[6] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[6] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[7])-0.6123724356957944*F_0l[6]; 
  F_0l_rx[1] = (-1.060660171779821*F_0l[5])-0.6123724356957944*F_0l[3]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[7]-0.6123724356957944*F_0c[6]; 
  F_0c_lx[1] = 1.060660171779821*F_0c[5]-0.6123724356957944*F_0c[3]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[7]+0.6123724356957944*F_0c[6]; 
  F_0c_rx[1] = 1.060660171779821*F_0c[5]+0.6123724356957944*F_0c[3]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[6]-1.060660171779821*F_0r[7]; 
  F_0r_lx[1] = 0.6123724356957944*F_0r[3]-1.060660171779821*F_0r[5]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[7])-0.6123724356957944*G_1l[6]; 
  G_1l_rx[1] = (-1.060660171779821*G_1l[5])-0.6123724356957944*G_1l[3]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[7]-0.6123724356957944*G_1c[6]; 
  G_1c_lx[1] = 1.060660171779821*G_1c[5]-0.6123724356957944*G_1c[3]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[7]+0.6123724356957944*G_1c[6]; 
  G_1c_rx[1] = 1.060660171779821*G_1c[5]+0.6123724356957944*G_1c[3]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[6]-1.060660171779821*G_1r[7]; 
  G_1r_lx[1] = 0.6123724356957944*G_1r[3]-1.060660171779821*G_1r[5]; 

  out_F_0[7] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[7] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.6123724356957944*F_0l[9]+0.3535533905932737*F_0l[8]; 
  F_0l_rx[1] = 0.6123724356957944*F_0l[11]+0.3535533905932737*F_0l[10]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[8]-0.6123724356957944*F_0c[9]; 
  F_0c_lx[1] = 0.3535533905932737*F_0c[10]-0.6123724356957944*F_0c[11]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[9]+0.3535533905932737*F_0c[8]; 
  F_0c_rx[1] = 0.6123724356957944*F_0c[11]+0.3535533905932737*F_0c[10]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[8]-0.6123724356957944*F_0r[9]; 
  F_0r_lx[1] = 0.3535533905932737*F_0r[10]-0.6123724356957944*F_0r[11]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[9]+0.3535533905932737*G_1l[8]; 
  G_1l_rx[1] = 0.6123724356957944*G_1l[11]+0.3535533905932737*G_1l[10]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[8]-0.6123724356957944*G_1c[9]; 
  G_1c_lx[1] = 0.3535533905932737*G_1c[10]-0.6123724356957944*G_1c[11]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[9]+0.3535533905932737*G_1c[8]; 
  G_1c_rx[1] = 0.6123724356957944*G_1c[11]+0.3535533905932737*G_1c[10]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[8]-0.6123724356957944*G_1r[9]; 
  G_1r_lx[1] = 0.3535533905932737*G_1r[10]-0.6123724356957944*G_1r[11]; 

  out_F_0[8] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[8] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[9])-0.6123724356957944*F_0l[8]; 
  F_0l_rx[1] = (-1.060660171779821*F_0l[11])-0.6123724356957944*F_0l[10]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[9]-0.6123724356957944*F_0c[8]; 
  F_0c_lx[1] = 1.060660171779821*F_0c[11]-0.6123724356957944*F_0c[10]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[9]+0.6123724356957944*F_0c[8]; 
  F_0c_rx[1] = 1.060660171779821*F_0c[11]+0.6123724356957944*F_0c[10]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[8]-1.060660171779821*F_0r[9]; 
  F_0r_lx[1] = 0.6123724356957944*F_0r[10]-1.060660171779821*F_0r[11]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[9])-0.6123724356957944*G_1l[8]; 
  G_1l_rx[1] = (-1.060660171779821*G_1l[11])-0.6123724356957944*G_1l[10]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[9]-0.6123724356957944*G_1c[8]; 
  G_1c_lx[1] = 1.060660171779821*G_1c[11]-0.6123724356957944*G_1c[10]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[9]+0.6123724356957944*G_1c[8]; 
  G_1c_rx[1] = 1.060660171779821*G_1c[11]+0.6123724356957944*G_1c[10]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[8]-1.060660171779821*G_1r[9]; 
  G_1r_lx[1] = 0.6123724356957944*G_1r[10]-1.060660171779821*G_1r[11]; 

  out_F_0[9] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[9] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.6123724356957944*F_0l[11]+0.3535533905932737*F_0l[10]; 
  F_0l_rx[1] = 0.6123724356957944*F_0l[9]+0.3535533905932737*F_0l[8]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[10]-0.6123724356957944*F_0c[11]; 
  F_0c_lx[1] = 0.3535533905932737*F_0c[8]-0.6123724356957944*F_0c[9]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[11]+0.3535533905932737*F_0c[10]; 
  F_0c_rx[1] = 0.6123724356957944*F_0c[9]+0.3535533905932737*F_0c[8]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[10]-0.6123724356957944*F_0r[11]; 
  F_0r_lx[1] = 0.3535533905932737*F_0r[8]-0.6123724356957944*F_0r[9]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[11]+0.3535533905932737*G_1l[10]; 
  G_1l_rx[1] = 0.6123724356957944*G_1l[9]+0.3535533905932737*G_1l[8]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[10]-0.6123724356957944*G_1c[11]; 
  G_1c_lx[1] = 0.3535533905932737*G_1c[8]-0.6123724356957944*G_1c[9]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[11]+0.3535533905932737*G_1c[10]; 
  G_1c_rx[1] = 0.6123724356957944*G_1c[9]+0.3535533905932737*G_1c[8]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[10]-0.6123724356957944*G_1r[11]; 
  G_1r_lx[1] = 0.3535533905932737*G_1r[8]-0.6123724356957944*G_1r[9]; 

  out_F_0[10] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[10] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[11])-0.6123724356957944*F_0l[10]; 
  F_0l_rx[1] = (-1.060660171779821*F_0l[9])-0.6123724356957944*F_0l[8]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[11]-0.6123724356957944*F_0c[10]; 
  F_0c_lx[1] = 1.060660171779821*F_0c[9]-0.6123724356957944*F_0c[8]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[11]+0.6123724356957944*F_0c[10]; 
  F_0c_rx[1] = 1.060660171779821*F_0c[9]+0.6123724356957944*F_0c[8]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[10]-1.060660171779821*F_0r[11]; 
  F_0r_lx[1] = 0.6123724356957944*F_0r[8]-1.060660171779821*F_0r[9]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[11])-0.6123724356957944*G_1l[10]; 
  G_1l_rx[1] = (-1.060660171779821*G_1l[9])-0.6123724356957944*G_1l[8]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[11]-0.6123724356957944*G_1c[10]; 
  G_1c_lx[1] = 1.060660171779821*G_1c[9]-0.6123724356957944*G_1c[8]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[11]+0.6123724356957944*G_1c[10]; 
  G_1c_rx[1] = 1.060660171779821*G_1c[9]+0.6123724356957944*G_1c[8]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[10]-1.060660171779821*G_1r[11]; 
  G_1r_lx[1] = 0.6123724356957944*G_1r[8]-1.060660171779821*G_1r[9]; 

  out_F_0[11] += 0.5*dx1*(F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[11] += 0.5*dx1*(G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

} 
