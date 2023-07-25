#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
#include <gkyl_basis_ser_2x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_2x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_pkpm_surfy_2x1v_ser_p2(const double *w, const double *dxv, 
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
  const double dx1 = 2.0/dxv[1]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *ul = &u_il[8]; 
  const double *uc = &u_ic[8]; 
  const double *ur = &u_ir[8]; 
  const double *bl = &bvarl[8]; 
  const double *bc = &bvarc[8]; 
  const double *br = &bvarr[8]; 
  // Get thermal velocity in direction of update for penalization vth^2 = 3.0*T_ii/m. 
  const double *vth_sql = &T_ijl[24]; 
  const double *vth_sqc = &T_ijc[24]; 
  const double *vth_sqr = &T_ijr[24]; 

  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[20]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[20]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[20]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[20]; 
  double alpha_l[20] = {0.0}; 
  double alpha_c[20] = {0.0}; 
  double alpha_r[20] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar; 
  alpha_l[2] = 1.414213562373095*bl[2]*wvpar; 
  alpha_l[3] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[4] = 1.414213562373095*bl[3]*wvpar; 
  alpha_l[5] = 0.408248290463863*bl[1]*dvpar; 
  alpha_l[6] = 0.408248290463863*bl[2]*dvpar; 
  alpha_l[7] = 1.414213562373095*bl[4]*wvpar; 
  alpha_l[8] = 1.414213562373095*bl[5]*wvpar; 
  alpha_l[10] = 0.408248290463863*bl[3]*dvpar; 
  alpha_l[11] = 1.414213562373095*bl[6]*wvpar; 
  alpha_l[12] = 1.414213562373095*bl[7]*wvpar; 
  alpha_l[13] = 0.408248290463863*bl[4]*dvpar; 
  alpha_l[14] = 0.408248290463863*bl[5]*dvpar; 
  alpha_l[17] = 0.408248290463863*bl[6]*dvpar; 
  alpha_l[18] = 0.408248290463863*bl[7]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar; 
  alpha_c[2] = 1.414213562373095*bc[2]*wvpar; 
  alpha_c[3] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[4] = 1.414213562373095*bc[3]*wvpar; 
  alpha_c[5] = 0.408248290463863*bc[1]*dvpar; 
  alpha_c[6] = 0.408248290463863*bc[2]*dvpar; 
  alpha_c[7] = 1.414213562373095*bc[4]*wvpar; 
  alpha_c[8] = 1.414213562373095*bc[5]*wvpar; 
  alpha_c[10] = 0.408248290463863*bc[3]*dvpar; 
  alpha_c[11] = 1.414213562373095*bc[6]*wvpar; 
  alpha_c[12] = 1.414213562373095*bc[7]*wvpar; 
  alpha_c[13] = 0.408248290463863*bc[4]*dvpar; 
  alpha_c[14] = 0.408248290463863*bc[5]*dvpar; 
  alpha_c[17] = 0.408248290463863*bc[6]*dvpar; 
  alpha_c[18] = 0.408248290463863*bc[7]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar; 
  alpha_r[2] = 1.414213562373095*br[2]*wvpar; 
  alpha_r[3] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[4] = 1.414213562373095*br[3]*wvpar; 
  alpha_r[5] = 0.408248290463863*br[1]*dvpar; 
  alpha_r[6] = 0.408248290463863*br[2]*dvpar; 
  alpha_r[7] = 1.414213562373095*br[4]*wvpar; 
  alpha_r[8] = 1.414213562373095*br[5]*wvpar; 
  alpha_r[10] = 0.408248290463863*br[3]*dvpar; 
  alpha_r[11] = 1.414213562373095*br[6]*wvpar; 
  alpha_r[12] = 1.414213562373095*br[7]*wvpar; 
  alpha_r[13] = 0.408248290463863*br[4]*dvpar; 
  alpha_r[14] = 0.408248290463863*br[5]*dvpar; 
  alpha_r[17] = 0.408248290463863*br[6]*dvpar; 
  alpha_r[18] = 0.408248290463863*br[7]*dvpar; 

  double alphaSurf_l[8] = {0.0}; 
  alphaSurf_l[0] = 0.3458741190809163*alpha_l[8]+0.3458741190809163*alpha_c[8]+0.4975526040028326*alpha_l[2]-0.4975526040028326*alpha_c[2]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.3458741190809163*alpha_l[12]+0.3458741190809163*alpha_c[12]+0.4975526040028326*alpha_l[4]-0.4975526040028326*alpha_c[4]+0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_c[1]; 
  alphaSurf_l[2] = 0.3458741190809163*alpha_l[14]+0.3458741190809163*alpha_c[14]+0.4975526040028326*alpha_l[6]-0.4975526040028326*alpha_c[6]+0.3535533905932737*alpha_l[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_l[3] = 0.3458741190809163*alpha_l[18]+0.3458741190809163*alpha_c[18]+0.4975526040028326*alpha_l[10]-0.4975526040028326*alpha_c[10]+0.3535533905932737*alpha_l[5]+0.3535533905932737*alpha_c[5]; 
  alphaSurf_l[4] = 0.4975526040028326*alpha_l[11]-0.4975526040028326*alpha_c[11]+0.3535533905932737*alpha_l[7]+0.3535533905932737*alpha_c[7]; 
  alphaSurf_l[6] = 0.4975526040028326*alpha_l[17]-0.4975526040028326*alpha_c[17]+0.3535533905932737*alpha_l[13]+0.3535533905932737*alpha_c[13]; 

  double alphaSurf_r[8] = {0.0}; 
  alphaSurf_r[0] = 0.3458741190809163*alpha_r[8]+0.3458741190809163*alpha_c[8]-0.4975526040028326*alpha_r[2]+0.4975526040028326*alpha_c[2]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = 0.3458741190809163*alpha_r[12]+0.3458741190809163*alpha_c[12]-0.4975526040028326*alpha_r[4]+0.4975526040028326*alpha_c[4]+0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_c[1]; 
  alphaSurf_r[2] = 0.3458741190809163*alpha_r[14]+0.3458741190809163*alpha_c[14]-0.4975526040028326*alpha_r[6]+0.4975526040028326*alpha_c[6]+0.3535533905932737*alpha_r[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_r[3] = 0.3458741190809163*alpha_r[18]+0.3458741190809163*alpha_c[18]-0.4975526040028326*alpha_r[10]+0.4975526040028326*alpha_c[10]+0.3535533905932737*alpha_r[5]+0.3535533905932737*alpha_c[5]; 
  alphaSurf_r[4] = (-0.4975526040028326*alpha_r[11])+0.4975526040028326*alpha_c[11]+0.3535533905932737*alpha_r[7]+0.3535533905932737*alpha_c[7]; 
  alphaSurf_r[6] = (-0.4975526040028326*alpha_r[17])+0.4975526040028326*alpha_c[17]+0.3535533905932737*alpha_r[13]+0.3535533905932737*alpha_c[13]; 

  double F_0_UpwindQuad_l[9] = {0.0};
  double F_0_UpwindQuad_r[9] = {0.0};
  double F_0_Upwind_l[8] = {0.0};
  double F_0_Upwind_r[8] = {0.0};
  double Ghat_F_0_l[8] = {0.0}; 
  double Ghat_F_0_r[8] = {0.0}; 
  double G_1_UpwindQuad_l[9] = {0.0};
  double G_1_UpwindQuad_r[9] = {0.0};
  double G_1_Upwind_l[8] = {0.0};
  double G_1_Upwind_r[8] = {0.0};
  double Ghat_G_1_l[8] = {0.0}; 
  double Ghat_G_1_r[8] = {0.0}; 

  if ((-0.5999999999999999*alphaSurf_l[6])+0.4472135954999579*alphaSurf_l[4]+0.9*alphaSurf_l[3]-0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(G_1c); 
  } 
  if ((-0.5999999999999999*alphaSurf_r[6])+0.4472135954999579*alphaSurf_r[4]+0.9*alphaSurf_r[3]-0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(G_1r); 
  } 
  if (0.4472135954999579*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(G_1c); 
  } 
  if (0.4472135954999579*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(G_1r); 
  } 
  if (0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[4]-0.9*alphaSurf_l[3]+0.6708203932499369*alphaSurf_l[2]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(G_1c); 
  } 
  if (0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[4]-0.9*alphaSurf_r[3]+0.6708203932499369*alphaSurf_r[2]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(G_1r); 
  } 
  if (0.75*alphaSurf_l[6]-0.5590169943749475*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(F_0l); 
    G_1_UpwindQuad_l[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(F_0c); 
    G_1_UpwindQuad_l[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(G_1c); 
  } 
  if (0.75*alphaSurf_r[6]-0.5590169943749475*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(F_0c); 
    G_1_UpwindQuad_r[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(F_0r); 
    G_1_UpwindQuad_r[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(G_1r); 
  } 
  if (0.5*alphaSurf_l[0]-0.5590169943749475*alphaSurf_l[4] > 0) { 
    F_0_UpwindQuad_l[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(F_0l); 
    G_1_UpwindQuad_l[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(F_0c); 
    G_1_UpwindQuad_l[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(G_1c); 
  } 
  if (0.5*alphaSurf_r[0]-0.5590169943749475*alphaSurf_r[4] > 0) { 
    F_0_UpwindQuad_r[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(F_0c); 
    G_1_UpwindQuad_r[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(F_0r); 
    G_1_UpwindQuad_r[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(G_1r); 
  } 
  if ((-0.75*alphaSurf_l[6])-0.5590169943749475*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(F_0l); 
    G_1_UpwindQuad_l[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(F_0c); 
    G_1_UpwindQuad_l[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(G_1c); 
  } 
  if ((-0.75*alphaSurf_r[6])-0.5590169943749475*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(F_0c); 
    G_1_UpwindQuad_r[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(F_0r); 
    G_1_UpwindQuad_r[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(G_1r); 
  } 
  if ((-0.5999999999999999*alphaSurf_l[6])+0.4472135954999579*alphaSurf_l[4]-0.9*alphaSurf_l[3]-0.6708203932499369*alphaSurf_l[2]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(F_0l); 
    G_1_UpwindQuad_l[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(F_0c); 
    G_1_UpwindQuad_l[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(G_1c); 
  } 
  if ((-0.5999999999999999*alphaSurf_r[6])+0.4472135954999579*alphaSurf_r[4]-0.9*alphaSurf_r[3]-0.6708203932499369*alphaSurf_r[2]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(F_0c); 
    G_1_UpwindQuad_r[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(F_0r); 
    G_1_UpwindQuad_r[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(G_1r); 
  } 
  if (0.4472135954999579*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(F_0l); 
    G_1_UpwindQuad_l[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(F_0c); 
    G_1_UpwindQuad_l[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(G_1c); 
  } 
  if (0.4472135954999579*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(F_0c); 
    G_1_UpwindQuad_r[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(F_0r); 
    G_1_UpwindQuad_r[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(G_1r); 
  } 
  if (0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[4]+0.9*alphaSurf_l[3]+0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(F_0l); 
    G_1_UpwindQuad_l[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(F_0c); 
    G_1_UpwindQuad_l[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(G_1c); 
  } 
  if (0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[4]+0.9*alphaSurf_r[3]+0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(F_0c); 
    G_1_UpwindQuad_r[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(F_0r); 
    G_1_UpwindQuad_r[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  ser_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  ser_3x_p2_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  ser_3x_p2_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  Ghat_F_0_l[0] = 0.5*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5*F_0_Upwind_l[4]*alphaSurf_l[4]+0.5*F_0_Upwind_l[3]*alphaSurf_l[3]+0.5*F_0_Upwind_l[2]*alphaSurf_l[2]+0.5*F_0_Upwind_l[1]*alphaSurf_l[1]+0.5*F_0_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_F_0_l[1] = 0.447213595499958*F_0_Upwind_l[3]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[6]+0.4472135954999579*F_0_Upwind_l[1]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[1]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[7]+0.5000000000000001*F_0_Upwind_l[4]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[4]*F_0_Upwind_l[6]+0.4472135954999579*alphaSurf_l[2]*F_0_Upwind_l[5]+0.5*F_0_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[3] = 0.4*alphaSurf_l[6]*F_0_Upwind_l[7]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[7]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[6]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[5]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[4] = 0.31943828249997*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5000000000000001*F_0_Upwind_l[2]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[2]*F_0_Upwind_l[6]+0.31943828249997*F_0_Upwind_l[4]*alphaSurf_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[4]+0.5*alphaSurf_l[0]*F_0_Upwind_l[4]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[1]*alphaSurf_l[1]; 
  Ghat_F_0_l[5] = 0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[7]+0.4472135954999579*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5*alphaSurf_l[0]*F_0_Upwind_l[5]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_F_0_l[6] = 0.4*alphaSurf_l[3]*F_0_Upwind_l[7]+0.4472135954999579*F_0_Upwind_l[5]*alphaSurf_l[6]+0.31943828249997*F_0_Upwind_l[4]*alphaSurf_l[6]+0.5*F_0_Upwind_l[0]*alphaSurf_l[6]+0.31943828249997*alphaSurf_l[4]*F_0_Upwind_l[6]+0.5*alphaSurf_l[0]*F_0_Upwind_l[6]+0.5000000000000001*F_0_Upwind_l[2]*alphaSurf_l[4]+0.5000000000000001*alphaSurf_l[2]*F_0_Upwind_l[4]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[7] = 0.4472135954999579*alphaSurf_l[4]*F_0_Upwind_l[7]+0.5*alphaSurf_l[0]*F_0_Upwind_l[7]+0.4*F_0_Upwind_l[3]*alphaSurf_l[6]+0.4*alphaSurf_l[3]*F_0_Upwind_l[6]+0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[5]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[3]; 
  Ghat_G_1_l[0] = 0.5*G_1_Upwind_l[6]*alphaSurf_l[6]+0.5*G_1_Upwind_l[4]*alphaSurf_l[4]+0.5*G_1_Upwind_l[3]*alphaSurf_l[3]+0.5*G_1_Upwind_l[2]*alphaSurf_l[2]+0.5*G_1_Upwind_l[1]*alphaSurf_l[1]+0.5*G_1_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_G_1_l[1] = 0.447213595499958*G_1_Upwind_l[3]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[3]*G_1_Upwind_l[6]+0.4472135954999579*G_1_Upwind_l[1]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[1]*G_1_Upwind_l[4]+0.5*G_1_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*G_1_Upwind_l[1]; 
  Ghat_G_1_l[2] = 0.447213595499958*alphaSurf_l[3]*G_1_Upwind_l[7]+0.5000000000000001*G_1_Upwind_l[4]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[4]*G_1_Upwind_l[6]+0.4472135954999579*alphaSurf_l[2]*G_1_Upwind_l[5]+0.5*G_1_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[3] = 0.4*alphaSurf_l[6]*G_1_Upwind_l[7]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[7]+0.447213595499958*G_1_Upwind_l[1]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[6]+0.4472135954999579*alphaSurf_l[3]*G_1_Upwind_l[5]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[3]*G_1_Upwind_l[4]+0.5*G_1_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[4] = 0.31943828249997*G_1_Upwind_l[6]*alphaSurf_l[6]+0.5000000000000001*G_1_Upwind_l[2]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[2]*G_1_Upwind_l[6]+0.31943828249997*G_1_Upwind_l[4]*alphaSurf_l[4]+0.5*G_1_Upwind_l[0]*alphaSurf_l[4]+0.5*alphaSurf_l[0]*G_1_Upwind_l[4]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*G_1_Upwind_l[1]*alphaSurf_l[1]; 
  Ghat_G_1_l[5] = 0.5000000000000001*alphaSurf_l[1]*G_1_Upwind_l[7]+0.4472135954999579*G_1_Upwind_l[6]*alphaSurf_l[6]+0.5*alphaSurf_l[0]*G_1_Upwind_l[5]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*G_1_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_G_1_l[6] = 0.4*alphaSurf_l[3]*G_1_Upwind_l[7]+0.4472135954999579*G_1_Upwind_l[5]*alphaSurf_l[6]+0.31943828249997*G_1_Upwind_l[4]*alphaSurf_l[6]+0.5*G_1_Upwind_l[0]*alphaSurf_l[6]+0.31943828249997*alphaSurf_l[4]*G_1_Upwind_l[6]+0.5*alphaSurf_l[0]*G_1_Upwind_l[6]+0.5000000000000001*G_1_Upwind_l[2]*alphaSurf_l[4]+0.5000000000000001*alphaSurf_l[2]*G_1_Upwind_l[4]+0.447213595499958*G_1_Upwind_l[1]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[7] = 0.4472135954999579*alphaSurf_l[4]*G_1_Upwind_l[7]+0.5*alphaSurf_l[0]*G_1_Upwind_l[7]+0.4*G_1_Upwind_l[3]*alphaSurf_l[6]+0.4*alphaSurf_l[3]*G_1_Upwind_l[6]+0.5000000000000001*alphaSurf_l[1]*G_1_Upwind_l[5]+0.447213595499958*G_1_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[3]; 

  Ghat_F_0_r[0] = 0.5*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5*F_0_Upwind_r[4]*alphaSurf_r[4]+0.5*F_0_Upwind_r[3]*alphaSurf_r[3]+0.5*F_0_Upwind_r[2]*alphaSurf_r[2]+0.5*F_0_Upwind_r[1]*alphaSurf_r[1]+0.5*F_0_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_F_0_r[1] = 0.447213595499958*F_0_Upwind_r[3]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[6]+0.4472135954999579*F_0_Upwind_r[1]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[1]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[7]+0.5000000000000001*F_0_Upwind_r[4]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[4]*F_0_Upwind_r[6]+0.4472135954999579*alphaSurf_r[2]*F_0_Upwind_r[5]+0.5*F_0_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[3] = 0.4*alphaSurf_r[6]*F_0_Upwind_r[7]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[7]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[6]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[5]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[4] = 0.31943828249997*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5000000000000001*F_0_Upwind_r[2]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[2]*F_0_Upwind_r[6]+0.31943828249997*F_0_Upwind_r[4]*alphaSurf_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[4]+0.5*alphaSurf_r[0]*F_0_Upwind_r[4]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[1]*alphaSurf_r[1]; 
  Ghat_F_0_r[5] = 0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[7]+0.4472135954999579*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5*alphaSurf_r[0]*F_0_Upwind_r[5]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_F_0_r[6] = 0.4*alphaSurf_r[3]*F_0_Upwind_r[7]+0.4472135954999579*F_0_Upwind_r[5]*alphaSurf_r[6]+0.31943828249997*F_0_Upwind_r[4]*alphaSurf_r[6]+0.5*F_0_Upwind_r[0]*alphaSurf_r[6]+0.31943828249997*alphaSurf_r[4]*F_0_Upwind_r[6]+0.5*alphaSurf_r[0]*F_0_Upwind_r[6]+0.5000000000000001*F_0_Upwind_r[2]*alphaSurf_r[4]+0.5000000000000001*alphaSurf_r[2]*F_0_Upwind_r[4]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[7] = 0.4472135954999579*alphaSurf_r[4]*F_0_Upwind_r[7]+0.5*alphaSurf_r[0]*F_0_Upwind_r[7]+0.4*F_0_Upwind_r[3]*alphaSurf_r[6]+0.4*alphaSurf_r[3]*F_0_Upwind_r[6]+0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[5]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[3]; 
  Ghat_G_1_r[0] = 0.5*G_1_Upwind_r[6]*alphaSurf_r[6]+0.5*G_1_Upwind_r[4]*alphaSurf_r[4]+0.5*G_1_Upwind_r[3]*alphaSurf_r[3]+0.5*G_1_Upwind_r[2]*alphaSurf_r[2]+0.5*G_1_Upwind_r[1]*alphaSurf_r[1]+0.5*G_1_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_G_1_r[1] = 0.447213595499958*G_1_Upwind_r[3]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[3]*G_1_Upwind_r[6]+0.4472135954999579*G_1_Upwind_r[1]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[1]*G_1_Upwind_r[4]+0.5*G_1_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*G_1_Upwind_r[1]; 
  Ghat_G_1_r[2] = 0.447213595499958*alphaSurf_r[3]*G_1_Upwind_r[7]+0.5000000000000001*G_1_Upwind_r[4]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[4]*G_1_Upwind_r[6]+0.4472135954999579*alphaSurf_r[2]*G_1_Upwind_r[5]+0.5*G_1_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[3] = 0.4*alphaSurf_r[6]*G_1_Upwind_r[7]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[7]+0.447213595499958*G_1_Upwind_r[1]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[6]+0.4472135954999579*alphaSurf_r[3]*G_1_Upwind_r[5]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[3]*G_1_Upwind_r[4]+0.5*G_1_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[4] = 0.31943828249997*G_1_Upwind_r[6]*alphaSurf_r[6]+0.5000000000000001*G_1_Upwind_r[2]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[2]*G_1_Upwind_r[6]+0.31943828249997*G_1_Upwind_r[4]*alphaSurf_r[4]+0.5*G_1_Upwind_r[0]*alphaSurf_r[4]+0.5*alphaSurf_r[0]*G_1_Upwind_r[4]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*G_1_Upwind_r[1]*alphaSurf_r[1]; 
  Ghat_G_1_r[5] = 0.5000000000000001*alphaSurf_r[1]*G_1_Upwind_r[7]+0.4472135954999579*G_1_Upwind_r[6]*alphaSurf_r[6]+0.5*alphaSurf_r[0]*G_1_Upwind_r[5]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*G_1_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_G_1_r[6] = 0.4*alphaSurf_r[3]*G_1_Upwind_r[7]+0.4472135954999579*G_1_Upwind_r[5]*alphaSurf_r[6]+0.31943828249997*G_1_Upwind_r[4]*alphaSurf_r[6]+0.5*G_1_Upwind_r[0]*alphaSurf_r[6]+0.31943828249997*alphaSurf_r[4]*G_1_Upwind_r[6]+0.5*alphaSurf_r[0]*G_1_Upwind_r[6]+0.5000000000000001*G_1_Upwind_r[2]*alphaSurf_r[4]+0.5000000000000001*alphaSurf_r[2]*G_1_Upwind_r[4]+0.447213595499958*G_1_Upwind_r[1]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[7] = 0.4472135954999579*alphaSurf_r[4]*G_1_Upwind_r[7]+0.5*alphaSurf_r[0]*G_1_Upwind_r[7]+0.4*G_1_Upwind_r[3]*alphaSurf_r[6]+0.4*alphaSurf_r[3]*G_1_Upwind_r[6]+0.5000000000000001*alphaSurf_r[1]*G_1_Upwind_r[5]+0.447213595499958*G_1_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[3]; 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_l[0]-0.7071067811865475*Ghat_F_0_r[0])*dx1; 
  out_F_0[1] += (0.7071067811865475*Ghat_F_0_l[1]-0.7071067811865475*Ghat_F_0_r[1])*dx1; 
  out_F_0[2] += -1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dx1; 
  out_F_0[3] += (0.7071067811865475*Ghat_F_0_l[2]-0.7071067811865475*Ghat_F_0_r[2])*dx1; 
  out_F_0[4] += -1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dx1; 
  out_F_0[5] += (0.7071067811865475*Ghat_F_0_l[3]-0.7071067811865475*Ghat_F_0_r[3])*dx1; 
  out_F_0[6] += -1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2])*dx1; 
  out_F_0[7] += (0.7071067811865475*Ghat_F_0_l[4]-0.7071067811865475*Ghat_F_0_r[4])*dx1; 
  out_F_0[8] += (1.58113883008419*Ghat_F_0_l[0]-1.58113883008419*Ghat_F_0_r[0])*dx1; 
  out_F_0[9] += (0.7071067811865475*Ghat_F_0_l[5]-0.7071067811865475*Ghat_F_0_r[5])*dx1; 
  out_F_0[10] += -1.224744871391589*(Ghat_F_0_r[3]+Ghat_F_0_l[3])*dx1; 
  out_F_0[11] += -1.224744871391589*(Ghat_F_0_r[4]+Ghat_F_0_l[4])*dx1; 
  out_F_0[12] += (1.58113883008419*Ghat_F_0_l[1]-1.58113883008419*Ghat_F_0_r[1])*dx1; 
  out_F_0[13] += (0.7071067811865475*Ghat_F_0_l[6]-0.7071067811865475*Ghat_F_0_r[6])*dx1; 
  out_F_0[14] += (1.58113883008419*Ghat_F_0_l[2]-1.58113883008419*Ghat_F_0_r[2])*dx1; 
  out_F_0[15] += (0.7071067811865475*Ghat_F_0_l[7]-0.7071067811865475*Ghat_F_0_r[7])*dx1; 
  out_F_0[16] += -1.224744871391589*(Ghat_F_0_r[5]+Ghat_F_0_l[5])*dx1; 
  out_F_0[17] += -1.224744871391589*(Ghat_F_0_r[6]+Ghat_F_0_l[6])*dx1; 
  out_F_0[18] += (1.58113883008419*Ghat_F_0_l[3]-1.58113883008419*Ghat_F_0_r[3])*dx1; 
  out_F_0[19] += -1.224744871391589*(Ghat_F_0_r[7]+Ghat_F_0_l[7])*dx1; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_l[0]-0.7071067811865475*Ghat_G_1_r[0])*dx1; 
  out_G_1[1] += (0.7071067811865475*Ghat_G_1_l[1]-0.7071067811865475*Ghat_G_1_r[1])*dx1; 
  out_G_1[2] += -1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dx1; 
  out_G_1[3] += (0.7071067811865475*Ghat_G_1_l[2]-0.7071067811865475*Ghat_G_1_r[2])*dx1; 
  out_G_1[4] += -1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dx1; 
  out_G_1[5] += (0.7071067811865475*Ghat_G_1_l[3]-0.7071067811865475*Ghat_G_1_r[3])*dx1; 
  out_G_1[6] += -1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dx1; 
  out_G_1[7] += (0.7071067811865475*Ghat_G_1_l[4]-0.7071067811865475*Ghat_G_1_r[4])*dx1; 
  out_G_1[8] += (1.58113883008419*Ghat_G_1_l[0]-1.58113883008419*Ghat_G_1_r[0])*dx1; 
  out_G_1[9] += (0.7071067811865475*Ghat_G_1_l[5]-0.7071067811865475*Ghat_G_1_r[5])*dx1; 
  out_G_1[10] += -1.224744871391589*(Ghat_G_1_r[3]+Ghat_G_1_l[3])*dx1; 
  out_G_1[11] += -1.224744871391589*(Ghat_G_1_r[4]+Ghat_G_1_l[4])*dx1; 
  out_G_1[12] += (1.58113883008419*Ghat_G_1_l[1]-1.58113883008419*Ghat_G_1_r[1])*dx1; 
  out_G_1[13] += (0.7071067811865475*Ghat_G_1_l[6]-0.7071067811865475*Ghat_G_1_r[6])*dx1; 
  out_G_1[14] += (1.58113883008419*Ghat_G_1_l[2]-1.58113883008419*Ghat_G_1_r[2])*dx1; 
  out_G_1[15] += (0.7071067811865475*Ghat_G_1_l[7]-0.7071067811865475*Ghat_G_1_r[7])*dx1; 
  out_G_1[16] += -1.224744871391589*(Ghat_G_1_r[5]+Ghat_G_1_l[5])*dx1; 
  out_G_1[17] += -1.224744871391589*(Ghat_G_1_r[6]+Ghat_G_1_l[6])*dx1; 
  out_G_1[18] += (1.58113883008419*Ghat_G_1_l[3]-1.58113883008419*Ghat_G_1_r[3])*dx1; 
  out_G_1[19] += -1.224744871391589*(Ghat_G_1_r[7]+Ghat_G_1_l[7])*dx1; 

  double F_0l_rx[3] = {0.0}; 
  double F_0c_lx[3] = {0.0}; 
  double F_0c_rx[3] = {0.0}; 
  double F_0r_lx[3] = {0.0}; 

  double G_1l_rx[3] = {0.0}; 
  double G_1c_lx[3] = {0.0}; 
  double G_1c_rx[3] = {0.0}; 
  double G_1r_lx[3] = {0.0}; 

  double avg_F_0_l[3] = {0.0}; 
  double avg_F_0_r[3] = {0.0}; 
  double avg_G_1_l[3] = {0.0}; 
  double avg_G_1_r[3] = {0.0}; 
  double avg_u_l[3] = {0.0}; 
  double avg_u_r[3] = {0.0}; 
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

  double max_speed_quad_l[3] = {0.0}; 
  double max_speed_quad_r[3] = {0.0}; 
  double max_speed_modal_l[3] = {0.0}; 
  double max_speed_modal_r[3] = {0.0}; 
  avg_u_l[0] = 0.7905694150420947*ul[5]+0.7905694150420947*uc[5]+0.6123724356957944*ul[2]-0.6123724356957944*uc[2]+0.3535533905932737*ul[0]+0.3535533905932737*uc[0]; 
  avg_u_l[1] = 0.7905694150420948*ul[7]+0.7905694150420948*uc[7]+0.6123724356957944*ul[3]-0.6123724356957944*uc[3]+0.3535533905932737*ul[1]+0.3535533905932737*uc[1]; 
  avg_u_l[2] = 0.6123724356957944*ul[6]-0.6123724356957944*uc[6]+0.3535533905932737*ul[4]+0.3535533905932737*uc[4]; 

  avg_u_r[0] = 0.7905694150420947*ur[5]+0.7905694150420947*uc[5]-0.6123724356957944*ur[2]+0.6123724356957944*uc[2]+0.3535533905932737*ur[0]+0.3535533905932737*uc[0]; 
  avg_u_r[1] = 0.7905694150420948*ur[7]+0.7905694150420948*uc[7]-0.6123724356957944*ur[3]+0.6123724356957944*uc[3]+0.3535533905932737*ur[1]+0.3535533905932737*uc[1]; 
  avg_u_r[2] = (-0.6123724356957944*ur[6])+0.6123724356957944*uc[6]+0.3535533905932737*ur[4]+0.3535533905932737*uc[4]; 

  ul_r = ser_2x_p2_surfx2_eval_quad_node_0_r(ul); 
  uc_l = ser_2x_p2_surfx2_eval_quad_node_0_l(uc); 
  uc_r = ser_2x_p2_surfx2_eval_quad_node_0_r(uc); 
  ur_l = ser_2x_p2_surfx2_eval_quad_node_0_l(ur); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_2x_p2_surfx2_eval_quad_node_0_r(vth_sql); 
  vth_sqc_l = ser_2x_p2_surfx2_eval_quad_node_0_l(vth_sqc); 
  vth_sqc_r = ser_2x_p2_surfx2_eval_quad_node_0_r(vth_sqc); 
  vth_sqr_l = ser_2x_p2_surfx2_eval_quad_node_0_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[0] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[0] = uQuad_r + vthQuad_r; 
  ul_r = ser_2x_p2_surfx2_eval_quad_node_1_r(ul); 
  uc_l = ser_2x_p2_surfx2_eval_quad_node_1_l(uc); 
  uc_r = ser_2x_p2_surfx2_eval_quad_node_1_r(uc); 
  ur_l = ser_2x_p2_surfx2_eval_quad_node_1_l(ur); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_2x_p2_surfx2_eval_quad_node_1_r(vth_sql); 
  vth_sqc_l = ser_2x_p2_surfx2_eval_quad_node_1_l(vth_sqc); 
  vth_sqc_r = ser_2x_p2_surfx2_eval_quad_node_1_r(vth_sqc); 
  vth_sqr_l = ser_2x_p2_surfx2_eval_quad_node_1_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[1] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[1] = uQuad_r + vthQuad_r; 
  ul_r = ser_2x_p2_surfx2_eval_quad_node_2_r(ul); 
  uc_l = ser_2x_p2_surfx2_eval_quad_node_2_l(uc); 
  uc_r = ser_2x_p2_surfx2_eval_quad_node_2_r(uc); 
  ur_l = ser_2x_p2_surfx2_eval_quad_node_2_l(ur); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_2x_p2_surfx2_eval_quad_node_2_r(vth_sql); 
  vth_sqc_l = ser_2x_p2_surfx2_eval_quad_node_2_l(vth_sqc); 
  vth_sqc_r = ser_2x_p2_surfx2_eval_quad_node_2_r(vth_sqc); 
  vth_sqr_l = ser_2x_p2_surfx2_eval_quad_node_2_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[2] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[2] = uQuad_r + vthQuad_r; 
  ser_2x_p2_upwind_quad_to_modal(max_speed_quad_l, max_speed_modal_l); 
  ser_2x_p2_upwind_quad_to_modal(max_speed_quad_r, max_speed_modal_r); 
  F_0l_rx[0] = 0.7905694150420947*F_0l[8]+0.6123724356957944*F_0l[2]+0.3535533905932737*F_0l[0]; 
  F_0l_rx[1] = 0.7905694150420948*F_0l[12]+0.6123724356957944*F_0l[4]+0.3535533905932737*F_0l[1]; 
  F_0l_rx[2] = 0.6123724356957944*F_0l[11]+0.3535533905932737*F_0l[7]; 

  F_0c_lx[0] = 0.7905694150420947*F_0c[8]-0.6123724356957944*F_0c[2]+0.3535533905932737*F_0c[0]; 
  F_0c_lx[1] = 0.7905694150420948*F_0c[12]-0.6123724356957944*F_0c[4]+0.3535533905932737*F_0c[1]; 
  F_0c_lx[2] = 0.3535533905932737*F_0c[7]-0.6123724356957944*F_0c[11]; 

  F_0c_rx[0] = 0.7905694150420947*F_0c[8]+0.6123724356957944*F_0c[2]+0.3535533905932737*F_0c[0]; 
  F_0c_rx[1] = 0.7905694150420948*F_0c[12]+0.6123724356957944*F_0c[4]+0.3535533905932737*F_0c[1]; 
  F_0c_rx[2] = 0.6123724356957944*F_0c[11]+0.3535533905932737*F_0c[7]; 

  F_0r_lx[0] = 0.7905694150420947*F_0r[8]-0.6123724356957944*F_0r[2]+0.3535533905932737*F_0r[0]; 
  F_0r_lx[1] = 0.7905694150420948*F_0r[12]-0.6123724356957944*F_0r[4]+0.3535533905932737*F_0r[1]; 
  F_0r_lx[2] = 0.3535533905932737*F_0r[7]-0.6123724356957944*F_0r[11]; 

  G_1l_rx[0] = 0.7905694150420947*G_1l[8]+0.6123724356957944*G_1l[2]+0.3535533905932737*G_1l[0]; 
  G_1l_rx[1] = 0.7905694150420948*G_1l[12]+0.6123724356957944*G_1l[4]+0.3535533905932737*G_1l[1]; 
  G_1l_rx[2] = 0.6123724356957944*G_1l[11]+0.3535533905932737*G_1l[7]; 

  G_1c_lx[0] = 0.7905694150420947*G_1c[8]-0.6123724356957944*G_1c[2]+0.3535533905932737*G_1c[0]; 
  G_1c_lx[1] = 0.7905694150420948*G_1c[12]-0.6123724356957944*G_1c[4]+0.3535533905932737*G_1c[1]; 
  G_1c_lx[2] = 0.3535533905932737*G_1c[7]-0.6123724356957944*G_1c[11]; 

  G_1c_rx[0] = 0.7905694150420947*G_1c[8]+0.6123724356957944*G_1c[2]+0.3535533905932737*G_1c[0]; 
  G_1c_rx[1] = 0.7905694150420948*G_1c[12]+0.6123724356957944*G_1c[4]+0.3535533905932737*G_1c[1]; 
  G_1c_rx[2] = 0.6123724356957944*G_1c[11]+0.3535533905932737*G_1c[7]; 

  G_1r_lx[0] = 0.7905694150420947*G_1r[8]-0.6123724356957944*G_1r[2]+0.3535533905932737*G_1r[0]; 
  G_1r_lx[1] = 0.7905694150420948*G_1r[12]-0.6123724356957944*G_1r[4]+0.3535533905932737*G_1r[1]; 
  G_1r_lx[2] = 0.3535533905932737*G_1r[7]-0.6123724356957944*G_1r[11]; 

  out_F_0[0] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[0] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.7905694150420948*F_0l[12]+0.6123724356957944*F_0l[4]+0.3535533905932737*F_0l[1]; 
  F_0l_rx[1] = 0.5477225575051661*F_0l[11]+0.7905694150420947*F_0l[8]+0.3162277660168379*F_0l[7]+0.6123724356957944*F_0l[2]+0.3535533905932737*F_0l[0]; 
  F_0l_rx[2] = 0.7071067811865475*F_0l[12]+0.5477225575051661*F_0l[4]+0.3162277660168379*F_0l[1]; 

  F_0c_lx[0] = 0.7905694150420948*F_0c[12]-0.6123724356957944*F_0c[4]+0.3535533905932737*F_0c[1]; 
  F_0c_lx[1] = (-0.5477225575051661*F_0c[11])+0.7905694150420947*F_0c[8]+0.3162277660168379*F_0c[7]-0.6123724356957944*F_0c[2]+0.3535533905932737*F_0c[0]; 
  F_0c_lx[2] = 0.7071067811865475*F_0c[12]-0.5477225575051661*F_0c[4]+0.3162277660168379*F_0c[1]; 

  F_0c_rx[0] = 0.7905694150420948*F_0c[12]+0.6123724356957944*F_0c[4]+0.3535533905932737*F_0c[1]; 
  F_0c_rx[1] = 0.5477225575051661*F_0c[11]+0.7905694150420947*F_0c[8]+0.3162277660168379*F_0c[7]+0.6123724356957944*F_0c[2]+0.3535533905932737*F_0c[0]; 
  F_0c_rx[2] = 0.7071067811865475*F_0c[12]+0.5477225575051661*F_0c[4]+0.3162277660168379*F_0c[1]; 

  F_0r_lx[0] = 0.7905694150420948*F_0r[12]-0.6123724356957944*F_0r[4]+0.3535533905932737*F_0r[1]; 
  F_0r_lx[1] = (-0.5477225575051661*F_0r[11])+0.7905694150420947*F_0r[8]+0.3162277660168379*F_0r[7]-0.6123724356957944*F_0r[2]+0.3535533905932737*F_0r[0]; 
  F_0r_lx[2] = 0.7071067811865475*F_0r[12]-0.5477225575051661*F_0r[4]+0.3162277660168379*F_0r[1]; 

  G_1l_rx[0] = 0.7905694150420948*G_1l[12]+0.6123724356957944*G_1l[4]+0.3535533905932737*G_1l[1]; 
  G_1l_rx[1] = 0.5477225575051661*G_1l[11]+0.7905694150420947*G_1l[8]+0.3162277660168379*G_1l[7]+0.6123724356957944*G_1l[2]+0.3535533905932737*G_1l[0]; 
  G_1l_rx[2] = 0.7071067811865475*G_1l[12]+0.5477225575051661*G_1l[4]+0.3162277660168379*G_1l[1]; 

  G_1c_lx[0] = 0.7905694150420948*G_1c[12]-0.6123724356957944*G_1c[4]+0.3535533905932737*G_1c[1]; 
  G_1c_lx[1] = (-0.5477225575051661*G_1c[11])+0.7905694150420947*G_1c[8]+0.3162277660168379*G_1c[7]-0.6123724356957944*G_1c[2]+0.3535533905932737*G_1c[0]; 
  G_1c_lx[2] = 0.7071067811865475*G_1c[12]-0.5477225575051661*G_1c[4]+0.3162277660168379*G_1c[1]; 

  G_1c_rx[0] = 0.7905694150420948*G_1c[12]+0.6123724356957944*G_1c[4]+0.3535533905932737*G_1c[1]; 
  G_1c_rx[1] = 0.5477225575051661*G_1c[11]+0.7905694150420947*G_1c[8]+0.3162277660168379*G_1c[7]+0.6123724356957944*G_1c[2]+0.3535533905932737*G_1c[0]; 
  G_1c_rx[2] = 0.7071067811865475*G_1c[12]+0.5477225575051661*G_1c[4]+0.3162277660168379*G_1c[1]; 

  G_1r_lx[0] = 0.7905694150420948*G_1r[12]-0.6123724356957944*G_1r[4]+0.3535533905932737*G_1r[1]; 
  G_1r_lx[1] = (-0.5477225575051661*G_1r[11])+0.7905694150420947*G_1r[8]+0.3162277660168379*G_1r[7]-0.6123724356957944*G_1r[2]+0.3535533905932737*G_1r[0]; 
  G_1r_lx[2] = 0.7071067811865475*G_1r[12]-0.5477225575051661*G_1r[4]+0.3162277660168379*G_1r[1]; 

  out_F_0[1] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[1] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.369306393762915*F_0l[8])-1.060660171779821*F_0l[2]-0.6123724356957944*F_0l[0]; 
  F_0l_rx[1] = (-1.369306393762915*F_0l[12])-1.060660171779821*F_0l[4]-0.6123724356957944*F_0l[1]; 
  F_0l_rx[2] = (-1.060660171779821*F_0l[11])-0.6123724356957944*F_0l[7]; 

  F_0c_lx[0] = (-1.369306393762915*F_0c[8])+1.060660171779821*F_0c[2]-0.6123724356957944*F_0c[0]; 
  F_0c_lx[1] = (-1.369306393762915*F_0c[12])+1.060660171779821*F_0c[4]-0.6123724356957944*F_0c[1]; 
  F_0c_lx[2] = 1.060660171779821*F_0c[11]-0.6123724356957944*F_0c[7]; 

  F_0c_rx[0] = 1.369306393762915*F_0c[8]+1.060660171779821*F_0c[2]+0.6123724356957944*F_0c[0]; 
  F_0c_rx[1] = 1.369306393762915*F_0c[12]+1.060660171779821*F_0c[4]+0.6123724356957944*F_0c[1]; 
  F_0c_rx[2] = 1.060660171779821*F_0c[11]+0.6123724356957944*F_0c[7]; 

  F_0r_lx[0] = 1.369306393762915*F_0r[8]-1.060660171779821*F_0r[2]+0.6123724356957944*F_0r[0]; 
  F_0r_lx[1] = 1.369306393762915*F_0r[12]-1.060660171779821*F_0r[4]+0.6123724356957944*F_0r[1]; 
  F_0r_lx[2] = 0.6123724356957944*F_0r[7]-1.060660171779821*F_0r[11]; 

  G_1l_rx[0] = (-1.369306393762915*G_1l[8])-1.060660171779821*G_1l[2]-0.6123724356957944*G_1l[0]; 
  G_1l_rx[1] = (-1.369306393762915*G_1l[12])-1.060660171779821*G_1l[4]-0.6123724356957944*G_1l[1]; 
  G_1l_rx[2] = (-1.060660171779821*G_1l[11])-0.6123724356957944*G_1l[7]; 

  G_1c_lx[0] = (-1.369306393762915*G_1c[8])+1.060660171779821*G_1c[2]-0.6123724356957944*G_1c[0]; 
  G_1c_lx[1] = (-1.369306393762915*G_1c[12])+1.060660171779821*G_1c[4]-0.6123724356957944*G_1c[1]; 
  G_1c_lx[2] = 1.060660171779821*G_1c[11]-0.6123724356957944*G_1c[7]; 

  G_1c_rx[0] = 1.369306393762915*G_1c[8]+1.060660171779821*G_1c[2]+0.6123724356957944*G_1c[0]; 
  G_1c_rx[1] = 1.369306393762915*G_1c[12]+1.060660171779821*G_1c[4]+0.6123724356957944*G_1c[1]; 
  G_1c_rx[2] = 1.060660171779821*G_1c[11]+0.6123724356957944*G_1c[7]; 

  G_1r_lx[0] = 1.369306393762915*G_1r[8]-1.060660171779821*G_1r[2]+0.6123724356957944*G_1r[0]; 
  G_1r_lx[1] = 1.369306393762915*G_1r[12]-1.060660171779821*G_1r[4]+0.6123724356957944*G_1r[1]; 
  G_1r_lx[2] = 0.6123724356957944*G_1r[7]-1.060660171779821*G_1r[11]; 

  out_F_0[2] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[2] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.7905694150420948*F_0l[14]+0.6123724356957944*F_0l[6]+0.3535533905932737*F_0l[3]; 
  F_0l_rx[1] = 0.7905694150420947*F_0l[18]+0.6123724356957944*F_0l[10]+0.3535533905932737*F_0l[5]; 
  F_0l_rx[2] = 0.6123724356957944*F_0l[17]+0.3535533905932737*F_0l[13]; 

  F_0c_lx[0] = 0.7905694150420948*F_0c[14]-0.6123724356957944*F_0c[6]+0.3535533905932737*F_0c[3]; 
  F_0c_lx[1] = 0.7905694150420947*F_0c[18]-0.6123724356957944*F_0c[10]+0.3535533905932737*F_0c[5]; 
  F_0c_lx[2] = 0.3535533905932737*F_0c[13]-0.6123724356957944*F_0c[17]; 

  F_0c_rx[0] = 0.7905694150420948*F_0c[14]+0.6123724356957944*F_0c[6]+0.3535533905932737*F_0c[3]; 
  F_0c_rx[1] = 0.7905694150420947*F_0c[18]+0.6123724356957944*F_0c[10]+0.3535533905932737*F_0c[5]; 
  F_0c_rx[2] = 0.6123724356957944*F_0c[17]+0.3535533905932737*F_0c[13]; 

  F_0r_lx[0] = 0.7905694150420948*F_0r[14]-0.6123724356957944*F_0r[6]+0.3535533905932737*F_0r[3]; 
  F_0r_lx[1] = 0.7905694150420947*F_0r[18]-0.6123724356957944*F_0r[10]+0.3535533905932737*F_0r[5]; 
  F_0r_lx[2] = 0.3535533905932737*F_0r[13]-0.6123724356957944*F_0r[17]; 

  G_1l_rx[0] = 0.7905694150420948*G_1l[14]+0.6123724356957944*G_1l[6]+0.3535533905932737*G_1l[3]; 
  G_1l_rx[1] = 0.7905694150420947*G_1l[18]+0.6123724356957944*G_1l[10]+0.3535533905932737*G_1l[5]; 
  G_1l_rx[2] = 0.6123724356957944*G_1l[17]+0.3535533905932737*G_1l[13]; 

  G_1c_lx[0] = 0.7905694150420948*G_1c[14]-0.6123724356957944*G_1c[6]+0.3535533905932737*G_1c[3]; 
  G_1c_lx[1] = 0.7905694150420947*G_1c[18]-0.6123724356957944*G_1c[10]+0.3535533905932737*G_1c[5]; 
  G_1c_lx[2] = 0.3535533905932737*G_1c[13]-0.6123724356957944*G_1c[17]; 

  G_1c_rx[0] = 0.7905694150420948*G_1c[14]+0.6123724356957944*G_1c[6]+0.3535533905932737*G_1c[3]; 
  G_1c_rx[1] = 0.7905694150420947*G_1c[18]+0.6123724356957944*G_1c[10]+0.3535533905932737*G_1c[5]; 
  G_1c_rx[2] = 0.6123724356957944*G_1c[17]+0.3535533905932737*G_1c[13]; 

  G_1r_lx[0] = 0.7905694150420948*G_1r[14]-0.6123724356957944*G_1r[6]+0.3535533905932737*G_1r[3]; 
  G_1r_lx[1] = 0.7905694150420947*G_1r[18]-0.6123724356957944*G_1r[10]+0.3535533905932737*G_1r[5]; 
  G_1r_lx[2] = 0.3535533905932737*G_1r[13]-0.6123724356957944*G_1r[17]; 

  out_F_0[3] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[3] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.369306393762915*F_0l[12])-1.060660171779821*F_0l[4]-0.6123724356957944*F_0l[1]; 
  F_0l_rx[1] = (-0.9486832980505138*F_0l[11])-1.369306393762915*F_0l[8]-0.5477225575051661*F_0l[7]-1.060660171779821*F_0l[2]-0.6123724356957944*F_0l[0]; 
  F_0l_rx[2] = (-1.224744871391589*F_0l[12])-0.9486832980505137*F_0l[4]-0.5477225575051661*F_0l[1]; 

  F_0c_lx[0] = (-1.369306393762915*F_0c[12])+1.060660171779821*F_0c[4]-0.6123724356957944*F_0c[1]; 
  F_0c_lx[1] = 0.9486832980505138*F_0c[11]-1.369306393762915*F_0c[8]-0.5477225575051661*F_0c[7]+1.060660171779821*F_0c[2]-0.6123724356957944*F_0c[0]; 
  F_0c_lx[2] = (-1.224744871391589*F_0c[12])+0.9486832980505137*F_0c[4]-0.5477225575051661*F_0c[1]; 

  F_0c_rx[0] = 1.369306393762915*F_0c[12]+1.060660171779821*F_0c[4]+0.6123724356957944*F_0c[1]; 
  F_0c_rx[1] = 0.9486832980505138*F_0c[11]+1.369306393762915*F_0c[8]+0.5477225575051661*F_0c[7]+1.060660171779821*F_0c[2]+0.6123724356957944*F_0c[0]; 
  F_0c_rx[2] = 1.224744871391589*F_0c[12]+0.9486832980505137*F_0c[4]+0.5477225575051661*F_0c[1]; 

  F_0r_lx[0] = 1.369306393762915*F_0r[12]-1.060660171779821*F_0r[4]+0.6123724356957944*F_0r[1]; 
  F_0r_lx[1] = (-0.9486832980505138*F_0r[11])+1.369306393762915*F_0r[8]+0.5477225575051661*F_0r[7]-1.060660171779821*F_0r[2]+0.6123724356957944*F_0r[0]; 
  F_0r_lx[2] = 1.224744871391589*F_0r[12]-0.9486832980505137*F_0r[4]+0.5477225575051661*F_0r[1]; 

  G_1l_rx[0] = (-1.369306393762915*G_1l[12])-1.060660171779821*G_1l[4]-0.6123724356957944*G_1l[1]; 
  G_1l_rx[1] = (-0.9486832980505138*G_1l[11])-1.369306393762915*G_1l[8]-0.5477225575051661*G_1l[7]-1.060660171779821*G_1l[2]-0.6123724356957944*G_1l[0]; 
  G_1l_rx[2] = (-1.224744871391589*G_1l[12])-0.9486832980505137*G_1l[4]-0.5477225575051661*G_1l[1]; 

  G_1c_lx[0] = (-1.369306393762915*G_1c[12])+1.060660171779821*G_1c[4]-0.6123724356957944*G_1c[1]; 
  G_1c_lx[1] = 0.9486832980505138*G_1c[11]-1.369306393762915*G_1c[8]-0.5477225575051661*G_1c[7]+1.060660171779821*G_1c[2]-0.6123724356957944*G_1c[0]; 
  G_1c_lx[2] = (-1.224744871391589*G_1c[12])+0.9486832980505137*G_1c[4]-0.5477225575051661*G_1c[1]; 

  G_1c_rx[0] = 1.369306393762915*G_1c[12]+1.060660171779821*G_1c[4]+0.6123724356957944*G_1c[1]; 
  G_1c_rx[1] = 0.9486832980505138*G_1c[11]+1.369306393762915*G_1c[8]+0.5477225575051661*G_1c[7]+1.060660171779821*G_1c[2]+0.6123724356957944*G_1c[0]; 
  G_1c_rx[2] = 1.224744871391589*G_1c[12]+0.9486832980505137*G_1c[4]+0.5477225575051661*G_1c[1]; 

  G_1r_lx[0] = 1.369306393762915*G_1r[12]-1.060660171779821*G_1r[4]+0.6123724356957944*G_1r[1]; 
  G_1r_lx[1] = (-0.9486832980505138*G_1r[11])+1.369306393762915*G_1r[8]+0.5477225575051661*G_1r[7]-1.060660171779821*G_1r[2]+0.6123724356957944*G_1r[0]; 
  G_1r_lx[2] = 1.224744871391589*G_1r[12]-0.9486832980505137*G_1r[4]+0.5477225575051661*G_1r[1]; 

  out_F_0[4] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[4] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.7905694150420947*F_0l[18]+0.6123724356957944*F_0l[10]+0.3535533905932737*F_0l[5]; 
  F_0l_rx[1] = 0.5477225575051661*F_0l[17]+0.7905694150420948*F_0l[14]+0.3162277660168379*F_0l[13]+0.6123724356957944*F_0l[6]+0.3535533905932737*F_0l[3]; 
  F_0l_rx[2] = 0.7071067811865475*F_0l[18]+0.5477225575051661*F_0l[10]+0.3162277660168379*F_0l[5]; 

  F_0c_lx[0] = 0.7905694150420947*F_0c[18]-0.6123724356957944*F_0c[10]+0.3535533905932737*F_0c[5]; 
  F_0c_lx[1] = (-0.5477225575051661*F_0c[17])+0.7905694150420948*F_0c[14]+0.3162277660168379*F_0c[13]-0.6123724356957944*F_0c[6]+0.3535533905932737*F_0c[3]; 
  F_0c_lx[2] = 0.7071067811865475*F_0c[18]-0.5477225575051661*F_0c[10]+0.3162277660168379*F_0c[5]; 

  F_0c_rx[0] = 0.7905694150420947*F_0c[18]+0.6123724356957944*F_0c[10]+0.3535533905932737*F_0c[5]; 
  F_0c_rx[1] = 0.5477225575051661*F_0c[17]+0.7905694150420948*F_0c[14]+0.3162277660168379*F_0c[13]+0.6123724356957944*F_0c[6]+0.3535533905932737*F_0c[3]; 
  F_0c_rx[2] = 0.7071067811865475*F_0c[18]+0.5477225575051661*F_0c[10]+0.3162277660168379*F_0c[5]; 

  F_0r_lx[0] = 0.7905694150420947*F_0r[18]-0.6123724356957944*F_0r[10]+0.3535533905932737*F_0r[5]; 
  F_0r_lx[1] = (-0.5477225575051661*F_0r[17])+0.7905694150420948*F_0r[14]+0.3162277660168379*F_0r[13]-0.6123724356957944*F_0r[6]+0.3535533905932737*F_0r[3]; 
  F_0r_lx[2] = 0.7071067811865475*F_0r[18]-0.5477225575051661*F_0r[10]+0.3162277660168379*F_0r[5]; 

  G_1l_rx[0] = 0.7905694150420947*G_1l[18]+0.6123724356957944*G_1l[10]+0.3535533905932737*G_1l[5]; 
  G_1l_rx[1] = 0.5477225575051661*G_1l[17]+0.7905694150420948*G_1l[14]+0.3162277660168379*G_1l[13]+0.6123724356957944*G_1l[6]+0.3535533905932737*G_1l[3]; 
  G_1l_rx[2] = 0.7071067811865475*G_1l[18]+0.5477225575051661*G_1l[10]+0.3162277660168379*G_1l[5]; 

  G_1c_lx[0] = 0.7905694150420947*G_1c[18]-0.6123724356957944*G_1c[10]+0.3535533905932737*G_1c[5]; 
  G_1c_lx[1] = (-0.5477225575051661*G_1c[17])+0.7905694150420948*G_1c[14]+0.3162277660168379*G_1c[13]-0.6123724356957944*G_1c[6]+0.3535533905932737*G_1c[3]; 
  G_1c_lx[2] = 0.7071067811865475*G_1c[18]-0.5477225575051661*G_1c[10]+0.3162277660168379*G_1c[5]; 

  G_1c_rx[0] = 0.7905694150420947*G_1c[18]+0.6123724356957944*G_1c[10]+0.3535533905932737*G_1c[5]; 
  G_1c_rx[1] = 0.5477225575051661*G_1c[17]+0.7905694150420948*G_1c[14]+0.3162277660168379*G_1c[13]+0.6123724356957944*G_1c[6]+0.3535533905932737*G_1c[3]; 
  G_1c_rx[2] = 0.7071067811865475*G_1c[18]+0.5477225575051661*G_1c[10]+0.3162277660168379*G_1c[5]; 

  G_1r_lx[0] = 0.7905694150420947*G_1r[18]-0.6123724356957944*G_1r[10]+0.3535533905932737*G_1r[5]; 
  G_1r_lx[1] = (-0.5477225575051661*G_1r[17])+0.7905694150420948*G_1r[14]+0.3162277660168379*G_1r[13]-0.6123724356957944*G_1r[6]+0.3535533905932737*G_1r[3]; 
  G_1r_lx[2] = 0.7071067811865475*G_1r[18]-0.5477225575051661*G_1r[10]+0.3162277660168379*G_1r[5]; 

  out_F_0[5] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[5] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.369306393762915*F_0l[14])-1.060660171779821*F_0l[6]-0.6123724356957944*F_0l[3]; 
  F_0l_rx[1] = (-1.369306393762915*F_0l[18])-1.060660171779821*F_0l[10]-0.6123724356957944*F_0l[5]; 
  F_0l_rx[2] = (-1.060660171779821*F_0l[17])-0.6123724356957944*F_0l[13]; 

  F_0c_lx[0] = (-1.369306393762915*F_0c[14])+1.060660171779821*F_0c[6]-0.6123724356957944*F_0c[3]; 
  F_0c_lx[1] = (-1.369306393762915*F_0c[18])+1.060660171779821*F_0c[10]-0.6123724356957944*F_0c[5]; 
  F_0c_lx[2] = 1.060660171779821*F_0c[17]-0.6123724356957944*F_0c[13]; 

  F_0c_rx[0] = 1.369306393762915*F_0c[14]+1.060660171779821*F_0c[6]+0.6123724356957944*F_0c[3]; 
  F_0c_rx[1] = 1.369306393762915*F_0c[18]+1.060660171779821*F_0c[10]+0.6123724356957944*F_0c[5]; 
  F_0c_rx[2] = 1.060660171779821*F_0c[17]+0.6123724356957944*F_0c[13]; 

  F_0r_lx[0] = 1.369306393762915*F_0r[14]-1.060660171779821*F_0r[6]+0.6123724356957944*F_0r[3]; 
  F_0r_lx[1] = 1.369306393762915*F_0r[18]-1.060660171779821*F_0r[10]+0.6123724356957944*F_0r[5]; 
  F_0r_lx[2] = 0.6123724356957944*F_0r[13]-1.060660171779821*F_0r[17]; 

  G_1l_rx[0] = (-1.369306393762915*G_1l[14])-1.060660171779821*G_1l[6]-0.6123724356957944*G_1l[3]; 
  G_1l_rx[1] = (-1.369306393762915*G_1l[18])-1.060660171779821*G_1l[10]-0.6123724356957944*G_1l[5]; 
  G_1l_rx[2] = (-1.060660171779821*G_1l[17])-0.6123724356957944*G_1l[13]; 

  G_1c_lx[0] = (-1.369306393762915*G_1c[14])+1.060660171779821*G_1c[6]-0.6123724356957944*G_1c[3]; 
  G_1c_lx[1] = (-1.369306393762915*G_1c[18])+1.060660171779821*G_1c[10]-0.6123724356957944*G_1c[5]; 
  G_1c_lx[2] = 1.060660171779821*G_1c[17]-0.6123724356957944*G_1c[13]; 

  G_1c_rx[0] = 1.369306393762915*G_1c[14]+1.060660171779821*G_1c[6]+0.6123724356957944*G_1c[3]; 
  G_1c_rx[1] = 1.369306393762915*G_1c[18]+1.060660171779821*G_1c[10]+0.6123724356957944*G_1c[5]; 
  G_1c_rx[2] = 1.060660171779821*G_1c[17]+0.6123724356957944*G_1c[13]; 

  G_1r_lx[0] = 1.369306393762915*G_1r[14]-1.060660171779821*G_1r[6]+0.6123724356957944*G_1r[3]; 
  G_1r_lx[1] = 1.369306393762915*G_1r[18]-1.060660171779821*G_1r[10]+0.6123724356957944*G_1r[5]; 
  G_1r_lx[2] = 0.6123724356957944*G_1r[13]-1.060660171779821*G_1r[17]; 

  out_F_0[6] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[6] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.6123724356957944*F_0l[11]+0.3535533905932737*F_0l[7]; 
  F_0l_rx[1] = 0.7071067811865475*F_0l[12]+0.5477225575051661*F_0l[4]+0.3162277660168379*F_0l[1]; 
  F_0l_rx[2] = 0.3912303982179757*F_0l[11]+0.7905694150420947*F_0l[8]+0.2258769757263128*F_0l[7]+0.6123724356957944*F_0l[2]+0.3535533905932737*F_0l[0]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[7]-0.6123724356957944*F_0c[11]; 
  F_0c_lx[1] = 0.7071067811865475*F_0c[12]-0.5477225575051661*F_0c[4]+0.3162277660168379*F_0c[1]; 
  F_0c_lx[2] = (-0.3912303982179757*F_0c[11])+0.7905694150420947*F_0c[8]+0.2258769757263128*F_0c[7]-0.6123724356957944*F_0c[2]+0.3535533905932737*F_0c[0]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[11]+0.3535533905932737*F_0c[7]; 
  F_0c_rx[1] = 0.7071067811865475*F_0c[12]+0.5477225575051661*F_0c[4]+0.3162277660168379*F_0c[1]; 
  F_0c_rx[2] = 0.3912303982179757*F_0c[11]+0.7905694150420947*F_0c[8]+0.2258769757263128*F_0c[7]+0.6123724356957944*F_0c[2]+0.3535533905932737*F_0c[0]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[7]-0.6123724356957944*F_0r[11]; 
  F_0r_lx[1] = 0.7071067811865475*F_0r[12]-0.5477225575051661*F_0r[4]+0.3162277660168379*F_0r[1]; 
  F_0r_lx[2] = (-0.3912303982179757*F_0r[11])+0.7905694150420947*F_0r[8]+0.2258769757263128*F_0r[7]-0.6123724356957944*F_0r[2]+0.3535533905932737*F_0r[0]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[11]+0.3535533905932737*G_1l[7]; 
  G_1l_rx[1] = 0.7071067811865475*G_1l[12]+0.5477225575051661*G_1l[4]+0.3162277660168379*G_1l[1]; 
  G_1l_rx[2] = 0.3912303982179757*G_1l[11]+0.7905694150420947*G_1l[8]+0.2258769757263128*G_1l[7]+0.6123724356957944*G_1l[2]+0.3535533905932737*G_1l[0]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[7]-0.6123724356957944*G_1c[11]; 
  G_1c_lx[1] = 0.7071067811865475*G_1c[12]-0.5477225575051661*G_1c[4]+0.3162277660168379*G_1c[1]; 
  G_1c_lx[2] = (-0.3912303982179757*G_1c[11])+0.7905694150420947*G_1c[8]+0.2258769757263128*G_1c[7]-0.6123724356957944*G_1c[2]+0.3535533905932737*G_1c[0]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[11]+0.3535533905932737*G_1c[7]; 
  G_1c_rx[1] = 0.7071067811865475*G_1c[12]+0.5477225575051661*G_1c[4]+0.3162277660168379*G_1c[1]; 
  G_1c_rx[2] = 0.3912303982179757*G_1c[11]+0.7905694150420947*G_1c[8]+0.2258769757263128*G_1c[7]+0.6123724356957944*G_1c[2]+0.3535533905932737*G_1c[0]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[7]-0.6123724356957944*G_1r[11]; 
  G_1r_lx[1] = 0.7071067811865475*G_1r[12]-0.5477225575051661*G_1r[4]+0.3162277660168379*G_1r[1]; 
  G_1r_lx[2] = (-0.3912303982179757*G_1r[11])+0.7905694150420947*G_1r[8]+0.2258769757263128*G_1r[7]-0.6123724356957944*G_1r[2]+0.3535533905932737*G_1r[0]; 

  out_F_0[7] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[7] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 1.767766952966368*F_0l[8]+1.369306393762915*F_0l[2]+0.7905694150420947*F_0l[0]; 
  F_0l_rx[1] = 1.767766952966369*F_0l[12]+1.369306393762915*F_0l[4]+0.7905694150420947*F_0l[1]; 
  F_0l_rx[2] = 1.369306393762915*F_0l[11]+0.7905694150420947*F_0l[7]; 

  F_0c_lx[0] = 1.767766952966368*F_0c[8]-1.369306393762915*F_0c[2]+0.7905694150420947*F_0c[0]; 
  F_0c_lx[1] = 1.767766952966369*F_0c[12]-1.369306393762915*F_0c[4]+0.7905694150420947*F_0c[1]; 
  F_0c_lx[2] = 0.7905694150420947*F_0c[7]-1.369306393762915*F_0c[11]; 

  F_0c_rx[0] = 1.767766952966368*F_0c[8]+1.369306393762915*F_0c[2]+0.7905694150420947*F_0c[0]; 
  F_0c_rx[1] = 1.767766952966369*F_0c[12]+1.369306393762915*F_0c[4]+0.7905694150420947*F_0c[1]; 
  F_0c_rx[2] = 1.369306393762915*F_0c[11]+0.7905694150420947*F_0c[7]; 

  F_0r_lx[0] = 1.767766952966368*F_0r[8]-1.369306393762915*F_0r[2]+0.7905694150420947*F_0r[0]; 
  F_0r_lx[1] = 1.767766952966369*F_0r[12]-1.369306393762915*F_0r[4]+0.7905694150420947*F_0r[1]; 
  F_0r_lx[2] = 0.7905694150420947*F_0r[7]-1.369306393762915*F_0r[11]; 

  G_1l_rx[0] = 1.767766952966368*G_1l[8]+1.369306393762915*G_1l[2]+0.7905694150420947*G_1l[0]; 
  G_1l_rx[1] = 1.767766952966369*G_1l[12]+1.369306393762915*G_1l[4]+0.7905694150420947*G_1l[1]; 
  G_1l_rx[2] = 1.369306393762915*G_1l[11]+0.7905694150420947*G_1l[7]; 

  G_1c_lx[0] = 1.767766952966368*G_1c[8]-1.369306393762915*G_1c[2]+0.7905694150420947*G_1c[0]; 
  G_1c_lx[1] = 1.767766952966369*G_1c[12]-1.369306393762915*G_1c[4]+0.7905694150420947*G_1c[1]; 
  G_1c_lx[2] = 0.7905694150420947*G_1c[7]-1.369306393762915*G_1c[11]; 

  G_1c_rx[0] = 1.767766952966368*G_1c[8]+1.369306393762915*G_1c[2]+0.7905694150420947*G_1c[0]; 
  G_1c_rx[1] = 1.767766952966369*G_1c[12]+1.369306393762915*G_1c[4]+0.7905694150420947*G_1c[1]; 
  G_1c_rx[2] = 1.369306393762915*G_1c[11]+0.7905694150420947*G_1c[7]; 

  G_1r_lx[0] = 1.767766952966368*G_1r[8]-1.369306393762915*G_1r[2]+0.7905694150420947*G_1r[0]; 
  G_1r_lx[1] = 1.767766952966369*G_1r[12]-1.369306393762915*G_1r[4]+0.7905694150420947*G_1r[1]; 
  G_1r_lx[2] = 0.7905694150420947*G_1r[7]-1.369306393762915*G_1r[11]; 

  out_F_0[8] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[8] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.6123724356957944*F_0l[16]+0.3535533905932737*F_0l[9]; 
  F_0l_rx[1] = 0.6123724356957944*F_0l[19]+0.3535533905932737*F_0l[15]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[9]-0.6123724356957944*F_0c[16]; 
  F_0c_lx[1] = 0.3535533905932737*F_0c[15]-0.6123724356957944*F_0c[19]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[16]+0.3535533905932737*F_0c[9]; 
  F_0c_rx[1] = 0.6123724356957944*F_0c[19]+0.3535533905932737*F_0c[15]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[9]-0.6123724356957944*F_0r[16]; 
  F_0r_lx[1] = 0.3535533905932737*F_0r[15]-0.6123724356957944*F_0r[19]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[16]+0.3535533905932737*G_1l[9]; 
  G_1l_rx[1] = 0.6123724356957944*G_1l[19]+0.3535533905932737*G_1l[15]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[9]-0.6123724356957944*G_1c[16]; 
  G_1c_lx[1] = 0.3535533905932737*G_1c[15]-0.6123724356957944*G_1c[19]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[16]+0.3535533905932737*G_1c[9]; 
  G_1c_rx[1] = 0.6123724356957944*G_1c[19]+0.3535533905932737*G_1c[15]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[9]-0.6123724356957944*G_1r[16]; 
  G_1r_lx[1] = 0.3535533905932737*G_1r[15]-0.6123724356957944*G_1r[19]; 

  out_F_0[9] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[9] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.369306393762915*F_0l[18])-1.060660171779821*F_0l[10]-0.6123724356957944*F_0l[5]; 
  F_0l_rx[1] = (-0.9486832980505137*F_0l[17])-1.369306393762915*F_0l[14]-0.5477225575051661*F_0l[13]-1.060660171779821*F_0l[6]-0.6123724356957944*F_0l[3]; 
  F_0l_rx[2] = (-1.224744871391589*F_0l[18])-0.9486832980505137*F_0l[10]-0.5477225575051661*F_0l[5]; 

  F_0c_lx[0] = (-1.369306393762915*F_0c[18])+1.060660171779821*F_0c[10]-0.6123724356957944*F_0c[5]; 
  F_0c_lx[1] = 0.9486832980505137*F_0c[17]-1.369306393762915*F_0c[14]-0.5477225575051661*F_0c[13]+1.060660171779821*F_0c[6]-0.6123724356957944*F_0c[3]; 
  F_0c_lx[2] = (-1.224744871391589*F_0c[18])+0.9486832980505137*F_0c[10]-0.5477225575051661*F_0c[5]; 

  F_0c_rx[0] = 1.369306393762915*F_0c[18]+1.060660171779821*F_0c[10]+0.6123724356957944*F_0c[5]; 
  F_0c_rx[1] = 0.9486832980505137*F_0c[17]+1.369306393762915*F_0c[14]+0.5477225575051661*F_0c[13]+1.060660171779821*F_0c[6]+0.6123724356957944*F_0c[3]; 
  F_0c_rx[2] = 1.224744871391589*F_0c[18]+0.9486832980505137*F_0c[10]+0.5477225575051661*F_0c[5]; 

  F_0r_lx[0] = 1.369306393762915*F_0r[18]-1.060660171779821*F_0r[10]+0.6123724356957944*F_0r[5]; 
  F_0r_lx[1] = (-0.9486832980505137*F_0r[17])+1.369306393762915*F_0r[14]+0.5477225575051661*F_0r[13]-1.060660171779821*F_0r[6]+0.6123724356957944*F_0r[3]; 
  F_0r_lx[2] = 1.224744871391589*F_0r[18]-0.9486832980505137*F_0r[10]+0.5477225575051661*F_0r[5]; 

  G_1l_rx[0] = (-1.369306393762915*G_1l[18])-1.060660171779821*G_1l[10]-0.6123724356957944*G_1l[5]; 
  G_1l_rx[1] = (-0.9486832980505137*G_1l[17])-1.369306393762915*G_1l[14]-0.5477225575051661*G_1l[13]-1.060660171779821*G_1l[6]-0.6123724356957944*G_1l[3]; 
  G_1l_rx[2] = (-1.224744871391589*G_1l[18])-0.9486832980505137*G_1l[10]-0.5477225575051661*G_1l[5]; 

  G_1c_lx[0] = (-1.369306393762915*G_1c[18])+1.060660171779821*G_1c[10]-0.6123724356957944*G_1c[5]; 
  G_1c_lx[1] = 0.9486832980505137*G_1c[17]-1.369306393762915*G_1c[14]-0.5477225575051661*G_1c[13]+1.060660171779821*G_1c[6]-0.6123724356957944*G_1c[3]; 
  G_1c_lx[2] = (-1.224744871391589*G_1c[18])+0.9486832980505137*G_1c[10]-0.5477225575051661*G_1c[5]; 

  G_1c_rx[0] = 1.369306393762915*G_1c[18]+1.060660171779821*G_1c[10]+0.6123724356957944*G_1c[5]; 
  G_1c_rx[1] = 0.9486832980505137*G_1c[17]+1.369306393762915*G_1c[14]+0.5477225575051661*G_1c[13]+1.060660171779821*G_1c[6]+0.6123724356957944*G_1c[3]; 
  G_1c_rx[2] = 1.224744871391589*G_1c[18]+0.9486832980505137*G_1c[10]+0.5477225575051661*G_1c[5]; 

  G_1r_lx[0] = 1.369306393762915*G_1r[18]-1.060660171779821*G_1r[10]+0.6123724356957944*G_1r[5]; 
  G_1r_lx[1] = (-0.9486832980505137*G_1r[17])+1.369306393762915*G_1r[14]+0.5477225575051661*G_1r[13]-1.060660171779821*G_1r[6]+0.6123724356957944*G_1r[3]; 
  G_1r_lx[2] = 1.224744871391589*G_1r[18]-0.9486832980505137*G_1r[10]+0.5477225575051661*G_1r[5]; 

  out_F_0[10] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[10] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[11])-0.6123724356957944*F_0l[7]; 
  F_0l_rx[1] = (-1.224744871391589*F_0l[12])-0.9486832980505138*F_0l[4]-0.5477225575051661*F_0l[1]; 
  F_0l_rx[2] = (-0.6776309271789384*F_0l[11])-1.369306393762915*F_0l[8]-0.3912303982179757*F_0l[7]-1.060660171779821*F_0l[2]-0.6123724356957944*F_0l[0]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[11]-0.6123724356957944*F_0c[7]; 
  F_0c_lx[1] = (-1.224744871391589*F_0c[12])+0.9486832980505138*F_0c[4]-0.5477225575051661*F_0c[1]; 
  F_0c_lx[2] = 0.6776309271789384*F_0c[11]-1.369306393762915*F_0c[8]-0.3912303982179757*F_0c[7]+1.060660171779821*F_0c[2]-0.6123724356957944*F_0c[0]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[11]+0.6123724356957944*F_0c[7]; 
  F_0c_rx[1] = 1.224744871391589*F_0c[12]+0.9486832980505138*F_0c[4]+0.5477225575051661*F_0c[1]; 
  F_0c_rx[2] = 0.6776309271789384*F_0c[11]+1.369306393762915*F_0c[8]+0.3912303982179757*F_0c[7]+1.060660171779821*F_0c[2]+0.6123724356957944*F_0c[0]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[7]-1.060660171779821*F_0r[11]; 
  F_0r_lx[1] = 1.224744871391589*F_0r[12]-0.9486832980505138*F_0r[4]+0.5477225575051661*F_0r[1]; 
  F_0r_lx[2] = (-0.6776309271789384*F_0r[11])+1.369306393762915*F_0r[8]+0.3912303982179757*F_0r[7]-1.060660171779821*F_0r[2]+0.6123724356957944*F_0r[0]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[11])-0.6123724356957944*G_1l[7]; 
  G_1l_rx[1] = (-1.224744871391589*G_1l[12])-0.9486832980505138*G_1l[4]-0.5477225575051661*G_1l[1]; 
  G_1l_rx[2] = (-0.6776309271789384*G_1l[11])-1.369306393762915*G_1l[8]-0.3912303982179757*G_1l[7]-1.060660171779821*G_1l[2]-0.6123724356957944*G_1l[0]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[11]-0.6123724356957944*G_1c[7]; 
  G_1c_lx[1] = (-1.224744871391589*G_1c[12])+0.9486832980505138*G_1c[4]-0.5477225575051661*G_1c[1]; 
  G_1c_lx[2] = 0.6776309271789384*G_1c[11]-1.369306393762915*G_1c[8]-0.3912303982179757*G_1c[7]+1.060660171779821*G_1c[2]-0.6123724356957944*G_1c[0]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[11]+0.6123724356957944*G_1c[7]; 
  G_1c_rx[1] = 1.224744871391589*G_1c[12]+0.9486832980505138*G_1c[4]+0.5477225575051661*G_1c[1]; 
  G_1c_rx[2] = 0.6776309271789384*G_1c[11]+1.369306393762915*G_1c[8]+0.3912303982179757*G_1c[7]+1.060660171779821*G_1c[2]+0.6123724356957944*G_1c[0]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[7]-1.060660171779821*G_1r[11]; 
  G_1r_lx[1] = 1.224744871391589*G_1r[12]-0.9486832980505138*G_1r[4]+0.5477225575051661*G_1r[1]; 
  G_1r_lx[2] = (-0.6776309271789384*G_1r[11])+1.369306393762915*G_1r[8]+0.3912303982179757*G_1r[7]-1.060660171779821*G_1r[2]+0.6123724356957944*G_1r[0]; 

  out_F_0[11] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[11] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 1.767766952966368*F_0l[12]+1.369306393762915*F_0l[4]+0.7905694150420948*F_0l[1]; 
  F_0l_rx[1] = 1.224744871391589*F_0l[11]+1.767766952966369*F_0l[8]+0.7071067811865475*F_0l[7]+1.369306393762915*F_0l[2]+0.7905694150420948*F_0l[0]; 
  F_0l_rx[2] = 1.58113883008419*F_0l[12]+1.224744871391589*F_0l[4]+0.7071067811865475*F_0l[1]; 

  F_0c_lx[0] = 1.767766952966368*F_0c[12]-1.369306393762915*F_0c[4]+0.7905694150420948*F_0c[1]; 
  F_0c_lx[1] = (-1.224744871391589*F_0c[11])+1.767766952966369*F_0c[8]+0.7071067811865475*F_0c[7]-1.369306393762915*F_0c[2]+0.7905694150420948*F_0c[0]; 
  F_0c_lx[2] = 1.58113883008419*F_0c[12]-1.224744871391589*F_0c[4]+0.7071067811865475*F_0c[1]; 

  F_0c_rx[0] = 1.767766952966368*F_0c[12]+1.369306393762915*F_0c[4]+0.7905694150420948*F_0c[1]; 
  F_0c_rx[1] = 1.224744871391589*F_0c[11]+1.767766952966369*F_0c[8]+0.7071067811865475*F_0c[7]+1.369306393762915*F_0c[2]+0.7905694150420948*F_0c[0]; 
  F_0c_rx[2] = 1.58113883008419*F_0c[12]+1.224744871391589*F_0c[4]+0.7071067811865475*F_0c[1]; 

  F_0r_lx[0] = 1.767766952966368*F_0r[12]-1.369306393762915*F_0r[4]+0.7905694150420948*F_0r[1]; 
  F_0r_lx[1] = (-1.224744871391589*F_0r[11])+1.767766952966369*F_0r[8]+0.7071067811865475*F_0r[7]-1.369306393762915*F_0r[2]+0.7905694150420948*F_0r[0]; 
  F_0r_lx[2] = 1.58113883008419*F_0r[12]-1.224744871391589*F_0r[4]+0.7071067811865475*F_0r[1]; 

  G_1l_rx[0] = 1.767766952966368*G_1l[12]+1.369306393762915*G_1l[4]+0.7905694150420948*G_1l[1]; 
  G_1l_rx[1] = 1.224744871391589*G_1l[11]+1.767766952966369*G_1l[8]+0.7071067811865475*G_1l[7]+1.369306393762915*G_1l[2]+0.7905694150420948*G_1l[0]; 
  G_1l_rx[2] = 1.58113883008419*G_1l[12]+1.224744871391589*G_1l[4]+0.7071067811865475*G_1l[1]; 

  G_1c_lx[0] = 1.767766952966368*G_1c[12]-1.369306393762915*G_1c[4]+0.7905694150420948*G_1c[1]; 
  G_1c_lx[1] = (-1.224744871391589*G_1c[11])+1.767766952966369*G_1c[8]+0.7071067811865475*G_1c[7]-1.369306393762915*G_1c[2]+0.7905694150420948*G_1c[0]; 
  G_1c_lx[2] = 1.58113883008419*G_1c[12]-1.224744871391589*G_1c[4]+0.7071067811865475*G_1c[1]; 

  G_1c_rx[0] = 1.767766952966368*G_1c[12]+1.369306393762915*G_1c[4]+0.7905694150420948*G_1c[1]; 
  G_1c_rx[1] = 1.224744871391589*G_1c[11]+1.767766952966369*G_1c[8]+0.7071067811865475*G_1c[7]+1.369306393762915*G_1c[2]+0.7905694150420948*G_1c[0]; 
  G_1c_rx[2] = 1.58113883008419*G_1c[12]+1.224744871391589*G_1c[4]+0.7071067811865475*G_1c[1]; 

  G_1r_lx[0] = 1.767766952966368*G_1r[12]-1.369306393762915*G_1r[4]+0.7905694150420948*G_1r[1]; 
  G_1r_lx[1] = (-1.224744871391589*G_1r[11])+1.767766952966369*G_1r[8]+0.7071067811865475*G_1r[7]-1.369306393762915*G_1r[2]+0.7905694150420948*G_1r[0]; 
  G_1r_lx[2] = 1.58113883008419*G_1r[12]-1.224744871391589*G_1r[4]+0.7071067811865475*G_1r[1]; 

  out_F_0[12] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[12] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.6123724356957944*F_0l[17]+0.3535533905932737*F_0l[13]; 
  F_0l_rx[1] = 0.7071067811865475*F_0l[18]+0.5477225575051661*F_0l[10]+0.3162277660168379*F_0l[5]; 
  F_0l_rx[2] = 0.3912303982179757*F_0l[17]+0.7905694150420947*F_0l[14]+0.2258769757263128*F_0l[13]+0.6123724356957944*F_0l[6]+0.3535533905932737*F_0l[3]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[13]-0.6123724356957944*F_0c[17]; 
  F_0c_lx[1] = 0.7071067811865475*F_0c[18]-0.5477225575051661*F_0c[10]+0.3162277660168379*F_0c[5]; 
  F_0c_lx[2] = (-0.3912303982179757*F_0c[17])+0.7905694150420947*F_0c[14]+0.2258769757263128*F_0c[13]-0.6123724356957944*F_0c[6]+0.3535533905932737*F_0c[3]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[17]+0.3535533905932737*F_0c[13]; 
  F_0c_rx[1] = 0.7071067811865475*F_0c[18]+0.5477225575051661*F_0c[10]+0.3162277660168379*F_0c[5]; 
  F_0c_rx[2] = 0.3912303982179757*F_0c[17]+0.7905694150420947*F_0c[14]+0.2258769757263128*F_0c[13]+0.6123724356957944*F_0c[6]+0.3535533905932737*F_0c[3]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[13]-0.6123724356957944*F_0r[17]; 
  F_0r_lx[1] = 0.7071067811865475*F_0r[18]-0.5477225575051661*F_0r[10]+0.3162277660168379*F_0r[5]; 
  F_0r_lx[2] = (-0.3912303982179757*F_0r[17])+0.7905694150420947*F_0r[14]+0.2258769757263128*F_0r[13]-0.6123724356957944*F_0r[6]+0.3535533905932737*F_0r[3]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[17]+0.3535533905932737*G_1l[13]; 
  G_1l_rx[1] = 0.7071067811865475*G_1l[18]+0.5477225575051661*G_1l[10]+0.3162277660168379*G_1l[5]; 
  G_1l_rx[2] = 0.3912303982179757*G_1l[17]+0.7905694150420947*G_1l[14]+0.2258769757263128*G_1l[13]+0.6123724356957944*G_1l[6]+0.3535533905932737*G_1l[3]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[13]-0.6123724356957944*G_1c[17]; 
  G_1c_lx[1] = 0.7071067811865475*G_1c[18]-0.5477225575051661*G_1c[10]+0.3162277660168379*G_1c[5]; 
  G_1c_lx[2] = (-0.3912303982179757*G_1c[17])+0.7905694150420947*G_1c[14]+0.2258769757263128*G_1c[13]-0.6123724356957944*G_1c[6]+0.3535533905932737*G_1c[3]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[17]+0.3535533905932737*G_1c[13]; 
  G_1c_rx[1] = 0.7071067811865475*G_1c[18]+0.5477225575051661*G_1c[10]+0.3162277660168379*G_1c[5]; 
  G_1c_rx[2] = 0.3912303982179757*G_1c[17]+0.7905694150420947*G_1c[14]+0.2258769757263128*G_1c[13]+0.6123724356957944*G_1c[6]+0.3535533905932737*G_1c[3]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[13]-0.6123724356957944*G_1r[17]; 
  G_1r_lx[1] = 0.7071067811865475*G_1r[18]-0.5477225575051661*G_1r[10]+0.3162277660168379*G_1r[5]; 
  G_1r_lx[2] = (-0.3912303982179757*G_1r[17])+0.7905694150420947*G_1r[14]+0.2258769757263128*G_1r[13]-0.6123724356957944*G_1r[6]+0.3535533905932737*G_1r[3]; 

  out_F_0[13] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[13] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 1.767766952966368*F_0l[14]+1.369306393762915*F_0l[6]+0.7905694150420948*F_0l[3]; 
  F_0l_rx[1] = 1.767766952966369*F_0l[18]+1.369306393762915*F_0l[10]+0.7905694150420948*F_0l[5]; 
  F_0l_rx[2] = 1.369306393762915*F_0l[17]+0.7905694150420947*F_0l[13]; 

  F_0c_lx[0] = 1.767766952966368*F_0c[14]-1.369306393762915*F_0c[6]+0.7905694150420948*F_0c[3]; 
  F_0c_lx[1] = 1.767766952966369*F_0c[18]-1.369306393762915*F_0c[10]+0.7905694150420948*F_0c[5]; 
  F_0c_lx[2] = 0.7905694150420947*F_0c[13]-1.369306393762915*F_0c[17]; 

  F_0c_rx[0] = 1.767766952966368*F_0c[14]+1.369306393762915*F_0c[6]+0.7905694150420948*F_0c[3]; 
  F_0c_rx[1] = 1.767766952966369*F_0c[18]+1.369306393762915*F_0c[10]+0.7905694150420948*F_0c[5]; 
  F_0c_rx[2] = 1.369306393762915*F_0c[17]+0.7905694150420947*F_0c[13]; 

  F_0r_lx[0] = 1.767766952966368*F_0r[14]-1.369306393762915*F_0r[6]+0.7905694150420948*F_0r[3]; 
  F_0r_lx[1] = 1.767766952966369*F_0r[18]-1.369306393762915*F_0r[10]+0.7905694150420948*F_0r[5]; 
  F_0r_lx[2] = 0.7905694150420947*F_0r[13]-1.369306393762915*F_0r[17]; 

  G_1l_rx[0] = 1.767766952966368*G_1l[14]+1.369306393762915*G_1l[6]+0.7905694150420948*G_1l[3]; 
  G_1l_rx[1] = 1.767766952966369*G_1l[18]+1.369306393762915*G_1l[10]+0.7905694150420948*G_1l[5]; 
  G_1l_rx[2] = 1.369306393762915*G_1l[17]+0.7905694150420947*G_1l[13]; 

  G_1c_lx[0] = 1.767766952966368*G_1c[14]-1.369306393762915*G_1c[6]+0.7905694150420948*G_1c[3]; 
  G_1c_lx[1] = 1.767766952966369*G_1c[18]-1.369306393762915*G_1c[10]+0.7905694150420948*G_1c[5]; 
  G_1c_lx[2] = 0.7905694150420947*G_1c[13]-1.369306393762915*G_1c[17]; 

  G_1c_rx[0] = 1.767766952966368*G_1c[14]+1.369306393762915*G_1c[6]+0.7905694150420948*G_1c[3]; 
  G_1c_rx[1] = 1.767766952966369*G_1c[18]+1.369306393762915*G_1c[10]+0.7905694150420948*G_1c[5]; 
  G_1c_rx[2] = 1.369306393762915*G_1c[17]+0.7905694150420947*G_1c[13]; 

  G_1r_lx[0] = 1.767766952966368*G_1r[14]-1.369306393762915*G_1r[6]+0.7905694150420948*G_1r[3]; 
  G_1r_lx[1] = 1.767766952966369*G_1r[18]-1.369306393762915*G_1r[10]+0.7905694150420948*G_1r[5]; 
  G_1r_lx[2] = 0.7905694150420947*G_1r[13]-1.369306393762915*G_1r[17]; 

  out_F_0[14] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[14] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.6123724356957944*F_0l[19]+0.3535533905932737*F_0l[15]; 
  F_0l_rx[1] = 0.6123724356957944*F_0l[16]+0.3535533905932737*F_0l[9]; 
  F_0l_rx[2] = 0.5477225575051661*F_0l[19]+0.3162277660168379*F_0l[15]; 

  F_0c_lx[0] = 0.3535533905932737*F_0c[15]-0.6123724356957944*F_0c[19]; 
  F_0c_lx[1] = 0.3535533905932737*F_0c[9]-0.6123724356957944*F_0c[16]; 
  F_0c_lx[2] = 0.3162277660168379*F_0c[15]-0.5477225575051661*F_0c[19]; 

  F_0c_rx[0] = 0.6123724356957944*F_0c[19]+0.3535533905932737*F_0c[15]; 
  F_0c_rx[1] = 0.6123724356957944*F_0c[16]+0.3535533905932737*F_0c[9]; 
  F_0c_rx[2] = 0.5477225575051661*F_0c[19]+0.3162277660168379*F_0c[15]; 

  F_0r_lx[0] = 0.3535533905932737*F_0r[15]-0.6123724356957944*F_0r[19]; 
  F_0r_lx[1] = 0.3535533905932737*F_0r[9]-0.6123724356957944*F_0r[16]; 
  F_0r_lx[2] = 0.3162277660168379*F_0r[15]-0.5477225575051661*F_0r[19]; 

  G_1l_rx[0] = 0.6123724356957944*G_1l[19]+0.3535533905932737*G_1l[15]; 
  G_1l_rx[1] = 0.6123724356957944*G_1l[16]+0.3535533905932737*G_1l[9]; 
  G_1l_rx[2] = 0.5477225575051661*G_1l[19]+0.3162277660168379*G_1l[15]; 

  G_1c_lx[0] = 0.3535533905932737*G_1c[15]-0.6123724356957944*G_1c[19]; 
  G_1c_lx[1] = 0.3535533905932737*G_1c[9]-0.6123724356957944*G_1c[16]; 
  G_1c_lx[2] = 0.3162277660168379*G_1c[15]-0.5477225575051661*G_1c[19]; 

  G_1c_rx[0] = 0.6123724356957944*G_1c[19]+0.3535533905932737*G_1c[15]; 
  G_1c_rx[1] = 0.6123724356957944*G_1c[16]+0.3535533905932737*G_1c[9]; 
  G_1c_rx[2] = 0.5477225575051661*G_1c[19]+0.3162277660168379*G_1c[15]; 

  G_1r_lx[0] = 0.3535533905932737*G_1r[15]-0.6123724356957944*G_1r[19]; 
  G_1r_lx[1] = 0.3535533905932737*G_1r[9]-0.6123724356957944*G_1r[16]; 
  G_1r_lx[2] = 0.3162277660168379*G_1r[15]-0.5477225575051661*G_1r[19]; 

  out_F_0[15] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[15] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[16])-0.6123724356957944*F_0l[9]; 
  F_0l_rx[1] = (-1.060660171779821*F_0l[19])-0.6123724356957944*F_0l[15]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[16]-0.6123724356957944*F_0c[9]; 
  F_0c_lx[1] = 1.060660171779821*F_0c[19]-0.6123724356957944*F_0c[15]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[16]+0.6123724356957944*F_0c[9]; 
  F_0c_rx[1] = 1.060660171779821*F_0c[19]+0.6123724356957944*F_0c[15]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[9]-1.060660171779821*F_0r[16]; 
  F_0r_lx[1] = 0.6123724356957944*F_0r[15]-1.060660171779821*F_0r[19]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[16])-0.6123724356957944*G_1l[9]; 
  G_1l_rx[1] = (-1.060660171779821*G_1l[19])-0.6123724356957944*G_1l[15]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[16]-0.6123724356957944*G_1c[9]; 
  G_1c_lx[1] = 1.060660171779821*G_1c[19]-0.6123724356957944*G_1c[15]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[16]+0.6123724356957944*G_1c[9]; 
  G_1c_rx[1] = 1.060660171779821*G_1c[19]+0.6123724356957944*G_1c[15]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[9]-1.060660171779821*G_1r[16]; 
  G_1r_lx[1] = 0.6123724356957944*G_1r[15]-1.060660171779821*G_1r[19]; 

  out_F_0[16] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[16] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[17])-0.6123724356957944*F_0l[13]; 
  F_0l_rx[1] = (-1.224744871391589*F_0l[18])-0.9486832980505137*F_0l[10]-0.5477225575051661*F_0l[5]; 
  F_0l_rx[2] = (-0.6776309271789384*F_0l[17])-1.369306393762915*F_0l[14]-0.3912303982179757*F_0l[13]-1.060660171779821*F_0l[6]-0.6123724356957944*F_0l[3]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[17]-0.6123724356957944*F_0c[13]; 
  F_0c_lx[1] = (-1.224744871391589*F_0c[18])+0.9486832980505137*F_0c[10]-0.5477225575051661*F_0c[5]; 
  F_0c_lx[2] = 0.6776309271789384*F_0c[17]-1.369306393762915*F_0c[14]-0.3912303982179757*F_0c[13]+1.060660171779821*F_0c[6]-0.6123724356957944*F_0c[3]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[17]+0.6123724356957944*F_0c[13]; 
  F_0c_rx[1] = 1.224744871391589*F_0c[18]+0.9486832980505137*F_0c[10]+0.5477225575051661*F_0c[5]; 
  F_0c_rx[2] = 0.6776309271789384*F_0c[17]+1.369306393762915*F_0c[14]+0.3912303982179757*F_0c[13]+1.060660171779821*F_0c[6]+0.6123724356957944*F_0c[3]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[13]-1.060660171779821*F_0r[17]; 
  F_0r_lx[1] = 1.224744871391589*F_0r[18]-0.9486832980505137*F_0r[10]+0.5477225575051661*F_0r[5]; 
  F_0r_lx[2] = (-0.6776309271789384*F_0r[17])+1.369306393762915*F_0r[14]+0.3912303982179757*F_0r[13]-1.060660171779821*F_0r[6]+0.6123724356957944*F_0r[3]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[17])-0.6123724356957944*G_1l[13]; 
  G_1l_rx[1] = (-1.224744871391589*G_1l[18])-0.9486832980505137*G_1l[10]-0.5477225575051661*G_1l[5]; 
  G_1l_rx[2] = (-0.6776309271789384*G_1l[17])-1.369306393762915*G_1l[14]-0.3912303982179757*G_1l[13]-1.060660171779821*G_1l[6]-0.6123724356957944*G_1l[3]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[17]-0.6123724356957944*G_1c[13]; 
  G_1c_lx[1] = (-1.224744871391589*G_1c[18])+0.9486832980505137*G_1c[10]-0.5477225575051661*G_1c[5]; 
  G_1c_lx[2] = 0.6776309271789384*G_1c[17]-1.369306393762915*G_1c[14]-0.3912303982179757*G_1c[13]+1.060660171779821*G_1c[6]-0.6123724356957944*G_1c[3]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[17]+0.6123724356957944*G_1c[13]; 
  G_1c_rx[1] = 1.224744871391589*G_1c[18]+0.9486832980505137*G_1c[10]+0.5477225575051661*G_1c[5]; 
  G_1c_rx[2] = 0.6776309271789384*G_1c[17]+1.369306393762915*G_1c[14]+0.3912303982179757*G_1c[13]+1.060660171779821*G_1c[6]+0.6123724356957944*G_1c[3]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[13]-1.060660171779821*G_1r[17]; 
  G_1r_lx[1] = 1.224744871391589*G_1r[18]-0.9486832980505137*G_1r[10]+0.5477225575051661*G_1r[5]; 
  G_1r_lx[2] = (-0.6776309271789384*G_1r[17])+1.369306393762915*G_1r[14]+0.3912303982179757*G_1r[13]-1.060660171779821*G_1r[6]+0.6123724356957944*G_1r[3]; 

  out_F_0[17] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[17] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 1.767766952966368*F_0l[18]+1.369306393762915*F_0l[10]+0.7905694150420947*F_0l[5]; 
  F_0l_rx[1] = 1.224744871391589*F_0l[17]+1.767766952966369*F_0l[14]+0.7071067811865475*F_0l[13]+1.369306393762915*F_0l[6]+0.7905694150420947*F_0l[3]; 
  F_0l_rx[2] = 1.58113883008419*F_0l[18]+1.224744871391589*F_0l[10]+0.7071067811865475*F_0l[5]; 

  F_0c_lx[0] = 1.767766952966368*F_0c[18]-1.369306393762915*F_0c[10]+0.7905694150420947*F_0c[5]; 
  F_0c_lx[1] = (-1.224744871391589*F_0c[17])+1.767766952966369*F_0c[14]+0.7071067811865475*F_0c[13]-1.369306393762915*F_0c[6]+0.7905694150420947*F_0c[3]; 
  F_0c_lx[2] = 1.58113883008419*F_0c[18]-1.224744871391589*F_0c[10]+0.7071067811865475*F_0c[5]; 

  F_0c_rx[0] = 1.767766952966368*F_0c[18]+1.369306393762915*F_0c[10]+0.7905694150420947*F_0c[5]; 
  F_0c_rx[1] = 1.224744871391589*F_0c[17]+1.767766952966369*F_0c[14]+0.7071067811865475*F_0c[13]+1.369306393762915*F_0c[6]+0.7905694150420947*F_0c[3]; 
  F_0c_rx[2] = 1.58113883008419*F_0c[18]+1.224744871391589*F_0c[10]+0.7071067811865475*F_0c[5]; 

  F_0r_lx[0] = 1.767766952966368*F_0r[18]-1.369306393762915*F_0r[10]+0.7905694150420947*F_0r[5]; 
  F_0r_lx[1] = (-1.224744871391589*F_0r[17])+1.767766952966369*F_0r[14]+0.7071067811865475*F_0r[13]-1.369306393762915*F_0r[6]+0.7905694150420947*F_0r[3]; 
  F_0r_lx[2] = 1.58113883008419*F_0r[18]-1.224744871391589*F_0r[10]+0.7071067811865475*F_0r[5]; 

  G_1l_rx[0] = 1.767766952966368*G_1l[18]+1.369306393762915*G_1l[10]+0.7905694150420947*G_1l[5]; 
  G_1l_rx[1] = 1.224744871391589*G_1l[17]+1.767766952966369*G_1l[14]+0.7071067811865475*G_1l[13]+1.369306393762915*G_1l[6]+0.7905694150420947*G_1l[3]; 
  G_1l_rx[2] = 1.58113883008419*G_1l[18]+1.224744871391589*G_1l[10]+0.7071067811865475*G_1l[5]; 

  G_1c_lx[0] = 1.767766952966368*G_1c[18]-1.369306393762915*G_1c[10]+0.7905694150420947*G_1c[5]; 
  G_1c_lx[1] = (-1.224744871391589*G_1c[17])+1.767766952966369*G_1c[14]+0.7071067811865475*G_1c[13]-1.369306393762915*G_1c[6]+0.7905694150420947*G_1c[3]; 
  G_1c_lx[2] = 1.58113883008419*G_1c[18]-1.224744871391589*G_1c[10]+0.7071067811865475*G_1c[5]; 

  G_1c_rx[0] = 1.767766952966368*G_1c[18]+1.369306393762915*G_1c[10]+0.7905694150420947*G_1c[5]; 
  G_1c_rx[1] = 1.224744871391589*G_1c[17]+1.767766952966369*G_1c[14]+0.7071067811865475*G_1c[13]+1.369306393762915*G_1c[6]+0.7905694150420947*G_1c[3]; 
  G_1c_rx[2] = 1.58113883008419*G_1c[18]+1.224744871391589*G_1c[10]+0.7071067811865475*G_1c[5]; 

  G_1r_lx[0] = 1.767766952966368*G_1r[18]-1.369306393762915*G_1r[10]+0.7905694150420947*G_1r[5]; 
  G_1r_lx[1] = (-1.224744871391589*G_1r[17])+1.767766952966369*G_1r[14]+0.7071067811865475*G_1r[13]-1.369306393762915*G_1r[6]+0.7905694150420947*G_1r[3]; 
  G_1r_lx[2] = 1.58113883008419*G_1r[18]-1.224744871391589*G_1r[10]+0.7071067811865475*G_1r[5]; 

  out_F_0[18] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[18] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-1.060660171779821*F_0l[19])-0.6123724356957944*F_0l[15]; 
  F_0l_rx[1] = (-1.060660171779821*F_0l[16])-0.6123724356957944*F_0l[9]; 
  F_0l_rx[2] = (-0.9486832980505137*F_0l[19])-0.5477225575051661*F_0l[15]; 

  F_0c_lx[0] = 1.060660171779821*F_0c[19]-0.6123724356957944*F_0c[15]; 
  F_0c_lx[1] = 1.060660171779821*F_0c[16]-0.6123724356957944*F_0c[9]; 
  F_0c_lx[2] = 0.9486832980505137*F_0c[19]-0.5477225575051661*F_0c[15]; 

  F_0c_rx[0] = 1.060660171779821*F_0c[19]+0.6123724356957944*F_0c[15]; 
  F_0c_rx[1] = 1.060660171779821*F_0c[16]+0.6123724356957944*F_0c[9]; 
  F_0c_rx[2] = 0.9486832980505137*F_0c[19]+0.5477225575051661*F_0c[15]; 

  F_0r_lx[0] = 0.6123724356957944*F_0r[15]-1.060660171779821*F_0r[19]; 
  F_0r_lx[1] = 0.6123724356957944*F_0r[9]-1.060660171779821*F_0r[16]; 
  F_0r_lx[2] = 0.5477225575051661*F_0r[15]-0.9486832980505137*F_0r[19]; 

  G_1l_rx[0] = (-1.060660171779821*G_1l[19])-0.6123724356957944*G_1l[15]; 
  G_1l_rx[1] = (-1.060660171779821*G_1l[16])-0.6123724356957944*G_1l[9]; 
  G_1l_rx[2] = (-0.9486832980505137*G_1l[19])-0.5477225575051661*G_1l[15]; 

  G_1c_lx[0] = 1.060660171779821*G_1c[19]-0.6123724356957944*G_1c[15]; 
  G_1c_lx[1] = 1.060660171779821*G_1c[16]-0.6123724356957944*G_1c[9]; 
  G_1c_lx[2] = 0.9486832980505137*G_1c[19]-0.5477225575051661*G_1c[15]; 

  G_1c_rx[0] = 1.060660171779821*G_1c[19]+0.6123724356957944*G_1c[15]; 
  G_1c_rx[1] = 1.060660171779821*G_1c[16]+0.6123724356957944*G_1c[9]; 
  G_1c_rx[2] = 0.9486832980505137*G_1c[19]+0.5477225575051661*G_1c[15]; 

  G_1r_lx[0] = 0.6123724356957944*G_1r[15]-1.060660171779821*G_1r[19]; 
  G_1r_lx[1] = 0.6123724356957944*G_1r[9]-1.060660171779821*G_1r[16]; 
  G_1r_lx[2] = 0.5477225575051661*G_1r[15]-0.9486832980505137*G_1r[19]; 

  out_F_0[19] += 0.5*dx1*(F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[19] += 0.5*dx1*(G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  return 0.;

} 
