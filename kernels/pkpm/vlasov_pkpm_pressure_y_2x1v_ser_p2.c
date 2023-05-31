#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_ser_3x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_pressure_y_2x1v_ser_p2(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                 Cell-center coordinates.
  // dxv[NDIM]:               Cell spacing.
  // bvarl/bvarc/bvarr:       Input magnetic field unit vector in left/center/right cells.
  // fl/fc/fr:                Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells.
  // out:                     Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[1]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double volFact = dxv[2]/2.0; 
  const double *bl = &bvarl[8]; 
  const double *bc = &bvarc[8]; 
  const double *br = &bvarr[8]; 
  const double *F_0l = &fl[0]; 
  const double *F_0c = &fc[0]; 
  const double *F_0r = &fr[0]; 
  double *out_pressure = &out[0]; 
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
  if ((-0.5999999999999999*alphaSurf_l[6])+0.4472135954999579*alphaSurf_l[4]+0.9*alphaSurf_l[3]-0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(F_0c); 
  } 
  if ((-0.5999999999999999*alphaSurf_r[6])+0.4472135954999579*alphaSurf_r[4]+0.9*alphaSurf_r[3]-0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(F_0r); 
  } 
  if (0.4472135954999579*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(F_0c); 
  } 
  if (0.4472135954999579*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(F_0r); 
  } 
  if (0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[4]-0.9*alphaSurf_l[3]+0.6708203932499369*alphaSurf_l[2]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(F_0c); 
  } 
  if (0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[4]-0.9*alphaSurf_r[3]+0.6708203932499369*alphaSurf_r[2]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(F_0r); 
  } 
  if (0.75*alphaSurf_l[6]-0.5590169943749475*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(F_0c); 
  } 
  if (0.75*alphaSurf_r[6]-0.5590169943749475*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(F_0r); 
  } 
  if (0.5*alphaSurf_l[0]-0.5590169943749475*alphaSurf_l[4] > 0) { 
    F_0_UpwindQuad_l[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(F_0c); 
  } 
  if (0.5*alphaSurf_r[0]-0.5590169943749475*alphaSurf_r[4] > 0) { 
    F_0_UpwindQuad_r[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(F_0r); 
  } 
  if ((-0.75*alphaSurf_l[6])-0.5590169943749475*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(F_0c); 
  } 
  if ((-0.75*alphaSurf_r[6])-0.5590169943749475*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(F_0r); 
  } 
  if ((-0.5999999999999999*alphaSurf_l[6])+0.4472135954999579*alphaSurf_l[4]-0.9*alphaSurf_l[3]-0.6708203932499369*alphaSurf_l[2]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(F_0c); 
  } 
  if ((-0.5999999999999999*alphaSurf_r[6])+0.4472135954999579*alphaSurf_r[4]-0.9*alphaSurf_r[3]-0.6708203932499369*alphaSurf_r[2]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(F_0r); 
  } 
  if (0.4472135954999579*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(F_0c); 
  } 
  if (0.4472135954999579*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(F_0r); 
  } 
  if (0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[4]+0.9*alphaSurf_l[3]+0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(F_0c); 
  } 
  if (0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[4]+0.9*alphaSurf_r[3]+0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(F_0r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  ser_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 

  Ghat_F_0_l[0] = 0.5*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5*F_0_Upwind_l[4]*alphaSurf_l[4]+0.5*F_0_Upwind_l[3]*alphaSurf_l[3]+0.5*F_0_Upwind_l[2]*alphaSurf_l[2]+0.5*F_0_Upwind_l[1]*alphaSurf_l[1]+0.5*F_0_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_F_0_l[1] = 0.447213595499958*F_0_Upwind_l[3]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[6]+0.4472135954999579*F_0_Upwind_l[1]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[1]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[7]+0.5000000000000001*F_0_Upwind_l[4]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[4]*F_0_Upwind_l[6]+0.4472135954999579*alphaSurf_l[2]*F_0_Upwind_l[5]+0.5*F_0_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[3] = 0.4*alphaSurf_l[6]*F_0_Upwind_l[7]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[7]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[6]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[5]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[4] = 0.31943828249997*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5000000000000001*F_0_Upwind_l[2]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[2]*F_0_Upwind_l[6]+0.31943828249997*F_0_Upwind_l[4]*alphaSurf_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[4]+0.5*alphaSurf_l[0]*F_0_Upwind_l[4]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[1]*alphaSurf_l[1]; 
  Ghat_F_0_l[5] = 0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[7]+0.4472135954999579*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5*alphaSurf_l[0]*F_0_Upwind_l[5]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_F_0_l[6] = 0.4*alphaSurf_l[3]*F_0_Upwind_l[7]+0.4472135954999579*F_0_Upwind_l[5]*alphaSurf_l[6]+0.31943828249997*F_0_Upwind_l[4]*alphaSurf_l[6]+0.5*F_0_Upwind_l[0]*alphaSurf_l[6]+0.31943828249997*alphaSurf_l[4]*F_0_Upwind_l[6]+0.5*alphaSurf_l[0]*F_0_Upwind_l[6]+0.5000000000000001*F_0_Upwind_l[2]*alphaSurf_l[4]+0.5000000000000001*alphaSurf_l[2]*F_0_Upwind_l[4]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[7] = 0.4472135954999579*alphaSurf_l[4]*F_0_Upwind_l[7]+0.5*alphaSurf_l[0]*F_0_Upwind_l[7]+0.4*F_0_Upwind_l[3]*alphaSurf_l[6]+0.4*alphaSurf_l[3]*F_0_Upwind_l[6]+0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[5]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[3]; 

  Ghat_F_0_r[0] = 0.5*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5*F_0_Upwind_r[4]*alphaSurf_r[4]+0.5*F_0_Upwind_r[3]*alphaSurf_r[3]+0.5*F_0_Upwind_r[2]*alphaSurf_r[2]+0.5*F_0_Upwind_r[1]*alphaSurf_r[1]+0.5*F_0_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_F_0_r[1] = 0.447213595499958*F_0_Upwind_r[3]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[6]+0.4472135954999579*F_0_Upwind_r[1]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[1]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[7]+0.5000000000000001*F_0_Upwind_r[4]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[4]*F_0_Upwind_r[6]+0.4472135954999579*alphaSurf_r[2]*F_0_Upwind_r[5]+0.5*F_0_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[3] = 0.4*alphaSurf_r[6]*F_0_Upwind_r[7]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[7]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[6]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[5]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[4] = 0.31943828249997*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5000000000000001*F_0_Upwind_r[2]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[2]*F_0_Upwind_r[6]+0.31943828249997*F_0_Upwind_r[4]*alphaSurf_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[4]+0.5*alphaSurf_r[0]*F_0_Upwind_r[4]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[1]*alphaSurf_r[1]; 
  Ghat_F_0_r[5] = 0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[7]+0.4472135954999579*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5*alphaSurf_r[0]*F_0_Upwind_r[5]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_F_0_r[6] = 0.4*alphaSurf_r[3]*F_0_Upwind_r[7]+0.4472135954999579*F_0_Upwind_r[5]*alphaSurf_r[6]+0.31943828249997*F_0_Upwind_r[4]*alphaSurf_r[6]+0.5*F_0_Upwind_r[0]*alphaSurf_r[6]+0.31943828249997*alphaSurf_r[4]*F_0_Upwind_r[6]+0.5*alphaSurf_r[0]*F_0_Upwind_r[6]+0.5000000000000001*F_0_Upwind_r[2]*alphaSurf_r[4]+0.5000000000000001*alphaSurf_r[2]*F_0_Upwind_r[4]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[7] = 0.4472135954999579*alphaSurf_r[4]*F_0_Upwind_r[7]+0.5*alphaSurf_r[0]*F_0_Upwind_r[7]+0.4*F_0_Upwind_r[3]*alphaSurf_r[6]+0.4*alphaSurf_r[3]*F_0_Upwind_r[6]+0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[5]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[3]; 

  out_pressure[0] += dx1*volFact*((Ghat_F_0_r[0]-1.0*Ghat_F_0_l[0])*wvpar+(0.2886751345948129*Ghat_F_0_r[2]-0.2886751345948129*Ghat_F_0_l[2])*dvpar); 
  out_pressure[1] += dx1*volFact*((Ghat_F_0_r[1]-1.0*Ghat_F_0_l[1])*wvpar+(0.2886751345948129*Ghat_F_0_r[3]-0.2886751345948129*Ghat_F_0_l[3])*dvpar); 
  out_pressure[2] += dx1*volFact*((1.732050807568877*(Ghat_F_0_r[0]+Ghat_F_0_l[0])-0.8660254037844386*(F_0c[18]*alpha_c[18]+F_0c[17]*alpha_c[17]+F_0c[14]*alpha_c[14]+F_0c[13]*alpha_c[13]+F_0c[12]*alpha_c[12]+F_0c[11]*alpha_c[11]+F_0c[10]*alpha_c[10]+F_0c[8]*alpha_c[8]+F_0c[7]*alpha_c[7]+F_0c[6]*alpha_c[6]+F_0c[5]*alpha_c[5]+F_0c[4]*alpha_c[4]+F_0c[3]*alpha_c[3]+F_0c[2]*alpha_c[2]+F_0c[1]*alpha_c[1]+F_0c[0]*alpha_c[0]))*wvpar+((-0.223606797749979*alpha_c[10]*F_0c[19])-0.2500000000000001*(F_0c[12]*alpha_c[18]+alpha_c[12]*F_0c[18]+F_0c[11]*alpha_c[17]+alpha_c[11]*F_0c[17])-0.223606797749979*(alpha_c[6]*F_0c[16]+alpha_c[5]*F_0c[15])-0.2500000000000001*(F_0c[8]*alpha_c[14]+alpha_c[8]*F_0c[14]+F_0c[7]*alpha_c[13]+alpha_c[7]*F_0c[13])-0.25*(F_0c[4]*alpha_c[10]+alpha_c[4]*F_0c[10])-0.223606797749979*alpha_c[3]*F_0c[9]-0.25*(F_0c[2]*alpha_c[6]+alpha_c[2]*F_0c[6]+F_0c[1]*alpha_c[5]+alpha_c[1]*F_0c[5]+F_0c[0]*alpha_c[3]+alpha_c[0]*F_0c[3])+0.5*(Ghat_F_0_r[2]+Ghat_F_0_l[2]))*dvpar); 
  out_pressure[3] += dx1*volFact*(((-0.8660254037844387*(F_0c[14]*alpha_c[18]+alpha_c[14]*F_0c[18]))-0.7745966692414833*(F_0c[10]*alpha_c[17]+alpha_c[10]*F_0c[17])-0.7745966692414834*(F_0c[5]*alpha_c[13]+alpha_c[5]*F_0c[13])-0.8660254037844387*(F_0c[8]*alpha_c[12]+alpha_c[8]*F_0c[12])-0.7745966692414834*(F_0c[4]*alpha_c[11]+alpha_c[4]*F_0c[11])-0.8660254037844386*(F_0c[6]*alpha_c[10]+alpha_c[6]*F_0c[10])-0.7745966692414833*(F_0c[1]*alpha_c[7]+alpha_c[1]*F_0c[7])-0.8660254037844386*(F_0c[3]*alpha_c[5]+alpha_c[3]*F_0c[5]+F_0c[2]*alpha_c[4]+alpha_c[2]*F_0c[4]+F_0c[0]*alpha_c[1])+1.732050807568877*(Ghat_F_0_r[1]+Ghat_F_0_l[1])-0.8660254037844386*alpha_c[0]*F_0c[1])*wvpar+(((-0.2*alpha_c[17])-0.223606797749979*alpha_c[6])*F_0c[19]-0.25*(F_0c[8]*alpha_c[18]+alpha_c[8]*F_0c[18])-0.223606797749979*(F_0c[4]*alpha_c[17]+alpha_c[4]*F_0c[17])-0.223606797749979*alpha_c[10]*F_0c[16]+((-0.2*alpha_c[13])-0.223606797749979*alpha_c[3])*F_0c[15]-0.25*(F_0c[12]*alpha_c[14]+alpha_c[12]*F_0c[14])-0.223606797749979*(F_0c[1]*alpha_c[13]+alpha_c[1]*F_0c[13]+F_0c[10]*alpha_c[11]+alpha_c[10]*F_0c[11])-0.25*(F_0c[2]*alpha_c[10]+alpha_c[2]*F_0c[10])-0.223606797749979*(alpha_c[5]*F_0c[9]+F_0c[5]*alpha_c[7]+alpha_c[5]*F_0c[7])-0.25*(F_0c[4]*alpha_c[6]+alpha_c[4]*F_0c[6]+F_0c[0]*alpha_c[5]+alpha_c[0]*F_0c[5]+F_0c[1]*alpha_c[3])+0.5*(Ghat_F_0_r[3]+Ghat_F_0_l[3])-0.25*alpha_c[1]*F_0c[3])*dvpar); 
  out_pressure[4] += dx1*volFact*((Ghat_F_0_r[4]-1.0*Ghat_F_0_l[4])*wvpar+(0.2886751345948129*Ghat_F_0_r[6]-0.2886751345948129*Ghat_F_0_l[6])*dvpar); 
  out_pressure[5] += dx1*volFact*(((-1.732050807568877*(F_0c[10]*alpha_c[18]+alpha_c[10]*F_0c[18]))-1.936491673103709*(F_0c[13]*alpha_c[17]+alpha_c[13]*F_0c[17])-1.732050807568877*(F_0c[6]*alpha_c[14]+alpha_c[6]*F_0c[14]+F_0c[4]*alpha_c[12]+alpha_c[4]*F_0c[12])-1.936491673103709*(F_0c[7]*alpha_c[11]+alpha_c[7]*F_0c[11]+F_0c[5]*alpha_c[10]+alpha_c[5]*F_0c[10])-1.732050807568877*(F_0c[2]*alpha_c[8]+alpha_c[2]*F_0c[8])-1.936491673103709*(F_0c[3]*alpha_c[6]+alpha_c[3]*F_0c[6]+F_0c[1]*alpha_c[4]+alpha_c[1]*F_0c[4]+F_0c[0]*alpha_c[2]+alpha_c[0]*F_0c[2])+2.23606797749979*Ghat_F_0_r[0]-2.23606797749979*Ghat_F_0_l[0])*wvpar+((-0.5*(alpha_c[5]*F_0c[19]+F_0c[4]*alpha_c[18]+alpha_c[4]*F_0c[18]))-0.4472135954999579*alpha_c[18]*F_0c[19]-0.5590169943749475*(F_0c[7]*alpha_c[17]+alpha_c[7]*F_0c[17])-0.5000000000000001*(alpha_c[3]*F_0c[16]+alpha_c[10]*F_0c[15]+F_0c[2]*alpha_c[14]+alpha_c[2]*F_0c[14])-0.4472135954999579*alpha_c[14]*F_0c[16]-0.5590169943749475*(F_0c[11]*alpha_c[13]+alpha_c[11]*F_0c[13])-0.5000000000000001*(F_0c[10]*alpha_c[12]+alpha_c[10]*F_0c[12])-0.5590169943749475*(F_0c[1]*alpha_c[10]+alpha_c[1]*F_0c[10])-0.5*(alpha_c[6]*F_0c[9]+F_0c[6]*alpha_c[8]+alpha_c[6]*F_0c[8])-0.5590169943749475*(F_0c[0]*alpha_c[6]+alpha_c[0]*F_0c[6]+F_0c[4]*alpha_c[5]+alpha_c[4]*F_0c[5]+F_0c[2]*alpha_c[3]+alpha_c[2]*F_0c[3])+0.6454972243679029*Ghat_F_0_r[2]-0.6454972243679029*Ghat_F_0_l[2])*dvpar); 
  out_pressure[6] += dx1*volFact*(((-0.7745966692414834*F_0c[18]*alpha_c[18])-0.8660254037844387*(F_0c[6]*alpha_c[17]+alpha_c[6]*F_0c[17])-0.5532833351724881*F_0c[17]*alpha_c[17]-0.8660254037844386*(F_0c[3]*alpha_c[13]+alpha_c[3]*F_0c[13])-0.5532833351724881*F_0c[13]*alpha_c[13]-0.7745966692414834*F_0c[12]*alpha_c[12]-0.8660254037844386*(F_0c[2]*alpha_c[11]+alpha_c[2]*F_0c[11])-0.5532833351724881*F_0c[11]*alpha_c[11]-0.7745966692414834*F_0c[10]*alpha_c[10]-0.8660254037844387*(F_0c[0]*alpha_c[7]+alpha_c[0]*F_0c[7])-0.5532833351724881*F_0c[7]*alpha_c[7]-0.7745966692414834*(F_0c[5]*alpha_c[5]+F_0c[4]*alpha_c[4])+1.732050807568877*(Ghat_F_0_r[4]+Ghat_F_0_l[4])-0.7745966692414834*F_0c[1]*alpha_c[1])*wvpar+((-0.2*alpha_c[10]*F_0c[19])-0.223606797749979*(F_0c[12]*alpha_c[18]+alpha_c[12]*F_0c[18])+((-0.223606797749979*F_0c[16])-0.159719141249985*F_0c[11]-0.2500000000000001*F_0c[2])*alpha_c[17]+((-0.159719141249985*alpha_c[11])-0.2500000000000001*alpha_c[2])*F_0c[17]-0.2*alpha_c[5]*F_0c[15]-0.159719141249985*(F_0c[7]*alpha_c[13]+alpha_c[7]*F_0c[13])-0.25*(F_0c[0]*alpha_c[13]+alpha_c[0]*F_0c[13]+F_0c[6]*alpha_c[11]+alpha_c[6]*F_0c[11])-0.223606797749979*F_0c[9]*alpha_c[13]-0.223606797749979*(F_0c[4]*alpha_c[10]+alpha_c[4]*F_0c[10])-0.2500000000000001*(F_0c[3]*alpha_c[7]+alpha_c[3]*F_0c[7])+0.5*(Ghat_F_0_r[6]+Ghat_F_0_l[6])-0.223606797749979*(F_0c[1]*alpha_c[5]+alpha_c[1]*F_0c[5]))*dvpar); 
  out_pressure[7] += dx1*volFact*(((-1.732050807568877*(F_0c[6]*alpha_c[18]+alpha_c[6]*F_0c[18]+F_0c[5]*alpha_c[17]+alpha_c[5]*F_0c[17]))-1.549193338482967*(F_0c[17]*alpha_c[18]+alpha_c[17]*F_0c[18])-1.732050807568877*(F_0c[10]*alpha_c[14]+alpha_c[10]*F_0c[14]+F_0c[10]*alpha_c[13]+alpha_c[10]*F_0c[13])-1.549193338482967*(F_0c[11]*alpha_c[12]+alpha_c[11]*F_0c[12])-1.732050807568877*(F_0c[2]*alpha_c[12]+alpha_c[2]*F_0c[12]+F_0c[1]*alpha_c[11]+alpha_c[1]*F_0c[11])-1.936491673103709*(F_0c[3]*alpha_c[10]+alpha_c[3]*F_0c[10])-1.732050807568877*(F_0c[4]*alpha_c[8]+alpha_c[4]*F_0c[8]+F_0c[4]*alpha_c[7]+alpha_c[4]*F_0c[7])-1.936491673103709*(F_0c[5]*alpha_c[6]+alpha_c[5]*F_0c[6]+F_0c[0]*alpha_c[4]+alpha_c[0]*F_0c[4]+F_0c[1]*alpha_c[2]+alpha_c[1]*F_0c[2])+2.23606797749979*Ghat_F_0_r[1]-2.23606797749979*Ghat_F_0_l[1])*wvpar+(((-0.4472135954999579*(alpha_c[14]+alpha_c[13]))-0.5000000000000001*alpha_c[3])*F_0c[19]+((-0.4472135954999579*(F_0c[16]+F_0c[11]))-0.5000000000000001*F_0c[2])*alpha_c[18]+((-0.4472135954999579*alpha_c[11])-0.5000000000000001*alpha_c[2])*F_0c[18]+((-0.4472135954999579*(F_0c[15]+F_0c[12]))-0.5000000000000001*F_0c[1])*alpha_c[17]+((-0.4472135954999579*alpha_c[12])-0.5000000000000001*alpha_c[1])*F_0c[17]-0.5*(alpha_c[5]*F_0c[16]+alpha_c[6]*F_0c[15]+F_0c[4]*alpha_c[14]+alpha_c[4]*F_0c[14]+F_0c[4]*alpha_c[13]+alpha_c[4]*F_0c[13]+F_0c[6]*alpha_c[12]+alpha_c[6]*F_0c[12]+F_0c[5]*alpha_c[11]+alpha_c[5]*F_0c[11])-0.5000000000000001*(F_0c[7]*alpha_c[10]+(alpha_c[8]+alpha_c[7])*F_0c[10])-0.5590169943749476*(F_0c[0]*alpha_c[10]+alpha_c[0]*F_0c[10]+F_0c[1]*alpha_c[6]+alpha_c[1]*F_0c[6]+F_0c[2]*alpha_c[5]+alpha_c[2]*F_0c[5]+F_0c[3]*alpha_c[4]+alpha_c[3]*F_0c[4])-0.5000000000000001*(F_0c[9]+F_0c[8])*alpha_c[10]+0.6454972243679028*Ghat_F_0_r[3]-0.6454972243679028*Ghat_F_0_l[3])*dvpar); 

} 
