#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_tensor_2x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_tensor_2x_p2_upwind_quad_to_modal.h> 
#include <gkyl_basis_ser_1x_p2_surfx1_eval_quad.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_1x1v_tensor_p2(const double *w, const double *dxv, 
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
  const double *G_1l = &fl[9]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[9]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[9]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[9]; 
  double alpha_l[9] = {0.0}; 
  double alpha_c[9] = {0.0}; 
  double alpha_r[9] = {0.0}; 
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
    F_0_UpwindQuad_l[0] = tensor_2x_p2_surfx1_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = tensor_2x_p2_surfx1_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = tensor_2x_p2_surfx1_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = tensor_2x_p2_surfx1_eval_quad_node_0_l(G_1c); 
  } 
  if (0.7071067811865475*alphaSurf_r[0]-0.9486832980505137*alphaSurf_r[1] > 0) { 
    F_0_UpwindQuad_r[0] = tensor_2x_p2_surfx1_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = tensor_2x_p2_surfx1_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = tensor_2x_p2_surfx1_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = tensor_2x_p2_surfx1_eval_quad_node_0_l(G_1r); 
  } 
  if (0.7071067811865475*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = tensor_2x_p2_surfx1_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = tensor_2x_p2_surfx1_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = tensor_2x_p2_surfx1_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = tensor_2x_p2_surfx1_eval_quad_node_1_l(G_1c); 
  } 
  if (0.7071067811865475*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = tensor_2x_p2_surfx1_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = tensor_2x_p2_surfx1_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = tensor_2x_p2_surfx1_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = tensor_2x_p2_surfx1_eval_quad_node_1_l(G_1r); 
  } 
  if (0.9486832980505137*alphaSurf_l[1]+0.7071067811865475*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = tensor_2x_p2_surfx1_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = tensor_2x_p2_surfx1_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = tensor_2x_p2_surfx1_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = tensor_2x_p2_surfx1_eval_quad_node_2_l(G_1c); 
  } 
  if (0.9486832980505137*alphaSurf_r[1]+0.7071067811865475*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = tensor_2x_p2_surfx1_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = tensor_2x_p2_surfx1_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = tensor_2x_p2_surfx1_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = tensor_2x_p2_surfx1_eval_quad_node_2_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  tensor_2x_p2_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  tensor_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  tensor_2x_p2_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

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
  out_F_0[8] += (1.58113883008419*Ghat_F_0_l[2]-1.58113883008419*Ghat_F_0_r[2])*dx1; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_l[0]-0.7071067811865475*Ghat_G_1_r[0])*dx1; 
  out_G_1[1] += -1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dx1; 
  out_G_1[2] += (0.7071067811865475*Ghat_G_1_l[1]-0.7071067811865475*Ghat_G_1_r[1])*dx1; 
  out_G_1[3] += -1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dx1; 
  out_G_1[4] += (1.58113883008419*Ghat_G_1_l[0]-1.58113883008419*Ghat_G_1_r[0])*dx1; 
  out_G_1[5] += (0.7071067811865475*Ghat_G_1_l[2]-0.7071067811865475*Ghat_G_1_r[2])*dx1; 
  out_G_1[6] += (1.58113883008419*Ghat_G_1_l[1]-1.58113883008419*Ghat_G_1_r[1])*dx1; 
  out_G_1[7] += -1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dx1; 
  out_G_1[8] += (1.58113883008419*Ghat_G_1_l[2]-1.58113883008419*Ghat_G_1_r[2])*dx1; 

  double F_0l_rx = 0.0; 
  double F_0c_lx = 0.0; 
  double F_0c_rx = 0.0; 
  double F_0r_lx = 0.0; 

  double G_1l_rx = 0.0; 
  double G_1c_lx = 0.0; 
  double G_1c_rx = 0.0; 
  double G_1r_lx = 0.0; 

  double ul_r = ser_1x_p2_surfx1_eval_quad_node_0_r(ul); 
  double uc_l = ser_1x_p2_surfx1_eval_quad_node_0_l(uc); 
  double uc_r = ser_1x_p2_surfx1_eval_quad_node_0_r(uc); 
  double ur_l = ser_1x_p2_surfx1_eval_quad_node_0_l(ur); 

  double vth_sql_r = ser_1x_p2_surfx1_eval_quad_node_0_r(vth_sql); 
  double vth_sqc_l = ser_1x_p2_surfx1_eval_quad_node_0_l(vth_sqc); 
  double vth_sqc_r = ser_1x_p2_surfx1_eval_quad_node_0_r(vth_sqc); 
  double vth_sqr_l = ser_1x_p2_surfx1_eval_quad_node_0_l(vth_sqr); 

  double avg_u_l = 0.5*(ul_r + uc_l); 
  double avg_u_r = 0.5*(uc_r + ur_l); 
  double max_u_l = fmax(fabs(ul_r), fabs(uc_l)); 
  double max_u_r = fmax(fabs(uc_r), fabs(ur_l)); 
  double max_vth_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  double max_vth_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  double max_speed_l = max_u_l + max_vth_l; 
  double max_speed_r = max_u_r + max_vth_r; 

  F_0l_rx = (-3.122502256758253e-16*F_0l[8])+1.118033988749895*F_0l[4]+0.8660254037844386*F_0l[1]+0.5*F_0l[0]; 
  F_0c_lx = (-3.122502256758253e-16*F_0c[8])+1.118033988749895*F_0c[4]-0.8660254037844386*F_0c[1]+0.5*F_0c[0]; 
  F_0c_rx = (-3.122502256758253e-16*F_0c[8])+1.118033988749895*F_0c[4]+0.8660254037844386*F_0c[1]+0.5*F_0c[0]; 
  F_0r_lx = (-3.122502256758253e-16*F_0r[8])+1.118033988749895*F_0r[4]-0.8660254037844386*F_0r[1]+0.5*F_0r[0]; 
  G_1l_rx = (-3.122502256758253e-16*G_1l[8])+1.118033988749895*G_1l[4]+0.8660254037844386*G_1l[1]+0.5*G_1l[0]; 
  G_1c_lx = (-3.122502256758253e-16*G_1c[8])+1.118033988749895*G_1c[4]-0.8660254037844386*G_1c[1]+0.5*G_1c[0]; 
  G_1c_rx = (-3.122502256758253e-16*G_1c[8])+1.118033988749895*G_1c[4]+0.8660254037844386*G_1c[1]+0.5*G_1c[0]; 
  G_1r_lx = (-3.122502256758253e-16*G_1r[8])+1.118033988749895*G_1r[4]-0.8660254037844386*G_1r[1]+0.5*G_1r[0]; 

  out_F_0[0] += 0.5*dx1*((F_0l_rx + F_0c_lx)*avg_u_l - (F_0c_rx + F_0r_lx)*avg_u_r - max_speed_l*(F_0c_lx - F_0l_rx) + max_speed_r*(F_0r_lx - F_0c_rx)); 
  out_G_1[0] += 0.5*dx1*((G_1l_rx + G_1c_lx)*avg_u_l - (G_1c_rx + G_1r_lx)*avg_u_r - max_speed_l*(G_1c_lx - G_1l_rx) + max_speed_r*(G_1r_lx - G_1c_rx)); 

  F_0l_rx = 5.408332555453774e-16*F_0l[8]-1.936491673103709*F_0l[4]-1.5*F_0l[1]-0.8660254037844386*F_0l[0]; 
  F_0c_lx = 5.408332555453774e-16*F_0c[8]-1.936491673103709*F_0c[4]+1.5*F_0c[1]-0.8660254037844386*F_0c[0]; 
  F_0c_rx = (-5.408332555453774e-16*F_0c[8])+1.936491673103709*F_0c[4]+1.5*F_0c[1]+0.8660254037844386*F_0c[0]; 
  F_0r_lx = (-5.408332555453774e-16*F_0r[8])+1.936491673103709*F_0r[4]-1.5*F_0r[1]+0.8660254037844386*F_0r[0]; 
  G_1l_rx = 5.408332555453774e-16*G_1l[8]-1.936491673103709*G_1l[4]-1.5*G_1l[1]-0.8660254037844386*G_1l[0]; 
  G_1c_lx = 5.408332555453774e-16*G_1c[8]-1.936491673103709*G_1c[4]+1.5*G_1c[1]-0.8660254037844386*G_1c[0]; 
  G_1c_rx = (-5.408332555453774e-16*G_1c[8])+1.936491673103709*G_1c[4]+1.5*G_1c[1]+0.8660254037844386*G_1c[0]; 
  G_1r_lx = (-5.408332555453774e-16*G_1r[8])+1.936491673103709*G_1r[4]-1.5*G_1r[1]+0.8660254037844386*G_1r[0]; 

  out_F_0[1] += 0.5*dx1*((F_0l_rx + F_0c_lx)*avg_u_l - (F_0c_rx + F_0r_lx)*avg_u_r - max_speed_l*(F_0c_lx - F_0l_rx) + max_speed_r*(F_0r_lx - F_0c_rx)); 
  out_G_1[1] += 0.5*dx1*((G_1l_rx + G_1c_lx)*avg_u_l - (G_1c_rx + G_1r_lx)*avg_u_r - max_speed_l*(G_1c_lx - G_1l_rx) + max_speed_r*(G_1r_lx - G_1c_rx)); 

  F_0l_rx = 1.118033988749895*F_0l[6]+0.8660254037844386*F_0l[3]+0.5*F_0l[2]; 
  F_0c_lx = 1.118033988749895*F_0c[6]-0.8660254037844386*F_0c[3]+0.5*F_0c[2]; 
  F_0c_rx = 1.118033988749895*F_0c[6]+0.8660254037844386*F_0c[3]+0.5*F_0c[2]; 
  F_0r_lx = 1.118033988749895*F_0r[6]-0.8660254037844386*F_0r[3]+0.5*F_0r[2]; 
  G_1l_rx = 1.118033988749895*G_1l[6]+0.8660254037844386*G_1l[3]+0.5*G_1l[2]; 
  G_1c_lx = 1.118033988749895*G_1c[6]-0.8660254037844386*G_1c[3]+0.5*G_1c[2]; 
  G_1c_rx = 1.118033988749895*G_1c[6]+0.8660254037844386*G_1c[3]+0.5*G_1c[2]; 
  G_1r_lx = 1.118033988749895*G_1r[6]-0.8660254037844386*G_1r[3]+0.5*G_1r[2]; 

  out_F_0[2] += 0.5*dx1*((F_0l_rx + F_0c_lx)*avg_u_l - (F_0c_rx + F_0r_lx)*avg_u_r - max_speed_l*(F_0c_lx - F_0l_rx) + max_speed_r*(F_0r_lx - F_0c_rx)); 
  out_G_1[2] += 0.5*dx1*((G_1l_rx + G_1c_lx)*avg_u_l - (G_1c_rx + G_1r_lx)*avg_u_r - max_speed_l*(G_1c_lx - G_1l_rx) + max_speed_r*(G_1r_lx - G_1c_rx)); 

  F_0l_rx = (-1.936491673103709*F_0l[6])-1.5*F_0l[3]-0.8660254037844386*F_0l[2]; 
  F_0c_lx = (-1.936491673103709*F_0c[6])+1.5*F_0c[3]-0.8660254037844386*F_0c[2]; 
  F_0c_rx = 1.936491673103709*F_0c[6]+1.5*F_0c[3]+0.8660254037844386*F_0c[2]; 
  F_0r_lx = 1.936491673103709*F_0r[6]-1.5*F_0r[3]+0.8660254037844386*F_0r[2]; 
  G_1l_rx = (-1.936491673103709*G_1l[6])-1.5*G_1l[3]-0.8660254037844386*G_1l[2]; 
  G_1c_lx = (-1.936491673103709*G_1c[6])+1.5*G_1c[3]-0.8660254037844386*G_1c[2]; 
  G_1c_rx = 1.936491673103709*G_1c[6]+1.5*G_1c[3]+0.8660254037844386*G_1c[2]; 
  G_1r_lx = 1.936491673103709*G_1r[6]-1.5*G_1r[3]+0.8660254037844386*G_1r[2]; 

  out_F_0[3] += 0.5*dx1*((F_0l_rx + F_0c_lx)*avg_u_l - (F_0c_rx + F_0r_lx)*avg_u_r - max_speed_l*(F_0c_lx - F_0l_rx) + max_speed_r*(F_0r_lx - F_0c_rx)); 
  out_G_1[3] += 0.5*dx1*((G_1l_rx + G_1c_lx)*avg_u_l - (G_1c_rx + G_1r_lx)*avg_u_r - max_speed_l*(G_1c_lx - G_1l_rx) + max_speed_r*(G_1r_lx - G_1c_rx)); 

  F_0l_rx = (-6.982127306007956e-16*F_0l[8])+2.5*F_0l[4]+1.936491673103709*F_0l[1]+1.118033988749895*F_0l[0]; 
  F_0c_lx = (-6.982127306007956e-16*F_0c[8])+2.5*F_0c[4]-1.936491673103709*F_0c[1]+1.118033988749895*F_0c[0]; 
  F_0c_rx = (-6.982127306007956e-16*F_0c[8])+2.5*F_0c[4]+1.936491673103709*F_0c[1]+1.118033988749895*F_0c[0]; 
  F_0r_lx = (-6.982127306007956e-16*F_0r[8])+2.5*F_0r[4]-1.936491673103709*F_0r[1]+1.118033988749895*F_0r[0]; 
  G_1l_rx = (-6.982127306007956e-16*G_1l[8])+2.5*G_1l[4]+1.936491673103709*G_1l[1]+1.118033988749895*G_1l[0]; 
  G_1c_lx = (-6.982127306007956e-16*G_1c[8])+2.5*G_1c[4]-1.936491673103709*G_1c[1]+1.118033988749895*G_1c[0]; 
  G_1c_rx = (-6.982127306007956e-16*G_1c[8])+2.5*G_1c[4]+1.936491673103709*G_1c[1]+1.118033988749895*G_1c[0]; 
  G_1r_lx = (-6.982127306007956e-16*G_1r[8])+2.5*G_1r[4]-1.936491673103709*G_1r[1]+1.118033988749895*G_1r[0]; 

  out_F_0[4] += 0.5*dx1*((F_0l_rx + F_0c_lx)*avg_u_l - (F_0c_rx + F_0r_lx)*avg_u_r - max_speed_l*(F_0c_lx - F_0l_rx) + max_speed_r*(F_0r_lx - F_0c_rx)); 
  out_G_1[4] += 0.5*dx1*((G_1l_rx + G_1c_lx)*avg_u_l - (G_1c_rx + G_1r_lx)*avg_u_r - max_speed_l*(G_1c_lx - G_1l_rx) + max_speed_r*(G_1r_lx - G_1c_rx)); 

  F_0l_rx = 1.118033988749895*F_0l[8]+0.8660254037844387*F_0l[7]+0.5*F_0l[5]; 
  F_0c_lx = 1.118033988749895*F_0c[8]-0.8660254037844387*F_0c[7]+0.5*F_0c[5]; 
  F_0c_rx = 1.118033988749895*F_0c[8]+0.8660254037844387*F_0c[7]+0.5*F_0c[5]; 
  F_0r_lx = 1.118033988749895*F_0r[8]-0.8660254037844387*F_0r[7]+0.5*F_0r[5]; 
  G_1l_rx = 1.118033988749895*G_1l[8]+0.8660254037844387*G_1l[7]+0.5*G_1l[5]; 
  G_1c_lx = 1.118033988749895*G_1c[8]-0.8660254037844387*G_1c[7]+0.5*G_1c[5]; 
  G_1c_rx = 1.118033988749895*G_1c[8]+0.8660254037844387*G_1c[7]+0.5*G_1c[5]; 
  G_1r_lx = 1.118033988749895*G_1r[8]-0.8660254037844387*G_1r[7]+0.5*G_1r[5]; 

  out_F_0[5] += 0.5*dx1*((F_0l_rx + F_0c_lx)*avg_u_l - (F_0c_rx + F_0r_lx)*avg_u_r - max_speed_l*(F_0c_lx - F_0l_rx) + max_speed_r*(F_0r_lx - F_0c_rx)); 
  out_G_1[5] += 0.5*dx1*((G_1l_rx + G_1c_lx)*avg_u_l - (G_1c_rx + G_1r_lx)*avg_u_r - max_speed_l*(G_1c_lx - G_1l_rx) + max_speed_r*(G_1r_lx - G_1c_rx)); 

  F_0l_rx = 2.5*F_0l[6]+1.936491673103709*F_0l[3]+1.118033988749895*F_0l[2]; 
  F_0c_lx = 2.5*F_0c[6]-1.936491673103709*F_0c[3]+1.118033988749895*F_0c[2]; 
  F_0c_rx = 2.5*F_0c[6]+1.936491673103709*F_0c[3]+1.118033988749895*F_0c[2]; 
  F_0r_lx = 2.5*F_0r[6]-1.936491673103709*F_0r[3]+1.118033988749895*F_0r[2]; 
  G_1l_rx = 2.5*G_1l[6]+1.936491673103709*G_1l[3]+1.118033988749895*G_1l[2]; 
  G_1c_lx = 2.5*G_1c[6]-1.936491673103709*G_1c[3]+1.118033988749895*G_1c[2]; 
  G_1c_rx = 2.5*G_1c[6]+1.936491673103709*G_1c[3]+1.118033988749895*G_1c[2]; 
  G_1r_lx = 2.5*G_1r[6]-1.936491673103709*G_1r[3]+1.118033988749895*G_1r[2]; 

  out_F_0[6] += 0.5*dx1*((F_0l_rx + F_0c_lx)*avg_u_l - (F_0c_rx + F_0r_lx)*avg_u_r - max_speed_l*(F_0c_lx - F_0l_rx) + max_speed_r*(F_0r_lx - F_0c_rx)); 
  out_G_1[6] += 0.5*dx1*((G_1l_rx + G_1c_lx)*avg_u_l - (G_1c_rx + G_1r_lx)*avg_u_r - max_speed_l*(G_1c_lx - G_1l_rx) + max_speed_r*(G_1r_lx - G_1c_rx)); 

  F_0l_rx = (-1.936491673103709*F_0l[8])-1.5*F_0l[7]-0.8660254037844387*F_0l[5]; 
  F_0c_lx = (-1.936491673103709*F_0c[8])+1.5*F_0c[7]-0.8660254037844387*F_0c[5]; 
  F_0c_rx = 1.936491673103709*F_0c[8]+1.5*F_0c[7]+0.8660254037844387*F_0c[5]; 
  F_0r_lx = 1.936491673103709*F_0r[8]-1.5*F_0r[7]+0.8660254037844387*F_0r[5]; 
  G_1l_rx = (-1.936491673103709*G_1l[8])-1.5*G_1l[7]-0.8660254037844387*G_1l[5]; 
  G_1c_lx = (-1.936491673103709*G_1c[8])+1.5*G_1c[7]-0.8660254037844387*G_1c[5]; 
  G_1c_rx = 1.936491673103709*G_1c[8]+1.5*G_1c[7]+0.8660254037844387*G_1c[5]; 
  G_1r_lx = 1.936491673103709*G_1r[8]-1.5*G_1r[7]+0.8660254037844387*G_1r[5]; 

  out_F_0[7] += 0.5*dx1*((F_0l_rx + F_0c_lx)*avg_u_l - (F_0c_rx + F_0r_lx)*avg_u_r - max_speed_l*(F_0c_lx - F_0l_rx) + max_speed_r*(F_0r_lx - F_0c_rx)); 
  out_G_1[7] += 0.5*dx1*((G_1l_rx + G_1c_lx)*avg_u_l - (G_1c_rx + G_1r_lx)*avg_u_r - max_speed_l*(G_1c_lx - G_1l_rx) + max_speed_r*(G_1r_lx - G_1c_rx)); 

  F_0l_rx = 2.5*F_0l[8]+1.936491673103709*F_0l[7]+1.118033988749895*F_0l[5]; 
  F_0c_lx = 2.5*F_0c[8]-1.936491673103709*F_0c[7]+1.118033988749895*F_0c[5]; 
  F_0c_rx = 2.5*F_0c[8]+1.936491673103709*F_0c[7]+1.118033988749895*F_0c[5]; 
  F_0r_lx = 2.5*F_0r[8]-1.936491673103709*F_0r[7]+1.118033988749895*F_0r[5]; 
  G_1l_rx = 2.5*G_1l[8]+1.936491673103709*G_1l[7]+1.118033988749895*G_1l[5]; 
  G_1c_lx = 2.5*G_1c[8]-1.936491673103709*G_1c[7]+1.118033988749895*G_1c[5]; 
  G_1c_rx = 2.5*G_1c[8]+1.936491673103709*G_1c[7]+1.118033988749895*G_1c[5]; 
  G_1r_lx = 2.5*G_1r[8]-1.936491673103709*G_1r[7]+1.118033988749895*G_1r[5]; 

  out_F_0[8] += 0.5*dx1*((F_0l_rx + F_0c_lx)*avg_u_l - (F_0c_rx + F_0r_lx)*avg_u_r - max_speed_l*(F_0c_lx - F_0l_rx) + max_speed_r*(F_0r_lx - F_0c_rx)); 
  out_G_1[8] += 0.5*dx1*((G_1l_rx + G_1c_lx)*avg_u_l - (G_1c_rx + G_1r_lx)*avg_u_r - max_speed_l*(G_1c_lx - G_1l_rx) + max_speed_r*(G_1r_lx - G_1c_rx)); 

} 
