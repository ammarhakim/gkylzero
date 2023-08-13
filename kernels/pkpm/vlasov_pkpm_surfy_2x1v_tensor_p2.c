#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_tensor_3x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_tensor_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfy_2x1v_tensor_p2(const double *w, const double *dxv, 
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
  const double dx1 = 2.0/dxv[1]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *bl = &bvar_l[9]; 
  const double *bc = &bvar_c[9]; 
  const double *br = &bvar_r[9]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[27]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[27]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[27]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[27]; 
  double alpha_l[27] = {0.0}; 
  double alpha_c[27] = {0.0}; 
  double alpha_r[27] = {0.0}; 
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
  alpha_l[20] = 1.414213562373095*bl[8]*wvpar; 
  alpha_l[23] = 0.408248290463863*bl[8]*dvpar; 

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
  alpha_c[20] = 1.414213562373095*bc[8]*wvpar; 
  alpha_c[23] = 0.408248290463863*bc[8]*dvpar; 

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
  alpha_r[20] = 1.414213562373095*br[8]*wvpar; 
  alpha_r[23] = 0.408248290463863*br[8]*dvpar; 

  double alphaSurf_l[9] = {0.0}; 
  alphaSurf_l[0] = 0.3458741190809163*alpha_l[8]+0.3458741190809163*alpha_c[8]+0.4975526040028326*alpha_l[2]-0.4975526040028326*alpha_c[2]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.3458741190809163*alpha_l[12]+0.3458741190809163*alpha_c[12]+0.4975526040028326*alpha_l[4]-0.4975526040028326*alpha_c[4]+0.3535533905932737*alpha_l[1]+0.3535533905932737*alpha_c[1]; 
  alphaSurf_l[2] = 0.3458741190809163*alpha_l[14]+0.3458741190809163*alpha_c[14]+0.4975526040028326*alpha_l[6]-0.4975526040028326*alpha_c[6]+0.3535533905932737*alpha_l[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_l[3] = 0.3458741190809163*alpha_l[18]+0.3458741190809163*alpha_c[18]+0.4975526040028326*alpha_l[10]-0.4975526040028326*alpha_c[10]+0.3535533905932737*alpha_l[5]+0.3535533905932737*alpha_c[5]; 
  alphaSurf_l[4] = 0.3458741190809163*alpha_l[20]+0.3458741190809163*alpha_c[20]+0.4975526040028326*alpha_l[11]-0.4975526040028326*alpha_c[11]+0.3535533905932737*alpha_l[7]+0.3535533905932737*alpha_c[7]; 
  alphaSurf_l[6] = 0.3458741190809163*alpha_l[23]+0.3458741190809163*alpha_c[23]+0.4975526040028326*alpha_l[17]-0.4975526040028326*alpha_c[17]+0.3535533905932737*alpha_l[13]+0.3535533905932737*alpha_c[13]; 

  double alphaSurf_r[9] = {0.0}; 
  alphaSurf_r[0] = 0.3458741190809163*alpha_r[8]+0.3458741190809163*alpha_c[8]-0.4975526040028326*alpha_r[2]+0.4975526040028326*alpha_c[2]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = 0.3458741190809163*alpha_r[12]+0.3458741190809163*alpha_c[12]-0.4975526040028326*alpha_r[4]+0.4975526040028326*alpha_c[4]+0.3535533905932737*alpha_r[1]+0.3535533905932737*alpha_c[1]; 
  alphaSurf_r[2] = 0.3458741190809163*alpha_r[14]+0.3458741190809163*alpha_c[14]-0.4975526040028326*alpha_r[6]+0.4975526040028326*alpha_c[6]+0.3535533905932737*alpha_r[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_r[3] = 0.3458741190809163*alpha_r[18]+0.3458741190809163*alpha_c[18]-0.4975526040028326*alpha_r[10]+0.4975526040028326*alpha_c[10]+0.3535533905932737*alpha_r[5]+0.3535533905932737*alpha_c[5]; 
  alphaSurf_r[4] = 0.3458741190809163*alpha_r[20]+0.3458741190809163*alpha_c[20]-0.4975526040028326*alpha_r[11]+0.4975526040028326*alpha_c[11]+0.3535533905932737*alpha_r[7]+0.3535533905932737*alpha_c[7]; 
  alphaSurf_r[6] = 0.3458741190809163*alpha_r[23]+0.3458741190809163*alpha_c[23]-0.4975526040028326*alpha_r[17]+0.4975526040028326*alpha_c[17]+0.3535533905932737*alpha_r[13]+0.3535533905932737*alpha_c[13]; 

  double F_0_UpwindQuad_l[9] = {0.0};
  double F_0_UpwindQuad_r[9] = {0.0};
  double F_0_Upwind_l[9] = {0.0};
  double F_0_Upwind_r[9] = {0.0};
  double Ghat_F_0_l[9] = {0.0}; 
  double Ghat_F_0_r[9] = {0.0}; 
  double G_1_UpwindQuad_l[9] = {0.0};
  double G_1_UpwindQuad_r[9] = {0.0};
  double G_1_Upwind_l[9] = {0.0};
  double G_1_Upwind_r[9] = {0.0};
  double Ghat_G_1_l[9] = {0.0}; 
  double Ghat_G_1_r[9] = {0.0}; 

  if ((-0.5999999999999999*alphaSurf_l[6])+0.4472135954999579*alphaSurf_l[4]+0.9*alphaSurf_l[3]-0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = tensor_3x_p2_surfx2_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = tensor_3x_p2_surfx2_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = tensor_3x_p2_surfx2_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = tensor_3x_p2_surfx2_eval_quad_node_0_l(G_1c); 
  } 
  if ((-0.5999999999999999*alphaSurf_r[6])+0.4472135954999579*alphaSurf_r[4]+0.9*alphaSurf_r[3]-0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = tensor_3x_p2_surfx2_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = tensor_3x_p2_surfx2_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = tensor_3x_p2_surfx2_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = tensor_3x_p2_surfx2_eval_quad_node_0_l(G_1r); 
  } 
  if (0.4472135954999579*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = tensor_3x_p2_surfx2_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = tensor_3x_p2_surfx2_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = tensor_3x_p2_surfx2_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = tensor_3x_p2_surfx2_eval_quad_node_1_l(G_1c); 
  } 
  if (0.4472135954999579*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = tensor_3x_p2_surfx2_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = tensor_3x_p2_surfx2_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = tensor_3x_p2_surfx2_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = tensor_3x_p2_surfx2_eval_quad_node_1_l(G_1r); 
  } 
  if (0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[4]-0.9*alphaSurf_l[3]+0.6708203932499369*alphaSurf_l[2]-0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = tensor_3x_p2_surfx2_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = tensor_3x_p2_surfx2_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = tensor_3x_p2_surfx2_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = tensor_3x_p2_surfx2_eval_quad_node_2_l(G_1c); 
  } 
  if (0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[4]-0.9*alphaSurf_r[3]+0.6708203932499369*alphaSurf_r[2]-0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = tensor_3x_p2_surfx2_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = tensor_3x_p2_surfx2_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = tensor_3x_p2_surfx2_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = tensor_3x_p2_surfx2_eval_quad_node_2_l(G_1r); 
  } 
  if (0.75*alphaSurf_l[6]-0.5590169943749475*alphaSurf_l[4]-0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[3] = tensor_3x_p2_surfx2_eval_quad_node_3_r(F_0l); 
    G_1_UpwindQuad_l[3] = tensor_3x_p2_surfx2_eval_quad_node_3_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[3] = tensor_3x_p2_surfx2_eval_quad_node_3_l(F_0c); 
    G_1_UpwindQuad_l[3] = tensor_3x_p2_surfx2_eval_quad_node_3_l(G_1c); 
  } 
  if (0.75*alphaSurf_r[6]-0.5590169943749475*alphaSurf_r[4]-0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[3] = tensor_3x_p2_surfx2_eval_quad_node_3_r(F_0c); 
    G_1_UpwindQuad_r[3] = tensor_3x_p2_surfx2_eval_quad_node_3_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[3] = tensor_3x_p2_surfx2_eval_quad_node_3_l(F_0r); 
    G_1_UpwindQuad_r[3] = tensor_3x_p2_surfx2_eval_quad_node_3_l(G_1r); 
  } 
  if (0.5*alphaSurf_l[0]-0.5590169943749475*alphaSurf_l[4] > 0) { 
    F_0_UpwindQuad_l[4] = tensor_3x_p2_surfx2_eval_quad_node_4_r(F_0l); 
    G_1_UpwindQuad_l[4] = tensor_3x_p2_surfx2_eval_quad_node_4_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[4] = tensor_3x_p2_surfx2_eval_quad_node_4_l(F_0c); 
    G_1_UpwindQuad_l[4] = tensor_3x_p2_surfx2_eval_quad_node_4_l(G_1c); 
  } 
  if (0.5*alphaSurf_r[0]-0.5590169943749475*alphaSurf_r[4] > 0) { 
    F_0_UpwindQuad_r[4] = tensor_3x_p2_surfx2_eval_quad_node_4_r(F_0c); 
    G_1_UpwindQuad_r[4] = tensor_3x_p2_surfx2_eval_quad_node_4_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[4] = tensor_3x_p2_surfx2_eval_quad_node_4_l(F_0r); 
    G_1_UpwindQuad_r[4] = tensor_3x_p2_surfx2_eval_quad_node_4_l(G_1r); 
  } 
  if ((-0.75*alphaSurf_l[6])-0.5590169943749475*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[2]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[5] = tensor_3x_p2_surfx2_eval_quad_node_5_r(F_0l); 
    G_1_UpwindQuad_l[5] = tensor_3x_p2_surfx2_eval_quad_node_5_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[5] = tensor_3x_p2_surfx2_eval_quad_node_5_l(F_0c); 
    G_1_UpwindQuad_l[5] = tensor_3x_p2_surfx2_eval_quad_node_5_l(G_1c); 
  } 
  if ((-0.75*alphaSurf_r[6])-0.5590169943749475*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[2]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[5] = tensor_3x_p2_surfx2_eval_quad_node_5_r(F_0c); 
    G_1_UpwindQuad_r[5] = tensor_3x_p2_surfx2_eval_quad_node_5_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[5] = tensor_3x_p2_surfx2_eval_quad_node_5_l(F_0r); 
    G_1_UpwindQuad_r[5] = tensor_3x_p2_surfx2_eval_quad_node_5_l(G_1r); 
  } 
  if ((-0.5999999999999999*alphaSurf_l[6])+0.4472135954999579*alphaSurf_l[4]-0.9*alphaSurf_l[3]-0.6708203932499369*alphaSurf_l[2]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[6] = tensor_3x_p2_surfx2_eval_quad_node_6_r(F_0l); 
    G_1_UpwindQuad_l[6] = tensor_3x_p2_surfx2_eval_quad_node_6_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[6] = tensor_3x_p2_surfx2_eval_quad_node_6_l(F_0c); 
    G_1_UpwindQuad_l[6] = tensor_3x_p2_surfx2_eval_quad_node_6_l(G_1c); 
  } 
  if ((-0.5999999999999999*alphaSurf_r[6])+0.4472135954999579*alphaSurf_r[4]-0.9*alphaSurf_r[3]-0.6708203932499369*alphaSurf_r[2]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[6] = tensor_3x_p2_surfx2_eval_quad_node_6_r(F_0c); 
    G_1_UpwindQuad_r[6] = tensor_3x_p2_surfx2_eval_quad_node_6_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[6] = tensor_3x_p2_surfx2_eval_quad_node_6_l(F_0r); 
    G_1_UpwindQuad_r[6] = tensor_3x_p2_surfx2_eval_quad_node_6_l(G_1r); 
  } 
  if (0.4472135954999579*alphaSurf_l[4]+0.6708203932499369*alphaSurf_l[1]+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[7] = tensor_3x_p2_surfx2_eval_quad_node_7_r(F_0l); 
    G_1_UpwindQuad_l[7] = tensor_3x_p2_surfx2_eval_quad_node_7_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[7] = tensor_3x_p2_surfx2_eval_quad_node_7_l(F_0c); 
    G_1_UpwindQuad_l[7] = tensor_3x_p2_surfx2_eval_quad_node_7_l(G_1c); 
  } 
  if (0.4472135954999579*alphaSurf_r[4]+0.6708203932499369*alphaSurf_r[1]+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[7] = tensor_3x_p2_surfx2_eval_quad_node_7_r(F_0c); 
    G_1_UpwindQuad_r[7] = tensor_3x_p2_surfx2_eval_quad_node_7_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[7] = tensor_3x_p2_surfx2_eval_quad_node_7_l(F_0r); 
    G_1_UpwindQuad_r[7] = tensor_3x_p2_surfx2_eval_quad_node_7_l(G_1r); 
  } 
  if (0.5999999999999999*alphaSurf_l[6]+0.4472135954999579*alphaSurf_l[4]+0.9*alphaSurf_l[3]+0.6708203932499369*(alphaSurf_l[2]+alphaSurf_l[1])+0.5*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[8] = tensor_3x_p2_surfx2_eval_quad_node_8_r(F_0l); 
    G_1_UpwindQuad_l[8] = tensor_3x_p2_surfx2_eval_quad_node_8_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[8] = tensor_3x_p2_surfx2_eval_quad_node_8_l(F_0c); 
    G_1_UpwindQuad_l[8] = tensor_3x_p2_surfx2_eval_quad_node_8_l(G_1c); 
  } 
  if (0.5999999999999999*alphaSurf_r[6]+0.4472135954999579*alphaSurf_r[4]+0.9*alphaSurf_r[3]+0.6708203932499369*(alphaSurf_r[2]+alphaSurf_r[1])+0.5*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[8] = tensor_3x_p2_surfx2_eval_quad_node_8_r(F_0c); 
    G_1_UpwindQuad_r[8] = tensor_3x_p2_surfx2_eval_quad_node_8_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[8] = tensor_3x_p2_surfx2_eval_quad_node_8_l(F_0r); 
    G_1_UpwindQuad_r[8] = tensor_3x_p2_surfx2_eval_quad_node_8_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  tensor_3x_p2_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  tensor_3x_p2_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  Ghat_F_0_l[0] = 0.5*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5*F_0_Upwind_l[4]*alphaSurf_l[4]+0.5*F_0_Upwind_l[3]*alphaSurf_l[3]+0.5*F_0_Upwind_l[2]*alphaSurf_l[2]+0.5*F_0_Upwind_l[1]*alphaSurf_l[1]+0.5*F_0_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_F_0_l[1] = 0.447213595499958*F_0_Upwind_l[3]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[6]+0.4472135954999579*F_0_Upwind_l[1]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[1]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.447213595499958*alphaSurf_l[6]*F_0_Upwind_l[8]+0.447213595499958*alphaSurf_l[3]*F_0_Upwind_l[7]+0.5000000000000001*F_0_Upwind_l[4]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[4]*F_0_Upwind_l[6]+0.4472135954999579*alphaSurf_l[2]*F_0_Upwind_l[5]+0.5*F_0_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[3] = 0.4*alphaSurf_l[3]*F_0_Upwind_l[8]+0.4*alphaSurf_l[6]*F_0_Upwind_l[7]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[7]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[6]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[5]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[3]*F_0_Upwind_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*F_0_Upwind_l[3]+0.5*F_0_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[4] = 0.31943828249997*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5000000000000001*F_0_Upwind_l[2]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[2]*F_0_Upwind_l[6]+0.31943828249997*F_0_Upwind_l[4]*alphaSurf_l[4]+0.5*F_0_Upwind_l[0]*alphaSurf_l[4]+0.5*alphaSurf_l[0]*F_0_Upwind_l[4]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[1]*alphaSurf_l[1]; 
  Ghat_F_0_l[5] = 0.5*alphaSurf_l[4]*F_0_Upwind_l[8]+0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[7]+0.4472135954999579*F_0_Upwind_l[6]*alphaSurf_l[6]+0.5*alphaSurf_l[0]*F_0_Upwind_l[5]+0.4472135954999579*F_0_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*F_0_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_F_0_l[6] = 0.2857142857142857*alphaSurf_l[6]*F_0_Upwind_l[8]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[8]+0.4*alphaSurf_l[3]*F_0_Upwind_l[7]+0.4472135954999579*F_0_Upwind_l[5]*alphaSurf_l[6]+0.31943828249997*F_0_Upwind_l[4]*alphaSurf_l[6]+0.5*F_0_Upwind_l[0]*alphaSurf_l[6]+0.31943828249997*alphaSurf_l[4]*F_0_Upwind_l[6]+0.5*alphaSurf_l[0]*F_0_Upwind_l[6]+0.5000000000000001*F_0_Upwind_l[2]*alphaSurf_l[4]+0.5000000000000001*alphaSurf_l[2]*F_0_Upwind_l[4]+0.447213595499958*F_0_Upwind_l[1]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[7] = 0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[8]+0.4472135954999579*alphaSurf_l[4]*F_0_Upwind_l[7]+0.5*alphaSurf_l[0]*F_0_Upwind_l[7]+0.4*F_0_Upwind_l[3]*alphaSurf_l[6]+0.4*alphaSurf_l[3]*F_0_Upwind_l[6]+0.5000000000000001*alphaSurf_l[1]*F_0_Upwind_l[5]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[8] = 0.31943828249997*alphaSurf_l[4]*F_0_Upwind_l[8]+0.5*alphaSurf_l[0]*F_0_Upwind_l[8]+0.447213595499958*alphaSurf_l[1]*F_0_Upwind_l[7]+0.2857142857142857*F_0_Upwind_l[6]*alphaSurf_l[6]+0.447213595499958*F_0_Upwind_l[2]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[2]*F_0_Upwind_l[6]+0.5*alphaSurf_l[4]*F_0_Upwind_l[5]+0.4*F_0_Upwind_l[3]*alphaSurf_l[3]; 
  Ghat_G_1_l[0] = 0.5*G_1_Upwind_l[6]*alphaSurf_l[6]+0.5*G_1_Upwind_l[4]*alphaSurf_l[4]+0.5*G_1_Upwind_l[3]*alphaSurf_l[3]+0.5*G_1_Upwind_l[2]*alphaSurf_l[2]+0.5*G_1_Upwind_l[1]*alphaSurf_l[1]+0.5*G_1_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_G_1_l[1] = 0.447213595499958*G_1_Upwind_l[3]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[3]*G_1_Upwind_l[6]+0.4472135954999579*G_1_Upwind_l[1]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[1]*G_1_Upwind_l[4]+0.5*G_1_Upwind_l[2]*alphaSurf_l[3]+0.5*alphaSurf_l[2]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_l[1]+0.5*alphaSurf_l[0]*G_1_Upwind_l[1]; 
  Ghat_G_1_l[2] = 0.447213595499958*alphaSurf_l[6]*G_1_Upwind_l[8]+0.447213595499958*alphaSurf_l[3]*G_1_Upwind_l[7]+0.5000000000000001*G_1_Upwind_l[4]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[4]*G_1_Upwind_l[6]+0.4472135954999579*alphaSurf_l[2]*G_1_Upwind_l[5]+0.5*G_1_Upwind_l[1]*alphaSurf_l[3]+0.5*alphaSurf_l[1]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[0]*alphaSurf_l[2]+0.5*alphaSurf_l[0]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[3] = 0.4*alphaSurf_l[3]*G_1_Upwind_l[8]+0.4*alphaSurf_l[6]*G_1_Upwind_l[7]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[7]+0.447213595499958*G_1_Upwind_l[1]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[6]+0.4472135954999579*alphaSurf_l[3]*G_1_Upwind_l[5]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[4]+0.4472135954999579*alphaSurf_l[3]*G_1_Upwind_l[4]+0.5*G_1_Upwind_l[0]*alphaSurf_l[3]+0.5*alphaSurf_l[0]*G_1_Upwind_l[3]+0.5*G_1_Upwind_l[1]*alphaSurf_l[2]+0.5*alphaSurf_l[1]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[4] = 0.31943828249997*G_1_Upwind_l[6]*alphaSurf_l[6]+0.5000000000000001*G_1_Upwind_l[2]*alphaSurf_l[6]+0.5000000000000001*alphaSurf_l[2]*G_1_Upwind_l[6]+0.31943828249997*G_1_Upwind_l[4]*alphaSurf_l[4]+0.5*G_1_Upwind_l[0]*alphaSurf_l[4]+0.5*alphaSurf_l[0]*G_1_Upwind_l[4]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*G_1_Upwind_l[1]*alphaSurf_l[1]; 
  Ghat_G_1_l[5] = 0.5*alphaSurf_l[4]*G_1_Upwind_l[8]+0.5000000000000001*alphaSurf_l[1]*G_1_Upwind_l[7]+0.4472135954999579*G_1_Upwind_l[6]*alphaSurf_l[6]+0.5*alphaSurf_l[0]*G_1_Upwind_l[5]+0.4472135954999579*G_1_Upwind_l[3]*alphaSurf_l[3]+0.4472135954999579*G_1_Upwind_l[2]*alphaSurf_l[2]; 
  Ghat_G_1_l[6] = 0.2857142857142857*alphaSurf_l[6]*G_1_Upwind_l[8]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[8]+0.4*alphaSurf_l[3]*G_1_Upwind_l[7]+0.4472135954999579*G_1_Upwind_l[5]*alphaSurf_l[6]+0.31943828249997*G_1_Upwind_l[4]*alphaSurf_l[6]+0.5*G_1_Upwind_l[0]*alphaSurf_l[6]+0.31943828249997*alphaSurf_l[4]*G_1_Upwind_l[6]+0.5*alphaSurf_l[0]*G_1_Upwind_l[6]+0.5000000000000001*G_1_Upwind_l[2]*alphaSurf_l[4]+0.5000000000000001*alphaSurf_l[2]*G_1_Upwind_l[4]+0.447213595499958*G_1_Upwind_l[1]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[7] = 0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[8]+0.4472135954999579*alphaSurf_l[4]*G_1_Upwind_l[7]+0.5*alphaSurf_l[0]*G_1_Upwind_l[7]+0.4*G_1_Upwind_l[3]*alphaSurf_l[6]+0.4*alphaSurf_l[3]*G_1_Upwind_l[6]+0.5000000000000001*alphaSurf_l[1]*G_1_Upwind_l[5]+0.447213595499958*G_1_Upwind_l[2]*alphaSurf_l[3]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[8] = 0.31943828249997*alphaSurf_l[4]*G_1_Upwind_l[8]+0.5*alphaSurf_l[0]*G_1_Upwind_l[8]+0.447213595499958*alphaSurf_l[1]*G_1_Upwind_l[7]+0.2857142857142857*G_1_Upwind_l[6]*alphaSurf_l[6]+0.447213595499958*G_1_Upwind_l[2]*alphaSurf_l[6]+0.447213595499958*alphaSurf_l[2]*G_1_Upwind_l[6]+0.5*alphaSurf_l[4]*G_1_Upwind_l[5]+0.4*G_1_Upwind_l[3]*alphaSurf_l[3]; 

  Ghat_F_0_r[0] = 0.5*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5*F_0_Upwind_r[4]*alphaSurf_r[4]+0.5*F_0_Upwind_r[3]*alphaSurf_r[3]+0.5*F_0_Upwind_r[2]*alphaSurf_r[2]+0.5*F_0_Upwind_r[1]*alphaSurf_r[1]+0.5*F_0_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_F_0_r[1] = 0.447213595499958*F_0_Upwind_r[3]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[6]+0.4472135954999579*F_0_Upwind_r[1]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[1]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.447213595499958*alphaSurf_r[6]*F_0_Upwind_r[8]+0.447213595499958*alphaSurf_r[3]*F_0_Upwind_r[7]+0.5000000000000001*F_0_Upwind_r[4]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[4]*F_0_Upwind_r[6]+0.4472135954999579*alphaSurf_r[2]*F_0_Upwind_r[5]+0.5*F_0_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[3] = 0.4*alphaSurf_r[3]*F_0_Upwind_r[8]+0.4*alphaSurf_r[6]*F_0_Upwind_r[7]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[7]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[6]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[5]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[3]*F_0_Upwind_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*F_0_Upwind_r[3]+0.5*F_0_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[4] = 0.31943828249997*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5000000000000001*F_0_Upwind_r[2]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[2]*F_0_Upwind_r[6]+0.31943828249997*F_0_Upwind_r[4]*alphaSurf_r[4]+0.5*F_0_Upwind_r[0]*alphaSurf_r[4]+0.5*alphaSurf_r[0]*F_0_Upwind_r[4]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[1]*alphaSurf_r[1]; 
  Ghat_F_0_r[5] = 0.5*alphaSurf_r[4]*F_0_Upwind_r[8]+0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[7]+0.4472135954999579*F_0_Upwind_r[6]*alphaSurf_r[6]+0.5*alphaSurf_r[0]*F_0_Upwind_r[5]+0.4472135954999579*F_0_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*F_0_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_F_0_r[6] = 0.2857142857142857*alphaSurf_r[6]*F_0_Upwind_r[8]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[8]+0.4*alphaSurf_r[3]*F_0_Upwind_r[7]+0.4472135954999579*F_0_Upwind_r[5]*alphaSurf_r[6]+0.31943828249997*F_0_Upwind_r[4]*alphaSurf_r[6]+0.5*F_0_Upwind_r[0]*alphaSurf_r[6]+0.31943828249997*alphaSurf_r[4]*F_0_Upwind_r[6]+0.5*alphaSurf_r[0]*F_0_Upwind_r[6]+0.5000000000000001*F_0_Upwind_r[2]*alphaSurf_r[4]+0.5000000000000001*alphaSurf_r[2]*F_0_Upwind_r[4]+0.447213595499958*F_0_Upwind_r[1]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[7] = 0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[8]+0.4472135954999579*alphaSurf_r[4]*F_0_Upwind_r[7]+0.5*alphaSurf_r[0]*F_0_Upwind_r[7]+0.4*F_0_Upwind_r[3]*alphaSurf_r[6]+0.4*alphaSurf_r[3]*F_0_Upwind_r[6]+0.5000000000000001*alphaSurf_r[1]*F_0_Upwind_r[5]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[8] = 0.31943828249997*alphaSurf_r[4]*F_0_Upwind_r[8]+0.5*alphaSurf_r[0]*F_0_Upwind_r[8]+0.447213595499958*alphaSurf_r[1]*F_0_Upwind_r[7]+0.2857142857142857*F_0_Upwind_r[6]*alphaSurf_r[6]+0.447213595499958*F_0_Upwind_r[2]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[2]*F_0_Upwind_r[6]+0.5*alphaSurf_r[4]*F_0_Upwind_r[5]+0.4*F_0_Upwind_r[3]*alphaSurf_r[3]; 
  Ghat_G_1_r[0] = 0.5*G_1_Upwind_r[6]*alphaSurf_r[6]+0.5*G_1_Upwind_r[4]*alphaSurf_r[4]+0.5*G_1_Upwind_r[3]*alphaSurf_r[3]+0.5*G_1_Upwind_r[2]*alphaSurf_r[2]+0.5*G_1_Upwind_r[1]*alphaSurf_r[1]+0.5*G_1_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_G_1_r[1] = 0.447213595499958*G_1_Upwind_r[3]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[3]*G_1_Upwind_r[6]+0.4472135954999579*G_1_Upwind_r[1]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[1]*G_1_Upwind_r[4]+0.5*G_1_Upwind_r[2]*alphaSurf_r[3]+0.5*alphaSurf_r[2]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_r[1]+0.5*alphaSurf_r[0]*G_1_Upwind_r[1]; 
  Ghat_G_1_r[2] = 0.447213595499958*alphaSurf_r[6]*G_1_Upwind_r[8]+0.447213595499958*alphaSurf_r[3]*G_1_Upwind_r[7]+0.5000000000000001*G_1_Upwind_r[4]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[4]*G_1_Upwind_r[6]+0.4472135954999579*alphaSurf_r[2]*G_1_Upwind_r[5]+0.5*G_1_Upwind_r[1]*alphaSurf_r[3]+0.5*alphaSurf_r[1]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[0]*alphaSurf_r[2]+0.5*alphaSurf_r[0]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[3] = 0.4*alphaSurf_r[3]*G_1_Upwind_r[8]+0.4*alphaSurf_r[6]*G_1_Upwind_r[7]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[7]+0.447213595499958*G_1_Upwind_r[1]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[6]+0.4472135954999579*alphaSurf_r[3]*G_1_Upwind_r[5]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[4]+0.4472135954999579*alphaSurf_r[3]*G_1_Upwind_r[4]+0.5*G_1_Upwind_r[0]*alphaSurf_r[3]+0.5*alphaSurf_r[0]*G_1_Upwind_r[3]+0.5*G_1_Upwind_r[1]*alphaSurf_r[2]+0.5*alphaSurf_r[1]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[4] = 0.31943828249997*G_1_Upwind_r[6]*alphaSurf_r[6]+0.5000000000000001*G_1_Upwind_r[2]*alphaSurf_r[6]+0.5000000000000001*alphaSurf_r[2]*G_1_Upwind_r[6]+0.31943828249997*G_1_Upwind_r[4]*alphaSurf_r[4]+0.5*G_1_Upwind_r[0]*alphaSurf_r[4]+0.5*alphaSurf_r[0]*G_1_Upwind_r[4]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*G_1_Upwind_r[1]*alphaSurf_r[1]; 
  Ghat_G_1_r[5] = 0.5*alphaSurf_r[4]*G_1_Upwind_r[8]+0.5000000000000001*alphaSurf_r[1]*G_1_Upwind_r[7]+0.4472135954999579*G_1_Upwind_r[6]*alphaSurf_r[6]+0.5*alphaSurf_r[0]*G_1_Upwind_r[5]+0.4472135954999579*G_1_Upwind_r[3]*alphaSurf_r[3]+0.4472135954999579*G_1_Upwind_r[2]*alphaSurf_r[2]; 
  Ghat_G_1_r[6] = 0.2857142857142857*alphaSurf_r[6]*G_1_Upwind_r[8]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[8]+0.4*alphaSurf_r[3]*G_1_Upwind_r[7]+0.4472135954999579*G_1_Upwind_r[5]*alphaSurf_r[6]+0.31943828249997*G_1_Upwind_r[4]*alphaSurf_r[6]+0.5*G_1_Upwind_r[0]*alphaSurf_r[6]+0.31943828249997*alphaSurf_r[4]*G_1_Upwind_r[6]+0.5*alphaSurf_r[0]*G_1_Upwind_r[6]+0.5000000000000001*G_1_Upwind_r[2]*alphaSurf_r[4]+0.5000000000000001*alphaSurf_r[2]*G_1_Upwind_r[4]+0.447213595499958*G_1_Upwind_r[1]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[7] = 0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[8]+0.4472135954999579*alphaSurf_r[4]*G_1_Upwind_r[7]+0.5*alphaSurf_r[0]*G_1_Upwind_r[7]+0.4*G_1_Upwind_r[3]*alphaSurf_r[6]+0.4*alphaSurf_r[3]*G_1_Upwind_r[6]+0.5000000000000001*alphaSurf_r[1]*G_1_Upwind_r[5]+0.447213595499958*G_1_Upwind_r[2]*alphaSurf_r[3]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[8] = 0.31943828249997*alphaSurf_r[4]*G_1_Upwind_r[8]+0.5*alphaSurf_r[0]*G_1_Upwind_r[8]+0.447213595499958*alphaSurf_r[1]*G_1_Upwind_r[7]+0.2857142857142857*G_1_Upwind_r[6]*alphaSurf_r[6]+0.447213595499958*G_1_Upwind_r[2]*alphaSurf_r[6]+0.447213595499958*alphaSurf_r[2]*G_1_Upwind_r[6]+0.5*alphaSurf_r[4]*G_1_Upwind_r[5]+0.4*G_1_Upwind_r[3]*alphaSurf_r[3]; 

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
  out_F_0[20] += (1.58113883008419*Ghat_F_0_l[4]-1.58113883008419*Ghat_F_0_r[4])*dx1; 
  out_F_0[21] += (0.7071067811865475*Ghat_F_0_l[8]-0.7071067811865475*Ghat_F_0_r[8])*dx1; 
  out_F_0[22] += (1.58113883008419*Ghat_F_0_l[5]-1.58113883008419*Ghat_F_0_r[5])*dx1; 
  out_F_0[23] += (1.58113883008419*Ghat_F_0_l[6]-1.58113883008419*Ghat_F_0_r[6])*dx1; 
  out_F_0[24] += -1.224744871391589*(Ghat_F_0_r[8]+Ghat_F_0_l[8])*dx1; 
  out_F_0[25] += (1.58113883008419*Ghat_F_0_l[7]-1.58113883008419*Ghat_F_0_r[7])*dx1; 
  out_F_0[26] += (1.58113883008419*Ghat_F_0_l[8]-1.58113883008419*Ghat_F_0_r[8])*dx1; 
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
  out_G_1[20] += (1.58113883008419*Ghat_G_1_l[4]-1.58113883008419*Ghat_G_1_r[4])*dx1; 
  out_G_1[21] += (0.7071067811865475*Ghat_G_1_l[8]-0.7071067811865475*Ghat_G_1_r[8])*dx1; 
  out_G_1[22] += (1.58113883008419*Ghat_G_1_l[5]-1.58113883008419*Ghat_G_1_r[5])*dx1; 
  out_G_1[23] += (1.58113883008419*Ghat_G_1_l[6]-1.58113883008419*Ghat_G_1_r[6])*dx1; 
  out_G_1[24] += -1.224744871391589*(Ghat_G_1_r[8]+Ghat_G_1_l[8])*dx1; 
  out_G_1[25] += (1.58113883008419*Ghat_G_1_l[7]-1.58113883008419*Ghat_G_1_r[7])*dx1; 
  out_G_1[26] += (1.58113883008419*Ghat_G_1_l[8]-1.58113883008419*Ghat_G_1_r[8])*dx1; 

  const double *u_surf_lr = &pkpm_prim_surf_l[33]; 
  const double *u_surf_cl = &pkpm_prim_surf_c[30]; 
  const double *u_surf_cr = &pkpm_prim_surf_c[33]; 
  const double *u_surf_rl = &pkpm_prim_surf_r[30]; 

  const double *pkpm_lax_l = &pkpm_lax[6]; 
  const double *pkpm_lax_r = &pkpm_lax[9]; 

  double F_0_lr[9] = {0.0}; 
  double F_0_cl[9] = {0.0}; 
  double F_0_cr[9] = {0.0}; 
  double F_0_rl[9] = {0.0}; 
  double Ghat_F_0_u_l[9] = {0.0}; 
  double Ghat_F_0_u_r[9] = {0.0}; 

  double G_1_lr[9] = {0.0}; 
  double G_1_cl[9] = {0.0}; 
  double G_1_cr[9] = {0.0}; 
  double G_1_rl[9] = {0.0}; 
  double Ghat_G_1_u_l[9] = {0.0}; 
  double Ghat_G_1_u_r[9] = {0.0}; 

  F_0_lr[0] = (-4.415885040048121e-16*F_0l[20])+1.58113883008419*F_0l[8]+1.224744871391589*F_0l[2]+0.7071067811865475*F_0l[0]; 
  F_0_lr[1] = 1.58113883008419*F_0l[12]+1.224744871391589*F_0l[4]+0.7071067811865475*F_0l[1]; 
  F_0_lr[2] = (-4.415885040048115e-16*F_0l[23])+1.58113883008419*F_0l[14]+1.224744871391589*F_0l[6]+0.7071067811865475*F_0l[3]; 
  F_0_lr[3] = 1.58113883008419*F_0l[18]+1.224744871391589*F_0l[10]+0.7071067811865475*F_0l[5]; 
  F_0_lr[4] = 1.58113883008419*F_0l[20]+1.224744871391589*F_0l[11]+0.7071067811865475*F_0l[7]; 
  F_0_lr[5] = 1.58113883008419*F_0l[22]+1.224744871391589*F_0l[16]+0.7071067811865475*F_0l[9]; 
  F_0_lr[6] = 1.58113883008419*F_0l[23]+1.224744871391589*F_0l[17]+0.7071067811865475*F_0l[13]; 
  F_0_lr[7] = 1.58113883008419*F_0l[25]+1.224744871391589*F_0l[19]+0.7071067811865475*F_0l[15]; 
  F_0_lr[8] = 1.58113883008419*F_0l[26]+1.224744871391589*F_0l[24]+0.7071067811865475*F_0l[21]; 

  F_0_cl[0] = (-4.415885040048121e-16*F_0c[20])+1.58113883008419*F_0c[8]-1.224744871391589*F_0c[2]+0.7071067811865475*F_0c[0]; 
  F_0_cl[1] = 1.58113883008419*F_0c[12]-1.224744871391589*F_0c[4]+0.7071067811865475*F_0c[1]; 
  F_0_cl[2] = (-4.415885040048115e-16*F_0c[23])+1.58113883008419*F_0c[14]-1.224744871391589*F_0c[6]+0.7071067811865475*F_0c[3]; 
  F_0_cl[3] = 1.58113883008419*F_0c[18]-1.224744871391589*F_0c[10]+0.7071067811865475*F_0c[5]; 
  F_0_cl[4] = 1.58113883008419*F_0c[20]-1.224744871391589*F_0c[11]+0.7071067811865475*F_0c[7]; 
  F_0_cl[5] = 1.58113883008419*F_0c[22]-1.224744871391589*F_0c[16]+0.7071067811865475*F_0c[9]; 
  F_0_cl[6] = 1.58113883008419*F_0c[23]-1.224744871391589*F_0c[17]+0.7071067811865475*F_0c[13]; 
  F_0_cl[7] = 1.58113883008419*F_0c[25]-1.224744871391589*F_0c[19]+0.7071067811865475*F_0c[15]; 
  F_0_cl[8] = 1.58113883008419*F_0c[26]-1.224744871391589*F_0c[24]+0.7071067811865475*F_0c[21]; 

  F_0_cr[0] = (-4.415885040048121e-16*F_0c[20])+1.58113883008419*F_0c[8]+1.224744871391589*F_0c[2]+0.7071067811865475*F_0c[0]; 
  F_0_cr[1] = 1.58113883008419*F_0c[12]+1.224744871391589*F_0c[4]+0.7071067811865475*F_0c[1]; 
  F_0_cr[2] = (-4.415885040048115e-16*F_0c[23])+1.58113883008419*F_0c[14]+1.224744871391589*F_0c[6]+0.7071067811865475*F_0c[3]; 
  F_0_cr[3] = 1.58113883008419*F_0c[18]+1.224744871391589*F_0c[10]+0.7071067811865475*F_0c[5]; 
  F_0_cr[4] = 1.58113883008419*F_0c[20]+1.224744871391589*F_0c[11]+0.7071067811865475*F_0c[7]; 
  F_0_cr[5] = 1.58113883008419*F_0c[22]+1.224744871391589*F_0c[16]+0.7071067811865475*F_0c[9]; 
  F_0_cr[6] = 1.58113883008419*F_0c[23]+1.224744871391589*F_0c[17]+0.7071067811865475*F_0c[13]; 
  F_0_cr[7] = 1.58113883008419*F_0c[25]+1.224744871391589*F_0c[19]+0.7071067811865475*F_0c[15]; 
  F_0_cr[8] = 1.58113883008419*F_0c[26]+1.224744871391589*F_0c[24]+0.7071067811865475*F_0c[21]; 

  F_0_rl[0] = (-4.415885040048121e-16*F_0r[20])+1.58113883008419*F_0r[8]-1.224744871391589*F_0r[2]+0.7071067811865475*F_0r[0]; 
  F_0_rl[1] = 1.58113883008419*F_0r[12]-1.224744871391589*F_0r[4]+0.7071067811865475*F_0r[1]; 
  F_0_rl[2] = (-4.415885040048115e-16*F_0r[23])+1.58113883008419*F_0r[14]-1.224744871391589*F_0r[6]+0.7071067811865475*F_0r[3]; 
  F_0_rl[3] = 1.58113883008419*F_0r[18]-1.224744871391589*F_0r[10]+0.7071067811865475*F_0r[5]; 
  F_0_rl[4] = 1.58113883008419*F_0r[20]-1.224744871391589*F_0r[11]+0.7071067811865475*F_0r[7]; 
  F_0_rl[5] = 1.58113883008419*F_0r[22]-1.224744871391589*F_0r[16]+0.7071067811865475*F_0r[9]; 
  F_0_rl[6] = 1.58113883008419*F_0r[23]-1.224744871391589*F_0r[17]+0.7071067811865475*F_0r[13]; 
  F_0_rl[7] = 1.58113883008419*F_0r[25]-1.224744871391589*F_0r[19]+0.7071067811865475*F_0r[15]; 
  F_0_rl[8] = 1.58113883008419*F_0r[26]-1.224744871391589*F_0r[24]+0.7071067811865475*F_0r[21]; 

  G_1_lr[0] = (-4.415885040048121e-16*G_1l[20])+1.58113883008419*G_1l[8]+1.224744871391589*G_1l[2]+0.7071067811865475*G_1l[0]; 
  G_1_lr[1] = 1.58113883008419*G_1l[12]+1.224744871391589*G_1l[4]+0.7071067811865475*G_1l[1]; 
  G_1_lr[2] = (-4.415885040048115e-16*G_1l[23])+1.58113883008419*G_1l[14]+1.224744871391589*G_1l[6]+0.7071067811865475*G_1l[3]; 
  G_1_lr[3] = 1.58113883008419*G_1l[18]+1.224744871391589*G_1l[10]+0.7071067811865475*G_1l[5]; 
  G_1_lr[4] = 1.58113883008419*G_1l[20]+1.224744871391589*G_1l[11]+0.7071067811865475*G_1l[7]; 
  G_1_lr[5] = 1.58113883008419*G_1l[22]+1.224744871391589*G_1l[16]+0.7071067811865475*G_1l[9]; 
  G_1_lr[6] = 1.58113883008419*G_1l[23]+1.224744871391589*G_1l[17]+0.7071067811865475*G_1l[13]; 
  G_1_lr[7] = 1.58113883008419*G_1l[25]+1.224744871391589*G_1l[19]+0.7071067811865475*G_1l[15]; 
  G_1_lr[8] = 1.58113883008419*G_1l[26]+1.224744871391589*G_1l[24]+0.7071067811865475*G_1l[21]; 

  G_1_cl[0] = (-4.415885040048121e-16*G_1c[20])+1.58113883008419*G_1c[8]-1.224744871391589*G_1c[2]+0.7071067811865475*G_1c[0]; 
  G_1_cl[1] = 1.58113883008419*G_1c[12]-1.224744871391589*G_1c[4]+0.7071067811865475*G_1c[1]; 
  G_1_cl[2] = (-4.415885040048115e-16*G_1c[23])+1.58113883008419*G_1c[14]-1.224744871391589*G_1c[6]+0.7071067811865475*G_1c[3]; 
  G_1_cl[3] = 1.58113883008419*G_1c[18]-1.224744871391589*G_1c[10]+0.7071067811865475*G_1c[5]; 
  G_1_cl[4] = 1.58113883008419*G_1c[20]-1.224744871391589*G_1c[11]+0.7071067811865475*G_1c[7]; 
  G_1_cl[5] = 1.58113883008419*G_1c[22]-1.224744871391589*G_1c[16]+0.7071067811865475*G_1c[9]; 
  G_1_cl[6] = 1.58113883008419*G_1c[23]-1.224744871391589*G_1c[17]+0.7071067811865475*G_1c[13]; 
  G_1_cl[7] = 1.58113883008419*G_1c[25]-1.224744871391589*G_1c[19]+0.7071067811865475*G_1c[15]; 
  G_1_cl[8] = 1.58113883008419*G_1c[26]-1.224744871391589*G_1c[24]+0.7071067811865475*G_1c[21]; 

  G_1_cr[0] = (-4.415885040048121e-16*G_1c[20])+1.58113883008419*G_1c[8]+1.224744871391589*G_1c[2]+0.7071067811865475*G_1c[0]; 
  G_1_cr[1] = 1.58113883008419*G_1c[12]+1.224744871391589*G_1c[4]+0.7071067811865475*G_1c[1]; 
  G_1_cr[2] = (-4.415885040048115e-16*G_1c[23])+1.58113883008419*G_1c[14]+1.224744871391589*G_1c[6]+0.7071067811865475*G_1c[3]; 
  G_1_cr[3] = 1.58113883008419*G_1c[18]+1.224744871391589*G_1c[10]+0.7071067811865475*G_1c[5]; 
  G_1_cr[4] = 1.58113883008419*G_1c[20]+1.224744871391589*G_1c[11]+0.7071067811865475*G_1c[7]; 
  G_1_cr[5] = 1.58113883008419*G_1c[22]+1.224744871391589*G_1c[16]+0.7071067811865475*G_1c[9]; 
  G_1_cr[6] = 1.58113883008419*G_1c[23]+1.224744871391589*G_1c[17]+0.7071067811865475*G_1c[13]; 
  G_1_cr[7] = 1.58113883008419*G_1c[25]+1.224744871391589*G_1c[19]+0.7071067811865475*G_1c[15]; 
  G_1_cr[8] = 1.58113883008419*G_1c[26]+1.224744871391589*G_1c[24]+0.7071067811865475*G_1c[21]; 

  G_1_rl[0] = (-4.415885040048121e-16*G_1r[20])+1.58113883008419*G_1r[8]-1.224744871391589*G_1r[2]+0.7071067811865475*G_1r[0]; 
  G_1_rl[1] = 1.58113883008419*G_1r[12]-1.224744871391589*G_1r[4]+0.7071067811865475*G_1r[1]; 
  G_1_rl[2] = (-4.415885040048115e-16*G_1r[23])+1.58113883008419*G_1r[14]-1.224744871391589*G_1r[6]+0.7071067811865475*G_1r[3]; 
  G_1_rl[3] = 1.58113883008419*G_1r[18]-1.224744871391589*G_1r[10]+0.7071067811865475*G_1r[5]; 
  G_1_rl[4] = 1.58113883008419*G_1r[20]-1.224744871391589*G_1r[11]+0.7071067811865475*G_1r[7]; 
  G_1_rl[5] = 1.58113883008419*G_1r[22]-1.224744871391589*G_1r[16]+0.7071067811865475*G_1r[9]; 
  G_1_rl[6] = 1.58113883008419*G_1r[23]-1.224744871391589*G_1r[17]+0.7071067811865475*G_1r[13]; 
  G_1_rl[7] = 1.58113883008419*G_1r[25]-1.224744871391589*G_1r[19]+0.7071067811865475*G_1r[15]; 
  G_1_rl[8] = 1.58113883008419*G_1r[26]-1.224744871391589*G_1r[24]+0.7071067811865475*G_1r[21]; 

  Ghat_F_0_u_l[0] = 0.1767766952966368*u_surf_lr[2]*F_0_lr[4]+0.1767766952966368*u_surf_cl[2]*F_0_lr[4]+0.3535533905932737*pkpm_lax_l[2]*F_0_lr[4]+0.1767766952966368*u_surf_lr[2]*F_0_cl[4]+0.1767766952966368*u_surf_cl[2]*F_0_cl[4]-0.3535533905932737*pkpm_lax_l[2]*F_0_cl[4]+0.1767766952966368*F_0_lr[1]*u_surf_lr[1]+0.1767766952966368*F_0_cl[1]*u_surf_lr[1]+0.1767766952966368*F_0_lr[1]*u_surf_cl[1]+0.1767766952966368*F_0_cl[1]*u_surf_cl[1]+0.3535533905932737*F_0_lr[1]*pkpm_lax_l[1]-0.3535533905932737*F_0_cl[1]*pkpm_lax_l[1]+0.1767766952966368*F_0_lr[0]*u_surf_lr[0]+0.1767766952966368*F_0_cl[0]*u_surf_lr[0]+0.1767766952966368*F_0_lr[0]*u_surf_cl[0]+0.1767766952966368*F_0_cl[0]*u_surf_cl[0]+0.3535533905932737*F_0_lr[0]*pkpm_lax_l[0]-0.3535533905932737*F_0_cl[0]*pkpm_lax_l[0]; 
  Ghat_F_0_u_l[1] = 0.1581138830084189*u_surf_lr[1]*F_0_lr[4]+0.1581138830084189*u_surf_cl[1]*F_0_lr[4]+0.3162277660168379*pkpm_lax_l[1]*F_0_lr[4]+0.1581138830084189*u_surf_lr[1]*F_0_cl[4]+0.1581138830084189*u_surf_cl[1]*F_0_cl[4]-0.3162277660168379*pkpm_lax_l[1]*F_0_cl[4]+0.1581138830084189*F_0_lr[1]*u_surf_lr[2]+0.1581138830084189*F_0_cl[1]*u_surf_lr[2]+0.1581138830084189*F_0_lr[1]*u_surf_cl[2]+0.1581138830084189*F_0_cl[1]*u_surf_cl[2]+0.3162277660168379*F_0_lr[1]*pkpm_lax_l[2]-0.3162277660168379*F_0_cl[1]*pkpm_lax_l[2]+0.1767766952966368*F_0_lr[0]*u_surf_lr[1]+0.1767766952966368*F_0_cl[0]*u_surf_lr[1]+0.1767766952966368*F_0_lr[0]*u_surf_cl[1]+0.1767766952966368*F_0_cl[0]*u_surf_cl[1]+0.3535533905932737*F_0_lr[0]*pkpm_lax_l[1]-0.3535533905932737*F_0_cl[0]*pkpm_lax_l[1]+0.1767766952966368*u_surf_lr[0]*F_0_lr[1]+0.1767766952966368*u_surf_cl[0]*F_0_lr[1]+0.3535533905932737*pkpm_lax_l[0]*F_0_lr[1]+0.1767766952966368*u_surf_lr[0]*F_0_cl[1]+0.1767766952966368*u_surf_cl[0]*F_0_cl[1]-0.3535533905932737*pkpm_lax_l[0]*F_0_cl[1]; 
  Ghat_F_0_u_l[2] = 0.1767766952966368*u_surf_lr[2]*F_0_lr[6]+0.1767766952966368*u_surf_cl[2]*F_0_lr[6]+0.3535533905932737*pkpm_lax_l[2]*F_0_lr[6]+0.1767766952966368*u_surf_lr[2]*F_0_cl[6]+0.1767766952966368*u_surf_cl[2]*F_0_cl[6]-0.3535533905932737*pkpm_lax_l[2]*F_0_cl[6]+0.1767766952966368*u_surf_lr[1]*F_0_lr[3]+0.1767766952966368*u_surf_cl[1]*F_0_lr[3]+0.3535533905932737*pkpm_lax_l[1]*F_0_lr[3]+0.1767766952966368*u_surf_lr[1]*F_0_cl[3]+0.1767766952966368*u_surf_cl[1]*F_0_cl[3]-0.3535533905932737*pkpm_lax_l[1]*F_0_cl[3]+0.1767766952966368*u_surf_lr[0]*F_0_lr[2]+0.1767766952966368*u_surf_cl[0]*F_0_lr[2]+0.3535533905932737*pkpm_lax_l[0]*F_0_lr[2]+0.1767766952966368*u_surf_lr[0]*F_0_cl[2]+0.1767766952966368*u_surf_cl[0]*F_0_cl[2]-0.3535533905932737*pkpm_lax_l[0]*F_0_cl[2]; 
  Ghat_F_0_u_l[3] = 0.1581138830084189*u_surf_lr[1]*F_0_lr[6]+0.1581138830084189*u_surf_cl[1]*F_0_lr[6]+0.3162277660168379*pkpm_lax_l[1]*F_0_lr[6]+0.1581138830084189*u_surf_lr[1]*F_0_cl[6]+0.1581138830084189*u_surf_cl[1]*F_0_cl[6]-0.3162277660168379*pkpm_lax_l[1]*F_0_cl[6]+0.1581138830084189*u_surf_lr[2]*F_0_lr[3]+0.1581138830084189*u_surf_cl[2]*F_0_lr[3]+0.3162277660168379*pkpm_lax_l[2]*F_0_lr[3]+0.1767766952966368*u_surf_lr[0]*F_0_lr[3]+0.1767766952966368*u_surf_cl[0]*F_0_lr[3]+0.3535533905932737*pkpm_lax_l[0]*F_0_lr[3]+0.1581138830084189*u_surf_lr[2]*F_0_cl[3]+0.1581138830084189*u_surf_cl[2]*F_0_cl[3]-0.3162277660168379*pkpm_lax_l[2]*F_0_cl[3]+0.1767766952966368*u_surf_lr[0]*F_0_cl[3]+0.1767766952966368*u_surf_cl[0]*F_0_cl[3]-0.3535533905932737*pkpm_lax_l[0]*F_0_cl[3]+0.1767766952966368*u_surf_lr[1]*F_0_lr[2]+0.1767766952966368*u_surf_cl[1]*F_0_lr[2]+0.3535533905932737*pkpm_lax_l[1]*F_0_lr[2]+0.1767766952966368*u_surf_lr[1]*F_0_cl[2]+0.1767766952966368*u_surf_cl[1]*F_0_cl[2]-0.3535533905932737*pkpm_lax_l[1]*F_0_cl[2]; 
  Ghat_F_0_u_l[4] = 0.1129384878631564*u_surf_lr[2]*F_0_lr[4]+0.1129384878631564*u_surf_cl[2]*F_0_lr[4]+0.2258769757263128*pkpm_lax_l[2]*F_0_lr[4]+0.1767766952966368*u_surf_lr[0]*F_0_lr[4]+0.1767766952966368*u_surf_cl[0]*F_0_lr[4]+0.3535533905932737*pkpm_lax_l[0]*F_0_lr[4]+0.1129384878631564*u_surf_lr[2]*F_0_cl[4]+0.1129384878631564*u_surf_cl[2]*F_0_cl[4]-0.2258769757263128*pkpm_lax_l[2]*F_0_cl[4]+0.1767766952966368*u_surf_lr[0]*F_0_cl[4]+0.1767766952966368*u_surf_cl[0]*F_0_cl[4]-0.3535533905932737*pkpm_lax_l[0]*F_0_cl[4]+0.1767766952966368*F_0_lr[0]*u_surf_lr[2]+0.1767766952966368*F_0_cl[0]*u_surf_lr[2]+0.1767766952966368*F_0_lr[0]*u_surf_cl[2]+0.1767766952966368*F_0_cl[0]*u_surf_cl[2]+0.3535533905932737*F_0_lr[0]*pkpm_lax_l[2]-0.3535533905932737*F_0_cl[0]*pkpm_lax_l[2]+0.1581138830084189*F_0_lr[1]*u_surf_lr[1]+0.1581138830084189*F_0_cl[1]*u_surf_lr[1]+0.1581138830084189*F_0_lr[1]*u_surf_cl[1]+0.1581138830084189*F_0_cl[1]*u_surf_cl[1]+0.3162277660168379*F_0_lr[1]*pkpm_lax_l[1]-0.3162277660168379*F_0_cl[1]*pkpm_lax_l[1]; 
  Ghat_F_0_u_l[5] = 0.1767766952966368*u_surf_lr[2]*F_0_lr[8]+0.1767766952966368*u_surf_cl[2]*F_0_lr[8]+0.3535533905932737*pkpm_lax_l[2]*F_0_lr[8]+0.1767766952966368*u_surf_lr[2]*F_0_cl[8]+0.1767766952966368*u_surf_cl[2]*F_0_cl[8]-0.3535533905932737*pkpm_lax_l[2]*F_0_cl[8]+0.1767766952966368*u_surf_lr[1]*F_0_lr[7]+0.1767766952966368*u_surf_cl[1]*F_0_lr[7]+0.3535533905932737*pkpm_lax_l[1]*F_0_lr[7]+0.1767766952966368*u_surf_lr[1]*F_0_cl[7]+0.1767766952966368*u_surf_cl[1]*F_0_cl[7]-0.3535533905932737*pkpm_lax_l[1]*F_0_cl[7]+0.1767766952966368*u_surf_lr[0]*F_0_lr[5]+0.1767766952966368*u_surf_cl[0]*F_0_lr[5]+0.3535533905932737*pkpm_lax_l[0]*F_0_lr[5]+0.1767766952966368*u_surf_lr[0]*F_0_cl[5]+0.1767766952966368*u_surf_cl[0]*F_0_cl[5]-0.3535533905932737*pkpm_lax_l[0]*F_0_cl[5]; 
  Ghat_F_0_u_l[6] = 0.1129384878631564*u_surf_lr[2]*F_0_lr[6]+0.1129384878631564*u_surf_cl[2]*F_0_lr[6]+0.2258769757263128*pkpm_lax_l[2]*F_0_lr[6]+0.1767766952966368*u_surf_lr[0]*F_0_lr[6]+0.1767766952966368*u_surf_cl[0]*F_0_lr[6]+0.3535533905932737*pkpm_lax_l[0]*F_0_lr[6]+0.1129384878631564*u_surf_lr[2]*F_0_cl[6]+0.1129384878631564*u_surf_cl[2]*F_0_cl[6]-0.2258769757263128*pkpm_lax_l[2]*F_0_cl[6]+0.1767766952966368*u_surf_lr[0]*F_0_cl[6]+0.1767766952966368*u_surf_cl[0]*F_0_cl[6]-0.3535533905932737*pkpm_lax_l[0]*F_0_cl[6]+0.1581138830084189*u_surf_lr[1]*F_0_lr[3]+0.1581138830084189*u_surf_cl[1]*F_0_lr[3]+0.3162277660168379*pkpm_lax_l[1]*F_0_lr[3]+0.1581138830084189*u_surf_lr[1]*F_0_cl[3]+0.1581138830084189*u_surf_cl[1]*F_0_cl[3]-0.3162277660168379*pkpm_lax_l[1]*F_0_cl[3]+0.1767766952966368*F_0_lr[2]*u_surf_lr[2]+0.1767766952966368*F_0_cl[2]*u_surf_lr[2]+0.1767766952966368*F_0_lr[2]*u_surf_cl[2]+0.1767766952966368*F_0_cl[2]*u_surf_cl[2]+0.3535533905932737*F_0_lr[2]*pkpm_lax_l[2]-0.3535533905932737*F_0_cl[2]*pkpm_lax_l[2]; 
  Ghat_F_0_u_l[7] = 0.1581138830084189*u_surf_lr[1]*F_0_lr[8]+0.1581138830084189*u_surf_cl[1]*F_0_lr[8]+0.3162277660168379*pkpm_lax_l[1]*F_0_lr[8]+0.1581138830084189*u_surf_lr[1]*F_0_cl[8]+0.1581138830084189*u_surf_cl[1]*F_0_cl[8]-0.3162277660168379*pkpm_lax_l[1]*F_0_cl[8]+0.1581138830084189*u_surf_lr[2]*F_0_lr[7]+0.1581138830084189*u_surf_cl[2]*F_0_lr[7]+0.3162277660168379*pkpm_lax_l[2]*F_0_lr[7]+0.1767766952966368*u_surf_lr[0]*F_0_lr[7]+0.1767766952966368*u_surf_cl[0]*F_0_lr[7]+0.3535533905932737*pkpm_lax_l[0]*F_0_lr[7]+0.1581138830084189*u_surf_lr[2]*F_0_cl[7]+0.1581138830084189*u_surf_cl[2]*F_0_cl[7]-0.3162277660168379*pkpm_lax_l[2]*F_0_cl[7]+0.1767766952966368*u_surf_lr[0]*F_0_cl[7]+0.1767766952966368*u_surf_cl[0]*F_0_cl[7]-0.3535533905932737*pkpm_lax_l[0]*F_0_cl[7]+0.1767766952966368*u_surf_lr[1]*F_0_lr[5]+0.1767766952966368*u_surf_cl[1]*F_0_lr[5]+0.3535533905932737*pkpm_lax_l[1]*F_0_lr[5]+0.1767766952966368*u_surf_lr[1]*F_0_cl[5]+0.1767766952966368*u_surf_cl[1]*F_0_cl[5]-0.3535533905932737*pkpm_lax_l[1]*F_0_cl[5]; 
  Ghat_F_0_u_l[8] = 0.1129384878631564*u_surf_lr[2]*F_0_lr[8]+0.1129384878631564*u_surf_cl[2]*F_0_lr[8]+0.2258769757263128*pkpm_lax_l[2]*F_0_lr[8]+0.1767766952966368*u_surf_lr[0]*F_0_lr[8]+0.1767766952966368*u_surf_cl[0]*F_0_lr[8]+0.3535533905932737*pkpm_lax_l[0]*F_0_lr[8]+0.1129384878631564*u_surf_lr[2]*F_0_cl[8]+0.1129384878631564*u_surf_cl[2]*F_0_cl[8]-0.2258769757263128*pkpm_lax_l[2]*F_0_cl[8]+0.1767766952966368*u_surf_lr[0]*F_0_cl[8]+0.1767766952966368*u_surf_cl[0]*F_0_cl[8]-0.3535533905932737*pkpm_lax_l[0]*F_0_cl[8]+0.1581138830084189*u_surf_lr[1]*F_0_lr[7]+0.1581138830084189*u_surf_cl[1]*F_0_lr[7]+0.3162277660168379*pkpm_lax_l[1]*F_0_lr[7]+0.1581138830084189*u_surf_lr[1]*F_0_cl[7]+0.1581138830084189*u_surf_cl[1]*F_0_cl[7]-0.3162277660168379*pkpm_lax_l[1]*F_0_cl[7]+0.1767766952966368*u_surf_lr[2]*F_0_lr[5]+0.1767766952966368*u_surf_cl[2]*F_0_lr[5]+0.3535533905932737*pkpm_lax_l[2]*F_0_lr[5]+0.1767766952966368*u_surf_lr[2]*F_0_cl[5]+0.1767766952966368*u_surf_cl[2]*F_0_cl[5]-0.3535533905932737*pkpm_lax_l[2]*F_0_cl[5]; 
  Ghat_G_1_u_l[0] = 0.1767766952966368*u_surf_lr[2]*G_1_lr[4]+0.1767766952966368*u_surf_cl[2]*G_1_lr[4]+0.3535533905932737*pkpm_lax_l[2]*G_1_lr[4]+0.1767766952966368*u_surf_lr[2]*G_1_cl[4]+0.1767766952966368*u_surf_cl[2]*G_1_cl[4]-0.3535533905932737*pkpm_lax_l[2]*G_1_cl[4]+0.1767766952966368*G_1_lr[1]*u_surf_lr[1]+0.1767766952966368*G_1_cl[1]*u_surf_lr[1]+0.1767766952966368*G_1_lr[1]*u_surf_cl[1]+0.1767766952966368*G_1_cl[1]*u_surf_cl[1]+0.3535533905932737*G_1_lr[1]*pkpm_lax_l[1]-0.3535533905932737*G_1_cl[1]*pkpm_lax_l[1]+0.1767766952966368*G_1_lr[0]*u_surf_lr[0]+0.1767766952966368*G_1_cl[0]*u_surf_lr[0]+0.1767766952966368*G_1_lr[0]*u_surf_cl[0]+0.1767766952966368*G_1_cl[0]*u_surf_cl[0]+0.3535533905932737*G_1_lr[0]*pkpm_lax_l[0]-0.3535533905932737*G_1_cl[0]*pkpm_lax_l[0]; 
  Ghat_G_1_u_l[1] = 0.1581138830084189*u_surf_lr[1]*G_1_lr[4]+0.1581138830084189*u_surf_cl[1]*G_1_lr[4]+0.3162277660168379*pkpm_lax_l[1]*G_1_lr[4]+0.1581138830084189*u_surf_lr[1]*G_1_cl[4]+0.1581138830084189*u_surf_cl[1]*G_1_cl[4]-0.3162277660168379*pkpm_lax_l[1]*G_1_cl[4]+0.1581138830084189*G_1_lr[1]*u_surf_lr[2]+0.1581138830084189*G_1_cl[1]*u_surf_lr[2]+0.1581138830084189*G_1_lr[1]*u_surf_cl[2]+0.1581138830084189*G_1_cl[1]*u_surf_cl[2]+0.3162277660168379*G_1_lr[1]*pkpm_lax_l[2]-0.3162277660168379*G_1_cl[1]*pkpm_lax_l[2]+0.1767766952966368*G_1_lr[0]*u_surf_lr[1]+0.1767766952966368*G_1_cl[0]*u_surf_lr[1]+0.1767766952966368*G_1_lr[0]*u_surf_cl[1]+0.1767766952966368*G_1_cl[0]*u_surf_cl[1]+0.3535533905932737*G_1_lr[0]*pkpm_lax_l[1]-0.3535533905932737*G_1_cl[0]*pkpm_lax_l[1]+0.1767766952966368*u_surf_lr[0]*G_1_lr[1]+0.1767766952966368*u_surf_cl[0]*G_1_lr[1]+0.3535533905932737*pkpm_lax_l[0]*G_1_lr[1]+0.1767766952966368*u_surf_lr[0]*G_1_cl[1]+0.1767766952966368*u_surf_cl[0]*G_1_cl[1]-0.3535533905932737*pkpm_lax_l[0]*G_1_cl[1]; 
  Ghat_G_1_u_l[2] = 0.1767766952966368*u_surf_lr[2]*G_1_lr[6]+0.1767766952966368*u_surf_cl[2]*G_1_lr[6]+0.3535533905932737*pkpm_lax_l[2]*G_1_lr[6]+0.1767766952966368*u_surf_lr[2]*G_1_cl[6]+0.1767766952966368*u_surf_cl[2]*G_1_cl[6]-0.3535533905932737*pkpm_lax_l[2]*G_1_cl[6]+0.1767766952966368*u_surf_lr[1]*G_1_lr[3]+0.1767766952966368*u_surf_cl[1]*G_1_lr[3]+0.3535533905932737*pkpm_lax_l[1]*G_1_lr[3]+0.1767766952966368*u_surf_lr[1]*G_1_cl[3]+0.1767766952966368*u_surf_cl[1]*G_1_cl[3]-0.3535533905932737*pkpm_lax_l[1]*G_1_cl[3]+0.1767766952966368*u_surf_lr[0]*G_1_lr[2]+0.1767766952966368*u_surf_cl[0]*G_1_lr[2]+0.3535533905932737*pkpm_lax_l[0]*G_1_lr[2]+0.1767766952966368*u_surf_lr[0]*G_1_cl[2]+0.1767766952966368*u_surf_cl[0]*G_1_cl[2]-0.3535533905932737*pkpm_lax_l[0]*G_1_cl[2]; 
  Ghat_G_1_u_l[3] = 0.1581138830084189*u_surf_lr[1]*G_1_lr[6]+0.1581138830084189*u_surf_cl[1]*G_1_lr[6]+0.3162277660168379*pkpm_lax_l[1]*G_1_lr[6]+0.1581138830084189*u_surf_lr[1]*G_1_cl[6]+0.1581138830084189*u_surf_cl[1]*G_1_cl[6]-0.3162277660168379*pkpm_lax_l[1]*G_1_cl[6]+0.1581138830084189*u_surf_lr[2]*G_1_lr[3]+0.1581138830084189*u_surf_cl[2]*G_1_lr[3]+0.3162277660168379*pkpm_lax_l[2]*G_1_lr[3]+0.1767766952966368*u_surf_lr[0]*G_1_lr[3]+0.1767766952966368*u_surf_cl[0]*G_1_lr[3]+0.3535533905932737*pkpm_lax_l[0]*G_1_lr[3]+0.1581138830084189*u_surf_lr[2]*G_1_cl[3]+0.1581138830084189*u_surf_cl[2]*G_1_cl[3]-0.3162277660168379*pkpm_lax_l[2]*G_1_cl[3]+0.1767766952966368*u_surf_lr[0]*G_1_cl[3]+0.1767766952966368*u_surf_cl[0]*G_1_cl[3]-0.3535533905932737*pkpm_lax_l[0]*G_1_cl[3]+0.1767766952966368*u_surf_lr[1]*G_1_lr[2]+0.1767766952966368*u_surf_cl[1]*G_1_lr[2]+0.3535533905932737*pkpm_lax_l[1]*G_1_lr[2]+0.1767766952966368*u_surf_lr[1]*G_1_cl[2]+0.1767766952966368*u_surf_cl[1]*G_1_cl[2]-0.3535533905932737*pkpm_lax_l[1]*G_1_cl[2]; 
  Ghat_G_1_u_l[4] = 0.1129384878631564*u_surf_lr[2]*G_1_lr[4]+0.1129384878631564*u_surf_cl[2]*G_1_lr[4]+0.2258769757263128*pkpm_lax_l[2]*G_1_lr[4]+0.1767766952966368*u_surf_lr[0]*G_1_lr[4]+0.1767766952966368*u_surf_cl[0]*G_1_lr[4]+0.3535533905932737*pkpm_lax_l[0]*G_1_lr[4]+0.1129384878631564*u_surf_lr[2]*G_1_cl[4]+0.1129384878631564*u_surf_cl[2]*G_1_cl[4]-0.2258769757263128*pkpm_lax_l[2]*G_1_cl[4]+0.1767766952966368*u_surf_lr[0]*G_1_cl[4]+0.1767766952966368*u_surf_cl[0]*G_1_cl[4]-0.3535533905932737*pkpm_lax_l[0]*G_1_cl[4]+0.1767766952966368*G_1_lr[0]*u_surf_lr[2]+0.1767766952966368*G_1_cl[0]*u_surf_lr[2]+0.1767766952966368*G_1_lr[0]*u_surf_cl[2]+0.1767766952966368*G_1_cl[0]*u_surf_cl[2]+0.3535533905932737*G_1_lr[0]*pkpm_lax_l[2]-0.3535533905932737*G_1_cl[0]*pkpm_lax_l[2]+0.1581138830084189*G_1_lr[1]*u_surf_lr[1]+0.1581138830084189*G_1_cl[1]*u_surf_lr[1]+0.1581138830084189*G_1_lr[1]*u_surf_cl[1]+0.1581138830084189*G_1_cl[1]*u_surf_cl[1]+0.3162277660168379*G_1_lr[1]*pkpm_lax_l[1]-0.3162277660168379*G_1_cl[1]*pkpm_lax_l[1]; 
  Ghat_G_1_u_l[5] = 0.1767766952966368*u_surf_lr[2]*G_1_lr[8]+0.1767766952966368*u_surf_cl[2]*G_1_lr[8]+0.3535533905932737*pkpm_lax_l[2]*G_1_lr[8]+0.1767766952966368*u_surf_lr[2]*G_1_cl[8]+0.1767766952966368*u_surf_cl[2]*G_1_cl[8]-0.3535533905932737*pkpm_lax_l[2]*G_1_cl[8]+0.1767766952966368*u_surf_lr[1]*G_1_lr[7]+0.1767766952966368*u_surf_cl[1]*G_1_lr[7]+0.3535533905932737*pkpm_lax_l[1]*G_1_lr[7]+0.1767766952966368*u_surf_lr[1]*G_1_cl[7]+0.1767766952966368*u_surf_cl[1]*G_1_cl[7]-0.3535533905932737*pkpm_lax_l[1]*G_1_cl[7]+0.1767766952966368*u_surf_lr[0]*G_1_lr[5]+0.1767766952966368*u_surf_cl[0]*G_1_lr[5]+0.3535533905932737*pkpm_lax_l[0]*G_1_lr[5]+0.1767766952966368*u_surf_lr[0]*G_1_cl[5]+0.1767766952966368*u_surf_cl[0]*G_1_cl[5]-0.3535533905932737*pkpm_lax_l[0]*G_1_cl[5]; 
  Ghat_G_1_u_l[6] = 0.1129384878631564*u_surf_lr[2]*G_1_lr[6]+0.1129384878631564*u_surf_cl[2]*G_1_lr[6]+0.2258769757263128*pkpm_lax_l[2]*G_1_lr[6]+0.1767766952966368*u_surf_lr[0]*G_1_lr[6]+0.1767766952966368*u_surf_cl[0]*G_1_lr[6]+0.3535533905932737*pkpm_lax_l[0]*G_1_lr[6]+0.1129384878631564*u_surf_lr[2]*G_1_cl[6]+0.1129384878631564*u_surf_cl[2]*G_1_cl[6]-0.2258769757263128*pkpm_lax_l[2]*G_1_cl[6]+0.1767766952966368*u_surf_lr[0]*G_1_cl[6]+0.1767766952966368*u_surf_cl[0]*G_1_cl[6]-0.3535533905932737*pkpm_lax_l[0]*G_1_cl[6]+0.1581138830084189*u_surf_lr[1]*G_1_lr[3]+0.1581138830084189*u_surf_cl[1]*G_1_lr[3]+0.3162277660168379*pkpm_lax_l[1]*G_1_lr[3]+0.1581138830084189*u_surf_lr[1]*G_1_cl[3]+0.1581138830084189*u_surf_cl[1]*G_1_cl[3]-0.3162277660168379*pkpm_lax_l[1]*G_1_cl[3]+0.1767766952966368*G_1_lr[2]*u_surf_lr[2]+0.1767766952966368*G_1_cl[2]*u_surf_lr[2]+0.1767766952966368*G_1_lr[2]*u_surf_cl[2]+0.1767766952966368*G_1_cl[2]*u_surf_cl[2]+0.3535533905932737*G_1_lr[2]*pkpm_lax_l[2]-0.3535533905932737*G_1_cl[2]*pkpm_lax_l[2]; 
  Ghat_G_1_u_l[7] = 0.1581138830084189*u_surf_lr[1]*G_1_lr[8]+0.1581138830084189*u_surf_cl[1]*G_1_lr[8]+0.3162277660168379*pkpm_lax_l[1]*G_1_lr[8]+0.1581138830084189*u_surf_lr[1]*G_1_cl[8]+0.1581138830084189*u_surf_cl[1]*G_1_cl[8]-0.3162277660168379*pkpm_lax_l[1]*G_1_cl[8]+0.1581138830084189*u_surf_lr[2]*G_1_lr[7]+0.1581138830084189*u_surf_cl[2]*G_1_lr[7]+0.3162277660168379*pkpm_lax_l[2]*G_1_lr[7]+0.1767766952966368*u_surf_lr[0]*G_1_lr[7]+0.1767766952966368*u_surf_cl[0]*G_1_lr[7]+0.3535533905932737*pkpm_lax_l[0]*G_1_lr[7]+0.1581138830084189*u_surf_lr[2]*G_1_cl[7]+0.1581138830084189*u_surf_cl[2]*G_1_cl[7]-0.3162277660168379*pkpm_lax_l[2]*G_1_cl[7]+0.1767766952966368*u_surf_lr[0]*G_1_cl[7]+0.1767766952966368*u_surf_cl[0]*G_1_cl[7]-0.3535533905932737*pkpm_lax_l[0]*G_1_cl[7]+0.1767766952966368*u_surf_lr[1]*G_1_lr[5]+0.1767766952966368*u_surf_cl[1]*G_1_lr[5]+0.3535533905932737*pkpm_lax_l[1]*G_1_lr[5]+0.1767766952966368*u_surf_lr[1]*G_1_cl[5]+0.1767766952966368*u_surf_cl[1]*G_1_cl[5]-0.3535533905932737*pkpm_lax_l[1]*G_1_cl[5]; 
  Ghat_G_1_u_l[8] = 0.1129384878631564*u_surf_lr[2]*G_1_lr[8]+0.1129384878631564*u_surf_cl[2]*G_1_lr[8]+0.2258769757263128*pkpm_lax_l[2]*G_1_lr[8]+0.1767766952966368*u_surf_lr[0]*G_1_lr[8]+0.1767766952966368*u_surf_cl[0]*G_1_lr[8]+0.3535533905932737*pkpm_lax_l[0]*G_1_lr[8]+0.1129384878631564*u_surf_lr[2]*G_1_cl[8]+0.1129384878631564*u_surf_cl[2]*G_1_cl[8]-0.2258769757263128*pkpm_lax_l[2]*G_1_cl[8]+0.1767766952966368*u_surf_lr[0]*G_1_cl[8]+0.1767766952966368*u_surf_cl[0]*G_1_cl[8]-0.3535533905932737*pkpm_lax_l[0]*G_1_cl[8]+0.1581138830084189*u_surf_lr[1]*G_1_lr[7]+0.1581138830084189*u_surf_cl[1]*G_1_lr[7]+0.3162277660168379*pkpm_lax_l[1]*G_1_lr[7]+0.1581138830084189*u_surf_lr[1]*G_1_cl[7]+0.1581138830084189*u_surf_cl[1]*G_1_cl[7]-0.3162277660168379*pkpm_lax_l[1]*G_1_cl[7]+0.1767766952966368*u_surf_lr[2]*G_1_lr[5]+0.1767766952966368*u_surf_cl[2]*G_1_lr[5]+0.3535533905932737*pkpm_lax_l[2]*G_1_lr[5]+0.1767766952966368*u_surf_lr[2]*G_1_cl[5]+0.1767766952966368*u_surf_cl[2]*G_1_cl[5]-0.3535533905932737*pkpm_lax_l[2]*G_1_cl[5]; 

  Ghat_F_0_u_r[0] = 0.1767766952966368*u_surf_rl[2]*F_0_rl[4]+0.1767766952966368*u_surf_cr[2]*F_0_rl[4]-0.3535533905932737*pkpm_lax_r[2]*F_0_rl[4]+0.1767766952966368*u_surf_rl[2]*F_0_cr[4]+0.1767766952966368*u_surf_cr[2]*F_0_cr[4]+0.3535533905932737*pkpm_lax_r[2]*F_0_cr[4]+0.1767766952966368*F_0_rl[1]*u_surf_rl[1]+0.1767766952966368*F_0_cr[1]*u_surf_rl[1]+0.1767766952966368*F_0_rl[1]*u_surf_cr[1]+0.1767766952966368*F_0_cr[1]*u_surf_cr[1]-0.3535533905932737*F_0_rl[1]*pkpm_lax_r[1]+0.3535533905932737*F_0_cr[1]*pkpm_lax_r[1]+0.1767766952966368*F_0_rl[0]*u_surf_rl[0]+0.1767766952966368*F_0_cr[0]*u_surf_rl[0]+0.1767766952966368*F_0_rl[0]*u_surf_cr[0]+0.1767766952966368*F_0_cr[0]*u_surf_cr[0]-0.3535533905932737*F_0_rl[0]*pkpm_lax_r[0]+0.3535533905932737*F_0_cr[0]*pkpm_lax_r[0]; 
  Ghat_F_0_u_r[1] = 0.1581138830084189*u_surf_rl[1]*F_0_rl[4]+0.1581138830084189*u_surf_cr[1]*F_0_rl[4]-0.3162277660168379*pkpm_lax_r[1]*F_0_rl[4]+0.1581138830084189*u_surf_rl[1]*F_0_cr[4]+0.1581138830084189*u_surf_cr[1]*F_0_cr[4]+0.3162277660168379*pkpm_lax_r[1]*F_0_cr[4]+0.1581138830084189*F_0_rl[1]*u_surf_rl[2]+0.1581138830084189*F_0_cr[1]*u_surf_rl[2]+0.1581138830084189*F_0_rl[1]*u_surf_cr[2]+0.1581138830084189*F_0_cr[1]*u_surf_cr[2]-0.3162277660168379*F_0_rl[1]*pkpm_lax_r[2]+0.3162277660168379*F_0_cr[1]*pkpm_lax_r[2]+0.1767766952966368*F_0_rl[0]*u_surf_rl[1]+0.1767766952966368*F_0_cr[0]*u_surf_rl[1]+0.1767766952966368*F_0_rl[0]*u_surf_cr[1]+0.1767766952966368*F_0_cr[0]*u_surf_cr[1]-0.3535533905932737*F_0_rl[0]*pkpm_lax_r[1]+0.3535533905932737*F_0_cr[0]*pkpm_lax_r[1]+0.1767766952966368*u_surf_rl[0]*F_0_rl[1]+0.1767766952966368*u_surf_cr[0]*F_0_rl[1]-0.3535533905932737*pkpm_lax_r[0]*F_0_rl[1]+0.1767766952966368*u_surf_rl[0]*F_0_cr[1]+0.1767766952966368*u_surf_cr[0]*F_0_cr[1]+0.3535533905932737*pkpm_lax_r[0]*F_0_cr[1]; 
  Ghat_F_0_u_r[2] = 0.1767766952966368*u_surf_rl[2]*F_0_rl[6]+0.1767766952966368*u_surf_cr[2]*F_0_rl[6]-0.3535533905932737*pkpm_lax_r[2]*F_0_rl[6]+0.1767766952966368*u_surf_rl[2]*F_0_cr[6]+0.1767766952966368*u_surf_cr[2]*F_0_cr[6]+0.3535533905932737*pkpm_lax_r[2]*F_0_cr[6]+0.1767766952966368*u_surf_rl[1]*F_0_rl[3]+0.1767766952966368*u_surf_cr[1]*F_0_rl[3]-0.3535533905932737*pkpm_lax_r[1]*F_0_rl[3]+0.1767766952966368*u_surf_rl[1]*F_0_cr[3]+0.1767766952966368*u_surf_cr[1]*F_0_cr[3]+0.3535533905932737*pkpm_lax_r[1]*F_0_cr[3]+0.1767766952966368*u_surf_rl[0]*F_0_rl[2]+0.1767766952966368*u_surf_cr[0]*F_0_rl[2]-0.3535533905932737*pkpm_lax_r[0]*F_0_rl[2]+0.1767766952966368*u_surf_rl[0]*F_0_cr[2]+0.1767766952966368*u_surf_cr[0]*F_0_cr[2]+0.3535533905932737*pkpm_lax_r[0]*F_0_cr[2]; 
  Ghat_F_0_u_r[3] = 0.1581138830084189*u_surf_rl[1]*F_0_rl[6]+0.1581138830084189*u_surf_cr[1]*F_0_rl[6]-0.3162277660168379*pkpm_lax_r[1]*F_0_rl[6]+0.1581138830084189*u_surf_rl[1]*F_0_cr[6]+0.1581138830084189*u_surf_cr[1]*F_0_cr[6]+0.3162277660168379*pkpm_lax_r[1]*F_0_cr[6]+0.1581138830084189*u_surf_rl[2]*F_0_rl[3]+0.1581138830084189*u_surf_cr[2]*F_0_rl[3]-0.3162277660168379*pkpm_lax_r[2]*F_0_rl[3]+0.1767766952966368*u_surf_rl[0]*F_0_rl[3]+0.1767766952966368*u_surf_cr[0]*F_0_rl[3]-0.3535533905932737*pkpm_lax_r[0]*F_0_rl[3]+0.1581138830084189*u_surf_rl[2]*F_0_cr[3]+0.1581138830084189*u_surf_cr[2]*F_0_cr[3]+0.3162277660168379*pkpm_lax_r[2]*F_0_cr[3]+0.1767766952966368*u_surf_rl[0]*F_0_cr[3]+0.1767766952966368*u_surf_cr[0]*F_0_cr[3]+0.3535533905932737*pkpm_lax_r[0]*F_0_cr[3]+0.1767766952966368*u_surf_rl[1]*F_0_rl[2]+0.1767766952966368*u_surf_cr[1]*F_0_rl[2]-0.3535533905932737*pkpm_lax_r[1]*F_0_rl[2]+0.1767766952966368*u_surf_rl[1]*F_0_cr[2]+0.1767766952966368*u_surf_cr[1]*F_0_cr[2]+0.3535533905932737*pkpm_lax_r[1]*F_0_cr[2]; 
  Ghat_F_0_u_r[4] = 0.1129384878631564*u_surf_rl[2]*F_0_rl[4]+0.1129384878631564*u_surf_cr[2]*F_0_rl[4]-0.2258769757263128*pkpm_lax_r[2]*F_0_rl[4]+0.1767766952966368*u_surf_rl[0]*F_0_rl[4]+0.1767766952966368*u_surf_cr[0]*F_0_rl[4]-0.3535533905932737*pkpm_lax_r[0]*F_0_rl[4]+0.1129384878631564*u_surf_rl[2]*F_0_cr[4]+0.1129384878631564*u_surf_cr[2]*F_0_cr[4]+0.2258769757263128*pkpm_lax_r[2]*F_0_cr[4]+0.1767766952966368*u_surf_rl[0]*F_0_cr[4]+0.1767766952966368*u_surf_cr[0]*F_0_cr[4]+0.3535533905932737*pkpm_lax_r[0]*F_0_cr[4]+0.1767766952966368*F_0_rl[0]*u_surf_rl[2]+0.1767766952966368*F_0_cr[0]*u_surf_rl[2]+0.1767766952966368*F_0_rl[0]*u_surf_cr[2]+0.1767766952966368*F_0_cr[0]*u_surf_cr[2]-0.3535533905932737*F_0_rl[0]*pkpm_lax_r[2]+0.3535533905932737*F_0_cr[0]*pkpm_lax_r[2]+0.1581138830084189*F_0_rl[1]*u_surf_rl[1]+0.1581138830084189*F_0_cr[1]*u_surf_rl[1]+0.1581138830084189*F_0_rl[1]*u_surf_cr[1]+0.1581138830084189*F_0_cr[1]*u_surf_cr[1]-0.3162277660168379*F_0_rl[1]*pkpm_lax_r[1]+0.3162277660168379*F_0_cr[1]*pkpm_lax_r[1]; 
  Ghat_F_0_u_r[5] = 0.1767766952966368*u_surf_rl[2]*F_0_rl[8]+0.1767766952966368*u_surf_cr[2]*F_0_rl[8]-0.3535533905932737*pkpm_lax_r[2]*F_0_rl[8]+0.1767766952966368*u_surf_rl[2]*F_0_cr[8]+0.1767766952966368*u_surf_cr[2]*F_0_cr[8]+0.3535533905932737*pkpm_lax_r[2]*F_0_cr[8]+0.1767766952966368*u_surf_rl[1]*F_0_rl[7]+0.1767766952966368*u_surf_cr[1]*F_0_rl[7]-0.3535533905932737*pkpm_lax_r[1]*F_0_rl[7]+0.1767766952966368*u_surf_rl[1]*F_0_cr[7]+0.1767766952966368*u_surf_cr[1]*F_0_cr[7]+0.3535533905932737*pkpm_lax_r[1]*F_0_cr[7]+0.1767766952966368*u_surf_rl[0]*F_0_rl[5]+0.1767766952966368*u_surf_cr[0]*F_0_rl[5]-0.3535533905932737*pkpm_lax_r[0]*F_0_rl[5]+0.1767766952966368*u_surf_rl[0]*F_0_cr[5]+0.1767766952966368*u_surf_cr[0]*F_0_cr[5]+0.3535533905932737*pkpm_lax_r[0]*F_0_cr[5]; 
  Ghat_F_0_u_r[6] = 0.1129384878631564*u_surf_rl[2]*F_0_rl[6]+0.1129384878631564*u_surf_cr[2]*F_0_rl[6]-0.2258769757263128*pkpm_lax_r[2]*F_0_rl[6]+0.1767766952966368*u_surf_rl[0]*F_0_rl[6]+0.1767766952966368*u_surf_cr[0]*F_0_rl[6]-0.3535533905932737*pkpm_lax_r[0]*F_0_rl[6]+0.1129384878631564*u_surf_rl[2]*F_0_cr[6]+0.1129384878631564*u_surf_cr[2]*F_0_cr[6]+0.2258769757263128*pkpm_lax_r[2]*F_0_cr[6]+0.1767766952966368*u_surf_rl[0]*F_0_cr[6]+0.1767766952966368*u_surf_cr[0]*F_0_cr[6]+0.3535533905932737*pkpm_lax_r[0]*F_0_cr[6]+0.1581138830084189*u_surf_rl[1]*F_0_rl[3]+0.1581138830084189*u_surf_cr[1]*F_0_rl[3]-0.3162277660168379*pkpm_lax_r[1]*F_0_rl[3]+0.1581138830084189*u_surf_rl[1]*F_0_cr[3]+0.1581138830084189*u_surf_cr[1]*F_0_cr[3]+0.3162277660168379*pkpm_lax_r[1]*F_0_cr[3]+0.1767766952966368*F_0_rl[2]*u_surf_rl[2]+0.1767766952966368*F_0_cr[2]*u_surf_rl[2]+0.1767766952966368*F_0_rl[2]*u_surf_cr[2]+0.1767766952966368*F_0_cr[2]*u_surf_cr[2]-0.3535533905932737*F_0_rl[2]*pkpm_lax_r[2]+0.3535533905932737*F_0_cr[2]*pkpm_lax_r[2]; 
  Ghat_F_0_u_r[7] = 0.1581138830084189*u_surf_rl[1]*F_0_rl[8]+0.1581138830084189*u_surf_cr[1]*F_0_rl[8]-0.3162277660168379*pkpm_lax_r[1]*F_0_rl[8]+0.1581138830084189*u_surf_rl[1]*F_0_cr[8]+0.1581138830084189*u_surf_cr[1]*F_0_cr[8]+0.3162277660168379*pkpm_lax_r[1]*F_0_cr[8]+0.1581138830084189*u_surf_rl[2]*F_0_rl[7]+0.1581138830084189*u_surf_cr[2]*F_0_rl[7]-0.3162277660168379*pkpm_lax_r[2]*F_0_rl[7]+0.1767766952966368*u_surf_rl[0]*F_0_rl[7]+0.1767766952966368*u_surf_cr[0]*F_0_rl[7]-0.3535533905932737*pkpm_lax_r[0]*F_0_rl[7]+0.1581138830084189*u_surf_rl[2]*F_0_cr[7]+0.1581138830084189*u_surf_cr[2]*F_0_cr[7]+0.3162277660168379*pkpm_lax_r[2]*F_0_cr[7]+0.1767766952966368*u_surf_rl[0]*F_0_cr[7]+0.1767766952966368*u_surf_cr[0]*F_0_cr[7]+0.3535533905932737*pkpm_lax_r[0]*F_0_cr[7]+0.1767766952966368*u_surf_rl[1]*F_0_rl[5]+0.1767766952966368*u_surf_cr[1]*F_0_rl[5]-0.3535533905932737*pkpm_lax_r[1]*F_0_rl[5]+0.1767766952966368*u_surf_rl[1]*F_0_cr[5]+0.1767766952966368*u_surf_cr[1]*F_0_cr[5]+0.3535533905932737*pkpm_lax_r[1]*F_0_cr[5]; 
  Ghat_F_0_u_r[8] = 0.1129384878631564*u_surf_rl[2]*F_0_rl[8]+0.1129384878631564*u_surf_cr[2]*F_0_rl[8]-0.2258769757263128*pkpm_lax_r[2]*F_0_rl[8]+0.1767766952966368*u_surf_rl[0]*F_0_rl[8]+0.1767766952966368*u_surf_cr[0]*F_0_rl[8]-0.3535533905932737*pkpm_lax_r[0]*F_0_rl[8]+0.1129384878631564*u_surf_rl[2]*F_0_cr[8]+0.1129384878631564*u_surf_cr[2]*F_0_cr[8]+0.2258769757263128*pkpm_lax_r[2]*F_0_cr[8]+0.1767766952966368*u_surf_rl[0]*F_0_cr[8]+0.1767766952966368*u_surf_cr[0]*F_0_cr[8]+0.3535533905932737*pkpm_lax_r[0]*F_0_cr[8]+0.1581138830084189*u_surf_rl[1]*F_0_rl[7]+0.1581138830084189*u_surf_cr[1]*F_0_rl[7]-0.3162277660168379*pkpm_lax_r[1]*F_0_rl[7]+0.1581138830084189*u_surf_rl[1]*F_0_cr[7]+0.1581138830084189*u_surf_cr[1]*F_0_cr[7]+0.3162277660168379*pkpm_lax_r[1]*F_0_cr[7]+0.1767766952966368*u_surf_rl[2]*F_0_rl[5]+0.1767766952966368*u_surf_cr[2]*F_0_rl[5]-0.3535533905932737*pkpm_lax_r[2]*F_0_rl[5]+0.1767766952966368*u_surf_rl[2]*F_0_cr[5]+0.1767766952966368*u_surf_cr[2]*F_0_cr[5]+0.3535533905932737*pkpm_lax_r[2]*F_0_cr[5]; 
  Ghat_G_1_u_r[0] = 0.1767766952966368*u_surf_rl[2]*G_1_rl[4]+0.1767766952966368*u_surf_cr[2]*G_1_rl[4]-0.3535533905932737*pkpm_lax_r[2]*G_1_rl[4]+0.1767766952966368*u_surf_rl[2]*G_1_cr[4]+0.1767766952966368*u_surf_cr[2]*G_1_cr[4]+0.3535533905932737*pkpm_lax_r[2]*G_1_cr[4]+0.1767766952966368*G_1_rl[1]*u_surf_rl[1]+0.1767766952966368*G_1_cr[1]*u_surf_rl[1]+0.1767766952966368*G_1_rl[1]*u_surf_cr[1]+0.1767766952966368*G_1_cr[1]*u_surf_cr[1]-0.3535533905932737*G_1_rl[1]*pkpm_lax_r[1]+0.3535533905932737*G_1_cr[1]*pkpm_lax_r[1]+0.1767766952966368*G_1_rl[0]*u_surf_rl[0]+0.1767766952966368*G_1_cr[0]*u_surf_rl[0]+0.1767766952966368*G_1_rl[0]*u_surf_cr[0]+0.1767766952966368*G_1_cr[0]*u_surf_cr[0]-0.3535533905932737*G_1_rl[0]*pkpm_lax_r[0]+0.3535533905932737*G_1_cr[0]*pkpm_lax_r[0]; 
  Ghat_G_1_u_r[1] = 0.1581138830084189*u_surf_rl[1]*G_1_rl[4]+0.1581138830084189*u_surf_cr[1]*G_1_rl[4]-0.3162277660168379*pkpm_lax_r[1]*G_1_rl[4]+0.1581138830084189*u_surf_rl[1]*G_1_cr[4]+0.1581138830084189*u_surf_cr[1]*G_1_cr[4]+0.3162277660168379*pkpm_lax_r[1]*G_1_cr[4]+0.1581138830084189*G_1_rl[1]*u_surf_rl[2]+0.1581138830084189*G_1_cr[1]*u_surf_rl[2]+0.1581138830084189*G_1_rl[1]*u_surf_cr[2]+0.1581138830084189*G_1_cr[1]*u_surf_cr[2]-0.3162277660168379*G_1_rl[1]*pkpm_lax_r[2]+0.3162277660168379*G_1_cr[1]*pkpm_lax_r[2]+0.1767766952966368*G_1_rl[0]*u_surf_rl[1]+0.1767766952966368*G_1_cr[0]*u_surf_rl[1]+0.1767766952966368*G_1_rl[0]*u_surf_cr[1]+0.1767766952966368*G_1_cr[0]*u_surf_cr[1]-0.3535533905932737*G_1_rl[0]*pkpm_lax_r[1]+0.3535533905932737*G_1_cr[0]*pkpm_lax_r[1]+0.1767766952966368*u_surf_rl[0]*G_1_rl[1]+0.1767766952966368*u_surf_cr[0]*G_1_rl[1]-0.3535533905932737*pkpm_lax_r[0]*G_1_rl[1]+0.1767766952966368*u_surf_rl[0]*G_1_cr[1]+0.1767766952966368*u_surf_cr[0]*G_1_cr[1]+0.3535533905932737*pkpm_lax_r[0]*G_1_cr[1]; 
  Ghat_G_1_u_r[2] = 0.1767766952966368*u_surf_rl[2]*G_1_rl[6]+0.1767766952966368*u_surf_cr[2]*G_1_rl[6]-0.3535533905932737*pkpm_lax_r[2]*G_1_rl[6]+0.1767766952966368*u_surf_rl[2]*G_1_cr[6]+0.1767766952966368*u_surf_cr[2]*G_1_cr[6]+0.3535533905932737*pkpm_lax_r[2]*G_1_cr[6]+0.1767766952966368*u_surf_rl[1]*G_1_rl[3]+0.1767766952966368*u_surf_cr[1]*G_1_rl[3]-0.3535533905932737*pkpm_lax_r[1]*G_1_rl[3]+0.1767766952966368*u_surf_rl[1]*G_1_cr[3]+0.1767766952966368*u_surf_cr[1]*G_1_cr[3]+0.3535533905932737*pkpm_lax_r[1]*G_1_cr[3]+0.1767766952966368*u_surf_rl[0]*G_1_rl[2]+0.1767766952966368*u_surf_cr[0]*G_1_rl[2]-0.3535533905932737*pkpm_lax_r[0]*G_1_rl[2]+0.1767766952966368*u_surf_rl[0]*G_1_cr[2]+0.1767766952966368*u_surf_cr[0]*G_1_cr[2]+0.3535533905932737*pkpm_lax_r[0]*G_1_cr[2]; 
  Ghat_G_1_u_r[3] = 0.1581138830084189*u_surf_rl[1]*G_1_rl[6]+0.1581138830084189*u_surf_cr[1]*G_1_rl[6]-0.3162277660168379*pkpm_lax_r[1]*G_1_rl[6]+0.1581138830084189*u_surf_rl[1]*G_1_cr[6]+0.1581138830084189*u_surf_cr[1]*G_1_cr[6]+0.3162277660168379*pkpm_lax_r[1]*G_1_cr[6]+0.1581138830084189*u_surf_rl[2]*G_1_rl[3]+0.1581138830084189*u_surf_cr[2]*G_1_rl[3]-0.3162277660168379*pkpm_lax_r[2]*G_1_rl[3]+0.1767766952966368*u_surf_rl[0]*G_1_rl[3]+0.1767766952966368*u_surf_cr[0]*G_1_rl[3]-0.3535533905932737*pkpm_lax_r[0]*G_1_rl[3]+0.1581138830084189*u_surf_rl[2]*G_1_cr[3]+0.1581138830084189*u_surf_cr[2]*G_1_cr[3]+0.3162277660168379*pkpm_lax_r[2]*G_1_cr[3]+0.1767766952966368*u_surf_rl[0]*G_1_cr[3]+0.1767766952966368*u_surf_cr[0]*G_1_cr[3]+0.3535533905932737*pkpm_lax_r[0]*G_1_cr[3]+0.1767766952966368*u_surf_rl[1]*G_1_rl[2]+0.1767766952966368*u_surf_cr[1]*G_1_rl[2]-0.3535533905932737*pkpm_lax_r[1]*G_1_rl[2]+0.1767766952966368*u_surf_rl[1]*G_1_cr[2]+0.1767766952966368*u_surf_cr[1]*G_1_cr[2]+0.3535533905932737*pkpm_lax_r[1]*G_1_cr[2]; 
  Ghat_G_1_u_r[4] = 0.1129384878631564*u_surf_rl[2]*G_1_rl[4]+0.1129384878631564*u_surf_cr[2]*G_1_rl[4]-0.2258769757263128*pkpm_lax_r[2]*G_1_rl[4]+0.1767766952966368*u_surf_rl[0]*G_1_rl[4]+0.1767766952966368*u_surf_cr[0]*G_1_rl[4]-0.3535533905932737*pkpm_lax_r[0]*G_1_rl[4]+0.1129384878631564*u_surf_rl[2]*G_1_cr[4]+0.1129384878631564*u_surf_cr[2]*G_1_cr[4]+0.2258769757263128*pkpm_lax_r[2]*G_1_cr[4]+0.1767766952966368*u_surf_rl[0]*G_1_cr[4]+0.1767766952966368*u_surf_cr[0]*G_1_cr[4]+0.3535533905932737*pkpm_lax_r[0]*G_1_cr[4]+0.1767766952966368*G_1_rl[0]*u_surf_rl[2]+0.1767766952966368*G_1_cr[0]*u_surf_rl[2]+0.1767766952966368*G_1_rl[0]*u_surf_cr[2]+0.1767766952966368*G_1_cr[0]*u_surf_cr[2]-0.3535533905932737*G_1_rl[0]*pkpm_lax_r[2]+0.3535533905932737*G_1_cr[0]*pkpm_lax_r[2]+0.1581138830084189*G_1_rl[1]*u_surf_rl[1]+0.1581138830084189*G_1_cr[1]*u_surf_rl[1]+0.1581138830084189*G_1_rl[1]*u_surf_cr[1]+0.1581138830084189*G_1_cr[1]*u_surf_cr[1]-0.3162277660168379*G_1_rl[1]*pkpm_lax_r[1]+0.3162277660168379*G_1_cr[1]*pkpm_lax_r[1]; 
  Ghat_G_1_u_r[5] = 0.1767766952966368*u_surf_rl[2]*G_1_rl[8]+0.1767766952966368*u_surf_cr[2]*G_1_rl[8]-0.3535533905932737*pkpm_lax_r[2]*G_1_rl[8]+0.1767766952966368*u_surf_rl[2]*G_1_cr[8]+0.1767766952966368*u_surf_cr[2]*G_1_cr[8]+0.3535533905932737*pkpm_lax_r[2]*G_1_cr[8]+0.1767766952966368*u_surf_rl[1]*G_1_rl[7]+0.1767766952966368*u_surf_cr[1]*G_1_rl[7]-0.3535533905932737*pkpm_lax_r[1]*G_1_rl[7]+0.1767766952966368*u_surf_rl[1]*G_1_cr[7]+0.1767766952966368*u_surf_cr[1]*G_1_cr[7]+0.3535533905932737*pkpm_lax_r[1]*G_1_cr[7]+0.1767766952966368*u_surf_rl[0]*G_1_rl[5]+0.1767766952966368*u_surf_cr[0]*G_1_rl[5]-0.3535533905932737*pkpm_lax_r[0]*G_1_rl[5]+0.1767766952966368*u_surf_rl[0]*G_1_cr[5]+0.1767766952966368*u_surf_cr[0]*G_1_cr[5]+0.3535533905932737*pkpm_lax_r[0]*G_1_cr[5]; 
  Ghat_G_1_u_r[6] = 0.1129384878631564*u_surf_rl[2]*G_1_rl[6]+0.1129384878631564*u_surf_cr[2]*G_1_rl[6]-0.2258769757263128*pkpm_lax_r[2]*G_1_rl[6]+0.1767766952966368*u_surf_rl[0]*G_1_rl[6]+0.1767766952966368*u_surf_cr[0]*G_1_rl[6]-0.3535533905932737*pkpm_lax_r[0]*G_1_rl[6]+0.1129384878631564*u_surf_rl[2]*G_1_cr[6]+0.1129384878631564*u_surf_cr[2]*G_1_cr[6]+0.2258769757263128*pkpm_lax_r[2]*G_1_cr[6]+0.1767766952966368*u_surf_rl[0]*G_1_cr[6]+0.1767766952966368*u_surf_cr[0]*G_1_cr[6]+0.3535533905932737*pkpm_lax_r[0]*G_1_cr[6]+0.1581138830084189*u_surf_rl[1]*G_1_rl[3]+0.1581138830084189*u_surf_cr[1]*G_1_rl[3]-0.3162277660168379*pkpm_lax_r[1]*G_1_rl[3]+0.1581138830084189*u_surf_rl[1]*G_1_cr[3]+0.1581138830084189*u_surf_cr[1]*G_1_cr[3]+0.3162277660168379*pkpm_lax_r[1]*G_1_cr[3]+0.1767766952966368*G_1_rl[2]*u_surf_rl[2]+0.1767766952966368*G_1_cr[2]*u_surf_rl[2]+0.1767766952966368*G_1_rl[2]*u_surf_cr[2]+0.1767766952966368*G_1_cr[2]*u_surf_cr[2]-0.3535533905932737*G_1_rl[2]*pkpm_lax_r[2]+0.3535533905932737*G_1_cr[2]*pkpm_lax_r[2]; 
  Ghat_G_1_u_r[7] = 0.1581138830084189*u_surf_rl[1]*G_1_rl[8]+0.1581138830084189*u_surf_cr[1]*G_1_rl[8]-0.3162277660168379*pkpm_lax_r[1]*G_1_rl[8]+0.1581138830084189*u_surf_rl[1]*G_1_cr[8]+0.1581138830084189*u_surf_cr[1]*G_1_cr[8]+0.3162277660168379*pkpm_lax_r[1]*G_1_cr[8]+0.1581138830084189*u_surf_rl[2]*G_1_rl[7]+0.1581138830084189*u_surf_cr[2]*G_1_rl[7]-0.3162277660168379*pkpm_lax_r[2]*G_1_rl[7]+0.1767766952966368*u_surf_rl[0]*G_1_rl[7]+0.1767766952966368*u_surf_cr[0]*G_1_rl[7]-0.3535533905932737*pkpm_lax_r[0]*G_1_rl[7]+0.1581138830084189*u_surf_rl[2]*G_1_cr[7]+0.1581138830084189*u_surf_cr[2]*G_1_cr[7]+0.3162277660168379*pkpm_lax_r[2]*G_1_cr[7]+0.1767766952966368*u_surf_rl[0]*G_1_cr[7]+0.1767766952966368*u_surf_cr[0]*G_1_cr[7]+0.3535533905932737*pkpm_lax_r[0]*G_1_cr[7]+0.1767766952966368*u_surf_rl[1]*G_1_rl[5]+0.1767766952966368*u_surf_cr[1]*G_1_rl[5]-0.3535533905932737*pkpm_lax_r[1]*G_1_rl[5]+0.1767766952966368*u_surf_rl[1]*G_1_cr[5]+0.1767766952966368*u_surf_cr[1]*G_1_cr[5]+0.3535533905932737*pkpm_lax_r[1]*G_1_cr[5]; 
  Ghat_G_1_u_r[8] = 0.1129384878631564*u_surf_rl[2]*G_1_rl[8]+0.1129384878631564*u_surf_cr[2]*G_1_rl[8]-0.2258769757263128*pkpm_lax_r[2]*G_1_rl[8]+0.1767766952966368*u_surf_rl[0]*G_1_rl[8]+0.1767766952966368*u_surf_cr[0]*G_1_rl[8]-0.3535533905932737*pkpm_lax_r[0]*G_1_rl[8]+0.1129384878631564*u_surf_rl[2]*G_1_cr[8]+0.1129384878631564*u_surf_cr[2]*G_1_cr[8]+0.2258769757263128*pkpm_lax_r[2]*G_1_cr[8]+0.1767766952966368*u_surf_rl[0]*G_1_cr[8]+0.1767766952966368*u_surf_cr[0]*G_1_cr[8]+0.3535533905932737*pkpm_lax_r[0]*G_1_cr[8]+0.1581138830084189*u_surf_rl[1]*G_1_rl[7]+0.1581138830084189*u_surf_cr[1]*G_1_rl[7]-0.3162277660168379*pkpm_lax_r[1]*G_1_rl[7]+0.1581138830084189*u_surf_rl[1]*G_1_cr[7]+0.1581138830084189*u_surf_cr[1]*G_1_cr[7]+0.3162277660168379*pkpm_lax_r[1]*G_1_cr[7]+0.1767766952966368*u_surf_rl[2]*G_1_rl[5]+0.1767766952966368*u_surf_cr[2]*G_1_rl[5]-0.3535533905932737*pkpm_lax_r[2]*G_1_rl[5]+0.1767766952966368*u_surf_rl[2]*G_1_cr[5]+0.1767766952966368*u_surf_cr[2]*G_1_cr[5]+0.3535533905932737*pkpm_lax_r[2]*G_1_cr[5]; 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_u_l[0]-0.7071067811865475*Ghat_F_0_u_r[0])*dx1; 
  out_F_0[1] += (0.7071067811865475*Ghat_F_0_u_l[1]-0.7071067811865475*Ghat_F_0_u_r[1])*dx1; 
  out_F_0[2] += -1.224744871391589*(Ghat_F_0_u_r[0]+Ghat_F_0_u_l[0])*dx1; 
  out_F_0[3] += (0.7071067811865475*Ghat_F_0_u_l[2]-0.7071067811865475*Ghat_F_0_u_r[2])*dx1; 
  out_F_0[4] += -1.224744871391589*(Ghat_F_0_u_r[1]+Ghat_F_0_u_l[1])*dx1; 
  out_F_0[5] += (0.7071067811865475*Ghat_F_0_u_l[3]-0.7071067811865475*Ghat_F_0_u_r[3])*dx1; 
  out_F_0[6] += -1.224744871391589*(Ghat_F_0_u_r[2]+Ghat_F_0_u_l[2])*dx1; 
  out_F_0[7] += (0.7071067811865475*Ghat_F_0_u_l[4]-0.7071067811865475*Ghat_F_0_u_r[4])*dx1; 
  out_F_0[8] += (1.58113883008419*Ghat_F_0_u_l[0]-1.58113883008419*Ghat_F_0_u_r[0])*dx1; 
  out_F_0[9] += (0.7071067811865475*Ghat_F_0_u_l[5]-0.7071067811865475*Ghat_F_0_u_r[5])*dx1; 
  out_F_0[10] += -1.224744871391589*(Ghat_F_0_u_r[3]+Ghat_F_0_u_l[3])*dx1; 
  out_F_0[11] += -1.224744871391589*(Ghat_F_0_u_r[4]+Ghat_F_0_u_l[4])*dx1; 
  out_F_0[12] += (1.58113883008419*Ghat_F_0_u_l[1]-1.58113883008419*Ghat_F_0_u_r[1])*dx1; 
  out_F_0[13] += (0.7071067811865475*Ghat_F_0_u_l[6]-0.7071067811865475*Ghat_F_0_u_r[6])*dx1; 
  out_F_0[14] += (1.58113883008419*Ghat_F_0_u_l[2]-1.58113883008419*Ghat_F_0_u_r[2])*dx1; 
  out_F_0[15] += (0.7071067811865475*Ghat_F_0_u_l[7]-0.7071067811865475*Ghat_F_0_u_r[7])*dx1; 
  out_F_0[16] += -1.224744871391589*(Ghat_F_0_u_r[5]+Ghat_F_0_u_l[5])*dx1; 
  out_F_0[17] += -1.224744871391589*(Ghat_F_0_u_r[6]+Ghat_F_0_u_l[6])*dx1; 
  out_F_0[18] += (1.58113883008419*Ghat_F_0_u_l[3]-1.58113883008419*Ghat_F_0_u_r[3])*dx1; 
  out_F_0[19] += -1.224744871391589*(Ghat_F_0_u_r[7]+Ghat_F_0_u_l[7])*dx1; 
  out_F_0[20] += (1.58113883008419*Ghat_F_0_u_l[4]-1.58113883008419*Ghat_F_0_u_r[4])*dx1; 
  out_F_0[21] += (0.7071067811865475*Ghat_F_0_u_l[8]-0.7071067811865475*Ghat_F_0_u_r[8])*dx1; 
  out_F_0[22] += (1.58113883008419*Ghat_F_0_u_l[5]-1.58113883008419*Ghat_F_0_u_r[5])*dx1; 
  out_F_0[23] += (1.58113883008419*Ghat_F_0_u_l[6]-1.58113883008419*Ghat_F_0_u_r[6])*dx1; 
  out_F_0[24] += -1.224744871391589*(Ghat_F_0_u_r[8]+Ghat_F_0_u_l[8])*dx1; 
  out_F_0[25] += (1.58113883008419*Ghat_F_0_u_l[7]-1.58113883008419*Ghat_F_0_u_r[7])*dx1; 
  out_F_0[26] += (1.58113883008419*Ghat_F_0_u_l[8]-1.58113883008419*Ghat_F_0_u_r[8])*dx1; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_u_l[0]-0.7071067811865475*Ghat_G_1_u_r[0])*dx1; 
  out_G_1[1] += (0.7071067811865475*Ghat_G_1_u_l[1]-0.7071067811865475*Ghat_G_1_u_r[1])*dx1; 
  out_G_1[2] += -1.224744871391589*(Ghat_G_1_u_r[0]+Ghat_G_1_u_l[0])*dx1; 
  out_G_1[3] += (0.7071067811865475*Ghat_G_1_u_l[2]-0.7071067811865475*Ghat_G_1_u_r[2])*dx1; 
  out_G_1[4] += -1.224744871391589*(Ghat_G_1_u_r[1]+Ghat_G_1_u_l[1])*dx1; 
  out_G_1[5] += (0.7071067811865475*Ghat_G_1_u_l[3]-0.7071067811865475*Ghat_G_1_u_r[3])*dx1; 
  out_G_1[6] += -1.224744871391589*(Ghat_G_1_u_r[2]+Ghat_G_1_u_l[2])*dx1; 
  out_G_1[7] += (0.7071067811865475*Ghat_G_1_u_l[4]-0.7071067811865475*Ghat_G_1_u_r[4])*dx1; 
  out_G_1[8] += (1.58113883008419*Ghat_G_1_u_l[0]-1.58113883008419*Ghat_G_1_u_r[0])*dx1; 
  out_G_1[9] += (0.7071067811865475*Ghat_G_1_u_l[5]-0.7071067811865475*Ghat_G_1_u_r[5])*dx1; 
  out_G_1[10] += -1.224744871391589*(Ghat_G_1_u_r[3]+Ghat_G_1_u_l[3])*dx1; 
  out_G_1[11] += -1.224744871391589*(Ghat_G_1_u_r[4]+Ghat_G_1_u_l[4])*dx1; 
  out_G_1[12] += (1.58113883008419*Ghat_G_1_u_l[1]-1.58113883008419*Ghat_G_1_u_r[1])*dx1; 
  out_G_1[13] += (0.7071067811865475*Ghat_G_1_u_l[6]-0.7071067811865475*Ghat_G_1_u_r[6])*dx1; 
  out_G_1[14] += (1.58113883008419*Ghat_G_1_u_l[2]-1.58113883008419*Ghat_G_1_u_r[2])*dx1; 
  out_G_1[15] += (0.7071067811865475*Ghat_G_1_u_l[7]-0.7071067811865475*Ghat_G_1_u_r[7])*dx1; 
  out_G_1[16] += -1.224744871391589*(Ghat_G_1_u_r[5]+Ghat_G_1_u_l[5])*dx1; 
  out_G_1[17] += -1.224744871391589*(Ghat_G_1_u_r[6]+Ghat_G_1_u_l[6])*dx1; 
  out_G_1[18] += (1.58113883008419*Ghat_G_1_u_l[3]-1.58113883008419*Ghat_G_1_u_r[3])*dx1; 
  out_G_1[19] += -1.224744871391589*(Ghat_G_1_u_r[7]+Ghat_G_1_u_l[7])*dx1; 
  out_G_1[20] += (1.58113883008419*Ghat_G_1_u_l[4]-1.58113883008419*Ghat_G_1_u_r[4])*dx1; 
  out_G_1[21] += (0.7071067811865475*Ghat_G_1_u_l[8]-0.7071067811865475*Ghat_G_1_u_r[8])*dx1; 
  out_G_1[22] += (1.58113883008419*Ghat_G_1_u_l[5]-1.58113883008419*Ghat_G_1_u_r[5])*dx1; 
  out_G_1[23] += (1.58113883008419*Ghat_G_1_u_l[6]-1.58113883008419*Ghat_G_1_u_r[6])*dx1; 
  out_G_1[24] += -1.224744871391589*(Ghat_G_1_u_r[8]+Ghat_G_1_u_l[8])*dx1; 
  out_G_1[25] += (1.58113883008419*Ghat_G_1_u_l[7]-1.58113883008419*Ghat_G_1_u_r[7])*dx1; 
  out_G_1[26] += (1.58113883008419*Ghat_G_1_u_l[8]-1.58113883008419*Ghat_G_1_u_r[8])*dx1; 

} 
