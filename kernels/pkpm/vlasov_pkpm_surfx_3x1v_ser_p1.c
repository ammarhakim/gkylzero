#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_hyb_3x1v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_hyb_3x1v_p1_upwind_quad_to_modal.h> 
#include <gkyl_basis_ser_3x_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_surfx_3x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *pkpm_priml, const double *pkpm_primc, const double *pkpm_primr, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                 Cell-center coordinates.
  // dxv[NDIM]:               Cell spacing.
  // bvarl/bvarc/bvarr:       Input magnetic field unit vector in left/center/right cells.
  // pkpm_priml/pkpm_primc/pkpm_primr: Input primitive variables in left/center/right cells.
  // fl/fc/fr:                Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells.
  // out:                     Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  const double *ul = &pkpm_priml[0]; 
  const double *uc = &pkpm_primc[0]; 
  const double *ur = &pkpm_primr[0]; 
  const double *bl = &bvarl[0]; 
  const double *bc = &bvarc[0]; 
  const double *br = &bvarr[0]; 
  // Get thermal velocity in direction of update for penalization vth^2 = 3.0*T_ii/m. 
  const double *vth_sql = &pkpm_priml[24]; 
  const double *vth_sqc = &pkpm_primc[24]; 
  const double *vth_sqr = &pkpm_primr[24]; 

  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[24]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[24]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[24]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[24]; 
  double alpha_l[24] = {0.0}; 
  double alpha_c[24] = {0.0}; 
  double alpha_r[24] = {0.0}; 
  alpha_l[0] = 1.414213562373095*bl[0]*wvpar; 
  alpha_l[1] = 1.414213562373095*bl[1]*wvpar; 
  alpha_l[2] = 1.414213562373095*bl[2]*wvpar; 
  alpha_l[3] = 1.414213562373095*bl[3]*wvpar; 
  alpha_l[4] = 0.408248290463863*bl[0]*dvpar; 
  alpha_l[5] = 1.414213562373095*bl[4]*wvpar; 
  alpha_l[6] = 1.414213562373095*bl[5]*wvpar; 
  alpha_l[7] = 1.414213562373095*bl[6]*wvpar; 
  alpha_l[8] = 0.408248290463863*bl[1]*dvpar; 
  alpha_l[9] = 0.408248290463863*bl[2]*dvpar; 
  alpha_l[10] = 0.408248290463863*bl[3]*dvpar; 
  alpha_l[11] = 1.414213562373095*bl[7]*wvpar; 
  alpha_l[12] = 0.408248290463863*bl[4]*dvpar; 
  alpha_l[13] = 0.408248290463863*bl[5]*dvpar; 
  alpha_l[14] = 0.408248290463863*bl[6]*dvpar; 
  alpha_l[15] = 0.408248290463863*bl[7]*dvpar; 

  alpha_c[0] = 1.414213562373095*bc[0]*wvpar; 
  alpha_c[1] = 1.414213562373095*bc[1]*wvpar; 
  alpha_c[2] = 1.414213562373095*bc[2]*wvpar; 
  alpha_c[3] = 1.414213562373095*bc[3]*wvpar; 
  alpha_c[4] = 0.408248290463863*bc[0]*dvpar; 
  alpha_c[5] = 1.414213562373095*bc[4]*wvpar; 
  alpha_c[6] = 1.414213562373095*bc[5]*wvpar; 
  alpha_c[7] = 1.414213562373095*bc[6]*wvpar; 
  alpha_c[8] = 0.408248290463863*bc[1]*dvpar; 
  alpha_c[9] = 0.408248290463863*bc[2]*dvpar; 
  alpha_c[10] = 0.408248290463863*bc[3]*dvpar; 
  alpha_c[11] = 1.414213562373095*bc[7]*wvpar; 
  alpha_c[12] = 0.408248290463863*bc[4]*dvpar; 
  alpha_c[13] = 0.408248290463863*bc[5]*dvpar; 
  alpha_c[14] = 0.408248290463863*bc[6]*dvpar; 
  alpha_c[15] = 0.408248290463863*bc[7]*dvpar; 

  alpha_r[0] = 1.414213562373095*br[0]*wvpar; 
  alpha_r[1] = 1.414213562373095*br[1]*wvpar; 
  alpha_r[2] = 1.414213562373095*br[2]*wvpar; 
  alpha_r[3] = 1.414213562373095*br[3]*wvpar; 
  alpha_r[4] = 0.408248290463863*br[0]*dvpar; 
  alpha_r[5] = 1.414213562373095*br[4]*wvpar; 
  alpha_r[6] = 1.414213562373095*br[5]*wvpar; 
  alpha_r[7] = 1.414213562373095*br[6]*wvpar; 
  alpha_r[8] = 0.408248290463863*br[1]*dvpar; 
  alpha_r[9] = 0.408248290463863*br[2]*dvpar; 
  alpha_r[10] = 0.408248290463863*br[3]*dvpar; 
  alpha_r[11] = 1.414213562373095*br[7]*wvpar; 
  alpha_r[12] = 0.408248290463863*br[4]*dvpar; 
  alpha_r[13] = 0.408248290463863*br[5]*dvpar; 
  alpha_r[14] = 0.408248290463863*br[6]*dvpar; 
  alpha_r[15] = 0.408248290463863*br[7]*dvpar; 

  double alphaSurf_l[12] = {0.0}; 
  alphaSurf_l[0] = 0.408248290463863*alpha_l[1]-0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_l[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_l[1] = 0.408248290463863*alpha_l[5]-0.408248290463863*alpha_c[5]+0.3535533905932737*alpha_l[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_l[2] = 0.408248290463863*alpha_l[6]-0.408248290463863*alpha_c[6]+0.3535533905932737*alpha_l[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_l[3] = 0.408248290463863*alpha_l[8]-0.408248290463863*alpha_c[8]+0.3535533905932737*alpha_l[4]+0.3535533905932737*alpha_c[4]; 
  alphaSurf_l[4] = 0.408248290463863*alpha_l[11]-0.408248290463863*alpha_c[11]+0.3535533905932737*alpha_l[7]+0.3535533905932737*alpha_c[7]; 
  alphaSurf_l[5] = 0.408248290463863*alpha_l[12]-0.408248290463863*alpha_c[12]+0.3535533905932737*alpha_l[9]+0.3535533905932737*alpha_c[9]; 
  alphaSurf_l[6] = 0.408248290463863*alpha_l[13]-0.408248290463863*alpha_c[13]+0.3535533905932737*alpha_l[10]+0.3535533905932737*alpha_c[10]; 
  alphaSurf_l[7] = 0.408248290463863*alpha_l[15]-0.408248290463863*alpha_c[15]+0.3535533905932737*alpha_l[14]+0.3535533905932737*alpha_c[14]; 

  double alphaSurf_r[12] = {0.0}; 
  alphaSurf_r[0] = (-0.408248290463863*alpha_r[1])+0.408248290463863*alpha_c[1]+0.3535533905932737*alpha_r[0]+0.3535533905932737*alpha_c[0]; 
  alphaSurf_r[1] = (-0.408248290463863*alpha_r[5])+0.408248290463863*alpha_c[5]+0.3535533905932737*alpha_r[2]+0.3535533905932737*alpha_c[2]; 
  alphaSurf_r[2] = (-0.408248290463863*alpha_r[6])+0.408248290463863*alpha_c[6]+0.3535533905932737*alpha_r[3]+0.3535533905932737*alpha_c[3]; 
  alphaSurf_r[3] = (-0.408248290463863*alpha_r[8])+0.408248290463863*alpha_c[8]+0.3535533905932737*alpha_r[4]+0.3535533905932737*alpha_c[4]; 
  alphaSurf_r[4] = (-0.408248290463863*alpha_r[11])+0.408248290463863*alpha_c[11]+0.3535533905932737*alpha_r[7]+0.3535533905932737*alpha_c[7]; 
  alphaSurf_r[5] = (-0.408248290463863*alpha_r[12])+0.408248290463863*alpha_c[12]+0.3535533905932737*alpha_r[9]+0.3535533905932737*alpha_c[9]; 
  alphaSurf_r[6] = (-0.408248290463863*alpha_r[13])+0.408248290463863*alpha_c[13]+0.3535533905932737*alpha_r[10]+0.3535533905932737*alpha_c[10]; 
  alphaSurf_r[7] = (-0.408248290463863*alpha_r[15])+0.408248290463863*alpha_c[15]+0.3535533905932737*alpha_r[14]+0.3535533905932737*alpha_c[14]; 

  double F_0_UpwindQuad_l[12] = {0.0};
  double F_0_UpwindQuad_r[12] = {0.0};
  double F_0_Upwind_l[12] = {0.0};
  double F_0_Upwind_r[12] = {0.0};
  double Ghat_F_0_l[12] = {0.0}; 
  double Ghat_F_0_r[12] = {0.0}; 
  double G_1_UpwindQuad_l[12] = {0.0};
  double G_1_UpwindQuad_r[12] = {0.0};
  double G_1_Upwind_l[12] = {0.0};
  double G_1_Upwind_r[12] = {0.0};
  double Ghat_G_1_l[12] = {0.0}; 
  double Ghat_G_1_r[12] = {0.0}; 

  if ((-0.4743416490252568*alphaSurf_l[7])+0.4743416490252568*(alphaSurf_l[6]+alphaSurf_l[5])+0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_r(F_0l); 
    G_1_UpwindQuad_l[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_l(F_0c); 
    G_1_UpwindQuad_l[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_l(G_1c); 
  } 
  if ((-0.4743416490252568*alphaSurf_r[7])+0.4743416490252568*(alphaSurf_r[6]+alphaSurf_r[5])+0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_r(F_0c); 
    G_1_UpwindQuad_r[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_l(F_0r); 
    G_1_UpwindQuad_r[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_l(G_1r); 
  } 
  if (0.3535533905932737*alphaSurf_l[4]-0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_r(F_0l); 
    G_1_UpwindQuad_l[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_l(F_0c); 
    G_1_UpwindQuad_l[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_l(G_1c); 
  } 
  if (0.3535533905932737*alphaSurf_r[4]-0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_r(F_0c); 
    G_1_UpwindQuad_r[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_l(F_0r); 
    G_1_UpwindQuad_r[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_l(G_1r); 
  } 
  if (0.4743416490252568*alphaSurf_l[7]-0.4743416490252568*(alphaSurf_l[6]+alphaSurf_l[5])+0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_r(F_0l); 
    G_1_UpwindQuad_l[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_l(F_0c); 
    G_1_UpwindQuad_l[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_l(G_1c); 
  } 
  if (0.4743416490252568*alphaSurf_r[7]-0.4743416490252568*(alphaSurf_r[6]+alphaSurf_r[5])+0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_r(F_0c); 
    G_1_UpwindQuad_r[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_l(F_0r); 
    G_1_UpwindQuad_r[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_l(G_1r); 
  } 
  if (0.4743416490252568*alphaSurf_l[7]-0.4743416490252568*alphaSurf_l[6]+0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[2]-0.3535533905932737*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_r(F_0l); 
    G_1_UpwindQuad_l[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_l(F_0c); 
    G_1_UpwindQuad_l[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_l(G_1c); 
  } 
  if (0.4743416490252568*alphaSurf_r[7]-0.4743416490252568*alphaSurf_r[6]+0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[2]-0.3535533905932737*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_r(F_0c); 
    G_1_UpwindQuad_r[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_l(F_0r); 
    G_1_UpwindQuad_r[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_l(G_1r); 
  } 
  if ((-0.3535533905932737*alphaSurf_l[4])+0.3535533905932737*alphaSurf_l[2]-0.3535533905932737*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_r(F_0l); 
    G_1_UpwindQuad_l[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_l(F_0c); 
    G_1_UpwindQuad_l[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_l(G_1c); 
  } 
  if ((-0.3535533905932737*alphaSurf_r[4])+0.3535533905932737*alphaSurf_r[2]-0.3535533905932737*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_r(F_0c); 
    G_1_UpwindQuad_r[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_l(F_0r); 
    G_1_UpwindQuad_r[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_l(G_1r); 
  } 
  if ((-0.4743416490252568*alphaSurf_l[7])+0.4743416490252568*alphaSurf_l[6]-0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[2]-0.3535533905932737*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_r(F_0l); 
    G_1_UpwindQuad_l[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_l(F_0c); 
    G_1_UpwindQuad_l[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_l(G_1c); 
  } 
  if ((-0.4743416490252568*alphaSurf_r[7])+0.4743416490252568*alphaSurf_r[6]-0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[2]-0.3535533905932737*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_r(F_0c); 
    G_1_UpwindQuad_r[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_l(F_0r); 
    G_1_UpwindQuad_r[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_l(G_1r); 
  } 
  if (0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6])-0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*alphaSurf_l[2]+0.3535533905932737*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_r(F_0l); 
    G_1_UpwindQuad_l[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_l(F_0c); 
    G_1_UpwindQuad_l[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_l(G_1c); 
  } 
  if (0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6])-0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*alphaSurf_r[2]+0.3535533905932737*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_r(F_0c); 
    G_1_UpwindQuad_r[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_l(F_0r); 
    G_1_UpwindQuad_r[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_l(G_1r); 
  } 
  if (0.3535533905932737*(alphaSurf_l[1]+alphaSurf_l[0])-0.3535533905932737*(alphaSurf_l[4]+alphaSurf_l[2]) > 0) { 
    F_0_UpwindQuad_l[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_r(F_0l); 
    G_1_UpwindQuad_l[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_l(F_0c); 
    G_1_UpwindQuad_l[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_l(G_1c); 
  } 
  if (0.3535533905932737*(alphaSurf_r[1]+alphaSurf_r[0])-0.3535533905932737*(alphaSurf_r[4]+alphaSurf_r[2]) > 0) { 
    F_0_UpwindQuad_r[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_r(F_0c); 
    G_1_UpwindQuad_r[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_l(F_0r); 
    G_1_UpwindQuad_r[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_l(G_1r); 
  } 
  if ((-0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6]))+0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*alphaSurf_l[2]+0.3535533905932737*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_r(F_0l); 
    G_1_UpwindQuad_l[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_l(F_0c); 
    G_1_UpwindQuad_l[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_l(G_1c); 
  } 
  if ((-0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6]))+0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*alphaSurf_r[2]+0.3535533905932737*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_r(F_0c); 
    G_1_UpwindQuad_r[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_l(F_0r); 
    G_1_UpwindQuad_r[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_l(G_1r); 
  } 
  if ((-0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6]+alphaSurf_l[5]))+0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_r(F_0l); 
    G_1_UpwindQuad_l[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_l(F_0c); 
    G_1_UpwindQuad_l[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_l(G_1c); 
  } 
  if ((-0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6]+alphaSurf_r[5]))+0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_r(F_0c); 
    G_1_UpwindQuad_r[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_l(F_0r); 
    G_1_UpwindQuad_r[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_l(G_1r); 
  } 
  if (0.3535533905932737*(alphaSurf_l[4]+alphaSurf_l[2]+alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_r(F_0l); 
    G_1_UpwindQuad_l[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_l(F_0c); 
    G_1_UpwindQuad_l[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_l(G_1c); 
  } 
  if (0.3535533905932737*(alphaSurf_r[4]+alphaSurf_r[2]+alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_r(F_0c); 
    G_1_UpwindQuad_r[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_l(F_0r); 
    G_1_UpwindQuad_r[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_l(G_1r); 
  } 
  if (0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6]+alphaSurf_l[5])+0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_r(F_0l); 
    G_1_UpwindQuad_l[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_r(G_1l); 
  } else { 
    F_0_UpwindQuad_l[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_l(F_0c); 
    G_1_UpwindQuad_l[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_l(G_1c); 
  } 
  if (0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6]+alphaSurf_r[5])+0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_r(F_0c); 
    G_1_UpwindQuad_r[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_r(G_1c); 
  } else { 
    F_0_UpwindQuad_r[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_l(F_0r); 
    G_1_UpwindQuad_r[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_l(G_1r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_3x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 
  hyb_3x1v_p1_xdir_upwind_quad_to_modal(G_1_UpwindQuad_l, G_1_Upwind_l); 
  hyb_3x1v_p1_xdir_upwind_quad_to_modal(G_1_UpwindQuad_r, G_1_Upwind_r); 

  Ghat_F_0_l[0] = 0.3535533905932737*F_0_Upwind_l[7]*alphaSurf_l[7]+0.3535533905932737*F_0_Upwind_l[6]*alphaSurf_l[6]+0.3535533905932737*F_0_Upwind_l[5]*alphaSurf_l[5]+0.3535533905932737*F_0_Upwind_l[4]*alphaSurf_l[4]+0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_l[3]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_l[2]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_l[1]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_F_0_l[1] = 0.3535533905932737*F_0_Upwind_l[6]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[6]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[3]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[2]*F_0_Upwind_l[4]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[1]; 
  Ghat_F_0_l[2] = 0.3535533905932737*F_0_Upwind_l[5]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[5]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[3]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[4]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[3] = 0.3162277660168379*alphaSurf_l[7]*F_0_Upwind_l[11]+0.3162277660168379*alphaSurf_l[6]*F_0_Upwind_l[10]+0.3162277660168379*alphaSurf_l[5]*F_0_Upwind_l[9]+0.3162277660168379*alphaSurf_l[3]*F_0_Upwind_l[8]+0.3535533905932737*F_0_Upwind_l[4]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[4]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[2]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[4] = 0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[3]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[5]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[5]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[4]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[2]; 
  Ghat_F_0_l[5] = 0.3162277660168379*alphaSurf_l[6]*F_0_Upwind_l[11]+0.3162277660168379*alphaSurf_l[7]*F_0_Upwind_l[10]+0.3162277660168379*alphaSurf_l[3]*F_0_Upwind_l[9]+0.3162277660168379*alphaSurf_l[5]*F_0_Upwind_l[8]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[2]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[4]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[4]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[6] = 0.3162277660168379*alphaSurf_l[5]*F_0_Upwind_l[11]+0.3162277660168379*alphaSurf_l[3]*F_0_Upwind_l[10]+0.3162277660168379*alphaSurf_l[7]*F_0_Upwind_l[9]+0.3162277660168379*alphaSurf_l[6]*F_0_Upwind_l[8]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[4]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[4]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[2]*F_0_Upwind_l[3]; 
  Ghat_F_0_l[7] = 0.3162277660168379*alphaSurf_l[3]*F_0_Upwind_l[11]+0.3162277660168379*alphaSurf_l[5]*F_0_Upwind_l[10]+0.3162277660168379*alphaSurf_l[6]*F_0_Upwind_l[9]+0.3162277660168379*alphaSurf_l[7]*F_0_Upwind_l[8]+0.3535533905932737*F_0_Upwind_l[0]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[7]+0.3535533905932737*F_0_Upwind_l[1]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[6]+0.3535533905932737*F_0_Upwind_l[2]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[2]*F_0_Upwind_l[5]+0.3535533905932737*F_0_Upwind_l[3]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[3]*F_0_Upwind_l[4]; 
  Ghat_F_0_l[8] = 0.3535533905932737*alphaSurf_l[4]*F_0_Upwind_l[11]+0.3535533905932737*alphaSurf_l[2]*F_0_Upwind_l[10]+0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[9]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[8]+0.3162277660168379*F_0_Upwind_l[7]*alphaSurf_l[7]+0.3162277660168379*F_0_Upwind_l[6]*alphaSurf_l[6]+0.3162277660168379*F_0_Upwind_l[5]*alphaSurf_l[5]+0.3162277660168379*F_0_Upwind_l[3]*alphaSurf_l[3]; 
  Ghat_F_0_l[9] = 0.3535533905932737*alphaSurf_l[2]*F_0_Upwind_l[11]+0.3535533905932737*alphaSurf_l[4]*F_0_Upwind_l[10]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[9]+0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[8]+0.3162277660168379*F_0_Upwind_l[6]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[6]*F_0_Upwind_l[7]+0.3162277660168379*F_0_Upwind_l[3]*alphaSurf_l[5]+0.3162277660168379*alphaSurf_l[3]*F_0_Upwind_l[5]; 
  Ghat_F_0_l[10] = 0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[11]+0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[10]+0.3535533905932737*alphaSurf_l[4]*F_0_Upwind_l[9]+0.3535533905932737*alphaSurf_l[2]*F_0_Upwind_l[8]+0.3162277660168379*F_0_Upwind_l[5]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[5]*F_0_Upwind_l[7]+0.3162277660168379*F_0_Upwind_l[3]*alphaSurf_l[6]+0.3162277660168379*alphaSurf_l[3]*F_0_Upwind_l[6]; 
  Ghat_F_0_l[11] = 0.3535533905932737*alphaSurf_l[0]*F_0_Upwind_l[11]+0.3535533905932737*alphaSurf_l[1]*F_0_Upwind_l[10]+0.3535533905932737*alphaSurf_l[2]*F_0_Upwind_l[9]+0.3535533905932737*alphaSurf_l[4]*F_0_Upwind_l[8]+0.3162277660168379*F_0_Upwind_l[3]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[3]*F_0_Upwind_l[7]+0.3162277660168379*F_0_Upwind_l[5]*alphaSurf_l[6]+0.3162277660168379*alphaSurf_l[5]*F_0_Upwind_l[6]; 
  Ghat_G_1_l[0] = 0.3535533905932737*G_1_Upwind_l[7]*alphaSurf_l[7]+0.3535533905932737*G_1_Upwind_l[6]*alphaSurf_l[6]+0.3535533905932737*G_1_Upwind_l[5]*alphaSurf_l[5]+0.3535533905932737*G_1_Upwind_l[4]*alphaSurf_l[4]+0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_l[3]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_l[2]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_l[1]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_l[0]; 
  Ghat_G_1_l[1] = 0.3535533905932737*G_1_Upwind_l[6]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[6]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[3]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[2]*G_1_Upwind_l[4]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[1]; 
  Ghat_G_1_l[2] = 0.3535533905932737*G_1_Upwind_l[5]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[5]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[3]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[4]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[3] = 0.3162277660168379*alphaSurf_l[7]*G_1_Upwind_l[11]+0.3162277660168379*alphaSurf_l[6]*G_1_Upwind_l[10]+0.3162277660168379*alphaSurf_l[5]*G_1_Upwind_l[9]+0.3162277660168379*alphaSurf_l[3]*G_1_Upwind_l[8]+0.3535533905932737*G_1_Upwind_l[4]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[4]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[2]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[4] = 0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[3]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[5]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[5]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[4]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_l[2]+0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[2]; 
  Ghat_G_1_l[5] = 0.3162277660168379*alphaSurf_l[6]*G_1_Upwind_l[11]+0.3162277660168379*alphaSurf_l[7]*G_1_Upwind_l[10]+0.3162277660168379*alphaSurf_l[3]*G_1_Upwind_l[9]+0.3162277660168379*alphaSurf_l[5]*G_1_Upwind_l[8]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[2]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[4]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[4]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[6] = 0.3162277660168379*alphaSurf_l[5]*G_1_Upwind_l[11]+0.3162277660168379*alphaSurf_l[3]*G_1_Upwind_l[10]+0.3162277660168379*alphaSurf_l[7]*G_1_Upwind_l[9]+0.3162277660168379*alphaSurf_l[6]*G_1_Upwind_l[8]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[4]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[4]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[2]*G_1_Upwind_l[3]; 
  Ghat_G_1_l[7] = 0.3162277660168379*alphaSurf_l[3]*G_1_Upwind_l[11]+0.3162277660168379*alphaSurf_l[5]*G_1_Upwind_l[10]+0.3162277660168379*alphaSurf_l[6]*G_1_Upwind_l[9]+0.3162277660168379*alphaSurf_l[7]*G_1_Upwind_l[8]+0.3535533905932737*G_1_Upwind_l[0]*alphaSurf_l[7]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[7]+0.3535533905932737*G_1_Upwind_l[1]*alphaSurf_l[6]+0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[6]+0.3535533905932737*G_1_Upwind_l[2]*alphaSurf_l[5]+0.3535533905932737*alphaSurf_l[2]*G_1_Upwind_l[5]+0.3535533905932737*G_1_Upwind_l[3]*alphaSurf_l[4]+0.3535533905932737*alphaSurf_l[3]*G_1_Upwind_l[4]; 
  Ghat_G_1_l[8] = 0.3535533905932737*alphaSurf_l[4]*G_1_Upwind_l[11]+0.3535533905932737*alphaSurf_l[2]*G_1_Upwind_l[10]+0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[9]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[8]+0.3162277660168379*G_1_Upwind_l[7]*alphaSurf_l[7]+0.3162277660168379*G_1_Upwind_l[6]*alphaSurf_l[6]+0.3162277660168379*G_1_Upwind_l[5]*alphaSurf_l[5]+0.3162277660168379*G_1_Upwind_l[3]*alphaSurf_l[3]; 
  Ghat_G_1_l[9] = 0.3535533905932737*alphaSurf_l[2]*G_1_Upwind_l[11]+0.3535533905932737*alphaSurf_l[4]*G_1_Upwind_l[10]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[9]+0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[8]+0.3162277660168379*G_1_Upwind_l[6]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[6]*G_1_Upwind_l[7]+0.3162277660168379*G_1_Upwind_l[3]*alphaSurf_l[5]+0.3162277660168379*alphaSurf_l[3]*G_1_Upwind_l[5]; 
  Ghat_G_1_l[10] = 0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[11]+0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[10]+0.3535533905932737*alphaSurf_l[4]*G_1_Upwind_l[9]+0.3535533905932737*alphaSurf_l[2]*G_1_Upwind_l[8]+0.3162277660168379*G_1_Upwind_l[5]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[5]*G_1_Upwind_l[7]+0.3162277660168379*G_1_Upwind_l[3]*alphaSurf_l[6]+0.3162277660168379*alphaSurf_l[3]*G_1_Upwind_l[6]; 
  Ghat_G_1_l[11] = 0.3535533905932737*alphaSurf_l[0]*G_1_Upwind_l[11]+0.3535533905932737*alphaSurf_l[1]*G_1_Upwind_l[10]+0.3535533905932737*alphaSurf_l[2]*G_1_Upwind_l[9]+0.3535533905932737*alphaSurf_l[4]*G_1_Upwind_l[8]+0.3162277660168379*G_1_Upwind_l[3]*alphaSurf_l[7]+0.3162277660168379*alphaSurf_l[3]*G_1_Upwind_l[7]+0.3162277660168379*G_1_Upwind_l[5]*alphaSurf_l[6]+0.3162277660168379*alphaSurf_l[5]*G_1_Upwind_l[6]; 

  Ghat_F_0_r[0] = 0.3535533905932737*F_0_Upwind_r[7]*alphaSurf_r[7]+0.3535533905932737*F_0_Upwind_r[6]*alphaSurf_r[6]+0.3535533905932737*F_0_Upwind_r[5]*alphaSurf_r[5]+0.3535533905932737*F_0_Upwind_r[4]*alphaSurf_r[4]+0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_r[3]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_r[2]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_r[1]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_F_0_r[1] = 0.3535533905932737*F_0_Upwind_r[6]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[6]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[3]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[2]*F_0_Upwind_r[4]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[1]; 
  Ghat_F_0_r[2] = 0.3535533905932737*F_0_Upwind_r[5]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[5]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[3]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[4]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[3] = 0.3162277660168379*alphaSurf_r[7]*F_0_Upwind_r[11]+0.3162277660168379*alphaSurf_r[6]*F_0_Upwind_r[10]+0.3162277660168379*alphaSurf_r[5]*F_0_Upwind_r[9]+0.3162277660168379*alphaSurf_r[3]*F_0_Upwind_r[8]+0.3535533905932737*F_0_Upwind_r[4]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[4]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[2]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[4] = 0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[3]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[5]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[5]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[4]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[2]; 
  Ghat_F_0_r[5] = 0.3162277660168379*alphaSurf_r[6]*F_0_Upwind_r[11]+0.3162277660168379*alphaSurf_r[7]*F_0_Upwind_r[10]+0.3162277660168379*alphaSurf_r[3]*F_0_Upwind_r[9]+0.3162277660168379*alphaSurf_r[5]*F_0_Upwind_r[8]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[2]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[4]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[4]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[6] = 0.3162277660168379*alphaSurf_r[5]*F_0_Upwind_r[11]+0.3162277660168379*alphaSurf_r[3]*F_0_Upwind_r[10]+0.3162277660168379*alphaSurf_r[7]*F_0_Upwind_r[9]+0.3162277660168379*alphaSurf_r[6]*F_0_Upwind_r[8]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[4]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[4]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[2]*F_0_Upwind_r[3]; 
  Ghat_F_0_r[7] = 0.3162277660168379*alphaSurf_r[3]*F_0_Upwind_r[11]+0.3162277660168379*alphaSurf_r[5]*F_0_Upwind_r[10]+0.3162277660168379*alphaSurf_r[6]*F_0_Upwind_r[9]+0.3162277660168379*alphaSurf_r[7]*F_0_Upwind_r[8]+0.3535533905932737*F_0_Upwind_r[0]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[7]+0.3535533905932737*F_0_Upwind_r[1]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[6]+0.3535533905932737*F_0_Upwind_r[2]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[2]*F_0_Upwind_r[5]+0.3535533905932737*F_0_Upwind_r[3]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[3]*F_0_Upwind_r[4]; 
  Ghat_F_0_r[8] = 0.3535533905932737*alphaSurf_r[4]*F_0_Upwind_r[11]+0.3535533905932737*alphaSurf_r[2]*F_0_Upwind_r[10]+0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[9]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[8]+0.3162277660168379*F_0_Upwind_r[7]*alphaSurf_r[7]+0.3162277660168379*F_0_Upwind_r[6]*alphaSurf_r[6]+0.3162277660168379*F_0_Upwind_r[5]*alphaSurf_r[5]+0.3162277660168379*F_0_Upwind_r[3]*alphaSurf_r[3]; 
  Ghat_F_0_r[9] = 0.3535533905932737*alphaSurf_r[2]*F_0_Upwind_r[11]+0.3535533905932737*alphaSurf_r[4]*F_0_Upwind_r[10]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[9]+0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[8]+0.3162277660168379*F_0_Upwind_r[6]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[6]*F_0_Upwind_r[7]+0.3162277660168379*F_0_Upwind_r[3]*alphaSurf_r[5]+0.3162277660168379*alphaSurf_r[3]*F_0_Upwind_r[5]; 
  Ghat_F_0_r[10] = 0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[11]+0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[10]+0.3535533905932737*alphaSurf_r[4]*F_0_Upwind_r[9]+0.3535533905932737*alphaSurf_r[2]*F_0_Upwind_r[8]+0.3162277660168379*F_0_Upwind_r[5]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[5]*F_0_Upwind_r[7]+0.3162277660168379*F_0_Upwind_r[3]*alphaSurf_r[6]+0.3162277660168379*alphaSurf_r[3]*F_0_Upwind_r[6]; 
  Ghat_F_0_r[11] = 0.3535533905932737*alphaSurf_r[0]*F_0_Upwind_r[11]+0.3535533905932737*alphaSurf_r[1]*F_0_Upwind_r[10]+0.3535533905932737*alphaSurf_r[2]*F_0_Upwind_r[9]+0.3535533905932737*alphaSurf_r[4]*F_0_Upwind_r[8]+0.3162277660168379*F_0_Upwind_r[3]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[3]*F_0_Upwind_r[7]+0.3162277660168379*F_0_Upwind_r[5]*alphaSurf_r[6]+0.3162277660168379*alphaSurf_r[5]*F_0_Upwind_r[6]; 
  Ghat_G_1_r[0] = 0.3535533905932737*G_1_Upwind_r[7]*alphaSurf_r[7]+0.3535533905932737*G_1_Upwind_r[6]*alphaSurf_r[6]+0.3535533905932737*G_1_Upwind_r[5]*alphaSurf_r[5]+0.3535533905932737*G_1_Upwind_r[4]*alphaSurf_r[4]+0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_r[3]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_r[2]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_r[1]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_r[0]; 
  Ghat_G_1_r[1] = 0.3535533905932737*G_1_Upwind_r[6]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[6]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[3]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[2]*G_1_Upwind_r[4]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[1]; 
  Ghat_G_1_r[2] = 0.3535533905932737*G_1_Upwind_r[5]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[5]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[3]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[4]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[3] = 0.3162277660168379*alphaSurf_r[7]*G_1_Upwind_r[11]+0.3162277660168379*alphaSurf_r[6]*G_1_Upwind_r[10]+0.3162277660168379*alphaSurf_r[5]*G_1_Upwind_r[9]+0.3162277660168379*alphaSurf_r[3]*G_1_Upwind_r[8]+0.3535533905932737*G_1_Upwind_r[4]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[4]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[2]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[4] = 0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[3]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[5]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[5]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[4]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_r[2]+0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[2]; 
  Ghat_G_1_r[5] = 0.3162277660168379*alphaSurf_r[6]*G_1_Upwind_r[11]+0.3162277660168379*alphaSurf_r[7]*G_1_Upwind_r[10]+0.3162277660168379*alphaSurf_r[3]*G_1_Upwind_r[9]+0.3162277660168379*alphaSurf_r[5]*G_1_Upwind_r[8]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[2]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[4]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[4]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[6] = 0.3162277660168379*alphaSurf_r[5]*G_1_Upwind_r[11]+0.3162277660168379*alphaSurf_r[3]*G_1_Upwind_r[10]+0.3162277660168379*alphaSurf_r[7]*G_1_Upwind_r[9]+0.3162277660168379*alphaSurf_r[6]*G_1_Upwind_r[8]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[4]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[4]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[2]*G_1_Upwind_r[3]; 
  Ghat_G_1_r[7] = 0.3162277660168379*alphaSurf_r[3]*G_1_Upwind_r[11]+0.3162277660168379*alphaSurf_r[5]*G_1_Upwind_r[10]+0.3162277660168379*alphaSurf_r[6]*G_1_Upwind_r[9]+0.3162277660168379*alphaSurf_r[7]*G_1_Upwind_r[8]+0.3535533905932737*G_1_Upwind_r[0]*alphaSurf_r[7]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[7]+0.3535533905932737*G_1_Upwind_r[1]*alphaSurf_r[6]+0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[6]+0.3535533905932737*G_1_Upwind_r[2]*alphaSurf_r[5]+0.3535533905932737*alphaSurf_r[2]*G_1_Upwind_r[5]+0.3535533905932737*G_1_Upwind_r[3]*alphaSurf_r[4]+0.3535533905932737*alphaSurf_r[3]*G_1_Upwind_r[4]; 
  Ghat_G_1_r[8] = 0.3535533905932737*alphaSurf_r[4]*G_1_Upwind_r[11]+0.3535533905932737*alphaSurf_r[2]*G_1_Upwind_r[10]+0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[9]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[8]+0.3162277660168379*G_1_Upwind_r[7]*alphaSurf_r[7]+0.3162277660168379*G_1_Upwind_r[6]*alphaSurf_r[6]+0.3162277660168379*G_1_Upwind_r[5]*alphaSurf_r[5]+0.3162277660168379*G_1_Upwind_r[3]*alphaSurf_r[3]; 
  Ghat_G_1_r[9] = 0.3535533905932737*alphaSurf_r[2]*G_1_Upwind_r[11]+0.3535533905932737*alphaSurf_r[4]*G_1_Upwind_r[10]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[9]+0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[8]+0.3162277660168379*G_1_Upwind_r[6]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[6]*G_1_Upwind_r[7]+0.3162277660168379*G_1_Upwind_r[3]*alphaSurf_r[5]+0.3162277660168379*alphaSurf_r[3]*G_1_Upwind_r[5]; 
  Ghat_G_1_r[10] = 0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[11]+0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[10]+0.3535533905932737*alphaSurf_r[4]*G_1_Upwind_r[9]+0.3535533905932737*alphaSurf_r[2]*G_1_Upwind_r[8]+0.3162277660168379*G_1_Upwind_r[5]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[5]*G_1_Upwind_r[7]+0.3162277660168379*G_1_Upwind_r[3]*alphaSurf_r[6]+0.3162277660168379*alphaSurf_r[3]*G_1_Upwind_r[6]; 
  Ghat_G_1_r[11] = 0.3535533905932737*alphaSurf_r[0]*G_1_Upwind_r[11]+0.3535533905932737*alphaSurf_r[1]*G_1_Upwind_r[10]+0.3535533905932737*alphaSurf_r[2]*G_1_Upwind_r[9]+0.3535533905932737*alphaSurf_r[4]*G_1_Upwind_r[8]+0.3162277660168379*G_1_Upwind_r[3]*alphaSurf_r[7]+0.3162277660168379*alphaSurf_r[3]*G_1_Upwind_r[7]+0.3162277660168379*G_1_Upwind_r[5]*alphaSurf_r[6]+0.3162277660168379*alphaSurf_r[5]*G_1_Upwind_r[6]; 

  out_F_0[0] += (0.7071067811865475*Ghat_F_0_l[0]-0.7071067811865475*Ghat_F_0_r[0])*dx1; 
  out_F_0[1] += -1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dx1; 
  out_F_0[2] += (0.7071067811865475*Ghat_F_0_l[1]-0.7071067811865475*Ghat_F_0_r[1])*dx1; 
  out_F_0[3] += (0.7071067811865475*Ghat_F_0_l[2]-0.7071067811865475*Ghat_F_0_r[2])*dx1; 
  out_F_0[4] += (0.7071067811865475*Ghat_F_0_l[3]-0.7071067811865475*Ghat_F_0_r[3])*dx1; 
  out_F_0[5] += -1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dx1; 
  out_F_0[6] += -1.224744871391589*(Ghat_F_0_r[2]+Ghat_F_0_l[2])*dx1; 
  out_F_0[7] += (0.7071067811865475*Ghat_F_0_l[4]-0.7071067811865475*Ghat_F_0_r[4])*dx1; 
  out_F_0[8] += -1.224744871391589*(Ghat_F_0_r[3]+Ghat_F_0_l[3])*dx1; 
  out_F_0[9] += (0.7071067811865475*Ghat_F_0_l[5]-0.7071067811865475*Ghat_F_0_r[5])*dx1; 
  out_F_0[10] += (0.7071067811865475*Ghat_F_0_l[6]-0.7071067811865475*Ghat_F_0_r[6])*dx1; 
  out_F_0[11] += -1.224744871391589*(Ghat_F_0_r[4]+Ghat_F_0_l[4])*dx1; 
  out_F_0[12] += -1.224744871391589*(Ghat_F_0_r[5]+Ghat_F_0_l[5])*dx1; 
  out_F_0[13] += -1.224744871391589*(Ghat_F_0_r[6]+Ghat_F_0_l[6])*dx1; 
  out_F_0[14] += (0.7071067811865475*Ghat_F_0_l[7]-0.7071067811865475*Ghat_F_0_r[7])*dx1; 
  out_F_0[15] += -1.224744871391589*(Ghat_F_0_r[7]+Ghat_F_0_l[7])*dx1; 
  out_F_0[16] += (0.7071067811865475*Ghat_F_0_l[8]-0.7071067811865475*Ghat_F_0_r[8])*dx1; 
  out_F_0[17] += -1.224744871391589*(Ghat_F_0_r[8]+Ghat_F_0_l[8])*dx1; 
  out_F_0[18] += (0.7071067811865475*Ghat_F_0_l[9]-0.7071067811865475*Ghat_F_0_r[9])*dx1; 
  out_F_0[19] += (0.7071067811865475*Ghat_F_0_l[10]-0.7071067811865475*Ghat_F_0_r[10])*dx1; 
  out_F_0[20] += -1.224744871391589*(Ghat_F_0_r[9]+Ghat_F_0_l[9])*dx1; 
  out_F_0[21] += -1.224744871391589*(Ghat_F_0_r[10]+Ghat_F_0_l[10])*dx1; 
  out_F_0[22] += (0.7071067811865475*Ghat_F_0_l[11]-0.7071067811865475*Ghat_F_0_r[11])*dx1; 
  out_F_0[23] += -1.224744871391589*(Ghat_F_0_r[11]+Ghat_F_0_l[11])*dx1; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_l[0]-0.7071067811865475*Ghat_G_1_r[0])*dx1; 
  out_G_1[1] += -1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dx1; 
  out_G_1[2] += (0.7071067811865475*Ghat_G_1_l[1]-0.7071067811865475*Ghat_G_1_r[1])*dx1; 
  out_G_1[3] += (0.7071067811865475*Ghat_G_1_l[2]-0.7071067811865475*Ghat_G_1_r[2])*dx1; 
  out_G_1[4] += (0.7071067811865475*Ghat_G_1_l[3]-0.7071067811865475*Ghat_G_1_r[3])*dx1; 
  out_G_1[5] += -1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dx1; 
  out_G_1[6] += -1.224744871391589*(Ghat_G_1_r[2]+Ghat_G_1_l[2])*dx1; 
  out_G_1[7] += (0.7071067811865475*Ghat_G_1_l[4]-0.7071067811865475*Ghat_G_1_r[4])*dx1; 
  out_G_1[8] += -1.224744871391589*(Ghat_G_1_r[3]+Ghat_G_1_l[3])*dx1; 
  out_G_1[9] += (0.7071067811865475*Ghat_G_1_l[5]-0.7071067811865475*Ghat_G_1_r[5])*dx1; 
  out_G_1[10] += (0.7071067811865475*Ghat_G_1_l[6]-0.7071067811865475*Ghat_G_1_r[6])*dx1; 
  out_G_1[11] += -1.224744871391589*(Ghat_G_1_r[4]+Ghat_G_1_l[4])*dx1; 
  out_G_1[12] += -1.224744871391589*(Ghat_G_1_r[5]+Ghat_G_1_l[5])*dx1; 
  out_G_1[13] += -1.224744871391589*(Ghat_G_1_r[6]+Ghat_G_1_l[6])*dx1; 
  out_G_1[14] += (0.7071067811865475*Ghat_G_1_l[7]-0.7071067811865475*Ghat_G_1_r[7])*dx1; 
  out_G_1[15] += -1.224744871391589*(Ghat_G_1_r[7]+Ghat_G_1_l[7])*dx1; 
  out_G_1[16] += (0.7071067811865475*Ghat_G_1_l[8]-0.7071067811865475*Ghat_G_1_r[8])*dx1; 
  out_G_1[17] += -1.224744871391589*(Ghat_G_1_r[8]+Ghat_G_1_l[8])*dx1; 
  out_G_1[18] += (0.7071067811865475*Ghat_G_1_l[9]-0.7071067811865475*Ghat_G_1_r[9])*dx1; 
  out_G_1[19] += (0.7071067811865475*Ghat_G_1_l[10]-0.7071067811865475*Ghat_G_1_r[10])*dx1; 
  out_G_1[20] += -1.224744871391589*(Ghat_G_1_r[9]+Ghat_G_1_l[9])*dx1; 
  out_G_1[21] += -1.224744871391589*(Ghat_G_1_r[10]+Ghat_G_1_l[10])*dx1; 
  out_G_1[22] += (0.7071067811865475*Ghat_G_1_l[11]-0.7071067811865475*Ghat_G_1_r[11])*dx1; 
  out_G_1[23] += -1.224744871391589*(Ghat_G_1_r[11]+Ghat_G_1_l[11])*dx1; 

  double F_0l_rx[4] = {0.0}; 
  double F_0c_lx[4] = {0.0}; 
  double F_0c_rx[4] = {0.0}; 
  double F_0r_lx[4] = {0.0}; 

  double G_1l_rx[4] = {0.0}; 
  double G_1c_lx[4] = {0.0}; 
  double G_1c_rx[4] = {0.0}; 
  double G_1r_lx[4] = {0.0}; 

  double avg_F_0_l[4] = {0.0}; 
  double avg_F_0_r[4] = {0.0}; 
  double avg_G_1_l[4] = {0.0}; 
  double avg_G_1_r[4] = {0.0}; 
  double avg_u_l[4] = {0.0}; 
  double avg_u_r[4] = {0.0}; 
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

  double max_speed_quad_l[4] = {0.0}; 
  double max_speed_quad_r[4] = {0.0}; 
  double max_speed_modal_l[4] = {0.0}; 
  double max_speed_modal_r[4] = {0.0}; 
  avg_u_l[0] = 0.6123724356957944*ul[1]-0.6123724356957944*uc[1]+0.3535533905932737*ul[0]+0.3535533905932737*uc[0]; 
  avg_u_l[1] = 0.6123724356957944*ul[4]-0.6123724356957944*uc[4]+0.3535533905932737*ul[2]+0.3535533905932737*uc[2]; 
  avg_u_l[2] = 0.6123724356957944*ul[5]-0.6123724356957944*uc[5]+0.3535533905932737*ul[3]+0.3535533905932737*uc[3]; 
  avg_u_l[3] = 0.6123724356957944*ul[7]-0.6123724356957944*uc[7]+0.3535533905932737*ul[6]+0.3535533905932737*uc[6]; 

  avg_u_r[0] = (-0.6123724356957944*ur[1])+0.6123724356957944*uc[1]+0.3535533905932737*ur[0]+0.3535533905932737*uc[0]; 
  avg_u_r[1] = (-0.6123724356957944*ur[4])+0.6123724356957944*uc[4]+0.3535533905932737*ur[2]+0.3535533905932737*uc[2]; 
  avg_u_r[2] = (-0.6123724356957944*ur[5])+0.6123724356957944*uc[5]+0.3535533905932737*ur[3]+0.3535533905932737*uc[3]; 
  avg_u_r[3] = (-0.6123724356957944*ur[7])+0.6123724356957944*uc[7]+0.3535533905932737*ur[6]+0.3535533905932737*uc[6]; 

  ul_r = ser_3x_p1_surfx1_eval_quad_node_0_r(ul); 
  uc_l = ser_3x_p1_surfx1_eval_quad_node_0_l(uc); 
  uc_r = ser_3x_p1_surfx1_eval_quad_node_0_r(uc); 
  ur_l = ser_3x_p1_surfx1_eval_quad_node_0_l(ur); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_3x_p1_surfx1_eval_quad_node_0_r(vth_sql); 
  vth_sqc_l = ser_3x_p1_surfx1_eval_quad_node_0_l(vth_sqc); 
  vth_sqc_r = ser_3x_p1_surfx1_eval_quad_node_0_r(vth_sqc); 
  vth_sqr_l = ser_3x_p1_surfx1_eval_quad_node_0_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[0] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[0] = uQuad_r + vthQuad_r; 
  ul_r = ser_3x_p1_surfx1_eval_quad_node_1_r(ul); 
  uc_l = ser_3x_p1_surfx1_eval_quad_node_1_l(uc); 
  uc_r = ser_3x_p1_surfx1_eval_quad_node_1_r(uc); 
  ur_l = ser_3x_p1_surfx1_eval_quad_node_1_l(ur); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_3x_p1_surfx1_eval_quad_node_1_r(vth_sql); 
  vth_sqc_l = ser_3x_p1_surfx1_eval_quad_node_1_l(vth_sqc); 
  vth_sqc_r = ser_3x_p1_surfx1_eval_quad_node_1_r(vth_sqc); 
  vth_sqr_l = ser_3x_p1_surfx1_eval_quad_node_1_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[1] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[1] = uQuad_r + vthQuad_r; 
  ul_r = ser_3x_p1_surfx1_eval_quad_node_2_r(ul); 
  uc_l = ser_3x_p1_surfx1_eval_quad_node_2_l(uc); 
  uc_r = ser_3x_p1_surfx1_eval_quad_node_2_r(uc); 
  ur_l = ser_3x_p1_surfx1_eval_quad_node_2_l(ur); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_3x_p1_surfx1_eval_quad_node_2_r(vth_sql); 
  vth_sqc_l = ser_3x_p1_surfx1_eval_quad_node_2_l(vth_sqc); 
  vth_sqc_r = ser_3x_p1_surfx1_eval_quad_node_2_r(vth_sqc); 
  vth_sqr_l = ser_3x_p1_surfx1_eval_quad_node_2_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[2] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[2] = uQuad_r + vthQuad_r; 
  ul_r = ser_3x_p1_surfx1_eval_quad_node_3_r(ul); 
  uc_l = ser_3x_p1_surfx1_eval_quad_node_3_l(uc); 
  uc_r = ser_3x_p1_surfx1_eval_quad_node_3_r(uc); 
  ur_l = ser_3x_p1_surfx1_eval_quad_node_3_l(ur); 
  uQuad_l = fmax(fabs(ul_r), fabs(uc_l)); 
  uQuad_r = fmax(fabs(uc_r), fabs(ur_l)); 
  vth_sql_r = ser_3x_p1_surfx1_eval_quad_node_3_r(vth_sql); 
  vth_sqc_l = ser_3x_p1_surfx1_eval_quad_node_3_l(vth_sqc); 
  vth_sqc_r = ser_3x_p1_surfx1_eval_quad_node_3_r(vth_sqc); 
  vth_sqr_l = ser_3x_p1_surfx1_eval_quad_node_3_l(vth_sqr); 
  vthQuad_l = fmax(sqrt(fabs(vth_sql_r)), sqrt(fabs(vth_sqc_l))); 
  vthQuad_r = fmax(sqrt(fabs(vth_sqc_r)), sqrt(fabs(vth_sqr_l))); 
  max_speed_quad_l[3] = uQuad_l + vthQuad_l; 
  max_speed_quad_r[3] = uQuad_r + vthQuad_r; 
  ser_3x_p1_upwind_quad_to_modal(max_speed_quad_l, max_speed_modal_l); 
  ser_3x_p1_upwind_quad_to_modal(max_speed_quad_r, max_speed_modal_r); 
  F_0l_rx[0] = 0.4330127018922193*F_0l[1]+0.25*F_0l[0]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[5]+0.25*F_0l[2]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[6]+0.25*F_0l[3]; 
  F_0l_rx[3] = 0.4330127018922193*F_0l[11]+0.25*F_0l[7]; 

  F_0c_lx[0] = 0.25*F_0c[0]-0.4330127018922193*F_0c[1]; 
  F_0c_lx[1] = 0.25*F_0c[2]-0.4330127018922193*F_0c[5]; 
  F_0c_lx[2] = 0.25*F_0c[3]-0.4330127018922193*F_0c[6]; 
  F_0c_lx[3] = 0.25*F_0c[7]-0.4330127018922193*F_0c[11]; 

  F_0c_rx[0] = 0.4330127018922193*F_0c[1]+0.25*F_0c[0]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[5]+0.25*F_0c[2]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[6]+0.25*F_0c[3]; 
  F_0c_rx[3] = 0.4330127018922193*F_0c[11]+0.25*F_0c[7]; 

  F_0r_lx[0] = 0.25*F_0r[0]-0.4330127018922193*F_0r[1]; 
  F_0r_lx[1] = 0.25*F_0r[2]-0.4330127018922193*F_0r[5]; 
  F_0r_lx[2] = 0.25*F_0r[3]-0.4330127018922193*F_0r[6]; 
  F_0r_lx[3] = 0.25*F_0r[7]-0.4330127018922193*F_0r[11]; 

  G_1l_rx[0] = 0.4330127018922193*G_1l[1]+0.25*G_1l[0]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[5]+0.25*G_1l[2]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[6]+0.25*G_1l[3]; 
  G_1l_rx[3] = 0.4330127018922193*G_1l[11]+0.25*G_1l[7]; 

  G_1c_lx[0] = 0.25*G_1c[0]-0.4330127018922193*G_1c[1]; 
  G_1c_lx[1] = 0.25*G_1c[2]-0.4330127018922193*G_1c[5]; 
  G_1c_lx[2] = 0.25*G_1c[3]-0.4330127018922193*G_1c[6]; 
  G_1c_lx[3] = 0.25*G_1c[7]-0.4330127018922193*G_1c[11]; 

  G_1c_rx[0] = 0.4330127018922193*G_1c[1]+0.25*G_1c[0]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[5]+0.25*G_1c[2]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[6]+0.25*G_1c[3]; 
  G_1c_rx[3] = 0.4330127018922193*G_1c[11]+0.25*G_1c[7]; 

  G_1r_lx[0] = 0.25*G_1r[0]-0.4330127018922193*G_1r[1]; 
  G_1r_lx[1] = 0.25*G_1r[2]-0.4330127018922193*G_1r[5]; 
  G_1r_lx[2] = 0.25*G_1r[3]-0.4330127018922193*G_1r[6]; 
  G_1r_lx[3] = 0.25*G_1r[7]-0.4330127018922193*G_1r[11]; 

  out_F_0[0] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[0] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[1])-0.4330127018922193*F_0l[0]; 
  F_0l_rx[1] = (-0.75*F_0l[5])-0.4330127018922193*F_0l[2]; 
  F_0l_rx[2] = (-0.75*F_0l[6])-0.4330127018922193*F_0l[3]; 
  F_0l_rx[3] = (-0.75*F_0l[11])-0.4330127018922193*F_0l[7]; 

  F_0c_lx[0] = 0.75*F_0c[1]-0.4330127018922193*F_0c[0]; 
  F_0c_lx[1] = 0.75*F_0c[5]-0.4330127018922193*F_0c[2]; 
  F_0c_lx[2] = 0.75*F_0c[6]-0.4330127018922193*F_0c[3]; 
  F_0c_lx[3] = 0.75*F_0c[11]-0.4330127018922193*F_0c[7]; 

  F_0c_rx[0] = 0.75*F_0c[1]+0.4330127018922193*F_0c[0]; 
  F_0c_rx[1] = 0.75*F_0c[5]+0.4330127018922193*F_0c[2]; 
  F_0c_rx[2] = 0.75*F_0c[6]+0.4330127018922193*F_0c[3]; 
  F_0c_rx[3] = 0.75*F_0c[11]+0.4330127018922193*F_0c[7]; 

  F_0r_lx[0] = 0.4330127018922193*F_0r[0]-0.75*F_0r[1]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[2]-0.75*F_0r[5]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[3]-0.75*F_0r[6]; 
  F_0r_lx[3] = 0.4330127018922193*F_0r[7]-0.75*F_0r[11]; 

  G_1l_rx[0] = (-0.75*G_1l[1])-0.4330127018922193*G_1l[0]; 
  G_1l_rx[1] = (-0.75*G_1l[5])-0.4330127018922193*G_1l[2]; 
  G_1l_rx[2] = (-0.75*G_1l[6])-0.4330127018922193*G_1l[3]; 
  G_1l_rx[3] = (-0.75*G_1l[11])-0.4330127018922193*G_1l[7]; 

  G_1c_lx[0] = 0.75*G_1c[1]-0.4330127018922193*G_1c[0]; 
  G_1c_lx[1] = 0.75*G_1c[5]-0.4330127018922193*G_1c[2]; 
  G_1c_lx[2] = 0.75*G_1c[6]-0.4330127018922193*G_1c[3]; 
  G_1c_lx[3] = 0.75*G_1c[11]-0.4330127018922193*G_1c[7]; 

  G_1c_rx[0] = 0.75*G_1c[1]+0.4330127018922193*G_1c[0]; 
  G_1c_rx[1] = 0.75*G_1c[5]+0.4330127018922193*G_1c[2]; 
  G_1c_rx[2] = 0.75*G_1c[6]+0.4330127018922193*G_1c[3]; 
  G_1c_rx[3] = 0.75*G_1c[11]+0.4330127018922193*G_1c[7]; 

  G_1r_lx[0] = 0.4330127018922193*G_1r[0]-0.75*G_1r[1]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[2]-0.75*G_1r[5]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[3]-0.75*G_1r[6]; 
  G_1r_lx[3] = 0.4330127018922193*G_1r[7]-0.75*G_1r[11]; 

  out_F_0[1] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[1] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922193*F_0l[5]+0.25*F_0l[2]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[1]+0.25*F_0l[0]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[11]+0.25*F_0l[7]; 
  F_0l_rx[3] = 0.4330127018922193*F_0l[6]+0.25*F_0l[3]; 

  F_0c_lx[0] = 0.25*F_0c[2]-0.4330127018922193*F_0c[5]; 
  F_0c_lx[1] = 0.25*F_0c[0]-0.4330127018922193*F_0c[1]; 
  F_0c_lx[2] = 0.25*F_0c[7]-0.4330127018922193*F_0c[11]; 
  F_0c_lx[3] = 0.25*F_0c[3]-0.4330127018922193*F_0c[6]; 

  F_0c_rx[0] = 0.4330127018922193*F_0c[5]+0.25*F_0c[2]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[1]+0.25*F_0c[0]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[11]+0.25*F_0c[7]; 
  F_0c_rx[3] = 0.4330127018922193*F_0c[6]+0.25*F_0c[3]; 

  F_0r_lx[0] = 0.25*F_0r[2]-0.4330127018922193*F_0r[5]; 
  F_0r_lx[1] = 0.25*F_0r[0]-0.4330127018922193*F_0r[1]; 
  F_0r_lx[2] = 0.25*F_0r[7]-0.4330127018922193*F_0r[11]; 
  F_0r_lx[3] = 0.25*F_0r[3]-0.4330127018922193*F_0r[6]; 

  G_1l_rx[0] = 0.4330127018922193*G_1l[5]+0.25*G_1l[2]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[1]+0.25*G_1l[0]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[11]+0.25*G_1l[7]; 
  G_1l_rx[3] = 0.4330127018922193*G_1l[6]+0.25*G_1l[3]; 

  G_1c_lx[0] = 0.25*G_1c[2]-0.4330127018922193*G_1c[5]; 
  G_1c_lx[1] = 0.25*G_1c[0]-0.4330127018922193*G_1c[1]; 
  G_1c_lx[2] = 0.25*G_1c[7]-0.4330127018922193*G_1c[11]; 
  G_1c_lx[3] = 0.25*G_1c[3]-0.4330127018922193*G_1c[6]; 

  G_1c_rx[0] = 0.4330127018922193*G_1c[5]+0.25*G_1c[2]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[1]+0.25*G_1c[0]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[11]+0.25*G_1c[7]; 
  G_1c_rx[3] = 0.4330127018922193*G_1c[6]+0.25*G_1c[3]; 

  G_1r_lx[0] = 0.25*G_1r[2]-0.4330127018922193*G_1r[5]; 
  G_1r_lx[1] = 0.25*G_1r[0]-0.4330127018922193*G_1r[1]; 
  G_1r_lx[2] = 0.25*G_1r[7]-0.4330127018922193*G_1r[11]; 
  G_1r_lx[3] = 0.25*G_1r[3]-0.4330127018922193*G_1r[6]; 

  out_F_0[2] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[2] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922193*F_0l[6]+0.25*F_0l[3]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[11]+0.25*F_0l[7]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[1]+0.25*F_0l[0]; 
  F_0l_rx[3] = 0.4330127018922193*F_0l[5]+0.25*F_0l[2]; 

  F_0c_lx[0] = 0.25*F_0c[3]-0.4330127018922193*F_0c[6]; 
  F_0c_lx[1] = 0.25*F_0c[7]-0.4330127018922193*F_0c[11]; 
  F_0c_lx[2] = 0.25*F_0c[0]-0.4330127018922193*F_0c[1]; 
  F_0c_lx[3] = 0.25*F_0c[2]-0.4330127018922193*F_0c[5]; 

  F_0c_rx[0] = 0.4330127018922193*F_0c[6]+0.25*F_0c[3]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[11]+0.25*F_0c[7]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[1]+0.25*F_0c[0]; 
  F_0c_rx[3] = 0.4330127018922193*F_0c[5]+0.25*F_0c[2]; 

  F_0r_lx[0] = 0.25*F_0r[3]-0.4330127018922193*F_0r[6]; 
  F_0r_lx[1] = 0.25*F_0r[7]-0.4330127018922193*F_0r[11]; 
  F_0r_lx[2] = 0.25*F_0r[0]-0.4330127018922193*F_0r[1]; 
  F_0r_lx[3] = 0.25*F_0r[2]-0.4330127018922193*F_0r[5]; 

  G_1l_rx[0] = 0.4330127018922193*G_1l[6]+0.25*G_1l[3]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[11]+0.25*G_1l[7]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[1]+0.25*G_1l[0]; 
  G_1l_rx[3] = 0.4330127018922193*G_1l[5]+0.25*G_1l[2]; 

  G_1c_lx[0] = 0.25*G_1c[3]-0.4330127018922193*G_1c[6]; 
  G_1c_lx[1] = 0.25*G_1c[7]-0.4330127018922193*G_1c[11]; 
  G_1c_lx[2] = 0.25*G_1c[0]-0.4330127018922193*G_1c[1]; 
  G_1c_lx[3] = 0.25*G_1c[2]-0.4330127018922193*G_1c[5]; 

  G_1c_rx[0] = 0.4330127018922193*G_1c[6]+0.25*G_1c[3]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[11]+0.25*G_1c[7]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[1]+0.25*G_1c[0]; 
  G_1c_rx[3] = 0.4330127018922193*G_1c[5]+0.25*G_1c[2]; 

  G_1r_lx[0] = 0.25*G_1r[3]-0.4330127018922193*G_1r[6]; 
  G_1r_lx[1] = 0.25*G_1r[7]-0.4330127018922193*G_1r[11]; 
  G_1r_lx[2] = 0.25*G_1r[0]-0.4330127018922193*G_1r[1]; 
  G_1r_lx[3] = 0.25*G_1r[2]-0.4330127018922193*G_1r[5]; 

  out_F_0[3] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[3] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922193*F_0l[8]+0.25*F_0l[4]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[12]+0.25*F_0l[9]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[13]+0.25*F_0l[10]; 
  F_0l_rx[3] = 0.4330127018922193*F_0l[15]+0.25*F_0l[14]; 

  F_0c_lx[0] = 0.25*F_0c[4]-0.4330127018922193*F_0c[8]; 
  F_0c_lx[1] = 0.25*F_0c[9]-0.4330127018922193*F_0c[12]; 
  F_0c_lx[2] = 0.25*F_0c[10]-0.4330127018922193*F_0c[13]; 
  F_0c_lx[3] = 0.25*F_0c[14]-0.4330127018922193*F_0c[15]; 

  F_0c_rx[0] = 0.4330127018922193*F_0c[8]+0.25*F_0c[4]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[12]+0.25*F_0c[9]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[13]+0.25*F_0c[10]; 
  F_0c_rx[3] = 0.4330127018922193*F_0c[15]+0.25*F_0c[14]; 

  F_0r_lx[0] = 0.25*F_0r[4]-0.4330127018922193*F_0r[8]; 
  F_0r_lx[1] = 0.25*F_0r[9]-0.4330127018922193*F_0r[12]; 
  F_0r_lx[2] = 0.25*F_0r[10]-0.4330127018922193*F_0r[13]; 
  F_0r_lx[3] = 0.25*F_0r[14]-0.4330127018922193*F_0r[15]; 

  G_1l_rx[0] = 0.4330127018922193*G_1l[8]+0.25*G_1l[4]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[12]+0.25*G_1l[9]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[13]+0.25*G_1l[10]; 
  G_1l_rx[3] = 0.4330127018922193*G_1l[15]+0.25*G_1l[14]; 

  G_1c_lx[0] = 0.25*G_1c[4]-0.4330127018922193*G_1c[8]; 
  G_1c_lx[1] = 0.25*G_1c[9]-0.4330127018922193*G_1c[12]; 
  G_1c_lx[2] = 0.25*G_1c[10]-0.4330127018922193*G_1c[13]; 
  G_1c_lx[3] = 0.25*G_1c[14]-0.4330127018922193*G_1c[15]; 

  G_1c_rx[0] = 0.4330127018922193*G_1c[8]+0.25*G_1c[4]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[12]+0.25*G_1c[9]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[13]+0.25*G_1c[10]; 
  G_1c_rx[3] = 0.4330127018922193*G_1c[15]+0.25*G_1c[14]; 

  G_1r_lx[0] = 0.25*G_1r[4]-0.4330127018922193*G_1r[8]; 
  G_1r_lx[1] = 0.25*G_1r[9]-0.4330127018922193*G_1r[12]; 
  G_1r_lx[2] = 0.25*G_1r[10]-0.4330127018922193*G_1r[13]; 
  G_1r_lx[3] = 0.25*G_1r[14]-0.4330127018922193*G_1r[15]; 

  out_F_0[4] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[4] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[5])-0.4330127018922193*F_0l[2]; 
  F_0l_rx[1] = (-0.75*F_0l[1])-0.4330127018922193*F_0l[0]; 
  F_0l_rx[2] = (-0.75*F_0l[11])-0.4330127018922193*F_0l[7]; 
  F_0l_rx[3] = (-0.75*F_0l[6])-0.4330127018922193*F_0l[3]; 

  F_0c_lx[0] = 0.75*F_0c[5]-0.4330127018922193*F_0c[2]; 
  F_0c_lx[1] = 0.75*F_0c[1]-0.4330127018922193*F_0c[0]; 
  F_0c_lx[2] = 0.75*F_0c[11]-0.4330127018922193*F_0c[7]; 
  F_0c_lx[3] = 0.75*F_0c[6]-0.4330127018922193*F_0c[3]; 

  F_0c_rx[0] = 0.75*F_0c[5]+0.4330127018922193*F_0c[2]; 
  F_0c_rx[1] = 0.75*F_0c[1]+0.4330127018922193*F_0c[0]; 
  F_0c_rx[2] = 0.75*F_0c[11]+0.4330127018922193*F_0c[7]; 
  F_0c_rx[3] = 0.75*F_0c[6]+0.4330127018922193*F_0c[3]; 

  F_0r_lx[0] = 0.4330127018922193*F_0r[2]-0.75*F_0r[5]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[0]-0.75*F_0r[1]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[7]-0.75*F_0r[11]; 
  F_0r_lx[3] = 0.4330127018922193*F_0r[3]-0.75*F_0r[6]; 

  G_1l_rx[0] = (-0.75*G_1l[5])-0.4330127018922193*G_1l[2]; 
  G_1l_rx[1] = (-0.75*G_1l[1])-0.4330127018922193*G_1l[0]; 
  G_1l_rx[2] = (-0.75*G_1l[11])-0.4330127018922193*G_1l[7]; 
  G_1l_rx[3] = (-0.75*G_1l[6])-0.4330127018922193*G_1l[3]; 

  G_1c_lx[0] = 0.75*G_1c[5]-0.4330127018922193*G_1c[2]; 
  G_1c_lx[1] = 0.75*G_1c[1]-0.4330127018922193*G_1c[0]; 
  G_1c_lx[2] = 0.75*G_1c[11]-0.4330127018922193*G_1c[7]; 
  G_1c_lx[3] = 0.75*G_1c[6]-0.4330127018922193*G_1c[3]; 

  G_1c_rx[0] = 0.75*G_1c[5]+0.4330127018922193*G_1c[2]; 
  G_1c_rx[1] = 0.75*G_1c[1]+0.4330127018922193*G_1c[0]; 
  G_1c_rx[2] = 0.75*G_1c[11]+0.4330127018922193*G_1c[7]; 
  G_1c_rx[3] = 0.75*G_1c[6]+0.4330127018922193*G_1c[3]; 

  G_1r_lx[0] = 0.4330127018922193*G_1r[2]-0.75*G_1r[5]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[0]-0.75*G_1r[1]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[7]-0.75*G_1r[11]; 
  G_1r_lx[3] = 0.4330127018922193*G_1r[3]-0.75*G_1r[6]; 

  out_F_0[5] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[5] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[6])-0.4330127018922193*F_0l[3]; 
  F_0l_rx[1] = (-0.75*F_0l[11])-0.4330127018922193*F_0l[7]; 
  F_0l_rx[2] = (-0.75*F_0l[1])-0.4330127018922193*F_0l[0]; 
  F_0l_rx[3] = (-0.75*F_0l[5])-0.4330127018922193*F_0l[2]; 

  F_0c_lx[0] = 0.75*F_0c[6]-0.4330127018922193*F_0c[3]; 
  F_0c_lx[1] = 0.75*F_0c[11]-0.4330127018922193*F_0c[7]; 
  F_0c_lx[2] = 0.75*F_0c[1]-0.4330127018922193*F_0c[0]; 
  F_0c_lx[3] = 0.75*F_0c[5]-0.4330127018922193*F_0c[2]; 

  F_0c_rx[0] = 0.75*F_0c[6]+0.4330127018922193*F_0c[3]; 
  F_0c_rx[1] = 0.75*F_0c[11]+0.4330127018922193*F_0c[7]; 
  F_0c_rx[2] = 0.75*F_0c[1]+0.4330127018922193*F_0c[0]; 
  F_0c_rx[3] = 0.75*F_0c[5]+0.4330127018922193*F_0c[2]; 

  F_0r_lx[0] = 0.4330127018922193*F_0r[3]-0.75*F_0r[6]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[7]-0.75*F_0r[11]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[0]-0.75*F_0r[1]; 
  F_0r_lx[3] = 0.4330127018922193*F_0r[2]-0.75*F_0r[5]; 

  G_1l_rx[0] = (-0.75*G_1l[6])-0.4330127018922193*G_1l[3]; 
  G_1l_rx[1] = (-0.75*G_1l[11])-0.4330127018922193*G_1l[7]; 
  G_1l_rx[2] = (-0.75*G_1l[1])-0.4330127018922193*G_1l[0]; 
  G_1l_rx[3] = (-0.75*G_1l[5])-0.4330127018922193*G_1l[2]; 

  G_1c_lx[0] = 0.75*G_1c[6]-0.4330127018922193*G_1c[3]; 
  G_1c_lx[1] = 0.75*G_1c[11]-0.4330127018922193*G_1c[7]; 
  G_1c_lx[2] = 0.75*G_1c[1]-0.4330127018922193*G_1c[0]; 
  G_1c_lx[3] = 0.75*G_1c[5]-0.4330127018922193*G_1c[2]; 

  G_1c_rx[0] = 0.75*G_1c[6]+0.4330127018922193*G_1c[3]; 
  G_1c_rx[1] = 0.75*G_1c[11]+0.4330127018922193*G_1c[7]; 
  G_1c_rx[2] = 0.75*G_1c[1]+0.4330127018922193*G_1c[0]; 
  G_1c_rx[3] = 0.75*G_1c[5]+0.4330127018922193*G_1c[2]; 

  G_1r_lx[0] = 0.4330127018922193*G_1r[3]-0.75*G_1r[6]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[7]-0.75*G_1r[11]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[0]-0.75*G_1r[1]; 
  G_1r_lx[3] = 0.4330127018922193*G_1r[2]-0.75*G_1r[5]; 

  out_F_0[6] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[6] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922193*F_0l[11]+0.25*F_0l[7]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[6]+0.25*F_0l[3]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[5]+0.25*F_0l[2]; 
  F_0l_rx[3] = 0.4330127018922193*F_0l[1]+0.25*F_0l[0]; 

  F_0c_lx[0] = 0.25*F_0c[7]-0.4330127018922193*F_0c[11]; 
  F_0c_lx[1] = 0.25*F_0c[3]-0.4330127018922193*F_0c[6]; 
  F_0c_lx[2] = 0.25*F_0c[2]-0.4330127018922193*F_0c[5]; 
  F_0c_lx[3] = 0.25*F_0c[0]-0.4330127018922193*F_0c[1]; 

  F_0c_rx[0] = 0.4330127018922193*F_0c[11]+0.25*F_0c[7]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[6]+0.25*F_0c[3]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[5]+0.25*F_0c[2]; 
  F_0c_rx[3] = 0.4330127018922193*F_0c[1]+0.25*F_0c[0]; 

  F_0r_lx[0] = 0.25*F_0r[7]-0.4330127018922193*F_0r[11]; 
  F_0r_lx[1] = 0.25*F_0r[3]-0.4330127018922193*F_0r[6]; 
  F_0r_lx[2] = 0.25*F_0r[2]-0.4330127018922193*F_0r[5]; 
  F_0r_lx[3] = 0.25*F_0r[0]-0.4330127018922193*F_0r[1]; 

  G_1l_rx[0] = 0.4330127018922193*G_1l[11]+0.25*G_1l[7]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[6]+0.25*G_1l[3]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[5]+0.25*G_1l[2]; 
  G_1l_rx[3] = 0.4330127018922193*G_1l[1]+0.25*G_1l[0]; 

  G_1c_lx[0] = 0.25*G_1c[7]-0.4330127018922193*G_1c[11]; 
  G_1c_lx[1] = 0.25*G_1c[3]-0.4330127018922193*G_1c[6]; 
  G_1c_lx[2] = 0.25*G_1c[2]-0.4330127018922193*G_1c[5]; 
  G_1c_lx[3] = 0.25*G_1c[0]-0.4330127018922193*G_1c[1]; 

  G_1c_rx[0] = 0.4330127018922193*G_1c[11]+0.25*G_1c[7]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[6]+0.25*G_1c[3]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[5]+0.25*G_1c[2]; 
  G_1c_rx[3] = 0.4330127018922193*G_1c[1]+0.25*G_1c[0]; 

  G_1r_lx[0] = 0.25*G_1r[7]-0.4330127018922193*G_1r[11]; 
  G_1r_lx[1] = 0.25*G_1r[3]-0.4330127018922193*G_1r[6]; 
  G_1r_lx[2] = 0.25*G_1r[2]-0.4330127018922193*G_1r[5]; 
  G_1r_lx[3] = 0.25*G_1r[0]-0.4330127018922193*G_1r[1]; 

  out_F_0[7] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[7] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[8])-0.4330127018922193*F_0l[4]; 
  F_0l_rx[1] = (-0.75*F_0l[12])-0.4330127018922193*F_0l[9]; 
  F_0l_rx[2] = (-0.75*F_0l[13])-0.4330127018922193*F_0l[10]; 
  F_0l_rx[3] = (-0.75*F_0l[15])-0.4330127018922193*F_0l[14]; 

  F_0c_lx[0] = 0.75*F_0c[8]-0.4330127018922193*F_0c[4]; 
  F_0c_lx[1] = 0.75*F_0c[12]-0.4330127018922193*F_0c[9]; 
  F_0c_lx[2] = 0.75*F_0c[13]-0.4330127018922193*F_0c[10]; 
  F_0c_lx[3] = 0.75*F_0c[15]-0.4330127018922193*F_0c[14]; 

  F_0c_rx[0] = 0.75*F_0c[8]+0.4330127018922193*F_0c[4]; 
  F_0c_rx[1] = 0.75*F_0c[12]+0.4330127018922193*F_0c[9]; 
  F_0c_rx[2] = 0.75*F_0c[13]+0.4330127018922193*F_0c[10]; 
  F_0c_rx[3] = 0.75*F_0c[15]+0.4330127018922193*F_0c[14]; 

  F_0r_lx[0] = 0.4330127018922193*F_0r[4]-0.75*F_0r[8]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[9]-0.75*F_0r[12]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[10]-0.75*F_0r[13]; 
  F_0r_lx[3] = 0.4330127018922193*F_0r[14]-0.75*F_0r[15]; 

  G_1l_rx[0] = (-0.75*G_1l[8])-0.4330127018922193*G_1l[4]; 
  G_1l_rx[1] = (-0.75*G_1l[12])-0.4330127018922193*G_1l[9]; 
  G_1l_rx[2] = (-0.75*G_1l[13])-0.4330127018922193*G_1l[10]; 
  G_1l_rx[3] = (-0.75*G_1l[15])-0.4330127018922193*G_1l[14]; 

  G_1c_lx[0] = 0.75*G_1c[8]-0.4330127018922193*G_1c[4]; 
  G_1c_lx[1] = 0.75*G_1c[12]-0.4330127018922193*G_1c[9]; 
  G_1c_lx[2] = 0.75*G_1c[13]-0.4330127018922193*G_1c[10]; 
  G_1c_lx[3] = 0.75*G_1c[15]-0.4330127018922193*G_1c[14]; 

  G_1c_rx[0] = 0.75*G_1c[8]+0.4330127018922193*G_1c[4]; 
  G_1c_rx[1] = 0.75*G_1c[12]+0.4330127018922193*G_1c[9]; 
  G_1c_rx[2] = 0.75*G_1c[13]+0.4330127018922193*G_1c[10]; 
  G_1c_rx[3] = 0.75*G_1c[15]+0.4330127018922193*G_1c[14]; 

  G_1r_lx[0] = 0.4330127018922193*G_1r[4]-0.75*G_1r[8]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[9]-0.75*G_1r[12]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[10]-0.75*G_1r[13]; 
  G_1r_lx[3] = 0.4330127018922193*G_1r[14]-0.75*G_1r[15]; 

  out_F_0[8] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[8] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922193*F_0l[12]+0.25*F_0l[9]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[8]+0.25*F_0l[4]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[15]+0.25*F_0l[14]; 
  F_0l_rx[3] = 0.4330127018922193*F_0l[13]+0.25*F_0l[10]; 

  F_0c_lx[0] = 0.25*F_0c[9]-0.4330127018922193*F_0c[12]; 
  F_0c_lx[1] = 0.25*F_0c[4]-0.4330127018922193*F_0c[8]; 
  F_0c_lx[2] = 0.25*F_0c[14]-0.4330127018922193*F_0c[15]; 
  F_0c_lx[3] = 0.25*F_0c[10]-0.4330127018922193*F_0c[13]; 

  F_0c_rx[0] = 0.4330127018922193*F_0c[12]+0.25*F_0c[9]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[8]+0.25*F_0c[4]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[15]+0.25*F_0c[14]; 
  F_0c_rx[3] = 0.4330127018922193*F_0c[13]+0.25*F_0c[10]; 

  F_0r_lx[0] = 0.25*F_0r[9]-0.4330127018922193*F_0r[12]; 
  F_0r_lx[1] = 0.25*F_0r[4]-0.4330127018922193*F_0r[8]; 
  F_0r_lx[2] = 0.25*F_0r[14]-0.4330127018922193*F_0r[15]; 
  F_0r_lx[3] = 0.25*F_0r[10]-0.4330127018922193*F_0r[13]; 

  G_1l_rx[0] = 0.4330127018922193*G_1l[12]+0.25*G_1l[9]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[8]+0.25*G_1l[4]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[15]+0.25*G_1l[14]; 
  G_1l_rx[3] = 0.4330127018922193*G_1l[13]+0.25*G_1l[10]; 

  G_1c_lx[0] = 0.25*G_1c[9]-0.4330127018922193*G_1c[12]; 
  G_1c_lx[1] = 0.25*G_1c[4]-0.4330127018922193*G_1c[8]; 
  G_1c_lx[2] = 0.25*G_1c[14]-0.4330127018922193*G_1c[15]; 
  G_1c_lx[3] = 0.25*G_1c[10]-0.4330127018922193*G_1c[13]; 

  G_1c_rx[0] = 0.4330127018922193*G_1c[12]+0.25*G_1c[9]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[8]+0.25*G_1c[4]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[15]+0.25*G_1c[14]; 
  G_1c_rx[3] = 0.4330127018922193*G_1c[13]+0.25*G_1c[10]; 

  G_1r_lx[0] = 0.25*G_1r[9]-0.4330127018922193*G_1r[12]; 
  G_1r_lx[1] = 0.25*G_1r[4]-0.4330127018922193*G_1r[8]; 
  G_1r_lx[2] = 0.25*G_1r[14]-0.4330127018922193*G_1r[15]; 
  G_1r_lx[3] = 0.25*G_1r[10]-0.4330127018922193*G_1r[13]; 

  out_F_0[9] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[9] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922193*F_0l[13]+0.25*F_0l[10]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[15]+0.25*F_0l[14]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[8]+0.25*F_0l[4]; 
  F_0l_rx[3] = 0.4330127018922193*F_0l[12]+0.25*F_0l[9]; 

  F_0c_lx[0] = 0.25*F_0c[10]-0.4330127018922193*F_0c[13]; 
  F_0c_lx[1] = 0.25*F_0c[14]-0.4330127018922193*F_0c[15]; 
  F_0c_lx[2] = 0.25*F_0c[4]-0.4330127018922193*F_0c[8]; 
  F_0c_lx[3] = 0.25*F_0c[9]-0.4330127018922193*F_0c[12]; 

  F_0c_rx[0] = 0.4330127018922193*F_0c[13]+0.25*F_0c[10]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[15]+0.25*F_0c[14]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[8]+0.25*F_0c[4]; 
  F_0c_rx[3] = 0.4330127018922193*F_0c[12]+0.25*F_0c[9]; 

  F_0r_lx[0] = 0.25*F_0r[10]-0.4330127018922193*F_0r[13]; 
  F_0r_lx[1] = 0.25*F_0r[14]-0.4330127018922193*F_0r[15]; 
  F_0r_lx[2] = 0.25*F_0r[4]-0.4330127018922193*F_0r[8]; 
  F_0r_lx[3] = 0.25*F_0r[9]-0.4330127018922193*F_0r[12]; 

  G_1l_rx[0] = 0.4330127018922193*G_1l[13]+0.25*G_1l[10]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[15]+0.25*G_1l[14]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[8]+0.25*G_1l[4]; 
  G_1l_rx[3] = 0.4330127018922193*G_1l[12]+0.25*G_1l[9]; 

  G_1c_lx[0] = 0.25*G_1c[10]-0.4330127018922193*G_1c[13]; 
  G_1c_lx[1] = 0.25*G_1c[14]-0.4330127018922193*G_1c[15]; 
  G_1c_lx[2] = 0.25*G_1c[4]-0.4330127018922193*G_1c[8]; 
  G_1c_lx[3] = 0.25*G_1c[9]-0.4330127018922193*G_1c[12]; 

  G_1c_rx[0] = 0.4330127018922193*G_1c[13]+0.25*G_1c[10]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[15]+0.25*G_1c[14]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[8]+0.25*G_1c[4]; 
  G_1c_rx[3] = 0.4330127018922193*G_1c[12]+0.25*G_1c[9]; 

  G_1r_lx[0] = 0.25*G_1r[10]-0.4330127018922193*G_1r[13]; 
  G_1r_lx[1] = 0.25*G_1r[14]-0.4330127018922193*G_1r[15]; 
  G_1r_lx[2] = 0.25*G_1r[4]-0.4330127018922193*G_1r[8]; 
  G_1r_lx[3] = 0.25*G_1r[9]-0.4330127018922193*G_1r[12]; 

  out_F_0[10] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[10] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[11])-0.4330127018922193*F_0l[7]; 
  F_0l_rx[1] = (-0.75*F_0l[6])-0.4330127018922193*F_0l[3]; 
  F_0l_rx[2] = (-0.75*F_0l[5])-0.4330127018922193*F_0l[2]; 
  F_0l_rx[3] = (-0.75*F_0l[1])-0.4330127018922193*F_0l[0]; 

  F_0c_lx[0] = 0.75*F_0c[11]-0.4330127018922193*F_0c[7]; 
  F_0c_lx[1] = 0.75*F_0c[6]-0.4330127018922193*F_0c[3]; 
  F_0c_lx[2] = 0.75*F_0c[5]-0.4330127018922193*F_0c[2]; 
  F_0c_lx[3] = 0.75*F_0c[1]-0.4330127018922193*F_0c[0]; 

  F_0c_rx[0] = 0.75*F_0c[11]+0.4330127018922193*F_0c[7]; 
  F_0c_rx[1] = 0.75*F_0c[6]+0.4330127018922193*F_0c[3]; 
  F_0c_rx[2] = 0.75*F_0c[5]+0.4330127018922193*F_0c[2]; 
  F_0c_rx[3] = 0.75*F_0c[1]+0.4330127018922193*F_0c[0]; 

  F_0r_lx[0] = 0.4330127018922193*F_0r[7]-0.75*F_0r[11]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[3]-0.75*F_0r[6]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[2]-0.75*F_0r[5]; 
  F_0r_lx[3] = 0.4330127018922193*F_0r[0]-0.75*F_0r[1]; 

  G_1l_rx[0] = (-0.75*G_1l[11])-0.4330127018922193*G_1l[7]; 
  G_1l_rx[1] = (-0.75*G_1l[6])-0.4330127018922193*G_1l[3]; 
  G_1l_rx[2] = (-0.75*G_1l[5])-0.4330127018922193*G_1l[2]; 
  G_1l_rx[3] = (-0.75*G_1l[1])-0.4330127018922193*G_1l[0]; 

  G_1c_lx[0] = 0.75*G_1c[11]-0.4330127018922193*G_1c[7]; 
  G_1c_lx[1] = 0.75*G_1c[6]-0.4330127018922193*G_1c[3]; 
  G_1c_lx[2] = 0.75*G_1c[5]-0.4330127018922193*G_1c[2]; 
  G_1c_lx[3] = 0.75*G_1c[1]-0.4330127018922193*G_1c[0]; 

  G_1c_rx[0] = 0.75*G_1c[11]+0.4330127018922193*G_1c[7]; 
  G_1c_rx[1] = 0.75*G_1c[6]+0.4330127018922193*G_1c[3]; 
  G_1c_rx[2] = 0.75*G_1c[5]+0.4330127018922193*G_1c[2]; 
  G_1c_rx[3] = 0.75*G_1c[1]+0.4330127018922193*G_1c[0]; 

  G_1r_lx[0] = 0.4330127018922193*G_1r[7]-0.75*G_1r[11]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[3]-0.75*G_1r[6]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[2]-0.75*G_1r[5]; 
  G_1r_lx[3] = 0.4330127018922193*G_1r[0]-0.75*G_1r[1]; 

  out_F_0[11] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[11] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[12])-0.4330127018922193*F_0l[9]; 
  F_0l_rx[1] = (-0.75*F_0l[8])-0.4330127018922193*F_0l[4]; 
  F_0l_rx[2] = (-0.75*F_0l[15])-0.4330127018922193*F_0l[14]; 
  F_0l_rx[3] = (-0.75*F_0l[13])-0.4330127018922193*F_0l[10]; 

  F_0c_lx[0] = 0.75*F_0c[12]-0.4330127018922193*F_0c[9]; 
  F_0c_lx[1] = 0.75*F_0c[8]-0.4330127018922193*F_0c[4]; 
  F_0c_lx[2] = 0.75*F_0c[15]-0.4330127018922193*F_0c[14]; 
  F_0c_lx[3] = 0.75*F_0c[13]-0.4330127018922193*F_0c[10]; 

  F_0c_rx[0] = 0.75*F_0c[12]+0.4330127018922193*F_0c[9]; 
  F_0c_rx[1] = 0.75*F_0c[8]+0.4330127018922193*F_0c[4]; 
  F_0c_rx[2] = 0.75*F_0c[15]+0.4330127018922193*F_0c[14]; 
  F_0c_rx[3] = 0.75*F_0c[13]+0.4330127018922193*F_0c[10]; 

  F_0r_lx[0] = 0.4330127018922193*F_0r[9]-0.75*F_0r[12]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[4]-0.75*F_0r[8]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[14]-0.75*F_0r[15]; 
  F_0r_lx[3] = 0.4330127018922193*F_0r[10]-0.75*F_0r[13]; 

  G_1l_rx[0] = (-0.75*G_1l[12])-0.4330127018922193*G_1l[9]; 
  G_1l_rx[1] = (-0.75*G_1l[8])-0.4330127018922193*G_1l[4]; 
  G_1l_rx[2] = (-0.75*G_1l[15])-0.4330127018922193*G_1l[14]; 
  G_1l_rx[3] = (-0.75*G_1l[13])-0.4330127018922193*G_1l[10]; 

  G_1c_lx[0] = 0.75*G_1c[12]-0.4330127018922193*G_1c[9]; 
  G_1c_lx[1] = 0.75*G_1c[8]-0.4330127018922193*G_1c[4]; 
  G_1c_lx[2] = 0.75*G_1c[15]-0.4330127018922193*G_1c[14]; 
  G_1c_lx[3] = 0.75*G_1c[13]-0.4330127018922193*G_1c[10]; 

  G_1c_rx[0] = 0.75*G_1c[12]+0.4330127018922193*G_1c[9]; 
  G_1c_rx[1] = 0.75*G_1c[8]+0.4330127018922193*G_1c[4]; 
  G_1c_rx[2] = 0.75*G_1c[15]+0.4330127018922193*G_1c[14]; 
  G_1c_rx[3] = 0.75*G_1c[13]+0.4330127018922193*G_1c[10]; 

  G_1r_lx[0] = 0.4330127018922193*G_1r[9]-0.75*G_1r[12]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[4]-0.75*G_1r[8]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[14]-0.75*G_1r[15]; 
  G_1r_lx[3] = 0.4330127018922193*G_1r[10]-0.75*G_1r[13]; 

  out_F_0[12] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[12] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[13])-0.4330127018922193*F_0l[10]; 
  F_0l_rx[1] = (-0.75*F_0l[15])-0.4330127018922193*F_0l[14]; 
  F_0l_rx[2] = (-0.75*F_0l[8])-0.4330127018922193*F_0l[4]; 
  F_0l_rx[3] = (-0.75*F_0l[12])-0.4330127018922193*F_0l[9]; 

  F_0c_lx[0] = 0.75*F_0c[13]-0.4330127018922193*F_0c[10]; 
  F_0c_lx[1] = 0.75*F_0c[15]-0.4330127018922193*F_0c[14]; 
  F_0c_lx[2] = 0.75*F_0c[8]-0.4330127018922193*F_0c[4]; 
  F_0c_lx[3] = 0.75*F_0c[12]-0.4330127018922193*F_0c[9]; 

  F_0c_rx[0] = 0.75*F_0c[13]+0.4330127018922193*F_0c[10]; 
  F_0c_rx[1] = 0.75*F_0c[15]+0.4330127018922193*F_0c[14]; 
  F_0c_rx[2] = 0.75*F_0c[8]+0.4330127018922193*F_0c[4]; 
  F_0c_rx[3] = 0.75*F_0c[12]+0.4330127018922193*F_0c[9]; 

  F_0r_lx[0] = 0.4330127018922193*F_0r[10]-0.75*F_0r[13]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[14]-0.75*F_0r[15]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[4]-0.75*F_0r[8]; 
  F_0r_lx[3] = 0.4330127018922193*F_0r[9]-0.75*F_0r[12]; 

  G_1l_rx[0] = (-0.75*G_1l[13])-0.4330127018922193*G_1l[10]; 
  G_1l_rx[1] = (-0.75*G_1l[15])-0.4330127018922193*G_1l[14]; 
  G_1l_rx[2] = (-0.75*G_1l[8])-0.4330127018922193*G_1l[4]; 
  G_1l_rx[3] = (-0.75*G_1l[12])-0.4330127018922193*G_1l[9]; 

  G_1c_lx[0] = 0.75*G_1c[13]-0.4330127018922193*G_1c[10]; 
  G_1c_lx[1] = 0.75*G_1c[15]-0.4330127018922193*G_1c[14]; 
  G_1c_lx[2] = 0.75*G_1c[8]-0.4330127018922193*G_1c[4]; 
  G_1c_lx[3] = 0.75*G_1c[12]-0.4330127018922193*G_1c[9]; 

  G_1c_rx[0] = 0.75*G_1c[13]+0.4330127018922193*G_1c[10]; 
  G_1c_rx[1] = 0.75*G_1c[15]+0.4330127018922193*G_1c[14]; 
  G_1c_rx[2] = 0.75*G_1c[8]+0.4330127018922193*G_1c[4]; 
  G_1c_rx[3] = 0.75*G_1c[12]+0.4330127018922193*G_1c[9]; 

  G_1r_lx[0] = 0.4330127018922193*G_1r[10]-0.75*G_1r[13]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[14]-0.75*G_1r[15]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[4]-0.75*G_1r[8]; 
  G_1r_lx[3] = 0.4330127018922193*G_1r[9]-0.75*G_1r[12]; 

  out_F_0[13] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[13] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922193*F_0l[15]+0.25*F_0l[14]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[13]+0.25*F_0l[10]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[12]+0.25*F_0l[9]; 
  F_0l_rx[3] = 0.4330127018922193*F_0l[8]+0.25*F_0l[4]; 

  F_0c_lx[0] = 0.25*F_0c[14]-0.4330127018922193*F_0c[15]; 
  F_0c_lx[1] = 0.25*F_0c[10]-0.4330127018922193*F_0c[13]; 
  F_0c_lx[2] = 0.25*F_0c[9]-0.4330127018922193*F_0c[12]; 
  F_0c_lx[3] = 0.25*F_0c[4]-0.4330127018922193*F_0c[8]; 

  F_0c_rx[0] = 0.4330127018922193*F_0c[15]+0.25*F_0c[14]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[13]+0.25*F_0c[10]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[12]+0.25*F_0c[9]; 
  F_0c_rx[3] = 0.4330127018922193*F_0c[8]+0.25*F_0c[4]; 

  F_0r_lx[0] = 0.25*F_0r[14]-0.4330127018922193*F_0r[15]; 
  F_0r_lx[1] = 0.25*F_0r[10]-0.4330127018922193*F_0r[13]; 
  F_0r_lx[2] = 0.25*F_0r[9]-0.4330127018922193*F_0r[12]; 
  F_0r_lx[3] = 0.25*F_0r[4]-0.4330127018922193*F_0r[8]; 

  G_1l_rx[0] = 0.4330127018922193*G_1l[15]+0.25*G_1l[14]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[13]+0.25*G_1l[10]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[12]+0.25*G_1l[9]; 
  G_1l_rx[3] = 0.4330127018922193*G_1l[8]+0.25*G_1l[4]; 

  G_1c_lx[0] = 0.25*G_1c[14]-0.4330127018922193*G_1c[15]; 
  G_1c_lx[1] = 0.25*G_1c[10]-0.4330127018922193*G_1c[13]; 
  G_1c_lx[2] = 0.25*G_1c[9]-0.4330127018922193*G_1c[12]; 
  G_1c_lx[3] = 0.25*G_1c[4]-0.4330127018922193*G_1c[8]; 

  G_1c_rx[0] = 0.4330127018922193*G_1c[15]+0.25*G_1c[14]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[13]+0.25*G_1c[10]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[12]+0.25*G_1c[9]; 
  G_1c_rx[3] = 0.4330127018922193*G_1c[8]+0.25*G_1c[4]; 

  G_1r_lx[0] = 0.25*G_1r[14]-0.4330127018922193*G_1r[15]; 
  G_1r_lx[1] = 0.25*G_1r[10]-0.4330127018922193*G_1r[13]; 
  G_1r_lx[2] = 0.25*G_1r[9]-0.4330127018922193*G_1r[12]; 
  G_1r_lx[3] = 0.25*G_1r[4]-0.4330127018922193*G_1r[8]; 

  out_F_0[14] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[14] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[15])-0.4330127018922193*F_0l[14]; 
  F_0l_rx[1] = (-0.75*F_0l[13])-0.4330127018922193*F_0l[10]; 
  F_0l_rx[2] = (-0.75*F_0l[12])-0.4330127018922193*F_0l[9]; 
  F_0l_rx[3] = (-0.75*F_0l[8])-0.4330127018922193*F_0l[4]; 

  F_0c_lx[0] = 0.75*F_0c[15]-0.4330127018922193*F_0c[14]; 
  F_0c_lx[1] = 0.75*F_0c[13]-0.4330127018922193*F_0c[10]; 
  F_0c_lx[2] = 0.75*F_0c[12]-0.4330127018922193*F_0c[9]; 
  F_0c_lx[3] = 0.75*F_0c[8]-0.4330127018922193*F_0c[4]; 

  F_0c_rx[0] = 0.75*F_0c[15]+0.4330127018922193*F_0c[14]; 
  F_0c_rx[1] = 0.75*F_0c[13]+0.4330127018922193*F_0c[10]; 
  F_0c_rx[2] = 0.75*F_0c[12]+0.4330127018922193*F_0c[9]; 
  F_0c_rx[3] = 0.75*F_0c[8]+0.4330127018922193*F_0c[4]; 

  F_0r_lx[0] = 0.4330127018922193*F_0r[14]-0.75*F_0r[15]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[10]-0.75*F_0r[13]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[9]-0.75*F_0r[12]; 
  F_0r_lx[3] = 0.4330127018922193*F_0r[4]-0.75*F_0r[8]; 

  G_1l_rx[0] = (-0.75*G_1l[15])-0.4330127018922193*G_1l[14]; 
  G_1l_rx[1] = (-0.75*G_1l[13])-0.4330127018922193*G_1l[10]; 
  G_1l_rx[2] = (-0.75*G_1l[12])-0.4330127018922193*G_1l[9]; 
  G_1l_rx[3] = (-0.75*G_1l[8])-0.4330127018922193*G_1l[4]; 

  G_1c_lx[0] = 0.75*G_1c[15]-0.4330127018922193*G_1c[14]; 
  G_1c_lx[1] = 0.75*G_1c[13]-0.4330127018922193*G_1c[10]; 
  G_1c_lx[2] = 0.75*G_1c[12]-0.4330127018922193*G_1c[9]; 
  G_1c_lx[3] = 0.75*G_1c[8]-0.4330127018922193*G_1c[4]; 

  G_1c_rx[0] = 0.75*G_1c[15]+0.4330127018922193*G_1c[14]; 
  G_1c_rx[1] = 0.75*G_1c[13]+0.4330127018922193*G_1c[10]; 
  G_1c_rx[2] = 0.75*G_1c[12]+0.4330127018922193*G_1c[9]; 
  G_1c_rx[3] = 0.75*G_1c[8]+0.4330127018922193*G_1c[4]; 

  G_1r_lx[0] = 0.4330127018922193*G_1r[14]-0.75*G_1r[15]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[10]-0.75*G_1r[13]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[9]-0.75*G_1r[12]; 
  G_1r_lx[3] = 0.4330127018922193*G_1r[4]-0.75*G_1r[8]; 

  out_F_0[15] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[15] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922194*F_0l[17]+0.25*F_0l[16]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[20]+0.2500000000000001*F_0l[18]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[21]+0.2500000000000001*F_0l[19]; 
  F_0l_rx[3] = 0.4330127018922194*F_0l[23]+0.25*F_0l[22]; 

  F_0c_lx[0] = 0.25*F_0c[16]-0.4330127018922194*F_0c[17]; 
  F_0c_lx[1] = 0.2500000000000001*F_0c[18]-0.4330127018922193*F_0c[20]; 
  F_0c_lx[2] = 0.2500000000000001*F_0c[19]-0.4330127018922193*F_0c[21]; 
  F_0c_lx[3] = 0.25*F_0c[22]-0.4330127018922194*F_0c[23]; 

  F_0c_rx[0] = 0.4330127018922194*F_0c[17]+0.25*F_0c[16]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[20]+0.2500000000000001*F_0c[18]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[21]+0.2500000000000001*F_0c[19]; 
  F_0c_rx[3] = 0.4330127018922194*F_0c[23]+0.25*F_0c[22]; 

  F_0r_lx[0] = 0.25*F_0r[16]-0.4330127018922194*F_0r[17]; 
  F_0r_lx[1] = 0.2500000000000001*F_0r[18]-0.4330127018922193*F_0r[20]; 
  F_0r_lx[2] = 0.2500000000000001*F_0r[19]-0.4330127018922193*F_0r[21]; 
  F_0r_lx[3] = 0.25*F_0r[22]-0.4330127018922194*F_0r[23]; 

  G_1l_rx[0] = 0.4330127018922194*G_1l[17]+0.25*G_1l[16]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[20]+0.2500000000000001*G_1l[18]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[21]+0.2500000000000001*G_1l[19]; 
  G_1l_rx[3] = 0.4330127018922194*G_1l[23]+0.25*G_1l[22]; 

  G_1c_lx[0] = 0.25*G_1c[16]-0.4330127018922194*G_1c[17]; 
  G_1c_lx[1] = 0.2500000000000001*G_1c[18]-0.4330127018922193*G_1c[20]; 
  G_1c_lx[2] = 0.2500000000000001*G_1c[19]-0.4330127018922193*G_1c[21]; 
  G_1c_lx[3] = 0.25*G_1c[22]-0.4330127018922194*G_1c[23]; 

  G_1c_rx[0] = 0.4330127018922194*G_1c[17]+0.25*G_1c[16]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[20]+0.2500000000000001*G_1c[18]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[21]+0.2500000000000001*G_1c[19]; 
  G_1c_rx[3] = 0.4330127018922194*G_1c[23]+0.25*G_1c[22]; 

  G_1r_lx[0] = 0.25*G_1r[16]-0.4330127018922194*G_1r[17]; 
  G_1r_lx[1] = 0.2500000000000001*G_1r[18]-0.4330127018922193*G_1r[20]; 
  G_1r_lx[2] = 0.2500000000000001*G_1r[19]-0.4330127018922193*G_1r[21]; 
  G_1r_lx[3] = 0.25*G_1r[22]-0.4330127018922194*G_1r[23]; 

  out_F_0[16] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[16] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[17])-0.4330127018922194*F_0l[16]; 
  F_0l_rx[1] = (-0.75*F_0l[20])-0.4330127018922193*F_0l[18]; 
  F_0l_rx[2] = (-0.75*F_0l[21])-0.4330127018922193*F_0l[19]; 
  F_0l_rx[3] = (-0.75*F_0l[23])-0.4330127018922194*F_0l[22]; 

  F_0c_lx[0] = 0.75*F_0c[17]-0.4330127018922194*F_0c[16]; 
  F_0c_lx[1] = 0.75*F_0c[20]-0.4330127018922193*F_0c[18]; 
  F_0c_lx[2] = 0.75*F_0c[21]-0.4330127018922193*F_0c[19]; 
  F_0c_lx[3] = 0.75*F_0c[23]-0.4330127018922194*F_0c[22]; 

  F_0c_rx[0] = 0.75*F_0c[17]+0.4330127018922194*F_0c[16]; 
  F_0c_rx[1] = 0.75*F_0c[20]+0.4330127018922193*F_0c[18]; 
  F_0c_rx[2] = 0.75*F_0c[21]+0.4330127018922193*F_0c[19]; 
  F_0c_rx[3] = 0.75*F_0c[23]+0.4330127018922194*F_0c[22]; 

  F_0r_lx[0] = 0.4330127018922194*F_0r[16]-0.75*F_0r[17]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[18]-0.75*F_0r[20]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[19]-0.75*F_0r[21]; 
  F_0r_lx[3] = 0.4330127018922194*F_0r[22]-0.75*F_0r[23]; 

  G_1l_rx[0] = (-0.75*G_1l[17])-0.4330127018922194*G_1l[16]; 
  G_1l_rx[1] = (-0.75*G_1l[20])-0.4330127018922193*G_1l[18]; 
  G_1l_rx[2] = (-0.75*G_1l[21])-0.4330127018922193*G_1l[19]; 
  G_1l_rx[3] = (-0.75*G_1l[23])-0.4330127018922194*G_1l[22]; 

  G_1c_lx[0] = 0.75*G_1c[17]-0.4330127018922194*G_1c[16]; 
  G_1c_lx[1] = 0.75*G_1c[20]-0.4330127018922193*G_1c[18]; 
  G_1c_lx[2] = 0.75*G_1c[21]-0.4330127018922193*G_1c[19]; 
  G_1c_lx[3] = 0.75*G_1c[23]-0.4330127018922194*G_1c[22]; 

  G_1c_rx[0] = 0.75*G_1c[17]+0.4330127018922194*G_1c[16]; 
  G_1c_rx[1] = 0.75*G_1c[20]+0.4330127018922193*G_1c[18]; 
  G_1c_rx[2] = 0.75*G_1c[21]+0.4330127018922193*G_1c[19]; 
  G_1c_rx[3] = 0.75*G_1c[23]+0.4330127018922194*G_1c[22]; 

  G_1r_lx[0] = 0.4330127018922194*G_1r[16]-0.75*G_1r[17]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[18]-0.75*G_1r[20]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[19]-0.75*G_1r[21]; 
  G_1r_lx[3] = 0.4330127018922194*G_1r[22]-0.75*G_1r[23]; 

  out_F_0[17] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[17] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922194*F_0l[20]+0.25*F_0l[18]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[17]+0.2500000000000001*F_0l[16]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[23]+0.2500000000000001*F_0l[22]; 
  F_0l_rx[3] = 0.4330127018922194*F_0l[21]+0.25*F_0l[19]; 

  F_0c_lx[0] = 0.25*F_0c[18]-0.4330127018922194*F_0c[20]; 
  F_0c_lx[1] = 0.2500000000000001*F_0c[16]-0.4330127018922193*F_0c[17]; 
  F_0c_lx[2] = 0.2500000000000001*F_0c[22]-0.4330127018922193*F_0c[23]; 
  F_0c_lx[3] = 0.25*F_0c[19]-0.4330127018922194*F_0c[21]; 

  F_0c_rx[0] = 0.4330127018922194*F_0c[20]+0.25*F_0c[18]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[17]+0.2500000000000001*F_0c[16]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[23]+0.2500000000000001*F_0c[22]; 
  F_0c_rx[3] = 0.4330127018922194*F_0c[21]+0.25*F_0c[19]; 

  F_0r_lx[0] = 0.25*F_0r[18]-0.4330127018922194*F_0r[20]; 
  F_0r_lx[1] = 0.2500000000000001*F_0r[16]-0.4330127018922193*F_0r[17]; 
  F_0r_lx[2] = 0.2500000000000001*F_0r[22]-0.4330127018922193*F_0r[23]; 
  F_0r_lx[3] = 0.25*F_0r[19]-0.4330127018922194*F_0r[21]; 

  G_1l_rx[0] = 0.4330127018922194*G_1l[20]+0.25*G_1l[18]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[17]+0.2500000000000001*G_1l[16]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[23]+0.2500000000000001*G_1l[22]; 
  G_1l_rx[3] = 0.4330127018922194*G_1l[21]+0.25*G_1l[19]; 

  G_1c_lx[0] = 0.25*G_1c[18]-0.4330127018922194*G_1c[20]; 
  G_1c_lx[1] = 0.2500000000000001*G_1c[16]-0.4330127018922193*G_1c[17]; 
  G_1c_lx[2] = 0.2500000000000001*G_1c[22]-0.4330127018922193*G_1c[23]; 
  G_1c_lx[3] = 0.25*G_1c[19]-0.4330127018922194*G_1c[21]; 

  G_1c_rx[0] = 0.4330127018922194*G_1c[20]+0.25*G_1c[18]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[17]+0.2500000000000001*G_1c[16]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[23]+0.2500000000000001*G_1c[22]; 
  G_1c_rx[3] = 0.4330127018922194*G_1c[21]+0.25*G_1c[19]; 

  G_1r_lx[0] = 0.25*G_1r[18]-0.4330127018922194*G_1r[20]; 
  G_1r_lx[1] = 0.2500000000000001*G_1r[16]-0.4330127018922193*G_1r[17]; 
  G_1r_lx[2] = 0.2500000000000001*G_1r[22]-0.4330127018922193*G_1r[23]; 
  G_1r_lx[3] = 0.25*G_1r[19]-0.4330127018922194*G_1r[21]; 

  out_F_0[18] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[18] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922194*F_0l[21]+0.25*F_0l[19]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[23]+0.2500000000000001*F_0l[22]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[17]+0.2500000000000001*F_0l[16]; 
  F_0l_rx[3] = 0.4330127018922194*F_0l[20]+0.25*F_0l[18]; 

  F_0c_lx[0] = 0.25*F_0c[19]-0.4330127018922194*F_0c[21]; 
  F_0c_lx[1] = 0.2500000000000001*F_0c[22]-0.4330127018922193*F_0c[23]; 
  F_0c_lx[2] = 0.2500000000000001*F_0c[16]-0.4330127018922193*F_0c[17]; 
  F_0c_lx[3] = 0.25*F_0c[18]-0.4330127018922194*F_0c[20]; 

  F_0c_rx[0] = 0.4330127018922194*F_0c[21]+0.25*F_0c[19]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[23]+0.2500000000000001*F_0c[22]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[17]+0.2500000000000001*F_0c[16]; 
  F_0c_rx[3] = 0.4330127018922194*F_0c[20]+0.25*F_0c[18]; 

  F_0r_lx[0] = 0.25*F_0r[19]-0.4330127018922194*F_0r[21]; 
  F_0r_lx[1] = 0.2500000000000001*F_0r[22]-0.4330127018922193*F_0r[23]; 
  F_0r_lx[2] = 0.2500000000000001*F_0r[16]-0.4330127018922193*F_0r[17]; 
  F_0r_lx[3] = 0.25*F_0r[18]-0.4330127018922194*F_0r[20]; 

  G_1l_rx[0] = 0.4330127018922194*G_1l[21]+0.25*G_1l[19]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[23]+0.2500000000000001*G_1l[22]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[17]+0.2500000000000001*G_1l[16]; 
  G_1l_rx[3] = 0.4330127018922194*G_1l[20]+0.25*G_1l[18]; 

  G_1c_lx[0] = 0.25*G_1c[19]-0.4330127018922194*G_1c[21]; 
  G_1c_lx[1] = 0.2500000000000001*G_1c[22]-0.4330127018922193*G_1c[23]; 
  G_1c_lx[2] = 0.2500000000000001*G_1c[16]-0.4330127018922193*G_1c[17]; 
  G_1c_lx[3] = 0.25*G_1c[18]-0.4330127018922194*G_1c[20]; 

  G_1c_rx[0] = 0.4330127018922194*G_1c[21]+0.25*G_1c[19]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[23]+0.2500000000000001*G_1c[22]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[17]+0.2500000000000001*G_1c[16]; 
  G_1c_rx[3] = 0.4330127018922194*G_1c[20]+0.25*G_1c[18]; 

  G_1r_lx[0] = 0.25*G_1r[19]-0.4330127018922194*G_1r[21]; 
  G_1r_lx[1] = 0.2500000000000001*G_1r[22]-0.4330127018922193*G_1r[23]; 
  G_1r_lx[2] = 0.2500000000000001*G_1r[16]-0.4330127018922193*G_1r[17]; 
  G_1r_lx[3] = 0.25*G_1r[18]-0.4330127018922194*G_1r[20]; 

  out_F_0[19] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[19] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[20])-0.4330127018922194*F_0l[18]; 
  F_0l_rx[1] = (-0.75*F_0l[17])-0.4330127018922193*F_0l[16]; 
  F_0l_rx[2] = (-0.75*F_0l[23])-0.4330127018922193*F_0l[22]; 
  F_0l_rx[3] = (-0.75*F_0l[21])-0.4330127018922194*F_0l[19]; 

  F_0c_lx[0] = 0.75*F_0c[20]-0.4330127018922194*F_0c[18]; 
  F_0c_lx[1] = 0.75*F_0c[17]-0.4330127018922193*F_0c[16]; 
  F_0c_lx[2] = 0.75*F_0c[23]-0.4330127018922193*F_0c[22]; 
  F_0c_lx[3] = 0.75*F_0c[21]-0.4330127018922194*F_0c[19]; 

  F_0c_rx[0] = 0.75*F_0c[20]+0.4330127018922194*F_0c[18]; 
  F_0c_rx[1] = 0.75*F_0c[17]+0.4330127018922193*F_0c[16]; 
  F_0c_rx[2] = 0.75*F_0c[23]+0.4330127018922193*F_0c[22]; 
  F_0c_rx[3] = 0.75*F_0c[21]+0.4330127018922194*F_0c[19]; 

  F_0r_lx[0] = 0.4330127018922194*F_0r[18]-0.75*F_0r[20]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[16]-0.75*F_0r[17]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[22]-0.75*F_0r[23]; 
  F_0r_lx[3] = 0.4330127018922194*F_0r[19]-0.75*F_0r[21]; 

  G_1l_rx[0] = (-0.75*G_1l[20])-0.4330127018922194*G_1l[18]; 
  G_1l_rx[1] = (-0.75*G_1l[17])-0.4330127018922193*G_1l[16]; 
  G_1l_rx[2] = (-0.75*G_1l[23])-0.4330127018922193*G_1l[22]; 
  G_1l_rx[3] = (-0.75*G_1l[21])-0.4330127018922194*G_1l[19]; 

  G_1c_lx[0] = 0.75*G_1c[20]-0.4330127018922194*G_1c[18]; 
  G_1c_lx[1] = 0.75*G_1c[17]-0.4330127018922193*G_1c[16]; 
  G_1c_lx[2] = 0.75*G_1c[23]-0.4330127018922193*G_1c[22]; 
  G_1c_lx[3] = 0.75*G_1c[21]-0.4330127018922194*G_1c[19]; 

  G_1c_rx[0] = 0.75*G_1c[20]+0.4330127018922194*G_1c[18]; 
  G_1c_rx[1] = 0.75*G_1c[17]+0.4330127018922193*G_1c[16]; 
  G_1c_rx[2] = 0.75*G_1c[23]+0.4330127018922193*G_1c[22]; 
  G_1c_rx[3] = 0.75*G_1c[21]+0.4330127018922194*G_1c[19]; 

  G_1r_lx[0] = 0.4330127018922194*G_1r[18]-0.75*G_1r[20]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[16]-0.75*G_1r[17]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[22]-0.75*G_1r[23]; 
  G_1r_lx[3] = 0.4330127018922194*G_1r[19]-0.75*G_1r[21]; 

  out_F_0[20] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[20] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[21])-0.4330127018922194*F_0l[19]; 
  F_0l_rx[1] = (-0.75*F_0l[23])-0.4330127018922193*F_0l[22]; 
  F_0l_rx[2] = (-0.75*F_0l[17])-0.4330127018922193*F_0l[16]; 
  F_0l_rx[3] = (-0.75*F_0l[20])-0.4330127018922194*F_0l[18]; 

  F_0c_lx[0] = 0.75*F_0c[21]-0.4330127018922194*F_0c[19]; 
  F_0c_lx[1] = 0.75*F_0c[23]-0.4330127018922193*F_0c[22]; 
  F_0c_lx[2] = 0.75*F_0c[17]-0.4330127018922193*F_0c[16]; 
  F_0c_lx[3] = 0.75*F_0c[20]-0.4330127018922194*F_0c[18]; 

  F_0c_rx[0] = 0.75*F_0c[21]+0.4330127018922194*F_0c[19]; 
  F_0c_rx[1] = 0.75*F_0c[23]+0.4330127018922193*F_0c[22]; 
  F_0c_rx[2] = 0.75*F_0c[17]+0.4330127018922193*F_0c[16]; 
  F_0c_rx[3] = 0.75*F_0c[20]+0.4330127018922194*F_0c[18]; 

  F_0r_lx[0] = 0.4330127018922194*F_0r[19]-0.75*F_0r[21]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[22]-0.75*F_0r[23]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[16]-0.75*F_0r[17]; 
  F_0r_lx[3] = 0.4330127018922194*F_0r[18]-0.75*F_0r[20]; 

  G_1l_rx[0] = (-0.75*G_1l[21])-0.4330127018922194*G_1l[19]; 
  G_1l_rx[1] = (-0.75*G_1l[23])-0.4330127018922193*G_1l[22]; 
  G_1l_rx[2] = (-0.75*G_1l[17])-0.4330127018922193*G_1l[16]; 
  G_1l_rx[3] = (-0.75*G_1l[20])-0.4330127018922194*G_1l[18]; 

  G_1c_lx[0] = 0.75*G_1c[21]-0.4330127018922194*G_1c[19]; 
  G_1c_lx[1] = 0.75*G_1c[23]-0.4330127018922193*G_1c[22]; 
  G_1c_lx[2] = 0.75*G_1c[17]-0.4330127018922193*G_1c[16]; 
  G_1c_lx[3] = 0.75*G_1c[20]-0.4330127018922194*G_1c[18]; 

  G_1c_rx[0] = 0.75*G_1c[21]+0.4330127018922194*G_1c[19]; 
  G_1c_rx[1] = 0.75*G_1c[23]+0.4330127018922193*G_1c[22]; 
  G_1c_rx[2] = 0.75*G_1c[17]+0.4330127018922193*G_1c[16]; 
  G_1c_rx[3] = 0.75*G_1c[20]+0.4330127018922194*G_1c[18]; 

  G_1r_lx[0] = 0.4330127018922194*G_1r[19]-0.75*G_1r[21]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[22]-0.75*G_1r[23]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[16]-0.75*G_1r[17]; 
  G_1r_lx[3] = 0.4330127018922194*G_1r[18]-0.75*G_1r[20]; 

  out_F_0[21] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[21] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = 0.4330127018922194*F_0l[23]+0.25*F_0l[22]; 
  F_0l_rx[1] = 0.4330127018922193*F_0l[21]+0.2500000000000001*F_0l[19]; 
  F_0l_rx[2] = 0.4330127018922193*F_0l[20]+0.2500000000000001*F_0l[18]; 
  F_0l_rx[3] = 0.4330127018922194*F_0l[17]+0.25*F_0l[16]; 

  F_0c_lx[0] = 0.25*F_0c[22]-0.4330127018922194*F_0c[23]; 
  F_0c_lx[1] = 0.2500000000000001*F_0c[19]-0.4330127018922193*F_0c[21]; 
  F_0c_lx[2] = 0.2500000000000001*F_0c[18]-0.4330127018922193*F_0c[20]; 
  F_0c_lx[3] = 0.25*F_0c[16]-0.4330127018922194*F_0c[17]; 

  F_0c_rx[0] = 0.4330127018922194*F_0c[23]+0.25*F_0c[22]; 
  F_0c_rx[1] = 0.4330127018922193*F_0c[21]+0.2500000000000001*F_0c[19]; 
  F_0c_rx[2] = 0.4330127018922193*F_0c[20]+0.2500000000000001*F_0c[18]; 
  F_0c_rx[3] = 0.4330127018922194*F_0c[17]+0.25*F_0c[16]; 

  F_0r_lx[0] = 0.25*F_0r[22]-0.4330127018922194*F_0r[23]; 
  F_0r_lx[1] = 0.2500000000000001*F_0r[19]-0.4330127018922193*F_0r[21]; 
  F_0r_lx[2] = 0.2500000000000001*F_0r[18]-0.4330127018922193*F_0r[20]; 
  F_0r_lx[3] = 0.25*F_0r[16]-0.4330127018922194*F_0r[17]; 

  G_1l_rx[0] = 0.4330127018922194*G_1l[23]+0.25*G_1l[22]; 
  G_1l_rx[1] = 0.4330127018922193*G_1l[21]+0.2500000000000001*G_1l[19]; 
  G_1l_rx[2] = 0.4330127018922193*G_1l[20]+0.2500000000000001*G_1l[18]; 
  G_1l_rx[3] = 0.4330127018922194*G_1l[17]+0.25*G_1l[16]; 

  G_1c_lx[0] = 0.25*G_1c[22]-0.4330127018922194*G_1c[23]; 
  G_1c_lx[1] = 0.2500000000000001*G_1c[19]-0.4330127018922193*G_1c[21]; 
  G_1c_lx[2] = 0.2500000000000001*G_1c[18]-0.4330127018922193*G_1c[20]; 
  G_1c_lx[3] = 0.25*G_1c[16]-0.4330127018922194*G_1c[17]; 

  G_1c_rx[0] = 0.4330127018922194*G_1c[23]+0.25*G_1c[22]; 
  G_1c_rx[1] = 0.4330127018922193*G_1c[21]+0.2500000000000001*G_1c[19]; 
  G_1c_rx[2] = 0.4330127018922193*G_1c[20]+0.2500000000000001*G_1c[18]; 
  G_1c_rx[3] = 0.4330127018922194*G_1c[17]+0.25*G_1c[16]; 

  G_1r_lx[0] = 0.25*G_1r[22]-0.4330127018922194*G_1r[23]; 
  G_1r_lx[1] = 0.2500000000000001*G_1r[19]-0.4330127018922193*G_1r[21]; 
  G_1r_lx[2] = 0.2500000000000001*G_1r[18]-0.4330127018922193*G_1r[20]; 
  G_1r_lx[3] = 0.25*G_1r[16]-0.4330127018922194*G_1r[17]; 

  out_F_0[22] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[22] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

  F_0l_rx[0] = (-0.75*F_0l[23])-0.4330127018922194*F_0l[22]; 
  F_0l_rx[1] = (-0.75*F_0l[21])-0.4330127018922193*F_0l[19]; 
  F_0l_rx[2] = (-0.75*F_0l[20])-0.4330127018922193*F_0l[18]; 
  F_0l_rx[3] = (-0.75*F_0l[17])-0.4330127018922194*F_0l[16]; 

  F_0c_lx[0] = 0.75*F_0c[23]-0.4330127018922194*F_0c[22]; 
  F_0c_lx[1] = 0.75*F_0c[21]-0.4330127018922193*F_0c[19]; 
  F_0c_lx[2] = 0.75*F_0c[20]-0.4330127018922193*F_0c[18]; 
  F_0c_lx[3] = 0.75*F_0c[17]-0.4330127018922194*F_0c[16]; 

  F_0c_rx[0] = 0.75*F_0c[23]+0.4330127018922194*F_0c[22]; 
  F_0c_rx[1] = 0.75*F_0c[21]+0.4330127018922193*F_0c[19]; 
  F_0c_rx[2] = 0.75*F_0c[20]+0.4330127018922193*F_0c[18]; 
  F_0c_rx[3] = 0.75*F_0c[17]+0.4330127018922194*F_0c[16]; 

  F_0r_lx[0] = 0.4330127018922194*F_0r[22]-0.75*F_0r[23]; 
  F_0r_lx[1] = 0.4330127018922193*F_0r[19]-0.75*F_0r[21]; 
  F_0r_lx[2] = 0.4330127018922193*F_0r[18]-0.75*F_0r[20]; 
  F_0r_lx[3] = 0.4330127018922194*F_0r[16]-0.75*F_0r[17]; 

  G_1l_rx[0] = (-0.75*G_1l[23])-0.4330127018922194*G_1l[22]; 
  G_1l_rx[1] = (-0.75*G_1l[21])-0.4330127018922193*G_1l[19]; 
  G_1l_rx[2] = (-0.75*G_1l[20])-0.4330127018922193*G_1l[18]; 
  G_1l_rx[3] = (-0.75*G_1l[17])-0.4330127018922194*G_1l[16]; 

  G_1c_lx[0] = 0.75*G_1c[23]-0.4330127018922194*G_1c[22]; 
  G_1c_lx[1] = 0.75*G_1c[21]-0.4330127018922193*G_1c[19]; 
  G_1c_lx[2] = 0.75*G_1c[20]-0.4330127018922193*G_1c[18]; 
  G_1c_lx[3] = 0.75*G_1c[17]-0.4330127018922194*G_1c[16]; 

  G_1c_rx[0] = 0.75*G_1c[23]+0.4330127018922194*G_1c[22]; 
  G_1c_rx[1] = 0.75*G_1c[21]+0.4330127018922193*G_1c[19]; 
  G_1c_rx[2] = 0.75*G_1c[20]+0.4330127018922193*G_1c[18]; 
  G_1c_rx[3] = 0.75*G_1c[17]+0.4330127018922194*G_1c[16]; 

  G_1r_lx[0] = 0.4330127018922194*G_1r[22]-0.75*G_1r[23]; 
  G_1r_lx[1] = 0.4330127018922193*G_1r[19]-0.75*G_1r[21]; 
  G_1r_lx[2] = 0.4330127018922193*G_1r[18]-0.75*G_1r[20]; 
  G_1r_lx[3] = 0.4330127018922194*G_1r[16]-0.75*G_1r[17]; 

  out_F_0[23] += 0.5*dx1*(F_0r_lx[3]*max_speed_modal_r[3]-1.0*F_0c_rx[3]*max_speed_modal_r[3]+F_0l_rx[3]*max_speed_modal_l[3]-1.0*F_0c_lx[3]*max_speed_modal_l[3]-1.0*F_0r_lx[3]*avg_u_r[3]-1.0*F_0c_rx[3]*avg_u_r[3]+F_0l_rx[3]*avg_u_l[3]+F_0c_lx[3]*avg_u_l[3]+F_0r_lx[2]*max_speed_modal_r[2]-1.0*F_0c_rx[2]*max_speed_modal_r[2]+F_0l_rx[2]*max_speed_modal_l[2]-1.0*F_0c_lx[2]*max_speed_modal_l[2]-1.0*F_0r_lx[2]*avg_u_r[2]-1.0*F_0c_rx[2]*avg_u_r[2]+F_0l_rx[2]*avg_u_l[2]+F_0c_lx[2]*avg_u_l[2]+F_0r_lx[1]*max_speed_modal_r[1]-1.0*F_0c_rx[1]*max_speed_modal_r[1]+F_0l_rx[1]*max_speed_modal_l[1]-1.0*F_0c_lx[1]*max_speed_modal_l[1]-1.0*F_0r_lx[1]*avg_u_r[1]-1.0*F_0c_rx[1]*avg_u_r[1]+F_0l_rx[1]*avg_u_l[1]+F_0c_lx[1]*avg_u_l[1]+F_0r_lx[0]*max_speed_modal_r[0]-1.0*F_0c_rx[0]*max_speed_modal_r[0]+F_0l_rx[0]*max_speed_modal_l[0]-1.0*F_0c_lx[0]*max_speed_modal_l[0]-1.0*F_0r_lx[0]*avg_u_r[0]-1.0*F_0c_rx[0]*avg_u_r[0]+F_0l_rx[0]*avg_u_l[0]+F_0c_lx[0]*avg_u_l[0]); 
  out_G_1[23] += 0.5*dx1*(G_1r_lx[3]*max_speed_modal_r[3]-1.0*G_1c_rx[3]*max_speed_modal_r[3]+G_1l_rx[3]*max_speed_modal_l[3]-1.0*G_1c_lx[3]*max_speed_modal_l[3]-1.0*G_1r_lx[3]*avg_u_r[3]-1.0*G_1c_rx[3]*avg_u_r[3]+G_1l_rx[3]*avg_u_l[3]+G_1c_lx[3]*avg_u_l[3]+G_1r_lx[2]*max_speed_modal_r[2]-1.0*G_1c_rx[2]*max_speed_modal_r[2]+G_1l_rx[2]*max_speed_modal_l[2]-1.0*G_1c_lx[2]*max_speed_modal_l[2]-1.0*G_1r_lx[2]*avg_u_r[2]-1.0*G_1c_rx[2]*avg_u_r[2]+G_1l_rx[2]*avg_u_l[2]+G_1c_lx[2]*avg_u_l[2]+G_1r_lx[1]*max_speed_modal_r[1]-1.0*G_1c_rx[1]*max_speed_modal_r[1]+G_1l_rx[1]*max_speed_modal_l[1]-1.0*G_1c_lx[1]*max_speed_modal_l[1]-1.0*G_1r_lx[1]*avg_u_r[1]-1.0*G_1c_rx[1]*avg_u_r[1]+G_1l_rx[1]*avg_u_l[1]+G_1c_lx[1]*avg_u_l[1]+G_1r_lx[0]*max_speed_modal_r[0]-1.0*G_1c_rx[0]*max_speed_modal_r[0]+G_1l_rx[0]*max_speed_modal_l[0]-1.0*G_1c_lx[0]*max_speed_modal_l[0]-1.0*G_1r_lx[0]*avg_u_r[0]-1.0*G_1c_rx[0]*avg_u_r[0]+G_1l_rx[0]*avg_u_l[0]+G_1c_lx[0]*avg_u_l[0]); 

} 
