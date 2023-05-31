#include <gkyl_vlasov_pkpm_kernels.h> 
#include <gkyl_basis_hyb_3x1v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_hyb_3x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void vlasov_pkpm_pressure_x_3x1v_ser_p1(const double *w, const double *dxv, 
     const double *bvarl, const double *bvarc, const double *bvarr, 
     const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:                 Cell-center coordinates.
  // dxv[NDIM]:               Cell spacing.
  // bvarl/bvarc/bvarr:       Input magnetic field unit vector in left/center/right cells.
  // fl/fc/fr:                Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells.
  // out:                     Incremented distribution function in center cell.
  const double dx1 = 2.0/dxv[0]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  const double volFact = dxv[3]/2.0; 
  const double *bl = &bvarl[0]; 
  const double *bc = &bvarc[0]; 
  const double *br = &bvarr[0]; 
  const double *F_0l = &fl[0]; 
  const double *F_0c = &fc[0]; 
  const double *F_0r = &fr[0]; 
  double *out_pressure = &out[0]; 
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
  if ((-0.4743416490252568*alphaSurf_l[7])+0.4743416490252568*(alphaSurf_l[6]+alphaSurf_l[5])+0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_l(F_0c); 
  } 
  if ((-0.4743416490252568*alphaSurf_r[7])+0.4743416490252568*(alphaSurf_r[6]+alphaSurf_r[5])+0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[0] = hyb_3x1v_p1_surfx1_eval_quad_node_0_l(F_0r); 
  } 
  if (0.3535533905932737*alphaSurf_l[4]-0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_l(F_0c); 
  } 
  if (0.3535533905932737*alphaSurf_r[4]-0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[1] = hyb_3x1v_p1_surfx1_eval_quad_node_1_l(F_0r); 
  } 
  if (0.4743416490252568*alphaSurf_l[7]-0.4743416490252568*(alphaSurf_l[6]+alphaSurf_l[5])+0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1])+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_l(F_0c); 
  } 
  if (0.4743416490252568*alphaSurf_r[7]-0.4743416490252568*(alphaSurf_r[6]+alphaSurf_r[5])+0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1])+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[2] = hyb_3x1v_p1_surfx1_eval_quad_node_2_l(F_0r); 
  } 
  if (0.4743416490252568*alphaSurf_l[7]-0.4743416490252568*alphaSurf_l[6]+0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[2]-0.3535533905932737*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_l(F_0c); 
  } 
  if (0.4743416490252568*alphaSurf_r[7]-0.4743416490252568*alphaSurf_r[6]+0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[2]-0.3535533905932737*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[3] = hyb_3x1v_p1_surfx1_eval_quad_node_3_l(F_0r); 
  } 
  if ((-0.3535533905932737*alphaSurf_l[4])+0.3535533905932737*alphaSurf_l[2]-0.3535533905932737*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_l(F_0c); 
  } 
  if ((-0.3535533905932737*alphaSurf_r[4])+0.3535533905932737*alphaSurf_r[2]-0.3535533905932737*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[4] = hyb_3x1v_p1_surfx1_eval_quad_node_4_l(F_0r); 
  } 
  if ((-0.4743416490252568*alphaSurf_l[7])+0.4743416490252568*alphaSurf_l[6]-0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*alphaSurf_l[2]-0.3535533905932737*alphaSurf_l[1]+0.3535533905932737*alphaSurf_l[0] > 0) { 
    F_0_UpwindQuad_l[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_l(F_0c); 
  } 
  if ((-0.4743416490252568*alphaSurf_r[7])+0.4743416490252568*alphaSurf_r[6]-0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*alphaSurf_r[2]-0.3535533905932737*alphaSurf_r[1]+0.3535533905932737*alphaSurf_r[0] > 0) { 
    F_0_UpwindQuad_r[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[5] = hyb_3x1v_p1_surfx1_eval_quad_node_5_l(F_0r); 
  } 
  if (0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6])-0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*alphaSurf_l[2]+0.3535533905932737*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_l(F_0c); 
  } 
  if (0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6])-0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*alphaSurf_r[2]+0.3535533905932737*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[6] = hyb_3x1v_p1_surfx1_eval_quad_node_6_l(F_0r); 
  } 
  if (0.3535533905932737*(alphaSurf_l[1]+alphaSurf_l[0])-0.3535533905932737*(alphaSurf_l[4]+alphaSurf_l[2]) > 0) { 
    F_0_UpwindQuad_l[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_l(F_0c); 
  } 
  if (0.3535533905932737*(alphaSurf_r[1]+alphaSurf_r[0])-0.3535533905932737*(alphaSurf_r[4]+alphaSurf_r[2]) > 0) { 
    F_0_UpwindQuad_r[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[7] = hyb_3x1v_p1_surfx1_eval_quad_node_7_l(F_0r); 
  } 
  if ((-0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6]))+0.4743416490252568*alphaSurf_l[5]-0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]-0.3535533905932737*alphaSurf_l[2]+0.3535533905932737*(alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_l(F_0c); 
  } 
  if ((-0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6]))+0.4743416490252568*alphaSurf_r[5]-0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]-0.3535533905932737*alphaSurf_r[2]+0.3535533905932737*(alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[8] = hyb_3x1v_p1_surfx1_eval_quad_node_8_l(F_0r); 
  } 
  if ((-0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6]+alphaSurf_l[5]))+0.3535533905932737*alphaSurf_l[4]-0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_l(F_0c); 
  } 
  if ((-0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6]+alphaSurf_r[5]))+0.3535533905932737*alphaSurf_r[4]-0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[9] = hyb_3x1v_p1_surfx1_eval_quad_node_9_l(F_0r); 
  } 
  if (0.3535533905932737*(alphaSurf_l[4]+alphaSurf_l[2]+alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_l(F_0c); 
  } 
  if (0.3535533905932737*(alphaSurf_r[4]+alphaSurf_r[2]+alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[10] = hyb_3x1v_p1_surfx1_eval_quad_node_10_l(F_0r); 
  } 
  if (0.4743416490252568*(alphaSurf_l[7]+alphaSurf_l[6]+alphaSurf_l[5])+0.3535533905932737*alphaSurf_l[4]+0.4743416490252568*alphaSurf_l[3]+0.3535533905932737*(alphaSurf_l[2]+alphaSurf_l[1]+alphaSurf_l[0]) > 0) { 
    F_0_UpwindQuad_l[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_r(F_0l); 
  } else { 
    F_0_UpwindQuad_l[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_l(F_0c); 
  } 
  if (0.4743416490252568*(alphaSurf_r[7]+alphaSurf_r[6]+alphaSurf_r[5])+0.3535533905932737*alphaSurf_r[4]+0.4743416490252568*alphaSurf_r[3]+0.3535533905932737*(alphaSurf_r[2]+alphaSurf_r[1]+alphaSurf_r[0]) > 0) { 
    F_0_UpwindQuad_r[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_r(F_0c); 
  } else { 
    F_0_UpwindQuad_r[11] = hyb_3x1v_p1_surfx1_eval_quad_node_11_l(F_0r); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_l, F_0_Upwind_l); 
  hyb_3x1v_p1_xdir_upwind_quad_to_modal(F_0_UpwindQuad_r, F_0_Upwind_r); 

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

  out_pressure[0] += dx1*volFact*((Ghat_F_0_r[0]-1.0*Ghat_F_0_l[0])*wvpar+(0.2886751345948129*Ghat_F_0_r[3]-0.2886751345948129*Ghat_F_0_l[3])*dvpar); 
  out_pressure[1] += dx1*volFact*((1.732050807568877*(Ghat_F_0_r[0]+Ghat_F_0_l[0])-0.6123724356957944*(F_0c[15]*alpha_c[15]+F_0c[14]*alpha_c[14]+F_0c[13]*alpha_c[13]+F_0c[12]*alpha_c[12]+F_0c[11]*alpha_c[11]+F_0c[10]*alpha_c[10]+F_0c[9]*alpha_c[9]+F_0c[8]*alpha_c[8]+F_0c[7]*alpha_c[7]+F_0c[6]*alpha_c[6]+F_0c[5]*alpha_c[5]+F_0c[4]*alpha_c[4]+F_0c[3]*alpha_c[3]+F_0c[2]*alpha_c[2]+F_0c[1]*alpha_c[1]+F_0c[0]*alpha_c[0]))*wvpar+((-0.1581138830084189*alpha_c[15]*F_0c[23])-0.1581138830084189*(alpha_c[14]*F_0c[22]+alpha_c[13]*F_0c[21]+alpha_c[12]*F_0c[20])-0.1581138830084189*(alpha_c[10]*F_0c[19]+alpha_c[9]*F_0c[18]+alpha_c[8]*F_0c[17])-0.1581138830084189*alpha_c[4]*F_0c[16]-0.1767766952966368*(F_0c[11]*alpha_c[15]+alpha_c[11]*F_0c[15]+F_0c[7]*alpha_c[14]+alpha_c[7]*F_0c[14]+F_0c[6]*alpha_c[13]+alpha_c[6]*F_0c[13]+F_0c[5]*alpha_c[12]+alpha_c[5]*F_0c[12]+F_0c[3]*alpha_c[10]+alpha_c[3]*F_0c[10]+F_0c[2]*alpha_c[9]+alpha_c[2]*F_0c[9]+F_0c[1]*alpha_c[8]+alpha_c[1]*F_0c[8]+F_0c[0]*alpha_c[4]+alpha_c[0]*F_0c[4])+0.5*(Ghat_F_0_r[3]+Ghat_F_0_l[3]))*dvpar); 
  out_pressure[2] += dx1*volFact*((Ghat_F_0_r[1]-1.0*Ghat_F_0_l[1])*wvpar+(0.2886751345948129*Ghat_F_0_r[5]-0.2886751345948129*Ghat_F_0_l[5])*dvpar); 
  out_pressure[3] += dx1*volFact*((Ghat_F_0_r[2]-1.0*Ghat_F_0_l[2])*wvpar+(0.2886751345948129*Ghat_F_0_r[6]-0.2886751345948129*Ghat_F_0_l[6])*dvpar); 
  out_pressure[4] += dx1*volFact*((1.732050807568877*(Ghat_F_0_r[1]+Ghat_F_0_l[1])-0.6123724356957944*(F_0c[13]*alpha_c[15]+alpha_c[13]*F_0c[15]+F_0c[10]*alpha_c[14]+alpha_c[10]*F_0c[14]+F_0c[8]*alpha_c[12]+alpha_c[8]*F_0c[12]+F_0c[6]*alpha_c[11]+alpha_c[6]*F_0c[11]+F_0c[4]*alpha_c[9]+alpha_c[4]*F_0c[9]+F_0c[3]*alpha_c[7]+alpha_c[3]*F_0c[7]+F_0c[1]*alpha_c[5]+alpha_c[1]*F_0c[5]+F_0c[0]*alpha_c[2]+alpha_c[0]*F_0c[2]))*wvpar+((-0.1581138830084189*alpha_c[13]*F_0c[23])-0.1581138830084189*(alpha_c[10]*F_0c[22]+alpha_c[15]*F_0c[21]+alpha_c[8]*F_0c[20])-0.1581138830084189*(alpha_c[14]*F_0c[19]+alpha_c[4]*F_0c[18]+alpha_c[12]*F_0c[17])-0.1581138830084189*alpha_c[9]*F_0c[16]-0.1767766952966368*(F_0c[6]*alpha_c[15]+alpha_c[6]*F_0c[15]+F_0c[3]*alpha_c[14]+alpha_c[3]*F_0c[14]+F_0c[11]*alpha_c[13]+alpha_c[11]*F_0c[13]+F_0c[1]*alpha_c[12]+alpha_c[1]*F_0c[12]+F_0c[7]*alpha_c[10]+alpha_c[7]*F_0c[10]+F_0c[0]*alpha_c[9]+alpha_c[0]*F_0c[9]+F_0c[5]*alpha_c[8]+alpha_c[5]*F_0c[8])+0.5*(Ghat_F_0_r[5]+Ghat_F_0_l[5])-0.1767766952966368*(F_0c[2]*alpha_c[4]+alpha_c[2]*F_0c[4]))*dvpar); 
  out_pressure[5] += dx1*volFact*((1.732050807568877*(Ghat_F_0_r[2]+Ghat_F_0_l[2])-0.6123724356957944*(F_0c[12]*alpha_c[15]+alpha_c[12]*F_0c[15]+F_0c[9]*alpha_c[14]+alpha_c[9]*F_0c[14]+F_0c[8]*alpha_c[13]+alpha_c[8]*F_0c[13]+F_0c[5]*alpha_c[11]+alpha_c[5]*F_0c[11]+F_0c[4]*alpha_c[10]+alpha_c[4]*F_0c[10]+F_0c[2]*alpha_c[7]+alpha_c[2]*F_0c[7]+F_0c[1]*alpha_c[6]+alpha_c[1]*F_0c[6]+F_0c[0]*alpha_c[3]+alpha_c[0]*F_0c[3]))*wvpar+((-0.1581138830084189*alpha_c[12]*F_0c[23])-0.1581138830084189*(alpha_c[9]*F_0c[22]+alpha_c[8]*F_0c[21]+alpha_c[15]*F_0c[20])-0.1581138830084189*(alpha_c[4]*F_0c[19]+alpha_c[14]*F_0c[18]+alpha_c[13]*F_0c[17])-0.1581138830084189*alpha_c[10]*F_0c[16]-0.1767766952966368*(F_0c[5]*alpha_c[15]+alpha_c[5]*F_0c[15]+F_0c[2]*alpha_c[14]+alpha_c[2]*F_0c[14]+F_0c[1]*alpha_c[13]+alpha_c[1]*F_0c[13]+F_0c[11]*alpha_c[12]+alpha_c[11]*F_0c[12]+F_0c[0]*alpha_c[10]+alpha_c[0]*F_0c[10]+F_0c[7]*alpha_c[9]+alpha_c[7]*F_0c[9]+F_0c[6]*alpha_c[8]+alpha_c[6]*F_0c[8])+0.5*(Ghat_F_0_r[6]+Ghat_F_0_l[6])-0.1767766952966368*(F_0c[3]*alpha_c[4]+alpha_c[3]*F_0c[4]))*dvpar); 
  out_pressure[6] += dx1*volFact*((Ghat_F_0_r[4]-1.0*Ghat_F_0_l[4])*wvpar+(0.2886751345948129*Ghat_F_0_r[7]-0.2886751345948129*Ghat_F_0_l[7])*dvpar); 
  out_pressure[7] += dx1*volFact*(((-0.6123724356957944*(F_0c[8]*alpha_c[15]+alpha_c[8]*F_0c[15]+F_0c[4]*alpha_c[14]+alpha_c[4]*F_0c[14]+F_0c[12]*alpha_c[13]+alpha_c[12]*F_0c[13]+F_0c[1]*alpha_c[11]+alpha_c[1]*F_0c[11]+F_0c[9]*alpha_c[10]+alpha_c[9]*F_0c[10]+F_0c[0]*alpha_c[7]+alpha_c[0]*F_0c[7]+F_0c[5]*alpha_c[6]+alpha_c[5]*F_0c[6]))+1.732050807568877*(Ghat_F_0_r[4]+Ghat_F_0_l[4])-0.6123724356957944*(F_0c[2]*alpha_c[3]+alpha_c[2]*F_0c[3]))*wvpar+((-0.1581138830084189*alpha_c[8]*F_0c[23])-0.1581138830084189*(alpha_c[4]*F_0c[22]+alpha_c[12]*F_0c[21]+alpha_c[13]*F_0c[20])-0.1581138830084189*(alpha_c[9]*F_0c[19]+alpha_c[10]*F_0c[18]+alpha_c[15]*F_0c[17])-0.1581138830084189*alpha_c[14]*F_0c[16]-0.1767766952966368*(F_0c[1]*alpha_c[15]+alpha_c[1]*F_0c[15]+F_0c[0]*alpha_c[14]+alpha_c[0]*F_0c[14]+F_0c[5]*alpha_c[13]+alpha_c[5]*F_0c[13]+F_0c[6]*alpha_c[12]+alpha_c[6]*F_0c[12]+F_0c[8]*alpha_c[11]+alpha_c[8]*F_0c[11]+F_0c[2]*alpha_c[10]+alpha_c[2]*F_0c[10]+F_0c[3]*alpha_c[9]+alpha_c[3]*F_0c[9]+F_0c[4]*alpha_c[7])+0.5*(Ghat_F_0_r[7]+Ghat_F_0_l[7])-0.1767766952966368*alpha_c[4]*F_0c[7])*dvpar); 

} 
