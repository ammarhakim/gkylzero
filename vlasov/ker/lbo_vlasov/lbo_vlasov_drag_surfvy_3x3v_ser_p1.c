#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_hyb_3x3v_p1_surfx5_eval_quad.h> 
#include <gkyl_basis_hyb_3x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[6]: cell-center coordinates. 
  // dxv[6]: cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[32]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fl/fc/fr: distribution function in cells 
  // out: incremented distribution function in cell 
  double rdv2 = 2.0/dxv[4]; 

  const double *sumNuUy = &nuPrimMomsSum[8]; 

  double alphaDrSurf_l[64] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*w[4]-1.0*nuSum[0]*dxv[4]-2.0*sumNuUy[0]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*w[4]-1.0*nuSum[1]*dxv[4]-2.0*sumNuUy[1]; 
  alphaDrSurf_l[2] = 2.0*nuSum[2]*w[4]-1.0*nuSum[2]*dxv[4]-2.0*sumNuUy[2]; 
  alphaDrSurf_l[3] = 2.0*nuSum[3]*w[4]-1.0*nuSum[3]*dxv[4]-2.0*sumNuUy[3]; 
  alphaDrSurf_l[6] = 2.0*nuSum[4]*w[4]-2.0*sumNuUy[4]-1.0*dxv[4]*nuSum[4]; 
  alphaDrSurf_l[7] = (-2.0*sumNuUy[5])+2.0*w[4]*nuSum[5]-1.0*dxv[4]*nuSum[5]; 
  alphaDrSurf_l[8] = (-2.0*sumNuUy[6])+2.0*w[4]*nuSum[6]-1.0*dxv[4]*nuSum[6]; 
  alphaDrSurf_l[16] = (-2.0*sumNuUy[7])+2.0*w[4]*nuSum[7]-1.0*dxv[4]*nuSum[7]; 

  double alphaDrSurf_r[64] = {0.0}; 
  alphaDrSurf_r[0] = 2.0*nuSum[0]*w[4]+nuSum[0]*dxv[4]-2.0*sumNuUy[0]; 
  alphaDrSurf_r[1] = 2.0*nuSum[1]*w[4]+nuSum[1]*dxv[4]-2.0*sumNuUy[1]; 
  alphaDrSurf_r[2] = 2.0*nuSum[2]*w[4]+nuSum[2]*dxv[4]-2.0*sumNuUy[2]; 
  alphaDrSurf_r[3] = 2.0*nuSum[3]*w[4]+nuSum[3]*dxv[4]-2.0*sumNuUy[3]; 
  alphaDrSurf_r[6] = 2.0*nuSum[4]*w[4]-2.0*sumNuUy[4]+dxv[4]*nuSum[4]; 
  alphaDrSurf_r[7] = (-2.0*sumNuUy[5])+2.0*w[4]*nuSum[5]+dxv[4]*nuSum[5]; 
  alphaDrSurf_r[8] = (-2.0*sumNuUy[6])+2.0*w[4]*nuSum[6]+dxv[4]*nuSum[6]; 
  alphaDrSurf_r[16] = (-2.0*sumNuUy[7])+2.0*w[4]*nuSum[7]+dxv[4]*nuSum[7]; 

  double fUpwindQuad_l[72] = {0.0};
  double fUpwindQuad_r[72] = {0.0};
  double fUpwind_l[64] = {0.0};
  double fUpwind_r[64] = {0.0};
  double Ghat_l[64] = {0.0}; 
  double Ghat_r[64] = {0.0}; 

  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = hyb_3x3v_p1_surfx5_eval_quad_node_0_r(fl); 
    fUpwindQuad_l[1] = hyb_3x3v_p1_surfx5_eval_quad_node_1_r(fl); 
    fUpwindQuad_l[2] = hyb_3x3v_p1_surfx5_eval_quad_node_2_r(fl); 
    fUpwindQuad_l[3] = hyb_3x3v_p1_surfx5_eval_quad_node_3_r(fl); 
    fUpwindQuad_l[4] = hyb_3x3v_p1_surfx5_eval_quad_node_4_r(fl); 
    fUpwindQuad_l[5] = hyb_3x3v_p1_surfx5_eval_quad_node_5_r(fl); 
    fUpwindQuad_l[6] = hyb_3x3v_p1_surfx5_eval_quad_node_6_r(fl); 
    fUpwindQuad_l[7] = hyb_3x3v_p1_surfx5_eval_quad_node_7_r(fl); 
    fUpwindQuad_l[8] = hyb_3x3v_p1_surfx5_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[0] = hyb_3x3v_p1_surfx5_eval_quad_node_0_l(fc); 
    fUpwindQuad_l[1] = hyb_3x3v_p1_surfx5_eval_quad_node_1_l(fc); 
    fUpwindQuad_l[2] = hyb_3x3v_p1_surfx5_eval_quad_node_2_l(fc); 
    fUpwindQuad_l[3] = hyb_3x3v_p1_surfx5_eval_quad_node_3_l(fc); 
    fUpwindQuad_l[4] = hyb_3x3v_p1_surfx5_eval_quad_node_4_l(fc); 
    fUpwindQuad_l[5] = hyb_3x3v_p1_surfx5_eval_quad_node_5_l(fc); 
    fUpwindQuad_l[6] = hyb_3x3v_p1_surfx5_eval_quad_node_6_l(fc); 
    fUpwindQuad_l[7] = hyb_3x3v_p1_surfx5_eval_quad_node_7_l(fc); 
    fUpwindQuad_l[8] = hyb_3x3v_p1_surfx5_eval_quad_node_8_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = hyb_3x3v_p1_surfx5_eval_quad_node_0_r(fc); 
    fUpwindQuad_r[1] = hyb_3x3v_p1_surfx5_eval_quad_node_1_r(fc); 
    fUpwindQuad_r[2] = hyb_3x3v_p1_surfx5_eval_quad_node_2_r(fc); 
    fUpwindQuad_r[3] = hyb_3x3v_p1_surfx5_eval_quad_node_3_r(fc); 
    fUpwindQuad_r[4] = hyb_3x3v_p1_surfx5_eval_quad_node_4_r(fc); 
    fUpwindQuad_r[5] = hyb_3x3v_p1_surfx5_eval_quad_node_5_r(fc); 
    fUpwindQuad_r[6] = hyb_3x3v_p1_surfx5_eval_quad_node_6_r(fc); 
    fUpwindQuad_r[7] = hyb_3x3v_p1_surfx5_eval_quad_node_7_r(fc); 
    fUpwindQuad_r[8] = hyb_3x3v_p1_surfx5_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[0] = hyb_3x3v_p1_surfx5_eval_quad_node_0_l(fr); 
    fUpwindQuad_r[1] = hyb_3x3v_p1_surfx5_eval_quad_node_1_l(fr); 
    fUpwindQuad_r[2] = hyb_3x3v_p1_surfx5_eval_quad_node_2_l(fr); 
    fUpwindQuad_r[3] = hyb_3x3v_p1_surfx5_eval_quad_node_3_l(fr); 
    fUpwindQuad_r[4] = hyb_3x3v_p1_surfx5_eval_quad_node_4_l(fr); 
    fUpwindQuad_r[5] = hyb_3x3v_p1_surfx5_eval_quad_node_5_l(fr); 
    fUpwindQuad_r[6] = hyb_3x3v_p1_surfx5_eval_quad_node_6_l(fr); 
    fUpwindQuad_r[7] = hyb_3x3v_p1_surfx5_eval_quad_node_7_l(fr); 
    fUpwindQuad_r[8] = hyb_3x3v_p1_surfx5_eval_quad_node_8_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[9] = hyb_3x3v_p1_surfx5_eval_quad_node_9_r(fl); 
    fUpwindQuad_l[10] = hyb_3x3v_p1_surfx5_eval_quad_node_10_r(fl); 
    fUpwindQuad_l[11] = hyb_3x3v_p1_surfx5_eval_quad_node_11_r(fl); 
    fUpwindQuad_l[12] = hyb_3x3v_p1_surfx5_eval_quad_node_12_r(fl); 
    fUpwindQuad_l[13] = hyb_3x3v_p1_surfx5_eval_quad_node_13_r(fl); 
    fUpwindQuad_l[14] = hyb_3x3v_p1_surfx5_eval_quad_node_14_r(fl); 
    fUpwindQuad_l[15] = hyb_3x3v_p1_surfx5_eval_quad_node_15_r(fl); 
    fUpwindQuad_l[16] = hyb_3x3v_p1_surfx5_eval_quad_node_16_r(fl); 
    fUpwindQuad_l[17] = hyb_3x3v_p1_surfx5_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_l[9] = hyb_3x3v_p1_surfx5_eval_quad_node_9_l(fc); 
    fUpwindQuad_l[10] = hyb_3x3v_p1_surfx5_eval_quad_node_10_l(fc); 
    fUpwindQuad_l[11] = hyb_3x3v_p1_surfx5_eval_quad_node_11_l(fc); 
    fUpwindQuad_l[12] = hyb_3x3v_p1_surfx5_eval_quad_node_12_l(fc); 
    fUpwindQuad_l[13] = hyb_3x3v_p1_surfx5_eval_quad_node_13_l(fc); 
    fUpwindQuad_l[14] = hyb_3x3v_p1_surfx5_eval_quad_node_14_l(fc); 
    fUpwindQuad_l[15] = hyb_3x3v_p1_surfx5_eval_quad_node_15_l(fc); 
    fUpwindQuad_l[16] = hyb_3x3v_p1_surfx5_eval_quad_node_16_l(fc); 
    fUpwindQuad_l[17] = hyb_3x3v_p1_surfx5_eval_quad_node_17_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[9] = hyb_3x3v_p1_surfx5_eval_quad_node_9_r(fc); 
    fUpwindQuad_r[10] = hyb_3x3v_p1_surfx5_eval_quad_node_10_r(fc); 
    fUpwindQuad_r[11] = hyb_3x3v_p1_surfx5_eval_quad_node_11_r(fc); 
    fUpwindQuad_r[12] = hyb_3x3v_p1_surfx5_eval_quad_node_12_r(fc); 
    fUpwindQuad_r[13] = hyb_3x3v_p1_surfx5_eval_quad_node_13_r(fc); 
    fUpwindQuad_r[14] = hyb_3x3v_p1_surfx5_eval_quad_node_14_r(fc); 
    fUpwindQuad_r[15] = hyb_3x3v_p1_surfx5_eval_quad_node_15_r(fc); 
    fUpwindQuad_r[16] = hyb_3x3v_p1_surfx5_eval_quad_node_16_r(fc); 
    fUpwindQuad_r[17] = hyb_3x3v_p1_surfx5_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_r[9] = hyb_3x3v_p1_surfx5_eval_quad_node_9_l(fr); 
    fUpwindQuad_r[10] = hyb_3x3v_p1_surfx5_eval_quad_node_10_l(fr); 
    fUpwindQuad_r[11] = hyb_3x3v_p1_surfx5_eval_quad_node_11_l(fr); 
    fUpwindQuad_r[12] = hyb_3x3v_p1_surfx5_eval_quad_node_12_l(fr); 
    fUpwindQuad_r[13] = hyb_3x3v_p1_surfx5_eval_quad_node_13_l(fr); 
    fUpwindQuad_r[14] = hyb_3x3v_p1_surfx5_eval_quad_node_14_l(fr); 
    fUpwindQuad_r[15] = hyb_3x3v_p1_surfx5_eval_quad_node_15_l(fr); 
    fUpwindQuad_r[16] = hyb_3x3v_p1_surfx5_eval_quad_node_16_l(fr); 
    fUpwindQuad_r[17] = hyb_3x3v_p1_surfx5_eval_quad_node_17_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[18] = hyb_3x3v_p1_surfx5_eval_quad_node_18_r(fl); 
    fUpwindQuad_l[19] = hyb_3x3v_p1_surfx5_eval_quad_node_19_r(fl); 
    fUpwindQuad_l[20] = hyb_3x3v_p1_surfx5_eval_quad_node_20_r(fl); 
    fUpwindQuad_l[21] = hyb_3x3v_p1_surfx5_eval_quad_node_21_r(fl); 
    fUpwindQuad_l[22] = hyb_3x3v_p1_surfx5_eval_quad_node_22_r(fl); 
    fUpwindQuad_l[23] = hyb_3x3v_p1_surfx5_eval_quad_node_23_r(fl); 
    fUpwindQuad_l[24] = hyb_3x3v_p1_surfx5_eval_quad_node_24_r(fl); 
    fUpwindQuad_l[25] = hyb_3x3v_p1_surfx5_eval_quad_node_25_r(fl); 
    fUpwindQuad_l[26] = hyb_3x3v_p1_surfx5_eval_quad_node_26_r(fl); 
  } else { 
    fUpwindQuad_l[18] = hyb_3x3v_p1_surfx5_eval_quad_node_18_l(fc); 
    fUpwindQuad_l[19] = hyb_3x3v_p1_surfx5_eval_quad_node_19_l(fc); 
    fUpwindQuad_l[20] = hyb_3x3v_p1_surfx5_eval_quad_node_20_l(fc); 
    fUpwindQuad_l[21] = hyb_3x3v_p1_surfx5_eval_quad_node_21_l(fc); 
    fUpwindQuad_l[22] = hyb_3x3v_p1_surfx5_eval_quad_node_22_l(fc); 
    fUpwindQuad_l[23] = hyb_3x3v_p1_surfx5_eval_quad_node_23_l(fc); 
    fUpwindQuad_l[24] = hyb_3x3v_p1_surfx5_eval_quad_node_24_l(fc); 
    fUpwindQuad_l[25] = hyb_3x3v_p1_surfx5_eval_quad_node_25_l(fc); 
    fUpwindQuad_l[26] = hyb_3x3v_p1_surfx5_eval_quad_node_26_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[18] = hyb_3x3v_p1_surfx5_eval_quad_node_18_r(fc); 
    fUpwindQuad_r[19] = hyb_3x3v_p1_surfx5_eval_quad_node_19_r(fc); 
    fUpwindQuad_r[20] = hyb_3x3v_p1_surfx5_eval_quad_node_20_r(fc); 
    fUpwindQuad_r[21] = hyb_3x3v_p1_surfx5_eval_quad_node_21_r(fc); 
    fUpwindQuad_r[22] = hyb_3x3v_p1_surfx5_eval_quad_node_22_r(fc); 
    fUpwindQuad_r[23] = hyb_3x3v_p1_surfx5_eval_quad_node_23_r(fc); 
    fUpwindQuad_r[24] = hyb_3x3v_p1_surfx5_eval_quad_node_24_r(fc); 
    fUpwindQuad_r[25] = hyb_3x3v_p1_surfx5_eval_quad_node_25_r(fc); 
    fUpwindQuad_r[26] = hyb_3x3v_p1_surfx5_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_r[18] = hyb_3x3v_p1_surfx5_eval_quad_node_18_l(fr); 
    fUpwindQuad_r[19] = hyb_3x3v_p1_surfx5_eval_quad_node_19_l(fr); 
    fUpwindQuad_r[20] = hyb_3x3v_p1_surfx5_eval_quad_node_20_l(fr); 
    fUpwindQuad_r[21] = hyb_3x3v_p1_surfx5_eval_quad_node_21_l(fr); 
    fUpwindQuad_r[22] = hyb_3x3v_p1_surfx5_eval_quad_node_22_l(fr); 
    fUpwindQuad_r[23] = hyb_3x3v_p1_surfx5_eval_quad_node_23_l(fr); 
    fUpwindQuad_r[24] = hyb_3x3v_p1_surfx5_eval_quad_node_24_l(fr); 
    fUpwindQuad_r[25] = hyb_3x3v_p1_surfx5_eval_quad_node_25_l(fr); 
    fUpwindQuad_r[26] = hyb_3x3v_p1_surfx5_eval_quad_node_26_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[27] = hyb_3x3v_p1_surfx5_eval_quad_node_27_r(fl); 
    fUpwindQuad_l[28] = hyb_3x3v_p1_surfx5_eval_quad_node_28_r(fl); 
    fUpwindQuad_l[29] = hyb_3x3v_p1_surfx5_eval_quad_node_29_r(fl); 
    fUpwindQuad_l[30] = hyb_3x3v_p1_surfx5_eval_quad_node_30_r(fl); 
    fUpwindQuad_l[31] = hyb_3x3v_p1_surfx5_eval_quad_node_31_r(fl); 
    fUpwindQuad_l[32] = hyb_3x3v_p1_surfx5_eval_quad_node_32_r(fl); 
    fUpwindQuad_l[33] = hyb_3x3v_p1_surfx5_eval_quad_node_33_r(fl); 
    fUpwindQuad_l[34] = hyb_3x3v_p1_surfx5_eval_quad_node_34_r(fl); 
    fUpwindQuad_l[35] = hyb_3x3v_p1_surfx5_eval_quad_node_35_r(fl); 
  } else { 
    fUpwindQuad_l[27] = hyb_3x3v_p1_surfx5_eval_quad_node_27_l(fc); 
    fUpwindQuad_l[28] = hyb_3x3v_p1_surfx5_eval_quad_node_28_l(fc); 
    fUpwindQuad_l[29] = hyb_3x3v_p1_surfx5_eval_quad_node_29_l(fc); 
    fUpwindQuad_l[30] = hyb_3x3v_p1_surfx5_eval_quad_node_30_l(fc); 
    fUpwindQuad_l[31] = hyb_3x3v_p1_surfx5_eval_quad_node_31_l(fc); 
    fUpwindQuad_l[32] = hyb_3x3v_p1_surfx5_eval_quad_node_32_l(fc); 
    fUpwindQuad_l[33] = hyb_3x3v_p1_surfx5_eval_quad_node_33_l(fc); 
    fUpwindQuad_l[34] = hyb_3x3v_p1_surfx5_eval_quad_node_34_l(fc); 
    fUpwindQuad_l[35] = hyb_3x3v_p1_surfx5_eval_quad_node_35_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[27] = hyb_3x3v_p1_surfx5_eval_quad_node_27_r(fc); 
    fUpwindQuad_r[28] = hyb_3x3v_p1_surfx5_eval_quad_node_28_r(fc); 
    fUpwindQuad_r[29] = hyb_3x3v_p1_surfx5_eval_quad_node_29_r(fc); 
    fUpwindQuad_r[30] = hyb_3x3v_p1_surfx5_eval_quad_node_30_r(fc); 
    fUpwindQuad_r[31] = hyb_3x3v_p1_surfx5_eval_quad_node_31_r(fc); 
    fUpwindQuad_r[32] = hyb_3x3v_p1_surfx5_eval_quad_node_32_r(fc); 
    fUpwindQuad_r[33] = hyb_3x3v_p1_surfx5_eval_quad_node_33_r(fc); 
    fUpwindQuad_r[34] = hyb_3x3v_p1_surfx5_eval_quad_node_34_r(fc); 
    fUpwindQuad_r[35] = hyb_3x3v_p1_surfx5_eval_quad_node_35_r(fc); 
  } else { 
    fUpwindQuad_r[27] = hyb_3x3v_p1_surfx5_eval_quad_node_27_l(fr); 
    fUpwindQuad_r[28] = hyb_3x3v_p1_surfx5_eval_quad_node_28_l(fr); 
    fUpwindQuad_r[29] = hyb_3x3v_p1_surfx5_eval_quad_node_29_l(fr); 
    fUpwindQuad_r[30] = hyb_3x3v_p1_surfx5_eval_quad_node_30_l(fr); 
    fUpwindQuad_r[31] = hyb_3x3v_p1_surfx5_eval_quad_node_31_l(fr); 
    fUpwindQuad_r[32] = hyb_3x3v_p1_surfx5_eval_quad_node_32_l(fr); 
    fUpwindQuad_r[33] = hyb_3x3v_p1_surfx5_eval_quad_node_33_l(fr); 
    fUpwindQuad_r[34] = hyb_3x3v_p1_surfx5_eval_quad_node_34_l(fr); 
    fUpwindQuad_r[35] = hyb_3x3v_p1_surfx5_eval_quad_node_35_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[36] = hyb_3x3v_p1_surfx5_eval_quad_node_36_r(fl); 
    fUpwindQuad_l[37] = hyb_3x3v_p1_surfx5_eval_quad_node_37_r(fl); 
    fUpwindQuad_l[38] = hyb_3x3v_p1_surfx5_eval_quad_node_38_r(fl); 
    fUpwindQuad_l[39] = hyb_3x3v_p1_surfx5_eval_quad_node_39_r(fl); 
    fUpwindQuad_l[40] = hyb_3x3v_p1_surfx5_eval_quad_node_40_r(fl); 
    fUpwindQuad_l[41] = hyb_3x3v_p1_surfx5_eval_quad_node_41_r(fl); 
    fUpwindQuad_l[42] = hyb_3x3v_p1_surfx5_eval_quad_node_42_r(fl); 
    fUpwindQuad_l[43] = hyb_3x3v_p1_surfx5_eval_quad_node_43_r(fl); 
    fUpwindQuad_l[44] = hyb_3x3v_p1_surfx5_eval_quad_node_44_r(fl); 
  } else { 
    fUpwindQuad_l[36] = hyb_3x3v_p1_surfx5_eval_quad_node_36_l(fc); 
    fUpwindQuad_l[37] = hyb_3x3v_p1_surfx5_eval_quad_node_37_l(fc); 
    fUpwindQuad_l[38] = hyb_3x3v_p1_surfx5_eval_quad_node_38_l(fc); 
    fUpwindQuad_l[39] = hyb_3x3v_p1_surfx5_eval_quad_node_39_l(fc); 
    fUpwindQuad_l[40] = hyb_3x3v_p1_surfx5_eval_quad_node_40_l(fc); 
    fUpwindQuad_l[41] = hyb_3x3v_p1_surfx5_eval_quad_node_41_l(fc); 
    fUpwindQuad_l[42] = hyb_3x3v_p1_surfx5_eval_quad_node_42_l(fc); 
    fUpwindQuad_l[43] = hyb_3x3v_p1_surfx5_eval_quad_node_43_l(fc); 
    fUpwindQuad_l[44] = hyb_3x3v_p1_surfx5_eval_quad_node_44_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[36] = hyb_3x3v_p1_surfx5_eval_quad_node_36_r(fc); 
    fUpwindQuad_r[37] = hyb_3x3v_p1_surfx5_eval_quad_node_37_r(fc); 
    fUpwindQuad_r[38] = hyb_3x3v_p1_surfx5_eval_quad_node_38_r(fc); 
    fUpwindQuad_r[39] = hyb_3x3v_p1_surfx5_eval_quad_node_39_r(fc); 
    fUpwindQuad_r[40] = hyb_3x3v_p1_surfx5_eval_quad_node_40_r(fc); 
    fUpwindQuad_r[41] = hyb_3x3v_p1_surfx5_eval_quad_node_41_r(fc); 
    fUpwindQuad_r[42] = hyb_3x3v_p1_surfx5_eval_quad_node_42_r(fc); 
    fUpwindQuad_r[43] = hyb_3x3v_p1_surfx5_eval_quad_node_43_r(fc); 
    fUpwindQuad_r[44] = hyb_3x3v_p1_surfx5_eval_quad_node_44_r(fc); 
  } else { 
    fUpwindQuad_r[36] = hyb_3x3v_p1_surfx5_eval_quad_node_36_l(fr); 
    fUpwindQuad_r[37] = hyb_3x3v_p1_surfx5_eval_quad_node_37_l(fr); 
    fUpwindQuad_r[38] = hyb_3x3v_p1_surfx5_eval_quad_node_38_l(fr); 
    fUpwindQuad_r[39] = hyb_3x3v_p1_surfx5_eval_quad_node_39_l(fr); 
    fUpwindQuad_r[40] = hyb_3x3v_p1_surfx5_eval_quad_node_40_l(fr); 
    fUpwindQuad_r[41] = hyb_3x3v_p1_surfx5_eval_quad_node_41_l(fr); 
    fUpwindQuad_r[42] = hyb_3x3v_p1_surfx5_eval_quad_node_42_l(fr); 
    fUpwindQuad_r[43] = hyb_3x3v_p1_surfx5_eval_quad_node_43_l(fr); 
    fUpwindQuad_r[44] = hyb_3x3v_p1_surfx5_eval_quad_node_44_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[45] = hyb_3x3v_p1_surfx5_eval_quad_node_45_r(fl); 
    fUpwindQuad_l[46] = hyb_3x3v_p1_surfx5_eval_quad_node_46_r(fl); 
    fUpwindQuad_l[47] = hyb_3x3v_p1_surfx5_eval_quad_node_47_r(fl); 
    fUpwindQuad_l[48] = hyb_3x3v_p1_surfx5_eval_quad_node_48_r(fl); 
    fUpwindQuad_l[49] = hyb_3x3v_p1_surfx5_eval_quad_node_49_r(fl); 
    fUpwindQuad_l[50] = hyb_3x3v_p1_surfx5_eval_quad_node_50_r(fl); 
    fUpwindQuad_l[51] = hyb_3x3v_p1_surfx5_eval_quad_node_51_r(fl); 
    fUpwindQuad_l[52] = hyb_3x3v_p1_surfx5_eval_quad_node_52_r(fl); 
    fUpwindQuad_l[53] = hyb_3x3v_p1_surfx5_eval_quad_node_53_r(fl); 
  } else { 
    fUpwindQuad_l[45] = hyb_3x3v_p1_surfx5_eval_quad_node_45_l(fc); 
    fUpwindQuad_l[46] = hyb_3x3v_p1_surfx5_eval_quad_node_46_l(fc); 
    fUpwindQuad_l[47] = hyb_3x3v_p1_surfx5_eval_quad_node_47_l(fc); 
    fUpwindQuad_l[48] = hyb_3x3v_p1_surfx5_eval_quad_node_48_l(fc); 
    fUpwindQuad_l[49] = hyb_3x3v_p1_surfx5_eval_quad_node_49_l(fc); 
    fUpwindQuad_l[50] = hyb_3x3v_p1_surfx5_eval_quad_node_50_l(fc); 
    fUpwindQuad_l[51] = hyb_3x3v_p1_surfx5_eval_quad_node_51_l(fc); 
    fUpwindQuad_l[52] = hyb_3x3v_p1_surfx5_eval_quad_node_52_l(fc); 
    fUpwindQuad_l[53] = hyb_3x3v_p1_surfx5_eval_quad_node_53_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[45] = hyb_3x3v_p1_surfx5_eval_quad_node_45_r(fc); 
    fUpwindQuad_r[46] = hyb_3x3v_p1_surfx5_eval_quad_node_46_r(fc); 
    fUpwindQuad_r[47] = hyb_3x3v_p1_surfx5_eval_quad_node_47_r(fc); 
    fUpwindQuad_r[48] = hyb_3x3v_p1_surfx5_eval_quad_node_48_r(fc); 
    fUpwindQuad_r[49] = hyb_3x3v_p1_surfx5_eval_quad_node_49_r(fc); 
    fUpwindQuad_r[50] = hyb_3x3v_p1_surfx5_eval_quad_node_50_r(fc); 
    fUpwindQuad_r[51] = hyb_3x3v_p1_surfx5_eval_quad_node_51_r(fc); 
    fUpwindQuad_r[52] = hyb_3x3v_p1_surfx5_eval_quad_node_52_r(fc); 
    fUpwindQuad_r[53] = hyb_3x3v_p1_surfx5_eval_quad_node_53_r(fc); 
  } else { 
    fUpwindQuad_r[45] = hyb_3x3v_p1_surfx5_eval_quad_node_45_l(fr); 
    fUpwindQuad_r[46] = hyb_3x3v_p1_surfx5_eval_quad_node_46_l(fr); 
    fUpwindQuad_r[47] = hyb_3x3v_p1_surfx5_eval_quad_node_47_l(fr); 
    fUpwindQuad_r[48] = hyb_3x3v_p1_surfx5_eval_quad_node_48_l(fr); 
    fUpwindQuad_r[49] = hyb_3x3v_p1_surfx5_eval_quad_node_49_l(fr); 
    fUpwindQuad_r[50] = hyb_3x3v_p1_surfx5_eval_quad_node_50_l(fr); 
    fUpwindQuad_r[51] = hyb_3x3v_p1_surfx5_eval_quad_node_51_l(fr); 
    fUpwindQuad_r[52] = hyb_3x3v_p1_surfx5_eval_quad_node_52_l(fr); 
    fUpwindQuad_r[53] = hyb_3x3v_p1_surfx5_eval_quad_node_53_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[54] = hyb_3x3v_p1_surfx5_eval_quad_node_54_r(fl); 
    fUpwindQuad_l[55] = hyb_3x3v_p1_surfx5_eval_quad_node_55_r(fl); 
    fUpwindQuad_l[56] = hyb_3x3v_p1_surfx5_eval_quad_node_56_r(fl); 
    fUpwindQuad_l[57] = hyb_3x3v_p1_surfx5_eval_quad_node_57_r(fl); 
    fUpwindQuad_l[58] = hyb_3x3v_p1_surfx5_eval_quad_node_58_r(fl); 
    fUpwindQuad_l[59] = hyb_3x3v_p1_surfx5_eval_quad_node_59_r(fl); 
    fUpwindQuad_l[60] = hyb_3x3v_p1_surfx5_eval_quad_node_60_r(fl); 
    fUpwindQuad_l[61] = hyb_3x3v_p1_surfx5_eval_quad_node_61_r(fl); 
    fUpwindQuad_l[62] = hyb_3x3v_p1_surfx5_eval_quad_node_62_r(fl); 
  } else { 
    fUpwindQuad_l[54] = hyb_3x3v_p1_surfx5_eval_quad_node_54_l(fc); 
    fUpwindQuad_l[55] = hyb_3x3v_p1_surfx5_eval_quad_node_55_l(fc); 
    fUpwindQuad_l[56] = hyb_3x3v_p1_surfx5_eval_quad_node_56_l(fc); 
    fUpwindQuad_l[57] = hyb_3x3v_p1_surfx5_eval_quad_node_57_l(fc); 
    fUpwindQuad_l[58] = hyb_3x3v_p1_surfx5_eval_quad_node_58_l(fc); 
    fUpwindQuad_l[59] = hyb_3x3v_p1_surfx5_eval_quad_node_59_l(fc); 
    fUpwindQuad_l[60] = hyb_3x3v_p1_surfx5_eval_quad_node_60_l(fc); 
    fUpwindQuad_l[61] = hyb_3x3v_p1_surfx5_eval_quad_node_61_l(fc); 
    fUpwindQuad_l[62] = hyb_3x3v_p1_surfx5_eval_quad_node_62_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[54] = hyb_3x3v_p1_surfx5_eval_quad_node_54_r(fc); 
    fUpwindQuad_r[55] = hyb_3x3v_p1_surfx5_eval_quad_node_55_r(fc); 
    fUpwindQuad_r[56] = hyb_3x3v_p1_surfx5_eval_quad_node_56_r(fc); 
    fUpwindQuad_r[57] = hyb_3x3v_p1_surfx5_eval_quad_node_57_r(fc); 
    fUpwindQuad_r[58] = hyb_3x3v_p1_surfx5_eval_quad_node_58_r(fc); 
    fUpwindQuad_r[59] = hyb_3x3v_p1_surfx5_eval_quad_node_59_r(fc); 
    fUpwindQuad_r[60] = hyb_3x3v_p1_surfx5_eval_quad_node_60_r(fc); 
    fUpwindQuad_r[61] = hyb_3x3v_p1_surfx5_eval_quad_node_61_r(fc); 
    fUpwindQuad_r[62] = hyb_3x3v_p1_surfx5_eval_quad_node_62_r(fc); 
  } else { 
    fUpwindQuad_r[54] = hyb_3x3v_p1_surfx5_eval_quad_node_54_l(fr); 
    fUpwindQuad_r[55] = hyb_3x3v_p1_surfx5_eval_quad_node_55_l(fr); 
    fUpwindQuad_r[56] = hyb_3x3v_p1_surfx5_eval_quad_node_56_l(fr); 
    fUpwindQuad_r[57] = hyb_3x3v_p1_surfx5_eval_quad_node_57_l(fr); 
    fUpwindQuad_r[58] = hyb_3x3v_p1_surfx5_eval_quad_node_58_l(fr); 
    fUpwindQuad_r[59] = hyb_3x3v_p1_surfx5_eval_quad_node_59_l(fr); 
    fUpwindQuad_r[60] = hyb_3x3v_p1_surfx5_eval_quad_node_60_l(fr); 
    fUpwindQuad_r[61] = hyb_3x3v_p1_surfx5_eval_quad_node_61_l(fr); 
    fUpwindQuad_r[62] = hyb_3x3v_p1_surfx5_eval_quad_node_62_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[63] = hyb_3x3v_p1_surfx5_eval_quad_node_63_r(fl); 
    fUpwindQuad_l[64] = hyb_3x3v_p1_surfx5_eval_quad_node_64_r(fl); 
    fUpwindQuad_l[65] = hyb_3x3v_p1_surfx5_eval_quad_node_65_r(fl); 
    fUpwindQuad_l[66] = hyb_3x3v_p1_surfx5_eval_quad_node_66_r(fl); 
    fUpwindQuad_l[67] = hyb_3x3v_p1_surfx5_eval_quad_node_67_r(fl); 
    fUpwindQuad_l[68] = hyb_3x3v_p1_surfx5_eval_quad_node_68_r(fl); 
    fUpwindQuad_l[69] = hyb_3x3v_p1_surfx5_eval_quad_node_69_r(fl); 
    fUpwindQuad_l[70] = hyb_3x3v_p1_surfx5_eval_quad_node_70_r(fl); 
    fUpwindQuad_l[71] = hyb_3x3v_p1_surfx5_eval_quad_node_71_r(fl); 
  } else { 
    fUpwindQuad_l[63] = hyb_3x3v_p1_surfx5_eval_quad_node_63_l(fc); 
    fUpwindQuad_l[64] = hyb_3x3v_p1_surfx5_eval_quad_node_64_l(fc); 
    fUpwindQuad_l[65] = hyb_3x3v_p1_surfx5_eval_quad_node_65_l(fc); 
    fUpwindQuad_l[66] = hyb_3x3v_p1_surfx5_eval_quad_node_66_l(fc); 
    fUpwindQuad_l[67] = hyb_3x3v_p1_surfx5_eval_quad_node_67_l(fc); 
    fUpwindQuad_l[68] = hyb_3x3v_p1_surfx5_eval_quad_node_68_l(fc); 
    fUpwindQuad_l[69] = hyb_3x3v_p1_surfx5_eval_quad_node_69_l(fc); 
    fUpwindQuad_l[70] = hyb_3x3v_p1_surfx5_eval_quad_node_70_l(fc); 
    fUpwindQuad_l[71] = hyb_3x3v_p1_surfx5_eval_quad_node_71_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[63] = hyb_3x3v_p1_surfx5_eval_quad_node_63_r(fc); 
    fUpwindQuad_r[64] = hyb_3x3v_p1_surfx5_eval_quad_node_64_r(fc); 
    fUpwindQuad_r[65] = hyb_3x3v_p1_surfx5_eval_quad_node_65_r(fc); 
    fUpwindQuad_r[66] = hyb_3x3v_p1_surfx5_eval_quad_node_66_r(fc); 
    fUpwindQuad_r[67] = hyb_3x3v_p1_surfx5_eval_quad_node_67_r(fc); 
    fUpwindQuad_r[68] = hyb_3x3v_p1_surfx5_eval_quad_node_68_r(fc); 
    fUpwindQuad_r[69] = hyb_3x3v_p1_surfx5_eval_quad_node_69_r(fc); 
    fUpwindQuad_r[70] = hyb_3x3v_p1_surfx5_eval_quad_node_70_r(fc); 
    fUpwindQuad_r[71] = hyb_3x3v_p1_surfx5_eval_quad_node_71_r(fc); 
  } else { 
    fUpwindQuad_r[63] = hyb_3x3v_p1_surfx5_eval_quad_node_63_l(fr); 
    fUpwindQuad_r[64] = hyb_3x3v_p1_surfx5_eval_quad_node_64_l(fr); 
    fUpwindQuad_r[65] = hyb_3x3v_p1_surfx5_eval_quad_node_65_l(fr); 
    fUpwindQuad_r[66] = hyb_3x3v_p1_surfx5_eval_quad_node_66_l(fr); 
    fUpwindQuad_r[67] = hyb_3x3v_p1_surfx5_eval_quad_node_67_l(fr); 
    fUpwindQuad_r[68] = hyb_3x3v_p1_surfx5_eval_quad_node_68_l(fr); 
    fUpwindQuad_r[69] = hyb_3x3v_p1_surfx5_eval_quad_node_69_l(fr); 
    fUpwindQuad_r[70] = hyb_3x3v_p1_surfx5_eval_quad_node_70_l(fr); 
    fUpwindQuad_r[71] = hyb_3x3v_p1_surfx5_eval_quad_node_71_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_3x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_3x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.1767766952966368*(alphaDrSurf_l[16]*fUpwind_l[16]+alphaDrSurf_l[8]*fUpwind_l[8]+alphaDrSurf_l[7]*fUpwind_l[7]+alphaDrSurf_l[6]*fUpwind_l[6]+alphaDrSurf_l[3]*fUpwind_l[3]+alphaDrSurf_l[2]*fUpwind_l[2]+alphaDrSurf_l[1]*fUpwind_l[1]+alphaDrSurf_l[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.1767766952966368*(alphaDrSurf_l[8]*fUpwind_l[16]+fUpwind_l[8]*alphaDrSurf_l[16]+alphaDrSurf_l[3]*fUpwind_l[7]+fUpwind_l[3]*alphaDrSurf_l[7]+alphaDrSurf_l[2]*fUpwind_l[6]+fUpwind_l[2]*alphaDrSurf_l[6]+alphaDrSurf_l[0]*fUpwind_l[1]+fUpwind_l[0]*alphaDrSurf_l[1]); 
  Ghat_l[2] = 0.1767766952966368*(alphaDrSurf_l[7]*fUpwind_l[16]+fUpwind_l[7]*alphaDrSurf_l[16]+alphaDrSurf_l[3]*fUpwind_l[8]+fUpwind_l[3]*alphaDrSurf_l[8]+alphaDrSurf_l[1]*fUpwind_l[6]+fUpwind_l[1]*alphaDrSurf_l[6]+alphaDrSurf_l[0]*fUpwind_l[2]+fUpwind_l[0]*alphaDrSurf_l[2]); 
  Ghat_l[3] = 0.1767766952966368*(alphaDrSurf_l[6]*fUpwind_l[16]+fUpwind_l[6]*alphaDrSurf_l[16]+alphaDrSurf_l[2]*fUpwind_l[8]+fUpwind_l[2]*alphaDrSurf_l[8]+alphaDrSurf_l[1]*fUpwind_l[7]+fUpwind_l[1]*alphaDrSurf_l[7]+alphaDrSurf_l[0]*fUpwind_l[3]+fUpwind_l[0]*alphaDrSurf_l[3]); 
  Ghat_l[4] = 0.1767766952966368*(alphaDrSurf_l[16]*fUpwind_l[26]+alphaDrSurf_l[8]*fUpwind_l[19]+alphaDrSurf_l[7]*fUpwind_l[18]+alphaDrSurf_l[6]*fUpwind_l[17]+alphaDrSurf_l[3]*fUpwind_l[11]+alphaDrSurf_l[2]*fUpwind_l[10]+alphaDrSurf_l[1]*fUpwind_l[9]+alphaDrSurf_l[0]*fUpwind_l[4]); 
  Ghat_l[5] = 0.1767766952966368*(alphaDrSurf_l[16]*fUpwind_l[27]+alphaDrSurf_l[8]*fUpwind_l[22]+alphaDrSurf_l[7]*fUpwind_l[21]+alphaDrSurf_l[6]*fUpwind_l[20]+alphaDrSurf_l[3]*fUpwind_l[14]+alphaDrSurf_l[2]*fUpwind_l[13]+alphaDrSurf_l[1]*fUpwind_l[12]+alphaDrSurf_l[0]*fUpwind_l[5]); 
  Ghat_l[6] = 0.1767766952966368*(alphaDrSurf_l[3]*fUpwind_l[16]+fUpwind_l[3]*alphaDrSurf_l[16]+alphaDrSurf_l[7]*fUpwind_l[8]+fUpwind_l[7]*alphaDrSurf_l[8]+alphaDrSurf_l[0]*fUpwind_l[6]+fUpwind_l[0]*alphaDrSurf_l[6]+alphaDrSurf_l[1]*fUpwind_l[2]+fUpwind_l[1]*alphaDrSurf_l[2]); 
  Ghat_l[7] = 0.1767766952966368*(alphaDrSurf_l[2]*fUpwind_l[16]+fUpwind_l[2]*alphaDrSurf_l[16]+alphaDrSurf_l[6]*fUpwind_l[8]+fUpwind_l[6]*alphaDrSurf_l[8]+alphaDrSurf_l[0]*fUpwind_l[7]+fUpwind_l[0]*alphaDrSurf_l[7]+alphaDrSurf_l[1]*fUpwind_l[3]+fUpwind_l[1]*alphaDrSurf_l[3]); 
  Ghat_l[8] = 0.1767766952966368*(alphaDrSurf_l[1]*fUpwind_l[16]+fUpwind_l[1]*alphaDrSurf_l[16]+alphaDrSurf_l[0]*fUpwind_l[8]+fUpwind_l[0]*alphaDrSurf_l[8]+alphaDrSurf_l[6]*fUpwind_l[7]+fUpwind_l[6]*alphaDrSurf_l[7]+alphaDrSurf_l[2]*fUpwind_l[3]+fUpwind_l[2]*alphaDrSurf_l[3]); 
  Ghat_l[9] = 0.1767766952966368*(alphaDrSurf_l[8]*fUpwind_l[26]+alphaDrSurf_l[16]*fUpwind_l[19]+alphaDrSurf_l[3]*fUpwind_l[18]+alphaDrSurf_l[2]*fUpwind_l[17]+alphaDrSurf_l[7]*fUpwind_l[11]+alphaDrSurf_l[6]*fUpwind_l[10]+alphaDrSurf_l[0]*fUpwind_l[9]+alphaDrSurf_l[1]*fUpwind_l[4]); 
  Ghat_l[10] = 0.1767766952966368*(alphaDrSurf_l[7]*fUpwind_l[26]+alphaDrSurf_l[3]*fUpwind_l[19]+alphaDrSurf_l[16]*fUpwind_l[18]+alphaDrSurf_l[1]*fUpwind_l[17]+alphaDrSurf_l[8]*fUpwind_l[11]+alphaDrSurf_l[0]*fUpwind_l[10]+alphaDrSurf_l[6]*fUpwind_l[9]+alphaDrSurf_l[2]*fUpwind_l[4]); 
  Ghat_l[11] = 0.1767766952966368*(alphaDrSurf_l[6]*fUpwind_l[26]+alphaDrSurf_l[2]*fUpwind_l[19]+alphaDrSurf_l[1]*fUpwind_l[18]+alphaDrSurf_l[16]*fUpwind_l[17]+alphaDrSurf_l[0]*fUpwind_l[11]+alphaDrSurf_l[8]*fUpwind_l[10]+alphaDrSurf_l[7]*fUpwind_l[9]+alphaDrSurf_l[3]*fUpwind_l[4]); 
  Ghat_l[12] = 0.1767766952966368*(alphaDrSurf_l[8]*fUpwind_l[27]+alphaDrSurf_l[16]*fUpwind_l[22]+alphaDrSurf_l[3]*fUpwind_l[21]+alphaDrSurf_l[2]*fUpwind_l[20]+alphaDrSurf_l[7]*fUpwind_l[14]+alphaDrSurf_l[6]*fUpwind_l[13]+alphaDrSurf_l[0]*fUpwind_l[12]+alphaDrSurf_l[1]*fUpwind_l[5]); 
  Ghat_l[13] = 0.1767766952966368*(alphaDrSurf_l[7]*fUpwind_l[27]+alphaDrSurf_l[3]*fUpwind_l[22]+alphaDrSurf_l[16]*fUpwind_l[21]+alphaDrSurf_l[1]*fUpwind_l[20]+alphaDrSurf_l[8]*fUpwind_l[14]+alphaDrSurf_l[0]*fUpwind_l[13]+alphaDrSurf_l[6]*fUpwind_l[12]+alphaDrSurf_l[2]*fUpwind_l[5]); 
  Ghat_l[14] = 0.1767766952966368*(alphaDrSurf_l[6]*fUpwind_l[27]+alphaDrSurf_l[2]*fUpwind_l[22]+alphaDrSurf_l[1]*fUpwind_l[21]+alphaDrSurf_l[16]*fUpwind_l[20]+alphaDrSurf_l[0]*fUpwind_l[14]+alphaDrSurf_l[8]*fUpwind_l[13]+alphaDrSurf_l[7]*fUpwind_l[12]+alphaDrSurf_l[3]*fUpwind_l[5]); 
  Ghat_l[15] = 0.1767766952966368*(alphaDrSurf_l[16]*fUpwind_l[31]+alphaDrSurf_l[8]*fUpwind_l[30]+alphaDrSurf_l[7]*fUpwind_l[29]+alphaDrSurf_l[6]*fUpwind_l[28]+alphaDrSurf_l[3]*fUpwind_l[25]+alphaDrSurf_l[2]*fUpwind_l[24]+alphaDrSurf_l[1]*fUpwind_l[23]+alphaDrSurf_l[0]*fUpwind_l[15]); 
  Ghat_l[16] = 0.1767766952966368*(alphaDrSurf_l[0]*fUpwind_l[16]+fUpwind_l[0]*alphaDrSurf_l[16]+alphaDrSurf_l[1]*fUpwind_l[8]+fUpwind_l[1]*alphaDrSurf_l[8]+alphaDrSurf_l[2]*fUpwind_l[7]+fUpwind_l[2]*alphaDrSurf_l[7]+alphaDrSurf_l[3]*fUpwind_l[6]+fUpwind_l[3]*alphaDrSurf_l[6]); 
  Ghat_l[17] = 0.1767766952966368*(alphaDrSurf_l[3]*fUpwind_l[26]+alphaDrSurf_l[7]*fUpwind_l[19]+alphaDrSurf_l[8]*fUpwind_l[18]+alphaDrSurf_l[0]*fUpwind_l[17]+fUpwind_l[11]*alphaDrSurf_l[16]+alphaDrSurf_l[1]*fUpwind_l[10]+alphaDrSurf_l[2]*fUpwind_l[9]+fUpwind_l[4]*alphaDrSurf_l[6]); 
  Ghat_l[18] = 0.1767766952966368*(alphaDrSurf_l[2]*fUpwind_l[26]+alphaDrSurf_l[6]*fUpwind_l[19]+alphaDrSurf_l[0]*fUpwind_l[18]+alphaDrSurf_l[8]*fUpwind_l[17]+fUpwind_l[10]*alphaDrSurf_l[16]+alphaDrSurf_l[1]*fUpwind_l[11]+alphaDrSurf_l[3]*fUpwind_l[9]+fUpwind_l[4]*alphaDrSurf_l[7]); 
  Ghat_l[19] = 0.1767766952966368*(alphaDrSurf_l[1]*fUpwind_l[26]+alphaDrSurf_l[0]*fUpwind_l[19]+alphaDrSurf_l[6]*fUpwind_l[18]+alphaDrSurf_l[7]*fUpwind_l[17]+fUpwind_l[9]*alphaDrSurf_l[16]+alphaDrSurf_l[2]*fUpwind_l[11]+alphaDrSurf_l[3]*fUpwind_l[10]+fUpwind_l[4]*alphaDrSurf_l[8]); 
  Ghat_l[20] = 0.1767766952966368*(alphaDrSurf_l[3]*fUpwind_l[27]+alphaDrSurf_l[7]*fUpwind_l[22]+alphaDrSurf_l[8]*fUpwind_l[21]+alphaDrSurf_l[0]*fUpwind_l[20]+fUpwind_l[14]*alphaDrSurf_l[16]+alphaDrSurf_l[1]*fUpwind_l[13]+alphaDrSurf_l[2]*fUpwind_l[12]+fUpwind_l[5]*alphaDrSurf_l[6]); 
  Ghat_l[21] = 0.1767766952966368*(alphaDrSurf_l[2]*fUpwind_l[27]+alphaDrSurf_l[6]*fUpwind_l[22]+alphaDrSurf_l[0]*fUpwind_l[21]+alphaDrSurf_l[8]*fUpwind_l[20]+fUpwind_l[13]*alphaDrSurf_l[16]+alphaDrSurf_l[1]*fUpwind_l[14]+alphaDrSurf_l[3]*fUpwind_l[12]+fUpwind_l[5]*alphaDrSurf_l[7]); 
  Ghat_l[22] = 0.1767766952966368*(alphaDrSurf_l[1]*fUpwind_l[27]+alphaDrSurf_l[0]*fUpwind_l[22]+alphaDrSurf_l[6]*fUpwind_l[21]+alphaDrSurf_l[7]*fUpwind_l[20]+fUpwind_l[12]*alphaDrSurf_l[16]+alphaDrSurf_l[2]*fUpwind_l[14]+alphaDrSurf_l[3]*fUpwind_l[13]+fUpwind_l[5]*alphaDrSurf_l[8]); 
  Ghat_l[23] = 0.1767766952966368*(alphaDrSurf_l[8]*fUpwind_l[31]+alphaDrSurf_l[16]*fUpwind_l[30]+alphaDrSurf_l[3]*fUpwind_l[29]+alphaDrSurf_l[2]*fUpwind_l[28]+alphaDrSurf_l[7]*fUpwind_l[25]+alphaDrSurf_l[6]*fUpwind_l[24]+alphaDrSurf_l[0]*fUpwind_l[23]+alphaDrSurf_l[1]*fUpwind_l[15]); 
  Ghat_l[24] = 0.1767766952966368*(alphaDrSurf_l[7]*fUpwind_l[31]+alphaDrSurf_l[3]*fUpwind_l[30]+alphaDrSurf_l[16]*fUpwind_l[29]+alphaDrSurf_l[1]*fUpwind_l[28]+alphaDrSurf_l[8]*fUpwind_l[25]+alphaDrSurf_l[0]*fUpwind_l[24]+alphaDrSurf_l[6]*fUpwind_l[23]+alphaDrSurf_l[2]*fUpwind_l[15]); 
  Ghat_l[25] = 0.1767766952966368*(alphaDrSurf_l[6]*fUpwind_l[31]+alphaDrSurf_l[2]*fUpwind_l[30]+alphaDrSurf_l[1]*fUpwind_l[29]+alphaDrSurf_l[16]*fUpwind_l[28]+alphaDrSurf_l[0]*fUpwind_l[25]+alphaDrSurf_l[8]*fUpwind_l[24]+alphaDrSurf_l[7]*fUpwind_l[23]+alphaDrSurf_l[3]*fUpwind_l[15]); 
  Ghat_l[26] = 0.1767766952966368*(alphaDrSurf_l[0]*fUpwind_l[26]+alphaDrSurf_l[1]*fUpwind_l[19]+alphaDrSurf_l[2]*fUpwind_l[18]+alphaDrSurf_l[3]*fUpwind_l[17]+fUpwind_l[4]*alphaDrSurf_l[16]+alphaDrSurf_l[6]*fUpwind_l[11]+alphaDrSurf_l[7]*fUpwind_l[10]+alphaDrSurf_l[8]*fUpwind_l[9]); 
  Ghat_l[27] = 0.1767766952966368*(alphaDrSurf_l[0]*fUpwind_l[27]+alphaDrSurf_l[1]*fUpwind_l[22]+alphaDrSurf_l[2]*fUpwind_l[21]+alphaDrSurf_l[3]*fUpwind_l[20]+fUpwind_l[5]*alphaDrSurf_l[16]+alphaDrSurf_l[6]*fUpwind_l[14]+alphaDrSurf_l[7]*fUpwind_l[13]+alphaDrSurf_l[8]*fUpwind_l[12]); 
  Ghat_l[28] = 0.1767766952966368*(alphaDrSurf_l[3]*fUpwind_l[31]+alphaDrSurf_l[7]*fUpwind_l[30]+alphaDrSurf_l[8]*fUpwind_l[29]+alphaDrSurf_l[0]*fUpwind_l[28]+alphaDrSurf_l[16]*fUpwind_l[25]+alphaDrSurf_l[1]*fUpwind_l[24]+alphaDrSurf_l[2]*fUpwind_l[23]+alphaDrSurf_l[6]*fUpwind_l[15]); 
  Ghat_l[29] = 0.1767766952966368*(alphaDrSurf_l[2]*fUpwind_l[31]+alphaDrSurf_l[6]*fUpwind_l[30]+alphaDrSurf_l[0]*fUpwind_l[29]+alphaDrSurf_l[8]*fUpwind_l[28]+alphaDrSurf_l[1]*fUpwind_l[25]+alphaDrSurf_l[16]*fUpwind_l[24]+alphaDrSurf_l[3]*fUpwind_l[23]+alphaDrSurf_l[7]*fUpwind_l[15]); 
  Ghat_l[30] = 0.1767766952966368*(alphaDrSurf_l[1]*fUpwind_l[31]+alphaDrSurf_l[0]*fUpwind_l[30]+alphaDrSurf_l[6]*fUpwind_l[29]+alphaDrSurf_l[7]*fUpwind_l[28]+alphaDrSurf_l[2]*fUpwind_l[25]+alphaDrSurf_l[3]*fUpwind_l[24]+alphaDrSurf_l[16]*fUpwind_l[23]+alphaDrSurf_l[8]*fUpwind_l[15]); 
  Ghat_l[31] = 0.1767766952966368*(alphaDrSurf_l[0]*fUpwind_l[31]+alphaDrSurf_l[1]*fUpwind_l[30]+alphaDrSurf_l[2]*fUpwind_l[29]+alphaDrSurf_l[3]*fUpwind_l[28]+alphaDrSurf_l[6]*fUpwind_l[25]+alphaDrSurf_l[7]*fUpwind_l[24]+alphaDrSurf_l[8]*fUpwind_l[23]+fUpwind_l[15]*alphaDrSurf_l[16]); 
  Ghat_l[32] = 0.01178511301977579*(15.0*alphaDrSurf_l[16]*fUpwind_l[43]+15.0*(alphaDrSurf_l[8]*fUpwind_l[39]+alphaDrSurf_l[7]*fUpwind_l[38]+alphaDrSurf_l[6]*fUpwind_l[37])+15.0*(alphaDrSurf_l[3]*fUpwind_l[35]+alphaDrSurf_l[2]*fUpwind_l[34]+alphaDrSurf_l[1]*fUpwind_l[33])+15.0*alphaDrSurf_l[0]*fUpwind_l[32]); 
  Ghat_l[33] = 0.01178511301977579*(15.0*alphaDrSurf_l[8]*fUpwind_l[43]+15.0*(alphaDrSurf_l[16]*fUpwind_l[39]+alphaDrSurf_l[3]*fUpwind_l[38]+alphaDrSurf_l[2]*fUpwind_l[37])+15.0*(alphaDrSurf_l[7]*fUpwind_l[35]+alphaDrSurf_l[6]*fUpwind_l[34]+alphaDrSurf_l[0]*fUpwind_l[33])+15.0*alphaDrSurf_l[1]*fUpwind_l[32]); 
  Ghat_l[34] = 0.01178511301977579*(15.0*alphaDrSurf_l[7]*fUpwind_l[43]+15.0*(alphaDrSurf_l[3]*fUpwind_l[39]+alphaDrSurf_l[16]*fUpwind_l[38]+alphaDrSurf_l[1]*fUpwind_l[37])+15.0*(alphaDrSurf_l[8]*fUpwind_l[35]+alphaDrSurf_l[0]*fUpwind_l[34]+alphaDrSurf_l[6]*fUpwind_l[33])+15.0*alphaDrSurf_l[2]*fUpwind_l[32]); 
  Ghat_l[35] = 0.01178511301977579*(15.0*alphaDrSurf_l[6]*fUpwind_l[43]+15.0*(alphaDrSurf_l[2]*fUpwind_l[39]+alphaDrSurf_l[1]*fUpwind_l[38]+alphaDrSurf_l[16]*fUpwind_l[37])+15.0*(alphaDrSurf_l[0]*fUpwind_l[35]+alphaDrSurf_l[8]*fUpwind_l[34]+alphaDrSurf_l[7]*fUpwind_l[33])+15.0*alphaDrSurf_l[3]*fUpwind_l[32]); 
  Ghat_l[36] = 0.01178511301977579*(15.0*alphaDrSurf_l[16]*fUpwind_l[47]+15.0*(alphaDrSurf_l[8]*fUpwind_l[46]+alphaDrSurf_l[7]*fUpwind_l[45]+alphaDrSurf_l[6]*fUpwind_l[44])+15.0*(alphaDrSurf_l[3]*fUpwind_l[42]+alphaDrSurf_l[2]*fUpwind_l[41]+alphaDrSurf_l[1]*fUpwind_l[40])+15.0*alphaDrSurf_l[0]*fUpwind_l[36]); 
  Ghat_l[37] = 0.01178511301977579*(15.0*alphaDrSurf_l[3]*fUpwind_l[43]+15.0*(alphaDrSurf_l[7]*fUpwind_l[39]+alphaDrSurf_l[8]*fUpwind_l[38]+alphaDrSurf_l[0]*fUpwind_l[37])+15.0*(alphaDrSurf_l[16]*fUpwind_l[35]+alphaDrSurf_l[1]*fUpwind_l[34]+alphaDrSurf_l[2]*fUpwind_l[33])+15.0*alphaDrSurf_l[6]*fUpwind_l[32]); 
  Ghat_l[38] = 0.01178511301977579*(15.0*alphaDrSurf_l[2]*fUpwind_l[43]+15.0*(alphaDrSurf_l[6]*fUpwind_l[39]+alphaDrSurf_l[0]*fUpwind_l[38]+alphaDrSurf_l[8]*fUpwind_l[37])+15.0*(alphaDrSurf_l[1]*fUpwind_l[35]+alphaDrSurf_l[16]*fUpwind_l[34]+alphaDrSurf_l[3]*fUpwind_l[33])+15.0*alphaDrSurf_l[7]*fUpwind_l[32]); 
  Ghat_l[39] = 0.01178511301977579*(15.0*alphaDrSurf_l[1]*fUpwind_l[43]+15.0*(alphaDrSurf_l[0]*fUpwind_l[39]+alphaDrSurf_l[6]*fUpwind_l[38]+alphaDrSurf_l[7]*fUpwind_l[37])+15.0*(alphaDrSurf_l[2]*fUpwind_l[35]+alphaDrSurf_l[3]*fUpwind_l[34]+alphaDrSurf_l[16]*fUpwind_l[33])+15.0*alphaDrSurf_l[8]*fUpwind_l[32]); 
  Ghat_l[40] = 0.01178511301977579*(15.0*alphaDrSurf_l[8]*fUpwind_l[47]+15.0*(alphaDrSurf_l[16]*fUpwind_l[46]+alphaDrSurf_l[3]*fUpwind_l[45]+alphaDrSurf_l[2]*fUpwind_l[44])+15.0*(alphaDrSurf_l[7]*fUpwind_l[42]+alphaDrSurf_l[6]*fUpwind_l[41]+alphaDrSurf_l[0]*fUpwind_l[40])+15.0*alphaDrSurf_l[1]*fUpwind_l[36]); 
  Ghat_l[41] = 0.01178511301977579*(15.0*alphaDrSurf_l[7]*fUpwind_l[47]+15.0*(alphaDrSurf_l[3]*fUpwind_l[46]+alphaDrSurf_l[16]*fUpwind_l[45]+alphaDrSurf_l[1]*fUpwind_l[44])+15.0*(alphaDrSurf_l[8]*fUpwind_l[42]+alphaDrSurf_l[0]*fUpwind_l[41]+alphaDrSurf_l[6]*fUpwind_l[40])+15.0*alphaDrSurf_l[2]*fUpwind_l[36]); 
  Ghat_l[42] = 0.01178511301977579*(15.0*alphaDrSurf_l[6]*fUpwind_l[47]+15.0*(alphaDrSurf_l[2]*fUpwind_l[46]+alphaDrSurf_l[1]*fUpwind_l[45]+alphaDrSurf_l[16]*fUpwind_l[44])+15.0*(alphaDrSurf_l[0]*fUpwind_l[42]+alphaDrSurf_l[8]*fUpwind_l[41]+alphaDrSurf_l[7]*fUpwind_l[40])+15.0*alphaDrSurf_l[3]*fUpwind_l[36]); 
  Ghat_l[43] = 0.01178511301977579*(15.0*alphaDrSurf_l[0]*fUpwind_l[43]+15.0*(alphaDrSurf_l[1]*fUpwind_l[39]+alphaDrSurf_l[2]*fUpwind_l[38]+alphaDrSurf_l[3]*fUpwind_l[37])+15.0*(alphaDrSurf_l[6]*fUpwind_l[35]+alphaDrSurf_l[7]*fUpwind_l[34]+alphaDrSurf_l[8]*fUpwind_l[33])+15.0*alphaDrSurf_l[16]*fUpwind_l[32]); 
  Ghat_l[44] = 0.01178511301977579*(15.0*alphaDrSurf_l[3]*fUpwind_l[47]+15.0*(alphaDrSurf_l[7]*fUpwind_l[46]+alphaDrSurf_l[8]*fUpwind_l[45]+alphaDrSurf_l[0]*fUpwind_l[44])+15.0*(alphaDrSurf_l[16]*fUpwind_l[42]+alphaDrSurf_l[1]*fUpwind_l[41]+alphaDrSurf_l[2]*fUpwind_l[40])+15.0*alphaDrSurf_l[6]*fUpwind_l[36]); 
  Ghat_l[45] = 0.01178511301977579*(15.0*alphaDrSurf_l[2]*fUpwind_l[47]+15.0*(alphaDrSurf_l[6]*fUpwind_l[46]+alphaDrSurf_l[0]*fUpwind_l[45]+alphaDrSurf_l[8]*fUpwind_l[44])+15.0*(alphaDrSurf_l[1]*fUpwind_l[42]+alphaDrSurf_l[16]*fUpwind_l[41]+alphaDrSurf_l[3]*fUpwind_l[40])+15.0*alphaDrSurf_l[7]*fUpwind_l[36]); 
  Ghat_l[46] = 0.01178511301977579*(15.0*alphaDrSurf_l[1]*fUpwind_l[47]+15.0*(alphaDrSurf_l[0]*fUpwind_l[46]+alphaDrSurf_l[6]*fUpwind_l[45]+alphaDrSurf_l[7]*fUpwind_l[44])+15.0*(alphaDrSurf_l[2]*fUpwind_l[42]+alphaDrSurf_l[3]*fUpwind_l[41]+alphaDrSurf_l[16]*fUpwind_l[40])+15.0*alphaDrSurf_l[8]*fUpwind_l[36]); 
  Ghat_l[47] = 0.01178511301977579*(15.0*alphaDrSurf_l[0]*fUpwind_l[47]+15.0*(alphaDrSurf_l[1]*fUpwind_l[46]+alphaDrSurf_l[2]*fUpwind_l[45]+alphaDrSurf_l[3]*fUpwind_l[44])+15.0*(alphaDrSurf_l[6]*fUpwind_l[42]+alphaDrSurf_l[7]*fUpwind_l[41]+alphaDrSurf_l[8]*fUpwind_l[40])+15.0*alphaDrSurf_l[16]*fUpwind_l[36]); 
  Ghat_l[48] = 0.01178511301977579*(15.0*alphaDrSurf_l[16]*fUpwind_l[59]+15.0*(alphaDrSurf_l[8]*fUpwind_l[55]+alphaDrSurf_l[7]*fUpwind_l[54]+alphaDrSurf_l[6]*fUpwind_l[53])+15.0*(alphaDrSurf_l[3]*fUpwind_l[51]+alphaDrSurf_l[2]*fUpwind_l[50]+alphaDrSurf_l[1]*fUpwind_l[49])+15.0*alphaDrSurf_l[0]*fUpwind_l[48]); 
  Ghat_l[49] = 0.01178511301977579*(15.0*alphaDrSurf_l[8]*fUpwind_l[59]+15.0*(alphaDrSurf_l[16]*fUpwind_l[55]+alphaDrSurf_l[3]*fUpwind_l[54]+alphaDrSurf_l[2]*fUpwind_l[53])+15.0*(alphaDrSurf_l[7]*fUpwind_l[51]+alphaDrSurf_l[6]*fUpwind_l[50]+alphaDrSurf_l[0]*fUpwind_l[49])+15.0*alphaDrSurf_l[1]*fUpwind_l[48]); 
  Ghat_l[50] = 0.01178511301977579*(15.0*alphaDrSurf_l[7]*fUpwind_l[59]+15.0*(alphaDrSurf_l[3]*fUpwind_l[55]+alphaDrSurf_l[16]*fUpwind_l[54]+alphaDrSurf_l[1]*fUpwind_l[53])+15.0*(alphaDrSurf_l[8]*fUpwind_l[51]+alphaDrSurf_l[0]*fUpwind_l[50]+alphaDrSurf_l[6]*fUpwind_l[49])+15.0*alphaDrSurf_l[2]*fUpwind_l[48]); 
  Ghat_l[51] = 0.01178511301977579*(15.0*alphaDrSurf_l[6]*fUpwind_l[59]+15.0*(alphaDrSurf_l[2]*fUpwind_l[55]+alphaDrSurf_l[1]*fUpwind_l[54]+alphaDrSurf_l[16]*fUpwind_l[53])+15.0*(alphaDrSurf_l[0]*fUpwind_l[51]+alphaDrSurf_l[8]*fUpwind_l[50]+alphaDrSurf_l[7]*fUpwind_l[49])+15.0*alphaDrSurf_l[3]*fUpwind_l[48]); 
  Ghat_l[52] = 0.01178511301977579*(15.0*alphaDrSurf_l[16]*fUpwind_l[63]+15.0*(alphaDrSurf_l[8]*fUpwind_l[62]+alphaDrSurf_l[7]*fUpwind_l[61]+alphaDrSurf_l[6]*fUpwind_l[60])+15.0*(alphaDrSurf_l[3]*fUpwind_l[58]+alphaDrSurf_l[2]*fUpwind_l[57]+alphaDrSurf_l[1]*fUpwind_l[56])+15.0*alphaDrSurf_l[0]*fUpwind_l[52]); 
  Ghat_l[53] = 0.01178511301977579*(15.0*alphaDrSurf_l[3]*fUpwind_l[59]+15.0*(alphaDrSurf_l[7]*fUpwind_l[55]+alphaDrSurf_l[8]*fUpwind_l[54]+alphaDrSurf_l[0]*fUpwind_l[53])+15.0*(alphaDrSurf_l[16]*fUpwind_l[51]+alphaDrSurf_l[1]*fUpwind_l[50]+alphaDrSurf_l[2]*fUpwind_l[49])+15.0*alphaDrSurf_l[6]*fUpwind_l[48]); 
  Ghat_l[54] = 0.01178511301977579*(15.0*alphaDrSurf_l[2]*fUpwind_l[59]+15.0*(alphaDrSurf_l[6]*fUpwind_l[55]+alphaDrSurf_l[0]*fUpwind_l[54]+alphaDrSurf_l[8]*fUpwind_l[53])+15.0*(alphaDrSurf_l[1]*fUpwind_l[51]+alphaDrSurf_l[16]*fUpwind_l[50]+alphaDrSurf_l[3]*fUpwind_l[49])+15.0*alphaDrSurf_l[7]*fUpwind_l[48]); 
  Ghat_l[55] = 0.01178511301977579*(15.0*alphaDrSurf_l[1]*fUpwind_l[59]+15.0*(alphaDrSurf_l[0]*fUpwind_l[55]+alphaDrSurf_l[6]*fUpwind_l[54]+alphaDrSurf_l[7]*fUpwind_l[53])+15.0*(alphaDrSurf_l[2]*fUpwind_l[51]+alphaDrSurf_l[3]*fUpwind_l[50]+alphaDrSurf_l[16]*fUpwind_l[49])+15.0*alphaDrSurf_l[8]*fUpwind_l[48]); 
  Ghat_l[56] = 0.01178511301977579*(15.0*alphaDrSurf_l[8]*fUpwind_l[63]+15.0*(alphaDrSurf_l[16]*fUpwind_l[62]+alphaDrSurf_l[3]*fUpwind_l[61]+alphaDrSurf_l[2]*fUpwind_l[60])+15.0*(alphaDrSurf_l[7]*fUpwind_l[58]+alphaDrSurf_l[6]*fUpwind_l[57]+alphaDrSurf_l[0]*fUpwind_l[56])+15.0*alphaDrSurf_l[1]*fUpwind_l[52]); 
  Ghat_l[57] = 0.01178511301977579*(15.0*alphaDrSurf_l[7]*fUpwind_l[63]+15.0*(alphaDrSurf_l[3]*fUpwind_l[62]+alphaDrSurf_l[16]*fUpwind_l[61]+alphaDrSurf_l[1]*fUpwind_l[60])+15.0*(alphaDrSurf_l[8]*fUpwind_l[58]+alphaDrSurf_l[0]*fUpwind_l[57]+alphaDrSurf_l[6]*fUpwind_l[56])+15.0*alphaDrSurf_l[2]*fUpwind_l[52]); 
  Ghat_l[58] = 0.01178511301977579*(15.0*alphaDrSurf_l[6]*fUpwind_l[63]+15.0*(alphaDrSurf_l[2]*fUpwind_l[62]+alphaDrSurf_l[1]*fUpwind_l[61]+alphaDrSurf_l[16]*fUpwind_l[60])+15.0*(alphaDrSurf_l[0]*fUpwind_l[58]+alphaDrSurf_l[8]*fUpwind_l[57]+alphaDrSurf_l[7]*fUpwind_l[56])+15.0*alphaDrSurf_l[3]*fUpwind_l[52]); 
  Ghat_l[59] = 0.01178511301977579*(15.0*alphaDrSurf_l[0]*fUpwind_l[59]+15.0*(alphaDrSurf_l[1]*fUpwind_l[55]+alphaDrSurf_l[2]*fUpwind_l[54]+alphaDrSurf_l[3]*fUpwind_l[53])+15.0*(alphaDrSurf_l[6]*fUpwind_l[51]+alphaDrSurf_l[7]*fUpwind_l[50]+alphaDrSurf_l[8]*fUpwind_l[49])+15.0*alphaDrSurf_l[16]*fUpwind_l[48]); 
  Ghat_l[60] = 0.01178511301977579*(15.0*alphaDrSurf_l[3]*fUpwind_l[63]+15.0*(alphaDrSurf_l[7]*fUpwind_l[62]+alphaDrSurf_l[8]*fUpwind_l[61]+alphaDrSurf_l[0]*fUpwind_l[60])+15.0*(alphaDrSurf_l[16]*fUpwind_l[58]+alphaDrSurf_l[1]*fUpwind_l[57]+alphaDrSurf_l[2]*fUpwind_l[56])+15.0*alphaDrSurf_l[6]*fUpwind_l[52]); 
  Ghat_l[61] = 0.01178511301977579*(15.0*alphaDrSurf_l[2]*fUpwind_l[63]+15.0*(alphaDrSurf_l[6]*fUpwind_l[62]+alphaDrSurf_l[0]*fUpwind_l[61]+alphaDrSurf_l[8]*fUpwind_l[60])+15.0*(alphaDrSurf_l[1]*fUpwind_l[58]+alphaDrSurf_l[16]*fUpwind_l[57]+alphaDrSurf_l[3]*fUpwind_l[56])+15.0*alphaDrSurf_l[7]*fUpwind_l[52]); 
  Ghat_l[62] = 0.01178511301977579*(15.0*alphaDrSurf_l[1]*fUpwind_l[63]+15.0*(alphaDrSurf_l[0]*fUpwind_l[62]+alphaDrSurf_l[6]*fUpwind_l[61]+alphaDrSurf_l[7]*fUpwind_l[60])+15.0*(alphaDrSurf_l[2]*fUpwind_l[58]+alphaDrSurf_l[3]*fUpwind_l[57]+alphaDrSurf_l[16]*fUpwind_l[56])+15.0*alphaDrSurf_l[8]*fUpwind_l[52]); 
  Ghat_l[63] = 0.01178511301977579*(15.0*alphaDrSurf_l[0]*fUpwind_l[63]+15.0*(alphaDrSurf_l[1]*fUpwind_l[62]+alphaDrSurf_l[2]*fUpwind_l[61]+alphaDrSurf_l[3]*fUpwind_l[60])+15.0*(alphaDrSurf_l[6]*fUpwind_l[58]+alphaDrSurf_l[7]*fUpwind_l[57]+alphaDrSurf_l[8]*fUpwind_l[56])+15.0*alphaDrSurf_l[16]*fUpwind_l[52]); 

  Ghat_r[0] = 0.1767766952966368*(alphaDrSurf_r[16]*fUpwind_r[16]+alphaDrSurf_r[8]*fUpwind_r[8]+alphaDrSurf_r[7]*fUpwind_r[7]+alphaDrSurf_r[6]*fUpwind_r[6]+alphaDrSurf_r[3]*fUpwind_r[3]+alphaDrSurf_r[2]*fUpwind_r[2]+alphaDrSurf_r[1]*fUpwind_r[1]+alphaDrSurf_r[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.1767766952966368*(alphaDrSurf_r[8]*fUpwind_r[16]+fUpwind_r[8]*alphaDrSurf_r[16]+alphaDrSurf_r[3]*fUpwind_r[7]+fUpwind_r[3]*alphaDrSurf_r[7]+alphaDrSurf_r[2]*fUpwind_r[6]+fUpwind_r[2]*alphaDrSurf_r[6]+alphaDrSurf_r[0]*fUpwind_r[1]+fUpwind_r[0]*alphaDrSurf_r[1]); 
  Ghat_r[2] = 0.1767766952966368*(alphaDrSurf_r[7]*fUpwind_r[16]+fUpwind_r[7]*alphaDrSurf_r[16]+alphaDrSurf_r[3]*fUpwind_r[8]+fUpwind_r[3]*alphaDrSurf_r[8]+alphaDrSurf_r[1]*fUpwind_r[6]+fUpwind_r[1]*alphaDrSurf_r[6]+alphaDrSurf_r[0]*fUpwind_r[2]+fUpwind_r[0]*alphaDrSurf_r[2]); 
  Ghat_r[3] = 0.1767766952966368*(alphaDrSurf_r[6]*fUpwind_r[16]+fUpwind_r[6]*alphaDrSurf_r[16]+alphaDrSurf_r[2]*fUpwind_r[8]+fUpwind_r[2]*alphaDrSurf_r[8]+alphaDrSurf_r[1]*fUpwind_r[7]+fUpwind_r[1]*alphaDrSurf_r[7]+alphaDrSurf_r[0]*fUpwind_r[3]+fUpwind_r[0]*alphaDrSurf_r[3]); 
  Ghat_r[4] = 0.1767766952966368*(alphaDrSurf_r[16]*fUpwind_r[26]+alphaDrSurf_r[8]*fUpwind_r[19]+alphaDrSurf_r[7]*fUpwind_r[18]+alphaDrSurf_r[6]*fUpwind_r[17]+alphaDrSurf_r[3]*fUpwind_r[11]+alphaDrSurf_r[2]*fUpwind_r[10]+alphaDrSurf_r[1]*fUpwind_r[9]+alphaDrSurf_r[0]*fUpwind_r[4]); 
  Ghat_r[5] = 0.1767766952966368*(alphaDrSurf_r[16]*fUpwind_r[27]+alphaDrSurf_r[8]*fUpwind_r[22]+alphaDrSurf_r[7]*fUpwind_r[21]+alphaDrSurf_r[6]*fUpwind_r[20]+alphaDrSurf_r[3]*fUpwind_r[14]+alphaDrSurf_r[2]*fUpwind_r[13]+alphaDrSurf_r[1]*fUpwind_r[12]+alphaDrSurf_r[0]*fUpwind_r[5]); 
  Ghat_r[6] = 0.1767766952966368*(alphaDrSurf_r[3]*fUpwind_r[16]+fUpwind_r[3]*alphaDrSurf_r[16]+alphaDrSurf_r[7]*fUpwind_r[8]+fUpwind_r[7]*alphaDrSurf_r[8]+alphaDrSurf_r[0]*fUpwind_r[6]+fUpwind_r[0]*alphaDrSurf_r[6]+alphaDrSurf_r[1]*fUpwind_r[2]+fUpwind_r[1]*alphaDrSurf_r[2]); 
  Ghat_r[7] = 0.1767766952966368*(alphaDrSurf_r[2]*fUpwind_r[16]+fUpwind_r[2]*alphaDrSurf_r[16]+alphaDrSurf_r[6]*fUpwind_r[8]+fUpwind_r[6]*alphaDrSurf_r[8]+alphaDrSurf_r[0]*fUpwind_r[7]+fUpwind_r[0]*alphaDrSurf_r[7]+alphaDrSurf_r[1]*fUpwind_r[3]+fUpwind_r[1]*alphaDrSurf_r[3]); 
  Ghat_r[8] = 0.1767766952966368*(alphaDrSurf_r[1]*fUpwind_r[16]+fUpwind_r[1]*alphaDrSurf_r[16]+alphaDrSurf_r[0]*fUpwind_r[8]+fUpwind_r[0]*alphaDrSurf_r[8]+alphaDrSurf_r[6]*fUpwind_r[7]+fUpwind_r[6]*alphaDrSurf_r[7]+alphaDrSurf_r[2]*fUpwind_r[3]+fUpwind_r[2]*alphaDrSurf_r[3]); 
  Ghat_r[9] = 0.1767766952966368*(alphaDrSurf_r[8]*fUpwind_r[26]+alphaDrSurf_r[16]*fUpwind_r[19]+alphaDrSurf_r[3]*fUpwind_r[18]+alphaDrSurf_r[2]*fUpwind_r[17]+alphaDrSurf_r[7]*fUpwind_r[11]+alphaDrSurf_r[6]*fUpwind_r[10]+alphaDrSurf_r[0]*fUpwind_r[9]+alphaDrSurf_r[1]*fUpwind_r[4]); 
  Ghat_r[10] = 0.1767766952966368*(alphaDrSurf_r[7]*fUpwind_r[26]+alphaDrSurf_r[3]*fUpwind_r[19]+alphaDrSurf_r[16]*fUpwind_r[18]+alphaDrSurf_r[1]*fUpwind_r[17]+alphaDrSurf_r[8]*fUpwind_r[11]+alphaDrSurf_r[0]*fUpwind_r[10]+alphaDrSurf_r[6]*fUpwind_r[9]+alphaDrSurf_r[2]*fUpwind_r[4]); 
  Ghat_r[11] = 0.1767766952966368*(alphaDrSurf_r[6]*fUpwind_r[26]+alphaDrSurf_r[2]*fUpwind_r[19]+alphaDrSurf_r[1]*fUpwind_r[18]+alphaDrSurf_r[16]*fUpwind_r[17]+alphaDrSurf_r[0]*fUpwind_r[11]+alphaDrSurf_r[8]*fUpwind_r[10]+alphaDrSurf_r[7]*fUpwind_r[9]+alphaDrSurf_r[3]*fUpwind_r[4]); 
  Ghat_r[12] = 0.1767766952966368*(alphaDrSurf_r[8]*fUpwind_r[27]+alphaDrSurf_r[16]*fUpwind_r[22]+alphaDrSurf_r[3]*fUpwind_r[21]+alphaDrSurf_r[2]*fUpwind_r[20]+alphaDrSurf_r[7]*fUpwind_r[14]+alphaDrSurf_r[6]*fUpwind_r[13]+alphaDrSurf_r[0]*fUpwind_r[12]+alphaDrSurf_r[1]*fUpwind_r[5]); 
  Ghat_r[13] = 0.1767766952966368*(alphaDrSurf_r[7]*fUpwind_r[27]+alphaDrSurf_r[3]*fUpwind_r[22]+alphaDrSurf_r[16]*fUpwind_r[21]+alphaDrSurf_r[1]*fUpwind_r[20]+alphaDrSurf_r[8]*fUpwind_r[14]+alphaDrSurf_r[0]*fUpwind_r[13]+alphaDrSurf_r[6]*fUpwind_r[12]+alphaDrSurf_r[2]*fUpwind_r[5]); 
  Ghat_r[14] = 0.1767766952966368*(alphaDrSurf_r[6]*fUpwind_r[27]+alphaDrSurf_r[2]*fUpwind_r[22]+alphaDrSurf_r[1]*fUpwind_r[21]+alphaDrSurf_r[16]*fUpwind_r[20]+alphaDrSurf_r[0]*fUpwind_r[14]+alphaDrSurf_r[8]*fUpwind_r[13]+alphaDrSurf_r[7]*fUpwind_r[12]+alphaDrSurf_r[3]*fUpwind_r[5]); 
  Ghat_r[15] = 0.1767766952966368*(alphaDrSurf_r[16]*fUpwind_r[31]+alphaDrSurf_r[8]*fUpwind_r[30]+alphaDrSurf_r[7]*fUpwind_r[29]+alphaDrSurf_r[6]*fUpwind_r[28]+alphaDrSurf_r[3]*fUpwind_r[25]+alphaDrSurf_r[2]*fUpwind_r[24]+alphaDrSurf_r[1]*fUpwind_r[23]+alphaDrSurf_r[0]*fUpwind_r[15]); 
  Ghat_r[16] = 0.1767766952966368*(alphaDrSurf_r[0]*fUpwind_r[16]+fUpwind_r[0]*alphaDrSurf_r[16]+alphaDrSurf_r[1]*fUpwind_r[8]+fUpwind_r[1]*alphaDrSurf_r[8]+alphaDrSurf_r[2]*fUpwind_r[7]+fUpwind_r[2]*alphaDrSurf_r[7]+alphaDrSurf_r[3]*fUpwind_r[6]+fUpwind_r[3]*alphaDrSurf_r[6]); 
  Ghat_r[17] = 0.1767766952966368*(alphaDrSurf_r[3]*fUpwind_r[26]+alphaDrSurf_r[7]*fUpwind_r[19]+alphaDrSurf_r[8]*fUpwind_r[18]+alphaDrSurf_r[0]*fUpwind_r[17]+fUpwind_r[11]*alphaDrSurf_r[16]+alphaDrSurf_r[1]*fUpwind_r[10]+alphaDrSurf_r[2]*fUpwind_r[9]+fUpwind_r[4]*alphaDrSurf_r[6]); 
  Ghat_r[18] = 0.1767766952966368*(alphaDrSurf_r[2]*fUpwind_r[26]+alphaDrSurf_r[6]*fUpwind_r[19]+alphaDrSurf_r[0]*fUpwind_r[18]+alphaDrSurf_r[8]*fUpwind_r[17]+fUpwind_r[10]*alphaDrSurf_r[16]+alphaDrSurf_r[1]*fUpwind_r[11]+alphaDrSurf_r[3]*fUpwind_r[9]+fUpwind_r[4]*alphaDrSurf_r[7]); 
  Ghat_r[19] = 0.1767766952966368*(alphaDrSurf_r[1]*fUpwind_r[26]+alphaDrSurf_r[0]*fUpwind_r[19]+alphaDrSurf_r[6]*fUpwind_r[18]+alphaDrSurf_r[7]*fUpwind_r[17]+fUpwind_r[9]*alphaDrSurf_r[16]+alphaDrSurf_r[2]*fUpwind_r[11]+alphaDrSurf_r[3]*fUpwind_r[10]+fUpwind_r[4]*alphaDrSurf_r[8]); 
  Ghat_r[20] = 0.1767766952966368*(alphaDrSurf_r[3]*fUpwind_r[27]+alphaDrSurf_r[7]*fUpwind_r[22]+alphaDrSurf_r[8]*fUpwind_r[21]+alphaDrSurf_r[0]*fUpwind_r[20]+fUpwind_r[14]*alphaDrSurf_r[16]+alphaDrSurf_r[1]*fUpwind_r[13]+alphaDrSurf_r[2]*fUpwind_r[12]+fUpwind_r[5]*alphaDrSurf_r[6]); 
  Ghat_r[21] = 0.1767766952966368*(alphaDrSurf_r[2]*fUpwind_r[27]+alphaDrSurf_r[6]*fUpwind_r[22]+alphaDrSurf_r[0]*fUpwind_r[21]+alphaDrSurf_r[8]*fUpwind_r[20]+fUpwind_r[13]*alphaDrSurf_r[16]+alphaDrSurf_r[1]*fUpwind_r[14]+alphaDrSurf_r[3]*fUpwind_r[12]+fUpwind_r[5]*alphaDrSurf_r[7]); 
  Ghat_r[22] = 0.1767766952966368*(alphaDrSurf_r[1]*fUpwind_r[27]+alphaDrSurf_r[0]*fUpwind_r[22]+alphaDrSurf_r[6]*fUpwind_r[21]+alphaDrSurf_r[7]*fUpwind_r[20]+fUpwind_r[12]*alphaDrSurf_r[16]+alphaDrSurf_r[2]*fUpwind_r[14]+alphaDrSurf_r[3]*fUpwind_r[13]+fUpwind_r[5]*alphaDrSurf_r[8]); 
  Ghat_r[23] = 0.1767766952966368*(alphaDrSurf_r[8]*fUpwind_r[31]+alphaDrSurf_r[16]*fUpwind_r[30]+alphaDrSurf_r[3]*fUpwind_r[29]+alphaDrSurf_r[2]*fUpwind_r[28]+alphaDrSurf_r[7]*fUpwind_r[25]+alphaDrSurf_r[6]*fUpwind_r[24]+alphaDrSurf_r[0]*fUpwind_r[23]+alphaDrSurf_r[1]*fUpwind_r[15]); 
  Ghat_r[24] = 0.1767766952966368*(alphaDrSurf_r[7]*fUpwind_r[31]+alphaDrSurf_r[3]*fUpwind_r[30]+alphaDrSurf_r[16]*fUpwind_r[29]+alphaDrSurf_r[1]*fUpwind_r[28]+alphaDrSurf_r[8]*fUpwind_r[25]+alphaDrSurf_r[0]*fUpwind_r[24]+alphaDrSurf_r[6]*fUpwind_r[23]+alphaDrSurf_r[2]*fUpwind_r[15]); 
  Ghat_r[25] = 0.1767766952966368*(alphaDrSurf_r[6]*fUpwind_r[31]+alphaDrSurf_r[2]*fUpwind_r[30]+alphaDrSurf_r[1]*fUpwind_r[29]+alphaDrSurf_r[16]*fUpwind_r[28]+alphaDrSurf_r[0]*fUpwind_r[25]+alphaDrSurf_r[8]*fUpwind_r[24]+alphaDrSurf_r[7]*fUpwind_r[23]+alphaDrSurf_r[3]*fUpwind_r[15]); 
  Ghat_r[26] = 0.1767766952966368*(alphaDrSurf_r[0]*fUpwind_r[26]+alphaDrSurf_r[1]*fUpwind_r[19]+alphaDrSurf_r[2]*fUpwind_r[18]+alphaDrSurf_r[3]*fUpwind_r[17]+fUpwind_r[4]*alphaDrSurf_r[16]+alphaDrSurf_r[6]*fUpwind_r[11]+alphaDrSurf_r[7]*fUpwind_r[10]+alphaDrSurf_r[8]*fUpwind_r[9]); 
  Ghat_r[27] = 0.1767766952966368*(alphaDrSurf_r[0]*fUpwind_r[27]+alphaDrSurf_r[1]*fUpwind_r[22]+alphaDrSurf_r[2]*fUpwind_r[21]+alphaDrSurf_r[3]*fUpwind_r[20]+fUpwind_r[5]*alphaDrSurf_r[16]+alphaDrSurf_r[6]*fUpwind_r[14]+alphaDrSurf_r[7]*fUpwind_r[13]+alphaDrSurf_r[8]*fUpwind_r[12]); 
  Ghat_r[28] = 0.1767766952966368*(alphaDrSurf_r[3]*fUpwind_r[31]+alphaDrSurf_r[7]*fUpwind_r[30]+alphaDrSurf_r[8]*fUpwind_r[29]+alphaDrSurf_r[0]*fUpwind_r[28]+alphaDrSurf_r[16]*fUpwind_r[25]+alphaDrSurf_r[1]*fUpwind_r[24]+alphaDrSurf_r[2]*fUpwind_r[23]+alphaDrSurf_r[6]*fUpwind_r[15]); 
  Ghat_r[29] = 0.1767766952966368*(alphaDrSurf_r[2]*fUpwind_r[31]+alphaDrSurf_r[6]*fUpwind_r[30]+alphaDrSurf_r[0]*fUpwind_r[29]+alphaDrSurf_r[8]*fUpwind_r[28]+alphaDrSurf_r[1]*fUpwind_r[25]+alphaDrSurf_r[16]*fUpwind_r[24]+alphaDrSurf_r[3]*fUpwind_r[23]+alphaDrSurf_r[7]*fUpwind_r[15]); 
  Ghat_r[30] = 0.1767766952966368*(alphaDrSurf_r[1]*fUpwind_r[31]+alphaDrSurf_r[0]*fUpwind_r[30]+alphaDrSurf_r[6]*fUpwind_r[29]+alphaDrSurf_r[7]*fUpwind_r[28]+alphaDrSurf_r[2]*fUpwind_r[25]+alphaDrSurf_r[3]*fUpwind_r[24]+alphaDrSurf_r[16]*fUpwind_r[23]+alphaDrSurf_r[8]*fUpwind_r[15]); 
  Ghat_r[31] = 0.1767766952966368*(alphaDrSurf_r[0]*fUpwind_r[31]+alphaDrSurf_r[1]*fUpwind_r[30]+alphaDrSurf_r[2]*fUpwind_r[29]+alphaDrSurf_r[3]*fUpwind_r[28]+alphaDrSurf_r[6]*fUpwind_r[25]+alphaDrSurf_r[7]*fUpwind_r[24]+alphaDrSurf_r[8]*fUpwind_r[23]+fUpwind_r[15]*alphaDrSurf_r[16]); 
  Ghat_r[32] = 0.01178511301977579*(15.0*alphaDrSurf_r[16]*fUpwind_r[43]+15.0*(alphaDrSurf_r[8]*fUpwind_r[39]+alphaDrSurf_r[7]*fUpwind_r[38]+alphaDrSurf_r[6]*fUpwind_r[37])+15.0*(alphaDrSurf_r[3]*fUpwind_r[35]+alphaDrSurf_r[2]*fUpwind_r[34]+alphaDrSurf_r[1]*fUpwind_r[33])+15.0*alphaDrSurf_r[0]*fUpwind_r[32]); 
  Ghat_r[33] = 0.01178511301977579*(15.0*alphaDrSurf_r[8]*fUpwind_r[43]+15.0*(alphaDrSurf_r[16]*fUpwind_r[39]+alphaDrSurf_r[3]*fUpwind_r[38]+alphaDrSurf_r[2]*fUpwind_r[37])+15.0*(alphaDrSurf_r[7]*fUpwind_r[35]+alphaDrSurf_r[6]*fUpwind_r[34]+alphaDrSurf_r[0]*fUpwind_r[33])+15.0*alphaDrSurf_r[1]*fUpwind_r[32]); 
  Ghat_r[34] = 0.01178511301977579*(15.0*alphaDrSurf_r[7]*fUpwind_r[43]+15.0*(alphaDrSurf_r[3]*fUpwind_r[39]+alphaDrSurf_r[16]*fUpwind_r[38]+alphaDrSurf_r[1]*fUpwind_r[37])+15.0*(alphaDrSurf_r[8]*fUpwind_r[35]+alphaDrSurf_r[0]*fUpwind_r[34]+alphaDrSurf_r[6]*fUpwind_r[33])+15.0*alphaDrSurf_r[2]*fUpwind_r[32]); 
  Ghat_r[35] = 0.01178511301977579*(15.0*alphaDrSurf_r[6]*fUpwind_r[43]+15.0*(alphaDrSurf_r[2]*fUpwind_r[39]+alphaDrSurf_r[1]*fUpwind_r[38]+alphaDrSurf_r[16]*fUpwind_r[37])+15.0*(alphaDrSurf_r[0]*fUpwind_r[35]+alphaDrSurf_r[8]*fUpwind_r[34]+alphaDrSurf_r[7]*fUpwind_r[33])+15.0*alphaDrSurf_r[3]*fUpwind_r[32]); 
  Ghat_r[36] = 0.01178511301977579*(15.0*alphaDrSurf_r[16]*fUpwind_r[47]+15.0*(alphaDrSurf_r[8]*fUpwind_r[46]+alphaDrSurf_r[7]*fUpwind_r[45]+alphaDrSurf_r[6]*fUpwind_r[44])+15.0*(alphaDrSurf_r[3]*fUpwind_r[42]+alphaDrSurf_r[2]*fUpwind_r[41]+alphaDrSurf_r[1]*fUpwind_r[40])+15.0*alphaDrSurf_r[0]*fUpwind_r[36]); 
  Ghat_r[37] = 0.01178511301977579*(15.0*alphaDrSurf_r[3]*fUpwind_r[43]+15.0*(alphaDrSurf_r[7]*fUpwind_r[39]+alphaDrSurf_r[8]*fUpwind_r[38]+alphaDrSurf_r[0]*fUpwind_r[37])+15.0*(alphaDrSurf_r[16]*fUpwind_r[35]+alphaDrSurf_r[1]*fUpwind_r[34]+alphaDrSurf_r[2]*fUpwind_r[33])+15.0*alphaDrSurf_r[6]*fUpwind_r[32]); 
  Ghat_r[38] = 0.01178511301977579*(15.0*alphaDrSurf_r[2]*fUpwind_r[43]+15.0*(alphaDrSurf_r[6]*fUpwind_r[39]+alphaDrSurf_r[0]*fUpwind_r[38]+alphaDrSurf_r[8]*fUpwind_r[37])+15.0*(alphaDrSurf_r[1]*fUpwind_r[35]+alphaDrSurf_r[16]*fUpwind_r[34]+alphaDrSurf_r[3]*fUpwind_r[33])+15.0*alphaDrSurf_r[7]*fUpwind_r[32]); 
  Ghat_r[39] = 0.01178511301977579*(15.0*alphaDrSurf_r[1]*fUpwind_r[43]+15.0*(alphaDrSurf_r[0]*fUpwind_r[39]+alphaDrSurf_r[6]*fUpwind_r[38]+alphaDrSurf_r[7]*fUpwind_r[37])+15.0*(alphaDrSurf_r[2]*fUpwind_r[35]+alphaDrSurf_r[3]*fUpwind_r[34]+alphaDrSurf_r[16]*fUpwind_r[33])+15.0*alphaDrSurf_r[8]*fUpwind_r[32]); 
  Ghat_r[40] = 0.01178511301977579*(15.0*alphaDrSurf_r[8]*fUpwind_r[47]+15.0*(alphaDrSurf_r[16]*fUpwind_r[46]+alphaDrSurf_r[3]*fUpwind_r[45]+alphaDrSurf_r[2]*fUpwind_r[44])+15.0*(alphaDrSurf_r[7]*fUpwind_r[42]+alphaDrSurf_r[6]*fUpwind_r[41]+alphaDrSurf_r[0]*fUpwind_r[40])+15.0*alphaDrSurf_r[1]*fUpwind_r[36]); 
  Ghat_r[41] = 0.01178511301977579*(15.0*alphaDrSurf_r[7]*fUpwind_r[47]+15.0*(alphaDrSurf_r[3]*fUpwind_r[46]+alphaDrSurf_r[16]*fUpwind_r[45]+alphaDrSurf_r[1]*fUpwind_r[44])+15.0*(alphaDrSurf_r[8]*fUpwind_r[42]+alphaDrSurf_r[0]*fUpwind_r[41]+alphaDrSurf_r[6]*fUpwind_r[40])+15.0*alphaDrSurf_r[2]*fUpwind_r[36]); 
  Ghat_r[42] = 0.01178511301977579*(15.0*alphaDrSurf_r[6]*fUpwind_r[47]+15.0*(alphaDrSurf_r[2]*fUpwind_r[46]+alphaDrSurf_r[1]*fUpwind_r[45]+alphaDrSurf_r[16]*fUpwind_r[44])+15.0*(alphaDrSurf_r[0]*fUpwind_r[42]+alphaDrSurf_r[8]*fUpwind_r[41]+alphaDrSurf_r[7]*fUpwind_r[40])+15.0*alphaDrSurf_r[3]*fUpwind_r[36]); 
  Ghat_r[43] = 0.01178511301977579*(15.0*alphaDrSurf_r[0]*fUpwind_r[43]+15.0*(alphaDrSurf_r[1]*fUpwind_r[39]+alphaDrSurf_r[2]*fUpwind_r[38]+alphaDrSurf_r[3]*fUpwind_r[37])+15.0*(alphaDrSurf_r[6]*fUpwind_r[35]+alphaDrSurf_r[7]*fUpwind_r[34]+alphaDrSurf_r[8]*fUpwind_r[33])+15.0*alphaDrSurf_r[16]*fUpwind_r[32]); 
  Ghat_r[44] = 0.01178511301977579*(15.0*alphaDrSurf_r[3]*fUpwind_r[47]+15.0*(alphaDrSurf_r[7]*fUpwind_r[46]+alphaDrSurf_r[8]*fUpwind_r[45]+alphaDrSurf_r[0]*fUpwind_r[44])+15.0*(alphaDrSurf_r[16]*fUpwind_r[42]+alphaDrSurf_r[1]*fUpwind_r[41]+alphaDrSurf_r[2]*fUpwind_r[40])+15.0*alphaDrSurf_r[6]*fUpwind_r[36]); 
  Ghat_r[45] = 0.01178511301977579*(15.0*alphaDrSurf_r[2]*fUpwind_r[47]+15.0*(alphaDrSurf_r[6]*fUpwind_r[46]+alphaDrSurf_r[0]*fUpwind_r[45]+alphaDrSurf_r[8]*fUpwind_r[44])+15.0*(alphaDrSurf_r[1]*fUpwind_r[42]+alphaDrSurf_r[16]*fUpwind_r[41]+alphaDrSurf_r[3]*fUpwind_r[40])+15.0*alphaDrSurf_r[7]*fUpwind_r[36]); 
  Ghat_r[46] = 0.01178511301977579*(15.0*alphaDrSurf_r[1]*fUpwind_r[47]+15.0*(alphaDrSurf_r[0]*fUpwind_r[46]+alphaDrSurf_r[6]*fUpwind_r[45]+alphaDrSurf_r[7]*fUpwind_r[44])+15.0*(alphaDrSurf_r[2]*fUpwind_r[42]+alphaDrSurf_r[3]*fUpwind_r[41]+alphaDrSurf_r[16]*fUpwind_r[40])+15.0*alphaDrSurf_r[8]*fUpwind_r[36]); 
  Ghat_r[47] = 0.01178511301977579*(15.0*alphaDrSurf_r[0]*fUpwind_r[47]+15.0*(alphaDrSurf_r[1]*fUpwind_r[46]+alphaDrSurf_r[2]*fUpwind_r[45]+alphaDrSurf_r[3]*fUpwind_r[44])+15.0*(alphaDrSurf_r[6]*fUpwind_r[42]+alphaDrSurf_r[7]*fUpwind_r[41]+alphaDrSurf_r[8]*fUpwind_r[40])+15.0*alphaDrSurf_r[16]*fUpwind_r[36]); 
  Ghat_r[48] = 0.01178511301977579*(15.0*alphaDrSurf_r[16]*fUpwind_r[59]+15.0*(alphaDrSurf_r[8]*fUpwind_r[55]+alphaDrSurf_r[7]*fUpwind_r[54]+alphaDrSurf_r[6]*fUpwind_r[53])+15.0*(alphaDrSurf_r[3]*fUpwind_r[51]+alphaDrSurf_r[2]*fUpwind_r[50]+alphaDrSurf_r[1]*fUpwind_r[49])+15.0*alphaDrSurf_r[0]*fUpwind_r[48]); 
  Ghat_r[49] = 0.01178511301977579*(15.0*alphaDrSurf_r[8]*fUpwind_r[59]+15.0*(alphaDrSurf_r[16]*fUpwind_r[55]+alphaDrSurf_r[3]*fUpwind_r[54]+alphaDrSurf_r[2]*fUpwind_r[53])+15.0*(alphaDrSurf_r[7]*fUpwind_r[51]+alphaDrSurf_r[6]*fUpwind_r[50]+alphaDrSurf_r[0]*fUpwind_r[49])+15.0*alphaDrSurf_r[1]*fUpwind_r[48]); 
  Ghat_r[50] = 0.01178511301977579*(15.0*alphaDrSurf_r[7]*fUpwind_r[59]+15.0*(alphaDrSurf_r[3]*fUpwind_r[55]+alphaDrSurf_r[16]*fUpwind_r[54]+alphaDrSurf_r[1]*fUpwind_r[53])+15.0*(alphaDrSurf_r[8]*fUpwind_r[51]+alphaDrSurf_r[0]*fUpwind_r[50]+alphaDrSurf_r[6]*fUpwind_r[49])+15.0*alphaDrSurf_r[2]*fUpwind_r[48]); 
  Ghat_r[51] = 0.01178511301977579*(15.0*alphaDrSurf_r[6]*fUpwind_r[59]+15.0*(alphaDrSurf_r[2]*fUpwind_r[55]+alphaDrSurf_r[1]*fUpwind_r[54]+alphaDrSurf_r[16]*fUpwind_r[53])+15.0*(alphaDrSurf_r[0]*fUpwind_r[51]+alphaDrSurf_r[8]*fUpwind_r[50]+alphaDrSurf_r[7]*fUpwind_r[49])+15.0*alphaDrSurf_r[3]*fUpwind_r[48]); 
  Ghat_r[52] = 0.01178511301977579*(15.0*alphaDrSurf_r[16]*fUpwind_r[63]+15.0*(alphaDrSurf_r[8]*fUpwind_r[62]+alphaDrSurf_r[7]*fUpwind_r[61]+alphaDrSurf_r[6]*fUpwind_r[60])+15.0*(alphaDrSurf_r[3]*fUpwind_r[58]+alphaDrSurf_r[2]*fUpwind_r[57]+alphaDrSurf_r[1]*fUpwind_r[56])+15.0*alphaDrSurf_r[0]*fUpwind_r[52]); 
  Ghat_r[53] = 0.01178511301977579*(15.0*alphaDrSurf_r[3]*fUpwind_r[59]+15.0*(alphaDrSurf_r[7]*fUpwind_r[55]+alphaDrSurf_r[8]*fUpwind_r[54]+alphaDrSurf_r[0]*fUpwind_r[53])+15.0*(alphaDrSurf_r[16]*fUpwind_r[51]+alphaDrSurf_r[1]*fUpwind_r[50]+alphaDrSurf_r[2]*fUpwind_r[49])+15.0*alphaDrSurf_r[6]*fUpwind_r[48]); 
  Ghat_r[54] = 0.01178511301977579*(15.0*alphaDrSurf_r[2]*fUpwind_r[59]+15.0*(alphaDrSurf_r[6]*fUpwind_r[55]+alphaDrSurf_r[0]*fUpwind_r[54]+alphaDrSurf_r[8]*fUpwind_r[53])+15.0*(alphaDrSurf_r[1]*fUpwind_r[51]+alphaDrSurf_r[16]*fUpwind_r[50]+alphaDrSurf_r[3]*fUpwind_r[49])+15.0*alphaDrSurf_r[7]*fUpwind_r[48]); 
  Ghat_r[55] = 0.01178511301977579*(15.0*alphaDrSurf_r[1]*fUpwind_r[59]+15.0*(alphaDrSurf_r[0]*fUpwind_r[55]+alphaDrSurf_r[6]*fUpwind_r[54]+alphaDrSurf_r[7]*fUpwind_r[53])+15.0*(alphaDrSurf_r[2]*fUpwind_r[51]+alphaDrSurf_r[3]*fUpwind_r[50]+alphaDrSurf_r[16]*fUpwind_r[49])+15.0*alphaDrSurf_r[8]*fUpwind_r[48]); 
  Ghat_r[56] = 0.01178511301977579*(15.0*alphaDrSurf_r[8]*fUpwind_r[63]+15.0*(alphaDrSurf_r[16]*fUpwind_r[62]+alphaDrSurf_r[3]*fUpwind_r[61]+alphaDrSurf_r[2]*fUpwind_r[60])+15.0*(alphaDrSurf_r[7]*fUpwind_r[58]+alphaDrSurf_r[6]*fUpwind_r[57]+alphaDrSurf_r[0]*fUpwind_r[56])+15.0*alphaDrSurf_r[1]*fUpwind_r[52]); 
  Ghat_r[57] = 0.01178511301977579*(15.0*alphaDrSurf_r[7]*fUpwind_r[63]+15.0*(alphaDrSurf_r[3]*fUpwind_r[62]+alphaDrSurf_r[16]*fUpwind_r[61]+alphaDrSurf_r[1]*fUpwind_r[60])+15.0*(alphaDrSurf_r[8]*fUpwind_r[58]+alphaDrSurf_r[0]*fUpwind_r[57]+alphaDrSurf_r[6]*fUpwind_r[56])+15.0*alphaDrSurf_r[2]*fUpwind_r[52]); 
  Ghat_r[58] = 0.01178511301977579*(15.0*alphaDrSurf_r[6]*fUpwind_r[63]+15.0*(alphaDrSurf_r[2]*fUpwind_r[62]+alphaDrSurf_r[1]*fUpwind_r[61]+alphaDrSurf_r[16]*fUpwind_r[60])+15.0*(alphaDrSurf_r[0]*fUpwind_r[58]+alphaDrSurf_r[8]*fUpwind_r[57]+alphaDrSurf_r[7]*fUpwind_r[56])+15.0*alphaDrSurf_r[3]*fUpwind_r[52]); 
  Ghat_r[59] = 0.01178511301977579*(15.0*alphaDrSurf_r[0]*fUpwind_r[59]+15.0*(alphaDrSurf_r[1]*fUpwind_r[55]+alphaDrSurf_r[2]*fUpwind_r[54]+alphaDrSurf_r[3]*fUpwind_r[53])+15.0*(alphaDrSurf_r[6]*fUpwind_r[51]+alphaDrSurf_r[7]*fUpwind_r[50]+alphaDrSurf_r[8]*fUpwind_r[49])+15.0*alphaDrSurf_r[16]*fUpwind_r[48]); 
  Ghat_r[60] = 0.01178511301977579*(15.0*alphaDrSurf_r[3]*fUpwind_r[63]+15.0*(alphaDrSurf_r[7]*fUpwind_r[62]+alphaDrSurf_r[8]*fUpwind_r[61]+alphaDrSurf_r[0]*fUpwind_r[60])+15.0*(alphaDrSurf_r[16]*fUpwind_r[58]+alphaDrSurf_r[1]*fUpwind_r[57]+alphaDrSurf_r[2]*fUpwind_r[56])+15.0*alphaDrSurf_r[6]*fUpwind_r[52]); 
  Ghat_r[61] = 0.01178511301977579*(15.0*alphaDrSurf_r[2]*fUpwind_r[63]+15.0*(alphaDrSurf_r[6]*fUpwind_r[62]+alphaDrSurf_r[0]*fUpwind_r[61]+alphaDrSurf_r[8]*fUpwind_r[60])+15.0*(alphaDrSurf_r[1]*fUpwind_r[58]+alphaDrSurf_r[16]*fUpwind_r[57]+alphaDrSurf_r[3]*fUpwind_r[56])+15.0*alphaDrSurf_r[7]*fUpwind_r[52]); 
  Ghat_r[62] = 0.01178511301977579*(15.0*alphaDrSurf_r[1]*fUpwind_r[63]+15.0*(alphaDrSurf_r[0]*fUpwind_r[62]+alphaDrSurf_r[6]*fUpwind_r[61]+alphaDrSurf_r[7]*fUpwind_r[60])+15.0*(alphaDrSurf_r[2]*fUpwind_r[58]+alphaDrSurf_r[3]*fUpwind_r[57]+alphaDrSurf_r[16]*fUpwind_r[56])+15.0*alphaDrSurf_r[8]*fUpwind_r[52]); 
  Ghat_r[63] = 0.01178511301977579*(15.0*alphaDrSurf_r[0]*fUpwind_r[63]+15.0*(alphaDrSurf_r[1]*fUpwind_r[62]+alphaDrSurf_r[2]*fUpwind_r[61]+alphaDrSurf_r[3]*fUpwind_r[60])+15.0*(alphaDrSurf_r[6]*fUpwind_r[58]+alphaDrSurf_r[7]*fUpwind_r[57]+alphaDrSurf_r[8]*fUpwind_r[56])+15.0*alphaDrSurf_r[16]*fUpwind_r[52]); 

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[5] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[9] += 0.7071067811865475*Ghat_r[8]*rdv2-0.7071067811865475*Ghat_l[8]*rdv2; 
  out[10] += 0.7071067811865475*Ghat_r[9]*rdv2-0.7071067811865475*Ghat_l[9]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_r[10]*rdv2-0.7071067811865475*Ghat_l[10]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_r[11]*rdv2-0.7071067811865475*Ghat_l[11]*rdv2; 
  out[13] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[14] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[15] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[16] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[17] += 0.7071067811865475*Ghat_r[12]*rdv2-0.7071067811865475*Ghat_l[12]*rdv2; 
  out[18] += 0.7071067811865475*Ghat_r[13]*rdv2-0.7071067811865475*Ghat_l[13]*rdv2; 
  out[19] += 0.7071067811865475*Ghat_r[14]*rdv2-0.7071067811865475*Ghat_l[14]*rdv2; 
  out[20] += 0.7071067811865475*Ghat_r[15]*rdv2-0.7071067811865475*Ghat_l[15]*rdv2; 
  out[21] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[22] += 0.7071067811865475*Ghat_r[16]*rdv2-0.7071067811865475*Ghat_l[16]*rdv2; 
  out[23] += 0.7071067811865475*Ghat_r[17]*rdv2-0.7071067811865475*Ghat_l[17]*rdv2; 
  out[24] += 0.7071067811865475*Ghat_r[18]*rdv2-0.7071067811865475*Ghat_l[18]*rdv2; 
  out[25] += 0.7071067811865475*Ghat_r[19]*rdv2-0.7071067811865475*Ghat_l[19]*rdv2; 
  out[26] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[27] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[28] += 1.224744871391589*Ghat_r[8]*rdv2+1.224744871391589*Ghat_l[8]*rdv2; 
  out[29] += 1.224744871391589*Ghat_r[9]*rdv2+1.224744871391589*Ghat_l[9]*rdv2; 
  out[30] += 1.224744871391589*Ghat_r[10]*rdv2+1.224744871391589*Ghat_l[10]*rdv2; 
  out[31] += 1.224744871391589*Ghat_r[11]*rdv2+1.224744871391589*Ghat_l[11]*rdv2; 
  out[32] += 0.7071067811865475*Ghat_r[20]*rdv2-0.7071067811865475*Ghat_l[20]*rdv2; 
  out[33] += 0.7071067811865475*Ghat_r[21]*rdv2-0.7071067811865475*Ghat_l[21]*rdv2; 
  out[34] += 0.7071067811865475*Ghat_r[22]*rdv2-0.7071067811865475*Ghat_l[22]*rdv2; 
  out[35] += 0.7071067811865475*Ghat_r[23]*rdv2-0.7071067811865475*Ghat_l[23]*rdv2; 
  out[36] += 0.7071067811865475*Ghat_r[24]*rdv2-0.7071067811865475*Ghat_l[24]*rdv2; 
  out[37] += 0.7071067811865475*Ghat_r[25]*rdv2-0.7071067811865475*Ghat_l[25]*rdv2; 
  out[38] += 1.224744871391589*Ghat_r[12]*rdv2+1.224744871391589*Ghat_l[12]*rdv2; 
  out[39] += 1.224744871391589*Ghat_r[13]*rdv2+1.224744871391589*Ghat_l[13]*rdv2; 
  out[40] += 1.224744871391589*Ghat_r[14]*rdv2+1.224744871391589*Ghat_l[14]*rdv2; 
  out[41] += 1.224744871391589*Ghat_r[15]*rdv2+1.224744871391589*Ghat_l[15]*rdv2; 
  out[42] += 0.7071067811865475*Ghat_r[26]*rdv2-0.7071067811865475*Ghat_l[26]*rdv2; 
  out[43] += 1.224744871391589*Ghat_r[16]*rdv2+1.224744871391589*Ghat_l[16]*rdv2; 
  out[44] += 1.224744871391589*Ghat_r[17]*rdv2+1.224744871391589*Ghat_l[17]*rdv2; 
  out[45] += 1.224744871391589*Ghat_r[18]*rdv2+1.224744871391589*Ghat_l[18]*rdv2; 
  out[46] += 1.224744871391589*Ghat_r[19]*rdv2+1.224744871391589*Ghat_l[19]*rdv2; 
  out[47] += 0.7071067811865475*Ghat_r[27]*rdv2-0.7071067811865475*Ghat_l[27]*rdv2; 
  out[48] += 0.7071067811865475*Ghat_r[28]*rdv2-0.7071067811865475*Ghat_l[28]*rdv2; 
  out[49] += 0.7071067811865475*Ghat_r[29]*rdv2-0.7071067811865475*Ghat_l[29]*rdv2; 
  out[50] += 0.7071067811865475*Ghat_r[30]*rdv2-0.7071067811865475*Ghat_l[30]*rdv2; 
  out[51] += 1.224744871391589*Ghat_r[20]*rdv2+1.224744871391589*Ghat_l[20]*rdv2; 
  out[52] += 1.224744871391589*Ghat_r[21]*rdv2+1.224744871391589*Ghat_l[21]*rdv2; 
  out[53] += 1.224744871391589*Ghat_r[22]*rdv2+1.224744871391589*Ghat_l[22]*rdv2; 
  out[54] += 1.224744871391589*Ghat_r[23]*rdv2+1.224744871391589*Ghat_l[23]*rdv2; 
  out[55] += 1.224744871391589*Ghat_r[24]*rdv2+1.224744871391589*Ghat_l[24]*rdv2; 
  out[56] += 1.224744871391589*Ghat_r[25]*rdv2+1.224744871391589*Ghat_l[25]*rdv2; 
  out[57] += 1.224744871391589*Ghat_r[26]*rdv2+1.224744871391589*Ghat_l[26]*rdv2; 
  out[58] += 0.7071067811865475*Ghat_r[31]*rdv2-0.7071067811865475*Ghat_l[31]*rdv2; 
  out[59] += 1.224744871391589*Ghat_r[27]*rdv2+1.224744871391589*Ghat_l[27]*rdv2; 
  out[60] += 1.224744871391589*Ghat_r[28]*rdv2+1.224744871391589*Ghat_l[28]*rdv2; 
  out[61] += 1.224744871391589*Ghat_r[29]*rdv2+1.224744871391589*Ghat_l[29]*rdv2; 
  out[62] += 1.224744871391589*Ghat_r[30]*rdv2+1.224744871391589*Ghat_l[30]*rdv2; 
  out[63] += 1.224744871391589*Ghat_r[31]*rdv2+1.224744871391589*Ghat_l[31]*rdv2; 
  out[64] += 0.7071067811865475*Ghat_r[32]*rdv2-0.7071067811865475*Ghat_l[32]*rdv2; 
  out[65] += 0.7071067811865475*Ghat_r[33]*rdv2-0.7071067811865475*Ghat_l[33]*rdv2; 
  out[66] += 0.7071067811865475*Ghat_r[34]*rdv2-0.7071067811865475*Ghat_l[34]*rdv2; 
  out[67] += 0.7071067811865475*Ghat_r[35]*rdv2-0.7071067811865475*Ghat_l[35]*rdv2; 
  out[68] += 1.224744871391589*Ghat_r[32]*rdv2+1.224744871391589*Ghat_l[32]*rdv2; 
  out[69] += 0.7071067811865475*Ghat_r[36]*rdv2-0.7071067811865475*Ghat_l[36]*rdv2; 
  out[70] += 0.7071067811865475*Ghat_r[37]*rdv2-0.7071067811865475*Ghat_l[37]*rdv2; 
  out[71] += 0.7071067811865475*Ghat_r[38]*rdv2-0.7071067811865475*Ghat_l[38]*rdv2; 
  out[72] += 0.7071067811865475*Ghat_r[39]*rdv2-0.7071067811865475*Ghat_l[39]*rdv2; 
  out[73] += 1.224744871391589*Ghat_r[33]*rdv2+1.224744871391589*Ghat_l[33]*rdv2; 
  out[74] += 1.224744871391589*Ghat_r[34]*rdv2+1.224744871391589*Ghat_l[34]*rdv2; 
  out[75] += 1.224744871391589*Ghat_r[35]*rdv2+1.224744871391589*Ghat_l[35]*rdv2; 
  out[76] += 0.7071067811865475*Ghat_r[40]*rdv2-0.7071067811865475*Ghat_l[40]*rdv2; 
  out[77] += 0.7071067811865475*Ghat_r[41]*rdv2-0.7071067811865475*Ghat_l[41]*rdv2; 
  out[78] += 0.7071067811865475*Ghat_r[42]*rdv2-0.7071067811865475*Ghat_l[42]*rdv2; 
  out[79] += 1.224744871391589*Ghat_r[36]*rdv2+1.224744871391589*Ghat_l[36]*rdv2; 
  out[80] += 0.7071067811865475*Ghat_r[43]*rdv2-0.7071067811865475*Ghat_l[43]*rdv2; 
  out[81] += 1.224744871391589*Ghat_r[37]*rdv2+1.224744871391589*Ghat_l[37]*rdv2; 
  out[82] += 1.224744871391589*Ghat_r[38]*rdv2+1.224744871391589*Ghat_l[38]*rdv2; 
  out[83] += 1.224744871391589*Ghat_r[39]*rdv2+1.224744871391589*Ghat_l[39]*rdv2; 
  out[84] += 0.7071067811865475*Ghat_r[44]*rdv2-0.7071067811865475*Ghat_l[44]*rdv2; 
  out[85] += 0.7071067811865475*Ghat_r[45]*rdv2-0.7071067811865475*Ghat_l[45]*rdv2; 
  out[86] += 0.7071067811865475*Ghat_r[46]*rdv2-0.7071067811865475*Ghat_l[46]*rdv2; 
  out[87] += 1.224744871391589*Ghat_r[40]*rdv2+1.224744871391589*Ghat_l[40]*rdv2; 
  out[88] += 1.224744871391589*Ghat_r[41]*rdv2+1.224744871391589*Ghat_l[41]*rdv2; 
  out[89] += 1.224744871391589*Ghat_r[42]*rdv2+1.224744871391589*Ghat_l[42]*rdv2; 
  out[90] += 1.224744871391589*Ghat_r[43]*rdv2+1.224744871391589*Ghat_l[43]*rdv2; 
  out[91] += 0.7071067811865475*Ghat_r[47]*rdv2-0.7071067811865475*Ghat_l[47]*rdv2; 
  out[92] += 1.224744871391589*Ghat_r[44]*rdv2+1.224744871391589*Ghat_l[44]*rdv2; 
  out[93] += 1.224744871391589*Ghat_r[45]*rdv2+1.224744871391589*Ghat_l[45]*rdv2; 
  out[94] += 1.224744871391589*Ghat_r[46]*rdv2+1.224744871391589*Ghat_l[46]*rdv2; 
  out[95] += 1.224744871391589*Ghat_r[47]*rdv2+1.224744871391589*Ghat_l[47]*rdv2; 
  out[96] += 1.58113883008419*Ghat_r[0]*rdv2-1.58113883008419*Ghat_l[0]*rdv2; 
  out[97] += 1.58113883008419*Ghat_r[1]*rdv2-1.58113883008419*Ghat_l[1]*rdv2; 
  out[98] += 1.58113883008419*Ghat_r[2]*rdv2-1.58113883008419*Ghat_l[2]*rdv2; 
  out[99] += 1.58113883008419*Ghat_r[3]*rdv2-1.58113883008419*Ghat_l[3]*rdv2; 
  out[100] += 1.58113883008419*Ghat_r[4]*rdv2-1.58113883008419*Ghat_l[4]*rdv2; 
  out[101] += 1.58113883008419*Ghat_r[5]*rdv2-1.58113883008419*Ghat_l[5]*rdv2; 
  out[102] += 1.58113883008419*Ghat_r[6]*rdv2-1.58113883008419*Ghat_l[6]*rdv2; 
  out[103] += 1.58113883008419*Ghat_r[7]*rdv2-1.58113883008419*Ghat_l[7]*rdv2; 
  out[104] += 1.58113883008419*Ghat_r[8]*rdv2-1.58113883008419*Ghat_l[8]*rdv2; 
  out[105] += 1.58113883008419*Ghat_r[9]*rdv2-1.58113883008419*Ghat_l[9]*rdv2; 
  out[106] += 1.58113883008419*Ghat_r[10]*rdv2-1.58113883008419*Ghat_l[10]*rdv2; 
  out[107] += 1.58113883008419*Ghat_r[11]*rdv2-1.58113883008419*Ghat_l[11]*rdv2; 
  out[108] += 1.58113883008419*Ghat_r[12]*rdv2-1.58113883008419*Ghat_l[12]*rdv2; 
  out[109] += 1.58113883008419*Ghat_r[13]*rdv2-1.58113883008419*Ghat_l[13]*rdv2; 
  out[110] += 1.58113883008419*Ghat_r[14]*rdv2-1.58113883008419*Ghat_l[14]*rdv2; 
  out[111] += 1.58113883008419*Ghat_r[15]*rdv2-1.58113883008419*Ghat_l[15]*rdv2; 
  out[112] += 1.58113883008419*Ghat_r[16]*rdv2-1.58113883008419*Ghat_l[16]*rdv2; 
  out[113] += 1.58113883008419*Ghat_r[17]*rdv2-1.58113883008419*Ghat_l[17]*rdv2; 
  out[114] += 1.58113883008419*Ghat_r[18]*rdv2-1.58113883008419*Ghat_l[18]*rdv2; 
  out[115] += 1.58113883008419*Ghat_r[19]*rdv2-1.58113883008419*Ghat_l[19]*rdv2; 
  out[116] += 1.58113883008419*Ghat_r[20]*rdv2-1.58113883008419*Ghat_l[20]*rdv2; 
  out[117] += 1.58113883008419*Ghat_r[21]*rdv2-1.58113883008419*Ghat_l[21]*rdv2; 
  out[118] += 1.58113883008419*Ghat_r[22]*rdv2-1.58113883008419*Ghat_l[22]*rdv2; 
  out[119] += 1.58113883008419*Ghat_r[23]*rdv2-1.58113883008419*Ghat_l[23]*rdv2; 
  out[120] += 1.58113883008419*Ghat_r[24]*rdv2-1.58113883008419*Ghat_l[24]*rdv2; 
  out[121] += 1.58113883008419*Ghat_r[25]*rdv2-1.58113883008419*Ghat_l[25]*rdv2; 
  out[122] += 1.58113883008419*Ghat_r[26]*rdv2-1.58113883008419*Ghat_l[26]*rdv2; 
  out[123] += 1.58113883008419*Ghat_r[27]*rdv2-1.58113883008419*Ghat_l[27]*rdv2; 
  out[124] += 1.58113883008419*Ghat_r[28]*rdv2-1.58113883008419*Ghat_l[28]*rdv2; 
  out[125] += 1.58113883008419*Ghat_r[29]*rdv2-1.58113883008419*Ghat_l[29]*rdv2; 
  out[126] += 1.58113883008419*Ghat_r[30]*rdv2-1.58113883008419*Ghat_l[30]*rdv2; 
  out[127] += 1.58113883008419*Ghat_r[31]*rdv2-1.58113883008419*Ghat_l[31]*rdv2; 
  out[128] += 0.7071067811865475*Ghat_r[48]*rdv2-0.7071067811865475*Ghat_l[48]*rdv2; 
  out[129] += 0.7071067811865475*Ghat_r[49]*rdv2-0.7071067811865475*Ghat_l[49]*rdv2; 
  out[130] += 0.7071067811865475*Ghat_r[50]*rdv2-0.7071067811865475*Ghat_l[50]*rdv2; 
  out[131] += 0.7071067811865475*Ghat_r[51]*rdv2-0.7071067811865475*Ghat_l[51]*rdv2; 
  out[132] += 0.7071067811865475*Ghat_r[52]*rdv2-0.7071067811865475*Ghat_l[52]*rdv2; 
  out[133] += 1.224744871391589*Ghat_r[48]*rdv2+1.224744871391589*Ghat_l[48]*rdv2; 
  out[134] += 0.7071067811865475*Ghat_r[53]*rdv2-0.7071067811865475*Ghat_l[53]*rdv2; 
  out[135] += 0.7071067811865475*Ghat_r[54]*rdv2-0.7071067811865475*Ghat_l[54]*rdv2; 
  out[136] += 0.7071067811865475*Ghat_r[55]*rdv2-0.7071067811865475*Ghat_l[55]*rdv2; 
  out[137] += 0.7071067811865475*Ghat_r[56]*rdv2-0.7071067811865475*Ghat_l[56]*rdv2; 
  out[138] += 0.7071067811865475*Ghat_r[57]*rdv2-0.7071067811865475*Ghat_l[57]*rdv2; 
  out[139] += 0.7071067811865475*Ghat_r[58]*rdv2-0.7071067811865475*Ghat_l[58]*rdv2; 
  out[140] += 1.224744871391589*Ghat_r[49]*rdv2+1.224744871391589*Ghat_l[49]*rdv2; 
  out[141] += 1.224744871391589*Ghat_r[50]*rdv2+1.224744871391589*Ghat_l[50]*rdv2; 
  out[142] += 1.224744871391589*Ghat_r[51]*rdv2+1.224744871391589*Ghat_l[51]*rdv2; 
  out[143] += 1.224744871391589*Ghat_r[52]*rdv2+1.224744871391589*Ghat_l[52]*rdv2; 
  out[144] += 0.7071067811865475*Ghat_r[59]*rdv2-0.7071067811865475*Ghat_l[59]*rdv2; 
  out[145] += 0.7071067811865475*Ghat_r[60]*rdv2-0.7071067811865475*Ghat_l[60]*rdv2; 
  out[146] += 0.7071067811865475*Ghat_r[61]*rdv2-0.7071067811865475*Ghat_l[61]*rdv2; 
  out[147] += 0.7071067811865475*Ghat_r[62]*rdv2-0.7071067811865475*Ghat_l[62]*rdv2; 
  out[148] += 1.224744871391589*Ghat_r[53]*rdv2+1.224744871391589*Ghat_l[53]*rdv2; 
  out[149] += 1.224744871391589*Ghat_r[54]*rdv2+1.224744871391589*Ghat_l[54]*rdv2; 
  out[150] += 1.224744871391589*Ghat_r[55]*rdv2+1.224744871391589*Ghat_l[55]*rdv2; 
  out[151] += 1.224744871391589*Ghat_r[56]*rdv2+1.224744871391589*Ghat_l[56]*rdv2; 
  out[152] += 1.224744871391589*Ghat_r[57]*rdv2+1.224744871391589*Ghat_l[57]*rdv2; 
  out[153] += 1.224744871391589*Ghat_r[58]*rdv2+1.224744871391589*Ghat_l[58]*rdv2; 
  out[154] += 0.7071067811865475*Ghat_r[63]*rdv2-0.7071067811865475*Ghat_l[63]*rdv2; 
  out[155] += 1.224744871391589*Ghat_r[59]*rdv2+1.224744871391589*Ghat_l[59]*rdv2; 
  out[156] += 1.224744871391589*Ghat_r[60]*rdv2+1.224744871391589*Ghat_l[60]*rdv2; 
  out[157] += 1.224744871391589*Ghat_r[61]*rdv2+1.224744871391589*Ghat_l[61]*rdv2; 
  out[158] += 1.224744871391589*Ghat_r[62]*rdv2+1.224744871391589*Ghat_l[62]*rdv2; 
  out[159] += 1.224744871391589*Ghat_r[63]*rdv2+1.224744871391589*Ghat_l[63]*rdv2; 

  return 0.;

} 
