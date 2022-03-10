#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_6x_p1_surfx6_eval_quad.h> 
#include <gkyl_basis_ser_6x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvz_3x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[6]:         cell-center coordinates. 
  // dxv[6]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[24]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[8]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[5]; 

  const double *sumNuUz = &nuUSum[16]; 

  double alphaDrSurf_l[32] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*w[5]-1.0*nuSum[0]*dxv[5]-2.0*sumNuUz[0]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*w[5]-1.0*nuSum[1]*dxv[5]-2.0*sumNuUz[1]; 
  alphaDrSurf_l[2] = 2.0*nuSum[2]*w[5]-1.0*nuSum[2]*dxv[5]-2.0*sumNuUz[2]; 
  alphaDrSurf_l[3] = 2.0*nuSum[3]*w[5]-1.0*nuSum[3]*dxv[5]-2.0*sumNuUz[3]; 
  alphaDrSurf_l[6] = 2.0*nuSum[4]*w[5]-1.0*nuSum[4]*dxv[5]-2.0*sumNuUz[4]; 
  alphaDrSurf_l[7] = 2.0*nuSum[5]*w[5]-2.0*sumNuUz[5]-1.0*dxv[5]*nuSum[5]; 
  alphaDrSurf_l[8] = (-2.0*sumNuUz[6])+2.0*w[5]*nuSum[6]-1.0*dxv[5]*nuSum[6]; 
  alphaDrSurf_l[16] = (-2.0*sumNuUz[7])+2.0*w[5]*nuSum[7]-1.0*dxv[5]*nuSum[7]; 

  double alphaDrSurf_r[32] = {0.0}; 
  alphaDrSurf_r[0] = 2.0*nuSum[0]*w[5]+nuSum[0]*dxv[5]-2.0*sumNuUz[0]; 
  alphaDrSurf_r[1] = 2.0*nuSum[1]*w[5]+nuSum[1]*dxv[5]-2.0*sumNuUz[1]; 
  alphaDrSurf_r[2] = 2.0*nuSum[2]*w[5]+nuSum[2]*dxv[5]-2.0*sumNuUz[2]; 
  alphaDrSurf_r[3] = 2.0*nuSum[3]*w[5]+nuSum[3]*dxv[5]-2.0*sumNuUz[3]; 
  alphaDrSurf_r[6] = 2.0*nuSum[4]*w[5]+nuSum[4]*dxv[5]-2.0*sumNuUz[4]; 
  alphaDrSurf_r[7] = 2.0*nuSum[5]*w[5]-2.0*sumNuUz[5]+dxv[5]*nuSum[5]; 
  alphaDrSurf_r[8] = (-2.0*sumNuUz[6])+2.0*w[5]*nuSum[6]+dxv[5]*nuSum[6]; 
  alphaDrSurf_r[16] = (-2.0*sumNuUz[7])+2.0*w[5]*nuSum[7]+dxv[5]*nuSum[7]; 

  double fUpwindQuad_l[32] = {0.0};
  double fUpwindQuad_r[32] = {0.0};
  double fUpwind_l[32] = {0.0};
  double fUpwind_r[32] = {0.0};
  double Ghat_l[32] = {0.0}; 
  double Ghat_r[32] = {0.0}; 

  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_6x_p1_surfx6_eval_quad_node_0_r(fl); 
    fUpwindQuad_l[1] = ser_6x_p1_surfx6_eval_quad_node_1_r(fl); 
    fUpwindQuad_l[2] = ser_6x_p1_surfx6_eval_quad_node_2_r(fl); 
    fUpwindQuad_l[3] = ser_6x_p1_surfx6_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_6x_p1_surfx6_eval_quad_node_0_l(fc); 
    fUpwindQuad_l[1] = ser_6x_p1_surfx6_eval_quad_node_1_l(fc); 
    fUpwindQuad_l[2] = ser_6x_p1_surfx6_eval_quad_node_2_l(fc); 
    fUpwindQuad_l[3] = ser_6x_p1_surfx6_eval_quad_node_3_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_6x_p1_surfx6_eval_quad_node_0_r(fc); 
    fUpwindQuad_r[1] = ser_6x_p1_surfx6_eval_quad_node_1_r(fc); 
    fUpwindQuad_r[2] = ser_6x_p1_surfx6_eval_quad_node_2_r(fc); 
    fUpwindQuad_r[3] = ser_6x_p1_surfx6_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_6x_p1_surfx6_eval_quad_node_0_l(fr); 
    fUpwindQuad_r[1] = ser_6x_p1_surfx6_eval_quad_node_1_l(fr); 
    fUpwindQuad_r[2] = ser_6x_p1_surfx6_eval_quad_node_2_l(fr); 
    fUpwindQuad_r[3] = ser_6x_p1_surfx6_eval_quad_node_3_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[4] = ser_6x_p1_surfx6_eval_quad_node_4_r(fl); 
    fUpwindQuad_l[5] = ser_6x_p1_surfx6_eval_quad_node_5_r(fl); 
    fUpwindQuad_l[6] = ser_6x_p1_surfx6_eval_quad_node_6_r(fl); 
    fUpwindQuad_l[7] = ser_6x_p1_surfx6_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_6x_p1_surfx6_eval_quad_node_4_l(fc); 
    fUpwindQuad_l[5] = ser_6x_p1_surfx6_eval_quad_node_5_l(fc); 
    fUpwindQuad_l[6] = ser_6x_p1_surfx6_eval_quad_node_6_l(fc); 
    fUpwindQuad_l[7] = ser_6x_p1_surfx6_eval_quad_node_7_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[4] = ser_6x_p1_surfx6_eval_quad_node_4_r(fc); 
    fUpwindQuad_r[5] = ser_6x_p1_surfx6_eval_quad_node_5_r(fc); 
    fUpwindQuad_r[6] = ser_6x_p1_surfx6_eval_quad_node_6_r(fc); 
    fUpwindQuad_r[7] = ser_6x_p1_surfx6_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_6x_p1_surfx6_eval_quad_node_4_l(fr); 
    fUpwindQuad_r[5] = ser_6x_p1_surfx6_eval_quad_node_5_l(fr); 
    fUpwindQuad_r[6] = ser_6x_p1_surfx6_eval_quad_node_6_l(fr); 
    fUpwindQuad_r[7] = ser_6x_p1_surfx6_eval_quad_node_7_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = ser_6x_p1_surfx6_eval_quad_node_8_r(fl); 
    fUpwindQuad_l[9] = ser_6x_p1_surfx6_eval_quad_node_9_r(fl); 
    fUpwindQuad_l[10] = ser_6x_p1_surfx6_eval_quad_node_10_r(fl); 
    fUpwindQuad_l[11] = ser_6x_p1_surfx6_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_6x_p1_surfx6_eval_quad_node_8_l(fc); 
    fUpwindQuad_l[9] = ser_6x_p1_surfx6_eval_quad_node_9_l(fc); 
    fUpwindQuad_l[10] = ser_6x_p1_surfx6_eval_quad_node_10_l(fc); 
    fUpwindQuad_l[11] = ser_6x_p1_surfx6_eval_quad_node_11_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = ser_6x_p1_surfx6_eval_quad_node_8_r(fc); 
    fUpwindQuad_r[9] = ser_6x_p1_surfx6_eval_quad_node_9_r(fc); 
    fUpwindQuad_r[10] = ser_6x_p1_surfx6_eval_quad_node_10_r(fc); 
    fUpwindQuad_r[11] = ser_6x_p1_surfx6_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_6x_p1_surfx6_eval_quad_node_8_l(fr); 
    fUpwindQuad_r[9] = ser_6x_p1_surfx6_eval_quad_node_9_l(fr); 
    fUpwindQuad_r[10] = ser_6x_p1_surfx6_eval_quad_node_10_l(fr); 
    fUpwindQuad_r[11] = ser_6x_p1_surfx6_eval_quad_node_11_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[12] = ser_6x_p1_surfx6_eval_quad_node_12_r(fl); 
    fUpwindQuad_l[13] = ser_6x_p1_surfx6_eval_quad_node_13_r(fl); 
    fUpwindQuad_l[14] = ser_6x_p1_surfx6_eval_quad_node_14_r(fl); 
    fUpwindQuad_l[15] = ser_6x_p1_surfx6_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_l[12] = ser_6x_p1_surfx6_eval_quad_node_12_l(fc); 
    fUpwindQuad_l[13] = ser_6x_p1_surfx6_eval_quad_node_13_l(fc); 
    fUpwindQuad_l[14] = ser_6x_p1_surfx6_eval_quad_node_14_l(fc); 
    fUpwindQuad_l[15] = ser_6x_p1_surfx6_eval_quad_node_15_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[12] = ser_6x_p1_surfx6_eval_quad_node_12_r(fc); 
    fUpwindQuad_r[13] = ser_6x_p1_surfx6_eval_quad_node_13_r(fc); 
    fUpwindQuad_r[14] = ser_6x_p1_surfx6_eval_quad_node_14_r(fc); 
    fUpwindQuad_r[15] = ser_6x_p1_surfx6_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_r[12] = ser_6x_p1_surfx6_eval_quad_node_12_l(fr); 
    fUpwindQuad_r[13] = ser_6x_p1_surfx6_eval_quad_node_13_l(fr); 
    fUpwindQuad_r[14] = ser_6x_p1_surfx6_eval_quad_node_14_l(fr); 
    fUpwindQuad_r[15] = ser_6x_p1_surfx6_eval_quad_node_15_l(fr); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[16] = ser_6x_p1_surfx6_eval_quad_node_16_r(fl); 
    fUpwindQuad_l[17] = ser_6x_p1_surfx6_eval_quad_node_17_r(fl); 
    fUpwindQuad_l[18] = ser_6x_p1_surfx6_eval_quad_node_18_r(fl); 
    fUpwindQuad_l[19] = ser_6x_p1_surfx6_eval_quad_node_19_r(fl); 
  } else { 
    fUpwindQuad_l[16] = ser_6x_p1_surfx6_eval_quad_node_16_l(fc); 
    fUpwindQuad_l[17] = ser_6x_p1_surfx6_eval_quad_node_17_l(fc); 
    fUpwindQuad_l[18] = ser_6x_p1_surfx6_eval_quad_node_18_l(fc); 
    fUpwindQuad_l[19] = ser_6x_p1_surfx6_eval_quad_node_19_l(fc); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[16] = ser_6x_p1_surfx6_eval_quad_node_16_r(fc); 
    fUpwindQuad_r[17] = ser_6x_p1_surfx6_eval_quad_node_17_r(fc); 
    fUpwindQuad_r[18] = ser_6x_p1_surfx6_eval_quad_node_18_r(fc); 
    fUpwindQuad_r[19] = ser_6x_p1_surfx6_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_r[16] = ser_6x_p1_surfx6_eval_quad_node_16_l(fr); 
    fUpwindQuad_r[17] = ser_6x_p1_surfx6_eval_quad_node_17_l(fr); 
    fUpwindQuad_r[18] = ser_6x_p1_surfx6_eval_quad_node_18_l(fr); 
    fUpwindQuad_r[19] = ser_6x_p1_surfx6_eval_quad_node_19_l(fr); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[20] = ser_6x_p1_surfx6_eval_quad_node_20_r(fl); 
    fUpwindQuad_l[21] = ser_6x_p1_surfx6_eval_quad_node_21_r(fl); 
    fUpwindQuad_l[22] = ser_6x_p1_surfx6_eval_quad_node_22_r(fl); 
    fUpwindQuad_l[23] = ser_6x_p1_surfx6_eval_quad_node_23_r(fl); 
  } else { 
    fUpwindQuad_l[20] = ser_6x_p1_surfx6_eval_quad_node_20_l(fc); 
    fUpwindQuad_l[21] = ser_6x_p1_surfx6_eval_quad_node_21_l(fc); 
    fUpwindQuad_l[22] = ser_6x_p1_surfx6_eval_quad_node_22_l(fc); 
    fUpwindQuad_l[23] = ser_6x_p1_surfx6_eval_quad_node_23_l(fc); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[20] = ser_6x_p1_surfx6_eval_quad_node_20_r(fc); 
    fUpwindQuad_r[21] = ser_6x_p1_surfx6_eval_quad_node_21_r(fc); 
    fUpwindQuad_r[22] = ser_6x_p1_surfx6_eval_quad_node_22_r(fc); 
    fUpwindQuad_r[23] = ser_6x_p1_surfx6_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_r[20] = ser_6x_p1_surfx6_eval_quad_node_20_l(fr); 
    fUpwindQuad_r[21] = ser_6x_p1_surfx6_eval_quad_node_21_l(fr); 
    fUpwindQuad_r[22] = ser_6x_p1_surfx6_eval_quad_node_22_l(fr); 
    fUpwindQuad_r[23] = ser_6x_p1_surfx6_eval_quad_node_23_l(fr); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[24] = ser_6x_p1_surfx6_eval_quad_node_24_r(fl); 
    fUpwindQuad_l[25] = ser_6x_p1_surfx6_eval_quad_node_25_r(fl); 
    fUpwindQuad_l[26] = ser_6x_p1_surfx6_eval_quad_node_26_r(fl); 
    fUpwindQuad_l[27] = ser_6x_p1_surfx6_eval_quad_node_27_r(fl); 
  } else { 
    fUpwindQuad_l[24] = ser_6x_p1_surfx6_eval_quad_node_24_l(fc); 
    fUpwindQuad_l[25] = ser_6x_p1_surfx6_eval_quad_node_25_l(fc); 
    fUpwindQuad_l[26] = ser_6x_p1_surfx6_eval_quad_node_26_l(fc); 
    fUpwindQuad_l[27] = ser_6x_p1_surfx6_eval_quad_node_27_l(fc); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[24] = ser_6x_p1_surfx6_eval_quad_node_24_r(fc); 
    fUpwindQuad_r[25] = ser_6x_p1_surfx6_eval_quad_node_25_r(fc); 
    fUpwindQuad_r[26] = ser_6x_p1_surfx6_eval_quad_node_26_r(fc); 
    fUpwindQuad_r[27] = ser_6x_p1_surfx6_eval_quad_node_27_r(fc); 
  } else { 
    fUpwindQuad_r[24] = ser_6x_p1_surfx6_eval_quad_node_24_l(fr); 
    fUpwindQuad_r[25] = ser_6x_p1_surfx6_eval_quad_node_25_l(fr); 
    fUpwindQuad_r[26] = ser_6x_p1_surfx6_eval_quad_node_26_l(fr); 
    fUpwindQuad_r[27] = ser_6x_p1_surfx6_eval_quad_node_27_l(fr); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[28] = ser_6x_p1_surfx6_eval_quad_node_28_r(fl); 
    fUpwindQuad_l[29] = ser_6x_p1_surfx6_eval_quad_node_29_r(fl); 
    fUpwindQuad_l[30] = ser_6x_p1_surfx6_eval_quad_node_30_r(fl); 
    fUpwindQuad_l[31] = ser_6x_p1_surfx6_eval_quad_node_31_r(fl); 
  } else { 
    fUpwindQuad_l[28] = ser_6x_p1_surfx6_eval_quad_node_28_l(fc); 
    fUpwindQuad_l[29] = ser_6x_p1_surfx6_eval_quad_node_29_l(fc); 
    fUpwindQuad_l[30] = ser_6x_p1_surfx6_eval_quad_node_30_l(fc); 
    fUpwindQuad_l[31] = ser_6x_p1_surfx6_eval_quad_node_31_l(fc); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[28] = ser_6x_p1_surfx6_eval_quad_node_28_r(fc); 
    fUpwindQuad_r[29] = ser_6x_p1_surfx6_eval_quad_node_29_r(fc); 
    fUpwindQuad_r[30] = ser_6x_p1_surfx6_eval_quad_node_30_r(fc); 
    fUpwindQuad_r[31] = ser_6x_p1_surfx6_eval_quad_node_31_r(fc); 
  } else { 
    fUpwindQuad_r[28] = ser_6x_p1_surfx6_eval_quad_node_28_l(fr); 
    fUpwindQuad_r[29] = ser_6x_p1_surfx6_eval_quad_node_29_l(fr); 
    fUpwindQuad_r[30] = ser_6x_p1_surfx6_eval_quad_node_30_l(fr); 
    fUpwindQuad_r[31] = ser_6x_p1_surfx6_eval_quad_node_31_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

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

  out[0] += 0.7071067811865475*Ghat_r[0]*rdv2-0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_r[1]*rdv2-0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_r[2]*rdv2-0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_r[3]*rdv2-0.7071067811865475*Ghat_l[3]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_r[4]*rdv2-0.7071067811865475*Ghat_l[4]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_r[5]*rdv2-0.7071067811865475*Ghat_l[5]*rdv2; 
  out[6] += 1.224744871391589*Ghat_r[0]*rdv2+1.224744871391589*Ghat_l[0]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_r[6]*rdv2-0.7071067811865475*Ghat_l[6]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_r[7]*rdv2-0.7071067811865475*Ghat_l[7]*rdv2; 
  out[9] += 0.7071067811865475*Ghat_r[8]*rdv2-0.7071067811865475*Ghat_l[8]*rdv2; 
  out[10] += 0.7071067811865475*Ghat_r[9]*rdv2-0.7071067811865475*Ghat_l[9]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_r[10]*rdv2-0.7071067811865475*Ghat_l[10]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_r[11]*rdv2-0.7071067811865475*Ghat_l[11]*rdv2; 
  out[13] += 0.7071067811865475*Ghat_r[12]*rdv2-0.7071067811865475*Ghat_l[12]*rdv2; 
  out[14] += 0.7071067811865475*Ghat_r[13]*rdv2-0.7071067811865475*Ghat_l[13]*rdv2; 
  out[15] += 0.7071067811865475*Ghat_r[14]*rdv2-0.7071067811865475*Ghat_l[14]*rdv2; 
  out[16] += 0.7071067811865475*Ghat_r[15]*rdv2-0.7071067811865475*Ghat_l[15]*rdv2; 
  out[17] += 1.224744871391589*Ghat_r[1]*rdv2+1.224744871391589*Ghat_l[1]*rdv2; 
  out[18] += 1.224744871391589*Ghat_r[2]*rdv2+1.224744871391589*Ghat_l[2]*rdv2; 
  out[19] += 1.224744871391589*Ghat_r[3]*rdv2+1.224744871391589*Ghat_l[3]*rdv2; 
  out[20] += 1.224744871391589*Ghat_r[4]*rdv2+1.224744871391589*Ghat_l[4]*rdv2; 
  out[21] += 1.224744871391589*Ghat_r[5]*rdv2+1.224744871391589*Ghat_l[5]*rdv2; 
  out[22] += 0.7071067811865475*Ghat_r[16]*rdv2-0.7071067811865475*Ghat_l[16]*rdv2; 
  out[23] += 0.7071067811865475*Ghat_r[17]*rdv2-0.7071067811865475*Ghat_l[17]*rdv2; 
  out[24] += 0.7071067811865475*Ghat_r[18]*rdv2-0.7071067811865475*Ghat_l[18]*rdv2; 
  out[25] += 0.7071067811865475*Ghat_r[19]*rdv2-0.7071067811865475*Ghat_l[19]*rdv2; 
  out[26] += 0.7071067811865475*Ghat_r[20]*rdv2-0.7071067811865475*Ghat_l[20]*rdv2; 
  out[27] += 0.7071067811865475*Ghat_r[21]*rdv2-0.7071067811865475*Ghat_l[21]*rdv2; 
  out[28] += 0.7071067811865475*Ghat_r[22]*rdv2-0.7071067811865475*Ghat_l[22]*rdv2; 
  out[29] += 0.7071067811865475*Ghat_r[23]*rdv2-0.7071067811865475*Ghat_l[23]*rdv2; 
  out[30] += 0.7071067811865475*Ghat_r[24]*rdv2-0.7071067811865475*Ghat_l[24]*rdv2; 
  out[31] += 0.7071067811865475*Ghat_r[25]*rdv2-0.7071067811865475*Ghat_l[25]*rdv2; 
  out[32] += 1.224744871391589*Ghat_r[6]*rdv2+1.224744871391589*Ghat_l[6]*rdv2; 
  out[33] += 1.224744871391589*Ghat_r[7]*rdv2+1.224744871391589*Ghat_l[7]*rdv2; 
  out[34] += 1.224744871391589*Ghat_r[8]*rdv2+1.224744871391589*Ghat_l[8]*rdv2; 
  out[35] += 1.224744871391589*Ghat_r[9]*rdv2+1.224744871391589*Ghat_l[9]*rdv2; 
  out[36] += 1.224744871391589*Ghat_r[10]*rdv2+1.224744871391589*Ghat_l[10]*rdv2; 
  out[37] += 1.224744871391589*Ghat_r[11]*rdv2+1.224744871391589*Ghat_l[11]*rdv2; 
  out[38] += 1.224744871391589*Ghat_r[12]*rdv2+1.224744871391589*Ghat_l[12]*rdv2; 
  out[39] += 1.224744871391589*Ghat_r[13]*rdv2+1.224744871391589*Ghat_l[13]*rdv2; 
  out[40] += 1.224744871391589*Ghat_r[14]*rdv2+1.224744871391589*Ghat_l[14]*rdv2; 
  out[41] += 1.224744871391589*Ghat_r[15]*rdv2+1.224744871391589*Ghat_l[15]*rdv2; 
  out[42] += 0.7071067811865475*Ghat_r[26]*rdv2-0.7071067811865475*Ghat_l[26]*rdv2; 
  out[43] += 0.7071067811865475*Ghat_r[27]*rdv2-0.7071067811865475*Ghat_l[27]*rdv2; 
  out[44] += 0.7071067811865475*Ghat_r[28]*rdv2-0.7071067811865475*Ghat_l[28]*rdv2; 
  out[45] += 0.7071067811865475*Ghat_r[29]*rdv2-0.7071067811865475*Ghat_l[29]*rdv2; 
  out[46] += 0.7071067811865475*Ghat_r[30]*rdv2-0.7071067811865475*Ghat_l[30]*rdv2; 
  out[47] += 1.224744871391589*Ghat_r[16]*rdv2+1.224744871391589*Ghat_l[16]*rdv2; 
  out[48] += 1.224744871391589*Ghat_r[17]*rdv2+1.224744871391589*Ghat_l[17]*rdv2; 
  out[49] += 1.224744871391589*Ghat_r[18]*rdv2+1.224744871391589*Ghat_l[18]*rdv2; 
  out[50] += 1.224744871391589*Ghat_r[19]*rdv2+1.224744871391589*Ghat_l[19]*rdv2; 
  out[51] += 1.224744871391589*Ghat_r[20]*rdv2+1.224744871391589*Ghat_l[20]*rdv2; 
  out[52] += 1.224744871391589*Ghat_r[21]*rdv2+1.224744871391589*Ghat_l[21]*rdv2; 
  out[53] += 1.224744871391589*Ghat_r[22]*rdv2+1.224744871391589*Ghat_l[22]*rdv2; 
  out[54] += 1.224744871391589*Ghat_r[23]*rdv2+1.224744871391589*Ghat_l[23]*rdv2; 
  out[55] += 1.224744871391589*Ghat_r[24]*rdv2+1.224744871391589*Ghat_l[24]*rdv2; 
  out[56] += 1.224744871391589*Ghat_r[25]*rdv2+1.224744871391589*Ghat_l[25]*rdv2; 
  out[57] += 0.7071067811865475*Ghat_r[31]*rdv2-0.7071067811865475*Ghat_l[31]*rdv2; 
  out[58] += 1.224744871391589*Ghat_r[26]*rdv2+1.224744871391589*Ghat_l[26]*rdv2; 
  out[59] += 1.224744871391589*Ghat_r[27]*rdv2+1.224744871391589*Ghat_l[27]*rdv2; 
  out[60] += 1.224744871391589*Ghat_r[28]*rdv2+1.224744871391589*Ghat_l[28]*rdv2; 
  out[61] += 1.224744871391589*Ghat_r[29]*rdv2+1.224744871391589*Ghat_l[29]*rdv2; 
  out[62] += 1.224744871391589*Ghat_r[30]*rdv2+1.224744871391589*Ghat_l[30]*rdv2; 
  out[63] += 1.224744871391589*Ghat_r[31]*rdv2+1.224744871391589*Ghat_l[31]*rdv2; 
} 
