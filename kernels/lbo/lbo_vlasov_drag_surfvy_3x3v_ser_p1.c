#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_6x_p1_surfx5_quad.h> 
#include <gkyl_basis_ser_6x_p1_upwind.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[6]:         cell-center coordinates. 
  // dxv[6]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[24]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[8]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[4]; 

  const double *sumNuUy = &nuUSum[8]; 

  double alphaDrSurf_l[32] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*w[4]-1.0*nuSum[0]*dxv[4]-2.0*sumNuUy[0]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*w[4]-1.0*nuSum[1]*dxv[4]-2.0*sumNuUy[1]; 
  alphaDrSurf_l[2] = 2.0*nuSum[2]*w[4]-1.0*nuSum[2]*dxv[4]-2.0*sumNuUy[2]; 
  alphaDrSurf_l[3] = 2.0*nuSum[3]*w[4]-1.0*nuSum[3]*dxv[4]-2.0*sumNuUy[3]; 
  alphaDrSurf_l[6] = 2.0*nuSum[4]*w[4]-2.0*sumNuUy[4]-1.0*dxv[4]*nuSum[4]; 
  alphaDrSurf_l[7] = (-2.0*sumNuUy[5])+2.0*w[4]*nuSum[5]-1.0*dxv[4]*nuSum[5]; 
  alphaDrSurf_l[8] = (-2.0*sumNuUy[6])+2.0*w[4]*nuSum[6]-1.0*dxv[4]*nuSum[6]; 
  alphaDrSurf_l[16] = (-2.0*sumNuUy[7])+2.0*w[4]*nuSum[7]-1.0*dxv[4]*nuSum[7]; 

  double alphaDrSurf_r[32] = {0.0}; 
  alphaDrSurf_r[0] = 2.0*nuSum[0]*w[4]+nuSum[0]*dxv[4]-2.0*sumNuUy[0]; 
  alphaDrSurf_r[1] = 2.0*nuSum[1]*w[4]+nuSum[1]*dxv[4]-2.0*sumNuUy[1]; 
  alphaDrSurf_r[2] = 2.0*nuSum[2]*w[4]+nuSum[2]*dxv[4]-2.0*sumNuUy[2]; 
  alphaDrSurf_r[3] = 2.0*nuSum[3]*w[4]+nuSum[3]*dxv[4]-2.0*sumNuUy[3]; 
  alphaDrSurf_r[6] = 2.0*nuSum[4]*w[4]-2.0*sumNuUy[4]+dxv[4]*nuSum[4]; 
  alphaDrSurf_r[7] = (-2.0*sumNuUy[5])+2.0*w[4]*nuSum[5]+dxv[4]*nuSum[5]; 
  alphaDrSurf_r[8] = (-2.0*sumNuUy[6])+2.0*w[4]*nuSum[6]+dxv[4]*nuSum[6]; 
  alphaDrSurf_r[16] = (-2.0*sumNuUy[7])+2.0*w[4]*nuSum[7]+dxv[4]*nuSum[7]; 

  double fUpwindQuad_l[32] = {0.0};
  double fUpwindQuad_r[32] = {0.0};
  double fUpwind_l[32] = {0.0};
  double fUpwind_r[32] = {0.0};
  double Gdrag_l[32] = {0.0}; 
  double Gdrag_r[32] = {0.0}; 

  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_6x_p1_surfx5_quad_0_r(fl); 
    fUpwindQuad_l[8] = ser_6x_p1_surfx5_quad_8_r(fl); 
    fUpwindQuad_l[16] = ser_6x_p1_surfx5_quad_16_r(fl); 
    fUpwindQuad_l[24] = ser_6x_p1_surfx5_quad_24_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_6x_p1_surfx5_quad_0_l(fc); 
    fUpwindQuad_l[8] = ser_6x_p1_surfx5_quad_8_l(fc); 
    fUpwindQuad_l[16] = ser_6x_p1_surfx5_quad_16_l(fc); 
    fUpwindQuad_l[24] = ser_6x_p1_surfx5_quad_24_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_6x_p1_surfx5_quad_0_r(fc); 
    fUpwindQuad_r[8] = ser_6x_p1_surfx5_quad_8_r(fc); 
    fUpwindQuad_r[16] = ser_6x_p1_surfx5_quad_16_r(fc); 
    fUpwindQuad_r[24] = ser_6x_p1_surfx5_quad_24_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_6x_p1_surfx5_quad_0_l(fr); 
    fUpwindQuad_r[8] = ser_6x_p1_surfx5_quad_8_l(fr); 
    fUpwindQuad_r[16] = ser_6x_p1_surfx5_quad_16_l(fr); 
    fUpwindQuad_r[24] = ser_6x_p1_surfx5_quad_24_l(fr); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_6x_p1_surfx5_quad_1_r(fl); 
    fUpwindQuad_l[9] = ser_6x_p1_surfx5_quad_9_r(fl); 
    fUpwindQuad_l[17] = ser_6x_p1_surfx5_quad_17_r(fl); 
    fUpwindQuad_l[25] = ser_6x_p1_surfx5_quad_25_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_6x_p1_surfx5_quad_1_l(fc); 
    fUpwindQuad_l[9] = ser_6x_p1_surfx5_quad_9_l(fc); 
    fUpwindQuad_l[17] = ser_6x_p1_surfx5_quad_17_l(fc); 
    fUpwindQuad_l[25] = ser_6x_p1_surfx5_quad_25_l(fc); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_6x_p1_surfx5_quad_1_r(fc); 
    fUpwindQuad_r[9] = ser_6x_p1_surfx5_quad_9_r(fc); 
    fUpwindQuad_r[17] = ser_6x_p1_surfx5_quad_17_r(fc); 
    fUpwindQuad_r[25] = ser_6x_p1_surfx5_quad_25_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_6x_p1_surfx5_quad_1_l(fr); 
    fUpwindQuad_r[9] = ser_6x_p1_surfx5_quad_9_l(fr); 
    fUpwindQuad_r[17] = ser_6x_p1_surfx5_quad_17_l(fr); 
    fUpwindQuad_r[25] = ser_6x_p1_surfx5_quad_25_l(fr); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_6x_p1_surfx5_quad_2_r(fl); 
    fUpwindQuad_l[10] = ser_6x_p1_surfx5_quad_10_r(fl); 
    fUpwindQuad_l[18] = ser_6x_p1_surfx5_quad_18_r(fl); 
    fUpwindQuad_l[26] = ser_6x_p1_surfx5_quad_26_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_6x_p1_surfx5_quad_2_l(fc); 
    fUpwindQuad_l[10] = ser_6x_p1_surfx5_quad_10_l(fc); 
    fUpwindQuad_l[18] = ser_6x_p1_surfx5_quad_18_l(fc); 
    fUpwindQuad_l[26] = ser_6x_p1_surfx5_quad_26_l(fc); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_6x_p1_surfx5_quad_2_r(fc); 
    fUpwindQuad_r[10] = ser_6x_p1_surfx5_quad_10_r(fc); 
    fUpwindQuad_r[18] = ser_6x_p1_surfx5_quad_18_r(fc); 
    fUpwindQuad_r[26] = ser_6x_p1_surfx5_quad_26_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_6x_p1_surfx5_quad_2_l(fr); 
    fUpwindQuad_r[10] = ser_6x_p1_surfx5_quad_10_l(fr); 
    fUpwindQuad_r[18] = ser_6x_p1_surfx5_quad_18_l(fr); 
    fUpwindQuad_r[26] = ser_6x_p1_surfx5_quad_26_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]-alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_6x_p1_surfx5_quad_3_r(fl); 
    fUpwindQuad_l[11] = ser_6x_p1_surfx5_quad_11_r(fl); 
    fUpwindQuad_l[19] = ser_6x_p1_surfx5_quad_19_r(fl); 
    fUpwindQuad_l[27] = ser_6x_p1_surfx5_quad_27_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_6x_p1_surfx5_quad_3_l(fc); 
    fUpwindQuad_l[11] = ser_6x_p1_surfx5_quad_11_l(fc); 
    fUpwindQuad_l[19] = ser_6x_p1_surfx5_quad_19_l(fc); 
    fUpwindQuad_l[27] = ser_6x_p1_surfx5_quad_27_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]-alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_6x_p1_surfx5_quad_3_r(fc); 
    fUpwindQuad_r[11] = ser_6x_p1_surfx5_quad_11_r(fc); 
    fUpwindQuad_r[19] = ser_6x_p1_surfx5_quad_19_r(fc); 
    fUpwindQuad_r[27] = ser_6x_p1_surfx5_quad_27_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_6x_p1_surfx5_quad_3_l(fr); 
    fUpwindQuad_r[11] = ser_6x_p1_surfx5_quad_11_l(fr); 
    fUpwindQuad_r[19] = ser_6x_p1_surfx5_quad_19_l(fr); 
    fUpwindQuad_r[27] = ser_6x_p1_surfx5_quad_27_l(fr); 
  } 
  if (alphaDrSurf_l[16]-alphaDrSurf_l[8]-alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[4] = ser_6x_p1_surfx5_quad_4_r(fl); 
    fUpwindQuad_l[12] = ser_6x_p1_surfx5_quad_12_r(fl); 
    fUpwindQuad_l[20] = ser_6x_p1_surfx5_quad_20_r(fl); 
    fUpwindQuad_l[28] = ser_6x_p1_surfx5_quad_28_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_6x_p1_surfx5_quad_4_l(fc); 
    fUpwindQuad_l[12] = ser_6x_p1_surfx5_quad_12_l(fc); 
    fUpwindQuad_l[20] = ser_6x_p1_surfx5_quad_20_l(fc); 
    fUpwindQuad_l[28] = ser_6x_p1_surfx5_quad_28_l(fc); 
  } 
  if (alphaDrSurf_r[16]-alphaDrSurf_r[8]-alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[4] = ser_6x_p1_surfx5_quad_4_r(fc); 
    fUpwindQuad_r[12] = ser_6x_p1_surfx5_quad_12_r(fc); 
    fUpwindQuad_r[20] = ser_6x_p1_surfx5_quad_20_r(fc); 
    fUpwindQuad_r[28] = ser_6x_p1_surfx5_quad_28_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_6x_p1_surfx5_quad_4_l(fr); 
    fUpwindQuad_r[12] = ser_6x_p1_surfx5_quad_12_l(fr); 
    fUpwindQuad_r[20] = ser_6x_p1_surfx5_quad_20_l(fr); 
    fUpwindQuad_r[28] = ser_6x_p1_surfx5_quad_28_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])-alphaDrSurf_l[8]+alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_6x_p1_surfx5_quad_5_r(fl); 
    fUpwindQuad_l[13] = ser_6x_p1_surfx5_quad_13_r(fl); 
    fUpwindQuad_l[21] = ser_6x_p1_surfx5_quad_21_r(fl); 
    fUpwindQuad_l[29] = ser_6x_p1_surfx5_quad_29_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_6x_p1_surfx5_quad_5_l(fc); 
    fUpwindQuad_l[13] = ser_6x_p1_surfx5_quad_13_l(fc); 
    fUpwindQuad_l[21] = ser_6x_p1_surfx5_quad_21_l(fc); 
    fUpwindQuad_l[29] = ser_6x_p1_surfx5_quad_29_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])-alphaDrSurf_r[8]+alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_6x_p1_surfx5_quad_5_r(fc); 
    fUpwindQuad_r[13] = ser_6x_p1_surfx5_quad_13_r(fc); 
    fUpwindQuad_r[21] = ser_6x_p1_surfx5_quad_21_r(fc); 
    fUpwindQuad_r[29] = ser_6x_p1_surfx5_quad_29_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_6x_p1_surfx5_quad_5_l(fr); 
    fUpwindQuad_r[13] = ser_6x_p1_surfx5_quad_13_l(fr); 
    fUpwindQuad_r[21] = ser_6x_p1_surfx5_quad_21_l(fr); 
    fUpwindQuad_r[29] = ser_6x_p1_surfx5_quad_29_l(fr); 
  } 
  if ((-alphaDrSurf_l[16])+alphaDrSurf_l[8]-alphaDrSurf_l[7]-alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_6x_p1_surfx5_quad_6_r(fl); 
    fUpwindQuad_l[14] = ser_6x_p1_surfx5_quad_14_r(fl); 
    fUpwindQuad_l[22] = ser_6x_p1_surfx5_quad_22_r(fl); 
    fUpwindQuad_l[30] = ser_6x_p1_surfx5_quad_30_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_6x_p1_surfx5_quad_6_l(fc); 
    fUpwindQuad_l[14] = ser_6x_p1_surfx5_quad_14_l(fc); 
    fUpwindQuad_l[22] = ser_6x_p1_surfx5_quad_22_l(fc); 
    fUpwindQuad_l[30] = ser_6x_p1_surfx5_quad_30_l(fc); 
  } 
  if ((-alphaDrSurf_r[16])+alphaDrSurf_r[8]-alphaDrSurf_r[7]-alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_6x_p1_surfx5_quad_6_r(fc); 
    fUpwindQuad_r[14] = ser_6x_p1_surfx5_quad_14_r(fc); 
    fUpwindQuad_r[22] = ser_6x_p1_surfx5_quad_22_r(fc); 
    fUpwindQuad_r[30] = ser_6x_p1_surfx5_quad_30_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_6x_p1_surfx5_quad_6_l(fr); 
    fUpwindQuad_r[14] = ser_6x_p1_surfx5_quad_14_l(fr); 
    fUpwindQuad_r[22] = ser_6x_p1_surfx5_quad_22_l(fr); 
    fUpwindQuad_r[30] = ser_6x_p1_surfx5_quad_30_l(fr); 
  } 
  if (alphaDrSurf_l[16]+alphaDrSurf_l[8]+alphaDrSurf_l[7]+alphaDrSurf_l[6]+alphaDrSurf_l[3]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[7] = ser_6x_p1_surfx5_quad_7_r(fl); 
    fUpwindQuad_l[15] = ser_6x_p1_surfx5_quad_15_r(fl); 
    fUpwindQuad_l[23] = ser_6x_p1_surfx5_quad_23_r(fl); 
    fUpwindQuad_l[31] = ser_6x_p1_surfx5_quad_31_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_6x_p1_surfx5_quad_7_l(fc); 
    fUpwindQuad_l[15] = ser_6x_p1_surfx5_quad_15_l(fc); 
    fUpwindQuad_l[23] = ser_6x_p1_surfx5_quad_23_l(fc); 
    fUpwindQuad_l[31] = ser_6x_p1_surfx5_quad_31_l(fc); 
  } 
  if (alphaDrSurf_r[16]+alphaDrSurf_r[8]+alphaDrSurf_r[7]+alphaDrSurf_r[6]+alphaDrSurf_r[3]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[7] = ser_6x_p1_surfx5_quad_7_r(fc); 
    fUpwindQuad_r[15] = ser_6x_p1_surfx5_quad_15_r(fc); 
    fUpwindQuad_r[23] = ser_6x_p1_surfx5_quad_23_r(fc); 
    fUpwindQuad_r[31] = ser_6x_p1_surfx5_quad_31_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_6x_p1_surfx5_quad_7_l(fr); 
    fUpwindQuad_r[15] = ser_6x_p1_surfx5_quad_15_l(fr); 
    fUpwindQuad_r[23] = ser_6x_p1_surfx5_quad_23_l(fr); 
    fUpwindQuad_r[31] = ser_6x_p1_surfx5_quad_31_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_6x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_6x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Gdrag_l[0] = 0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[16]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[8]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[7]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[6]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[3]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[2]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[1]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[8]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[3]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[6]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[1]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[7]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[3]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[6]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[2]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Gdrag_l[3] = 0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[6]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[3]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[3]; 
  Gdrag_l[4] = 0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[17]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[9]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[4]; 
  Gdrag_l[5] = 0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[20]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[12]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[5]; 
  Gdrag_l[6] = 0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[3]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[7]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[6]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[2]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Gdrag_l[7] = 0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[6]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[3]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[3]; 
  Gdrag_l[8] = 0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[6]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[3]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[3]; 
  Gdrag_l[9] = 0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[17]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[9]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[4]; 
  Gdrag_l[10] = 0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[17]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[9]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[4]; 
  Gdrag_l[11] = 0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[17]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[9]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[4]; 
  Gdrag_l[12] = 0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[20]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[12]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[5]; 
  Gdrag_l[13] = 0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[20]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[12]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[5]; 
  Gdrag_l[14] = 0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[20]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[12]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[5]; 
  Gdrag_l[15] = 0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[15]; 
  Gdrag_l[16] = 0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[16]+0.1767766952966368*fUpwind_l[0]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[8]+0.1767766952966368*fUpwind_l[1]*alphaDrSurf_l[8]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[7]+0.1767766952966368*fUpwind_l[2]*alphaDrSurf_l[7]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[6]+0.1767766952966368*fUpwind_l[3]*alphaDrSurf_l[6]; 
  Gdrag_l[17] = 0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[11]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[4]*alphaDrSurf_l[6]; 
  Gdrag_l[18] = 0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[10]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[9]+0.1767766952966368*fUpwind_l[4]*alphaDrSurf_l[7]; 
  Gdrag_l[19] = 0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[9]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[10]+0.1767766952966368*fUpwind_l[4]*alphaDrSurf_l[8]; 
  Gdrag_l[20] = 0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[14]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[5]*alphaDrSurf_l[6]; 
  Gdrag_l[21] = 0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[13]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[12]+0.1767766952966368*fUpwind_l[5]*alphaDrSurf_l[7]; 
  Gdrag_l[22] = 0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[12]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[13]+0.1767766952966368*fUpwind_l[5]*alphaDrSurf_l[8]; 
  Gdrag_l[23] = 0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[15]; 
  Gdrag_l[24] = 0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[15]; 
  Gdrag_l[25] = 0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[15]; 
  Gdrag_l[26] = 0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[26]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[19]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[18]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[17]+0.1767766952966368*fUpwind_l[4]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[11]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[10]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[9]; 
  Gdrag_l[27] = 0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[27]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[22]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[21]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[20]+0.1767766952966368*fUpwind_l[5]*alphaDrSurf_l[16]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[14]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[13]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[12]; 
  Gdrag_l[28] = 0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[15]; 
  Gdrag_l[29] = 0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[15]; 
  Gdrag_l[30] = 0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[16]*fUpwind_l[23]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[15]; 
  Gdrag_l[31] = 0.1767766952966368*alphaDrSurf_l[0]*fUpwind_l[31]+0.1767766952966368*alphaDrSurf_l[1]*fUpwind_l[30]+0.1767766952966368*alphaDrSurf_l[2]*fUpwind_l[29]+0.1767766952966368*alphaDrSurf_l[3]*fUpwind_l[28]+0.1767766952966368*alphaDrSurf_l[6]*fUpwind_l[25]+0.1767766952966368*alphaDrSurf_l[7]*fUpwind_l[24]+0.1767766952966368*alphaDrSurf_l[8]*fUpwind_l[23]+0.1767766952966368*fUpwind_l[15]*alphaDrSurf_l[16]; 

  Gdrag_r[0] = 0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[16]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[8]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[7]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[6]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[3]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[2]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[1]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[8]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[3]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[6]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[1]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[7]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[3]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[6]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[2]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Gdrag_r[3] = 0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[6]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[3]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[3]; 
  Gdrag_r[4] = 0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[17]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[9]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[4]; 
  Gdrag_r[5] = 0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[20]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[12]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[5]; 
  Gdrag_r[6] = 0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[3]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[7]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[6]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[2]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Gdrag_r[7] = 0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[6]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[3]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[3]; 
  Gdrag_r[8] = 0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[6]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[3]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[3]; 
  Gdrag_r[9] = 0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[17]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[9]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[4]; 
  Gdrag_r[10] = 0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[17]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[9]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[4]; 
  Gdrag_r[11] = 0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[17]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[9]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[4]; 
  Gdrag_r[12] = 0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[20]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[12]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[5]; 
  Gdrag_r[13] = 0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[20]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[12]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[5]; 
  Gdrag_r[14] = 0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[20]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[12]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[5]; 
  Gdrag_r[15] = 0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[15]; 
  Gdrag_r[16] = 0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[16]+0.1767766952966368*fUpwind_r[0]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[8]+0.1767766952966368*fUpwind_r[1]*alphaDrSurf_r[8]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[7]+0.1767766952966368*fUpwind_r[2]*alphaDrSurf_r[7]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[6]+0.1767766952966368*fUpwind_r[3]*alphaDrSurf_r[6]; 
  Gdrag_r[17] = 0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[11]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[4]*alphaDrSurf_r[6]; 
  Gdrag_r[18] = 0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[10]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[9]+0.1767766952966368*fUpwind_r[4]*alphaDrSurf_r[7]; 
  Gdrag_r[19] = 0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[9]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[10]+0.1767766952966368*fUpwind_r[4]*alphaDrSurf_r[8]; 
  Gdrag_r[20] = 0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[14]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[5]*alphaDrSurf_r[6]; 
  Gdrag_r[21] = 0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[13]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[12]+0.1767766952966368*fUpwind_r[5]*alphaDrSurf_r[7]; 
  Gdrag_r[22] = 0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[12]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[13]+0.1767766952966368*fUpwind_r[5]*alphaDrSurf_r[8]; 
  Gdrag_r[23] = 0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[15]; 
  Gdrag_r[24] = 0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[15]; 
  Gdrag_r[25] = 0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[15]; 
  Gdrag_r[26] = 0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[26]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[19]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[18]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[17]+0.1767766952966368*fUpwind_r[4]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[11]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[10]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[9]; 
  Gdrag_r[27] = 0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[27]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[22]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[21]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[20]+0.1767766952966368*fUpwind_r[5]*alphaDrSurf_r[16]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[14]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[13]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[12]; 
  Gdrag_r[28] = 0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[15]; 
  Gdrag_r[29] = 0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[15]; 
  Gdrag_r[30] = 0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[16]*fUpwind_r[23]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[15]; 
  Gdrag_r[31] = 0.1767766952966368*alphaDrSurf_r[0]*fUpwind_r[31]+0.1767766952966368*alphaDrSurf_r[1]*fUpwind_r[30]+0.1767766952966368*alphaDrSurf_r[2]*fUpwind_r[29]+0.1767766952966368*alphaDrSurf_r[3]*fUpwind_r[28]+0.1767766952966368*alphaDrSurf_r[6]*fUpwind_r[25]+0.1767766952966368*alphaDrSurf_r[7]*fUpwind_r[24]+0.1767766952966368*alphaDrSurf_r[8]*fUpwind_r[23]+0.1767766952966368*fUpwind_r[15]*alphaDrSurf_r[16]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[4] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[5] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[6] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[8] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[9] += 0.7071067811865475*Gdrag_r[8]*rdv2-0.7071067811865475*Gdrag_l[8]*rdv2; 
  out[10] += 0.7071067811865475*Gdrag_r[9]*rdv2-0.7071067811865475*Gdrag_l[9]*rdv2; 
  out[11] += 0.7071067811865475*Gdrag_r[10]*rdv2-0.7071067811865475*Gdrag_l[10]*rdv2; 
  out[12] += 0.7071067811865475*Gdrag_r[11]*rdv2-0.7071067811865475*Gdrag_l[11]*rdv2; 
  out[13] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[14] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[15] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[16] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[17] += 0.7071067811865475*Gdrag_r[12]*rdv2-0.7071067811865475*Gdrag_l[12]*rdv2; 
  out[18] += 0.7071067811865475*Gdrag_r[13]*rdv2-0.7071067811865475*Gdrag_l[13]*rdv2; 
  out[19] += 0.7071067811865475*Gdrag_r[14]*rdv2-0.7071067811865475*Gdrag_l[14]*rdv2; 
  out[20] += 0.7071067811865475*Gdrag_r[15]*rdv2-0.7071067811865475*Gdrag_l[15]*rdv2; 
  out[21] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[22] += 0.7071067811865475*Gdrag_r[16]*rdv2-0.7071067811865475*Gdrag_l[16]*rdv2; 
  out[23] += 0.7071067811865475*Gdrag_r[17]*rdv2-0.7071067811865475*Gdrag_l[17]*rdv2; 
  out[24] += 0.7071067811865475*Gdrag_r[18]*rdv2-0.7071067811865475*Gdrag_l[18]*rdv2; 
  out[25] += 0.7071067811865475*Gdrag_r[19]*rdv2-0.7071067811865475*Gdrag_l[19]*rdv2; 
  out[26] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[27] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
  out[28] += 1.224744871391589*Gdrag_r[8]*rdv2+1.224744871391589*Gdrag_l[8]*rdv2; 
  out[29] += 1.224744871391589*Gdrag_r[9]*rdv2+1.224744871391589*Gdrag_l[9]*rdv2; 
  out[30] += 1.224744871391589*Gdrag_r[10]*rdv2+1.224744871391589*Gdrag_l[10]*rdv2; 
  out[31] += 1.224744871391589*Gdrag_r[11]*rdv2+1.224744871391589*Gdrag_l[11]*rdv2; 
  out[32] += 0.7071067811865475*Gdrag_r[20]*rdv2-0.7071067811865475*Gdrag_l[20]*rdv2; 
  out[33] += 0.7071067811865475*Gdrag_r[21]*rdv2-0.7071067811865475*Gdrag_l[21]*rdv2; 
  out[34] += 0.7071067811865475*Gdrag_r[22]*rdv2-0.7071067811865475*Gdrag_l[22]*rdv2; 
  out[35] += 0.7071067811865475*Gdrag_r[23]*rdv2-0.7071067811865475*Gdrag_l[23]*rdv2; 
  out[36] += 0.7071067811865475*Gdrag_r[24]*rdv2-0.7071067811865475*Gdrag_l[24]*rdv2; 
  out[37] += 0.7071067811865475*Gdrag_r[25]*rdv2-0.7071067811865475*Gdrag_l[25]*rdv2; 
  out[38] += 1.224744871391589*Gdrag_r[12]*rdv2+1.224744871391589*Gdrag_l[12]*rdv2; 
  out[39] += 1.224744871391589*Gdrag_r[13]*rdv2+1.224744871391589*Gdrag_l[13]*rdv2; 
  out[40] += 1.224744871391589*Gdrag_r[14]*rdv2+1.224744871391589*Gdrag_l[14]*rdv2; 
  out[41] += 1.224744871391589*Gdrag_r[15]*rdv2+1.224744871391589*Gdrag_l[15]*rdv2; 
  out[42] += 0.7071067811865475*Gdrag_r[26]*rdv2-0.7071067811865475*Gdrag_l[26]*rdv2; 
  out[43] += 1.224744871391589*Gdrag_r[16]*rdv2+1.224744871391589*Gdrag_l[16]*rdv2; 
  out[44] += 1.224744871391589*Gdrag_r[17]*rdv2+1.224744871391589*Gdrag_l[17]*rdv2; 
  out[45] += 1.224744871391589*Gdrag_r[18]*rdv2+1.224744871391589*Gdrag_l[18]*rdv2; 
  out[46] += 1.224744871391589*Gdrag_r[19]*rdv2+1.224744871391589*Gdrag_l[19]*rdv2; 
  out[47] += 0.7071067811865475*Gdrag_r[27]*rdv2-0.7071067811865475*Gdrag_l[27]*rdv2; 
  out[48] += 0.7071067811865475*Gdrag_r[28]*rdv2-0.7071067811865475*Gdrag_l[28]*rdv2; 
  out[49] += 0.7071067811865475*Gdrag_r[29]*rdv2-0.7071067811865475*Gdrag_l[29]*rdv2; 
  out[50] += 0.7071067811865475*Gdrag_r[30]*rdv2-0.7071067811865475*Gdrag_l[30]*rdv2; 
  out[51] += 1.224744871391589*Gdrag_r[20]*rdv2+1.224744871391589*Gdrag_l[20]*rdv2; 
  out[52] += 1.224744871391589*Gdrag_r[21]*rdv2+1.224744871391589*Gdrag_l[21]*rdv2; 
  out[53] += 1.224744871391589*Gdrag_r[22]*rdv2+1.224744871391589*Gdrag_l[22]*rdv2; 
  out[54] += 1.224744871391589*Gdrag_r[23]*rdv2+1.224744871391589*Gdrag_l[23]*rdv2; 
  out[55] += 1.224744871391589*Gdrag_r[24]*rdv2+1.224744871391589*Gdrag_l[24]*rdv2; 
  out[56] += 1.224744871391589*Gdrag_r[25]*rdv2+1.224744871391589*Gdrag_l[25]*rdv2; 
  out[57] += 1.224744871391589*Gdrag_r[26]*rdv2+1.224744871391589*Gdrag_l[26]*rdv2; 
  out[58] += 0.7071067811865475*Gdrag_r[31]*rdv2-0.7071067811865475*Gdrag_l[31]*rdv2; 
  out[59] += 1.224744871391589*Gdrag_r[27]*rdv2+1.224744871391589*Gdrag_l[27]*rdv2; 
  out[60] += 1.224744871391589*Gdrag_r[28]*rdv2+1.224744871391589*Gdrag_l[28]*rdv2; 
  out[61] += 1.224744871391589*Gdrag_r[29]*rdv2+1.224744871391589*Gdrag_l[29]*rdv2; 
  out[62] += 1.224744871391589*Gdrag_r[30]*rdv2+1.224744871391589*Gdrag_l[30]*rdv2; 
  out[63] += 1.224744871391589*Gdrag_r[31]*rdv2+1.224744871391589*Gdrag_l[31]*rdv2; 
} 
