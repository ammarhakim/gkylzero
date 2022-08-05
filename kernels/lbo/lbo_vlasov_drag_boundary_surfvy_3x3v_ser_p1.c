#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_6x_p1_surfx5_eval_quad.h> 
#include <gkyl_basis_ser_6x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void lbo_vlasov_drag_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[6]:         Cell-center coordinates. 
  // dxv[6]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[24]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[8]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[4]; 

  const double *sumNuUy = &nuUSum[8]; 

  double alphaDrSurf[32] = {0.0}; 
  double fUpwindQuad[32] = {0.0};
  double fUpwind[32] = {0.0};;
  double Ghat[32] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[4]+dxv[4])-2.0*sumNuUy[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[4]+dxv[4])-2.0*sumNuUy[1]; 
  alphaDrSurf[2] = nuSum[2]*(2.0*w[4]+dxv[4])-2.0*sumNuUy[2]; 
  alphaDrSurf[3] = nuSum[3]*(2.0*w[4]+dxv[4])-2.0*sumNuUy[3]; 
  alphaDrSurf[6] = 2.0*nuSum[4]*w[4]-2.0*sumNuUy[4]+dxv[4]*nuSum[4]; 
  alphaDrSurf[7] = (2.0*w[4]+dxv[4])*nuSum[5]-2.0*sumNuUy[5]; 
  alphaDrSurf[8] = (2.0*w[4]+dxv[4])*nuSum[6]-2.0*sumNuUy[6]; 
  alphaDrSurf[16] = (2.0*w[4]+dxv[4])*nuSum[7]-2.0*sumNuUy[7]; 

  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_r(fSkin); 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_l(fEdge); 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_l(fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_r(fSkin); 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_l(fEdge); 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_l(fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_r(fSkin); 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_l(fEdge); 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_l(fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_r(fSkin); 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_l(fEdge); 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_l(fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_r(fSkin); 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_r(fSkin); 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_r(fSkin); 
  } else { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_l(fEdge); 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_l(fEdge); 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_l(fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_r(fSkin); 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_r(fSkin); 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_r(fSkin); 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_l(fEdge); 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_l(fEdge); 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_l(fEdge); 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_l(fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_r(fSkin); 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_r(fSkin); 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_r(fSkin); 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_l(fEdge); 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_l(fEdge); 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_l(fEdge); 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_l(fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_r(fSkin); 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_r(fSkin); 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_r(fSkin); 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_r(fSkin); 
  } else { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_l(fEdge); 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_l(fEdge); 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_l(fEdge); 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*(alphaDrSurf[16]*fUpwind[16]+alphaDrSurf[8]*fUpwind[8]+alphaDrSurf[7]*fUpwind[7]+alphaDrSurf[6]*fUpwind[6]+alphaDrSurf[3]*fUpwind[3]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.1767766952966368*(alphaDrSurf[8]*fUpwind[16]+fUpwind[8]*alphaDrSurf[16]+alphaDrSurf[3]*fUpwind[7]+fUpwind[3]*alphaDrSurf[7]+alphaDrSurf[2]*fUpwind[6]+fUpwind[2]*alphaDrSurf[6]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.1767766952966368*(alphaDrSurf[7]*fUpwind[16]+fUpwind[7]*alphaDrSurf[16]+alphaDrSurf[3]*fUpwind[8]+fUpwind[3]*alphaDrSurf[8]+alphaDrSurf[1]*fUpwind[6]+fUpwind[1]*alphaDrSurf[6]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.1767766952966368*(alphaDrSurf[6]*fUpwind[16]+fUpwind[6]*alphaDrSurf[16]+alphaDrSurf[2]*fUpwind[8]+fUpwind[2]*alphaDrSurf[8]+alphaDrSurf[1]*fUpwind[7]+fUpwind[1]*alphaDrSurf[7]+alphaDrSurf[0]*fUpwind[3]+fUpwind[0]*alphaDrSurf[3]); 
  Ghat[4] = 0.1767766952966368*(alphaDrSurf[16]*fUpwind[26]+alphaDrSurf[8]*fUpwind[19]+alphaDrSurf[7]*fUpwind[18]+alphaDrSurf[6]*fUpwind[17]+alphaDrSurf[3]*fUpwind[11]+alphaDrSurf[2]*fUpwind[10]+alphaDrSurf[1]*fUpwind[9]+alphaDrSurf[0]*fUpwind[4]); 
  Ghat[5] = 0.1767766952966368*(alphaDrSurf[16]*fUpwind[27]+alphaDrSurf[8]*fUpwind[22]+alphaDrSurf[7]*fUpwind[21]+alphaDrSurf[6]*fUpwind[20]+alphaDrSurf[3]*fUpwind[14]+alphaDrSurf[2]*fUpwind[13]+alphaDrSurf[1]*fUpwind[12]+alphaDrSurf[0]*fUpwind[5]); 
  Ghat[6] = 0.1767766952966368*(alphaDrSurf[3]*fUpwind[16]+fUpwind[3]*alphaDrSurf[16]+alphaDrSurf[7]*fUpwind[8]+fUpwind[7]*alphaDrSurf[8]+alphaDrSurf[0]*fUpwind[6]+fUpwind[0]*alphaDrSurf[6]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[7] = 0.1767766952966368*(alphaDrSurf[2]*fUpwind[16]+fUpwind[2]*alphaDrSurf[16]+alphaDrSurf[6]*fUpwind[8]+fUpwind[6]*alphaDrSurf[8]+alphaDrSurf[0]*fUpwind[7]+fUpwind[0]*alphaDrSurf[7]+alphaDrSurf[1]*fUpwind[3]+fUpwind[1]*alphaDrSurf[3]); 
  Ghat[8] = 0.1767766952966368*(alphaDrSurf[1]*fUpwind[16]+fUpwind[1]*alphaDrSurf[16]+alphaDrSurf[0]*fUpwind[8]+fUpwind[0]*alphaDrSurf[8]+alphaDrSurf[6]*fUpwind[7]+fUpwind[6]*alphaDrSurf[7]+alphaDrSurf[2]*fUpwind[3]+fUpwind[2]*alphaDrSurf[3]); 
  Ghat[9] = 0.1767766952966368*(alphaDrSurf[8]*fUpwind[26]+alphaDrSurf[16]*fUpwind[19]+alphaDrSurf[3]*fUpwind[18]+alphaDrSurf[2]*fUpwind[17]+alphaDrSurf[7]*fUpwind[11]+alphaDrSurf[6]*fUpwind[10]+alphaDrSurf[0]*fUpwind[9]+alphaDrSurf[1]*fUpwind[4]); 
  Ghat[10] = 0.1767766952966368*(alphaDrSurf[7]*fUpwind[26]+alphaDrSurf[3]*fUpwind[19]+alphaDrSurf[16]*fUpwind[18]+alphaDrSurf[1]*fUpwind[17]+alphaDrSurf[8]*fUpwind[11]+alphaDrSurf[0]*fUpwind[10]+alphaDrSurf[6]*fUpwind[9]+alphaDrSurf[2]*fUpwind[4]); 
  Ghat[11] = 0.1767766952966368*(alphaDrSurf[6]*fUpwind[26]+alphaDrSurf[2]*fUpwind[19]+alphaDrSurf[1]*fUpwind[18]+alphaDrSurf[16]*fUpwind[17]+alphaDrSurf[0]*fUpwind[11]+alphaDrSurf[8]*fUpwind[10]+alphaDrSurf[7]*fUpwind[9]+alphaDrSurf[3]*fUpwind[4]); 
  Ghat[12] = 0.1767766952966368*(alphaDrSurf[8]*fUpwind[27]+alphaDrSurf[16]*fUpwind[22]+alphaDrSurf[3]*fUpwind[21]+alphaDrSurf[2]*fUpwind[20]+alphaDrSurf[7]*fUpwind[14]+alphaDrSurf[6]*fUpwind[13]+alphaDrSurf[0]*fUpwind[12]+alphaDrSurf[1]*fUpwind[5]); 
  Ghat[13] = 0.1767766952966368*(alphaDrSurf[7]*fUpwind[27]+alphaDrSurf[3]*fUpwind[22]+alphaDrSurf[16]*fUpwind[21]+alphaDrSurf[1]*fUpwind[20]+alphaDrSurf[8]*fUpwind[14]+alphaDrSurf[0]*fUpwind[13]+alphaDrSurf[6]*fUpwind[12]+alphaDrSurf[2]*fUpwind[5]); 
  Ghat[14] = 0.1767766952966368*(alphaDrSurf[6]*fUpwind[27]+alphaDrSurf[2]*fUpwind[22]+alphaDrSurf[1]*fUpwind[21]+alphaDrSurf[16]*fUpwind[20]+alphaDrSurf[0]*fUpwind[14]+alphaDrSurf[8]*fUpwind[13]+alphaDrSurf[7]*fUpwind[12]+alphaDrSurf[3]*fUpwind[5]); 
  Ghat[15] = 0.1767766952966368*(alphaDrSurf[16]*fUpwind[31]+alphaDrSurf[8]*fUpwind[30]+alphaDrSurf[7]*fUpwind[29]+alphaDrSurf[6]*fUpwind[28]+alphaDrSurf[3]*fUpwind[25]+alphaDrSurf[2]*fUpwind[24]+alphaDrSurf[1]*fUpwind[23]+alphaDrSurf[0]*fUpwind[15]); 
  Ghat[16] = 0.1767766952966368*(alphaDrSurf[0]*fUpwind[16]+fUpwind[0]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[8]+fUpwind[1]*alphaDrSurf[8]+alphaDrSurf[2]*fUpwind[7]+fUpwind[2]*alphaDrSurf[7]+alphaDrSurf[3]*fUpwind[6]+fUpwind[3]*alphaDrSurf[6]); 
  Ghat[17] = 0.1767766952966368*(alphaDrSurf[3]*fUpwind[26]+alphaDrSurf[7]*fUpwind[19]+alphaDrSurf[8]*fUpwind[18]+alphaDrSurf[0]*fUpwind[17]+fUpwind[11]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[10]+alphaDrSurf[2]*fUpwind[9]+fUpwind[4]*alphaDrSurf[6]); 
  Ghat[18] = 0.1767766952966368*(alphaDrSurf[2]*fUpwind[26]+alphaDrSurf[6]*fUpwind[19]+alphaDrSurf[0]*fUpwind[18]+alphaDrSurf[8]*fUpwind[17]+fUpwind[10]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[11]+alphaDrSurf[3]*fUpwind[9]+fUpwind[4]*alphaDrSurf[7]); 
  Ghat[19] = 0.1767766952966368*(alphaDrSurf[1]*fUpwind[26]+alphaDrSurf[0]*fUpwind[19]+alphaDrSurf[6]*fUpwind[18]+alphaDrSurf[7]*fUpwind[17]+fUpwind[9]*alphaDrSurf[16]+alphaDrSurf[2]*fUpwind[11]+alphaDrSurf[3]*fUpwind[10]+fUpwind[4]*alphaDrSurf[8]); 
  Ghat[20] = 0.1767766952966368*(alphaDrSurf[3]*fUpwind[27]+alphaDrSurf[7]*fUpwind[22]+alphaDrSurf[8]*fUpwind[21]+alphaDrSurf[0]*fUpwind[20]+fUpwind[14]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[13]+alphaDrSurf[2]*fUpwind[12]+fUpwind[5]*alphaDrSurf[6]); 
  Ghat[21] = 0.1767766952966368*(alphaDrSurf[2]*fUpwind[27]+alphaDrSurf[6]*fUpwind[22]+alphaDrSurf[0]*fUpwind[21]+alphaDrSurf[8]*fUpwind[20]+fUpwind[13]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[14]+alphaDrSurf[3]*fUpwind[12]+fUpwind[5]*alphaDrSurf[7]); 
  Ghat[22] = 0.1767766952966368*(alphaDrSurf[1]*fUpwind[27]+alphaDrSurf[0]*fUpwind[22]+alphaDrSurf[6]*fUpwind[21]+alphaDrSurf[7]*fUpwind[20]+fUpwind[12]*alphaDrSurf[16]+alphaDrSurf[2]*fUpwind[14]+alphaDrSurf[3]*fUpwind[13]+fUpwind[5]*alphaDrSurf[8]); 
  Ghat[23] = 0.1767766952966368*(alphaDrSurf[8]*fUpwind[31]+alphaDrSurf[16]*fUpwind[30]+alphaDrSurf[3]*fUpwind[29]+alphaDrSurf[2]*fUpwind[28]+alphaDrSurf[7]*fUpwind[25]+alphaDrSurf[6]*fUpwind[24]+alphaDrSurf[0]*fUpwind[23]+alphaDrSurf[1]*fUpwind[15]); 
  Ghat[24] = 0.1767766952966368*(alphaDrSurf[7]*fUpwind[31]+alphaDrSurf[3]*fUpwind[30]+alphaDrSurf[16]*fUpwind[29]+alphaDrSurf[1]*fUpwind[28]+alphaDrSurf[8]*fUpwind[25]+alphaDrSurf[0]*fUpwind[24]+alphaDrSurf[6]*fUpwind[23]+alphaDrSurf[2]*fUpwind[15]); 
  Ghat[25] = 0.1767766952966368*(alphaDrSurf[6]*fUpwind[31]+alphaDrSurf[2]*fUpwind[30]+alphaDrSurf[1]*fUpwind[29]+alphaDrSurf[16]*fUpwind[28]+alphaDrSurf[0]*fUpwind[25]+alphaDrSurf[8]*fUpwind[24]+alphaDrSurf[7]*fUpwind[23]+alphaDrSurf[3]*fUpwind[15]); 
  Ghat[26] = 0.1767766952966368*(alphaDrSurf[0]*fUpwind[26]+alphaDrSurf[1]*fUpwind[19]+alphaDrSurf[2]*fUpwind[18]+alphaDrSurf[3]*fUpwind[17]+fUpwind[4]*alphaDrSurf[16]+alphaDrSurf[6]*fUpwind[11]+alphaDrSurf[7]*fUpwind[10]+alphaDrSurf[8]*fUpwind[9]); 
  Ghat[27] = 0.1767766952966368*(alphaDrSurf[0]*fUpwind[27]+alphaDrSurf[1]*fUpwind[22]+alphaDrSurf[2]*fUpwind[21]+alphaDrSurf[3]*fUpwind[20]+fUpwind[5]*alphaDrSurf[16]+alphaDrSurf[6]*fUpwind[14]+alphaDrSurf[7]*fUpwind[13]+alphaDrSurf[8]*fUpwind[12]); 
  Ghat[28] = 0.1767766952966368*(alphaDrSurf[3]*fUpwind[31]+alphaDrSurf[7]*fUpwind[30]+alphaDrSurf[8]*fUpwind[29]+alphaDrSurf[0]*fUpwind[28]+alphaDrSurf[16]*fUpwind[25]+alphaDrSurf[1]*fUpwind[24]+alphaDrSurf[2]*fUpwind[23]+alphaDrSurf[6]*fUpwind[15]); 
  Ghat[29] = 0.1767766952966368*(alphaDrSurf[2]*fUpwind[31]+alphaDrSurf[6]*fUpwind[30]+alphaDrSurf[0]*fUpwind[29]+alphaDrSurf[8]*fUpwind[28]+alphaDrSurf[1]*fUpwind[25]+alphaDrSurf[16]*fUpwind[24]+alphaDrSurf[3]*fUpwind[23]+alphaDrSurf[7]*fUpwind[15]); 
  Ghat[30] = 0.1767766952966368*(alphaDrSurf[1]*fUpwind[31]+alphaDrSurf[0]*fUpwind[30]+alphaDrSurf[6]*fUpwind[29]+alphaDrSurf[7]*fUpwind[28]+alphaDrSurf[2]*fUpwind[25]+alphaDrSurf[3]*fUpwind[24]+alphaDrSurf[16]*fUpwind[23]+alphaDrSurf[8]*fUpwind[15]); 
  Ghat[31] = 0.1767766952966368*(alphaDrSurf[0]*fUpwind[31]+alphaDrSurf[1]*fUpwind[30]+alphaDrSurf[2]*fUpwind[29]+alphaDrSurf[3]*fUpwind[28]+alphaDrSurf[6]*fUpwind[25]+alphaDrSurf[7]*fUpwind[24]+alphaDrSurf[8]*fUpwind[23]+fUpwind[15]*alphaDrSurf[16]); 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[5] += 1.224744871391589*Ghat[0]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[10] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[13] += 1.224744871391589*Ghat[1]*rdv2; 
  out[14] += 1.224744871391589*Ghat[2]*rdv2; 
  out[15] += 1.224744871391589*Ghat[3]*rdv2; 
  out[16] += 1.224744871391589*Ghat[4]*rdv2; 
  out[17] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[18] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[19] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[20] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[21] += 1.224744871391589*Ghat[5]*rdv2; 
  out[22] += 0.7071067811865475*Ghat[16]*rdv2; 
  out[23] += 0.7071067811865475*Ghat[17]*rdv2; 
  out[24] += 0.7071067811865475*Ghat[18]*rdv2; 
  out[25] += 0.7071067811865475*Ghat[19]*rdv2; 
  out[26] += 1.224744871391589*Ghat[6]*rdv2; 
  out[27] += 1.224744871391589*Ghat[7]*rdv2; 
  out[28] += 1.224744871391589*Ghat[8]*rdv2; 
  out[29] += 1.224744871391589*Ghat[9]*rdv2; 
  out[30] += 1.224744871391589*Ghat[10]*rdv2; 
  out[31] += 1.224744871391589*Ghat[11]*rdv2; 
  out[32] += 0.7071067811865475*Ghat[20]*rdv2; 
  out[33] += 0.7071067811865475*Ghat[21]*rdv2; 
  out[34] += 0.7071067811865475*Ghat[22]*rdv2; 
  out[35] += 0.7071067811865475*Ghat[23]*rdv2; 
  out[36] += 0.7071067811865475*Ghat[24]*rdv2; 
  out[37] += 0.7071067811865475*Ghat[25]*rdv2; 
  out[38] += 1.224744871391589*Ghat[12]*rdv2; 
  out[39] += 1.224744871391589*Ghat[13]*rdv2; 
  out[40] += 1.224744871391589*Ghat[14]*rdv2; 
  out[41] += 1.224744871391589*Ghat[15]*rdv2; 
  out[42] += 0.7071067811865475*Ghat[26]*rdv2; 
  out[43] += 1.224744871391589*Ghat[16]*rdv2; 
  out[44] += 1.224744871391589*Ghat[17]*rdv2; 
  out[45] += 1.224744871391589*Ghat[18]*rdv2; 
  out[46] += 1.224744871391589*Ghat[19]*rdv2; 
  out[47] += 0.7071067811865475*Ghat[27]*rdv2; 
  out[48] += 0.7071067811865475*Ghat[28]*rdv2; 
  out[49] += 0.7071067811865475*Ghat[29]*rdv2; 
  out[50] += 0.7071067811865475*Ghat[30]*rdv2; 
  out[51] += 1.224744871391589*Ghat[20]*rdv2; 
  out[52] += 1.224744871391589*Ghat[21]*rdv2; 
  out[53] += 1.224744871391589*Ghat[22]*rdv2; 
  out[54] += 1.224744871391589*Ghat[23]*rdv2; 
  out[55] += 1.224744871391589*Ghat[24]*rdv2; 
  out[56] += 1.224744871391589*Ghat[25]*rdv2; 
  out[57] += 1.224744871391589*Ghat[26]*rdv2; 
  out[58] += 0.7071067811865475*Ghat[31]*rdv2; 
  out[59] += 1.224744871391589*Ghat[27]*rdv2; 
  out[60] += 1.224744871391589*Ghat[28]*rdv2; 
  out[61] += 1.224744871391589*Ghat[29]*rdv2; 
  out[62] += 1.224744871391589*Ghat[30]*rdv2; 
  out[63] += 1.224744871391589*Ghat[31]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[4]-1.0*dxv[4])-2.0*sumNuUy[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[4]-1.0*dxv[4])-2.0*sumNuUy[1]; 
  alphaDrSurf[2] = nuSum[2]*(2.0*w[4]-1.0*dxv[4])-2.0*sumNuUy[2]; 
  alphaDrSurf[3] = nuSum[3]*(2.0*w[4]-1.0*dxv[4])-2.0*sumNuUy[3]; 
  alphaDrSurf[6] = 2.0*nuSum[4]*w[4]-2.0*sumNuUy[4]-1.0*dxv[4]*nuSum[4]; 
  alphaDrSurf[7] = (2.0*w[4]-1.0*dxv[4])*nuSum[5]-2.0*sumNuUy[5]; 
  alphaDrSurf[8] = (2.0*w[4]-1.0*dxv[4])*nuSum[6]-2.0*sumNuUy[6]; 
  alphaDrSurf[16] = (2.0*w[4]-1.0*dxv[4])*nuSum[7]-2.0*sumNuUy[7]; 

  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_r(fEdge); 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_6x_p1_surfx5_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_6x_p1_surfx5_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = ser_6x_p1_surfx5_eval_quad_node_2_l(fSkin); 
    fUpwindQuad[3] = ser_6x_p1_surfx5_eval_quad_node_3_l(fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_r(fEdge); 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_6x_p1_surfx5_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = ser_6x_p1_surfx5_eval_quad_node_5_l(fSkin); 
    fUpwindQuad[6] = ser_6x_p1_surfx5_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = ser_6x_p1_surfx5_eval_quad_node_7_l(fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_r(fEdge); 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_6x_p1_surfx5_eval_quad_node_8_l(fSkin); 
    fUpwindQuad[9] = ser_6x_p1_surfx5_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = ser_6x_p1_surfx5_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = ser_6x_p1_surfx5_eval_quad_node_11_l(fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_r(fEdge); 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_6x_p1_surfx5_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = ser_6x_p1_surfx5_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = ser_6x_p1_surfx5_eval_quad_node_14_l(fSkin); 
    fUpwindQuad[15] = ser_6x_p1_surfx5_eval_quad_node_15_l(fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_r(fEdge); 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_r(fEdge); 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_r(fEdge); 
  } else { 
    fUpwindQuad[16] = ser_6x_p1_surfx5_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = ser_6x_p1_surfx5_eval_quad_node_17_l(fSkin); 
    fUpwindQuad[18] = ser_6x_p1_surfx5_eval_quad_node_18_l(fSkin); 
    fUpwindQuad[19] = ser_6x_p1_surfx5_eval_quad_node_19_l(fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_r(fEdge); 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_r(fEdge); 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_r(fEdge); 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[20] = ser_6x_p1_surfx5_eval_quad_node_20_l(fSkin); 
    fUpwindQuad[21] = ser_6x_p1_surfx5_eval_quad_node_21_l(fSkin); 
    fUpwindQuad[22] = ser_6x_p1_surfx5_eval_quad_node_22_l(fSkin); 
    fUpwindQuad[23] = ser_6x_p1_surfx5_eval_quad_node_23_l(fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_r(fEdge); 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_r(fEdge); 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_r(fEdge); 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_6x_p1_surfx5_eval_quad_node_24_l(fSkin); 
    fUpwindQuad[25] = ser_6x_p1_surfx5_eval_quad_node_25_l(fSkin); 
    fUpwindQuad[26] = ser_6x_p1_surfx5_eval_quad_node_26_l(fSkin); 
    fUpwindQuad[27] = ser_6x_p1_surfx5_eval_quad_node_27_l(fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_r(fEdge); 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_r(fEdge); 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_r(fEdge); 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_r(fEdge); 
  } else { 
    fUpwindQuad[28] = ser_6x_p1_surfx5_eval_quad_node_28_l(fSkin); 
    fUpwindQuad[29] = ser_6x_p1_surfx5_eval_quad_node_29_l(fSkin); 
    fUpwindQuad[30] = ser_6x_p1_surfx5_eval_quad_node_30_l(fSkin); 
    fUpwindQuad[31] = ser_6x_p1_surfx5_eval_quad_node_31_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_6x_p1_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.1767766952966368*(alphaDrSurf[16]*fUpwind[16]+alphaDrSurf[8]*fUpwind[8]+alphaDrSurf[7]*fUpwind[7]+alphaDrSurf[6]*fUpwind[6]+alphaDrSurf[3]*fUpwind[3]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = 0.1767766952966368*(alphaDrSurf[8]*fUpwind[16]+fUpwind[8]*alphaDrSurf[16]+alphaDrSurf[3]*fUpwind[7]+fUpwind[3]*alphaDrSurf[7]+alphaDrSurf[2]*fUpwind[6]+fUpwind[2]*alphaDrSurf[6]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = 0.1767766952966368*(alphaDrSurf[7]*fUpwind[16]+fUpwind[7]*alphaDrSurf[16]+alphaDrSurf[3]*fUpwind[8]+fUpwind[3]*alphaDrSurf[8]+alphaDrSurf[1]*fUpwind[6]+fUpwind[1]*alphaDrSurf[6]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[3] = 0.1767766952966368*(alphaDrSurf[6]*fUpwind[16]+fUpwind[6]*alphaDrSurf[16]+alphaDrSurf[2]*fUpwind[8]+fUpwind[2]*alphaDrSurf[8]+alphaDrSurf[1]*fUpwind[7]+fUpwind[1]*alphaDrSurf[7]+alphaDrSurf[0]*fUpwind[3]+fUpwind[0]*alphaDrSurf[3]); 
  Ghat[4] = 0.1767766952966368*(alphaDrSurf[16]*fUpwind[26]+alphaDrSurf[8]*fUpwind[19]+alphaDrSurf[7]*fUpwind[18]+alphaDrSurf[6]*fUpwind[17]+alphaDrSurf[3]*fUpwind[11]+alphaDrSurf[2]*fUpwind[10]+alphaDrSurf[1]*fUpwind[9]+alphaDrSurf[0]*fUpwind[4]); 
  Ghat[5] = 0.1767766952966368*(alphaDrSurf[16]*fUpwind[27]+alphaDrSurf[8]*fUpwind[22]+alphaDrSurf[7]*fUpwind[21]+alphaDrSurf[6]*fUpwind[20]+alphaDrSurf[3]*fUpwind[14]+alphaDrSurf[2]*fUpwind[13]+alphaDrSurf[1]*fUpwind[12]+alphaDrSurf[0]*fUpwind[5]); 
  Ghat[6] = 0.1767766952966368*(alphaDrSurf[3]*fUpwind[16]+fUpwind[3]*alphaDrSurf[16]+alphaDrSurf[7]*fUpwind[8]+fUpwind[7]*alphaDrSurf[8]+alphaDrSurf[0]*fUpwind[6]+fUpwind[0]*alphaDrSurf[6]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[7] = 0.1767766952966368*(alphaDrSurf[2]*fUpwind[16]+fUpwind[2]*alphaDrSurf[16]+alphaDrSurf[6]*fUpwind[8]+fUpwind[6]*alphaDrSurf[8]+alphaDrSurf[0]*fUpwind[7]+fUpwind[0]*alphaDrSurf[7]+alphaDrSurf[1]*fUpwind[3]+fUpwind[1]*alphaDrSurf[3]); 
  Ghat[8] = 0.1767766952966368*(alphaDrSurf[1]*fUpwind[16]+fUpwind[1]*alphaDrSurf[16]+alphaDrSurf[0]*fUpwind[8]+fUpwind[0]*alphaDrSurf[8]+alphaDrSurf[6]*fUpwind[7]+fUpwind[6]*alphaDrSurf[7]+alphaDrSurf[2]*fUpwind[3]+fUpwind[2]*alphaDrSurf[3]); 
  Ghat[9] = 0.1767766952966368*(alphaDrSurf[8]*fUpwind[26]+alphaDrSurf[16]*fUpwind[19]+alphaDrSurf[3]*fUpwind[18]+alphaDrSurf[2]*fUpwind[17]+alphaDrSurf[7]*fUpwind[11]+alphaDrSurf[6]*fUpwind[10]+alphaDrSurf[0]*fUpwind[9]+alphaDrSurf[1]*fUpwind[4]); 
  Ghat[10] = 0.1767766952966368*(alphaDrSurf[7]*fUpwind[26]+alphaDrSurf[3]*fUpwind[19]+alphaDrSurf[16]*fUpwind[18]+alphaDrSurf[1]*fUpwind[17]+alphaDrSurf[8]*fUpwind[11]+alphaDrSurf[0]*fUpwind[10]+alphaDrSurf[6]*fUpwind[9]+alphaDrSurf[2]*fUpwind[4]); 
  Ghat[11] = 0.1767766952966368*(alphaDrSurf[6]*fUpwind[26]+alphaDrSurf[2]*fUpwind[19]+alphaDrSurf[1]*fUpwind[18]+alphaDrSurf[16]*fUpwind[17]+alphaDrSurf[0]*fUpwind[11]+alphaDrSurf[8]*fUpwind[10]+alphaDrSurf[7]*fUpwind[9]+alphaDrSurf[3]*fUpwind[4]); 
  Ghat[12] = 0.1767766952966368*(alphaDrSurf[8]*fUpwind[27]+alphaDrSurf[16]*fUpwind[22]+alphaDrSurf[3]*fUpwind[21]+alphaDrSurf[2]*fUpwind[20]+alphaDrSurf[7]*fUpwind[14]+alphaDrSurf[6]*fUpwind[13]+alphaDrSurf[0]*fUpwind[12]+alphaDrSurf[1]*fUpwind[5]); 
  Ghat[13] = 0.1767766952966368*(alphaDrSurf[7]*fUpwind[27]+alphaDrSurf[3]*fUpwind[22]+alphaDrSurf[16]*fUpwind[21]+alphaDrSurf[1]*fUpwind[20]+alphaDrSurf[8]*fUpwind[14]+alphaDrSurf[0]*fUpwind[13]+alphaDrSurf[6]*fUpwind[12]+alphaDrSurf[2]*fUpwind[5]); 
  Ghat[14] = 0.1767766952966368*(alphaDrSurf[6]*fUpwind[27]+alphaDrSurf[2]*fUpwind[22]+alphaDrSurf[1]*fUpwind[21]+alphaDrSurf[16]*fUpwind[20]+alphaDrSurf[0]*fUpwind[14]+alphaDrSurf[8]*fUpwind[13]+alphaDrSurf[7]*fUpwind[12]+alphaDrSurf[3]*fUpwind[5]); 
  Ghat[15] = 0.1767766952966368*(alphaDrSurf[16]*fUpwind[31]+alphaDrSurf[8]*fUpwind[30]+alphaDrSurf[7]*fUpwind[29]+alphaDrSurf[6]*fUpwind[28]+alphaDrSurf[3]*fUpwind[25]+alphaDrSurf[2]*fUpwind[24]+alphaDrSurf[1]*fUpwind[23]+alphaDrSurf[0]*fUpwind[15]); 
  Ghat[16] = 0.1767766952966368*(alphaDrSurf[0]*fUpwind[16]+fUpwind[0]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[8]+fUpwind[1]*alphaDrSurf[8]+alphaDrSurf[2]*fUpwind[7]+fUpwind[2]*alphaDrSurf[7]+alphaDrSurf[3]*fUpwind[6]+fUpwind[3]*alphaDrSurf[6]); 
  Ghat[17] = 0.1767766952966368*(alphaDrSurf[3]*fUpwind[26]+alphaDrSurf[7]*fUpwind[19]+alphaDrSurf[8]*fUpwind[18]+alphaDrSurf[0]*fUpwind[17]+fUpwind[11]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[10]+alphaDrSurf[2]*fUpwind[9]+fUpwind[4]*alphaDrSurf[6]); 
  Ghat[18] = 0.1767766952966368*(alphaDrSurf[2]*fUpwind[26]+alphaDrSurf[6]*fUpwind[19]+alphaDrSurf[0]*fUpwind[18]+alphaDrSurf[8]*fUpwind[17]+fUpwind[10]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[11]+alphaDrSurf[3]*fUpwind[9]+fUpwind[4]*alphaDrSurf[7]); 
  Ghat[19] = 0.1767766952966368*(alphaDrSurf[1]*fUpwind[26]+alphaDrSurf[0]*fUpwind[19]+alphaDrSurf[6]*fUpwind[18]+alphaDrSurf[7]*fUpwind[17]+fUpwind[9]*alphaDrSurf[16]+alphaDrSurf[2]*fUpwind[11]+alphaDrSurf[3]*fUpwind[10]+fUpwind[4]*alphaDrSurf[8]); 
  Ghat[20] = 0.1767766952966368*(alphaDrSurf[3]*fUpwind[27]+alphaDrSurf[7]*fUpwind[22]+alphaDrSurf[8]*fUpwind[21]+alphaDrSurf[0]*fUpwind[20]+fUpwind[14]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[13]+alphaDrSurf[2]*fUpwind[12]+fUpwind[5]*alphaDrSurf[6]); 
  Ghat[21] = 0.1767766952966368*(alphaDrSurf[2]*fUpwind[27]+alphaDrSurf[6]*fUpwind[22]+alphaDrSurf[0]*fUpwind[21]+alphaDrSurf[8]*fUpwind[20]+fUpwind[13]*alphaDrSurf[16]+alphaDrSurf[1]*fUpwind[14]+alphaDrSurf[3]*fUpwind[12]+fUpwind[5]*alphaDrSurf[7]); 
  Ghat[22] = 0.1767766952966368*(alphaDrSurf[1]*fUpwind[27]+alphaDrSurf[0]*fUpwind[22]+alphaDrSurf[6]*fUpwind[21]+alphaDrSurf[7]*fUpwind[20]+fUpwind[12]*alphaDrSurf[16]+alphaDrSurf[2]*fUpwind[14]+alphaDrSurf[3]*fUpwind[13]+fUpwind[5]*alphaDrSurf[8]); 
  Ghat[23] = 0.1767766952966368*(alphaDrSurf[8]*fUpwind[31]+alphaDrSurf[16]*fUpwind[30]+alphaDrSurf[3]*fUpwind[29]+alphaDrSurf[2]*fUpwind[28]+alphaDrSurf[7]*fUpwind[25]+alphaDrSurf[6]*fUpwind[24]+alphaDrSurf[0]*fUpwind[23]+alphaDrSurf[1]*fUpwind[15]); 
  Ghat[24] = 0.1767766952966368*(alphaDrSurf[7]*fUpwind[31]+alphaDrSurf[3]*fUpwind[30]+alphaDrSurf[16]*fUpwind[29]+alphaDrSurf[1]*fUpwind[28]+alphaDrSurf[8]*fUpwind[25]+alphaDrSurf[0]*fUpwind[24]+alphaDrSurf[6]*fUpwind[23]+alphaDrSurf[2]*fUpwind[15]); 
  Ghat[25] = 0.1767766952966368*(alphaDrSurf[6]*fUpwind[31]+alphaDrSurf[2]*fUpwind[30]+alphaDrSurf[1]*fUpwind[29]+alphaDrSurf[16]*fUpwind[28]+alphaDrSurf[0]*fUpwind[25]+alphaDrSurf[8]*fUpwind[24]+alphaDrSurf[7]*fUpwind[23]+alphaDrSurf[3]*fUpwind[15]); 
  Ghat[26] = 0.1767766952966368*(alphaDrSurf[0]*fUpwind[26]+alphaDrSurf[1]*fUpwind[19]+alphaDrSurf[2]*fUpwind[18]+alphaDrSurf[3]*fUpwind[17]+fUpwind[4]*alphaDrSurf[16]+alphaDrSurf[6]*fUpwind[11]+alphaDrSurf[7]*fUpwind[10]+alphaDrSurf[8]*fUpwind[9]); 
  Ghat[27] = 0.1767766952966368*(alphaDrSurf[0]*fUpwind[27]+alphaDrSurf[1]*fUpwind[22]+alphaDrSurf[2]*fUpwind[21]+alphaDrSurf[3]*fUpwind[20]+fUpwind[5]*alphaDrSurf[16]+alphaDrSurf[6]*fUpwind[14]+alphaDrSurf[7]*fUpwind[13]+alphaDrSurf[8]*fUpwind[12]); 
  Ghat[28] = 0.1767766952966368*(alphaDrSurf[3]*fUpwind[31]+alphaDrSurf[7]*fUpwind[30]+alphaDrSurf[8]*fUpwind[29]+alphaDrSurf[0]*fUpwind[28]+alphaDrSurf[16]*fUpwind[25]+alphaDrSurf[1]*fUpwind[24]+alphaDrSurf[2]*fUpwind[23]+alphaDrSurf[6]*fUpwind[15]); 
  Ghat[29] = 0.1767766952966368*(alphaDrSurf[2]*fUpwind[31]+alphaDrSurf[6]*fUpwind[30]+alphaDrSurf[0]*fUpwind[29]+alphaDrSurf[8]*fUpwind[28]+alphaDrSurf[1]*fUpwind[25]+alphaDrSurf[16]*fUpwind[24]+alphaDrSurf[3]*fUpwind[23]+alphaDrSurf[7]*fUpwind[15]); 
  Ghat[30] = 0.1767766952966368*(alphaDrSurf[1]*fUpwind[31]+alphaDrSurf[0]*fUpwind[30]+alphaDrSurf[6]*fUpwind[29]+alphaDrSurf[7]*fUpwind[28]+alphaDrSurf[2]*fUpwind[25]+alphaDrSurf[3]*fUpwind[24]+alphaDrSurf[16]*fUpwind[23]+alphaDrSurf[8]*fUpwind[15]); 
  Ghat[31] = 0.1767766952966368*(alphaDrSurf[0]*fUpwind[31]+alphaDrSurf[1]*fUpwind[30]+alphaDrSurf[2]*fUpwind[29]+alphaDrSurf[3]*fUpwind[28]+alphaDrSurf[6]*fUpwind[25]+alphaDrSurf[7]*fUpwind[24]+alphaDrSurf[8]*fUpwind[23]+fUpwind[15]*alphaDrSurf[16]); 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[5] += 1.224744871391589*Ghat[0]*rdv2; 
  out[6] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[9] += -0.7071067811865475*Ghat[8]*rdv2; 
  out[10] += -0.7071067811865475*Ghat[9]*rdv2; 
  out[11] += -0.7071067811865475*Ghat[10]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[11]*rdv2; 
  out[13] += 1.224744871391589*Ghat[1]*rdv2; 
  out[14] += 1.224744871391589*Ghat[2]*rdv2; 
  out[15] += 1.224744871391589*Ghat[3]*rdv2; 
  out[16] += 1.224744871391589*Ghat[4]*rdv2; 
  out[17] += -0.7071067811865475*Ghat[12]*rdv2; 
  out[18] += -0.7071067811865475*Ghat[13]*rdv2; 
  out[19] += -0.7071067811865475*Ghat[14]*rdv2; 
  out[20] += -0.7071067811865475*Ghat[15]*rdv2; 
  out[21] += 1.224744871391589*Ghat[5]*rdv2; 
  out[22] += -0.7071067811865475*Ghat[16]*rdv2; 
  out[23] += -0.7071067811865475*Ghat[17]*rdv2; 
  out[24] += -0.7071067811865475*Ghat[18]*rdv2; 
  out[25] += -0.7071067811865475*Ghat[19]*rdv2; 
  out[26] += 1.224744871391589*Ghat[6]*rdv2; 
  out[27] += 1.224744871391589*Ghat[7]*rdv2; 
  out[28] += 1.224744871391589*Ghat[8]*rdv2; 
  out[29] += 1.224744871391589*Ghat[9]*rdv2; 
  out[30] += 1.224744871391589*Ghat[10]*rdv2; 
  out[31] += 1.224744871391589*Ghat[11]*rdv2; 
  out[32] += -0.7071067811865475*Ghat[20]*rdv2; 
  out[33] += -0.7071067811865475*Ghat[21]*rdv2; 
  out[34] += -0.7071067811865475*Ghat[22]*rdv2; 
  out[35] += -0.7071067811865475*Ghat[23]*rdv2; 
  out[36] += -0.7071067811865475*Ghat[24]*rdv2; 
  out[37] += -0.7071067811865475*Ghat[25]*rdv2; 
  out[38] += 1.224744871391589*Ghat[12]*rdv2; 
  out[39] += 1.224744871391589*Ghat[13]*rdv2; 
  out[40] += 1.224744871391589*Ghat[14]*rdv2; 
  out[41] += 1.224744871391589*Ghat[15]*rdv2; 
  out[42] += -0.7071067811865475*Ghat[26]*rdv2; 
  out[43] += 1.224744871391589*Ghat[16]*rdv2; 
  out[44] += 1.224744871391589*Ghat[17]*rdv2; 
  out[45] += 1.224744871391589*Ghat[18]*rdv2; 
  out[46] += 1.224744871391589*Ghat[19]*rdv2; 
  out[47] += -0.7071067811865475*Ghat[27]*rdv2; 
  out[48] += -0.7071067811865475*Ghat[28]*rdv2; 
  out[49] += -0.7071067811865475*Ghat[29]*rdv2; 
  out[50] += -0.7071067811865475*Ghat[30]*rdv2; 
  out[51] += 1.224744871391589*Ghat[20]*rdv2; 
  out[52] += 1.224744871391589*Ghat[21]*rdv2; 
  out[53] += 1.224744871391589*Ghat[22]*rdv2; 
  out[54] += 1.224744871391589*Ghat[23]*rdv2; 
  out[55] += 1.224744871391589*Ghat[24]*rdv2; 
  out[56] += 1.224744871391589*Ghat[25]*rdv2; 
  out[57] += 1.224744871391589*Ghat[26]*rdv2; 
  out[58] += -0.7071067811865475*Ghat[31]*rdv2; 
  out[59] += 1.224744871391589*Ghat[27]*rdv2; 
  out[60] += 1.224744871391589*Ghat[28]*rdv2; 
  out[61] += 1.224744871391589*Ghat[29]*rdv2; 
  out[62] += 1.224744871391589*Ghat[30]*rdv2; 
  out[63] += 1.224744871391589*Ghat[31]*rdv2; 

  } 
} 
