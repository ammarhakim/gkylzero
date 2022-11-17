#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_boundary_surfvpar_3x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[5]:     cell-center coordinates. 
  // dxv[5]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[40]: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 

  const double *nuUSum = nuPrimMomsSum;

  double rdv2 = 2.0/dxv[3]; 


  double alphaDrSurf[48] = {0.0}; 
  double fUpwindQuad[81] = {0.0};
  double fUpwind[48] = {0.0};;
  double Ghat[48] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[3]+1.414213562373095*dxv[3])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[3]+1.414213562373095*dxv[3])-2.828427124746191*nuUSum[1]); 
  alphaDrSurf[2] = 0.5*(nuSum[2]*(2.828427124746191*w[3]+1.414213562373095*dxv[3])-2.828427124746191*nuUSum[2]); 
  alphaDrSurf[3] = 0.5*(2.828427124746191*nuSum[3]*w[3]-2.828427124746191*nuUSum[3]+1.414213562373095*dxv[3]*nuSum[3]); 
  alphaDrSurf[5] = -0.5*(2.828427124746191*nuUSum[4]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[4]); 
  alphaDrSurf[6] = -0.5*(2.828427124746191*nuUSum[5]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[5]); 
  alphaDrSurf[7] = -0.5*(2.828427124746191*nuUSum[6]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[6]); 
  alphaDrSurf[11] = -0.5*(2.828427124746191*nuUSum[7]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[7]); 
  alphaDrSurf[12] = -0.5*(2.828427124746191*nuUSum[8]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[8]); 
  alphaDrSurf[13] = -0.5*(2.828427124746191*nuUSum[9]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[9]); 
  alphaDrSurf[15] = -0.5*(2.828427124746191*nuUSum[10]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[10]); 
  alphaDrSurf[19] = -0.5*(2.828427124746191*nuUSum[11]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[11]); 
  alphaDrSurf[20] = -0.5*(2.828427124746191*nuUSum[12]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[12]); 
  alphaDrSurf[21] = -0.5*(2.828427124746191*nuUSum[13]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[13]); 
  alphaDrSurf[22] = -0.5*(2.828427124746191*nuUSum[14]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[14]); 
  alphaDrSurf[23] = -0.5*(2.828427124746191*nuUSum[15]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[15]); 
  alphaDrSurf[24] = -0.5*(2.828427124746191*nuUSum[16]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[16]); 
  alphaDrSurf[32] = -0.5*(2.828427124746191*nuUSum[17]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[17]); 
  alphaDrSurf[33] = -0.5*(2.828427124746191*nuUSum[18]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[18]); 
  alphaDrSurf[34] = -0.5*(2.828427124746191*nuUSum[19]+((-2.828427124746191*w[3])-1.414213562373095*dxv[3])*nuSum[19]); 

  if (0.4024922359499623*alphaDrSurf[34]+0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx4_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_5x_p2_surfx4_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = ser_5x_p2_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx4_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_5x_p2_surfx4_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = ser_5x_p2_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]+0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_5x_p2_surfx4_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = ser_5x_p2_surfx4_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = ser_5x_p2_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_5x_p2_surfx4_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = ser_5x_p2_surfx4_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = ser_5x_p2_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]+0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_5x_p2_surfx4_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = ser_5x_p2_surfx4_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = ser_5x_p2_surfx4_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_5x_p2_surfx4_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = ser_5x_p2_surfx4_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = ser_5x_p2_surfx4_eval_quad_node_8_l(fEdge); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[34])+0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_5x_p2_surfx4_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = ser_5x_p2_surfx4_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = ser_5x_p2_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_5x_p2_surfx4_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = ser_5x_p2_surfx4_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = ser_5x_p2_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[34])+0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_5x_p2_surfx4_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = ser_5x_p2_surfx4_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = ser_5x_p2_surfx4_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_5x_p2_surfx4_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = ser_5x_p2_surfx4_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = ser_5x_p2_surfx4_eval_quad_node_14_l(fEdge); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[34])+0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_5x_p2_surfx4_eval_quad_node_15_r(fSkin); 
    fUpwindQuad[16] = ser_5x_p2_surfx4_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = ser_5x_p2_surfx4_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_5x_p2_surfx4_eval_quad_node_15_l(fEdge); 
    fUpwindQuad[16] = ser_5x_p2_surfx4_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = ser_5x_p2_surfx4_eval_quad_node_17_l(fEdge); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]-0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_5x_p2_surfx4_eval_quad_node_18_r(fSkin); 
    fUpwindQuad[19] = ser_5x_p2_surfx4_eval_quad_node_19_r(fSkin); 
    fUpwindQuad[20] = ser_5x_p2_surfx4_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_5x_p2_surfx4_eval_quad_node_18_l(fEdge); 
    fUpwindQuad[19] = ser_5x_p2_surfx4_eval_quad_node_19_l(fEdge); 
    fUpwindQuad[20] = ser_5x_p2_surfx4_eval_quad_node_20_l(fEdge); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]-0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_5x_p2_surfx4_eval_quad_node_21_r(fSkin); 
    fUpwindQuad[22] = ser_5x_p2_surfx4_eval_quad_node_22_r(fSkin); 
    fUpwindQuad[23] = ser_5x_p2_surfx4_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[21] = ser_5x_p2_surfx4_eval_quad_node_21_l(fEdge); 
    fUpwindQuad[22] = ser_5x_p2_surfx4_eval_quad_node_22_l(fEdge); 
    fUpwindQuad[23] = ser_5x_p2_surfx4_eval_quad_node_23_l(fEdge); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]-0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_5x_p2_surfx4_eval_quad_node_24_r(fSkin); 
    fUpwindQuad[25] = ser_5x_p2_surfx4_eval_quad_node_25_r(fSkin); 
    fUpwindQuad[26] = ser_5x_p2_surfx4_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_5x_p2_surfx4_eval_quad_node_24_l(fEdge); 
    fUpwindQuad[25] = ser_5x_p2_surfx4_eval_quad_node_25_l(fEdge); 
    fUpwindQuad[26] = ser_5x_p2_surfx4_eval_quad_node_26_l(fEdge); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[33])-0.3*alphaDrSurf[23]+0.375*alphaDrSurf[22]-0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[6]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[27] = ser_5x_p2_surfx4_eval_quad_node_27_r(fSkin); 
    fUpwindQuad[28] = ser_5x_p2_surfx4_eval_quad_node_28_r(fSkin); 
    fUpwindQuad[29] = ser_5x_p2_surfx4_eval_quad_node_29_r(fSkin); 
  } else { 
    fUpwindQuad[27] = ser_5x_p2_surfx4_eval_quad_node_27_l(fEdge); 
    fUpwindQuad[28] = ser_5x_p2_surfx4_eval_quad_node_28_l(fEdge); 
    fUpwindQuad[29] = ser_5x_p2_surfx4_eval_quad_node_29_l(fEdge); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[33])-0.3*alphaDrSurf[23]+0.375*alphaDrSurf[22]-0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[6]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[30] = ser_5x_p2_surfx4_eval_quad_node_30_r(fSkin); 
    fUpwindQuad[31] = ser_5x_p2_surfx4_eval_quad_node_31_r(fSkin); 
    fUpwindQuad[32] = ser_5x_p2_surfx4_eval_quad_node_32_r(fSkin); 
  } else { 
    fUpwindQuad[30] = ser_5x_p2_surfx4_eval_quad_node_30_l(fEdge); 
    fUpwindQuad[31] = ser_5x_p2_surfx4_eval_quad_node_31_l(fEdge); 
    fUpwindQuad[32] = ser_5x_p2_surfx4_eval_quad_node_32_l(fEdge); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[33])-0.3*alphaDrSurf[23]+0.375*alphaDrSurf[22]-0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[6]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[33] = ser_5x_p2_surfx4_eval_quad_node_33_r(fSkin); 
    fUpwindQuad[34] = ser_5x_p2_surfx4_eval_quad_node_34_r(fSkin); 
    fUpwindQuad[35] = ser_5x_p2_surfx4_eval_quad_node_35_r(fSkin); 
  } else { 
    fUpwindQuad[33] = ser_5x_p2_surfx4_eval_quad_node_33_l(fEdge); 
    fUpwindQuad[34] = ser_5x_p2_surfx4_eval_quad_node_34_l(fEdge); 
    fUpwindQuad[35] = ser_5x_p2_surfx4_eval_quad_node_35_l(fEdge); 
  } 
  if (0.375*alphaDrSurf[23]+0.375*alphaDrSurf[20]-0.2795084971874732*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[36] = ser_5x_p2_surfx4_eval_quad_node_36_r(fSkin); 
    fUpwindQuad[37] = ser_5x_p2_surfx4_eval_quad_node_37_r(fSkin); 
    fUpwindQuad[38] = ser_5x_p2_surfx4_eval_quad_node_38_r(fSkin); 
  } else { 
    fUpwindQuad[36] = ser_5x_p2_surfx4_eval_quad_node_36_l(fEdge); 
    fUpwindQuad[37] = ser_5x_p2_surfx4_eval_quad_node_37_l(fEdge); 
    fUpwindQuad[38] = ser_5x_p2_surfx4_eval_quad_node_38_l(fEdge); 
  } 
  if (0.375*alphaDrSurf[23]+0.375*alphaDrSurf[20]-0.2795084971874732*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[39] = ser_5x_p2_surfx4_eval_quad_node_39_r(fSkin); 
    fUpwindQuad[40] = ser_5x_p2_surfx4_eval_quad_node_40_r(fSkin); 
    fUpwindQuad[41] = ser_5x_p2_surfx4_eval_quad_node_41_r(fSkin); 
  } else { 
    fUpwindQuad[39] = ser_5x_p2_surfx4_eval_quad_node_39_l(fEdge); 
    fUpwindQuad[40] = ser_5x_p2_surfx4_eval_quad_node_40_l(fEdge); 
    fUpwindQuad[41] = ser_5x_p2_surfx4_eval_quad_node_41_l(fEdge); 
  } 
  if (0.375*alphaDrSurf[23]+0.375*alphaDrSurf[20]-0.2795084971874732*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[42] = ser_5x_p2_surfx4_eval_quad_node_42_r(fSkin); 
    fUpwindQuad[43] = ser_5x_p2_surfx4_eval_quad_node_43_r(fSkin); 
    fUpwindQuad[44] = ser_5x_p2_surfx4_eval_quad_node_44_r(fSkin); 
  } else { 
    fUpwindQuad[42] = ser_5x_p2_surfx4_eval_quad_node_42_l(fEdge); 
    fUpwindQuad[43] = ser_5x_p2_surfx4_eval_quad_node_43_l(fEdge); 
    fUpwindQuad[44] = ser_5x_p2_surfx4_eval_quad_node_44_l(fEdge); 
  } 
  if (0.5031152949374518*alphaDrSurf[33]-0.3*alphaDrSurf[23]-0.375*alphaDrSurf[22]+0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[6]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[45] = ser_5x_p2_surfx4_eval_quad_node_45_r(fSkin); 
    fUpwindQuad[46] = ser_5x_p2_surfx4_eval_quad_node_46_r(fSkin); 
    fUpwindQuad[47] = ser_5x_p2_surfx4_eval_quad_node_47_r(fSkin); 
  } else { 
    fUpwindQuad[45] = ser_5x_p2_surfx4_eval_quad_node_45_l(fEdge); 
    fUpwindQuad[46] = ser_5x_p2_surfx4_eval_quad_node_46_l(fEdge); 
    fUpwindQuad[47] = ser_5x_p2_surfx4_eval_quad_node_47_l(fEdge); 
  } 
  if (0.5031152949374518*alphaDrSurf[33]-0.3*alphaDrSurf[23]-0.375*alphaDrSurf[22]+0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[6]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[48] = ser_5x_p2_surfx4_eval_quad_node_48_r(fSkin); 
    fUpwindQuad[49] = ser_5x_p2_surfx4_eval_quad_node_49_r(fSkin); 
    fUpwindQuad[50] = ser_5x_p2_surfx4_eval_quad_node_50_r(fSkin); 
  } else { 
    fUpwindQuad[48] = ser_5x_p2_surfx4_eval_quad_node_48_l(fEdge); 
    fUpwindQuad[49] = ser_5x_p2_surfx4_eval_quad_node_49_l(fEdge); 
    fUpwindQuad[50] = ser_5x_p2_surfx4_eval_quad_node_50_l(fEdge); 
  } 
  if (0.5031152949374518*alphaDrSurf[33]-0.3*alphaDrSurf[23]-0.375*alphaDrSurf[22]+0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[6]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[51] = ser_5x_p2_surfx4_eval_quad_node_51_r(fSkin); 
    fUpwindQuad[52] = ser_5x_p2_surfx4_eval_quad_node_52_r(fSkin); 
    fUpwindQuad[53] = ser_5x_p2_surfx4_eval_quad_node_53_r(fSkin); 
  } else { 
    fUpwindQuad[51] = ser_5x_p2_surfx4_eval_quad_node_51_l(fEdge); 
    fUpwindQuad[52] = ser_5x_p2_surfx4_eval_quad_node_52_l(fEdge); 
    fUpwindQuad[53] = ser_5x_p2_surfx4_eval_quad_node_53_l(fEdge); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])+0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[54] = ser_5x_p2_surfx4_eval_quad_node_54_r(fSkin); 
    fUpwindQuad[55] = ser_5x_p2_surfx4_eval_quad_node_55_r(fSkin); 
    fUpwindQuad[56] = ser_5x_p2_surfx4_eval_quad_node_56_r(fSkin); 
  } else { 
    fUpwindQuad[54] = ser_5x_p2_surfx4_eval_quad_node_54_l(fEdge); 
    fUpwindQuad[55] = ser_5x_p2_surfx4_eval_quad_node_55_l(fEdge); 
    fUpwindQuad[56] = ser_5x_p2_surfx4_eval_quad_node_56_l(fEdge); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])+0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[57] = ser_5x_p2_surfx4_eval_quad_node_57_r(fSkin); 
    fUpwindQuad[58] = ser_5x_p2_surfx4_eval_quad_node_58_r(fSkin); 
    fUpwindQuad[59] = ser_5x_p2_surfx4_eval_quad_node_59_r(fSkin); 
  } else { 
    fUpwindQuad[57] = ser_5x_p2_surfx4_eval_quad_node_57_l(fEdge); 
    fUpwindQuad[58] = ser_5x_p2_surfx4_eval_quad_node_58_l(fEdge); 
    fUpwindQuad[59] = ser_5x_p2_surfx4_eval_quad_node_59_l(fEdge); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])+0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[60] = ser_5x_p2_surfx4_eval_quad_node_60_r(fSkin); 
    fUpwindQuad[61] = ser_5x_p2_surfx4_eval_quad_node_61_r(fSkin); 
    fUpwindQuad[62] = ser_5x_p2_surfx4_eval_quad_node_62_r(fSkin); 
  } else { 
    fUpwindQuad[60] = ser_5x_p2_surfx4_eval_quad_node_60_l(fEdge); 
    fUpwindQuad[61] = ser_5x_p2_surfx4_eval_quad_node_61_l(fEdge); 
    fUpwindQuad[62] = ser_5x_p2_surfx4_eval_quad_node_62_l(fEdge); 
  } 
  if (0.5031152949374518*alphaDrSurf[34]-0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[63] = ser_5x_p2_surfx4_eval_quad_node_63_r(fSkin); 
    fUpwindQuad[64] = ser_5x_p2_surfx4_eval_quad_node_64_r(fSkin); 
    fUpwindQuad[65] = ser_5x_p2_surfx4_eval_quad_node_65_r(fSkin); 
  } else { 
    fUpwindQuad[63] = ser_5x_p2_surfx4_eval_quad_node_63_l(fEdge); 
    fUpwindQuad[64] = ser_5x_p2_surfx4_eval_quad_node_64_l(fEdge); 
    fUpwindQuad[65] = ser_5x_p2_surfx4_eval_quad_node_65_l(fEdge); 
  } 
  if (0.5031152949374518*alphaDrSurf[34]-0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[66] = ser_5x_p2_surfx4_eval_quad_node_66_r(fSkin); 
    fUpwindQuad[67] = ser_5x_p2_surfx4_eval_quad_node_67_r(fSkin); 
    fUpwindQuad[68] = ser_5x_p2_surfx4_eval_quad_node_68_r(fSkin); 
  } else { 
    fUpwindQuad[66] = ser_5x_p2_surfx4_eval_quad_node_66_l(fEdge); 
    fUpwindQuad[67] = ser_5x_p2_surfx4_eval_quad_node_67_l(fEdge); 
    fUpwindQuad[68] = ser_5x_p2_surfx4_eval_quad_node_68_l(fEdge); 
  } 
  if (0.5031152949374518*alphaDrSurf[34]-0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[69] = ser_5x_p2_surfx4_eval_quad_node_69_r(fSkin); 
    fUpwindQuad[70] = ser_5x_p2_surfx4_eval_quad_node_70_r(fSkin); 
    fUpwindQuad[71] = ser_5x_p2_surfx4_eval_quad_node_71_r(fSkin); 
  } else { 
    fUpwindQuad[69] = ser_5x_p2_surfx4_eval_quad_node_69_l(fEdge); 
    fUpwindQuad[70] = ser_5x_p2_surfx4_eval_quad_node_70_l(fEdge); 
    fUpwindQuad[71] = ser_5x_p2_surfx4_eval_quad_node_71_l(fEdge); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])-0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[72] = ser_5x_p2_surfx4_eval_quad_node_72_r(fSkin); 
    fUpwindQuad[73] = ser_5x_p2_surfx4_eval_quad_node_73_r(fSkin); 
    fUpwindQuad[74] = ser_5x_p2_surfx4_eval_quad_node_74_r(fSkin); 
  } else { 
    fUpwindQuad[72] = ser_5x_p2_surfx4_eval_quad_node_72_l(fEdge); 
    fUpwindQuad[73] = ser_5x_p2_surfx4_eval_quad_node_73_l(fEdge); 
    fUpwindQuad[74] = ser_5x_p2_surfx4_eval_quad_node_74_l(fEdge); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])-0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[75] = ser_5x_p2_surfx4_eval_quad_node_75_r(fSkin); 
    fUpwindQuad[76] = ser_5x_p2_surfx4_eval_quad_node_76_r(fSkin); 
    fUpwindQuad[77] = ser_5x_p2_surfx4_eval_quad_node_77_r(fSkin); 
  } else { 
    fUpwindQuad[75] = ser_5x_p2_surfx4_eval_quad_node_75_l(fEdge); 
    fUpwindQuad[76] = ser_5x_p2_surfx4_eval_quad_node_76_l(fEdge); 
    fUpwindQuad[77] = ser_5x_p2_surfx4_eval_quad_node_77_l(fEdge); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])-0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[78] = ser_5x_p2_surfx4_eval_quad_node_78_r(fSkin); 
    fUpwindQuad[79] = ser_5x_p2_surfx4_eval_quad_node_79_r(fSkin); 
    fUpwindQuad[80] = ser_5x_p2_surfx4_eval_quad_node_80_r(fSkin); 
  } else { 
    fUpwindQuad[78] = ser_5x_p2_surfx4_eval_quad_node_78_l(fEdge); 
    fUpwindQuad[79] = ser_5x_p2_surfx4_eval_quad_node_79_l(fEdge); 
    fUpwindQuad[80] = ser_5x_p2_surfx4_eval_quad_node_80_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alphaDrSurf[34]*fUpwind[34]+0.25*alphaDrSurf[33]*fUpwind[33]+0.25*alphaDrSurf[32]*fUpwind[32]+0.25*alphaDrSurf[24]*fUpwind[24]+0.25*alphaDrSurf[23]*fUpwind[23]+0.25*alphaDrSurf[22]*fUpwind[22]+0.25*alphaDrSurf[21]*fUpwind[21]+0.25*alphaDrSurf[20]*fUpwind[20]+0.25*alphaDrSurf[19]*fUpwind[19]+0.25*alphaDrSurf[15]*fUpwind[15]+0.25*alphaDrSurf[13]*fUpwind[13]+0.25*alphaDrSurf[12]*fUpwind[12]+0.25*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[7]*fUpwind[7]+0.25*alphaDrSurf[6]*fUpwind[6]+0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[3]*fUpwind[3]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.2500000000000001*alphaDrSurf[24]*fUpwind[34]+0.2500000000000001*fUpwind[24]*alphaDrSurf[34]+0.2500000000000001*alphaDrSurf[22]*fUpwind[33]+0.2500000000000001*fUpwind[22]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[15]*fUpwind[32]+0.223606797749979*fUpwind[15]*alphaDrSurf[32]+0.2500000000000001*alphaDrSurf[13]*fUpwind[23]+0.2500000000000001*fUpwind[13]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[6]*fUpwind[21]+0.223606797749979*fUpwind[6]*alphaDrSurf[21]+0.2500000000000001*alphaDrSurf[12]*fUpwind[20]+0.2500000000000001*fUpwind[12]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[19]+0.223606797749979*fUpwind[5]*alphaDrSurf[19]+0.25*alphaDrSurf[7]*fUpwind[15]+0.25*fUpwind[7]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[1]*fUpwind[11]+0.223606797749979*fUpwind[1]*alphaDrSurf[11]+0.25*alphaDrSurf[3]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[6]+0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.2500000000000001*alphaDrSurf[23]*fUpwind[34]+0.2500000000000001*fUpwind[23]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[15]*fUpwind[33]+0.223606797749979*fUpwind[15]*alphaDrSurf[33]+0.2500000000000001*alphaDrSurf[21]*fUpwind[32]+0.2500000000000001*fUpwind[21]*alphaDrSurf[32]+0.2500000000000001*alphaDrSurf[13]*fUpwind[24]+0.2500000000000001*fUpwind[13]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[7]*fUpwind[22]+0.223606797749979*fUpwind[7]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[5]*fUpwind[20]+0.223606797749979*fUpwind[5]*alphaDrSurf[20]+0.2500000000000001*alphaDrSurf[11]*fUpwind[19]+0.2500000000000001*fUpwind[11]*alphaDrSurf[19]+0.25*alphaDrSurf[6]*fUpwind[15]+0.25*fUpwind[6]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[2]*fUpwind[12]+0.223606797749979*fUpwind[2]*alphaDrSurf[12]+0.25*alphaDrSurf[3]*fUpwind[7]+0.25*fUpwind[3]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.223606797749979*alphaDrSurf[15]*fUpwind[34]+0.223606797749979*fUpwind[15]*alphaDrSurf[34]+0.2500000000000001*alphaDrSurf[20]*fUpwind[33]+0.2500000000000001*fUpwind[20]*alphaDrSurf[33]+0.2500000000000001*alphaDrSurf[19]*fUpwind[32]+0.2500000000000001*fUpwind[19]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[7]*fUpwind[24]+0.223606797749979*fUpwind[7]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[6]*fUpwind[23]+0.223606797749979*fUpwind[6]*alphaDrSurf[23]+0.2500000000000001*alphaDrSurf[12]*fUpwind[22]+0.2500000000000001*fUpwind[12]*alphaDrSurf[22]+0.2500000000000001*alphaDrSurf[11]*fUpwind[21]+0.2500000000000001*fUpwind[11]*alphaDrSurf[21]+0.25*alphaDrSurf[5]*fUpwind[15]+0.25*fUpwind[5]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[3]*fUpwind[13]+0.223606797749979*fUpwind[3]*alphaDrSurf[13]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*fUpwind[2]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*fUpwind[1]*alphaDrSurf[6]+0.25*alphaDrSurf[0]*fUpwind[3]+0.25*fUpwind[0]*alphaDrSurf[3]; 
  Ghat[4] = 0.2500000000000001*alphaDrSurf[34]*fUpwind[46]+0.2500000000000001*alphaDrSurf[33]*fUpwind[45]+0.2500000000000001*alphaDrSurf[32]*fUpwind[44]+0.2500000000000001*alphaDrSurf[24]*fUpwind[40]+0.2500000000000001*alphaDrSurf[23]*fUpwind[39]+0.2500000000000001*alphaDrSurf[22]*fUpwind[38]+0.2500000000000001*alphaDrSurf[21]*fUpwind[37]+0.2500000000000001*alphaDrSurf[20]*fUpwind[36]+0.2500000000000001*alphaDrSurf[19]*fUpwind[35]+0.25*alphaDrSurf[15]*fUpwind[31]+0.2500000000000001*alphaDrSurf[13]*fUpwind[27]+0.2500000000000001*alphaDrSurf[12]*fUpwind[26]+0.2500000000000001*alphaDrSurf[11]*fUpwind[25]+0.25*alphaDrSurf[7]*fUpwind[18]+0.25*alphaDrSurf[6]*fUpwind[17]+0.25*alphaDrSurf[5]*fUpwind[16]+0.25*alphaDrSurf[3]*fUpwind[10]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  Ghat[5] = 0.25*alphaDrSurf[13]*fUpwind[34]+0.25*fUpwind[13]*alphaDrSurf[34]+0.2*alphaDrSurf[32]*fUpwind[33]+0.223606797749979*alphaDrSurf[7]*fUpwind[33]+0.2*fUpwind[32]*alphaDrSurf[33]+0.223606797749979*fUpwind[7]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[6]*fUpwind[32]+0.223606797749979*fUpwind[6]*alphaDrSurf[32]+0.25*alphaDrSurf[23]*fUpwind[24]+0.25*fUpwind[23]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[15]*fUpwind[22]+0.223606797749979*fUpwind[15]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[15]*fUpwind[21]+0.223606797749979*fUpwind[15]*alphaDrSurf[21]+0.2*alphaDrSurf[19]*fUpwind[20]+0.223606797749979*alphaDrSurf[2]*fUpwind[20]+0.2*fUpwind[19]*alphaDrSurf[20]+0.223606797749979*fUpwind[2]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[1]*fUpwind[19]+0.223606797749979*fUpwind[1]*alphaDrSurf[19]+0.25*alphaDrSurf[3]*fUpwind[15]+0.25*fUpwind[3]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[5]*fUpwind[12]+0.223606797749979*fUpwind[5]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[11]+0.223606797749979*fUpwind[5]*alphaDrSurf[11]+0.25*alphaDrSurf[6]*fUpwind[7]+0.25*fUpwind[6]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[6] = 0.2*alphaDrSurf[32]*fUpwind[34]+0.223606797749979*alphaDrSurf[7]*fUpwind[34]+0.2*fUpwind[32]*alphaDrSurf[34]+0.223606797749979*fUpwind[7]*alphaDrSurf[34]+0.25*alphaDrSurf[12]*fUpwind[33]+0.25*fUpwind[12]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[5]*fUpwind[32]+0.223606797749979*fUpwind[5]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[15]*fUpwind[24]+0.223606797749979*fUpwind[15]*alphaDrSurf[24]+0.2*alphaDrSurf[21]*fUpwind[23]+0.223606797749979*alphaDrSurf[3]*fUpwind[23]+0.2*fUpwind[21]*alphaDrSurf[23]+0.223606797749979*fUpwind[3]*alphaDrSurf[23]+0.25*alphaDrSurf[20]*fUpwind[22]+0.25*fUpwind[20]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[1]*fUpwind[21]+0.223606797749979*fUpwind[1]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[15]*fUpwind[19]+0.223606797749979*fUpwind[15]*alphaDrSurf[19]+0.25*alphaDrSurf[2]*fUpwind[15]+0.25*fUpwind[2]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[6]*fUpwind[13]+0.223606797749979*fUpwind[6]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[6]*fUpwind[11]+0.223606797749979*fUpwind[6]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*fUpwind[5]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*fUpwind[0]*alphaDrSurf[6]+0.25*alphaDrSurf[1]*fUpwind[3]+0.25*fUpwind[1]*alphaDrSurf[3]; 
  Ghat[7] = 0.2*alphaDrSurf[33]*fUpwind[34]+0.223606797749979*alphaDrSurf[6]*fUpwind[34]+0.2*fUpwind[33]*alphaDrSurf[34]+0.223606797749979*fUpwind[6]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[5]*fUpwind[33]+0.223606797749979*fUpwind[5]*alphaDrSurf[33]+0.25*alphaDrSurf[11]*fUpwind[32]+0.25*fUpwind[11]*alphaDrSurf[32]+0.2*alphaDrSurf[22]*fUpwind[24]+0.223606797749979*alphaDrSurf[3]*fUpwind[24]+0.2*fUpwind[22]*alphaDrSurf[24]+0.223606797749979*fUpwind[3]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[15]*fUpwind[23]+0.223606797749979*fUpwind[15]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[2]*fUpwind[22]+0.223606797749979*fUpwind[2]*alphaDrSurf[22]+0.25*alphaDrSurf[19]*fUpwind[21]+0.25*fUpwind[19]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[15]*fUpwind[20]+0.223606797749979*fUpwind[15]*alphaDrSurf[20]+0.25*alphaDrSurf[1]*fUpwind[15]+0.25*fUpwind[1]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[7]*fUpwind[13]+0.223606797749979*fUpwind[7]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[7]*fUpwind[12]+0.223606797749979*fUpwind[7]*alphaDrSurf[12]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*fUpwind[0]*alphaDrSurf[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*fUpwind[5]*alphaDrSurf[6]+0.25*alphaDrSurf[2]*fUpwind[3]+0.25*fUpwind[2]*alphaDrSurf[3]; 
  Ghat[8] = 0.25*alphaDrSurf[24]*fUpwind[46]+0.25*alphaDrSurf[22]*fUpwind[45]+0.223606797749979*alphaDrSurf[15]*fUpwind[44]+0.25*alphaDrSurf[34]*fUpwind[40]+0.25*alphaDrSurf[13]*fUpwind[39]+0.25*alphaDrSurf[33]*fUpwind[38]+0.223606797749979*alphaDrSurf[6]*fUpwind[37]+0.25*alphaDrSurf[12]*fUpwind[36]+0.223606797749979*alphaDrSurf[5]*fUpwind[35]+0.223606797749979*fUpwind[31]*alphaDrSurf[32]+0.25*alphaDrSurf[7]*fUpwind[31]+0.25*alphaDrSurf[23]*fUpwind[27]+0.25*alphaDrSurf[20]*fUpwind[26]+0.223606797749979*alphaDrSurf[1]*fUpwind[25]+0.223606797749979*fUpwind[17]*alphaDrSurf[21]+0.223606797749979*fUpwind[16]*alphaDrSurf[19]+0.25*alphaDrSurf[15]*fUpwind[18]+0.25*alphaDrSurf[3]*fUpwind[17]+0.25*alphaDrSurf[2]*fUpwind[16]+0.223606797749979*fUpwind[8]*alphaDrSurf[11]+0.25*alphaDrSurf[6]*fUpwind[10]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  Ghat[9] = 0.25*alphaDrSurf[23]*fUpwind[46]+0.223606797749979*alphaDrSurf[15]*fUpwind[45]+0.25*alphaDrSurf[21]*fUpwind[44]+0.25*alphaDrSurf[13]*fUpwind[40]+0.25*alphaDrSurf[34]*fUpwind[39]+0.223606797749979*alphaDrSurf[7]*fUpwind[38]+0.25*alphaDrSurf[32]*fUpwind[37]+0.223606797749979*alphaDrSurf[5]*fUpwind[36]+0.25*alphaDrSurf[11]*fUpwind[35]+0.223606797749979*fUpwind[31]*alphaDrSurf[33]+0.25*alphaDrSurf[6]*fUpwind[31]+0.25*alphaDrSurf[24]*fUpwind[27]+0.223606797749979*alphaDrSurf[2]*fUpwind[26]+0.25*alphaDrSurf[19]*fUpwind[25]+0.223606797749979*fUpwind[18]*alphaDrSurf[22]+0.223606797749979*fUpwind[16]*alphaDrSurf[20]+0.25*alphaDrSurf[3]*fUpwind[18]+0.25*alphaDrSurf[15]*fUpwind[17]+0.25*alphaDrSurf[1]*fUpwind[16]+0.223606797749979*fUpwind[9]*alphaDrSurf[12]+0.25*alphaDrSurf[7]*fUpwind[10]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  Ghat[10] = 0.223606797749979*alphaDrSurf[15]*fUpwind[46]+0.25*alphaDrSurf[20]*fUpwind[45]+0.25*alphaDrSurf[19]*fUpwind[44]+0.223606797749979*alphaDrSurf[7]*fUpwind[40]+0.223606797749979*alphaDrSurf[6]*fUpwind[39]+0.25*alphaDrSurf[12]*fUpwind[38]+0.25*alphaDrSurf[11]*fUpwind[37]+0.25*alphaDrSurf[33]*fUpwind[36]+0.25*alphaDrSurf[32]*fUpwind[35]+0.223606797749979*fUpwind[31]*alphaDrSurf[34]+0.25*alphaDrSurf[5]*fUpwind[31]+0.223606797749979*alphaDrSurf[3]*fUpwind[27]+0.25*alphaDrSurf[22]*fUpwind[26]+0.25*alphaDrSurf[21]*fUpwind[25]+0.223606797749979*fUpwind[18]*alphaDrSurf[24]+0.223606797749979*fUpwind[17]*alphaDrSurf[23]+0.25*alphaDrSurf[2]*fUpwind[18]+0.25*alphaDrSurf[1]*fUpwind[17]+0.25*alphaDrSurf[15]*fUpwind[16]+0.223606797749979*fUpwind[10]*alphaDrSurf[13]+0.25*alphaDrSurf[0]*fUpwind[10]+0.25*alphaDrSurf[7]*fUpwind[9]+0.25*alphaDrSurf[6]*fUpwind[8]+0.25*alphaDrSurf[3]*fUpwind[4]; 
  Ghat[11] = 0.223606797749979*alphaDrSurf[34]*fUpwind[34]+0.223606797749979*alphaDrSurf[33]*fUpwind[33]+0.159719141249985*alphaDrSurf[32]*fUpwind[32]+0.25*alphaDrSurf[7]*fUpwind[32]+0.25*fUpwind[7]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[23]*fUpwind[23]+0.159719141249985*alphaDrSurf[21]*fUpwind[21]+0.2500000000000001*alphaDrSurf[3]*fUpwind[21]+0.2500000000000001*fUpwind[3]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[20]*fUpwind[20]+0.159719141249985*alphaDrSurf[19]*fUpwind[19]+0.2500000000000001*alphaDrSurf[2]*fUpwind[19]+0.2500000000000001*fUpwind[2]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[15]*fUpwind[15]+0.159719141249985*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[0]*fUpwind[11]+0.25*fUpwind[0]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[6]*fUpwind[6]+0.223606797749979*alphaDrSurf[5]*fUpwind[5]+0.223606797749979*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[12] = 0.223606797749979*alphaDrSurf[34]*fUpwind[34]+0.159719141249985*alphaDrSurf[33]*fUpwind[33]+0.25*alphaDrSurf[6]*fUpwind[33]+0.25*fUpwind[6]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[32]*fUpwind[32]+0.223606797749979*alphaDrSurf[24]*fUpwind[24]+0.159719141249985*alphaDrSurf[22]*fUpwind[22]+0.2500000000000001*alphaDrSurf[3]*fUpwind[22]+0.2500000000000001*fUpwind[3]*alphaDrSurf[22]+0.159719141249985*alphaDrSurf[20]*fUpwind[20]+0.2500000000000001*alphaDrSurf[1]*fUpwind[20]+0.2500000000000001*fUpwind[1]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[19]*fUpwind[19]+0.223606797749979*alphaDrSurf[15]*fUpwind[15]+0.159719141249985*alphaDrSurf[12]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[12]+0.25*fUpwind[0]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[7]*fUpwind[7]+0.223606797749979*alphaDrSurf[5]*fUpwind[5]+0.223606797749979*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[13] = 0.159719141249985*alphaDrSurf[34]*fUpwind[34]+0.25*alphaDrSurf[5]*fUpwind[34]+0.25*fUpwind[5]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[33]*fUpwind[33]+0.223606797749979*alphaDrSurf[32]*fUpwind[32]+0.159719141249985*alphaDrSurf[24]*fUpwind[24]+0.2500000000000001*alphaDrSurf[2]*fUpwind[24]+0.2500000000000001*fUpwind[2]*alphaDrSurf[24]+0.159719141249985*alphaDrSurf[23]*fUpwind[23]+0.2500000000000001*alphaDrSurf[1]*fUpwind[23]+0.2500000000000001*fUpwind[1]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[22]*fUpwind[22]+0.223606797749979*alphaDrSurf[21]*fUpwind[21]+0.223606797749979*alphaDrSurf[15]*fUpwind[15]+0.159719141249985*alphaDrSurf[13]*fUpwind[13]+0.25*alphaDrSurf[0]*fUpwind[13]+0.25*fUpwind[0]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[7]*fUpwind[7]+0.223606797749979*alphaDrSurf[6]*fUpwind[6]+0.223606797749979*alphaDrSurf[3]*fUpwind[3]; 
  Ghat[14] = 0.2500000000000001*alphaDrSurf[15]*fUpwind[47]+0.25*alphaDrSurf[7]*fUpwind[43]+0.25*alphaDrSurf[6]*fUpwind[42]+0.25*alphaDrSurf[5]*fUpwind[41]+0.2500000000000001*alphaDrSurf[3]*fUpwind[30]+0.2500000000000001*alphaDrSurf[2]*fUpwind[29]+0.2500000000000001*alphaDrSurf[1]*fUpwind[28]+0.25*alphaDrSurf[0]*fUpwind[14]; 
  Ghat[15] = 0.2*alphaDrSurf[22]*fUpwind[34]+0.2*alphaDrSurf[21]*fUpwind[34]+0.223606797749979*alphaDrSurf[3]*fUpwind[34]+0.2*fUpwind[22]*alphaDrSurf[34]+0.2*fUpwind[21]*alphaDrSurf[34]+0.223606797749979*fUpwind[3]*alphaDrSurf[34]+0.2*alphaDrSurf[24]*fUpwind[33]+0.2*alphaDrSurf[19]*fUpwind[33]+0.223606797749979*alphaDrSurf[2]*fUpwind[33]+0.2*fUpwind[24]*alphaDrSurf[33]+0.2*fUpwind[19]*alphaDrSurf[33]+0.223606797749979*fUpwind[2]*alphaDrSurf[33]+0.2*alphaDrSurf[23]*fUpwind[32]+0.2*alphaDrSurf[20]*fUpwind[32]+0.223606797749979*alphaDrSurf[1]*fUpwind[32]+0.2*fUpwind[23]*alphaDrSurf[32]+0.2*fUpwind[20]*alphaDrSurf[32]+0.223606797749979*fUpwind[1]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[6]*fUpwind[24]+0.223606797749979*fUpwind[6]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[7]*fUpwind[23]+0.223606797749979*fUpwind[7]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[5]*fUpwind[22]+0.223606797749979*fUpwind[5]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[5]*fUpwind[21]+0.223606797749979*fUpwind[5]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[7]*fUpwind[20]+0.223606797749979*fUpwind[7]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[6]*fUpwind[19]+0.223606797749979*fUpwind[6]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[13]*fUpwind[15]+0.223606797749979*alphaDrSurf[12]*fUpwind[15]+0.223606797749979*alphaDrSurf[11]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[15]+0.223606797749979*fUpwind[13]*alphaDrSurf[15]+0.223606797749979*fUpwind[12]*alphaDrSurf[15]+0.223606797749979*fUpwind[11]*alphaDrSurf[15]+0.25*fUpwind[0]*alphaDrSurf[15]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*fUpwind[1]*alphaDrSurf[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[2]*alphaDrSurf[6]+0.25*alphaDrSurf[3]*fUpwind[5]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  Ghat[16] = 0.2500000000000001*alphaDrSurf[13]*fUpwind[46]+0.2*alphaDrSurf[32]*fUpwind[45]+0.223606797749979*alphaDrSurf[7]*fUpwind[45]+0.2*alphaDrSurf[33]*fUpwind[44]+0.223606797749979*alphaDrSurf[6]*fUpwind[44]+0.2500000000000001*alphaDrSurf[23]*fUpwind[40]+0.2500000000000001*alphaDrSurf[24]*fUpwind[39]+0.223606797749979*alphaDrSurf[15]*fUpwind[38]+0.223606797749979*alphaDrSurf[15]*fUpwind[37]+0.2*alphaDrSurf[19]*fUpwind[36]+0.223606797749979*alphaDrSurf[2]*fUpwind[36]+0.2*alphaDrSurf[20]*fUpwind[35]+0.223606797749979*alphaDrSurf[1]*fUpwind[35]+0.2500000000000001*fUpwind[27]*alphaDrSurf[34]+0.223606797749979*fUpwind[18]*alphaDrSurf[33]+0.223606797749979*fUpwind[17]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[22]*fUpwind[31]+0.223606797749979*alphaDrSurf[21]*fUpwind[31]+0.25*alphaDrSurf[3]*fUpwind[31]+0.223606797749979*alphaDrSurf[5]*fUpwind[26]+0.223606797749979*alphaDrSurf[5]*fUpwind[25]+0.223606797749979*fUpwind[9]*alphaDrSurf[20]+0.223606797749979*fUpwind[8]*alphaDrSurf[19]+0.25*alphaDrSurf[6]*fUpwind[18]+0.25*alphaDrSurf[7]*fUpwind[17]+0.223606797749979*alphaDrSurf[12]*fUpwind[16]+0.223606797749979*alphaDrSurf[11]*fUpwind[16]+0.25*alphaDrSurf[0]*fUpwind[16]+0.25*fUpwind[10]*alphaDrSurf[15]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  Ghat[17] = 0.2*alphaDrSurf[32]*fUpwind[46]+0.223606797749979*alphaDrSurf[7]*fUpwind[46]+0.2500000000000001*alphaDrSurf[12]*fUpwind[45]+0.2*alphaDrSurf[34]*fUpwind[44]+0.223606797749979*alphaDrSurf[5]*fUpwind[44]+0.223606797749979*alphaDrSurf[15]*fUpwind[40]+0.2*alphaDrSurf[21]*fUpwind[39]+0.223606797749979*alphaDrSurf[3]*fUpwind[39]+0.2500000000000001*alphaDrSurf[20]*fUpwind[38]+0.2*alphaDrSurf[23]*fUpwind[37]+0.223606797749979*alphaDrSurf[1]*fUpwind[37]+0.2500000000000001*alphaDrSurf[22]*fUpwind[36]+0.223606797749979*alphaDrSurf[15]*fUpwind[35]+0.223606797749979*fUpwind[18]*alphaDrSurf[34]+0.2500000000000001*fUpwind[26]*alphaDrSurf[33]+0.223606797749979*fUpwind[16]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[24]*fUpwind[31]+0.223606797749979*alphaDrSurf[19]*fUpwind[31]+0.25*alphaDrSurf[2]*fUpwind[31]+0.223606797749979*alphaDrSurf[6]*fUpwind[27]+0.223606797749979*alphaDrSurf[6]*fUpwind[25]+0.223606797749979*fUpwind[10]*alphaDrSurf[23]+0.223606797749979*fUpwind[8]*alphaDrSurf[21]+0.25*alphaDrSurf[5]*fUpwind[18]+0.223606797749979*alphaDrSurf[13]*fUpwind[17]+0.223606797749979*alphaDrSurf[11]*fUpwind[17]+0.25*alphaDrSurf[0]*fUpwind[17]+0.25*alphaDrSurf[7]*fUpwind[16]+0.25*fUpwind[9]*alphaDrSurf[15]+0.25*alphaDrSurf[1]*fUpwind[10]+0.25*alphaDrSurf[3]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[6]; 
  Ghat[18] = 0.2*alphaDrSurf[33]*fUpwind[46]+0.223606797749979*alphaDrSurf[6]*fUpwind[46]+0.2*alphaDrSurf[34]*fUpwind[45]+0.223606797749979*alphaDrSurf[5]*fUpwind[45]+0.2500000000000001*alphaDrSurf[11]*fUpwind[44]+0.2*alphaDrSurf[22]*fUpwind[40]+0.223606797749979*alphaDrSurf[3]*fUpwind[40]+0.223606797749979*alphaDrSurf[15]*fUpwind[39]+0.2*alphaDrSurf[24]*fUpwind[38]+0.223606797749979*alphaDrSurf[2]*fUpwind[38]+0.2500000000000001*alphaDrSurf[19]*fUpwind[37]+0.223606797749979*alphaDrSurf[15]*fUpwind[36]+0.2500000000000001*alphaDrSurf[21]*fUpwind[35]+0.223606797749979*fUpwind[17]*alphaDrSurf[34]+0.223606797749979*fUpwind[16]*alphaDrSurf[33]+0.2500000000000001*fUpwind[25]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[23]*fUpwind[31]+0.223606797749979*alphaDrSurf[20]*fUpwind[31]+0.25*alphaDrSurf[1]*fUpwind[31]+0.223606797749979*alphaDrSurf[7]*fUpwind[27]+0.223606797749979*alphaDrSurf[7]*fUpwind[26]+0.223606797749979*fUpwind[10]*alphaDrSurf[24]+0.223606797749979*fUpwind[9]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[13]*fUpwind[18]+0.223606797749979*alphaDrSurf[12]*fUpwind[18]+0.25*alphaDrSurf[0]*fUpwind[18]+0.25*alphaDrSurf[5]*fUpwind[17]+0.25*alphaDrSurf[6]*fUpwind[16]+0.25*fUpwind[8]*alphaDrSurf[15]+0.25*alphaDrSurf[2]*fUpwind[10]+0.25*alphaDrSurf[3]*fUpwind[9]+0.25*fUpwind[4]*alphaDrSurf[7]; 
  Ghat[19] = 0.223606797749979*alphaDrSurf[23]*fUpwind[34]+0.223606797749979*fUpwind[23]*alphaDrSurf[34]+0.2*alphaDrSurf[15]*fUpwind[33]+0.2*fUpwind[15]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[22]*fUpwind[32]+0.159719141249985*alphaDrSurf[21]*fUpwind[32]+0.2500000000000001*alphaDrSurf[3]*fUpwind[32]+0.223606797749979*fUpwind[22]*alphaDrSurf[32]+0.159719141249985*fUpwind[21]*alphaDrSurf[32]+0.2500000000000001*fUpwind[3]*alphaDrSurf[32]+0.25*alphaDrSurf[7]*fUpwind[21]+0.25*fUpwind[7]*alphaDrSurf[21]+0.2*alphaDrSurf[5]*fUpwind[20]+0.2*fUpwind[5]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[12]*fUpwind[19]+0.159719141249985*alphaDrSurf[11]*fUpwind[19]+0.25*alphaDrSurf[0]*fUpwind[19]+0.223606797749979*fUpwind[12]*alphaDrSurf[19]+0.159719141249985*fUpwind[11]*alphaDrSurf[19]+0.25*fUpwind[0]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[6]*fUpwind[15]+0.223606797749979*fUpwind[6]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[2]*fUpwind[11]+0.2500000000000001*fUpwind[2]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[5]+0.223606797749979*fUpwind[1]*alphaDrSurf[5]; 
  Ghat[20] = 0.223606797749979*alphaDrSurf[24]*fUpwind[34]+0.223606797749979*fUpwind[24]*alphaDrSurf[34]+0.159719141249985*alphaDrSurf[22]*fUpwind[33]+0.223606797749979*alphaDrSurf[21]*fUpwind[33]+0.2500000000000001*alphaDrSurf[3]*fUpwind[33]+0.159719141249985*fUpwind[22]*alphaDrSurf[33]+0.223606797749979*fUpwind[21]*alphaDrSurf[33]+0.2500000000000001*fUpwind[3]*alphaDrSurf[33]+0.2*alphaDrSurf[15]*fUpwind[32]+0.2*fUpwind[15]*alphaDrSurf[32]+0.25*alphaDrSurf[6]*fUpwind[22]+0.25*fUpwind[6]*alphaDrSurf[22]+0.159719141249985*alphaDrSurf[12]*fUpwind[20]+0.223606797749979*alphaDrSurf[11]*fUpwind[20]+0.25*alphaDrSurf[0]*fUpwind[20]+0.159719141249985*fUpwind[12]*alphaDrSurf[20]+0.223606797749979*fUpwind[11]*alphaDrSurf[20]+0.25*fUpwind[0]*alphaDrSurf[20]+0.2*alphaDrSurf[5]*fUpwind[19]+0.2*fUpwind[5]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[7]*fUpwind[15]+0.223606797749979*fUpwind[7]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[1]*fUpwind[12]+0.2500000000000001*fUpwind[1]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[5]+0.223606797749979*fUpwind[2]*alphaDrSurf[5]; 
  Ghat[21] = 0.2*alphaDrSurf[15]*fUpwind[34]+0.2*fUpwind[15]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[20]*fUpwind[33]+0.223606797749979*fUpwind[20]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[24]*fUpwind[32]+0.159719141249985*alphaDrSurf[19]*fUpwind[32]+0.2500000000000001*alphaDrSurf[2]*fUpwind[32]+0.223606797749979*fUpwind[24]*alphaDrSurf[32]+0.159719141249985*fUpwind[19]*alphaDrSurf[32]+0.2500000000000001*fUpwind[2]*alphaDrSurf[32]+0.2*alphaDrSurf[6]*fUpwind[23]+0.2*fUpwind[6]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[13]*fUpwind[21]+0.159719141249985*alphaDrSurf[11]*fUpwind[21]+0.25*alphaDrSurf[0]*fUpwind[21]+0.223606797749979*fUpwind[13]*alphaDrSurf[21]+0.159719141249985*fUpwind[11]*alphaDrSurf[21]+0.25*fUpwind[0]*alphaDrSurf[21]+0.25*alphaDrSurf[7]*fUpwind[19]+0.25*fUpwind[7]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[5]*fUpwind[15]+0.223606797749979*fUpwind[5]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[3]*fUpwind[11]+0.2500000000000001*fUpwind[3]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[6]+0.223606797749979*fUpwind[1]*alphaDrSurf[6]; 
  Ghat[22] = 0.2*alphaDrSurf[15]*fUpwind[34]+0.2*fUpwind[15]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[23]*fUpwind[33]+0.159719141249985*alphaDrSurf[20]*fUpwind[33]+0.2500000000000001*alphaDrSurf[1]*fUpwind[33]+0.223606797749979*fUpwind[23]*alphaDrSurf[33]+0.159719141249985*fUpwind[20]*alphaDrSurf[33]+0.2500000000000001*fUpwind[1]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[19]*fUpwind[32]+0.223606797749979*fUpwind[19]*alphaDrSurf[32]+0.2*alphaDrSurf[7]*fUpwind[24]+0.2*fUpwind[7]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[13]*fUpwind[22]+0.159719141249985*alphaDrSurf[12]*fUpwind[22]+0.25*alphaDrSurf[0]*fUpwind[22]+0.223606797749979*fUpwind[13]*alphaDrSurf[22]+0.159719141249985*fUpwind[12]*alphaDrSurf[22]+0.25*fUpwind[0]*alphaDrSurf[22]+0.25*alphaDrSurf[6]*fUpwind[20]+0.25*fUpwind[6]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[15]+0.223606797749979*fUpwind[5]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[3]*fUpwind[12]+0.2500000000000001*fUpwind[3]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[7]+0.223606797749979*fUpwind[2]*alphaDrSurf[7]; 
  Ghat[23] = 0.159719141249985*alphaDrSurf[24]*fUpwind[34]+0.223606797749979*alphaDrSurf[19]*fUpwind[34]+0.2500000000000001*alphaDrSurf[2]*fUpwind[34]+0.159719141249985*fUpwind[24]*alphaDrSurf[34]+0.223606797749979*fUpwind[19]*alphaDrSurf[34]+0.2500000000000001*fUpwind[2]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[22]*fUpwind[33]+0.223606797749979*fUpwind[22]*alphaDrSurf[33]+0.2*alphaDrSurf[15]*fUpwind[32]+0.2*fUpwind[15]*alphaDrSurf[32]+0.25*alphaDrSurf[5]*fUpwind[24]+0.25*fUpwind[5]*alphaDrSurf[24]+0.159719141249985*alphaDrSurf[13]*fUpwind[23]+0.223606797749979*alphaDrSurf[11]*fUpwind[23]+0.25*alphaDrSurf[0]*fUpwind[23]+0.159719141249985*fUpwind[13]*alphaDrSurf[23]+0.223606797749979*fUpwind[11]*alphaDrSurf[23]+0.25*fUpwind[0]*alphaDrSurf[23]+0.2*alphaDrSurf[6]*fUpwind[21]+0.2*fUpwind[6]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[7]*fUpwind[15]+0.223606797749979*fUpwind[7]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[1]*fUpwind[13]+0.2500000000000001*fUpwind[1]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[3]*fUpwind[6]+0.223606797749979*fUpwind[3]*alphaDrSurf[6]; 
  Ghat[24] = 0.159719141249985*alphaDrSurf[23]*fUpwind[34]+0.223606797749979*alphaDrSurf[20]*fUpwind[34]+0.2500000000000001*alphaDrSurf[1]*fUpwind[34]+0.159719141249985*fUpwind[23]*alphaDrSurf[34]+0.223606797749979*fUpwind[20]*alphaDrSurf[34]+0.2500000000000001*fUpwind[1]*alphaDrSurf[34]+0.2*alphaDrSurf[15]*fUpwind[33]+0.2*fUpwind[15]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[21]*fUpwind[32]+0.223606797749979*fUpwind[21]*alphaDrSurf[32]+0.159719141249985*alphaDrSurf[13]*fUpwind[24]+0.223606797749979*alphaDrSurf[12]*fUpwind[24]+0.25*alphaDrSurf[0]*fUpwind[24]+0.159719141249985*fUpwind[13]*alphaDrSurf[24]+0.223606797749979*fUpwind[12]*alphaDrSurf[24]+0.25*fUpwind[0]*alphaDrSurf[24]+0.25*alphaDrSurf[5]*fUpwind[23]+0.25*fUpwind[5]*alphaDrSurf[23]+0.2*alphaDrSurf[7]*fUpwind[22]+0.2*fUpwind[7]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[6]*fUpwind[15]+0.223606797749979*fUpwind[6]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[2]*fUpwind[13]+0.2500000000000001*fUpwind[2]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[3]*fUpwind[7]+0.223606797749979*fUpwind[3]*alphaDrSurf[7]; 
  Ghat[25] = 0.223606797749979*alphaDrSurf[34]*fUpwind[46]+0.223606797749979*alphaDrSurf[33]*fUpwind[45]+0.159719141249985*alphaDrSurf[32]*fUpwind[44]+0.25*alphaDrSurf[7]*fUpwind[44]+0.223606797749979*alphaDrSurf[23]*fUpwind[39]+0.159719141249985*alphaDrSurf[21]*fUpwind[37]+0.2500000000000001*alphaDrSurf[3]*fUpwind[37]+0.223606797749979*alphaDrSurf[20]*fUpwind[36]+0.159719141249985*alphaDrSurf[19]*fUpwind[35]+0.2500000000000001*alphaDrSurf[2]*fUpwind[35]+0.2500000000000001*fUpwind[18]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[15]*fUpwind[31]+0.159719141249985*alphaDrSurf[11]*fUpwind[25]+0.25*alphaDrSurf[0]*fUpwind[25]+0.25*fUpwind[10]*alphaDrSurf[21]+0.25*fUpwind[9]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[6]*fUpwind[17]+0.223606797749979*alphaDrSurf[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[8]; 
  Ghat[26] = 0.223606797749979*alphaDrSurf[34]*fUpwind[46]+0.159719141249985*alphaDrSurf[33]*fUpwind[45]+0.25*alphaDrSurf[6]*fUpwind[45]+0.223606797749979*alphaDrSurf[32]*fUpwind[44]+0.223606797749979*alphaDrSurf[24]*fUpwind[40]+0.159719141249985*alphaDrSurf[22]*fUpwind[38]+0.2500000000000001*alphaDrSurf[3]*fUpwind[38]+0.159719141249985*alphaDrSurf[20]*fUpwind[36]+0.2500000000000001*alphaDrSurf[1]*fUpwind[36]+0.223606797749979*alphaDrSurf[19]*fUpwind[35]+0.2500000000000001*fUpwind[17]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[15]*fUpwind[31]+0.159719141249985*alphaDrSurf[12]*fUpwind[26]+0.25*alphaDrSurf[0]*fUpwind[26]+0.25*fUpwind[10]*alphaDrSurf[22]+0.25*fUpwind[8]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[7]*fUpwind[18]+0.223606797749979*alphaDrSurf[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[9]; 
  Ghat[27] = 0.159719141249985*alphaDrSurf[34]*fUpwind[46]+0.25*alphaDrSurf[5]*fUpwind[46]+0.223606797749979*alphaDrSurf[33]*fUpwind[45]+0.223606797749979*alphaDrSurf[32]*fUpwind[44]+0.159719141249985*alphaDrSurf[24]*fUpwind[40]+0.2500000000000001*alphaDrSurf[2]*fUpwind[40]+0.159719141249985*alphaDrSurf[23]*fUpwind[39]+0.2500000000000001*alphaDrSurf[1]*fUpwind[39]+0.223606797749979*alphaDrSurf[22]*fUpwind[38]+0.223606797749979*alphaDrSurf[21]*fUpwind[37]+0.2500000000000001*fUpwind[16]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[15]*fUpwind[31]+0.159719141249985*alphaDrSurf[13]*fUpwind[27]+0.25*alphaDrSurf[0]*fUpwind[27]+0.25*fUpwind[9]*alphaDrSurf[24]+0.25*fUpwind[8]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[7]*fUpwind[18]+0.223606797749979*alphaDrSurf[6]*fUpwind[17]+0.2500000000000001*fUpwind[4]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[3]*fUpwind[10]; 
  Ghat[28] = 0.223606797749979*alphaDrSurf[32]*fUpwind[47]+0.25*alphaDrSurf[7]*fUpwind[47]+0.2500000000000001*alphaDrSurf[15]*fUpwind[43]+0.223606797749979*alphaDrSurf[21]*fUpwind[42]+0.2500000000000001*alphaDrSurf[3]*fUpwind[42]+0.223606797749979*alphaDrSurf[19]*fUpwind[41]+0.2500000000000001*alphaDrSurf[2]*fUpwind[41]+0.25*alphaDrSurf[6]*fUpwind[30]+0.25*alphaDrSurf[5]*fUpwind[29]+0.223606797749979*alphaDrSurf[11]*fUpwind[28]+0.25*alphaDrSurf[0]*fUpwind[28]+0.2500000000000001*alphaDrSurf[1]*fUpwind[14]; 
  Ghat[29] = 0.223606797749979*alphaDrSurf[33]*fUpwind[47]+0.25*alphaDrSurf[6]*fUpwind[47]+0.223606797749979*alphaDrSurf[22]*fUpwind[43]+0.2500000000000001*alphaDrSurf[3]*fUpwind[43]+0.2500000000000001*alphaDrSurf[15]*fUpwind[42]+0.223606797749979*alphaDrSurf[20]*fUpwind[41]+0.2500000000000001*alphaDrSurf[1]*fUpwind[41]+0.25*alphaDrSurf[7]*fUpwind[30]+0.223606797749979*alphaDrSurf[12]*fUpwind[29]+0.25*alphaDrSurf[0]*fUpwind[29]+0.25*alphaDrSurf[5]*fUpwind[28]+0.2500000000000001*alphaDrSurf[2]*fUpwind[14]; 
  Ghat[30] = 0.223606797749979*alphaDrSurf[34]*fUpwind[47]+0.25*alphaDrSurf[5]*fUpwind[47]+0.223606797749979*alphaDrSurf[24]*fUpwind[43]+0.2500000000000001*alphaDrSurf[2]*fUpwind[43]+0.223606797749979*alphaDrSurf[23]*fUpwind[42]+0.2500000000000001*alphaDrSurf[1]*fUpwind[42]+0.2500000000000001*alphaDrSurf[15]*fUpwind[41]+0.223606797749979*alphaDrSurf[13]*fUpwind[30]+0.25*alphaDrSurf[0]*fUpwind[30]+0.25*alphaDrSurf[7]*fUpwind[29]+0.25*alphaDrSurf[6]*fUpwind[28]+0.2500000000000001*alphaDrSurf[3]*fUpwind[14]; 
  Ghat[31] = 0.2*alphaDrSurf[22]*fUpwind[46]+0.2*alphaDrSurf[21]*fUpwind[46]+0.223606797749979*alphaDrSurf[3]*fUpwind[46]+0.2*alphaDrSurf[24]*fUpwind[45]+0.2*alphaDrSurf[19]*fUpwind[45]+0.223606797749979*alphaDrSurf[2]*fUpwind[45]+0.2*alphaDrSurf[23]*fUpwind[44]+0.2*alphaDrSurf[20]*fUpwind[44]+0.223606797749979*alphaDrSurf[1]*fUpwind[44]+0.2*alphaDrSurf[33]*fUpwind[40]+0.223606797749979*alphaDrSurf[6]*fUpwind[40]+0.2*alphaDrSurf[32]*fUpwind[39]+0.223606797749979*alphaDrSurf[7]*fUpwind[39]+0.2*alphaDrSurf[34]*fUpwind[38]+0.223606797749979*alphaDrSurf[5]*fUpwind[38]+0.2*alphaDrSurf[34]*fUpwind[37]+0.223606797749979*alphaDrSurf[5]*fUpwind[37]+0.2*alphaDrSurf[32]*fUpwind[36]+0.223606797749979*alphaDrSurf[7]*fUpwind[36]+0.2*alphaDrSurf[33]*fUpwind[35]+0.223606797749979*alphaDrSurf[6]*fUpwind[35]+0.223606797749979*fUpwind[10]*alphaDrSurf[34]+0.223606797749979*fUpwind[9]*alphaDrSurf[33]+0.223606797749979*fUpwind[8]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[13]*fUpwind[31]+0.223606797749979*alphaDrSurf[12]*fUpwind[31]+0.223606797749979*alphaDrSurf[11]*fUpwind[31]+0.25*alphaDrSurf[0]*fUpwind[31]+0.223606797749979*alphaDrSurf[15]*fUpwind[27]+0.223606797749979*alphaDrSurf[15]*fUpwind[26]+0.223606797749979*alphaDrSurf[15]*fUpwind[25]+0.223606797749979*fUpwind[17]*alphaDrSurf[24]+0.223606797749979*fUpwind[18]*alphaDrSurf[23]+0.223606797749979*fUpwind[16]*alphaDrSurf[22]+0.223606797749979*fUpwind[16]*alphaDrSurf[21]+0.223606797749979*fUpwind[18]*alphaDrSurf[20]+0.223606797749979*fUpwind[17]*alphaDrSurf[19]+0.25*alphaDrSurf[1]*fUpwind[18]+0.25*alphaDrSurf[2]*fUpwind[17]+0.25*alphaDrSurf[3]*fUpwind[16]+0.25*fUpwind[4]*alphaDrSurf[15]+0.25*alphaDrSurf[5]*fUpwind[10]+0.25*alphaDrSurf[6]*fUpwind[9]+0.25*alphaDrSurf[7]*fUpwind[8]; 
  Ghat[32] = 0.1788854381999831*alphaDrSurf[33]*fUpwind[34]+0.2*alphaDrSurf[6]*fUpwind[34]+0.1788854381999831*fUpwind[33]*alphaDrSurf[34]+0.2*fUpwind[6]*alphaDrSurf[34]+0.2*alphaDrSurf[5]*fUpwind[33]+0.2*fUpwind[5]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[13]*fUpwind[32]+0.223606797749979*alphaDrSurf[12]*fUpwind[32]+0.159719141249985*alphaDrSurf[11]*fUpwind[32]+0.25*alphaDrSurf[0]*fUpwind[32]+0.223606797749979*fUpwind[13]*alphaDrSurf[32]+0.223606797749979*fUpwind[12]*alphaDrSurf[32]+0.159719141249985*fUpwind[11]*alphaDrSurf[32]+0.25*fUpwind[0]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[21]*fUpwind[24]+0.223606797749979*fUpwind[21]*alphaDrSurf[24]+0.2*alphaDrSurf[15]*fUpwind[23]+0.2*fUpwind[15]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[19]*fUpwind[22]+0.223606797749979*fUpwind[19]*alphaDrSurf[22]+0.159719141249985*alphaDrSurf[19]*fUpwind[21]+0.2500000000000001*alphaDrSurf[2]*fUpwind[21]+0.159719141249985*fUpwind[19]*alphaDrSurf[21]+0.2500000000000001*fUpwind[2]*alphaDrSurf[21]+0.2*alphaDrSurf[15]*fUpwind[20]+0.2*fUpwind[15]*alphaDrSurf[20]+0.2500000000000001*alphaDrSurf[3]*fUpwind[19]+0.2500000000000001*fUpwind[3]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[15]+0.223606797749979*fUpwind[1]*alphaDrSurf[15]+0.25*alphaDrSurf[7]*fUpwind[11]+0.25*fUpwind[7]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[6]+0.223606797749979*fUpwind[5]*alphaDrSurf[6]; 
  Ghat[33] = 0.1788854381999831*alphaDrSurf[32]*fUpwind[34]+0.2*alphaDrSurf[7]*fUpwind[34]+0.1788854381999831*fUpwind[32]*alphaDrSurf[34]+0.2*fUpwind[7]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[13]*fUpwind[33]+0.159719141249985*alphaDrSurf[12]*fUpwind[33]+0.223606797749979*alphaDrSurf[11]*fUpwind[33]+0.25*alphaDrSurf[0]*fUpwind[33]+0.223606797749979*fUpwind[13]*alphaDrSurf[33]+0.159719141249985*fUpwind[12]*alphaDrSurf[33]+0.223606797749979*fUpwind[11]*alphaDrSurf[33]+0.25*fUpwind[0]*alphaDrSurf[33]+0.2*alphaDrSurf[5]*fUpwind[32]+0.2*fUpwind[5]*alphaDrSurf[32]+0.2*alphaDrSurf[15]*fUpwind[24]+0.2*fUpwind[15]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[22]*fUpwind[23]+0.223606797749979*fUpwind[22]*alphaDrSurf[23]+0.159719141249985*alphaDrSurf[20]*fUpwind[22]+0.2500000000000001*alphaDrSurf[1]*fUpwind[22]+0.159719141249985*fUpwind[20]*alphaDrSurf[22]+0.2500000000000001*fUpwind[1]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[20]*fUpwind[21]+0.223606797749979*fUpwind[20]*alphaDrSurf[21]+0.2500000000000001*alphaDrSurf[3]*fUpwind[20]+0.2500000000000001*fUpwind[3]*alphaDrSurf[20]+0.2*alphaDrSurf[15]*fUpwind[19]+0.2*fUpwind[15]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[2]*fUpwind[15]+0.223606797749979*fUpwind[2]*alphaDrSurf[15]+0.25*alphaDrSurf[6]*fUpwind[12]+0.25*fUpwind[6]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[7]+0.223606797749979*fUpwind[5]*alphaDrSurf[7]; 
  Ghat[34] = 0.159719141249985*alphaDrSurf[13]*fUpwind[34]+0.223606797749979*alphaDrSurf[12]*fUpwind[34]+0.223606797749979*alphaDrSurf[11]*fUpwind[34]+0.25*alphaDrSurf[0]*fUpwind[34]+0.159719141249985*fUpwind[13]*alphaDrSurf[34]+0.223606797749979*fUpwind[12]*alphaDrSurf[34]+0.223606797749979*fUpwind[11]*alphaDrSurf[34]+0.25*fUpwind[0]*alphaDrSurf[34]+0.1788854381999831*alphaDrSurf[32]*fUpwind[33]+0.2*alphaDrSurf[7]*fUpwind[33]+0.1788854381999831*fUpwind[32]*alphaDrSurf[33]+0.2*fUpwind[7]*alphaDrSurf[33]+0.2*alphaDrSurf[6]*fUpwind[32]+0.2*fUpwind[6]*alphaDrSurf[32]+0.159719141249985*alphaDrSurf[23]*fUpwind[24]+0.223606797749979*alphaDrSurf[20]*fUpwind[24]+0.2500000000000001*alphaDrSurf[1]*fUpwind[24]+0.159719141249985*fUpwind[23]*alphaDrSurf[24]+0.223606797749979*fUpwind[20]*alphaDrSurf[24]+0.2500000000000001*fUpwind[1]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[19]*fUpwind[23]+0.2500000000000001*alphaDrSurf[2]*fUpwind[23]+0.223606797749979*fUpwind[19]*alphaDrSurf[23]+0.2500000000000001*fUpwind[2]*alphaDrSurf[23]+0.2*alphaDrSurf[15]*fUpwind[22]+0.2*fUpwind[15]*alphaDrSurf[22]+0.2*alphaDrSurf[15]*fUpwind[21]+0.2*fUpwind[15]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[3]*fUpwind[15]+0.223606797749979*fUpwind[3]*alphaDrSurf[15]+0.25*alphaDrSurf[5]*fUpwind[13]+0.25*fUpwind[5]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[6]*fUpwind[7]+0.223606797749979*fUpwind[6]*alphaDrSurf[7]; 
  Ghat[35] = 0.223606797749979*alphaDrSurf[23]*fUpwind[46]+0.2*alphaDrSurf[15]*fUpwind[45]+0.223606797749979*alphaDrSurf[22]*fUpwind[44]+0.159719141249985*alphaDrSurf[21]*fUpwind[44]+0.2500000000000001*alphaDrSurf[3]*fUpwind[44]+0.223606797749979*alphaDrSurf[34]*fUpwind[39]+0.223606797749979*alphaDrSurf[32]*fUpwind[38]+0.159719141249985*alphaDrSurf[32]*fUpwind[37]+0.25*alphaDrSurf[7]*fUpwind[37]+0.2*alphaDrSurf[5]*fUpwind[36]+0.223606797749979*alphaDrSurf[12]*fUpwind[35]+0.159719141249985*alphaDrSurf[11]*fUpwind[35]+0.25*alphaDrSurf[0]*fUpwind[35]+0.2*fUpwind[31]*alphaDrSurf[33]+0.25*fUpwind[10]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[6]*fUpwind[31]+0.223606797749979*alphaDrSurf[19]*fUpwind[26]+0.159719141249985*alphaDrSurf[19]*fUpwind[25]+0.2500000000000001*alphaDrSurf[2]*fUpwind[25]+0.2500000000000001*fUpwind[18]*alphaDrSurf[21]+0.2*fUpwind[16]*alphaDrSurf[20]+0.2500000000000001*fUpwind[4]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[15]*fUpwind[17]+0.223606797749979*alphaDrSurf[1]*fUpwind[16]+0.25*fUpwind[9]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[8]; 
  Ghat[36] = 0.223606797749979*alphaDrSurf[24]*fUpwind[46]+0.159719141249985*alphaDrSurf[22]*fUpwind[45]+0.223606797749979*alphaDrSurf[21]*fUpwind[45]+0.2500000000000001*alphaDrSurf[3]*fUpwind[45]+0.2*alphaDrSurf[15]*fUpwind[44]+0.223606797749979*alphaDrSurf[34]*fUpwind[40]+0.159719141249985*alphaDrSurf[33]*fUpwind[38]+0.25*alphaDrSurf[6]*fUpwind[38]+0.223606797749979*alphaDrSurf[33]*fUpwind[37]+0.159719141249985*alphaDrSurf[12]*fUpwind[36]+0.223606797749979*alphaDrSurf[11]*fUpwind[36]+0.25*alphaDrSurf[0]*fUpwind[36]+0.2*alphaDrSurf[5]*fUpwind[35]+0.25*fUpwind[10]*alphaDrSurf[33]+0.2*fUpwind[31]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[7]*fUpwind[31]+0.159719141249985*alphaDrSurf[20]*fUpwind[26]+0.2500000000000001*alphaDrSurf[1]*fUpwind[26]+0.223606797749979*alphaDrSurf[20]*fUpwind[25]+0.2500000000000001*fUpwind[17]*alphaDrSurf[22]+0.2500000000000001*fUpwind[4]*alphaDrSurf[20]+0.2*fUpwind[16]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[15]*fUpwind[18]+0.223606797749979*alphaDrSurf[2]*fUpwind[16]+0.25*fUpwind[8]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[9]; 
  Ghat[37] = 0.2*alphaDrSurf[15]*fUpwind[46]+0.223606797749979*alphaDrSurf[20]*fUpwind[45]+0.223606797749979*alphaDrSurf[24]*fUpwind[44]+0.159719141249985*alphaDrSurf[19]*fUpwind[44]+0.2500000000000001*alphaDrSurf[2]*fUpwind[44]+0.223606797749979*alphaDrSurf[32]*fUpwind[40]+0.2*alphaDrSurf[6]*fUpwind[39]+0.223606797749979*alphaDrSurf[13]*fUpwind[37]+0.159719141249985*alphaDrSurf[11]*fUpwind[37]+0.25*alphaDrSurf[0]*fUpwind[37]+0.223606797749979*alphaDrSurf[33]*fUpwind[36]+0.159719141249985*alphaDrSurf[32]*fUpwind[35]+0.25*alphaDrSurf[7]*fUpwind[35]+0.2*fUpwind[31]*alphaDrSurf[34]+0.25*fUpwind[9]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[5]*fUpwind[31]+0.223606797749979*alphaDrSurf[21]*fUpwind[27]+0.159719141249985*alphaDrSurf[21]*fUpwind[25]+0.2500000000000001*alphaDrSurf[3]*fUpwind[25]+0.2*fUpwind[17]*alphaDrSurf[23]+0.2500000000000001*fUpwind[4]*alphaDrSurf[21]+0.2500000000000001*fUpwind[18]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[17]+0.223606797749979*alphaDrSurf[15]*fUpwind[16]+0.25*fUpwind[10]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[6]*fUpwind[8]; 
  Ghat[38] = 0.2*alphaDrSurf[15]*fUpwind[46]+0.223606797749979*alphaDrSurf[23]*fUpwind[45]+0.159719141249985*alphaDrSurf[20]*fUpwind[45]+0.2500000000000001*alphaDrSurf[1]*fUpwind[45]+0.223606797749979*alphaDrSurf[19]*fUpwind[44]+0.2*alphaDrSurf[7]*fUpwind[40]+0.223606797749979*alphaDrSurf[33]*fUpwind[39]+0.223606797749979*alphaDrSurf[13]*fUpwind[38]+0.159719141249985*alphaDrSurf[12]*fUpwind[38]+0.25*alphaDrSurf[0]*fUpwind[38]+0.159719141249985*alphaDrSurf[33]*fUpwind[36]+0.25*alphaDrSurf[6]*fUpwind[36]+0.223606797749979*alphaDrSurf[32]*fUpwind[35]+0.2*fUpwind[31]*alphaDrSurf[34]+0.25*fUpwind[8]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[5]*fUpwind[31]+0.223606797749979*alphaDrSurf[22]*fUpwind[27]+0.159719141249985*alphaDrSurf[22]*fUpwind[26]+0.2500000000000001*alphaDrSurf[3]*fUpwind[26]+0.2*fUpwind[18]*alphaDrSurf[24]+0.2500000000000001*fUpwind[4]*alphaDrSurf[22]+0.2500000000000001*fUpwind[17]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[2]*fUpwind[18]+0.223606797749979*alphaDrSurf[15]*fUpwind[16]+0.25*fUpwind[10]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[7]*fUpwind[9]; 
  Ghat[39] = 0.159719141249985*alphaDrSurf[24]*fUpwind[46]+0.223606797749979*alphaDrSurf[19]*fUpwind[46]+0.2500000000000001*alphaDrSurf[2]*fUpwind[46]+0.223606797749979*alphaDrSurf[22]*fUpwind[45]+0.2*alphaDrSurf[15]*fUpwind[44]+0.159719141249985*alphaDrSurf[34]*fUpwind[40]+0.25*alphaDrSurf[5]*fUpwind[40]+0.159719141249985*alphaDrSurf[13]*fUpwind[39]+0.223606797749979*alphaDrSurf[11]*fUpwind[39]+0.25*alphaDrSurf[0]*fUpwind[39]+0.223606797749979*alphaDrSurf[33]*fUpwind[38]+0.2*alphaDrSurf[6]*fUpwind[37]+0.223606797749979*alphaDrSurf[34]*fUpwind[35]+0.25*fUpwind[9]*alphaDrSurf[34]+0.2*fUpwind[31]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[7]*fUpwind[31]+0.159719141249985*alphaDrSurf[23]*fUpwind[27]+0.2500000000000001*alphaDrSurf[1]*fUpwind[27]+0.223606797749979*alphaDrSurf[23]*fUpwind[25]+0.2500000000000001*fUpwind[16]*alphaDrSurf[24]+0.2500000000000001*fUpwind[4]*alphaDrSurf[23]+0.2*fUpwind[17]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[15]*fUpwind[18]+0.223606797749979*alphaDrSurf[3]*fUpwind[17]+0.25*fUpwind[8]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[6]*fUpwind[10]; 
  Ghat[40] = 0.159719141249985*alphaDrSurf[23]*fUpwind[46]+0.223606797749979*alphaDrSurf[20]*fUpwind[46]+0.2500000000000001*alphaDrSurf[1]*fUpwind[46]+0.2*alphaDrSurf[15]*fUpwind[45]+0.223606797749979*alphaDrSurf[21]*fUpwind[44]+0.159719141249985*alphaDrSurf[13]*fUpwind[40]+0.223606797749979*alphaDrSurf[12]*fUpwind[40]+0.25*alphaDrSurf[0]*fUpwind[40]+0.159719141249985*alphaDrSurf[34]*fUpwind[39]+0.25*alphaDrSurf[5]*fUpwind[39]+0.2*alphaDrSurf[7]*fUpwind[38]+0.223606797749979*alphaDrSurf[32]*fUpwind[37]+0.223606797749979*alphaDrSurf[34]*fUpwind[36]+0.25*fUpwind[8]*alphaDrSurf[34]+0.2*fUpwind[31]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[6]*fUpwind[31]+0.159719141249985*alphaDrSurf[24]*fUpwind[27]+0.2500000000000001*alphaDrSurf[2]*fUpwind[27]+0.223606797749979*alphaDrSurf[24]*fUpwind[26]+0.2500000000000001*fUpwind[4]*alphaDrSurf[24]+0.2500000000000001*fUpwind[16]*alphaDrSurf[23]+0.2*fUpwind[18]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[3]*fUpwind[18]+0.223606797749979*alphaDrSurf[15]*fUpwind[17]+0.25*fUpwind[9]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[7]*fUpwind[10]; 
  Ghat[41] = 0.223606797749979*alphaDrSurf[22]*fUpwind[47]+0.223606797749979*alphaDrSurf[21]*fUpwind[47]+0.2500000000000001*alphaDrSurf[3]*fUpwind[47]+0.223606797749979*alphaDrSurf[33]*fUpwind[43]+0.25*alphaDrSurf[6]*fUpwind[43]+0.223606797749979*alphaDrSurf[32]*fUpwind[42]+0.25*alphaDrSurf[7]*fUpwind[42]+0.223606797749979*alphaDrSurf[12]*fUpwind[41]+0.223606797749979*alphaDrSurf[11]*fUpwind[41]+0.25*alphaDrSurf[0]*fUpwind[41]+0.2500000000000001*alphaDrSurf[15]*fUpwind[30]+0.223606797749979*alphaDrSurf[20]*fUpwind[29]+0.2500000000000001*alphaDrSurf[1]*fUpwind[29]+0.223606797749979*alphaDrSurf[19]*fUpwind[28]+0.2500000000000001*alphaDrSurf[2]*fUpwind[28]+0.25*alphaDrSurf[5]*fUpwind[14]; 
  Ghat[42] = 0.223606797749979*alphaDrSurf[24]*fUpwind[47]+0.223606797749979*alphaDrSurf[19]*fUpwind[47]+0.2500000000000001*alphaDrSurf[2]*fUpwind[47]+0.223606797749979*alphaDrSurf[34]*fUpwind[43]+0.25*alphaDrSurf[5]*fUpwind[43]+0.223606797749979*alphaDrSurf[13]*fUpwind[42]+0.223606797749979*alphaDrSurf[11]*fUpwind[42]+0.25*alphaDrSurf[0]*fUpwind[42]+0.223606797749979*alphaDrSurf[32]*fUpwind[41]+0.25*alphaDrSurf[7]*fUpwind[41]+0.223606797749979*alphaDrSurf[23]*fUpwind[30]+0.2500000000000001*alphaDrSurf[1]*fUpwind[30]+0.2500000000000001*alphaDrSurf[15]*fUpwind[29]+0.223606797749979*alphaDrSurf[21]*fUpwind[28]+0.2500000000000001*alphaDrSurf[3]*fUpwind[28]+0.25*alphaDrSurf[6]*fUpwind[14]; 
  Ghat[43] = 0.223606797749979*alphaDrSurf[23]*fUpwind[47]+0.223606797749979*alphaDrSurf[20]*fUpwind[47]+0.2500000000000001*alphaDrSurf[1]*fUpwind[47]+0.223606797749979*alphaDrSurf[13]*fUpwind[43]+0.223606797749979*alphaDrSurf[12]*fUpwind[43]+0.25*alphaDrSurf[0]*fUpwind[43]+0.223606797749979*alphaDrSurf[34]*fUpwind[42]+0.25*alphaDrSurf[5]*fUpwind[42]+0.223606797749979*alphaDrSurf[33]*fUpwind[41]+0.25*alphaDrSurf[6]*fUpwind[41]+0.223606797749979*alphaDrSurf[24]*fUpwind[30]+0.2500000000000001*alphaDrSurf[2]*fUpwind[30]+0.223606797749979*alphaDrSurf[22]*fUpwind[29]+0.2500000000000001*alphaDrSurf[3]*fUpwind[29]+0.2500000000000001*alphaDrSurf[15]*fUpwind[28]+0.25*alphaDrSurf[7]*fUpwind[14]; 
  Ghat[44] = 0.1788854381999831*alphaDrSurf[33]*fUpwind[46]+0.2*alphaDrSurf[6]*fUpwind[46]+0.1788854381999831*alphaDrSurf[34]*fUpwind[45]+0.2*alphaDrSurf[5]*fUpwind[45]+0.223606797749979*alphaDrSurf[13]*fUpwind[44]+0.223606797749979*alphaDrSurf[12]*fUpwind[44]+0.159719141249985*alphaDrSurf[11]*fUpwind[44]+0.25*alphaDrSurf[0]*fUpwind[44]+0.223606797749979*alphaDrSurf[21]*fUpwind[40]+0.2*alphaDrSurf[15]*fUpwind[39]+0.223606797749979*alphaDrSurf[19]*fUpwind[38]+0.223606797749979*alphaDrSurf[24]*fUpwind[37]+0.159719141249985*alphaDrSurf[19]*fUpwind[37]+0.2500000000000001*alphaDrSurf[2]*fUpwind[37]+0.2*alphaDrSurf[15]*fUpwind[36]+0.223606797749979*alphaDrSurf[22]*fUpwind[35]+0.159719141249985*alphaDrSurf[21]*fUpwind[35]+0.2500000000000001*alphaDrSurf[3]*fUpwind[35]+0.2*fUpwind[17]*alphaDrSurf[34]+0.2*fUpwind[16]*alphaDrSurf[33]+0.223606797749979*fUpwind[27]*alphaDrSurf[32]+0.223606797749979*fUpwind[26]*alphaDrSurf[32]+0.159719141249985*fUpwind[25]*alphaDrSurf[32]+0.2500000000000001*fUpwind[4]*alphaDrSurf[32]+0.2*alphaDrSurf[23]*fUpwind[31]+0.2*alphaDrSurf[20]*fUpwind[31]+0.223606797749979*alphaDrSurf[1]*fUpwind[31]+0.25*alphaDrSurf[7]*fUpwind[25]+0.25*fUpwind[9]*alphaDrSurf[21]+0.25*fUpwind[10]*alphaDrSurf[19]+0.2500000000000001*alphaDrSurf[11]*fUpwind[18]+0.223606797749979*alphaDrSurf[5]*fUpwind[17]+0.223606797749979*alphaDrSurf[6]*fUpwind[16]+0.223606797749979*fUpwind[8]*alphaDrSurf[15]; 
  Ghat[45] = 0.1788854381999831*alphaDrSurf[32]*fUpwind[46]+0.2*alphaDrSurf[7]*fUpwind[46]+0.223606797749979*alphaDrSurf[13]*fUpwind[45]+0.159719141249985*alphaDrSurf[12]*fUpwind[45]+0.223606797749979*alphaDrSurf[11]*fUpwind[45]+0.25*alphaDrSurf[0]*fUpwind[45]+0.1788854381999831*alphaDrSurf[34]*fUpwind[44]+0.2*alphaDrSurf[5]*fUpwind[44]+0.2*alphaDrSurf[15]*fUpwind[40]+0.223606797749979*alphaDrSurf[22]*fUpwind[39]+0.223606797749979*alphaDrSurf[23]*fUpwind[38]+0.159719141249985*alphaDrSurf[20]*fUpwind[38]+0.2500000000000001*alphaDrSurf[1]*fUpwind[38]+0.223606797749979*alphaDrSurf[20]*fUpwind[37]+0.159719141249985*alphaDrSurf[22]*fUpwind[36]+0.223606797749979*alphaDrSurf[21]*fUpwind[36]+0.2500000000000001*alphaDrSurf[3]*fUpwind[36]+0.2*alphaDrSurf[15]*fUpwind[35]+0.2*fUpwind[18]*alphaDrSurf[34]+0.223606797749979*fUpwind[27]*alphaDrSurf[33]+0.159719141249985*fUpwind[26]*alphaDrSurf[33]+0.223606797749979*fUpwind[25]*alphaDrSurf[33]+0.2500000000000001*fUpwind[4]*alphaDrSurf[33]+0.2*fUpwind[16]*alphaDrSurf[32]+0.2*alphaDrSurf[24]*fUpwind[31]+0.2*alphaDrSurf[19]*fUpwind[31]+0.223606797749979*alphaDrSurf[2]*fUpwind[31]+0.25*alphaDrSurf[6]*fUpwind[26]+0.25*fUpwind[8]*alphaDrSurf[22]+0.25*fUpwind[10]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[18]+0.2500000000000001*alphaDrSurf[12]*fUpwind[17]+0.223606797749979*alphaDrSurf[7]*fUpwind[16]+0.223606797749979*fUpwind[9]*alphaDrSurf[15]; 
  Ghat[46] = 0.159719141249985*alphaDrSurf[13]*fUpwind[46]+0.223606797749979*alphaDrSurf[12]*fUpwind[46]+0.223606797749979*alphaDrSurf[11]*fUpwind[46]+0.25*alphaDrSurf[0]*fUpwind[46]+0.1788854381999831*alphaDrSurf[32]*fUpwind[45]+0.2*alphaDrSurf[7]*fUpwind[45]+0.1788854381999831*alphaDrSurf[33]*fUpwind[44]+0.2*alphaDrSurf[6]*fUpwind[44]+0.159719141249985*alphaDrSurf[23]*fUpwind[40]+0.223606797749979*alphaDrSurf[20]*fUpwind[40]+0.2500000000000001*alphaDrSurf[1]*fUpwind[40]+0.159719141249985*alphaDrSurf[24]*fUpwind[39]+0.223606797749979*alphaDrSurf[19]*fUpwind[39]+0.2500000000000001*alphaDrSurf[2]*fUpwind[39]+0.2*alphaDrSurf[15]*fUpwind[38]+0.2*alphaDrSurf[15]*fUpwind[37]+0.223606797749979*alphaDrSurf[24]*fUpwind[36]+0.223606797749979*alphaDrSurf[23]*fUpwind[35]+0.159719141249985*fUpwind[27]*alphaDrSurf[34]+0.223606797749979*fUpwind[26]*alphaDrSurf[34]+0.223606797749979*fUpwind[25]*alphaDrSurf[34]+0.2500000000000001*fUpwind[4]*alphaDrSurf[34]+0.2*fUpwind[18]*alphaDrSurf[33]+0.2*fUpwind[17]*alphaDrSurf[32]+0.2*alphaDrSurf[22]*fUpwind[31]+0.2*alphaDrSurf[21]*fUpwind[31]+0.223606797749979*alphaDrSurf[3]*fUpwind[31]+0.25*alphaDrSurf[5]*fUpwind[27]+0.25*fUpwind[8]*alphaDrSurf[24]+0.25*fUpwind[9]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[6]*fUpwind[18]+0.223606797749979*alphaDrSurf[7]*fUpwind[17]+0.2500000000000001*alphaDrSurf[13]*fUpwind[16]+0.223606797749979*fUpwind[10]*alphaDrSurf[15]; 
  Ghat[47] = 0.223606797749979*alphaDrSurf[13]*fUpwind[47]+0.223606797749979*alphaDrSurf[12]*fUpwind[47]+0.223606797749979*alphaDrSurf[11]*fUpwind[47]+0.25*alphaDrSurf[0]*fUpwind[47]+0.223606797749979*alphaDrSurf[23]*fUpwind[43]+0.223606797749979*alphaDrSurf[20]*fUpwind[43]+0.2500000000000001*alphaDrSurf[1]*fUpwind[43]+0.223606797749979*alphaDrSurf[24]*fUpwind[42]+0.223606797749979*alphaDrSurf[19]*fUpwind[42]+0.2500000000000001*alphaDrSurf[2]*fUpwind[42]+0.223606797749979*alphaDrSurf[22]*fUpwind[41]+0.223606797749979*alphaDrSurf[21]*fUpwind[41]+0.2500000000000001*alphaDrSurf[3]*fUpwind[41]+0.223606797749979*fUpwind[30]*alphaDrSurf[34]+0.223606797749979*fUpwind[29]*alphaDrSurf[33]+0.223606797749979*fUpwind[28]*alphaDrSurf[32]+0.25*alphaDrSurf[5]*fUpwind[30]+0.25*alphaDrSurf[6]*fUpwind[29]+0.25*alphaDrSurf[7]*fUpwind[28]+0.2500000000000001*fUpwind[14]*alphaDrSurf[15]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[9] += 1.224744871391589*Ghat[1]*rdv2; 
  out[10] += 1.224744871391589*Ghat[2]*rdv2; 
  out[11] += 1.224744871391589*Ghat[3]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[15] += 1.224744871391589*Ghat[4]*rdv2; 
  out[16] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[17] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[18] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[19] += 1.58113883008419*Ghat[0]*rdv2; 
  out[20] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[21] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[22] += 1.224744871391589*Ghat[5]*rdv2; 
  out[23] += 1.224744871391589*Ghat[6]*rdv2; 
  out[24] += 1.224744871391589*Ghat[7]*rdv2; 
  out[25] += 0.7071067811865475*Ghat[16]*rdv2; 
  out[26] += 0.7071067811865475*Ghat[17]*rdv2; 
  out[27] += 0.7071067811865475*Ghat[18]*rdv2; 
  out[28] += 1.224744871391589*Ghat[8]*rdv2; 
  out[29] += 1.224744871391589*Ghat[9]*rdv2; 
  out[30] += 1.224744871391589*Ghat[10]*rdv2; 
  out[31] += 0.7071067811865475*Ghat[19]*rdv2; 
  out[32] += 0.7071067811865475*Ghat[20]*rdv2; 
  out[33] += 0.7071067811865475*Ghat[21]*rdv2; 
  out[34] += 0.7071067811865475*Ghat[22]*rdv2; 
  out[35] += 0.7071067811865475*Ghat[23]*rdv2; 
  out[36] += 0.7071067811865475*Ghat[24]*rdv2; 
  out[37] += 1.224744871391589*Ghat[11]*rdv2; 
  out[38] += 1.224744871391589*Ghat[12]*rdv2; 
  out[39] += 1.224744871391589*Ghat[13]*rdv2; 
  out[40] += 1.58113883008419*Ghat[1]*rdv2; 
  out[41] += 1.58113883008419*Ghat[2]*rdv2; 
  out[42] += 1.58113883008419*Ghat[3]*rdv2; 
  out[43] += 0.7071067811865475*Ghat[25]*rdv2; 
  out[44] += 0.7071067811865475*Ghat[26]*rdv2; 
  out[45] += 0.7071067811865475*Ghat[27]*rdv2; 
  out[46] += 1.58113883008419*Ghat[4]*rdv2; 
  out[47] += 0.7071067811865475*Ghat[28]*rdv2; 
  out[48] += 0.7071067811865475*Ghat[29]*rdv2; 
  out[49] += 0.7071067811865475*Ghat[30]*rdv2; 
  out[50] += 1.224744871391589*Ghat[14]*rdv2; 
  out[51] += 1.224744871391589*Ghat[15]*rdv2; 
  out[52] += 0.7071067811865475*Ghat[31]*rdv2; 
  out[53] += 1.224744871391589*Ghat[16]*rdv2; 
  out[54] += 1.224744871391589*Ghat[17]*rdv2; 
  out[55] += 1.224744871391589*Ghat[18]*rdv2; 
  out[56] += 0.7071067811865475*Ghat[32]*rdv2; 
  out[57] += 0.7071067811865475*Ghat[33]*rdv2; 
  out[58] += 0.7071067811865475*Ghat[34]*rdv2; 
  out[59] += 1.224744871391589*Ghat[19]*rdv2; 
  out[60] += 1.224744871391589*Ghat[20]*rdv2; 
  out[61] += 1.224744871391589*Ghat[21]*rdv2; 
  out[62] += 1.224744871391589*Ghat[22]*rdv2; 
  out[63] += 1.224744871391589*Ghat[23]*rdv2; 
  out[64] += 1.224744871391589*Ghat[24]*rdv2; 
  out[65] += 1.58113883008419*Ghat[5]*rdv2; 
  out[66] += 1.58113883008419*Ghat[6]*rdv2; 
  out[67] += 1.58113883008419*Ghat[7]*rdv2; 
  out[68] += 0.7071067811865475*Ghat[35]*rdv2; 
  out[69] += 0.7071067811865475*Ghat[36]*rdv2; 
  out[70] += 0.7071067811865475*Ghat[37]*rdv2; 
  out[71] += 0.7071067811865475*Ghat[38]*rdv2; 
  out[72] += 0.7071067811865475*Ghat[39]*rdv2; 
  out[73] += 0.7071067811865475*Ghat[40]*rdv2; 
  out[74] += 1.224744871391589*Ghat[25]*rdv2; 
  out[75] += 1.224744871391589*Ghat[26]*rdv2; 
  out[76] += 1.224744871391589*Ghat[27]*rdv2; 
  out[77] += 1.58113883008419*Ghat[8]*rdv2; 
  out[78] += 1.58113883008419*Ghat[9]*rdv2; 
  out[79] += 1.58113883008419*Ghat[10]*rdv2; 
  out[80] += 0.7071067811865475*Ghat[41]*rdv2; 
  out[81] += 0.7071067811865475*Ghat[42]*rdv2; 
  out[82] += 0.7071067811865475*Ghat[43]*rdv2; 
  out[83] += 1.224744871391589*Ghat[28]*rdv2; 
  out[84] += 1.224744871391589*Ghat[29]*rdv2; 
  out[85] += 1.224744871391589*Ghat[30]*rdv2; 
  out[86] += 1.224744871391589*Ghat[31]*rdv2; 
  out[87] += 1.224744871391589*Ghat[32]*rdv2; 
  out[88] += 1.224744871391589*Ghat[33]*rdv2; 
  out[89] += 1.224744871391589*Ghat[34]*rdv2; 
  out[90] += 1.58113883008419*Ghat[15]*rdv2; 
  out[91] += 0.7071067811865475*Ghat[44]*rdv2; 
  out[92] += 0.7071067811865475*Ghat[45]*rdv2; 
  out[93] += 0.7071067811865475*Ghat[46]*rdv2; 
  out[94] += 1.224744871391589*Ghat[35]*rdv2; 
  out[95] += 1.224744871391589*Ghat[36]*rdv2; 
  out[96] += 1.224744871391589*Ghat[37]*rdv2; 
  out[97] += 1.224744871391589*Ghat[38]*rdv2; 
  out[98] += 1.224744871391589*Ghat[39]*rdv2; 
  out[99] += 1.224744871391589*Ghat[40]*rdv2; 
  out[100] += 1.58113883008419*Ghat[16]*rdv2; 
  out[101] += 1.58113883008419*Ghat[17]*rdv2; 
  out[102] += 1.58113883008419*Ghat[18]*rdv2; 
  out[103] += 0.7071067811865475*Ghat[47]*rdv2; 
  out[104] += 1.224744871391589*Ghat[41]*rdv2; 
  out[105] += 1.224744871391589*Ghat[42]*rdv2; 
  out[106] += 1.224744871391589*Ghat[43]*rdv2; 
  out[107] += 1.224744871391589*Ghat[44]*rdv2; 
  out[108] += 1.224744871391589*Ghat[45]*rdv2; 
  out[109] += 1.224744871391589*Ghat[46]*rdv2; 
  out[110] += 1.58113883008419*Ghat[31]*rdv2; 
  out[111] += 1.224744871391589*Ghat[47]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.828427124746191*w[3]-1.414213562373095*dxv[3])-2.828427124746191*nuUSum[0]); 
  alphaDrSurf[1] = 0.5*(nuSum[1]*(2.828427124746191*w[3]-1.414213562373095*dxv[3])-2.828427124746191*nuUSum[1]); 
  alphaDrSurf[2] = 0.5*(nuSum[2]*(2.828427124746191*w[3]-1.414213562373095*dxv[3])-2.828427124746191*nuUSum[2]); 
  alphaDrSurf[3] = 0.5*(2.828427124746191*nuSum[3]*w[3]-2.828427124746191*nuUSum[3]-1.414213562373095*dxv[3]*nuSum[3]); 
  alphaDrSurf[5] = -0.5*(2.828427124746191*nuUSum[4]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[4]); 
  alphaDrSurf[6] = -0.5*(2.828427124746191*nuUSum[5]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[5]); 
  alphaDrSurf[7] = -0.5*(2.828427124746191*nuUSum[6]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[6]); 
  alphaDrSurf[11] = -0.5*(2.828427124746191*nuUSum[7]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[7]); 
  alphaDrSurf[12] = -0.5*(2.828427124746191*nuUSum[8]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[8]); 
  alphaDrSurf[13] = -0.5*(2.828427124746191*nuUSum[9]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[9]); 
  alphaDrSurf[15] = -0.5*(2.828427124746191*nuUSum[10]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[10]); 
  alphaDrSurf[19] = -0.5*(2.828427124746191*nuUSum[11]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[11]); 
  alphaDrSurf[20] = -0.5*(2.828427124746191*nuUSum[12]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[12]); 
  alphaDrSurf[21] = -0.5*(2.828427124746191*nuUSum[13]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[13]); 
  alphaDrSurf[22] = -0.5*(2.828427124746191*nuUSum[14]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[14]); 
  alphaDrSurf[23] = -0.5*(2.828427124746191*nuUSum[15]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[15]); 
  alphaDrSurf[24] = -0.5*(2.828427124746191*nuUSum[16]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[16]); 
  alphaDrSurf[32] = -0.5*(2.828427124746191*nuUSum[17]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[17]); 
  alphaDrSurf[33] = -0.5*(2.828427124746191*nuUSum[18]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[18]); 
  alphaDrSurf[34] = -0.5*(2.828427124746191*nuUSum[19]+(1.414213562373095*dxv[3]-2.828427124746191*w[3])*nuSum[19]); 

  if (0.4024922359499623*alphaDrSurf[34]+0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx4_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_5x_p2_surfx4_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = ser_5x_p2_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx4_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_5x_p2_surfx4_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = ser_5x_p2_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]+0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_5x_p2_surfx4_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = ser_5x_p2_surfx4_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = ser_5x_p2_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_5x_p2_surfx4_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = ser_5x_p2_surfx4_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = ser_5x_p2_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]+0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_5x_p2_surfx4_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = ser_5x_p2_surfx4_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = ser_5x_p2_surfx4_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_5x_p2_surfx4_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = ser_5x_p2_surfx4_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = ser_5x_p2_surfx4_eval_quad_node_8_l(fSkin); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[34])+0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_5x_p2_surfx4_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = ser_5x_p2_surfx4_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = ser_5x_p2_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_5x_p2_surfx4_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = ser_5x_p2_surfx4_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = ser_5x_p2_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[34])+0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_5x_p2_surfx4_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = ser_5x_p2_surfx4_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = ser_5x_p2_surfx4_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_5x_p2_surfx4_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = ser_5x_p2_surfx4_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = ser_5x_p2_surfx4_eval_quad_node_14_l(fSkin); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[34])+0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_5x_p2_surfx4_eval_quad_node_15_r(fEdge); 
    fUpwindQuad[16] = ser_5x_p2_surfx4_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = ser_5x_p2_surfx4_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_5x_p2_surfx4_eval_quad_node_15_l(fSkin); 
    fUpwindQuad[16] = ser_5x_p2_surfx4_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = ser_5x_p2_surfx4_eval_quad_node_17_l(fSkin); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]-0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_5x_p2_surfx4_eval_quad_node_18_r(fEdge); 
    fUpwindQuad[19] = ser_5x_p2_surfx4_eval_quad_node_19_r(fEdge); 
    fUpwindQuad[20] = ser_5x_p2_surfx4_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_5x_p2_surfx4_eval_quad_node_18_l(fSkin); 
    fUpwindQuad[19] = ser_5x_p2_surfx4_eval_quad_node_19_l(fSkin); 
    fUpwindQuad[20] = ser_5x_p2_surfx4_eval_quad_node_20_l(fSkin); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]-0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_5x_p2_surfx4_eval_quad_node_21_r(fEdge); 
    fUpwindQuad[22] = ser_5x_p2_surfx4_eval_quad_node_22_r(fEdge); 
    fUpwindQuad[23] = ser_5x_p2_surfx4_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[21] = ser_5x_p2_surfx4_eval_quad_node_21_l(fSkin); 
    fUpwindQuad[22] = ser_5x_p2_surfx4_eval_quad_node_22_l(fSkin); 
    fUpwindQuad[23] = ser_5x_p2_surfx4_eval_quad_node_23_l(fSkin); 
  } 
  if (0.4024922359499623*alphaDrSurf[34]-0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]-0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]-0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]+0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_5x_p2_surfx4_eval_quad_node_24_r(fEdge); 
    fUpwindQuad[25] = ser_5x_p2_surfx4_eval_quad_node_25_r(fEdge); 
    fUpwindQuad[26] = ser_5x_p2_surfx4_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_5x_p2_surfx4_eval_quad_node_24_l(fSkin); 
    fUpwindQuad[25] = ser_5x_p2_surfx4_eval_quad_node_25_l(fSkin); 
    fUpwindQuad[26] = ser_5x_p2_surfx4_eval_quad_node_26_l(fSkin); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[33])-0.3*alphaDrSurf[23]+0.375*alphaDrSurf[22]-0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[6]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[27] = ser_5x_p2_surfx4_eval_quad_node_27_r(fEdge); 
    fUpwindQuad[28] = ser_5x_p2_surfx4_eval_quad_node_28_r(fEdge); 
    fUpwindQuad[29] = ser_5x_p2_surfx4_eval_quad_node_29_r(fEdge); 
  } else { 
    fUpwindQuad[27] = ser_5x_p2_surfx4_eval_quad_node_27_l(fSkin); 
    fUpwindQuad[28] = ser_5x_p2_surfx4_eval_quad_node_28_l(fSkin); 
    fUpwindQuad[29] = ser_5x_p2_surfx4_eval_quad_node_29_l(fSkin); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[33])-0.3*alphaDrSurf[23]+0.375*alphaDrSurf[22]-0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[6]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[30] = ser_5x_p2_surfx4_eval_quad_node_30_r(fEdge); 
    fUpwindQuad[31] = ser_5x_p2_surfx4_eval_quad_node_31_r(fEdge); 
    fUpwindQuad[32] = ser_5x_p2_surfx4_eval_quad_node_32_r(fEdge); 
  } else { 
    fUpwindQuad[30] = ser_5x_p2_surfx4_eval_quad_node_30_l(fSkin); 
    fUpwindQuad[31] = ser_5x_p2_surfx4_eval_quad_node_31_l(fSkin); 
    fUpwindQuad[32] = ser_5x_p2_surfx4_eval_quad_node_32_l(fSkin); 
  } 
  if ((-0.5031152949374518*alphaDrSurf[33])-0.3*alphaDrSurf[23]+0.375*alphaDrSurf[22]-0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[6]-0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[33] = ser_5x_p2_surfx4_eval_quad_node_33_r(fEdge); 
    fUpwindQuad[34] = ser_5x_p2_surfx4_eval_quad_node_34_r(fEdge); 
    fUpwindQuad[35] = ser_5x_p2_surfx4_eval_quad_node_35_r(fEdge); 
  } else { 
    fUpwindQuad[33] = ser_5x_p2_surfx4_eval_quad_node_33_l(fSkin); 
    fUpwindQuad[34] = ser_5x_p2_surfx4_eval_quad_node_34_l(fSkin); 
    fUpwindQuad[35] = ser_5x_p2_surfx4_eval_quad_node_35_l(fSkin); 
  } 
  if (0.375*alphaDrSurf[23]+0.375*alphaDrSurf[20]-0.2795084971874732*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[36] = ser_5x_p2_surfx4_eval_quad_node_36_r(fEdge); 
    fUpwindQuad[37] = ser_5x_p2_surfx4_eval_quad_node_37_r(fEdge); 
    fUpwindQuad[38] = ser_5x_p2_surfx4_eval_quad_node_38_r(fEdge); 
  } else { 
    fUpwindQuad[36] = ser_5x_p2_surfx4_eval_quad_node_36_l(fSkin); 
    fUpwindQuad[37] = ser_5x_p2_surfx4_eval_quad_node_37_l(fSkin); 
    fUpwindQuad[38] = ser_5x_p2_surfx4_eval_quad_node_38_l(fSkin); 
  } 
  if (0.375*alphaDrSurf[23]+0.375*alphaDrSurf[20]-0.2795084971874732*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[39] = ser_5x_p2_surfx4_eval_quad_node_39_r(fEdge); 
    fUpwindQuad[40] = ser_5x_p2_surfx4_eval_quad_node_40_r(fEdge); 
    fUpwindQuad[41] = ser_5x_p2_surfx4_eval_quad_node_41_r(fEdge); 
  } else { 
    fUpwindQuad[39] = ser_5x_p2_surfx4_eval_quad_node_39_l(fSkin); 
    fUpwindQuad[40] = ser_5x_p2_surfx4_eval_quad_node_40_l(fSkin); 
    fUpwindQuad[41] = ser_5x_p2_surfx4_eval_quad_node_41_l(fSkin); 
  } 
  if (0.375*alphaDrSurf[23]+0.375*alphaDrSurf[20]-0.2795084971874732*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[42] = ser_5x_p2_surfx4_eval_quad_node_42_r(fEdge); 
    fUpwindQuad[43] = ser_5x_p2_surfx4_eval_quad_node_43_r(fEdge); 
    fUpwindQuad[44] = ser_5x_p2_surfx4_eval_quad_node_44_r(fEdge); 
  } else { 
    fUpwindQuad[42] = ser_5x_p2_surfx4_eval_quad_node_42_l(fSkin); 
    fUpwindQuad[43] = ser_5x_p2_surfx4_eval_quad_node_43_l(fSkin); 
    fUpwindQuad[44] = ser_5x_p2_surfx4_eval_quad_node_44_l(fSkin); 
  } 
  if (0.5031152949374518*alphaDrSurf[33]-0.3*alphaDrSurf[23]-0.375*alphaDrSurf[22]+0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[6]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[45] = ser_5x_p2_surfx4_eval_quad_node_45_r(fEdge); 
    fUpwindQuad[46] = ser_5x_p2_surfx4_eval_quad_node_46_r(fEdge); 
    fUpwindQuad[47] = ser_5x_p2_surfx4_eval_quad_node_47_r(fEdge); 
  } else { 
    fUpwindQuad[45] = ser_5x_p2_surfx4_eval_quad_node_45_l(fSkin); 
    fUpwindQuad[46] = ser_5x_p2_surfx4_eval_quad_node_46_l(fSkin); 
    fUpwindQuad[47] = ser_5x_p2_surfx4_eval_quad_node_47_l(fSkin); 
  } 
  if (0.5031152949374518*alphaDrSurf[33]-0.3*alphaDrSurf[23]-0.375*alphaDrSurf[22]+0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[6]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[48] = ser_5x_p2_surfx4_eval_quad_node_48_r(fEdge); 
    fUpwindQuad[49] = ser_5x_p2_surfx4_eval_quad_node_49_r(fEdge); 
    fUpwindQuad[50] = ser_5x_p2_surfx4_eval_quad_node_50_r(fEdge); 
  } else { 
    fUpwindQuad[48] = ser_5x_p2_surfx4_eval_quad_node_48_l(fSkin); 
    fUpwindQuad[49] = ser_5x_p2_surfx4_eval_quad_node_49_l(fSkin); 
    fUpwindQuad[50] = ser_5x_p2_surfx4_eval_quad_node_50_l(fSkin); 
  } 
  if (0.5031152949374518*alphaDrSurf[33]-0.3*alphaDrSurf[23]-0.375*alphaDrSurf[22]+0.3*alphaDrSurf[21]+0.375*alphaDrSurf[20]+0.2236067977499786*alphaDrSurf[13]-0.2795084971874732*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[6]+0.3354101966249678*alphaDrSurf[3]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[51] = ser_5x_p2_surfx4_eval_quad_node_51_r(fEdge); 
    fUpwindQuad[52] = ser_5x_p2_surfx4_eval_quad_node_52_r(fEdge); 
    fUpwindQuad[53] = ser_5x_p2_surfx4_eval_quad_node_53_r(fEdge); 
  } else { 
    fUpwindQuad[51] = ser_5x_p2_surfx4_eval_quad_node_51_l(fSkin); 
    fUpwindQuad[52] = ser_5x_p2_surfx4_eval_quad_node_52_l(fSkin); 
    fUpwindQuad[53] = ser_5x_p2_surfx4_eval_quad_node_53_l(fSkin); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])+0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[54] = ser_5x_p2_surfx4_eval_quad_node_54_r(fEdge); 
    fUpwindQuad[55] = ser_5x_p2_surfx4_eval_quad_node_55_r(fEdge); 
    fUpwindQuad[56] = ser_5x_p2_surfx4_eval_quad_node_56_r(fEdge); 
  } else { 
    fUpwindQuad[54] = ser_5x_p2_surfx4_eval_quad_node_54_l(fSkin); 
    fUpwindQuad[55] = ser_5x_p2_surfx4_eval_quad_node_55_l(fSkin); 
    fUpwindQuad[56] = ser_5x_p2_surfx4_eval_quad_node_56_l(fSkin); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])+0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[57] = ser_5x_p2_surfx4_eval_quad_node_57_r(fEdge); 
    fUpwindQuad[58] = ser_5x_p2_surfx4_eval_quad_node_58_r(fEdge); 
    fUpwindQuad[59] = ser_5x_p2_surfx4_eval_quad_node_59_r(fEdge); 
  } else { 
    fUpwindQuad[57] = ser_5x_p2_surfx4_eval_quad_node_57_l(fSkin); 
    fUpwindQuad[58] = ser_5x_p2_surfx4_eval_quad_node_58_l(fSkin); 
    fUpwindQuad[59] = ser_5x_p2_surfx4_eval_quad_node_59_l(fSkin); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])+0.4024922359499623*alphaDrSurf[33]-0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]-0.3*alphaDrSurf[22]-0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]+0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[7]+0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]-0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[60] = ser_5x_p2_surfx4_eval_quad_node_60_r(fEdge); 
    fUpwindQuad[61] = ser_5x_p2_surfx4_eval_quad_node_61_r(fEdge); 
    fUpwindQuad[62] = ser_5x_p2_surfx4_eval_quad_node_62_r(fEdge); 
  } else { 
    fUpwindQuad[60] = ser_5x_p2_surfx4_eval_quad_node_60_l(fSkin); 
    fUpwindQuad[61] = ser_5x_p2_surfx4_eval_quad_node_61_l(fSkin); 
    fUpwindQuad[62] = ser_5x_p2_surfx4_eval_quad_node_62_l(fSkin); 
  } 
  if (0.5031152949374518*alphaDrSurf[34]-0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[63] = ser_5x_p2_surfx4_eval_quad_node_63_r(fEdge); 
    fUpwindQuad[64] = ser_5x_p2_surfx4_eval_quad_node_64_r(fEdge); 
    fUpwindQuad[65] = ser_5x_p2_surfx4_eval_quad_node_65_r(fEdge); 
  } else { 
    fUpwindQuad[63] = ser_5x_p2_surfx4_eval_quad_node_63_l(fSkin); 
    fUpwindQuad[64] = ser_5x_p2_surfx4_eval_quad_node_64_l(fSkin); 
    fUpwindQuad[65] = ser_5x_p2_surfx4_eval_quad_node_65_l(fSkin); 
  } 
  if (0.5031152949374518*alphaDrSurf[34]-0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[66] = ser_5x_p2_surfx4_eval_quad_node_66_r(fEdge); 
    fUpwindQuad[67] = ser_5x_p2_surfx4_eval_quad_node_67_r(fEdge); 
    fUpwindQuad[68] = ser_5x_p2_surfx4_eval_quad_node_68_r(fEdge); 
  } else { 
    fUpwindQuad[66] = ser_5x_p2_surfx4_eval_quad_node_66_l(fSkin); 
    fUpwindQuad[67] = ser_5x_p2_surfx4_eval_quad_node_67_l(fSkin); 
    fUpwindQuad[68] = ser_5x_p2_surfx4_eval_quad_node_68_l(fSkin); 
  } 
  if (0.5031152949374518*alphaDrSurf[34]-0.375*alphaDrSurf[24]+0.375*alphaDrSurf[23]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.2795084971874732*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[69] = ser_5x_p2_surfx4_eval_quad_node_69_r(fEdge); 
    fUpwindQuad[70] = ser_5x_p2_surfx4_eval_quad_node_70_r(fEdge); 
    fUpwindQuad[71] = ser_5x_p2_surfx4_eval_quad_node_71_r(fEdge); 
  } else { 
    fUpwindQuad[69] = ser_5x_p2_surfx4_eval_quad_node_69_l(fSkin); 
    fUpwindQuad[70] = ser_5x_p2_surfx4_eval_quad_node_70_l(fSkin); 
    fUpwindQuad[71] = ser_5x_p2_surfx4_eval_quad_node_71_l(fSkin); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])-0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[72] = ser_5x_p2_surfx4_eval_quad_node_72_r(fEdge); 
    fUpwindQuad[73] = ser_5x_p2_surfx4_eval_quad_node_73_r(fEdge); 
    fUpwindQuad[74] = ser_5x_p2_surfx4_eval_quad_node_74_r(fEdge); 
  } else { 
    fUpwindQuad[72] = ser_5x_p2_surfx4_eval_quad_node_72_l(fSkin); 
    fUpwindQuad[73] = ser_5x_p2_surfx4_eval_quad_node_73_l(fSkin); 
    fUpwindQuad[74] = ser_5x_p2_surfx4_eval_quad_node_74_l(fSkin); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])-0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[75] = ser_5x_p2_surfx4_eval_quad_node_75_r(fEdge); 
    fUpwindQuad[76] = ser_5x_p2_surfx4_eval_quad_node_76_r(fEdge); 
    fUpwindQuad[77] = ser_5x_p2_surfx4_eval_quad_node_77_r(fEdge); 
  } else { 
    fUpwindQuad[75] = ser_5x_p2_surfx4_eval_quad_node_75_l(fSkin); 
    fUpwindQuad[76] = ser_5x_p2_surfx4_eval_quad_node_76_l(fSkin); 
    fUpwindQuad[77] = ser_5x_p2_surfx4_eval_quad_node_77_l(fSkin); 
  } 
  if ((-0.4024922359499623*alphaDrSurf[34])-0.4024922359499623*alphaDrSurf[33]+0.4024922359499623*alphaDrSurf[32]+0.3*alphaDrSurf[24]-0.3*alphaDrSurf[23]+0.3*alphaDrSurf[22]+0.3*alphaDrSurf[21]-0.3*alphaDrSurf[20]+0.3*alphaDrSurf[19]-0.603738353924943*alphaDrSurf[15]+0.2236067977499786*alphaDrSurf[13]+0.2236067977499786*alphaDrSurf[12]+0.2236067977499786*alphaDrSurf[11]+0.45*alphaDrSurf[7]-0.45*alphaDrSurf[6]-0.45*alphaDrSurf[5]+0.3354101966249678*alphaDrSurf[3]+0.3354101966249678*alphaDrSurf[2]-0.3354101966249678*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[78] = ser_5x_p2_surfx4_eval_quad_node_78_r(fEdge); 
    fUpwindQuad[79] = ser_5x_p2_surfx4_eval_quad_node_79_r(fEdge); 
    fUpwindQuad[80] = ser_5x_p2_surfx4_eval_quad_node_80_r(fEdge); 
  } else { 
    fUpwindQuad[78] = ser_5x_p2_surfx4_eval_quad_node_78_l(fSkin); 
    fUpwindQuad[79] = ser_5x_p2_surfx4_eval_quad_node_79_l(fSkin); 
    fUpwindQuad[80] = ser_5x_p2_surfx4_eval_quad_node_80_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_5x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.25*alphaDrSurf[34]*fUpwind[34]+0.25*alphaDrSurf[33]*fUpwind[33]+0.25*alphaDrSurf[32]*fUpwind[32]+0.25*alphaDrSurf[24]*fUpwind[24]+0.25*alphaDrSurf[23]*fUpwind[23]+0.25*alphaDrSurf[22]*fUpwind[22]+0.25*alphaDrSurf[21]*fUpwind[21]+0.25*alphaDrSurf[20]*fUpwind[20]+0.25*alphaDrSurf[19]*fUpwind[19]+0.25*alphaDrSurf[15]*fUpwind[15]+0.25*alphaDrSurf[13]*fUpwind[13]+0.25*alphaDrSurf[12]*fUpwind[12]+0.25*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[7]*fUpwind[7]+0.25*alphaDrSurf[6]*fUpwind[6]+0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[3]*fUpwind[3]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.2500000000000001*alphaDrSurf[24]*fUpwind[34]+0.2500000000000001*fUpwind[24]*alphaDrSurf[34]+0.2500000000000001*alphaDrSurf[22]*fUpwind[33]+0.2500000000000001*fUpwind[22]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[15]*fUpwind[32]+0.223606797749979*fUpwind[15]*alphaDrSurf[32]+0.2500000000000001*alphaDrSurf[13]*fUpwind[23]+0.2500000000000001*fUpwind[13]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[6]*fUpwind[21]+0.223606797749979*fUpwind[6]*alphaDrSurf[21]+0.2500000000000001*alphaDrSurf[12]*fUpwind[20]+0.2500000000000001*fUpwind[12]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[19]+0.223606797749979*fUpwind[5]*alphaDrSurf[19]+0.25*alphaDrSurf[7]*fUpwind[15]+0.25*fUpwind[7]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[1]*fUpwind[11]+0.223606797749979*fUpwind[1]*alphaDrSurf[11]+0.25*alphaDrSurf[3]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[6]+0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.2500000000000001*alphaDrSurf[23]*fUpwind[34]+0.2500000000000001*fUpwind[23]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[15]*fUpwind[33]+0.223606797749979*fUpwind[15]*alphaDrSurf[33]+0.2500000000000001*alphaDrSurf[21]*fUpwind[32]+0.2500000000000001*fUpwind[21]*alphaDrSurf[32]+0.2500000000000001*alphaDrSurf[13]*fUpwind[24]+0.2500000000000001*fUpwind[13]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[7]*fUpwind[22]+0.223606797749979*fUpwind[7]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[5]*fUpwind[20]+0.223606797749979*fUpwind[5]*alphaDrSurf[20]+0.2500000000000001*alphaDrSurf[11]*fUpwind[19]+0.2500000000000001*fUpwind[11]*alphaDrSurf[19]+0.25*alphaDrSurf[6]*fUpwind[15]+0.25*fUpwind[6]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[2]*fUpwind[12]+0.223606797749979*fUpwind[2]*alphaDrSurf[12]+0.25*alphaDrSurf[3]*fUpwind[7]+0.25*fUpwind[3]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.223606797749979*alphaDrSurf[15]*fUpwind[34]+0.223606797749979*fUpwind[15]*alphaDrSurf[34]+0.2500000000000001*alphaDrSurf[20]*fUpwind[33]+0.2500000000000001*fUpwind[20]*alphaDrSurf[33]+0.2500000000000001*alphaDrSurf[19]*fUpwind[32]+0.2500000000000001*fUpwind[19]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[7]*fUpwind[24]+0.223606797749979*fUpwind[7]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[6]*fUpwind[23]+0.223606797749979*fUpwind[6]*alphaDrSurf[23]+0.2500000000000001*alphaDrSurf[12]*fUpwind[22]+0.2500000000000001*fUpwind[12]*alphaDrSurf[22]+0.2500000000000001*alphaDrSurf[11]*fUpwind[21]+0.2500000000000001*fUpwind[11]*alphaDrSurf[21]+0.25*alphaDrSurf[5]*fUpwind[15]+0.25*fUpwind[5]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[3]*fUpwind[13]+0.223606797749979*fUpwind[3]*alphaDrSurf[13]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*fUpwind[2]*alphaDrSurf[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*fUpwind[1]*alphaDrSurf[6]+0.25*alphaDrSurf[0]*fUpwind[3]+0.25*fUpwind[0]*alphaDrSurf[3]; 
  Ghat[4] = 0.2500000000000001*alphaDrSurf[34]*fUpwind[46]+0.2500000000000001*alphaDrSurf[33]*fUpwind[45]+0.2500000000000001*alphaDrSurf[32]*fUpwind[44]+0.2500000000000001*alphaDrSurf[24]*fUpwind[40]+0.2500000000000001*alphaDrSurf[23]*fUpwind[39]+0.2500000000000001*alphaDrSurf[22]*fUpwind[38]+0.2500000000000001*alphaDrSurf[21]*fUpwind[37]+0.2500000000000001*alphaDrSurf[20]*fUpwind[36]+0.2500000000000001*alphaDrSurf[19]*fUpwind[35]+0.25*alphaDrSurf[15]*fUpwind[31]+0.2500000000000001*alphaDrSurf[13]*fUpwind[27]+0.2500000000000001*alphaDrSurf[12]*fUpwind[26]+0.2500000000000001*alphaDrSurf[11]*fUpwind[25]+0.25*alphaDrSurf[7]*fUpwind[18]+0.25*alphaDrSurf[6]*fUpwind[17]+0.25*alphaDrSurf[5]*fUpwind[16]+0.25*alphaDrSurf[3]*fUpwind[10]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  Ghat[5] = 0.25*alphaDrSurf[13]*fUpwind[34]+0.25*fUpwind[13]*alphaDrSurf[34]+0.2*alphaDrSurf[32]*fUpwind[33]+0.223606797749979*alphaDrSurf[7]*fUpwind[33]+0.2*fUpwind[32]*alphaDrSurf[33]+0.223606797749979*fUpwind[7]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[6]*fUpwind[32]+0.223606797749979*fUpwind[6]*alphaDrSurf[32]+0.25*alphaDrSurf[23]*fUpwind[24]+0.25*fUpwind[23]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[15]*fUpwind[22]+0.223606797749979*fUpwind[15]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[15]*fUpwind[21]+0.223606797749979*fUpwind[15]*alphaDrSurf[21]+0.2*alphaDrSurf[19]*fUpwind[20]+0.223606797749979*alphaDrSurf[2]*fUpwind[20]+0.2*fUpwind[19]*alphaDrSurf[20]+0.223606797749979*fUpwind[2]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[1]*fUpwind[19]+0.223606797749979*fUpwind[1]*alphaDrSurf[19]+0.25*alphaDrSurf[3]*fUpwind[15]+0.25*fUpwind[3]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[5]*fUpwind[12]+0.223606797749979*fUpwind[5]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[11]+0.223606797749979*fUpwind[5]*alphaDrSurf[11]+0.25*alphaDrSurf[6]*fUpwind[7]+0.25*fUpwind[6]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[6] = 0.2*alphaDrSurf[32]*fUpwind[34]+0.223606797749979*alphaDrSurf[7]*fUpwind[34]+0.2*fUpwind[32]*alphaDrSurf[34]+0.223606797749979*fUpwind[7]*alphaDrSurf[34]+0.25*alphaDrSurf[12]*fUpwind[33]+0.25*fUpwind[12]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[5]*fUpwind[32]+0.223606797749979*fUpwind[5]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[15]*fUpwind[24]+0.223606797749979*fUpwind[15]*alphaDrSurf[24]+0.2*alphaDrSurf[21]*fUpwind[23]+0.223606797749979*alphaDrSurf[3]*fUpwind[23]+0.2*fUpwind[21]*alphaDrSurf[23]+0.223606797749979*fUpwind[3]*alphaDrSurf[23]+0.25*alphaDrSurf[20]*fUpwind[22]+0.25*fUpwind[20]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[1]*fUpwind[21]+0.223606797749979*fUpwind[1]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[15]*fUpwind[19]+0.223606797749979*fUpwind[15]*alphaDrSurf[19]+0.25*alphaDrSurf[2]*fUpwind[15]+0.25*fUpwind[2]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[6]*fUpwind[13]+0.223606797749979*fUpwind[6]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[6]*fUpwind[11]+0.223606797749979*fUpwind[6]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*fUpwind[5]*alphaDrSurf[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*fUpwind[0]*alphaDrSurf[6]+0.25*alphaDrSurf[1]*fUpwind[3]+0.25*fUpwind[1]*alphaDrSurf[3]; 
  Ghat[7] = 0.2*alphaDrSurf[33]*fUpwind[34]+0.223606797749979*alphaDrSurf[6]*fUpwind[34]+0.2*fUpwind[33]*alphaDrSurf[34]+0.223606797749979*fUpwind[6]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[5]*fUpwind[33]+0.223606797749979*fUpwind[5]*alphaDrSurf[33]+0.25*alphaDrSurf[11]*fUpwind[32]+0.25*fUpwind[11]*alphaDrSurf[32]+0.2*alphaDrSurf[22]*fUpwind[24]+0.223606797749979*alphaDrSurf[3]*fUpwind[24]+0.2*fUpwind[22]*alphaDrSurf[24]+0.223606797749979*fUpwind[3]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[15]*fUpwind[23]+0.223606797749979*fUpwind[15]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[2]*fUpwind[22]+0.223606797749979*fUpwind[2]*alphaDrSurf[22]+0.25*alphaDrSurf[19]*fUpwind[21]+0.25*fUpwind[19]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[15]*fUpwind[20]+0.223606797749979*fUpwind[15]*alphaDrSurf[20]+0.25*alphaDrSurf[1]*fUpwind[15]+0.25*fUpwind[1]*alphaDrSurf[15]+0.223606797749979*alphaDrSurf[7]*fUpwind[13]+0.223606797749979*fUpwind[7]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[7]*fUpwind[12]+0.223606797749979*fUpwind[7]*alphaDrSurf[12]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*fUpwind[0]*alphaDrSurf[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*fUpwind[5]*alphaDrSurf[6]+0.25*alphaDrSurf[2]*fUpwind[3]+0.25*fUpwind[2]*alphaDrSurf[3]; 
  Ghat[8] = 0.25*alphaDrSurf[24]*fUpwind[46]+0.25*alphaDrSurf[22]*fUpwind[45]+0.223606797749979*alphaDrSurf[15]*fUpwind[44]+0.25*alphaDrSurf[34]*fUpwind[40]+0.25*alphaDrSurf[13]*fUpwind[39]+0.25*alphaDrSurf[33]*fUpwind[38]+0.223606797749979*alphaDrSurf[6]*fUpwind[37]+0.25*alphaDrSurf[12]*fUpwind[36]+0.223606797749979*alphaDrSurf[5]*fUpwind[35]+0.223606797749979*fUpwind[31]*alphaDrSurf[32]+0.25*alphaDrSurf[7]*fUpwind[31]+0.25*alphaDrSurf[23]*fUpwind[27]+0.25*alphaDrSurf[20]*fUpwind[26]+0.223606797749979*alphaDrSurf[1]*fUpwind[25]+0.223606797749979*fUpwind[17]*alphaDrSurf[21]+0.223606797749979*fUpwind[16]*alphaDrSurf[19]+0.25*alphaDrSurf[15]*fUpwind[18]+0.25*alphaDrSurf[3]*fUpwind[17]+0.25*alphaDrSurf[2]*fUpwind[16]+0.223606797749979*fUpwind[8]*alphaDrSurf[11]+0.25*alphaDrSurf[6]*fUpwind[10]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  Ghat[9] = 0.25*alphaDrSurf[23]*fUpwind[46]+0.223606797749979*alphaDrSurf[15]*fUpwind[45]+0.25*alphaDrSurf[21]*fUpwind[44]+0.25*alphaDrSurf[13]*fUpwind[40]+0.25*alphaDrSurf[34]*fUpwind[39]+0.223606797749979*alphaDrSurf[7]*fUpwind[38]+0.25*alphaDrSurf[32]*fUpwind[37]+0.223606797749979*alphaDrSurf[5]*fUpwind[36]+0.25*alphaDrSurf[11]*fUpwind[35]+0.223606797749979*fUpwind[31]*alphaDrSurf[33]+0.25*alphaDrSurf[6]*fUpwind[31]+0.25*alphaDrSurf[24]*fUpwind[27]+0.223606797749979*alphaDrSurf[2]*fUpwind[26]+0.25*alphaDrSurf[19]*fUpwind[25]+0.223606797749979*fUpwind[18]*alphaDrSurf[22]+0.223606797749979*fUpwind[16]*alphaDrSurf[20]+0.25*alphaDrSurf[3]*fUpwind[18]+0.25*alphaDrSurf[15]*fUpwind[17]+0.25*alphaDrSurf[1]*fUpwind[16]+0.223606797749979*fUpwind[9]*alphaDrSurf[12]+0.25*alphaDrSurf[7]*fUpwind[10]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  Ghat[10] = 0.223606797749979*alphaDrSurf[15]*fUpwind[46]+0.25*alphaDrSurf[20]*fUpwind[45]+0.25*alphaDrSurf[19]*fUpwind[44]+0.223606797749979*alphaDrSurf[7]*fUpwind[40]+0.223606797749979*alphaDrSurf[6]*fUpwind[39]+0.25*alphaDrSurf[12]*fUpwind[38]+0.25*alphaDrSurf[11]*fUpwind[37]+0.25*alphaDrSurf[33]*fUpwind[36]+0.25*alphaDrSurf[32]*fUpwind[35]+0.223606797749979*fUpwind[31]*alphaDrSurf[34]+0.25*alphaDrSurf[5]*fUpwind[31]+0.223606797749979*alphaDrSurf[3]*fUpwind[27]+0.25*alphaDrSurf[22]*fUpwind[26]+0.25*alphaDrSurf[21]*fUpwind[25]+0.223606797749979*fUpwind[18]*alphaDrSurf[24]+0.223606797749979*fUpwind[17]*alphaDrSurf[23]+0.25*alphaDrSurf[2]*fUpwind[18]+0.25*alphaDrSurf[1]*fUpwind[17]+0.25*alphaDrSurf[15]*fUpwind[16]+0.223606797749979*fUpwind[10]*alphaDrSurf[13]+0.25*alphaDrSurf[0]*fUpwind[10]+0.25*alphaDrSurf[7]*fUpwind[9]+0.25*alphaDrSurf[6]*fUpwind[8]+0.25*alphaDrSurf[3]*fUpwind[4]; 
  Ghat[11] = 0.223606797749979*alphaDrSurf[34]*fUpwind[34]+0.223606797749979*alphaDrSurf[33]*fUpwind[33]+0.159719141249985*alphaDrSurf[32]*fUpwind[32]+0.25*alphaDrSurf[7]*fUpwind[32]+0.25*fUpwind[7]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[23]*fUpwind[23]+0.159719141249985*alphaDrSurf[21]*fUpwind[21]+0.2500000000000001*alphaDrSurf[3]*fUpwind[21]+0.2500000000000001*fUpwind[3]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[20]*fUpwind[20]+0.159719141249985*alphaDrSurf[19]*fUpwind[19]+0.2500000000000001*alphaDrSurf[2]*fUpwind[19]+0.2500000000000001*fUpwind[2]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[15]*fUpwind[15]+0.159719141249985*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[0]*fUpwind[11]+0.25*fUpwind[0]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[6]*fUpwind[6]+0.223606797749979*alphaDrSurf[5]*fUpwind[5]+0.223606797749979*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[12] = 0.223606797749979*alphaDrSurf[34]*fUpwind[34]+0.159719141249985*alphaDrSurf[33]*fUpwind[33]+0.25*alphaDrSurf[6]*fUpwind[33]+0.25*fUpwind[6]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[32]*fUpwind[32]+0.223606797749979*alphaDrSurf[24]*fUpwind[24]+0.159719141249985*alphaDrSurf[22]*fUpwind[22]+0.2500000000000001*alphaDrSurf[3]*fUpwind[22]+0.2500000000000001*fUpwind[3]*alphaDrSurf[22]+0.159719141249985*alphaDrSurf[20]*fUpwind[20]+0.2500000000000001*alphaDrSurf[1]*fUpwind[20]+0.2500000000000001*fUpwind[1]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[19]*fUpwind[19]+0.223606797749979*alphaDrSurf[15]*fUpwind[15]+0.159719141249985*alphaDrSurf[12]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[12]+0.25*fUpwind[0]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[7]*fUpwind[7]+0.223606797749979*alphaDrSurf[5]*fUpwind[5]+0.223606797749979*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[13] = 0.159719141249985*alphaDrSurf[34]*fUpwind[34]+0.25*alphaDrSurf[5]*fUpwind[34]+0.25*fUpwind[5]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[33]*fUpwind[33]+0.223606797749979*alphaDrSurf[32]*fUpwind[32]+0.159719141249985*alphaDrSurf[24]*fUpwind[24]+0.2500000000000001*alphaDrSurf[2]*fUpwind[24]+0.2500000000000001*fUpwind[2]*alphaDrSurf[24]+0.159719141249985*alphaDrSurf[23]*fUpwind[23]+0.2500000000000001*alphaDrSurf[1]*fUpwind[23]+0.2500000000000001*fUpwind[1]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[22]*fUpwind[22]+0.223606797749979*alphaDrSurf[21]*fUpwind[21]+0.223606797749979*alphaDrSurf[15]*fUpwind[15]+0.159719141249985*alphaDrSurf[13]*fUpwind[13]+0.25*alphaDrSurf[0]*fUpwind[13]+0.25*fUpwind[0]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[7]*fUpwind[7]+0.223606797749979*alphaDrSurf[6]*fUpwind[6]+0.223606797749979*alphaDrSurf[3]*fUpwind[3]; 
  Ghat[14] = 0.2500000000000001*alphaDrSurf[15]*fUpwind[47]+0.25*alphaDrSurf[7]*fUpwind[43]+0.25*alphaDrSurf[6]*fUpwind[42]+0.25*alphaDrSurf[5]*fUpwind[41]+0.2500000000000001*alphaDrSurf[3]*fUpwind[30]+0.2500000000000001*alphaDrSurf[2]*fUpwind[29]+0.2500000000000001*alphaDrSurf[1]*fUpwind[28]+0.25*alphaDrSurf[0]*fUpwind[14]; 
  Ghat[15] = 0.2*alphaDrSurf[22]*fUpwind[34]+0.2*alphaDrSurf[21]*fUpwind[34]+0.223606797749979*alphaDrSurf[3]*fUpwind[34]+0.2*fUpwind[22]*alphaDrSurf[34]+0.2*fUpwind[21]*alphaDrSurf[34]+0.223606797749979*fUpwind[3]*alphaDrSurf[34]+0.2*alphaDrSurf[24]*fUpwind[33]+0.2*alphaDrSurf[19]*fUpwind[33]+0.223606797749979*alphaDrSurf[2]*fUpwind[33]+0.2*fUpwind[24]*alphaDrSurf[33]+0.2*fUpwind[19]*alphaDrSurf[33]+0.223606797749979*fUpwind[2]*alphaDrSurf[33]+0.2*alphaDrSurf[23]*fUpwind[32]+0.2*alphaDrSurf[20]*fUpwind[32]+0.223606797749979*alphaDrSurf[1]*fUpwind[32]+0.2*fUpwind[23]*alphaDrSurf[32]+0.2*fUpwind[20]*alphaDrSurf[32]+0.223606797749979*fUpwind[1]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[6]*fUpwind[24]+0.223606797749979*fUpwind[6]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[7]*fUpwind[23]+0.223606797749979*fUpwind[7]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[5]*fUpwind[22]+0.223606797749979*fUpwind[5]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[5]*fUpwind[21]+0.223606797749979*fUpwind[5]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[7]*fUpwind[20]+0.223606797749979*fUpwind[7]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[6]*fUpwind[19]+0.223606797749979*fUpwind[6]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[13]*fUpwind[15]+0.223606797749979*alphaDrSurf[12]*fUpwind[15]+0.223606797749979*alphaDrSurf[11]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[15]+0.223606797749979*fUpwind[13]*alphaDrSurf[15]+0.223606797749979*fUpwind[12]*alphaDrSurf[15]+0.223606797749979*fUpwind[11]*alphaDrSurf[15]+0.25*fUpwind[0]*alphaDrSurf[15]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*fUpwind[1]*alphaDrSurf[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[2]*alphaDrSurf[6]+0.25*alphaDrSurf[3]*fUpwind[5]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  Ghat[16] = 0.2500000000000001*alphaDrSurf[13]*fUpwind[46]+0.2*alphaDrSurf[32]*fUpwind[45]+0.223606797749979*alphaDrSurf[7]*fUpwind[45]+0.2*alphaDrSurf[33]*fUpwind[44]+0.223606797749979*alphaDrSurf[6]*fUpwind[44]+0.2500000000000001*alphaDrSurf[23]*fUpwind[40]+0.2500000000000001*alphaDrSurf[24]*fUpwind[39]+0.223606797749979*alphaDrSurf[15]*fUpwind[38]+0.223606797749979*alphaDrSurf[15]*fUpwind[37]+0.2*alphaDrSurf[19]*fUpwind[36]+0.223606797749979*alphaDrSurf[2]*fUpwind[36]+0.2*alphaDrSurf[20]*fUpwind[35]+0.223606797749979*alphaDrSurf[1]*fUpwind[35]+0.2500000000000001*fUpwind[27]*alphaDrSurf[34]+0.223606797749979*fUpwind[18]*alphaDrSurf[33]+0.223606797749979*fUpwind[17]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[22]*fUpwind[31]+0.223606797749979*alphaDrSurf[21]*fUpwind[31]+0.25*alphaDrSurf[3]*fUpwind[31]+0.223606797749979*alphaDrSurf[5]*fUpwind[26]+0.223606797749979*alphaDrSurf[5]*fUpwind[25]+0.223606797749979*fUpwind[9]*alphaDrSurf[20]+0.223606797749979*fUpwind[8]*alphaDrSurf[19]+0.25*alphaDrSurf[6]*fUpwind[18]+0.25*alphaDrSurf[7]*fUpwind[17]+0.223606797749979*alphaDrSurf[12]*fUpwind[16]+0.223606797749979*alphaDrSurf[11]*fUpwind[16]+0.25*alphaDrSurf[0]*fUpwind[16]+0.25*fUpwind[10]*alphaDrSurf[15]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  Ghat[17] = 0.2*alphaDrSurf[32]*fUpwind[46]+0.223606797749979*alphaDrSurf[7]*fUpwind[46]+0.2500000000000001*alphaDrSurf[12]*fUpwind[45]+0.2*alphaDrSurf[34]*fUpwind[44]+0.223606797749979*alphaDrSurf[5]*fUpwind[44]+0.223606797749979*alphaDrSurf[15]*fUpwind[40]+0.2*alphaDrSurf[21]*fUpwind[39]+0.223606797749979*alphaDrSurf[3]*fUpwind[39]+0.2500000000000001*alphaDrSurf[20]*fUpwind[38]+0.2*alphaDrSurf[23]*fUpwind[37]+0.223606797749979*alphaDrSurf[1]*fUpwind[37]+0.2500000000000001*alphaDrSurf[22]*fUpwind[36]+0.223606797749979*alphaDrSurf[15]*fUpwind[35]+0.223606797749979*fUpwind[18]*alphaDrSurf[34]+0.2500000000000001*fUpwind[26]*alphaDrSurf[33]+0.223606797749979*fUpwind[16]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[24]*fUpwind[31]+0.223606797749979*alphaDrSurf[19]*fUpwind[31]+0.25*alphaDrSurf[2]*fUpwind[31]+0.223606797749979*alphaDrSurf[6]*fUpwind[27]+0.223606797749979*alphaDrSurf[6]*fUpwind[25]+0.223606797749979*fUpwind[10]*alphaDrSurf[23]+0.223606797749979*fUpwind[8]*alphaDrSurf[21]+0.25*alphaDrSurf[5]*fUpwind[18]+0.223606797749979*alphaDrSurf[13]*fUpwind[17]+0.223606797749979*alphaDrSurf[11]*fUpwind[17]+0.25*alphaDrSurf[0]*fUpwind[17]+0.25*alphaDrSurf[7]*fUpwind[16]+0.25*fUpwind[9]*alphaDrSurf[15]+0.25*alphaDrSurf[1]*fUpwind[10]+0.25*alphaDrSurf[3]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[6]; 
  Ghat[18] = 0.2*alphaDrSurf[33]*fUpwind[46]+0.223606797749979*alphaDrSurf[6]*fUpwind[46]+0.2*alphaDrSurf[34]*fUpwind[45]+0.223606797749979*alphaDrSurf[5]*fUpwind[45]+0.2500000000000001*alphaDrSurf[11]*fUpwind[44]+0.2*alphaDrSurf[22]*fUpwind[40]+0.223606797749979*alphaDrSurf[3]*fUpwind[40]+0.223606797749979*alphaDrSurf[15]*fUpwind[39]+0.2*alphaDrSurf[24]*fUpwind[38]+0.223606797749979*alphaDrSurf[2]*fUpwind[38]+0.2500000000000001*alphaDrSurf[19]*fUpwind[37]+0.223606797749979*alphaDrSurf[15]*fUpwind[36]+0.2500000000000001*alphaDrSurf[21]*fUpwind[35]+0.223606797749979*fUpwind[17]*alphaDrSurf[34]+0.223606797749979*fUpwind[16]*alphaDrSurf[33]+0.2500000000000001*fUpwind[25]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[23]*fUpwind[31]+0.223606797749979*alphaDrSurf[20]*fUpwind[31]+0.25*alphaDrSurf[1]*fUpwind[31]+0.223606797749979*alphaDrSurf[7]*fUpwind[27]+0.223606797749979*alphaDrSurf[7]*fUpwind[26]+0.223606797749979*fUpwind[10]*alphaDrSurf[24]+0.223606797749979*fUpwind[9]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[13]*fUpwind[18]+0.223606797749979*alphaDrSurf[12]*fUpwind[18]+0.25*alphaDrSurf[0]*fUpwind[18]+0.25*alphaDrSurf[5]*fUpwind[17]+0.25*alphaDrSurf[6]*fUpwind[16]+0.25*fUpwind[8]*alphaDrSurf[15]+0.25*alphaDrSurf[2]*fUpwind[10]+0.25*alphaDrSurf[3]*fUpwind[9]+0.25*fUpwind[4]*alphaDrSurf[7]; 
  Ghat[19] = 0.223606797749979*alphaDrSurf[23]*fUpwind[34]+0.223606797749979*fUpwind[23]*alphaDrSurf[34]+0.2*alphaDrSurf[15]*fUpwind[33]+0.2*fUpwind[15]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[22]*fUpwind[32]+0.159719141249985*alphaDrSurf[21]*fUpwind[32]+0.2500000000000001*alphaDrSurf[3]*fUpwind[32]+0.223606797749979*fUpwind[22]*alphaDrSurf[32]+0.159719141249985*fUpwind[21]*alphaDrSurf[32]+0.2500000000000001*fUpwind[3]*alphaDrSurf[32]+0.25*alphaDrSurf[7]*fUpwind[21]+0.25*fUpwind[7]*alphaDrSurf[21]+0.2*alphaDrSurf[5]*fUpwind[20]+0.2*fUpwind[5]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[12]*fUpwind[19]+0.159719141249985*alphaDrSurf[11]*fUpwind[19]+0.25*alphaDrSurf[0]*fUpwind[19]+0.223606797749979*fUpwind[12]*alphaDrSurf[19]+0.159719141249985*fUpwind[11]*alphaDrSurf[19]+0.25*fUpwind[0]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[6]*fUpwind[15]+0.223606797749979*fUpwind[6]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[2]*fUpwind[11]+0.2500000000000001*fUpwind[2]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[5]+0.223606797749979*fUpwind[1]*alphaDrSurf[5]; 
  Ghat[20] = 0.223606797749979*alphaDrSurf[24]*fUpwind[34]+0.223606797749979*fUpwind[24]*alphaDrSurf[34]+0.159719141249985*alphaDrSurf[22]*fUpwind[33]+0.223606797749979*alphaDrSurf[21]*fUpwind[33]+0.2500000000000001*alphaDrSurf[3]*fUpwind[33]+0.159719141249985*fUpwind[22]*alphaDrSurf[33]+0.223606797749979*fUpwind[21]*alphaDrSurf[33]+0.2500000000000001*fUpwind[3]*alphaDrSurf[33]+0.2*alphaDrSurf[15]*fUpwind[32]+0.2*fUpwind[15]*alphaDrSurf[32]+0.25*alphaDrSurf[6]*fUpwind[22]+0.25*fUpwind[6]*alphaDrSurf[22]+0.159719141249985*alphaDrSurf[12]*fUpwind[20]+0.223606797749979*alphaDrSurf[11]*fUpwind[20]+0.25*alphaDrSurf[0]*fUpwind[20]+0.159719141249985*fUpwind[12]*alphaDrSurf[20]+0.223606797749979*fUpwind[11]*alphaDrSurf[20]+0.25*fUpwind[0]*alphaDrSurf[20]+0.2*alphaDrSurf[5]*fUpwind[19]+0.2*fUpwind[5]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[7]*fUpwind[15]+0.223606797749979*fUpwind[7]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[1]*fUpwind[12]+0.2500000000000001*fUpwind[1]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[5]+0.223606797749979*fUpwind[2]*alphaDrSurf[5]; 
  Ghat[21] = 0.2*alphaDrSurf[15]*fUpwind[34]+0.2*fUpwind[15]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[20]*fUpwind[33]+0.223606797749979*fUpwind[20]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[24]*fUpwind[32]+0.159719141249985*alphaDrSurf[19]*fUpwind[32]+0.2500000000000001*alphaDrSurf[2]*fUpwind[32]+0.223606797749979*fUpwind[24]*alphaDrSurf[32]+0.159719141249985*fUpwind[19]*alphaDrSurf[32]+0.2500000000000001*fUpwind[2]*alphaDrSurf[32]+0.2*alphaDrSurf[6]*fUpwind[23]+0.2*fUpwind[6]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[13]*fUpwind[21]+0.159719141249985*alphaDrSurf[11]*fUpwind[21]+0.25*alphaDrSurf[0]*fUpwind[21]+0.223606797749979*fUpwind[13]*alphaDrSurf[21]+0.159719141249985*fUpwind[11]*alphaDrSurf[21]+0.25*fUpwind[0]*alphaDrSurf[21]+0.25*alphaDrSurf[7]*fUpwind[19]+0.25*fUpwind[7]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[5]*fUpwind[15]+0.223606797749979*fUpwind[5]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[3]*fUpwind[11]+0.2500000000000001*fUpwind[3]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[6]+0.223606797749979*fUpwind[1]*alphaDrSurf[6]; 
  Ghat[22] = 0.2*alphaDrSurf[15]*fUpwind[34]+0.2*fUpwind[15]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[23]*fUpwind[33]+0.159719141249985*alphaDrSurf[20]*fUpwind[33]+0.2500000000000001*alphaDrSurf[1]*fUpwind[33]+0.223606797749979*fUpwind[23]*alphaDrSurf[33]+0.159719141249985*fUpwind[20]*alphaDrSurf[33]+0.2500000000000001*fUpwind[1]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[19]*fUpwind[32]+0.223606797749979*fUpwind[19]*alphaDrSurf[32]+0.2*alphaDrSurf[7]*fUpwind[24]+0.2*fUpwind[7]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[13]*fUpwind[22]+0.159719141249985*alphaDrSurf[12]*fUpwind[22]+0.25*alphaDrSurf[0]*fUpwind[22]+0.223606797749979*fUpwind[13]*alphaDrSurf[22]+0.159719141249985*fUpwind[12]*alphaDrSurf[22]+0.25*fUpwind[0]*alphaDrSurf[22]+0.25*alphaDrSurf[6]*fUpwind[20]+0.25*fUpwind[6]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[15]+0.223606797749979*fUpwind[5]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[3]*fUpwind[12]+0.2500000000000001*fUpwind[3]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[7]+0.223606797749979*fUpwind[2]*alphaDrSurf[7]; 
  Ghat[23] = 0.159719141249985*alphaDrSurf[24]*fUpwind[34]+0.223606797749979*alphaDrSurf[19]*fUpwind[34]+0.2500000000000001*alphaDrSurf[2]*fUpwind[34]+0.159719141249985*fUpwind[24]*alphaDrSurf[34]+0.223606797749979*fUpwind[19]*alphaDrSurf[34]+0.2500000000000001*fUpwind[2]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[22]*fUpwind[33]+0.223606797749979*fUpwind[22]*alphaDrSurf[33]+0.2*alphaDrSurf[15]*fUpwind[32]+0.2*fUpwind[15]*alphaDrSurf[32]+0.25*alphaDrSurf[5]*fUpwind[24]+0.25*fUpwind[5]*alphaDrSurf[24]+0.159719141249985*alphaDrSurf[13]*fUpwind[23]+0.223606797749979*alphaDrSurf[11]*fUpwind[23]+0.25*alphaDrSurf[0]*fUpwind[23]+0.159719141249985*fUpwind[13]*alphaDrSurf[23]+0.223606797749979*fUpwind[11]*alphaDrSurf[23]+0.25*fUpwind[0]*alphaDrSurf[23]+0.2*alphaDrSurf[6]*fUpwind[21]+0.2*fUpwind[6]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[7]*fUpwind[15]+0.223606797749979*fUpwind[7]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[1]*fUpwind[13]+0.2500000000000001*fUpwind[1]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[3]*fUpwind[6]+0.223606797749979*fUpwind[3]*alphaDrSurf[6]; 
  Ghat[24] = 0.159719141249985*alphaDrSurf[23]*fUpwind[34]+0.223606797749979*alphaDrSurf[20]*fUpwind[34]+0.2500000000000001*alphaDrSurf[1]*fUpwind[34]+0.159719141249985*fUpwind[23]*alphaDrSurf[34]+0.223606797749979*fUpwind[20]*alphaDrSurf[34]+0.2500000000000001*fUpwind[1]*alphaDrSurf[34]+0.2*alphaDrSurf[15]*fUpwind[33]+0.2*fUpwind[15]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[21]*fUpwind[32]+0.223606797749979*fUpwind[21]*alphaDrSurf[32]+0.159719141249985*alphaDrSurf[13]*fUpwind[24]+0.223606797749979*alphaDrSurf[12]*fUpwind[24]+0.25*alphaDrSurf[0]*fUpwind[24]+0.159719141249985*fUpwind[13]*alphaDrSurf[24]+0.223606797749979*fUpwind[12]*alphaDrSurf[24]+0.25*fUpwind[0]*alphaDrSurf[24]+0.25*alphaDrSurf[5]*fUpwind[23]+0.25*fUpwind[5]*alphaDrSurf[23]+0.2*alphaDrSurf[7]*fUpwind[22]+0.2*fUpwind[7]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[6]*fUpwind[15]+0.223606797749979*fUpwind[6]*alphaDrSurf[15]+0.2500000000000001*alphaDrSurf[2]*fUpwind[13]+0.2500000000000001*fUpwind[2]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[3]*fUpwind[7]+0.223606797749979*fUpwind[3]*alphaDrSurf[7]; 
  Ghat[25] = 0.223606797749979*alphaDrSurf[34]*fUpwind[46]+0.223606797749979*alphaDrSurf[33]*fUpwind[45]+0.159719141249985*alphaDrSurf[32]*fUpwind[44]+0.25*alphaDrSurf[7]*fUpwind[44]+0.223606797749979*alphaDrSurf[23]*fUpwind[39]+0.159719141249985*alphaDrSurf[21]*fUpwind[37]+0.2500000000000001*alphaDrSurf[3]*fUpwind[37]+0.223606797749979*alphaDrSurf[20]*fUpwind[36]+0.159719141249985*alphaDrSurf[19]*fUpwind[35]+0.2500000000000001*alphaDrSurf[2]*fUpwind[35]+0.2500000000000001*fUpwind[18]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[15]*fUpwind[31]+0.159719141249985*alphaDrSurf[11]*fUpwind[25]+0.25*alphaDrSurf[0]*fUpwind[25]+0.25*fUpwind[10]*alphaDrSurf[21]+0.25*fUpwind[9]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[6]*fUpwind[17]+0.223606797749979*alphaDrSurf[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[8]; 
  Ghat[26] = 0.223606797749979*alphaDrSurf[34]*fUpwind[46]+0.159719141249985*alphaDrSurf[33]*fUpwind[45]+0.25*alphaDrSurf[6]*fUpwind[45]+0.223606797749979*alphaDrSurf[32]*fUpwind[44]+0.223606797749979*alphaDrSurf[24]*fUpwind[40]+0.159719141249985*alphaDrSurf[22]*fUpwind[38]+0.2500000000000001*alphaDrSurf[3]*fUpwind[38]+0.159719141249985*alphaDrSurf[20]*fUpwind[36]+0.2500000000000001*alphaDrSurf[1]*fUpwind[36]+0.223606797749979*alphaDrSurf[19]*fUpwind[35]+0.2500000000000001*fUpwind[17]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[15]*fUpwind[31]+0.159719141249985*alphaDrSurf[12]*fUpwind[26]+0.25*alphaDrSurf[0]*fUpwind[26]+0.25*fUpwind[10]*alphaDrSurf[22]+0.25*fUpwind[8]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[7]*fUpwind[18]+0.223606797749979*alphaDrSurf[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[9]; 
  Ghat[27] = 0.159719141249985*alphaDrSurf[34]*fUpwind[46]+0.25*alphaDrSurf[5]*fUpwind[46]+0.223606797749979*alphaDrSurf[33]*fUpwind[45]+0.223606797749979*alphaDrSurf[32]*fUpwind[44]+0.159719141249985*alphaDrSurf[24]*fUpwind[40]+0.2500000000000001*alphaDrSurf[2]*fUpwind[40]+0.159719141249985*alphaDrSurf[23]*fUpwind[39]+0.2500000000000001*alphaDrSurf[1]*fUpwind[39]+0.223606797749979*alphaDrSurf[22]*fUpwind[38]+0.223606797749979*alphaDrSurf[21]*fUpwind[37]+0.2500000000000001*fUpwind[16]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[15]*fUpwind[31]+0.159719141249985*alphaDrSurf[13]*fUpwind[27]+0.25*alphaDrSurf[0]*fUpwind[27]+0.25*fUpwind[9]*alphaDrSurf[24]+0.25*fUpwind[8]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[7]*fUpwind[18]+0.223606797749979*alphaDrSurf[6]*fUpwind[17]+0.2500000000000001*fUpwind[4]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[3]*fUpwind[10]; 
  Ghat[28] = 0.223606797749979*alphaDrSurf[32]*fUpwind[47]+0.25*alphaDrSurf[7]*fUpwind[47]+0.2500000000000001*alphaDrSurf[15]*fUpwind[43]+0.223606797749979*alphaDrSurf[21]*fUpwind[42]+0.2500000000000001*alphaDrSurf[3]*fUpwind[42]+0.223606797749979*alphaDrSurf[19]*fUpwind[41]+0.2500000000000001*alphaDrSurf[2]*fUpwind[41]+0.25*alphaDrSurf[6]*fUpwind[30]+0.25*alphaDrSurf[5]*fUpwind[29]+0.223606797749979*alphaDrSurf[11]*fUpwind[28]+0.25*alphaDrSurf[0]*fUpwind[28]+0.2500000000000001*alphaDrSurf[1]*fUpwind[14]; 
  Ghat[29] = 0.223606797749979*alphaDrSurf[33]*fUpwind[47]+0.25*alphaDrSurf[6]*fUpwind[47]+0.223606797749979*alphaDrSurf[22]*fUpwind[43]+0.2500000000000001*alphaDrSurf[3]*fUpwind[43]+0.2500000000000001*alphaDrSurf[15]*fUpwind[42]+0.223606797749979*alphaDrSurf[20]*fUpwind[41]+0.2500000000000001*alphaDrSurf[1]*fUpwind[41]+0.25*alphaDrSurf[7]*fUpwind[30]+0.223606797749979*alphaDrSurf[12]*fUpwind[29]+0.25*alphaDrSurf[0]*fUpwind[29]+0.25*alphaDrSurf[5]*fUpwind[28]+0.2500000000000001*alphaDrSurf[2]*fUpwind[14]; 
  Ghat[30] = 0.223606797749979*alphaDrSurf[34]*fUpwind[47]+0.25*alphaDrSurf[5]*fUpwind[47]+0.223606797749979*alphaDrSurf[24]*fUpwind[43]+0.2500000000000001*alphaDrSurf[2]*fUpwind[43]+0.223606797749979*alphaDrSurf[23]*fUpwind[42]+0.2500000000000001*alphaDrSurf[1]*fUpwind[42]+0.2500000000000001*alphaDrSurf[15]*fUpwind[41]+0.223606797749979*alphaDrSurf[13]*fUpwind[30]+0.25*alphaDrSurf[0]*fUpwind[30]+0.25*alphaDrSurf[7]*fUpwind[29]+0.25*alphaDrSurf[6]*fUpwind[28]+0.2500000000000001*alphaDrSurf[3]*fUpwind[14]; 
  Ghat[31] = 0.2*alphaDrSurf[22]*fUpwind[46]+0.2*alphaDrSurf[21]*fUpwind[46]+0.223606797749979*alphaDrSurf[3]*fUpwind[46]+0.2*alphaDrSurf[24]*fUpwind[45]+0.2*alphaDrSurf[19]*fUpwind[45]+0.223606797749979*alphaDrSurf[2]*fUpwind[45]+0.2*alphaDrSurf[23]*fUpwind[44]+0.2*alphaDrSurf[20]*fUpwind[44]+0.223606797749979*alphaDrSurf[1]*fUpwind[44]+0.2*alphaDrSurf[33]*fUpwind[40]+0.223606797749979*alphaDrSurf[6]*fUpwind[40]+0.2*alphaDrSurf[32]*fUpwind[39]+0.223606797749979*alphaDrSurf[7]*fUpwind[39]+0.2*alphaDrSurf[34]*fUpwind[38]+0.223606797749979*alphaDrSurf[5]*fUpwind[38]+0.2*alphaDrSurf[34]*fUpwind[37]+0.223606797749979*alphaDrSurf[5]*fUpwind[37]+0.2*alphaDrSurf[32]*fUpwind[36]+0.223606797749979*alphaDrSurf[7]*fUpwind[36]+0.2*alphaDrSurf[33]*fUpwind[35]+0.223606797749979*alphaDrSurf[6]*fUpwind[35]+0.223606797749979*fUpwind[10]*alphaDrSurf[34]+0.223606797749979*fUpwind[9]*alphaDrSurf[33]+0.223606797749979*fUpwind[8]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[13]*fUpwind[31]+0.223606797749979*alphaDrSurf[12]*fUpwind[31]+0.223606797749979*alphaDrSurf[11]*fUpwind[31]+0.25*alphaDrSurf[0]*fUpwind[31]+0.223606797749979*alphaDrSurf[15]*fUpwind[27]+0.223606797749979*alphaDrSurf[15]*fUpwind[26]+0.223606797749979*alphaDrSurf[15]*fUpwind[25]+0.223606797749979*fUpwind[17]*alphaDrSurf[24]+0.223606797749979*fUpwind[18]*alphaDrSurf[23]+0.223606797749979*fUpwind[16]*alphaDrSurf[22]+0.223606797749979*fUpwind[16]*alphaDrSurf[21]+0.223606797749979*fUpwind[18]*alphaDrSurf[20]+0.223606797749979*fUpwind[17]*alphaDrSurf[19]+0.25*alphaDrSurf[1]*fUpwind[18]+0.25*alphaDrSurf[2]*fUpwind[17]+0.25*alphaDrSurf[3]*fUpwind[16]+0.25*fUpwind[4]*alphaDrSurf[15]+0.25*alphaDrSurf[5]*fUpwind[10]+0.25*alphaDrSurf[6]*fUpwind[9]+0.25*alphaDrSurf[7]*fUpwind[8]; 
  Ghat[32] = 0.1788854381999831*alphaDrSurf[33]*fUpwind[34]+0.2*alphaDrSurf[6]*fUpwind[34]+0.1788854381999831*fUpwind[33]*alphaDrSurf[34]+0.2*fUpwind[6]*alphaDrSurf[34]+0.2*alphaDrSurf[5]*fUpwind[33]+0.2*fUpwind[5]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[13]*fUpwind[32]+0.223606797749979*alphaDrSurf[12]*fUpwind[32]+0.159719141249985*alphaDrSurf[11]*fUpwind[32]+0.25*alphaDrSurf[0]*fUpwind[32]+0.223606797749979*fUpwind[13]*alphaDrSurf[32]+0.223606797749979*fUpwind[12]*alphaDrSurf[32]+0.159719141249985*fUpwind[11]*alphaDrSurf[32]+0.25*fUpwind[0]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[21]*fUpwind[24]+0.223606797749979*fUpwind[21]*alphaDrSurf[24]+0.2*alphaDrSurf[15]*fUpwind[23]+0.2*fUpwind[15]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[19]*fUpwind[22]+0.223606797749979*fUpwind[19]*alphaDrSurf[22]+0.159719141249985*alphaDrSurf[19]*fUpwind[21]+0.2500000000000001*alphaDrSurf[2]*fUpwind[21]+0.159719141249985*fUpwind[19]*alphaDrSurf[21]+0.2500000000000001*fUpwind[2]*alphaDrSurf[21]+0.2*alphaDrSurf[15]*fUpwind[20]+0.2*fUpwind[15]*alphaDrSurf[20]+0.2500000000000001*alphaDrSurf[3]*fUpwind[19]+0.2500000000000001*fUpwind[3]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[15]+0.223606797749979*fUpwind[1]*alphaDrSurf[15]+0.25*alphaDrSurf[7]*fUpwind[11]+0.25*fUpwind[7]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[6]+0.223606797749979*fUpwind[5]*alphaDrSurf[6]; 
  Ghat[33] = 0.1788854381999831*alphaDrSurf[32]*fUpwind[34]+0.2*alphaDrSurf[7]*fUpwind[34]+0.1788854381999831*fUpwind[32]*alphaDrSurf[34]+0.2*fUpwind[7]*alphaDrSurf[34]+0.223606797749979*alphaDrSurf[13]*fUpwind[33]+0.159719141249985*alphaDrSurf[12]*fUpwind[33]+0.223606797749979*alphaDrSurf[11]*fUpwind[33]+0.25*alphaDrSurf[0]*fUpwind[33]+0.223606797749979*fUpwind[13]*alphaDrSurf[33]+0.159719141249985*fUpwind[12]*alphaDrSurf[33]+0.223606797749979*fUpwind[11]*alphaDrSurf[33]+0.25*fUpwind[0]*alphaDrSurf[33]+0.2*alphaDrSurf[5]*fUpwind[32]+0.2*fUpwind[5]*alphaDrSurf[32]+0.2*alphaDrSurf[15]*fUpwind[24]+0.2*fUpwind[15]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[22]*fUpwind[23]+0.223606797749979*fUpwind[22]*alphaDrSurf[23]+0.159719141249985*alphaDrSurf[20]*fUpwind[22]+0.2500000000000001*alphaDrSurf[1]*fUpwind[22]+0.159719141249985*fUpwind[20]*alphaDrSurf[22]+0.2500000000000001*fUpwind[1]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[20]*fUpwind[21]+0.223606797749979*fUpwind[20]*alphaDrSurf[21]+0.2500000000000001*alphaDrSurf[3]*fUpwind[20]+0.2500000000000001*fUpwind[3]*alphaDrSurf[20]+0.2*alphaDrSurf[15]*fUpwind[19]+0.2*fUpwind[15]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[2]*fUpwind[15]+0.223606797749979*fUpwind[2]*alphaDrSurf[15]+0.25*alphaDrSurf[6]*fUpwind[12]+0.25*fUpwind[6]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[7]+0.223606797749979*fUpwind[5]*alphaDrSurf[7]; 
  Ghat[34] = 0.159719141249985*alphaDrSurf[13]*fUpwind[34]+0.223606797749979*alphaDrSurf[12]*fUpwind[34]+0.223606797749979*alphaDrSurf[11]*fUpwind[34]+0.25*alphaDrSurf[0]*fUpwind[34]+0.159719141249985*fUpwind[13]*alphaDrSurf[34]+0.223606797749979*fUpwind[12]*alphaDrSurf[34]+0.223606797749979*fUpwind[11]*alphaDrSurf[34]+0.25*fUpwind[0]*alphaDrSurf[34]+0.1788854381999831*alphaDrSurf[32]*fUpwind[33]+0.2*alphaDrSurf[7]*fUpwind[33]+0.1788854381999831*fUpwind[32]*alphaDrSurf[33]+0.2*fUpwind[7]*alphaDrSurf[33]+0.2*alphaDrSurf[6]*fUpwind[32]+0.2*fUpwind[6]*alphaDrSurf[32]+0.159719141249985*alphaDrSurf[23]*fUpwind[24]+0.223606797749979*alphaDrSurf[20]*fUpwind[24]+0.2500000000000001*alphaDrSurf[1]*fUpwind[24]+0.159719141249985*fUpwind[23]*alphaDrSurf[24]+0.223606797749979*fUpwind[20]*alphaDrSurf[24]+0.2500000000000001*fUpwind[1]*alphaDrSurf[24]+0.223606797749979*alphaDrSurf[19]*fUpwind[23]+0.2500000000000001*alphaDrSurf[2]*fUpwind[23]+0.223606797749979*fUpwind[19]*alphaDrSurf[23]+0.2500000000000001*fUpwind[2]*alphaDrSurf[23]+0.2*alphaDrSurf[15]*fUpwind[22]+0.2*fUpwind[15]*alphaDrSurf[22]+0.2*alphaDrSurf[15]*fUpwind[21]+0.2*fUpwind[15]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[3]*fUpwind[15]+0.223606797749979*fUpwind[3]*alphaDrSurf[15]+0.25*alphaDrSurf[5]*fUpwind[13]+0.25*fUpwind[5]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[6]*fUpwind[7]+0.223606797749979*fUpwind[6]*alphaDrSurf[7]; 
  Ghat[35] = 0.223606797749979*alphaDrSurf[23]*fUpwind[46]+0.2*alphaDrSurf[15]*fUpwind[45]+0.223606797749979*alphaDrSurf[22]*fUpwind[44]+0.159719141249985*alphaDrSurf[21]*fUpwind[44]+0.2500000000000001*alphaDrSurf[3]*fUpwind[44]+0.223606797749979*alphaDrSurf[34]*fUpwind[39]+0.223606797749979*alphaDrSurf[32]*fUpwind[38]+0.159719141249985*alphaDrSurf[32]*fUpwind[37]+0.25*alphaDrSurf[7]*fUpwind[37]+0.2*alphaDrSurf[5]*fUpwind[36]+0.223606797749979*alphaDrSurf[12]*fUpwind[35]+0.159719141249985*alphaDrSurf[11]*fUpwind[35]+0.25*alphaDrSurf[0]*fUpwind[35]+0.2*fUpwind[31]*alphaDrSurf[33]+0.25*fUpwind[10]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[6]*fUpwind[31]+0.223606797749979*alphaDrSurf[19]*fUpwind[26]+0.159719141249985*alphaDrSurf[19]*fUpwind[25]+0.2500000000000001*alphaDrSurf[2]*fUpwind[25]+0.2500000000000001*fUpwind[18]*alphaDrSurf[21]+0.2*fUpwind[16]*alphaDrSurf[20]+0.2500000000000001*fUpwind[4]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[15]*fUpwind[17]+0.223606797749979*alphaDrSurf[1]*fUpwind[16]+0.25*fUpwind[9]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[8]; 
  Ghat[36] = 0.223606797749979*alphaDrSurf[24]*fUpwind[46]+0.159719141249985*alphaDrSurf[22]*fUpwind[45]+0.223606797749979*alphaDrSurf[21]*fUpwind[45]+0.2500000000000001*alphaDrSurf[3]*fUpwind[45]+0.2*alphaDrSurf[15]*fUpwind[44]+0.223606797749979*alphaDrSurf[34]*fUpwind[40]+0.159719141249985*alphaDrSurf[33]*fUpwind[38]+0.25*alphaDrSurf[6]*fUpwind[38]+0.223606797749979*alphaDrSurf[33]*fUpwind[37]+0.159719141249985*alphaDrSurf[12]*fUpwind[36]+0.223606797749979*alphaDrSurf[11]*fUpwind[36]+0.25*alphaDrSurf[0]*fUpwind[36]+0.2*alphaDrSurf[5]*fUpwind[35]+0.25*fUpwind[10]*alphaDrSurf[33]+0.2*fUpwind[31]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[7]*fUpwind[31]+0.159719141249985*alphaDrSurf[20]*fUpwind[26]+0.2500000000000001*alphaDrSurf[1]*fUpwind[26]+0.223606797749979*alphaDrSurf[20]*fUpwind[25]+0.2500000000000001*fUpwind[17]*alphaDrSurf[22]+0.2500000000000001*fUpwind[4]*alphaDrSurf[20]+0.2*fUpwind[16]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[15]*fUpwind[18]+0.223606797749979*alphaDrSurf[2]*fUpwind[16]+0.25*fUpwind[8]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[9]; 
  Ghat[37] = 0.2*alphaDrSurf[15]*fUpwind[46]+0.223606797749979*alphaDrSurf[20]*fUpwind[45]+0.223606797749979*alphaDrSurf[24]*fUpwind[44]+0.159719141249985*alphaDrSurf[19]*fUpwind[44]+0.2500000000000001*alphaDrSurf[2]*fUpwind[44]+0.223606797749979*alphaDrSurf[32]*fUpwind[40]+0.2*alphaDrSurf[6]*fUpwind[39]+0.223606797749979*alphaDrSurf[13]*fUpwind[37]+0.159719141249985*alphaDrSurf[11]*fUpwind[37]+0.25*alphaDrSurf[0]*fUpwind[37]+0.223606797749979*alphaDrSurf[33]*fUpwind[36]+0.159719141249985*alphaDrSurf[32]*fUpwind[35]+0.25*alphaDrSurf[7]*fUpwind[35]+0.2*fUpwind[31]*alphaDrSurf[34]+0.25*fUpwind[9]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[5]*fUpwind[31]+0.223606797749979*alphaDrSurf[21]*fUpwind[27]+0.159719141249985*alphaDrSurf[21]*fUpwind[25]+0.2500000000000001*alphaDrSurf[3]*fUpwind[25]+0.2*fUpwind[17]*alphaDrSurf[23]+0.2500000000000001*fUpwind[4]*alphaDrSurf[21]+0.2500000000000001*fUpwind[18]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[17]+0.223606797749979*alphaDrSurf[15]*fUpwind[16]+0.25*fUpwind[10]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[6]*fUpwind[8]; 
  Ghat[38] = 0.2*alphaDrSurf[15]*fUpwind[46]+0.223606797749979*alphaDrSurf[23]*fUpwind[45]+0.159719141249985*alphaDrSurf[20]*fUpwind[45]+0.2500000000000001*alphaDrSurf[1]*fUpwind[45]+0.223606797749979*alphaDrSurf[19]*fUpwind[44]+0.2*alphaDrSurf[7]*fUpwind[40]+0.223606797749979*alphaDrSurf[33]*fUpwind[39]+0.223606797749979*alphaDrSurf[13]*fUpwind[38]+0.159719141249985*alphaDrSurf[12]*fUpwind[38]+0.25*alphaDrSurf[0]*fUpwind[38]+0.159719141249985*alphaDrSurf[33]*fUpwind[36]+0.25*alphaDrSurf[6]*fUpwind[36]+0.223606797749979*alphaDrSurf[32]*fUpwind[35]+0.2*fUpwind[31]*alphaDrSurf[34]+0.25*fUpwind[8]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[5]*fUpwind[31]+0.223606797749979*alphaDrSurf[22]*fUpwind[27]+0.159719141249985*alphaDrSurf[22]*fUpwind[26]+0.2500000000000001*alphaDrSurf[3]*fUpwind[26]+0.2*fUpwind[18]*alphaDrSurf[24]+0.2500000000000001*fUpwind[4]*alphaDrSurf[22]+0.2500000000000001*fUpwind[17]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[2]*fUpwind[18]+0.223606797749979*alphaDrSurf[15]*fUpwind[16]+0.25*fUpwind[10]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[7]*fUpwind[9]; 
  Ghat[39] = 0.159719141249985*alphaDrSurf[24]*fUpwind[46]+0.223606797749979*alphaDrSurf[19]*fUpwind[46]+0.2500000000000001*alphaDrSurf[2]*fUpwind[46]+0.223606797749979*alphaDrSurf[22]*fUpwind[45]+0.2*alphaDrSurf[15]*fUpwind[44]+0.159719141249985*alphaDrSurf[34]*fUpwind[40]+0.25*alphaDrSurf[5]*fUpwind[40]+0.159719141249985*alphaDrSurf[13]*fUpwind[39]+0.223606797749979*alphaDrSurf[11]*fUpwind[39]+0.25*alphaDrSurf[0]*fUpwind[39]+0.223606797749979*alphaDrSurf[33]*fUpwind[38]+0.2*alphaDrSurf[6]*fUpwind[37]+0.223606797749979*alphaDrSurf[34]*fUpwind[35]+0.25*fUpwind[9]*alphaDrSurf[34]+0.2*fUpwind[31]*alphaDrSurf[32]+0.223606797749979*alphaDrSurf[7]*fUpwind[31]+0.159719141249985*alphaDrSurf[23]*fUpwind[27]+0.2500000000000001*alphaDrSurf[1]*fUpwind[27]+0.223606797749979*alphaDrSurf[23]*fUpwind[25]+0.2500000000000001*fUpwind[16]*alphaDrSurf[24]+0.2500000000000001*fUpwind[4]*alphaDrSurf[23]+0.2*fUpwind[17]*alphaDrSurf[21]+0.223606797749979*alphaDrSurf[15]*fUpwind[18]+0.223606797749979*alphaDrSurf[3]*fUpwind[17]+0.25*fUpwind[8]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[6]*fUpwind[10]; 
  Ghat[40] = 0.159719141249985*alphaDrSurf[23]*fUpwind[46]+0.223606797749979*alphaDrSurf[20]*fUpwind[46]+0.2500000000000001*alphaDrSurf[1]*fUpwind[46]+0.2*alphaDrSurf[15]*fUpwind[45]+0.223606797749979*alphaDrSurf[21]*fUpwind[44]+0.159719141249985*alphaDrSurf[13]*fUpwind[40]+0.223606797749979*alphaDrSurf[12]*fUpwind[40]+0.25*alphaDrSurf[0]*fUpwind[40]+0.159719141249985*alphaDrSurf[34]*fUpwind[39]+0.25*alphaDrSurf[5]*fUpwind[39]+0.2*alphaDrSurf[7]*fUpwind[38]+0.223606797749979*alphaDrSurf[32]*fUpwind[37]+0.223606797749979*alphaDrSurf[34]*fUpwind[36]+0.25*fUpwind[8]*alphaDrSurf[34]+0.2*fUpwind[31]*alphaDrSurf[33]+0.223606797749979*alphaDrSurf[6]*fUpwind[31]+0.159719141249985*alphaDrSurf[24]*fUpwind[27]+0.2500000000000001*alphaDrSurf[2]*fUpwind[27]+0.223606797749979*alphaDrSurf[24]*fUpwind[26]+0.2500000000000001*fUpwind[4]*alphaDrSurf[24]+0.2500000000000001*fUpwind[16]*alphaDrSurf[23]+0.2*fUpwind[18]*alphaDrSurf[22]+0.223606797749979*alphaDrSurf[3]*fUpwind[18]+0.223606797749979*alphaDrSurf[15]*fUpwind[17]+0.25*fUpwind[9]*alphaDrSurf[13]+0.223606797749979*alphaDrSurf[7]*fUpwind[10]; 
  Ghat[41] = 0.223606797749979*alphaDrSurf[22]*fUpwind[47]+0.223606797749979*alphaDrSurf[21]*fUpwind[47]+0.2500000000000001*alphaDrSurf[3]*fUpwind[47]+0.223606797749979*alphaDrSurf[33]*fUpwind[43]+0.25*alphaDrSurf[6]*fUpwind[43]+0.223606797749979*alphaDrSurf[32]*fUpwind[42]+0.25*alphaDrSurf[7]*fUpwind[42]+0.223606797749979*alphaDrSurf[12]*fUpwind[41]+0.223606797749979*alphaDrSurf[11]*fUpwind[41]+0.25*alphaDrSurf[0]*fUpwind[41]+0.2500000000000001*alphaDrSurf[15]*fUpwind[30]+0.223606797749979*alphaDrSurf[20]*fUpwind[29]+0.2500000000000001*alphaDrSurf[1]*fUpwind[29]+0.223606797749979*alphaDrSurf[19]*fUpwind[28]+0.2500000000000001*alphaDrSurf[2]*fUpwind[28]+0.25*alphaDrSurf[5]*fUpwind[14]; 
  Ghat[42] = 0.223606797749979*alphaDrSurf[24]*fUpwind[47]+0.223606797749979*alphaDrSurf[19]*fUpwind[47]+0.2500000000000001*alphaDrSurf[2]*fUpwind[47]+0.223606797749979*alphaDrSurf[34]*fUpwind[43]+0.25*alphaDrSurf[5]*fUpwind[43]+0.223606797749979*alphaDrSurf[13]*fUpwind[42]+0.223606797749979*alphaDrSurf[11]*fUpwind[42]+0.25*alphaDrSurf[0]*fUpwind[42]+0.223606797749979*alphaDrSurf[32]*fUpwind[41]+0.25*alphaDrSurf[7]*fUpwind[41]+0.223606797749979*alphaDrSurf[23]*fUpwind[30]+0.2500000000000001*alphaDrSurf[1]*fUpwind[30]+0.2500000000000001*alphaDrSurf[15]*fUpwind[29]+0.223606797749979*alphaDrSurf[21]*fUpwind[28]+0.2500000000000001*alphaDrSurf[3]*fUpwind[28]+0.25*alphaDrSurf[6]*fUpwind[14]; 
  Ghat[43] = 0.223606797749979*alphaDrSurf[23]*fUpwind[47]+0.223606797749979*alphaDrSurf[20]*fUpwind[47]+0.2500000000000001*alphaDrSurf[1]*fUpwind[47]+0.223606797749979*alphaDrSurf[13]*fUpwind[43]+0.223606797749979*alphaDrSurf[12]*fUpwind[43]+0.25*alphaDrSurf[0]*fUpwind[43]+0.223606797749979*alphaDrSurf[34]*fUpwind[42]+0.25*alphaDrSurf[5]*fUpwind[42]+0.223606797749979*alphaDrSurf[33]*fUpwind[41]+0.25*alphaDrSurf[6]*fUpwind[41]+0.223606797749979*alphaDrSurf[24]*fUpwind[30]+0.2500000000000001*alphaDrSurf[2]*fUpwind[30]+0.223606797749979*alphaDrSurf[22]*fUpwind[29]+0.2500000000000001*alphaDrSurf[3]*fUpwind[29]+0.2500000000000001*alphaDrSurf[15]*fUpwind[28]+0.25*alphaDrSurf[7]*fUpwind[14]; 
  Ghat[44] = 0.1788854381999831*alphaDrSurf[33]*fUpwind[46]+0.2*alphaDrSurf[6]*fUpwind[46]+0.1788854381999831*alphaDrSurf[34]*fUpwind[45]+0.2*alphaDrSurf[5]*fUpwind[45]+0.223606797749979*alphaDrSurf[13]*fUpwind[44]+0.223606797749979*alphaDrSurf[12]*fUpwind[44]+0.159719141249985*alphaDrSurf[11]*fUpwind[44]+0.25*alphaDrSurf[0]*fUpwind[44]+0.223606797749979*alphaDrSurf[21]*fUpwind[40]+0.2*alphaDrSurf[15]*fUpwind[39]+0.223606797749979*alphaDrSurf[19]*fUpwind[38]+0.223606797749979*alphaDrSurf[24]*fUpwind[37]+0.159719141249985*alphaDrSurf[19]*fUpwind[37]+0.2500000000000001*alphaDrSurf[2]*fUpwind[37]+0.2*alphaDrSurf[15]*fUpwind[36]+0.223606797749979*alphaDrSurf[22]*fUpwind[35]+0.159719141249985*alphaDrSurf[21]*fUpwind[35]+0.2500000000000001*alphaDrSurf[3]*fUpwind[35]+0.2*fUpwind[17]*alphaDrSurf[34]+0.2*fUpwind[16]*alphaDrSurf[33]+0.223606797749979*fUpwind[27]*alphaDrSurf[32]+0.223606797749979*fUpwind[26]*alphaDrSurf[32]+0.159719141249985*fUpwind[25]*alphaDrSurf[32]+0.2500000000000001*fUpwind[4]*alphaDrSurf[32]+0.2*alphaDrSurf[23]*fUpwind[31]+0.2*alphaDrSurf[20]*fUpwind[31]+0.223606797749979*alphaDrSurf[1]*fUpwind[31]+0.25*alphaDrSurf[7]*fUpwind[25]+0.25*fUpwind[9]*alphaDrSurf[21]+0.25*fUpwind[10]*alphaDrSurf[19]+0.2500000000000001*alphaDrSurf[11]*fUpwind[18]+0.223606797749979*alphaDrSurf[5]*fUpwind[17]+0.223606797749979*alphaDrSurf[6]*fUpwind[16]+0.223606797749979*fUpwind[8]*alphaDrSurf[15]; 
  Ghat[45] = 0.1788854381999831*alphaDrSurf[32]*fUpwind[46]+0.2*alphaDrSurf[7]*fUpwind[46]+0.223606797749979*alphaDrSurf[13]*fUpwind[45]+0.159719141249985*alphaDrSurf[12]*fUpwind[45]+0.223606797749979*alphaDrSurf[11]*fUpwind[45]+0.25*alphaDrSurf[0]*fUpwind[45]+0.1788854381999831*alphaDrSurf[34]*fUpwind[44]+0.2*alphaDrSurf[5]*fUpwind[44]+0.2*alphaDrSurf[15]*fUpwind[40]+0.223606797749979*alphaDrSurf[22]*fUpwind[39]+0.223606797749979*alphaDrSurf[23]*fUpwind[38]+0.159719141249985*alphaDrSurf[20]*fUpwind[38]+0.2500000000000001*alphaDrSurf[1]*fUpwind[38]+0.223606797749979*alphaDrSurf[20]*fUpwind[37]+0.159719141249985*alphaDrSurf[22]*fUpwind[36]+0.223606797749979*alphaDrSurf[21]*fUpwind[36]+0.2500000000000001*alphaDrSurf[3]*fUpwind[36]+0.2*alphaDrSurf[15]*fUpwind[35]+0.2*fUpwind[18]*alphaDrSurf[34]+0.223606797749979*fUpwind[27]*alphaDrSurf[33]+0.159719141249985*fUpwind[26]*alphaDrSurf[33]+0.223606797749979*fUpwind[25]*alphaDrSurf[33]+0.2500000000000001*fUpwind[4]*alphaDrSurf[33]+0.2*fUpwind[16]*alphaDrSurf[32]+0.2*alphaDrSurf[24]*fUpwind[31]+0.2*alphaDrSurf[19]*fUpwind[31]+0.223606797749979*alphaDrSurf[2]*fUpwind[31]+0.25*alphaDrSurf[6]*fUpwind[26]+0.25*fUpwind[8]*alphaDrSurf[22]+0.25*fUpwind[10]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[18]+0.2500000000000001*alphaDrSurf[12]*fUpwind[17]+0.223606797749979*alphaDrSurf[7]*fUpwind[16]+0.223606797749979*fUpwind[9]*alphaDrSurf[15]; 
  Ghat[46] = 0.159719141249985*alphaDrSurf[13]*fUpwind[46]+0.223606797749979*alphaDrSurf[12]*fUpwind[46]+0.223606797749979*alphaDrSurf[11]*fUpwind[46]+0.25*alphaDrSurf[0]*fUpwind[46]+0.1788854381999831*alphaDrSurf[32]*fUpwind[45]+0.2*alphaDrSurf[7]*fUpwind[45]+0.1788854381999831*alphaDrSurf[33]*fUpwind[44]+0.2*alphaDrSurf[6]*fUpwind[44]+0.159719141249985*alphaDrSurf[23]*fUpwind[40]+0.223606797749979*alphaDrSurf[20]*fUpwind[40]+0.2500000000000001*alphaDrSurf[1]*fUpwind[40]+0.159719141249985*alphaDrSurf[24]*fUpwind[39]+0.223606797749979*alphaDrSurf[19]*fUpwind[39]+0.2500000000000001*alphaDrSurf[2]*fUpwind[39]+0.2*alphaDrSurf[15]*fUpwind[38]+0.2*alphaDrSurf[15]*fUpwind[37]+0.223606797749979*alphaDrSurf[24]*fUpwind[36]+0.223606797749979*alphaDrSurf[23]*fUpwind[35]+0.159719141249985*fUpwind[27]*alphaDrSurf[34]+0.223606797749979*fUpwind[26]*alphaDrSurf[34]+0.223606797749979*fUpwind[25]*alphaDrSurf[34]+0.2500000000000001*fUpwind[4]*alphaDrSurf[34]+0.2*fUpwind[18]*alphaDrSurf[33]+0.2*fUpwind[17]*alphaDrSurf[32]+0.2*alphaDrSurf[22]*fUpwind[31]+0.2*alphaDrSurf[21]*fUpwind[31]+0.223606797749979*alphaDrSurf[3]*fUpwind[31]+0.25*alphaDrSurf[5]*fUpwind[27]+0.25*fUpwind[8]*alphaDrSurf[24]+0.25*fUpwind[9]*alphaDrSurf[23]+0.223606797749979*alphaDrSurf[6]*fUpwind[18]+0.223606797749979*alphaDrSurf[7]*fUpwind[17]+0.2500000000000001*alphaDrSurf[13]*fUpwind[16]+0.223606797749979*fUpwind[10]*alphaDrSurf[15]; 
  Ghat[47] = 0.223606797749979*alphaDrSurf[13]*fUpwind[47]+0.223606797749979*alphaDrSurf[12]*fUpwind[47]+0.223606797749979*alphaDrSurf[11]*fUpwind[47]+0.25*alphaDrSurf[0]*fUpwind[47]+0.223606797749979*alphaDrSurf[23]*fUpwind[43]+0.223606797749979*alphaDrSurf[20]*fUpwind[43]+0.2500000000000001*alphaDrSurf[1]*fUpwind[43]+0.223606797749979*alphaDrSurf[24]*fUpwind[42]+0.223606797749979*alphaDrSurf[19]*fUpwind[42]+0.2500000000000001*alphaDrSurf[2]*fUpwind[42]+0.223606797749979*alphaDrSurf[22]*fUpwind[41]+0.223606797749979*alphaDrSurf[21]*fUpwind[41]+0.2500000000000001*alphaDrSurf[3]*fUpwind[41]+0.223606797749979*fUpwind[30]*alphaDrSurf[34]+0.223606797749979*fUpwind[29]*alphaDrSurf[33]+0.223606797749979*fUpwind[28]*alphaDrSurf[32]+0.25*alphaDrSurf[5]*fUpwind[30]+0.25*alphaDrSurf[6]*fUpwind[29]+0.25*alphaDrSurf[7]*fUpwind[28]+0.2500000000000001*fUpwind[14]*alphaDrSurf[15]; 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[9] += 1.224744871391589*Ghat[1]*rdv2; 
  out[10] += 1.224744871391589*Ghat[2]*rdv2; 
  out[11] += 1.224744871391589*Ghat[3]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += -0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += -0.7071067811865475*Ghat[10]*rdv2; 
  out[15] += 1.224744871391589*Ghat[4]*rdv2; 
  out[16] += -0.7071067811865475*Ghat[11]*rdv2; 
  out[17] += -0.7071067811865475*Ghat[12]*rdv2; 
  out[18] += -0.7071067811865475*Ghat[13]*rdv2; 
  out[19] += -1.58113883008419*Ghat[0]*rdv2; 
  out[20] += -0.7071067811865475*Ghat[14]*rdv2; 
  out[21] += -0.7071067811865475*Ghat[15]*rdv2; 
  out[22] += 1.224744871391589*Ghat[5]*rdv2; 
  out[23] += 1.224744871391589*Ghat[6]*rdv2; 
  out[24] += 1.224744871391589*Ghat[7]*rdv2; 
  out[25] += -0.7071067811865475*Ghat[16]*rdv2; 
  out[26] += -0.7071067811865475*Ghat[17]*rdv2; 
  out[27] += -0.7071067811865475*Ghat[18]*rdv2; 
  out[28] += 1.224744871391589*Ghat[8]*rdv2; 
  out[29] += 1.224744871391589*Ghat[9]*rdv2; 
  out[30] += 1.224744871391589*Ghat[10]*rdv2; 
  out[31] += -0.7071067811865475*Ghat[19]*rdv2; 
  out[32] += -0.7071067811865475*Ghat[20]*rdv2; 
  out[33] += -0.7071067811865475*Ghat[21]*rdv2; 
  out[34] += -0.7071067811865475*Ghat[22]*rdv2; 
  out[35] += -0.7071067811865475*Ghat[23]*rdv2; 
  out[36] += -0.7071067811865475*Ghat[24]*rdv2; 
  out[37] += 1.224744871391589*Ghat[11]*rdv2; 
  out[38] += 1.224744871391589*Ghat[12]*rdv2; 
  out[39] += 1.224744871391589*Ghat[13]*rdv2; 
  out[40] += -1.58113883008419*Ghat[1]*rdv2; 
  out[41] += -1.58113883008419*Ghat[2]*rdv2; 
  out[42] += -1.58113883008419*Ghat[3]*rdv2; 
  out[43] += -0.7071067811865475*Ghat[25]*rdv2; 
  out[44] += -0.7071067811865475*Ghat[26]*rdv2; 
  out[45] += -0.7071067811865475*Ghat[27]*rdv2; 
  out[46] += -1.58113883008419*Ghat[4]*rdv2; 
  out[47] += -0.7071067811865475*Ghat[28]*rdv2; 
  out[48] += -0.7071067811865475*Ghat[29]*rdv2; 
  out[49] += -0.7071067811865475*Ghat[30]*rdv2; 
  out[50] += 1.224744871391589*Ghat[14]*rdv2; 
  out[51] += 1.224744871391589*Ghat[15]*rdv2; 
  out[52] += -0.7071067811865475*Ghat[31]*rdv2; 
  out[53] += 1.224744871391589*Ghat[16]*rdv2; 
  out[54] += 1.224744871391589*Ghat[17]*rdv2; 
  out[55] += 1.224744871391589*Ghat[18]*rdv2; 
  out[56] += -0.7071067811865475*Ghat[32]*rdv2; 
  out[57] += -0.7071067811865475*Ghat[33]*rdv2; 
  out[58] += -0.7071067811865475*Ghat[34]*rdv2; 
  out[59] += 1.224744871391589*Ghat[19]*rdv2; 
  out[60] += 1.224744871391589*Ghat[20]*rdv2; 
  out[61] += 1.224744871391589*Ghat[21]*rdv2; 
  out[62] += 1.224744871391589*Ghat[22]*rdv2; 
  out[63] += 1.224744871391589*Ghat[23]*rdv2; 
  out[64] += 1.224744871391589*Ghat[24]*rdv2; 
  out[65] += -1.58113883008419*Ghat[5]*rdv2; 
  out[66] += -1.58113883008419*Ghat[6]*rdv2; 
  out[67] += -1.58113883008419*Ghat[7]*rdv2; 
  out[68] += -0.7071067811865475*Ghat[35]*rdv2; 
  out[69] += -0.7071067811865475*Ghat[36]*rdv2; 
  out[70] += -0.7071067811865475*Ghat[37]*rdv2; 
  out[71] += -0.7071067811865475*Ghat[38]*rdv2; 
  out[72] += -0.7071067811865475*Ghat[39]*rdv2; 
  out[73] += -0.7071067811865475*Ghat[40]*rdv2; 
  out[74] += 1.224744871391589*Ghat[25]*rdv2; 
  out[75] += 1.224744871391589*Ghat[26]*rdv2; 
  out[76] += 1.224744871391589*Ghat[27]*rdv2; 
  out[77] += -1.58113883008419*Ghat[8]*rdv2; 
  out[78] += -1.58113883008419*Ghat[9]*rdv2; 
  out[79] += -1.58113883008419*Ghat[10]*rdv2; 
  out[80] += -0.7071067811865475*Ghat[41]*rdv2; 
  out[81] += -0.7071067811865475*Ghat[42]*rdv2; 
  out[82] += -0.7071067811865475*Ghat[43]*rdv2; 
  out[83] += 1.224744871391589*Ghat[28]*rdv2; 
  out[84] += 1.224744871391589*Ghat[29]*rdv2; 
  out[85] += 1.224744871391589*Ghat[30]*rdv2; 
  out[86] += 1.224744871391589*Ghat[31]*rdv2; 
  out[87] += 1.224744871391589*Ghat[32]*rdv2; 
  out[88] += 1.224744871391589*Ghat[33]*rdv2; 
  out[89] += 1.224744871391589*Ghat[34]*rdv2; 
  out[90] += -1.58113883008419*Ghat[15]*rdv2; 
  out[91] += -0.7071067811865475*Ghat[44]*rdv2; 
  out[92] += -0.7071067811865475*Ghat[45]*rdv2; 
  out[93] += -0.7071067811865475*Ghat[46]*rdv2; 
  out[94] += 1.224744871391589*Ghat[35]*rdv2; 
  out[95] += 1.224744871391589*Ghat[36]*rdv2; 
  out[96] += 1.224744871391589*Ghat[37]*rdv2; 
  out[97] += 1.224744871391589*Ghat[38]*rdv2; 
  out[98] += 1.224744871391589*Ghat[39]*rdv2; 
  out[99] += 1.224744871391589*Ghat[40]*rdv2; 
  out[100] += -1.58113883008419*Ghat[16]*rdv2; 
  out[101] += -1.58113883008419*Ghat[17]*rdv2; 
  out[102] += -1.58113883008419*Ghat[18]*rdv2; 
  out[103] += -0.7071067811865475*Ghat[47]*rdv2; 
  out[104] += 1.224744871391589*Ghat[41]*rdv2; 
  out[105] += 1.224744871391589*Ghat[42]*rdv2; 
  out[106] += 1.224744871391589*Ghat[43]*rdv2; 
  out[107] += 1.224744871391589*Ghat[44]*rdv2; 
  out[108] += 1.224744871391589*Ghat[45]*rdv2; 
  out[109] += 1.224744871391589*Ghat[46]*rdv2; 
  out[110] += -1.58113883008419*Ghat[31]*rdv2; 
  out[111] += 1.224744871391589*Ghat[47]*rdv2; 

  } 
} 
