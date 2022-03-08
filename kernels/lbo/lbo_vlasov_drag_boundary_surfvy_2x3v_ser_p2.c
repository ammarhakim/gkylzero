#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx4_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind.h> 
GKYL_CU_DH void lbo_vlasov_drag_boundary_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         Cell-center coordinates. 
  // dxv[5]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[24]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[8]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *sumNuUy = &nuUSum[8]; 

  double alphaDrSurf[48] = {0.0}; 
  double fUpwindQuad[81] = {0.0};
  double fUpwind[48] = {0.0};;
  double drag_incr[48] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[1]; 
  alphaDrSurf[2] = nuSum[2]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[2]; 
  alphaDrSurf[5] = 2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]; 
  alphaDrSurf[11] = (2.0*w[3]+dxv[3])*nuSum[4]-2.0*sumNuUy[4]; 
  alphaDrSurf[12] = (2.0*w[3]+dxv[3])*nuSum[5]-2.0*sumNuUy[5]; 
  alphaDrSurf[19] = (2.0*w[3]+dxv[3])*nuSum[6]-2.0*sumNuUy[6]; 
  alphaDrSurf[20] = (2.0*w[3]+dxv[3])*nuSum[7]-2.0*sumNuUy[7]; 

  if ((-0.2999999999999998*alphaDrSurf[20])-0.2999999999999999*alphaDrSurf[19]+0.223606797749979*(alphaDrSurf[12]+alphaDrSurf[11])+0.45*alphaDrSurf[5]-0.3354101966249685*(alphaDrSurf[2]+alphaDrSurf[1])+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx4_quad_0_r(fSkin); 
    fUpwindQuad[9] = ser_5x_p2_surfx4_quad_9_r(fSkin); 
    fUpwindQuad[18] = ser_5x_p2_surfx4_quad_18_r(fSkin); 
    fUpwindQuad[27] = ser_5x_p2_surfx4_quad_27_r(fSkin); 
    fUpwindQuad[36] = ser_5x_p2_surfx4_quad_36_r(fSkin); 
    fUpwindQuad[45] = ser_5x_p2_surfx4_quad_45_r(fSkin); 
    fUpwindQuad[54] = ser_5x_p2_surfx4_quad_54_r(fSkin); 
    fUpwindQuad[63] = ser_5x_p2_surfx4_quad_63_r(fSkin); 
    fUpwindQuad[72] = ser_5x_p2_surfx4_quad_72_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx4_quad_0_l(fEdge); 
    fUpwindQuad[9] = ser_5x_p2_surfx4_quad_9_l(fEdge); 
    fUpwindQuad[18] = ser_5x_p2_surfx4_quad_18_l(fEdge); 
    fUpwindQuad[27] = ser_5x_p2_surfx4_quad_27_l(fEdge); 
    fUpwindQuad[36] = ser_5x_p2_surfx4_quad_36_l(fEdge); 
    fUpwindQuad[45] = ser_5x_p2_surfx4_quad_45_l(fEdge); 
    fUpwindQuad[54] = ser_5x_p2_surfx4_quad_54_l(fEdge); 
    fUpwindQuad[63] = ser_5x_p2_surfx4_quad_63_l(fEdge); 
    fUpwindQuad[72] = ser_5x_p2_surfx4_quad_72_l(fEdge); 
  } 
  if (0.375*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[12]-0.2795084971874737*alphaDrSurf[11]-0.3354101966249685*alphaDrSurf[2]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_5x_p2_surfx4_quad_1_r(fSkin); 
    fUpwindQuad[10] = ser_5x_p2_surfx4_quad_10_r(fSkin); 
    fUpwindQuad[19] = ser_5x_p2_surfx4_quad_19_r(fSkin); 
    fUpwindQuad[28] = ser_5x_p2_surfx4_quad_28_r(fSkin); 
    fUpwindQuad[37] = ser_5x_p2_surfx4_quad_37_r(fSkin); 
    fUpwindQuad[46] = ser_5x_p2_surfx4_quad_46_r(fSkin); 
    fUpwindQuad[55] = ser_5x_p2_surfx4_quad_55_r(fSkin); 
    fUpwindQuad[64] = ser_5x_p2_surfx4_quad_64_r(fSkin); 
    fUpwindQuad[73] = ser_5x_p2_surfx4_quad_73_r(fSkin); 
  } else { 
    fUpwindQuad[1] = ser_5x_p2_surfx4_quad_1_l(fEdge); 
    fUpwindQuad[10] = ser_5x_p2_surfx4_quad_10_l(fEdge); 
    fUpwindQuad[19] = ser_5x_p2_surfx4_quad_19_l(fEdge); 
    fUpwindQuad[28] = ser_5x_p2_surfx4_quad_28_l(fEdge); 
    fUpwindQuad[37] = ser_5x_p2_surfx4_quad_37_l(fEdge); 
    fUpwindQuad[46] = ser_5x_p2_surfx4_quad_46_l(fEdge); 
    fUpwindQuad[55] = ser_5x_p2_surfx4_quad_55_l(fEdge); 
    fUpwindQuad[64] = ser_5x_p2_surfx4_quad_64_l(fEdge); 
    fUpwindQuad[73] = ser_5x_p2_surfx4_quad_73_l(fEdge); 
  } 
  if (0.2999999999999998*alphaDrSurf[20]-0.2999999999999999*alphaDrSurf[19]+0.223606797749979*(alphaDrSurf[12]+alphaDrSurf[11])-0.45*alphaDrSurf[5]-0.3354101966249685*alphaDrSurf[2]+0.3354101966249685*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_5x_p2_surfx4_quad_2_r(fSkin); 
    fUpwindQuad[11] = ser_5x_p2_surfx4_quad_11_r(fSkin); 
    fUpwindQuad[20] = ser_5x_p2_surfx4_quad_20_r(fSkin); 
    fUpwindQuad[29] = ser_5x_p2_surfx4_quad_29_r(fSkin); 
    fUpwindQuad[38] = ser_5x_p2_surfx4_quad_38_r(fSkin); 
    fUpwindQuad[47] = ser_5x_p2_surfx4_quad_47_r(fSkin); 
    fUpwindQuad[56] = ser_5x_p2_surfx4_quad_56_r(fSkin); 
    fUpwindQuad[65] = ser_5x_p2_surfx4_quad_65_r(fSkin); 
    fUpwindQuad[74] = ser_5x_p2_surfx4_quad_74_r(fSkin); 
  } else { 
    fUpwindQuad[2] = ser_5x_p2_surfx4_quad_2_l(fEdge); 
    fUpwindQuad[11] = ser_5x_p2_surfx4_quad_11_l(fEdge); 
    fUpwindQuad[20] = ser_5x_p2_surfx4_quad_20_l(fEdge); 
    fUpwindQuad[29] = ser_5x_p2_surfx4_quad_29_l(fEdge); 
    fUpwindQuad[38] = ser_5x_p2_surfx4_quad_38_l(fEdge); 
    fUpwindQuad[47] = ser_5x_p2_surfx4_quad_47_l(fEdge); 
    fUpwindQuad[56] = ser_5x_p2_surfx4_quad_56_l(fEdge); 
    fUpwindQuad[65] = ser_5x_p2_surfx4_quad_65_l(fEdge); 
    fUpwindQuad[74] = ser_5x_p2_surfx4_quad_74_l(fEdge); 
  } 
  if (0.375*alphaDrSurf[20]-0.2795084971874737*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[11]-0.3354101966249685*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_5x_p2_surfx4_quad_3_r(fSkin); 
    fUpwindQuad[12] = ser_5x_p2_surfx4_quad_12_r(fSkin); 
    fUpwindQuad[21] = ser_5x_p2_surfx4_quad_21_r(fSkin); 
    fUpwindQuad[30] = ser_5x_p2_surfx4_quad_30_r(fSkin); 
    fUpwindQuad[39] = ser_5x_p2_surfx4_quad_39_r(fSkin); 
    fUpwindQuad[48] = ser_5x_p2_surfx4_quad_48_r(fSkin); 
    fUpwindQuad[57] = ser_5x_p2_surfx4_quad_57_r(fSkin); 
    fUpwindQuad[66] = ser_5x_p2_surfx4_quad_66_r(fSkin); 
    fUpwindQuad[75] = ser_5x_p2_surfx4_quad_75_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_5x_p2_surfx4_quad_3_l(fEdge); 
    fUpwindQuad[12] = ser_5x_p2_surfx4_quad_12_l(fEdge); 
    fUpwindQuad[21] = ser_5x_p2_surfx4_quad_21_l(fEdge); 
    fUpwindQuad[30] = ser_5x_p2_surfx4_quad_30_l(fEdge); 
    fUpwindQuad[39] = ser_5x_p2_surfx4_quad_39_l(fEdge); 
    fUpwindQuad[48] = ser_5x_p2_surfx4_quad_48_l(fEdge); 
    fUpwindQuad[57] = ser_5x_p2_surfx4_quad_57_l(fEdge); 
    fUpwindQuad[66] = ser_5x_p2_surfx4_quad_66_l(fEdge); 
    fUpwindQuad[75] = ser_5x_p2_surfx4_quad_75_l(fEdge); 
  } 
  if (0.25*alphaDrSurf[0]-0.2795084971874737*(alphaDrSurf[12]+alphaDrSurf[11]) < 0) { 
    fUpwindQuad[4] = ser_5x_p2_surfx4_quad_4_r(fSkin); 
    fUpwindQuad[13] = ser_5x_p2_surfx4_quad_13_r(fSkin); 
    fUpwindQuad[22] = ser_5x_p2_surfx4_quad_22_r(fSkin); 
    fUpwindQuad[31] = ser_5x_p2_surfx4_quad_31_r(fSkin); 
    fUpwindQuad[40] = ser_5x_p2_surfx4_quad_40_r(fSkin); 
    fUpwindQuad[49] = ser_5x_p2_surfx4_quad_49_r(fSkin); 
    fUpwindQuad[58] = ser_5x_p2_surfx4_quad_58_r(fSkin); 
    fUpwindQuad[67] = ser_5x_p2_surfx4_quad_67_r(fSkin); 
    fUpwindQuad[76] = ser_5x_p2_surfx4_quad_76_r(fSkin); 
  } else { 
    fUpwindQuad[4] = ser_5x_p2_surfx4_quad_4_l(fEdge); 
    fUpwindQuad[13] = ser_5x_p2_surfx4_quad_13_l(fEdge); 
    fUpwindQuad[22] = ser_5x_p2_surfx4_quad_22_l(fEdge); 
    fUpwindQuad[31] = ser_5x_p2_surfx4_quad_31_l(fEdge); 
    fUpwindQuad[40] = ser_5x_p2_surfx4_quad_40_l(fEdge); 
    fUpwindQuad[49] = ser_5x_p2_surfx4_quad_49_l(fEdge); 
    fUpwindQuad[58] = ser_5x_p2_surfx4_quad_58_l(fEdge); 
    fUpwindQuad[67] = ser_5x_p2_surfx4_quad_67_l(fEdge); 
    fUpwindQuad[76] = ser_5x_p2_surfx4_quad_76_l(fEdge); 
  } 
  if ((-0.375*alphaDrSurf[20])-0.2795084971874737*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[11]+0.3354101966249685*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_5x_p2_surfx4_quad_5_r(fSkin); 
    fUpwindQuad[14] = ser_5x_p2_surfx4_quad_14_r(fSkin); 
    fUpwindQuad[23] = ser_5x_p2_surfx4_quad_23_r(fSkin); 
    fUpwindQuad[32] = ser_5x_p2_surfx4_quad_32_r(fSkin); 
    fUpwindQuad[41] = ser_5x_p2_surfx4_quad_41_r(fSkin); 
    fUpwindQuad[50] = ser_5x_p2_surfx4_quad_50_r(fSkin); 
    fUpwindQuad[59] = ser_5x_p2_surfx4_quad_59_r(fSkin); 
    fUpwindQuad[68] = ser_5x_p2_surfx4_quad_68_r(fSkin); 
    fUpwindQuad[77] = ser_5x_p2_surfx4_quad_77_r(fSkin); 
  } else { 
    fUpwindQuad[5] = ser_5x_p2_surfx4_quad_5_l(fEdge); 
    fUpwindQuad[14] = ser_5x_p2_surfx4_quad_14_l(fEdge); 
    fUpwindQuad[23] = ser_5x_p2_surfx4_quad_23_l(fEdge); 
    fUpwindQuad[32] = ser_5x_p2_surfx4_quad_32_l(fEdge); 
    fUpwindQuad[41] = ser_5x_p2_surfx4_quad_41_l(fEdge); 
    fUpwindQuad[50] = ser_5x_p2_surfx4_quad_50_l(fEdge); 
    fUpwindQuad[59] = ser_5x_p2_surfx4_quad_59_l(fEdge); 
    fUpwindQuad[68] = ser_5x_p2_surfx4_quad_68_l(fEdge); 
    fUpwindQuad[77] = ser_5x_p2_surfx4_quad_77_l(fEdge); 
  } 
  if ((-0.2999999999999998*alphaDrSurf[20])+0.2999999999999999*alphaDrSurf[19]+0.223606797749979*(alphaDrSurf[12]+alphaDrSurf[11])-0.45*alphaDrSurf[5]+0.3354101966249685*alphaDrSurf[2]-0.3354101966249685*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_5x_p2_surfx4_quad_6_r(fSkin); 
    fUpwindQuad[15] = ser_5x_p2_surfx4_quad_15_r(fSkin); 
    fUpwindQuad[24] = ser_5x_p2_surfx4_quad_24_r(fSkin); 
    fUpwindQuad[33] = ser_5x_p2_surfx4_quad_33_r(fSkin); 
    fUpwindQuad[42] = ser_5x_p2_surfx4_quad_42_r(fSkin); 
    fUpwindQuad[51] = ser_5x_p2_surfx4_quad_51_r(fSkin); 
    fUpwindQuad[60] = ser_5x_p2_surfx4_quad_60_r(fSkin); 
    fUpwindQuad[69] = ser_5x_p2_surfx4_quad_69_r(fSkin); 
    fUpwindQuad[78] = ser_5x_p2_surfx4_quad_78_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_5x_p2_surfx4_quad_6_l(fEdge); 
    fUpwindQuad[15] = ser_5x_p2_surfx4_quad_15_l(fEdge); 
    fUpwindQuad[24] = ser_5x_p2_surfx4_quad_24_l(fEdge); 
    fUpwindQuad[33] = ser_5x_p2_surfx4_quad_33_l(fEdge); 
    fUpwindQuad[42] = ser_5x_p2_surfx4_quad_42_l(fEdge); 
    fUpwindQuad[51] = ser_5x_p2_surfx4_quad_51_l(fEdge); 
    fUpwindQuad[60] = ser_5x_p2_surfx4_quad_60_l(fEdge); 
    fUpwindQuad[69] = ser_5x_p2_surfx4_quad_69_l(fEdge); 
    fUpwindQuad[78] = ser_5x_p2_surfx4_quad_78_l(fEdge); 
  } 
  if ((-0.375*alphaDrSurf[19])+0.223606797749979*alphaDrSurf[12]-0.2795084971874737*alphaDrSurf[11]+0.3354101966249685*alphaDrSurf[2]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_5x_p2_surfx4_quad_7_r(fSkin); 
    fUpwindQuad[16] = ser_5x_p2_surfx4_quad_16_r(fSkin); 
    fUpwindQuad[25] = ser_5x_p2_surfx4_quad_25_r(fSkin); 
    fUpwindQuad[34] = ser_5x_p2_surfx4_quad_34_r(fSkin); 
    fUpwindQuad[43] = ser_5x_p2_surfx4_quad_43_r(fSkin); 
    fUpwindQuad[52] = ser_5x_p2_surfx4_quad_52_r(fSkin); 
    fUpwindQuad[61] = ser_5x_p2_surfx4_quad_61_r(fSkin); 
    fUpwindQuad[70] = ser_5x_p2_surfx4_quad_70_r(fSkin); 
    fUpwindQuad[79] = ser_5x_p2_surfx4_quad_79_r(fSkin); 
  } else { 
    fUpwindQuad[7] = ser_5x_p2_surfx4_quad_7_l(fEdge); 
    fUpwindQuad[16] = ser_5x_p2_surfx4_quad_16_l(fEdge); 
    fUpwindQuad[25] = ser_5x_p2_surfx4_quad_25_l(fEdge); 
    fUpwindQuad[34] = ser_5x_p2_surfx4_quad_34_l(fEdge); 
    fUpwindQuad[43] = ser_5x_p2_surfx4_quad_43_l(fEdge); 
    fUpwindQuad[52] = ser_5x_p2_surfx4_quad_52_l(fEdge); 
    fUpwindQuad[61] = ser_5x_p2_surfx4_quad_61_l(fEdge); 
    fUpwindQuad[70] = ser_5x_p2_surfx4_quad_70_l(fEdge); 
    fUpwindQuad[79] = ser_5x_p2_surfx4_quad_79_l(fEdge); 
  } 
  if (0.2999999999999998*alphaDrSurf[20]+0.2999999999999999*alphaDrSurf[19]+0.223606797749979*(alphaDrSurf[12]+alphaDrSurf[11])+0.45*alphaDrSurf[5]+0.3354101966249685*(alphaDrSurf[2]+alphaDrSurf[1])+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_5x_p2_surfx4_quad_8_r(fSkin); 
    fUpwindQuad[17] = ser_5x_p2_surfx4_quad_17_r(fSkin); 
    fUpwindQuad[26] = ser_5x_p2_surfx4_quad_26_r(fSkin); 
    fUpwindQuad[35] = ser_5x_p2_surfx4_quad_35_r(fSkin); 
    fUpwindQuad[44] = ser_5x_p2_surfx4_quad_44_r(fSkin); 
    fUpwindQuad[53] = ser_5x_p2_surfx4_quad_53_r(fSkin); 
    fUpwindQuad[62] = ser_5x_p2_surfx4_quad_62_r(fSkin); 
    fUpwindQuad[71] = ser_5x_p2_surfx4_quad_71_r(fSkin); 
    fUpwindQuad[80] = ser_5x_p2_surfx4_quad_80_r(fSkin); 
  } else { 
    fUpwindQuad[8] = ser_5x_p2_surfx4_quad_8_l(fEdge); 
    fUpwindQuad[17] = ser_5x_p2_surfx4_quad_17_l(fEdge); 
    fUpwindQuad[26] = ser_5x_p2_surfx4_quad_26_l(fEdge); 
    fUpwindQuad[35] = ser_5x_p2_surfx4_quad_35_l(fEdge); 
    fUpwindQuad[44] = ser_5x_p2_surfx4_quad_44_l(fEdge); 
    fUpwindQuad[53] = ser_5x_p2_surfx4_quad_53_l(fEdge); 
    fUpwindQuad[62] = ser_5x_p2_surfx4_quad_62_l(fEdge); 
    fUpwindQuad[71] = ser_5x_p2_surfx4_quad_71_l(fEdge); 
    fUpwindQuad[80] = ser_5x_p2_surfx4_quad_80_l(fEdge); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p2_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.25*alphaDrSurf[20]*fUpwind[20]+0.25*alphaDrSurf[19]*fUpwind[19]+0.25*alphaDrSurf[12]*fUpwind[12]+0.25*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.2500000000000001*alphaDrSurf[12]*fUpwind[20]+0.2500000000000001*fUpwind[12]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[19]+0.223606797749979*fUpwind[5]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[11]+0.223606797749979*fUpwind[1]*alphaDrSurf[11]+0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.223606797749979*alphaDrSurf[5]*fUpwind[20]+0.223606797749979*fUpwind[5]*alphaDrSurf[20]+0.2500000000000001*alphaDrSurf[11]*fUpwind[19]+0.2500000000000001*fUpwind[11]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[2]*fUpwind[12]+0.223606797749979*fUpwind[2]*alphaDrSurf[12]+0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.2500000000000001*alphaDrSurf[20]*fUpwind[33]+0.2500000000000001*alphaDrSurf[19]*fUpwind[32]+0.2500000000000001*alphaDrSurf[12]*fUpwind[22]+0.2500000000000001*alphaDrSurf[11]*fUpwind[21]+0.25*alphaDrSurf[5]*fUpwind[15]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.2500000000000001*alphaDrSurf[20]*fUpwind[36]+0.2500000000000001*alphaDrSurf[19]*fUpwind[35]+0.2500000000000001*alphaDrSurf[12]*fUpwind[26]+0.2500000000000001*alphaDrSurf[11]*fUpwind[25]+0.25*alphaDrSurf[5]*fUpwind[16]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  drag_incr[5] = 0.2*alphaDrSurf[19]*fUpwind[20]+0.223606797749979*alphaDrSurf[2]*fUpwind[20]+0.2*fUpwind[19]*alphaDrSurf[20]+0.223606797749979*fUpwind[2]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[1]*fUpwind[19]+0.223606797749979*fUpwind[1]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[5]*fUpwind[12]+0.223606797749979*fUpwind[5]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[11]+0.223606797749979*fUpwind[5]*alphaDrSurf[11]+0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[6] = 0.25*alphaDrSurf[12]*fUpwind[33]+0.223606797749979*alphaDrSurf[5]*fUpwind[32]+0.25*alphaDrSurf[20]*fUpwind[22]+0.223606797749979*alphaDrSurf[1]*fUpwind[21]+0.223606797749979*fUpwind[15]*alphaDrSurf[19]+0.25*alphaDrSurf[2]*fUpwind[15]+0.223606797749979*fUpwind[6]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[7] = 0.223606797749979*alphaDrSurf[5]*fUpwind[33]+0.25*alphaDrSurf[11]*fUpwind[32]+0.223606797749979*alphaDrSurf[2]*fUpwind[22]+0.25*alphaDrSurf[19]*fUpwind[21]+0.223606797749979*fUpwind[15]*alphaDrSurf[20]+0.25*alphaDrSurf[1]*fUpwind[15]+0.223606797749979*fUpwind[7]*alphaDrSurf[12]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*alphaDrSurf[2]*fUpwind[3]; 
  drag_incr[8] = 0.25*alphaDrSurf[12]*fUpwind[36]+0.223606797749979*alphaDrSurf[5]*fUpwind[35]+0.25*alphaDrSurf[20]*fUpwind[26]+0.223606797749979*alphaDrSurf[1]*fUpwind[25]+0.223606797749979*fUpwind[16]*alphaDrSurf[19]+0.25*alphaDrSurf[2]*fUpwind[16]+0.223606797749979*fUpwind[8]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[9] = 0.223606797749979*alphaDrSurf[5]*fUpwind[36]+0.25*alphaDrSurf[11]*fUpwind[35]+0.223606797749979*alphaDrSurf[2]*fUpwind[26]+0.25*alphaDrSurf[19]*fUpwind[25]+0.223606797749979*fUpwind[16]*alphaDrSurf[20]+0.25*alphaDrSurf[1]*fUpwind[16]+0.223606797749979*fUpwind[9]*alphaDrSurf[12]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  drag_incr[10] = 0.25*alphaDrSurf[20]*fUpwind[45]+0.25*alphaDrSurf[19]*fUpwind[44]+0.25*alphaDrSurf[12]*fUpwind[38]+0.25*alphaDrSurf[11]*fUpwind[37]+0.25*alphaDrSurf[5]*fUpwind[31]+0.25*alphaDrSurf[2]*fUpwind[18]+0.25*alphaDrSurf[1]*fUpwind[17]+0.25*alphaDrSurf[0]*fUpwind[10]; 
  drag_incr[11] = 0.223606797749979*alphaDrSurf[20]*fUpwind[20]+0.159719141249985*alphaDrSurf[19]*fUpwind[19]+0.2500000000000001*alphaDrSurf[2]*fUpwind[19]+0.2500000000000001*fUpwind[2]*alphaDrSurf[19]+0.159719141249985*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[0]*fUpwind[11]+0.25*fUpwind[0]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[5]+0.223606797749979*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[12] = 0.159719141249985*alphaDrSurf[20]*fUpwind[20]+0.2500000000000001*alphaDrSurf[1]*fUpwind[20]+0.2500000000000001*fUpwind[1]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[19]*fUpwind[19]+0.159719141249985*alphaDrSurf[12]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[12]+0.25*fUpwind[0]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[5]+0.223606797749979*alphaDrSurf[2]*fUpwind[2]; 
  drag_incr[13] = 0.25*alphaDrSurf[5]*fUpwind[34]+0.2500000000000001*alphaDrSurf[2]*fUpwind[24]+0.2500000000000001*alphaDrSurf[1]*fUpwind[23]+0.25*alphaDrSurf[0]*fUpwind[13]; 
  drag_incr[14] = 0.25*alphaDrSurf[5]*fUpwind[41]+0.2500000000000001*alphaDrSurf[2]*fUpwind[29]+0.2500000000000001*alphaDrSurf[1]*fUpwind[28]+0.25*alphaDrSurf[0]*fUpwind[14]; 
  drag_incr[15] = 0.2*alphaDrSurf[19]*fUpwind[33]+0.223606797749979*alphaDrSurf[2]*fUpwind[33]+0.2*alphaDrSurf[20]*fUpwind[32]+0.223606797749979*alphaDrSurf[1]*fUpwind[32]+0.223606797749979*alphaDrSurf[5]*fUpwind[22]+0.223606797749979*alphaDrSurf[5]*fUpwind[21]+0.223606797749979*fUpwind[7]*alphaDrSurf[20]+0.223606797749979*fUpwind[6]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[12]*fUpwind[15]+0.223606797749979*alphaDrSurf[11]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[15]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  drag_incr[16] = 0.2*alphaDrSurf[19]*fUpwind[36]+0.223606797749979*alphaDrSurf[2]*fUpwind[36]+0.2*alphaDrSurf[20]*fUpwind[35]+0.223606797749979*alphaDrSurf[1]*fUpwind[35]+0.223606797749979*alphaDrSurf[5]*fUpwind[26]+0.223606797749979*alphaDrSurf[5]*fUpwind[25]+0.223606797749979*fUpwind[9]*alphaDrSurf[20]+0.223606797749979*fUpwind[8]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[12]*fUpwind[16]+0.223606797749979*alphaDrSurf[11]*fUpwind[16]+0.25*alphaDrSurf[0]*fUpwind[16]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  drag_incr[17] = 0.2500000000000001*alphaDrSurf[12]*fUpwind[45]+0.223606797749979*alphaDrSurf[5]*fUpwind[44]+0.2500000000000001*alphaDrSurf[20]*fUpwind[38]+0.223606797749979*alphaDrSurf[1]*fUpwind[37]+0.223606797749979*alphaDrSurf[19]*fUpwind[31]+0.25*alphaDrSurf[2]*fUpwind[31]+0.25*alphaDrSurf[5]*fUpwind[18]+0.223606797749979*alphaDrSurf[11]*fUpwind[17]+0.25*alphaDrSurf[0]*fUpwind[17]+0.25*alphaDrSurf[1]*fUpwind[10]; 
  drag_incr[18] = 0.223606797749979*alphaDrSurf[5]*fUpwind[45]+0.2500000000000001*alphaDrSurf[11]*fUpwind[44]+0.223606797749979*alphaDrSurf[2]*fUpwind[38]+0.2500000000000001*alphaDrSurf[19]*fUpwind[37]+0.223606797749979*alphaDrSurf[20]*fUpwind[31]+0.25*alphaDrSurf[1]*fUpwind[31]+0.223606797749979*alphaDrSurf[12]*fUpwind[18]+0.25*alphaDrSurf[0]*fUpwind[18]+0.25*alphaDrSurf[5]*fUpwind[17]+0.25*alphaDrSurf[2]*fUpwind[10]; 
  drag_incr[19] = 0.2*alphaDrSurf[5]*fUpwind[20]+0.2*fUpwind[5]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[12]*fUpwind[19]+0.159719141249985*alphaDrSurf[11]*fUpwind[19]+0.25*alphaDrSurf[0]*fUpwind[19]+0.223606797749979*fUpwind[12]*alphaDrSurf[19]+0.159719141249985*fUpwind[11]*alphaDrSurf[19]+0.25*fUpwind[0]*alphaDrSurf[19]+0.2500000000000001*alphaDrSurf[2]*fUpwind[11]+0.2500000000000001*fUpwind[2]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[5]+0.223606797749979*fUpwind[1]*alphaDrSurf[5]; 
  drag_incr[20] = 0.159719141249985*alphaDrSurf[12]*fUpwind[20]+0.223606797749979*alphaDrSurf[11]*fUpwind[20]+0.25*alphaDrSurf[0]*fUpwind[20]+0.159719141249985*fUpwind[12]*alphaDrSurf[20]+0.223606797749979*fUpwind[11]*alphaDrSurf[20]+0.25*fUpwind[0]*alphaDrSurf[20]+0.2*alphaDrSurf[5]*fUpwind[19]+0.2*fUpwind[5]*alphaDrSurf[19]+0.2500000000000001*alphaDrSurf[1]*fUpwind[12]+0.2500000000000001*fUpwind[1]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[5]+0.223606797749979*fUpwind[2]*alphaDrSurf[5]; 
  drag_incr[21] = 0.223606797749979*alphaDrSurf[20]*fUpwind[33]+0.159719141249985*alphaDrSurf[19]*fUpwind[32]+0.2500000000000001*alphaDrSurf[2]*fUpwind[32]+0.159719141249985*alphaDrSurf[11]*fUpwind[21]+0.25*alphaDrSurf[0]*fUpwind[21]+0.25*fUpwind[7]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[5]*fUpwind[15]+0.2500000000000001*fUpwind[3]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[6]; 
  drag_incr[22] = 0.159719141249985*alphaDrSurf[20]*fUpwind[33]+0.2500000000000001*alphaDrSurf[1]*fUpwind[33]+0.223606797749979*alphaDrSurf[19]*fUpwind[32]+0.159719141249985*alphaDrSurf[12]*fUpwind[22]+0.25*alphaDrSurf[0]*fUpwind[22]+0.25*fUpwind[6]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[15]+0.2500000000000001*fUpwind[3]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[7]; 
  drag_incr[23] = 0.223606797749979*alphaDrSurf[19]*fUpwind[34]+0.2500000000000001*alphaDrSurf[2]*fUpwind[34]+0.25*alphaDrSurf[5]*fUpwind[24]+0.223606797749979*alphaDrSurf[11]*fUpwind[23]+0.25*alphaDrSurf[0]*fUpwind[23]+0.2500000000000001*alphaDrSurf[1]*fUpwind[13]; 
  drag_incr[24] = 0.223606797749979*alphaDrSurf[20]*fUpwind[34]+0.2500000000000001*alphaDrSurf[1]*fUpwind[34]+0.223606797749979*alphaDrSurf[12]*fUpwind[24]+0.25*alphaDrSurf[0]*fUpwind[24]+0.25*alphaDrSurf[5]*fUpwind[23]+0.2500000000000001*alphaDrSurf[2]*fUpwind[13]; 
  drag_incr[25] = 0.223606797749979*alphaDrSurf[20]*fUpwind[36]+0.159719141249985*alphaDrSurf[19]*fUpwind[35]+0.2500000000000001*alphaDrSurf[2]*fUpwind[35]+0.159719141249985*alphaDrSurf[11]*fUpwind[25]+0.25*alphaDrSurf[0]*fUpwind[25]+0.25*fUpwind[9]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[8]; 
  drag_incr[26] = 0.159719141249985*alphaDrSurf[20]*fUpwind[36]+0.2500000000000001*alphaDrSurf[1]*fUpwind[36]+0.223606797749979*alphaDrSurf[19]*fUpwind[35]+0.159719141249985*alphaDrSurf[12]*fUpwind[26]+0.25*alphaDrSurf[0]*fUpwind[26]+0.25*fUpwind[8]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[9]; 
  drag_incr[27] = 0.25*alphaDrSurf[5]*fUpwind[46]+0.2500000000000001*alphaDrSurf[2]*fUpwind[40]+0.2500000000000001*alphaDrSurf[1]*fUpwind[39]+0.25*alphaDrSurf[0]*fUpwind[27]; 
  drag_incr[28] = 0.223606797749979*alphaDrSurf[19]*fUpwind[41]+0.2500000000000001*alphaDrSurf[2]*fUpwind[41]+0.25*alphaDrSurf[5]*fUpwind[29]+0.223606797749979*alphaDrSurf[11]*fUpwind[28]+0.25*alphaDrSurf[0]*fUpwind[28]+0.2500000000000001*alphaDrSurf[1]*fUpwind[14]; 
  drag_incr[29] = 0.223606797749979*alphaDrSurf[20]*fUpwind[41]+0.2500000000000001*alphaDrSurf[1]*fUpwind[41]+0.223606797749979*alphaDrSurf[12]*fUpwind[29]+0.25*alphaDrSurf[0]*fUpwind[29]+0.25*alphaDrSurf[5]*fUpwind[28]+0.2500000000000001*alphaDrSurf[2]*fUpwind[14]; 
  drag_incr[30] = 0.25*alphaDrSurf[5]*fUpwind[47]+0.2500000000000001*alphaDrSurf[2]*fUpwind[43]+0.2500000000000001*alphaDrSurf[1]*fUpwind[42]+0.25*alphaDrSurf[0]*fUpwind[30]; 
  drag_incr[31] = 0.2*alphaDrSurf[19]*fUpwind[45]+0.223606797749979*alphaDrSurf[2]*fUpwind[45]+0.2*alphaDrSurf[20]*fUpwind[44]+0.223606797749979*alphaDrSurf[1]*fUpwind[44]+0.223606797749979*alphaDrSurf[5]*fUpwind[38]+0.223606797749979*alphaDrSurf[5]*fUpwind[37]+0.223606797749979*alphaDrSurf[12]*fUpwind[31]+0.223606797749979*alphaDrSurf[11]*fUpwind[31]+0.25*alphaDrSurf[0]*fUpwind[31]+0.223606797749979*fUpwind[18]*alphaDrSurf[20]+0.223606797749979*fUpwind[17]*alphaDrSurf[19]+0.25*alphaDrSurf[1]*fUpwind[18]+0.25*alphaDrSurf[2]*fUpwind[17]+0.25*alphaDrSurf[5]*fUpwind[10]; 
  drag_incr[32] = 0.2*alphaDrSurf[5]*fUpwind[33]+0.223606797749979*alphaDrSurf[12]*fUpwind[32]+0.159719141249985*alphaDrSurf[11]*fUpwind[32]+0.25*alphaDrSurf[0]*fUpwind[32]+0.223606797749979*alphaDrSurf[19]*fUpwind[22]+0.159719141249985*alphaDrSurf[19]*fUpwind[21]+0.2500000000000001*alphaDrSurf[2]*fUpwind[21]+0.2*fUpwind[15]*alphaDrSurf[20]+0.2500000000000001*fUpwind[3]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[15]+0.25*fUpwind[7]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[6]; 
  drag_incr[33] = 0.159719141249985*alphaDrSurf[12]*fUpwind[33]+0.223606797749979*alphaDrSurf[11]*fUpwind[33]+0.25*alphaDrSurf[0]*fUpwind[33]+0.2*alphaDrSurf[5]*fUpwind[32]+0.159719141249985*alphaDrSurf[20]*fUpwind[22]+0.2500000000000001*alphaDrSurf[1]*fUpwind[22]+0.223606797749979*alphaDrSurf[20]*fUpwind[21]+0.2500000000000001*fUpwind[3]*alphaDrSurf[20]+0.2*fUpwind[15]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[2]*fUpwind[15]+0.25*fUpwind[6]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[7]; 
  drag_incr[34] = 0.223606797749979*alphaDrSurf[12]*fUpwind[34]+0.223606797749979*alphaDrSurf[11]*fUpwind[34]+0.25*alphaDrSurf[0]*fUpwind[34]+0.223606797749979*alphaDrSurf[20]*fUpwind[24]+0.2500000000000001*alphaDrSurf[1]*fUpwind[24]+0.223606797749979*alphaDrSurf[19]*fUpwind[23]+0.2500000000000001*alphaDrSurf[2]*fUpwind[23]+0.25*alphaDrSurf[5]*fUpwind[13]; 
  drag_incr[35] = 0.2*alphaDrSurf[5]*fUpwind[36]+0.223606797749979*alphaDrSurf[12]*fUpwind[35]+0.159719141249985*alphaDrSurf[11]*fUpwind[35]+0.25*alphaDrSurf[0]*fUpwind[35]+0.223606797749979*alphaDrSurf[19]*fUpwind[26]+0.159719141249985*alphaDrSurf[19]*fUpwind[25]+0.2500000000000001*alphaDrSurf[2]*fUpwind[25]+0.2*fUpwind[16]*alphaDrSurf[20]+0.2500000000000001*fUpwind[4]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[16]+0.25*fUpwind[9]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[8]; 
  drag_incr[36] = 0.159719141249985*alphaDrSurf[12]*fUpwind[36]+0.223606797749979*alphaDrSurf[11]*fUpwind[36]+0.25*alphaDrSurf[0]*fUpwind[36]+0.2*alphaDrSurf[5]*fUpwind[35]+0.159719141249985*alphaDrSurf[20]*fUpwind[26]+0.2500000000000001*alphaDrSurf[1]*fUpwind[26]+0.223606797749979*alphaDrSurf[20]*fUpwind[25]+0.2500000000000001*fUpwind[4]*alphaDrSurf[20]+0.2*fUpwind[16]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[2]*fUpwind[16]+0.25*fUpwind[8]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[9]; 
  drag_incr[37] = 0.223606797749979*alphaDrSurf[20]*fUpwind[45]+0.159719141249985*alphaDrSurf[19]*fUpwind[44]+0.2500000000000001*alphaDrSurf[2]*fUpwind[44]+0.159719141249985*alphaDrSurf[11]*fUpwind[37]+0.25*alphaDrSurf[0]*fUpwind[37]+0.223606797749979*alphaDrSurf[5]*fUpwind[31]+0.2500000000000001*fUpwind[18]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[17]+0.25*fUpwind[10]*alphaDrSurf[11]; 
  drag_incr[38] = 0.159719141249985*alphaDrSurf[20]*fUpwind[45]+0.2500000000000001*alphaDrSurf[1]*fUpwind[45]+0.223606797749979*alphaDrSurf[19]*fUpwind[44]+0.159719141249985*alphaDrSurf[12]*fUpwind[38]+0.25*alphaDrSurf[0]*fUpwind[38]+0.223606797749979*alphaDrSurf[5]*fUpwind[31]+0.2500000000000001*fUpwind[17]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[2]*fUpwind[18]+0.25*fUpwind[10]*alphaDrSurf[12]; 
  drag_incr[39] = 0.223606797749979*alphaDrSurf[19]*fUpwind[46]+0.2500000000000001*alphaDrSurf[2]*fUpwind[46]+0.25*alphaDrSurf[5]*fUpwind[40]+0.223606797749979*alphaDrSurf[11]*fUpwind[39]+0.25*alphaDrSurf[0]*fUpwind[39]+0.2500000000000001*alphaDrSurf[1]*fUpwind[27]; 
  drag_incr[40] = 0.223606797749979*alphaDrSurf[20]*fUpwind[46]+0.2500000000000001*alphaDrSurf[1]*fUpwind[46]+0.223606797749979*alphaDrSurf[12]*fUpwind[40]+0.25*alphaDrSurf[0]*fUpwind[40]+0.25*alphaDrSurf[5]*fUpwind[39]+0.2500000000000001*alphaDrSurf[2]*fUpwind[27]; 
  drag_incr[41] = 0.223606797749979*alphaDrSurf[12]*fUpwind[41]+0.223606797749979*alphaDrSurf[11]*fUpwind[41]+0.25*alphaDrSurf[0]*fUpwind[41]+0.223606797749979*alphaDrSurf[20]*fUpwind[29]+0.2500000000000001*alphaDrSurf[1]*fUpwind[29]+0.223606797749979*alphaDrSurf[19]*fUpwind[28]+0.2500000000000001*alphaDrSurf[2]*fUpwind[28]+0.25*alphaDrSurf[5]*fUpwind[14]; 
  drag_incr[42] = 0.223606797749979*alphaDrSurf[19]*fUpwind[47]+0.2500000000000001*alphaDrSurf[2]*fUpwind[47]+0.25*alphaDrSurf[5]*fUpwind[43]+0.223606797749979*alphaDrSurf[11]*fUpwind[42]+0.25*alphaDrSurf[0]*fUpwind[42]+0.2500000000000001*alphaDrSurf[1]*fUpwind[30]; 
  drag_incr[43] = 0.223606797749979*alphaDrSurf[20]*fUpwind[47]+0.2500000000000001*alphaDrSurf[1]*fUpwind[47]+0.223606797749979*alphaDrSurf[12]*fUpwind[43]+0.25*alphaDrSurf[0]*fUpwind[43]+0.25*alphaDrSurf[5]*fUpwind[42]+0.2500000000000001*alphaDrSurf[2]*fUpwind[30]; 
  drag_incr[44] = 0.2*alphaDrSurf[5]*fUpwind[45]+0.223606797749979*alphaDrSurf[12]*fUpwind[44]+0.159719141249985*alphaDrSurf[11]*fUpwind[44]+0.25*alphaDrSurf[0]*fUpwind[44]+0.223606797749979*alphaDrSurf[19]*fUpwind[38]+0.159719141249985*alphaDrSurf[19]*fUpwind[37]+0.2500000000000001*alphaDrSurf[2]*fUpwind[37]+0.2*alphaDrSurf[20]*fUpwind[31]+0.223606797749979*alphaDrSurf[1]*fUpwind[31]+0.25*fUpwind[10]*alphaDrSurf[19]+0.2500000000000001*alphaDrSurf[11]*fUpwind[18]+0.223606797749979*alphaDrSurf[5]*fUpwind[17]; 
  drag_incr[45] = 0.159719141249985*alphaDrSurf[12]*fUpwind[45]+0.223606797749979*alphaDrSurf[11]*fUpwind[45]+0.25*alphaDrSurf[0]*fUpwind[45]+0.2*alphaDrSurf[5]*fUpwind[44]+0.159719141249985*alphaDrSurf[20]*fUpwind[38]+0.2500000000000001*alphaDrSurf[1]*fUpwind[38]+0.223606797749979*alphaDrSurf[20]*fUpwind[37]+0.2*alphaDrSurf[19]*fUpwind[31]+0.223606797749979*alphaDrSurf[2]*fUpwind[31]+0.25*fUpwind[10]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[18]+0.2500000000000001*alphaDrSurf[12]*fUpwind[17]; 
  drag_incr[46] = 0.223606797749979*alphaDrSurf[12]*fUpwind[46]+0.223606797749979*alphaDrSurf[11]*fUpwind[46]+0.25*alphaDrSurf[0]*fUpwind[46]+0.223606797749979*alphaDrSurf[20]*fUpwind[40]+0.2500000000000001*alphaDrSurf[1]*fUpwind[40]+0.223606797749979*alphaDrSurf[19]*fUpwind[39]+0.2500000000000001*alphaDrSurf[2]*fUpwind[39]+0.25*alphaDrSurf[5]*fUpwind[27]; 
  drag_incr[47] = 0.223606797749979*alphaDrSurf[12]*fUpwind[47]+0.223606797749979*alphaDrSurf[11]*fUpwind[47]+0.25*alphaDrSurf[0]*fUpwind[47]+0.223606797749979*alphaDrSurf[20]*fUpwind[43]+0.2500000000000001*alphaDrSurf[1]*fUpwind[43]+0.223606797749979*alphaDrSurf[19]*fUpwind[42]+0.2500000000000001*alphaDrSurf[2]*fUpwind[42]+0.25*alphaDrSurf[5]*fUpwind[30]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[11] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[12] += 0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += 0.7071067811865475*drag_incr[9]*rdv2; 
  out[14] += 0.7071067811865475*drag_incr[10]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[16] += 0.7071067811865475*drag_incr[11]*rdv2; 
  out[17] += 0.7071067811865475*drag_incr[12]*rdv2; 
  out[18] += 0.7071067811865475*drag_incr[13]*rdv2; 
  out[19] += 1.58113883008419*drag_incr[0]*rdv2; 
  out[20] += 0.7071067811865475*drag_incr[14]*rdv2; 
  out[21] += 0.7071067811865475*drag_incr[15]*rdv2; 
  out[22] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[23] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[24] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[25] += 0.7071067811865475*drag_incr[16]*rdv2; 
  out[26] += 0.7071067811865475*drag_incr[17]*rdv2; 
  out[27] += 0.7071067811865475*drag_incr[18]*rdv2; 
  out[28] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[31] += 0.7071067811865475*drag_incr[19]*rdv2; 
  out[32] += 0.7071067811865475*drag_incr[20]*rdv2; 
  out[33] += 0.7071067811865475*drag_incr[21]*rdv2; 
  out[34] += 0.7071067811865475*drag_incr[22]*rdv2; 
  out[35] += 0.7071067811865475*drag_incr[23]*rdv2; 
  out[36] += 0.7071067811865475*drag_incr[24]*rdv2; 
  out[37] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[38] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[39] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[40] += 1.58113883008419*drag_incr[1]*rdv2; 
  out[41] += 1.58113883008419*drag_incr[2]*rdv2; 
  out[42] += 1.58113883008419*drag_incr[3]*rdv2; 
  out[43] += 0.7071067811865475*drag_incr[25]*rdv2; 
  out[44] += 0.7071067811865475*drag_incr[26]*rdv2; 
  out[45] += 0.7071067811865475*drag_incr[27]*rdv2; 
  out[46] += 1.58113883008419*drag_incr[4]*rdv2; 
  out[47] += 0.7071067811865475*drag_incr[28]*rdv2; 
  out[48] += 0.7071067811865475*drag_incr[29]*rdv2; 
  out[49] += 0.7071067811865475*drag_incr[30]*rdv2; 
  out[50] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[51] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[52] += 0.7071067811865475*drag_incr[31]*rdv2; 
  out[53] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[54] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[55] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[56] += 0.7071067811865475*drag_incr[32]*rdv2; 
  out[57] += 0.7071067811865475*drag_incr[33]*rdv2; 
  out[58] += 0.7071067811865475*drag_incr[34]*rdv2; 
  out[59] += 1.224744871391589*drag_incr[19]*rdv2; 
  out[60] += 1.224744871391589*drag_incr[20]*rdv2; 
  out[61] += 1.224744871391589*drag_incr[21]*rdv2; 
  out[62] += 1.224744871391589*drag_incr[22]*rdv2; 
  out[63] += 1.224744871391589*drag_incr[23]*rdv2; 
  out[64] += 1.224744871391589*drag_incr[24]*rdv2; 
  out[65] += 1.58113883008419*drag_incr[5]*rdv2; 
  out[66] += 1.58113883008419*drag_incr[6]*rdv2; 
  out[67] += 1.58113883008419*drag_incr[7]*rdv2; 
  out[68] += 0.7071067811865475*drag_incr[35]*rdv2; 
  out[69] += 0.7071067811865475*drag_incr[36]*rdv2; 
  out[70] += 0.7071067811865475*drag_incr[37]*rdv2; 
  out[71] += 0.7071067811865475*drag_incr[38]*rdv2; 
  out[72] += 0.7071067811865475*drag_incr[39]*rdv2; 
  out[73] += 0.7071067811865475*drag_incr[40]*rdv2; 
  out[74] += 1.224744871391589*drag_incr[25]*rdv2; 
  out[75] += 1.224744871391589*drag_incr[26]*rdv2; 
  out[76] += 1.224744871391589*drag_incr[27]*rdv2; 
  out[77] += 1.58113883008419*drag_incr[8]*rdv2; 
  out[78] += 1.58113883008419*drag_incr[9]*rdv2; 
  out[79] += 1.58113883008419*drag_incr[10]*rdv2; 
  out[80] += 0.7071067811865475*drag_incr[41]*rdv2; 
  out[81] += 0.7071067811865475*drag_incr[42]*rdv2; 
  out[82] += 0.7071067811865475*drag_incr[43]*rdv2; 
  out[83] += 1.224744871391589*drag_incr[28]*rdv2; 
  out[84] += 1.224744871391589*drag_incr[29]*rdv2; 
  out[85] += 1.224744871391589*drag_incr[30]*rdv2; 
  out[86] += 1.224744871391589*drag_incr[31]*rdv2; 
  out[87] += 1.224744871391589*drag_incr[32]*rdv2; 
  out[88] += 1.224744871391589*drag_incr[33]*rdv2; 
  out[89] += 1.224744871391589*drag_incr[34]*rdv2; 
  out[90] += 1.58113883008419*drag_incr[15]*rdv2; 
  out[91] += 0.7071067811865475*drag_incr[44]*rdv2; 
  out[92] += 0.7071067811865475*drag_incr[45]*rdv2; 
  out[93] += 0.7071067811865475*drag_incr[46]*rdv2; 
  out[94] += 1.224744871391589*drag_incr[35]*rdv2; 
  out[95] += 1.224744871391589*drag_incr[36]*rdv2; 
  out[96] += 1.224744871391589*drag_incr[37]*rdv2; 
  out[97] += 1.224744871391589*drag_incr[38]*rdv2; 
  out[98] += 1.224744871391589*drag_incr[39]*rdv2; 
  out[99] += 1.224744871391589*drag_incr[40]*rdv2; 
  out[100] += 1.58113883008419*drag_incr[16]*rdv2; 
  out[101] += 1.58113883008419*drag_incr[17]*rdv2; 
  out[102] += 1.58113883008419*drag_incr[18]*rdv2; 
  out[103] += 0.7071067811865475*drag_incr[47]*rdv2; 
  out[104] += 1.224744871391589*drag_incr[41]*rdv2; 
  out[105] += 1.224744871391589*drag_incr[42]*rdv2; 
  out[106] += 1.224744871391589*drag_incr[43]*rdv2; 
  out[107] += 1.224744871391589*drag_incr[44]*rdv2; 
  out[108] += 1.224744871391589*drag_incr[45]*rdv2; 
  out[109] += 1.224744871391589*drag_incr[46]*rdv2; 
  out[110] += 1.58113883008419*drag_incr[31]*rdv2; 
  out[111] += 1.224744871391589*drag_incr[47]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[0]; 
  alphaDrSurf[1] = nuSum[1]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[1]; 
  alphaDrSurf[2] = nuSum[2]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[2]; 
  alphaDrSurf[5] = 2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]-1.0*dxv[3]*nuSum[3]; 
  alphaDrSurf[11] = (2.0*w[3]-1.0*dxv[3])*nuSum[4]-2.0*sumNuUy[4]; 
  alphaDrSurf[12] = (2.0*w[3]-1.0*dxv[3])*nuSum[5]-2.0*sumNuUy[5]; 
  alphaDrSurf[19] = (2.0*w[3]-1.0*dxv[3])*nuSum[6]-2.0*sumNuUy[6]; 
  alphaDrSurf[20] = (2.0*w[3]-1.0*dxv[3])*nuSum[7]-2.0*sumNuUy[7]; 

  if ((-0.2999999999999998*alphaDrSurf[20])-0.2999999999999999*alphaDrSurf[19]+0.223606797749979*(alphaDrSurf[12]+alphaDrSurf[11])+0.45*alphaDrSurf[5]-0.3354101966249685*(alphaDrSurf[2]+alphaDrSurf[1])+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_5x_p2_surfx4_quad_0_r(fEdge); 
    fUpwindQuad[9] = ser_5x_p2_surfx4_quad_9_r(fEdge); 
    fUpwindQuad[18] = ser_5x_p2_surfx4_quad_18_r(fEdge); 
    fUpwindQuad[27] = ser_5x_p2_surfx4_quad_27_r(fEdge); 
    fUpwindQuad[36] = ser_5x_p2_surfx4_quad_36_r(fEdge); 
    fUpwindQuad[45] = ser_5x_p2_surfx4_quad_45_r(fEdge); 
    fUpwindQuad[54] = ser_5x_p2_surfx4_quad_54_r(fEdge); 
    fUpwindQuad[63] = ser_5x_p2_surfx4_quad_63_r(fEdge); 
    fUpwindQuad[72] = ser_5x_p2_surfx4_quad_72_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_5x_p2_surfx4_quad_0_l(fSkin); 
    fUpwindQuad[9] = ser_5x_p2_surfx4_quad_9_l(fSkin); 
    fUpwindQuad[18] = ser_5x_p2_surfx4_quad_18_l(fSkin); 
    fUpwindQuad[27] = ser_5x_p2_surfx4_quad_27_l(fSkin); 
    fUpwindQuad[36] = ser_5x_p2_surfx4_quad_36_l(fSkin); 
    fUpwindQuad[45] = ser_5x_p2_surfx4_quad_45_l(fSkin); 
    fUpwindQuad[54] = ser_5x_p2_surfx4_quad_54_l(fSkin); 
    fUpwindQuad[63] = ser_5x_p2_surfx4_quad_63_l(fSkin); 
    fUpwindQuad[72] = ser_5x_p2_surfx4_quad_72_l(fSkin); 
  } 
  if (0.375*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[12]-0.2795084971874737*alphaDrSurf[11]-0.3354101966249685*alphaDrSurf[2]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_5x_p2_surfx4_quad_1_r(fEdge); 
    fUpwindQuad[10] = ser_5x_p2_surfx4_quad_10_r(fEdge); 
    fUpwindQuad[19] = ser_5x_p2_surfx4_quad_19_r(fEdge); 
    fUpwindQuad[28] = ser_5x_p2_surfx4_quad_28_r(fEdge); 
    fUpwindQuad[37] = ser_5x_p2_surfx4_quad_37_r(fEdge); 
    fUpwindQuad[46] = ser_5x_p2_surfx4_quad_46_r(fEdge); 
    fUpwindQuad[55] = ser_5x_p2_surfx4_quad_55_r(fEdge); 
    fUpwindQuad[64] = ser_5x_p2_surfx4_quad_64_r(fEdge); 
    fUpwindQuad[73] = ser_5x_p2_surfx4_quad_73_r(fEdge); 
  } else { 
    fUpwindQuad[1] = ser_5x_p2_surfx4_quad_1_l(fSkin); 
    fUpwindQuad[10] = ser_5x_p2_surfx4_quad_10_l(fSkin); 
    fUpwindQuad[19] = ser_5x_p2_surfx4_quad_19_l(fSkin); 
    fUpwindQuad[28] = ser_5x_p2_surfx4_quad_28_l(fSkin); 
    fUpwindQuad[37] = ser_5x_p2_surfx4_quad_37_l(fSkin); 
    fUpwindQuad[46] = ser_5x_p2_surfx4_quad_46_l(fSkin); 
    fUpwindQuad[55] = ser_5x_p2_surfx4_quad_55_l(fSkin); 
    fUpwindQuad[64] = ser_5x_p2_surfx4_quad_64_l(fSkin); 
    fUpwindQuad[73] = ser_5x_p2_surfx4_quad_73_l(fSkin); 
  } 
  if (0.2999999999999998*alphaDrSurf[20]-0.2999999999999999*alphaDrSurf[19]+0.223606797749979*(alphaDrSurf[12]+alphaDrSurf[11])-0.45*alphaDrSurf[5]-0.3354101966249685*alphaDrSurf[2]+0.3354101966249685*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_5x_p2_surfx4_quad_2_r(fEdge); 
    fUpwindQuad[11] = ser_5x_p2_surfx4_quad_11_r(fEdge); 
    fUpwindQuad[20] = ser_5x_p2_surfx4_quad_20_r(fEdge); 
    fUpwindQuad[29] = ser_5x_p2_surfx4_quad_29_r(fEdge); 
    fUpwindQuad[38] = ser_5x_p2_surfx4_quad_38_r(fEdge); 
    fUpwindQuad[47] = ser_5x_p2_surfx4_quad_47_r(fEdge); 
    fUpwindQuad[56] = ser_5x_p2_surfx4_quad_56_r(fEdge); 
    fUpwindQuad[65] = ser_5x_p2_surfx4_quad_65_r(fEdge); 
    fUpwindQuad[74] = ser_5x_p2_surfx4_quad_74_r(fEdge); 
  } else { 
    fUpwindQuad[2] = ser_5x_p2_surfx4_quad_2_l(fSkin); 
    fUpwindQuad[11] = ser_5x_p2_surfx4_quad_11_l(fSkin); 
    fUpwindQuad[20] = ser_5x_p2_surfx4_quad_20_l(fSkin); 
    fUpwindQuad[29] = ser_5x_p2_surfx4_quad_29_l(fSkin); 
    fUpwindQuad[38] = ser_5x_p2_surfx4_quad_38_l(fSkin); 
    fUpwindQuad[47] = ser_5x_p2_surfx4_quad_47_l(fSkin); 
    fUpwindQuad[56] = ser_5x_p2_surfx4_quad_56_l(fSkin); 
    fUpwindQuad[65] = ser_5x_p2_surfx4_quad_65_l(fSkin); 
    fUpwindQuad[74] = ser_5x_p2_surfx4_quad_74_l(fSkin); 
  } 
  if (0.375*alphaDrSurf[20]-0.2795084971874737*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[11]-0.3354101966249685*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_5x_p2_surfx4_quad_3_r(fEdge); 
    fUpwindQuad[12] = ser_5x_p2_surfx4_quad_12_r(fEdge); 
    fUpwindQuad[21] = ser_5x_p2_surfx4_quad_21_r(fEdge); 
    fUpwindQuad[30] = ser_5x_p2_surfx4_quad_30_r(fEdge); 
    fUpwindQuad[39] = ser_5x_p2_surfx4_quad_39_r(fEdge); 
    fUpwindQuad[48] = ser_5x_p2_surfx4_quad_48_r(fEdge); 
    fUpwindQuad[57] = ser_5x_p2_surfx4_quad_57_r(fEdge); 
    fUpwindQuad[66] = ser_5x_p2_surfx4_quad_66_r(fEdge); 
    fUpwindQuad[75] = ser_5x_p2_surfx4_quad_75_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_5x_p2_surfx4_quad_3_l(fSkin); 
    fUpwindQuad[12] = ser_5x_p2_surfx4_quad_12_l(fSkin); 
    fUpwindQuad[21] = ser_5x_p2_surfx4_quad_21_l(fSkin); 
    fUpwindQuad[30] = ser_5x_p2_surfx4_quad_30_l(fSkin); 
    fUpwindQuad[39] = ser_5x_p2_surfx4_quad_39_l(fSkin); 
    fUpwindQuad[48] = ser_5x_p2_surfx4_quad_48_l(fSkin); 
    fUpwindQuad[57] = ser_5x_p2_surfx4_quad_57_l(fSkin); 
    fUpwindQuad[66] = ser_5x_p2_surfx4_quad_66_l(fSkin); 
    fUpwindQuad[75] = ser_5x_p2_surfx4_quad_75_l(fSkin); 
  } 
  if (0.25*alphaDrSurf[0]-0.2795084971874737*(alphaDrSurf[12]+alphaDrSurf[11]) < 0) { 
    fUpwindQuad[4] = ser_5x_p2_surfx4_quad_4_r(fEdge); 
    fUpwindQuad[13] = ser_5x_p2_surfx4_quad_13_r(fEdge); 
    fUpwindQuad[22] = ser_5x_p2_surfx4_quad_22_r(fEdge); 
    fUpwindQuad[31] = ser_5x_p2_surfx4_quad_31_r(fEdge); 
    fUpwindQuad[40] = ser_5x_p2_surfx4_quad_40_r(fEdge); 
    fUpwindQuad[49] = ser_5x_p2_surfx4_quad_49_r(fEdge); 
    fUpwindQuad[58] = ser_5x_p2_surfx4_quad_58_r(fEdge); 
    fUpwindQuad[67] = ser_5x_p2_surfx4_quad_67_r(fEdge); 
    fUpwindQuad[76] = ser_5x_p2_surfx4_quad_76_r(fEdge); 
  } else { 
    fUpwindQuad[4] = ser_5x_p2_surfx4_quad_4_l(fSkin); 
    fUpwindQuad[13] = ser_5x_p2_surfx4_quad_13_l(fSkin); 
    fUpwindQuad[22] = ser_5x_p2_surfx4_quad_22_l(fSkin); 
    fUpwindQuad[31] = ser_5x_p2_surfx4_quad_31_l(fSkin); 
    fUpwindQuad[40] = ser_5x_p2_surfx4_quad_40_l(fSkin); 
    fUpwindQuad[49] = ser_5x_p2_surfx4_quad_49_l(fSkin); 
    fUpwindQuad[58] = ser_5x_p2_surfx4_quad_58_l(fSkin); 
    fUpwindQuad[67] = ser_5x_p2_surfx4_quad_67_l(fSkin); 
    fUpwindQuad[76] = ser_5x_p2_surfx4_quad_76_l(fSkin); 
  } 
  if ((-0.375*alphaDrSurf[20])-0.2795084971874737*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[11]+0.3354101966249685*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_5x_p2_surfx4_quad_5_r(fEdge); 
    fUpwindQuad[14] = ser_5x_p2_surfx4_quad_14_r(fEdge); 
    fUpwindQuad[23] = ser_5x_p2_surfx4_quad_23_r(fEdge); 
    fUpwindQuad[32] = ser_5x_p2_surfx4_quad_32_r(fEdge); 
    fUpwindQuad[41] = ser_5x_p2_surfx4_quad_41_r(fEdge); 
    fUpwindQuad[50] = ser_5x_p2_surfx4_quad_50_r(fEdge); 
    fUpwindQuad[59] = ser_5x_p2_surfx4_quad_59_r(fEdge); 
    fUpwindQuad[68] = ser_5x_p2_surfx4_quad_68_r(fEdge); 
    fUpwindQuad[77] = ser_5x_p2_surfx4_quad_77_r(fEdge); 
  } else { 
    fUpwindQuad[5] = ser_5x_p2_surfx4_quad_5_l(fSkin); 
    fUpwindQuad[14] = ser_5x_p2_surfx4_quad_14_l(fSkin); 
    fUpwindQuad[23] = ser_5x_p2_surfx4_quad_23_l(fSkin); 
    fUpwindQuad[32] = ser_5x_p2_surfx4_quad_32_l(fSkin); 
    fUpwindQuad[41] = ser_5x_p2_surfx4_quad_41_l(fSkin); 
    fUpwindQuad[50] = ser_5x_p2_surfx4_quad_50_l(fSkin); 
    fUpwindQuad[59] = ser_5x_p2_surfx4_quad_59_l(fSkin); 
    fUpwindQuad[68] = ser_5x_p2_surfx4_quad_68_l(fSkin); 
    fUpwindQuad[77] = ser_5x_p2_surfx4_quad_77_l(fSkin); 
  } 
  if ((-0.2999999999999998*alphaDrSurf[20])+0.2999999999999999*alphaDrSurf[19]+0.223606797749979*(alphaDrSurf[12]+alphaDrSurf[11])-0.45*alphaDrSurf[5]+0.3354101966249685*alphaDrSurf[2]-0.3354101966249685*alphaDrSurf[1]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_5x_p2_surfx4_quad_6_r(fEdge); 
    fUpwindQuad[15] = ser_5x_p2_surfx4_quad_15_r(fEdge); 
    fUpwindQuad[24] = ser_5x_p2_surfx4_quad_24_r(fEdge); 
    fUpwindQuad[33] = ser_5x_p2_surfx4_quad_33_r(fEdge); 
    fUpwindQuad[42] = ser_5x_p2_surfx4_quad_42_r(fEdge); 
    fUpwindQuad[51] = ser_5x_p2_surfx4_quad_51_r(fEdge); 
    fUpwindQuad[60] = ser_5x_p2_surfx4_quad_60_r(fEdge); 
    fUpwindQuad[69] = ser_5x_p2_surfx4_quad_69_r(fEdge); 
    fUpwindQuad[78] = ser_5x_p2_surfx4_quad_78_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_5x_p2_surfx4_quad_6_l(fSkin); 
    fUpwindQuad[15] = ser_5x_p2_surfx4_quad_15_l(fSkin); 
    fUpwindQuad[24] = ser_5x_p2_surfx4_quad_24_l(fSkin); 
    fUpwindQuad[33] = ser_5x_p2_surfx4_quad_33_l(fSkin); 
    fUpwindQuad[42] = ser_5x_p2_surfx4_quad_42_l(fSkin); 
    fUpwindQuad[51] = ser_5x_p2_surfx4_quad_51_l(fSkin); 
    fUpwindQuad[60] = ser_5x_p2_surfx4_quad_60_l(fSkin); 
    fUpwindQuad[69] = ser_5x_p2_surfx4_quad_69_l(fSkin); 
    fUpwindQuad[78] = ser_5x_p2_surfx4_quad_78_l(fSkin); 
  } 
  if ((-0.375*alphaDrSurf[19])+0.223606797749979*alphaDrSurf[12]-0.2795084971874737*alphaDrSurf[11]+0.3354101966249685*alphaDrSurf[2]+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_5x_p2_surfx4_quad_7_r(fEdge); 
    fUpwindQuad[16] = ser_5x_p2_surfx4_quad_16_r(fEdge); 
    fUpwindQuad[25] = ser_5x_p2_surfx4_quad_25_r(fEdge); 
    fUpwindQuad[34] = ser_5x_p2_surfx4_quad_34_r(fEdge); 
    fUpwindQuad[43] = ser_5x_p2_surfx4_quad_43_r(fEdge); 
    fUpwindQuad[52] = ser_5x_p2_surfx4_quad_52_r(fEdge); 
    fUpwindQuad[61] = ser_5x_p2_surfx4_quad_61_r(fEdge); 
    fUpwindQuad[70] = ser_5x_p2_surfx4_quad_70_r(fEdge); 
    fUpwindQuad[79] = ser_5x_p2_surfx4_quad_79_r(fEdge); 
  } else { 
    fUpwindQuad[7] = ser_5x_p2_surfx4_quad_7_l(fSkin); 
    fUpwindQuad[16] = ser_5x_p2_surfx4_quad_16_l(fSkin); 
    fUpwindQuad[25] = ser_5x_p2_surfx4_quad_25_l(fSkin); 
    fUpwindQuad[34] = ser_5x_p2_surfx4_quad_34_l(fSkin); 
    fUpwindQuad[43] = ser_5x_p2_surfx4_quad_43_l(fSkin); 
    fUpwindQuad[52] = ser_5x_p2_surfx4_quad_52_l(fSkin); 
    fUpwindQuad[61] = ser_5x_p2_surfx4_quad_61_l(fSkin); 
    fUpwindQuad[70] = ser_5x_p2_surfx4_quad_70_l(fSkin); 
    fUpwindQuad[79] = ser_5x_p2_surfx4_quad_79_l(fSkin); 
  } 
  if (0.2999999999999998*alphaDrSurf[20]+0.2999999999999999*alphaDrSurf[19]+0.223606797749979*(alphaDrSurf[12]+alphaDrSurf[11])+0.45*alphaDrSurf[5]+0.3354101966249685*(alphaDrSurf[2]+alphaDrSurf[1])+0.25*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_5x_p2_surfx4_quad_8_r(fEdge); 
    fUpwindQuad[17] = ser_5x_p2_surfx4_quad_17_r(fEdge); 
    fUpwindQuad[26] = ser_5x_p2_surfx4_quad_26_r(fEdge); 
    fUpwindQuad[35] = ser_5x_p2_surfx4_quad_35_r(fEdge); 
    fUpwindQuad[44] = ser_5x_p2_surfx4_quad_44_r(fEdge); 
    fUpwindQuad[53] = ser_5x_p2_surfx4_quad_53_r(fEdge); 
    fUpwindQuad[62] = ser_5x_p2_surfx4_quad_62_r(fEdge); 
    fUpwindQuad[71] = ser_5x_p2_surfx4_quad_71_r(fEdge); 
    fUpwindQuad[80] = ser_5x_p2_surfx4_quad_80_r(fEdge); 
  } else { 
    fUpwindQuad[8] = ser_5x_p2_surfx4_quad_8_l(fSkin); 
    fUpwindQuad[17] = ser_5x_p2_surfx4_quad_17_l(fSkin); 
    fUpwindQuad[26] = ser_5x_p2_surfx4_quad_26_l(fSkin); 
    fUpwindQuad[35] = ser_5x_p2_surfx4_quad_35_l(fSkin); 
    fUpwindQuad[44] = ser_5x_p2_surfx4_quad_44_l(fSkin); 
    fUpwindQuad[53] = ser_5x_p2_surfx4_quad_53_l(fSkin); 
    fUpwindQuad[62] = ser_5x_p2_surfx4_quad_62_l(fSkin); 
    fUpwindQuad[71] = ser_5x_p2_surfx4_quad_71_l(fSkin); 
    fUpwindQuad[80] = ser_5x_p2_surfx4_quad_80_l(fSkin); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p2_upwind(fUpwindQuad, fUpwind); 

  drag_incr[0] = 0.25*alphaDrSurf[20]*fUpwind[20]+0.25*alphaDrSurf[19]*fUpwind[19]+0.25*alphaDrSurf[12]*fUpwind[12]+0.25*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.2500000000000001*alphaDrSurf[12]*fUpwind[20]+0.2500000000000001*fUpwind[12]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[19]+0.223606797749979*fUpwind[5]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[11]+0.223606797749979*fUpwind[1]*alphaDrSurf[11]+0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.223606797749979*alphaDrSurf[5]*fUpwind[20]+0.223606797749979*fUpwind[5]*alphaDrSurf[20]+0.2500000000000001*alphaDrSurf[11]*fUpwind[19]+0.2500000000000001*fUpwind[11]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[2]*fUpwind[12]+0.223606797749979*fUpwind[2]*alphaDrSurf[12]+0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.2500000000000001*alphaDrSurf[20]*fUpwind[33]+0.2500000000000001*alphaDrSurf[19]*fUpwind[32]+0.2500000000000001*alphaDrSurf[12]*fUpwind[22]+0.2500000000000001*alphaDrSurf[11]*fUpwind[21]+0.25*alphaDrSurf[5]*fUpwind[15]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.2500000000000001*alphaDrSurf[20]*fUpwind[36]+0.2500000000000001*alphaDrSurf[19]*fUpwind[35]+0.2500000000000001*alphaDrSurf[12]*fUpwind[26]+0.2500000000000001*alphaDrSurf[11]*fUpwind[25]+0.25*alphaDrSurf[5]*fUpwind[16]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  drag_incr[5] = 0.2*alphaDrSurf[19]*fUpwind[20]+0.223606797749979*alphaDrSurf[2]*fUpwind[20]+0.2*fUpwind[19]*alphaDrSurf[20]+0.223606797749979*fUpwind[2]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[1]*fUpwind[19]+0.223606797749979*fUpwind[1]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[5]*fUpwind[12]+0.223606797749979*fUpwind[5]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[11]+0.223606797749979*fUpwind[5]*alphaDrSurf[11]+0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[6] = 0.25*alphaDrSurf[12]*fUpwind[33]+0.223606797749979*alphaDrSurf[5]*fUpwind[32]+0.25*alphaDrSurf[20]*fUpwind[22]+0.223606797749979*alphaDrSurf[1]*fUpwind[21]+0.223606797749979*fUpwind[15]*alphaDrSurf[19]+0.25*alphaDrSurf[2]*fUpwind[15]+0.223606797749979*fUpwind[6]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[7] = 0.223606797749979*alphaDrSurf[5]*fUpwind[33]+0.25*alphaDrSurf[11]*fUpwind[32]+0.223606797749979*alphaDrSurf[2]*fUpwind[22]+0.25*alphaDrSurf[19]*fUpwind[21]+0.223606797749979*fUpwind[15]*alphaDrSurf[20]+0.25*alphaDrSurf[1]*fUpwind[15]+0.223606797749979*fUpwind[7]*alphaDrSurf[12]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*alphaDrSurf[2]*fUpwind[3]; 
  drag_incr[8] = 0.25*alphaDrSurf[12]*fUpwind[36]+0.223606797749979*alphaDrSurf[5]*fUpwind[35]+0.25*alphaDrSurf[20]*fUpwind[26]+0.223606797749979*alphaDrSurf[1]*fUpwind[25]+0.223606797749979*fUpwind[16]*alphaDrSurf[19]+0.25*alphaDrSurf[2]*fUpwind[16]+0.223606797749979*fUpwind[8]*alphaDrSurf[11]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[9] = 0.223606797749979*alphaDrSurf[5]*fUpwind[36]+0.25*alphaDrSurf[11]*fUpwind[35]+0.223606797749979*alphaDrSurf[2]*fUpwind[26]+0.25*alphaDrSurf[19]*fUpwind[25]+0.223606797749979*fUpwind[16]*alphaDrSurf[20]+0.25*alphaDrSurf[1]*fUpwind[16]+0.223606797749979*fUpwind[9]*alphaDrSurf[12]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  drag_incr[10] = 0.25*alphaDrSurf[20]*fUpwind[45]+0.25*alphaDrSurf[19]*fUpwind[44]+0.25*alphaDrSurf[12]*fUpwind[38]+0.25*alphaDrSurf[11]*fUpwind[37]+0.25*alphaDrSurf[5]*fUpwind[31]+0.25*alphaDrSurf[2]*fUpwind[18]+0.25*alphaDrSurf[1]*fUpwind[17]+0.25*alphaDrSurf[0]*fUpwind[10]; 
  drag_incr[11] = 0.223606797749979*alphaDrSurf[20]*fUpwind[20]+0.159719141249985*alphaDrSurf[19]*fUpwind[19]+0.2500000000000001*alphaDrSurf[2]*fUpwind[19]+0.2500000000000001*fUpwind[2]*alphaDrSurf[19]+0.159719141249985*alphaDrSurf[11]*fUpwind[11]+0.25*alphaDrSurf[0]*fUpwind[11]+0.25*fUpwind[0]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[5]+0.223606797749979*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[12] = 0.159719141249985*alphaDrSurf[20]*fUpwind[20]+0.2500000000000001*alphaDrSurf[1]*fUpwind[20]+0.2500000000000001*fUpwind[1]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[19]*fUpwind[19]+0.159719141249985*alphaDrSurf[12]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[12]+0.25*fUpwind[0]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[5]+0.223606797749979*alphaDrSurf[2]*fUpwind[2]; 
  drag_incr[13] = 0.25*alphaDrSurf[5]*fUpwind[34]+0.2500000000000001*alphaDrSurf[2]*fUpwind[24]+0.2500000000000001*alphaDrSurf[1]*fUpwind[23]+0.25*alphaDrSurf[0]*fUpwind[13]; 
  drag_incr[14] = 0.25*alphaDrSurf[5]*fUpwind[41]+0.2500000000000001*alphaDrSurf[2]*fUpwind[29]+0.2500000000000001*alphaDrSurf[1]*fUpwind[28]+0.25*alphaDrSurf[0]*fUpwind[14]; 
  drag_incr[15] = 0.2*alphaDrSurf[19]*fUpwind[33]+0.223606797749979*alphaDrSurf[2]*fUpwind[33]+0.2*alphaDrSurf[20]*fUpwind[32]+0.223606797749979*alphaDrSurf[1]*fUpwind[32]+0.223606797749979*alphaDrSurf[5]*fUpwind[22]+0.223606797749979*alphaDrSurf[5]*fUpwind[21]+0.223606797749979*fUpwind[7]*alphaDrSurf[20]+0.223606797749979*fUpwind[6]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[12]*fUpwind[15]+0.223606797749979*alphaDrSurf[11]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[15]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  drag_incr[16] = 0.2*alphaDrSurf[19]*fUpwind[36]+0.223606797749979*alphaDrSurf[2]*fUpwind[36]+0.2*alphaDrSurf[20]*fUpwind[35]+0.223606797749979*alphaDrSurf[1]*fUpwind[35]+0.223606797749979*alphaDrSurf[5]*fUpwind[26]+0.223606797749979*alphaDrSurf[5]*fUpwind[25]+0.223606797749979*fUpwind[9]*alphaDrSurf[20]+0.223606797749979*fUpwind[8]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[12]*fUpwind[16]+0.223606797749979*alphaDrSurf[11]*fUpwind[16]+0.25*alphaDrSurf[0]*fUpwind[16]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  drag_incr[17] = 0.2500000000000001*alphaDrSurf[12]*fUpwind[45]+0.223606797749979*alphaDrSurf[5]*fUpwind[44]+0.2500000000000001*alphaDrSurf[20]*fUpwind[38]+0.223606797749979*alphaDrSurf[1]*fUpwind[37]+0.223606797749979*alphaDrSurf[19]*fUpwind[31]+0.25*alphaDrSurf[2]*fUpwind[31]+0.25*alphaDrSurf[5]*fUpwind[18]+0.223606797749979*alphaDrSurf[11]*fUpwind[17]+0.25*alphaDrSurf[0]*fUpwind[17]+0.25*alphaDrSurf[1]*fUpwind[10]; 
  drag_incr[18] = 0.223606797749979*alphaDrSurf[5]*fUpwind[45]+0.2500000000000001*alphaDrSurf[11]*fUpwind[44]+0.223606797749979*alphaDrSurf[2]*fUpwind[38]+0.2500000000000001*alphaDrSurf[19]*fUpwind[37]+0.223606797749979*alphaDrSurf[20]*fUpwind[31]+0.25*alphaDrSurf[1]*fUpwind[31]+0.223606797749979*alphaDrSurf[12]*fUpwind[18]+0.25*alphaDrSurf[0]*fUpwind[18]+0.25*alphaDrSurf[5]*fUpwind[17]+0.25*alphaDrSurf[2]*fUpwind[10]; 
  drag_incr[19] = 0.2*alphaDrSurf[5]*fUpwind[20]+0.2*fUpwind[5]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[12]*fUpwind[19]+0.159719141249985*alphaDrSurf[11]*fUpwind[19]+0.25*alphaDrSurf[0]*fUpwind[19]+0.223606797749979*fUpwind[12]*alphaDrSurf[19]+0.159719141249985*fUpwind[11]*alphaDrSurf[19]+0.25*fUpwind[0]*alphaDrSurf[19]+0.2500000000000001*alphaDrSurf[2]*fUpwind[11]+0.2500000000000001*fUpwind[2]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[5]+0.223606797749979*fUpwind[1]*alphaDrSurf[5]; 
  drag_incr[20] = 0.159719141249985*alphaDrSurf[12]*fUpwind[20]+0.223606797749979*alphaDrSurf[11]*fUpwind[20]+0.25*alphaDrSurf[0]*fUpwind[20]+0.159719141249985*fUpwind[12]*alphaDrSurf[20]+0.223606797749979*fUpwind[11]*alphaDrSurf[20]+0.25*fUpwind[0]*alphaDrSurf[20]+0.2*alphaDrSurf[5]*fUpwind[19]+0.2*fUpwind[5]*alphaDrSurf[19]+0.2500000000000001*alphaDrSurf[1]*fUpwind[12]+0.2500000000000001*fUpwind[1]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[5]+0.223606797749979*fUpwind[2]*alphaDrSurf[5]; 
  drag_incr[21] = 0.223606797749979*alphaDrSurf[20]*fUpwind[33]+0.159719141249985*alphaDrSurf[19]*fUpwind[32]+0.2500000000000001*alphaDrSurf[2]*fUpwind[32]+0.159719141249985*alphaDrSurf[11]*fUpwind[21]+0.25*alphaDrSurf[0]*fUpwind[21]+0.25*fUpwind[7]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[5]*fUpwind[15]+0.2500000000000001*fUpwind[3]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[6]; 
  drag_incr[22] = 0.159719141249985*alphaDrSurf[20]*fUpwind[33]+0.2500000000000001*alphaDrSurf[1]*fUpwind[33]+0.223606797749979*alphaDrSurf[19]*fUpwind[32]+0.159719141249985*alphaDrSurf[12]*fUpwind[22]+0.25*alphaDrSurf[0]*fUpwind[22]+0.25*fUpwind[6]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[15]+0.2500000000000001*fUpwind[3]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[7]; 
  drag_incr[23] = 0.223606797749979*alphaDrSurf[19]*fUpwind[34]+0.2500000000000001*alphaDrSurf[2]*fUpwind[34]+0.25*alphaDrSurf[5]*fUpwind[24]+0.223606797749979*alphaDrSurf[11]*fUpwind[23]+0.25*alphaDrSurf[0]*fUpwind[23]+0.2500000000000001*alphaDrSurf[1]*fUpwind[13]; 
  drag_incr[24] = 0.223606797749979*alphaDrSurf[20]*fUpwind[34]+0.2500000000000001*alphaDrSurf[1]*fUpwind[34]+0.223606797749979*alphaDrSurf[12]*fUpwind[24]+0.25*alphaDrSurf[0]*fUpwind[24]+0.25*alphaDrSurf[5]*fUpwind[23]+0.2500000000000001*alphaDrSurf[2]*fUpwind[13]; 
  drag_incr[25] = 0.223606797749979*alphaDrSurf[20]*fUpwind[36]+0.159719141249985*alphaDrSurf[19]*fUpwind[35]+0.2500000000000001*alphaDrSurf[2]*fUpwind[35]+0.159719141249985*alphaDrSurf[11]*fUpwind[25]+0.25*alphaDrSurf[0]*fUpwind[25]+0.25*fUpwind[9]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[1]*fUpwind[8]; 
  drag_incr[26] = 0.159719141249985*alphaDrSurf[20]*fUpwind[36]+0.2500000000000001*alphaDrSurf[1]*fUpwind[36]+0.223606797749979*alphaDrSurf[19]*fUpwind[35]+0.159719141249985*alphaDrSurf[12]*fUpwind[26]+0.25*alphaDrSurf[0]*fUpwind[26]+0.25*fUpwind[8]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[16]+0.2500000000000001*fUpwind[4]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[2]*fUpwind[9]; 
  drag_incr[27] = 0.25*alphaDrSurf[5]*fUpwind[46]+0.2500000000000001*alphaDrSurf[2]*fUpwind[40]+0.2500000000000001*alphaDrSurf[1]*fUpwind[39]+0.25*alphaDrSurf[0]*fUpwind[27]; 
  drag_incr[28] = 0.223606797749979*alphaDrSurf[19]*fUpwind[41]+0.2500000000000001*alphaDrSurf[2]*fUpwind[41]+0.25*alphaDrSurf[5]*fUpwind[29]+0.223606797749979*alphaDrSurf[11]*fUpwind[28]+0.25*alphaDrSurf[0]*fUpwind[28]+0.2500000000000001*alphaDrSurf[1]*fUpwind[14]; 
  drag_incr[29] = 0.223606797749979*alphaDrSurf[20]*fUpwind[41]+0.2500000000000001*alphaDrSurf[1]*fUpwind[41]+0.223606797749979*alphaDrSurf[12]*fUpwind[29]+0.25*alphaDrSurf[0]*fUpwind[29]+0.25*alphaDrSurf[5]*fUpwind[28]+0.2500000000000001*alphaDrSurf[2]*fUpwind[14]; 
  drag_incr[30] = 0.25*alphaDrSurf[5]*fUpwind[47]+0.2500000000000001*alphaDrSurf[2]*fUpwind[43]+0.2500000000000001*alphaDrSurf[1]*fUpwind[42]+0.25*alphaDrSurf[0]*fUpwind[30]; 
  drag_incr[31] = 0.2*alphaDrSurf[19]*fUpwind[45]+0.223606797749979*alphaDrSurf[2]*fUpwind[45]+0.2*alphaDrSurf[20]*fUpwind[44]+0.223606797749979*alphaDrSurf[1]*fUpwind[44]+0.223606797749979*alphaDrSurf[5]*fUpwind[38]+0.223606797749979*alphaDrSurf[5]*fUpwind[37]+0.223606797749979*alphaDrSurf[12]*fUpwind[31]+0.223606797749979*alphaDrSurf[11]*fUpwind[31]+0.25*alphaDrSurf[0]*fUpwind[31]+0.223606797749979*fUpwind[18]*alphaDrSurf[20]+0.223606797749979*fUpwind[17]*alphaDrSurf[19]+0.25*alphaDrSurf[1]*fUpwind[18]+0.25*alphaDrSurf[2]*fUpwind[17]+0.25*alphaDrSurf[5]*fUpwind[10]; 
  drag_incr[32] = 0.2*alphaDrSurf[5]*fUpwind[33]+0.223606797749979*alphaDrSurf[12]*fUpwind[32]+0.159719141249985*alphaDrSurf[11]*fUpwind[32]+0.25*alphaDrSurf[0]*fUpwind[32]+0.223606797749979*alphaDrSurf[19]*fUpwind[22]+0.159719141249985*alphaDrSurf[19]*fUpwind[21]+0.2500000000000001*alphaDrSurf[2]*fUpwind[21]+0.2*fUpwind[15]*alphaDrSurf[20]+0.2500000000000001*fUpwind[3]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[15]+0.25*fUpwind[7]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[6]; 
  drag_incr[33] = 0.159719141249985*alphaDrSurf[12]*fUpwind[33]+0.223606797749979*alphaDrSurf[11]*fUpwind[33]+0.25*alphaDrSurf[0]*fUpwind[33]+0.2*alphaDrSurf[5]*fUpwind[32]+0.159719141249985*alphaDrSurf[20]*fUpwind[22]+0.2500000000000001*alphaDrSurf[1]*fUpwind[22]+0.223606797749979*alphaDrSurf[20]*fUpwind[21]+0.2500000000000001*fUpwind[3]*alphaDrSurf[20]+0.2*fUpwind[15]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[2]*fUpwind[15]+0.25*fUpwind[6]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[7]; 
  drag_incr[34] = 0.223606797749979*alphaDrSurf[12]*fUpwind[34]+0.223606797749979*alphaDrSurf[11]*fUpwind[34]+0.25*alphaDrSurf[0]*fUpwind[34]+0.223606797749979*alphaDrSurf[20]*fUpwind[24]+0.2500000000000001*alphaDrSurf[1]*fUpwind[24]+0.223606797749979*alphaDrSurf[19]*fUpwind[23]+0.2500000000000001*alphaDrSurf[2]*fUpwind[23]+0.25*alphaDrSurf[5]*fUpwind[13]; 
  drag_incr[35] = 0.2*alphaDrSurf[5]*fUpwind[36]+0.223606797749979*alphaDrSurf[12]*fUpwind[35]+0.159719141249985*alphaDrSurf[11]*fUpwind[35]+0.25*alphaDrSurf[0]*fUpwind[35]+0.223606797749979*alphaDrSurf[19]*fUpwind[26]+0.159719141249985*alphaDrSurf[19]*fUpwind[25]+0.2500000000000001*alphaDrSurf[2]*fUpwind[25]+0.2*fUpwind[16]*alphaDrSurf[20]+0.2500000000000001*fUpwind[4]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[16]+0.25*fUpwind[9]*alphaDrSurf[11]+0.223606797749979*alphaDrSurf[5]*fUpwind[8]; 
  drag_incr[36] = 0.159719141249985*alphaDrSurf[12]*fUpwind[36]+0.223606797749979*alphaDrSurf[11]*fUpwind[36]+0.25*alphaDrSurf[0]*fUpwind[36]+0.2*alphaDrSurf[5]*fUpwind[35]+0.159719141249985*alphaDrSurf[20]*fUpwind[26]+0.2500000000000001*alphaDrSurf[1]*fUpwind[26]+0.223606797749979*alphaDrSurf[20]*fUpwind[25]+0.2500000000000001*fUpwind[4]*alphaDrSurf[20]+0.2*fUpwind[16]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[2]*fUpwind[16]+0.25*fUpwind[8]*alphaDrSurf[12]+0.223606797749979*alphaDrSurf[5]*fUpwind[9]; 
  drag_incr[37] = 0.223606797749979*alphaDrSurf[20]*fUpwind[45]+0.159719141249985*alphaDrSurf[19]*fUpwind[44]+0.2500000000000001*alphaDrSurf[2]*fUpwind[44]+0.159719141249985*alphaDrSurf[11]*fUpwind[37]+0.25*alphaDrSurf[0]*fUpwind[37]+0.223606797749979*alphaDrSurf[5]*fUpwind[31]+0.2500000000000001*fUpwind[18]*alphaDrSurf[19]+0.223606797749979*alphaDrSurf[1]*fUpwind[17]+0.25*fUpwind[10]*alphaDrSurf[11]; 
  drag_incr[38] = 0.159719141249985*alphaDrSurf[20]*fUpwind[45]+0.2500000000000001*alphaDrSurf[1]*fUpwind[45]+0.223606797749979*alphaDrSurf[19]*fUpwind[44]+0.159719141249985*alphaDrSurf[12]*fUpwind[38]+0.25*alphaDrSurf[0]*fUpwind[38]+0.223606797749979*alphaDrSurf[5]*fUpwind[31]+0.2500000000000001*fUpwind[17]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[2]*fUpwind[18]+0.25*fUpwind[10]*alphaDrSurf[12]; 
  drag_incr[39] = 0.223606797749979*alphaDrSurf[19]*fUpwind[46]+0.2500000000000001*alphaDrSurf[2]*fUpwind[46]+0.25*alphaDrSurf[5]*fUpwind[40]+0.223606797749979*alphaDrSurf[11]*fUpwind[39]+0.25*alphaDrSurf[0]*fUpwind[39]+0.2500000000000001*alphaDrSurf[1]*fUpwind[27]; 
  drag_incr[40] = 0.223606797749979*alphaDrSurf[20]*fUpwind[46]+0.2500000000000001*alphaDrSurf[1]*fUpwind[46]+0.223606797749979*alphaDrSurf[12]*fUpwind[40]+0.25*alphaDrSurf[0]*fUpwind[40]+0.25*alphaDrSurf[5]*fUpwind[39]+0.2500000000000001*alphaDrSurf[2]*fUpwind[27]; 
  drag_incr[41] = 0.223606797749979*alphaDrSurf[12]*fUpwind[41]+0.223606797749979*alphaDrSurf[11]*fUpwind[41]+0.25*alphaDrSurf[0]*fUpwind[41]+0.223606797749979*alphaDrSurf[20]*fUpwind[29]+0.2500000000000001*alphaDrSurf[1]*fUpwind[29]+0.223606797749979*alphaDrSurf[19]*fUpwind[28]+0.2500000000000001*alphaDrSurf[2]*fUpwind[28]+0.25*alphaDrSurf[5]*fUpwind[14]; 
  drag_incr[42] = 0.223606797749979*alphaDrSurf[19]*fUpwind[47]+0.2500000000000001*alphaDrSurf[2]*fUpwind[47]+0.25*alphaDrSurf[5]*fUpwind[43]+0.223606797749979*alphaDrSurf[11]*fUpwind[42]+0.25*alphaDrSurf[0]*fUpwind[42]+0.2500000000000001*alphaDrSurf[1]*fUpwind[30]; 
  drag_incr[43] = 0.223606797749979*alphaDrSurf[20]*fUpwind[47]+0.2500000000000001*alphaDrSurf[1]*fUpwind[47]+0.223606797749979*alphaDrSurf[12]*fUpwind[43]+0.25*alphaDrSurf[0]*fUpwind[43]+0.25*alphaDrSurf[5]*fUpwind[42]+0.2500000000000001*alphaDrSurf[2]*fUpwind[30]; 
  drag_incr[44] = 0.2*alphaDrSurf[5]*fUpwind[45]+0.223606797749979*alphaDrSurf[12]*fUpwind[44]+0.159719141249985*alphaDrSurf[11]*fUpwind[44]+0.25*alphaDrSurf[0]*fUpwind[44]+0.223606797749979*alphaDrSurf[19]*fUpwind[38]+0.159719141249985*alphaDrSurf[19]*fUpwind[37]+0.2500000000000001*alphaDrSurf[2]*fUpwind[37]+0.2*alphaDrSurf[20]*fUpwind[31]+0.223606797749979*alphaDrSurf[1]*fUpwind[31]+0.25*fUpwind[10]*alphaDrSurf[19]+0.2500000000000001*alphaDrSurf[11]*fUpwind[18]+0.223606797749979*alphaDrSurf[5]*fUpwind[17]; 
  drag_incr[45] = 0.159719141249985*alphaDrSurf[12]*fUpwind[45]+0.223606797749979*alphaDrSurf[11]*fUpwind[45]+0.25*alphaDrSurf[0]*fUpwind[45]+0.2*alphaDrSurf[5]*fUpwind[44]+0.159719141249985*alphaDrSurf[20]*fUpwind[38]+0.2500000000000001*alphaDrSurf[1]*fUpwind[38]+0.223606797749979*alphaDrSurf[20]*fUpwind[37]+0.2*alphaDrSurf[19]*fUpwind[31]+0.223606797749979*alphaDrSurf[2]*fUpwind[31]+0.25*fUpwind[10]*alphaDrSurf[20]+0.223606797749979*alphaDrSurf[5]*fUpwind[18]+0.2500000000000001*alphaDrSurf[12]*fUpwind[17]; 
  drag_incr[46] = 0.223606797749979*alphaDrSurf[12]*fUpwind[46]+0.223606797749979*alphaDrSurf[11]*fUpwind[46]+0.25*alphaDrSurf[0]*fUpwind[46]+0.223606797749979*alphaDrSurf[20]*fUpwind[40]+0.2500000000000001*alphaDrSurf[1]*fUpwind[40]+0.223606797749979*alphaDrSurf[19]*fUpwind[39]+0.2500000000000001*alphaDrSurf[2]*fUpwind[39]+0.25*alphaDrSurf[5]*fUpwind[27]; 
  drag_incr[47] = 0.223606797749979*alphaDrSurf[12]*fUpwind[47]+0.223606797749979*alphaDrSurf[11]*fUpwind[47]+0.25*alphaDrSurf[0]*fUpwind[47]+0.223606797749979*alphaDrSurf[20]*fUpwind[43]+0.2500000000000001*alphaDrSurf[1]*fUpwind[43]+0.223606797749979*alphaDrSurf[19]*fUpwind[42]+0.2500000000000001*alphaDrSurf[2]*fUpwind[42]+0.25*alphaDrSurf[5]*fUpwind[30]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[11] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[12] += -0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += -0.7071067811865475*drag_incr[9]*rdv2; 
  out[14] += -0.7071067811865475*drag_incr[10]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[16] += -0.7071067811865475*drag_incr[11]*rdv2; 
  out[17] += -0.7071067811865475*drag_incr[12]*rdv2; 
  out[18] += -0.7071067811865475*drag_incr[13]*rdv2; 
  out[19] += -1.58113883008419*drag_incr[0]*rdv2; 
  out[20] += -0.7071067811865475*drag_incr[14]*rdv2; 
  out[21] += -0.7071067811865475*drag_incr[15]*rdv2; 
  out[22] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[23] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[24] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[25] += -0.7071067811865475*drag_incr[16]*rdv2; 
  out[26] += -0.7071067811865475*drag_incr[17]*rdv2; 
  out[27] += -0.7071067811865475*drag_incr[18]*rdv2; 
  out[28] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[31] += -0.7071067811865475*drag_incr[19]*rdv2; 
  out[32] += -0.7071067811865475*drag_incr[20]*rdv2; 
  out[33] += -0.7071067811865475*drag_incr[21]*rdv2; 
  out[34] += -0.7071067811865475*drag_incr[22]*rdv2; 
  out[35] += -0.7071067811865475*drag_incr[23]*rdv2; 
  out[36] += -0.7071067811865475*drag_incr[24]*rdv2; 
  out[37] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[38] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[39] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[40] += -1.58113883008419*drag_incr[1]*rdv2; 
  out[41] += -1.58113883008419*drag_incr[2]*rdv2; 
  out[42] += -1.58113883008419*drag_incr[3]*rdv2; 
  out[43] += -0.7071067811865475*drag_incr[25]*rdv2; 
  out[44] += -0.7071067811865475*drag_incr[26]*rdv2; 
  out[45] += -0.7071067811865475*drag_incr[27]*rdv2; 
  out[46] += -1.58113883008419*drag_incr[4]*rdv2; 
  out[47] += -0.7071067811865475*drag_incr[28]*rdv2; 
  out[48] += -0.7071067811865475*drag_incr[29]*rdv2; 
  out[49] += -0.7071067811865475*drag_incr[30]*rdv2; 
  out[50] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[51] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[52] += -0.7071067811865475*drag_incr[31]*rdv2; 
  out[53] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[54] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[55] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[56] += -0.7071067811865475*drag_incr[32]*rdv2; 
  out[57] += -0.7071067811865475*drag_incr[33]*rdv2; 
  out[58] += -0.7071067811865475*drag_incr[34]*rdv2; 
  out[59] += 1.224744871391589*drag_incr[19]*rdv2; 
  out[60] += 1.224744871391589*drag_incr[20]*rdv2; 
  out[61] += 1.224744871391589*drag_incr[21]*rdv2; 
  out[62] += 1.224744871391589*drag_incr[22]*rdv2; 
  out[63] += 1.224744871391589*drag_incr[23]*rdv2; 
  out[64] += 1.224744871391589*drag_incr[24]*rdv2; 
  out[65] += -1.58113883008419*drag_incr[5]*rdv2; 
  out[66] += -1.58113883008419*drag_incr[6]*rdv2; 
  out[67] += -1.58113883008419*drag_incr[7]*rdv2; 
  out[68] += -0.7071067811865475*drag_incr[35]*rdv2; 
  out[69] += -0.7071067811865475*drag_incr[36]*rdv2; 
  out[70] += -0.7071067811865475*drag_incr[37]*rdv2; 
  out[71] += -0.7071067811865475*drag_incr[38]*rdv2; 
  out[72] += -0.7071067811865475*drag_incr[39]*rdv2; 
  out[73] += -0.7071067811865475*drag_incr[40]*rdv2; 
  out[74] += 1.224744871391589*drag_incr[25]*rdv2; 
  out[75] += 1.224744871391589*drag_incr[26]*rdv2; 
  out[76] += 1.224744871391589*drag_incr[27]*rdv2; 
  out[77] += -1.58113883008419*drag_incr[8]*rdv2; 
  out[78] += -1.58113883008419*drag_incr[9]*rdv2; 
  out[79] += -1.58113883008419*drag_incr[10]*rdv2; 
  out[80] += -0.7071067811865475*drag_incr[41]*rdv2; 
  out[81] += -0.7071067811865475*drag_incr[42]*rdv2; 
  out[82] += -0.7071067811865475*drag_incr[43]*rdv2; 
  out[83] += 1.224744871391589*drag_incr[28]*rdv2; 
  out[84] += 1.224744871391589*drag_incr[29]*rdv2; 
  out[85] += 1.224744871391589*drag_incr[30]*rdv2; 
  out[86] += 1.224744871391589*drag_incr[31]*rdv2; 
  out[87] += 1.224744871391589*drag_incr[32]*rdv2; 
  out[88] += 1.224744871391589*drag_incr[33]*rdv2; 
  out[89] += 1.224744871391589*drag_incr[34]*rdv2; 
  out[90] += -1.58113883008419*drag_incr[15]*rdv2; 
  out[91] += -0.7071067811865475*drag_incr[44]*rdv2; 
  out[92] += -0.7071067811865475*drag_incr[45]*rdv2; 
  out[93] += -0.7071067811865475*drag_incr[46]*rdv2; 
  out[94] += 1.224744871391589*drag_incr[35]*rdv2; 
  out[95] += 1.224744871391589*drag_incr[36]*rdv2; 
  out[96] += 1.224744871391589*drag_incr[37]*rdv2; 
  out[97] += 1.224744871391589*drag_incr[38]*rdv2; 
  out[98] += 1.224744871391589*drag_incr[39]*rdv2; 
  out[99] += 1.224744871391589*drag_incr[40]*rdv2; 
  out[100] += -1.58113883008419*drag_incr[16]*rdv2; 
  out[101] += -1.58113883008419*drag_incr[17]*rdv2; 
  out[102] += -1.58113883008419*drag_incr[18]*rdv2; 
  out[103] += -0.7071067811865475*drag_incr[47]*rdv2; 
  out[104] += 1.224744871391589*drag_incr[41]*rdv2; 
  out[105] += 1.224744871391589*drag_incr[42]*rdv2; 
  out[106] += 1.224744871391589*drag_incr[43]*rdv2; 
  out[107] += 1.224744871391589*drag_incr[44]*rdv2; 
  out[108] += 1.224744871391589*drag_incr[45]*rdv2; 
  out[109] += 1.224744871391589*drag_incr[46]*rdv2; 
  out[110] += -1.58113883008419*drag_incr[31]*rdv2; 
  out[111] += 1.224744871391589*drag_incr[47]*rdv2; 

  } 
} 
