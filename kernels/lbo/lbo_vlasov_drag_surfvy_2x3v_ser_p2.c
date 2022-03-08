#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_5x_p2_surfx4_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvy_2x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         cell-center coordinates. 
  // dxv[5]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[24]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[8]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *sumNuUy = &nuUSum[8]; 

  double alphaDrSurf_l[48] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*w[3]-1.0*nuSum[0]*dxv[3]-2.0*sumNuUy[0]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*w[3]-1.0*nuSum[1]*dxv[3]-2.0*sumNuUy[1]; 
  alphaDrSurf_l[2] = 2.0*nuSum[2]*w[3]-1.0*nuSum[2]*dxv[3]-2.0*sumNuUy[2]; 
  alphaDrSurf_l[5] = 2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]-1.0*dxv[3]*nuSum[3]; 
  alphaDrSurf_l[11] = (-2.0*sumNuUy[4])+2.0*w[3]*nuSum[4]-1.0*dxv[3]*nuSum[4]; 
  alphaDrSurf_l[12] = (-2.0*sumNuUy[5])+2.0*w[3]*nuSum[5]-1.0*dxv[3]*nuSum[5]; 
  alphaDrSurf_l[19] = (-2.0*sumNuUy[6])+2.0*w[3]*nuSum[6]-1.0*dxv[3]*nuSum[6]; 
  alphaDrSurf_l[20] = (-2.0*sumNuUy[7])+2.0*w[3]*nuSum[7]-1.0*dxv[3]*nuSum[7]; 

  double alphaDrSurf_r[48] = {0.0}; 
  alphaDrSurf_r[0] = 2.0*nuSum[0]*w[3]+nuSum[0]*dxv[3]-2.0*sumNuUy[0]; 
  alphaDrSurf_r[1] = 2.0*nuSum[1]*w[3]+nuSum[1]*dxv[3]-2.0*sumNuUy[1]; 
  alphaDrSurf_r[2] = 2.0*nuSum[2]*w[3]+nuSum[2]*dxv[3]-2.0*sumNuUy[2]; 
  alphaDrSurf_r[5] = 2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]; 
  alphaDrSurf_r[11] = (-2.0*sumNuUy[4])+2.0*w[3]*nuSum[4]+dxv[3]*nuSum[4]; 
  alphaDrSurf_r[12] = (-2.0*sumNuUy[5])+2.0*w[3]*nuSum[5]+dxv[3]*nuSum[5]; 
  alphaDrSurf_r[19] = (-2.0*sumNuUy[6])+2.0*w[3]*nuSum[6]+dxv[3]*nuSum[6]; 
  alphaDrSurf_r[20] = (-2.0*sumNuUy[7])+2.0*w[3]*nuSum[7]+dxv[3]*nuSum[7]; 

  double fUpwindQuad_l[81] = {0.0};
  double fUpwindQuad_r[81] = {0.0};
  double fUpwind_l[48] = {0.0};
  double fUpwind_r[48] = {0.0};
  double Gdrag_l[48] = {0.0}; 
  double Gdrag_r[48] = {0.0}; 

  if ((-0.2999999999999998*alphaDrSurf_l[20])-0.2999999999999999*alphaDrSurf_l[19]+0.223606797749979*(alphaDrSurf_l[12]+alphaDrSurf_l[11])+0.45*alphaDrSurf_l[5]-0.3354101966249685*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.25*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_5x_p2_surfx4_quad_0_r(fl); 
    fUpwindQuad_l[9] = ser_5x_p2_surfx4_quad_9_r(fl); 
    fUpwindQuad_l[18] = ser_5x_p2_surfx4_quad_18_r(fl); 
    fUpwindQuad_l[27] = ser_5x_p2_surfx4_quad_27_r(fl); 
    fUpwindQuad_l[36] = ser_5x_p2_surfx4_quad_36_r(fl); 
    fUpwindQuad_l[45] = ser_5x_p2_surfx4_quad_45_r(fl); 
    fUpwindQuad_l[54] = ser_5x_p2_surfx4_quad_54_r(fl); 
    fUpwindQuad_l[63] = ser_5x_p2_surfx4_quad_63_r(fl); 
    fUpwindQuad_l[72] = ser_5x_p2_surfx4_quad_72_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_5x_p2_surfx4_quad_0_l(fc); 
    fUpwindQuad_l[9] = ser_5x_p2_surfx4_quad_9_l(fc); 
    fUpwindQuad_l[18] = ser_5x_p2_surfx4_quad_18_l(fc); 
    fUpwindQuad_l[27] = ser_5x_p2_surfx4_quad_27_l(fc); 
    fUpwindQuad_l[36] = ser_5x_p2_surfx4_quad_36_l(fc); 
    fUpwindQuad_l[45] = ser_5x_p2_surfx4_quad_45_l(fc); 
    fUpwindQuad_l[54] = ser_5x_p2_surfx4_quad_54_l(fc); 
    fUpwindQuad_l[63] = ser_5x_p2_surfx4_quad_63_l(fc); 
    fUpwindQuad_l[72] = ser_5x_p2_surfx4_quad_72_l(fc); 
  } 
  if ((-0.2999999999999998*alphaDrSurf_r[20])-0.2999999999999999*alphaDrSurf_r[19]+0.223606797749979*(alphaDrSurf_r[12]+alphaDrSurf_r[11])+0.45*alphaDrSurf_r[5]-0.3354101966249685*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.25*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_5x_p2_surfx4_quad_0_r(fc); 
    fUpwindQuad_r[9] = ser_5x_p2_surfx4_quad_9_r(fc); 
    fUpwindQuad_r[18] = ser_5x_p2_surfx4_quad_18_r(fc); 
    fUpwindQuad_r[27] = ser_5x_p2_surfx4_quad_27_r(fc); 
    fUpwindQuad_r[36] = ser_5x_p2_surfx4_quad_36_r(fc); 
    fUpwindQuad_r[45] = ser_5x_p2_surfx4_quad_45_r(fc); 
    fUpwindQuad_r[54] = ser_5x_p2_surfx4_quad_54_r(fc); 
    fUpwindQuad_r[63] = ser_5x_p2_surfx4_quad_63_r(fc); 
    fUpwindQuad_r[72] = ser_5x_p2_surfx4_quad_72_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_5x_p2_surfx4_quad_0_l(fr); 
    fUpwindQuad_r[9] = ser_5x_p2_surfx4_quad_9_l(fr); 
    fUpwindQuad_r[18] = ser_5x_p2_surfx4_quad_18_l(fr); 
    fUpwindQuad_r[27] = ser_5x_p2_surfx4_quad_27_l(fr); 
    fUpwindQuad_r[36] = ser_5x_p2_surfx4_quad_36_l(fr); 
    fUpwindQuad_r[45] = ser_5x_p2_surfx4_quad_45_l(fr); 
    fUpwindQuad_r[54] = ser_5x_p2_surfx4_quad_54_l(fr); 
    fUpwindQuad_r[63] = ser_5x_p2_surfx4_quad_63_l(fr); 
    fUpwindQuad_r[72] = ser_5x_p2_surfx4_quad_72_l(fr); 
  } 
  if (0.375*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[12]-0.2795084971874737*alphaDrSurf_l[11]-0.3354101966249685*alphaDrSurf_l[2]+0.25*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_5x_p2_surfx4_quad_1_r(fl); 
    fUpwindQuad_l[10] = ser_5x_p2_surfx4_quad_10_r(fl); 
    fUpwindQuad_l[19] = ser_5x_p2_surfx4_quad_19_r(fl); 
    fUpwindQuad_l[28] = ser_5x_p2_surfx4_quad_28_r(fl); 
    fUpwindQuad_l[37] = ser_5x_p2_surfx4_quad_37_r(fl); 
    fUpwindQuad_l[46] = ser_5x_p2_surfx4_quad_46_r(fl); 
    fUpwindQuad_l[55] = ser_5x_p2_surfx4_quad_55_r(fl); 
    fUpwindQuad_l[64] = ser_5x_p2_surfx4_quad_64_r(fl); 
    fUpwindQuad_l[73] = ser_5x_p2_surfx4_quad_73_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_5x_p2_surfx4_quad_1_l(fc); 
    fUpwindQuad_l[10] = ser_5x_p2_surfx4_quad_10_l(fc); 
    fUpwindQuad_l[19] = ser_5x_p2_surfx4_quad_19_l(fc); 
    fUpwindQuad_l[28] = ser_5x_p2_surfx4_quad_28_l(fc); 
    fUpwindQuad_l[37] = ser_5x_p2_surfx4_quad_37_l(fc); 
    fUpwindQuad_l[46] = ser_5x_p2_surfx4_quad_46_l(fc); 
    fUpwindQuad_l[55] = ser_5x_p2_surfx4_quad_55_l(fc); 
    fUpwindQuad_l[64] = ser_5x_p2_surfx4_quad_64_l(fc); 
    fUpwindQuad_l[73] = ser_5x_p2_surfx4_quad_73_l(fc); 
  } 
  if (0.375*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[12]-0.2795084971874737*alphaDrSurf_r[11]-0.3354101966249685*alphaDrSurf_r[2]+0.25*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_5x_p2_surfx4_quad_1_r(fc); 
    fUpwindQuad_r[10] = ser_5x_p2_surfx4_quad_10_r(fc); 
    fUpwindQuad_r[19] = ser_5x_p2_surfx4_quad_19_r(fc); 
    fUpwindQuad_r[28] = ser_5x_p2_surfx4_quad_28_r(fc); 
    fUpwindQuad_r[37] = ser_5x_p2_surfx4_quad_37_r(fc); 
    fUpwindQuad_r[46] = ser_5x_p2_surfx4_quad_46_r(fc); 
    fUpwindQuad_r[55] = ser_5x_p2_surfx4_quad_55_r(fc); 
    fUpwindQuad_r[64] = ser_5x_p2_surfx4_quad_64_r(fc); 
    fUpwindQuad_r[73] = ser_5x_p2_surfx4_quad_73_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_5x_p2_surfx4_quad_1_l(fr); 
    fUpwindQuad_r[10] = ser_5x_p2_surfx4_quad_10_l(fr); 
    fUpwindQuad_r[19] = ser_5x_p2_surfx4_quad_19_l(fr); 
    fUpwindQuad_r[28] = ser_5x_p2_surfx4_quad_28_l(fr); 
    fUpwindQuad_r[37] = ser_5x_p2_surfx4_quad_37_l(fr); 
    fUpwindQuad_r[46] = ser_5x_p2_surfx4_quad_46_l(fr); 
    fUpwindQuad_r[55] = ser_5x_p2_surfx4_quad_55_l(fr); 
    fUpwindQuad_r[64] = ser_5x_p2_surfx4_quad_64_l(fr); 
    fUpwindQuad_r[73] = ser_5x_p2_surfx4_quad_73_l(fr); 
  } 
  if (0.2999999999999998*alphaDrSurf_l[20]-0.2999999999999999*alphaDrSurf_l[19]+0.223606797749979*(alphaDrSurf_l[12]+alphaDrSurf_l[11])-0.45*alphaDrSurf_l[5]-0.3354101966249685*alphaDrSurf_l[2]+0.3354101966249685*alphaDrSurf_l[1]+0.25*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_5x_p2_surfx4_quad_2_r(fl); 
    fUpwindQuad_l[11] = ser_5x_p2_surfx4_quad_11_r(fl); 
    fUpwindQuad_l[20] = ser_5x_p2_surfx4_quad_20_r(fl); 
    fUpwindQuad_l[29] = ser_5x_p2_surfx4_quad_29_r(fl); 
    fUpwindQuad_l[38] = ser_5x_p2_surfx4_quad_38_r(fl); 
    fUpwindQuad_l[47] = ser_5x_p2_surfx4_quad_47_r(fl); 
    fUpwindQuad_l[56] = ser_5x_p2_surfx4_quad_56_r(fl); 
    fUpwindQuad_l[65] = ser_5x_p2_surfx4_quad_65_r(fl); 
    fUpwindQuad_l[74] = ser_5x_p2_surfx4_quad_74_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_5x_p2_surfx4_quad_2_l(fc); 
    fUpwindQuad_l[11] = ser_5x_p2_surfx4_quad_11_l(fc); 
    fUpwindQuad_l[20] = ser_5x_p2_surfx4_quad_20_l(fc); 
    fUpwindQuad_l[29] = ser_5x_p2_surfx4_quad_29_l(fc); 
    fUpwindQuad_l[38] = ser_5x_p2_surfx4_quad_38_l(fc); 
    fUpwindQuad_l[47] = ser_5x_p2_surfx4_quad_47_l(fc); 
    fUpwindQuad_l[56] = ser_5x_p2_surfx4_quad_56_l(fc); 
    fUpwindQuad_l[65] = ser_5x_p2_surfx4_quad_65_l(fc); 
    fUpwindQuad_l[74] = ser_5x_p2_surfx4_quad_74_l(fc); 
  } 
  if (0.2999999999999998*alphaDrSurf_r[20]-0.2999999999999999*alphaDrSurf_r[19]+0.223606797749979*(alphaDrSurf_r[12]+alphaDrSurf_r[11])-0.45*alphaDrSurf_r[5]-0.3354101966249685*alphaDrSurf_r[2]+0.3354101966249685*alphaDrSurf_r[1]+0.25*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_5x_p2_surfx4_quad_2_r(fc); 
    fUpwindQuad_r[11] = ser_5x_p2_surfx4_quad_11_r(fc); 
    fUpwindQuad_r[20] = ser_5x_p2_surfx4_quad_20_r(fc); 
    fUpwindQuad_r[29] = ser_5x_p2_surfx4_quad_29_r(fc); 
    fUpwindQuad_r[38] = ser_5x_p2_surfx4_quad_38_r(fc); 
    fUpwindQuad_r[47] = ser_5x_p2_surfx4_quad_47_r(fc); 
    fUpwindQuad_r[56] = ser_5x_p2_surfx4_quad_56_r(fc); 
    fUpwindQuad_r[65] = ser_5x_p2_surfx4_quad_65_r(fc); 
    fUpwindQuad_r[74] = ser_5x_p2_surfx4_quad_74_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_5x_p2_surfx4_quad_2_l(fr); 
    fUpwindQuad_r[11] = ser_5x_p2_surfx4_quad_11_l(fr); 
    fUpwindQuad_r[20] = ser_5x_p2_surfx4_quad_20_l(fr); 
    fUpwindQuad_r[29] = ser_5x_p2_surfx4_quad_29_l(fr); 
    fUpwindQuad_r[38] = ser_5x_p2_surfx4_quad_38_l(fr); 
    fUpwindQuad_r[47] = ser_5x_p2_surfx4_quad_47_l(fr); 
    fUpwindQuad_r[56] = ser_5x_p2_surfx4_quad_56_l(fr); 
    fUpwindQuad_r[65] = ser_5x_p2_surfx4_quad_65_l(fr); 
    fUpwindQuad_r[74] = ser_5x_p2_surfx4_quad_74_l(fr); 
  } 
  if (0.375*alphaDrSurf_l[20]-0.2795084971874737*alphaDrSurf_l[12]+0.223606797749979*alphaDrSurf_l[11]-0.3354101966249685*alphaDrSurf_l[1]+0.25*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_5x_p2_surfx4_quad_3_r(fl); 
    fUpwindQuad_l[12] = ser_5x_p2_surfx4_quad_12_r(fl); 
    fUpwindQuad_l[21] = ser_5x_p2_surfx4_quad_21_r(fl); 
    fUpwindQuad_l[30] = ser_5x_p2_surfx4_quad_30_r(fl); 
    fUpwindQuad_l[39] = ser_5x_p2_surfx4_quad_39_r(fl); 
    fUpwindQuad_l[48] = ser_5x_p2_surfx4_quad_48_r(fl); 
    fUpwindQuad_l[57] = ser_5x_p2_surfx4_quad_57_r(fl); 
    fUpwindQuad_l[66] = ser_5x_p2_surfx4_quad_66_r(fl); 
    fUpwindQuad_l[75] = ser_5x_p2_surfx4_quad_75_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_5x_p2_surfx4_quad_3_l(fc); 
    fUpwindQuad_l[12] = ser_5x_p2_surfx4_quad_12_l(fc); 
    fUpwindQuad_l[21] = ser_5x_p2_surfx4_quad_21_l(fc); 
    fUpwindQuad_l[30] = ser_5x_p2_surfx4_quad_30_l(fc); 
    fUpwindQuad_l[39] = ser_5x_p2_surfx4_quad_39_l(fc); 
    fUpwindQuad_l[48] = ser_5x_p2_surfx4_quad_48_l(fc); 
    fUpwindQuad_l[57] = ser_5x_p2_surfx4_quad_57_l(fc); 
    fUpwindQuad_l[66] = ser_5x_p2_surfx4_quad_66_l(fc); 
    fUpwindQuad_l[75] = ser_5x_p2_surfx4_quad_75_l(fc); 
  } 
  if (0.375*alphaDrSurf_r[20]-0.2795084971874737*alphaDrSurf_r[12]+0.223606797749979*alphaDrSurf_r[11]-0.3354101966249685*alphaDrSurf_r[1]+0.25*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_5x_p2_surfx4_quad_3_r(fc); 
    fUpwindQuad_r[12] = ser_5x_p2_surfx4_quad_12_r(fc); 
    fUpwindQuad_r[21] = ser_5x_p2_surfx4_quad_21_r(fc); 
    fUpwindQuad_r[30] = ser_5x_p2_surfx4_quad_30_r(fc); 
    fUpwindQuad_r[39] = ser_5x_p2_surfx4_quad_39_r(fc); 
    fUpwindQuad_r[48] = ser_5x_p2_surfx4_quad_48_r(fc); 
    fUpwindQuad_r[57] = ser_5x_p2_surfx4_quad_57_r(fc); 
    fUpwindQuad_r[66] = ser_5x_p2_surfx4_quad_66_r(fc); 
    fUpwindQuad_r[75] = ser_5x_p2_surfx4_quad_75_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_5x_p2_surfx4_quad_3_l(fr); 
    fUpwindQuad_r[12] = ser_5x_p2_surfx4_quad_12_l(fr); 
    fUpwindQuad_r[21] = ser_5x_p2_surfx4_quad_21_l(fr); 
    fUpwindQuad_r[30] = ser_5x_p2_surfx4_quad_30_l(fr); 
    fUpwindQuad_r[39] = ser_5x_p2_surfx4_quad_39_l(fr); 
    fUpwindQuad_r[48] = ser_5x_p2_surfx4_quad_48_l(fr); 
    fUpwindQuad_r[57] = ser_5x_p2_surfx4_quad_57_l(fr); 
    fUpwindQuad_r[66] = ser_5x_p2_surfx4_quad_66_l(fr); 
    fUpwindQuad_r[75] = ser_5x_p2_surfx4_quad_75_l(fr); 
  } 
  if (0.25*alphaDrSurf_l[0]-0.2795084971874737*(alphaDrSurf_l[12]+alphaDrSurf_l[11]) < 0) { 
    fUpwindQuad_l[4] = ser_5x_p2_surfx4_quad_4_r(fl); 
    fUpwindQuad_l[13] = ser_5x_p2_surfx4_quad_13_r(fl); 
    fUpwindQuad_l[22] = ser_5x_p2_surfx4_quad_22_r(fl); 
    fUpwindQuad_l[31] = ser_5x_p2_surfx4_quad_31_r(fl); 
    fUpwindQuad_l[40] = ser_5x_p2_surfx4_quad_40_r(fl); 
    fUpwindQuad_l[49] = ser_5x_p2_surfx4_quad_49_r(fl); 
    fUpwindQuad_l[58] = ser_5x_p2_surfx4_quad_58_r(fl); 
    fUpwindQuad_l[67] = ser_5x_p2_surfx4_quad_67_r(fl); 
    fUpwindQuad_l[76] = ser_5x_p2_surfx4_quad_76_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_5x_p2_surfx4_quad_4_l(fc); 
    fUpwindQuad_l[13] = ser_5x_p2_surfx4_quad_13_l(fc); 
    fUpwindQuad_l[22] = ser_5x_p2_surfx4_quad_22_l(fc); 
    fUpwindQuad_l[31] = ser_5x_p2_surfx4_quad_31_l(fc); 
    fUpwindQuad_l[40] = ser_5x_p2_surfx4_quad_40_l(fc); 
    fUpwindQuad_l[49] = ser_5x_p2_surfx4_quad_49_l(fc); 
    fUpwindQuad_l[58] = ser_5x_p2_surfx4_quad_58_l(fc); 
    fUpwindQuad_l[67] = ser_5x_p2_surfx4_quad_67_l(fc); 
    fUpwindQuad_l[76] = ser_5x_p2_surfx4_quad_76_l(fc); 
  } 
  if (0.25*alphaDrSurf_r[0]-0.2795084971874737*(alphaDrSurf_r[12]+alphaDrSurf_r[11]) < 0) { 
    fUpwindQuad_r[4] = ser_5x_p2_surfx4_quad_4_r(fc); 
    fUpwindQuad_r[13] = ser_5x_p2_surfx4_quad_13_r(fc); 
    fUpwindQuad_r[22] = ser_5x_p2_surfx4_quad_22_r(fc); 
    fUpwindQuad_r[31] = ser_5x_p2_surfx4_quad_31_r(fc); 
    fUpwindQuad_r[40] = ser_5x_p2_surfx4_quad_40_r(fc); 
    fUpwindQuad_r[49] = ser_5x_p2_surfx4_quad_49_r(fc); 
    fUpwindQuad_r[58] = ser_5x_p2_surfx4_quad_58_r(fc); 
    fUpwindQuad_r[67] = ser_5x_p2_surfx4_quad_67_r(fc); 
    fUpwindQuad_r[76] = ser_5x_p2_surfx4_quad_76_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_5x_p2_surfx4_quad_4_l(fr); 
    fUpwindQuad_r[13] = ser_5x_p2_surfx4_quad_13_l(fr); 
    fUpwindQuad_r[22] = ser_5x_p2_surfx4_quad_22_l(fr); 
    fUpwindQuad_r[31] = ser_5x_p2_surfx4_quad_31_l(fr); 
    fUpwindQuad_r[40] = ser_5x_p2_surfx4_quad_40_l(fr); 
    fUpwindQuad_r[49] = ser_5x_p2_surfx4_quad_49_l(fr); 
    fUpwindQuad_r[58] = ser_5x_p2_surfx4_quad_58_l(fr); 
    fUpwindQuad_r[67] = ser_5x_p2_surfx4_quad_67_l(fr); 
    fUpwindQuad_r[76] = ser_5x_p2_surfx4_quad_76_l(fr); 
  } 
  if ((-0.375*alphaDrSurf_l[20])-0.2795084971874737*alphaDrSurf_l[12]+0.223606797749979*alphaDrSurf_l[11]+0.3354101966249685*alphaDrSurf_l[1]+0.25*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_5x_p2_surfx4_quad_5_r(fl); 
    fUpwindQuad_l[14] = ser_5x_p2_surfx4_quad_14_r(fl); 
    fUpwindQuad_l[23] = ser_5x_p2_surfx4_quad_23_r(fl); 
    fUpwindQuad_l[32] = ser_5x_p2_surfx4_quad_32_r(fl); 
    fUpwindQuad_l[41] = ser_5x_p2_surfx4_quad_41_r(fl); 
    fUpwindQuad_l[50] = ser_5x_p2_surfx4_quad_50_r(fl); 
    fUpwindQuad_l[59] = ser_5x_p2_surfx4_quad_59_r(fl); 
    fUpwindQuad_l[68] = ser_5x_p2_surfx4_quad_68_r(fl); 
    fUpwindQuad_l[77] = ser_5x_p2_surfx4_quad_77_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_5x_p2_surfx4_quad_5_l(fc); 
    fUpwindQuad_l[14] = ser_5x_p2_surfx4_quad_14_l(fc); 
    fUpwindQuad_l[23] = ser_5x_p2_surfx4_quad_23_l(fc); 
    fUpwindQuad_l[32] = ser_5x_p2_surfx4_quad_32_l(fc); 
    fUpwindQuad_l[41] = ser_5x_p2_surfx4_quad_41_l(fc); 
    fUpwindQuad_l[50] = ser_5x_p2_surfx4_quad_50_l(fc); 
    fUpwindQuad_l[59] = ser_5x_p2_surfx4_quad_59_l(fc); 
    fUpwindQuad_l[68] = ser_5x_p2_surfx4_quad_68_l(fc); 
    fUpwindQuad_l[77] = ser_5x_p2_surfx4_quad_77_l(fc); 
  } 
  if ((-0.375*alphaDrSurf_r[20])-0.2795084971874737*alphaDrSurf_r[12]+0.223606797749979*alphaDrSurf_r[11]+0.3354101966249685*alphaDrSurf_r[1]+0.25*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_5x_p2_surfx4_quad_5_r(fc); 
    fUpwindQuad_r[14] = ser_5x_p2_surfx4_quad_14_r(fc); 
    fUpwindQuad_r[23] = ser_5x_p2_surfx4_quad_23_r(fc); 
    fUpwindQuad_r[32] = ser_5x_p2_surfx4_quad_32_r(fc); 
    fUpwindQuad_r[41] = ser_5x_p2_surfx4_quad_41_r(fc); 
    fUpwindQuad_r[50] = ser_5x_p2_surfx4_quad_50_r(fc); 
    fUpwindQuad_r[59] = ser_5x_p2_surfx4_quad_59_r(fc); 
    fUpwindQuad_r[68] = ser_5x_p2_surfx4_quad_68_r(fc); 
    fUpwindQuad_r[77] = ser_5x_p2_surfx4_quad_77_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_5x_p2_surfx4_quad_5_l(fr); 
    fUpwindQuad_r[14] = ser_5x_p2_surfx4_quad_14_l(fr); 
    fUpwindQuad_r[23] = ser_5x_p2_surfx4_quad_23_l(fr); 
    fUpwindQuad_r[32] = ser_5x_p2_surfx4_quad_32_l(fr); 
    fUpwindQuad_r[41] = ser_5x_p2_surfx4_quad_41_l(fr); 
    fUpwindQuad_r[50] = ser_5x_p2_surfx4_quad_50_l(fr); 
    fUpwindQuad_r[59] = ser_5x_p2_surfx4_quad_59_l(fr); 
    fUpwindQuad_r[68] = ser_5x_p2_surfx4_quad_68_l(fr); 
    fUpwindQuad_r[77] = ser_5x_p2_surfx4_quad_77_l(fr); 
  } 
  if ((-0.2999999999999998*alphaDrSurf_l[20])+0.2999999999999999*alphaDrSurf_l[19]+0.223606797749979*(alphaDrSurf_l[12]+alphaDrSurf_l[11])-0.45*alphaDrSurf_l[5]+0.3354101966249685*alphaDrSurf_l[2]-0.3354101966249685*alphaDrSurf_l[1]+0.25*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_5x_p2_surfx4_quad_6_r(fl); 
    fUpwindQuad_l[15] = ser_5x_p2_surfx4_quad_15_r(fl); 
    fUpwindQuad_l[24] = ser_5x_p2_surfx4_quad_24_r(fl); 
    fUpwindQuad_l[33] = ser_5x_p2_surfx4_quad_33_r(fl); 
    fUpwindQuad_l[42] = ser_5x_p2_surfx4_quad_42_r(fl); 
    fUpwindQuad_l[51] = ser_5x_p2_surfx4_quad_51_r(fl); 
    fUpwindQuad_l[60] = ser_5x_p2_surfx4_quad_60_r(fl); 
    fUpwindQuad_l[69] = ser_5x_p2_surfx4_quad_69_r(fl); 
    fUpwindQuad_l[78] = ser_5x_p2_surfx4_quad_78_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_5x_p2_surfx4_quad_6_l(fc); 
    fUpwindQuad_l[15] = ser_5x_p2_surfx4_quad_15_l(fc); 
    fUpwindQuad_l[24] = ser_5x_p2_surfx4_quad_24_l(fc); 
    fUpwindQuad_l[33] = ser_5x_p2_surfx4_quad_33_l(fc); 
    fUpwindQuad_l[42] = ser_5x_p2_surfx4_quad_42_l(fc); 
    fUpwindQuad_l[51] = ser_5x_p2_surfx4_quad_51_l(fc); 
    fUpwindQuad_l[60] = ser_5x_p2_surfx4_quad_60_l(fc); 
    fUpwindQuad_l[69] = ser_5x_p2_surfx4_quad_69_l(fc); 
    fUpwindQuad_l[78] = ser_5x_p2_surfx4_quad_78_l(fc); 
  } 
  if ((-0.2999999999999998*alphaDrSurf_r[20])+0.2999999999999999*alphaDrSurf_r[19]+0.223606797749979*(alphaDrSurf_r[12]+alphaDrSurf_r[11])-0.45*alphaDrSurf_r[5]+0.3354101966249685*alphaDrSurf_r[2]-0.3354101966249685*alphaDrSurf_r[1]+0.25*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_5x_p2_surfx4_quad_6_r(fc); 
    fUpwindQuad_r[15] = ser_5x_p2_surfx4_quad_15_r(fc); 
    fUpwindQuad_r[24] = ser_5x_p2_surfx4_quad_24_r(fc); 
    fUpwindQuad_r[33] = ser_5x_p2_surfx4_quad_33_r(fc); 
    fUpwindQuad_r[42] = ser_5x_p2_surfx4_quad_42_r(fc); 
    fUpwindQuad_r[51] = ser_5x_p2_surfx4_quad_51_r(fc); 
    fUpwindQuad_r[60] = ser_5x_p2_surfx4_quad_60_r(fc); 
    fUpwindQuad_r[69] = ser_5x_p2_surfx4_quad_69_r(fc); 
    fUpwindQuad_r[78] = ser_5x_p2_surfx4_quad_78_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_5x_p2_surfx4_quad_6_l(fr); 
    fUpwindQuad_r[15] = ser_5x_p2_surfx4_quad_15_l(fr); 
    fUpwindQuad_r[24] = ser_5x_p2_surfx4_quad_24_l(fr); 
    fUpwindQuad_r[33] = ser_5x_p2_surfx4_quad_33_l(fr); 
    fUpwindQuad_r[42] = ser_5x_p2_surfx4_quad_42_l(fr); 
    fUpwindQuad_r[51] = ser_5x_p2_surfx4_quad_51_l(fr); 
    fUpwindQuad_r[60] = ser_5x_p2_surfx4_quad_60_l(fr); 
    fUpwindQuad_r[69] = ser_5x_p2_surfx4_quad_69_l(fr); 
    fUpwindQuad_r[78] = ser_5x_p2_surfx4_quad_78_l(fr); 
  } 
  if ((-0.375*alphaDrSurf_l[19])+0.223606797749979*alphaDrSurf_l[12]-0.2795084971874737*alphaDrSurf_l[11]+0.3354101966249685*alphaDrSurf_l[2]+0.25*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[7] = ser_5x_p2_surfx4_quad_7_r(fl); 
    fUpwindQuad_l[16] = ser_5x_p2_surfx4_quad_16_r(fl); 
    fUpwindQuad_l[25] = ser_5x_p2_surfx4_quad_25_r(fl); 
    fUpwindQuad_l[34] = ser_5x_p2_surfx4_quad_34_r(fl); 
    fUpwindQuad_l[43] = ser_5x_p2_surfx4_quad_43_r(fl); 
    fUpwindQuad_l[52] = ser_5x_p2_surfx4_quad_52_r(fl); 
    fUpwindQuad_l[61] = ser_5x_p2_surfx4_quad_61_r(fl); 
    fUpwindQuad_l[70] = ser_5x_p2_surfx4_quad_70_r(fl); 
    fUpwindQuad_l[79] = ser_5x_p2_surfx4_quad_79_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_5x_p2_surfx4_quad_7_l(fc); 
    fUpwindQuad_l[16] = ser_5x_p2_surfx4_quad_16_l(fc); 
    fUpwindQuad_l[25] = ser_5x_p2_surfx4_quad_25_l(fc); 
    fUpwindQuad_l[34] = ser_5x_p2_surfx4_quad_34_l(fc); 
    fUpwindQuad_l[43] = ser_5x_p2_surfx4_quad_43_l(fc); 
    fUpwindQuad_l[52] = ser_5x_p2_surfx4_quad_52_l(fc); 
    fUpwindQuad_l[61] = ser_5x_p2_surfx4_quad_61_l(fc); 
    fUpwindQuad_l[70] = ser_5x_p2_surfx4_quad_70_l(fc); 
    fUpwindQuad_l[79] = ser_5x_p2_surfx4_quad_79_l(fc); 
  } 
  if ((-0.375*alphaDrSurf_r[19])+0.223606797749979*alphaDrSurf_r[12]-0.2795084971874737*alphaDrSurf_r[11]+0.3354101966249685*alphaDrSurf_r[2]+0.25*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[7] = ser_5x_p2_surfx4_quad_7_r(fc); 
    fUpwindQuad_r[16] = ser_5x_p2_surfx4_quad_16_r(fc); 
    fUpwindQuad_r[25] = ser_5x_p2_surfx4_quad_25_r(fc); 
    fUpwindQuad_r[34] = ser_5x_p2_surfx4_quad_34_r(fc); 
    fUpwindQuad_r[43] = ser_5x_p2_surfx4_quad_43_r(fc); 
    fUpwindQuad_r[52] = ser_5x_p2_surfx4_quad_52_r(fc); 
    fUpwindQuad_r[61] = ser_5x_p2_surfx4_quad_61_r(fc); 
    fUpwindQuad_r[70] = ser_5x_p2_surfx4_quad_70_r(fc); 
    fUpwindQuad_r[79] = ser_5x_p2_surfx4_quad_79_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_5x_p2_surfx4_quad_7_l(fr); 
    fUpwindQuad_r[16] = ser_5x_p2_surfx4_quad_16_l(fr); 
    fUpwindQuad_r[25] = ser_5x_p2_surfx4_quad_25_l(fr); 
    fUpwindQuad_r[34] = ser_5x_p2_surfx4_quad_34_l(fr); 
    fUpwindQuad_r[43] = ser_5x_p2_surfx4_quad_43_l(fr); 
    fUpwindQuad_r[52] = ser_5x_p2_surfx4_quad_52_l(fr); 
    fUpwindQuad_r[61] = ser_5x_p2_surfx4_quad_61_l(fr); 
    fUpwindQuad_r[70] = ser_5x_p2_surfx4_quad_70_l(fr); 
    fUpwindQuad_r[79] = ser_5x_p2_surfx4_quad_79_l(fr); 
  } 
  if (0.2999999999999998*alphaDrSurf_l[20]+0.2999999999999999*alphaDrSurf_l[19]+0.223606797749979*(alphaDrSurf_l[12]+alphaDrSurf_l[11])+0.45*alphaDrSurf_l[5]+0.3354101966249685*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.25*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = ser_5x_p2_surfx4_quad_8_r(fl); 
    fUpwindQuad_l[17] = ser_5x_p2_surfx4_quad_17_r(fl); 
    fUpwindQuad_l[26] = ser_5x_p2_surfx4_quad_26_r(fl); 
    fUpwindQuad_l[35] = ser_5x_p2_surfx4_quad_35_r(fl); 
    fUpwindQuad_l[44] = ser_5x_p2_surfx4_quad_44_r(fl); 
    fUpwindQuad_l[53] = ser_5x_p2_surfx4_quad_53_r(fl); 
    fUpwindQuad_l[62] = ser_5x_p2_surfx4_quad_62_r(fl); 
    fUpwindQuad_l[71] = ser_5x_p2_surfx4_quad_71_r(fl); 
    fUpwindQuad_l[80] = ser_5x_p2_surfx4_quad_80_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_5x_p2_surfx4_quad_8_l(fc); 
    fUpwindQuad_l[17] = ser_5x_p2_surfx4_quad_17_l(fc); 
    fUpwindQuad_l[26] = ser_5x_p2_surfx4_quad_26_l(fc); 
    fUpwindQuad_l[35] = ser_5x_p2_surfx4_quad_35_l(fc); 
    fUpwindQuad_l[44] = ser_5x_p2_surfx4_quad_44_l(fc); 
    fUpwindQuad_l[53] = ser_5x_p2_surfx4_quad_53_l(fc); 
    fUpwindQuad_l[62] = ser_5x_p2_surfx4_quad_62_l(fc); 
    fUpwindQuad_l[71] = ser_5x_p2_surfx4_quad_71_l(fc); 
    fUpwindQuad_l[80] = ser_5x_p2_surfx4_quad_80_l(fc); 
  } 
  if (0.2999999999999998*alphaDrSurf_r[20]+0.2999999999999999*alphaDrSurf_r[19]+0.223606797749979*(alphaDrSurf_r[12]+alphaDrSurf_r[11])+0.45*alphaDrSurf_r[5]+0.3354101966249685*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.25*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = ser_5x_p2_surfx4_quad_8_r(fc); 
    fUpwindQuad_r[17] = ser_5x_p2_surfx4_quad_17_r(fc); 
    fUpwindQuad_r[26] = ser_5x_p2_surfx4_quad_26_r(fc); 
    fUpwindQuad_r[35] = ser_5x_p2_surfx4_quad_35_r(fc); 
    fUpwindQuad_r[44] = ser_5x_p2_surfx4_quad_44_r(fc); 
    fUpwindQuad_r[53] = ser_5x_p2_surfx4_quad_53_r(fc); 
    fUpwindQuad_r[62] = ser_5x_p2_surfx4_quad_62_r(fc); 
    fUpwindQuad_r[71] = ser_5x_p2_surfx4_quad_71_r(fc); 
    fUpwindQuad_r[80] = ser_5x_p2_surfx4_quad_80_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_5x_p2_surfx4_quad_8_l(fr); 
    fUpwindQuad_r[17] = ser_5x_p2_surfx4_quad_17_l(fr); 
    fUpwindQuad_r[26] = ser_5x_p2_surfx4_quad_26_l(fr); 
    fUpwindQuad_r[35] = ser_5x_p2_surfx4_quad_35_l(fr); 
    fUpwindQuad_r[44] = ser_5x_p2_surfx4_quad_44_l(fr); 
    fUpwindQuad_r[53] = ser_5x_p2_surfx4_quad_53_l(fr); 
    fUpwindQuad_r[62] = ser_5x_p2_surfx4_quad_62_l(fr); 
    fUpwindQuad_r[71] = ser_5x_p2_surfx4_quad_71_l(fr); 
    fUpwindQuad_r[80] = ser_5x_p2_surfx4_quad_80_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_5x_p2_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_5x_p2_upwind(fUpwindQuad_r, fUpwind_r); 

  Gdrag_l[0] = 0.25*alphaDrSurf_l[20]*fUpwind_l[20]+0.25*alphaDrSurf_l[19]*fUpwind_l[19]+0.25*alphaDrSurf_l[12]*fUpwind_l[12]+0.25*alphaDrSurf_l[11]*fUpwind_l[11]+0.25*alphaDrSurf_l[5]*fUpwind_l[5]+0.25*alphaDrSurf_l[2]*fUpwind_l[2]+0.25*alphaDrSurf_l[1]*fUpwind_l[1]+0.25*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.2500000000000001*alphaDrSurf_l[12]*fUpwind_l[20]+0.2500000000000001*fUpwind_l[12]*alphaDrSurf_l[20]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[19]+0.223606797749979*fUpwind_l[5]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[11]+0.223606797749979*fUpwind_l[1]*alphaDrSurf_l[11]+0.25*alphaDrSurf_l[2]*fUpwind_l[5]+0.25*fUpwind_l[2]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[0]*fUpwind_l[1]+0.25*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[20]+0.223606797749979*fUpwind_l[5]*alphaDrSurf_l[20]+0.2500000000000001*alphaDrSurf_l[11]*fUpwind_l[19]+0.2500000000000001*fUpwind_l[11]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[12]+0.223606797749979*fUpwind_l[2]*alphaDrSurf_l[12]+0.25*alphaDrSurf_l[1]*fUpwind_l[5]+0.25*fUpwind_l[1]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[0]*fUpwind_l[2]+0.25*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Gdrag_l[3] = 0.2500000000000001*alphaDrSurf_l[20]*fUpwind_l[33]+0.2500000000000001*alphaDrSurf_l[19]*fUpwind_l[32]+0.2500000000000001*alphaDrSurf_l[12]*fUpwind_l[22]+0.2500000000000001*alphaDrSurf_l[11]*fUpwind_l[21]+0.25*alphaDrSurf_l[5]*fUpwind_l[15]+0.25*alphaDrSurf_l[2]*fUpwind_l[7]+0.25*alphaDrSurf_l[1]*fUpwind_l[6]+0.25*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Gdrag_l[4] = 0.2500000000000001*alphaDrSurf_l[20]*fUpwind_l[36]+0.2500000000000001*alphaDrSurf_l[19]*fUpwind_l[35]+0.2500000000000001*alphaDrSurf_l[12]*fUpwind_l[26]+0.2500000000000001*alphaDrSurf_l[11]*fUpwind_l[25]+0.25*alphaDrSurf_l[5]*fUpwind_l[16]+0.25*alphaDrSurf_l[2]*fUpwind_l[9]+0.25*alphaDrSurf_l[1]*fUpwind_l[8]+0.25*alphaDrSurf_l[0]*fUpwind_l[4]; 
  Gdrag_l[5] = 0.2*alphaDrSurf_l[19]*fUpwind_l[20]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[20]+0.2*fUpwind_l[19]*alphaDrSurf_l[20]+0.223606797749979*fUpwind_l[2]*alphaDrSurf_l[20]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[19]+0.223606797749979*fUpwind_l[1]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[12]+0.223606797749979*fUpwind_l[5]*alphaDrSurf_l[12]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[11]+0.223606797749979*fUpwind_l[5]*alphaDrSurf_l[11]+0.25*alphaDrSurf_l[0]*fUpwind_l[5]+0.25*fUpwind_l[0]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[1]*fUpwind_l[2]+0.25*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Gdrag_l[6] = 0.25*alphaDrSurf_l[12]*fUpwind_l[33]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[32]+0.25*alphaDrSurf_l[20]*fUpwind_l[22]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[21]+0.223606797749979*fUpwind_l[15]*alphaDrSurf_l[19]+0.25*alphaDrSurf_l[2]*fUpwind_l[15]+0.223606797749979*fUpwind_l[6]*alphaDrSurf_l[11]+0.25*alphaDrSurf_l[5]*fUpwind_l[7]+0.25*alphaDrSurf_l[0]*fUpwind_l[6]+0.25*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Gdrag_l[7] = 0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[33]+0.25*alphaDrSurf_l[11]*fUpwind_l[32]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[22]+0.25*alphaDrSurf_l[19]*fUpwind_l[21]+0.223606797749979*fUpwind_l[15]*alphaDrSurf_l[20]+0.25*alphaDrSurf_l[1]*fUpwind_l[15]+0.223606797749979*fUpwind_l[7]*alphaDrSurf_l[12]+0.25*alphaDrSurf_l[0]*fUpwind_l[7]+0.25*alphaDrSurf_l[5]*fUpwind_l[6]+0.25*alphaDrSurf_l[2]*fUpwind_l[3]; 
  Gdrag_l[8] = 0.25*alphaDrSurf_l[12]*fUpwind_l[36]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[35]+0.25*alphaDrSurf_l[20]*fUpwind_l[26]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[25]+0.223606797749979*fUpwind_l[16]*alphaDrSurf_l[19]+0.25*alphaDrSurf_l[2]*fUpwind_l[16]+0.223606797749979*fUpwind_l[8]*alphaDrSurf_l[11]+0.25*alphaDrSurf_l[5]*fUpwind_l[9]+0.25*alphaDrSurf_l[0]*fUpwind_l[8]+0.25*alphaDrSurf_l[1]*fUpwind_l[4]; 
  Gdrag_l[9] = 0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[36]+0.25*alphaDrSurf_l[11]*fUpwind_l[35]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[26]+0.25*alphaDrSurf_l[19]*fUpwind_l[25]+0.223606797749979*fUpwind_l[16]*alphaDrSurf_l[20]+0.25*alphaDrSurf_l[1]*fUpwind_l[16]+0.223606797749979*fUpwind_l[9]*alphaDrSurf_l[12]+0.25*alphaDrSurf_l[0]*fUpwind_l[9]+0.25*alphaDrSurf_l[5]*fUpwind_l[8]+0.25*alphaDrSurf_l[2]*fUpwind_l[4]; 
  Gdrag_l[10] = 0.25*alphaDrSurf_l[20]*fUpwind_l[45]+0.25*alphaDrSurf_l[19]*fUpwind_l[44]+0.25*alphaDrSurf_l[12]*fUpwind_l[38]+0.25*alphaDrSurf_l[11]*fUpwind_l[37]+0.25*alphaDrSurf_l[5]*fUpwind_l[31]+0.25*alphaDrSurf_l[2]*fUpwind_l[18]+0.25*alphaDrSurf_l[1]*fUpwind_l[17]+0.25*alphaDrSurf_l[0]*fUpwind_l[10]; 
  Gdrag_l[11] = 0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[20]+0.159719141249985*alphaDrSurf_l[19]*fUpwind_l[19]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[19]+0.2500000000000001*fUpwind_l[2]*alphaDrSurf_l[19]+0.159719141249985*alphaDrSurf_l[11]*fUpwind_l[11]+0.25*alphaDrSurf_l[0]*fUpwind_l[11]+0.25*fUpwind_l[0]*alphaDrSurf_l[11]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[5]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Gdrag_l[12] = 0.159719141249985*alphaDrSurf_l[20]*fUpwind_l[20]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[20]+0.2500000000000001*fUpwind_l[1]*alphaDrSurf_l[20]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[19]+0.159719141249985*alphaDrSurf_l[12]*fUpwind_l[12]+0.25*alphaDrSurf_l[0]*fUpwind_l[12]+0.25*fUpwind_l[0]*alphaDrSurf_l[12]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[5]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[2]; 
  Gdrag_l[13] = 0.25*alphaDrSurf_l[5]*fUpwind_l[34]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[24]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[23]+0.25*alphaDrSurf_l[0]*fUpwind_l[13]; 
  Gdrag_l[14] = 0.25*alphaDrSurf_l[5]*fUpwind_l[41]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[29]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[28]+0.25*alphaDrSurf_l[0]*fUpwind_l[14]; 
  Gdrag_l[15] = 0.2*alphaDrSurf_l[19]*fUpwind_l[33]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[33]+0.2*alphaDrSurf_l[20]*fUpwind_l[32]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[32]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[22]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[21]+0.223606797749979*fUpwind_l[7]*alphaDrSurf_l[20]+0.223606797749979*fUpwind_l[6]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[15]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[15]+0.25*alphaDrSurf_l[0]*fUpwind_l[15]+0.25*alphaDrSurf_l[1]*fUpwind_l[7]+0.25*alphaDrSurf_l[2]*fUpwind_l[6]+0.25*fUpwind_l[3]*alphaDrSurf_l[5]; 
  Gdrag_l[16] = 0.2*alphaDrSurf_l[19]*fUpwind_l[36]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[36]+0.2*alphaDrSurf_l[20]*fUpwind_l[35]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[35]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[26]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[25]+0.223606797749979*fUpwind_l[9]*alphaDrSurf_l[20]+0.223606797749979*fUpwind_l[8]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[16]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[16]+0.25*alphaDrSurf_l[0]*fUpwind_l[16]+0.25*alphaDrSurf_l[1]*fUpwind_l[9]+0.25*alphaDrSurf_l[2]*fUpwind_l[8]+0.25*fUpwind_l[4]*alphaDrSurf_l[5]; 
  Gdrag_l[17] = 0.2500000000000001*alphaDrSurf_l[12]*fUpwind_l[45]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[44]+0.2500000000000001*alphaDrSurf_l[20]*fUpwind_l[38]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[37]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[31]+0.25*alphaDrSurf_l[2]*fUpwind_l[31]+0.25*alphaDrSurf_l[5]*fUpwind_l[18]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[17]+0.25*alphaDrSurf_l[0]*fUpwind_l[17]+0.25*alphaDrSurf_l[1]*fUpwind_l[10]; 
  Gdrag_l[18] = 0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[45]+0.2500000000000001*alphaDrSurf_l[11]*fUpwind_l[44]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[38]+0.2500000000000001*alphaDrSurf_l[19]*fUpwind_l[37]+0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[31]+0.25*alphaDrSurf_l[1]*fUpwind_l[31]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[18]+0.25*alphaDrSurf_l[0]*fUpwind_l[18]+0.25*alphaDrSurf_l[5]*fUpwind_l[17]+0.25*alphaDrSurf_l[2]*fUpwind_l[10]; 
  Gdrag_l[19] = 0.2*alphaDrSurf_l[5]*fUpwind_l[20]+0.2*fUpwind_l[5]*alphaDrSurf_l[20]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[19]+0.159719141249985*alphaDrSurf_l[11]*fUpwind_l[19]+0.25*alphaDrSurf_l[0]*fUpwind_l[19]+0.223606797749979*fUpwind_l[12]*alphaDrSurf_l[19]+0.159719141249985*fUpwind_l[11]*alphaDrSurf_l[19]+0.25*fUpwind_l[0]*alphaDrSurf_l[19]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[11]+0.2500000000000001*fUpwind_l[2]*alphaDrSurf_l[11]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[5]+0.223606797749979*fUpwind_l[1]*alphaDrSurf_l[5]; 
  Gdrag_l[20] = 0.159719141249985*alphaDrSurf_l[12]*fUpwind_l[20]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[20]+0.25*alphaDrSurf_l[0]*fUpwind_l[20]+0.159719141249985*fUpwind_l[12]*alphaDrSurf_l[20]+0.223606797749979*fUpwind_l[11]*alphaDrSurf_l[20]+0.25*fUpwind_l[0]*alphaDrSurf_l[20]+0.2*alphaDrSurf_l[5]*fUpwind_l[19]+0.2*fUpwind_l[5]*alphaDrSurf_l[19]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[12]+0.2500000000000001*fUpwind_l[1]*alphaDrSurf_l[12]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[5]+0.223606797749979*fUpwind_l[2]*alphaDrSurf_l[5]; 
  Gdrag_l[21] = 0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[33]+0.159719141249985*alphaDrSurf_l[19]*fUpwind_l[32]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[32]+0.159719141249985*alphaDrSurf_l[11]*fUpwind_l[21]+0.25*alphaDrSurf_l[0]*fUpwind_l[21]+0.25*fUpwind_l[7]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[15]+0.2500000000000001*fUpwind_l[3]*alphaDrSurf_l[11]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[6]; 
  Gdrag_l[22] = 0.159719141249985*alphaDrSurf_l[20]*fUpwind_l[33]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[33]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[32]+0.159719141249985*alphaDrSurf_l[12]*fUpwind_l[22]+0.25*alphaDrSurf_l[0]*fUpwind_l[22]+0.25*fUpwind_l[6]*alphaDrSurf_l[20]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[15]+0.2500000000000001*fUpwind_l[3]*alphaDrSurf_l[12]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[7]; 
  Gdrag_l[23] = 0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[34]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[34]+0.25*alphaDrSurf_l[5]*fUpwind_l[24]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[23]+0.25*alphaDrSurf_l[0]*fUpwind_l[23]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[13]; 
  Gdrag_l[24] = 0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[34]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[34]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[24]+0.25*alphaDrSurf_l[0]*fUpwind_l[24]+0.25*alphaDrSurf_l[5]*fUpwind_l[23]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[13]; 
  Gdrag_l[25] = 0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[36]+0.159719141249985*alphaDrSurf_l[19]*fUpwind_l[35]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[35]+0.159719141249985*alphaDrSurf_l[11]*fUpwind_l[25]+0.25*alphaDrSurf_l[0]*fUpwind_l[25]+0.25*fUpwind_l[9]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[16]+0.2500000000000001*fUpwind_l[4]*alphaDrSurf_l[11]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[8]; 
  Gdrag_l[26] = 0.159719141249985*alphaDrSurf_l[20]*fUpwind_l[36]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[36]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[35]+0.159719141249985*alphaDrSurf_l[12]*fUpwind_l[26]+0.25*alphaDrSurf_l[0]*fUpwind_l[26]+0.25*fUpwind_l[8]*alphaDrSurf_l[20]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[16]+0.2500000000000001*fUpwind_l[4]*alphaDrSurf_l[12]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[9]; 
  Gdrag_l[27] = 0.25*alphaDrSurf_l[5]*fUpwind_l[46]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[40]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[39]+0.25*alphaDrSurf_l[0]*fUpwind_l[27]; 
  Gdrag_l[28] = 0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[41]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[41]+0.25*alphaDrSurf_l[5]*fUpwind_l[29]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[28]+0.25*alphaDrSurf_l[0]*fUpwind_l[28]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[14]; 
  Gdrag_l[29] = 0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[41]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[41]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[29]+0.25*alphaDrSurf_l[0]*fUpwind_l[29]+0.25*alphaDrSurf_l[5]*fUpwind_l[28]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[14]; 
  Gdrag_l[30] = 0.25*alphaDrSurf_l[5]*fUpwind_l[47]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[43]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[42]+0.25*alphaDrSurf_l[0]*fUpwind_l[30]; 
  Gdrag_l[31] = 0.2*alphaDrSurf_l[19]*fUpwind_l[45]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[45]+0.2*alphaDrSurf_l[20]*fUpwind_l[44]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[44]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[38]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[37]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[31]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[31]+0.25*alphaDrSurf_l[0]*fUpwind_l[31]+0.223606797749979*fUpwind_l[18]*alphaDrSurf_l[20]+0.223606797749979*fUpwind_l[17]*alphaDrSurf_l[19]+0.25*alphaDrSurf_l[1]*fUpwind_l[18]+0.25*alphaDrSurf_l[2]*fUpwind_l[17]+0.25*alphaDrSurf_l[5]*fUpwind_l[10]; 
  Gdrag_l[32] = 0.2*alphaDrSurf_l[5]*fUpwind_l[33]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[32]+0.159719141249985*alphaDrSurf_l[11]*fUpwind_l[32]+0.25*alphaDrSurf_l[0]*fUpwind_l[32]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[22]+0.159719141249985*alphaDrSurf_l[19]*fUpwind_l[21]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[21]+0.2*fUpwind_l[15]*alphaDrSurf_l[20]+0.2500000000000001*fUpwind_l[3]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[15]+0.25*fUpwind_l[7]*alphaDrSurf_l[11]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[6]; 
  Gdrag_l[33] = 0.159719141249985*alphaDrSurf_l[12]*fUpwind_l[33]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[33]+0.25*alphaDrSurf_l[0]*fUpwind_l[33]+0.2*alphaDrSurf_l[5]*fUpwind_l[32]+0.159719141249985*alphaDrSurf_l[20]*fUpwind_l[22]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[22]+0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[21]+0.2500000000000001*fUpwind_l[3]*alphaDrSurf_l[20]+0.2*fUpwind_l[15]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[15]+0.25*fUpwind_l[6]*alphaDrSurf_l[12]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[7]; 
  Gdrag_l[34] = 0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[34]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[34]+0.25*alphaDrSurf_l[0]*fUpwind_l[34]+0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[24]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[24]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[23]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[23]+0.25*alphaDrSurf_l[5]*fUpwind_l[13]; 
  Gdrag_l[35] = 0.2*alphaDrSurf_l[5]*fUpwind_l[36]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[35]+0.159719141249985*alphaDrSurf_l[11]*fUpwind_l[35]+0.25*alphaDrSurf_l[0]*fUpwind_l[35]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[26]+0.159719141249985*alphaDrSurf_l[19]*fUpwind_l[25]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[25]+0.2*fUpwind_l[16]*alphaDrSurf_l[20]+0.2500000000000001*fUpwind_l[4]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[16]+0.25*fUpwind_l[9]*alphaDrSurf_l[11]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[8]; 
  Gdrag_l[36] = 0.159719141249985*alphaDrSurf_l[12]*fUpwind_l[36]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[36]+0.25*alphaDrSurf_l[0]*fUpwind_l[36]+0.2*alphaDrSurf_l[5]*fUpwind_l[35]+0.159719141249985*alphaDrSurf_l[20]*fUpwind_l[26]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[26]+0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[25]+0.2500000000000001*fUpwind_l[4]*alphaDrSurf_l[20]+0.2*fUpwind_l[16]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[16]+0.25*fUpwind_l[8]*alphaDrSurf_l[12]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[9]; 
  Gdrag_l[37] = 0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[45]+0.159719141249985*alphaDrSurf_l[19]*fUpwind_l[44]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[44]+0.159719141249985*alphaDrSurf_l[11]*fUpwind_l[37]+0.25*alphaDrSurf_l[0]*fUpwind_l[37]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[31]+0.2500000000000001*fUpwind_l[18]*alphaDrSurf_l[19]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[17]+0.25*fUpwind_l[10]*alphaDrSurf_l[11]; 
  Gdrag_l[38] = 0.159719141249985*alphaDrSurf_l[20]*fUpwind_l[45]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[45]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[44]+0.159719141249985*alphaDrSurf_l[12]*fUpwind_l[38]+0.25*alphaDrSurf_l[0]*fUpwind_l[38]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[31]+0.2500000000000001*fUpwind_l[17]*alphaDrSurf_l[20]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[18]+0.25*fUpwind_l[10]*alphaDrSurf_l[12]; 
  Gdrag_l[39] = 0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[46]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[46]+0.25*alphaDrSurf_l[5]*fUpwind_l[40]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[39]+0.25*alphaDrSurf_l[0]*fUpwind_l[39]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[27]; 
  Gdrag_l[40] = 0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[46]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[46]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[40]+0.25*alphaDrSurf_l[0]*fUpwind_l[40]+0.25*alphaDrSurf_l[5]*fUpwind_l[39]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[27]; 
  Gdrag_l[41] = 0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[41]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[41]+0.25*alphaDrSurf_l[0]*fUpwind_l[41]+0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[29]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[29]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[28]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[28]+0.25*alphaDrSurf_l[5]*fUpwind_l[14]; 
  Gdrag_l[42] = 0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[47]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[47]+0.25*alphaDrSurf_l[5]*fUpwind_l[43]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[42]+0.25*alphaDrSurf_l[0]*fUpwind_l[42]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[30]; 
  Gdrag_l[43] = 0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[47]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[47]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[43]+0.25*alphaDrSurf_l[0]*fUpwind_l[43]+0.25*alphaDrSurf_l[5]*fUpwind_l[42]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[30]; 
  Gdrag_l[44] = 0.2*alphaDrSurf_l[5]*fUpwind_l[45]+0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[44]+0.159719141249985*alphaDrSurf_l[11]*fUpwind_l[44]+0.25*alphaDrSurf_l[0]*fUpwind_l[44]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[38]+0.159719141249985*alphaDrSurf_l[19]*fUpwind_l[37]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[37]+0.2*alphaDrSurf_l[20]*fUpwind_l[31]+0.223606797749979*alphaDrSurf_l[1]*fUpwind_l[31]+0.25*fUpwind_l[10]*alphaDrSurf_l[19]+0.2500000000000001*alphaDrSurf_l[11]*fUpwind_l[18]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[17]; 
  Gdrag_l[45] = 0.159719141249985*alphaDrSurf_l[12]*fUpwind_l[45]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[45]+0.25*alphaDrSurf_l[0]*fUpwind_l[45]+0.2*alphaDrSurf_l[5]*fUpwind_l[44]+0.159719141249985*alphaDrSurf_l[20]*fUpwind_l[38]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[38]+0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[37]+0.2*alphaDrSurf_l[19]*fUpwind_l[31]+0.223606797749979*alphaDrSurf_l[2]*fUpwind_l[31]+0.25*fUpwind_l[10]*alphaDrSurf_l[20]+0.223606797749979*alphaDrSurf_l[5]*fUpwind_l[18]+0.2500000000000001*alphaDrSurf_l[12]*fUpwind_l[17]; 
  Gdrag_l[46] = 0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[46]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[46]+0.25*alphaDrSurf_l[0]*fUpwind_l[46]+0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[40]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[40]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[39]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[39]+0.25*alphaDrSurf_l[5]*fUpwind_l[27]; 
  Gdrag_l[47] = 0.223606797749979*alphaDrSurf_l[12]*fUpwind_l[47]+0.223606797749979*alphaDrSurf_l[11]*fUpwind_l[47]+0.25*alphaDrSurf_l[0]*fUpwind_l[47]+0.223606797749979*alphaDrSurf_l[20]*fUpwind_l[43]+0.2500000000000001*alphaDrSurf_l[1]*fUpwind_l[43]+0.223606797749979*alphaDrSurf_l[19]*fUpwind_l[42]+0.2500000000000001*alphaDrSurf_l[2]*fUpwind_l[42]+0.25*alphaDrSurf_l[5]*fUpwind_l[30]; 

  Gdrag_r[0] = 0.25*alphaDrSurf_r[20]*fUpwind_r[20]+0.25*alphaDrSurf_r[19]*fUpwind_r[19]+0.25*alphaDrSurf_r[12]*fUpwind_r[12]+0.25*alphaDrSurf_r[11]*fUpwind_r[11]+0.25*alphaDrSurf_r[5]*fUpwind_r[5]+0.25*alphaDrSurf_r[2]*fUpwind_r[2]+0.25*alphaDrSurf_r[1]*fUpwind_r[1]+0.25*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.2500000000000001*alphaDrSurf_r[12]*fUpwind_r[20]+0.2500000000000001*fUpwind_r[12]*alphaDrSurf_r[20]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[19]+0.223606797749979*fUpwind_r[5]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[11]+0.223606797749979*fUpwind_r[1]*alphaDrSurf_r[11]+0.25*alphaDrSurf_r[2]*fUpwind_r[5]+0.25*fUpwind_r[2]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[0]*fUpwind_r[1]+0.25*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[20]+0.223606797749979*fUpwind_r[5]*alphaDrSurf_r[20]+0.2500000000000001*alphaDrSurf_r[11]*fUpwind_r[19]+0.2500000000000001*fUpwind_r[11]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[12]+0.223606797749979*fUpwind_r[2]*alphaDrSurf_r[12]+0.25*alphaDrSurf_r[1]*fUpwind_r[5]+0.25*fUpwind_r[1]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[0]*fUpwind_r[2]+0.25*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Gdrag_r[3] = 0.2500000000000001*alphaDrSurf_r[20]*fUpwind_r[33]+0.2500000000000001*alphaDrSurf_r[19]*fUpwind_r[32]+0.2500000000000001*alphaDrSurf_r[12]*fUpwind_r[22]+0.2500000000000001*alphaDrSurf_r[11]*fUpwind_r[21]+0.25*alphaDrSurf_r[5]*fUpwind_r[15]+0.25*alphaDrSurf_r[2]*fUpwind_r[7]+0.25*alphaDrSurf_r[1]*fUpwind_r[6]+0.25*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Gdrag_r[4] = 0.2500000000000001*alphaDrSurf_r[20]*fUpwind_r[36]+0.2500000000000001*alphaDrSurf_r[19]*fUpwind_r[35]+0.2500000000000001*alphaDrSurf_r[12]*fUpwind_r[26]+0.2500000000000001*alphaDrSurf_r[11]*fUpwind_r[25]+0.25*alphaDrSurf_r[5]*fUpwind_r[16]+0.25*alphaDrSurf_r[2]*fUpwind_r[9]+0.25*alphaDrSurf_r[1]*fUpwind_r[8]+0.25*alphaDrSurf_r[0]*fUpwind_r[4]; 
  Gdrag_r[5] = 0.2*alphaDrSurf_r[19]*fUpwind_r[20]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[20]+0.2*fUpwind_r[19]*alphaDrSurf_r[20]+0.223606797749979*fUpwind_r[2]*alphaDrSurf_r[20]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[19]+0.223606797749979*fUpwind_r[1]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[12]+0.223606797749979*fUpwind_r[5]*alphaDrSurf_r[12]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[11]+0.223606797749979*fUpwind_r[5]*alphaDrSurf_r[11]+0.25*alphaDrSurf_r[0]*fUpwind_r[5]+0.25*fUpwind_r[0]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[1]*fUpwind_r[2]+0.25*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Gdrag_r[6] = 0.25*alphaDrSurf_r[12]*fUpwind_r[33]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[32]+0.25*alphaDrSurf_r[20]*fUpwind_r[22]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[21]+0.223606797749979*fUpwind_r[15]*alphaDrSurf_r[19]+0.25*alphaDrSurf_r[2]*fUpwind_r[15]+0.223606797749979*fUpwind_r[6]*alphaDrSurf_r[11]+0.25*alphaDrSurf_r[5]*fUpwind_r[7]+0.25*alphaDrSurf_r[0]*fUpwind_r[6]+0.25*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Gdrag_r[7] = 0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[33]+0.25*alphaDrSurf_r[11]*fUpwind_r[32]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[22]+0.25*alphaDrSurf_r[19]*fUpwind_r[21]+0.223606797749979*fUpwind_r[15]*alphaDrSurf_r[20]+0.25*alphaDrSurf_r[1]*fUpwind_r[15]+0.223606797749979*fUpwind_r[7]*alphaDrSurf_r[12]+0.25*alphaDrSurf_r[0]*fUpwind_r[7]+0.25*alphaDrSurf_r[5]*fUpwind_r[6]+0.25*alphaDrSurf_r[2]*fUpwind_r[3]; 
  Gdrag_r[8] = 0.25*alphaDrSurf_r[12]*fUpwind_r[36]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[35]+0.25*alphaDrSurf_r[20]*fUpwind_r[26]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[25]+0.223606797749979*fUpwind_r[16]*alphaDrSurf_r[19]+0.25*alphaDrSurf_r[2]*fUpwind_r[16]+0.223606797749979*fUpwind_r[8]*alphaDrSurf_r[11]+0.25*alphaDrSurf_r[5]*fUpwind_r[9]+0.25*alphaDrSurf_r[0]*fUpwind_r[8]+0.25*alphaDrSurf_r[1]*fUpwind_r[4]; 
  Gdrag_r[9] = 0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[36]+0.25*alphaDrSurf_r[11]*fUpwind_r[35]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[26]+0.25*alphaDrSurf_r[19]*fUpwind_r[25]+0.223606797749979*fUpwind_r[16]*alphaDrSurf_r[20]+0.25*alphaDrSurf_r[1]*fUpwind_r[16]+0.223606797749979*fUpwind_r[9]*alphaDrSurf_r[12]+0.25*alphaDrSurf_r[0]*fUpwind_r[9]+0.25*alphaDrSurf_r[5]*fUpwind_r[8]+0.25*alphaDrSurf_r[2]*fUpwind_r[4]; 
  Gdrag_r[10] = 0.25*alphaDrSurf_r[20]*fUpwind_r[45]+0.25*alphaDrSurf_r[19]*fUpwind_r[44]+0.25*alphaDrSurf_r[12]*fUpwind_r[38]+0.25*alphaDrSurf_r[11]*fUpwind_r[37]+0.25*alphaDrSurf_r[5]*fUpwind_r[31]+0.25*alphaDrSurf_r[2]*fUpwind_r[18]+0.25*alphaDrSurf_r[1]*fUpwind_r[17]+0.25*alphaDrSurf_r[0]*fUpwind_r[10]; 
  Gdrag_r[11] = 0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[20]+0.159719141249985*alphaDrSurf_r[19]*fUpwind_r[19]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[19]+0.2500000000000001*fUpwind_r[2]*alphaDrSurf_r[19]+0.159719141249985*alphaDrSurf_r[11]*fUpwind_r[11]+0.25*alphaDrSurf_r[0]*fUpwind_r[11]+0.25*fUpwind_r[0]*alphaDrSurf_r[11]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[5]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Gdrag_r[12] = 0.159719141249985*alphaDrSurf_r[20]*fUpwind_r[20]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[20]+0.2500000000000001*fUpwind_r[1]*alphaDrSurf_r[20]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[19]+0.159719141249985*alphaDrSurf_r[12]*fUpwind_r[12]+0.25*alphaDrSurf_r[0]*fUpwind_r[12]+0.25*fUpwind_r[0]*alphaDrSurf_r[12]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[5]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[2]; 
  Gdrag_r[13] = 0.25*alphaDrSurf_r[5]*fUpwind_r[34]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[24]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[23]+0.25*alphaDrSurf_r[0]*fUpwind_r[13]; 
  Gdrag_r[14] = 0.25*alphaDrSurf_r[5]*fUpwind_r[41]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[29]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[28]+0.25*alphaDrSurf_r[0]*fUpwind_r[14]; 
  Gdrag_r[15] = 0.2*alphaDrSurf_r[19]*fUpwind_r[33]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[33]+0.2*alphaDrSurf_r[20]*fUpwind_r[32]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[32]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[22]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[21]+0.223606797749979*fUpwind_r[7]*alphaDrSurf_r[20]+0.223606797749979*fUpwind_r[6]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[15]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[15]+0.25*alphaDrSurf_r[0]*fUpwind_r[15]+0.25*alphaDrSurf_r[1]*fUpwind_r[7]+0.25*alphaDrSurf_r[2]*fUpwind_r[6]+0.25*fUpwind_r[3]*alphaDrSurf_r[5]; 
  Gdrag_r[16] = 0.2*alphaDrSurf_r[19]*fUpwind_r[36]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[36]+0.2*alphaDrSurf_r[20]*fUpwind_r[35]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[35]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[26]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[25]+0.223606797749979*fUpwind_r[9]*alphaDrSurf_r[20]+0.223606797749979*fUpwind_r[8]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[16]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[16]+0.25*alphaDrSurf_r[0]*fUpwind_r[16]+0.25*alphaDrSurf_r[1]*fUpwind_r[9]+0.25*alphaDrSurf_r[2]*fUpwind_r[8]+0.25*fUpwind_r[4]*alphaDrSurf_r[5]; 
  Gdrag_r[17] = 0.2500000000000001*alphaDrSurf_r[12]*fUpwind_r[45]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[44]+0.2500000000000001*alphaDrSurf_r[20]*fUpwind_r[38]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[37]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[31]+0.25*alphaDrSurf_r[2]*fUpwind_r[31]+0.25*alphaDrSurf_r[5]*fUpwind_r[18]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[17]+0.25*alphaDrSurf_r[0]*fUpwind_r[17]+0.25*alphaDrSurf_r[1]*fUpwind_r[10]; 
  Gdrag_r[18] = 0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[45]+0.2500000000000001*alphaDrSurf_r[11]*fUpwind_r[44]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[38]+0.2500000000000001*alphaDrSurf_r[19]*fUpwind_r[37]+0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[31]+0.25*alphaDrSurf_r[1]*fUpwind_r[31]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[18]+0.25*alphaDrSurf_r[0]*fUpwind_r[18]+0.25*alphaDrSurf_r[5]*fUpwind_r[17]+0.25*alphaDrSurf_r[2]*fUpwind_r[10]; 
  Gdrag_r[19] = 0.2*alphaDrSurf_r[5]*fUpwind_r[20]+0.2*fUpwind_r[5]*alphaDrSurf_r[20]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[19]+0.159719141249985*alphaDrSurf_r[11]*fUpwind_r[19]+0.25*alphaDrSurf_r[0]*fUpwind_r[19]+0.223606797749979*fUpwind_r[12]*alphaDrSurf_r[19]+0.159719141249985*fUpwind_r[11]*alphaDrSurf_r[19]+0.25*fUpwind_r[0]*alphaDrSurf_r[19]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[11]+0.2500000000000001*fUpwind_r[2]*alphaDrSurf_r[11]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[5]+0.223606797749979*fUpwind_r[1]*alphaDrSurf_r[5]; 
  Gdrag_r[20] = 0.159719141249985*alphaDrSurf_r[12]*fUpwind_r[20]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[20]+0.25*alphaDrSurf_r[0]*fUpwind_r[20]+0.159719141249985*fUpwind_r[12]*alphaDrSurf_r[20]+0.223606797749979*fUpwind_r[11]*alphaDrSurf_r[20]+0.25*fUpwind_r[0]*alphaDrSurf_r[20]+0.2*alphaDrSurf_r[5]*fUpwind_r[19]+0.2*fUpwind_r[5]*alphaDrSurf_r[19]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[12]+0.2500000000000001*fUpwind_r[1]*alphaDrSurf_r[12]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[5]+0.223606797749979*fUpwind_r[2]*alphaDrSurf_r[5]; 
  Gdrag_r[21] = 0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[33]+0.159719141249985*alphaDrSurf_r[19]*fUpwind_r[32]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[32]+0.159719141249985*alphaDrSurf_r[11]*fUpwind_r[21]+0.25*alphaDrSurf_r[0]*fUpwind_r[21]+0.25*fUpwind_r[7]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[15]+0.2500000000000001*fUpwind_r[3]*alphaDrSurf_r[11]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[6]; 
  Gdrag_r[22] = 0.159719141249985*alphaDrSurf_r[20]*fUpwind_r[33]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[33]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[32]+0.159719141249985*alphaDrSurf_r[12]*fUpwind_r[22]+0.25*alphaDrSurf_r[0]*fUpwind_r[22]+0.25*fUpwind_r[6]*alphaDrSurf_r[20]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[15]+0.2500000000000001*fUpwind_r[3]*alphaDrSurf_r[12]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[7]; 
  Gdrag_r[23] = 0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[34]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[34]+0.25*alphaDrSurf_r[5]*fUpwind_r[24]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[23]+0.25*alphaDrSurf_r[0]*fUpwind_r[23]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[13]; 
  Gdrag_r[24] = 0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[34]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[34]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[24]+0.25*alphaDrSurf_r[0]*fUpwind_r[24]+0.25*alphaDrSurf_r[5]*fUpwind_r[23]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[13]; 
  Gdrag_r[25] = 0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[36]+0.159719141249985*alphaDrSurf_r[19]*fUpwind_r[35]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[35]+0.159719141249985*alphaDrSurf_r[11]*fUpwind_r[25]+0.25*alphaDrSurf_r[0]*fUpwind_r[25]+0.25*fUpwind_r[9]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[16]+0.2500000000000001*fUpwind_r[4]*alphaDrSurf_r[11]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[8]; 
  Gdrag_r[26] = 0.159719141249985*alphaDrSurf_r[20]*fUpwind_r[36]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[36]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[35]+0.159719141249985*alphaDrSurf_r[12]*fUpwind_r[26]+0.25*alphaDrSurf_r[0]*fUpwind_r[26]+0.25*fUpwind_r[8]*alphaDrSurf_r[20]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[16]+0.2500000000000001*fUpwind_r[4]*alphaDrSurf_r[12]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[9]; 
  Gdrag_r[27] = 0.25*alphaDrSurf_r[5]*fUpwind_r[46]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[40]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[39]+0.25*alphaDrSurf_r[0]*fUpwind_r[27]; 
  Gdrag_r[28] = 0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[41]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[41]+0.25*alphaDrSurf_r[5]*fUpwind_r[29]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[28]+0.25*alphaDrSurf_r[0]*fUpwind_r[28]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[14]; 
  Gdrag_r[29] = 0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[41]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[41]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[29]+0.25*alphaDrSurf_r[0]*fUpwind_r[29]+0.25*alphaDrSurf_r[5]*fUpwind_r[28]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[14]; 
  Gdrag_r[30] = 0.25*alphaDrSurf_r[5]*fUpwind_r[47]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[43]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[42]+0.25*alphaDrSurf_r[0]*fUpwind_r[30]; 
  Gdrag_r[31] = 0.2*alphaDrSurf_r[19]*fUpwind_r[45]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[45]+0.2*alphaDrSurf_r[20]*fUpwind_r[44]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[44]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[38]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[37]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[31]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[31]+0.25*alphaDrSurf_r[0]*fUpwind_r[31]+0.223606797749979*fUpwind_r[18]*alphaDrSurf_r[20]+0.223606797749979*fUpwind_r[17]*alphaDrSurf_r[19]+0.25*alphaDrSurf_r[1]*fUpwind_r[18]+0.25*alphaDrSurf_r[2]*fUpwind_r[17]+0.25*alphaDrSurf_r[5]*fUpwind_r[10]; 
  Gdrag_r[32] = 0.2*alphaDrSurf_r[5]*fUpwind_r[33]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[32]+0.159719141249985*alphaDrSurf_r[11]*fUpwind_r[32]+0.25*alphaDrSurf_r[0]*fUpwind_r[32]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[22]+0.159719141249985*alphaDrSurf_r[19]*fUpwind_r[21]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[21]+0.2*fUpwind_r[15]*alphaDrSurf_r[20]+0.2500000000000001*fUpwind_r[3]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[15]+0.25*fUpwind_r[7]*alphaDrSurf_r[11]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[6]; 
  Gdrag_r[33] = 0.159719141249985*alphaDrSurf_r[12]*fUpwind_r[33]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[33]+0.25*alphaDrSurf_r[0]*fUpwind_r[33]+0.2*alphaDrSurf_r[5]*fUpwind_r[32]+0.159719141249985*alphaDrSurf_r[20]*fUpwind_r[22]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[22]+0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[21]+0.2500000000000001*fUpwind_r[3]*alphaDrSurf_r[20]+0.2*fUpwind_r[15]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[15]+0.25*fUpwind_r[6]*alphaDrSurf_r[12]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[7]; 
  Gdrag_r[34] = 0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[34]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[34]+0.25*alphaDrSurf_r[0]*fUpwind_r[34]+0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[24]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[24]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[23]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[23]+0.25*alphaDrSurf_r[5]*fUpwind_r[13]; 
  Gdrag_r[35] = 0.2*alphaDrSurf_r[5]*fUpwind_r[36]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[35]+0.159719141249985*alphaDrSurf_r[11]*fUpwind_r[35]+0.25*alphaDrSurf_r[0]*fUpwind_r[35]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[26]+0.159719141249985*alphaDrSurf_r[19]*fUpwind_r[25]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[25]+0.2*fUpwind_r[16]*alphaDrSurf_r[20]+0.2500000000000001*fUpwind_r[4]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[16]+0.25*fUpwind_r[9]*alphaDrSurf_r[11]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[8]; 
  Gdrag_r[36] = 0.159719141249985*alphaDrSurf_r[12]*fUpwind_r[36]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[36]+0.25*alphaDrSurf_r[0]*fUpwind_r[36]+0.2*alphaDrSurf_r[5]*fUpwind_r[35]+0.159719141249985*alphaDrSurf_r[20]*fUpwind_r[26]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[26]+0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[25]+0.2500000000000001*fUpwind_r[4]*alphaDrSurf_r[20]+0.2*fUpwind_r[16]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[16]+0.25*fUpwind_r[8]*alphaDrSurf_r[12]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[9]; 
  Gdrag_r[37] = 0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[45]+0.159719141249985*alphaDrSurf_r[19]*fUpwind_r[44]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[44]+0.159719141249985*alphaDrSurf_r[11]*fUpwind_r[37]+0.25*alphaDrSurf_r[0]*fUpwind_r[37]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[31]+0.2500000000000001*fUpwind_r[18]*alphaDrSurf_r[19]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[17]+0.25*fUpwind_r[10]*alphaDrSurf_r[11]; 
  Gdrag_r[38] = 0.159719141249985*alphaDrSurf_r[20]*fUpwind_r[45]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[45]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[44]+0.159719141249985*alphaDrSurf_r[12]*fUpwind_r[38]+0.25*alphaDrSurf_r[0]*fUpwind_r[38]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[31]+0.2500000000000001*fUpwind_r[17]*alphaDrSurf_r[20]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[18]+0.25*fUpwind_r[10]*alphaDrSurf_r[12]; 
  Gdrag_r[39] = 0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[46]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[46]+0.25*alphaDrSurf_r[5]*fUpwind_r[40]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[39]+0.25*alphaDrSurf_r[0]*fUpwind_r[39]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[27]; 
  Gdrag_r[40] = 0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[46]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[46]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[40]+0.25*alphaDrSurf_r[0]*fUpwind_r[40]+0.25*alphaDrSurf_r[5]*fUpwind_r[39]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[27]; 
  Gdrag_r[41] = 0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[41]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[41]+0.25*alphaDrSurf_r[0]*fUpwind_r[41]+0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[29]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[29]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[28]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[28]+0.25*alphaDrSurf_r[5]*fUpwind_r[14]; 
  Gdrag_r[42] = 0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[47]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[47]+0.25*alphaDrSurf_r[5]*fUpwind_r[43]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[42]+0.25*alphaDrSurf_r[0]*fUpwind_r[42]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[30]; 
  Gdrag_r[43] = 0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[47]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[47]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[43]+0.25*alphaDrSurf_r[0]*fUpwind_r[43]+0.25*alphaDrSurf_r[5]*fUpwind_r[42]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[30]; 
  Gdrag_r[44] = 0.2*alphaDrSurf_r[5]*fUpwind_r[45]+0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[44]+0.159719141249985*alphaDrSurf_r[11]*fUpwind_r[44]+0.25*alphaDrSurf_r[0]*fUpwind_r[44]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[38]+0.159719141249985*alphaDrSurf_r[19]*fUpwind_r[37]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[37]+0.2*alphaDrSurf_r[20]*fUpwind_r[31]+0.223606797749979*alphaDrSurf_r[1]*fUpwind_r[31]+0.25*fUpwind_r[10]*alphaDrSurf_r[19]+0.2500000000000001*alphaDrSurf_r[11]*fUpwind_r[18]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[17]; 
  Gdrag_r[45] = 0.159719141249985*alphaDrSurf_r[12]*fUpwind_r[45]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[45]+0.25*alphaDrSurf_r[0]*fUpwind_r[45]+0.2*alphaDrSurf_r[5]*fUpwind_r[44]+0.159719141249985*alphaDrSurf_r[20]*fUpwind_r[38]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[38]+0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[37]+0.2*alphaDrSurf_r[19]*fUpwind_r[31]+0.223606797749979*alphaDrSurf_r[2]*fUpwind_r[31]+0.25*fUpwind_r[10]*alphaDrSurf_r[20]+0.223606797749979*alphaDrSurf_r[5]*fUpwind_r[18]+0.2500000000000001*alphaDrSurf_r[12]*fUpwind_r[17]; 
  Gdrag_r[46] = 0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[46]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[46]+0.25*alphaDrSurf_r[0]*fUpwind_r[46]+0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[40]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[40]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[39]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[39]+0.25*alphaDrSurf_r[5]*fUpwind_r[27]; 
  Gdrag_r[47] = 0.223606797749979*alphaDrSurf_r[12]*fUpwind_r[47]+0.223606797749979*alphaDrSurf_r[11]*fUpwind_r[47]+0.25*alphaDrSurf_r[0]*fUpwind_r[47]+0.223606797749979*alphaDrSurf_r[20]*fUpwind_r[43]+0.2500000000000001*alphaDrSurf_r[1]*fUpwind_r[43]+0.223606797749979*alphaDrSurf_r[19]*fUpwind_r[42]+0.2500000000000001*alphaDrSurf_r[2]*fUpwind_r[42]+0.25*alphaDrSurf_r[5]*fUpwind_r[30]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[4] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[8] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[9] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[10] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[11] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[12] += 0.7071067811865475*Gdrag_r[8]*rdv2-0.7071067811865475*Gdrag_l[8]*rdv2; 
  out[13] += 0.7071067811865475*Gdrag_r[9]*rdv2-0.7071067811865475*Gdrag_l[9]*rdv2; 
  out[14] += 0.7071067811865475*Gdrag_r[10]*rdv2-0.7071067811865475*Gdrag_l[10]*rdv2; 
  out[15] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[16] += 0.7071067811865475*Gdrag_r[11]*rdv2-0.7071067811865475*Gdrag_l[11]*rdv2; 
  out[17] += 0.7071067811865475*Gdrag_r[12]*rdv2-0.7071067811865475*Gdrag_l[12]*rdv2; 
  out[18] += 0.7071067811865475*Gdrag_r[13]*rdv2-0.7071067811865475*Gdrag_l[13]*rdv2; 
  out[19] += 1.58113883008419*Gdrag_r[0]*rdv2-1.58113883008419*Gdrag_l[0]*rdv2; 
  out[20] += 0.7071067811865475*Gdrag_r[14]*rdv2-0.7071067811865475*Gdrag_l[14]*rdv2; 
  out[21] += 0.7071067811865475*Gdrag_r[15]*rdv2-0.7071067811865475*Gdrag_l[15]*rdv2; 
  out[22] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[23] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[24] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
  out[25] += 0.7071067811865475*Gdrag_r[16]*rdv2-0.7071067811865475*Gdrag_l[16]*rdv2; 
  out[26] += 0.7071067811865475*Gdrag_r[17]*rdv2-0.7071067811865475*Gdrag_l[17]*rdv2; 
  out[27] += 0.7071067811865475*Gdrag_r[18]*rdv2-0.7071067811865475*Gdrag_l[18]*rdv2; 
  out[28] += 1.224744871391589*Gdrag_r[8]*rdv2+1.224744871391589*Gdrag_l[8]*rdv2; 
  out[29] += 1.224744871391589*Gdrag_r[9]*rdv2+1.224744871391589*Gdrag_l[9]*rdv2; 
  out[30] += 1.224744871391589*Gdrag_r[10]*rdv2+1.224744871391589*Gdrag_l[10]*rdv2; 
  out[31] += 0.7071067811865475*Gdrag_r[19]*rdv2-0.7071067811865475*Gdrag_l[19]*rdv2; 
  out[32] += 0.7071067811865475*Gdrag_r[20]*rdv2-0.7071067811865475*Gdrag_l[20]*rdv2; 
  out[33] += 0.7071067811865475*Gdrag_r[21]*rdv2-0.7071067811865475*Gdrag_l[21]*rdv2; 
  out[34] += 0.7071067811865475*Gdrag_r[22]*rdv2-0.7071067811865475*Gdrag_l[22]*rdv2; 
  out[35] += 0.7071067811865475*Gdrag_r[23]*rdv2-0.7071067811865475*Gdrag_l[23]*rdv2; 
  out[36] += 0.7071067811865475*Gdrag_r[24]*rdv2-0.7071067811865475*Gdrag_l[24]*rdv2; 
  out[37] += 1.224744871391589*Gdrag_r[11]*rdv2+1.224744871391589*Gdrag_l[11]*rdv2; 
  out[38] += 1.224744871391589*Gdrag_r[12]*rdv2+1.224744871391589*Gdrag_l[12]*rdv2; 
  out[39] += 1.224744871391589*Gdrag_r[13]*rdv2+1.224744871391589*Gdrag_l[13]*rdv2; 
  out[40] += 1.58113883008419*Gdrag_r[1]*rdv2-1.58113883008419*Gdrag_l[1]*rdv2; 
  out[41] += 1.58113883008419*Gdrag_r[2]*rdv2-1.58113883008419*Gdrag_l[2]*rdv2; 
  out[42] += 1.58113883008419*Gdrag_r[3]*rdv2-1.58113883008419*Gdrag_l[3]*rdv2; 
  out[43] += 0.7071067811865475*Gdrag_r[25]*rdv2-0.7071067811865475*Gdrag_l[25]*rdv2; 
  out[44] += 0.7071067811865475*Gdrag_r[26]*rdv2-0.7071067811865475*Gdrag_l[26]*rdv2; 
  out[45] += 0.7071067811865475*Gdrag_r[27]*rdv2-0.7071067811865475*Gdrag_l[27]*rdv2; 
  out[46] += 1.58113883008419*Gdrag_r[4]*rdv2-1.58113883008419*Gdrag_l[4]*rdv2; 
  out[47] += 0.7071067811865475*Gdrag_r[28]*rdv2-0.7071067811865475*Gdrag_l[28]*rdv2; 
  out[48] += 0.7071067811865475*Gdrag_r[29]*rdv2-0.7071067811865475*Gdrag_l[29]*rdv2; 
  out[49] += 0.7071067811865475*Gdrag_r[30]*rdv2-0.7071067811865475*Gdrag_l[30]*rdv2; 
  out[50] += 1.224744871391589*Gdrag_r[14]*rdv2+1.224744871391589*Gdrag_l[14]*rdv2; 
  out[51] += 1.224744871391589*Gdrag_r[15]*rdv2+1.224744871391589*Gdrag_l[15]*rdv2; 
  out[52] += 0.7071067811865475*Gdrag_r[31]*rdv2-0.7071067811865475*Gdrag_l[31]*rdv2; 
  out[53] += 1.224744871391589*Gdrag_r[16]*rdv2+1.224744871391589*Gdrag_l[16]*rdv2; 
  out[54] += 1.224744871391589*Gdrag_r[17]*rdv2+1.224744871391589*Gdrag_l[17]*rdv2; 
  out[55] += 1.224744871391589*Gdrag_r[18]*rdv2+1.224744871391589*Gdrag_l[18]*rdv2; 
  out[56] += 0.7071067811865475*Gdrag_r[32]*rdv2-0.7071067811865475*Gdrag_l[32]*rdv2; 
  out[57] += 0.7071067811865475*Gdrag_r[33]*rdv2-0.7071067811865475*Gdrag_l[33]*rdv2; 
  out[58] += 0.7071067811865475*Gdrag_r[34]*rdv2-0.7071067811865475*Gdrag_l[34]*rdv2; 
  out[59] += 1.224744871391589*Gdrag_r[19]*rdv2+1.224744871391589*Gdrag_l[19]*rdv2; 
  out[60] += 1.224744871391589*Gdrag_r[20]*rdv2+1.224744871391589*Gdrag_l[20]*rdv2; 
  out[61] += 1.224744871391589*Gdrag_r[21]*rdv2+1.224744871391589*Gdrag_l[21]*rdv2; 
  out[62] += 1.224744871391589*Gdrag_r[22]*rdv2+1.224744871391589*Gdrag_l[22]*rdv2; 
  out[63] += 1.224744871391589*Gdrag_r[23]*rdv2+1.224744871391589*Gdrag_l[23]*rdv2; 
  out[64] += 1.224744871391589*Gdrag_r[24]*rdv2+1.224744871391589*Gdrag_l[24]*rdv2; 
  out[65] += 1.58113883008419*Gdrag_r[5]*rdv2-1.58113883008419*Gdrag_l[5]*rdv2; 
  out[66] += 1.58113883008419*Gdrag_r[6]*rdv2-1.58113883008419*Gdrag_l[6]*rdv2; 
  out[67] += 1.58113883008419*Gdrag_r[7]*rdv2-1.58113883008419*Gdrag_l[7]*rdv2; 
  out[68] += 0.7071067811865475*Gdrag_r[35]*rdv2-0.7071067811865475*Gdrag_l[35]*rdv2; 
  out[69] += 0.7071067811865475*Gdrag_r[36]*rdv2-0.7071067811865475*Gdrag_l[36]*rdv2; 
  out[70] += 0.7071067811865475*Gdrag_r[37]*rdv2-0.7071067811865475*Gdrag_l[37]*rdv2; 
  out[71] += 0.7071067811865475*Gdrag_r[38]*rdv2-0.7071067811865475*Gdrag_l[38]*rdv2; 
  out[72] += 0.7071067811865475*Gdrag_r[39]*rdv2-0.7071067811865475*Gdrag_l[39]*rdv2; 
  out[73] += 0.7071067811865475*Gdrag_r[40]*rdv2-0.7071067811865475*Gdrag_l[40]*rdv2; 
  out[74] += 1.224744871391589*Gdrag_r[25]*rdv2+1.224744871391589*Gdrag_l[25]*rdv2; 
  out[75] += 1.224744871391589*Gdrag_r[26]*rdv2+1.224744871391589*Gdrag_l[26]*rdv2; 
  out[76] += 1.224744871391589*Gdrag_r[27]*rdv2+1.224744871391589*Gdrag_l[27]*rdv2; 
  out[77] += 1.58113883008419*Gdrag_r[8]*rdv2-1.58113883008419*Gdrag_l[8]*rdv2; 
  out[78] += 1.58113883008419*Gdrag_r[9]*rdv2-1.58113883008419*Gdrag_l[9]*rdv2; 
  out[79] += 1.58113883008419*Gdrag_r[10]*rdv2-1.58113883008419*Gdrag_l[10]*rdv2; 
  out[80] += 0.7071067811865475*Gdrag_r[41]*rdv2-0.7071067811865475*Gdrag_l[41]*rdv2; 
  out[81] += 0.7071067811865475*Gdrag_r[42]*rdv2-0.7071067811865475*Gdrag_l[42]*rdv2; 
  out[82] += 0.7071067811865475*Gdrag_r[43]*rdv2-0.7071067811865475*Gdrag_l[43]*rdv2; 
  out[83] += 1.224744871391589*Gdrag_r[28]*rdv2+1.224744871391589*Gdrag_l[28]*rdv2; 
  out[84] += 1.224744871391589*Gdrag_r[29]*rdv2+1.224744871391589*Gdrag_l[29]*rdv2; 
  out[85] += 1.224744871391589*Gdrag_r[30]*rdv2+1.224744871391589*Gdrag_l[30]*rdv2; 
  out[86] += 1.224744871391589*Gdrag_r[31]*rdv2+1.224744871391589*Gdrag_l[31]*rdv2; 
  out[87] += 1.224744871391589*Gdrag_r[32]*rdv2+1.224744871391589*Gdrag_l[32]*rdv2; 
  out[88] += 1.224744871391589*Gdrag_r[33]*rdv2+1.224744871391589*Gdrag_l[33]*rdv2; 
  out[89] += 1.224744871391589*Gdrag_r[34]*rdv2+1.224744871391589*Gdrag_l[34]*rdv2; 
  out[90] += 1.58113883008419*Gdrag_r[15]*rdv2-1.58113883008419*Gdrag_l[15]*rdv2; 
  out[91] += 0.7071067811865475*Gdrag_r[44]*rdv2-0.7071067811865475*Gdrag_l[44]*rdv2; 
  out[92] += 0.7071067811865475*Gdrag_r[45]*rdv2-0.7071067811865475*Gdrag_l[45]*rdv2; 
  out[93] += 0.7071067811865475*Gdrag_r[46]*rdv2-0.7071067811865475*Gdrag_l[46]*rdv2; 
  out[94] += 1.224744871391589*Gdrag_r[35]*rdv2+1.224744871391589*Gdrag_l[35]*rdv2; 
  out[95] += 1.224744871391589*Gdrag_r[36]*rdv2+1.224744871391589*Gdrag_l[36]*rdv2; 
  out[96] += 1.224744871391589*Gdrag_r[37]*rdv2+1.224744871391589*Gdrag_l[37]*rdv2; 
  out[97] += 1.224744871391589*Gdrag_r[38]*rdv2+1.224744871391589*Gdrag_l[38]*rdv2; 
  out[98] += 1.224744871391589*Gdrag_r[39]*rdv2+1.224744871391589*Gdrag_l[39]*rdv2; 
  out[99] += 1.224744871391589*Gdrag_r[40]*rdv2+1.224744871391589*Gdrag_l[40]*rdv2; 
  out[100] += 1.58113883008419*Gdrag_r[16]*rdv2-1.58113883008419*Gdrag_l[16]*rdv2; 
  out[101] += 1.58113883008419*Gdrag_r[17]*rdv2-1.58113883008419*Gdrag_l[17]*rdv2; 
  out[102] += 1.58113883008419*Gdrag_r[18]*rdv2-1.58113883008419*Gdrag_l[18]*rdv2; 
  out[103] += 0.7071067811865475*Gdrag_r[47]*rdv2-0.7071067811865475*Gdrag_l[47]*rdv2; 
  out[104] += 1.224744871391589*Gdrag_r[41]*rdv2+1.224744871391589*Gdrag_l[41]*rdv2; 
  out[105] += 1.224744871391589*Gdrag_r[42]*rdv2+1.224744871391589*Gdrag_l[42]*rdv2; 
  out[106] += 1.224744871391589*Gdrag_r[43]*rdv2+1.224744871391589*Gdrag_l[43]*rdv2; 
  out[107] += 1.224744871391589*Gdrag_r[44]*rdv2+1.224744871391589*Gdrag_l[44]*rdv2; 
  out[108] += 1.224744871391589*Gdrag_r[45]*rdv2+1.224744871391589*Gdrag_l[45]*rdv2; 
  out[109] += 1.224744871391589*Gdrag_r[46]*rdv2+1.224744871391589*Gdrag_l[46]*rdv2; 
  out[110] += 1.58113883008419*Gdrag_r[31]*rdv2-1.58113883008419*Gdrag_l[31]*rdv2; 
  out[111] += 1.224744871391589*Gdrag_r[47]*rdv2+1.224744871391589*Gdrag_l[47]*rdv2; 
} 
