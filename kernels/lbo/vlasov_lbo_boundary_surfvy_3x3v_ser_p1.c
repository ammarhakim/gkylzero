#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_3x3v_p1_surfvy_quad.h> 
GKYL_CU_DH void vlasov_lbo_boundary_surfvy_3x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
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
  double drag_incr[32] = {0.0}; 

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
    fUpwindQuad[0] = ser_3x3v_p1_surfvy_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_3x3v_p1_surfvy_quad_0(-1, fEdge); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_3x3v_p1_surfvy_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_3x3v_p1_surfvy_quad_1(-1, fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_3x3v_p1_surfvy_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_3x3v_p1_surfvy_quad_2(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_3x3v_p1_surfvy_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_3x3v_p1_surfvy_quad_3(-1, fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_3x3v_p1_surfvy_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_3x3v_p1_surfvy_quad_4(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_3x3v_p1_surfvy_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_3x3v_p1_surfvy_quad_5(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_3x3v_p1_surfvy_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_3x3v_p1_surfvy_quad_6(-1, fEdge); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_3x3v_p1_surfvy_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_3x3v_p1_surfvy_quad_7(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_3x3v_p1_surfvy_quad_8(1, fSkin); 
  } else { 

    fUpwindQuad[8] = ser_3x3v_p1_surfvy_quad_8(-1, fEdge); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_3x3v_p1_surfvy_quad_9(1, fSkin); 
  } else { 

    fUpwindQuad[9] = ser_3x3v_p1_surfvy_quad_9(-1, fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = ser_3x3v_p1_surfvy_quad_10(1, fSkin); 
  } else { 

    fUpwindQuad[10] = ser_3x3v_p1_surfvy_quad_10(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_3x3v_p1_surfvy_quad_11(1, fSkin); 
  } else { 

    fUpwindQuad[11] = ser_3x3v_p1_surfvy_quad_11(-1, fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_3x3v_p1_surfvy_quad_12(1, fSkin); 
  } else { 

    fUpwindQuad[12] = ser_3x3v_p1_surfvy_quad_12(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[13] = ser_3x3v_p1_surfvy_quad_13(1, fSkin); 
  } else { 

    fUpwindQuad[13] = ser_3x3v_p1_surfvy_quad_13(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_3x3v_p1_surfvy_quad_14(1, fSkin); 
  } else { 

    fUpwindQuad[14] = ser_3x3v_p1_surfvy_quad_14(-1, fEdge); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_3x3v_p1_surfvy_quad_15(1, fSkin); 
  } else { 

    fUpwindQuad[15] = ser_3x3v_p1_surfvy_quad_15(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[16] = ser_3x3v_p1_surfvy_quad_16(1, fSkin); 
  } else { 

    fUpwindQuad[16] = ser_3x3v_p1_surfvy_quad_16(-1, fEdge); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[17] = ser_3x3v_p1_surfvy_quad_17(1, fSkin); 
  } else { 

    fUpwindQuad[17] = ser_3x3v_p1_surfvy_quad_17(-1, fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_3x3v_p1_surfvy_quad_18(1, fSkin); 
  } else { 

    fUpwindQuad[18] = ser_3x3v_p1_surfvy_quad_18(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[19] = ser_3x3v_p1_surfvy_quad_19(1, fSkin); 
  } else { 

    fUpwindQuad[19] = ser_3x3v_p1_surfvy_quad_19(-1, fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[20] = ser_3x3v_p1_surfvy_quad_20(1, fSkin); 
  } else { 

    fUpwindQuad[20] = ser_3x3v_p1_surfvy_quad_20(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_3x3v_p1_surfvy_quad_21(1, fSkin); 
  } else { 

    fUpwindQuad[21] = ser_3x3v_p1_surfvy_quad_21(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[22] = ser_3x3v_p1_surfvy_quad_22(1, fSkin); 
  } else { 

    fUpwindQuad[22] = ser_3x3v_p1_surfvy_quad_22(-1, fEdge); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[23] = ser_3x3v_p1_surfvy_quad_23(1, fSkin); 
  } else { 

    fUpwindQuad[23] = ser_3x3v_p1_surfvy_quad_23(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_3x3v_p1_surfvy_quad_24(1, fSkin); 
  } else { 

    fUpwindQuad[24] = ser_3x3v_p1_surfvy_quad_24(-1, fEdge); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[25] = ser_3x3v_p1_surfvy_quad_25(1, fSkin); 
  } else { 

    fUpwindQuad[25] = ser_3x3v_p1_surfvy_quad_25(-1, fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[26] = ser_3x3v_p1_surfvy_quad_26(1, fSkin); 
  } else { 

    fUpwindQuad[26] = ser_3x3v_p1_surfvy_quad_26(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[27] = ser_3x3v_p1_surfvy_quad_27(1, fSkin); 
  } else { 

    fUpwindQuad[27] = ser_3x3v_p1_surfvy_quad_27(-1, fEdge); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[28] = ser_3x3v_p1_surfvy_quad_28(1, fSkin); 
  } else { 

    fUpwindQuad[28] = ser_3x3v_p1_surfvy_quad_28(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[29] = ser_3x3v_p1_surfvy_quad_29(1, fSkin); 
  } else { 

    fUpwindQuad[29] = ser_3x3v_p1_surfvy_quad_29(-1, fEdge); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[30] = ser_3x3v_p1_surfvy_quad_30(1, fSkin); 
  } else { 

    fUpwindQuad[30] = ser_3x3v_p1_surfvy_quad_30(-1, fEdge); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[31] = ser_3x3v_p1_surfvy_quad_31(1, fSkin); 
  } else { 

    fUpwindQuad[31] = ser_3x3v_p1_surfvy_quad_31(-1, fEdge); 
  } 

  fUpwind[0] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[6] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[7] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[12] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[13] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[14] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[15] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[16] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[17] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[18] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[19] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[20] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[21] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[22] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[23] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[24] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[25] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[26] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[27] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[28] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[29] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[30] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[31] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  drag_incr[0] = 0.1767766952966368*alphaDrSurf[16]*fUpwind[16]+0.1767766952966368*alphaDrSurf[8]*fUpwind[8]+0.1767766952966368*alphaDrSurf[7]*fUpwind[7]+0.1767766952966368*alphaDrSurf[6]*fUpwind[6]+0.1767766952966368*alphaDrSurf[3]*fUpwind[3]+0.1767766952966368*alphaDrSurf[2]*fUpwind[2]+0.1767766952966368*alphaDrSurf[1]*fUpwind[1]+0.1767766952966368*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.1767766952966368*alphaDrSurf[8]*fUpwind[16]+0.1767766952966368*fUpwind[8]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[3]*fUpwind[7]+0.1767766952966368*fUpwind[3]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[2]*fUpwind[6]+0.1767766952966368*fUpwind[2]*alphaDrSurf[6]+0.1767766952966368*alphaDrSurf[0]*fUpwind[1]+0.1767766952966368*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.1767766952966368*alphaDrSurf[7]*fUpwind[16]+0.1767766952966368*fUpwind[7]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[3]*fUpwind[8]+0.1767766952966368*fUpwind[3]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[1]*fUpwind[6]+0.1767766952966368*fUpwind[1]*alphaDrSurf[6]+0.1767766952966368*alphaDrSurf[0]*fUpwind[2]+0.1767766952966368*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.1767766952966368*alphaDrSurf[6]*fUpwind[16]+0.1767766952966368*fUpwind[6]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[2]*fUpwind[8]+0.1767766952966368*fUpwind[2]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[1]*fUpwind[7]+0.1767766952966368*fUpwind[1]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[0]*fUpwind[3]+0.1767766952966368*fUpwind[0]*alphaDrSurf[3]; 
  drag_incr[4] = 0.1767766952966368*alphaDrSurf[16]*fUpwind[26]+0.1767766952966368*alphaDrSurf[8]*fUpwind[19]+0.1767766952966368*alphaDrSurf[7]*fUpwind[18]+0.1767766952966368*alphaDrSurf[6]*fUpwind[17]+0.1767766952966368*alphaDrSurf[3]*fUpwind[11]+0.1767766952966368*alphaDrSurf[2]*fUpwind[10]+0.1767766952966368*alphaDrSurf[1]*fUpwind[9]+0.1767766952966368*alphaDrSurf[0]*fUpwind[4]; 
  drag_incr[5] = 0.1767766952966368*alphaDrSurf[16]*fUpwind[27]+0.1767766952966368*alphaDrSurf[8]*fUpwind[22]+0.1767766952966368*alphaDrSurf[7]*fUpwind[21]+0.1767766952966368*alphaDrSurf[6]*fUpwind[20]+0.1767766952966368*alphaDrSurf[3]*fUpwind[14]+0.1767766952966368*alphaDrSurf[2]*fUpwind[13]+0.1767766952966368*alphaDrSurf[1]*fUpwind[12]+0.1767766952966368*alphaDrSurf[0]*fUpwind[5]; 
  drag_incr[6] = 0.1767766952966368*alphaDrSurf[3]*fUpwind[16]+0.1767766952966368*fUpwind[3]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[7]*fUpwind[8]+0.1767766952966368*fUpwind[7]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[0]*fUpwind[6]+0.1767766952966368*fUpwind[0]*alphaDrSurf[6]+0.1767766952966368*alphaDrSurf[1]*fUpwind[2]+0.1767766952966368*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[7] = 0.1767766952966368*alphaDrSurf[2]*fUpwind[16]+0.1767766952966368*fUpwind[2]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[6]*fUpwind[8]+0.1767766952966368*fUpwind[6]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[0]*fUpwind[7]+0.1767766952966368*fUpwind[0]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[1]*fUpwind[3]+0.1767766952966368*fUpwind[1]*alphaDrSurf[3]; 
  drag_incr[8] = 0.1767766952966368*alphaDrSurf[1]*fUpwind[16]+0.1767766952966368*fUpwind[1]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[0]*fUpwind[8]+0.1767766952966368*fUpwind[0]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[6]*fUpwind[7]+0.1767766952966368*fUpwind[6]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[2]*fUpwind[3]+0.1767766952966368*fUpwind[2]*alphaDrSurf[3]; 
  drag_incr[9] = 0.1767766952966368*alphaDrSurf[8]*fUpwind[26]+0.1767766952966368*alphaDrSurf[16]*fUpwind[19]+0.1767766952966368*alphaDrSurf[3]*fUpwind[18]+0.1767766952966368*alphaDrSurf[2]*fUpwind[17]+0.1767766952966368*alphaDrSurf[7]*fUpwind[11]+0.1767766952966368*alphaDrSurf[6]*fUpwind[10]+0.1767766952966368*alphaDrSurf[0]*fUpwind[9]+0.1767766952966368*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[10] = 0.1767766952966368*alphaDrSurf[7]*fUpwind[26]+0.1767766952966368*alphaDrSurf[3]*fUpwind[19]+0.1767766952966368*alphaDrSurf[16]*fUpwind[18]+0.1767766952966368*alphaDrSurf[1]*fUpwind[17]+0.1767766952966368*alphaDrSurf[8]*fUpwind[11]+0.1767766952966368*alphaDrSurf[0]*fUpwind[10]+0.1767766952966368*alphaDrSurf[6]*fUpwind[9]+0.1767766952966368*alphaDrSurf[2]*fUpwind[4]; 
  drag_incr[11] = 0.1767766952966368*alphaDrSurf[6]*fUpwind[26]+0.1767766952966368*alphaDrSurf[2]*fUpwind[19]+0.1767766952966368*alphaDrSurf[1]*fUpwind[18]+0.1767766952966368*alphaDrSurf[16]*fUpwind[17]+0.1767766952966368*alphaDrSurf[0]*fUpwind[11]+0.1767766952966368*alphaDrSurf[8]*fUpwind[10]+0.1767766952966368*alphaDrSurf[7]*fUpwind[9]+0.1767766952966368*alphaDrSurf[3]*fUpwind[4]; 
  drag_incr[12] = 0.1767766952966368*alphaDrSurf[8]*fUpwind[27]+0.1767766952966368*alphaDrSurf[16]*fUpwind[22]+0.1767766952966368*alphaDrSurf[3]*fUpwind[21]+0.1767766952966368*alphaDrSurf[2]*fUpwind[20]+0.1767766952966368*alphaDrSurf[7]*fUpwind[14]+0.1767766952966368*alphaDrSurf[6]*fUpwind[13]+0.1767766952966368*alphaDrSurf[0]*fUpwind[12]+0.1767766952966368*alphaDrSurf[1]*fUpwind[5]; 
  drag_incr[13] = 0.1767766952966368*alphaDrSurf[7]*fUpwind[27]+0.1767766952966368*alphaDrSurf[3]*fUpwind[22]+0.1767766952966368*alphaDrSurf[16]*fUpwind[21]+0.1767766952966368*alphaDrSurf[1]*fUpwind[20]+0.1767766952966368*alphaDrSurf[8]*fUpwind[14]+0.1767766952966368*alphaDrSurf[0]*fUpwind[13]+0.1767766952966368*alphaDrSurf[6]*fUpwind[12]+0.1767766952966368*alphaDrSurf[2]*fUpwind[5]; 
  drag_incr[14] = 0.1767766952966368*alphaDrSurf[6]*fUpwind[27]+0.1767766952966368*alphaDrSurf[2]*fUpwind[22]+0.1767766952966368*alphaDrSurf[1]*fUpwind[21]+0.1767766952966368*alphaDrSurf[16]*fUpwind[20]+0.1767766952966368*alphaDrSurf[0]*fUpwind[14]+0.1767766952966368*alphaDrSurf[8]*fUpwind[13]+0.1767766952966368*alphaDrSurf[7]*fUpwind[12]+0.1767766952966368*alphaDrSurf[3]*fUpwind[5]; 
  drag_incr[15] = 0.1767766952966368*alphaDrSurf[16]*fUpwind[31]+0.1767766952966368*alphaDrSurf[8]*fUpwind[30]+0.1767766952966368*alphaDrSurf[7]*fUpwind[29]+0.1767766952966368*alphaDrSurf[6]*fUpwind[28]+0.1767766952966368*alphaDrSurf[3]*fUpwind[25]+0.1767766952966368*alphaDrSurf[2]*fUpwind[24]+0.1767766952966368*alphaDrSurf[1]*fUpwind[23]+0.1767766952966368*alphaDrSurf[0]*fUpwind[15]; 
  drag_incr[16] = 0.1767766952966368*alphaDrSurf[0]*fUpwind[16]+0.1767766952966368*fUpwind[0]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[8]+0.1767766952966368*fUpwind[1]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[2]*fUpwind[7]+0.1767766952966368*fUpwind[2]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[3]*fUpwind[6]+0.1767766952966368*fUpwind[3]*alphaDrSurf[6]; 
  drag_incr[17] = 0.1767766952966368*alphaDrSurf[3]*fUpwind[26]+0.1767766952966368*alphaDrSurf[7]*fUpwind[19]+0.1767766952966368*alphaDrSurf[8]*fUpwind[18]+0.1767766952966368*alphaDrSurf[0]*fUpwind[17]+0.1767766952966368*fUpwind[11]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[10]+0.1767766952966368*alphaDrSurf[2]*fUpwind[9]+0.1767766952966368*fUpwind[4]*alphaDrSurf[6]; 
  drag_incr[18] = 0.1767766952966368*alphaDrSurf[2]*fUpwind[26]+0.1767766952966368*alphaDrSurf[6]*fUpwind[19]+0.1767766952966368*alphaDrSurf[0]*fUpwind[18]+0.1767766952966368*alphaDrSurf[8]*fUpwind[17]+0.1767766952966368*fUpwind[10]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[11]+0.1767766952966368*alphaDrSurf[3]*fUpwind[9]+0.1767766952966368*fUpwind[4]*alphaDrSurf[7]; 
  drag_incr[19] = 0.1767766952966368*alphaDrSurf[1]*fUpwind[26]+0.1767766952966368*alphaDrSurf[0]*fUpwind[19]+0.1767766952966368*alphaDrSurf[6]*fUpwind[18]+0.1767766952966368*alphaDrSurf[7]*fUpwind[17]+0.1767766952966368*fUpwind[9]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[2]*fUpwind[11]+0.1767766952966368*alphaDrSurf[3]*fUpwind[10]+0.1767766952966368*fUpwind[4]*alphaDrSurf[8]; 
  drag_incr[20] = 0.1767766952966368*alphaDrSurf[3]*fUpwind[27]+0.1767766952966368*alphaDrSurf[7]*fUpwind[22]+0.1767766952966368*alphaDrSurf[8]*fUpwind[21]+0.1767766952966368*alphaDrSurf[0]*fUpwind[20]+0.1767766952966368*fUpwind[14]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[13]+0.1767766952966368*alphaDrSurf[2]*fUpwind[12]+0.1767766952966368*fUpwind[5]*alphaDrSurf[6]; 
  drag_incr[21] = 0.1767766952966368*alphaDrSurf[2]*fUpwind[27]+0.1767766952966368*alphaDrSurf[6]*fUpwind[22]+0.1767766952966368*alphaDrSurf[0]*fUpwind[21]+0.1767766952966368*alphaDrSurf[8]*fUpwind[20]+0.1767766952966368*fUpwind[13]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[14]+0.1767766952966368*alphaDrSurf[3]*fUpwind[12]+0.1767766952966368*fUpwind[5]*alphaDrSurf[7]; 
  drag_incr[22] = 0.1767766952966368*alphaDrSurf[1]*fUpwind[27]+0.1767766952966368*alphaDrSurf[0]*fUpwind[22]+0.1767766952966368*alphaDrSurf[6]*fUpwind[21]+0.1767766952966368*alphaDrSurf[7]*fUpwind[20]+0.1767766952966368*fUpwind[12]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[2]*fUpwind[14]+0.1767766952966368*alphaDrSurf[3]*fUpwind[13]+0.1767766952966368*fUpwind[5]*alphaDrSurf[8]; 
  drag_incr[23] = 0.1767766952966368*alphaDrSurf[8]*fUpwind[31]+0.1767766952966368*alphaDrSurf[16]*fUpwind[30]+0.1767766952966368*alphaDrSurf[3]*fUpwind[29]+0.1767766952966368*alphaDrSurf[2]*fUpwind[28]+0.1767766952966368*alphaDrSurf[7]*fUpwind[25]+0.1767766952966368*alphaDrSurf[6]*fUpwind[24]+0.1767766952966368*alphaDrSurf[0]*fUpwind[23]+0.1767766952966368*alphaDrSurf[1]*fUpwind[15]; 
  drag_incr[24] = 0.1767766952966368*alphaDrSurf[7]*fUpwind[31]+0.1767766952966368*alphaDrSurf[3]*fUpwind[30]+0.1767766952966368*alphaDrSurf[16]*fUpwind[29]+0.1767766952966368*alphaDrSurf[1]*fUpwind[28]+0.1767766952966368*alphaDrSurf[8]*fUpwind[25]+0.1767766952966368*alphaDrSurf[0]*fUpwind[24]+0.1767766952966368*alphaDrSurf[6]*fUpwind[23]+0.1767766952966368*alphaDrSurf[2]*fUpwind[15]; 
  drag_incr[25] = 0.1767766952966368*alphaDrSurf[6]*fUpwind[31]+0.1767766952966368*alphaDrSurf[2]*fUpwind[30]+0.1767766952966368*alphaDrSurf[1]*fUpwind[29]+0.1767766952966368*alphaDrSurf[16]*fUpwind[28]+0.1767766952966368*alphaDrSurf[0]*fUpwind[25]+0.1767766952966368*alphaDrSurf[8]*fUpwind[24]+0.1767766952966368*alphaDrSurf[7]*fUpwind[23]+0.1767766952966368*alphaDrSurf[3]*fUpwind[15]; 
  drag_incr[26] = 0.1767766952966368*alphaDrSurf[0]*fUpwind[26]+0.1767766952966368*alphaDrSurf[1]*fUpwind[19]+0.1767766952966368*alphaDrSurf[2]*fUpwind[18]+0.1767766952966368*alphaDrSurf[3]*fUpwind[17]+0.1767766952966368*fUpwind[4]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[6]*fUpwind[11]+0.1767766952966368*alphaDrSurf[7]*fUpwind[10]+0.1767766952966368*alphaDrSurf[8]*fUpwind[9]; 
  drag_incr[27] = 0.1767766952966368*alphaDrSurf[0]*fUpwind[27]+0.1767766952966368*alphaDrSurf[1]*fUpwind[22]+0.1767766952966368*alphaDrSurf[2]*fUpwind[21]+0.1767766952966368*alphaDrSurf[3]*fUpwind[20]+0.1767766952966368*fUpwind[5]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[6]*fUpwind[14]+0.1767766952966368*alphaDrSurf[7]*fUpwind[13]+0.1767766952966368*alphaDrSurf[8]*fUpwind[12]; 
  drag_incr[28] = 0.1767766952966368*alphaDrSurf[3]*fUpwind[31]+0.1767766952966368*alphaDrSurf[7]*fUpwind[30]+0.1767766952966368*alphaDrSurf[8]*fUpwind[29]+0.1767766952966368*alphaDrSurf[0]*fUpwind[28]+0.1767766952966368*alphaDrSurf[16]*fUpwind[25]+0.1767766952966368*alphaDrSurf[1]*fUpwind[24]+0.1767766952966368*alphaDrSurf[2]*fUpwind[23]+0.1767766952966368*alphaDrSurf[6]*fUpwind[15]; 
  drag_incr[29] = 0.1767766952966368*alphaDrSurf[2]*fUpwind[31]+0.1767766952966368*alphaDrSurf[6]*fUpwind[30]+0.1767766952966368*alphaDrSurf[0]*fUpwind[29]+0.1767766952966368*alphaDrSurf[8]*fUpwind[28]+0.1767766952966368*alphaDrSurf[1]*fUpwind[25]+0.1767766952966368*alphaDrSurf[16]*fUpwind[24]+0.1767766952966368*alphaDrSurf[3]*fUpwind[23]+0.1767766952966368*alphaDrSurf[7]*fUpwind[15]; 
  drag_incr[30] = 0.1767766952966368*alphaDrSurf[1]*fUpwind[31]+0.1767766952966368*alphaDrSurf[0]*fUpwind[30]+0.1767766952966368*alphaDrSurf[6]*fUpwind[29]+0.1767766952966368*alphaDrSurf[7]*fUpwind[28]+0.1767766952966368*alphaDrSurf[2]*fUpwind[25]+0.1767766952966368*alphaDrSurf[3]*fUpwind[24]+0.1767766952966368*alphaDrSurf[16]*fUpwind[23]+0.1767766952966368*alphaDrSurf[8]*fUpwind[15]; 
  drag_incr[31] = 0.1767766952966368*alphaDrSurf[0]*fUpwind[31]+0.1767766952966368*alphaDrSurf[1]*fUpwind[30]+0.1767766952966368*alphaDrSurf[2]*fUpwind[29]+0.1767766952966368*alphaDrSurf[3]*fUpwind[28]+0.1767766952966368*alphaDrSurf[6]*fUpwind[25]+0.1767766952966368*alphaDrSurf[7]*fUpwind[24]+0.1767766952966368*alphaDrSurf[8]*fUpwind[23]+0.1767766952966368*fUpwind[15]*alphaDrSurf[16]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[9] += 0.7071067811865475*drag_incr[8]*rdv2; 
  out[10] += 0.7071067811865475*drag_incr[9]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[10]*rdv2; 
  out[12] += 0.7071067811865475*drag_incr[11]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[16] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[17] += 0.7071067811865475*drag_incr[12]*rdv2; 
  out[18] += 0.7071067811865475*drag_incr[13]*rdv2; 
  out[19] += 0.7071067811865475*drag_incr[14]*rdv2; 
  out[20] += 0.7071067811865475*drag_incr[15]*rdv2; 
  out[21] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[22] += 0.7071067811865475*drag_incr[16]*rdv2; 
  out[23] += 0.7071067811865475*drag_incr[17]*rdv2; 
  out[24] += 0.7071067811865475*drag_incr[18]*rdv2; 
  out[25] += 0.7071067811865475*drag_incr[19]*rdv2; 
  out[26] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[27] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[28] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[32] += 0.7071067811865475*drag_incr[20]*rdv2; 
  out[33] += 0.7071067811865475*drag_incr[21]*rdv2; 
  out[34] += 0.7071067811865475*drag_incr[22]*rdv2; 
  out[35] += 0.7071067811865475*drag_incr[23]*rdv2; 
  out[36] += 0.7071067811865475*drag_incr[24]*rdv2; 
  out[37] += 0.7071067811865475*drag_incr[25]*rdv2; 
  out[38] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[39] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[40] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[41] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[42] += 0.7071067811865475*drag_incr[26]*rdv2; 
  out[43] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[44] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[45] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[46] += 1.224744871391589*drag_incr[19]*rdv2; 
  out[47] += 0.7071067811865475*drag_incr[27]*rdv2; 
  out[48] += 0.7071067811865475*drag_incr[28]*rdv2; 
  out[49] += 0.7071067811865475*drag_incr[29]*rdv2; 
  out[50] += 0.7071067811865475*drag_incr[30]*rdv2; 
  out[51] += 1.224744871391589*drag_incr[20]*rdv2; 
  out[52] += 1.224744871391589*drag_incr[21]*rdv2; 
  out[53] += 1.224744871391589*drag_incr[22]*rdv2; 
  out[54] += 1.224744871391589*drag_incr[23]*rdv2; 
  out[55] += 1.224744871391589*drag_incr[24]*rdv2; 
  out[56] += 1.224744871391589*drag_incr[25]*rdv2; 
  out[57] += 1.224744871391589*drag_incr[26]*rdv2; 
  out[58] += 0.7071067811865475*drag_incr[31]*rdv2; 
  out[59] += 1.224744871391589*drag_incr[27]*rdv2; 
  out[60] += 1.224744871391589*drag_incr[28]*rdv2; 
  out[61] += 1.224744871391589*drag_incr[29]*rdv2; 
  out[62] += 1.224744871391589*drag_incr[30]*rdv2; 
  out[63] += 1.224744871391589*drag_incr[31]*rdv2; 

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
    fUpwindQuad[0] = ser_3x3v_p1_surfvy_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_3x3v_p1_surfvy_quad_0(-1, fSkin); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_3x3v_p1_surfvy_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_3x3v_p1_surfvy_quad_1(-1, fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_3x3v_p1_surfvy_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_3x3v_p1_surfvy_quad_2(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_3x3v_p1_surfvy_quad_3(1, fEdge); 
  } else { 
    fUpwindQuad[3] = ser_3x3v_p1_surfvy_quad_3(-1, fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_3x3v_p1_surfvy_quad_4(1, fEdge); 
  } else { 
    fUpwindQuad[4] = ser_3x3v_p1_surfvy_quad_4(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_3x3v_p1_surfvy_quad_5(1, fEdge); 
  } else { 
    fUpwindQuad[5] = ser_3x3v_p1_surfvy_quad_5(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_3x3v_p1_surfvy_quad_6(1, fEdge); 
  } else { 
    fUpwindQuad[6] = ser_3x3v_p1_surfvy_quad_6(-1, fSkin); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_3x3v_p1_surfvy_quad_7(1, fEdge); 
  } else { 
    fUpwindQuad[7] = ser_3x3v_p1_surfvy_quad_7(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_3x3v_p1_surfvy_quad_8(1, fEdge); 
  } else { 
    fUpwindQuad[8] = ser_3x3v_p1_surfvy_quad_8(-1, fSkin); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_3x3v_p1_surfvy_quad_9(1, fEdge); 
  } else { 
    fUpwindQuad[9] = ser_3x3v_p1_surfvy_quad_9(-1, fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = ser_3x3v_p1_surfvy_quad_10(1, fEdge); 
  } else { 
    fUpwindQuad[10] = ser_3x3v_p1_surfvy_quad_10(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_3x3v_p1_surfvy_quad_11(1, fEdge); 
  } else { 
    fUpwindQuad[11] = ser_3x3v_p1_surfvy_quad_11(-1, fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_3x3v_p1_surfvy_quad_12(1, fEdge); 
  } else { 
    fUpwindQuad[12] = ser_3x3v_p1_surfvy_quad_12(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[13] = ser_3x3v_p1_surfvy_quad_13(1, fEdge); 
  } else { 
    fUpwindQuad[13] = ser_3x3v_p1_surfvy_quad_13(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_3x3v_p1_surfvy_quad_14(1, fEdge); 
  } else { 
    fUpwindQuad[14] = ser_3x3v_p1_surfvy_quad_14(-1, fSkin); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_3x3v_p1_surfvy_quad_15(1, fEdge); 
  } else { 
    fUpwindQuad[15] = ser_3x3v_p1_surfvy_quad_15(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[16] = ser_3x3v_p1_surfvy_quad_16(1, fEdge); 
  } else { 
    fUpwindQuad[16] = ser_3x3v_p1_surfvy_quad_16(-1, fSkin); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[17] = ser_3x3v_p1_surfvy_quad_17(1, fEdge); 
  } else { 
    fUpwindQuad[17] = ser_3x3v_p1_surfvy_quad_17(-1, fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_3x3v_p1_surfvy_quad_18(1, fEdge); 
  } else { 
    fUpwindQuad[18] = ser_3x3v_p1_surfvy_quad_18(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[19] = ser_3x3v_p1_surfvy_quad_19(1, fEdge); 
  } else { 
    fUpwindQuad[19] = ser_3x3v_p1_surfvy_quad_19(-1, fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[20] = ser_3x3v_p1_surfvy_quad_20(1, fEdge); 
  } else { 
    fUpwindQuad[20] = ser_3x3v_p1_surfvy_quad_20(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_3x3v_p1_surfvy_quad_21(1, fEdge); 
  } else { 
    fUpwindQuad[21] = ser_3x3v_p1_surfvy_quad_21(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[22] = ser_3x3v_p1_surfvy_quad_22(1, fEdge); 
  } else { 
    fUpwindQuad[22] = ser_3x3v_p1_surfvy_quad_22(-1, fSkin); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[23] = ser_3x3v_p1_surfvy_quad_23(1, fEdge); 
  } else { 
    fUpwindQuad[23] = ser_3x3v_p1_surfvy_quad_23(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_3x3v_p1_surfvy_quad_24(1, fEdge); 
  } else { 
    fUpwindQuad[24] = ser_3x3v_p1_surfvy_quad_24(-1, fSkin); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[25] = ser_3x3v_p1_surfvy_quad_25(1, fEdge); 
  } else { 
    fUpwindQuad[25] = ser_3x3v_p1_surfvy_quad_25(-1, fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[26] = ser_3x3v_p1_surfvy_quad_26(1, fEdge); 
  } else { 
    fUpwindQuad[26] = ser_3x3v_p1_surfvy_quad_26(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]-alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[27] = ser_3x3v_p1_surfvy_quad_27(1, fEdge); 
  } else { 
    fUpwindQuad[27] = ser_3x3v_p1_surfvy_quad_27(-1, fSkin); 
  } 
  if (alphaDrSurf[16]-alphaDrSurf[8]-alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[28] = ser_3x3v_p1_surfvy_quad_28(1, fEdge); 
  } else { 
    fUpwindQuad[28] = ser_3x3v_p1_surfvy_quad_28(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])-alphaDrSurf[8]+alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[29] = ser_3x3v_p1_surfvy_quad_29(1, fEdge); 
  } else { 
    fUpwindQuad[29] = ser_3x3v_p1_surfvy_quad_29(-1, fSkin); 
  } 
  if ((-alphaDrSurf[16])+alphaDrSurf[8]-alphaDrSurf[7]-alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[30] = ser_3x3v_p1_surfvy_quad_30(1, fEdge); 
  } else { 
    fUpwindQuad[30] = ser_3x3v_p1_surfvy_quad_30(-1, fSkin); 
  } 
  if (alphaDrSurf[16]+alphaDrSurf[8]+alphaDrSurf[7]+alphaDrSurf[6]+alphaDrSurf[3]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[31] = ser_3x3v_p1_surfvy_quad_31(1, fEdge); 
  } else { 
    fUpwindQuad[31] = ser_3x3v_p1_surfvy_quad_31(-1, fSkin); 
  } 

  fUpwind[0] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[6] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[7] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[12] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[13] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[14] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24])+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[15] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[16] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[17] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[18] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*fUpwindQuad[16]+fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[19] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16])+fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[20] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[21] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]-1.0*fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[22] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]-1.0*(fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[23] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*fUpwindQuad[28]+fUpwindQuad[27]-1.0*fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]-1.0*fUpwindQuad[17]+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[24] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28])+fUpwindQuad[27]+fUpwindQuad[26]-1.0*(fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]-1.0*(fUpwindQuad[19]+fUpwindQuad[18])+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[25] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]+fUpwindQuad[29]+fUpwindQuad[28]-1.0*(fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]+fUpwindQuad[24]+fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]+fUpwindQuad[20])+fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]+fUpwindQuad[16]-1.0*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[26] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]+fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[27] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*fUpwindQuad[24]+fUpwindQuad[23]-1.0*(fUpwindQuad[22]+fUpwindQuad[21])+fUpwindQuad[20]-1.0*fUpwindQuad[19]+fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[28] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]+fUpwindQuad[27]-1.0*(fUpwindQuad[26]+fUpwindQuad[25])+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]+fUpwindQuad[21]-1.0*(fUpwindQuad[20]+fUpwindQuad[19])+fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[29] = 0.1767766952966368*(fUpwindQuad[31]-1.0*fUpwindQuad[30]+fUpwindQuad[29]-1.0*(fUpwindQuad[28]+fUpwindQuad[27])+fUpwindQuad[26]-1.0*fUpwindQuad[25]+fUpwindQuad[24]-1.0*fUpwindQuad[23]+fUpwindQuad[22]-1.0*fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]-1.0*fUpwindQuad[18]+fUpwindQuad[17]-1.0*(fUpwindQuad[16]+fUpwindQuad[15])+fUpwindQuad[14]-1.0*fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[30] = 0.1767766952966368*(fUpwindQuad[31]+fUpwindQuad[30]-1.0*(fUpwindQuad[29]+fUpwindQuad[28]+fUpwindQuad[27]+fUpwindQuad[26])+fUpwindQuad[25]+fUpwindQuad[24]-1.0*(fUpwindQuad[23]+fUpwindQuad[22])+fUpwindQuad[21]+fUpwindQuad[20]+fUpwindQuad[19]+fUpwindQuad[18]-1.0*(fUpwindQuad[17]+fUpwindQuad[16]+fUpwindQuad[15]+fUpwindQuad[14])+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[31] = 0.1767766952966368*(fUpwindQuad[31]-1.0*(fUpwindQuad[30]+fUpwindQuad[29])+fUpwindQuad[28]-1.0*fUpwindQuad[27]+fUpwindQuad[26]+fUpwindQuad[25]-1.0*(fUpwindQuad[24]+fUpwindQuad[23])+fUpwindQuad[22]+fUpwindQuad[21]-1.0*fUpwindQuad[20]+fUpwindQuad[19]-1.0*(fUpwindQuad[18]+fUpwindQuad[17])+fUpwindQuad[16]-1.0*fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  drag_incr[0] = 0.1767766952966368*alphaDrSurf[16]*fUpwind[16]+0.1767766952966368*alphaDrSurf[8]*fUpwind[8]+0.1767766952966368*alphaDrSurf[7]*fUpwind[7]+0.1767766952966368*alphaDrSurf[6]*fUpwind[6]+0.1767766952966368*alphaDrSurf[3]*fUpwind[3]+0.1767766952966368*alphaDrSurf[2]*fUpwind[2]+0.1767766952966368*alphaDrSurf[1]*fUpwind[1]+0.1767766952966368*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.1767766952966368*alphaDrSurf[8]*fUpwind[16]+0.1767766952966368*fUpwind[8]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[3]*fUpwind[7]+0.1767766952966368*fUpwind[3]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[2]*fUpwind[6]+0.1767766952966368*fUpwind[2]*alphaDrSurf[6]+0.1767766952966368*alphaDrSurf[0]*fUpwind[1]+0.1767766952966368*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.1767766952966368*alphaDrSurf[7]*fUpwind[16]+0.1767766952966368*fUpwind[7]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[3]*fUpwind[8]+0.1767766952966368*fUpwind[3]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[1]*fUpwind[6]+0.1767766952966368*fUpwind[1]*alphaDrSurf[6]+0.1767766952966368*alphaDrSurf[0]*fUpwind[2]+0.1767766952966368*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.1767766952966368*alphaDrSurf[6]*fUpwind[16]+0.1767766952966368*fUpwind[6]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[2]*fUpwind[8]+0.1767766952966368*fUpwind[2]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[1]*fUpwind[7]+0.1767766952966368*fUpwind[1]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[0]*fUpwind[3]+0.1767766952966368*fUpwind[0]*alphaDrSurf[3]; 
  drag_incr[4] = 0.1767766952966368*alphaDrSurf[16]*fUpwind[26]+0.1767766952966368*alphaDrSurf[8]*fUpwind[19]+0.1767766952966368*alphaDrSurf[7]*fUpwind[18]+0.1767766952966368*alphaDrSurf[6]*fUpwind[17]+0.1767766952966368*alphaDrSurf[3]*fUpwind[11]+0.1767766952966368*alphaDrSurf[2]*fUpwind[10]+0.1767766952966368*alphaDrSurf[1]*fUpwind[9]+0.1767766952966368*alphaDrSurf[0]*fUpwind[4]; 
  drag_incr[5] = 0.1767766952966368*alphaDrSurf[16]*fUpwind[27]+0.1767766952966368*alphaDrSurf[8]*fUpwind[22]+0.1767766952966368*alphaDrSurf[7]*fUpwind[21]+0.1767766952966368*alphaDrSurf[6]*fUpwind[20]+0.1767766952966368*alphaDrSurf[3]*fUpwind[14]+0.1767766952966368*alphaDrSurf[2]*fUpwind[13]+0.1767766952966368*alphaDrSurf[1]*fUpwind[12]+0.1767766952966368*alphaDrSurf[0]*fUpwind[5]; 
  drag_incr[6] = 0.1767766952966368*alphaDrSurf[3]*fUpwind[16]+0.1767766952966368*fUpwind[3]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[7]*fUpwind[8]+0.1767766952966368*fUpwind[7]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[0]*fUpwind[6]+0.1767766952966368*fUpwind[0]*alphaDrSurf[6]+0.1767766952966368*alphaDrSurf[1]*fUpwind[2]+0.1767766952966368*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[7] = 0.1767766952966368*alphaDrSurf[2]*fUpwind[16]+0.1767766952966368*fUpwind[2]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[6]*fUpwind[8]+0.1767766952966368*fUpwind[6]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[0]*fUpwind[7]+0.1767766952966368*fUpwind[0]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[1]*fUpwind[3]+0.1767766952966368*fUpwind[1]*alphaDrSurf[3]; 
  drag_incr[8] = 0.1767766952966368*alphaDrSurf[1]*fUpwind[16]+0.1767766952966368*fUpwind[1]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[0]*fUpwind[8]+0.1767766952966368*fUpwind[0]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[6]*fUpwind[7]+0.1767766952966368*fUpwind[6]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[2]*fUpwind[3]+0.1767766952966368*fUpwind[2]*alphaDrSurf[3]; 
  drag_incr[9] = 0.1767766952966368*alphaDrSurf[8]*fUpwind[26]+0.1767766952966368*alphaDrSurf[16]*fUpwind[19]+0.1767766952966368*alphaDrSurf[3]*fUpwind[18]+0.1767766952966368*alphaDrSurf[2]*fUpwind[17]+0.1767766952966368*alphaDrSurf[7]*fUpwind[11]+0.1767766952966368*alphaDrSurf[6]*fUpwind[10]+0.1767766952966368*alphaDrSurf[0]*fUpwind[9]+0.1767766952966368*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[10] = 0.1767766952966368*alphaDrSurf[7]*fUpwind[26]+0.1767766952966368*alphaDrSurf[3]*fUpwind[19]+0.1767766952966368*alphaDrSurf[16]*fUpwind[18]+0.1767766952966368*alphaDrSurf[1]*fUpwind[17]+0.1767766952966368*alphaDrSurf[8]*fUpwind[11]+0.1767766952966368*alphaDrSurf[0]*fUpwind[10]+0.1767766952966368*alphaDrSurf[6]*fUpwind[9]+0.1767766952966368*alphaDrSurf[2]*fUpwind[4]; 
  drag_incr[11] = 0.1767766952966368*alphaDrSurf[6]*fUpwind[26]+0.1767766952966368*alphaDrSurf[2]*fUpwind[19]+0.1767766952966368*alphaDrSurf[1]*fUpwind[18]+0.1767766952966368*alphaDrSurf[16]*fUpwind[17]+0.1767766952966368*alphaDrSurf[0]*fUpwind[11]+0.1767766952966368*alphaDrSurf[8]*fUpwind[10]+0.1767766952966368*alphaDrSurf[7]*fUpwind[9]+0.1767766952966368*alphaDrSurf[3]*fUpwind[4]; 
  drag_incr[12] = 0.1767766952966368*alphaDrSurf[8]*fUpwind[27]+0.1767766952966368*alphaDrSurf[16]*fUpwind[22]+0.1767766952966368*alphaDrSurf[3]*fUpwind[21]+0.1767766952966368*alphaDrSurf[2]*fUpwind[20]+0.1767766952966368*alphaDrSurf[7]*fUpwind[14]+0.1767766952966368*alphaDrSurf[6]*fUpwind[13]+0.1767766952966368*alphaDrSurf[0]*fUpwind[12]+0.1767766952966368*alphaDrSurf[1]*fUpwind[5]; 
  drag_incr[13] = 0.1767766952966368*alphaDrSurf[7]*fUpwind[27]+0.1767766952966368*alphaDrSurf[3]*fUpwind[22]+0.1767766952966368*alphaDrSurf[16]*fUpwind[21]+0.1767766952966368*alphaDrSurf[1]*fUpwind[20]+0.1767766952966368*alphaDrSurf[8]*fUpwind[14]+0.1767766952966368*alphaDrSurf[0]*fUpwind[13]+0.1767766952966368*alphaDrSurf[6]*fUpwind[12]+0.1767766952966368*alphaDrSurf[2]*fUpwind[5]; 
  drag_incr[14] = 0.1767766952966368*alphaDrSurf[6]*fUpwind[27]+0.1767766952966368*alphaDrSurf[2]*fUpwind[22]+0.1767766952966368*alphaDrSurf[1]*fUpwind[21]+0.1767766952966368*alphaDrSurf[16]*fUpwind[20]+0.1767766952966368*alphaDrSurf[0]*fUpwind[14]+0.1767766952966368*alphaDrSurf[8]*fUpwind[13]+0.1767766952966368*alphaDrSurf[7]*fUpwind[12]+0.1767766952966368*alphaDrSurf[3]*fUpwind[5]; 
  drag_incr[15] = 0.1767766952966368*alphaDrSurf[16]*fUpwind[31]+0.1767766952966368*alphaDrSurf[8]*fUpwind[30]+0.1767766952966368*alphaDrSurf[7]*fUpwind[29]+0.1767766952966368*alphaDrSurf[6]*fUpwind[28]+0.1767766952966368*alphaDrSurf[3]*fUpwind[25]+0.1767766952966368*alphaDrSurf[2]*fUpwind[24]+0.1767766952966368*alphaDrSurf[1]*fUpwind[23]+0.1767766952966368*alphaDrSurf[0]*fUpwind[15]; 
  drag_incr[16] = 0.1767766952966368*alphaDrSurf[0]*fUpwind[16]+0.1767766952966368*fUpwind[0]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[8]+0.1767766952966368*fUpwind[1]*alphaDrSurf[8]+0.1767766952966368*alphaDrSurf[2]*fUpwind[7]+0.1767766952966368*fUpwind[2]*alphaDrSurf[7]+0.1767766952966368*alphaDrSurf[3]*fUpwind[6]+0.1767766952966368*fUpwind[3]*alphaDrSurf[6]; 
  drag_incr[17] = 0.1767766952966368*alphaDrSurf[3]*fUpwind[26]+0.1767766952966368*alphaDrSurf[7]*fUpwind[19]+0.1767766952966368*alphaDrSurf[8]*fUpwind[18]+0.1767766952966368*alphaDrSurf[0]*fUpwind[17]+0.1767766952966368*fUpwind[11]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[10]+0.1767766952966368*alphaDrSurf[2]*fUpwind[9]+0.1767766952966368*fUpwind[4]*alphaDrSurf[6]; 
  drag_incr[18] = 0.1767766952966368*alphaDrSurf[2]*fUpwind[26]+0.1767766952966368*alphaDrSurf[6]*fUpwind[19]+0.1767766952966368*alphaDrSurf[0]*fUpwind[18]+0.1767766952966368*alphaDrSurf[8]*fUpwind[17]+0.1767766952966368*fUpwind[10]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[11]+0.1767766952966368*alphaDrSurf[3]*fUpwind[9]+0.1767766952966368*fUpwind[4]*alphaDrSurf[7]; 
  drag_incr[19] = 0.1767766952966368*alphaDrSurf[1]*fUpwind[26]+0.1767766952966368*alphaDrSurf[0]*fUpwind[19]+0.1767766952966368*alphaDrSurf[6]*fUpwind[18]+0.1767766952966368*alphaDrSurf[7]*fUpwind[17]+0.1767766952966368*fUpwind[9]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[2]*fUpwind[11]+0.1767766952966368*alphaDrSurf[3]*fUpwind[10]+0.1767766952966368*fUpwind[4]*alphaDrSurf[8]; 
  drag_incr[20] = 0.1767766952966368*alphaDrSurf[3]*fUpwind[27]+0.1767766952966368*alphaDrSurf[7]*fUpwind[22]+0.1767766952966368*alphaDrSurf[8]*fUpwind[21]+0.1767766952966368*alphaDrSurf[0]*fUpwind[20]+0.1767766952966368*fUpwind[14]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[13]+0.1767766952966368*alphaDrSurf[2]*fUpwind[12]+0.1767766952966368*fUpwind[5]*alphaDrSurf[6]; 
  drag_incr[21] = 0.1767766952966368*alphaDrSurf[2]*fUpwind[27]+0.1767766952966368*alphaDrSurf[6]*fUpwind[22]+0.1767766952966368*alphaDrSurf[0]*fUpwind[21]+0.1767766952966368*alphaDrSurf[8]*fUpwind[20]+0.1767766952966368*fUpwind[13]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[1]*fUpwind[14]+0.1767766952966368*alphaDrSurf[3]*fUpwind[12]+0.1767766952966368*fUpwind[5]*alphaDrSurf[7]; 
  drag_incr[22] = 0.1767766952966368*alphaDrSurf[1]*fUpwind[27]+0.1767766952966368*alphaDrSurf[0]*fUpwind[22]+0.1767766952966368*alphaDrSurf[6]*fUpwind[21]+0.1767766952966368*alphaDrSurf[7]*fUpwind[20]+0.1767766952966368*fUpwind[12]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[2]*fUpwind[14]+0.1767766952966368*alphaDrSurf[3]*fUpwind[13]+0.1767766952966368*fUpwind[5]*alphaDrSurf[8]; 
  drag_incr[23] = 0.1767766952966368*alphaDrSurf[8]*fUpwind[31]+0.1767766952966368*alphaDrSurf[16]*fUpwind[30]+0.1767766952966368*alphaDrSurf[3]*fUpwind[29]+0.1767766952966368*alphaDrSurf[2]*fUpwind[28]+0.1767766952966368*alphaDrSurf[7]*fUpwind[25]+0.1767766952966368*alphaDrSurf[6]*fUpwind[24]+0.1767766952966368*alphaDrSurf[0]*fUpwind[23]+0.1767766952966368*alphaDrSurf[1]*fUpwind[15]; 
  drag_incr[24] = 0.1767766952966368*alphaDrSurf[7]*fUpwind[31]+0.1767766952966368*alphaDrSurf[3]*fUpwind[30]+0.1767766952966368*alphaDrSurf[16]*fUpwind[29]+0.1767766952966368*alphaDrSurf[1]*fUpwind[28]+0.1767766952966368*alphaDrSurf[8]*fUpwind[25]+0.1767766952966368*alphaDrSurf[0]*fUpwind[24]+0.1767766952966368*alphaDrSurf[6]*fUpwind[23]+0.1767766952966368*alphaDrSurf[2]*fUpwind[15]; 
  drag_incr[25] = 0.1767766952966368*alphaDrSurf[6]*fUpwind[31]+0.1767766952966368*alphaDrSurf[2]*fUpwind[30]+0.1767766952966368*alphaDrSurf[1]*fUpwind[29]+0.1767766952966368*alphaDrSurf[16]*fUpwind[28]+0.1767766952966368*alphaDrSurf[0]*fUpwind[25]+0.1767766952966368*alphaDrSurf[8]*fUpwind[24]+0.1767766952966368*alphaDrSurf[7]*fUpwind[23]+0.1767766952966368*alphaDrSurf[3]*fUpwind[15]; 
  drag_incr[26] = 0.1767766952966368*alphaDrSurf[0]*fUpwind[26]+0.1767766952966368*alphaDrSurf[1]*fUpwind[19]+0.1767766952966368*alphaDrSurf[2]*fUpwind[18]+0.1767766952966368*alphaDrSurf[3]*fUpwind[17]+0.1767766952966368*fUpwind[4]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[6]*fUpwind[11]+0.1767766952966368*alphaDrSurf[7]*fUpwind[10]+0.1767766952966368*alphaDrSurf[8]*fUpwind[9]; 
  drag_incr[27] = 0.1767766952966368*alphaDrSurf[0]*fUpwind[27]+0.1767766952966368*alphaDrSurf[1]*fUpwind[22]+0.1767766952966368*alphaDrSurf[2]*fUpwind[21]+0.1767766952966368*alphaDrSurf[3]*fUpwind[20]+0.1767766952966368*fUpwind[5]*alphaDrSurf[16]+0.1767766952966368*alphaDrSurf[6]*fUpwind[14]+0.1767766952966368*alphaDrSurf[7]*fUpwind[13]+0.1767766952966368*alphaDrSurf[8]*fUpwind[12]; 
  drag_incr[28] = 0.1767766952966368*alphaDrSurf[3]*fUpwind[31]+0.1767766952966368*alphaDrSurf[7]*fUpwind[30]+0.1767766952966368*alphaDrSurf[8]*fUpwind[29]+0.1767766952966368*alphaDrSurf[0]*fUpwind[28]+0.1767766952966368*alphaDrSurf[16]*fUpwind[25]+0.1767766952966368*alphaDrSurf[1]*fUpwind[24]+0.1767766952966368*alphaDrSurf[2]*fUpwind[23]+0.1767766952966368*alphaDrSurf[6]*fUpwind[15]; 
  drag_incr[29] = 0.1767766952966368*alphaDrSurf[2]*fUpwind[31]+0.1767766952966368*alphaDrSurf[6]*fUpwind[30]+0.1767766952966368*alphaDrSurf[0]*fUpwind[29]+0.1767766952966368*alphaDrSurf[8]*fUpwind[28]+0.1767766952966368*alphaDrSurf[1]*fUpwind[25]+0.1767766952966368*alphaDrSurf[16]*fUpwind[24]+0.1767766952966368*alphaDrSurf[3]*fUpwind[23]+0.1767766952966368*alphaDrSurf[7]*fUpwind[15]; 
  drag_incr[30] = 0.1767766952966368*alphaDrSurf[1]*fUpwind[31]+0.1767766952966368*alphaDrSurf[0]*fUpwind[30]+0.1767766952966368*alphaDrSurf[6]*fUpwind[29]+0.1767766952966368*alphaDrSurf[7]*fUpwind[28]+0.1767766952966368*alphaDrSurf[2]*fUpwind[25]+0.1767766952966368*alphaDrSurf[3]*fUpwind[24]+0.1767766952966368*alphaDrSurf[16]*fUpwind[23]+0.1767766952966368*alphaDrSurf[8]*fUpwind[15]; 
  drag_incr[31] = 0.1767766952966368*alphaDrSurf[0]*fUpwind[31]+0.1767766952966368*alphaDrSurf[1]*fUpwind[30]+0.1767766952966368*alphaDrSurf[2]*fUpwind[29]+0.1767766952966368*alphaDrSurf[3]*fUpwind[28]+0.1767766952966368*alphaDrSurf[6]*fUpwind[25]+0.1767766952966368*alphaDrSurf[7]*fUpwind[24]+0.1767766952966368*alphaDrSurf[8]*fUpwind[23]+0.1767766952966368*fUpwind[15]*alphaDrSurf[16]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[6] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[9] += -0.7071067811865475*drag_incr[8]*rdv2; 
  out[10] += -0.7071067811865475*drag_incr[9]*rdv2; 
  out[11] += -0.7071067811865475*drag_incr[10]*rdv2; 
  out[12] += -0.7071067811865475*drag_incr[11]*rdv2; 
  out[13] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[14] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[16] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[17] += -0.7071067811865475*drag_incr[12]*rdv2; 
  out[18] += -0.7071067811865475*drag_incr[13]*rdv2; 
  out[19] += -0.7071067811865475*drag_incr[14]*rdv2; 
  out[20] += -0.7071067811865475*drag_incr[15]*rdv2; 
  out[21] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[22] += -0.7071067811865475*drag_incr[16]*rdv2; 
  out[23] += -0.7071067811865475*drag_incr[17]*rdv2; 
  out[24] += -0.7071067811865475*drag_incr[18]*rdv2; 
  out[25] += -0.7071067811865475*drag_incr[19]*rdv2; 
  out[26] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[27] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[28] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[30] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[32] += -0.7071067811865475*drag_incr[20]*rdv2; 
  out[33] += -0.7071067811865475*drag_incr[21]*rdv2; 
  out[34] += -0.7071067811865475*drag_incr[22]*rdv2; 
  out[35] += -0.7071067811865475*drag_incr[23]*rdv2; 
  out[36] += -0.7071067811865475*drag_incr[24]*rdv2; 
  out[37] += -0.7071067811865475*drag_incr[25]*rdv2; 
  out[38] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[39] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[40] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[41] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[42] += -0.7071067811865475*drag_incr[26]*rdv2; 
  out[43] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[44] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[45] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[46] += 1.224744871391589*drag_incr[19]*rdv2; 
  out[47] += -0.7071067811865475*drag_incr[27]*rdv2; 
  out[48] += -0.7071067811865475*drag_incr[28]*rdv2; 
  out[49] += -0.7071067811865475*drag_incr[29]*rdv2; 
  out[50] += -0.7071067811865475*drag_incr[30]*rdv2; 
  out[51] += 1.224744871391589*drag_incr[20]*rdv2; 
  out[52] += 1.224744871391589*drag_incr[21]*rdv2; 
  out[53] += 1.224744871391589*drag_incr[22]*rdv2; 
  out[54] += 1.224744871391589*drag_incr[23]*rdv2; 
  out[55] += 1.224744871391589*drag_incr[24]*rdv2; 
  out[56] += 1.224744871391589*drag_incr[25]*rdv2; 
  out[57] += 1.224744871391589*drag_incr[26]*rdv2; 
  out[58] += -0.7071067811865475*drag_incr[31]*rdv2; 
  out[59] += 1.224744871391589*drag_incr[27]*rdv2; 
  out[60] += 1.224744871391589*drag_incr[28]*rdv2; 
  out[61] += 1.224744871391589*drag_incr[29]*rdv2; 
  out[62] += 1.224744871391589*drag_incr[30]*rdv2; 
  out[63] += 1.224744871391589*drag_incr[31]*rdv2; 

  } 
} 
