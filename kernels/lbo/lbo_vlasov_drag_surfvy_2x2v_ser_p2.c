#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x2v_p2_surfvy_quad.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         cell-center coordinates. 
  // dxv[4]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[16]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[8]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  const double *sumNuUy = &nuUSum[8]; 

  double alphaDrSurf_l[20] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[3]-0.7071067811865475*nuSum[0]*dxv[3]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_l[1] = 1.414213562373095*nuSum[1]*w[3]-0.7071067811865475*nuSum[1]*dxv[3]-1.414213562373095*sumNuUy[1]; 
  alphaDrSurf_l[2] = 1.414213562373095*nuSum[2]*w[3]-0.7071067811865475*nuSum[2]*dxv[3]-1.414213562373095*sumNuUy[2]; 
  alphaDrSurf_l[4] = 1.414213562373095*nuSum[3]*w[3]-1.414213562373095*sumNuUy[3]-0.7071067811865475*dxv[3]*nuSum[3]; 
  alphaDrSurf_l[7] = (-1.414213562373095*sumNuUy[4])+1.414213562373095*w[3]*nuSum[4]-0.7071067811865475*dxv[3]*nuSum[4]; 
  alphaDrSurf_l[8] = (-1.414213562373095*sumNuUy[5])+1.414213562373095*w[3]*nuSum[5]-0.7071067811865475*dxv[3]*nuSum[5]; 
  alphaDrSurf_l[11] = (-1.414213562373095*sumNuUy[6])+1.414213562373095*w[3]*nuSum[6]-0.7071067811865475*dxv[3]*nuSum[6]; 
  alphaDrSurf_l[12] = (-1.414213562373095*sumNuUy[7])+1.414213562373095*w[3]*nuSum[7]-0.7071067811865475*dxv[3]*nuSum[7]; 

  double alphaDrSurf_r[20] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[3]+0.7071067811865475*nuSum[0]*dxv[3]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*nuSum[1]*w[3]+0.7071067811865475*nuSum[1]*dxv[3]-1.414213562373095*sumNuUy[1]; 
  alphaDrSurf_r[2] = 1.414213562373095*nuSum[2]*w[3]+0.7071067811865475*nuSum[2]*dxv[3]-1.414213562373095*sumNuUy[2]; 
  alphaDrSurf_r[4] = 1.414213562373095*nuSum[3]*w[3]-1.414213562373095*sumNuUy[3]+0.7071067811865475*dxv[3]*nuSum[3]; 
  alphaDrSurf_r[7] = (-1.414213562373095*sumNuUy[4])+1.414213562373095*w[3]*nuSum[4]+0.7071067811865475*dxv[3]*nuSum[4]; 
  alphaDrSurf_r[8] = (-1.414213562373095*sumNuUy[5])+1.414213562373095*w[3]*nuSum[5]+0.7071067811865475*dxv[3]*nuSum[5]; 
  alphaDrSurf_r[11] = (-1.414213562373095*sumNuUy[6])+1.414213562373095*w[3]*nuSum[6]+0.7071067811865475*dxv[3]*nuSum[6]; 
  alphaDrSurf_r[12] = (-1.414213562373095*sumNuUy[7])+1.414213562373095*w[3]*nuSum[7]+0.7071067811865475*dxv[3]*nuSum[7]; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};;
  double fUpwind_r[20] = {0.0};
  double Gdrag_l[20] = {0.0}; 
  double Gdrag_r[20] = {0.0}; 

  if ((-0.4242640687119281*alphaDrSurf_l[12])-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_2x2v_p2_surfvy_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x2v_p2_surfvy_quad_0(-1, fc); 
  } 
  if (0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_2x2v_p2_surfvy_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x2v_p2_surfvy_quad_1(-1, fc); 
  } 
  if (0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*alphaDrSurf_l[2]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_2x2v_p2_surfvy_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x2v_p2_surfvy_quad_2(-1, fc); 
  } 
  if (0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_2x2v_p2_surfvy_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = ser_2x2v_p2_surfvy_quad_3(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*(alphaDrSurf_l[8]+alphaDrSurf_l[7]) < 0) { 
    fUpwindQuad_l[4] = ser_2x2v_p2_surfvy_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = ser_2x2v_p2_surfvy_quad_4(-1, fc); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_l[12])-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_2x2v_p2_surfvy_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = ser_2x2v_p2_surfvy_quad_5(-1, fc); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_l[12])+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*alphaDrSurf_l[2]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_2x2v_p2_surfvy_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = ser_2x2v_p2_surfvy_quad_6(-1, fc); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_l[11])+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[7] = ser_2x2v_p2_surfvy_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = ser_2x2v_p2_surfvy_quad_7(-1, fc); 
  } 
  if (0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = ser_2x2v_p2_surfvy_quad_8(1, fl); 
  } else { 
    fUpwindQuad_l[8] = ser_2x2v_p2_surfvy_quad_8(-1, fc); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_l[12])-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[9] = ser_2x2v_p2_surfvy_quad_9(1, fl); 
  } else { 
    fUpwindQuad_l[9] = ser_2x2v_p2_surfvy_quad_9(-1, fc); 
  } 
  if (0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[10] = ser_2x2v_p2_surfvy_quad_10(1, fl); 
  } else { 
    fUpwindQuad_l[10] = ser_2x2v_p2_surfvy_quad_10(-1, fc); 
  } 
  if (0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*alphaDrSurf_l[2]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[11] = ser_2x2v_p2_surfvy_quad_11(1, fl); 
  } else { 
    fUpwindQuad_l[11] = ser_2x2v_p2_surfvy_quad_11(-1, fc); 
  } 
  if (0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[12] = ser_2x2v_p2_surfvy_quad_12(1, fl); 
  } else { 
    fUpwindQuad_l[12] = ser_2x2v_p2_surfvy_quad_12(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*(alphaDrSurf_l[8]+alphaDrSurf_l[7]) < 0) { 
    fUpwindQuad_l[13] = ser_2x2v_p2_surfvy_quad_13(1, fl); 
  } else { 
    fUpwindQuad_l[13] = ser_2x2v_p2_surfvy_quad_13(-1, fc); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_l[12])-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[14] = ser_2x2v_p2_surfvy_quad_14(1, fl); 
  } else { 
    fUpwindQuad_l[14] = ser_2x2v_p2_surfvy_quad_14(-1, fc); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_l[12])+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*alphaDrSurf_l[2]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[15] = ser_2x2v_p2_surfvy_quad_15(1, fl); 
  } else { 
    fUpwindQuad_l[15] = ser_2x2v_p2_surfvy_quad_15(-1, fc); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_l[11])+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[16] = ser_2x2v_p2_surfvy_quad_16(1, fl); 
  } else { 
    fUpwindQuad_l[16] = ser_2x2v_p2_surfvy_quad_16(-1, fc); 
  } 
  if (0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[17] = ser_2x2v_p2_surfvy_quad_17(1, fl); 
  } else { 
    fUpwindQuad_l[17] = ser_2x2v_p2_surfvy_quad_17(-1, fc); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_l[12])-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[18] = ser_2x2v_p2_surfvy_quad_18(1, fl); 
  } else { 
    fUpwindQuad_l[18] = ser_2x2v_p2_surfvy_quad_18(-1, fc); 
  } 
  if (0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[19] = ser_2x2v_p2_surfvy_quad_19(1, fl); 
  } else { 
    fUpwindQuad_l[19] = ser_2x2v_p2_surfvy_quad_19(-1, fc); 
  } 
  if (0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*alphaDrSurf_l[2]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[20] = ser_2x2v_p2_surfvy_quad_20(1, fl); 
  } else { 
    fUpwindQuad_l[20] = ser_2x2v_p2_surfvy_quad_20(-1, fc); 
  } 
  if (0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[21] = ser_2x2v_p2_surfvy_quad_21(1, fl); 
  } else { 
    fUpwindQuad_l[21] = ser_2x2v_p2_surfvy_quad_21(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*(alphaDrSurf_l[8]+alphaDrSurf_l[7]) < 0) { 
    fUpwindQuad_l[22] = ser_2x2v_p2_surfvy_quad_22(1, fl); 
  } else { 
    fUpwindQuad_l[22] = ser_2x2v_p2_surfvy_quad_22(-1, fc); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_l[12])-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[23] = ser_2x2v_p2_surfvy_quad_23(1, fl); 
  } else { 
    fUpwindQuad_l[23] = ser_2x2v_p2_surfvy_quad_23(-1, fc); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_l[12])+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*alphaDrSurf_l[2]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[24] = ser_2x2v_p2_surfvy_quad_24(1, fl); 
  } else { 
    fUpwindQuad_l[24] = ser_2x2v_p2_surfvy_quad_24(-1, fc); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_l[11])+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[25] = ser_2x2v_p2_surfvy_quad_25(1, fl); 
  } else { 
    fUpwindQuad_l[25] = ser_2x2v_p2_surfvy_quad_25(-1, fc); 
  } 
  if (0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[26] = ser_2x2v_p2_surfvy_quad_26(1, fl); 
  } else { 
    fUpwindQuad_l[26] = ser_2x2v_p2_surfvy_quad_26(-1, fc); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_r[12])-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_2x2v_p2_surfvy_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x2v_p2_surfvy_quad_0(-1, fr); 
  } 
  if (0.5303300858899104*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_2x2v_p2_surfvy_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x2v_p2_surfvy_quad_1(-1, fr); 
  } 
  if (0.4242640687119281*alphaDrSurf_r[12]-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*alphaDrSurf_r[2]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_2x2v_p2_surfvy_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x2v_p2_surfvy_quad_2(-1, fr); 
  } 
  if (0.5303300858899104*alphaDrSurf_r[12]-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_2x2v_p2_surfvy_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = ser_2x2v_p2_surfvy_quad_3(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*(alphaDrSurf_r[8]+alphaDrSurf_r[7]) < 0) { 
    fUpwindQuad_r[4] = ser_2x2v_p2_surfvy_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = ser_2x2v_p2_surfvy_quad_4(-1, fr); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_r[12])-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_2x2v_p2_surfvy_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = ser_2x2v_p2_surfvy_quad_5(-1, fr); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_r[12])+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*alphaDrSurf_r[2]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_2x2v_p2_surfvy_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = ser_2x2v_p2_surfvy_quad_6(-1, fr); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_r[11])+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[7] = ser_2x2v_p2_surfvy_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = ser_2x2v_p2_surfvy_quad_7(-1, fr); 
  } 
  if (0.4242640687119281*alphaDrSurf_r[12]+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = ser_2x2v_p2_surfvy_quad_8(1, fc); 
  } else { 
    fUpwindQuad_r[8] = ser_2x2v_p2_surfvy_quad_8(-1, fr); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_r[12])-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[9] = ser_2x2v_p2_surfvy_quad_9(1, fc); 
  } else { 
    fUpwindQuad_r[9] = ser_2x2v_p2_surfvy_quad_9(-1, fr); 
  } 
  if (0.5303300858899104*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[10] = ser_2x2v_p2_surfvy_quad_10(1, fc); 
  } else { 
    fUpwindQuad_r[10] = ser_2x2v_p2_surfvy_quad_10(-1, fr); 
  } 
  if (0.4242640687119281*alphaDrSurf_r[12]-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*alphaDrSurf_r[2]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[11] = ser_2x2v_p2_surfvy_quad_11(1, fc); 
  } else { 
    fUpwindQuad_r[11] = ser_2x2v_p2_surfvy_quad_11(-1, fr); 
  } 
  if (0.5303300858899104*alphaDrSurf_r[12]-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[12] = ser_2x2v_p2_surfvy_quad_12(1, fc); 
  } else { 
    fUpwindQuad_r[12] = ser_2x2v_p2_surfvy_quad_12(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*(alphaDrSurf_r[8]+alphaDrSurf_r[7]) < 0) { 
    fUpwindQuad_r[13] = ser_2x2v_p2_surfvy_quad_13(1, fc); 
  } else { 
    fUpwindQuad_r[13] = ser_2x2v_p2_surfvy_quad_13(-1, fr); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_r[12])-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[14] = ser_2x2v_p2_surfvy_quad_14(1, fc); 
  } else { 
    fUpwindQuad_r[14] = ser_2x2v_p2_surfvy_quad_14(-1, fr); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_r[12])+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*alphaDrSurf_r[2]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[15] = ser_2x2v_p2_surfvy_quad_15(1, fc); 
  } else { 
    fUpwindQuad_r[15] = ser_2x2v_p2_surfvy_quad_15(-1, fr); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_r[11])+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[16] = ser_2x2v_p2_surfvy_quad_16(1, fc); 
  } else { 
    fUpwindQuad_r[16] = ser_2x2v_p2_surfvy_quad_16(-1, fr); 
  } 
  if (0.4242640687119281*alphaDrSurf_r[12]+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[17] = ser_2x2v_p2_surfvy_quad_17(1, fc); 
  } else { 
    fUpwindQuad_r[17] = ser_2x2v_p2_surfvy_quad_17(-1, fr); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_r[12])-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[18] = ser_2x2v_p2_surfvy_quad_18(1, fc); 
  } else { 
    fUpwindQuad_r[18] = ser_2x2v_p2_surfvy_quad_18(-1, fr); 
  } 
  if (0.5303300858899104*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[19] = ser_2x2v_p2_surfvy_quad_19(1, fc); 
  } else { 
    fUpwindQuad_r[19] = ser_2x2v_p2_surfvy_quad_19(-1, fr); 
  } 
  if (0.4242640687119281*alphaDrSurf_r[12]-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*alphaDrSurf_r[2]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[20] = ser_2x2v_p2_surfvy_quad_20(1, fc); 
  } else { 
    fUpwindQuad_r[20] = ser_2x2v_p2_surfvy_quad_20(-1, fr); 
  } 
  if (0.5303300858899104*alphaDrSurf_r[12]-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[21] = ser_2x2v_p2_surfvy_quad_21(1, fc); 
  } else { 
    fUpwindQuad_r[21] = ser_2x2v_p2_surfvy_quad_21(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*(alphaDrSurf_r[8]+alphaDrSurf_r[7]) < 0) { 
    fUpwindQuad_r[22] = ser_2x2v_p2_surfvy_quad_22(1, fc); 
  } else { 
    fUpwindQuad_r[22] = ser_2x2v_p2_surfvy_quad_22(-1, fr); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_r[12])-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[23] = ser_2x2v_p2_surfvy_quad_23(1, fc); 
  } else { 
    fUpwindQuad_r[23] = ser_2x2v_p2_surfvy_quad_23(-1, fr); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_r[12])+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*alphaDrSurf_r[2]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[24] = ser_2x2v_p2_surfvy_quad_24(1, fc); 
  } else { 
    fUpwindQuad_r[24] = ser_2x2v_p2_surfvy_quad_24(-1, fr); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_r[11])+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[25] = ser_2x2v_p2_surfvy_quad_25(1, fc); 
  } else { 
    fUpwindQuad_r[25] = ser_2x2v_p2_surfvy_quad_25(-1, fr); 
  } 
  if (0.4242640687119281*alphaDrSurf_r[12]+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[26] = ser_2x2v_p2_surfvy_quad_26(1, fc); 
  } else { 
    fUpwindQuad_r[26] = ser_2x2v_p2_surfvy_quad_26(-1, fr); 
  } 
  fUpwind_l[0] = 0.06062300936098657*fUpwindQuad_l[26]+0.09699681497757856*fUpwindQuad_l[25]+0.06062300936098657*fUpwindQuad_l[24]+0.09699681497757856*fUpwindQuad_l[23]+0.1551949039641257*fUpwindQuad_l[22]+0.09699681497757856*fUpwindQuad_l[21]+0.06062300936098657*fUpwindQuad_l[20]+0.09699681497757856*fUpwindQuad_l[19]+0.06062300936098657*fUpwindQuad_l[18]+0.09699681497757856*fUpwindQuad_l[17]+0.1551949039641257*fUpwindQuad_l[16]+0.09699681497757856*fUpwindQuad_l[15]+0.1551949039641257*fUpwindQuad_l[14]+0.2483118463426013*fUpwindQuad_l[13]+0.1551949039641257*fUpwindQuad_l[12]+0.09699681497757856*fUpwindQuad_l[11]+0.1551949039641257*fUpwindQuad_l[10]+0.09699681497757856*fUpwindQuad_l[9]+0.06062300936098657*fUpwindQuad_l[8]+0.09699681497757856*fUpwindQuad_l[7]+0.06062300936098657*fUpwindQuad_l[6]+0.09699681497757856*fUpwindQuad_l[5]+0.1551949039641257*fUpwindQuad_l[4]+0.09699681497757856*fUpwindQuad_l[3]+0.06062300936098657*fUpwindQuad_l[2]+0.09699681497757856*fUpwindQuad_l[1]+0.06062300936098657*fUpwindQuad_l[0]; 
  fUpwind_l[1] = 0.08133430195906327*fUpwindQuad_l[26]-0.08133430195906327*fUpwindQuad_l[24]+0.1301348831345013*fUpwindQuad_l[23]-0.1301348831345013*fUpwindQuad_l[21]+0.08133430195906327*fUpwindQuad_l[20]-0.08133430195906327*fUpwindQuad_l[18]+0.1301348831345013*fUpwindQuad_l[17]-0.1301348831345013*fUpwindQuad_l[15]+0.2082158130152021*fUpwindQuad_l[14]-0.2082158130152021*fUpwindQuad_l[12]+0.1301348831345013*fUpwindQuad_l[11]-0.1301348831345013*fUpwindQuad_l[9]+0.08133430195906327*fUpwindQuad_l[8]-0.08133430195906327*fUpwindQuad_l[6]+0.1301348831345013*fUpwindQuad_l[5]-0.1301348831345013*fUpwindQuad_l[3]+0.08133430195906327*fUpwindQuad_l[2]-0.08133430195906327*fUpwindQuad_l[0]; 
  fUpwind_l[2] = 0.08133430195906327*fUpwindQuad_l[26]+0.1301348831345013*fUpwindQuad_l[25]+0.08133430195906327*fUpwindQuad_l[24]-0.08133430195906327*fUpwindQuad_l[20]-0.1301348831345013*fUpwindQuad_l[19]-0.08133430195906327*fUpwindQuad_l[18]+0.1301348831345013*fUpwindQuad_l[17]+0.2082158130152021*fUpwindQuad_l[16]+0.1301348831345013*fUpwindQuad_l[15]-0.1301348831345013*fUpwindQuad_l[11]-0.2082158130152021*fUpwindQuad_l[10]-0.1301348831345013*fUpwindQuad_l[9]+0.08133430195906327*fUpwindQuad_l[8]+0.1301348831345013*fUpwindQuad_l[7]+0.08133430195906327*fUpwindQuad_l[6]-0.08133430195906327*fUpwindQuad_l[2]-0.1301348831345013*fUpwindQuad_l[1]-0.08133430195906327*fUpwindQuad_l[0]; 
  fUpwind_l[3] = 0.08133430195906327*fUpwindQuad_l[26]+0.1301348831345013*fUpwindQuad_l[25]+0.08133430195906327*fUpwindQuad_l[24]+0.1301348831345013*fUpwindQuad_l[23]+0.2082158130152021*fUpwindQuad_l[22]+0.1301348831345013*fUpwindQuad_l[21]+0.08133430195906327*fUpwindQuad_l[20]+0.1301348831345013*fUpwindQuad_l[19]+0.08133430195906327*fUpwindQuad_l[18]-0.08133430195906327*fUpwindQuad_l[8]-0.1301348831345013*fUpwindQuad_l[7]-0.08133430195906327*fUpwindQuad_l[6]-0.1301348831345013*fUpwindQuad_l[5]-0.2082158130152021*fUpwindQuad_l[4]-0.1301348831345013*fUpwindQuad_l[3]-0.08133430195906327*fUpwindQuad_l[2]-0.1301348831345013*fUpwindQuad_l[1]-0.08133430195906327*fUpwindQuad_l[0]; 
  fUpwind_l[4] = 0.1091214168497758*fUpwindQuad_l[26]-0.1091214168497758*fUpwindQuad_l[24]-0.1091214168497758*fUpwindQuad_l[20]+0.1091214168497758*fUpwindQuad_l[18]+0.1745942669596414*fUpwindQuad_l[17]-0.1745942669596414*fUpwindQuad_l[15]-0.1745942669596414*fUpwindQuad_l[11]+0.1745942669596414*fUpwindQuad_l[9]+0.1091214168497758*fUpwindQuad_l[8]-0.1091214168497758*fUpwindQuad_l[6]-0.1091214168497758*fUpwindQuad_l[2]+0.1091214168497758*fUpwindQuad_l[0]; 
  fUpwind_l[5] = 0.1091214168497758*fUpwindQuad_l[26]-0.1091214168497758*fUpwindQuad_l[24]+0.1745942669596414*fUpwindQuad_l[23]-0.1745942669596414*fUpwindQuad_l[21]+0.1091214168497758*fUpwindQuad_l[20]-0.1091214168497758*fUpwindQuad_l[18]-0.1091214168497758*fUpwindQuad_l[8]+0.1091214168497758*fUpwindQuad_l[6]-0.1745942669596414*fUpwindQuad_l[5]+0.1745942669596414*fUpwindQuad_l[3]-0.1091214168497758*fUpwindQuad_l[2]+0.1091214168497758*fUpwindQuad_l[0]; 
  fUpwind_l[6] = 0.1091214168497758*fUpwindQuad_l[26]+0.1745942669596414*fUpwindQuad_l[25]+0.1091214168497758*fUpwindQuad_l[24]-0.1091214168497758*fUpwindQuad_l[20]-0.1745942669596414*fUpwindQuad_l[19]-0.1091214168497758*fUpwindQuad_l[18]-0.1091214168497758*fUpwindQuad_l[8]-0.1745942669596414*fUpwindQuad_l[7]-0.1091214168497758*fUpwindQuad_l[6]+0.1091214168497758*fUpwindQuad_l[2]+0.1745942669596414*fUpwindQuad_l[1]+0.1091214168497758*fUpwindQuad_l[0]; 
  fUpwind_l[7] = 0.05422286797270884*fUpwindQuad_l[26]-0.1084457359454177*fUpwindQuad_l[25]+0.05422286797270884*fUpwindQuad_l[24]+0.08675658875633419*fUpwindQuad_l[23]-0.1735131775126684*fUpwindQuad_l[22]+0.08675658875633419*fUpwindQuad_l[21]+0.05422286797270884*fUpwindQuad_l[20]-0.1084457359454177*fUpwindQuad_l[19]+0.05422286797270884*fUpwindQuad_l[18]+0.08675658875633419*fUpwindQuad_l[17]-0.1735131775126684*fUpwindQuad_l[16]+0.08675658875633419*fUpwindQuad_l[15]+0.1388105420101347*fUpwindQuad_l[14]-0.2776210840202695*fUpwindQuad_l[13]+0.1388105420101347*fUpwindQuad_l[12]+0.08675658875633419*fUpwindQuad_l[11]-0.1735131775126684*fUpwindQuad_l[10]+0.08675658875633419*fUpwindQuad_l[9]+0.05422286797270884*fUpwindQuad_l[8]-0.1084457359454177*fUpwindQuad_l[7]+0.05422286797270884*fUpwindQuad_l[6]+0.08675658875633419*fUpwindQuad_l[5]-0.1735131775126684*fUpwindQuad_l[4]+0.08675658875633419*fUpwindQuad_l[3]+0.05422286797270884*fUpwindQuad_l[2]-0.1084457359454177*fUpwindQuad_l[1]+0.05422286797270884*fUpwindQuad_l[0]; 
  fUpwind_l[8] = 0.05422286797270884*fUpwindQuad_l[26]+0.08675658875633419*fUpwindQuad_l[25]+0.05422286797270884*fUpwindQuad_l[24]-0.1084457359454177*fUpwindQuad_l[23]-0.1735131775126684*fUpwindQuad_l[22]-0.1084457359454177*fUpwindQuad_l[21]+0.05422286797270884*fUpwindQuad_l[20]+0.08675658875633419*fUpwindQuad_l[19]+0.05422286797270884*fUpwindQuad_l[18]+0.08675658875633419*fUpwindQuad_l[17]+0.1388105420101347*fUpwindQuad_l[16]+0.08675658875633419*fUpwindQuad_l[15]-0.1735131775126684*fUpwindQuad_l[14]-0.2776210840202695*fUpwindQuad_l[13]-0.1735131775126684*fUpwindQuad_l[12]+0.08675658875633419*fUpwindQuad_l[11]+0.1388105420101347*fUpwindQuad_l[10]+0.08675658875633419*fUpwindQuad_l[9]+0.05422286797270884*fUpwindQuad_l[8]+0.08675658875633419*fUpwindQuad_l[7]+0.05422286797270884*fUpwindQuad_l[6]-0.1084457359454177*fUpwindQuad_l[5]-0.1735131775126684*fUpwindQuad_l[4]-0.1084457359454177*fUpwindQuad_l[3]+0.05422286797270884*fUpwindQuad_l[2]+0.08675658875633419*fUpwindQuad_l[1]+0.05422286797270884*fUpwindQuad_l[0]; 
  fUpwind_l[9] = 0.05422286797270884*fUpwindQuad_l[26]+0.08675658875633419*fUpwindQuad_l[25]+0.05422286797270884*fUpwindQuad_l[24]+0.08675658875633419*fUpwindQuad_l[23]+0.1388105420101347*fUpwindQuad_l[22]+0.08675658875633419*fUpwindQuad_l[21]+0.05422286797270884*fUpwindQuad_l[20]+0.08675658875633419*fUpwindQuad_l[19]+0.05422286797270884*fUpwindQuad_l[18]-0.1084457359454177*fUpwindQuad_l[17]-0.1735131775126684*fUpwindQuad_l[16]-0.1084457359454177*fUpwindQuad_l[15]-0.1735131775126684*fUpwindQuad_l[14]-0.2776210840202695*fUpwindQuad_l[13]-0.1735131775126684*fUpwindQuad_l[12]-0.1084457359454177*fUpwindQuad_l[11]-0.1735131775126684*fUpwindQuad_l[10]-0.1084457359454177*fUpwindQuad_l[9]+0.05422286797270884*fUpwindQuad_l[8]+0.08675658875633419*fUpwindQuad_l[7]+0.05422286797270884*fUpwindQuad_l[6]+0.08675658875633419*fUpwindQuad_l[5]+0.1388105420101347*fUpwindQuad_l[4]+0.08675658875633419*fUpwindQuad_l[3]+0.05422286797270884*fUpwindQuad_l[2]+0.08675658875633419*fUpwindQuad_l[1]+0.05422286797270884*fUpwindQuad_l[0]; 
  fUpwind_l[10] = 0.1464017435263139*fUpwindQuad_l[26]-0.1464017435263139*fUpwindQuad_l[24]-0.1464017435263139*fUpwindQuad_l[20]+0.1464017435263139*fUpwindQuad_l[18]-0.1464017435263139*fUpwindQuad_l[8]+0.1464017435263139*fUpwindQuad_l[6]+0.1464017435263139*fUpwindQuad_l[2]-0.1464017435263139*fUpwindQuad_l[0]; 
  fUpwind_l[11] = 0.07274761123318395*fUpwindQuad_l[26]-0.1454952224663679*fUpwindQuad_l[25]+0.07274761123318395*fUpwindQuad_l[24]-0.07274761123318395*fUpwindQuad_l[20]+0.1454952224663679*fUpwindQuad_l[19]-0.07274761123318395*fUpwindQuad_l[18]+0.1163961779730944*fUpwindQuad_l[17]-0.2327923559461888*fUpwindQuad_l[16]+0.1163961779730944*fUpwindQuad_l[15]-0.1163961779730944*fUpwindQuad_l[11]+0.2327923559461888*fUpwindQuad_l[10]-0.1163961779730944*fUpwindQuad_l[9]+0.07274761123318395*fUpwindQuad_l[8]-0.1454952224663679*fUpwindQuad_l[7]+0.07274761123318395*fUpwindQuad_l[6]-0.07274761123318395*fUpwindQuad_l[2]+0.1454952224663679*fUpwindQuad_l[1]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[12] = 0.07274761123318395*fUpwindQuad_l[26]-0.07274761123318395*fUpwindQuad_l[24]-0.1454952224663679*fUpwindQuad_l[23]+0.1454952224663679*fUpwindQuad_l[21]+0.07274761123318395*fUpwindQuad_l[20]-0.07274761123318395*fUpwindQuad_l[18]+0.1163961779730944*fUpwindQuad_l[17]-0.1163961779730944*fUpwindQuad_l[15]-0.2327923559461888*fUpwindQuad_l[14]+0.2327923559461888*fUpwindQuad_l[12]+0.1163961779730944*fUpwindQuad_l[11]-0.1163961779730944*fUpwindQuad_l[9]+0.07274761123318395*fUpwindQuad_l[8]-0.07274761123318395*fUpwindQuad_l[6]-0.1454952224663679*fUpwindQuad_l[5]+0.1454952224663679*fUpwindQuad_l[3]+0.07274761123318395*fUpwindQuad_l[2]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[13] = 0.07274761123318395*fUpwindQuad_l[26]-0.1454952224663679*fUpwindQuad_l[25]+0.07274761123318395*fUpwindQuad_l[24]+0.1163961779730944*fUpwindQuad_l[23]-0.2327923559461888*fUpwindQuad_l[22]+0.1163961779730944*fUpwindQuad_l[21]+0.07274761123318395*fUpwindQuad_l[20]-0.1454952224663679*fUpwindQuad_l[19]+0.07274761123318395*fUpwindQuad_l[18]-0.07274761123318395*fUpwindQuad_l[8]+0.1454952224663679*fUpwindQuad_l[7]-0.07274761123318395*fUpwindQuad_l[6]-0.1163961779730944*fUpwindQuad_l[5]+0.2327923559461888*fUpwindQuad_l[4]-0.1163961779730944*fUpwindQuad_l[3]-0.07274761123318395*fUpwindQuad_l[2]+0.1454952224663679*fUpwindQuad_l[1]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[14] = 0.07274761123318395*fUpwindQuad_l[26]+0.1163961779730944*fUpwindQuad_l[25]+0.07274761123318395*fUpwindQuad_l[24]-0.1454952224663679*fUpwindQuad_l[23]-0.2327923559461888*fUpwindQuad_l[22]-0.1454952224663679*fUpwindQuad_l[21]+0.07274761123318395*fUpwindQuad_l[20]+0.1163961779730944*fUpwindQuad_l[19]+0.07274761123318395*fUpwindQuad_l[18]-0.07274761123318395*fUpwindQuad_l[8]-0.1163961779730944*fUpwindQuad_l[7]-0.07274761123318395*fUpwindQuad_l[6]+0.1454952224663679*fUpwindQuad_l[5]+0.2327923559461888*fUpwindQuad_l[4]+0.1454952224663679*fUpwindQuad_l[3]-0.07274761123318395*fUpwindQuad_l[2]-0.1163961779730944*fUpwindQuad_l[1]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[15] = 0.07274761123318395*fUpwindQuad_l[26]-0.07274761123318395*fUpwindQuad_l[24]+0.1163961779730944*fUpwindQuad_l[23]-0.1163961779730944*fUpwindQuad_l[21]+0.07274761123318395*fUpwindQuad_l[20]-0.07274761123318395*fUpwindQuad_l[18]-0.1454952224663679*fUpwindQuad_l[17]+0.1454952224663679*fUpwindQuad_l[15]-0.2327923559461888*fUpwindQuad_l[14]+0.2327923559461888*fUpwindQuad_l[12]-0.1454952224663679*fUpwindQuad_l[11]+0.1454952224663679*fUpwindQuad_l[9]+0.07274761123318395*fUpwindQuad_l[8]-0.07274761123318395*fUpwindQuad_l[6]+0.1163961779730944*fUpwindQuad_l[5]-0.1163961779730944*fUpwindQuad_l[3]+0.07274761123318395*fUpwindQuad_l[2]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[16] = 0.07274761123318395*fUpwindQuad_l[26]+0.1163961779730944*fUpwindQuad_l[25]+0.07274761123318395*fUpwindQuad_l[24]-0.07274761123318395*fUpwindQuad_l[20]-0.1163961779730944*fUpwindQuad_l[19]-0.07274761123318395*fUpwindQuad_l[18]-0.1454952224663679*fUpwindQuad_l[17]-0.2327923559461888*fUpwindQuad_l[16]-0.1454952224663679*fUpwindQuad_l[15]+0.1454952224663679*fUpwindQuad_l[11]+0.2327923559461888*fUpwindQuad_l[10]+0.1454952224663679*fUpwindQuad_l[9]+0.07274761123318395*fUpwindQuad_l[8]+0.1163961779730944*fUpwindQuad_l[7]+0.07274761123318395*fUpwindQuad_l[6]-0.07274761123318395*fUpwindQuad_l[2]-0.1163961779730944*fUpwindQuad_l[1]-0.07274761123318395*fUpwindQuad_l[0]; 
  fUpwind_l[17] = 0.09760116235087592*fUpwindQuad_l[26]-0.1952023247017519*fUpwindQuad_l[25]+0.09760116235087592*fUpwindQuad_l[24]-0.09760116235087592*fUpwindQuad_l[20]+0.1952023247017519*fUpwindQuad_l[19]-0.09760116235087592*fUpwindQuad_l[18]-0.09760116235087592*fUpwindQuad_l[8]+0.1952023247017519*fUpwindQuad_l[7]-0.09760116235087592*fUpwindQuad_l[6]+0.09760116235087592*fUpwindQuad_l[2]-0.1952023247017519*fUpwindQuad_l[1]+0.09760116235087592*fUpwindQuad_l[0]; 
  fUpwind_l[18] = 0.09760116235087592*fUpwindQuad_l[26]-0.09760116235087592*fUpwindQuad_l[24]-0.1952023247017519*fUpwindQuad_l[23]+0.1952023247017519*fUpwindQuad_l[21]+0.09760116235087592*fUpwindQuad_l[20]-0.09760116235087592*fUpwindQuad_l[18]-0.09760116235087592*fUpwindQuad_l[8]+0.09760116235087592*fUpwindQuad_l[6]+0.1952023247017519*fUpwindQuad_l[5]-0.1952023247017519*fUpwindQuad_l[3]-0.09760116235087592*fUpwindQuad_l[2]+0.09760116235087592*fUpwindQuad_l[0]; 
  fUpwind_l[19] = 0.09760116235087592*fUpwindQuad_l[26]-0.09760116235087592*fUpwindQuad_l[24]-0.09760116235087592*fUpwindQuad_l[20]+0.09760116235087592*fUpwindQuad_l[18]-0.1952023247017519*fUpwindQuad_l[17]+0.1952023247017519*fUpwindQuad_l[15]+0.1952023247017519*fUpwindQuad_l[11]-0.1952023247017519*fUpwindQuad_l[9]+0.09760116235087592*fUpwindQuad_l[8]-0.09760116235087592*fUpwindQuad_l[6]-0.09760116235087592*fUpwindQuad_l[2]+0.09760116235087592*fUpwindQuad_l[0]; 

  fUpwind_r[0] = 0.06062300936098657*fUpwindQuad_r[26]+0.09699681497757856*fUpwindQuad_r[25]+0.06062300936098657*fUpwindQuad_r[24]+0.09699681497757856*fUpwindQuad_r[23]+0.1551949039641257*fUpwindQuad_r[22]+0.09699681497757856*fUpwindQuad_r[21]+0.06062300936098657*fUpwindQuad_r[20]+0.09699681497757856*fUpwindQuad_r[19]+0.06062300936098657*fUpwindQuad_r[18]+0.09699681497757856*fUpwindQuad_r[17]+0.1551949039641257*fUpwindQuad_r[16]+0.09699681497757856*fUpwindQuad_r[15]+0.1551949039641257*fUpwindQuad_r[14]+0.2483118463426013*fUpwindQuad_r[13]+0.1551949039641257*fUpwindQuad_r[12]+0.09699681497757856*fUpwindQuad_r[11]+0.1551949039641257*fUpwindQuad_r[10]+0.09699681497757856*fUpwindQuad_r[9]+0.06062300936098657*fUpwindQuad_r[8]+0.09699681497757856*fUpwindQuad_r[7]+0.06062300936098657*fUpwindQuad_r[6]+0.09699681497757856*fUpwindQuad_r[5]+0.1551949039641257*fUpwindQuad_r[4]+0.09699681497757856*fUpwindQuad_r[3]+0.06062300936098657*fUpwindQuad_r[2]+0.09699681497757856*fUpwindQuad_r[1]+0.06062300936098657*fUpwindQuad_r[0]; 
  fUpwind_r[1] = 0.08133430195906327*fUpwindQuad_r[26]-0.08133430195906327*fUpwindQuad_r[24]+0.1301348831345013*fUpwindQuad_r[23]-0.1301348831345013*fUpwindQuad_r[21]+0.08133430195906327*fUpwindQuad_r[20]-0.08133430195906327*fUpwindQuad_r[18]+0.1301348831345013*fUpwindQuad_r[17]-0.1301348831345013*fUpwindQuad_r[15]+0.2082158130152021*fUpwindQuad_r[14]-0.2082158130152021*fUpwindQuad_r[12]+0.1301348831345013*fUpwindQuad_r[11]-0.1301348831345013*fUpwindQuad_r[9]+0.08133430195906327*fUpwindQuad_r[8]-0.08133430195906327*fUpwindQuad_r[6]+0.1301348831345013*fUpwindQuad_r[5]-0.1301348831345013*fUpwindQuad_r[3]+0.08133430195906327*fUpwindQuad_r[2]-0.08133430195906327*fUpwindQuad_r[0]; 
  fUpwind_r[2] = 0.08133430195906327*fUpwindQuad_r[26]+0.1301348831345013*fUpwindQuad_r[25]+0.08133430195906327*fUpwindQuad_r[24]-0.08133430195906327*fUpwindQuad_r[20]-0.1301348831345013*fUpwindQuad_r[19]-0.08133430195906327*fUpwindQuad_r[18]+0.1301348831345013*fUpwindQuad_r[17]+0.2082158130152021*fUpwindQuad_r[16]+0.1301348831345013*fUpwindQuad_r[15]-0.1301348831345013*fUpwindQuad_r[11]-0.2082158130152021*fUpwindQuad_r[10]-0.1301348831345013*fUpwindQuad_r[9]+0.08133430195906327*fUpwindQuad_r[8]+0.1301348831345013*fUpwindQuad_r[7]+0.08133430195906327*fUpwindQuad_r[6]-0.08133430195906327*fUpwindQuad_r[2]-0.1301348831345013*fUpwindQuad_r[1]-0.08133430195906327*fUpwindQuad_r[0]; 
  fUpwind_r[3] = 0.08133430195906327*fUpwindQuad_r[26]+0.1301348831345013*fUpwindQuad_r[25]+0.08133430195906327*fUpwindQuad_r[24]+0.1301348831345013*fUpwindQuad_r[23]+0.2082158130152021*fUpwindQuad_r[22]+0.1301348831345013*fUpwindQuad_r[21]+0.08133430195906327*fUpwindQuad_r[20]+0.1301348831345013*fUpwindQuad_r[19]+0.08133430195906327*fUpwindQuad_r[18]-0.08133430195906327*fUpwindQuad_r[8]-0.1301348831345013*fUpwindQuad_r[7]-0.08133430195906327*fUpwindQuad_r[6]-0.1301348831345013*fUpwindQuad_r[5]-0.2082158130152021*fUpwindQuad_r[4]-0.1301348831345013*fUpwindQuad_r[3]-0.08133430195906327*fUpwindQuad_r[2]-0.1301348831345013*fUpwindQuad_r[1]-0.08133430195906327*fUpwindQuad_r[0]; 
  fUpwind_r[4] = 0.1091214168497758*fUpwindQuad_r[26]-0.1091214168497758*fUpwindQuad_r[24]-0.1091214168497758*fUpwindQuad_r[20]+0.1091214168497758*fUpwindQuad_r[18]+0.1745942669596414*fUpwindQuad_r[17]-0.1745942669596414*fUpwindQuad_r[15]-0.1745942669596414*fUpwindQuad_r[11]+0.1745942669596414*fUpwindQuad_r[9]+0.1091214168497758*fUpwindQuad_r[8]-0.1091214168497758*fUpwindQuad_r[6]-0.1091214168497758*fUpwindQuad_r[2]+0.1091214168497758*fUpwindQuad_r[0]; 
  fUpwind_r[5] = 0.1091214168497758*fUpwindQuad_r[26]-0.1091214168497758*fUpwindQuad_r[24]+0.1745942669596414*fUpwindQuad_r[23]-0.1745942669596414*fUpwindQuad_r[21]+0.1091214168497758*fUpwindQuad_r[20]-0.1091214168497758*fUpwindQuad_r[18]-0.1091214168497758*fUpwindQuad_r[8]+0.1091214168497758*fUpwindQuad_r[6]-0.1745942669596414*fUpwindQuad_r[5]+0.1745942669596414*fUpwindQuad_r[3]-0.1091214168497758*fUpwindQuad_r[2]+0.1091214168497758*fUpwindQuad_r[0]; 
  fUpwind_r[6] = 0.1091214168497758*fUpwindQuad_r[26]+0.1745942669596414*fUpwindQuad_r[25]+0.1091214168497758*fUpwindQuad_r[24]-0.1091214168497758*fUpwindQuad_r[20]-0.1745942669596414*fUpwindQuad_r[19]-0.1091214168497758*fUpwindQuad_r[18]-0.1091214168497758*fUpwindQuad_r[8]-0.1745942669596414*fUpwindQuad_r[7]-0.1091214168497758*fUpwindQuad_r[6]+0.1091214168497758*fUpwindQuad_r[2]+0.1745942669596414*fUpwindQuad_r[1]+0.1091214168497758*fUpwindQuad_r[0]; 
  fUpwind_r[7] = 0.05422286797270884*fUpwindQuad_r[26]-0.1084457359454177*fUpwindQuad_r[25]+0.05422286797270884*fUpwindQuad_r[24]+0.08675658875633419*fUpwindQuad_r[23]-0.1735131775126684*fUpwindQuad_r[22]+0.08675658875633419*fUpwindQuad_r[21]+0.05422286797270884*fUpwindQuad_r[20]-0.1084457359454177*fUpwindQuad_r[19]+0.05422286797270884*fUpwindQuad_r[18]+0.08675658875633419*fUpwindQuad_r[17]-0.1735131775126684*fUpwindQuad_r[16]+0.08675658875633419*fUpwindQuad_r[15]+0.1388105420101347*fUpwindQuad_r[14]-0.2776210840202695*fUpwindQuad_r[13]+0.1388105420101347*fUpwindQuad_r[12]+0.08675658875633419*fUpwindQuad_r[11]-0.1735131775126684*fUpwindQuad_r[10]+0.08675658875633419*fUpwindQuad_r[9]+0.05422286797270884*fUpwindQuad_r[8]-0.1084457359454177*fUpwindQuad_r[7]+0.05422286797270884*fUpwindQuad_r[6]+0.08675658875633419*fUpwindQuad_r[5]-0.1735131775126684*fUpwindQuad_r[4]+0.08675658875633419*fUpwindQuad_r[3]+0.05422286797270884*fUpwindQuad_r[2]-0.1084457359454177*fUpwindQuad_r[1]+0.05422286797270884*fUpwindQuad_r[0]; 
  fUpwind_r[8] = 0.05422286797270884*fUpwindQuad_r[26]+0.08675658875633419*fUpwindQuad_r[25]+0.05422286797270884*fUpwindQuad_r[24]-0.1084457359454177*fUpwindQuad_r[23]-0.1735131775126684*fUpwindQuad_r[22]-0.1084457359454177*fUpwindQuad_r[21]+0.05422286797270884*fUpwindQuad_r[20]+0.08675658875633419*fUpwindQuad_r[19]+0.05422286797270884*fUpwindQuad_r[18]+0.08675658875633419*fUpwindQuad_r[17]+0.1388105420101347*fUpwindQuad_r[16]+0.08675658875633419*fUpwindQuad_r[15]-0.1735131775126684*fUpwindQuad_r[14]-0.2776210840202695*fUpwindQuad_r[13]-0.1735131775126684*fUpwindQuad_r[12]+0.08675658875633419*fUpwindQuad_r[11]+0.1388105420101347*fUpwindQuad_r[10]+0.08675658875633419*fUpwindQuad_r[9]+0.05422286797270884*fUpwindQuad_r[8]+0.08675658875633419*fUpwindQuad_r[7]+0.05422286797270884*fUpwindQuad_r[6]-0.1084457359454177*fUpwindQuad_r[5]-0.1735131775126684*fUpwindQuad_r[4]-0.1084457359454177*fUpwindQuad_r[3]+0.05422286797270884*fUpwindQuad_r[2]+0.08675658875633419*fUpwindQuad_r[1]+0.05422286797270884*fUpwindQuad_r[0]; 
  fUpwind_r[9] = 0.05422286797270884*fUpwindQuad_r[26]+0.08675658875633419*fUpwindQuad_r[25]+0.05422286797270884*fUpwindQuad_r[24]+0.08675658875633419*fUpwindQuad_r[23]+0.1388105420101347*fUpwindQuad_r[22]+0.08675658875633419*fUpwindQuad_r[21]+0.05422286797270884*fUpwindQuad_r[20]+0.08675658875633419*fUpwindQuad_r[19]+0.05422286797270884*fUpwindQuad_r[18]-0.1084457359454177*fUpwindQuad_r[17]-0.1735131775126684*fUpwindQuad_r[16]-0.1084457359454177*fUpwindQuad_r[15]-0.1735131775126684*fUpwindQuad_r[14]-0.2776210840202695*fUpwindQuad_r[13]-0.1735131775126684*fUpwindQuad_r[12]-0.1084457359454177*fUpwindQuad_r[11]-0.1735131775126684*fUpwindQuad_r[10]-0.1084457359454177*fUpwindQuad_r[9]+0.05422286797270884*fUpwindQuad_r[8]+0.08675658875633419*fUpwindQuad_r[7]+0.05422286797270884*fUpwindQuad_r[6]+0.08675658875633419*fUpwindQuad_r[5]+0.1388105420101347*fUpwindQuad_r[4]+0.08675658875633419*fUpwindQuad_r[3]+0.05422286797270884*fUpwindQuad_r[2]+0.08675658875633419*fUpwindQuad_r[1]+0.05422286797270884*fUpwindQuad_r[0]; 
  fUpwind_r[10] = 0.1464017435263139*fUpwindQuad_r[26]-0.1464017435263139*fUpwindQuad_r[24]-0.1464017435263139*fUpwindQuad_r[20]+0.1464017435263139*fUpwindQuad_r[18]-0.1464017435263139*fUpwindQuad_r[8]+0.1464017435263139*fUpwindQuad_r[6]+0.1464017435263139*fUpwindQuad_r[2]-0.1464017435263139*fUpwindQuad_r[0]; 
  fUpwind_r[11] = 0.07274761123318395*fUpwindQuad_r[26]-0.1454952224663679*fUpwindQuad_r[25]+0.07274761123318395*fUpwindQuad_r[24]-0.07274761123318395*fUpwindQuad_r[20]+0.1454952224663679*fUpwindQuad_r[19]-0.07274761123318395*fUpwindQuad_r[18]+0.1163961779730944*fUpwindQuad_r[17]-0.2327923559461888*fUpwindQuad_r[16]+0.1163961779730944*fUpwindQuad_r[15]-0.1163961779730944*fUpwindQuad_r[11]+0.2327923559461888*fUpwindQuad_r[10]-0.1163961779730944*fUpwindQuad_r[9]+0.07274761123318395*fUpwindQuad_r[8]-0.1454952224663679*fUpwindQuad_r[7]+0.07274761123318395*fUpwindQuad_r[6]-0.07274761123318395*fUpwindQuad_r[2]+0.1454952224663679*fUpwindQuad_r[1]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[12] = 0.07274761123318395*fUpwindQuad_r[26]-0.07274761123318395*fUpwindQuad_r[24]-0.1454952224663679*fUpwindQuad_r[23]+0.1454952224663679*fUpwindQuad_r[21]+0.07274761123318395*fUpwindQuad_r[20]-0.07274761123318395*fUpwindQuad_r[18]+0.1163961779730944*fUpwindQuad_r[17]-0.1163961779730944*fUpwindQuad_r[15]-0.2327923559461888*fUpwindQuad_r[14]+0.2327923559461888*fUpwindQuad_r[12]+0.1163961779730944*fUpwindQuad_r[11]-0.1163961779730944*fUpwindQuad_r[9]+0.07274761123318395*fUpwindQuad_r[8]-0.07274761123318395*fUpwindQuad_r[6]-0.1454952224663679*fUpwindQuad_r[5]+0.1454952224663679*fUpwindQuad_r[3]+0.07274761123318395*fUpwindQuad_r[2]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[13] = 0.07274761123318395*fUpwindQuad_r[26]-0.1454952224663679*fUpwindQuad_r[25]+0.07274761123318395*fUpwindQuad_r[24]+0.1163961779730944*fUpwindQuad_r[23]-0.2327923559461888*fUpwindQuad_r[22]+0.1163961779730944*fUpwindQuad_r[21]+0.07274761123318395*fUpwindQuad_r[20]-0.1454952224663679*fUpwindQuad_r[19]+0.07274761123318395*fUpwindQuad_r[18]-0.07274761123318395*fUpwindQuad_r[8]+0.1454952224663679*fUpwindQuad_r[7]-0.07274761123318395*fUpwindQuad_r[6]-0.1163961779730944*fUpwindQuad_r[5]+0.2327923559461888*fUpwindQuad_r[4]-0.1163961779730944*fUpwindQuad_r[3]-0.07274761123318395*fUpwindQuad_r[2]+0.1454952224663679*fUpwindQuad_r[1]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[14] = 0.07274761123318395*fUpwindQuad_r[26]+0.1163961779730944*fUpwindQuad_r[25]+0.07274761123318395*fUpwindQuad_r[24]-0.1454952224663679*fUpwindQuad_r[23]-0.2327923559461888*fUpwindQuad_r[22]-0.1454952224663679*fUpwindQuad_r[21]+0.07274761123318395*fUpwindQuad_r[20]+0.1163961779730944*fUpwindQuad_r[19]+0.07274761123318395*fUpwindQuad_r[18]-0.07274761123318395*fUpwindQuad_r[8]-0.1163961779730944*fUpwindQuad_r[7]-0.07274761123318395*fUpwindQuad_r[6]+0.1454952224663679*fUpwindQuad_r[5]+0.2327923559461888*fUpwindQuad_r[4]+0.1454952224663679*fUpwindQuad_r[3]-0.07274761123318395*fUpwindQuad_r[2]-0.1163961779730944*fUpwindQuad_r[1]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[15] = 0.07274761123318395*fUpwindQuad_r[26]-0.07274761123318395*fUpwindQuad_r[24]+0.1163961779730944*fUpwindQuad_r[23]-0.1163961779730944*fUpwindQuad_r[21]+0.07274761123318395*fUpwindQuad_r[20]-0.07274761123318395*fUpwindQuad_r[18]-0.1454952224663679*fUpwindQuad_r[17]+0.1454952224663679*fUpwindQuad_r[15]-0.2327923559461888*fUpwindQuad_r[14]+0.2327923559461888*fUpwindQuad_r[12]-0.1454952224663679*fUpwindQuad_r[11]+0.1454952224663679*fUpwindQuad_r[9]+0.07274761123318395*fUpwindQuad_r[8]-0.07274761123318395*fUpwindQuad_r[6]+0.1163961779730944*fUpwindQuad_r[5]-0.1163961779730944*fUpwindQuad_r[3]+0.07274761123318395*fUpwindQuad_r[2]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[16] = 0.07274761123318395*fUpwindQuad_r[26]+0.1163961779730944*fUpwindQuad_r[25]+0.07274761123318395*fUpwindQuad_r[24]-0.07274761123318395*fUpwindQuad_r[20]-0.1163961779730944*fUpwindQuad_r[19]-0.07274761123318395*fUpwindQuad_r[18]-0.1454952224663679*fUpwindQuad_r[17]-0.2327923559461888*fUpwindQuad_r[16]-0.1454952224663679*fUpwindQuad_r[15]+0.1454952224663679*fUpwindQuad_r[11]+0.2327923559461888*fUpwindQuad_r[10]+0.1454952224663679*fUpwindQuad_r[9]+0.07274761123318395*fUpwindQuad_r[8]+0.1163961779730944*fUpwindQuad_r[7]+0.07274761123318395*fUpwindQuad_r[6]-0.07274761123318395*fUpwindQuad_r[2]-0.1163961779730944*fUpwindQuad_r[1]-0.07274761123318395*fUpwindQuad_r[0]; 
  fUpwind_r[17] = 0.09760116235087592*fUpwindQuad_r[26]-0.1952023247017519*fUpwindQuad_r[25]+0.09760116235087592*fUpwindQuad_r[24]-0.09760116235087592*fUpwindQuad_r[20]+0.1952023247017519*fUpwindQuad_r[19]-0.09760116235087592*fUpwindQuad_r[18]-0.09760116235087592*fUpwindQuad_r[8]+0.1952023247017519*fUpwindQuad_r[7]-0.09760116235087592*fUpwindQuad_r[6]+0.09760116235087592*fUpwindQuad_r[2]-0.1952023247017519*fUpwindQuad_r[1]+0.09760116235087592*fUpwindQuad_r[0]; 
  fUpwind_r[18] = 0.09760116235087592*fUpwindQuad_r[26]-0.09760116235087592*fUpwindQuad_r[24]-0.1952023247017519*fUpwindQuad_r[23]+0.1952023247017519*fUpwindQuad_r[21]+0.09760116235087592*fUpwindQuad_r[20]-0.09760116235087592*fUpwindQuad_r[18]-0.09760116235087592*fUpwindQuad_r[8]+0.09760116235087592*fUpwindQuad_r[6]+0.1952023247017519*fUpwindQuad_r[5]-0.1952023247017519*fUpwindQuad_r[3]-0.09760116235087592*fUpwindQuad_r[2]+0.09760116235087592*fUpwindQuad_r[0]; 
  fUpwind_r[19] = 0.09760116235087592*fUpwindQuad_r[26]-0.09760116235087592*fUpwindQuad_r[24]-0.09760116235087592*fUpwindQuad_r[20]+0.09760116235087592*fUpwindQuad_r[18]-0.1952023247017519*fUpwindQuad_r[17]+0.1952023247017519*fUpwindQuad_r[15]+0.1952023247017519*fUpwindQuad_r[11]-0.1952023247017519*fUpwindQuad_r[9]+0.09760116235087592*fUpwindQuad_r[8]-0.09760116235087592*fUpwindQuad_r[6]-0.09760116235087592*fUpwindQuad_r[2]+0.09760116235087592*fUpwindQuad_r[0]; 

  Gdrag_l[0] = 0.3535533905932737*alphaDrSurf_l[12]*fUpwind_l[12]+0.3535533905932737*alphaDrSurf_l[11]*fUpwind_l[11]+0.3535533905932737*alphaDrSurf_l[8]*fUpwind_l[8]+0.3535533905932737*alphaDrSurf_l[7]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[2]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.3535533905932737*alphaDrSurf_l[8]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[8]*alphaDrSurf_l[12]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[4]*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[1]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[4]*alphaDrSurf_l[12]+0.3535533905932737*alphaDrSurf_l[7]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[7]*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[2]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[2]*alphaDrSurf_l[8]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Gdrag_l[3] = 0.3535533905932737*alphaDrSurf_l[12]*fUpwind_l[18]+0.3535533905932737*alphaDrSurf_l[11]*fUpwind_l[17]+0.3535533905932737*alphaDrSurf_l[8]*fUpwind_l[14]+0.3535533905932737*alphaDrSurf_l[7]*fUpwind_l[13]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[10]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Gdrag_l[4] = 0.2828427124746191*alphaDrSurf_l[11]*fUpwind_l[12]+0.3162277660168379*alphaDrSurf_l[2]*fUpwind_l[12]+0.2828427124746191*fUpwind_l[11]*alphaDrSurf_l[12]+0.3162277660168379*fUpwind_l[2]*alphaDrSurf_l[12]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[1]*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[4]*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[4]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Gdrag_l[5] = 0.3535533905932737*alphaDrSurf_l[8]*fUpwind_l[18]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[17]+0.3535533905932737*alphaDrSurf_l[12]*fUpwind_l[14]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[10]*alphaDrSurf_l[11]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[5]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Gdrag_l[6] = 0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[18]+0.3535533905932737*alphaDrSurf_l[7]*fUpwind_l[17]+0.3162277660168379*alphaDrSurf_l[2]*fUpwind_l[14]+0.3535533905932737*alphaDrSurf_l[11]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[10]*alphaDrSurf_l[12]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[6]*alphaDrSurf_l[8]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[3]; 
  Gdrag_l[7] = 0.3162277660168379*alphaDrSurf_l[12]*fUpwind_l[12]+0.2258769757263128*alphaDrSurf_l[11]*fUpwind_l[11]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[11]+0.2258769757263128*alphaDrSurf_l[7]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[7]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[4]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Gdrag_l[8] = 0.2258769757263128*alphaDrSurf_l[12]*fUpwind_l[12]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[12]+0.3162277660168379*alphaDrSurf_l[11]*fUpwind_l[11]+0.2258769757263128*alphaDrSurf_l[8]*fUpwind_l[8]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[4]+0.3162277660168379*alphaDrSurf_l[2]*fUpwind_l[2]; 
  Gdrag_l[9] = 0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[19]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[16]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[15]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[9]; 
  Gdrag_l[10] = 0.282842712474619*alphaDrSurf_l[11]*fUpwind_l[18]+0.3162277660168379*alphaDrSurf_l[2]*fUpwind_l[18]+0.282842712474619*alphaDrSurf_l[12]*fUpwind_l[17]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[17]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[14]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[6]*alphaDrSurf_l[12]+0.3162277660168379*fUpwind_l[5]*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]*fUpwind_l[10]+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[10]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[10]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[4]; 
  Gdrag_l[11] = 0.2828427124746191*alphaDrSurf_l[4]*fUpwind_l[12]+0.2828427124746191*fUpwind_l[4]*alphaDrSurf_l[12]+0.3162277660168379*alphaDrSurf_l[8]*fUpwind_l[11]+0.2258769757263128*alphaDrSurf_l[7]*fUpwind_l[11]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[8]*alphaDrSurf_l[11]+0.2258769757263128*fUpwind_l[7]*alphaDrSurf_l[11]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[11]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[7]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[1]*alphaDrSurf_l[4]; 
  Gdrag_l[12] = 0.2258769757263128*alphaDrSurf_l[8]*fUpwind_l[12]+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[12]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[12]+0.2258769757263128*fUpwind_l[8]*alphaDrSurf_l[12]+0.3162277660168379*fUpwind_l[7]*alphaDrSurf_l[12]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[12]+0.2828427124746191*alphaDrSurf_l[4]*fUpwind_l[11]+0.2828427124746191*fUpwind_l[4]*alphaDrSurf_l[11]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[2]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[2]*alphaDrSurf_l[4]; 
  Gdrag_l[13] = 0.3162277660168379*alphaDrSurf_l[12]*fUpwind_l[18]+0.2258769757263128*alphaDrSurf_l[11]*fUpwind_l[17]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[17]+0.2258769757263128*alphaDrSurf_l[7]*fUpwind_l[13]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[6]*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[7]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[5]; 
  Gdrag_l[14] = 0.2258769757263128*alphaDrSurf_l[12]*fUpwind_l[18]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[18]+0.3162277660168379*alphaDrSurf_l[11]*fUpwind_l[17]+0.2258769757263128*alphaDrSurf_l[8]*fUpwind_l[14]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[5]*alphaDrSurf_l[12]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[2]*fUpwind_l[6]; 
  Gdrag_l[15] = 0.3162277660168379*alphaDrSurf_l[11]*fUpwind_l[19]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[19]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[16]+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[15]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[15]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[9]; 
  Gdrag_l[16] = 0.3162277660168379*alphaDrSurf_l[12]*fUpwind_l[19]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[19]+0.3162277660168379*alphaDrSurf_l[8]*fUpwind_l[16]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[16]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[15]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[9]; 
  Gdrag_l[17] = 0.2828427124746191*alphaDrSurf_l[4]*fUpwind_l[18]+0.3162277660168379*alphaDrSurf_l[8]*fUpwind_l[17]+0.2258769757263128*alphaDrSurf_l[7]*fUpwind_l[17]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[17]+0.3162277660168379*alphaDrSurf_l[11]*fUpwind_l[14]+0.2258769757263128*alphaDrSurf_l[11]*fUpwind_l[13]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[13]+0.282842712474619*fUpwind_l[10]*alphaDrSurf_l[12]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[6]*alphaDrSurf_l[7]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[5]; 
  Gdrag_l[18] = 0.2258769757263128*alphaDrSurf_l[8]*fUpwind_l[18]+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[18]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[18]+0.2828427124746191*alphaDrSurf_l[4]*fUpwind_l[17]+0.2258769757263128*alphaDrSurf_l[12]*fUpwind_l[14]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[14]+0.3162277660168379*alphaDrSurf_l[12]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[12]+0.282842712474619*fUpwind_l[10]*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[2]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[5]*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[4]*fUpwind_l[6]; 
  Gdrag_l[19] = 0.3162277660168379*alphaDrSurf_l[8]*fUpwind_l[19]+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[19]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[19]+0.3162277660168379*alphaDrSurf_l[12]*fUpwind_l[16]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[16]+0.3162277660168379*alphaDrSurf_l[11]*fUpwind_l[15]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[15]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[9]; 

  Gdrag_r[0] = 0.3535533905932737*alphaDrSurf_r[12]*fUpwind_r[12]+0.3535533905932737*alphaDrSurf_r[11]*fUpwind_r[11]+0.3535533905932737*alphaDrSurf_r[8]*fUpwind_r[8]+0.3535533905932737*alphaDrSurf_r[7]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[2]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[1]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.3535533905932737*alphaDrSurf_r[8]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[8]*alphaDrSurf_r[12]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[4]*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[1]*alphaDrSurf_r[7]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[4]*alphaDrSurf_r[12]+0.3535533905932737*alphaDrSurf_r[7]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[7]*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[2]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[2]*alphaDrSurf_r[8]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Gdrag_r[3] = 0.3535533905932737*alphaDrSurf_r[12]*fUpwind_r[18]+0.3535533905932737*alphaDrSurf_r[11]*fUpwind_r[17]+0.3535533905932737*alphaDrSurf_r[8]*fUpwind_r[14]+0.3535533905932737*alphaDrSurf_r[7]*fUpwind_r[13]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[10]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Gdrag_r[4] = 0.2828427124746191*alphaDrSurf_r[11]*fUpwind_r[12]+0.3162277660168379*alphaDrSurf_r[2]*fUpwind_r[12]+0.2828427124746191*fUpwind_r[11]*alphaDrSurf_r[12]+0.3162277660168379*fUpwind_r[2]*alphaDrSurf_r[12]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[1]*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[4]*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[4]*alphaDrSurf_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Gdrag_r[5] = 0.3535533905932737*alphaDrSurf_r[8]*fUpwind_r[18]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[17]+0.3535533905932737*alphaDrSurf_r[12]*fUpwind_r[14]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[10]*alphaDrSurf_r[11]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[5]*alphaDrSurf_r[7]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Gdrag_r[6] = 0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[18]+0.3535533905932737*alphaDrSurf_r[7]*fUpwind_r[17]+0.3162277660168379*alphaDrSurf_r[2]*fUpwind_r[14]+0.3535533905932737*alphaDrSurf_r[11]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[10]*alphaDrSurf_r[12]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[6]*alphaDrSurf_r[8]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[3]; 
  Gdrag_r[7] = 0.3162277660168379*alphaDrSurf_r[12]*fUpwind_r[12]+0.2258769757263128*alphaDrSurf_r[11]*fUpwind_r[11]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[2]*alphaDrSurf_r[11]+0.2258769757263128*alphaDrSurf_r[7]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[7]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[4]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Gdrag_r[8] = 0.2258769757263128*alphaDrSurf_r[12]*fUpwind_r[12]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[1]*alphaDrSurf_r[12]+0.3162277660168379*alphaDrSurf_r[11]*fUpwind_r[11]+0.2258769757263128*alphaDrSurf_r[8]*fUpwind_r[8]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[4]+0.3162277660168379*alphaDrSurf_r[2]*fUpwind_r[2]; 
  Gdrag_r[9] = 0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[19]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[16]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[15]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[9]; 
  Gdrag_r[10] = 0.282842712474619*alphaDrSurf_r[11]*fUpwind_r[18]+0.3162277660168379*alphaDrSurf_r[2]*fUpwind_r[18]+0.282842712474619*alphaDrSurf_r[12]*fUpwind_r[17]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[17]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[14]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[6]*alphaDrSurf_r[12]+0.3162277660168379*fUpwind_r[5]*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[8]*fUpwind_r[10]+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[10]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[10]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alphaDrSurf_r[4]; 
  Gdrag_r[11] = 0.2828427124746191*alphaDrSurf_r[4]*fUpwind_r[12]+0.2828427124746191*fUpwind_r[4]*alphaDrSurf_r[12]+0.3162277660168379*alphaDrSurf_r[8]*fUpwind_r[11]+0.2258769757263128*alphaDrSurf_r[7]*fUpwind_r[11]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[8]*alphaDrSurf_r[11]+0.2258769757263128*fUpwind_r[7]*alphaDrSurf_r[11]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[11]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alphaDrSurf_r[7]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[1]*alphaDrSurf_r[4]; 
  Gdrag_r[12] = 0.2258769757263128*alphaDrSurf_r[8]*fUpwind_r[12]+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[12]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[12]+0.2258769757263128*fUpwind_r[8]*alphaDrSurf_r[12]+0.3162277660168379*fUpwind_r[7]*alphaDrSurf_r[12]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[12]+0.2828427124746191*alphaDrSurf_r[4]*fUpwind_r[11]+0.2828427124746191*fUpwind_r[4]*alphaDrSurf_r[11]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[1]*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[2]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[2]*alphaDrSurf_r[4]; 
  Gdrag_r[13] = 0.3162277660168379*alphaDrSurf_r[12]*fUpwind_r[18]+0.2258769757263128*alphaDrSurf_r[11]*fUpwind_r[17]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[17]+0.2258769757263128*alphaDrSurf_r[7]*fUpwind_r[13]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[6]*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[3]*alphaDrSurf_r[7]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[5]; 
  Gdrag_r[14] = 0.2258769757263128*alphaDrSurf_r[12]*fUpwind_r[18]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[18]+0.3162277660168379*alphaDrSurf_r[11]*fUpwind_r[17]+0.2258769757263128*alphaDrSurf_r[8]*fUpwind_r[14]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[5]*alphaDrSurf_r[12]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[3]*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[2]*fUpwind_r[6]; 
  Gdrag_r[15] = 0.3162277660168379*alphaDrSurf_r[11]*fUpwind_r[19]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[19]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[16]+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[15]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[15]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[9]; 
  Gdrag_r[16] = 0.3162277660168379*alphaDrSurf_r[12]*fUpwind_r[19]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[19]+0.3162277660168379*alphaDrSurf_r[8]*fUpwind_r[16]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[16]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[15]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[9]; 
  Gdrag_r[17] = 0.2828427124746191*alphaDrSurf_r[4]*fUpwind_r[18]+0.3162277660168379*alphaDrSurf_r[8]*fUpwind_r[17]+0.2258769757263128*alphaDrSurf_r[7]*fUpwind_r[17]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[17]+0.3162277660168379*alphaDrSurf_r[11]*fUpwind_r[14]+0.2258769757263128*alphaDrSurf_r[11]*fUpwind_r[13]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[13]+0.282842712474619*fUpwind_r[10]*alphaDrSurf_r[12]+0.3535533905932737*fUpwind_r[3]*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[6]*alphaDrSurf_r[7]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[5]; 
  Gdrag_r[18] = 0.2258769757263128*alphaDrSurf_r[8]*fUpwind_r[18]+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[18]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[18]+0.2828427124746191*alphaDrSurf_r[4]*fUpwind_r[17]+0.2258769757263128*alphaDrSurf_r[12]*fUpwind_r[14]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[14]+0.3162277660168379*alphaDrSurf_r[12]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[3]*alphaDrSurf_r[12]+0.282842712474619*fUpwind_r[10]*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[2]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[5]*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[4]*fUpwind_r[6]; 
  Gdrag_r[19] = 0.3162277660168379*alphaDrSurf_r[8]*fUpwind_r[19]+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[19]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[19]+0.3162277660168379*alphaDrSurf_r[12]*fUpwind_r[16]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[16]+0.3162277660168379*alphaDrSurf_r[11]*fUpwind_r[15]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[15]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[9]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[4] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[8] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[9] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[10] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[12] += 0.7071067811865475*Gdrag_r[8]*rdv2-0.7071067811865475*Gdrag_l[8]*rdv2; 
  out[13] += 0.7071067811865475*Gdrag_r[9]*rdv2-0.7071067811865475*Gdrag_l[9]*rdv2; 
  out[14] += 1.58113883008419*Gdrag_r[0]*rdv2-1.58113883008419*Gdrag_l[0]*rdv2; 
  out[15] += 0.7071067811865475*Gdrag_r[10]*rdv2-0.7071067811865475*Gdrag_l[10]*rdv2; 
  out[16] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[17] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[18] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[19] += 0.7071067811865475*Gdrag_r[11]*rdv2-0.7071067811865475*Gdrag_l[11]*rdv2; 
  out[20] += 0.7071067811865475*Gdrag_r[12]*rdv2-0.7071067811865475*Gdrag_l[12]*rdv2; 
  out[21] += 0.7071067811865475*Gdrag_r[13]*rdv2-0.7071067811865475*Gdrag_l[13]*rdv2; 
  out[22] += 0.7071067811865475*Gdrag_r[14]*rdv2-0.7071067811865475*Gdrag_l[14]*rdv2; 
  out[23] += 0.7071067811865475*Gdrag_r[15]*rdv2-0.7071067811865475*Gdrag_l[15]*rdv2; 
  out[24] += 0.7071067811865475*Gdrag_r[16]*rdv2-0.7071067811865475*Gdrag_l[16]*rdv2; 
  out[25] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
  out[26] += 1.224744871391589*Gdrag_r[8]*rdv2+1.224744871391589*Gdrag_l[8]*rdv2; 
  out[27] += 1.224744871391589*Gdrag_r[9]*rdv2+1.224744871391589*Gdrag_l[9]*rdv2; 
  out[28] += 1.58113883008419*Gdrag_r[1]*rdv2-1.58113883008419*Gdrag_l[1]*rdv2; 
  out[29] += 1.58113883008419*Gdrag_r[2]*rdv2-1.58113883008419*Gdrag_l[2]*rdv2; 
  out[30] += 1.58113883008419*Gdrag_r[3]*rdv2-1.58113883008419*Gdrag_l[3]*rdv2; 
  out[31] += 1.224744871391589*Gdrag_r[10]*rdv2+1.224744871391589*Gdrag_l[10]*rdv2; 
  out[32] += 0.7071067811865475*Gdrag_r[17]*rdv2-0.7071067811865475*Gdrag_l[17]*rdv2; 
  out[33] += 0.7071067811865475*Gdrag_r[18]*rdv2-0.7071067811865475*Gdrag_l[18]*rdv2; 
  out[34] += 0.7071067811865475*Gdrag_r[19]*rdv2-0.7071067811865475*Gdrag_l[19]*rdv2; 
  out[35] += 1.224744871391589*Gdrag_r[11]*rdv2+1.224744871391589*Gdrag_l[11]*rdv2; 
  out[36] += 1.224744871391589*Gdrag_r[12]*rdv2+1.224744871391589*Gdrag_l[12]*rdv2; 
  out[37] += 1.224744871391589*Gdrag_r[13]*rdv2+1.224744871391589*Gdrag_l[13]*rdv2; 
  out[38] += 1.224744871391589*Gdrag_r[14]*rdv2+1.224744871391589*Gdrag_l[14]*rdv2; 
  out[39] += 1.224744871391589*Gdrag_r[15]*rdv2+1.224744871391589*Gdrag_l[15]*rdv2; 
  out[40] += 1.224744871391589*Gdrag_r[16]*rdv2+1.224744871391589*Gdrag_l[16]*rdv2; 
  out[41] += 1.58113883008419*Gdrag_r[4]*rdv2-1.58113883008419*Gdrag_l[4]*rdv2; 
  out[42] += 1.58113883008419*Gdrag_r[5]*rdv2-1.58113883008419*Gdrag_l[5]*rdv2; 
  out[43] += 1.58113883008419*Gdrag_r[6]*rdv2-1.58113883008419*Gdrag_l[6]*rdv2; 
  out[44] += 1.224744871391589*Gdrag_r[17]*rdv2+1.224744871391589*Gdrag_l[17]*rdv2; 
  out[45] += 1.224744871391589*Gdrag_r[18]*rdv2+1.224744871391589*Gdrag_l[18]*rdv2; 
  out[46] += 1.224744871391589*Gdrag_r[19]*rdv2+1.224744871391589*Gdrag_l[19]*rdv2; 
  out[47] += 1.58113883008419*Gdrag_r[10]*rdv2-1.58113883008419*Gdrag_l[10]*rdv2; 
} 
