#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_tensor_2x2v_p2_surfvy_quad.h> 
GKYL_CU_DH void vlasov_lbo_surfvy_2x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         cell-center coordinates. 
  // dxv[4]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[18]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[9]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  const double *sumNuUy = &nuUSum[9]; 

  double alphaDrSurf_l[27] = {0.0}; 
  alphaDrSurf_l[0] = 2.828427124746191*nuSum[0]*w[3]-1.414213562373095*nuSum[0]*dxv[3]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_l[1] = -1.414213562373095*sumNuUy[1]; 
  alphaDrSurf_l[2] = -1.414213562373095*sumNuUy[2]; 
  alphaDrSurf_l[4] = -1.414213562373095*sumNuUy[3]; 
  alphaDrSurf_l[7] = -1.414213562373095*sumNuUy[4]; 
  alphaDrSurf_l[8] = -1.414213562373095*sumNuUy[5]; 
  alphaDrSurf_l[11] = -1.414213562373095*sumNuUy[6]; 
  alphaDrSurf_l[12] = -1.414213562373095*sumNuUy[7]; 
  alphaDrSurf_l[20] = -1.414213562373095*sumNuUy[8]; 

  double alphaDrSurf_r[27] = {0.0}; 
  alphaDrSurf_r[0] = 2.828427124746191*nuSum[0]*w[3]+1.414213562373095*nuSum[0]*dxv[3]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_r[1] = -1.414213562373095*sumNuUy[1]; 
  alphaDrSurf_r[2] = -1.414213562373095*sumNuUy[2]; 
  alphaDrSurf_r[4] = -1.414213562373095*sumNuUy[3]; 
  alphaDrSurf_r[7] = -1.414213562373095*sumNuUy[4]; 
  alphaDrSurf_r[8] = -1.414213562373095*sumNuUy[5]; 
  alphaDrSurf_r[11] = -1.414213562373095*sumNuUy[6]; 
  alphaDrSurf_r[12] = -1.414213562373095*sumNuUy[7]; 
  alphaDrSurf_r[20] = -1.414213562373095*sumNuUy[8]; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[27] = {0.0};;
  double fUpwind_r[27] = {0.0};
  double Ghat_l[27] = {0.0}; 
  double Ghat_r[27] = {0.0}; 
  double Gdiff_l[27] = {0.0}; 
  double Gdiff_r[27] = {0.0}; 
  double Gdiff2_l[27] = {0.0}; 
  double Gdiff2_r[27] = {0.0}; 

  if (0.2828427124746191*alphaDrSurf_l[20]-0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = tensor_2x2v_p2_surfvy_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = tensor_2x2v_p2_surfvy_quad_0(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])+0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = tensor_2x2v_p2_surfvy_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = tensor_2x2v_p2_surfvy_quad_1(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]+0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*alphaDrSurf_l[2]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = tensor_2x2v_p2_surfvy_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = tensor_2x2v_p2_surfvy_quad_2(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])+0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = tensor_2x2v_p2_surfvy_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = tensor_2x2v_p2_surfvy_quad_3(-1, fc); 
  } 
  if (0.441941738241592*alphaDrSurf_l[20]-0.3952847075210473*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[4] = tensor_2x2v_p2_surfvy_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = tensor_2x2v_p2_surfvy_quad_4(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])-0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = tensor_2x2v_p2_surfvy_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = tensor_2x2v_p2_surfvy_quad_5(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]-0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*alphaDrSurf_l[2]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[6] = tensor_2x2v_p2_surfvy_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = tensor_2x2v_p2_surfvy_quad_6(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])-0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = tensor_2x2v_p2_surfvy_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = tensor_2x2v_p2_surfvy_quad_7(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]+0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[8] = tensor_2x2v_p2_surfvy_quad_8(1, fl); 
  } else { 
    fUpwindQuad_l[8] = tensor_2x2v_p2_surfvy_quad_8(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]-0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[9] = tensor_2x2v_p2_surfvy_quad_9(1, fl); 
  } else { 
    fUpwindQuad_l[9] = tensor_2x2v_p2_surfvy_quad_9(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])+0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[10] = tensor_2x2v_p2_surfvy_quad_10(1, fl); 
  } else { 
    fUpwindQuad_l[10] = tensor_2x2v_p2_surfvy_quad_10(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]+0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*alphaDrSurf_l[2]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[11] = tensor_2x2v_p2_surfvy_quad_11(1, fl); 
  } else { 
    fUpwindQuad_l[11] = tensor_2x2v_p2_surfvy_quad_11(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])+0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[12] = tensor_2x2v_p2_surfvy_quad_12(1, fl); 
  } else { 
    fUpwindQuad_l[12] = tensor_2x2v_p2_surfvy_quad_12(-1, fc); 
  } 
  if (0.441941738241592*alphaDrSurf_l[20]-0.3952847075210473*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[13] = tensor_2x2v_p2_surfvy_quad_13(1, fl); 
  } else { 
    fUpwindQuad_l[13] = tensor_2x2v_p2_surfvy_quad_13(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])-0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[14] = tensor_2x2v_p2_surfvy_quad_14(1, fl); 
  } else { 
    fUpwindQuad_l[14] = tensor_2x2v_p2_surfvy_quad_14(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]-0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*alphaDrSurf_l[2]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[15] = tensor_2x2v_p2_surfvy_quad_15(1, fl); 
  } else { 
    fUpwindQuad_l[15] = tensor_2x2v_p2_surfvy_quad_15(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])-0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[16] = tensor_2x2v_p2_surfvy_quad_16(1, fl); 
  } else { 
    fUpwindQuad_l[16] = tensor_2x2v_p2_surfvy_quad_16(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]+0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[17] = tensor_2x2v_p2_surfvy_quad_17(1, fl); 
  } else { 
    fUpwindQuad_l[17] = tensor_2x2v_p2_surfvy_quad_17(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]-0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[18] = tensor_2x2v_p2_surfvy_quad_18(1, fl); 
  } else { 
    fUpwindQuad_l[18] = tensor_2x2v_p2_surfvy_quad_18(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])+0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[19] = tensor_2x2v_p2_surfvy_quad_19(1, fl); 
  } else { 
    fUpwindQuad_l[19] = tensor_2x2v_p2_surfvy_quad_19(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]+0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*alphaDrSurf_l[2]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[20] = tensor_2x2v_p2_surfvy_quad_20(1, fl); 
  } else { 
    fUpwindQuad_l[20] = tensor_2x2v_p2_surfvy_quad_20(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])+0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[21] = tensor_2x2v_p2_surfvy_quad_21(1, fl); 
  } else { 
    fUpwindQuad_l[21] = tensor_2x2v_p2_surfvy_quad_21(-1, fc); 
  } 
  if (0.441941738241592*alphaDrSurf_l[20]-0.3952847075210473*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[22] = tensor_2x2v_p2_surfvy_quad_22(1, fl); 
  } else { 
    fUpwindQuad_l[22] = tensor_2x2v_p2_surfvy_quad_22(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])-0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[23] = tensor_2x2v_p2_surfvy_quad_23(1, fl); 
  } else { 
    fUpwindQuad_l[23] = tensor_2x2v_p2_surfvy_quad_23(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]-0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*alphaDrSurf_l[2]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[24] = tensor_2x2v_p2_surfvy_quad_24(1, fl); 
  } else { 
    fUpwindQuad_l[24] = tensor_2x2v_p2_surfvy_quad_24(-1, fc); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_l[20])-0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[25] = tensor_2x2v_p2_surfvy_quad_25(1, fl); 
  } else { 
    fUpwindQuad_l[25] = tensor_2x2v_p2_surfvy_quad_25(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_l[20]+0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[26] = tensor_2x2v_p2_surfvy_quad_26(1, fl); 
  } else { 
    fUpwindQuad_l[26] = tensor_2x2v_p2_surfvy_quad_26(-1, fc); 
  } 
  if (0.2828427124746191*alphaDrSurf_r[20]-0.4242640687119281*alphaDrSurf_r[12]-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = tensor_2x2v_p2_surfvy_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = tensor_2x2v_p2_surfvy_quad_0(-1, fr); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_r[20])+0.5303300858899104*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = tensor_2x2v_p2_surfvy_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = tensor_2x2v_p2_surfvy_quad_1(-1, fr); 
  } 
  if (0.2828427124746191*alphaDrSurf_r[20]+0.4242640687119281*alphaDrSurf_r[12]-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*alphaDrSurf_r[2]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = tensor_2x2v_p2_surfvy_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = tensor_2x2v_p2_surfvy_quad_2(-1, fr); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_r[20])+0.5303300858899104*alphaDrSurf_r[12]-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = tensor_2x2v_p2_surfvy_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = tensor_2x2v_p2_surfvy_quad_3(-1, fr); 
  } 
  if (0.441941738241592*alphaDrSurf_r[20]-0.3952847075210473*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[4] = tensor_2x2v_p2_surfvy_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = tensor_2x2v_p2_surfvy_quad_4(-1, fr); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_r[20])-0.5303300858899104*alphaDrSurf_r[12]-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = tensor_2x2v_p2_surfvy_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = tensor_2x2v_p2_surfvy_quad_5(-1, fr); 
  } 
  if (0.2828427124746191*alphaDrSurf_r[20]-0.4242640687119281*alphaDrSurf_r[12]+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*alphaDrSurf_r[2]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[6] = tensor_2x2v_p2_surfvy_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = tensor_2x2v_p2_surfvy_quad_6(-1, fr); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_r[20])-0.5303300858899104*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = tensor_2x2v_p2_surfvy_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = tensor_2x2v_p2_surfvy_quad_7(-1, fr); 
  } 
  if (0.2828427124746191*alphaDrSurf_r[20]+0.4242640687119281*alphaDrSurf_r[12]+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[8] = tensor_2x2v_p2_surfvy_quad_8(1, fc); 
  } else { 
    fUpwindQuad_r[8] = tensor_2x2v_p2_surfvy_quad_8(-1, fr); 
  } 
  if (0.2828427124746191*alphaDrSurf_r[20]-0.4242640687119281*alphaDrSurf_r[12]-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[9] = tensor_2x2v_p2_surfvy_quad_9(1, fc); 
  } else { 
    fUpwindQuad_r[9] = tensor_2x2v_p2_surfvy_quad_9(-1, fr); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_r[20])+0.5303300858899104*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[10] = tensor_2x2v_p2_surfvy_quad_10(1, fc); 
  } else { 
    fUpwindQuad_r[10] = tensor_2x2v_p2_surfvy_quad_10(-1, fr); 
  } 
  if (0.2828427124746191*alphaDrSurf_r[20]+0.4242640687119281*alphaDrSurf_r[12]-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*alphaDrSurf_r[2]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[11] = tensor_2x2v_p2_surfvy_quad_11(1, fc); 
  } else { 
    fUpwindQuad_r[11] = tensor_2x2v_p2_surfvy_quad_11(-1, fr); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_r[20])+0.5303300858899104*alphaDrSurf_r[12]-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[12] = tensor_2x2v_p2_surfvy_quad_12(1, fc); 
  } else { 
    fUpwindQuad_r[12] = tensor_2x2v_p2_surfvy_quad_12(-1, fr); 
  } 
  if (0.441941738241592*alphaDrSurf_r[20]-0.3952847075210473*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[13] = tensor_2x2v_p2_surfvy_quad_13(1, fc); 
  } else { 
    fUpwindQuad_r[13] = tensor_2x2v_p2_surfvy_quad_13(-1, fr); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_r[20])-0.5303300858899104*alphaDrSurf_r[12]-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[14] = tensor_2x2v_p2_surfvy_quad_14(1, fc); 
  } else { 
    fUpwindQuad_r[14] = tensor_2x2v_p2_surfvy_quad_14(-1, fr); 
  } 
  if (0.2828427124746191*alphaDrSurf_r[20]-0.4242640687119281*alphaDrSurf_r[12]+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*alphaDrSurf_r[2]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[15] = tensor_2x2v_p2_surfvy_quad_15(1, fc); 
  } else { 
    fUpwindQuad_r[15] = tensor_2x2v_p2_surfvy_quad_15(-1, fr); 
  } 
  if ((-0.3535533905932737*alphaDrSurf_r[20])-0.5303300858899104*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[16] = tensor_2x2v_p2_surfvy_quad_16(1, fc); 
  } else { 
    fUpwindQuad_r[16] = tensor_2x2v_p2_surfvy_quad_16(-1, fr); 
  } 
  if (0.2828427124746191*alphaDrSurf_r[20]+0.4242640687119281*alphaDrSurf_r[12]+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[17] = tensor_2x2v_p2_surfvy_quad_17(1, fc); 
  } else { 
    fUpwindQuad_r[17] = tensor_2x2v