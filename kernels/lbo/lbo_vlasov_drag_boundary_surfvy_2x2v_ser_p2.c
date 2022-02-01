#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x2v_p2_surfvy_quad.h> 
GKYL_CU_DH void lbo_vlasov_drag_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[16]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[8]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *sumNuUy = &nuUSum[8]; 

  double alphaDrSurf[20] = {0.0}; 
  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};;
  double drag_incr[20] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]); 
  alphaDrSurf[7] = -0.7071067811865475*(2.0*sumNuUy[4]+((-2.0*w[3])-1.0*dxv[3])*nuSum[4]); 
  alphaDrSurf[8] = -0.7071067811865475*(2.0*sumNuUy[5]+((-2.0*w[3])-1.0*dxv[3])*nuSum[5]); 
  alphaDrSurf[11] = -0.7071067811865475*(2.0*sumNuUy[6]+((-2.0*w[3])-1.0*dxv[3])*nuSum[6]); 
  alphaDrSurf[12] = -0.7071067811865475*(2.0*sumNuUy[7]+((-2.0*w[3])-1.0*dxv[3])*nuSum[7]); 

  if ((-0.4242640687119281*alphaDrSurf[12])-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_2x2v_p2_surfvy_quad_0(1, fSkin); 
  } else { 
    fUpwindQuad[0] = ser_2x2v_p2_surfvy_quad_0(-1, fEdge); 
  } 
  if (0.5303300858899104*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x2v_p2_surfvy_quad_1(1, fSkin); 
  } else { 
    fUpwindQuad[1] = ser_2x2v_p2_surfvy_quad_1(-1, fEdge); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_2x2v_p2_surfvy_quad_2(1, fSkin); 
  } else { 
    fUpwindQuad[2] = ser_2x2v_p2_surfvy_quad_2(-1, fEdge); 
  } 
  if (0.5303300858899104*alphaDrSurf[12]-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_2x2v_p2_surfvy_quad_3(1, fSkin); 
  } else { 
    fUpwindQuad[3] = ser_2x2v_p2_surfvy_quad_3(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*(alphaDrSurf[8]+alphaDrSurf[7]) < 0) { 
    fUpwindQuad[4] = ser_2x2v_p2_surfvy_quad_4(1, fSkin); 
  } else { 
    fUpwindQuad[4] = ser_2x2v_p2_surfvy_quad_4(-1, fEdge); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[12])-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_2x2v_p2_surfvy_quad_5(1, fSkin); 
  } else { 
    fUpwindQuad[5] = ser_2x2v_p2_surfvy_quad_5(-1, fEdge); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_2x2v_p2_surfvy_quad_6(1, fSkin); 
  } else { 
    fUpwindQuad[6] = ser_2x2v_p2_surfvy_quad_6(-1, fEdge); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[11])+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_2x2v_p2_surfvy_quad_7(1, fSkin); 
  } else { 
    fUpwindQuad[7] = ser_2x2v_p2_surfvy_quad_7(-1, fEdge); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_2x2v_p2_surfvy_quad_8(1, fSkin); 
  } else { 
    fUpwindQuad[8] = ser_2x2v_p2_surfvy_quad_8(-1, fEdge); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_2x2v_p2_surfvy_quad_9(1, fSkin); 
  } else { 
    fUpwindQuad[9] = ser_2x2v_p2_surfvy_quad_9(-1, fEdge); 
  } 
  if (0.5303300858899104*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = ser_2x2v_p2_surfvy_quad_10(1, fSkin); 
  } else { 
    fUpwindQuad[10] = ser_2x2v_p2_surfvy_quad_10(-1, fEdge); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_2x2v_p2_surfvy_quad_11(1, fSkin); 
  } else { 
    fUpwindQuad[11] = ser_2x2v_p2_surfvy_quad_11(-1, fEdge); 
  } 
  if (0.5303300858899104*alphaDrSurf[12]-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_2x2v_p2_surfvy_quad_12(1, fSkin); 
  } else { 
    fUpwindQuad[12] = ser_2x2v_p2_surfvy_quad_12(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*(alphaDrSurf[8]+alphaDrSurf[7]) < 0) { 
    fUpwindQuad[13] = ser_2x2v_p2_surfvy_quad_13(1, fSkin); 
  } else { 
    fUpwindQuad[13] = ser_2x2v_p2_surfvy_quad_13(-1, fEdge); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[12])-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_2x2v_p2_surfvy_quad_14(1, fSkin); 
  } else { 
    fUpwindQuad[14] = ser_2x2v_p2_surfvy_quad_14(-1, fEdge); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_2x2v_p2_surfvy_quad_15(1, fSkin); 
  } else { 
    fUpwindQuad[15] = ser_2x2v_p2_surfvy_quad_15(-1, fEdge); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[11])+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[16] = ser_2x2v_p2_surfvy_quad_16(1, fSkin); 
  } else { 
    fUpwindQuad[16] = ser_2x2v_p2_surfvy_quad_16(-1, fEdge); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[17] = ser_2x2v_p2_surfvy_quad_17(1, fSkin); 
  } else { 
    fUpwindQuad[17] = ser_2x2v_p2_surfvy_quad_17(-1, fEdge); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_2x2v_p2_surfvy_quad_18(1, fSkin); 
  } else { 
    fUpwindQuad[18] = ser_2x2v_p2_surfvy_quad_18(-1, fEdge); 
  } 
  if (0.5303300858899104*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[19] = ser_2x2v_p2_surfvy_quad_19(1, fSkin); 
  } else { 
    fUpwindQuad[19] = ser_2x2v_p2_surfvy_quad_19(-1, fEdge); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[20] = ser_2x2v_p2_surfvy_quad_20(1, fSkin); 
  } else { 
    fUpwindQuad[20] = ser_2x2v_p2_surfvy_quad_20(-1, fEdge); 
  } 
  if (0.5303300858899104*alphaDrSurf[12]-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_2x2v_p2_surfvy_quad_21(1, fSkin); 
  } else { 
    fUpwindQuad[21] = ser_2x2v_p2_surfvy_quad_21(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*(alphaDrSurf[8]+alphaDrSurf[7]) < 0) { 
    fUpwindQuad[22] = ser_2x2v_p2_surfvy_quad_22(1, fSkin); 
  } else { 
    fUpwindQuad[22] = ser_2x2v_p2_surfvy_quad_22(-1, fEdge); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[12])-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[23] = ser_2x2v_p2_surfvy_quad_23(1, fSkin); 
  } else { 
    fUpwindQuad[23] = ser_2x2v_p2_surfvy_quad_23(-1, fEdge); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_2x2v_p2_surfvy_quad_24(1, fSkin); 
  } else { 
    fUpwindQuad[24] = ser_2x2v_p2_surfvy_quad_24(-1, fEdge); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[11])+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[25] = ser_2x2v_p2_surfvy_quad_25(1, fSkin); 
  } else { 
    fUpwindQuad[25] = ser_2x2v_p2_surfvy_quad_25(-1, fEdge); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[26] = ser_2x2v_p2_surfvy_quad_26(1, fSkin); 
  } else { 
    fUpwindQuad[26] = ser_2x2v_p2_surfvy_quad_26(-1, fEdge); 
  } 

  fUpwind[0] = 0.06062300936098657*fUpwindQuad[26]+0.09699681497757856*fUpwindQuad[25]+0.06062300936098657*fUpwindQuad[24]+0.09699681497757856*fUpwindQuad[23]+0.1551949039641257*fUpwindQuad[22]+0.09699681497757856*fUpwindQuad[21]+0.06062300936098657*fUpwindQuad[20]+0.09699681497757856*fUpwindQuad[19]+0.06062300936098657*fUpwindQuad[18]+0.09699681497757856*fUpwindQuad[17]+0.1551949039641257*fUpwindQuad[16]+0.09699681497757856*fUpwindQuad[15]+0.1551949039641257*fUpwindQuad[14]+0.2483118463426013*fUpwindQuad[13]+0.1551949039641257*fUpwindQuad[12]+0.09699681497757856*fUpwindQuad[11]+0.1551949039641257*fUpwindQuad[10]+0.09699681497757856*fUpwindQuad[9]+0.06062300936098657*fUpwindQuad[8]+0.09699681497757856*fUpwindQuad[7]+0.06062300936098657*fUpwindQuad[6]+0.09699681497757856*fUpwindQuad[5]+0.1551949039641257*fUpwindQuad[4]+0.09699681497757856*fUpwindQuad[3]+0.06062300936098657*fUpwindQuad[2]+0.09699681497757856*fUpwindQuad[1]+0.06062300936098657*fUpwindQuad[0]; 
  fUpwind[1] = 0.08133430195906327*fUpwindQuad[26]-0.08133430195906327*fUpwindQuad[24]+0.1301348831345013*fUpwindQuad[23]-0.1301348831345013*fUpwindQuad[21]+0.08133430195906327*fUpwindQuad[20]-0.08133430195906327*fUpwindQuad[18]+0.1301348831345013*fUpwindQuad[17]-0.1301348831345013*fUpwindQuad[15]+0.2082158130152021*fUpwindQuad[14]-0.2082158130152021*fUpwindQuad[12]+0.1301348831345013*fUpwindQuad[11]-0.1301348831345013*fUpwindQuad[9]+0.08133430195906327*fUpwindQuad[8]-0.08133430195906327*fUpwindQuad[6]+0.1301348831345013*fUpwindQuad[5]-0.1301348831345013*fUpwindQuad[3]+0.08133430195906327*fUpwindQuad[2]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[2] = 0.08133430195906327*fUpwindQuad[26]+0.1301348831345013*fUpwindQuad[25]+0.08133430195906327*fUpwindQuad[24]-0.08133430195906327*fUpwindQuad[20]-0.1301348831345013*fUpwindQuad[19]-0.08133430195906327*fUpwindQuad[18]+0.1301348831345013*fUpwindQuad[17]+0.2082158130152021*fUpwindQuad[16]+0.1301348831345013*fUpwindQuad[15]-0.1301348831345013*fUpwindQuad[11]-0.2082158130152021*fUpwindQuad[10]-0.1301348831345013*fUpwindQuad[9]+0.08133430195906327*fUpwindQuad[8]+0.1301348831345013*fUpwindQuad[7]+0.08133430195906327*fUpwindQuad[6]-0.08133430195906327*fUpwindQuad[2]-0.1301348831345013*fUpwindQuad[1]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[3] = 0.08133430195906327*fUpwindQuad[26]+0.1301348831345013*fUpwindQuad[25]+0.08133430195906327*fUpwindQuad[24]+0.1301348831345013*fUpwindQuad[23]+0.2082158130152021*fUpwindQuad[22]+0.1301348831345013*fUpwindQuad[21]+0.08133430195906327*fUpwindQuad[20]+0.1301348831345013*fUpwindQuad[19]+0.08133430195906327*fUpwindQuad[18]-0.08133430195906327*fUpwindQuad[8]-0.1301348831345013*fUpwindQuad[7]-0.08133430195906327*fUpwindQuad[6]-0.1301348831345013*fUpwindQuad[5]-0.2082158130152021*fUpwindQuad[4]-0.1301348831345013*fUpwindQuad[3]-0.08133430195906327*fUpwindQuad[2]-0.1301348831345013*fUpwindQuad[1]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[4] = 0.1091214168497758*fUpwindQuad[26]-0.1091214168497758*fUpwindQuad[24]-0.1091214168497758*fUpwindQuad[20]+0.1091214168497758*fUpwindQuad[18]+0.1745942669596414*fUpwindQuad[17]-0.1745942669596414*fUpwindQuad[15]-0.1745942669596414*fUpwindQuad[11]+0.1745942669596414*fUpwindQuad[9]+0.1091214168497758*fUpwindQuad[8]-0.1091214168497758*fUpwindQuad[6]-0.1091214168497758*fUpwindQuad[2]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[5] = 0.1091214168497758*fUpwindQuad[26]-0.1091214168497758*fUpwindQuad[24]+0.1745942669596414*fUpwindQuad[23]-0.1745942669596414*fUpwindQuad[21]+0.1091214168497758*fUpwindQuad[20]-0.1091214168497758*fUpwindQuad[18]-0.1091214168497758*fUpwindQuad[8]+0.1091214168497758*fUpwindQuad[6]-0.1745942669596414*fUpwindQuad[5]+0.1745942669596414*fUpwindQuad[3]-0.1091214168497758*fUpwindQuad[2]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[6] = 0.1091214168497758*fUpwindQuad[26]+0.1745942669596414*fUpwindQuad[25]+0.1091214168497758*fUpwindQuad[24]-0.1091214168497758*fUpwindQuad[20]-0.1745942669596414*fUpwindQuad[19]-0.1091214168497758*fUpwindQuad[18]-0.1091214168497758*fUpwindQuad[8]-0.1745942669596414*fUpwindQuad[7]-0.1091214168497758*fUpwindQuad[6]+0.1091214168497758*fUpwindQuad[2]+0.1745942669596414*fUpwindQuad[1]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[7] = 0.05422286797270884*fUpwindQuad[26]-0.1084457359454177*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]+0.08675658875633419*fUpwindQuad[23]-0.1735131775126684*fUpwindQuad[22]+0.08675658875633419*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]-0.1084457359454177*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]+0.08675658875633419*fUpwindQuad[17]-0.1735131775126684*fUpwindQuad[16]+0.08675658875633419*fUpwindQuad[15]+0.1388105420101347*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]+0.1388105420101347*fUpwindQuad[12]+0.08675658875633419*fUpwindQuad[11]-0.1735131775126684*fUpwindQuad[10]+0.08675658875633419*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]-0.1084457359454177*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]+0.08675658875633419*fUpwindQuad[5]-0.1735131775126684*fUpwindQuad[4]+0.08675658875633419*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]-0.1084457359454177*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[8] = 0.05422286797270884*fUpwindQuad[26]+0.08675658875633419*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]-0.1084457359454177*fUpwindQuad[23]-0.1735131775126684*fUpwindQuad[22]-0.1084457359454177*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]+0.08675658875633419*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]+0.08675658875633419*fUpwindQuad[17]+0.1388105420101347*fUpwindQuad[16]+0.08675658875633419*fUpwindQuad[15]-0.1735131775126684*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]-0.1735131775126684*fUpwindQuad[12]+0.08675658875633419*fUpwindQuad[11]+0.1388105420101347*fUpwindQuad[10]+0.08675658875633419*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]+0.08675658875633419*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]-0.1084457359454177*fUpwindQuad[5]-0.1735131775126684*fUpwindQuad[4]-0.1084457359454177*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]+0.08675658875633419*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[9] = 0.05422286797270884*fUpwindQuad[26]+0.08675658875633419*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]+0.08675658875633419*fUpwindQuad[23]+0.1388105420101347*fUpwindQuad[22]+0.08675658875633419*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]+0.08675658875633419*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]-0.1084457359454177*fUpwindQuad[17]-0.1735131775126684*fUpwindQuad[16]-0.1084457359454177*fUpwindQuad[15]-0.1735131775126684*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]-0.1735131775126684*fUpwindQuad[12]-0.1084457359454177*fUpwindQuad[11]-0.1735131775126684*fUpwindQuad[10]-0.1084457359454177*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]+0.08675658875633419*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]+0.08675658875633419*fUpwindQuad[5]+0.1388105420101347*fUpwindQuad[4]+0.08675658875633419*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]+0.08675658875633419*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[10] = 0.1464017435263139*fUpwindQuad[26]-0.1464017435263139*fUpwindQuad[24]-0.1464017435263139*fUpwindQuad[20]+0.1464017435263139*fUpwindQuad[18]-0.1464017435263139*fUpwindQuad[8]+0.1464017435263139*fUpwindQuad[6]+0.1464017435263139*fUpwindQuad[2]-0.1464017435263139*fUpwindQuad[0]; 
  fUpwind[11] = 0.07274761123318395*fUpwindQuad[26]-0.1454952224663679*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.07274761123318395*fUpwindQuad[20]+0.1454952224663679*fUpwindQuad[19]-0.07274761123318395*fUpwindQuad[18]+0.1163961779730944*fUpwindQuad[17]-0.2327923559461888*fUpwindQuad[16]+0.1163961779730944*fUpwindQuad[15]-0.1163961779730944*fUpwindQuad[11]+0.2327923559461888*fUpwindQuad[10]-0.1163961779730944*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.1454952224663679*fUpwindQuad[7]+0.07274761123318395*fUpwindQuad[6]-0.07274761123318395*fUpwindQuad[2]+0.1454952224663679*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[12] = 0.07274761123318395*fUpwindQuad[26]-0.07274761123318395*fUpwindQuad[24]-0.1454952224663679*fUpwindQuad[23]+0.1454952224663679*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.07274761123318395*fUpwindQuad[18]+0.1163961779730944*fUpwindQuad[17]-0.1163961779730944*fUpwindQuad[15]-0.2327923559461888*fUpwindQuad[14]+0.2327923559461888*fUpwindQuad[12]+0.1163961779730944*fUpwindQuad[11]-0.1163961779730944*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.07274761123318395*fUpwindQuad[6]-0.1454952224663679*fUpwindQuad[5]+0.1454952224663679*fUpwindQuad[3]+0.07274761123318395*fUpwindQuad[2]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[13] = 0.07274761123318395*fUpwindQuad[26]-0.1454952224663679*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]+0.1163961779730944*fUpwindQuad[23]-0.2327923559461888*fUpwindQuad[22]+0.1163961779730944*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.1454952224663679*fUpwindQuad[19]+0.07274761123318395*fUpwindQuad[18]-0.07274761123318395*fUpwindQuad[8]+0.1454952224663679*fUpwindQuad[7]-0.07274761123318395*fUpwindQuad[6]-0.1163961779730944*fUpwindQuad[5]+0.2327923559461888*fUpwindQuad[4]-0.1163961779730944*fUpwindQuad[3]-0.07274761123318395*fUpwindQuad[2]+0.1454952224663679*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[14] = 0.07274761123318395*fUpwindQuad[26]+0.1163961779730944*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.1454952224663679*fUpwindQuad[23]-0.2327923559461888*fUpwindQuad[22]-0.1454952224663679*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]+0.1163961779730944*fUpwindQuad[19]+0.07274761123318395*fUpwindQuad[18]-0.07274761123318395*fUpwindQuad[8]-0.1163961779730944*fUpwindQuad[7]-0.07274761123318395*fUpwindQuad[6]+0.1454952224663679*fUpwindQuad[5]+0.2327923559461888*fUpwindQuad[4]+0.1454952224663679*fUpwindQuad[3]-0.07274761123318395*fUpwindQuad[2]-0.1163961779730944*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[15] = 0.07274761123318395*fUpwindQuad[26]-0.07274761123318395*fUpwindQuad[24]+0.1163961779730944*fUpwindQuad[23]-0.1163961779730944*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.07274761123318395*fUpwindQuad[18]-0.1454952224663679*fUpwindQuad[17]+0.1454952224663679*fUpwindQuad[15]-0.2327923559461888*fUpwindQuad[14]+0.2327923559461888*fUpwindQuad[12]-0.1454952224663679*fUpwindQuad[11]+0.1454952224663679*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.07274761123318395*fUpwindQuad[6]+0.1163961779730944*fUpwindQuad[5]-0.1163961779730944*fUpwindQuad[3]+0.07274761123318395*fUpwindQuad[2]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[16] = 0.07274761123318395*fUpwindQuad[26]+0.1163961779730944*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.07274761123318395*fUpwindQuad[20]-0.1163961779730944*fUpwindQuad[19]-0.07274761123318395*fUpwindQuad[18]-0.1454952224663679*fUpwindQuad[17]-0.2327923559461888*fUpwindQuad[16]-0.1454952224663679*fUpwindQuad[15]+0.1454952224663679*fUpwindQuad[11]+0.2327923559461888*fUpwindQuad[10]+0.1454952224663679*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]+0.1163961779730944*fUpwindQuad[7]+0.07274761123318395*fUpwindQuad[6]-0.07274761123318395*fUpwindQuad[2]-0.1163961779730944*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[17] = 0.09760116235087592*fUpwindQuad[26]-0.1952023247017519*fUpwindQuad[25]+0.09760116235087592*fUpwindQuad[24]-0.09760116235087592*fUpwindQuad[20]+0.1952023247017519*fUpwindQuad[19]-0.09760116235087592*fUpwindQuad[18]-0.09760116235087592*fUpwindQuad[8]+0.1952023247017519*fUpwindQuad[7]-0.09760116235087592*fUpwindQuad[6]+0.09760116235087592*fUpwindQuad[2]-0.1952023247017519*fUpwindQuad[1]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[18] = 0.09760116235087592*fUpwindQuad[26]-0.09760116235087592*fUpwindQuad[24]-0.1952023247017519*fUpwindQuad[23]+0.1952023247017519*fUpwindQuad[21]+0.09760116235087592*fUpwindQuad[20]-0.09760116235087592*fUpwindQuad[18]-0.09760116235087592*fUpwindQuad[8]+0.09760116235087592*fUpwindQuad[6]+0.1952023247017519*fUpwindQuad[5]-0.1952023247017519*fUpwindQuad[3]-0.09760116235087592*fUpwindQuad[2]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[19] = 0.09760116235087592*fUpwindQuad[26]-0.09760116235087592*fUpwindQuad[24]-0.09760116235087592*fUpwindQuad[20]+0.09760116235087592*fUpwindQuad[18]-0.1952023247017519*fUpwindQuad[17]+0.1952023247017519*fUpwindQuad[15]+0.1952023247017519*fUpwindQuad[11]-0.1952023247017519*fUpwindQuad[9]+0.09760116235087592*fUpwindQuad[8]-0.09760116235087592*fUpwindQuad[6]-0.09760116235087592*fUpwindQuad[2]+0.09760116235087592*fUpwindQuad[0]; 

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3535533905932737*alphaDrSurf[8]*fUpwind[12]+0.3535533905932737*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[4]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3162277660168379*alphaDrSurf[4]*fUpwind[12]+0.3162277660168379*fUpwind[4]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*fUpwind[7]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.2828427124746191*alphaDrSurf[11]*fUpwind[12]+0.3162277660168379*alphaDrSurf[2]*fUpwind[12]+0.2828427124746191*fUpwind[11]*alphaDrSurf[12]+0.3162277660168379*fUpwind[2]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[1]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[5] = 0.3535533905932737*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[4]*fUpwind[17]+0.3535533905932737*alphaDrSurf[12]*fUpwind[14]+0.3162277660168379*alphaDrSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3162277660168379*alphaDrSurf[4]*fUpwind[18]+0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.3162277660168379*alphaDrSurf[2]*fUpwind[14]+0.3535533905932737*alphaDrSurf[11]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  drag_incr[7] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[12]+0.2258769757263128*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[11]+0.2258769757263128*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[8] = 0.2258769757263128*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[11]*fUpwind[11]+0.2258769757263128*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[2]*fUpwind[2]; 
  drag_incr[9] = 0.3535533905932737*alphaDrSurf[4]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  drag_incr[10] = 0.282842712474619*alphaDrSurf[11]*fUpwind[18]+0.3162277660168379*alphaDrSurf[2]*fUpwind[18]+0.282842712474619*alphaDrSurf[12]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.3162277660168379*alphaDrSurf[4]*fUpwind[14]+0.3162277660168379*alphaDrSurf[4]*fUpwind[13]+0.3162277660168379*fUpwind[6]*alphaDrSurf[12]+0.3162277660168379*fUpwind[5]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]*fUpwind[10]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 
  drag_incr[11] = 0.2828427124746191*alphaDrSurf[4]*fUpwind[12]+0.2828427124746191*fUpwind[4]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[8]*fUpwind[11]+0.2258769757263128*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3162277660168379*fUpwind[8]*alphaDrSurf[11]+0.2258769757263128*fUpwind[7]*alphaDrSurf[11]+0.3535533905932737*fUpwind[0]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[4]+0.3162277660168379*fUpwind[1]*alphaDrSurf[4]; 
  drag_incr[12] = 0.2258769757263128*alphaDrSurf[8]*fUpwind[12]+0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.2258769757263128*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*fUpwind[7]*alphaDrSurf[12]+0.3535533905932737*fUpwind[0]*alphaDrSurf[12]+0.2828427124746191*alphaDrSurf[4]*fUpwind[11]+0.2828427124746191*fUpwind[4]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[2]*fUpwind[4]+0.3162277660168379*fUpwind[2]*alphaDrSurf[4]; 
  drag_incr[13] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[18]+0.2258769757263128*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[2]*fUpwind[17]+0.2258769757263128*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[6]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[5]; 
  drag_incr[14] = 0.2258769757263128*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3162277660168379*alphaDrSurf[11]*fUpwind[17]+0.2258769757263128*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]+0.3535533905932737*fUpwind[5]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[2]*fUpwind[6]; 
  drag_incr[15] = 0.3162277660168379*alphaDrSurf[11]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[19]+0.3535533905932737*alphaDrSurf[4]*fUpwind[16]+0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  drag_incr[16] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3162277660168379*alphaDrSurf[8]*fUpwind[16]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]+0.3535533905932737*alphaDrSurf[4]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[9]; 
  drag_incr[17] = 0.2828427124746191*alphaDrSurf[4]*fUpwind[18]+0.3162277660168379*alphaDrSurf[8]*fUpwind[17]+0.2258769757263128*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[11]*fUpwind[14]+0.2258769757263128*alphaDrSurf[11]*fUpwind[13]+0.3535533905932737*alphaDrSurf[2]*fUpwind[13]+0.282842712474619*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*fUpwind[3]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[5]; 
  drag_incr[18] = 0.2258769757263128*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.2828427124746191*alphaDrSurf[4]*fUpwind[17]+0.2258769757263128*alphaDrSurf[12]*fUpwind[14]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]+0.3162277660168379*alphaDrSurf[12]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[12]+0.282842712474619*fUpwind[10]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[6]; 
  drag_incr[19] = 0.3162277660168379*alphaDrSurf[8]*fUpwind[19]+0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3162277660168379*alphaDrSurf[12]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]+0.3162277660168379*alphaDrSurf[11]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[15]+0.3535533905932737*alphaDrSurf[4]*fUpwind[9]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += 0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += 0.7071067811865475*drag_incr[9]*rdv2; 
  out[14] += 1.58113883008419*drag_incr[0]*rdv2; 
  out[15] += 0.7071067811865475*drag_incr[10]*rdv2; 
  out[16] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += 0.7071067811865475*drag_incr[11]*rdv2; 
  out[20] += 0.7071067811865475*drag_incr[12]*rdv2; 
  out[21] += 0.7071067811865475*drag_incr[13]*rdv2; 
  out[22] += 0.7071067811865475*drag_incr[14]*rdv2; 
  out[23] += 0.7071067811865475*drag_incr[15]*rdv2; 
  out[24] += 0.7071067811865475*drag_incr[16]*rdv2; 
  out[25] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[26] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[27] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[28] += 1.58113883008419*drag_incr[1]*rdv2; 
  out[29] += 1.58113883008419*drag_incr[2]*rdv2; 
  out[30] += 1.58113883008419*drag_incr[3]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[32] += 0.7071067811865475*drag_incr[17]*rdv2; 
  out[33] += 0.7071067811865475*drag_incr[18]*rdv2; 
  out[34] += 0.7071067811865475*drag_incr[19]*rdv2; 
  out[35] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[36] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[37] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[38] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[39] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[40] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[41] += 1.58113883008419*drag_incr[4]*rdv2; 
  out[42] += 1.58113883008419*drag_incr[5]*rdv2; 
  out[43] += 1.58113883008419*drag_incr[6]*rdv2; 
  out[44] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[45] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[46] += 1.224744871391589*drag_incr[19]*rdv2; 
  out[47] += 1.58113883008419*drag_incr[10]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]-1.0*dxv[3]*nuSum[3]); 
  alphaDrSurf[7] = -0.7071067811865475*(2.0*sumNuUy[4]+(dxv[3]-2.0*w[3])*nuSum[4]); 
  alphaDrSurf[8] = -0.7071067811865475*(2.0*sumNuUy[5]+(dxv[3]-2.0*w[3])*nuSum[5]); 
  alphaDrSurf[11] = -0.7071067811865475*(2.0*sumNuUy[6]+(dxv[3]-2.0*w[3])*nuSum[6]); 
  alphaDrSurf[12] = -0.7071067811865475*(2.0*sumNuUy[7]+(dxv[3]-2.0*w[3])*nuSum[7]); 

  if ((-0.4242640687119281*alphaDrSurf[12])-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_2x2v_p2_surfvy_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x2v_p2_surfvy_quad_0(-1, fSkin); 
  } 
  if (0.5303300858899104*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x2v_p2_surfvy_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x2v_p2_surfvy_quad_1(-1, fSkin); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_2x2v_p2_surfvy_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_2x2v_p2_surfvy_quad_2(-1, fSkin); 
  } 
  if (0.5303300858899104*alphaDrSurf[12]-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_2x2v_p2_surfvy_quad_3(1, fEdge); 
  } else { 
    fUpwindQuad[3] = ser_2x2v_p2_surfvy_quad_3(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*(alphaDrSurf[8]+alphaDrSurf[7]) < 0) { 
    fUpwindQuad[4] = ser_2x2v_p2_surfvy_quad_4(1, fEdge); 
  } else { 
    fUpwindQuad[4] = ser_2x2v_p2_surfvy_quad_4(-1, fSkin); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[12])-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_2x2v_p2_surfvy_quad_5(1, fEdge); 
  } else { 
    fUpwindQuad[5] = ser_2x2v_p2_surfvy_quad_5(-1, fSkin); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_2x2v_p2_surfvy_quad_6(1, fEdge); 
  } else { 
    fUpwindQuad[6] = ser_2x2v_p2_surfvy_quad_6(-1, fSkin); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[11])+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_2x2v_p2_surfvy_quad_7(1, fEdge); 
  } else { 
    fUpwindQuad[7] = ser_2x2v_p2_surfvy_quad_7(-1, fSkin); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_2x2v_p2_surfvy_quad_8(1, fEdge); 
  } else { 
    fUpwindQuad[8] = ser_2x2v_p2_surfvy_quad_8(-1, fSkin); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_2x2v_p2_surfvy_quad_9(1, fEdge); 
  } else { 
    fUpwindQuad[9] = ser_2x2v_p2_surfvy_quad_9(-1, fSkin); 
  } 
  if (0.5303300858899104*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[10] = ser_2x2v_p2_surfvy_quad_10(1, fEdge); 
  } else { 
    fUpwindQuad[10] = ser_2x2v_p2_surfvy_quad_10(-1, fSkin); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_2x2v_p2_surfvy_quad_11(1, fEdge); 
  } else { 
    fUpwindQuad[11] = ser_2x2v_p2_surfvy_quad_11(-1, fSkin); 
  } 
  if (0.5303300858899104*alphaDrSurf[12]-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_2x2v_p2_surfvy_quad_12(1, fEdge); 
  } else { 
    fUpwindQuad[12] = ser_2x2v_p2_surfvy_quad_12(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*(alphaDrSurf[8]+alphaDrSurf[7]) < 0) { 
    fUpwindQuad[13] = ser_2x2v_p2_surfvy_quad_13(1, fEdge); 
  } else { 
    fUpwindQuad[13] = ser_2x2v_p2_surfvy_quad_13(-1, fSkin); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[12])-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_2x2v_p2_surfvy_quad_14(1, fEdge); 
  } else { 
    fUpwindQuad[14] = ser_2x2v_p2_surfvy_quad_14(-1, fSkin); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_2x2v_p2_surfvy_quad_15(1, fEdge); 
  } else { 
    fUpwindQuad[15] = ser_2x2v_p2_surfvy_quad_15(-1, fSkin); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[11])+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[16] = ser_2x2v_p2_surfvy_quad_16(1, fEdge); 
  } else { 
    fUpwindQuad[16] = ser_2x2v_p2_surfvy_quad_16(-1, fSkin); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[17] = ser_2x2v_p2_surfvy_quad_17(1, fEdge); 
  } else { 
    fUpwindQuad[17] = ser_2x2v_p2_surfvy_quad_17(-1, fSkin); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_2x2v_p2_surfvy_quad_18(1, fEdge); 
  } else { 
    fUpwindQuad[18] = ser_2x2v_p2_surfvy_quad_18(-1, fSkin); 
  } 
  if (0.5303300858899104*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[19] = ser_2x2v_p2_surfvy_quad_19(1, fEdge); 
  } else { 
    fUpwindQuad[19] = ser_2x2v_p2_surfvy_quad_19(-1, fSkin); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]-0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[20] = ser_2x2v_p2_surfvy_quad_20(1, fEdge); 
  } else { 
    fUpwindQuad[20] = ser_2x2v_p2_surfvy_quad_20(-1, fSkin); 
  } 
  if (0.5303300858899104*alphaDrSurf[12]-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_2x2v_p2_surfvy_quad_21(1, fEdge); 
  } else { 
    fUpwindQuad[21] = ser_2x2v_p2_surfvy_quad_21(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*(alphaDrSurf[8]+alphaDrSurf[7]) < 0) { 
    fUpwindQuad[22] = ser_2x2v_p2_surfvy_quad_22(1, fEdge); 
  } else { 
    fUpwindQuad[22] = ser_2x2v_p2_surfvy_quad_22(-1, fSkin); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[12])-0.3952847075210473*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[23] = ser_2x2v_p2_surfvy_quad_23(1, fEdge); 
  } else { 
    fUpwindQuad[23] = ser_2x2v_p2_surfvy_quad_23(-1, fSkin); 
  } 
  if ((-0.4242640687119281*alphaDrSurf[12])+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])-0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_2x2v_p2_surfvy_quad_24(1, fEdge); 
  } else { 
    fUpwindQuad[24] = ser_2x2v_p2_surfvy_quad_24(-1, fSkin); 
  } 
  if ((-0.5303300858899104*alphaDrSurf[11])+0.3162277660168379*alphaDrSurf[8]-0.3952847075210473*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[25] = ser_2x2v_p2_surfvy_quad_25(1, fEdge); 
  } else { 
    fUpwindQuad[25] = ser_2x2v_p2_surfvy_quad_25(-1, fSkin); 
  } 
  if (0.4242640687119281*alphaDrSurf[12]+0.4242640687119285*alphaDrSurf[11]+0.3162277660168379*(alphaDrSurf[8]+alphaDrSurf[7])+0.6363961030678926*alphaDrSurf[4]+0.4743416490252568*(alphaDrSurf[2]+alphaDrSurf[1])+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[26] = ser_2x2v_p2_surfvy_quad_26(1, fEdge); 
  } else { 
    fUpwindQuad[26] = ser_2x2v_p2_surfvy_quad_26(-1, fSkin); 
  } 

  fUpwind[0] = 0.06062300936098657*fUpwindQuad[26]+0.09699681497757856*fUpwindQuad[25]+0.06062300936098657*fUpwindQuad[24]+0.09699681497757856*fUpwindQuad[23]+0.1551949039641257*fUpwindQuad[22]+0.09699681497757856*fUpwindQuad[21]+0.06062300936098657*fUpwindQuad[20]+0.09699681497757856*fUpwindQuad[19]+0.06062300936098657*fUpwindQuad[18]+0.09699681497757856*fUpwindQuad[17]+0.1551949039641257*fUpwindQuad[16]+0.09699681497757856*fUpwindQuad[15]+0.1551949039641257*fUpwindQuad[14]+0.2483118463426013*fUpwindQuad[13]+0.1551949039641257*fUpwindQuad[12]+0.09699681497757856*fUpwindQuad[11]+0.1551949039641257*fUpwindQuad[10]+0.09699681497757856*fUpwindQuad[9]+0.06062300936098657*fUpwindQuad[8]+0.09699681497757856*fUpwindQuad[7]+0.06062300936098657*fUpwindQuad[6]+0.09699681497757856*fUpwindQuad[5]+0.1551949039641257*fUpwindQuad[4]+0.09699681497757856*fUpwindQuad[3]+0.06062300936098657*fUpwindQuad[2]+0.09699681497757856*fUpwindQuad[1]+0.06062300936098657*fUpwindQuad[0]; 
  fUpwind[1] = 0.08133430195906327*fUpwindQuad[26]-0.08133430195906327*fUpwindQuad[24]+0.1301348831345013*fUpwindQuad[23]-0.1301348831345013*fUpwindQuad[21]+0.08133430195906327*fUpwindQuad[20]-0.08133430195906327*fUpwindQuad[18]+0.1301348831345013*fUpwindQuad[17]-0.1301348831345013*fUpwindQuad[15]+0.2082158130152021*fUpwindQuad[14]-0.2082158130152021*fUpwindQuad[12]+0.1301348831345013*fUpwindQuad[11]-0.1301348831345013*fUpwindQuad[9]+0.08133430195906327*fUpwindQuad[8]-0.08133430195906327*fUpwindQuad[6]+0.1301348831345013*fUpwindQuad[5]-0.1301348831345013*fUpwindQuad[3]+0.08133430195906327*fUpwindQuad[2]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[2] = 0.08133430195906327*fUpwindQuad[26]+0.1301348831345013*fUpwindQuad[25]+0.08133430195906327*fUpwindQuad[24]-0.08133430195906327*fUpwindQuad[20]-0.1301348831345013*fUpwindQuad[19]-0.08133430195906327*fUpwindQuad[18]+0.1301348831345013*fUpwindQuad[17]+0.2082158130152021*fUpwindQuad[16]+0.1301348831345013*fUpwindQuad[15]-0.1301348831345013*fUpwindQuad[11]-0.2082158130152021*fUpwindQuad[10]-0.1301348831345013*fUpwindQuad[9]+0.08133430195906327*fUpwindQuad[8]+0.1301348831345013*fUpwindQuad[7]+0.08133430195906327*fUpwindQuad[6]-0.08133430195906327*fUpwindQuad[2]-0.1301348831345013*fUpwindQuad[1]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[3] = 0.08133430195906327*fUpwindQuad[26]+0.1301348831345013*fUpwindQuad[25]+0.08133430195906327*fUpwindQuad[24]+0.1301348831345013*fUpwindQuad[23]+0.2082158130152021*fUpwindQuad[22]+0.1301348831345013*fUpwindQuad[21]+0.08133430195906327*fUpwindQuad[20]+0.1301348831345013*fUpwindQuad[19]+0.08133430195906327*fUpwindQuad[18]-0.08133430195906327*fUpwindQuad[8]-0.1301348831345013*fUpwindQuad[7]-0.08133430195906327*fUpwindQuad[6]-0.1301348831345013*fUpwindQuad[5]-0.2082158130152021*fUpwindQuad[4]-0.1301348831345013*fUpwindQuad[3]-0.08133430195906327*fUpwindQuad[2]-0.1301348831345013*fUpwindQuad[1]-0.08133430195906327*fUpwindQuad[0]; 
  fUpwind[4] = 0.1091214168497758*fUpwindQuad[26]-0.1091214168497758*fUpwindQuad[24]-0.1091214168497758*fUpwindQuad[20]+0.1091214168497758*fUpwindQuad[18]+0.1745942669596414*fUpwindQuad[17]-0.1745942669596414*fUpwindQuad[15]-0.1745942669596414*fUpwindQuad[11]+0.1745942669596414*fUpwindQuad[9]+0.1091214168497758*fUpwindQuad[8]-0.1091214168497758*fUpwindQuad[6]-0.1091214168497758*fUpwindQuad[2]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[5] = 0.1091214168497758*fUpwindQuad[26]-0.1091214168497758*fUpwindQuad[24]+0.1745942669596414*fUpwindQuad[23]-0.1745942669596414*fUpwindQuad[21]+0.1091214168497758*fUpwindQuad[20]-0.1091214168497758*fUpwindQuad[18]-0.1091214168497758*fUpwindQuad[8]+0.1091214168497758*fUpwindQuad[6]-0.1745942669596414*fUpwindQuad[5]+0.1745942669596414*fUpwindQuad[3]-0.1091214168497758*fUpwindQuad[2]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[6] = 0.1091214168497758*fUpwindQuad[26]+0.1745942669596414*fUpwindQuad[25]+0.1091214168497758*fUpwindQuad[24]-0.1091214168497758*fUpwindQuad[20]-0.1745942669596414*fUpwindQuad[19]-0.1091214168497758*fUpwindQuad[18]-0.1091214168497758*fUpwindQuad[8]-0.1745942669596414*fUpwindQuad[7]-0.1091214168497758*fUpwindQuad[6]+0.1091214168497758*fUpwindQuad[2]+0.1745942669596414*fUpwindQuad[1]+0.1091214168497758*fUpwindQuad[0]; 
  fUpwind[7] = 0.05422286797270884*fUpwindQuad[26]-0.1084457359454177*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]+0.08675658875633419*fUpwindQuad[23]-0.1735131775126684*fUpwindQuad[22]+0.08675658875633419*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]-0.1084457359454177*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]+0.08675658875633419*fUpwindQuad[17]-0.1735131775126684*fUpwindQuad[16]+0.08675658875633419*fUpwindQuad[15]+0.1388105420101347*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]+0.1388105420101347*fUpwindQuad[12]+0.08675658875633419*fUpwindQuad[11]-0.1735131775126684*fUpwindQuad[10]+0.08675658875633419*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]-0.1084457359454177*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]+0.08675658875633419*fUpwindQuad[5]-0.1735131775126684*fUpwindQuad[4]+0.08675658875633419*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]-0.1084457359454177*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[8] = 0.05422286797270884*fUpwindQuad[26]+0.08675658875633419*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]-0.1084457359454177*fUpwindQuad[23]-0.1735131775126684*fUpwindQuad[22]-0.1084457359454177*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]+0.08675658875633419*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]+0.08675658875633419*fUpwindQuad[17]+0.1388105420101347*fUpwindQuad[16]+0.08675658875633419*fUpwindQuad[15]-0.1735131775126684*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]-0.1735131775126684*fUpwindQuad[12]+0.08675658875633419*fUpwindQuad[11]+0.1388105420101347*fUpwindQuad[10]+0.08675658875633419*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]+0.08675658875633419*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]-0.1084457359454177*fUpwindQuad[5]-0.1735131775126684*fUpwindQuad[4]-0.1084457359454177*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]+0.08675658875633419*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[9] = 0.05422286797270884*fUpwindQuad[26]+0.08675658875633419*fUpwindQuad[25]+0.05422286797270884*fUpwindQuad[24]+0.08675658875633419*fUpwindQuad[23]+0.1388105420101347*fUpwindQuad[22]+0.08675658875633419*fUpwindQuad[21]+0.05422286797270884*fUpwindQuad[20]+0.08675658875633419*fUpwindQuad[19]+0.05422286797270884*fUpwindQuad[18]-0.1084457359454177*fUpwindQuad[17]-0.1735131775126684*fUpwindQuad[16]-0.1084457359454177*fUpwindQuad[15]-0.1735131775126684*fUpwindQuad[14]-0.2776210840202695*fUpwindQuad[13]-0.1735131775126684*fUpwindQuad[12]-0.1084457359454177*fUpwindQuad[11]-0.1735131775126684*fUpwindQuad[10]-0.1084457359454177*fUpwindQuad[9]+0.05422286797270884*fUpwindQuad[8]+0.08675658875633419*fUpwindQuad[7]+0.05422286797270884*fUpwindQuad[6]+0.08675658875633419*fUpwindQuad[5]+0.1388105420101347*fUpwindQuad[4]+0.08675658875633419*fUpwindQuad[3]+0.05422286797270884*fUpwindQuad[2]+0.08675658875633419*fUpwindQuad[1]+0.05422286797270884*fUpwindQuad[0]; 
  fUpwind[10] = 0.1464017435263139*fUpwindQuad[26]-0.1464017435263139*fUpwindQuad[24]-0.1464017435263139*fUpwindQuad[20]+0.1464017435263139*fUpwindQuad[18]-0.1464017435263139*fUpwindQuad[8]+0.1464017435263139*fUpwindQuad[6]+0.1464017435263139*fUpwindQuad[2]-0.1464017435263139*fUpwindQuad[0]; 
  fUpwind[11] = 0.07274761123318395*fUpwindQuad[26]-0.1454952224663679*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.07274761123318395*fUpwindQuad[20]+0.1454952224663679*fUpwindQuad[19]-0.07274761123318395*fUpwindQuad[18]+0.1163961779730944*fUpwindQuad[17]-0.2327923559461888*fUpwindQuad[16]+0.1163961779730944*fUpwindQuad[15]-0.1163961779730944*fUpwindQuad[11]+0.2327923559461888*fUpwindQuad[10]-0.1163961779730944*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.1454952224663679*fUpwindQuad[7]+0.07274761123318395*fUpwindQuad[6]-0.07274761123318395*fUpwindQuad[2]+0.1454952224663679*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[12] = 0.07274761123318395*fUpwindQuad[26]-0.07274761123318395*fUpwindQuad[24]-0.1454952224663679*fUpwindQuad[23]+0.1454952224663679*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.07274761123318395*fUpwindQuad[18]+0.1163961779730944*fUpwindQuad[17]-0.1163961779730944*fUpwindQuad[15]-0.2327923559461888*fUpwindQuad[14]+0.2327923559461888*fUpwindQuad[12]+0.1163961779730944*fUpwindQuad[11]-0.1163961779730944*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.07274761123318395*fUpwindQuad[6]-0.1454952224663679*fUpwindQuad[5]+0.1454952224663679*fUpwindQuad[3]+0.07274761123318395*fUpwindQuad[2]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[13] = 0.07274761123318395*fUpwindQuad[26]-0.1454952224663679*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]+0.1163961779730944*fUpwindQuad[23]-0.2327923559461888*fUpwindQuad[22]+0.1163961779730944*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.1454952224663679*fUpwindQuad[19]+0.07274761123318395*fUpwindQuad[18]-0.07274761123318395*fUpwindQuad[8]+0.1454952224663679*fUpwindQuad[7]-0.07274761123318395*fUpwindQuad[6]-0.1163961779730944*fUpwindQuad[5]+0.2327923559461888*fUpwindQuad[4]-0.1163961779730944*fUpwindQuad[3]-0.07274761123318395*fUpwindQuad[2]+0.1454952224663679*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[14] = 0.07274761123318395*fUpwindQuad[26]+0.1163961779730944*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.1454952224663679*fUpwindQuad[23]-0.2327923559461888*fUpwindQuad[22]-0.1454952224663679*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]+0.1163961779730944*fUpwindQuad[19]+0.07274761123318395*fUpwindQuad[18]-0.07274761123318395*fUpwindQuad[8]-0.1163961779730944*fUpwindQuad[7]-0.07274761123318395*fUpwindQuad[6]+0.1454952224663679*fUpwindQuad[5]+0.2327923559461888*fUpwindQuad[4]+0.1454952224663679*fUpwindQuad[3]-0.07274761123318395*fUpwindQuad[2]-0.1163961779730944*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[15] = 0.07274761123318395*fUpwindQuad[26]-0.07274761123318395*fUpwindQuad[24]+0.1163961779730944*fUpwindQuad[23]-0.1163961779730944*fUpwindQuad[21]+0.07274761123318395*fUpwindQuad[20]-0.07274761123318395*fUpwindQuad[18]-0.1454952224663679*fUpwindQuad[17]+0.1454952224663679*fUpwindQuad[15]-0.2327923559461888*fUpwindQuad[14]+0.2327923559461888*fUpwindQuad[12]-0.1454952224663679*fUpwindQuad[11]+0.1454952224663679*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]-0.07274761123318395*fUpwindQuad[6]+0.1163961779730944*fUpwindQuad[5]-0.1163961779730944*fUpwindQuad[3]+0.07274761123318395*fUpwindQuad[2]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[16] = 0.07274761123318395*fUpwindQuad[26]+0.1163961779730944*fUpwindQuad[25]+0.07274761123318395*fUpwindQuad[24]-0.07274761123318395*fUpwindQuad[20]-0.1163961779730944*fUpwindQuad[19]-0.07274761123318395*fUpwindQuad[18]-0.1454952224663679*fUpwindQuad[17]-0.2327923559461888*fUpwindQuad[16]-0.1454952224663679*fUpwindQuad[15]+0.1454952224663679*fUpwindQuad[11]+0.2327923559461888*fUpwindQuad[10]+0.1454952224663679*fUpwindQuad[9]+0.07274761123318395*fUpwindQuad[8]+0.1163961779730944*fUpwindQuad[7]+0.07274761123318395*fUpwindQuad[6]-0.07274761123318395*fUpwindQuad[2]-0.1163961779730944*fUpwindQuad[1]-0.07274761123318395*fUpwindQuad[0]; 
  fUpwind[17] = 0.09760116235087592*fUpwindQuad[26]-0.1952023247017519*fUpwindQuad[25]+0.09760116235087592*fUpwindQuad[24]-0.09760116235087592*fUpwindQuad[20]+0.1952023247017519*fUpwindQuad[19]-0.09760116235087592*fUpwindQuad[18]-0.09760116235087592*fUpwindQuad[8]+0.1952023247017519*fUpwindQuad[7]-0.09760116235087592*fUpwindQuad[6]+0.09760116235087592*fUpwindQuad[2]-0.1952023247017519*fUpwindQuad[1]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[18] = 0.09760116235087592*fUpwindQuad[26]-0.09760116235087592*fUpwindQuad[24]-0.1952023247017519*fUpwindQuad[23]+0.1952023247017519*fUpwindQuad[21]+0.09760116235087592*fUpwindQuad[20]-0.09760116235087592*fUpwindQuad[18]-0.09760116235087592*fUpwindQuad[8]+0.09760116235087592*fUpwindQuad[6]+0.1952023247017519*fUpwindQuad[5]-0.1952023247017519*fUpwindQuad[3]-0.09760116235087592*fUpwindQuad[2]+0.09760116235087592*fUpwindQuad[0]; 
  fUpwind[19] = 0.09760116235087592*fUpwindQuad[26]-0.09760116235087592*fUpwindQuad[24]-0.09760116235087592*fUpwindQuad[20]+0.09760116235087592*fUpwindQuad[18]-0.1952023247017519*fUpwindQuad[17]+0.1952023247017519*fUpwindQuad[15]+0.1952023247017519*fUpwindQuad[11]-0.1952023247017519*fUpwindQuad[9]+0.09760116235087592*fUpwindQuad[8]-0.09760116235087592*fUpwindQuad[6]-0.09760116235087592*fUpwindQuad[2]+0.09760116235087592*fUpwindQuad[0]; 

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3535533905932737*alphaDrSurf[8]*fUpwind[12]+0.3535533905932737*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[4]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3162277660168379*alphaDrSurf[4]*fUpwind[12]+0.3162277660168379*fUpwind[4]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*fUpwind[7]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.2828427124746191*alphaDrSurf[11]*fUpwind[12]+0.3162277660168379*alphaDrSurf[2]*fUpwind[12]+0.2828427124746191*fUpwind[11]*alphaDrSurf[12]+0.3162277660168379*fUpwind[2]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[1]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  drag_incr[5] = 0.3535533905932737*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[4]*fUpwind[17]+0.3535533905932737*alphaDrSurf[12]*fUpwind[14]+0.3162277660168379*alphaDrSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3162277660168379*alphaDrSurf[4]*fUpwind[18]+0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.3162277660168379*alphaDrSurf[2]*fUpwind[14]+0.3535533905932737*alphaDrSurf[11]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  drag_incr[7] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[12]+0.2258769757263128*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[11]+0.2258769757263128*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[8] = 0.2258769757263128*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[11]*fUpwind[11]+0.2258769757263128*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[2]*fUpwind[2]; 
  drag_incr[9] = 0.3535533905932737*alphaDrSurf[4]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  drag_incr[10] = 0.282842712474619*alphaDrSurf[11]*fUpwind[18]+0.3162277660168379*alphaDrSurf[2]*fUpwind[18]+0.282842712474619*alphaDrSurf[12]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.3162277660168379*alphaDrSurf[4]*fUpwind[14]+0.3162277660168379*alphaDrSurf[4]*fUpwind[13]+0.3162277660168379*fUpwind[6]*alphaDrSurf[12]+0.3162277660168379*fUpwind[5]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]*fUpwind[10]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 
  drag_incr[11] = 0.2828427124746191*alphaDrSurf[4]*fUpwind[12]+0.2828427124746191*fUpwind[4]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[8]*fUpwind[11]+0.2258769757263128*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3162277660168379*fUpwind[8]*alphaDrSurf[11]+0.2258769757263128*fUpwind[7]*alphaDrSurf[11]+0.3535533905932737*fUpwind[0]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[4]+0.3162277660168379*fUpwind[1]*alphaDrSurf[4]; 
  drag_incr[12] = 0.2258769757263128*alphaDrSurf[8]*fUpwind[12]+0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.2258769757263128*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*fUpwind[7]*alphaDrSurf[12]+0.3535533905932737*fUpwind[0]*alphaDrSurf[12]+0.2828427124746191*alphaDrSurf[4]*fUpwind[11]+0.2828427124746191*fUpwind[4]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[2]*fUpwind[4]+0.3162277660168379*fUpwind[2]*alphaDrSurf[4]; 
  drag_incr[13] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[18]+0.2258769757263128*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[2]*fUpwind[17]+0.2258769757263128*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[6]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[5]; 
  drag_incr[14] = 0.2258769757263128*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3162277660168379*alphaDrSurf[11]*fUpwind[17]+0.2258769757263128*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]+0.3535533905932737*fUpwind[5]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[2]*fUpwind[6]; 
  drag_incr[15] = 0.3162277660168379*alphaDrSurf[11]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[19]+0.3535533905932737*alphaDrSurf[4]*fUpwind[16]+0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  drag_incr[16] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3162277660168379*alphaDrSurf[8]*fUpwind[16]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]+0.3535533905932737*alphaDrSurf[4]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[9]; 
  drag_incr[17] = 0.2828427124746191*alphaDrSurf[4]*fUpwind[18]+0.3162277660168379*alphaDrSurf[8]*fUpwind[17]+0.2258769757263128*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[11]*fUpwind[14]+0.2258769757263128*alphaDrSurf[11]*fUpwind[13]+0.3535533905932737*alphaDrSurf[2]*fUpwind[13]+0.282842712474619*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*fUpwind[3]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[5]; 
  drag_incr[18] = 0.2258769757263128*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.2828427124746191*alphaDrSurf[4]*fUpwind[17]+0.2258769757263128*alphaDrSurf[12]*fUpwind[14]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]+0.3162277660168379*alphaDrSurf[12]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[12]+0.282842712474619*fUpwind[10]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[6]; 
  drag_incr[19] = 0.3162277660168379*alphaDrSurf[8]*fUpwind[19]+0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3162277660168379*alphaDrSurf[12]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]+0.3162277660168379*alphaDrSurf[11]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[15]+0.3535533905932737*alphaDrSurf[4]*fUpwind[9]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[3] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[4] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[5] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[6] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[7] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[8] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[10] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[11] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += -0.7071067811865475*drag_incr[8]*rdv2; 
  out[13] += -0.7071067811865475*drag_incr[9]*rdv2; 
  out[14] += -1.58113883008419*drag_incr[0]*rdv2; 
  out[15] += -0.7071067811865475*drag_incr[10]*rdv2; 
  out[16] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[17] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += -0.7071067811865475*drag_incr[11]*rdv2; 
  out[20] += -0.7071067811865475*drag_incr[12]*rdv2; 
  out[21] += -0.7071067811865475*drag_incr[13]*rdv2; 
  out[22] += -0.7071067811865475*drag_incr[14]*rdv2; 
  out[23] += -0.7071067811865475*drag_incr[15]*rdv2; 
  out[24] += -0.7071067811865475*drag_incr[16]*rdv2; 
  out[25] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[26] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[27] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[28] += -1.58113883008419*drag_incr[1]*rdv2; 
  out[29] += -1.58113883008419*drag_incr[2]*rdv2; 
  out[30] += -1.58113883008419*drag_incr[3]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[32] += -0.7071067811865475*drag_incr[17]*rdv2; 
  out[33] += -0.7071067811865475*drag_incr[18]*rdv2; 
  out[34] += -0.7071067811865475*drag_incr[19]*rdv2; 
  out[35] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[36] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[37] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[38] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[39] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[40] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[41] += -1.58113883008419*drag_incr[4]*rdv2; 
  out[42] += -1.58113883008419*drag_incr[5]*rdv2; 
  out[43] += -1.58113883008419*drag_incr[6]*rdv2; 
  out[44] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[45] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[46] += 1.224744871391589*drag_incr[19]*rdv2; 
  out[47] += -1.58113883008419*drag_incr[10]*rdv2; 

  } 
} 
