#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x3v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_lbo_boundary_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[9]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf[20] = {0.0}; 
  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};;
  double drag_incr[20] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[1]+dxv[1])-2.0*sumNuUx[0]; 
  alphaDrSurf[1] = 2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]+dxv[1]*nuSum[1]; 
  alphaDrSurf[7] = (2.0*w[1]+dxv[1])*nuSum[2]-2.0*sumNuUx[2]; 

  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x3v_p2_surfvx_quad_0(1, fSkin); 
  } else { 
    fUpwindQuad[0] = ser_1x3v_p2_surfvx_quad_0(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[1] = ser_1x3v_p2_surfvx_quad_1(1, fSkin); 
  } else { 
    fUpwindQuad[1] = ser_1x3v_p2_surfvx_quad_1(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x3v_p2_surfvx_quad_2(1, fSkin); 
  } else { 
    fUpwindQuad[2] = ser_1x3v_p2_surfvx_quad_2(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_1x3v_p2_surfvx_quad_3(1, fSkin); 
  } else { 
    fUpwindQuad[3] = ser_1x3v_p2_surfvx_quad_3(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[4] = ser_1x3v_p2_surfvx_quad_4(1, fSkin); 
  } else { 
    fUpwindQuad[4] = ser_1x3v_p2_surfvx_quad_4(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_1x3v_p2_surfvx_quad_5(1, fSkin); 
  } else { 
    fUpwindQuad[5] = ser_1x3v_p2_surfvx_quad_5(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_1x3v_p2_surfvx_quad_6(1, fSkin); 
  } else { 
    fUpwindQuad[6] = ser_1x3v_p2_surfvx_quad_6(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[7] = ser_1x3v_p2_surfvx_quad_7(1, fSkin); 
  } else { 
    fUpwindQuad[7] = ser_1x3v_p2_surfvx_quad_7(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_1x3v_p2_surfvx_quad_8(1, fSkin); 
  } else { 
    fUpwindQuad[8] = ser_1x3v_p2_surfvx_quad_8(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_1x3v_p2_surfvx_quad_9(1, fSkin); 
  } else { 
    fUpwindQuad[9] = ser_1x3v_p2_surfvx_quad_9(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[10] = ser_1x3v_p2_surfvx_quad_10(1, fSkin); 
  } else { 
    fUpwindQuad[10] = ser_1x3v_p2_surfvx_quad_10(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_1x3v_p2_surfvx_quad_11(1, fSkin); 
  } else { 
    fUpwindQuad[11] = ser_1x3v_p2_surfvx_quad_11(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_1x3v_p2_surfvx_quad_12(1, fSkin); 
  } else { 
    fUpwindQuad[12] = ser_1x3v_p2_surfvx_quad_12(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[13] = ser_1x3v_p2_surfvx_quad_13(1, fSkin); 
  } else { 
    fUpwindQuad[13] = ser_1x3v_p2_surfvx_quad_13(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_1x3v_p2_surfvx_quad_14(1, fSkin); 
  } else { 
    fUpwindQuad[14] = ser_1x3v_p2_surfvx_quad_14(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_1x3v_p2_surfvx_quad_15(1, fSkin); 
  } else { 
    fUpwindQuad[15] = ser_1x3v_p2_surfvx_quad_15(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[16] = ser_1x3v_p2_surfvx_quad_16(1, fSkin); 
  } else { 
    fUpwindQuad[16] = ser_1x3v_p2_surfvx_quad_16(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[17] = ser_1x3v_p2_surfvx_quad_17(1, fSkin); 
  } else { 
    fUpwindQuad[17] = ser_1x3v_p2_surfvx_quad_17(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_1x3v_p2_surfvx_quad_18(1, fSkin); 
  } else { 
    fUpwindQuad[18] = ser_1x3v_p2_surfvx_quad_18(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[19] = ser_1x3v_p2_surfvx_quad_19(1, fSkin); 
  } else { 
    fUpwindQuad[19] = ser_1x3v_p2_surfvx_quad_19(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[20] = ser_1x3v_p2_surfvx_quad_20(1, fSkin); 
  } else { 
    fUpwindQuad[20] = ser_1x3v_p2_surfvx_quad_20(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_1x3v_p2_surfvx_quad_21(1, fSkin); 
  } else { 
    fUpwindQuad[21] = ser_1x3v_p2_surfvx_quad_21(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[22] = ser_1x3v_p2_surfvx_quad_22(1, fSkin); 
  } else { 
    fUpwindQuad[22] = ser_1x3v_p2_surfvx_quad_22(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[23] = ser_1x3v_p2_surfvx_quad_23(1, fSkin); 
  } else { 
    fUpwindQuad[23] = ser_1x3v_p2_surfvx_quad_23(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_1x3v_p2_surfvx_quad_24(1, fSkin); 
  } else { 
    fUpwindQuad[24] = ser_1x3v_p2_surfvx_quad_24(-1, fEdge); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[25] = ser_1x3v_p2_surfvx_quad_25(1, fSkin); 
  } else { 
    fUpwindQuad[25] = ser_1x3v_p2_surfvx_quad_25(-1, fEdge); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[26] = ser_1x3v_p2_surfvx_quad_26(1, fSkin); 
  } else { 
    fUpwindQuad[26] = ser_1x3v_p2_surfvx_quad_26(-1, fEdge); 
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

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]; 
  drag_incr[5] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]; 
  drag_incr[7] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[8] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]; 
  drag_incr[9] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  drag_incr[10] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]; 
  drag_incr[11] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[12] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]; 
  drag_incr[13] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[5]; 
  drag_incr[14] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]; 
  drag_incr[15] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  drag_incr[16] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]; 
  drag_incr[17] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]; 
  drag_incr[18] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]; 
  drag_incr[19] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]; 

  out[0] += 0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += 0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[3] += 0.7071067811865475*drag_incr[2]*rdv2; 
  out[4] += 0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[6] += 0.7071067811865475*drag_incr[4]*rdv2; 
  out[7] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[8] += 0.7071067811865475*drag_incr[5]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[10] += 0.7071067811865475*drag_incr[6]*rdv2; 
  out[11] += 0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += 1.58113883008419*drag_incr[0]*rdv2; 
  out[13] += 0.7071067811865475*drag_incr[8]*rdv2; 
  out[14] += 0.7071067811865475*drag_incr[9]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[16] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[17] += 0.7071067811865475*drag_incr[10]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[20] += 1.58113883008419*drag_incr[1]*rdv2; 
  out[21] += 0.7071067811865475*drag_incr[11]*rdv2; 
  out[22] += 1.58113883008419*drag_incr[2]*rdv2; 
  out[23] += 0.7071067811865475*drag_incr[12]*rdv2; 
  out[24] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[25] += 0.7071067811865475*drag_incr[13]*rdv2; 
  out[26] += 1.58113883008419*drag_incr[3]*rdv2; 
  out[27] += 0.7071067811865475*drag_incr[14]*rdv2; 
  out[28] += 0.7071067811865475*drag_incr[15]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[30] += 0.7071067811865475*drag_incr[16]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[32] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[33] += 1.58113883008419*drag_incr[4]*rdv2; 
  out[34] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[35] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[36] += 1.58113883008419*drag_incr[5]*rdv2; 
  out[37] += 0.7071067811865475*drag_incr[17]*rdv2; 
  out[38] += 1.58113883008419*drag_incr[6]*rdv2; 
  out[39] += 0.7071067811865475*drag_incr[18]*rdv2; 
  out[40] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[41] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[42] += 0.7071067811865475*drag_incr[19]*rdv2; 
  out[43] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[44] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[45] += 1.58113883008419*drag_incr[10]*rdv2; 
  out[46] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[47] += 1.224744871391589*drag_incr[19]*rdv2; 

  } else { 

  alphaDrSurf[0] = nuSum[0]*(2.0*w[1]-1.0*dxv[1])-2.0*sumNuUx[0]; 
  alphaDrSurf[1] = 2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]-1.0*dxv[1]*nuSum[1]; 
  alphaDrSurf[7] = (2.0*w[1]-1.0*dxv[1])*nuSum[2]-2.0*sumNuUx[2]; 

  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x3v_p2_surfvx_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_1x3v_p2_surfvx_quad_0(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[1] = ser_1x3v_p2_surfvx_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_1x3v_p2_surfvx_quad_1(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x3v_p2_surfvx_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_1x3v_p2_surfvx_quad_2(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_1x3v_p2_surfvx_quad_3(1, fEdge); 
  } else { 
    fUpwindQuad[3] = ser_1x3v_p2_surfvx_quad_3(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[4] = ser_1x3v_p2_surfvx_quad_4(1, fEdge); 
  } else { 
    fUpwindQuad[4] = ser_1x3v_p2_surfvx_quad_4(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_1x3v_p2_surfvx_quad_5(1, fEdge); 
  } else { 
    fUpwindQuad[5] = ser_1x3v_p2_surfvx_quad_5(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_1x3v_p2_surfvx_quad_6(1, fEdge); 
  } else { 
    fUpwindQuad[6] = ser_1x3v_p2_surfvx_quad_6(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[7] = ser_1x3v_p2_surfvx_quad_7(1, fEdge); 
  } else { 
    fUpwindQuad[7] = ser_1x3v_p2_surfvx_quad_7(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[8] = ser_1x3v_p2_surfvx_quad_8(1, fEdge); 
  } else { 
    fUpwindQuad[8] = ser_1x3v_p2_surfvx_quad_8(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_1x3v_p2_surfvx_quad_9(1, fEdge); 
  } else { 
    fUpwindQuad[9] = ser_1x3v_p2_surfvx_quad_9(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[10] = ser_1x3v_p2_surfvx_quad_10(1, fEdge); 
  } else { 
    fUpwindQuad[10] = ser_1x3v_p2_surfvx_quad_10(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[11] = ser_1x3v_p2_surfvx_quad_11(1, fEdge); 
  } else { 
    fUpwindQuad[11] = ser_1x3v_p2_surfvx_quad_11(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_1x3v_p2_surfvx_quad_12(1, fEdge); 
  } else { 
    fUpwindQuad[12] = ser_1x3v_p2_surfvx_quad_12(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[13] = ser_1x3v_p2_surfvx_quad_13(1, fEdge); 
  } else { 
    fUpwindQuad[13] = ser_1x3v_p2_surfvx_quad_13(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[14] = ser_1x3v_p2_surfvx_quad_14(1, fEdge); 
  } else { 
    fUpwindQuad[14] = ser_1x3v_p2_surfvx_quad_14(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_1x3v_p2_surfvx_quad_15(1, fEdge); 
  } else { 
    fUpwindQuad[15] = ser_1x3v_p2_surfvx_quad_15(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[16] = ser_1x3v_p2_surfvx_quad_16(1, fEdge); 
  } else { 
    fUpwindQuad[16] = ser_1x3v_p2_surfvx_quad_16(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[17] = ser_1x3v_p2_surfvx_quad_17(1, fEdge); 
  } else { 
    fUpwindQuad[17] = ser_1x3v_p2_surfvx_quad_17(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_1x3v_p2_surfvx_quad_18(1, fEdge); 
  } else { 
    fUpwindQuad[18] = ser_1x3v_p2_surfvx_quad_18(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[19] = ser_1x3v_p2_surfvx_quad_19(1, fEdge); 
  } else { 
    fUpwindQuad[19] = ser_1x3v_p2_surfvx_quad_19(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[20] = ser_1x3v_p2_surfvx_quad_20(1, fEdge); 
  } else { 
    fUpwindQuad[20] = ser_1x3v_p2_surfvx_quad_20(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_1x3v_p2_surfvx_quad_21(1, fEdge); 
  } else { 
    fUpwindQuad[21] = ser_1x3v_p2_surfvx_quad_21(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[22] = ser_1x3v_p2_surfvx_quad_22(1, fEdge); 
  } else { 
    fUpwindQuad[22] = ser_1x3v_p2_surfvx_quad_22(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[23] = ser_1x3v_p2_surfvx_quad_23(1, fEdge); 
  } else { 
    fUpwindQuad[23] = ser_1x3v_p2_surfvx_quad_23(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_1x3v_p2_surfvx_quad_24(1, fEdge); 
  } else { 
    fUpwindQuad[24] = ser_1x3v_p2_surfvx_quad_24(-1, fSkin); 
  } 
  if (0.3535533905932737*alphaDrSurf[0]-0.3952847075210473*alphaDrSurf[7] < 0) { 
    fUpwindQuad[25] = ser_1x3v_p2_surfvx_quad_25(1, fEdge); 
  } else { 
    fUpwindQuad[25] = ser_1x3v_p2_surfvx_quad_25(-1, fSkin); 
  } 
  if (0.3162277660168379*alphaDrSurf[7]+0.4743416490252568*alphaDrSurf[1]+0.3535533905932737*alphaDrSurf[0] < 0) { 
    fUpwindQuad[26] = ser_1x3v_p2_surfvx_quad_26(1, fEdge); 
  } else { 
    fUpwindQuad[26] = ser_1x3v_p2_surfvx_quad_26(-1, fSkin); 
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

  drag_incr[0] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  drag_incr[1] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  drag_incr[2] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]; 
  drag_incr[3] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  drag_incr[4] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]; 
  drag_incr[5] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  drag_incr[6] = 0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]; 
  drag_incr[7] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  drag_incr[8] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]; 
  drag_incr[9] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  drag_incr[10] = 0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]; 
  drag_incr[11] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[4]; 
  drag_incr[12] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]; 
  drag_incr[13] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[5]; 
  drag_incr[14] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]; 
  drag_incr[15] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  drag_incr[16] = 0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]; 
  drag_incr[17] = 0.2258769757263128*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]; 
  drag_incr[18] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]; 
  drag_incr[19] = 0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]; 

  out[0] += -0.7071067811865475*drag_incr[0]*rdv2; 
  out[1] += -0.7071067811865475*drag_incr[1]*rdv2; 
  out[2] += 1.224744871391589*drag_incr[0]*rdv2; 
  out[3] += -0.7071067811865475*drag_incr[2]*rdv2; 
  out[4] += -0.7071067811865475*drag_incr[3]*rdv2; 
  out[5] += 1.224744871391589*drag_incr[1]*rdv2; 
  out[6] += -0.7071067811865475*drag_incr[4]*rdv2; 
  out[7] += 1.224744871391589*drag_incr[2]*rdv2; 
  out[8] += -0.7071067811865475*drag_incr[5]*rdv2; 
  out[9] += 1.224744871391589*drag_incr[3]*rdv2; 
  out[10] += -0.7071067811865475*drag_incr[6]*rdv2; 
  out[11] += -0.7071067811865475*drag_incr[7]*rdv2; 
  out[12] += -1.58113883008419*drag_incr[0]*rdv2; 
  out[13] += -0.7071067811865475*drag_incr[8]*rdv2; 
  out[14] += -0.7071067811865475*drag_incr[9]*rdv2; 
  out[15] += 1.224744871391589*drag_incr[4]*rdv2; 
  out[16] += 1.224744871391589*drag_incr[5]*rdv2; 
  out[17] += -0.7071067811865475*drag_incr[10]*rdv2; 
  out[18] += 1.224744871391589*drag_incr[6]*rdv2; 
  out[19] += 1.224744871391589*drag_incr[7]*rdv2; 
  out[20] += -1.58113883008419*drag_incr[1]*rdv2; 
  out[21] += -0.7071067811865475*drag_incr[11]*rdv2; 
  out[22] += -1.58113883008419*drag_incr[2]*rdv2; 
  out[23] += -0.7071067811865475*drag_incr[12]*rdv2; 
  out[24] += 1.224744871391589*drag_incr[8]*rdv2; 
  out[25] += -0.7071067811865475*drag_incr[13]*rdv2; 
  out[26] += -1.58113883008419*drag_incr[3]*rdv2; 
  out[27] += -0.7071067811865475*drag_incr[14]*rdv2; 
  out[28] += -0.7071067811865475*drag_incr[15]*rdv2; 
  out[29] += 1.224744871391589*drag_incr[9]*rdv2; 
  out[30] += -0.7071067811865475*drag_incr[16]*rdv2; 
  out[31] += 1.224744871391589*drag_incr[10]*rdv2; 
  out[32] += 1.224744871391589*drag_incr[11]*rdv2; 
  out[33] += -1.58113883008419*drag_incr[4]*rdv2; 
  out[34] += 1.224744871391589*drag_incr[12]*rdv2; 
  out[35] += 1.224744871391589*drag_incr[13]*rdv2; 
  out[36] += -1.58113883008419*drag_incr[5]*rdv2; 
  out[37] += -0.7071067811865475*drag_incr[17]*rdv2; 
  out[38] += -1.58113883008419*drag_incr[6]*rdv2; 
  out[39] += -0.7071067811865475*drag_incr[18]*rdv2; 
  out[40] += 1.224744871391589*drag_incr[14]*rdv2; 
  out[41] += 1.224744871391589*drag_incr[15]*rdv2; 
  out[42] += -0.7071067811865475*drag_incr[19]*rdv2; 
  out[43] += 1.224744871391589*drag_incr[16]*rdv2; 
  out[44] += 1.224744871391589*drag_incr[17]*rdv2; 
  out[45] += -1.58113883008419*drag_incr[10]*rdv2; 
  out[46] += 1.224744871391589*drag_incr[18]*rdv2; 
  out[47] += 1.224744871391589*drag_incr[19]*rdv2; 

  } 
} 
