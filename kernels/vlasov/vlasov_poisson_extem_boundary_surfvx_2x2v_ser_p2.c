#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x2v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_poisson_extem_boundary_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out) 
{ 
  // w:           Cell-center coordinates.
  // dxv[NDIM]:   Cell spacing.
  // fac_phi:     potential (scaled by appropriate factors).
  // vecA:        vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // edge:        Determines if the update is for the left edge (-1) or right edge (+1).
  // fSkin/fEdge: Input Distribution function in skin cell/last edge cell 
  // out:         Output distribution function in skin cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[8]; 
  double alpha[20] = {0.0}; 

  alpha[0] = (-2.449489742783178*A0[2]*dx11*wv2)+2.449489742783178*A1[1]*dx10*wv2-2.449489742783178*phi[1]*dx10; 
  alpha[1] = (-2.449489742783178*A0[3]*dx11*wv2)+5.477225575051662*A1[4]*dx10*wv2-5.477225575051662*phi[4]*dx10; 
  alpha[2] = (-5.477225575051662*A0[5]*dx11*wv2)+2.449489742783178*A1[3]*dx10*wv2-2.449489742783178*phi[3]*dx10; 
  alpha[3] = 0.7071067811865475*A1[1]*dv2*dx10-0.7071067811865475*A0[2]*dv2*dx11; 
  alpha[4] = (-5.477225575051662*A0[7]*dx11*wv2)+5.477225575051662*A1[6]*dx10*wv2-5.477225575051662*phi[6]*dx10; 
  alpha[5] = 1.58113883008419*A1[4]*dv2*dx10-0.7071067811865475*A0[3]*dv2*dx11; 
  alpha[6] = 0.7071067811865475*A1[3]*dv2*dx10-1.58113883008419*A0[5]*dv2*dx11; 
  alpha[7] = -2.449489742783178*A0[6]*dx11*wv2; 
  alpha[8] = 2.449489742783178*A1[7]*dx10*wv2-2.449489742783178*phi[7]*dx10; 
  alpha[10] = 1.58113883008419*A1[6]*dv2*dx10-1.58113883008419*A0[7]*dv2*dx11; 
  alpha[13] = -0.7071067811865475*A0[6]*dv2*dx11; 
  alpha[14] = 0.7071067811865475*A1[7]*dv2*dx10; 

  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};
  double Ghat[20] = {0.0}; 

  if (edge == -1) { 

  if ((-0.4242640687119285*(alpha[14]+alpha[13]))-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[0] = ser_2x2v_p2_surfvx_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_2x2v_p2_surfvx_quad_0(-1, fEdge); 
  } 
  if ((-0.4242640687119285*alpha[14])+0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]-0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[1] = ser_2x2v_p2_surfvx_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_2x2v_p2_surfvx_quad_1(-1, fEdge); 
  } 
  if ((-0.4242640687119285*(alpha[14]+alpha[13]))+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[2] = ser_2x2v_p2_surfvx_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_2x2v_p2_surfvx_quad_2(-1, fEdge); 
  } 
  if (0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[3] = ser_2x2v_p2_surfvx_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_2x2v_p2_surfvx_quad_3(-1, fEdge); 
  } 
  if (0.5303300858899104*(alpha[14]+alpha[13])-0.3952847075210473*(alpha[8]+alpha[7])-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[4] = ser_2x2v_p2_surfvx_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_2x2v_p2_surfvx_quad_4(-1, fEdge); 
  } 
  if (0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[5] = ser_2x2v_p2_surfvx_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_2x2v_p2_surfvx_quad_5(-1, fEdge); 
  } 
  if ((-0.4242640687119285*(alpha[14]+alpha[13]))+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[6] = ser_2x2v_p2_surfvx_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_2x2v_p2_surfvx_quad_6(-1, fEdge); 
  } 
  if ((-0.4242640687119285*alpha[14])+0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[7] = ser_2x2v_p2_surfvx_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_2x2v_p2_surfvx_quad_7(-1, fEdge); 
  } 
  if ((-0.4242640687119285*(alpha[14]+alpha[13]))-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[8] = ser_2x2v_p2_surfvx_quad_8(1, fSkin); 
  } else { 

    fUpwindQuad[8] = ser_2x2v_p2_surfvx_quad_8(-1, fEdge); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[9] = ser_2x2v_p2_surfvx_quad_9(1, fSkin); 
  } else { 

    fUpwindQuad[9] = ser_2x2v_p2_surfvx_quad_9(-1, fEdge); 
  } 
  if (0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[10] = ser_2x2v_p2_surfvx_quad_10(1, fSkin); 
  } else { 

    fUpwindQuad[10] = ser_2x2v_p2_surfvx_quad_10(-1, fEdge); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[11] = ser_2x2v_p2_surfvx_quad_11(1, fSkin); 
  } else { 

    fUpwindQuad[11] = ser_2x2v_p2_surfvx_quad_11(-1, fEdge); 
  } 
  if ((-0.3952847075210473*alpha[8])+0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[12] = ser_2x2v_p2_surfvx_quad_12(1, fSkin); 
  } else { 

    fUpwindQuad[12] = ser_2x2v_p2_surfvx_quad_12(-1, fEdge); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*(alpha[8]+alpha[7]) > 0) { 

    fUpwindQuad[13] = ser_2x2v_p2_surfvx_quad_13(1, fSkin); 
  } else { 

    fUpwindQuad[13] = ser_2x2v_p2_surfvx_quad_13(-1, fEdge); 
  } 
  if ((-0.3952847075210473*alpha[8])+0.3162277660168379*alpha[7]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[14] = ser_2x2v_p2_surfvx_quad_14(1, fSkin); 
  } else { 

    fUpwindQuad[14] = ser_2x2v_p2_surfvx_quad_14(-1, fEdge); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[15] = ser_2x2v_p2_surfvx_quad_15(1, fSkin); 
  } else { 

    fUpwindQuad[15] = ser_2x2v_p2_surfvx_quad_15(-1, fEdge); 
  } 
  if (0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[16] = ser_2x2v_p2_surfvx_quad_16(1, fSkin); 
  } else { 

    fUpwindQuad[16] = ser_2x2v_p2_surfvx_quad_16(-1, fEdge); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[17] = ser_2x2v_p2_surfvx_quad_17(1, fSkin); 
  } else { 

    fUpwindQuad[17] = ser_2x2v_p2_surfvx_quad_17(-1, fEdge); 
  } 
  if (0.4242640687119285*(alpha[14]+alpha[13])+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[18] = ser_2x2v_p2_surfvx_quad_18(1, fSkin); 
  } else { 

    fUpwindQuad[18] = ser_2x2v_p2_surfvx_quad_18(-1, fEdge); 
  } 
  if (0.4242640687119285*alpha[14]-0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[19] = ser_2x2v_p2_surfvx_quad_19(1, fSkin); 
  } else { 

    fUpwindQuad[19] = ser_2x2v_p2_surfvx_quad_19(-1, fEdge); 
  } 
  if (0.4242640687119285*(alpha[14]+alpha[13])-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[20] = ser_2x2v_p2_surfvx_quad_20(1, fSkin); 
  } else { 

    fUpwindQuad[20] = ser_2x2v_p2_surfvx_quad_20(-1, fEdge); 
  } 
  if ((-0.5303300858899104*alpha[14])+0.4242640687119285*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[21] = ser_2x2v_p2_surfvx_quad_21(1, fSkin); 
  } else { 

    fUpwindQuad[21] = ser_2x2v_p2_surfvx_quad_21(-1, fEdge); 
  } 
  if ((-0.5303300858899104*(alpha[14]+alpha[13]))-0.3952847075210473*(alpha[8]+alpha[7])+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[22] = ser_2x2v_p2_surfvx_quad_22(1, fSkin); 
  } else { 

    fUpwindQuad[22] = ser_2x2v_p2_surfvx_quad_22(-1, fEdge); 
  } 
  if ((-0.5303300858899104*alpha[14])+0.4242640687119285*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]+0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[23] = ser_2x2v_p2_surfvx_quad_23(1, fSkin); 
  } else { 

    fUpwindQuad[23] = ser_2x2v_p2_surfvx_quad_23(-1, fEdge); 
  } 
  if (0.4242640687119285*(alpha[14]+alpha[13])-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[24] = ser_2x2v_p2_surfvx_quad_24(1, fSkin); 
  } else { 

    fUpwindQuad[24] = ser_2x2v_p2_surfvx_quad_24(-1, fEdge); 
  } 
  if (0.4242640687119285*alpha[14]-0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]+0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[25] = ser_2x2v_p2_surfvx_quad_25(1, fSkin); 
  } else { 

    fUpwindQuad[25] = ser_2x2v_p2_surfvx_quad_25(-1, fEdge); 
  } 
  if (0.4242640687119285*(alpha[14]+alpha[13])+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[26] = ser_2x2v_p2_surfvx_quad_26(1, fSkin); 
  } else { 

    fUpwindQuad[26] = ser_2x2v_p2_surfvx_quad_26(-1, fEdge); 
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

  Ghat[0] += 0.3535533905932737*(alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[10]*fUpwind[10]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*alpha[14]*fUpwind[18]+0.3162277660168379*alpha[10]*fUpwind[17]+0.3162277660168379*(alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])+0.3535533905932737*alpha[8]*fUpwind[12]+0.3162277660168379*alpha[4]*fUpwind[11]+0.3535533905932737*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3162277660168379*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7])+0.3535533905932737*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3162277660168379*alpha[10]*fUpwind[18]+0.3535533905932737*alpha[13]*fUpwind[17]+0.3162277660168379*(alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[4]*fUpwind[12])+0.3535533905932737*(alpha[7]*fUpwind[11]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3162277660168379*(alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8])+0.3535533905932737*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3162277660168379*alpha[10]*fUpwind[19]+0.3162277660168379*(alpha[6]*fUpwind[16]+alpha[5]*fUpwind[15])+0.3535533905932737*(alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3162277660168379*alpha[3]*fUpwind[9]+0.3535533905932737*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.3162277660168379*(alpha[6]*fUpwind[18]+alpha[5]*fUpwind[17])+0.3162277660168379*(alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[2]*fUpwind[12]+alpha[1]*fUpwind[11])+0.3535533905932737*(alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10])+0.3162277660168379*(alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7])+0.3535533905932737*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3162277660168379*alpha[6]*fUpwind[19]+0.3535533905932737*alpha[8]*fUpwind[18]+0.3162277660168379*alpha[4]*fUpwind[17]+0.3162277660168379*alpha[10]*fUpwind[16]+(0.2828427124746191*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[15]+0.3535533905932737*fUpwind[12]*alpha[14]+0.3162277660168379*(alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[10]*fUpwind[11])+0.3535533905932737*(alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10])+0.3162277660168379*(alpha[5]*(fUpwind[9]+fUpwind[7])+fUpwind[5]*alpha[7])+0.3535533905932737*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] += 0.3162277660168379*(alpha[5]*fUpwind[19]+alpha[4]*fUpwind[18])+0.3535533905932737*alpha[7]*fUpwind[17]+0.2828427124746191*alpha[14]*fUpwind[16]+0.3162277660168379*(alpha[3]*fUpwind[16]+alpha[10]*fUpwind[15]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14])+0.3535533905932737*fUpwind[11]*alpha[13]+0.3162277660168379*alpha[10]*fUpwind[12]+0.3535533905932737*(alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10])+0.3162277660168379*(alpha[6]*(fUpwind[9]+fUpwind[8])+fUpwind[6]*alpha[8])+0.3535533905932737*(alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] += 0.3535533905932737*alpha[6]*fUpwind[17]+0.2258769757263128*alpha[13]*fUpwind[13]+0.3535533905932737*(alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[2]*fUpwind[11])+0.3162277660168379*alpha[10]*fUpwind[10]+0.2258769757263128*alpha[7]*fUpwind[7]+0.3535533905932737*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7])+0.3162277660168379*(alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[1]*fUpwind[1]); 
  Ghat[8] += 0.3535533905932737*alpha[5]*fUpwind[18]+0.2258769757263128*alpha[14]*fUpwind[14]+0.3535533905932737*(alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[1]*fUpwind[12])+0.3162277660168379*alpha[10]*fUpwind[10]+0.2258769757263128*alpha[8]*fUpwind[8]+0.3535533905932737*(alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8])+0.3162277660168379*(alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]); 
  Ghat[9] += 0.3535533905932737*(alpha[4]*fUpwind[19]+alpha[2]*fUpwind[16]+alpha[1]*fUpwind[15])+0.3162277660168379*(alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[10]*fUpwind[10])+0.3535533905932737*alpha[0]*fUpwind[9]+0.3162277660168379*(alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[3]*fUpwind[3]); 
  Ghat[10] += 0.282842712474619*(alpha[14]+alpha[13])*fUpwind[19]+0.3162277660168379*(alpha[3]*fUpwind[19]+alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17])+0.3162277660168379*(alpha[5]*fUpwind[16]+alpha[6]*fUpwind[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[6]*fUpwind[12]+alpha[5]*fUpwind[11])+(0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[10]+0.3162277660168379*(fUpwind[9]+fUpwind[8]+fUpwind[7])*alpha[10]+0.3535533905932737*(fUpwind[0]*alpha[10]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.282842712474619*alpha[10]*fUpwind[18]+(0.3162277660168379*alpha[14]+0.2258769757263128*alpha[13])*fUpwind[17]+0.3535533905932737*(alpha[3]*fUpwind[17]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13])+0.2828427124746191*alpha[4]*fUpwind[12]+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[11]+0.3162277660168379*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3535533905932737*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[12] += (0.2258769757263128*alpha[14]+0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])*fUpwind[18]+0.282842712474619*alpha[10]*fUpwind[17]+0.3535533905932737*(alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14])+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[12]+0.2828427124746191*alpha[4]*fUpwind[11]+0.3162277660168379*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3535533905932737*(alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[13] += 0.282842712474619*alpha[10]*fUpwind[19]+0.3535533905932737*alpha[2]*fUpwind[17]+0.2828427124746191*alpha[5]*fUpwind[15]+(0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[13]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[7])*alpha[13]+0.3535533905932737*(fUpwind[0]*alpha[13]+alpha[6]*fUpwind[11])+0.3162277660168379*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[14] += 0.282842712474619*alpha[10]*fUpwind[19]+0.3535533905932737*alpha[1]*fUpwind[18]+0.2828427124746191*alpha[6]*fUpwind[16]+(0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind[14]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[8])*alpha[14]+0.3535533905932737*(fUpwind[0]*alpha[14]+alpha[5]*fUpwind[12])+0.3162277660168379*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]); 
  Ghat[15] += 0.3535533905932737*alpha[2]*fUpwind[19]+0.3162277660168379*alpha[14]*fUpwind[18]+0.282842712474619*alpha[10]*fUpwind[17]+0.3535533905932737*alpha[4]*fUpwind[16]+(0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[15]+0.2828427124746191*(alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])+0.3162277660168379*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3535533905932737*alpha[1]*fUpwind[9]+0.3162277660168379*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[16] += 0.3535533905932737*alpha[1]*fUpwind[19]+0.282842712474619*alpha[10]*fUpwind[18]+0.3162277660168379*(alpha[13]*fUpwind[17]+alpha[8]*fUpwind[16])+0.3535533905932737*(alpha[0]*fUpwind[16]+alpha[4]*fUpwind[15])+0.2828427124746191*(alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14])+0.3162277660168379*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3535533905932737*alpha[2]*fUpwind[9]+0.3162277660168379*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] += 0.2828427124746191*(alpha[5]*fUpwind[19]+alpha[4]*fUpwind[18])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[17]+0.3162277660168379*alpha[13]*fUpwind[16]+0.282842712474619*alpha[10]*fUpwind[15]+0.3162277660168379*fUpwind[11]*alpha[14]+0.3535533905932737*alpha[2]*fUpwind[13]+(0.2258769757263128*fUpwind[11]+0.3535533905932737*fUpwind[2])*alpha[13]+0.282842712474619*alpha[10]*fUpwind[12]+0.3535533905932737*alpha[3]*fUpwind[11]+0.3162277660168379*(alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10])+0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7])+0.3162277660168379*(alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[18] += 0.2828427124746191*alpha[6]*fUpwind[19]+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[18]+0.2828427124746191*alpha[4]*fUpwind[17]+0.282842712474619*alpha[10]*fUpwind[16]+0.3162277660168379*alpha[14]*fUpwind[15]+0.3535533905932737*alpha[1]*fUpwind[14]+(0.2258769757263128*fUpwind[12]+0.3535533905932737*fUpwind[1])*alpha[14]+fUpwind[12]*(0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])+0.282842712474619*alpha[10]*fUpwind[11]+0.3162277660168379*(alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10])+0.3535533905932737*(alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8])+0.3162277660168379*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[19] += (0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[19]+0.2828427124746191*(alpha[6]*fUpwind[18]+alpha[5]*fUpwind[17])+0.3535533905932737*(alpha[1]*fUpwind[16]+alpha[2]*fUpwind[15])+0.282842712474619*(alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13])+0.3162277660168379*(alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10])+0.3535533905932737*alpha[4]*fUpwind[9]+0.3162277660168379*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 

  out[0] += -0.7071067811865475*Ghat[0]*dv10; 
  out[1] += -0.7071067811865475*Ghat[1]*dv10; 
  out[2] += -0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += -0.7071067811865475*Ghat[3]*dv10; 
  out[5] += -0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -1.224744871391589*Ghat[1]*dv10; 
  out[7] += -1.224744871391589*Ghat[2]*dv10; 
  out[8] += -0.7071067811865475*Ghat[5]*dv10; 
  out[9] += -0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += -0.7071067811865475*Ghat[7]*dv10; 
  out[12] += -0.7071067811865475*Ghat[8]*dv10; 
  out[13] += -1.58113883008419*Ghat[0]*dv10; 
  out[14] += -0.7071067811865475*Ghat[9]*dv10; 
  out[15] += -1.224744871391589*Ghat[4]*dv10; 
  out[16] += -0.7071067811865475*Ghat[10]*dv10; 
  out[17] += -1.224744871391589*Ghat[5]*dv10; 
  out[18] += -1.224744871391589*Ghat[6]*dv10; 
  out[19] += -0.7071067811865475*Ghat[11]*dv10; 
  out[20] += -0.7071067811865475*Ghat[12]*dv10; 
  out[21] += -1.224744871391589*Ghat[7]*dv10; 
  out[22] += -1.224744871391589*Ghat[8]*dv10; 
  out[23] += -1.58113883008419*Ghat[1]*dv10; 
  out[24] += -1.58113883008419*Ghat[2]*dv10; 
  out[25] += -0.7071067811865475*Ghat[13]*dv10; 
  out[26] += -0.7071067811865475*Ghat[14]*dv10; 
  out[27] += -1.58113883008419*Ghat[3]*dv10; 
  out[28] += -0.7071067811865475*Ghat[15]*dv10; 
  out[29] += -0.7071067811865475*Ghat[16]*dv10; 
  out[30] += -1.224744871391589*Ghat[9]*dv10; 
  out[31] += -1.224744871391589*Ghat[10]*dv10; 
  out[32] += -1.224744871391589*Ghat[11]*dv10; 
  out[33] += -1.224744871391589*Ghat[12]*dv10; 
  out[34] += -1.58113883008419*Ghat[4]*dv10; 
  out[35] += -0.7071067811865475*Ghat[17]*dv10; 
  out[36] += -0.7071067811865475*Ghat[18]*dv10; 
  out[37] += -1.224744871391589*Ghat[13]*dv10; 
  out[38] += -1.224744871391589*Ghat[14]*dv10; 
  out[39] += -1.58113883008419*Ghat[5]*dv10; 
  out[40] += -1.58113883008419*Ghat[6]*dv10; 
  out[41] += -0.7071067811865475*Ghat[19]*dv10; 
  out[42] += -1.224744871391589*Ghat[15]*dv10; 
  out[43] += -1.224744871391589*Ghat[16]*dv10; 
  out[44] += -1.224744871391589*Ghat[17]*dv10; 
  out[45] += -1.224744871391589*Ghat[18]*dv10; 
  out[46] += -1.58113883008419*Ghat[10]*dv10; 
  out[47] += -1.224744871391589*Ghat[19]*dv10; 

  } else { 

  if ((-0.4242640687119285*(alpha[14]+alpha[13]))-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[0] = ser_2x2v_p2_surfvx_quad_0(1, fEdge); 
  } else { 

    fUpwindQuad[0] = ser_2x2v_p2_surfvx_quad_0(-1, fSkin); 
  } 
  if ((-0.4242640687119285*alpha[14])+0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]-0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[1] = ser_2x2v_p2_surfvx_quad_1(1, fEdge); 
  } else { 

    fUpwindQuad[1] = ser_2x2v_p2_surfvx_quad_1(-1, fSkin); 
  } 
  if ((-0.4242640687119285*(alpha[14]+alpha[13]))+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])-0.4743416490252568*(alpha[3]+alpha[2])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[2] = ser_2x2v_p2_surfvx_quad_2(1, fEdge); 
  } else { 

    fUpwindQuad[2] = ser_2x2v_p2_surfvx_quad_2(-1, fSkin); 
  } 
  if (0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]-0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[3] = ser_2x2v_p2_surfvx_quad_3(1, fEdge); 
  } else { 

    fUpwindQuad[3] = ser_2x2v_p2_surfvx_quad_3(-1, fSkin); 
  } 
  if (0.5303300858899104*(alpha[14]+alpha[13])-0.3952847075210473*(alpha[8]+alpha[7])-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[4] = ser_2x2v_p2_surfvx_quad_4(1, fEdge); 
  } else { 

    fUpwindQuad[4] = ser_2x2v_p2_surfvx_quad_4(-1, fSkin); 
  } 
  if (0.5303300858899104*alpha[14]-0.4242640687119285*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[5] = ser_2x2v_p2_surfvx_quad_5(1, fEdge); 
  } else { 

    fUpwindQuad[5] = ser_2x2v_p2_surfvx_quad_5(-1, fSkin); 
  } 
  if ((-0.4242640687119285*(alpha[14]+alpha[13]))+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[6] = ser_2x2v_p2_surfvx_quad_6(1, fEdge); 
  } else { 

    fUpwindQuad[6] = ser_2x2v_p2_surfvx_quad_6(-1, fSkin); 
  } 
  if ((-0.4242640687119285*alpha[14])+0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]-0.4743416490252568*alpha[3]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[7] = ser_2x2v_p2_surfvx_quad_7(1, fEdge); 
  } else { 

    fUpwindQuad[7] = ser_2x2v_p2_surfvx_quad_7(-1, fSkin); 
  } 
  if ((-0.4242640687119285*(alpha[14]+alpha[13]))-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]-0.4743416490252568*alpha[3]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[8] = ser_2x2v_p2_surfvx_quad_8(1, fEdge); 
  } else { 

    fUpwindQuad[8] = ser_2x2v_p2_surfvx_quad_8(-1, fSkin); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[9] = ser_2x2v_p2_surfvx_quad_9(1, fEdge); 
  } else { 

    fUpwindQuad[9] = ser_2x2v_p2_surfvx_quad_9(-1, fSkin); 
  } 
  if (0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[10] = ser_2x2v_p2_surfvx_quad_10(1, fEdge); 
  } else { 

    fUpwindQuad[10] = ser_2x2v_p2_surfvx_quad_10(-1, fSkin); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[11] = ser_2x2v_p2_surfvx_quad_11(1, fEdge); 
  } else { 

    fUpwindQuad[11] = ser_2x2v_p2_surfvx_quad_11(-1, fSkin); 
  } 
  if ((-0.3952847075210473*alpha[8])+0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[12] = ser_2x2v_p2_surfvx_quad_12(1, fEdge); 
  } else { 

    fUpwindQuad[12] = ser_2x2v_p2_surfvx_quad_12(-1, fSkin); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*(alpha[8]+alpha[7]) > 0) { 

    fUpwindQuad[13] = ser_2x2v_p2_surfvx_quad_13(1, fEdge); 
  } else { 

    fUpwindQuad[13] = ser_2x2v_p2_surfvx_quad_13(-1, fSkin); 
  } 
  if ((-0.3952847075210473*alpha[8])+0.3162277660168379*alpha[7]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[14] = ser_2x2v_p2_surfvx_quad_14(1, fEdge); 
  } else { 

    fUpwindQuad[14] = ser_2x2v_p2_surfvx_quad_14(-1, fSkin); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[15] = ser_2x2v_p2_surfvx_quad_15(1, fEdge); 
  } else { 

    fUpwindQuad[15] = ser_2x2v_p2_surfvx_quad_15(-1, fSkin); 
  } 
  if (0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[16] = ser_2x2v_p2_surfvx_quad_16(1, fEdge); 
  } else { 

    fUpwindQuad[16] = ser_2x2v_p2_surfvx_quad_16(-1, fSkin); 
  } 
  if (0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[17] = ser_2x2v_p2_surfvx_quad_17(1, fEdge); 
  } else { 

    fUpwindQuad[17] = ser_2x2v_p2_surfvx_quad_17(-1, fSkin); 
  } 
  if (0.4242640687119285*(alpha[14]+alpha[13])+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*(alpha[6]+alpha[5])+0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[18] = ser_2x2v_p2_surfvx_quad_18(1, fEdge); 
  } else { 

    fUpwindQuad[18] = ser_2x2v_p2_surfvx_quad_18(-1, fSkin); 
  } 
  if (0.4242640687119285*alpha[14]-0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]-0.6363961030678926*alpha[6]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[19] = ser_2x2v_p2_surfvx_quad_19(1, fEdge); 
  } else { 

    fUpwindQuad[19] = ser_2x2v_p2_surfvx_quad_19(-1, fSkin); 
  } 
  if (0.4242640687119285*(alpha[14]+alpha[13])-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])-0.6363961030678926*alpha[6]+0.6363961030678926*alpha[5]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[20] = ser_2x2v_p2_surfvx_quad_20(1, fEdge); 
  } else { 

    fUpwindQuad[20] = ser_2x2v_p2_surfvx_quad_20(-1, fSkin); 
  } 
  if ((-0.5303300858899104*alpha[14])+0.4242640687119285*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]-0.6363961030678926*alpha[5]+0.4743416490252568*alpha[3]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[21] = ser_2x2v_p2_surfvx_quad_21(1, fEdge); 
  } else { 

    fUpwindQuad[21] = ser_2x2v_p2_surfvx_quad_21(-1, fSkin); 
  } 
  if ((-0.5303300858899104*(alpha[14]+alpha[13]))-0.3952847075210473*(alpha[8]+alpha[7])+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[22] = ser_2x2v_p2_surfvx_quad_22(1, fEdge); 
  } else { 

    fUpwindQuad[22] = ser_2x2v_p2_surfvx_quad_22(-1, fSkin); 
  } 
  if ((-0.5303300858899104*alpha[14])+0.4242640687119285*alpha[13]-0.3952847075210473*alpha[8]+0.3162277660168379*alpha[7]+0.6363961030678926*alpha[5]+0.4743416490252568*(alpha[3]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[23] = ser_2x2v_p2_surfvx_quad_23(1, fEdge); 
  } else { 

    fUpwindQuad[23] = ser_2x2v_p2_surfvx_quad_23(-1, fSkin); 
  } 
  if (0.4242640687119285*(alpha[14]+alpha[13])-0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*alpha[6]-0.6363961030678926*(alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[24] = ser_2x2v_p2_surfvx_quad_24(1, fEdge); 
  } else { 

    fUpwindQuad[24] = ser_2x2v_p2_surfvx_quad_24(-1, fSkin); 
  } 
  if (0.4242640687119285*alpha[14]-0.5303300858899104*alpha[13]+0.3162277660168379*alpha[8]-0.3952847075210473*alpha[7]+0.6363961030678926*alpha[6]+0.4743416490252568*(alpha[3]+alpha[2])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[25] = ser_2x2v_p2_surfvx_quad_25(1, fEdge); 
  } else { 

    fUpwindQuad[25] = ser_2x2v_p2_surfvx_quad_25(-1, fSkin); 
  } 
  if (0.4242640687119285*(alpha[14]+alpha[13])+0.853814968245462*alpha[10]+0.3162277660168379*(alpha[8]+alpha[7])+0.6363961030678926*(alpha[6]+alpha[5]+alpha[4])+0.4743416490252568*(alpha[3]+alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad[26] = ser_2x2v_p2_surfvx_quad_26(1, fEdge); 
  } else { 

    fUpwindQuad[26] = ser_2x2v_p2_surfvx_quad_26(-1, fSkin); 
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

  Ghat[0] += 0.3535533905932737*(alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[10]*fUpwind[10]+alpha[8]*fUpwind[8]+alpha[7]*fUpwind[7]+alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[3]*fUpwind[3]+alpha[2]*fUpwind[2]+alpha[1]*fUpwind[1]+alpha[0]*fUpwind[0]); 
  Ghat[1] += 0.3535533905932737*alpha[14]*fUpwind[18]+0.3162277660168379*alpha[10]*fUpwind[17]+0.3162277660168379*(alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])+0.3535533905932737*alpha[8]*fUpwind[12]+0.3162277660168379*alpha[4]*fUpwind[11]+0.3535533905932737*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3162277660168379*(alpha[1]*fUpwind[7]+fUpwind[1]*alpha[7])+0.3535533905932737*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]+alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]+alpha[0]*fUpwind[1]+fUpwind[0]*alpha[1]); 
  Ghat[2] += 0.3162277660168379*alpha[10]*fUpwind[18]+0.3535533905932737*alpha[13]*fUpwind[17]+0.3162277660168379*(alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14]+alpha[4]*fUpwind[12])+0.3535533905932737*(alpha[7]*fUpwind[11]+alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3162277660168379*(alpha[2]*fUpwind[8]+fUpwind[2]*alpha[8])+0.3535533905932737*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]+alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]+alpha[0]*fUpwind[2]+fUpwind[0]*alpha[2]); 
  Ghat[3] += 0.3162277660168379*alpha[10]*fUpwind[19]+0.3162277660168379*(alpha[6]*fUpwind[16]+alpha[5]*fUpwind[15])+0.3535533905932737*(alpha[8]*fUpwind[14]+fUpwind[8]*alpha[14]+alpha[7]*fUpwind[13]+fUpwind[7]*alpha[13]+alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3162277660168379*alpha[3]*fUpwind[9]+0.3535533905932737*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]+alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]+alpha[0]*fUpwind[3]+fUpwind[0]*alpha[3]); 
  Ghat[4] += 0.3162277660168379*(alpha[6]*fUpwind[18]+alpha[5]*fUpwind[17])+0.3162277660168379*(alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13]+alpha[2]*fUpwind[12]+alpha[1]*fUpwind[11])+0.3535533905932737*(alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10])+0.3162277660168379*(alpha[4]*fUpwind[8]+fUpwind[4]*alpha[8]+alpha[4]*fUpwind[7]+fUpwind[4]*alpha[7])+0.3535533905932737*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]+alpha[0]*fUpwind[4]+fUpwind[0]*alpha[4]+alpha[1]*fUpwind[2]+fUpwind[1]*alpha[2]); 
  Ghat[5] += 0.3162277660168379*alpha[6]*fUpwind[19]+0.3535533905932737*alpha[8]*fUpwind[18]+0.3162277660168379*alpha[4]*fUpwind[17]+0.3162277660168379*alpha[10]*fUpwind[16]+(0.2828427124746191*alpha[13]+0.3162277660168379*alpha[3])*fUpwind[15]+0.3535533905932737*fUpwind[12]*alpha[14]+0.3162277660168379*(alpha[1]*fUpwind[13]+fUpwind[1]*alpha[13]+alpha[10]*fUpwind[11])+0.3535533905932737*(alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10])+0.3162277660168379*(alpha[5]*(fUpwind[9]+fUpwind[7])+fUpwind[5]*alpha[7])+0.3535533905932737*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]+alpha[0]*fUpwind[5]+fUpwind[0]*alpha[5]+alpha[1]*fUpwind[3]+fUpwind[1]*alpha[3]); 
  Ghat[6] += 0.3162277660168379*(alpha[5]*fUpwind[19]+alpha[4]*fUpwind[18])+0.3535533905932737*alpha[7]*fUpwind[17]+0.2828427124746191*alpha[14]*fUpwind[16]+0.3162277660168379*(alpha[3]*fUpwind[16]+alpha[10]*fUpwind[15]+alpha[2]*fUpwind[14]+fUpwind[2]*alpha[14])+0.3535533905932737*fUpwind[11]*alpha[13]+0.3162277660168379*alpha[10]*fUpwind[12]+0.3535533905932737*(alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10])+0.3162277660168379*(alpha[6]*(fUpwind[9]+fUpwind[8])+fUpwind[6]*alpha[8])+0.3535533905932737*(alpha[0]*fUpwind[6]+fUpwind[0]*alpha[6]+alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]+alpha[2]*fUpwind[3]+fUpwind[2]*alpha[3]); 
  Ghat[7] += 0.3535533905932737*alpha[6]*fUpwind[17]+0.2258769757263128*alpha[13]*fUpwind[13]+0.3535533905932737*(alpha[3]*fUpwind[13]+fUpwind[3]*alpha[13]+alpha[2]*fUpwind[11])+0.3162277660168379*alpha[10]*fUpwind[10]+0.2258769757263128*alpha[7]*fUpwind[7]+0.3535533905932737*(alpha[0]*fUpwind[7]+fUpwind[0]*alpha[7])+0.3162277660168379*(alpha[5]*fUpwind[5]+alpha[4]*fUpwind[4]+alpha[1]*fUpwind[1]); 
  Ghat[8] += 0.3535533905932737*alpha[5]*fUpwind[18]+0.2258769757263128*alpha[14]*fUpwind[14]+0.3535533905932737*(alpha[3]*fUpwind[14]+fUpwind[3]*alpha[14]+alpha[1]*fUpwind[12])+0.3162277660168379*alpha[10]*fUpwind[10]+0.2258769757263128*alpha[8]*fUpwind[8]+0.3535533905932737*(alpha[0]*fUpwind[8]+fUpwind[0]*alpha[8])+0.3162277660168379*(alpha[6]*fUpwind[6]+alpha[4]*fUpwind[4]+alpha[2]*fUpwind[2]); 
  Ghat[9] += 0.3535533905932737*(alpha[4]*fUpwind[19]+alpha[2]*fUpwind[16]+alpha[1]*fUpwind[15])+0.3162277660168379*(alpha[14]*fUpwind[14]+alpha[13]*fUpwind[13]+alpha[10]*fUpwind[10])+0.3535533905932737*alpha[0]*fUpwind[9]+0.3162277660168379*(alpha[6]*fUpwind[6]+alpha[5]*fUpwind[5]+alpha[3]*fUpwind[3]); 
  Ghat[10] += 0.282842712474619*(alpha[14]+alpha[13])*fUpwind[19]+0.3162277660168379*(alpha[3]*fUpwind[19]+alpha[2]*fUpwind[18]+alpha[1]*fUpwind[17])+0.3162277660168379*(alpha[5]*fUpwind[16]+alpha[6]*fUpwind[15]+alpha[4]*fUpwind[14]+fUpwind[4]*alpha[14]+alpha[4]*fUpwind[13]+fUpwind[4]*alpha[13]+alpha[6]*fUpwind[12]+alpha[5]*fUpwind[11])+(0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[10]+0.3162277660168379*(fUpwind[9]+fUpwind[8]+fUpwind[7])*alpha[10]+0.3535533905932737*(fUpwind[0]*alpha[10]+alpha[1]*fUpwind[6]+fUpwind[1]*alpha[6]+alpha[2]*fUpwind[5]+fUpwind[2]*alpha[5]+alpha[3]*fUpwind[4]+fUpwind[3]*alpha[4]); 
  Ghat[11] += 0.282842712474619*alpha[10]*fUpwind[18]+(0.3162277660168379*alpha[14]+0.2258769757263128*alpha[13])*fUpwind[17]+0.3535533905932737*(alpha[3]*fUpwind[17]+alpha[6]*fUpwind[13]+fUpwind[6]*alpha[13])+0.2828427124746191*alpha[4]*fUpwind[12]+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[11]+0.3162277660168379*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3535533905932737*(alpha[2]*fUpwind[7]+fUpwind[2]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind[4]+fUpwind[1]*alpha[4]); 
  Ghat[12] += (0.2258769757263128*alpha[14]+0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])*fUpwind[18]+0.282842712474619*alpha[10]*fUpwind[17]+0.3535533905932737*(alpha[5]*fUpwind[14]+fUpwind[5]*alpha[14])+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[12]+0.2828427124746191*alpha[4]*fUpwind[11]+0.3162277660168379*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3535533905932737*(alpha[1]*fUpwind[8]+fUpwind[1]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind[4]+fUpwind[2]*alpha[4]); 
  Ghat[13] += 0.282842712474619*alpha[10]*fUpwind[19]+0.3535533905932737*alpha[2]*fUpwind[17]+0.2828427124746191*alpha[5]*fUpwind[15]+(0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[13]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[7])*alpha[13]+0.3535533905932737*(fUpwind[0]*alpha[13]+alpha[6]*fUpwind[11])+0.3162277660168379*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind[7]+fUpwind[3]*alpha[7])+0.3162277660168379*(alpha[1]*fUpwind[5]+fUpwind[1]*alpha[5]); 
  Ghat[14] += 0.282842712474619*alpha[10]*fUpwind[19]+0.3535533905932737*alpha[1]*fUpwind[18]+0.2828427124746191*alpha[6]*fUpwind[16]+(0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind[14]+(0.3162277660168379*fUpwind[9]+0.2258769757263128*fUpwind[8])*alpha[14]+0.3535533905932737*(fUpwind[0]*alpha[14]+alpha[5]*fUpwind[12])+0.3162277660168379*(alpha[4]*fUpwind[10]+fUpwind[4]*alpha[10])+0.3535533905932737*(alpha[3]*fUpwind[8]+fUpwind[3]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind[6]+fUpwind[2]*alpha[6]); 
  Ghat[15] += 0.3535533905932737*alpha[2]*fUpwind[19]+0.3162277660168379*alpha[14]*fUpwind[18]+0.282842712474619*alpha[10]*fUpwind[17]+0.3535533905932737*alpha[4]*fUpwind[16]+(0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[15]+0.2828427124746191*(alpha[5]*fUpwind[13]+fUpwind[5]*alpha[13])+0.3162277660168379*(alpha[6]*fUpwind[10]+fUpwind[6]*alpha[10])+0.3535533905932737*alpha[1]*fUpwind[9]+0.3162277660168379*(alpha[3]*fUpwind[5]+fUpwind[3]*alpha[5]); 
  Ghat[16] += 0.3535533905932737*alpha[1]*fUpwind[19]+0.282842712474619*alpha[10]*fUpwind[18]+0.3162277660168379*(alpha[13]*fUpwind[17]+alpha[8]*fUpwind[16])+0.3535533905932737*(alpha[0]*fUpwind[16]+alpha[4]*fUpwind[15])+0.2828427124746191*(alpha[6]*fUpwind[14]+fUpwind[6]*alpha[14])+0.3162277660168379*(alpha[5]*fUpwind[10]+fUpwind[5]*alpha[10])+0.3535533905932737*alpha[2]*fUpwind[9]+0.3162277660168379*(alpha[3]*fUpwind[6]+fUpwind[3]*alpha[6]); 
  Ghat[17] += 0.2828427124746191*(alpha[5]*fUpwind[19]+alpha[4]*fUpwind[18])+(0.3162277660168379*alpha[8]+0.2258769757263128*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[17]+0.3162277660168379*alpha[13]*fUpwind[16]+0.282842712474619*alpha[10]*fUpwind[15]+0.3162277660168379*fUpwind[11]*alpha[14]+0.3535533905932737*alpha[2]*fUpwind[13]+(0.2258769757263128*fUpwind[11]+0.3535533905932737*fUpwind[2])*alpha[13]+0.282842712474619*alpha[10]*fUpwind[12]+0.3535533905932737*alpha[3]*fUpwind[11]+0.3162277660168379*(alpha[1]*fUpwind[10]+fUpwind[1]*alpha[10])+0.3535533905932737*(alpha[6]*fUpwind[7]+fUpwind[6]*alpha[7])+0.3162277660168379*(alpha[4]*fUpwind[5]+fUpwind[4]*alpha[5]); 
  Ghat[18] += 0.2828427124746191*alpha[6]*fUpwind[19]+(0.2258769757263128*alpha[8]+0.3162277660168379*alpha[7]+0.3535533905932737*alpha[0])*fUpwind[18]+0.2828427124746191*alpha[4]*fUpwind[17]+0.282842712474619*alpha[10]*fUpwind[16]+0.3162277660168379*alpha[14]*fUpwind[15]+0.3535533905932737*alpha[1]*fUpwind[14]+(0.2258769757263128*fUpwind[12]+0.3535533905932737*fUpwind[1])*alpha[14]+fUpwind[12]*(0.3162277660168379*alpha[13]+0.3535533905932737*alpha[3])+0.282842712474619*alpha[10]*fUpwind[11]+0.3162277660168379*(alpha[2]*fUpwind[10]+fUpwind[2]*alpha[10])+0.3535533905932737*(alpha[5]*fUpwind[8]+fUpwind[5]*alpha[8])+0.3162277660168379*(alpha[4]*fUpwind[6]+fUpwind[4]*alpha[6]); 
  Ghat[19] += (0.3162277660168379*(alpha[8]+alpha[7])+0.3535533905932737*alpha[0])*fUpwind[19]+0.2828427124746191*(alpha[6]*fUpwind[18]+alpha[5]*fUpwind[17])+0.3535533905932737*(alpha[1]*fUpwind[16]+alpha[2]*fUpwind[15])+0.282842712474619*(alpha[10]*fUpwind[14]+fUpwind[10]*alpha[14]+alpha[10]*fUpwind[13]+fUpwind[10]*alpha[13])+0.3162277660168379*(alpha[3]*fUpwind[10]+fUpwind[3]*alpha[10])+0.3535533905932737*alpha[4]*fUpwind[9]+0.3162277660168379*(alpha[5]*fUpwind[6]+fUpwind[5]*alpha[6]); 

  out[0] += 0.7071067811865475*Ghat[0]*dv10; 
  out[1] += 0.7071067811865475*Ghat[1]*dv10; 
  out[2] += 0.7071067811865475*Ghat[2]*dv10; 
  out[3] += -1.224744871391589*Ghat[0]*dv10; 
  out[4] += 0.7071067811865475*Ghat[3]*dv10; 
  out[5] += 0.7071067811865475*Ghat[4]*dv10; 
  out[6] += -1.224744871391589*Ghat[1]*dv10; 
  out[7] += -1.224744871391589*Ghat[2]*dv10; 
  out[8] += 0.7071067811865475*Ghat[5]*dv10; 
  out[9] += 0.7071067811865475*Ghat[6]*dv10; 
  out[10] += -1.224744871391589*Ghat[3]*dv10; 
  out[11] += 0.7071067811865475*Ghat[7]*dv10; 
  out[12] += 0.7071067811865475*Ghat[8]*dv10; 
  out[13] += 1.58113883008419*Ghat[0]*dv10; 
  out[14] += 0.7071067811865475*Ghat[9]*dv10; 
  out[15] += -1.224744871391589*Ghat[4]*dv10; 
  out[16] += 0.7071067811865475*Ghat[10]*dv10; 
  out[17] += -1.224744871391589*Ghat[5]*dv10; 
  out[18] += -1.224744871391589*Ghat[6]*dv10; 
  out[19] += 0.7071067811865475*Ghat[11]*dv10; 
  out[20] += 0.7071067811865475*Ghat[12]*dv10; 
  out[21] += -1.224744871391589*Ghat[7]*dv10; 
  out[22] += -1.224744871391589*Ghat[8]*dv10; 
  out[23] += 1.58113883008419*Ghat[1]*dv10; 
  out[24] += 1.58113883008419*Ghat[2]*dv10; 
  out[25] += 0.7071067811865475*Ghat[13]*dv10; 
  out[26] += 0.7071067811865475*Ghat[14]*dv10; 
  out[27] += 1.58113883008419*Ghat[3]*dv10; 
  out[28] += 0.7071067811865475*Ghat[15]*dv10; 
  out[29] += 0.7071067811865475*Ghat[16]*dv10; 
  out[30] += -1.224744871391589*Ghat[9]*dv10; 
  out[31] += -1.224744871391589*Ghat[10]*dv10; 
  out[32] += -1.224744871391589*Ghat[11]*dv10; 
  out[33] += -1.224744871391589*Ghat[12]*dv10; 
  out[34] += 1.58113883008419*Ghat[4]*dv10; 
  out[35] += 0.7071067811865475*Ghat[17]*dv10; 
  out[36] += 0.7071067811865475*Ghat[18]*dv10; 
  out[37] += -1.224744871391589*Ghat[13]*dv10; 
  out[38] += -1.224744871391589*Ghat[14]*dv10; 
  out[39] += 1.58113883008419*Ghat[5]*dv10; 
  out[40] += 1.58113883008419*Ghat[6]*dv10; 
  out[41] += 0.7071067811865475*Ghat[19]*dv10; 
  out[42] += -1.224744871391589*Ghat[15]*dv10; 
  out[43] += -1.224744871391589*Ghat[16]*dv10; 
  out[44] += -1.224744871391589*Ghat[17]*dv10; 
  out[45] += -1.224744871391589*Ghat[18]*dv10; 
  out[46] += 1.58113883008419*Ghat[10]*dv10; 
  out[47] += -1.224744871391589*Ghat[19]*dv10; 

  } 
} 
