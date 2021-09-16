#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x2v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_poisson_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  double alpha[20] = {0.0}; 

  alpha[0] = -2.449489742783178*phi[1]*dx10; 
  alpha[1] = -5.477225575051662*phi[4]*dx10; 
  alpha[2] = -2.449489742783178*phi[3]*dx10; 
  alpha[4] = -5.477225575051662*phi[6]*dx10; 
  alpha[8] = -2.449489742783178*phi[7]*dx10; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};;
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 

  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[0] = ser_2x2v_p2_surfvx_quad_0(1, fl); 
    fUpwindQuad_r[0] = ser_2x2v_p2_surfvx_quad_0(1, fc); 
  } else { 

    fUpwindQuad_l[0] = ser_2x2v_p2_surfvx_quad_0(-1, fc); 
    fUpwindQuad_r[0] = ser_2x2v_p2_surfvx_quad_0(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[1] = ser_2x2v_p2_surfvx_quad_1(1, fl); 
    fUpwindQuad_r[1] = ser_2x2v_p2_surfvx_quad_1(1, fc); 
  } else { 

    fUpwindQuad_l[1] = ser_2x2v_p2_surfvx_quad_1(-1, fc); 
    fUpwindQuad_r[1] = ser_2x2v_p2_surfvx_quad_1(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[2] = ser_2x2v_p2_surfvx_quad_2(1, fl); 
    fUpwindQuad_r[2] = ser_2x2v_p2_surfvx_quad_2(1, fc); 
  } else { 

    fUpwindQuad_l[2] = ser_2x2v_p2_surfvx_quad_2(-1, fc); 
    fUpwindQuad_r[2] = ser_2x2v_p2_surfvx_quad_2(-1, fr); 
  } 
  if ((-0.3952847075210473*alpha[8])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[3] = ser_2x2v_p2_surfvx_quad_3(1, fl); 
    fUpwindQuad_r[3] = ser_2x2v_p2_surfvx_quad_3(1, fc); 
  } else { 

    fUpwindQuad_l[3] = ser_2x2v_p2_surfvx_quad_3(-1, fc); 
    fUpwindQuad_r[3] = ser_2x2v_p2_surfvx_quad_3(-1, fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*alpha[8] > 0) { 

    fUpwindQuad_l[4] = ser_2x2v_p2_surfvx_quad_4(1, fl); 
    fUpwindQuad_r[4] = ser_2x2v_p2_surfvx_quad_4(1, fc); 
  } else { 

    fUpwindQuad_l[4] = ser_2x2v_p2_surfvx_quad_4(-1, fc); 
    fUpwindQuad_r[4] = ser_2x2v_p2_surfvx_quad_4(-1, fr); 
  } 
  if ((-0.3952847075210473*alpha[8])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[5] = ser_2x2v_p2_surfvx_quad_5(1, fl); 
    fUpwindQuad_r[5] = ser_2x2v_p2_surfvx_quad_5(1, fc); 
  } else { 

    fUpwindQuad_l[5] = ser_2x2v_p2_surfvx_quad_5(-1, fc); 
    fUpwindQuad_r[5] = ser_2x2v_p2_surfvx_quad_5(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[6] = ser_2x2v_p2_surfvx_quad_6(1, fl); 
    fUpwindQuad_r[6] = ser_2x2v_p2_surfvx_quad_6(1, fc); 
  } else { 

    fUpwindQuad_l[6] = ser_2x2v_p2_surfvx_quad_6(-1, fc); 
    fUpwindQuad_r[6] = ser_2x2v_p2_surfvx_quad_6(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[7] = ser_2x2v_p2_surfvx_quad_7(1, fl); 
    fUpwindQuad_r[7] = ser_2x2v_p2_surfvx_quad_7(1, fc); 
  } else { 

    fUpwindQuad_l[7] = ser_2x2v_p2_surfvx_quad_7(-1, fc); 
    fUpwindQuad_r[7] = ser_2x2v_p2_surfvx_quad_7(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[8] = ser_2x2v_p2_surfvx_quad_8(1, fl); 
    fUpwindQuad_r[8] = ser_2x2v_p2_surfvx_quad_8(1, fc); 
  } else { 

    fUpwindQuad_l[8] = ser_2x2v_p2_surfvx_quad_8(-1, fc); 
    fUpwindQuad_r[8] = ser_2x2v_p2_surfvx_quad_8(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[9] = ser_2x2v_p2_surfvx_quad_9(1, fl); 
    fUpwindQuad_r[9] = ser_2x2v_p2_surfvx_quad_9(1, fc); 
  } else { 

    fUpwindQuad_l[9] = ser_2x2v_p2_surfvx_quad_9(-1, fc); 
    fUpwindQuad_r[9] = ser_2x2v_p2_surfvx_quad_9(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[10] = ser_2x2v_p2_surfvx_quad_10(1, fl); 
    fUpwindQuad_r[10] = ser_2x2v_p2_surfvx_quad_10(1, fc); 
  } else { 

    fUpwindQuad_l[10] = ser_2x2v_p2_surfvx_quad_10(-1, fc); 
    fUpwindQuad_r[10] = ser_2x2v_p2_surfvx_quad_10(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[11] = ser_2x2v_p2_surfvx_quad_11(1, fl); 
    fUpwindQuad_r[11] = ser_2x2v_p2_surfvx_quad_11(1, fc); 
  } else { 

    fUpwindQuad_l[11] = ser_2x2v_p2_surfvx_quad_11(-1, fc); 
    fUpwindQuad_r[11] = ser_2x2v_p2_surfvx_quad_11(-1, fr); 
  } 
  if ((-0.3952847075210473*alpha[8])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[12] = ser_2x2v_p2_surfvx_quad_12(1, fl); 
    fUpwindQuad_r[12] = ser_2x2v_p2_surfvx_quad_12(1, fc); 
  } else { 

    fUpwindQuad_l[12] = ser_2x2v_p2_surfvx_quad_12(-1, fc); 
    fUpwindQuad_r[12] = ser_2x2v_p2_surfvx_quad_12(-1, fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*alpha[8] > 0) { 

    fUpwindQuad_l[13] = ser_2x2v_p2_surfvx_quad_13(1, fl); 
    fUpwindQuad_r[13] = ser_2x2v_p2_surfvx_quad_13(1, fc); 
  } else { 

    fUpwindQuad_l[13] = ser_2x2v_p2_surfvx_quad_13(-1, fc); 
    fUpwindQuad_r[13] = ser_2x2v_p2_surfvx_quad_13(-1, fr); 
  } 
  if ((-0.3952847075210473*alpha[8])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[14] = ser_2x2v_p2_surfvx_quad_14(1, fl); 
    fUpwindQuad_r[14] = ser_2x2v_p2_surfvx_quad_14(1, fc); 
  } else { 

    fUpwindQuad_l[14] = ser_2x2v_p2_surfvx_quad_14(-1, fc); 
    fUpwindQuad_r[14] = ser_2x2v_p2_surfvx_quad_14(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[15] = ser_2x2v_p2_surfvx_quad_15(1, fl); 
    fUpwindQuad_r[15] = ser_2x2v_p2_surfvx_quad_15(1, fc); 
  } else { 

    fUpwindQuad_l[15] = ser_2x2v_p2_surfvx_quad_15(-1, fc); 
    fUpwindQuad_r[15] = ser_2x2v_p2_surfvx_quad_15(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[16] = ser_2x2v_p2_surfvx_quad_16(1, fl); 
    fUpwindQuad_r[16] = ser_2x2v_p2_surfvx_quad_16(1, fc); 
  } else { 

    fUpwindQuad_l[16] = ser_2x2v_p2_surfvx_quad_16(-1, fc); 
    fUpwindQuad_r[16] = ser_2x2v_p2_surfvx_quad_16(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[17] = ser_2x2v_p2_surfvx_quad_17(1, fl); 
    fUpwindQuad_r[17] = ser_2x2v_p2_surfvx_quad_17(1, fc); 
  } else { 

    fUpwindQuad_l[17] = ser_2x2v_p2_surfvx_quad_17(-1, fc); 
    fUpwindQuad_r[17] = ser_2x2v_p2_surfvx_quad_17(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[18] = ser_2x2v_p2_surfvx_quad_18(1, fl); 
    fUpwindQuad_r[18] = ser_2x2v_p2_surfvx_quad_18(1, fc); 
  } else { 

    fUpwindQuad_l[18] = ser_2x2v_p2_surfvx_quad_18(-1, fc); 
    fUpwindQuad_r[18] = ser_2x2v_p2_surfvx_quad_18(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[19] = ser_2x2v_p2_surfvx_quad_19(1, fl); 
    fUpwindQuad_r[19] = ser_2x2v_p2_surfvx_quad_19(1, fc); 
  } else { 

    fUpwindQuad_l[19] = ser_2x2v_p2_surfvx_quad_19(-1, fc); 
    fUpwindQuad_r[19] = ser_2x2v_p2_surfvx_quad_19(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]-0.4743416490252568*alpha[2]+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[20] = ser_2x2v_p2_surfvx_quad_20(1, fl); 
    fUpwindQuad_r[20] = ser_2x2v_p2_surfvx_quad_20(1, fc); 
  } else { 

    fUpwindQuad_l[20] = ser_2x2v_p2_surfvx_quad_20(-1, fc); 
    fUpwindQuad_r[20] = ser_2x2v_p2_surfvx_quad_20(-1, fr); 
  } 
  if ((-0.3952847075210473*alpha[8])-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[21] = ser_2x2v_p2_surfvx_quad_21(1, fl); 
    fUpwindQuad_r[21] = ser_2x2v_p2_surfvx_quad_21(1, fc); 
  } else { 

    fUpwindQuad_l[21] = ser_2x2v_p2_surfvx_quad_21(-1, fc); 
    fUpwindQuad_r[21] = ser_2x2v_p2_surfvx_quad_21(-1, fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3952847075210473*alpha[8] > 0) { 

    fUpwindQuad_l[22] = ser_2x2v_p2_surfvx_quad_22(1, fl); 
    fUpwindQuad_r[22] = ser_2x2v_p2_surfvx_quad_22(1, fc); 
  } else { 

    fUpwindQuad_l[22] = ser_2x2v_p2_surfvx_quad_22(-1, fc); 
    fUpwindQuad_r[22] = ser_2x2v_p2_surfvx_quad_22(-1, fr); 
  } 
  if ((-0.3952847075210473*alpha[8])+0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[23] = ser_2x2v_p2_surfvx_quad_23(1, fl); 
    fUpwindQuad_r[23] = ser_2x2v_p2_surfvx_quad_23(1, fc); 
  } else { 

    fUpwindQuad_l[23] = ser_2x2v_p2_surfvx_quad_23(-1, fc); 
    fUpwindQuad_r[23] = ser_2x2v_p2_surfvx_quad_23(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[24] = ser_2x2v_p2_surfvx_quad_24(1, fl); 
    fUpwindQuad_r[24] = ser_2x2v_p2_surfvx_quad_24(1, fc); 
  } else { 

    fUpwindQuad_l[24] = ser_2x2v_p2_surfvx_quad_24(-1, fc); 
    fUpwindQuad_r[24] = ser_2x2v_p2_surfvx_quad_24(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]+0.4743416490252568*alpha[2]+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[25] = ser_2x2v_p2_surfvx_quad_25(1, fl); 
    fUpwindQuad_r[25] = ser_2x2v_p2_surfvx_quad_25(1, fc); 
  } else { 

    fUpwindQuad_l[25] = ser_2x2v_p2_surfvx_quad_25(-1, fc); 
    fUpwindQuad_r[25] = ser_2x2v_p2_surfvx_quad_25(-1, fr); 
  } 
  if (0.3162277660168379*alpha[8]+0.6363961030678926*alpha[4]+0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 

    fUpwindQuad_l[26] = ser_2x2v_p2_surfvx_quad_26(1, fl); 
    fUpwindQuad_r[26] = ser_2x2v_p2_surfvx_quad_26(1, fc); 
  } else { 

    fUpwindQuad_l[26] = ser_2x2v_p2_surfvx_quad_26(-1, fc); 
    fUpwindQuad_r[26] = ser_2x2v_p2_surfvx_quad_26(-1, fr); 
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

  Ghat_l[0] += 0.3535533905932737*(alpha[8]*fUpwind_l[8]+alpha[4]*fUpwind_l[4]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.3535533905932737*alpha[8]*fUpwind_l[12]+0.3162277660168379*alpha[4]*fUpwind_l[11]+0.3162277660168379*alpha[1]*fUpwind_l[7]+0.3535533905932737*(alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.3162277660168379*alpha[4]*fUpwind_l[12]+0.3162277660168379*(alpha[2]*fUpwind_l[8]+fUpwind_l[2]*alpha[8])+0.3535533905932737*(alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.3535533905932737*(alpha[8]*fUpwind_l[14]+alpha[4]*fUpwind_l[10]+alpha[2]*fUpwind_l[6]+alpha[1]*fUpwind_l[5]+alpha[0]*fUpwind_l[3]); 
  Ghat_l[4] += 0.3162277660168379*(alpha[2]*fUpwind_l[12]+alpha[1]*fUpwind_l[11])+0.3162277660168379*(alpha[4]*fUpwind_l[8]+fUpwind_l[4]*alpha[8]+alpha[4]*fUpwind_l[7])+0.3535533905932737*(alpha[0]*fUpwind_l[4]+fUpwind_l[0]*alpha[4]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[5] += 0.3535533905932737*alpha[8]*fUpwind_l[18]+0.3162277660168379*alpha[4]*fUpwind_l[17]+0.3162277660168379*alpha[1]*fUpwind_l[13]+0.3535533905932737*(alpha[2]*fUpwind_l[10]+alpha[4]*fUpwind_l[6]+alpha[0]*fUpwind_l[5]+alpha[1]*fUpwind_l[3]); 
  Ghat_l[6] += 0.3162277660168379*alpha[4]*fUpwind_l[18]+0.3162277660168379*alpha[2]*fUpwind_l[14]+0.3535533905932737*alpha[1]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[6]*alpha[8]+0.3535533905932737*(alpha[0]*fUpwind_l[6]+alpha[4]*fUpwind_l[5]+alpha[2]*fUpwind_l[3]); 
  Ghat_l[7] += 0.3535533905932737*(alpha[2]*fUpwind_l[11]+alpha[0]*fUpwind_l[7])+0.3162277660168379*(alpha[4]*fUpwind_l[4]+alpha[1]*fUpwind_l[1]); 
  Ghat_l[8] += 0.3535533905932737*alpha[1]*fUpwind_l[12]+0.2258769757263128*alpha[8]*fUpwind_l[8]+0.3535533905932737*(alpha[0]*fUpwind_l[8]+fUpwind_l[0]*alpha[8])+0.3162277660168379*(alpha[4]*fUpwind_l[4]+alpha[2]*fUpwind_l[2]); 
  Ghat_l[9] += 0.3535533905932737*(alpha[4]*fUpwind_l[19]+alpha[2]*fUpwind_l[16]+alpha[1]*fUpwind_l[15]+alpha[0]*fUpwind_l[9]); 
  Ghat_l[10] += 0.3162277660168379*(alpha[2]*fUpwind_l[18]+alpha[1]*fUpwind_l[17])+0.3162277660168379*alpha[4]*(fUpwind_l[14]+fUpwind_l[13])+0.3162277660168379*alpha[8]*fUpwind_l[10]+0.3535533905932737*(alpha[0]*fUpwind_l[10]+alpha[1]*fUpwind_l[6]+alpha[2]*fUpwind_l[5]+fUpwind_l[3]*alpha[4]); 
  Ghat_l[11] += 0.2828427124746191*alpha[4]*fUpwind_l[12]+0.3162277660168379*alpha[8]*fUpwind_l[11]+0.3535533905932737*(alpha[0]*fUpwind_l[11]+alpha[2]*fUpwind_l[7])+0.3162277660168379*(alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]); 
  Ghat_l[12] += (0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind_l[12]+0.2828427124746191*alpha[4]*fUpwind_l[11]+0.3535533905932737*(alpha[1]*fUpwind_l[8]+fUpwind_l[1]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]); 
  Ghat_l[13] += 0.3535533905932737*(alpha[2]*fUpwind_l[17]+alpha[0]*fUpwind_l[13])+0.3162277660168379*(alpha[4]*fUpwind_l[10]+alpha[1]*fUpwind_l[5]); 
  Ghat_l[14] += 0.3535533905932737*alpha[1]*fUpwind_l[18]+(0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind_l[14]+0.3162277660168379*alpha[4]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[3]*alpha[8]+0.3162277660168379*alpha[2]*fUpwind_l[6]; 
  Ghat_l[15] += 0.3535533905932737*(alpha[2]*fUpwind_l[19]+alpha[4]*fUpwind_l[16]+alpha[0]*fUpwind_l[15]+alpha[1]*fUpwind_l[9]); 
  Ghat_l[16] += 0.3535533905932737*alpha[1]*fUpwind_l[19]+0.3162277660168379*alpha[8]*fUpwind_l[16]+0.3535533905932737*(alpha[0]*fUpwind_l[16]+alpha[4]*fUpwind_l[15]+alpha[2]*fUpwind_l[9]); 
  Ghat_l[17] += 0.2828427124746191*alpha[4]*fUpwind_l[18]+0.3162277660168379*alpha[8]*fUpwind_l[17]+0.3535533905932737*(alpha[0]*fUpwind_l[17]+alpha[2]*fUpwind_l[13])+0.3162277660168379*(alpha[1]*fUpwind_l[10]+alpha[4]*fUpwind_l[5]); 
  Ghat_l[18] += (0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind_l[18]+0.2828427124746191*alpha[4]*fUpwind_l[17]+0.3535533905932737*alpha[1]*fUpwind_l[14]+0.3162277660168379*alpha[2]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[5]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_l[6]; 
  Ghat_l[19] += 0.3162277660168379*alpha[8]*fUpwind_l[19]+0.3535533905932737*(alpha[0]*fUpwind_l[19]+alpha[1]*fUpwind_l[16]+alpha[2]*fUpwind_l[15]+alpha[4]*fUpwind_l[9]); 

  Ghat_r[0] += 0.3535533905932737*(alpha[8]*fUpwind_r[8]+alpha[4]*fUpwind_r[4]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.3535533905932737*alpha[8]*fUpwind_r[12]+0.3162277660168379*alpha[4]*fUpwind_r[11]+0.3162277660168379*alpha[1]*fUpwind_r[7]+0.3535533905932737*(alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.3162277660168379*alpha[4]*fUpwind_r[12]+0.3162277660168379*(alpha[2]*fUpwind_r[8]+fUpwind_r[2]*alpha[8])+0.3535533905932737*(alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.3535533905932737*(alpha[8]*fUpwind_r[14]+alpha[4]*fUpwind_r[10]+alpha[2]*fUpwind_r[6]+alpha[1]*fUpwind_r[5]+alpha[0]*fUpwind_r[3]); 
  Ghat_r[4] += 0.3162277660168379*(alpha[2]*fUpwind_r[12]+alpha[1]*fUpwind_r[11])+0.3162277660168379*(alpha[4]*fUpwind_r[8]+fUpwind_r[4]*alpha[8]+alpha[4]*fUpwind_r[7])+0.3535533905932737*(alpha[0]*fUpwind_r[4]+fUpwind_r[0]*alpha[4]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[5] += 0.3535533905932737*alpha[8]*fUpwind_r[18]+0.3162277660168379*alpha[4]*fUpwind_r[17]+0.3162277660168379*alpha[1]*fUpwind_r[13]+0.3535533905932737*(alpha[2]*fUpwind_r[10]+alpha[4]*fUpwind_r[6]+alpha[0]*fUpwind_r[5]+alpha[1]*fUpwind_r[3]); 
  Ghat_r[6] += 0.3162277660168379*alpha[4]*fUpwind_r[18]+0.3162277660168379*alpha[2]*fUpwind_r[14]+0.3535533905932737*alpha[1]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[6]*alpha[8]+0.3535533905932737*(alpha[0]*fUpwind_r[6]+alpha[4]*fUpwind_r[5]+alpha[2]*fUpwind_r[3]); 
  Ghat_r[7] += 0.3535533905932737*(alpha[2]*fUpwind_r[11]+alpha[0]*fUpwind_r[7])+0.3162277660168379*(alpha[4]*fUpwind_r[4]+alpha[1]*fUpwind_r[1]); 
  Ghat_r[8] += 0.3535533905932737*alpha[1]*fUpwind_r[12]+0.2258769757263128*alpha[8]*fUpwind_r[8]+0.3535533905932737*(alpha[0]*fUpwind_r[8]+fUpwind_r[0]*alpha[8])+0.3162277660168379*(alpha[4]*fUpwind_r[4]+alpha[2]*fUpwind_r[2]); 
  Ghat_r[9] += 0.3535533905932737*(alpha[4]*fUpwind_r[19]+alpha[2]*fUpwind_r[16]+alpha[1]*fUpwind_r[15]+alpha[0]*fUpwind_r[9]); 
  Ghat_r[10] += 0.3162277660168379*(alpha[2]*fUpwind_r[18]+alpha[1]*fUpwind_r[17])+0.3162277660168379*alpha[4]*(fUpwind_r[14]+fUpwind_r[13])+0.3162277660168379*alpha[8]*fUpwind_r[10]+0.3535533905932737*(alpha[0]*fUpwind_r[10]+alpha[1]*fUpwind_r[6]+alpha[2]*fUpwind_r[5]+fUpwind_r[3]*alpha[4]); 
  Ghat_r[11] += 0.2828427124746191*alpha[4]*fUpwind_r[12]+0.3162277660168379*alpha[8]*fUpwind_r[11]+0.3535533905932737*(alpha[0]*fUpwind_r[11]+alpha[2]*fUpwind_r[7])+0.3162277660168379*(alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]); 
  Ghat_r[12] += (0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind_r[12]+0.2828427124746191*alpha[4]*fUpwind_r[11]+0.3535533905932737*(alpha[1]*fUpwind_r[8]+fUpwind_r[1]*alpha[8])+0.3162277660168379*(alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]); 
  Ghat_r[13] += 0.3535533905932737*(alpha[2]*fUpwind_r[17]+alpha[0]*fUpwind_r[13])+0.3162277660168379*(alpha[4]*fUpwind_r[10]+alpha[1]*fUpwind_r[5]); 
  Ghat_r[14] += 0.3535533905932737*alpha[1]*fUpwind_r[18]+(0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind_r[14]+0.3162277660168379*alpha[4]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[3]*alpha[8]+0.3162277660168379*alpha[2]*fUpwind_r[6]; 
  Ghat_r[15] += 0.3535533905932737*(alpha[2]*fUpwind_r[19]+alpha[4]*fUpwind_r[16]+alpha[0]*fUpwind_r[15]+alpha[1]*fUpwind_r[9]); 
  Ghat_r[16] += 0.3535533905932737*alpha[1]*fUpwind_r[19]+0.3162277660168379*alpha[8]*fUpwind_r[16]+0.3535533905932737*(alpha[0]*fUpwind_r[16]+alpha[4]*fUpwind_r[15]+alpha[2]*fUpwind_r[9]); 
  Ghat_r[17] += 0.2828427124746191*alpha[4]*fUpwind_r[18]+0.3162277660168379*alpha[8]*fUpwind_r[17]+0.3535533905932737*(alpha[0]*fUpwind_r[17]+alpha[2]*fUpwind_r[13])+0.3162277660168379*(alpha[1]*fUpwind_r[10]+alpha[4]*fUpwind_r[5]); 
  Ghat_r[18] += (0.2258769757263128*alpha[8]+0.3535533905932737*alpha[0])*fUpwind_r[18]+0.2828427124746191*alpha[4]*fUpwind_r[17]+0.3535533905932737*alpha[1]*fUpwind_r[14]+0.3162277660168379*alpha[2]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[5]*alpha[8]+0.3162277660168379*alpha[4]*fUpwind_r[6]; 
  Ghat_r[19] += 0.3162277660168379*alpha[8]*fUpwind_r[19]+0.3535533905932737*(alpha[0]*fUpwind_r[19]+alpha[1]*fUpwind_r[16]+alpha[2]*fUpwind_r[15]+alpha[4]*fUpwind_r[9]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[6] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[13] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv10; 
  out[14] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[15] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[16] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[17] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[19] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[20] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[21] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 
  out[22] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10; 
  out[23] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv10; 
  out[24] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv10; 
  out[25] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[26] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[27] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv10; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[29] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*dv10; 
  out[30] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10; 
  out[31] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10; 
  out[32] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10; 
  out[33] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10; 
  out[34] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv10; 
  out[35] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*dv10; 
  out[36] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*dv10; 
  out[37] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10; 
  out[38] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10; 
  out[39] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv10; 
  out[40] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv10; 
  out[41] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*dv10; 
  out[42] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10; 
  out[43] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*dv10; 
  out[44] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*dv10; 
  out[45] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*dv10; 
  out[46] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*dv10; 
  out[47] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*dv10; 

} 
