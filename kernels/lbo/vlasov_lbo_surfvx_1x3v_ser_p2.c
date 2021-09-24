#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x3v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_lbo_surfvx_1x3v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         cell-center coordinates. 
  // dxv[4]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[9]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf_l[20] = {0.0}; 
  alphaDrSurf_l[0] = 2.828427124746191*nuSum[0]*w[1]-1.414213562373095*nuSum[0]*dxv[1]-2.0*sumNuUx[0]; 
  alphaDrSurf_l[1] = -2.0*sumNuUx[1]; 
  alphaDrSurf_l[7] = -2.0*sumNuUx[2]; 

  double alphaDrSurf_r[20] = {0.0}; 
  alphaDrSurf_r[0] = 2.828427124746191*nuSum[0]*w[1]+1.414213562373095*nuSum[0]*dxv[1]-2.0*sumNuUx[0]; 
  alphaDrSurf_r[1] = -2.0*sumNuUx[1]; 
  alphaDrSurf_r[7] = -2.0*sumNuUx[2]; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};;
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 
  double Gdiff_l[20] = {0.0}; 
  double Gdiff_r[20] = {0.0}; 
  double Gdiff2_l[20] = {0.0}; 
  double Gdiff2_r[20] = {0.0}; 

  if (0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_1x3v_p2_surfvx_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_1x3v_p2_surfvx_quad_0(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*alphaDrSurf_l[7] > 0) { 
    fUpwindQuad_l[1] = ser_1x3v_p2_surfvx_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_1x3v_p2_surfvx_quad_1(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_1x3v_p2_surfvx_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_1x3v_p2_surfvx_quad_2(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_1x3v_p2_surfvx_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = ser_1x3v_p2_surfvx_quad_3(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*alphaDrSurf_l[7] > 0) { 
    fUpwindQuad_l[4] = ser_1x3v_p2_surfvx_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = ser_1x3v_p2_surfvx_quad_4(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = ser_1x3v_p2_surfvx_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = ser_1x3v_p2_surfvx_quad_5(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[6] = ser_1x3v_p2_surfvx_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = ser_1x3v_p2_surfvx_quad_6(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*alphaDrSurf_l[7] > 0) { 
    fUpwindQuad_l[7] = ser_1x3v_p2_surfvx_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = ser_1x3v_p2_surfvx_quad_7(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[8] = ser_1x3v_p2_surfvx_quad_8(1, fl); 
  } else { 
    fUpwindQuad_l[8] = ser_1x3v_p2_surfvx_quad_8(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[9] = ser_1x3v_p2_surfvx_quad_9(1, fl); 
  } else { 
    fUpwindQuad_l[9] = ser_1x3v_p2_surfvx_quad_9(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*alphaDrSurf_l[7] > 0) { 
    fUpwindQuad_l[10] = ser_1x3v_p2_surfvx_quad_10(1, fl); 
  } else { 
    fUpwindQuad_l[10] = ser_1x3v_p2_surfvx_quad_10(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[11] = ser_1x3v_p2_surfvx_quad_11(1, fl); 
  } else { 
    fUpwindQuad_l[11] = ser_1x3v_p2_surfvx_quad_11(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[12] = ser_1x3v_p2_surfvx_quad_12(1, fl); 
  } else { 
    fUpwindQuad_l[12] = ser_1x3v_p2_surfvx_quad_12(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*alphaDrSurf_l[7] > 0) { 
    fUpwindQuad_l[13] = ser_1x3v_p2_surfvx_quad_13(1, fl); 
  } else { 
    fUpwindQuad_l[13] = ser_1x3v_p2_surfvx_quad_13(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[14] = ser_1x3v_p2_surfvx_quad_14(1, fl); 
  } else { 
    fUpwindQuad_l[14] = ser_1x3v_p2_surfvx_quad_14(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[15] = ser_1x3v_p2_surfvx_quad_15(1, fl); 
  } else { 
    fUpwindQuad_l[15] = ser_1x3v_p2_surfvx_quad_15(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*alphaDrSurf_l[7] > 0) { 
    fUpwindQuad_l[16] = ser_1x3v_p2_surfvx_quad_16(1, fl); 
  } else { 
    fUpwindQuad_l[16] = ser_1x3v_p2_surfvx_quad_16(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[17] = ser_1x3v_p2_surfvx_quad_17(1, fl); 
  } else { 
    fUpwindQuad_l[17] = ser_1x3v_p2_surfvx_quad_17(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[18] = ser_1x3v_p2_surfvx_quad_18(1, fl); 
  } else { 
    fUpwindQuad_l[18] = ser_1x3v_p2_surfvx_quad_18(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*alphaDrSurf_l[7] > 0) { 
    fUpwindQuad_l[19] = ser_1x3v_p2_surfvx_quad_19(1, fl); 
  } else { 
    fUpwindQuad_l[19] = ser_1x3v_p2_surfvx_quad_19(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[20] = ser_1x3v_p2_surfvx_quad_20(1, fl); 
  } else { 
    fUpwindQuad_l[20] = ser_1x3v_p2_surfvx_quad_20(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[21] = ser_1x3v_p2_surfvx_quad_21(1, fl); 
  } else { 
    fUpwindQuad_l[21] = ser_1x3v_p2_surfvx_quad_21(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*alphaDrSurf_l[7] > 0) { 
    fUpwindQuad_l[22] = ser_1x3v_p2_surfvx_quad_22(1, fl); 
  } else { 
    fUpwindQuad_l[22] = ser_1x3v_p2_surfvx_quad_22(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[23] = ser_1x3v_p2_surfvx_quad_23(1, fl); 
  } else { 
    fUpwindQuad_l[23] = ser_1x3v_p2_surfvx_quad_23(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[24] = ser_1x3v_p2_surfvx_quad_24(1, fl); 
  } else { 
    fUpwindQuad_l[24] = ser_1x3v_p2_surfvx_quad_24(-1, fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*alphaDrSurf_l[7] > 0) { 
    fUpwindQuad_l[25] = ser_1x3v_p2_surfvx_quad_25(1, fl); 
  } else { 
    fUpwindQuad_l[25] = ser_1x3v_p2_surfvx_quad_25(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[26] = ser_1x3v_p2_surfvx_quad_26(1, fl); 
  } else { 
    fUpwindQuad_l[26] = ser_1x3v_p2_surfvx_quad_26(-1, fc); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_1x3v_p2_surfvx_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_1x3v_p2_surfvx_quad_0(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*alphaDrSurf_r[7] > 0) { 
    fUpwindQuad_r[1] = ser_1x3v_p2_surfvx_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_1x3v_p2_surfvx_quad_1(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_1x3v_p2_surfvx_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_1x3v_p2_surfvx_quad_2(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_1x3v_p2_surfvx_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = ser_1x3v_p2_surfvx_quad_3(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*alphaDrSurf_r[7] > 0) { 
    fUpwindQuad_r[4] = ser_1x3v_p2_surfvx_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = ser_1x3v_p2_surfvx_quad_4(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = ser_1x3v_p2_surfvx_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = ser_1x3v_p2_surfvx_quad_5(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[6] = ser_1x3v_p2_surfvx_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = ser_1x3v_p2_surfvx_quad_6(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*alphaDrSurf_r[7] > 0) { 
    fUpwindQuad_r[7] = ser_1x3v_p2_surfvx_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = ser_1x3v_p2_surfvx_quad_7(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[8] = ser_1x3v_p2_surfvx_quad_8(1, fc); 
  } else { 
    fUpwindQuad_r[8] = ser_1x3v_p2_surfvx_quad_8(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[9] = ser_1x3v_p2_surfvx_quad_9(1, fc); 
  } else { 
    fUpwindQuad_r[9] = ser_1x3v_p2_surfvx_quad_9(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*alphaDrSurf_r[7] > 0) { 
    fUpwindQuad_r[10] = ser_1x3v_p2_surfvx_quad_10(1, fc); 
  } else { 
    fUpwindQuad_r[10] = ser_1x3v_p2_surfvx_quad_10(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[11] = ser_1x3v_p2_surfvx_quad_11(1, fc); 
  } else { 
    fUpwindQuad_r[11] = ser_1x3v_p2_surfvx_quad_11(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[12] = ser_1x3v_p2_surfvx_quad_12(1, fc); 
  } else { 
    fUpwindQuad_r[12] = ser_1x3v_p2_surfvx_quad_12(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*alphaDrSurf_r[7] > 0) { 
    fUpwindQuad_r[13] = ser_1x3v_p2_surfvx_quad_13(1, fc); 
  } else { 
    fUpwindQuad_r[13] = ser_1x3v_p2_surfvx_quad_13(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[14] = ser_1x3v_p2_surfvx_quad_14(1, fc); 
  } else { 
    fUpwindQuad_r[14] = ser_1x3v_p2_surfvx_quad_14(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[15] = ser_1x3v_p2_surfvx_quad_15(1, fc); 
  } else { 
    fUpwindQuad_r[15] = ser_1x3v_p2_surfvx_quad_15(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*alphaDrSurf_r[7] > 0) { 
    fUpwindQuad_r[16] = ser_1x3v_p2_surfvx_quad_16(1, fc); 
  } else { 
    fUpwindQuad_r[16] = ser_1x3v_p2_surfvx_quad_16(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[17] = ser_1x3v_p2_surfvx_quad_17(1, fc); 
  } else { 
    fUpwindQuad_r[17] = ser_1x3v_p2_surfvx_quad_17(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[18] = ser_1x3v_p2_surfvx_quad_18(1, fc); 
  } else { 
    fUpwindQuad_r[18] = ser_1x3v_p2_surfvx_quad_18(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*alphaDrSurf_r[7] > 0) { 
    fUpwindQuad_r[19] = ser_1x3v_p2_surfvx_quad_19(1, fc); 
  } else { 
    fUpwindQuad_r[19] = ser_1x3v_p2_surfvx_quad_19(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[20] = ser_1x3v_p2_surfvx_quad_20(1, fc); 
  } else { 
    fUpwindQuad_r[20] = ser_1x3v_p2_surfvx_quad_20(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[21] = ser_1x3v_p2_surfvx_quad_21(1, fc); 
  } else { 
    fUpwindQuad_r[21] = ser_1x3v_p2_surfvx_quad_21(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*alphaDrSurf_r[7] > 0) { 
    fUpwindQuad_r[22] = ser_1x3v_p2_surfvx_quad_22(1, fc); 
  } else { 
    fUpwindQuad_r[22] = ser_1x3v_p2_surfvx_quad_22(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[23] = ser_1x3v_p2_surfvx_quad_23(1, fc); 
  } else { 
    fUpwindQuad_r[23] = ser_1x3v_p2_surfvx_quad_23(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[24] = ser_1x3v_p2_surfvx_quad_24(1, fc); 
  } else { 
    fUpwindQuad_r[24] = ser_1x3v_p2_surfvx_quad_24(-1, fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*alphaDrSurf_r[7] > 0) { 
    fUpwindQuad_r[25] = ser_1x3v_p2_surfvx_quad_25(1, fc); 
  } else { 
    fUpwindQuad_r[25] = ser_1x3v_p2_surfvx_quad_25(-1, fr); 
  } 
  if (0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[26] = ser_1x3v_p2_surfvx_quad_26(1, fc); 
  } else { 
    fUpwindQuad_r[26] = ser_1x3v_p2_surfvx_quad_26(-1, fr); 
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

  Gdiff2_l[0] = 0.2445699350390395*nuVtSqSum[1]*fl[20]+0.2445699350390395*nuVtSqSum[1]*fc[20]+0.3518228202874282*nuVtSqSum[2]*fl[19]-0.3518228202874282*nuVtSqSum[2]*fc[19]+0.2445699350390395*nuVtSqSum[0]*fl[12]+0.2445699350390395*nuVtSqSum[0]*fc[12]+0.25*nuVtSqSum[2]*fl[11]+0.25*nuVtSqSum[2]*fc[11]+0.3518228202874282*nuVtSqSum[1]*fl[5]-0.3518228202874282*nuVtSqSum[1]*fc[5]+0.3518228202874282*nuVtSqSum[0]*fl[2]-0.3518228202874282*nuVtSqSum[0]*fc[2]+0.25*fl[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fl[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_l[1] = 0.21875*nuVtSqSum[2]*fl[20]+0.2445699350390395*nuVtSqSum[0]*fl[20]+0.21875*nuVtSqSum[2]*fc[20]+0.2445699350390395*nuVtSqSum[0]*fc[20]+0.3146798968793526*nuVtSqSum[1]*fl[19]-0.3146798968793526*nuVtSqSum[1]*fc[19]+0.2445699350390395*nuVtSqSum[1]*fl[12]+0.2445699350390395*nuVtSqSum[1]*fc[12]+0.223606797749979*nuVtSqSum[1]*fl[11]+0.223606797749979*nuVtSqSum[1]*fc[11]+0.3146798968793526*nuVtSqSum[2]*fl[5]+0.3518228202874282*nuVtSqSum[0]*fl[5]-0.3146798968793526*nuVtSqSum[2]*fc[5]-0.3518228202874282*nuVtSqSum[0]*fc[5]+0.223606797749979*fl[1]*nuVtSqSum[2]+0.223606797749979*fc[1]*nuVtSqSum[2]+0.3518228202874282*nuVtSqSum[1]*fl[2]-0.3518228202874282*nuVtSqSum[1]*fc[2]+0.25*fl[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fl[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_l[2] = 0.2445699350390395*nuVtSqSum[1]*fl[33]+0.2445699350390395*nuVtSqSum[1]*fc[33]+0.3518228202874282*nuVtSqSum[2]*fl[32]-0.3518228202874282*nuVtSqSum[2]*fc[32]+0.2445699350390395*nuVtSqSum[0]*fl[22]+0.2445699350390395*nuVtSqSum[0]*fc[22]+0.25*nuVtSqSum[2]*fl[21]+0.25*nuVtSqSum[2]*fc[21]+0.3518228202874282*nuVtSqSum[1]*fl[15]-0.3518228202874282*nuVtSqSum[1]*fc[15]+0.3518228202874282*nuVtSqSum[0]*fl[7]-0.3518228202874282*nuVtSqSum[0]*fc[7]+0.25*nuVtSqSum[1]*fl[6]+0.25*nuVtSqSum[1]*fc[6]+0.25*nuVtSqSum[0]*fl[3]+0.25*nuVtSqSum[0]*fc[3]; 
  Gdiff2_l[3] = 0.2445699350390395*nuVtSqSum[1]*fl[36]+0.2445699350390395*nuVtSqSum[1]*fc[36]+0.3518228202874282*nuVtSqSum[2]*fl[35]-0.3518228202874282*nuVtSqSum[2]*fc[35]+0.2445699350390395*nuVtSqSum[0]*fl[26]+0.2445699350390395*nuVtSqSum[0]*fc[26]+0.25*nuVtSqSum[2]*fl[25]+0.25*nuVtSqSum[2]*fc[25]+0.3518228202874282*nuVtSqSum[1]*fl[16]-0.3518228202874282*nuVtSqSum[1]*fc[16]+0.3518228202874282*nuVtSqSum[0]*fl[9]-0.3518228202874282*nuVtSqSum[0]*fc[9]+0.25*nuVtSqSum[1]*fl[8]+0.25*nuVtSqSum[1]*fc[8]+0.25*nuVtSqSum[0]*fl[4]+0.25*nuVtSqSum[0]*fc[4]; 
  Gdiff2_l[4] = 0.21875*nuVtSqSum[2]*fl[33]+0.2445699350390395*nuVtSqSum[0]*fl[33]+0.21875*nuVtSqSum[2]*fc[33]+0.2445699350390395*nuVtSqSum[0]*fc[33]+0.3146798968793526*nuVtSqSum[1]*fl[32]-0.3146798968793526*nuVtSqSum[1]*fc[32]+0.2445699350390395*nuVtSqSum[1]*fl[22]+0.2445699350390395*nuVtSqSum[1]*fc[22]+0.223606797749979*nuVtSqSum[1]*fl[21]+0.223606797749979*nuVtSqSum[1]*fc[21]+0.3146798968793526*nuVtSqSum[2]*fl[15]+0.3518228202874282*nuVtSqSum[0]*fl[15]-0.3146798968793526*nuVtSqSum[2]*fc[15]-0.3518228202874282*nuVtSqSum[0]*fc[15]+0.3518228202874282*nuVtSqSum[1]*fl[7]-0.3518228202874282*nuVtSqSum[1]*fc[7]+0.223606797749979*nuVtSqSum[2]*fl[6]+0.25*nuVtSqSum[0]*fl[6]+0.223606797749979*nuVtSqSum[2]*fc[6]+0.25*nuVtSqSum[0]*fc[6]+0.25*nuVtSqSum[1]*fl[3]+0.25*nuVtSqSum[1]*fc[3]; 
  Gdiff2_l[5] = 0.21875*nuVtSqSum[2]*fl[36]+0.2445699350390395*nuVtSqSum[0]*fl[36]+0.21875*nuVtSqSum[2]*fc[36]+0.2445699350390395*nuVtSqSum[0]*fc[36]+0.3146798968793526*nuVtSqSum[1]*fl[35]-0.3146798968793526*nuVtSqSum[1]*fc[35]+0.2445699350390395*nuVtSqSum[1]*fl[26]+0.2445699350390395*nuVtSqSum[1]*fc[26]+0.223606797749979*nuVtSqSum[1]*fl[25]+0.223606797749979*nuVtSqSum[1]*fc[25]+0.3146798968793526*nuVtSqSum[2]*fl[16]+0.3518228202874282*nuVtSqSum[0]*fl[16]-0.3146798968793526*nuVtSqSum[2]*fc[16]-0.3518228202874282*nuVtSqSum[0]*fc[16]+0.3518228202874282*nuVtSqSum[1]*fl[9]-0.3518228202874282*nuVtSqSum[1]*fc[9]+0.223606797749979*nuVtSqSum[2]*fl[8]+0.25*nuVtSqSum[0]*fl[8]+0.223606797749979*nuVtSqSum[2]*fc[8]+0.25*nuVtSqSum[0]*fc[8]+0.25*nuVtSqSum[1]*fl[4]+0.25*nuVtSqSum[1]*fc[4]; 
  Gdiff2_l[6] = 0.2445699350390395*nuVtSqSum[1]*fl[45]+0.2445699350390395*nuVtSqSum[1]*fc[45]+0.3518228202874282*nuVtSqSum[2]*fl[44]-0.3518228202874282*nuVtSqSum[2]*fc[44]+0.2445699350390395*nuVtSqSum[0]*fl[38]+0.2445699350390395*nuVtSqSum[0]*fc[38]+0.25*nuVtSqSum[2]*fl[37]+0.25*nuVtSqSum[2]*fc[37]+0.3518228202874282*nuVtSqSum[1]*fl[31]-0.3518228202874282*nuVtSqSum[1]*fc[31]+0.3518228202874282*nuVtSqSum[0]*fl[18]-0.3518228202874282*nuVtSqSum[0]*fc[18]+0.25*nuVtSqSum[1]*fl[17]+0.25*nuVtSqSum[1]*fc[17]+0.25*nuVtSqSum[0]*fl[10]+0.25*nuVtSqSum[0]*fc[10]; 
  Gdiff2_l[7] = 0.21875*nuVtSqSum[1]*fl[20]+0.21875*nuVtSqSum[1]*fc[20]+0.2247713549138233*nuVtSqSum[2]*fl[19]+0.3518228202874282*nuVtSqSum[0]*fl[19]-0.2247713549138233*nuVtSqSum[2]*fc[19]-0.3518228202874282*nuVtSqSum[0]*fc[19]+0.2445699350390395*nuVtSqSum[2]*fl[12]+0.2445699350390395*nuVtSqSum[2]*fc[12]+0.159719141249985*nuVtSqSum[2]*fl[11]+0.25*nuVtSqSum[0]*fl[11]+0.159719141249985*nuVtSqSum[2]*fc[11]+0.25*nuVtSqSum[0]*fc[11]+0.3146798968793526*nuVtSqSum[1]*fl[5]-0.3146798968793526*nuVtSqSum[1]*fc[5]+0.3518228202874282*fl[2]*nuVtSqSum[2]-0.3518228202874282*fc[2]*nuVtSqSum[2]+0.25*fl[0]*nuVtSqSum[2]+0.25*fc[0]*nuVtSqSum[2]+0.223606797749979*fl[1]*nuVtSqSum[1]+0.223606797749979*fc[1]*nuVtSqSum[1]; 
  Gdiff2_l[8] = 0.3518228202874282*nuVtSqSum[1]*fl[34]-0.3518228202874282*nuVtSqSum[1]*fc[34]+0.3518228202874282*nuVtSqSum[0]*fl[24]-0.3518228202874282*nuVtSqSum[0]*fc[24]+0.25*nuVtSqSum[1]*fl[23]+0.25*nuVtSqSum[1]*fc[23]+0.25*nuVtSqSum[0]*fl[13]+0.25*nuVtSqSum[0]*fc[13]; 
  Gdiff2_l[9] = 0.3518228202874282*nuVtSqSum[1]*fl[41]-0.3518228202874282*nuVtSqSum[1]*fc[41]+0.3518228202874282*nuVtSqSum[0]*fl[29]-0.3518228202874282*nuVtSqSum[0]*fc[29]+0.25*nuVtSqSum[1]*fl[28]+0.25*nuVtSqSum[1]*fc[28]+0.25*nuVtSqSum[0]*fl[14]+0.25*nuVtSqSum[0]*fc[14]; 
  Gdiff2_l[10] = 0.21875*nuVtSqSum[2]*fl[45]+0.2445699350390395*nuVtSqSum[0]*fl[45]+0.21875*nuVtSqSum[2]*fc[45]+0.2445699350390395*nuVtSqSum[0]*fc[45]+0.3146798968793526*nuVtSqSum[1]*fl[44]-0.3146798968793526*nuVtSqSum[1]*fc[44]+0.2445699350390395*nuVtSqSum[1]*fl[38]+0.2445699350390395*nuVtSqSum[1]*fc[38]+0.223606797749979*nuVtSqSum[1]*fl[37]+0.223606797749979*nuVtSqSum[1]*fc[37]+0.3146798968793526*nuVtSqSum[2]*fl[31]+0.3518228202874282*nuVtSqSum[0]*fl[31]-0.3146798968793526*nuVtSqSum[2]*fc[31]-0.3518228202874282*nuVtSqSum[0]*fc[31]+0.3518228202874282*nuVtSqSum[1]*fl[18]-0.3518228202874282*nuVtSqSum[1]*fc[18]+0.223606797749979*nuVtSqSum[2]*fl[17]+0.25*nuVtSqSum[0]*fl[17]+0.223606797749979*nuVtSqSum[2]*fc[17]+0.25*nuVtSqSum[0]*fc[17]+0.25*nuVtSqSum[1]*fl[10]+0.25*nuVtSqSum[1]*fc[10]; 
  Gdiff2_l[11] = 0.21875*nuVtSqSum[1]*fl[33]+0.21875*nuVtSqSum[1]*fc[33]+0.2247713549138233*nuVtSqSum[2]*fl[32]+0.3518228202874282*nuVtSqSum[0]*fl[32]-0.2247713549138233*nuVtSqSum[2]*fc[32]-0.3518228202874282*nuVtSqSum[0]*fc[32]+0.2445699350390395*nuVtSqSum[2]*fl[22]+0.2445699350390395*nuVtSqSum[2]*fc[22]+0.159719141249985*nuVtSqSum[2]*fl[21]+0.25*nuVtSqSum[0]*fl[21]+0.159719141249985*nuVtSqSum[2]*fc[21]+0.25*nuVtSqSum[0]*fc[21]+0.3146798968793526*nuVtSqSum[1]*fl[15]-0.3146798968793526*nuVtSqSum[1]*fc[15]+0.3518228202874282*nuVtSqSum[2]*fl[7]-0.3518228202874282*nuVtSqSum[2]*fc[7]+0.223606797749979*nuVtSqSum[1]*fl[6]+0.223606797749979*nuVtSqSum[1]*fc[6]+0.25*nuVtSqSum[2]*fl[3]+0.25*nuVtSqSum[2]*fc[3]; 
  Gdiff2_l[12] = 0.3146798968793526*nuVtSqSum[2]*fl[34]+0.3518228202874282*nuVtSqSum[0]*fl[34]-0.3146798968793526*nuVtSqSum[2]*fc[34]-0.3518228202874282*nuVtSqSum[0]*fc[34]+0.3518228202874282*nuVtSqSum[1]*fl[24]-0.3518228202874282*nuVtSqSum[1]*fc[24]+0.223606797749979*nuVtSqSum[2]*fl[23]+0.25*nuVtSqSum[0]*fl[23]+0.223606797749979*nuVtSqSum[2]*fc[23]+0.25*nuVtSqSum[0]*fc[23]+0.25*nuVtSqSum[1]*fl[13]+0.25*nuVtSqSum[1]*fc[13]; 
  Gdiff2_l[13] = 0.21875*nuVtSqSum[1]*fl[36]+0.21875*nuVtSqSum[1]*fc[36]+0.2247713549138233*nuVtSqSum[2]*fl[35]+0.3518228202874282*nuVtSqSum[0]*fl[35]-0.2247713549138233*nuVtSqSum[2]*fc[35]-0.3518228202874282*nuVtSqSum[0]*fc[35]+0.2445699350390395*nuVtSqSum[2]*fl[26]+0.2445699350390395*nuVtSqSum[2]*fc[26]+0.159719141249985*nuVtSqSum[2]*fl[25]+0.25*nuVtSqSum[0]*fl[25]+0.159719141249985*nuVtSqSum[2]*fc[25]+0.25*nuVtSqSum[0]*fc[25]+0.3146798968793526*nuVtSqSum[1]*fl[16]-0.3146798968793526*nuVtSqSum[1]*fc[16]+0.3518228202874282*nuVtSqSum[2]*fl[9]-0.3518228202874282*nuVtSqSum[2]*fc[9]+0.223606797749979*nuVtSqSum[1]*fl[8]+0.223606797749979*nuVtSqSum[1]*fc[8]+0.25*nuVtSqSum[2]*fl[4]+0.25*nuVtSqSum[2]*fc[4]; 
  Gdiff2_l[14] = 0.3518228202874282*nuVtSqSum[1]*fl[46]-0.3518228202874282*nuVtSqSum[1]*fc[46]+0.3518228202874282*nuVtSqSum[0]*fl[40]-0.3518228202874282*nuVtSqSum[0]*fc[40]+0.25*nuVtSqSum[1]*fl[39]+0.25*nuVtSqSum[1]*fc[39]+0.25*nuVtSqSum[0]*fl[27]+0.25*nuVtSqSum[0]*fc[27]; 
  Gdiff2_l[15] = 0.3146798968793526*nuVtSqSum[2]*fl[41]+0.3518228202874282*nuVtSqSum[0]*fl[41]-0.3146798968793526*nuVtSqSum[2]*fc[41]-0.3518228202874282*nuVtSqSum[0]*fc[41]+0.3518228202874282*nuVtSqSum[1]*fl[29]-0.3518228202874282*nuVtSqSum[1]*fc[29]+0.223606797749979*nuVtSqSum[2]*fl[28]+0.25*nuVtSqSum[0]*fl[28]+0.223606797749979*nuVtSqSum[2]*fc[28]+0.25*nuVtSqSum[0]*fc[28]+0.25*nuVtSqSum[1]*fl[14]+0.25*nuVtSqSum[1]*fc[14]; 
  Gdiff2_l[16] = 0.3518228202874282*nuVtSqSum[1]*fl[47]-0.3518228202874282*nuVtSqSum[1]*fc[47]+0.3518228202874282*nuVtSqSum[0]*fl[43]-0.3518228202874282*nuVtSqSum[0]*fc[43]+0.25*nuVtSqSum[1]*fl[42]+0.25*nuVtSqSum[1]*fc[42]+0.25*nuVtSqSum[0]*fl[30]+0.25*nuVtSqSum[0]*fc[30]; 
  Gdiff2_l[17] = 0.21875*nuVtSqSum[1]*fl[45]+0.21875*nuVtSqSum[1]*fc[45]+0.2247713549138233*nuVtSqSum[2]*fl[44]+0.3518228202874282*nuVtSqSum[0]*fl[44]-0.2247713549138233*nuVtSqSum[2]*fc[44]-0.3518228202874282*nuVtSqSum[0]*fc[44]+0.2445699350390395*nuVtSqSum[2]*fl[38]+0.2445699350390395*nuVtSqSum[2]*fc[38]+0.159719141249985*nuVtSqSum[2]*fl[37]+0.25*nuVtSqSum[0]*fl[37]+0.159719141249985*nuVtSqSum[2]*fc[37]+0.25*nuVtSqSum[0]*fc[37]+0.3146798968793526*nuVtSqSum[1]*fl[31]-0.3146798968793526*nuVtSqSum[1]*fc[31]+0.3518228202874282*nuVtSqSum[2]*fl[18]-0.3518228202874282*nuVtSqSum[2]*fc[18]+0.223606797749979*nuVtSqSum[1]*fl[17]+0.223606797749979*nuVtSqSum[1]*fc[17]+0.25*nuVtSqSum[2]*fl[10]+0.25*nuVtSqSum[2]*fc[10]; 
  Gdiff2_l[18] = 0.3146798968793526*nuVtSqSum[2]*fl[46]+0.3518228202874282*nuVtSqSum[0]*fl[46]-0.3146798968793526*nuVtSqSum[2]*fc[46]-0.3518228202874282*nuVtSqSum[0]*fc[46]+0.3518228202874282*nuVtSqSum[1]*fl[40]-0.3518228202874282*nuVtSqSum[1]*fc[40]+0.223606797749979*nuVtSqSum[2]*fl[39]+0.25*nuVtSqSum[0]*fl[39]+0.223606797749979*nuVtSqSum[2]*fc[39]+0.25*nuVtSqSum[0]*fc[39]+0.25*nuVtSqSum[1]*fl[27]+0.25*nuVtSqSum[1]*fc[27]; 
  Gdiff2_l[19] = 0.3146798968793526*nuVtSqSum[2]*fl[47]+0.3518228202874282*nuVtSqSum[0]*fl[47]-0.3146798968793526*nuVtSqSum[2]*fc[47]-0.3518228202874282*nuVtSqSum[0]*fc[47]+0.3518228202874282*nuVtSqSum[1]*fl[43]-0.3518228202874282*nuVtSqSum[1]*fc[43]+0.223606797749979*nuVtSqSum[2]*fl[42]+0.25*nuVtSqSum[0]*fl[42]+0.223606797749979*nuVtSqSum[2]*fc[42]+0.25*nuVtSqSum[0]*fc[42]+0.25*nuVtSqSum[1]*fl[30]+0.25*nuVtSqSum[1]*fc[30]; 

  Gdiff2_r[0] = 0.2445699350390395*nuVtSqSum[1]*fr[20]+0.2445699350390395*nuVtSqSum[1]*fc[20]-0.3518228202874282*nuVtSqSum[2]*fr[19]+0.3518228202874282*nuVtSqSum[2]*fc[19]+0.2445699350390395*nuVtSqSum[0]*fr[12]+0.2445699350390395*nuVtSqSum[0]*fc[12]+0.25*nuVtSqSum[2]*fr[11]+0.25*nuVtSqSum[2]*fc[11]-0.3518228202874282*nuVtSqSum[1]*fr[5]+0.3518228202874282*nuVtSqSum[1]*fc[5]-0.3518228202874282*nuVtSqSum[0]*fr[2]+0.3518228202874282*nuVtSqSum[0]*fc[2]+0.25*fr[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fr[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_r[1] = 0.21875*nuVtSqSum[2]*fr[20]+0.2445699350390395*nuVtSqSum[0]*fr[20]+0.21875*nuVtSqSum[2]*fc[20]+0.2445699350390395*nuVtSqSum[0]*fc[20]-0.3146798968793526*nuVtSqSum[1]*fr[19]+0.3146798968793526*nuVtSqSum[1]*fc[19]+0.2445699350390395*nuVtSqSum[1]*fr[12]+0.2445699350390395*nuVtSqSum[1]*fc[12]+0.223606797749979*nuVtSqSum[1]*fr[11]+0.223606797749979*nuVtSqSum[1]*fc[11]-0.3146798968793526*nuVtSqSum[2]*fr[5]-0.3518228202874282*nuVtSqSum[0]*fr[5]+0.3146798968793526*nuVtSqSum[2]*fc[5]+0.3518228202874282*nuVtSqSum[0]*fc[5]+0.223606797749979*fr[1]*nuVtSqSum[2]+0.223606797749979*fc[1]*nuVtSqSum[2]-0.3518228202874282*nuVtSqSum[1]*fr[2]+0.3518228202874282*nuVtSqSum[1]*fc[2]+0.25*fr[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fr[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_r[2] = 0.2445699350390395*nuVtSqSum[1]*fr[33]+0.2445699350390395*nuVtSqSum[1]*fc[33]-0.3518228202874282*nuVtSqSum[2]*fr[32]+0.3518228202874282*nuVtSqSum[2]*fc[32]+0.2445699350390395*nuVtSqSum[0]*fr[22]+0.2445699350390395*nuVtSqSum[0]*fc[22]+0.25*nuVtSqSum[2]*fr[21]+0.25*nuVtSqSum[2]*fc[21]-0.3518228202874282*nuVtSqSum[1]*fr[15]+0.3518228202874282*nuVtSqSum[1]*fc[15]-0.3518228202874282*nuVtSqSum[0]*fr[7]+0.3518228202874282*nuVtSqSum[0]*fc[7]+0.25*nuVtSqSum[1]*fr[6]+0.25*nuVtSqSum[1]*fc[6]+0.25*nuVtSqSum[0]*fr[3]+0.25*nuVtSqSum[0]*fc[3]; 
  Gdiff2_r[3] = 0.2445699350390395*nuVtSqSum[1]*fr[36]+0.2445699350390395*nuVtSqSum[1]*fc[36]-0.3518228202874282*nuVtSqSum[2]*fr[35]+0.3518228202874282*nuVtSqSum[2]*fc[35]+0.2445699350390395*nuVtSqSum[0]*fr[26]+0.2445699350390395*nuVtSqSum[0]*fc[26]+0.25*nuVtSqSum[2]*fr[25]+0.25*nuVtSqSum[2]*fc[25]-0.3518228202874282*nuVtSqSum[1]*fr[16]+0.3518228202874282*nuVtSqSum[1]*fc[16]-0.3518228202874282*nuVtSqSum[0]*fr[9]+0.3518228202874282*nuVtSqSum[0]*fc[9]+0.25*nuVtSqSum[1]*fr[8]+0.25*nuVtSqSum[1]*fc[8]+0.25*nuVtSqSum[0]*fr[4]+0.25*nuVtSqSum[0]*fc[4]; 
  Gdiff2_r[4] = 0.21875*nuVtSqSum[2]*fr[33]+0.2445699350390395*nuVtSqSum[0]*fr[33]+0.21875*nuVtSqSum[2]*fc[33]+0.2445699350390395*nuVtSqSum[0]*fc[33]-0.3146798968793526*nuVtSqSum[1]*fr[32]+0.3146798968793526*nuVtSqSum[1]*fc[32]+0.2445699350390395*nuVtSqSum[1]*fr[22]+0.2445699350390395*nuVtSqSum[1]*fc[22]+0.223606797749979*nuVtSqSum[1]*fr[21]+0.223606797749979*nuVtSqSum[1]*fc[21]-0.3146798968793526*nuVtSqSum[2]*fr[15]-0.3518228202874282*nuVtSqSum[0]*fr[15]+0.3146798968793526*nuVtSqSum[2]*fc[15]+0.3518228202874282*nuVtSqSum[0]*fc[15]-0.3518228202874282*nuVtSqSum[1]*fr[7]+0.3518228202874282*nuVtSqSum[1]*fc[7]+0.223606797749979*nuVtSqSum[2]*fr[6]+0.25*nuVtSqSum[0]*fr[6]+0.223606797749979*nuVtSqSum[2]*fc[6]+0.25*nuVtSqSum[0]*fc[6]+0.25*nuVtSqSum[1]*fr[3]+0.25*nuVtSqSum[1]*fc[3]; 
  Gdiff2_r[5] = 0.21875*nuVtSqSum[2]*fr[36]+0.2445699350390395*nuVtSqSum[0]*fr[36]+0.21875*nuVtSqSum[2]*fc[36]+0.2445699350390395*nuVtSqSum[0]*fc[36]-0.3146798968793526*nuVtSqSum[1]*fr[35]+0.3146798968793526*nuVtSqSum[1]*fc[35]+0.2445699350390395*nuVtSqSum[1]*fr[26]+0.2445699350390395*nuVtSqSum[1]*fc[26]+0.223606797749979*nuVtSqSum[1]*fr[25]+0.223606797749979*nuVtSqSum[1]*fc[25]-0.3146798968793526*nuVtSqSum[2]*fr[16]-0.3518228202874282*nuVtSqSum[0]*fr[16]+0.3146798968793526*nuVtSqSum[2]*fc[16]+0.3518228202874282*nuVtSqSum[0]*fc[16]-0.3518228202874282*nuVtSqSum[1]*fr[9]+0.3518228202874282*nuVtSqSum[1]*fc[9]+0.223606797749979*nuVtSqSum[2]*fr[8]+0.25*nuVtSqSum[0]*fr[8]+0.223606797749979*nuVtSqSum[2]*fc[8]+0.25*nuVtSqSum[0]*fc[8]+0.25*nuVtSqSum[1]*fr[4]+0.25*nuVtSqSum[1]*fc[4]; 
  Gdiff2_r[6] = 0.2445699350390395*nuVtSqSum[1]*fr[45]+0.2445699350390395*nuVtSqSum[1]*fc[45]-0.3518228202874282*nuVtSqSum[2]*fr[44]+0.3518228202874282*nuVtSqSum[2]*fc[44]+0.2445699350390395*nuVtSqSum[0]*fr[38]+0.2445699350390395*nuVtSqSum[0]*fc[38]+0.25*nuVtSqSum[2]*fr[37]+0.25*nuVtSqSum[2]*fc[37]-0.3518228202874282*nuVtSqSum[1]*fr[31]+0.3518228202874282*nuVtSqSum[1]*fc[31]-0.3518228202874282*nuVtSqSum[0]*fr[18]+0.3518228202874282*nuVtSqSum[0]*fc[18]+0.25*nuVtSqSum[1]*fr[17]+0.25*nuVtSqSum[1]*fc[17]+0.25*nuVtSqSum[0]*fr[10]+0.25*nuVtSqSum[0]*fc[10]; 
  Gdiff2_r[7] = 0.21875*nuVtSqSum[1]*fr[20]+0.21875*nuVtSqSum[1]*fc[20]-0.2247713549138233*nuVtSqSum[2]*fr[19]-0.3518228202874282*nuVtSqSum[0]*fr[19]+0.2247713549138233*nuVtSqSum[2]*fc[19]+0.3518228202874282*nuVtSqSum[0]*fc[19]+0.2445699350390395*nuVtSqSum[2]*fr[12]+0.2445699350390395*nuVtSqSum[2]*fc[12]+0.159719141249985*nuVtSqSum[2]*fr[11]+0.25*nuVtSqSum[0]*fr[11]+0.159719141249985*nuVtSqSum[2]*fc[11]+0.25*nuVtSqSum[0]*fc[11]-0.3146798968793526*nuVtSqSum[1]*fr[5]+0.3146798968793526*nuVtSqSum[1]*fc[5]-0.3518228202874282*fr[2]*nuVtSqSum[2]+0.3518228202874282*fc[2]*nuVtSqSum[2]+0.25*fr[0]*nuVtSqSum[2]+0.25*fc[0]*nuVtSqSum[2]+0.223606797749979*fr[1]*nuVtSqSum[1]+0.223606797749979*fc[1]*nuVtSqSum[1]; 
  Gdiff2_r[8] = (-0.3518228202874282*nuVtSqSum[1]*fr[34])+0.3518228202874282*nuVtSqSum[1]*fc[34]-0.3518228202874282*nuVtSqSum[0]*fr[24]+0.3518228202874282*nuVtSqSum[0]*fc[24]+0.25*nuVtSqSum[1]*fr[23]+0.25*nuVtSqSum[1]*fc[23]+0.25*nuVtSqSum[0]*fr[13]+0.25*nuVtSqSum[0]*fc[13]; 
  Gdiff2_r[9] = (-0.3518228202874282*nuVtSqSum[1]*fr[41])+0.3518228202874282*nuVtSqSum[1]*fc[41]-0.3518228202874282*nuVtSqSum[0]*fr[29]+0.3518228202874282*nuVtSqSum[0]*fc[29]+0.25*nuVtSqSum[1]*fr[28]+0.25*nuVtSqSum[1]*fc[28]+0.25*nuVtSqSum[0]*fr[14]+0.25*nuVtSqSum[0]*fc[14]; 
  Gdiff2_r[10] = 0.21875*nuVtSqSum[2]*fr[45]+0.2445699350390395*nuVtSqSum[0]*fr[45]+0.21875*nuVtSqSum[2]*fc[45]+0.2445699350390395*nuVtSqSum[0]*fc[45]-0.3146798968793526*nuVtSqSum[1]*fr[44]+0.3146798968793526*nuVtSqSum[1]*fc[44]+0.2445699350390395*nuVtSqSum[1]*fr[38]+0.2445699350390395*nuVtSqSum[1]*fc[38]+0.223606797749979*nuVtSqSum[1]*fr[37]+0.223606797749979*nuVtSqSum[1]*fc[37]-0.3146798968793526*nuVtSqSum[2]*fr[31]-0.3518228202874282*nuVtSqSum[0]*fr[31]+0.3146798968793526*nuVtSqSum[2]*fc[31]+0.3518228202874282*nuVtSqSum[0]*fc[31]-0.3518228202874282*nuVtSqSum[1]*fr[18]+0.3518228202874282*nuVtSqSum[1]*fc[18]+0.223606797749979*nuVtSqSum[2]*fr[17]+0.25*nuVtSqSum[0]*fr[17]+0.223606797749979*nuVtSqSum[2]*fc[17]+0.25*nuVtSqSum[0]*fc[17]+0.25*nuVtSqSum[1]*fr[10]+0.25*nuVtSqSum[1]*fc[10]; 
  Gdiff2_r[11] = 0.21875*nuVtSqSum[1]*fr[33]+0.21875*nuVtSqSum[1]*fc[33]-0.2247713549138233*nuVtSqSum[2]*fr[32]-0.3518228202874282*nuVtSqSum[0]*fr[32]+0.2247713549138233*nuVtSqSum[2]*fc[32]+0.3518228202874282*nuVtSqSum[0]*fc[32]+0.2445699350390395*nuVtSqSum[2]*fr[22]+0.2445699350390395*nuVtSqSum[2]*fc[22]+0.159719141249985*nuVtSqSum[2]*fr[21]+0.25*nuVtSqSum[0]*fr[21]+0.159719141249985*nuVtSqSum[2]*fc[21]+0.25*nuVtSqSum[0]*fc[21]-0.3146798968793526*nuVtSqSum[1]*fr[15]+0.3146798968793526*nuVtSqSum[1]*fc[15]-0.3518228202874282*nuVtSqSum[2]*fr[7]+0.3518228202874282*nuVtSqSum[2]*fc[7]+0.223606797749979*nuVtSqSum[1]*fr[6]+0.223606797749979*nuVtSqSum[1]*fc[6]+0.25*nuVtSqSum[2]*fr[3]+0.25*nuVtSqSum[2]*fc[3]; 
  Gdiff2_r[12] = (-0.3146798968793526*nuVtSqSum[2]*fr[34])-0.3518228202874282*nuVtSqSum[0]*fr[34]+0.3146798968793526*nuVtSqSum[2]*fc[34]+0.3518228202874282*nuVtSqSum[0]*fc[34]-0.3518228202874282*nuVtSqSum[1]*fr[24]+0.3518228202874282*nuVtSqSum[1]*fc[24]+0.223606797749979*nuVtSqSum[2]*fr[23]+0.25*nuVtSqSum[0]*fr[23]+0.223606797749979*nuVtSqSum[2]*fc[23]+0.25*nuVtSqSum[0]*fc[23]+0.25*nuVtSqSum[1]*fr[13]+0.25*nuVtSqSum[1]*fc[13]; 
  Gdiff2_r[13] = 0.21875*nuVtSqSum[1]*fr[36]+0.21875*nuVtSqSum[1]*fc[36]-0.2247713549138233*nuVtSqSum[2]*fr[35]-0.3518228202874282*nuVtSqSum[0]*fr[35]+0.2247713549138233*nuVtSqSum[2]*fc[35]+0.3518228202874282*nuVtSqSum[0]*fc[35]+0.2445699350390395*nuVtSqSum[2]*fr[26]+0.2445699350390395*nuVtSqSum[2]*fc[26]+0.159719141249985*nuVtSqSum[2]*fr[25]+0.25*nuVtSqSum[0]*fr[25]+0.159719141249985*nuVtSqSum[2]*fc[25]+0.25*nuVtSqSum[0]*fc[25]-0.3146798968793526*nuVtSqSum[1]*fr[16]+0.3146798968793526*nuVtSqSum[1]*fc[16]-0.3518228202874282*nuVtSqSum[2]*fr[9]+0.3518228202874282*nuVtSqSum[2]*fc[9]+0.223606797749979*nuVtSqSum[1]*fr[8]+0.223606797749979*nuVtSqSum[1]*fc[8]+0.25*nuVtSqSum[2]*fr[4]+0.25*nuVtSqSum[2]*fc[4]; 
  Gdiff2_r[14] = (-0.3518228202874282*nuVtSqSum[1]*fr[46])+0.3518228202874282*nuVtSqSum[1]*fc[46]-0.3518228202874282*nuVtSqSum[0]*fr[40]+0.3518228202874282*nuVtSqSum[0]*fc[40]+0.25*nuVtSqSum[1]*fr[39]+0.25*nuVtSqSum[1]*fc[39]+0.25*nuVtSqSum[0]*fr[27]+0.25*nuVtSqSum[0]*fc[27]; 
  Gdiff2_r[15] = (-0.3146798968793526*nuVtSqSum[2]*fr[41])-0.3518228202874282*nuVtSqSum[0]*fr[41]+0.3146798968793526*nuVtSqSum[2]*fc[41]+0.3518228202874282*nuVtSqSum[0]*fc[41]-0.3518228202874282*nuVtSqSum[1]*fr[29]+0.3518228202874282*nuVtSqSum[1]*fc[29]+0.223606797749979*nuVtSqSum[2]*fr[28]+0.25*nuVtSqSum[0]*fr[28]+0.223606797749979*nuVtSqSum[2]*fc[28]+0.25*nuVtSqSum[0]*fc[28]+0.25*nuVtSqSum[1]*fr[14]+0.25*nuVtSqSum[1]*fc[14]; 
  Gdiff2_r[16] = (-0.3518228202874282*nuVtSqSum[1]*fr[47])+0.3518228202874282*nuVtSqSum[1]*fc[47]-0.3518228202874282*nuVtSqSum[0]*fr[43]+0.3518228202874282*nuVtSqSum[0]*fc[43]+0.25*nuVtSqSum[1]*fr[42]+0.25*nuVtSqSum[1]*fc[42]+0.25*nuVtSqSum[0]*fr[30]+0.25*nuVtSqSum[0]*fc[30]; 
  Gdiff2_r[17] = 0.21875*nuVtSqSum[1]*fr[45]+0.21875*nuVtSqSum[1]*fc[45]-0.2247713549138233*nuVtSqSum[2]*fr[44]-0.3518228202874282*nuVtSqSum[0]*fr[44]+0.2247713549138233*nuVtSqSum[2]*fc[44]+0.3518228202874282*nuVtSqSum[0]*fc[44]+0.2445699350390395*nuVtSqSum[2]*fr[38]+0.2445699350390395*nuVtSqSum[2]*fc[38]+0.159719141249985*nuVtSqSum[2]*fr[37]+0.25*nuVtSqSum[0]*fr[37]+0.159719141249985*nuVtSqSum[2]*fc[37]+0.25*nuVtSqSum[0]*fc[37]-0.3146798968793526*nuVtSqSum[1]*fr[31]+0.3146798968793526*nuVtSqSum[1]*fc[31]-0.3518228202874282*nuVtSqSum[2]*fr[18]+0.3518228202874282*nuVtSqSum[2]*fc[18]+0.223606797749979*nuVtSqSum[1]*fr[17]+0.223606797749979*nuVtSqSum[1]*fc[17]+0.25*nuVtSqSum[2]*fr[10]+0.25*nuVtSqSum[2]*fc[10]; 
  Gdiff2_r[18] = (-0.3146798968793526*nuVtSqSum[2]*fr[46])-0.3518228202874282*nuVtSqSum[0]*fr[46]+0.3146798968793526*nuVtSqSum[2]*fc[46]+0.3518228202874282*nuVtSqSum[0]*fc[46]-0.3518228202874282*nuVtSqSum[1]*fr[40]+0.3518228202874282*nuVtSqSum[1]*fc[40]+0.223606797749979*nuVtSqSum[2]*fr[39]+0.25*nuVtSqSum[0]*fr[39]+0.223606797749979*nuVtSqSum[2]*fc[39]+0.25*nuVtSqSum[0]*fc[39]+0.25*nuVtSqSum[1]*fr[27]+0.25*nuVtSqSum[1]*fc[27]; 
  Gdiff2_r[19] = (-0.3146798968793526*nuVtSqSum[2]*fr[47])-0.3518228202874282*nuVtSqSum[0]*fr[47]+0.3146798968793526*nuVtSqSum[2]*fc[47]+0.3518228202874282*nuVtSqSum[0]*fc[47]-0.3518228202874282*nuVtSqSum[1]*fr[43]+0.3518228202874282*nuVtSqSum[1]*fc[43]+0.223606797749979*nuVtSqSum[2]*fr[42]+0.25*nuVtSqSum[0]*fr[42]+0.223606797749979*nuVtSqSum[2]*fc[42]+0.25*nuVtSqSum[0]*fc[42]+0.25*nuVtSqSum[1]*fr[30]+0.25*nuVtSqSum[1]*fc[30]; 

  Gdiff_l[0] = (-0.6708203932499369*nuVtSqSum[1]*fl[20])+0.6708203932499369*nuVtSqSum[1]*fc[20]-1.190784930203603*nuVtSqSum[2]*fl[19]-1.190784930203603*nuVtSqSum[2]*fc[19]-0.6708203932499369*nuVtSqSum[0]*fl[12]+0.6708203932499369*nuVtSqSum[0]*fc[12]-0.9375*nuVtSqSum[2]*fl[11]+0.9375*nuVtSqSum[2]*fc[11]-1.190784930203603*nuVtSqSum[1]*fl[5]-1.190784930203603*nuVtSqSum[1]*fc[5]-1.190784930203603*nuVtSqSum[0]*fl[2]-1.190784930203603*nuVtSqSum[0]*fc[2]-0.9375*fl[1]*nuVtSqSum[1]+0.9375*fc[1]*nuVtSqSum[1]-0.9375*fl[0]*nuVtSqSum[0]+0.9375*fc[0]*nuVtSqSum[0]; 
  Gdiff_l[1] = (-0.5999999999999999*nuVtSqSum[2]*fl[20])-0.6708203932499369*nuVtSqSum[0]*fl[20]+0.5999999999999999*nuVtSqSum[2]*fc[20]+0.6708203932499369*nuVtSqSum[0]*fc[20]-1.06507042020704*nuVtSqSum[1]*fl[19]-1.06507042020704*nuVtSqSum[1]*fc[19]-0.6708203932499369*nuVtSqSum[1]*fl[12]+0.6708203932499369*nuVtSqSum[1]*fc[12]-0.8385254915624212*nuVtSqSum[1]*fl[11]+0.8385254915624212*nuVtSqSum[1]*fc[11]-1.06507042020704*nuVtSqSum[2]*fl[5]-1.190784930203603*nuVtSqSum[0]*fl[5]-1.06507042020704*nuVtSqSum[2]*fc[5]-1.190784930203603*nuVtSqSum[0]*fc[5]-0.8385254915624212*fl[1]*nuVtSqSum[2]+0.8385254915624212*fc[1]*nuVtSqSum[2]-1.190784930203603*nuVtSqSum[1]*fl[2]-1.190784930203603*nuVtSqSum[1]*fc[2]-0.9375*fl[0]*nuVtSqSum[1]+0.9375*fc[0]*nuVtSqSum[1]-0.9375*nuVtSqSum[0]*fl[1]+0.9375*nuVtSqSum[0]*fc[1]; 
  Gdiff_l[2] = (-0.6708203932499369*nuVtSqSum[1]*fl[33])+0.6708203932499369*nuVtSqSum[1]*fc[33]-1.190784930203603*nuVtSqSum[2]*fl[32]-1.190784930203603*nuVtSqSum[2]*fc[32]-0.6708203932499369*nuVtSqSum[0]*fl[22]+0.6708203932499369*nuVtSqSum[0]*fc[22]-0.9375000000000001*nuVtSqSum[2]*fl[21]+0.9375000000000001*nuVtSqSum[2]*fc[21]-1.190784930203603*nuVtSqSum[1]*fl[15]-1.190784930203603*nuVtSqSum[1]*fc[15]-1.190784930203603*nuVtSqSum[0]*fl[7]-1.190784930203603*nuVtSqSum[0]*fc[7]-0.9375*nuVtSqSum[1]*fl[6]+0.9375*nuVtSqSum[1]*fc[6]-0.9375*nuVtSqSum[0]*fl[3]+0.9375*nuVtSqSum[0]*fc[3]; 
  Gdiff_l[3] = (-0.6708203932499369*nuVtSqSum[1]*fl[36])+0.6708203932499369*nuVtSqSum[1]*fc[36]-1.190784930203603*nuVtSqSum[2]*fl[35]-1.190784930203603*nuVtSqSum[2]*fc[35]-0.6708203932499369*nuVtSqSum[0]*fl[26]+0.6708203932499369*nuVtSqSum[0]*fc[26]-0.9375000000000001*nuVtSqSum[2]*fl[25]+0.9375000000000001*nuVtSqSum[2]*fc[25]-1.190784930203603*nuVtSqSum[1]*fl[16]-1.190784930203603*nuVtSqSum[1]*fc[16]-1.190784930203603*nuVtSqSum[0]*fl[9]-1.190784930203603*nuVtSqSum[0]*fc[9]-0.9375*nuVtSqSum[1]*fl[8]+0.9375*nuVtSqSum[1]*fc[8]-0.9375*nuVtSqSum[0]*fl[4]+0.9375*nuVtSqSum[0]*fc[4]; 
  Gdiff_l[4] = (-0.6*nuVtSqSum[2]*fl[33])-0.6708203932499369*nuVtSqSum[0]*fl[33]+0.6*nuVtSqSum[2]*fc[33]+0.6708203932499369*nuVtSqSum[0]*fc[33]-1.06507042020704*nuVtSqSum[1]*fl[32]-1.06507042020704*nuVtSqSum[1]*fc[32]-0.6708203932499369*nuVtSqSum[1]*fl[22]+0.6708203932499369*nuVtSqSum[1]*fc[22]-0.8385254915624211*nuVtSqSum[1]*fl[21]+0.8385254915624211*nuVtSqSum[1]*fc[21]-1.06507042020704*nuVtSqSum[2]*fl[15]-1.190784930203603*nuVtSqSum[0]*fl[15]-1.06507042020704*nuVtSqSum[2]*fc[15]-1.190784930203603*nuVtSqSum[0]*fc[15]-1.190784930203603*nuVtSqSum[1]*fl[7]-1.190784930203603*nuVtSqSum[1]*fc[7]-0.8385254915624212*nuVtSqSum[2]*fl[6]-0.9375*nuVtSqSum[0]*fl[6]+0.8385254915624212*nuVtSqSum[2]*fc[6]+0.9375*nuVtSqSum[0]*fc[6]-0.9375*nuVtSqSum[1]*fl[3]+0.9375*nuVtSqSum[1]*fc[3]; 
  Gdiff_l[5] = (-0.6*nuVtSqSum[2]*fl[36])-0.6708203932499369*nuVtSqSum[0]*fl[36]+0.6*nuVtSqSum[2]*fc[36]+0.6708203932499369*nuVtSqSum[0]*fc[36]-1.06507042020704*nuVtSqSum[1]*fl[35]-1.06507042020704*nuVtSqSum[1]*fc[35]-0.6708203932499369*nuVtSqSum[1]*fl[26]+0.6708203932499369*nuVtSqSum[1]*fc[26]-0.8385254915624211*nuVtSqSum[1]*fl[25]+0.8385254915624211*nuVtSqSum[1]*fc[25]-1.06507042020704*nuVtSqSum[2]*fl[16]-1.190784930203603*nuVtSqSum[0]*fl[16]-1.06507042020704*nuVtSqSum[2]*fc[16]-1.190784930203603*nuVtSqSum[0]*fc[16]-1.190784930203603*nuVtSqSum[1]*fl[9]-1.190784930203603*nuVtSqSum[1]*fc[9]-0.8385254915624212*nuVtSqSum[2]*fl[8]-0.9375*nuVtSqSum[0]*fl[8]+0.8385254915624212*nuVtSqSum[2]*fc[8]+0.9375*nuVtSqSum[0]*fc[8]-0.9375*nuVtSqSum[1]*fl[4]+0.9375*nuVtSqSum[1]*fc[4]; 
  Gdiff_l[6] = (-0.6708203932499369*nuVtSqSum[1]*fl[45])+0.6708203932499369*nuVtSqSum[1]*fc[45]-1.190784930203603*nuVtSqSum[2]*fl[44]-1.190784930203603*nuVtSqSum[2]*fc[44]-0.6708203932499369*nuVtSqSum[0]*fl[38]+0.6708203932499369*nuVtSqSum[0]*fc[38]-0.9375*nuVtSqSum[2]*fl[37]+0.9375*nuVtSqSum[2]*fc[37]-1.190784930203603*nuVtSqSum[1]*fl[31]-1.190784930203603*nuVtSqSum[1]*fc[31]-1.190784930203603*nuVtSqSum[0]*fl[18]-1.190784930203603*nuVtSqSum[0]*fc[18]-0.9375*nuVtSqSum[1]*fl[17]+0.9375*nuVtSqSum[1]*fc[17]-0.9375*nuVtSqSum[0]*fl[10]+0.9375*nuVtSqSum[0]*fc[10]; 
  Gdiff_l[7] = (-0.5999999999999999*nuVtSqSum[1]*fl[20])+0.5999999999999999*nuVtSqSum[1]*fc[20]-0.7607645858621712*nuVtSqSum[2]*fl[19]-1.190784930203603*nuVtSqSum[0]*fl[19]-0.7607645858621712*nuVtSqSum[2]*fc[19]-1.190784930203603*nuVtSqSum[0]*fc[19]-0.6708203932499369*nuVtSqSum[2]*fl[12]+0.6708203932499369*nuVtSqSum[2]*fc[12]-0.5989467796874438*nuVtSqSum[2]*fl[11]-0.9375*nuVtSqSum[0]*fl[11]+0.5989467796874438*nuVtSqSum[2]*fc[11]+0.9375*nuVtSqSum[0]*fc[11]-1.06507042020704*nuVtSqSum[1]*fl[5]-1.06507042020704*nuVtSqSum[1]*fc[5]-1.190784930203603*fl[2]*nuVtSqSum[2]-1.190784930203603*fc[2]*nuVtSqSum[2]-0.9375*fl[0]*nuVtSqSum[2]+0.9375*fc[0]*nuVtSqSum[2]-0.8385254915624212*fl[1]*nuVtSqSum[1]+0.8385254915624212*fc[1]*nuVtSqSum[1]; 
  Gdiff_l[8] = (-1.190784930203603*nuVtSqSum[1]*fl[34])-1.190784930203603*nuVtSqSum[1]*fc[34]-1.190784930203603*nuVtSqSum[0]*fl[24]-1.190784930203603*nuVtSqSum[0]*fc[24]-0.9375000000000001*nuVtSqSum[1]*fl[23]+0.9375000000000001*nuVtSqSum[1]*fc[23]-0.9375*nuVtSqSum[0]*fl[13]+0.9375*nuVtSqSum[0]*fc[13]; 
  Gdiff_l[9] = (-1.190784930203603*nuVtSqSum[1]*fl[41])-1.190784930203603*nuVtSqSum[1]*fc[41]-1.190784930203603*nuVtSqSum[0]*fl[29]-1.190784930203603*nuVtSqSum[0]*fc[29]-0.9375000000000001*nuVtSqSum[1]*fl[28]+0.9375000000000001*nuVtSqSum[1]*fc[28]-0.9375*nuVtSqSum[0]*fl[14]+0.9375*nuVtSqSum[0]*fc[14]; 
  Gdiff_l[10] = (-0.5999999999999999*nuVtSqSum[2]*fl[45])-0.6708203932499369*nuVtSqSum[0]*fl[45]+0.5999999999999999*nuVtSqSum[2]*fc[45]+0.6708203932499369*nuVtSqSum[0]*fc[45]-1.06507042020704*nuVtSqSum[1]*fl[44]-1.06507042020704*nuVtSqSum[1]*fc[44]-0.6708203932499369*nuVtSqSum[1]*fl[38]+0.6708203932499369*nuVtSqSum[1]*fc[38]-0.8385254915624212*nuVtSqSum[1]*fl[37]+0.8385254915624212*nuVtSqSum[1]*fc[37]-1.06507042020704*nuVtSqSum[2]*fl[31]-1.190784930203603*nuVtSqSum[0]*fl[31]-1.06507042020704*nuVtSqSum[2]*fc[31]-1.190784930203603*nuVtSqSum[0]*fc[31]-1.190784930203603*nuVtSqSum[1]*fl[18]-1.190784930203603*nuVtSqSum[1]*fc[18]-0.8385254915624212*nuVtSqSum[2]*fl[17]-0.9375*nuVtSqSum[0]*fl[17]+0.8385254915624212*nuVtSqSum[2]*fc[17]+0.9375*nuVtSqSum[0]*fc[17]-0.9375*nuVtSqSum[1]*fl[10]+0.9375*nuVtSqSum[1]*fc[10]; 
  Gdiff_l[11] = (-0.5999999999999999*nuVtSqSum[1]*fl[33])+0.5999999999999999*nuVtSqSum[1]*fc[33]-0.7607645858621712*nuVtSqSum[2]*fl[32]-1.190784930203603*nuVtSqSum[0]*fl[32]-0.7607645858621712*nuVtSqSum[2]*fc[32]-1.190784930203603*nuVtSqSum[0]*fc[32]-0.6708203932499369*nuVtSqSum[2]*fl[22]+0.6708203932499369*nuVtSqSum[2]*fc[22]-0.5989467796874438*nuVtSqSum[2]*fl[21]-0.9375*nuVtSqSum[0]*fl[21]+0.5989467796874438*nuVtSqSum[2]*fc[21]+0.9375*nuVtSqSum[0]*fc[21]-1.06507042020704*nuVtSqSum[1]*fl[15]-1.06507042020704*nuVtSqSum[1]*fc[15]-1.190784930203603*nuVtSqSum[2]*fl[7]-1.190784930203603*nuVtSqSum[2]*fc[7]-0.8385254915624211*nuVtSqSum[1]*fl[6]+0.8385254915624211*nuVtSqSum[1]*fc[6]-0.9375000000000001*nuVtSqSum[2]*fl[3]+0.9375000000000001*nuVtSqSum[2]*fc[3]; 
  Gdiff_l[12] = (-1.06507042020704*nuVtSqSum[2]*fl[34])-1.190784930203603*nuVtSqSum[0]*fl[34]-1.06507042020704*nuVtSqSum[2]*fc[34]-1.190784930203603*nuVtSqSum[0]*fc[34]-1.190784930203603*nuVtSqSum[1]*fl[24]-1.190784930203603*nuVtSqSum[1]*fc[24]-0.8385254915624212*nuVtSqSum[2]*fl[23]-0.9375*nuVtSqSum[0]*fl[23]+0.8385254915624212*nuVtSqSum[2]*fc[23]+0.9375*nuVtSqSum[0]*fc[23]-0.9375000000000001*nuVtSqSum[1]*fl[13]+0.9375000000000001*nuVtSqSum[1]*fc[13]; 
  Gdiff_l[13] = (-0.5999999999999999*nuVtSqSum[1]*fl[36])+0.5999999999999999*nuVtSqSum[1]*fc[36]-0.7607645858621712*nuVtSqSum[2]*fl[35]-1.190784930203603*nuVtSqSum[0]*fl[35]-0.7607645858621712*nuVtSqSum[2]*fc[35]-1.190784930203603*nuVtSqSum[0]*fc[35]-0.6708203932499369*nuVtSqSum[2]*fl[26]+0.6708203932499369*nuVtSqSum[2]*fc[26]-0.5989467796874438*nuVtSqSum[2]*fl[25]-0.9375*nuVtSqSum[0]*fl[25]+0.5989467796874438*nuVtSqSum[2]*fc[25]+0.9375*nuVtSqSum[0]*fc[25]-1.06507042020704*nuVtSqSum[1]*fl[16]-1.06507042020704*nuVtSqSum[1]*fc[16]-1.190784930203603*nuVtSqSum[2]*fl[9]-1.190784930203603*nuVtSqSum[2]*fc[9]-0.8385254915624211*nuVtSqSum[1]*fl[8]+0.8385254915624211*nuVtSqSum[1]*fc[8]-0.9375000000000001*nuVtSqSum[2]*fl[4]+0.9375000000000001*nuVtSqSum[2]*fc[4]; 
  Gdiff_l[14] = (-1.190784930203603*nuVtSqSum[1]*fl[46])-1.190784930203603*nuVtSqSum[1]*fc[46]-1.190784930203603*nuVtSqSum[0]*fl[40]-1.190784930203603*nuVtSqSum[0]*fc[40]-0.9375000000000001*nuVtSqSum[1]*fl[39]+0.9375000000000001*nuVtSqSum[1]*fc[39]-0.9375*nuVtSqSum[0]*fl[27]+0.9375*nuVtSqSum[0]*fc[27]; 
  Gdiff_l[15] = (-1.06507042020704*nuVtSqSum[2]*fl[41])-1.190784930203603*nuVtSqSum[0]*fl[41]-1.06507042020704*nuVtSqSum[2]*fc[41]-1.190784930203603*nuVtSqSum[0]*fc[41]-1.190784930203603*nuVtSqSum[1]*fl[29]-1.190784930203603*nuVtSqSum[1]*fc[29]-0.8385254915624212*nuVtSqSum[2]*fl[28]-0.9375*nuVtSqSum[0]*fl[28]+0.8385254915624212*nuVtSqSum[2]*fc[28]+0.9375*nuVtSqSum[0]*fc[28]-0.9375000000000001*nuVtSqSum[1]*fl[14]+0.9375000000000001*nuVtSqSum[1]*fc[14]; 
  Gdiff_l[16] = (-1.190784930203603*nuVtSqSum[1]*fl[47])-1.190784930203603*nuVtSqSum[1]*fc[47]-1.190784930203603*nuVtSqSum[0]*fl[43]-1.190784930203603*nuVtSqSum[0]*fc[43]-0.9375000000000001*nuVtSqSum[1]*fl[42]+0.9375000000000001*nuVtSqSum[1]*fc[42]-0.9375*nuVtSqSum[0]*fl[30]+0.9375*nuVtSqSum[0]*fc[30]; 
  Gdiff_l[17] = (-0.5999999999999999*nuVtSqSum[1]*fl[45])+0.5999999999999999*nuVtSqSum[1]*fc[45]-0.7607645858621712*nuVtSqSum[2]*fl[44]-1.190784930203603*nuVtSqSum[0]*fl[44]-0.7607645858621712*nuVtSqSum[2]*fc[44]-1.190784930203603*nuVtSqSum[0]*fc[44]-0.6708203932499369*nuVtSqSum[2]*fl[38]+0.6708203932499369*nuVtSqSum[2]*fc[38]-0.5989467796874438*nuVtSqSum[2]*fl[37]-0.9375*nuVtSqSum[0]*fl[37]+0.5989467796874438*nuVtSqSum[2]*fc[37]+0.9375*nuVtSqSum[0]*fc[37]-1.06507042020704*nuVtSqSum[1]*fl[31]-1.06507042020704*nuVtSqSum[1]*fc[31]-1.190784930203603*nuVtSqSum[2]*fl[18]-1.190784930203603*nuVtSqSum[2]*fc[18]-0.8385254915624212*nuVtSqSum[1]*fl[17]+0.8385254915624212*nuVtSqSum[1]*fc[17]-0.9375*nuVtSqSum[2]*fl[10]+0.9375*nuVtSqSum[2]*fc[10]; 
  Gdiff_l[18] = (-1.06507042020704*nuVtSqSum[2]*fl[46])-1.190784930203603*nuVtSqSum[0]*fl[46]-1.06507042020704*nuVtSqSum[2]*fc[46]-1.190784930203603*nuVtSqSum[0]*fc[46]-1.190784930203603*nuVtSqSum[1]*fl[40]-1.190784930203603*nuVtSqSum[1]*fc[40]-0.8385254915624212*nuVtSqSum[2]*fl[39]-0.9375*nuVtSqSum[0]*fl[39]+0.8385254915624212*nuVtSqSum[2]*fc[39]+0.9375*nuVtSqSum[0]*fc[39]-0.9375000000000001*nuVtSqSum[1]*fl[27]+0.9375000000000001*nuVtSqSum[1]*fc[27]; 
  Gdiff_l[19] = (-1.06507042020704*nuVtSqSum[2]*fl[47])-1.190784930203603*nuVtSqSum[0]*fl[47]-1.06507042020704*nuVtSqSum[2]*fc[47]-1.190784930203603*nuVtSqSum[0]*fc[47]-1.190784930203603*nuVtSqSum[1]*fl[43]-1.190784930203603*nuVtSqSum[1]*fc[43]-0.8385254915624212*nuVtSqSum[2]*fl[42]-0.9375*nuVtSqSum[0]*fl[42]+0.8385254915624212*nuVtSqSum[2]*fc[42]+0.9375*nuVtSqSum[0]*fc[42]-0.9375000000000001*nuVtSqSum[1]*fl[30]+0.9375000000000001*nuVtSqSum[1]*fc[30]; 

  Gdiff_r[0] = 0.6708203932499369*nuVtSqSum[1]*fr[20]-0.6708203932499369*nuVtSqSum[1]*fc[20]-1.190784930203603*nuVtSqSum[2]*fr[19]-1.190784930203603*nuVtSqSum[2]*fc[19]+0.6708203932499369*nuVtSqSum[0]*fr[12]-0.6708203932499369*nuVtSqSum[0]*fc[12]+0.9375*nuVtSqSum[2]*fr[11]-0.9375*nuVtSqSum[2]*fc[11]-1.190784930203603*nuVtSqSum[1]*fr[5]-1.190784930203603*nuVtSqSum[1]*fc[5]-1.190784930203603*nuVtSqSum[0]*fr[2]-1.190784930203603*nuVtSqSum[0]*fc[2]+0.9375*fr[1]*nuVtSqSum[1]-0.9375*fc[1]*nuVtSqSum[1]+0.9375*fr[0]*nuVtSqSum[0]-0.9375*fc[0]*nuVtSqSum[0]; 
  Gdiff_r[1] = 0.5999999999999999*nuVtSqSum[2]*fr[20]+0.6708203932499369*nuVtSqSum[0]*fr[20]-0.5999999999999999*nuVtSqSum[2]*fc[20]-0.6708203932499369*nuVtSqSum[0]*fc[20]-1.06507042020704*nuVtSqSum[1]*fr[19]-1.06507042020704*nuVtSqSum[1]*fc[19]+0.6708203932499369*nuVtSqSum[1]*fr[12]-0.6708203932499369*nuVtSqSum[1]*fc[12]+0.8385254915624212*nuVtSqSum[1]*fr[11]-0.8385254915624212*nuVtSqSum[1]*fc[11]-1.06507042020704*nuVtSqSum[2]*fr[5]-1.190784930203603*nuVtSqSum[0]*fr[5]-1.06507042020704*nuVtSqSum[2]*fc[5]-1.190784930203603*nuVtSqSum[0]*fc[5]+0.8385254915624212*fr[1]*nuVtSqSum[2]-0.8385254915624212*fc[1]*nuVtSqSum[2]-1.190784930203603*nuVtSqSum[1]*fr[2]-1.190784930203603*nuVtSqSum[1]*fc[2]+0.9375*fr[0]*nuVtSqSum[1]-0.9375*fc[0]*nuVtSqSum[1]+0.9375*nuVtSqSum[0]*fr[1]-0.9375*nuVtSqSum[0]*fc[1]; 
  Gdiff_r[2] = 0.6708203932499369*nuVtSqSum[1]*fr[33]-0.6708203932499369*nuVtSqSum[1]*fc[33]-1.190784930203603*nuVtSqSum[2]*fr[32]-1.190784930203603*nuVtSqSum[2]*fc[32]+0.6708203932499369*nuVtSqSum[0]*fr[22]-0.6708203932499369*nuVtSqSum[0]*fc[22]+0.9375000000000001*nuVtSqSum[2]*fr[21]-0.9375000000000001*nuVtSqSum[2]*fc[21]-1.190784930203603*nuVtSqSum[1]*fr[15]-1.190784930203603*nuVtSqSum[1]*fc[15]-1.190784930203603*nuVtSqSum[0]*fr[7]-1.190784930203603*nuVtSqSum[0]*fc[7]+0.9375*nuVtSqSum[1]*fr[6]-0.9375*nuVtSqSum[1]*fc[6]+0.9375*nuVtSqSum[0]*fr[3]-0.9375*nuVtSqSum[0]*fc[3]; 
  Gdiff_r[3] = 0.6708203932499369*nuVtSqSum[1]*fr[36]-0.6708203932499369*nuVtSqSum[1]*fc[36]-1.190784930203603*nuVtSqSum[2]*fr[35]-1.190784930203603*nuVtSqSum[2]*fc[35]+0.6708203932499369*nuVtSqSum[0]*fr[26]-0.6708203932499369*nuVtSqSum[0]*fc[26]+0.9375000000000001*nuVtSqSum[2]*fr[25]-0.9375000000000001*nuVtSqSum[2]*fc[25]-1.190784930203603*nuVtSqSum[1]*fr[16]-1.190784930203603*nuVtSqSum[1]*fc[16]-1.190784930203603*nuVtSqSum[0]*fr[9]-1.190784930203603*nuVtSqSum[0]*fc[9]+0.9375*nuVtSqSum[1]*fr[8]-0.9375*nuVtSqSum[1]*fc[8]+0.9375*nuVtSqSum[0]*fr[4]-0.9375*nuVtSqSum[0]*fc[4]; 
  Gdiff_r[4] = 0.6*nuVtSqSum[2]*fr[33]+0.6708203932499369*nuVtSqSum[0]*fr[33]-0.6*nuVtSqSum[2]*fc[33]-0.6708203932499369*nuVtSqSum[0]*fc[33]-1.06507042020704*nuVtSqSum[1]*fr[32]-1.06507042020704*nuVtSqSum[1]*fc[32]+0.6708203932499369*nuVtSqSum[1]*fr[22]-0.6708203932499369*nuVtSqSum[1]*fc[22]+0.8385254915624211*nuVtSqSum[1]*fr[21]-0.8385254915624211*nuVtSqSum[1]*fc[21]-1.06507042020704*nuVtSqSum[2]*fr[15]-1.190784930203603*nuVtSqSum[0]*fr[15]-1.06507042020704*nuVtSqSum[2]*fc[15]-1.190784930203603*nuVtSqSum[0]*fc[15]-1.190784930203603*nuVtSqSum[1]*fr[7]-1.190784930203603*nuVtSqSum[1]*fc[7]+0.8385254915624212*nuVtSqSum[2]*fr[6]+0.9375*nuVtSqSum[0]*fr[6]-0.8385254915624212*nuVtSqSum[2]*fc[6]-0.9375*nuVtSqSum[0]*fc[6]+0.9375*nuVtSqSum[1]*fr[3]-0.9375*nuVtSqSum[1]*fc[3]; 
  Gdiff_r[5] = 0.6*nuVtSqSum[2]*fr[36]+0.6708203932499369*nuVtSqSum[0]*fr[36]-0.6*nuVtSqSum[2]*fc[36]-0.6708203932499369*nuVtSqSum[0]*fc[36]-1.06507042020704*nuVtSqSum[1]*fr[35]-1.06507042020704*nuVtSqSum[1]*fc[35]+0.6708203932499369*nuVtSqSum[1]*fr[26]-0.6708203932499369*nuVtSqSum[1]*fc[26]+0.8385254915624211*nuVtSqSum[1]*fr[25]-0.8385254915624211*nuVtSqSum[1]*fc[25]-1.06507042020704*nuVtSqSum[2]*fr[16]-1.190784930203603*nuVtSqSum[0]*fr[16]-1.06507042020704*nuVtSqSum[2]*fc[16]-1.190784930203603*nuVtSqSum[0]*fc[16]-1.190784930203603*nuVtSqSum[1]*fr[9]-1.190784930203603*nuVtSqSum[1]*fc[9]+0.8385254915624212*nuVtSqSum[2]*fr[8]+0.9375*nuVtSqSum[0]*fr[8]-0.8385254915624212*nuVtSqSum[2]*fc[8]-0.9375*nuVtSqSum[0]*fc[8]+0.9375*nuVtSqSum[1]*fr[4]-0.9375*nuVtSqSum[1]*fc[4]; 
  Gdiff_r[6] = 0.6708203932499369*nuVtSqSum[1]*fr[45]-0.6708203932499369*nuVtSqSum[1]*fc[45]-1.190784930203603*nuVtSqSum[2]*fr[44]-1.190784930203603*nuVtSqSum[2]*fc[44]+0.6708203932499369*nuVtSqSum[0]*fr[38]-0.6708203932499369*nuVtSqSum[0]*fc[38]+0.9375*nuVtSqSum[2]*fr[37]-0.9375*nuVtSqSum[2]*fc[37]-1.190784930203603*nuVtSqSum[1]*fr[31]-1.190784930203603*nuVtSqSum[1]*fc[31]-1.190784930203603*nuVtSqSum[0]*fr[18]-1.190784930203603*nuVtSqSum[0]*fc[18]+0.9375*nuVtSqSum[1]*fr[17]-0.9375*nuVtSqSum[1]*fc[17]+0.9375*nuVtSqSum[0]*fr[10]-0.9375*nuVtSqSum[0]*fc[10]; 
  Gdiff_r[7] = 0.5999999999999999*nuVtSqSum[1]*fr[20]-0.5999999999999999*nuVtSqSum[1]*fc[20]-0.7607645858621712*nuVtSqSum[2]*fr[19]-1.190784930203603*nuVtSqSum[0]*fr[19]-0.7607645858621712*nuVtSqSum[2]*fc[19]-1.190784930203603*nuVtSqSum[0]*fc[19]+0.6708203932499369*nuVtSqSum[2]*fr[12]-0.6708203932499369*nuVtSqSum[2]*fc[12]+0.5989467796874438*nuVtSqSum[2]*fr[11]+0.9375*nuVtSqSum[0]*fr[11]-0.5989467796874438*nuVtSqSum[2]*fc[11]-0.9375*nuVtSqSum[0]*fc[11]-1.06507042020704*nuVtSqSum[1]*fr[5]-1.06507042020704*nuVtSqSum[1]*fc[5]-1.190784930203603*fr[2]*nuVtSqSum[2]-1.190784930203603*fc[2]*nuVtSqSum[2]+0.9375*fr[0]*nuVtSqSum[2]-0.9375*fc[0]*nuVtSqSum[2]+0.8385254915624212*fr[1]*nuVtSqSum[1]-0.8385254915624212*fc[1]*nuVtSqSum[1]; 
  Gdiff_r[8] = (-1.190784930203603*nuVtSqSum[1]*fr[34])-1.190784930203603*nuVtSqSum[1]*fc[34]-1.190784930203603*nuVtSqSum[0]*fr[24]-1.190784930203603*nuVtSqSum[0]*fc[24]+0.9375000000000001*nuVtSqSum[1]*fr[23]-0.9375000000000001*nuVtSqSum[1]*fc[23]+0.9375*nuVtSqSum[0]*fr[13]-0.9375*nuVtSqSum[0]*fc[13]; 
  Gdiff_r[9] = (-1.190784930203603*nuVtSqSum[1]*fr[41])-1.190784930203603*nuVtSqSum[1]*fc[41]-1.190784930203603*nuVtSqSum[0]*fr[29]-1.190784930203603*nuVtSqSum[0]*fc[29]+0.9375000000000001*nuVtSqSum[1]*fr[28]-0.9375000000000001*nuVtSqSum[1]*fc[28]+0.9375*nuVtSqSum[0]*fr[14]-0.9375*nuVtSqSum[0]*fc[14]; 
  Gdiff_r[10] = 0.5999999999999999*nuVtSqSum[2]*fr[45]+0.6708203932499369*nuVtSqSum[0]*fr[45]-0.5999999999999999*nuVtSqSum[2]*fc[45]-0.6708203932499369*nuVtSqSum[0]*fc[45]-1.06507042020704*nuVtSqSum[1]*fr[44]-1.06507042020704*nuVtSqSum[1]*fc[44]+0.6708203932499369*nuVtSqSum[1]*fr[38]-0.6708203932499369*nuVtSqSum[1]*fc[38]+0.8385254915624212*nuVtSqSum[1]*fr[37]-0.8385254915624212*nuVtSqSum[1]*fc[37]-1.06507042020704*nuVtSqSum[2]*fr[31]-1.190784930203603*nuVtSqSum[0]*fr[31]-1.06507042020704*nuVtSqSum[2]*fc[31]-1.190784930203603*nuVtSqSum[0]*fc[31]-1.190784930203603*nuVtSqSum[1]*fr[18]-1.190784930203603*nuVtSqSum[1]*fc[18]+0.8385254915624212*nuVtSqSum[2]*fr[17]+0.9375*nuVtSqSum[0]*fr[17]-0.8385254915624212*nuVtSqSum[2]*fc[17]-0.9375*nuVtSqSum[0]*fc[17]+0.9375*nuVtSqSum[1]*fr[10]-0.9375*nuVtSqSum[1]*fc[10]; 
  Gdiff_r[11] = 0.5999999999999999*nuVtSqSum[1]*fr[33]-0.5999999999999999*nuVtSqSum[1]*fc[33]-0.7607645858621712*nuVtSqSum[2]*fr[32]-1.190784930203603*nuVtSqSum[0]*fr[32]-0.7607645858621712*nuVtSqSum[2]*fc[32]-1.190784930203603*nuVtSqSum[0]*fc[32]+0.6708203932499369*nuVtSqSum[2]*fr[22]-0.6708203932499369*nuVtSqSum[2]*fc[22]+0.5989467796874438*nuVtSqSum[2]*fr[21]+0.9375*nuVtSqSum[0]*fr[21]-0.5989467796874438*nuVtSqSum[2]*fc[21]-0.9375*nuVtSqSum[0]*fc[21]-1.06507042020704*nuVtSqSum[1]*fr[15]-1.06507042020704*nuVtSqSum[1]*fc[15]-1.190784930203603*nuVtSqSum[2]*fr[7]-1.190784930203603*nuVtSqSum[2]*fc[7]+0.8385254915624211*nuVtSqSum[1]*fr[6]-0.8385254915624211*nuVtSqSum[1]*fc[6]+0.9375000000000001*nuVtSqSum[2]*fr[3]-0.9375000000000001*nuVtSqSum[2]*fc[3]; 
  Gdiff_r[12] = (-1.06507042020704*nuVtSqSum[2]*fr[34])-1.190784930203603*nuVtSqSum[0]*fr[34]-1.06507042020704*nuVtSqSum[2]*fc[34]-1.190784930203603*nuVtSqSum[0]*fc[34]-1.190784930203603*nuVtSqSum[1]*fr[24]-1.190784930203603*nuVtSqSum[1]*fc[24]+0.8385254915624212*nuVtSqSum[2]*fr[23]+0.9375*nuVtSqSum[0]*fr[23]-0.8385254915624212*nuVtSqSum[2]*fc[23]-0.9375*nuVtSqSum[0]*fc[23]+0.9375000000000001*nuVtSqSum[1]*fr[13]-0.9375000000000001*nuVtSqSum[1]*fc[13]; 
  Gdiff_r[13] = 0.5999999999999999*nuVtSqSum[1]*fr[36]-0.5999999999999999*nuVtSqSum[1]*fc[36]-0.7607645858621712*nuVtSqSum[2]*fr[35]-1.190784930203603*nuVtSqSum[0]*fr[35]-0.7607645858621712*nuVtSqSum[2]*fc[35]-1.190784930203603*nuVtSqSum[0]*fc[35]+0.6708203932499369*nuVtSqSum[2]*fr[26]-0.6708203932499369*nuVtSqSum[2]*fc[26]+0.5989467796874438*nuVtSqSum[2]*fr[25]+0.9375*nuVtSqSum[0]*fr[25]-0.5989467796874438*nuVtSqSum[2]*fc[25]-0.9375*nuVtSqSum[0]*fc[25]-1.06507042020704*nuVtSqSum[1]*fr[16]-1.06507042020704*nuVtSqSum[1]*fc[16]-1.190784930203603*nuVtSqSum[2]*fr[9]-1.190784930203603*nuVtSqSum[2]*fc[9]+0.8385254915624211*nuVtSqSum[1]*fr[8]-0.8385254915624211*nuVtSqSum[1]*fc[8]+0.9375000000000001*nuVtSqSum[2]*fr[4]-0.9375000000000001*nuVtSqSum[2]*fc[4]; 
  Gdiff_r[14] = (-1.190784930203603*nuVtSqSum[1]*fr[46])-1.190784930203603*nuVtSqSum[1]*fc[46]-1.190784930203603*nuVtSqSum[0]*fr[40]-1.190784930203603*nuVtSqSum[0]*fc[40]+0.9375000000000001*nuVtSqSum[1]*fr[39]-0.9375000000000001*nuVtSqSum[1]*fc[39]+0.9375*nuVtSqSum[0]*fr[27]-0.9375*nuVtSqSum[0]*fc[27]; 
  Gdiff_r[15] = (-1.06507042020704*nuVtSqSum[2]*fr[41])-1.190784930203603*nuVtSqSum[0]*fr[41]-1.06507042020704*nuVtSqSum[2]*fc[41]-1.190784930203603*nuVtSqSum[0]*fc[41]-1.190784930203603*nuVtSqSum[1]*fr[29]-1.190784930203603*nuVtSqSum[1]*fc[29]+0.8385254915624212*nuVtSqSum[2]*fr[28]+0.9375*nuVtSqSum[0]*fr[28]-0.8385254915624212*nuVtSqSum[2]*fc[28]-0.9375*nuVtSqSum[0]*fc[28]+0.9375000000000001*nuVtSqSum[1]*fr[14]-0.9375000000000001*nuVtSqSum[1]*fc[14]; 
  Gdiff_r[16] = (-1.190784930203603*nuVtSqSum[1]*fr[47])-1.190784930203603*nuVtSqSum[1]*fc[47]-1.190784930203603*nuVtSqSum[0]*fr[43]-1.190784930203603*nuVtSqSum[0]*fc[43]+0.9375000000000001*nuVtSqSum[1]*fr[42]-0.9375000000000001*nuVtSqSum[1]*fc[42]+0.9375*nuVtSqSum[0]*fr[30]-0.9375*nuVtSqSum[0]*fc[30]; 
  Gdiff_r[17] = 0.5999999999999999*nuVtSqSum[1]*fr[45]-0.5999999999999999*nuVtSqSum[1]*fc[45]-0.7607645858621712*nuVtSqSum[2]*fr[44]-1.190784930203603*nuVtSqSum[0]*fr[44]-0.7607645858621712*nuVtSqSum[2]*fc[44]-1.190784930203603*nuVtSqSum[0]*fc[44]+0.6708203932499369*nuVtSqSum[2]*fr[38]-0.6708203932499369*nuVtSqSum[2]*fc[38]+0.5989467796874438*nuVtSqSum[2]*fr[37]+0.9375*nuVtSqSum[0]*fr[37]-0.5989467796874438*nuVtSqSum[2]*fc[37]-0.9375*nuVtSqSum[0]*fc[37]-1.06507042020704*nuVtSqSum[1]*fr[31]-1.06507042020704*nuVtSqSum[1]*fc[31]-1.190784930203603*nuVtSqSum[2]*fr[18]-1.190784930203603*nuVtSqSum[2]*fc[18]+0.8385254915624212*nuVtSqSum[1]*fr[17]-0.8385254915624212*nuVtSqSum[1]*fc[17]+0.9375*nuVtSqSum[2]*fr[10]-0.9375*nuVtSqSum[2]*fc[10]; 
  Gdiff_r[18] = (-1.06507042020704*nuVtSqSum[2]*fr[46])-1.190784930203603*nuVtSqSum[0]*fr[46]-1.06507042020704*nuVtSqSum[2]*fc[46]-1.190784930203603*nuVtSqSum[0]*fc[46]-1.190784930203603*nuVtSqSum[1]*fr[40]-1.190784930203603*nuVtSqSum[1]*fc[40]+0.8385254915624212*nuVtSqSum[2]*fr[39]+0.9375*nuVtSqSum[0]*fr[39]-0.8385254915624212*nuVtSqSum[2]*fc[39]-0.9375*nuVtSqSum[0]*fc[39]+0.9375000000000001*nuVtSqSum[1]*fr[27]-0.9375000000000001*nuVtSqSum[1]*fc[27]; 
  Gdiff_r[19] = (-1.06507042020704*nuVtSqSum[2]*fr[47])-1.190784930203603*nuVtSqSum[0]*fr[47]-1.06507042020704*nuVtSqSum[2]*fc[47]-1.190784930203603*nuVtSqSum[0]*fc[47]-1.190784930203603*nuVtSqSum[1]*fr[43]-1.190784930203603*nuVtSqSum[1]*fc[43]+0.8385254915624212*nuVtSqSum[2]*fr[42]+0.9375*nuVtSqSum[0]*fr[42]-0.8385254915624212*nuVtSqSum[2]*fc[42]-0.9375*nuVtSqSum[0]*fc[42]+0.9375000000000001*nuVtSqSum[1]*fr[30]-0.9375000000000001*nuVtSqSum[1]*fc[30]; 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.3535533905932737*alphaDrSurf_l[7]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[1]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.3535533905932737*alphaDrSurf_l[7]*fUpwind_l[11]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.3535533905932737*alphaDrSurf_l[7]*fUpwind_l[13]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Ghat_l[4] = Gdiff_l[4]*rdv2+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[4]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[2]; 
  Ghat_l[5] = Gdiff_l[5]*rdv2+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[5]*alphaDrSurf_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Ghat_l[6] = Gdiff_l[6]*rdv2+0.3535533905932737*alphaDrSurf_l[7]*fUpwind_l[17]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[10]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[6]; 
  Ghat_l[7] = Gdiff_l[7]*rdv2+0.2258769757263128*alphaDrSurf_l[7]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[7]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[8] = Gdiff_l[8]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[12]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[8]; 
  Ghat_l[9] = Gdiff_l[9]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[15]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[9]; 
  Ghat_l[10] = Gdiff_l[10]*rdv2+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[17]+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[10]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[10]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[6]; 
  Ghat_l[11] = Gdiff_l[11]*rdv2+0.2258769757263128*alphaDrSurf_l[7]*fUpwind_l[11]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[7]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[4]; 
  Ghat_l[12] = Gdiff_l[12]*rdv2+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[12]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[12]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[8]; 
  Ghat_l[13] = Gdiff_l[13]*rdv2+0.2258769757263128*alphaDrSurf_l[7]*fUpwind_l[13]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[7]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[5]; 
  Ghat_l[14] = Gdiff_l[14]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[18]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[14]; 
  Ghat_l[15] = Gdiff_l[15]*rdv2+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[15]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[15]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[9]; 
  Ghat_l[16] = Gdiff_l[16]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[19]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[16]; 
  Ghat_l[17] = Gdiff_l[17]*rdv2+0.2258769757263128*alphaDrSurf_l[7]*fUpwind_l[17]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[17]+0.3162277660168379*alphaDrSurf_l[1]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[6]*alphaDrSurf_l[7]; 
  Ghat_l[18] = Gdiff_l[18]*rdv2+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[18]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[18]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[14]; 
  Ghat_l[19] = Gdiff_l[19]*rdv2+0.3162277660168379*alphaDrSurf_l[7]*fUpwind_l[19]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[19]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[16]; 

  Ghat_r[0] = Gdiff_r[0]*rdv2+0.3535533905932737*alphaDrSurf_r[7]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[1]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = Gdiff_r[1]*rdv2+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[1]*alphaDrSurf_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = Gdiff_r[2]*rdv2+0.3535533905932737*alphaDrSurf_r[7]*fUpwind_r[11]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Ghat_r[3] = Gdiff_r[3]*rdv2+0.3535533905932737*alphaDrSurf_r[7]*fUpwind_r[13]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Ghat_r[4] = Gdiff_r[4]*rdv2+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[4]*alphaDrSurf_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[2]; 
  Ghat_r[5] = Gdiff_r[5]*rdv2+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[5]*alphaDrSurf_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Ghat_r[6] = Gdiff_r[6]*rdv2+0.3535533905932737*alphaDrSurf_r[7]*fUpwind_r[17]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[10]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[6]; 
  Ghat_r[7] = Gdiff_r[7]*rdv2+0.2258769757263128*alphaDrSurf_r[7]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[7]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[8] = Gdiff_r[8]*rdv2+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[12]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[8]; 
  Ghat_r[9] = Gdiff_r[9]*rdv2+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[15]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[9]; 
  Ghat_r[10] = Gdiff_r[10]*rdv2+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[17]+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[10]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[10]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[6]; 
  Ghat_r[11] = Gdiff_r[11]*rdv2+0.2258769757263128*alphaDrSurf_r[7]*fUpwind_r[11]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[2]*alphaDrSurf_r[7]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[4]; 
  Ghat_r[12] = Gdiff_r[12]*rdv2+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[12]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[12]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[8]; 
  Ghat_r[13] = Gdiff_r[13]*rdv2+0.2258769757263128*alphaDrSurf_r[7]*fUpwind_r[13]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[3]*alphaDrSurf_r[7]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[5]; 
  Ghat_r[14] = Gdiff_r[14]*rdv2+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[18]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[14]; 
  Ghat_r[15] = Gdiff_r[15]*rdv2+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[15]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[15]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[9]; 
  Ghat_r[16] = Gdiff_r[16]*rdv2+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[19]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[16]; 
  Ghat_r[17] = Gdiff_r[17]*rdv2+0.2258769757263128*alphaDrSurf_r[7]*fUpwind_r[17]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[17]+0.3162277660168379*alphaDrSurf_r[1]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[6]*alphaDrSurf_r[7]; 
  Ghat_r[18] = Gdiff_r[18]*rdv2+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[18]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[18]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[14]; 
  Ghat_r[19] = Gdiff_r[19]*rdv2+0.3162277660168379*alphaDrSurf_r[7]*fUpwind_r[19]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[19]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[16]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2-0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2-0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2-1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_l[2]*rdv2-0.7071067811865475*Ghat_r[2]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_l[3]*rdv2-0.7071067811865475*Ghat_r[3]*rdv2; 
  out[5] += 1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2-1.224744871391589*Ghat_l[1]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_l[4]*rdv2-0.7071067811865475*Ghat_r[4]*rdv2; 
  out[7] += 1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2-1.224744871391589*Ghat_l[2]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_l[5]*rdv2-0.7071067811865475*Ghat_r[5]*rdv2; 
  out[9] += 1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2-1.224744871391589*Ghat_l[3]*rdv2; 
  out[10] += 0.7071067811865475*Ghat_l[6]*rdv2-0.7071067811865475*Ghat_r[6]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_l[7]*rdv2-0.7071067811865475*Ghat_r[7]*rdv2; 
  out[12] += 4.743416490252569*Gdiff2_r[0]*rdvSq4+4.743416490252569*Gdiff2_l[0]*rdvSq4-1.58113883008419*Ghat_r[0]*rdv2+1.58113883008419*Ghat_l[0]*rdv2; 
  out[13] += 0.7071067811865475*Ghat_l[8]*rdv2-0.7071067811865475*Ghat_r[8]*rdv2; 
  out[14] += 0.7071067811865475*Ghat_l[9]*rdv2-0.7071067811865475*Ghat_r[9]*rdv2; 
  out[15] += 1.224744871391589*Gdiff2_r[4]*rdvSq4-1.224744871391589*Gdiff2_l[4]*rdvSq4-1.224744871391589*Ghat_r[4]*rdv2-1.224744871391589*Ghat_l[4]*rdv2; 
  out[16] += 1.224744871391589*Gdiff2_r[5]*rdvSq4-1.224744871391589*Gdiff2_l[5]*rdvSq4-1.224744871391589*Ghat_r[5]*rdv2-1.224744871391589*Ghat_l[5]*rdv2; 
  out[17] += 0.7071067811865475*Ghat_l[10]*rdv2-0.7071067811865475*Ghat_r[10]*rdv2; 
  out[18] += 1.224744871391589*Gdiff2_r[6]*rdvSq4-1.224744871391589*Gdiff2_l[6]*rdvSq4-1.224744871391589*Ghat_r[6]*rdv2-1.224744871391589*Ghat_l[6]*rdv2; 
  out[19] += 1.224744871391589*Gdiff2_r[7]*rdvSq4-1.224744871391589*Gdiff2_l[7]*rdvSq4-1.224744871391589*Ghat_r[7]*rdv2-1.224744871391589*Ghat_l[7]*rdv2; 
  out[20] += 4.743416490252569*Gdiff2_r[1]*rdvSq4+4.743416490252569*Gdiff2_l[1]*rdvSq4-1.58113883008419*Ghat_r[1]*rdv2+1.58113883008419*Ghat_l[1]*rdv2; 
  out[21] += 0.7071067811865475*Ghat_l[11]*rdv2-0.7071067811865475*Ghat_r[11]*rdv2; 
  out[22] += 4.743416490252569*Gdiff2_r[2]*rdvSq4+4.743416490252569*Gdiff2_l[2]*rdvSq4-1.58113883008419*Ghat_r[2]*rdv2+1.58113883008419*Ghat_l[2]*rdv2; 
  out[23] += 0.7071067811865475*Ghat_l[12]*rdv2-0.7071067811865475*Ghat_r[12]*rdv2; 
  out[24] += 1.224744871391589*Gdiff2_r[8]*rdvSq4-1.224744871391589*Gdiff2_l[8]*rdvSq4-1.224744871391589*Ghat_r[8]*rdv2-1.224744871391589*Ghat_l[8]*rdv2; 
  out[25] += 0.7071067811865475*Ghat_l[13]*rdv2-0.7071067811865475*Ghat_r[13]*rdv2; 
  out[26] += 4.743416490252569*Gdiff2_r[3]*rdvSq4+4.743416490252569*Gdiff2_l[3]*rdvSq4-1.58113883008419*Ghat_r[3]*rdv2+1.58113883008419*Ghat_l[3]*rdv2; 
  out[27] += 0.7071067811865475*Ghat_l[14]*rdv2-0.7071067811865475*Ghat_r[14]*rdv2; 
  out[28] += 0.7071067811865475*Ghat_l[15]*rdv2-0.7071067811865475*Ghat_r[15]*rdv2; 
  out[29] += 1.224744871391589*Gdiff2_r[9]*rdvSq4-1.224744871391589*Gdiff2_l[9]*rdvSq4-1.224744871391589*Ghat_r[9]*rdv2-1.224744871391589*Ghat_l[9]*rdv2; 
  out[30] += 0.7071067811865475*Ghat_l[16]*rdv2-0.7071067811865475*Ghat_r[16]*rdv2; 
  out[31] += 1.224744871391589*Gdiff2_r[10]*rdvSq4-1.224744871391589*Gdiff2_l[10]*rdvSq4-1.224744871391589*Ghat_r[10]*rdv2-1.224744871391589*Ghat_l[10]*rdv2; 
  out[32] += 1.224744871391589*Gdiff2_r[11]*rdvSq4-1.224744871391589*Gdiff2_l[11]*rdvSq4-1.224744871391589*Ghat_r[11]*rdv2-1.224744871391589*Ghat_l[11]*rdv2; 
  out[33] += 4.743416490252569*Gdiff2_r[4]*rdvSq4+4.743416490252569*Gdiff2_l[4]*rdvSq4-1.58113883008419*Ghat_r[4]*rdv2+1.58113883008419*Ghat_l[4]*rdv2; 
  out[34] += 1.224744871391589*Gdiff2_r[12]*rdvSq4-1.224744871391589*Gdiff2_l[12]*rdvSq4-1.224744871391589*Ghat_r[12]*rdv2-1.224744871391589*Ghat_l[12]*rdv2; 
  out[35] += 1.224744871391589*Gdiff2_r[13]*rdvSq4-1.224744871391589*Gdiff2_l[13]*rdvSq4-1.224744871391589*Ghat_r[13]*rdv2-1.224744871391589*Ghat_l[13]*rdv2; 
  out[36] += 4.743416490252569*Gdiff2_r[5]*rdvSq4+4.743416490252569*Gdiff2_l[5]*rdvSq4-1.58113883008419*Ghat_r[5]*rdv2+1.58113883008419*Ghat_l[5]*rdv2; 
  out[37] += 0.7071067811865475*Ghat_l[17]*rdv2-0.7071067811865475*Ghat_r[17]*rdv2; 
  out[38] += 4.743416490252569*Gdiff2_r[6]*rdvSq4+4.743416490252569*Gdiff2_l[6]*rdvSq4-1.58113883008419*Ghat_r[6]*rdv2+1.58113883008419*Ghat_l[6]*rdv2; 
  out[39] += 0.7071067811865475*Ghat_l[18]*rdv2-0.7071067811865475*Ghat_r[18]*rdv2; 
  out[40] += 1.224744871391589*Gdiff2_r[14]*rdvSq4-1.224744871391589*Gdiff2_l[14]*rdvSq4-1.224744871391589*Ghat_r[14]*rdv2-1.224744871391589*Ghat_l[14]*rdv2; 
  out[41] += 1.224744871391589*Gdiff2_r[15]*rdvSq4-1.224744871391589*Gdiff2_l[15]*rdvSq4-1.224744871391589*Ghat_r[15]*rdv2-1.224744871391589*Ghat_l[15]*rdv2; 
  out[42] += 0.7071067811865475*Ghat_l[19]*rdv2-0.7071067811865475*Ghat_r[19]*rdv2; 
  out[43] += 1.224744871391589*Gdiff2_r[16]*rdvSq4-1.224744871391589*Gdiff2_l[16]*rdvSq4-1.224744871391589*Ghat_r[16]*rdv2-1.224744871391589*Ghat_l[16]*rdv2; 
  out[44] += 1.224744871391589*Gdiff2_r[17]*rdvSq4-1.224744871391589*Gdiff2_l[17]*rdvSq4-1.224744871391589*Ghat_r[17]*rdv2-1.224744871391589*Ghat_l[17]*rdv2; 
  out[45] += 4.743416490252569*Gdiff2_r[10]*rdvSq4+4.743416490252569*Gdiff2_l[10]*rdvSq4-1.58113883008419*Ghat_r[10]*rdv2+1.58113883008419*Ghat_l[10]*rdv2; 
  out[46] += 1.224744871391589*Gdiff2_r[18]*rdvSq4-1.224744871391589*Gdiff2_l[18]*rdvSq4-1.224744871391589*Ghat_r[18]*rdv2-1.224744871391589*Ghat_l[18]*rdv2; 
  out[47] += 1.224744871391589*Gdiff2_r[19]*rdvSq4-1.224744871391589*Gdiff2_l[19]*rdvSq4-1.224744871391589*Ghat_r[19]*rdv2-1.224744871391589*Ghat_l[19]*rdv2; 
} 
