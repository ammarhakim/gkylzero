#include <gkyl_lbo_gyrokinetic_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx3_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_surfvpar_2x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:     cell-center coordinates. 
  // dxv[4]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[16]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:  distribution function in cells 
  // out:       incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  double alphaDrSurf_l[20] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[2]-0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*nuUSum[0]; 
  alphaDrSurf_l[1] = 1.414213562373095*nuSum[1]*w[2]-0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*nuUSum[1]; 
  alphaDrSurf_l[2] = 1.414213562373095*nuSum[2]*w[2]-1.414213562373095*nuUSum[2]-0.7071067811865475*dxv[2]*nuSum[2]; 
  alphaDrSurf_l[4] = (-1.414213562373095*nuUSum[3])+1.414213562373095*w[2]*nuSum[3]-0.7071067811865475*dxv[2]*nuSum[3]; 
  alphaDrSurf_l[7] = (-1.414213562373095*nuUSum[4])+1.414213562373095*w[2]*nuSum[4]-0.7071067811865475*dxv[2]*nuSum[4]; 
  alphaDrSurf_l[8] = (-1.414213562373095*nuUSum[5])+1.414213562373095*w[2]*nuSum[5]-0.7071067811865475*dxv[2]*nuSum[5]; 
  alphaDrSurf_l[11] = (-1.414213562373095*nuUSum[6])+1.414213562373095*w[2]*nuSum[6]-0.7071067811865475*dxv[2]*nuSum[6]; 
  alphaDrSurf_l[12] = (-1.414213562373095*nuUSum[7])+1.414213562373095*w[2]*nuSum[7]-0.7071067811865475*dxv[2]*nuSum[7]; 

  double alphaDrSurf_r[20] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[2]+0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*nuUSum[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*nuSum[1]*w[2]+0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*nuUSum[1]; 
  alphaDrSurf_r[2] = 1.414213562373095*nuSum[2]*w[2]-1.414213562373095*nuUSum[2]+0.7071067811865475*dxv[2]*nuSum[2]; 
  alphaDrSurf_r[4] = (-1.414213562373095*nuUSum[3])+1.414213562373095*w[2]*nuSum[3]+0.7071067811865475*dxv[2]*nuSum[3]; 
  alphaDrSurf_r[7] = (-1.414213562373095*nuUSum[4])+1.414213562373095*w[2]*nuSum[4]+0.7071067811865475*dxv[2]*nuSum[4]; 
  alphaDrSurf_r[8] = (-1.414213562373095*nuUSum[5])+1.414213562373095*w[2]*nuSum[5]+0.7071067811865475*dxv[2]*nuSum[5]; 
  alphaDrSurf_r[11] = (-1.414213562373095*nuUSum[6])+1.414213562373095*w[2]*nuSum[6]+0.7071067811865475*dxv[2]*nuSum[6]; 
  alphaDrSurf_r[12] = (-1.414213562373095*nuUSum[7])+1.414213562373095*w[2]*nuSum[7]+0.7071067811865475*dxv[2]*nuSum[7]; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};
  double fUpwind_r[20] = {0.0};
  double Gdrag_l[20] = {0.0}; 
  double Gdrag_r[20] = {0.0}; 

  if ((-0.4242640687119281*alphaDrSurf_l[12])-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx3_quad_0_r(fl); 
    fUpwindQuad_l[9] = ser_4x_p2_surfx3_quad_9_r(fl); 
    fUpwindQuad_l[18] = ser_4x_p2_surfx3_quad_18_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx3_quad_0_l(fc); 
    fUpwindQuad_l[9] = ser_4x_p2_surfx3_quad_9_l(fc); 
    fUpwindQuad_l[18] = ser_4x_p2_surfx3_quad_18_l(fc); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_r[12])-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_4x_p2_surfx3_quad_0_r(fc); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx3_quad_9_r(fc); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx3_quad_18_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_4x_p2_surfx3_quad_0_l(fr); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx3_quad_9_l(fr); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx3_quad_18_l(fr); 
  } 
  if (0.5303300858899104*alphaDrSurf_l[11]+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx3_quad_1_r(fl); 
    fUpwindQuad_l[10] = ser_4x_p2_surfx3_quad_10_r(fl); 
    fUpwindQuad_l[19] = ser_4x_p2_surfx3_quad_19_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx3_quad_1_l(fc); 
    fUpwindQuad_l[10] = ser_4x_p2_surfx3_quad_10_l(fc); 
    fUpwindQuad_l[19] = ser_4x_p2_surfx3_quad_19_l(fc); 
  } 
  if (0.5303300858899104*alphaDrSurf_r[11]+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_4x_p2_surfx3_quad_1_r(fc); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx3_quad_10_r(fc); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx3_quad_19_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_4x_p2_surfx3_quad_1_l(fr); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx3_quad_10_l(fr); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx3_quad_19_l(fr); 
  } 
  if (0.4242640687119281*alphaDrSurf_l[12]-0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]-0.4743416490252568*alphaDrSurf_l[2]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx3_quad_2_r(fl); 
    fUpwindQuad_l[11] = ser_4x_p2_surfx3_quad_11_r(fl); 
    fUpwindQuad_l[20] = ser_4x_p2_surfx3_quad_20_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx3_quad_2_l(fc); 
    fUpwindQuad_l[11] = ser_4x_p2_surfx3_quad_11_l(fc); 
    fUpwindQuad_l[20] = ser_4x_p2_surfx3_quad_20_l(fc); 
  } 
  if (0.4242640687119281*alphaDrSurf_r[12]-0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]-0.4743416490252568*alphaDrSurf_r[2]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_4x_p2_surfx3_quad_2_r(fc); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx3_quad_11_r(fc); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx3_quad_20_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_4x_p2_surfx3_quad_2_l(fr); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx3_quad_11_l(fr); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx3_quad_20_l(fr); 
  } 
  if (0.5303300858899104*alphaDrSurf_l[12]-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx3_quad_3_r(fl); 
    fUpwindQuad_l[12] = ser_4x_p2_surfx3_quad_12_r(fl); 
    fUpwindQuad_l[21] = ser_4x_p2_surfx3_quad_21_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx3_quad_3_l(fc); 
    fUpwindQuad_l[12] = ser_4x_p2_surfx3_quad_12_l(fc); 
    fUpwindQuad_l[21] = ser_4x_p2_surfx3_quad_21_l(fc); 
  } 
  if (0.5303300858899104*alphaDrSurf_r[12]-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_4x_p2_surfx3_quad_3_r(fc); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx3_quad_12_r(fc); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx3_quad_21_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_4x_p2_surfx3_quad_3_l(fr); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx3_quad_12_l(fr); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx3_quad_21_l(fr); 
  } 
  if (0.3535533905932737*alphaDrSurf_l[0]-0.3952847075210473*(alphaDrSurf_l[8]+alphaDrSurf_l[7]) < 0) { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx3_quad_4_r(fl); 
    fUpwindQuad_l[13] = ser_4x_p2_surfx3_quad_13_r(fl); 
    fUpwindQuad_l[22] = ser_4x_p2_surfx3_quad_22_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx3_quad_4_l(fc); 
    fUpwindQuad_l[13] = ser_4x_p2_surfx3_quad_13_l(fc); 
    fUpwindQuad_l[22] = ser_4x_p2_surfx3_quad_22_l(fc); 
  } 
  if (0.3535533905932737*alphaDrSurf_r[0]-0.3952847075210473*(alphaDrSurf_r[8]+alphaDrSurf_r[7]) < 0) { 
    fUpwindQuad_r[4] = ser_4x_p2_surfx3_quad_4_r(fc); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx3_quad_13_r(fc); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx3_quad_22_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_4x_p2_surfx3_quad_4_l(fr); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx3_quad_13_l(fr); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx3_quad_22_l(fr); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_l[12])-0.3952847075210473*alphaDrSurf_l[8]+0.3162277660168379*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx3_quad_5_r(fl); 
    fUpwindQuad_l[14] = ser_4x_p2_surfx3_quad_14_r(fl); 
    fUpwindQuad_l[23] = ser_4x_p2_surfx3_quad_23_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx3_quad_5_l(fc); 
    fUpwindQuad_l[14] = ser_4x_p2_surfx3_quad_14_l(fc); 
    fUpwindQuad_l[23] = ser_4x_p2_surfx3_quad_23_l(fc); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_r[12])-0.3952847075210473*alphaDrSurf_r[8]+0.3162277660168379*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_4x_p2_surfx3_quad_5_r(fc); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx3_quad_14_r(fc); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx3_quad_23_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_4x_p2_surfx3_quad_5_l(fr); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx3_quad_14_l(fr); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx3_quad_23_l(fr); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_l[12])+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])-0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*alphaDrSurf_l[2]-0.4743416490252568*alphaDrSurf_l[1]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx3_quad_6_r(fl); 
    fUpwindQuad_l[15] = ser_4x_p2_surfx3_quad_15_r(fl); 
    fUpwindQuad_l[24] = ser_4x_p2_surfx3_quad_24_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx3_quad_6_l(fc); 
    fUpwindQuad_l[15] = ser_4x_p2_surfx3_quad_15_l(fc); 
    fUpwindQuad_l[24] = ser_4x_p2_surfx3_quad_24_l(fc); 
  } 
  if ((-0.4242640687119281*alphaDrSurf_r[12])+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])-0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*alphaDrSurf_r[2]-0.4743416490252568*alphaDrSurf_r[1]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_4x_p2_surfx3_quad_6_r(fc); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx3_quad_15_r(fc); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx3_quad_24_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_4x_p2_surfx3_quad_6_l(fr); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx3_quad_15_l(fr); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx3_quad_24_l(fr); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_l[11])+0.3162277660168379*alphaDrSurf_l[8]-0.3952847075210473*alphaDrSurf_l[7]+0.4743416490252568*alphaDrSurf_l[2]+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx3_quad_7_r(fl); 
    fUpwindQuad_l[16] = ser_4x_p2_surfx3_quad_16_r(fl); 
    fUpwindQuad_l[25] = ser_4x_p2_surfx3_quad_25_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx3_quad_7_l(fc); 
    fUpwindQuad_l[16] = ser_4x_p2_surfx3_quad_16_l(fc); 
    fUpwindQuad_l[25] = ser_4x_p2_surfx3_quad_25_l(fc); 
  } 
  if ((-0.5303300858899104*alphaDrSurf_r[11])+0.3162277660168379*alphaDrSurf_r[8]-0.3952847075210473*alphaDrSurf_r[7]+0.4743416490252568*alphaDrSurf_r[2]+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[7] = ser_4x_p2_surfx3_quad_7_r(fc); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx3_quad_16_r(fc); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx3_quad_25_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_4x_p2_surfx3_quad_7_l(fr); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx3_quad_16_l(fr); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx3_quad_25_l(fr); 
  } 
  if (0.4242640687119281*alphaDrSurf_l[12]+0.4242640687119285*alphaDrSurf_l[11]+0.3162277660168379*(alphaDrSurf_l[8]+alphaDrSurf_l[7])+0.6363961030678926*alphaDrSurf_l[4]+0.4743416490252568*(alphaDrSurf_l[2]+alphaDrSurf_l[1])+0.3535533905932737*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx3_quad_8_r(fl); 
    fUpwindQuad_l[17] = ser_4x_p2_surfx3_quad_17_r(fl); 
    fUpwindQuad_l[26] = ser_4x_p2_surfx3_quad_26_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx3_quad_8_l(fc); 
    fUpwindQuad_l[17] = ser_4x_p2_surfx3_quad_17_l(fc); 
    fUpwindQuad_l[26] = ser_4x_p2_surfx3_quad_26_l(fc); 
  } 
  if (0.4242640687119281*alphaDrSurf_r[12]+0.4242640687119285*alphaDrSurf_r[11]+0.3162277660168379*(alphaDrSurf_r[8]+alphaDrSurf_r[7])+0.6363961030678926*alphaDrSurf_r[4]+0.4743416490252568*(alphaDrSurf_r[2]+alphaDrSurf_r[1])+0.3535533905932737*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = ser_4x_p2_surfx3_quad_8_r(fc); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx3_quad_17_r(fc); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx3_quad_26_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_4x_p2_surfx3_quad_8_l(fr); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx3_quad_17_l(fr); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx3_quad_26_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p2_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p2_upwind(fUpwindQuad_r, fUpwind_r); 

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
  out[3] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[5] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[6] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[7] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[8] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[9] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[10] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[12] += 0.7071067811865475*Gdrag_r[8]*rdv2-0.7071067811865475*Gdrag_l[8]*rdv2; 
  out[13] += 1.58113883008419*Gdrag_r[0]*rdv2-1.58113883008419*Gdrag_l[0]*rdv2; 
  out[14] += 0.7071067811865475*Gdrag_r[9]*rdv2-0.7071067811865475*Gdrag_l[9]*rdv2; 
  out[15] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[16] += 0.7071067811865475*Gdrag_r[10]*rdv2-0.7071067811865475*Gdrag_l[10]*rdv2; 
  out[17] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[18] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[19] += 0.7071067811865475*Gdrag_r[11]*rdv2-0.7071067811865475*Gdrag_l[11]*rdv2; 
  out[20] += 0.7071067811865475*Gdrag_r[12]*rdv2-0.7071067811865475*Gdrag_l[12]*rdv2; 
  out[21] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
  out[22] += 1.224744871391589*Gdrag_r[8]*rdv2+1.224744871391589*Gdrag_l[8]*rdv2; 
  out[23] += 1.58113883008419*Gdrag_r[1]*rdv2-1.58113883008419*Gdrag_l[1]*rdv2; 
  out[24] += 1.58113883008419*Gdrag_r[2]*rdv2-1.58113883008419*Gdrag_l[2]*rdv2; 
  out[25] += 0.7071067811865475*Gdrag_r[13]*rdv2-0.7071067811865475*Gdrag_l[13]*rdv2; 
  out[26] += 0.7071067811865475*Gdrag_r[14]*rdv2-0.7071067811865475*Gdrag_l[14]*rdv2; 
  out[27] += 1.58113883008419*Gdrag_r[3]*rdv2-1.58113883008419*Gdrag_l[3]*rdv2; 
  out[28] += 0.7071067811865475*Gdrag_r[15]*rdv2-0.7071067811865475*Gdrag_l[15]*rdv2; 
  out[29] += 0.7071067811865475*Gdrag_r[16]*rdv2-0.7071067811865475*Gdrag_l[16]*rdv2; 
  out[30] += 1.224744871391589*Gdrag_r[9]*rdv2+1.224744871391589*Gdrag_l[9]*rdv2; 
  out[31] += 1.224744871391589*Gdrag_r[10]*rdv2+1.224744871391589*Gdrag_l[10]*rdv2; 
  out[32] += 1.224744871391589*Gdrag_r[11]*rdv2+1.224744871391589*Gdrag_l[11]*rdv2; 
  out[33] += 1.224744871391589*Gdrag_r[12]*rdv2+1.224744871391589*Gdrag_l[12]*rdv2; 
  out[34] += 1.58113883008419*Gdrag_r[4]*rdv2-1.58113883008419*Gdrag_l[4]*rdv2; 
  out[35] += 0.7071067811865475*Gdrag_r[17]*rdv2-0.7071067811865475*Gdrag_l[17]*rdv2; 
  out[36] += 0.7071067811865475*Gdrag_r[18]*rdv2-0.7071067811865475*Gdrag_l[18]*rdv2; 
  out[37] += 1.224744871391589*Gdrag_r[13]*rdv2+1.224744871391589*Gdrag_l[13]*rdv2; 
  out[38] += 1.224744871391589*Gdrag_r[14]*rdv2+1.224744871391589*Gdrag_l[14]*rdv2; 
  out[39] += 1.58113883008419*Gdrag_r[5]*rdv2-1.58113883008419*Gdrag_l[5]*rdv2; 
  out[40] += 1.58113883008419*Gdrag_r[6]*rdv2-1.58113883008419*Gdrag_l[6]*rdv2; 
  out[41] += 0.7071067811865475*Gdrag_r[19]*rdv2-0.7071067811865475*Gdrag_l[19]*rdv2; 
  out[42] += 1.224744871391589*Gdrag_r[15]*rdv2+1.224744871391589*Gdrag_l[15]*rdv2; 
  out[43] += 1.224744871391589*Gdrag_r[16]*rdv2+1.224744871391589*Gdrag_l[16]*rdv2; 
  out[44] += 1.224744871391589*Gdrag_r[17]*rdv2+1.224744871391589*Gdrag_l[17]*rdv2; 
  out[45] += 1.224744871391589*Gdrag_r[18]*rdv2+1.224744871391589*Gdrag_l[18]*rdv2; 
  out[46] += 1.58113883008419*Gdrag_r[10]*rdv2-1.58113883008419*Gdrag_l[10]*rdv2; 
  out[47] += 1.224744871391589*Gdrag_r[19]*rdv2+1.224744871391589*Gdrag_l[19]*rdv2; 
} 
