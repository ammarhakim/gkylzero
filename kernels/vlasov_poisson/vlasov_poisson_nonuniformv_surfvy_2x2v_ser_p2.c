#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_nonuniformv_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *vcoord, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // vcoord:    Discrete (DG) velocity coordinate.
  // field:     potential (scaled by appropriate factors).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double rdv2 = 2/dxv[3]; 
  const double *phi = &field[0]; 
  const double rdx2 = 2/dxv[0]; 
  const double rdy2 = 2/dxv[1]; 
  double alpha[20] = {0.0}; 

  alpha[0] = -2.449489742783178*phi[2]*rdy2; 
  alpha[1] = -2.449489742783178*phi[3]*rdy2; 
  alpha[2] = -5.477225575051662*phi[5]*rdy2; 
  alpha[4] = -5.477225575051662*phi[7]*rdy2; 
  alpha[7] = -2.449489742783178*phi[6]*rdy2; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 

  if (0.3162277660168379*alpha[7]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx4_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = ser_4x_p2_surfx4_eval_quad_node_0_r(fc); 
    fUpwindQuad_l[1] = ser_4x_p2_surfx4_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = ser_4x_p2_surfx4_eval_quad_node_1_r(fc); 
    fUpwindQuad_l[2] = ser_4x_p2_surfx4_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = ser_4x_p2_surfx4_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx4_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = ser_4x_p2_surfx4_eval_quad_node_0_l(fr); 
    fUpwindQuad_l[1] = ser_4x_p2_surfx4_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = ser_4x_p2_surfx4_eval_quad_node_1_l(fr); 
    fUpwindQuad_l[2] = ser_4x_p2_surfx4_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = ser_4x_p2_surfx4_eval_quad_node_2_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx4_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = ser_4x_p2_surfx4_eval_quad_node_3_r(fc); 
    fUpwindQuad_l[4] = ser_4x_p2_surfx4_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = ser_4x_p2_surfx4_eval_quad_node_4_r(fc); 
    fUpwindQuad_l[5] = ser_4x_p2_surfx4_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = ser_4x_p2_surfx4_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx4_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = ser_4x_p2_surfx4_eval_quad_node_3_l(fr); 
    fUpwindQuad_l[4] = ser_4x_p2_surfx4_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = ser_4x_p2_surfx4_eval_quad_node_4_l(fr); 
    fUpwindQuad_l[5] = ser_4x_p2_surfx4_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = ser_4x_p2_surfx4_eval_quad_node_5_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]+0.6363961030678926*alpha[4]-0.4743416490252568*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx4_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = ser_4x_p2_surfx4_eval_quad_node_6_r(fc); 
    fUpwindQuad_l[7] = ser_4x_p2_surfx4_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = ser_4x_p2_surfx4_eval_quad_node_7_r(fc); 
    fUpwindQuad_l[8] = ser_4x_p2_surfx4_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = ser_4x_p2_surfx4_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx4_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = ser_4x_p2_surfx4_eval_quad_node_6_l(fr); 
    fUpwindQuad_l[7] = ser_4x_p2_surfx4_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = ser_4x_p2_surfx4_eval_quad_node_7_l(fr); 
    fUpwindQuad_l[8] = ser_4x_p2_surfx4_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = ser_4x_p2_surfx4_eval_quad_node_8_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx4_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx4_eval_quad_node_9_r(fc); 
    fUpwindQuad_l[10] = ser_4x_p2_surfx4_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx4_eval_quad_node_10_r(fc); 
    fUpwindQuad_l[11] = ser_4x_p2_surfx4_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx4_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx4_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = ser_4x_p2_surfx4_eval_quad_node_9_l(fr); 
    fUpwindQuad_l[10] = ser_4x_p2_surfx4_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = ser_4x_p2_surfx4_eval_quad_node_10_l(fr); 
    fUpwindQuad_l[11] = ser_4x_p2_surfx4_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = ser_4x_p2_surfx4_eval_quad_node_11_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx4_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx4_eval_quad_node_12_r(fc); 
    fUpwindQuad_l[13] = ser_4x_p2_surfx4_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx4_eval_quad_node_13_r(fc); 
    fUpwindQuad_l[14] = ser_4x_p2_surfx4_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx4_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx4_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = ser_4x_p2_surfx4_eval_quad_node_12_l(fr); 
    fUpwindQuad_l[13] = ser_4x_p2_surfx4_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = ser_4x_p2_surfx4_eval_quad_node_13_l(fr); 
    fUpwindQuad_l[14] = ser_4x_p2_surfx4_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = ser_4x_p2_surfx4_eval_quad_node_14_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx4_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx4_eval_quad_node_15_r(fc); 
    fUpwindQuad_l[16] = ser_4x_p2_surfx4_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx4_eval_quad_node_16_r(fc); 
    fUpwindQuad_l[17] = ser_4x_p2_surfx4_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx4_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx4_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = ser_4x_p2_surfx4_eval_quad_node_15_l(fr); 
    fUpwindQuad_l[16] = ser_4x_p2_surfx4_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = ser_4x_p2_surfx4_eval_quad_node_16_l(fr); 
    fUpwindQuad_l[17] = ser_4x_p2_surfx4_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = ser_4x_p2_surfx4_eval_quad_node_17_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx4_eval_quad_node_18_r(fl); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx4_eval_quad_node_18_r(fc); 
    fUpwindQuad_l[19] = ser_4x_p2_surfx4_eval_quad_node_19_r(fl); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx4_eval_quad_node_19_r(fc); 
    fUpwindQuad_l[20] = ser_4x_p2_surfx4_eval_quad_node_20_r(fl); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx4_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx4_eval_quad_node_18_l(fc); 
    fUpwindQuad_r[18] = ser_4x_p2_surfx4_eval_quad_node_18_l(fr); 
    fUpwindQuad_l[19] = ser_4x_p2_surfx4_eval_quad_node_19_l(fc); 
    fUpwindQuad_r[19] = ser_4x_p2_surfx4_eval_quad_node_19_l(fr); 
    fUpwindQuad_l[20] = ser_4x_p2_surfx4_eval_quad_node_20_l(fc); 
    fUpwindQuad_r[20] = ser_4x_p2_surfx4_eval_quad_node_20_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx4_eval_quad_node_21_r(fl); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx4_eval_quad_node_21_r(fc); 
    fUpwindQuad_l[22] = ser_4x_p2_surfx4_eval_quad_node_22_r(fl); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx4_eval_quad_node_22_r(fc); 
    fUpwindQuad_l[23] = ser_4x_p2_surfx4_eval_quad_node_23_r(fl); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx4_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx4_eval_quad_node_21_l(fc); 
    fUpwindQuad_r[21] = ser_4x_p2_surfx4_eval_quad_node_21_l(fr); 
    fUpwindQuad_l[22] = ser_4x_p2_surfx4_eval_quad_node_22_l(fc); 
    fUpwindQuad_r[22] = ser_4x_p2_surfx4_eval_quad_node_22_l(fr); 
    fUpwindQuad_l[23] = ser_4x_p2_surfx4_eval_quad_node_23_l(fc); 
    fUpwindQuad_r[23] = ser_4x_p2_surfx4_eval_quad_node_23_l(fr); 
  } 
  if (0.3162277660168379*alpha[7]-0.6363961030678926*alpha[4]+0.4743416490252568*alpha[2]-0.4743416490252568*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx4_eval_quad_node_24_r(fl); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx4_eval_quad_node_24_r(fc); 
    fUpwindQuad_l[25] = ser_4x_p2_surfx4_eval_quad_node_25_r(fl); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx4_eval_quad_node_25_r(fc); 
    fUpwindQuad_l[26] = ser_4x_p2_surfx4_eval_quad_node_26_r(fl); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx4_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx4_eval_quad_node_24_l(fc); 
    fUpwindQuad_r[24] = ser_4x_p2_surfx4_eval_quad_node_24_l(fr); 
    fUpwindQuad_l[25] = ser_4x_p2_surfx4_eval_quad_node_25_l(fc); 
    fUpwindQuad_r[25] = ser_4x_p2_surfx4_eval_quad_node_25_l(fr); 
    fUpwindQuad_l[26] = ser_4x_p2_surfx4_eval_quad_node_26_l(fc); 
    fUpwindQuad_r[26] = ser_4x_p2_surfx4_eval_quad_node_26_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[4]*fUpwind_l[4]+0.3535533905932737*alpha[2]*fUpwind_l[2]+0.3535533905932737*alpha[1]*fUpwind_l[1]+0.3535533905932737*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3162277660168379*alpha[4]*fUpwind_l[11]+0.3162277660168379*alpha[1]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[1]*alpha[7]+0.3535533905932737*alpha[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.3162277660168379*alpha[4]*fUpwind_l[12]+0.3535533905932737*alpha[7]*fUpwind_l[11]+0.3162277660168379*alpha[2]*fUpwind_l[8]+0.3535533905932737*alpha[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.3535533905932737*alpha[7]*fUpwind_l[13]+0.3535533905932737*alpha[4]*fUpwind_l[10]+0.3535533905932737*alpha[2]*fUpwind_l[6]+0.3535533905932737*alpha[1]*fUpwind_l[5]+0.3535533905932737*alpha[0]*fUpwind_l[3]; 
  Ghat_l[4] = 0.3162277660168379*alpha[2]*fUpwind_l[12]+0.3162277660168379*alpha[1]*fUpwind_l[11]+0.3162277660168379*alpha[4]*fUpwind_l[8]+0.3162277660168379*alpha[4]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[4]*alpha[7]+0.3535533905932737*alpha[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alpha[2]; 
  Ghat_l[5] = 0.3162277660168379*alpha[4]*fUpwind_l[17]+0.3162277660168379*alpha[1]*fUpwind_l[13]+0.3535533905932737*alpha[2]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_l[6]+0.3535533905932737*alpha[0]*fUpwind_l[5]+0.3535533905932737*alpha[1]*fUpwind_l[3]; 
  Ghat_l[6] = 0.3162277660168379*alpha[4]*fUpwind_l[18]+0.3535533905932737*alpha[7]*fUpwind_l[17]+0.3162277660168379*alpha[2]*fUpwind_l[14]+0.3535533905932737*alpha[1]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[6]+0.3535533905932737*alpha[4]*fUpwind_l[5]+0.3535533905932737*alpha[2]*fUpwind_l[3]; 
  Ghat_l[7] = 0.3535533905932737*alpha[2]*fUpwind_l[11]+0.2258769757263128*alpha[7]*fUpwind_l[7]+0.3535533905932737*alpha[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[1]*fUpwind_l[1]; 
  Ghat_l[8] = 0.3535533905932737*alpha[1]*fUpwind_l[12]+0.3535533905932737*alpha[0]*fUpwind_l[8]+0.3162277660168379*alpha[4]*fUpwind_l[4]+0.3162277660168379*alpha[2]*fUpwind_l[2]; 
  Ghat_l[9] = 0.3535533905932737*alpha[4]*fUpwind_l[19]+0.3535533905932737*alpha[2]*fUpwind_l[16]+0.3535533905932737*alpha[1]*fUpwind_l[15]+0.3535533905932737*alpha[0]*fUpwind_l[9]; 
  Ghat_l[10] = 0.3162277660168379*alpha[2]*fUpwind_l[18]+0.3162277660168379*alpha[1]*fUpwind_l[17]+0.3162277660168379*alpha[4]*fUpwind_l[14]+0.3162277660168379*alpha[4]*fUpwind_l[13]+0.3162277660168379*alpha[7]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[10]+0.3535533905932737*alpha[1]*fUpwind_l[6]+0.3535533905932737*alpha[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alpha[4]; 
  Ghat_l[11] = 0.2828427124746191*alpha[4]*fUpwind_l[12]+0.2258769757263128*alpha[7]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[11]+0.3535533905932737*alpha[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alpha[7]+0.3162277660168379*alpha[1]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[1]*alpha[4]; 
  Ghat_l[12] = 0.3162277660168379*alpha[7]*fUpwind_l[12]+0.3535533905932737*alpha[0]*fUpwind_l[12]+0.2828427124746191*alpha[4]*fUpwind_l[11]+0.3535533905932737*alpha[1]*fUpwind_l[8]+0.3162277660168379*alpha[2]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[2]*alpha[4]; 
  Ghat_l[13] = 0.3535533905932737*alpha[2]*fUpwind_l[17]+0.2258769757263128*alpha[7]*fUpwind_l[13]+0.3535533905932737*alpha[0]*fUpwind_l[13]+0.3162277660168379*alpha[4]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[3]*alpha[7]+0.3162277660168379*alpha[1]*fUpwind_l[5]; 
  Ghat_l[14] = 0.3535533905932737*alpha[1]*fUpwind_l[18]+0.3535533905932737*alpha[0]*fUpwind_l[14]+0.3162277660168379*alpha[4]*fUpwind_l[10]+0.3162277660168379*alpha[2]*fUpwind_l[6]; 
  Ghat_l[15] = 0.3535533905932737*alpha[2]*fUpwind_l[19]+0.3535533905932737*alpha[4]*fUpwind_l[16]+0.3162277660168379*alpha[7]*fUpwind_l[15]+0.3535533905932737*alpha[0]*fUpwind_l[15]+0.3535533905932737*alpha[1]*fUpwind_l[9]; 
  Ghat_l[16] = 0.3535533905932737*alpha[1]*fUpwind_l[19]+0.3535533905932737*alpha[0]*fUpwind_l[16]+0.3535533905932737*alpha[4]*fUpwind_l[15]+0.3535533905932737*alpha[2]*fUpwind_l[9]; 
  Ghat_l[17] = 0.2828427124746191*alpha[4]*fUpwind_l[18]+0.2258769757263128*alpha[7]*fUpwind_l[17]+0.3535533905932737*alpha[0]*fUpwind_l[17]+0.3535533905932737*alpha[2]*fUpwind_l[13]+0.3162277660168379*alpha[1]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_l[5]; 
  Ghat_l[18] = 0.3162277660168379*alpha[7]*fUpwind_l[18]+0.3535533905932737*alpha[0]*fUpwind_l[18]+0.2828427124746191*alpha[4]*fUpwind_l[17]+0.3535533905932737*alpha[1]*fUpwind_l[14]+0.3162277660168379*alpha[2]*fUpwind_l[10]+0.3162277660168379*alpha[4]*fUpwind_l[6]; 
  Ghat_l[19] = 0.3162277660168379*alpha[7]*fUpwind_l[19]+0.3535533905932737*alpha[0]*fUpwind_l[19]+0.3535533905932737*alpha[1]*fUpwind_l[16]+0.3535533905932737*alpha[2]*fUpwind_l[15]+0.3535533905932737*alpha[4]*fUpwind_l[9]; 

  Ghat_r[0] = 0.3535533905932737*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[4]*fUpwind_r[4]+0.3535533905932737*alpha[2]*fUpwind_r[2]+0.3535533905932737*alpha[1]*fUpwind_r[1]+0.3535533905932737*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3162277660168379*alpha[4]*fUpwind_r[11]+0.3162277660168379*alpha[1]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[1]*alpha[7]+0.3535533905932737*alpha[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.3162277660168379*alpha[4]*fUpwind_r[12]+0.3535533905932737*alpha[7]*fUpwind_r[11]+0.3162277660168379*alpha[2]*fUpwind_r[8]+0.3535533905932737*alpha[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alpha[4]+0.3535533905932737*alpha[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.3535533905932737*alpha[7]*fUpwind_r[13]+0.3535533905932737*alpha[4]*fUpwind_r[10]+0.3535533905932737*alpha[2]*fUpwind_r[6]+0.3535533905932737*alpha[1]*fUpwind_r[5]+0.3535533905932737*alpha[0]*fUpwind_r[3]; 
  Ghat_r[4] = 0.3162277660168379*alpha[2]*fUpwind_r[12]+0.3162277660168379*alpha[1]*fUpwind_r[11]+0.3162277660168379*alpha[4]*fUpwind_r[8]+0.3162277660168379*alpha[4]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[4]*alpha[7]+0.3535533905932737*alpha[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alpha[4]+0.3535533905932737*alpha[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alpha[2]; 
  Ghat_r[5] = 0.3162277660168379*alpha[4]*fUpwind_r[17]+0.3162277660168379*alpha[1]*fUpwind_r[13]+0.3535533905932737*alpha[2]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[5]*alpha[7]+0.3535533905932737*alpha[4]*fUpwind_r[6]+0.3535533905932737*alpha[0]*fUpwind_r[5]+0.3535533905932737*alpha[1]*fUpwind_r[3]; 
  Ghat_r[6] = 0.3162277660168379*alpha[4]*fUpwind_r[18]+0.3535533905932737*alpha[7]*fUpwind_r[17]+0.3162277660168379*alpha[2]*fUpwind_r[14]+0.3535533905932737*alpha[1]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[6]+0.3535533905932737*alpha[4]*fUpwind_r[5]+0.3535533905932737*alpha[2]*fUpwind_r[3]; 
  Ghat_r[7] = 0.3535533905932737*alpha[2]*fUpwind_r[11]+0.2258769757263128*alpha[7]*fUpwind_r[7]+0.3535533905932737*alpha[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[1]*fUpwind_r[1]; 
  Ghat_r[8] = 0.3535533905932737*alpha[1]*fUpwind_r[12]+0.3535533905932737*alpha[0]*fUpwind_r[8]+0.3162277660168379*alpha[4]*fUpwind_r[4]+0.3162277660168379*alpha[2]*fUpwind_r[2]; 
  Ghat_r[9] = 0.3535533905932737*alpha[4]*fUpwind_r[19]+0.3535533905932737*alpha[2]*fUpwind_r[16]+0.3535533905932737*alpha[1]*fUpwind_r[15]+0.3535533905932737*alpha[0]*fUpwind_r[9]; 
  Ghat_r[10] = 0.3162277660168379*alpha[2]*fUpwind_r[18]+0.3162277660168379*alpha[1]*fUpwind_r[17]+0.3162277660168379*alpha[4]*fUpwind_r[14]+0.3162277660168379*alpha[4]*fUpwind_r[13]+0.3162277660168379*alpha[7]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[10]+0.3535533905932737*alpha[1]*fUpwind_r[6]+0.3535533905932737*alpha[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alpha[4]; 
  Ghat_r[11] = 0.2828427124746191*alpha[4]*fUpwind_r[12]+0.2258769757263128*alpha[7]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[11]+0.3535533905932737*alpha[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alpha[7]+0.3162277660168379*alpha[1]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[1]*alpha[4]; 
  Ghat_r[12] = 0.3162277660168379*alpha[7]*fUpwind_r[12]+0.3535533905932737*alpha[0]*fUpwind_r[12]+0.2828427124746191*alpha[4]*fUpwind_r[11]+0.3535533905932737*alpha[1]*fUpwind_r[8]+0.3162277660168379*alpha[2]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[2]*alpha[4]; 
  Ghat_r[13] = 0.3535533905932737*alpha[2]*fUpwind_r[17]+0.2258769757263128*alpha[7]*fUpwind_r[13]+0.3535533905932737*alpha[0]*fUpwind_r[13]+0.3162277660168379*alpha[4]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[3]*alpha[7]+0.3162277660168379*alpha[1]*fUpwind_r[5]; 
  Ghat_r[14] = 0.3535533905932737*alpha[1]*fUpwind_r[18]+0.3535533905932737*alpha[0]*fUpwind_r[14]+0.3162277660168379*alpha[4]*fUpwind_r[10]+0.3162277660168379*alpha[2]*fUpwind_r[6]; 
  Ghat_r[15] = 0.3535533905932737*alpha[2]*fUpwind_r[19]+0.3535533905932737*alpha[4]*fUpwind_r[16]+0.3162277660168379*alpha[7]*fUpwind_r[15]+0.3535533905932737*alpha[0]*fUpwind_r[15]+0.3535533905932737*alpha[1]*fUpwind_r[9]; 
  Ghat_r[16] = 0.3535533905932737*alpha[1]*fUpwind_r[19]+0.3535533905932737*alpha[0]*fUpwind_r[16]+0.3535533905932737*alpha[4]*fUpwind_r[15]+0.3535533905932737*alpha[2]*fUpwind_r[9]; 
  Ghat_r[17] = 0.2828427124746191*alpha[4]*fUpwind_r[18]+0.2258769757263128*alpha[7]*fUpwind_r[17]+0.3535533905932737*alpha[0]*fUpwind_r[17]+0.3535533905932737*alpha[2]*fUpwind_r[13]+0.3162277660168379*alpha[1]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[6]*alpha[7]+0.3162277660168379*alpha[4]*fUpwind_r[5]; 
  Ghat_r[18] = 0.3162277660168379*alpha[7]*fUpwind_r[18]+0.3535533905932737*alpha[0]*fUpwind_r[18]+0.2828427124746191*alpha[4]*fUpwind_r[17]+0.3535533905932737*alpha[1]*fUpwind_r[14]+0.3162277660168379*alpha[2]*fUpwind_r[10]+0.3162277660168379*alpha[4]*fUpwind_r[6]; 
  Ghat_r[19] = 0.3162277660168379*alpha[7]*fUpwind_r[19]+0.3535533905932737*alpha[0]*fUpwind_r[19]+0.3535533905932737*alpha[1]*fUpwind_r[16]+0.3535533905932737*alpha[2]*fUpwind_r[15]+0.3535533905932737*alpha[4]*fUpwind_r[9]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*rdv2; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*rdv2; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*rdv2; 
  out[4] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*rdv2; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*rdv2; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*rdv2; 
  out[8] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[9] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdv2; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*rdv2; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*rdv2; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*rdv2; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*rdv2; 
  out[14] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*rdv2; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*rdv2; 
  out[16] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*rdv2; 
  out[17] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*rdv2; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*rdv2; 
  out[19] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*rdv2; 
  out[20] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*rdv2; 
  out[21] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*rdv2; 
  out[22] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*rdv2; 
  out[23] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*rdv2; 
  out[24] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*rdv2; 
  out[25] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*rdv2; 
  out[26] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*rdv2; 
  out[27] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*rdv2; 
  out[28] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*rdv2; 
  out[29] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*rdv2; 
  out[30] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*rdv2; 
  out[31] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*rdv2; 
  out[32] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*rdv2; 
  out[33] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*rdv2; 
  out[34] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*rdv2; 
  out[35] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*rdv2; 
  out[36] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*rdv2; 
  out[37] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*rdv2; 
  out[38] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*rdv2; 
  out[39] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*rdv2; 
  out[40] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*rdv2; 
  out[41] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*rdv2; 
  out[42] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*rdv2; 
  out[43] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*rdv2; 
  out[44] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*rdv2; 
  out[45] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*rdv2; 
  out[46] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*rdv2; 
  out[47] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*rdv2; 

  return 0.;

} 
