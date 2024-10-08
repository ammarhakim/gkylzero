#include <gkyl_vlasov_sr_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx3_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_sr_surfvx_2x2v_ser_p2(const double *w, const double *dxv, const double *gamma, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // gamma:     Particle Lorentz boost factor sqrt(1 + p^2).
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double *E0 = &qmem[0]; 
  double p1_over_gamma_l[3] = {0.0}; 
  double p1_over_gamma_r[3] = {0.0}; 
  p1_over_gamma_l[0] = 2.738612787525831*gamma[6]*dv11-2.121320343559642*gamma[3]*dv11+1.224744871391589*gamma[2]*dv11; 
  p1_over_gamma_l[1] = 2.738612787525831*gamma[5]*dv11-4.743416490252569*gamma[7]*dv11; 
  p1_over_gamma_r[0] = 2.738612787525831*gamma[6]*dv11+2.121320343559642*gamma[3]*dv11+1.224744871391589*gamma[2]*dv11; 
  p1_over_gamma_r[1] = 4.743416490252569*gamma[7]*dv11+2.738612787525831*gamma[5]*dv11; 

  const double *B2 = &qmem[40]; 

  double alpha_l[20] = {0.0}; 
  double alpha_r[20] = {0.0}; 

  alpha_l[0] = B2[0]*p1_over_gamma_l[0]+1.414213562373095*E0[0]; 
  alpha_l[1] = 1.414213562373095*E0[1]+p1_over_gamma_l[0]*B2[1]; 
  alpha_l[2] = 1.414213562373095*E0[2]+p1_over_gamma_l[0]*B2[2]; 
  alpha_l[3] = B2[0]*p1_over_gamma_l[1]; 
  alpha_l[4] = 1.414213562373095*E0[3]+p1_over_gamma_l[0]*B2[3]; 
  alpha_l[5] = B2[1]*p1_over_gamma_l[1]; 
  alpha_l[6] = p1_over_gamma_l[1]*B2[2]; 
  alpha_l[7] = 1.414213562373095*E0[4]+p1_over_gamma_l[0]*B2[4]; 
  alpha_l[8] = 1.414213562373095*E0[5]+p1_over_gamma_l[0]*B2[5]; 
  alpha_l[10] = p1_over_gamma_l[1]*B2[3]; 
  alpha_l[11] = 1.414213562373095*E0[6]+p1_over_gamma_l[0]*B2[6]; 
  alpha_l[12] = 1.414213562373095*E0[7]+p1_over_gamma_l[0]*B2[7]; 
  alpha_l[13] = 1.0*p1_over_gamma_l[1]*B2[4]; 
  alpha_l[14] = 1.0*p1_over_gamma_l[1]*B2[5]; 
  alpha_l[17] = 1.0*p1_over_gamma_l[1]*B2[6]; 
  alpha_l[18] = 1.0*p1_over_gamma_l[1]*B2[7]; 

  alpha_r[0] = B2[0]*p1_over_gamma_r[0]+1.414213562373095*E0[0]; 
  alpha_r[1] = 1.414213562373095*E0[1]+p1_over_gamma_r[0]*B2[1]; 
  alpha_r[2] = 1.414213562373095*E0[2]+p1_over_gamma_r[0]*B2[2]; 
  alpha_r[3] = B2[0]*p1_over_gamma_r[1]; 
  alpha_r[4] = 1.414213562373095*E0[3]+p1_over_gamma_r[0]*B2[3]; 
  alpha_r[5] = B2[1]*p1_over_gamma_r[1]; 
  alpha_r[6] = p1_over_gamma_r[1]*B2[2]; 
  alpha_r[7] = 1.414213562373095*E0[4]+p1_over_gamma_r[0]*B2[4]; 
  alpha_r[8] = 1.414213562373095*E0[5]+p1_over_gamma_r[0]*B2[5]; 
  alpha_r[10] = p1_over_gamma_r[1]*B2[3]; 
  alpha_r[11] = 1.414213562373095*E0[6]+p1_over_gamma_r[0]*B2[6]; 
  alpha_r[12] = 1.414213562373095*E0[7]+p1_over_gamma_r[0]*B2[7]; 
  alpha_r[13] = 1.0*p1_over_gamma_r[1]*B2[4]; 
  alpha_r[14] = 1.0*p1_over_gamma_r[1]*B2[5]; 
  alpha_r[17] = 1.0*p1_over_gamma_r[1]*B2[6]; 
  alpha_r[18] = 1.0*p1_over_gamma_r[1]*B2[7]; 

  double fUpwindQuad_l[27] = {0.0};
  double fUpwindQuad_r[27] = {0.0};
  double fUpwind_l[20] = {0.0};;
  double fUpwind_r[20] = {0.0};
  double Ghat_l[20] = {0.0}; 
  double Ghat_r[20] = {0.0}; 

  if (0.5692099788303082*(alpha_l[18]+alpha_l[17])-0.4242640687119285*(alpha_l[14]+alpha_l[13])-0.4242640687119281*alpha_l[12]-0.4242640687119285*alpha_l[11]-0.853814968245462*alpha_l[10]+0.3162277660168379*(alpha_l[8]+alpha_l[7])+0.6363961030678926*(alpha_l[6]+alpha_l[5]+alpha_l[4])-0.4743416490252568*(alpha_l[3]+alpha_l[2]+alpha_l[1])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx3_eval_quad_node_0_r(fl); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p2_surfx3_eval_quad_node_0_l(fc); 
  } 
  if (0.5692099788303082*(alpha_r[18]+alpha_r[17])-0.4242640687119285*(alpha_r[14]+alpha_r[13])-0.4242640687119281*alpha_r[12]-0.4242640687119285*alpha_r[11]-0.853814968245462*alpha_r[10]+0.3162277660168379*(alpha_r[8]+alpha_r[7])+0.6363961030678926*(alpha_r[6]+alpha_r[5]+alpha_r[4])-0.4743416490252568*(alpha_r[3]+alpha_r[2]+alpha_r[1])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_4x_p2_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_r[0] = ser_4x_p2_surfx3_eval_quad_node_0_l(fr); 
  } 
  if ((-0.4242640687119281*alpha_l[12])-0.4242640687119285*alpha_l[11]+0.3162277660168379*(alpha_l[8]+alpha_l[7])+0.6363961030678926*alpha_l[4]-0.4743416490252568*(alpha_l[2]+alpha_l[1])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx3_eval_quad_node_1_r(fl); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p2_surfx3_eval_quad_node_1_l(fc); 
  } 
  if ((-0.4242640687119281*alpha_r[12])-0.4242640687119285*alpha_r[11]+0.3162277660168379*(alpha_r[8]+alpha_r[7])+0.6363961030678926*alpha_r[4]-0.4743416490252568*(alpha_r[2]+alpha_r[1])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_4x_p2_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_r[1] = ser_4x_p2_surfx3_eval_quad_node_1_l(fr); 
  } 
  if ((-0.5692099788303082*(alpha_l[18]+alpha_l[17]))+0.4242640687119285*(alpha_l[14]+alpha_l[13])-0.4242640687119281*alpha_l[12]-0.4242640687119285*alpha_l[11]+0.853814968245462*alpha_l[10]+0.3162277660168379*(alpha_l[8]+alpha_l[7])-0.6363961030678926*(alpha_l[6]+alpha_l[5])+0.6363961030678926*alpha_l[4]+0.4743416490252568*alpha_l[3]-0.4743416490252568*(alpha_l[2]+alpha_l[1])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx3_eval_quad_node_2_r(fl); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p2_surfx3_eval_quad_node_2_l(fc); 
  } 
  if ((-0.5692099788303082*(alpha_r[18]+alpha_r[17]))+0.4242640687119285*(alpha_r[14]+alpha_r[13])-0.4242640687119281*alpha_r[12]-0.4242640687119285*alpha_r[11]+0.853814968245462*alpha_r[10]+0.3162277660168379*(alpha_r[8]+alpha_r[7])-0.6363961030678926*(alpha_r[6]+alpha_r[5])+0.6363961030678926*alpha_r[4]+0.4743416490252568*alpha_r[3]-0.4743416490252568*(alpha_r[2]+alpha_r[1])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_4x_p2_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_r[2] = ser_4x_p2_surfx3_eval_quad_node_2_l(fr); 
  } 
  if ((-0.711512473537885*alpha_l[18])+0.5303300858899104*alpha_l[14]-0.4242640687119285*alpha_l[13]+0.5303300858899104*alpha_l[12]-0.3952847075210473*alpha_l[8]+0.3162277660168379*alpha_l[7]+0.6363961030678926*alpha_l[5]-0.4743416490252568*(alpha_l[3]+alpha_l[1])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx3_eval_quad_node_3_r(fl); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p2_surfx3_eval_quad_node_3_l(fc); 
  } 
  if ((-0.711512473537885*alpha_r[18])+0.5303300858899104*alpha_r[14]-0.4242640687119285*alpha_r[13]+0.5303300858899104*alpha_r[12]-0.3952847075210473*alpha_r[8]+0.3162277660168379*alpha_r[7]+0.6363961030678926*alpha_r[5]-0.4743416490252568*(alpha_r[3]+alpha_r[1])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_4x_p2_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_r[3] = ser_4x_p2_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (0.5303300858899104*alpha_l[12]-0.3952847075210473*alpha_l[8]+0.3162277660168379*alpha_l[7]-0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx3_eval_quad_node_4_r(fl); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p2_surfx3_eval_quad_node_4_l(fc); 
  } 
  if (0.5303300858899104*alpha_r[12]-0.3952847075210473*alpha_r[8]+0.3162277660168379*alpha_r[7]-0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[4] = ser_4x_p2_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_r[4] = ser_4x_p2_surfx3_eval_quad_node_4_l(fr); 
  } 
  if (0.711512473537885*alpha_l[18]-0.5303300858899104*alpha_l[14]+0.4242640687119285*alpha_l[13]+0.5303300858899104*alpha_l[12]-0.3952847075210473*alpha_l[8]+0.3162277660168379*alpha_l[7]-0.6363961030678926*alpha_l[5]+0.4743416490252568*alpha_l[3]-0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx3_eval_quad_node_5_r(fl); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p2_surfx3_eval_quad_node_5_l(fc); 
  } 
  if (0.711512473537885*alpha_r[18]-0.5303300858899104*alpha_r[14]+0.4242640687119285*alpha_r[13]+0.5303300858899104*alpha_r[12]-0.3952847075210473*alpha_r[8]+0.3162277660168379*alpha_r[7]-0.6363961030678926*alpha_r[5]+0.4743416490252568*alpha_r[3]-0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[5] = ser_4x_p2_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_r[5] = ser_4x_p2_surfx3_eval_quad_node_5_l(fr); 
  } 
  if (0.5692099788303082*alpha_l[18]-0.5692099788303082*alpha_l[17]-0.4242640687119285*(alpha_l[14]+alpha_l[13])-0.4242640687119281*alpha_l[12]+0.4242640687119285*alpha_l[11]+0.853814968245462*alpha_l[10]+0.3162277660168379*(alpha_l[8]+alpha_l[7])-0.6363961030678926*alpha_l[6]+0.6363961030678926*alpha_l[5]-0.6363961030678926*alpha_l[4]-0.4743416490252568*alpha_l[3]+0.4743416490252568*alpha_l[2]-0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx3_eval_quad_node_6_r(fl); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p2_surfx3_eval_quad_node_6_l(fc); 
  } 
  if (0.5692099788303082*alpha_r[18]-0.5692099788303082*alpha_r[17]-0.4242640687119285*(alpha_r[14]+alpha_r[13])-0.4242640687119281*alpha_r[12]+0.4242640687119285*alpha_r[11]+0.853814968245462*alpha_r[10]+0.3162277660168379*(alpha_r[8]+alpha_r[7])-0.6363961030678926*alpha_r[6]+0.6363961030678926*alpha_r[5]-0.6363961030678926*alpha_r[4]-0.4743416490252568*alpha_r[3]+0.4743416490252568*alpha_r[2]-0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[6] = ser_4x_p2_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_r[6] = ser_4x_p2_surfx3_eval_quad_node_6_l(fr); 
  } 
  if ((-0.4242640687119281*alpha_l[12])+0.4242640687119285*alpha_l[11]+0.3162277660168379*(alpha_l[8]+alpha_l[7])-0.6363961030678926*alpha_l[4]+0.4743416490252568*alpha_l[2]-0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx3_eval_quad_node_7_r(fl); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p2_surfx3_eval_quad_node_7_l(fc); 
  } 
  if ((-0.4242640687119281*alpha_r[12])+0.4242640687119285*alpha_r[11]+0.3162277660168379*(alpha_r[8]+alpha_r[7])-0.6363961030678926*alpha_r[4]+0.4743416490252568*alpha_r[2]-0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[7] = ser_4x_p2_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_r[7] = ser_4x_p2_surfx3_eval_quad_node_7_l(fr); 
  } 
  if ((-0.5692099788303082*alpha_l[18])+0.5692099788303082*alpha_l[17]+0.4242640687119285*(alpha_l[14]+alpha_l[13])-0.4242640687119281*alpha_l[12]+0.4242640687119285*alpha_l[11]-0.853814968245462*alpha_l[10]+0.3162277660168379*(alpha_l[8]+alpha_l[7])+0.6363961030678926*alpha_l[6]-0.6363961030678926*(alpha_l[5]+alpha_l[4])+0.4743416490252568*(alpha_l[3]+alpha_l[2])-0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx3_eval_quad_node_8_r(fl); 
  } else { 
    fUpwindQuad_l[8] = ser_4x_p2_surfx3_eval_quad_node_8_l(fc); 
  } 
  if ((-0.5692099788303082*alpha_r[18])+0.5692099788303082*alpha_r[17]+0.4242640687119285*(alpha_r[14]+alpha_r[13])-0.4242640687119281*alpha_r[12]+0.4242640687119285*alpha_r[11]-0.853814968245462*alpha_r[10]+0.3162277660168379*(alpha_r[8]+alpha_r[7])+0.6363961030678926*alpha_r[6]-0.6363961030678926*(alpha_r[5]+alpha_r[4])+0.4743416490252568*(alpha_r[3]+alpha_r[2])-0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[8] = ser_4x_p2_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_r[8] = ser_4x_p2_surfx3_eval_quad_node_8_l(fr); 
  } 
  if ((-0.711512473537885*alpha_l[17])-0.4242640687119285*alpha_l[14]+0.5303300858899104*(alpha_l[13]+alpha_l[11])+0.3162277660168379*alpha_l[8]-0.3952847075210473*alpha_l[7]+0.6363961030678926*alpha_l[6]-0.4743416490252568*(alpha_l[3]+alpha_l[2])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx3_eval_quad_node_9_r(fl); 
  } else { 
    fUpwindQuad_l[9] = ser_4x_p2_surfx3_eval_quad_node_9_l(fc); 
  } 
  if ((-0.711512473537885*alpha_r[17])-0.4242640687119285*alpha_r[14]+0.5303300858899104*(alpha_r[13]+alpha_r[11])+0.3162277660168379*alpha_r[8]-0.3952847075210473*alpha_r[7]+0.6363961030678926*alpha_r[6]-0.4743416490252568*(alpha_r[3]+alpha_r[2])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[9] = ser_4x_p2_surfx3_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_r[9] = ser_4x_p2_surfx3_eval_quad_node_9_l(fr); 
  } 
  if (0.5303300858899104*alpha_l[11]+0.3162277660168379*alpha_l[8]-0.3952847075210473*alpha_l[7]-0.4743416490252568*alpha_l[2]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx3_eval_quad_node_10_r(fl); 
  } else { 
    fUpwindQuad_l[10] = ser_4x_p2_surfx3_eval_quad_node_10_l(fc); 
  } 
  if (0.5303300858899104*alpha_r[11]+0.3162277660168379*alpha_r[8]-0.3952847075210473*alpha_r[7]-0.4743416490252568*alpha_r[2]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[10] = ser_4x_p2_surfx3_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_r[10] = ser_4x_p2_surfx3_eval_quad_node_10_l(fr); 
  } 
  if (0.711512473537885*alpha_l[17]+0.4242640687119285*alpha_l[14]-0.5303300858899104*alpha_l[13]+0.5303300858899104*alpha_l[11]+0.3162277660168379*alpha_l[8]-0.3952847075210473*alpha_l[7]-0.6363961030678926*alpha_l[6]+0.4743416490252568*alpha_l[3]-0.4743416490252568*alpha_l[2]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx3_eval_quad_node_11_r(fl); 
  } else { 
    fUpwindQuad_l[11] = ser_4x_p2_surfx3_eval_quad_node_11_l(fc); 
  } 
  if (0.711512473537885*alpha_r[17]+0.4242640687119285*alpha_r[14]-0.5303300858899104*alpha_r[13]+0.5303300858899104*alpha_r[11]+0.3162277660168379*alpha_r[8]-0.3952847075210473*alpha_r[7]-0.6363961030678926*alpha_r[6]+0.4743416490252568*alpha_r[3]-0.4743416490252568*alpha_r[2]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[11] = ser_4x_p2_surfx3_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_r[11] = ser_4x_p2_surfx3_eval_quad_node_11_l(fr); 
  } 
  if (0.5303300858899104*(alpha_l[14]+alpha_l[13])-0.3952847075210473*(alpha_l[8]+alpha_l[7])-0.4743416490252568*alpha_l[3]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx3_eval_quad_node_12_r(fl); 
  } else { 
    fUpwindQuad_l[12] = ser_4x_p2_surfx3_eval_quad_node_12_l(fc); 
  } 
  if (0.5303300858899104*(alpha_r[14]+alpha_r[13])-0.3952847075210473*(alpha_r[8]+alpha_r[7])-0.4743416490252568*alpha_r[3]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[12] = ser_4x_p2_surfx3_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_r[12] = ser_4x_p2_surfx3_eval_quad_node_12_l(fr); 
  } 
  if (0.3535533905932737*alpha_l[0]-0.3952847075210473*(alpha_l[8]+alpha_l[7]) > 0) { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx3_eval_quad_node_13_r(fl); 
  } else { 
    fUpwindQuad_l[13] = ser_4x_p2_surfx3_eval_quad_node_13_l(fc); 
  } 
  if (0.3535533905932737*alpha_r[0]-0.3952847075210473*(alpha_r[8]+alpha_r[7]) > 0) { 
    fUpwindQuad_r[13] = ser_4x_p2_surfx3_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_r[13] = ser_4x_p2_surfx3_eval_quad_node_13_l(fr); 
  } 
  if ((-0.5303300858899104*(alpha_l[14]+alpha_l[13]))-0.3952847075210473*(alpha_l[8]+alpha_l[7])+0.4743416490252568*alpha_l[3]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx3_eval_quad_node_14_r(fl); 
  } else { 
    fUpwindQuad_l[14] = ser_4x_p2_surfx3_eval_quad_node_14_l(fc); 
  } 
  if ((-0.5303300858899104*(alpha_r[14]+alpha_r[13]))-0.3952847075210473*(alpha_r[8]+alpha_r[7])+0.4743416490252568*alpha_r[3]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[14] = ser_4x_p2_surfx3_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_r[14] = ser_4x_p2_surfx3_eval_quad_node_14_l(fr); 
  } 
  if (0.711512473537885*alpha_l[17]-0.4242640687119285*alpha_l[14]+0.5303300858899104*alpha_l[13]-0.5303300858899104*alpha_l[11]+0.3162277660168379*alpha_l[8]-0.3952847075210473*alpha_l[7]-0.6363961030678926*alpha_l[6]-0.4743416490252568*alpha_l[3]+0.4743416490252568*alpha_l[2]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx3_eval_quad_node_15_r(fl); 
  } else { 
    fUpwindQuad_l[15] = ser_4x_p2_surfx3_eval_quad_node_15_l(fc); 
  } 
  if (0.711512473537885*alpha_r[17]-0.4242640687119285*alpha_r[14]+0.5303300858899104*alpha_r[13]-0.5303300858899104*alpha_r[11]+0.3162277660168379*alpha_r[8]-0.3952847075210473*alpha_r[7]-0.6363961030678926*alpha_r[6]-0.4743416490252568*alpha_r[3]+0.4743416490252568*alpha_r[2]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[15] = ser_4x_p2_surfx3_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_r[15] = ser_4x_p2_surfx3_eval_quad_node_15_l(fr); 
  } 
  if ((-0.5303300858899104*alpha_l[11])+0.3162277660168379*alpha_l[8]-0.3952847075210473*alpha_l[7]+0.4743416490252568*alpha_l[2]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx3_eval_quad_node_16_r(fl); 
  } else { 
    fUpwindQuad_l[16] = ser_4x_p2_surfx3_eval_quad_node_16_l(fc); 
  } 
  if ((-0.5303300858899104*alpha_r[11])+0.3162277660168379*alpha_r[8]-0.3952847075210473*alpha_r[7]+0.4743416490252568*alpha_r[2]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[16] = ser_4x_p2_surfx3_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_r[16] = ser_4x_p2_surfx3_eval_quad_node_16_l(fr); 
  } 
  if ((-0.711512473537885*alpha_l[17])+0.4242640687119285*alpha_l[14]-0.5303300858899104*(alpha_l[13]+alpha_l[11])+0.3162277660168379*alpha_l[8]-0.3952847075210473*alpha_l[7]+0.6363961030678926*alpha_l[6]+0.4743416490252568*(alpha_l[3]+alpha_l[2])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx3_eval_quad_node_17_r(fl); 
  } else { 
    fUpwindQuad_l[17] = ser_4x_p2_surfx3_eval_quad_node_17_l(fc); 
  } 
  if ((-0.711512473537885*alpha_r[17])+0.4242640687119285*alpha_r[14]-0.5303300858899104*(alpha_r[13]+alpha_r[11])+0.3162277660168379*alpha_r[8]-0.3952847075210473*alpha_r[7]+0.6363961030678926*alpha_r[6]+0.4743416490252568*(alpha_r[3]+alpha_r[2])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[17] = ser_4x_p2_surfx3_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_r[17] = ser_4x_p2_surfx3_eval_quad_node_17_l(fr); 
  } 
  if ((-0.5692099788303082*alpha_l[18])+0.5692099788303082*alpha_l[17]-0.4242640687119285*(alpha_l[14]+alpha_l[13])+0.4242640687119281*alpha_l[12]-0.4242640687119285*alpha_l[11]+0.853814968245462*alpha_l[10]+0.3162277660168379*(alpha_l[8]+alpha_l[7])+0.6363961030678926*alpha_l[6]-0.6363961030678926*(alpha_l[5]+alpha_l[4])-0.4743416490252568*(alpha_l[3]+alpha_l[2])+0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx3_eval_quad_node_18_r(fl); 
  } else { 
    fUpwindQuad_l[18] = ser_4x_p2_surfx3_eval_quad_node_18_l(fc); 
  } 
  if ((-0.5692099788303082*alpha_r[18])+0.5692099788303082*alpha_r[17]-0.4242640687119285*(alpha_r[14]+alpha_r[13])+0.4242640687119281*alpha_r[12]-0.4242640687119285*alpha_r[11]+0.853814968245462*alpha_r[10]+0.3162277660168379*(alpha_r[8]+alpha_r[7])+0.6363961030678926*alpha_r[6]-0.6363961030678926*(alpha_r[5]+alpha_r[4])-0.4743416490252568*(alpha_r[3]+alpha_r[2])+0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[18] = ser_4x_p2_surfx3_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_r[18] = ser_4x_p2_surfx3_eval_quad_node_18_l(fr); 
  } 
  if (0.4242640687119281*alpha_l[12]-0.4242640687119285*alpha_l[11]+0.3162277660168379*(alpha_l[8]+alpha_l[7])-0.6363961030678926*alpha_l[4]-0.4743416490252568*alpha_l[2]+0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx3_eval_quad_node_19_r(fl); 
  } else { 
    fUpwindQuad_l[19] = ser_4x_p2_surfx3_eval_quad_node_19_l(fc); 
  } 
  if (0.4242640687119281*alpha_r[12]-0.4242640687119285*alpha_r[11]+0.3162277660168379*(alpha_r[8]+alpha_r[7])-0.6363961030678926*alpha_r[4]-0.4743416490252568*alpha_r[2]+0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[19] = ser_4x_p2_surfx3_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_r[19] = ser_4x_p2_surfx3_eval_quad_node_19_l(fr); 
  } 
  if (0.5692099788303082*alpha_l[18]-0.5692099788303082*alpha_l[17]+0.4242640687119285*(alpha_l[14]+alpha_l[13])+0.4242640687119281*alpha_l[12]-0.4242640687119285*alpha_l[11]-0.853814968245462*alpha_l[10]+0.3162277660168379*(alpha_l[8]+alpha_l[7])-0.6363961030678926*alpha_l[6]+0.6363961030678926*alpha_l[5]-0.6363961030678926*alpha_l[4]+0.4743416490252568*alpha_l[3]-0.4743416490252568*alpha_l[2]+0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx3_eval_quad_node_20_r(fl); 
  } else { 
    fUpwindQuad_l[20] = ser_4x_p2_surfx3_eval_quad_node_20_l(fc); 
  } 
  if (0.5692099788303082*alpha_r[18]-0.5692099788303082*alpha_r[17]+0.4242640687119285*(alpha_r[14]+alpha_r[13])+0.4242640687119281*alpha_r[12]-0.4242640687119285*alpha_r[11]-0.853814968245462*alpha_r[10]+0.3162277660168379*(alpha_r[8]+alpha_r[7])-0.6363961030678926*alpha_r[6]+0.6363961030678926*alpha_r[5]-0.6363961030678926*alpha_r[4]+0.4743416490252568*alpha_r[3]-0.4743416490252568*alpha_r[2]+0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[20] = ser_4x_p2_surfx3_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_r[20] = ser_4x_p2_surfx3_eval_quad_node_20_l(fr); 
  } 
  if (0.711512473537885*alpha_l[18]+0.5303300858899104*alpha_l[14]-0.4242640687119285*alpha_l[13]-0.5303300858899104*alpha_l[12]-0.3952847075210473*alpha_l[8]+0.3162277660168379*alpha_l[7]-0.6363961030678926*alpha_l[5]-0.4743416490252568*alpha_l[3]+0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx3_eval_quad_node_21_r(fl); 
  } else { 
    fUpwindQuad_l[21] = ser_4x_p2_surfx3_eval_quad_node_21_l(fc); 
  } 
  if (0.711512473537885*alpha_r[18]+0.5303300858899104*alpha_r[14]-0.4242640687119285*alpha_r[13]-0.5303300858899104*alpha_r[12]-0.3952847075210473*alpha_r[8]+0.3162277660168379*alpha_r[7]-0.6363961030678926*alpha_r[5]-0.4743416490252568*alpha_r[3]+0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[21] = ser_4x_p2_surfx3_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_r[21] = ser_4x_p2_surfx3_eval_quad_node_21_l(fr); 
  } 
  if ((-0.5303300858899104*alpha_l[12])-0.3952847075210473*alpha_l[8]+0.3162277660168379*alpha_l[7]+0.4743416490252568*alpha_l[1]+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx3_eval_quad_node_22_r(fl); 
  } else { 
    fUpwindQuad_l[22] = ser_4x_p2_surfx3_eval_quad_node_22_l(fc); 
  } 
  if ((-0.5303300858899104*alpha_r[12])-0.3952847075210473*alpha_r[8]+0.3162277660168379*alpha_r[7]+0.4743416490252568*alpha_r[1]+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[22] = ser_4x_p2_surfx3_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_r[22] = ser_4x_p2_surfx3_eval_quad_node_22_l(fr); 
  } 
  if ((-0.711512473537885*alpha_l[18])-0.5303300858899104*alpha_l[14]+0.4242640687119285*alpha_l[13]-0.5303300858899104*alpha_l[12]-0.3952847075210473*alpha_l[8]+0.3162277660168379*alpha_l[7]+0.6363961030678926*alpha_l[5]+0.4743416490252568*(alpha_l[3]+alpha_l[1])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx3_eval_quad_node_23_r(fl); 
  } else { 
    fUpwindQuad_l[23] = ser_4x_p2_surfx3_eval_quad_node_23_l(fc); 
  } 
  if ((-0.711512473537885*alpha_r[18])-0.5303300858899104*alpha_r[14]+0.4242640687119285*alpha_r[13]-0.5303300858899104*alpha_r[12]-0.3952847075210473*alpha_r[8]+0.3162277660168379*alpha_r[7]+0.6363961030678926*alpha_r[5]+0.4743416490252568*(alpha_r[3]+alpha_r[1])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[23] = ser_4x_p2_surfx3_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_r[23] = ser_4x_p2_surfx3_eval_quad_node_23_l(fr); 
  } 
  if ((-0.5692099788303082*(alpha_l[18]+alpha_l[17]))-0.4242640687119285*(alpha_l[14]+alpha_l[13])+0.4242640687119281*alpha_l[12]+0.4242640687119285*alpha_l[11]-0.853814968245462*alpha_l[10]+0.3162277660168379*(alpha_l[8]+alpha_l[7])-0.6363961030678926*(alpha_l[6]+alpha_l[5])+0.6363961030678926*alpha_l[4]-0.4743416490252568*alpha_l[3]+0.4743416490252568*(alpha_l[2]+alpha_l[1])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx3_eval_quad_node_24_r(fl); 
  } else { 
    fUpwindQuad_l[24] = ser_4x_p2_surfx3_eval_quad_node_24_l(fc); 
  } 
  if ((-0.5692099788303082*(alpha_r[18]+alpha_r[17]))-0.4242640687119285*(alpha_r[14]+alpha_r[13])+0.4242640687119281*alpha_r[12]+0.4242640687119285*alpha_r[11]-0.853814968245462*alpha_r[10]+0.3162277660168379*(alpha_r[8]+alpha_r[7])-0.6363961030678926*(alpha_r[6]+alpha_r[5])+0.6363961030678926*alpha_r[4]-0.4743416490252568*alpha_r[3]+0.4743416490252568*(alpha_r[2]+alpha_r[1])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[24] = ser_4x_p2_surfx3_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_r[24] = ser_4x_p2_surfx3_eval_quad_node_24_l(fr); 
  } 
  if (0.4242640687119281*alpha_l[12]+0.4242640687119285*alpha_l[11]+0.3162277660168379*(alpha_l[8]+alpha_l[7])+0.6363961030678926*alpha_l[4]+0.4743416490252568*(alpha_l[2]+alpha_l[1])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx3_eval_quad_node_25_r(fl); 
  } else { 
    fUpwindQuad_l[25] = ser_4x_p2_surfx3_eval_quad_node_25_l(fc); 
  } 
  if (0.4242640687119281*alpha_r[12]+0.4242640687119285*alpha_r[11]+0.3162277660168379*(alpha_r[8]+alpha_r[7])+0.6363961030678926*alpha_r[4]+0.4743416490252568*(alpha_r[2]+alpha_r[1])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[25] = ser_4x_p2_surfx3_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_r[25] = ser_4x_p2_surfx3_eval_quad_node_25_l(fr); 
  } 
  if (0.5692099788303082*(alpha_l[18]+alpha_l[17])+0.4242640687119285*(alpha_l[14]+alpha_l[13])+0.4242640687119281*alpha_l[12]+0.4242640687119285*alpha_l[11]+0.853814968245462*alpha_l[10]+0.3162277660168379*(alpha_l[8]+alpha_l[7])+0.6363961030678926*(alpha_l[6]+alpha_l[5]+alpha_l[4])+0.4743416490252568*(alpha_l[3]+alpha_l[2]+alpha_l[1])+0.3535533905932737*alpha_l[0] > 0) { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx3_eval_quad_node_26_r(fl); 
  } else { 
    fUpwindQuad_l[26] = ser_4x_p2_surfx3_eval_quad_node_26_l(fc); 
  } 
  if (0.5692099788303082*(alpha_r[18]+alpha_r[17])+0.4242640687119285*(alpha_r[14]+alpha_r[13])+0.4242640687119281*alpha_r[12]+0.4242640687119285*alpha_r[11]+0.853814968245462*alpha_r[10]+0.3162277660168379*(alpha_r[8]+alpha_r[7])+0.6363961030678926*(alpha_r[6]+alpha_r[5]+alpha_r[4])+0.4743416490252568*(alpha_r[3]+alpha_r[2]+alpha_r[1])+0.3535533905932737*alpha_r[0] > 0) { 
    fUpwindQuad_r[26] = ser_4x_p2_surfx3_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_r[26] = ser_4x_p2_surfx3_eval_quad_node_26_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 
  Ghat_l[0] = 0.3535533905932737*alpha_l[18]*fUpwind_l[18]+0.3535533905932737*alpha_l[17]*fUpwind_l[17]+0.3535533905932737*alpha_l[14]*fUpwind_l[14]+0.3535533905932737*alpha_l[13]*fUpwind_l[13]+0.3535533905932737*alpha_l[12]*fUpwind_l[12]+0.3535533905932737*alpha_l[11]*fUpwind_l[11]+0.3535533905932737*alpha_l[10]*fUpwind_l[10]+0.3535533905932737*alpha_l[8]*fUpwind_l[8]+0.3535533905932737*alpha_l[7]*fUpwind_l[7]+0.3535533905932737*alpha_l[6]*fUpwind_l[6]+0.3535533905932737*alpha_l[5]*fUpwind_l[5]+0.3535533905932737*alpha_l[4]*fUpwind_l[4]+0.3535533905932737*alpha_l[3]*fUpwind_l[3]+0.3535533905932737*alpha_l[2]*fUpwind_l[2]+0.3535533905932737*alpha_l[1]*fUpwind_l[1]+0.3535533905932737*alpha_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alpha_l[14]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[14]*alpha_l[18]+0.3162277660168379*alpha_l[10]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[10]*alpha_l[17]+0.3162277660168379*alpha_l[5]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[5]*alpha_l[13]+0.3535533905932737*alpha_l[8]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[8]*alpha_l[12]+0.3162277660168379*alpha_l[4]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[4]*alpha_l[11]+0.3535533905932737*alpha_l[6]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[6]*alpha_l[10]+0.3162277660168379*alpha_l[1]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[1]*alpha_l[7]+0.3535533905932737*alpha_l[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alpha_l[5]+0.3535533905932737*alpha_l[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alpha_l[4]+0.3535533905932737*alpha_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alpha_l[1]; 
  Ghat_l[2] = 0.3162277660168379*alpha_l[10]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[10]*alpha_l[18]+0.3535533905932737*alpha_l[13]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[13]*alpha_l[17]+0.3162277660168379*alpha_l[6]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[6]*alpha_l[14]+0.3162277660168379*alpha_l[4]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[4]*alpha_l[12]+0.3535533905932737*alpha_l[7]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[7]*alpha_l[11]+0.3535533905932737*alpha_l[5]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[5]*alpha_l[10]+0.3162277660168379*alpha_l[2]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[2]*alpha_l[8]+0.3535533905932737*alpha_l[3]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[3]*alpha_l[6]+0.3535533905932737*alpha_l[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alpha_l[4]+0.3535533905932737*alpha_l[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha_l[2]; 
  Ghat_l[3] = 0.3162277660168379*alpha_l[10]*fUpwind_l[19]+0.3535533905932737*alpha_l[12]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[12]*alpha_l[18]+0.3535533905932737*alpha_l[11]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[11]*alpha_l[17]+0.3162277660168379*alpha_l[6]*fUpwind_l[16]+0.3162277660168379*alpha_l[5]*fUpwind_l[15]+0.3535533905932737*alpha_l[8]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[8]*alpha_l[14]+0.3535533905932737*alpha_l[7]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[7]*alpha_l[13]+0.3535533905932737*alpha_l[4]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[4]*alpha_l[10]+0.3162277660168379*alpha_l[3]*fUpwind_l[9]+0.3535533905932737*alpha_l[2]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[2]*alpha_l[6]+0.3535533905932737*alpha_l[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alpha_l[5]+0.3535533905932737*alpha_l[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alpha_l[3]; 
  Ghat_l[4] = 0.2828427124746191*alpha_l[17]*fUpwind_l[18]+0.3162277660168379*alpha_l[6]*fUpwind_l[18]+0.2828427124746191*fUpwind_l[17]*alpha_l[18]+0.3162277660168379*fUpwind_l[6]*alpha_l[18]+0.3162277660168379*alpha_l[5]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[5]*alpha_l[17]+0.3162277660168379*alpha_l[10]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[10]*alpha_l[14]+0.3162277660168379*alpha_l[10]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[10]*alpha_l[13]+0.2828427124746191*alpha_l[11]*fUpwind_l[12]+0.3162277660168379*alpha_l[2]*fUpwind_l[12]+0.2828427124746191*fUpwind_l[11]*alpha_l[12]+0.3162277660168379*fUpwind_l[2]*alpha_l[12]+0.3162277660168379*alpha_l[1]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[1]*alpha_l[11]+0.3535533905932737*alpha_l[3]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[3]*alpha_l[10]+0.3162277660168379*alpha_l[4]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[4]*alpha_l[8]+0.3162277660168379*alpha_l[4]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[4]*alpha_l[7]+0.3535533905932737*alpha_l[5]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[5]*alpha_l[6]+0.3535533905932737*alpha_l[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alpha_l[4]+0.3535533905932737*alpha_l[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alpha_l[2]; 
  Ghat_l[5] = 0.2828427124746191*alpha_l[17]*fUpwind_l[19]+0.3162277660168379*alpha_l[6]*fUpwind_l[19]+0.3535533905932737*alpha_l[8]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[8]*alpha_l[18]+0.3162277660168379*alpha_l[4]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[4]*alpha_l[17]+0.3162277660168379*alpha_l[10]*fUpwind_l[16]+0.2828427124746191*alpha_l[13]*fUpwind_l[15]+0.3162277660168379*alpha_l[3]*fUpwind_l[15]+0.3535533905932737*alpha_l[12]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[12]*alpha_l[14]+0.3162277660168379*alpha_l[1]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[1]*alpha_l[13]+0.3162277660168379*alpha_l[10]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[10]*alpha_l[11]+0.3535533905932737*alpha_l[2]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[2]*alpha_l[10]+0.3162277660168379*alpha_l[5]*fUpwind_l[9]+0.3162277660168379*alpha_l[5]*fUpwind_l[7]+0.3162277660168379*fUpwind_l[5]*alpha_l[7]+0.3535533905932737*alpha_l[4]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[4]*alpha_l[6]+0.3535533905932737*alpha_l[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alpha_l[5]+0.3535533905932737*alpha_l[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alpha_l[3]; 
  Ghat_l[6] = 0.2828427124746191*alpha_l[18]*fUpwind_l[19]+0.3162277660168379*alpha_l[5]*fUpwind_l[19]+0.3162277660168379*alpha_l[4]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[4]*alpha_l[18]+0.3535533905932737*alpha_l[7]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[7]*alpha_l[17]+0.2828427124746191*alpha_l[14]*fUpwind_l[16]+0.3162277660168379*alpha_l[3]*fUpwind_l[16]+0.3162277660168379*alpha_l[10]*fUpwind_l[15]+0.3162277660168379*alpha_l[2]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[2]*alpha_l[14]+0.3535533905932737*alpha_l[11]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[11]*alpha_l[13]+0.3162277660168379*alpha_l[10]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[10]*alpha_l[12]+0.3535533905932737*alpha_l[1]*fUpwind_l[10]+0.3535533905932737*fUpwind_l[1]*alpha_l[10]+0.3162277660168379*alpha_l[6]*fUpwind_l[9]+0.3162277660168379*alpha_l[6]*fUpwind_l[8]+0.3162277660168379*fUpwind_l[6]*alpha_l[8]+0.3535533905932737*alpha_l[0]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[0]*alpha_l[6]+0.3535533905932737*alpha_l[4]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[4]*alpha_l[5]+0.3535533905932737*alpha_l[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alpha_l[3]; 
  Ghat_l[7] = 0.3162277660168379*alpha_l[18]*fUpwind_l[18]+0.2258769757263128*alpha_l[17]*fUpwind_l[17]+0.3535533905932737*alpha_l[6]*fUpwind_l[17]+0.3535533905932737*fUpwind_l[6]*alpha_l[17]+0.2258769757263128*alpha_l[13]*fUpwind_l[13]+0.3535533905932737*alpha_l[3]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[3]*alpha_l[13]+0.3162277660168379*alpha_l[12]*fUpwind_l[12]+0.2258769757263128*alpha_l[11]*fUpwind_l[11]+0.3535533905932737*alpha_l[2]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[2]*alpha_l[11]+0.3162277660168379*alpha_l[10]*fUpwind_l[10]+0.2258769757263128*alpha_l[7]*fUpwind_l[7]+0.3535533905932737*alpha_l[0]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[0]*alpha_l[7]+0.3162277660168379*alpha_l[5]*fUpwind_l[5]+0.3162277660168379*alpha_l[4]*fUpwind_l[4]+0.3162277660168379*alpha_l[1]*fUpwind_l[1]; 
  Ghat_l[8] = 0.2258769757263128*alpha_l[18]*fUpwind_l[18]+0.3535533905932737*alpha_l[5]*fUpwind_l[18]+0.3535533905932737*fUpwind_l[5]*alpha_l[18]+0.3162277660168379*alpha_l[17]*fUpwind_l[17]+0.2258769757263128*alpha_l[14]*fUpwind_l[14]+0.3535533905932737*alpha_l[3]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[3]*alpha_l[14]+0.2258769757263128*alpha_l[12]*fUpwind_l[12]+0.3535533905932737*alpha_l[1]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[1]*alpha_l[12]+0.3162277660168379*alpha_l[11]*fUpwind_l[11]+0.3162277660168379*alpha_l[10]*fUpwind_l[10]+0.2258769757263128*alpha_l[8]*fUpwind_l[8]+0.3535533905932737*alpha_l[0]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[0]*alpha_l[8]+0.3162277660168379*alpha_l[6]*fUpwind_l[6]+0.3162277660168379*alpha_l[4]*fUpwind_l[4]+0.3162277660168379*alpha_l[2]*fUpwind_l[2]; 
  Ghat_l[9] = 0.3535533905932737*alpha_l[4]*fUpwind_l[19]+0.3162277660168379*alpha_l[18]*fUpwind_l[18]+0.3162277660168379*alpha_l[17]*fUpwind_l[17]+0.3535533905932737*alpha_l[2]*fUpwind_l[16]+0.3535533905932737*alpha_l[1]*fUpwind_l[15]+0.3162277660168379*alpha_l[14]*fUpwind_l[14]+0.3162277660168379*alpha_l[13]*fUpwind_l[13]+0.3162277660168379*alpha_l[10]*fUpwind_l[10]+0.3535533905932737*alpha_l[0]*fUpwind_l[9]+0.3162277660168379*alpha_l[6]*fUpwind_l[6]+0.3162277660168379*alpha_l[5]*fUpwind_l[5]+0.3162277660168379*alpha_l[3]*fUpwind_l[3]; 
  Ghat_l[10] = 0.282842712474619*alpha_l[14]*fUpwind_l[19]+0.282842712474619*alpha_l[13]*fUpwind_l[19]+0.3162277660168379*alpha_l[3]*fUpwind_l[19]+0.282842712474619*alpha_l[11]*fUpwind_l[18]+0.3162277660168379*alpha_l[2]*fUpwind_l[18]+0.282842712474619*fUpwind_l[16]*alpha_l[18]+0.282842712474619*fUpwind_l[11]*alpha_l[18]+0.3162277660168379*fUpwind_l[2]*alpha_l[18]+0.282842712474619*alpha_l[12]*fUpwind_l[17]+0.3162277660168379*alpha_l[1]*fUpwind_l[17]+0.282842712474619*fUpwind_l[15]*alpha_l[17]+0.282842712474619*fUpwind_l[12]*alpha_l[17]+0.3162277660168379*fUpwind_l[1]*alpha_l[17]+0.3162277660168379*alpha_l[5]*fUpwind_l[16]+0.3162277660168379*alpha_l[6]*fUpwind_l[15]+0.3162277660168379*alpha_l[4]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[4]*alpha_l[14]+0.3162277660168379*alpha_l[4]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[4]*alpha_l[13]+0.3162277660168379*alpha_l[6]*fUpwind_l[12]+0.3162277660168379*fUpwind_l[6]*alpha_l[12]+0.3162277660168379*alpha_l[5]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[5]*alpha_l[11]+0.3162277660168379*alpha_l[8]*fUpwind_l[10]+0.3162277660168379*alpha_l[7]*fUpwind_l[10]+0.3535533905932737*alpha_l[0]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[9]*alpha_l[10]+0.3162277660168379*fUpwind_l[8]*alpha_l[10]+0.3162277660168379*fUpwind_l[7]*alpha_l[10]+0.3535533905932737*fUpwind_l[0]*alpha_l[10]+0.3535533905932737*alpha_l[1]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[1]*alpha_l[6]+0.3535533905932737*alpha_l[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alpha_l[5]+0.3535533905932737*alpha_l[3]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[3]*alpha_l[4]; 
  Ghat_l[11] = 0.282842712474619*alpha_l[10]*fUpwind_l[18]+0.282842712474619*fUpwind_l[10]*alpha_l[18]+0.3162277660168379*alpha_l[14]*fUpwind_l[17]+0.2258769757263128*alpha_l[13]*fUpwind_l[17]+0.3535533905932737*alpha_l[3]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[14]*alpha_l[17]+0.2258769757263128*fUpwind_l[13]*alpha_l[17]+0.3535533905932737*fUpwind_l[3]*alpha_l[17]+0.3535533905932737*alpha_l[6]*fUpwind_l[13]+0.3535533905932737*fUpwind_l[6]*alpha_l[13]+0.2828427124746191*alpha_l[4]*fUpwind_l[12]+0.2828427124746191*fUpwind_l[4]*alpha_l[12]+0.3162277660168379*alpha_l[8]*fUpwind_l[11]+0.2258769757263128*alpha_l[7]*fUpwind_l[11]+0.3535533905932737*alpha_l[0]*fUpwind_l[11]+0.3162277660168379*fUpwind_l[8]*alpha_l[11]+0.2258769757263128*fUpwind_l[7]*alpha_l[11]+0.3535533905932737*fUpwind_l[0]*alpha_l[11]+0.3162277660168379*alpha_l[5]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[5]*alpha_l[10]+0.3535533905932737*alpha_l[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[2]*alpha_l[7]+0.3162277660168379*alpha_l[1]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[1]*alpha_l[4]; 
  Ghat_l[12] = 0.2258769757263128*alpha_l[14]*fUpwind_l[18]+0.3162277660168379*alpha_l[13]*fUpwind_l[18]+0.3535533905932737*alpha_l[3]*fUpwind_l[18]+0.2258769757263128*fUpwind_l[14]*alpha_l[18]+0.3162277660168379*fUpwind_l[13]*alpha_l[18]+0.3535533905932737*fUpwind_l[3]*alpha_l[18]+0.282842712474619*alpha_l[10]*fUpwind_l[17]+0.282842712474619*fUpwind_l[10]*alpha_l[17]+0.3535533905932737*alpha_l[5]*fUpwind_l[14]+0.3535533905932737*fUpwind_l[5]*alpha_l[14]+0.2258769757263128*alpha_l[8]*fUpwind_l[12]+0.3162277660168379*alpha_l[7]*fUpwind_l[12]+0.3535533905932737*alpha_l[0]*fUpwind_l[12]+0.2258769757263128*fUpwind_l[8]*alpha_l[12]+0.3162277660168379*fUpwind_l[7]*alpha_l[12]+0.3535533905932737*fUpwind_l[0]*alpha_l[12]+0.2828427124746191*alpha_l[4]*fUpwind_l[11]+0.2828427124746191*fUpwind_l[4]*alpha_l[11]+0.3162277660168379*alpha_l[6]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[6]*alpha_l[10]+0.3535533905932737*alpha_l[1]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[1]*alpha_l[8]+0.3162277660168379*alpha_l[2]*fUpwind_l[4]+0.3162277660168379*fUpwind_l[2]*alpha_l[4]; 
  Ghat_l[13] = 0.282842712474619*alpha_l[10]*fUpwind_l[19]+0.3162277660168379*alpha_l[12]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[12]*alpha_l[18]+0.2258769757263128*alpha_l[11]*fUpwind_l[17]+0.3535533905932737*alpha_l[2]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[16]*alpha_l[17]+0.2258769757263128*fUpwind_l[11]*alpha_l[17]+0.3535533905932737*fUpwind_l[2]*alpha_l[17]+0.2828427124746191*alpha_l[5]*fUpwind_l[15]+0.2258769757263128*alpha_l[7]*fUpwind_l[13]+0.3535533905932737*alpha_l[0]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[9]*alpha_l[13]+0.2258769757263128*fUpwind_l[7]*alpha_l[13]+0.3535533905932737*fUpwind_l[0]*alpha_l[13]+0.3535533905932737*alpha_l[6]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[6]*alpha_l[11]+0.3162277660168379*alpha_l[4]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[4]*alpha_l[10]+0.3535533905932737*alpha_l[3]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[3]*alpha_l[7]+0.3162277660168379*alpha_l[1]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[1]*alpha_l[5]; 
  Ghat_l[14] = 0.282842712474619*alpha_l[10]*fUpwind_l[19]+0.2258769757263128*alpha_l[12]*fUpwind_l[18]+0.3535533905932737*alpha_l[1]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[15]*alpha_l[18]+0.2258769757263128*fUpwind_l[12]*alpha_l[18]+0.3535533905932737*fUpwind_l[1]*alpha_l[18]+0.3162277660168379*alpha_l[11]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[11]*alpha_l[17]+0.2828427124746191*alpha_l[6]*fUpwind_l[16]+0.2258769757263128*alpha_l[8]*fUpwind_l[14]+0.3535533905932737*alpha_l[0]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[9]*alpha_l[14]+0.2258769757263128*fUpwind_l[8]*alpha_l[14]+0.3535533905932737*fUpwind_l[0]*alpha_l[14]+0.3535533905932737*alpha_l[5]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[5]*alpha_l[12]+0.3162277660168379*alpha_l[4]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[4]*alpha_l[10]+0.3535533905932737*alpha_l[3]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[3]*alpha_l[8]+0.3162277660168379*alpha_l[2]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[2]*alpha_l[6]; 
  Ghat_l[15] = 0.3162277660168379*alpha_l[11]*fUpwind_l[19]+0.3535533905932737*alpha_l[2]*fUpwind_l[19]+0.3162277660168379*alpha_l[14]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[14]*alpha_l[18]+0.282842712474619*alpha_l[10]*fUpwind_l[17]+0.282842712474619*fUpwind_l[10]*alpha_l[17]+0.3535533905932737*alpha_l[4]*fUpwind_l[16]+0.3162277660168379*alpha_l[7]*fUpwind_l[15]+0.3535533905932737*alpha_l[0]*fUpwind_l[15]+0.2828427124746191*alpha_l[5]*fUpwind_l[13]+0.2828427124746191*fUpwind_l[5]*alpha_l[13]+0.3162277660168379*alpha_l[6]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[6]*alpha_l[10]+0.3535533905932737*alpha_l[1]*fUpwind_l[9]+0.3162277660168379*alpha_l[3]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[3]*alpha_l[5]; 
  Ghat_l[16] = 0.3162277660168379*alpha_l[12]*fUpwind_l[19]+0.3535533905932737*alpha_l[1]*fUpwind_l[19]+0.282842712474619*alpha_l[10]*fUpwind_l[18]+0.282842712474619*fUpwind_l[10]*alpha_l[18]+0.3162277660168379*alpha_l[13]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[13]*alpha_l[17]+0.3162277660168379*alpha_l[8]*fUpwind_l[16]+0.3535533905932737*alpha_l[0]*fUpwind_l[16]+0.3535533905932737*alpha_l[4]*fUpwind_l[15]+0.2828427124746191*alpha_l[6]*fUpwind_l[14]+0.2828427124746191*fUpwind_l[6]*alpha_l[14]+0.3162277660168379*alpha_l[5]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[5]*alpha_l[10]+0.3535533905932737*alpha_l[2]*fUpwind_l[9]+0.3162277660168379*alpha_l[3]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[3]*alpha_l[6]; 
  Ghat_l[17] = 0.2529822128134704*alpha_l[18]*fUpwind_l[19]+0.2828427124746191*alpha_l[5]*fUpwind_l[19]+0.2828427124746191*alpha_l[4]*fUpwind_l[18]+0.2828427124746191*fUpwind_l[4]*alpha_l[18]+0.3162277660168379*alpha_l[8]*fUpwind_l[17]+0.2258769757263128*alpha_l[7]*fUpwind_l[17]+0.3535533905932737*alpha_l[0]*fUpwind_l[17]+0.3162277660168379*fUpwind_l[9]*alpha_l[17]+0.3162277660168379*fUpwind_l[8]*alpha_l[17]+0.2258769757263128*fUpwind_l[7]*alpha_l[17]+0.3535533905932737*fUpwind_l[0]*alpha_l[17]+0.3162277660168379*alpha_l[13]*fUpwind_l[16]+0.282842712474619*alpha_l[10]*fUpwind_l[15]+0.3162277660168379*alpha_l[11]*fUpwind_l[14]+0.3162277660168379*fUpwind_l[11]*alpha_l[14]+0.2258769757263128*alpha_l[11]*fUpwind_l[13]+0.3535533905932737*alpha_l[2]*fUpwind_l[13]+0.2258769757263128*fUpwind_l[11]*alpha_l[13]+0.3535533905932737*fUpwind_l[2]*alpha_l[13]+0.282842712474619*alpha_l[10]*fUpwind_l[12]+0.282842712474619*fUpwind_l[10]*alpha_l[12]+0.3535533905932737*alpha_l[3]*fUpwind_l[11]+0.3535533905932737*fUpwind_l[3]*alpha_l[11]+0.3162277660168379*alpha_l[1]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[1]*alpha_l[10]+0.3535533905932737*alpha_l[6]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[6]*alpha_l[7]+0.3162277660168379*alpha_l[4]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[4]*alpha_l[5]; 
  Ghat_l[18] = 0.2529822128134704*alpha_l[17]*fUpwind_l[19]+0.2828427124746191*alpha_l[6]*fUpwind_l[19]+0.2258769757263128*alpha_l[8]*fUpwind_l[18]+0.3162277660168379*alpha_l[7]*fUpwind_l[18]+0.3535533905932737*alpha_l[0]*fUpwind_l[18]+0.3162277660168379*fUpwind_l[9]*alpha_l[18]+0.2258769757263128*fUpwind_l[8]*alpha_l[18]+0.3162277660168379*fUpwind_l[7]*alpha_l[18]+0.3535533905932737*fUpwind_l[0]*alpha_l[18]+0.2828427124746191*alpha_l[4]*fUpwind_l[17]+0.2828427124746191*fUpwind_l[4]*alpha_l[17]+0.282842712474619*alpha_l[10]*fUpwind_l[16]+0.3162277660168379*alpha_l[14]*fUpwind_l[15]+0.2258769757263128*alpha_l[12]*fUpwind_l[14]+0.3535533905932737*alpha_l[1]*fUpwind_l[14]+0.2258769757263128*fUpwind_l[12]*alpha_l[14]+0.3535533905932737*fUpwind_l[1]*alpha_l[14]+0.3162277660168379*alpha_l[12]*fUpwind_l[13]+0.3162277660168379*fUpwind_l[12]*alpha_l[13]+0.3535533905932737*alpha_l[3]*fUpwind_l[12]+0.3535533905932737*fUpwind_l[3]*alpha_l[12]+0.282842712474619*alpha_l[10]*fUpwind_l[11]+0.282842712474619*fUpwind_l[10]*alpha_l[11]+0.3162277660168379*alpha_l[2]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[2]*alpha_l[10]+0.3535533905932737*alpha_l[5]*fUpwind_l[8]+0.3535533905932737*fUpwind_l[5]*alpha_l[8]+0.3162277660168379*alpha_l[4]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[4]*alpha_l[6]; 
  Ghat_l[19] = 0.3162277660168379*alpha_l[8]*fUpwind_l[19]+0.3162277660168379*alpha_l[7]*fUpwind_l[19]+0.3535533905932737*alpha_l[0]*fUpwind_l[19]+0.2529822128134704*alpha_l[17]*fUpwind_l[18]+0.2828427124746191*alpha_l[6]*fUpwind_l[18]+0.2529822128134704*fUpwind_l[17]*alpha_l[18]+0.2828427124746191*fUpwind_l[6]*alpha_l[18]+0.2828427124746191*alpha_l[5]*fUpwind_l[17]+0.2828427124746191*fUpwind_l[5]*alpha_l[17]+0.3162277660168379*alpha_l[12]*fUpwind_l[16]+0.3535533905932737*alpha_l[1]*fUpwind_l[16]+0.3162277660168379*alpha_l[11]*fUpwind_l[15]+0.3535533905932737*alpha_l[2]*fUpwind_l[15]+0.282842712474619*alpha_l[10]*fUpwind_l[14]+0.282842712474619*fUpwind_l[10]*alpha_l[14]+0.282842712474619*alpha_l[10]*fUpwind_l[13]+0.282842712474619*fUpwind_l[10]*alpha_l[13]+0.3162277660168379*alpha_l[3]*fUpwind_l[10]+0.3162277660168379*fUpwind_l[3]*alpha_l[10]+0.3535533905932737*alpha_l[4]*fUpwind_l[9]+0.3162277660168379*alpha_l[5]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[5]*alpha_l[6]; 

  Ghat_r[0] = 0.3535533905932737*alpha_r[18]*fUpwind_r[18]+0.3535533905932737*alpha_r[17]*fUpwind_r[17]+0.3535533905932737*alpha_r[14]*fUpwind_r[14]+0.3535533905932737*alpha_r[13]*fUpwind_r[13]+0.3535533905932737*alpha_r[12]*fUpwind_r[12]+0.3535533905932737*alpha_r[11]*fUpwind_r[11]+0.3535533905932737*alpha_r[10]*fUpwind_r[10]+0.3535533905932737*alpha_r[8]*fUpwind_r[8]+0.3535533905932737*alpha_r[7]*fUpwind_r[7]+0.3535533905932737*alpha_r[6]*fUpwind_r[6]+0.3535533905932737*alpha_r[5]*fUpwind_r[5]+0.3535533905932737*alpha_r[4]*fUpwind_r[4]+0.3535533905932737*alpha_r[3]*fUpwind_r[3]+0.3535533905932737*alpha_r[2]*fUpwind_r[2]+0.3535533905932737*alpha_r[1]*fUpwind_r[1]+0.3535533905932737*alpha_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3535533905932737*alpha_r[14]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[14]*alpha_r[18]+0.3162277660168379*alpha_r[10]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[10]*alpha_r[17]+0.3162277660168379*alpha_r[5]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[5]*alpha_r[13]+0.3535533905932737*alpha_r[8]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[8]*alpha_r[12]+0.3162277660168379*alpha_r[4]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[4]*alpha_r[11]+0.3535533905932737*alpha_r[6]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[6]*alpha_r[10]+0.3162277660168379*alpha_r[1]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[1]*alpha_r[7]+0.3535533905932737*alpha_r[3]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alpha_r[5]+0.3535533905932737*alpha_r[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alpha_r[4]+0.3535533905932737*alpha_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alpha_r[1]; 
  Ghat_r[2] = 0.3162277660168379*alpha_r[10]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[10]*alpha_r[18]+0.3535533905932737*alpha_r[13]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[13]*alpha_r[17]+0.3162277660168379*alpha_r[6]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[6]*alpha_r[14]+0.3162277660168379*alpha_r[4]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[4]*alpha_r[12]+0.3535533905932737*alpha_r[7]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[7]*alpha_r[11]+0.3535533905932737*alpha_r[5]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[5]*alpha_r[10]+0.3162277660168379*alpha_r[2]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[2]*alpha_r[8]+0.3535533905932737*alpha_r[3]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[3]*alpha_r[6]+0.3535533905932737*alpha_r[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alpha_r[4]+0.3535533905932737*alpha_r[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha_r[2]; 
  Ghat_r[3] = 0.3162277660168379*alpha_r[10]*fUpwind_r[19]+0.3535533905932737*alpha_r[12]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[12]*alpha_r[18]+0.3535533905932737*alpha_r[11]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[11]*alpha_r[17]+0.3162277660168379*alpha_r[6]*fUpwind_r[16]+0.3162277660168379*alpha_r[5]*fUpwind_r[15]+0.3535533905932737*alpha_r[8]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[8]*alpha_r[14]+0.3535533905932737*alpha_r[7]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[7]*alpha_r[13]+0.3535533905932737*alpha_r[4]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[4]*alpha_r[10]+0.3162277660168379*alpha_r[3]*fUpwind_r[9]+0.3535533905932737*alpha_r[2]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[2]*alpha_r[6]+0.3535533905932737*alpha_r[1]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[1]*alpha_r[5]+0.3535533905932737*alpha_r[0]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[0]*alpha_r[3]; 
  Ghat_r[4] = 0.2828427124746191*alpha_r[17]*fUpwind_r[18]+0.3162277660168379*alpha_r[6]*fUpwind_r[18]+0.2828427124746191*fUpwind_r[17]*alpha_r[18]+0.3162277660168379*fUpwind_r[6]*alpha_r[18]+0.3162277660168379*alpha_r[5]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[5]*alpha_r[17]+0.3162277660168379*alpha_r[10]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[10]*alpha_r[14]+0.3162277660168379*alpha_r[10]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[10]*alpha_r[13]+0.2828427124746191*alpha_r[11]*fUpwind_r[12]+0.3162277660168379*alpha_r[2]*fUpwind_r[12]+0.2828427124746191*fUpwind_r[11]*alpha_r[12]+0.3162277660168379*fUpwind_r[2]*alpha_r[12]+0.3162277660168379*alpha_r[1]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[1]*alpha_r[11]+0.3535533905932737*alpha_r[3]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[3]*alpha_r[10]+0.3162277660168379*alpha_r[4]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[4]*alpha_r[8]+0.3162277660168379*alpha_r[4]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[4]*alpha_r[7]+0.3535533905932737*alpha_r[5]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[5]*alpha_r[6]+0.3535533905932737*alpha_r[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alpha_r[4]+0.3535533905932737*alpha_r[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alpha_r[2]; 
  Ghat_r[5] = 0.2828427124746191*alpha_r[17]*fUpwind_r[19]+0.3162277660168379*alpha_r[6]*fUpwind_r[19]+0.3535533905932737*alpha_r[8]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[8]*alpha_r[18]+0.3162277660168379*alpha_r[4]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[4]*alpha_r[17]+0.3162277660168379*alpha_r[10]*fUpwind_r[16]+0.2828427124746191*alpha_r[13]*fUpwind_r[15]+0.3162277660168379*alpha_r[3]*fUpwind_r[15]+0.3535533905932737*alpha_r[12]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[12]*alpha_r[14]+0.3162277660168379*alpha_r[1]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[1]*alpha_r[13]+0.3162277660168379*alpha_r[10]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[10]*alpha_r[11]+0.3535533905932737*alpha_r[2]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[2]*alpha_r[10]+0.3162277660168379*alpha_r[5]*fUpwind_r[9]+0.3162277660168379*alpha_r[5]*fUpwind_r[7]+0.3162277660168379*fUpwind_r[5]*alpha_r[7]+0.3535533905932737*alpha_r[4]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[4]*alpha_r[6]+0.3535533905932737*alpha_r[0]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[0]*alpha_r[5]+0.3535533905932737*alpha_r[1]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[1]*alpha_r[3]; 
  Ghat_r[6] = 0.2828427124746191*alpha_r[18]*fUpwind_r[19]+0.3162277660168379*alpha_r[5]*fUpwind_r[19]+0.3162277660168379*alpha_r[4]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[4]*alpha_r[18]+0.3535533905932737*alpha_r[7]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[7]*alpha_r[17]+0.2828427124746191*alpha_r[14]*fUpwind_r[16]+0.3162277660168379*alpha_r[3]*fUpwind_r[16]+0.3162277660168379*alpha_r[10]*fUpwind_r[15]+0.3162277660168379*alpha_r[2]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[2]*alpha_r[14]+0.3535533905932737*alpha_r[11]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[11]*alpha_r[13]+0.3162277660168379*alpha_r[10]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[10]*alpha_r[12]+0.3535533905932737*alpha_r[1]*fUpwind_r[10]+0.3535533905932737*fUpwind_r[1]*alpha_r[10]+0.3162277660168379*alpha_r[6]*fUpwind_r[9]+0.3162277660168379*alpha_r[6]*fUpwind_r[8]+0.3162277660168379*fUpwind_r[6]*alpha_r[8]+0.3535533905932737*alpha_r[0]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[0]*alpha_r[6]+0.3535533905932737*alpha_r[4]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[4]*alpha_r[5]+0.3535533905932737*alpha_r[2]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[2]*alpha_r[3]; 
  Ghat_r[7] = 0.3162277660168379*alpha_r[18]*fUpwind_r[18]+0.2258769757263128*alpha_r[17]*fUpwind_r[17]+0.3535533905932737*alpha_r[6]*fUpwind_r[17]+0.3535533905932737*fUpwind_r[6]*alpha_r[17]+0.2258769757263128*alpha_r[13]*fUpwind_r[13]+0.3535533905932737*alpha_r[3]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[3]*alpha_r[13]+0.3162277660168379*alpha_r[12]*fUpwind_r[12]+0.2258769757263128*alpha_r[11]*fUpwind_r[11]+0.3535533905932737*alpha_r[2]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[2]*alpha_r[11]+0.3162277660168379*alpha_r[10]*fUpwind_r[10]+0.2258769757263128*alpha_r[7]*fUpwind_r[7]+0.3535533905932737*alpha_r[0]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[0]*alpha_r[7]+0.3162277660168379*alpha_r[5]*fUpwind_r[5]+0.3162277660168379*alpha_r[4]*fUpwind_r[4]+0.3162277660168379*alpha_r[1]*fUpwind_r[1]; 
  Ghat_r[8] = 0.2258769757263128*alpha_r[18]*fUpwind_r[18]+0.3535533905932737*alpha_r[5]*fUpwind_r[18]+0.3535533905932737*fUpwind_r[5]*alpha_r[18]+0.3162277660168379*alpha_r[17]*fUpwind_r[17]+0.2258769757263128*alpha_r[14]*fUpwind_r[14]+0.3535533905932737*alpha_r[3]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[3]*alpha_r[14]+0.2258769757263128*alpha_r[12]*fUpwind_r[12]+0.3535533905932737*alpha_r[1]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[1]*alpha_r[12]+0.3162277660168379*alpha_r[11]*fUpwind_r[11]+0.3162277660168379*alpha_r[10]*fUpwind_r[10]+0.2258769757263128*alpha_r[8]*fUpwind_r[8]+0.3535533905932737*alpha_r[0]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[0]*alpha_r[8]+0.3162277660168379*alpha_r[6]*fUpwind_r[6]+0.3162277660168379*alpha_r[4]*fUpwind_r[4]+0.3162277660168379*alpha_r[2]*fUpwind_r[2]; 
  Ghat_r[9] = 0.3535533905932737*alpha_r[4]*fUpwind_r[19]+0.3162277660168379*alpha_r[18]*fUpwind_r[18]+0.3162277660168379*alpha_r[17]*fUpwind_r[17]+0.3535533905932737*alpha_r[2]*fUpwind_r[16]+0.3535533905932737*alpha_r[1]*fUpwind_r[15]+0.3162277660168379*alpha_r[14]*fUpwind_r[14]+0.3162277660168379*alpha_r[13]*fUpwind_r[13]+0.3162277660168379*alpha_r[10]*fUpwind_r[10]+0.3535533905932737*alpha_r[0]*fUpwind_r[9]+0.3162277660168379*alpha_r[6]*fUpwind_r[6]+0.3162277660168379*alpha_r[5]*fUpwind_r[5]+0.3162277660168379*alpha_r[3]*fUpwind_r[3]; 
  Ghat_r[10] = 0.282842712474619*alpha_r[14]*fUpwind_r[19]+0.282842712474619*alpha_r[13]*fUpwind_r[19]+0.3162277660168379*alpha_r[3]*fUpwind_r[19]+0.282842712474619*alpha_r[11]*fUpwind_r[18]+0.3162277660168379*alpha_r[2]*fUpwind_r[18]+0.282842712474619*fUpwind_r[16]*alpha_r[18]+0.282842712474619*fUpwind_r[11]*alpha_r[18]+0.3162277660168379*fUpwind_r[2]*alpha_r[18]+0.282842712474619*alpha_r[12]*fUpwind_r[17]+0.3162277660168379*alpha_r[1]*fUpwind_r[17]+0.282842712474619*fUpwind_r[15]*alpha_r[17]+0.282842712474619*fUpwind_r[12]*alpha_r[17]+0.3162277660168379*fUpwind_r[1]*alpha_r[17]+0.3162277660168379*alpha_r[5]*fUpwind_r[16]+0.3162277660168379*alpha_r[6]*fUpwind_r[15]+0.3162277660168379*alpha_r[4]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[4]*alpha_r[14]+0.3162277660168379*alpha_r[4]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[4]*alpha_r[13]+0.3162277660168379*alpha_r[6]*fUpwind_r[12]+0.3162277660168379*fUpwind_r[6]*alpha_r[12]+0.3162277660168379*alpha_r[5]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[5]*alpha_r[11]+0.3162277660168379*alpha_r[8]*fUpwind_r[10]+0.3162277660168379*alpha_r[7]*fUpwind_r[10]+0.3535533905932737*alpha_r[0]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[9]*alpha_r[10]+0.3162277660168379*fUpwind_r[8]*alpha_r[10]+0.3162277660168379*fUpwind_r[7]*alpha_r[10]+0.3535533905932737*fUpwind_r[0]*alpha_r[10]+0.3535533905932737*alpha_r[1]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[1]*alpha_r[6]+0.3535533905932737*alpha_r[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[2]*alpha_r[5]+0.3535533905932737*alpha_r[3]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[3]*alpha_r[4]; 
  Ghat_r[11] = 0.282842712474619*alpha_r[10]*fUpwind_r[18]+0.282842712474619*fUpwind_r[10]*alpha_r[18]+0.3162277660168379*alpha_r[14]*fUpwind_r[17]+0.2258769757263128*alpha_r[13]*fUpwind_r[17]+0.3535533905932737*alpha_r[3]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[14]*alpha_r[17]+0.2258769757263128*fUpwind_r[13]*alpha_r[17]+0.3535533905932737*fUpwind_r[3]*alpha_r[17]+0.3535533905932737*alpha_r[6]*fUpwind_r[13]+0.3535533905932737*fUpwind_r[6]*alpha_r[13]+0.2828427124746191*alpha_r[4]*fUpwind_r[12]+0.2828427124746191*fUpwind_r[4]*alpha_r[12]+0.3162277660168379*alpha_r[8]*fUpwind_r[11]+0.2258769757263128*alpha_r[7]*fUpwind_r[11]+0.3535533905932737*alpha_r[0]*fUpwind_r[11]+0.3162277660168379*fUpwind_r[8]*alpha_r[11]+0.2258769757263128*fUpwind_r[7]*alpha_r[11]+0.3535533905932737*fUpwind_r[0]*alpha_r[11]+0.3162277660168379*alpha_r[5]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[5]*alpha_r[10]+0.3535533905932737*alpha_r[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[2]*alpha_r[7]+0.3162277660168379*alpha_r[1]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[1]*alpha_r[4]; 
  Ghat_r[12] = 0.2258769757263128*alpha_r[14]*fUpwind_r[18]+0.3162277660168379*alpha_r[13]*fUpwind_r[18]+0.3535533905932737*alpha_r[3]*fUpwind_r[18]+0.2258769757263128*fUpwind_r[14]*alpha_r[18]+0.3162277660168379*fUpwind_r[13]*alpha_r[18]+0.3535533905932737*fUpwind_r[3]*alpha_r[18]+0.282842712474619*alpha_r[10]*fUpwind_r[17]+0.282842712474619*fUpwind_r[10]*alpha_r[17]+0.3535533905932737*alpha_r[5]*fUpwind_r[14]+0.3535533905932737*fUpwind_r[5]*alpha_r[14]+0.2258769757263128*alpha_r[8]*fUpwind_r[12]+0.3162277660168379*alpha_r[7]*fUpwind_r[12]+0.3535533905932737*alpha_r[0]*fUpwind_r[12]+0.2258769757263128*fUpwind_r[8]*alpha_r[12]+0.3162277660168379*fUpwind_r[7]*alpha_r[12]+0.3535533905932737*fUpwind_r[0]*alpha_r[12]+0.2828427124746191*alpha_r[4]*fUpwind_r[11]+0.2828427124746191*fUpwind_r[4]*alpha_r[11]+0.3162277660168379*alpha_r[6]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[6]*alpha_r[10]+0.3535533905932737*alpha_r[1]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[1]*alpha_r[8]+0.3162277660168379*alpha_r[2]*fUpwind_r[4]+0.3162277660168379*fUpwind_r[2]*alpha_r[4]; 
  Ghat_r[13] = 0.282842712474619*alpha_r[10]*fUpwind_r[19]+0.3162277660168379*alpha_r[12]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[12]*alpha_r[18]+0.2258769757263128*alpha_r[11]*fUpwind_r[17]+0.3535533905932737*alpha_r[2]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[16]*alpha_r[17]+0.2258769757263128*fUpwind_r[11]*alpha_r[17]+0.3535533905932737*fUpwind_r[2]*alpha_r[17]+0.2828427124746191*alpha_r[5]*fUpwind_r[15]+0.2258769757263128*alpha_r[7]*fUpwind_r[13]+0.3535533905932737*alpha_r[0]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[9]*alpha_r[13]+0.2258769757263128*fUpwind_r[7]*alpha_r[13]+0.3535533905932737*fUpwind_r[0]*alpha_r[13]+0.3535533905932737*alpha_r[6]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[6]*alpha_r[11]+0.3162277660168379*alpha_r[4]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[4]*alpha_r[10]+0.3535533905932737*alpha_r[3]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[3]*alpha_r[7]+0.3162277660168379*alpha_r[1]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[1]*alpha_r[5]; 
  Ghat_r[14] = 0.282842712474619*alpha_r[10]*fUpwind_r[19]+0.2258769757263128*alpha_r[12]*fUpwind_r[18]+0.3535533905932737*alpha_r[1]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[15]*alpha_r[18]+0.2258769757263128*fUpwind_r[12]*alpha_r[18]+0.3535533905932737*fUpwind_r[1]*alpha_r[18]+0.3162277660168379*alpha_r[11]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[11]*alpha_r[17]+0.2828427124746191*alpha_r[6]*fUpwind_r[16]+0.2258769757263128*alpha_r[8]*fUpwind_r[14]+0.3535533905932737*alpha_r[0]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[9]*alpha_r[14]+0.2258769757263128*fUpwind_r[8]*alpha_r[14]+0.3535533905932737*fUpwind_r[0]*alpha_r[14]+0.3535533905932737*alpha_r[5]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[5]*alpha_r[12]+0.3162277660168379*alpha_r[4]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[4]*alpha_r[10]+0.3535533905932737*alpha_r[3]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[3]*alpha_r[8]+0.3162277660168379*alpha_r[2]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[2]*alpha_r[6]; 
  Ghat_r[15] = 0.3162277660168379*alpha_r[11]*fUpwind_r[19]+0.3535533905932737*alpha_r[2]*fUpwind_r[19]+0.3162277660168379*alpha_r[14]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[14]*alpha_r[18]+0.282842712474619*alpha_r[10]*fUpwind_r[17]+0.282842712474619*fUpwind_r[10]*alpha_r[17]+0.3535533905932737*alpha_r[4]*fUpwind_r[16]+0.3162277660168379*alpha_r[7]*fUpwind_r[15]+0.3535533905932737*alpha_r[0]*fUpwind_r[15]+0.2828427124746191*alpha_r[5]*fUpwind_r[13]+0.2828427124746191*fUpwind_r[5]*alpha_r[13]+0.3162277660168379*alpha_r[6]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[6]*alpha_r[10]+0.3535533905932737*alpha_r[1]*fUpwind_r[9]+0.3162277660168379*alpha_r[3]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[3]*alpha_r[5]; 
  Ghat_r[16] = 0.3162277660168379*alpha_r[12]*fUpwind_r[19]+0.3535533905932737*alpha_r[1]*fUpwind_r[19]+0.282842712474619*alpha_r[10]*fUpwind_r[18]+0.282842712474619*fUpwind_r[10]*alpha_r[18]+0.3162277660168379*alpha_r[13]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[13]*alpha_r[17]+0.3162277660168379*alpha_r[8]*fUpwind_r[16]+0.3535533905932737*alpha_r[0]*fUpwind_r[16]+0.3535533905932737*alpha_r[4]*fUpwind_r[15]+0.2828427124746191*alpha_r[6]*fUpwind_r[14]+0.2828427124746191*fUpwind_r[6]*alpha_r[14]+0.3162277660168379*alpha_r[5]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[5]*alpha_r[10]+0.3535533905932737*alpha_r[2]*fUpwind_r[9]+0.3162277660168379*alpha_r[3]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[3]*alpha_r[6]; 
  Ghat_r[17] = 0.2529822128134704*alpha_r[18]*fUpwind_r[19]+0.2828427124746191*alpha_r[5]*fUpwind_r[19]+0.2828427124746191*alpha_r[4]*fUpwind_r[18]+0.2828427124746191*fUpwind_r[4]*alpha_r[18]+0.3162277660168379*alpha_r[8]*fUpwind_r[17]+0.2258769757263128*alpha_r[7]*fUpwind_r[17]+0.3535533905932737*alpha_r[0]*fUpwind_r[17]+0.3162277660168379*fUpwind_r[9]*alpha_r[17]+0.3162277660168379*fUpwind_r[8]*alpha_r[17]+0.2258769757263128*fUpwind_r[7]*alpha_r[17]+0.3535533905932737*fUpwind_r[0]*alpha_r[17]+0.3162277660168379*alpha_r[13]*fUpwind_r[16]+0.282842712474619*alpha_r[10]*fUpwind_r[15]+0.3162277660168379*alpha_r[11]*fUpwind_r[14]+0.3162277660168379*fUpwind_r[11]*alpha_r[14]+0.2258769757263128*alpha_r[11]*fUpwind_r[13]+0.3535533905932737*alpha_r[2]*fUpwind_r[13]+0.2258769757263128*fUpwind_r[11]*alpha_r[13]+0.3535533905932737*fUpwind_r[2]*alpha_r[13]+0.282842712474619*alpha_r[10]*fUpwind_r[12]+0.282842712474619*fUpwind_r[10]*alpha_r[12]+0.3535533905932737*alpha_r[3]*fUpwind_r[11]+0.3535533905932737*fUpwind_r[3]*alpha_r[11]+0.3162277660168379*alpha_r[1]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[1]*alpha_r[10]+0.3535533905932737*alpha_r[6]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[6]*alpha_r[7]+0.3162277660168379*alpha_r[4]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[4]*alpha_r[5]; 
  Ghat_r[18] = 0.2529822128134704*alpha_r[17]*fUpwind_r[19]+0.2828427124746191*alpha_r[6]*fUpwind_r[19]+0.2258769757263128*alpha_r[8]*fUpwind_r[18]+0.3162277660168379*alpha_r[7]*fUpwind_r[18]+0.3535533905932737*alpha_r[0]*fUpwind_r[18]+0.3162277660168379*fUpwind_r[9]*alpha_r[18]+0.2258769757263128*fUpwind_r[8]*alpha_r[18]+0.3162277660168379*fUpwind_r[7]*alpha_r[18]+0.3535533905932737*fUpwind_r[0]*alpha_r[18]+0.2828427124746191*alpha_r[4]*fUpwind_r[17]+0.2828427124746191*fUpwind_r[4]*alpha_r[17]+0.282842712474619*alpha_r[10]*fUpwind_r[16]+0.3162277660168379*alpha_r[14]*fUpwind_r[15]+0.2258769757263128*alpha_r[12]*fUpwind_r[14]+0.3535533905932737*alpha_r[1]*fUpwind_r[14]+0.2258769757263128*fUpwind_r[12]*alpha_r[14]+0.3535533905932737*fUpwind_r[1]*alpha_r[14]+0.3162277660168379*alpha_r[12]*fUpwind_r[13]+0.3162277660168379*fUpwind_r[12]*alpha_r[13]+0.3535533905932737*alpha_r[3]*fUpwind_r[12]+0.3535533905932737*fUpwind_r[3]*alpha_r[12]+0.282842712474619*alpha_r[10]*fUpwind_r[11]+0.282842712474619*fUpwind_r[10]*alpha_r[11]+0.3162277660168379*alpha_r[2]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[2]*alpha_r[10]+0.3535533905932737*alpha_r[5]*fUpwind_r[8]+0.3535533905932737*fUpwind_r[5]*alpha_r[8]+0.3162277660168379*alpha_r[4]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[4]*alpha_r[6]; 
  Ghat_r[19] = 0.3162277660168379*alpha_r[8]*fUpwind_r[19]+0.3162277660168379*alpha_r[7]*fUpwind_r[19]+0.3535533905932737*alpha_r[0]*fUpwind_r[19]+0.2529822128134704*alpha_r[17]*fUpwind_r[18]+0.2828427124746191*alpha_r[6]*fUpwind_r[18]+0.2529822128134704*fUpwind_r[17]*alpha_r[18]+0.2828427124746191*fUpwind_r[6]*alpha_r[18]+0.2828427124746191*alpha_r[5]*fUpwind_r[17]+0.2828427124746191*fUpwind_r[5]*alpha_r[17]+0.3162277660168379*alpha_r[12]*fUpwind_r[16]+0.3535533905932737*alpha_r[1]*fUpwind_r[16]+0.3162277660168379*alpha_r[11]*fUpwind_r[15]+0.3535533905932737*alpha_r[2]*fUpwind_r[15]+0.282842712474619*alpha_r[10]*fUpwind_r[14]+0.282842712474619*fUpwind_r[10]*alpha_r[14]+0.282842712474619*alpha_r[10]*fUpwind_r[13]+0.282842712474619*fUpwind_r[10]*alpha_r[13]+0.3162277660168379*alpha_r[3]*fUpwind_r[10]+0.3162277660168379*fUpwind_r[3]*alpha_r[10]+0.3535533905932737*alpha_r[4]*fUpwind_r[9]+0.3162277660168379*alpha_r[5]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[5]*alpha_r[6]; 

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

  return 0.;

} 
