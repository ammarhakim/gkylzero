#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double lbo_vlasov_drag_boundary_surfvy_2x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum[24]: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // fSkin/Edge: Distribution function in cells 
  // out: Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 

  const double *sumNuUy = &nuPrimMomsSum[8]; 

  double alphaDrSurf[20] = {0.0}; 
  double fUpwindQuad[27] = {0.0};
  double fUpwind[20] = {0.0};;
  double Ghat[20] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]); 
  alphaDrSurf[7] = -0.7071067811865475*(2.0*sumNuUy[4]+((-2.0*w[3])-1.0*dxv[3])*nuSum[4]); 
  alphaDrSurf[8] = -0.7071067811865475*(2.0*sumNuUy[5]+((-2.0*w[3])-1.0*dxv[3])*nuSum[5]); 
  alphaDrSurf[11] = -0.7071067811865475*(2.0*sumNuUy[6]+((-2.0*w[3])-1.0*dxv[3])*nuSum[6]); 
  alphaDrSurf[12] = -0.7071067811865475*(2.0*sumNuUy[7]+((-2.0*w[3])-1.0*dxv[3])*nuSum[7]); 

  if ((-0.4242640687119282*alphaDrSurf[12])-0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_r(fSkin); 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_r(fSkin); 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_r(fSkin); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_l(fEdge); 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_l(fEdge); 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_l(fEdge); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])-0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_r(fSkin); 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_r(fSkin); 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_r(fSkin); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_l(fEdge); 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_l(fEdge); 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_l(fEdge); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])-0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_r(fSkin); 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_r(fSkin); 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_r(fSkin); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_l(fEdge); 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_l(fEdge); 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_l(fEdge); 
  } 
  if (0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_r(fSkin); 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_r(fSkin); 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_r(fSkin); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_l(fEdge); 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_l(fEdge); 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_l(fEdge); 
  } 
  if (0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_r(fSkin); 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_r(fSkin); 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_r(fSkin); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_l(fEdge); 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_l(fEdge); 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_l(fEdge); 
  } 
  if (0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_r(fSkin); 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_r(fSkin); 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_r(fSkin); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_l(fEdge); 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_l(fEdge); 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_l(fEdge); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])+0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_r(fSkin); 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_r(fSkin); 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_r(fSkin); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_l(fEdge); 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_l(fEdge); 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_l(fEdge); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])+0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_r(fSkin); 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_r(fSkin); 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_r(fSkin); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_l(fEdge); 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_l(fEdge); 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_l(fEdge); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])+0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_r(fSkin); 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_r(fSkin); 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_r(fSkin); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_l(fEdge); 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_l(fEdge); 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_l(fEdge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alphaDrSurf[8]*fUpwind[12]+0.3535533905932737*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[4]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.3162277660168379*alphaDrSurf[4]*fUpwind[12]+0.3162277660168379*fUpwind[4]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*fUpwind[7]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.3535533905932737*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[4] = 0.2828427124746191*alphaDrSurf[11]*fUpwind[12]+0.3162277660168379*alphaDrSurf[2]*fUpwind[12]+0.2828427124746191*fUpwind[11]*alphaDrSurf[12]+0.3162277660168379*fUpwind[2]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[1]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[5] = 0.3535533905932737*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[4]*fUpwind[17]+0.3535533905932737*alphaDrSurf[12]*fUpwind[14]+0.3162277660168379*alphaDrSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[6] = 0.3162277660168379*alphaDrSurf[4]*fUpwind[18]+0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.3162277660168379*alphaDrSurf[2]*fUpwind[14]+0.3535533905932737*alphaDrSurf[11]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[7] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[12]+0.2258769757263128*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[11]+0.2258769757263128*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[8] = 0.2258769757263128*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[11]*fUpwind[11]+0.2258769757263128*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alphaDrSurf[4]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  Ghat[10] = 0.282842712474619*alphaDrSurf[11]*fUpwind[18]+0.3162277660168379*alphaDrSurf[2]*fUpwind[18]+0.282842712474619*alphaDrSurf[12]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.3162277660168379*alphaDrSurf[4]*fUpwind[14]+0.3162277660168379*alphaDrSurf[4]*fUpwind[13]+0.3162277660168379*fUpwind[6]*alphaDrSurf[12]+0.3162277660168379*fUpwind[5]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]*fUpwind[10]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 
  Ghat[11] = 0.2828427124746191*alphaDrSurf[4]*fUpwind[12]+0.2828427124746191*fUpwind[4]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[8]*fUpwind[11]+0.2258769757263128*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3162277660168379*fUpwind[8]*alphaDrSurf[11]+0.2258769757263128*fUpwind[7]*alphaDrSurf[11]+0.3535533905932737*fUpwind[0]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[4]+0.3162277660168379*fUpwind[1]*alphaDrSurf[4]; 
  Ghat[12] = 0.2258769757263128*alphaDrSurf[8]*fUpwind[12]+0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.2258769757263128*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*fUpwind[7]*alphaDrSurf[12]+0.3535533905932737*fUpwind[0]*alphaDrSurf[12]+0.2828427124746191*alphaDrSurf[4]*fUpwind[11]+0.2828427124746191*fUpwind[4]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[2]*fUpwind[4]+0.3162277660168379*fUpwind[2]*alphaDrSurf[4]; 
  Ghat[13] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[18]+0.2258769757263128*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[2]*fUpwind[17]+0.2258769757263128*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[6]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[5]; 
  Ghat[14] = 0.2258769757263128*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3162277660168379*alphaDrSurf[11]*fUpwind[17]+0.2258769757263128*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]+0.3535533905932737*fUpwind[5]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[2]*fUpwind[6]; 
  Ghat[15] = 0.3162277660168379*alphaDrSurf[11]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[19]+0.3535533905932737*alphaDrSurf[4]*fUpwind[16]+0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  Ghat[16] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3162277660168379*alphaDrSurf[8]*fUpwind[16]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]+0.3535533905932737*alphaDrSurf[4]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[9]; 
  Ghat[17] = 0.2828427124746191*alphaDrSurf[4]*fUpwind[18]+0.3162277660168379*alphaDrSurf[8]*fUpwind[17]+0.2258769757263128*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[11]*fUpwind[14]+0.2258769757263128*alphaDrSurf[11]*fUpwind[13]+0.3535533905932737*alphaDrSurf[2]*fUpwind[13]+0.282842712474619*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*fUpwind[3]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[5]; 
  Ghat[18] = 0.2258769757263128*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.2828427124746191*alphaDrSurf[4]*fUpwind[17]+0.2258769757263128*alphaDrSurf[12]*fUpwind[14]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]+0.3162277660168379*alphaDrSurf[12]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[12]+0.282842712474619*fUpwind[10]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[6]; 
  Ghat[19] = 0.3162277660168379*alphaDrSurf[8]*fUpwind[19]+0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3162277660168379*alphaDrSurf[12]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]+0.3162277660168379*alphaDrSurf[11]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[15]+0.3535533905932737*alphaDrSurf[4]*fUpwind[9]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += 1.58113883008419*Ghat[0]*rdv2; 
  out[15] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[16] += 1.224744871391589*Ghat[4]*rdv2; 
  out[17] += 1.224744871391589*Ghat[5]*rdv2; 
  out[18] += 1.224744871391589*Ghat[6]*rdv2; 
  out[19] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[20] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[21] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[22] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[23] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[24] += 0.7071067811865475*Ghat[16]*rdv2; 
  out[25] += 1.224744871391589*Ghat[7]*rdv2; 
  out[26] += 1.224744871391589*Ghat[8]*rdv2; 
  out[27] += 1.224744871391589*Ghat[9]*rdv2; 
  out[28] += 1.58113883008419*Ghat[1]*rdv2; 
  out[29] += 1.58113883008419*Ghat[2]*rdv2; 
  out[30] += 1.58113883008419*Ghat[3]*rdv2; 
  out[31] += 1.224744871391589*Ghat[10]*rdv2; 
  out[32] += 0.7071067811865475*Ghat[17]*rdv2; 
  out[33] += 0.7071067811865475*Ghat[18]*rdv2; 
  out[34] += 0.7071067811865475*Ghat[19]*rdv2; 
  out[35] += 1.224744871391589*Ghat[11]*rdv2; 
  out[36] += 1.224744871391589*Ghat[12]*rdv2; 
  out[37] += 1.224744871391589*Ghat[13]*rdv2; 
  out[38] += 1.224744871391589*Ghat[14]*rdv2; 
  out[39] += 1.224744871391589*Ghat[15]*rdv2; 
  out[40] += 1.224744871391589*Ghat[16]*rdv2; 
  out[41] += 1.58113883008419*Ghat[4]*rdv2; 
  out[42] += 1.58113883008419*Ghat[5]*rdv2; 
  out[43] += 1.58113883008419*Ghat[6]*rdv2; 
  out[44] += 1.224744871391589*Ghat[17]*rdv2; 
  out[45] += 1.224744871391589*Ghat[18]*rdv2; 
  out[46] += 1.224744871391589*Ghat[19]*rdv2; 
  out[47] += 1.58113883008419*Ghat[10]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]-1.0*dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]-1.0*dxv[3]*nuSum[3]); 
  alphaDrSurf[7] = -0.7071067811865475*(2.0*sumNuUy[4]+(dxv[3]-2.0*w[3])*nuSum[4]); 
  alphaDrSurf[8] = -0.7071067811865475*(2.0*sumNuUy[5]+(dxv[3]-2.0*w[3])*nuSum[5]); 
  alphaDrSurf[11] = -0.7071067811865475*(2.0*sumNuUy[6]+(dxv[3]-2.0*w[3])*nuSum[6]); 
  alphaDrSurf[12] = -0.7071067811865475*(2.0*sumNuUy[7]+(dxv[3]-2.0*w[3])*nuSum[7]); 

  if ((-0.4242640687119282*alphaDrSurf[12])-0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_r(fEdge); 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_r(fEdge); 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_r(fEdge); 
  } else { 
    fUpwindQuad[0] = ser_4x_p2_surfx4_eval_quad_node_0_l(fSkin); 
    fUpwindQuad[1] = ser_4x_p2_surfx4_eval_quad_node_1_l(fSkin); 
    fUpwindQuad[2] = ser_4x_p2_surfx4_eval_quad_node_2_l(fSkin); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])-0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_r(fEdge); 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_r(fEdge); 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_r(fEdge); 
  } else { 
    fUpwindQuad[3] = ser_4x_p2_surfx4_eval_quad_node_3_l(fSkin); 
    fUpwindQuad[4] = ser_4x_p2_surfx4_eval_quad_node_4_l(fSkin); 
    fUpwindQuad[5] = ser_4x_p2_surfx4_eval_quad_node_5_l(fSkin); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])-0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]+0.6363961030678927*alphaDrSurf[4]-0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_r(fEdge); 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_r(fEdge); 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_r(fEdge); 
  } else { 
    fUpwindQuad[6] = ser_4x_p2_surfx4_eval_quad_node_6_l(fSkin); 
    fUpwindQuad[7] = ser_4x_p2_surfx4_eval_quad_node_7_l(fSkin); 
    fUpwindQuad[8] = ser_4x_p2_surfx4_eval_quad_node_8_l(fSkin); 
  } 
  if (0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_r(fEdge); 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_r(fEdge); 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_r(fEdge); 
  } else { 
    fUpwindQuad[9] = ser_4x_p2_surfx4_eval_quad_node_9_l(fSkin); 
    fUpwindQuad[10] = ser_4x_p2_surfx4_eval_quad_node_10_l(fSkin); 
    fUpwindQuad[11] = ser_4x_p2_surfx4_eval_quad_node_11_l(fSkin); 
  } 
  if (0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_r(fEdge); 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_r(fEdge); 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_r(fEdge); 
  } else { 
    fUpwindQuad[12] = ser_4x_p2_surfx4_eval_quad_node_12_l(fSkin); 
    fUpwindQuad[13] = ser_4x_p2_surfx4_eval_quad_node_13_l(fSkin); 
    fUpwindQuad[14] = ser_4x_p2_surfx4_eval_quad_node_14_l(fSkin); 
  } 
  if (0.5303300858899102*alphaDrSurf[12]-0.3952847075210471*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_r(fEdge); 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_r(fEdge); 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_r(fEdge); 
  } else { 
    fUpwindQuad[15] = ser_4x_p2_surfx4_eval_quad_node_15_l(fSkin); 
    fUpwindQuad[16] = ser_4x_p2_surfx4_eval_quad_node_16_l(fSkin); 
    fUpwindQuad[17] = ser_4x_p2_surfx4_eval_quad_node_17_l(fSkin); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])+0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_r(fEdge); 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_r(fEdge); 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_r(fEdge); 
  } else { 
    fUpwindQuad[18] = ser_4x_p2_surfx4_eval_quad_node_18_l(fSkin); 
    fUpwindQuad[19] = ser_4x_p2_surfx4_eval_quad_node_19_l(fSkin); 
    fUpwindQuad[20] = ser_4x_p2_surfx4_eval_quad_node_20_l(fSkin); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])+0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_r(fEdge); 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_r(fEdge); 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_r(fEdge); 
  } else { 
    fUpwindQuad[21] = ser_4x_p2_surfx4_eval_quad_node_21_l(fSkin); 
    fUpwindQuad[22] = ser_4x_p2_surfx4_eval_quad_node_22_l(fSkin); 
    fUpwindQuad[23] = ser_4x_p2_surfx4_eval_quad_node_23_l(fSkin); 
  } 
  if ((-0.4242640687119282*alphaDrSurf[12])+0.4242640687119286*alphaDrSurf[11]+0.3162277660168378*alphaDrSurf[8]+0.3162277660168378*alphaDrSurf[7]-0.6363961030678927*alphaDrSurf[4]+0.4743416490252568*alphaDrSurf[2]-0.4743416490252568*alphaDrSurf[1]+0.3535533905932734*alphaDrSurf[0] < 0) { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_r(fEdge); 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_r(fEdge); 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_r(fEdge); 
  } else { 
    fUpwindQuad[24] = ser_4x_p2_surfx4_eval_quad_node_24_l(fSkin); 
    fUpwindQuad[25] = ser_4x_p2_surfx4_eval_quad_node_25_l(fSkin); 
    fUpwindQuad[26] = ser_4x_p2_surfx4_eval_quad_node_26_l(fSkin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  ser_4x_p2_upwind_quad_to_modal(fUpwindQuad, fUpwind); 

  Ghat[0] = 0.3535533905932737*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = 0.3535533905932737*alphaDrSurf[8]*fUpwind[12]+0.3535533905932737*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[4]*fUpwind[11]+0.3162277660168379*fUpwind[4]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[7]+0.3162277660168379*fUpwind[1]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = 0.3162277660168379*alphaDrSurf[4]*fUpwind[12]+0.3162277660168379*fUpwind[4]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*fUpwind[7]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[8]+0.3162277660168379*fUpwind[2]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = 0.3535533905932737*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[4] = 0.2828427124746191*alphaDrSurf[11]*fUpwind[12]+0.3162277660168379*alphaDrSurf[2]*fUpwind[12]+0.2828427124746191*fUpwind[11]*alphaDrSurf[12]+0.3162277660168379*fUpwind[2]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[1]*fUpwind[11]+0.3162277660168379*fUpwind[1]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[8]+0.3162277660168379*fUpwind[4]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[7]+0.3162277660168379*fUpwind[4]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[5] = 0.3535533905932737*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[4]*fUpwind[17]+0.3535533905932737*alphaDrSurf[12]*fUpwind[14]+0.3162277660168379*alphaDrSurf[1]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[10]+0.3162277660168379*fUpwind[5]*alphaDrSurf[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[6] = 0.3162277660168379*alphaDrSurf[4]*fUpwind[18]+0.3535533905932737*alphaDrSurf[7]*fUpwind[17]+0.3162277660168379*alphaDrSurf[2]*fUpwind[14]+0.3535533905932737*alphaDrSurf[11]*fUpwind[13]+0.3162277660168379*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[10]+0.3162277660168379*fUpwind[6]*alphaDrSurf[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[7] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[12]+0.2258769757263128*alphaDrSurf[11]*fUpwind[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[11]+0.3535533905932737*fUpwind[2]*alphaDrSurf[11]+0.2258769757263128*alphaDrSurf[7]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*fUpwind[0]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[1]*fUpwind[1]; 
  Ghat[8] = 0.2258769757263128*alphaDrSurf[12]*fUpwind[12]+0.3535533905932737*alphaDrSurf[1]*fUpwind[12]+0.3535533905932737*fUpwind[1]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[11]*fUpwind[11]+0.2258769757263128*alphaDrSurf[8]*fUpwind[8]+0.3535533905932737*alphaDrSurf[0]*fUpwind[8]+0.3535533905932737*fUpwind[0]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[4]+0.3162277660168379*alphaDrSurf[2]*fUpwind[2]; 
  Ghat[9] = 0.3535533905932737*alphaDrSurf[4]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[9]; 
  Ghat[10] = 0.282842712474619*alphaDrSurf[11]*fUpwind[18]+0.3162277660168379*alphaDrSurf[2]*fUpwind[18]+0.282842712474619*alphaDrSurf[12]*fUpwind[17]+0.3162277660168379*alphaDrSurf[1]*fUpwind[17]+0.3162277660168379*alphaDrSurf[4]*fUpwind[14]+0.3162277660168379*alphaDrSurf[4]*fUpwind[13]+0.3162277660168379*fUpwind[6]*alphaDrSurf[12]+0.3162277660168379*fUpwind[5]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[8]*fUpwind[10]+0.3162277660168379*alphaDrSurf[7]*fUpwind[10]+0.3535533905932737*alphaDrSurf[0]*fUpwind[10]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 
  Ghat[11] = 0.2828427124746191*alphaDrSurf[4]*fUpwind[12]+0.2828427124746191*fUpwind[4]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[8]*fUpwind[11]+0.2258769757263128*alphaDrSurf[7]*fUpwind[11]+0.3535533905932737*alphaDrSurf[0]*fUpwind[11]+0.3162277660168379*fUpwind[8]*alphaDrSurf[11]+0.2258769757263128*fUpwind[7]*alphaDrSurf[11]+0.3535533905932737*fUpwind[0]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*fUpwind[2]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[4]+0.3162277660168379*fUpwind[1]*alphaDrSurf[4]; 
  Ghat[12] = 0.2258769757263128*alphaDrSurf[8]*fUpwind[12]+0.3162277660168379*alphaDrSurf[7]*fUpwind[12]+0.3535533905932737*alphaDrSurf[0]*fUpwind[12]+0.2258769757263128*fUpwind[8]*alphaDrSurf[12]+0.3162277660168379*fUpwind[7]*alphaDrSurf[12]+0.3535533905932737*fUpwind[0]*alphaDrSurf[12]+0.2828427124746191*alphaDrSurf[4]*fUpwind[11]+0.2828427124746191*fUpwind[4]*alphaDrSurf[11]+0.3535533905932737*alphaDrSurf[1]*fUpwind[8]+0.3535533905932737*fUpwind[1]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[2]*fUpwind[4]+0.3162277660168379*fUpwind[2]*alphaDrSurf[4]; 
  Ghat[13] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[18]+0.2258769757263128*alphaDrSurf[11]*fUpwind[17]+0.3535533905932737*alphaDrSurf[2]*fUpwind[17]+0.2258769757263128*alphaDrSurf[7]*fUpwind[13]+0.3535533905932737*alphaDrSurf[0]*fUpwind[13]+0.3535533905932737*fUpwind[6]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[1]*fUpwind[5]; 
  Ghat[14] = 0.2258769757263128*alphaDrSurf[12]*fUpwind[18]+0.3535533905932737*alphaDrSurf[1]*fUpwind[18]+0.3162277660168379*alphaDrSurf[11]*fUpwind[17]+0.2258769757263128*alphaDrSurf[8]*fUpwind[14]+0.3535533905932737*alphaDrSurf[0]*fUpwind[14]+0.3535533905932737*fUpwind[5]*alphaDrSurf[12]+0.3162277660168379*alphaDrSurf[4]*fUpwind[10]+0.3535533905932737*fUpwind[3]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[2]*fUpwind[6]; 
  Ghat[15] = 0.3162277660168379*alphaDrSurf[11]*fUpwind[19]+0.3535533905932737*alphaDrSurf[2]*fUpwind[19]+0.3535533905932737*alphaDrSurf[4]*fUpwind[16]+0.3162277660168379*alphaDrSurf[7]*fUpwind[15]+0.3535533905932737*alphaDrSurf[0]*fUpwind[15]+0.3535533905932737*alphaDrSurf[1]*fUpwind[9]; 
  Ghat[16] = 0.3162277660168379*alphaDrSurf[12]*fUpwind[19]+0.3535533905932737*alphaDrSurf[1]*fUpwind[19]+0.3162277660168379*alphaDrSurf[8]*fUpwind[16]+0.3535533905932737*alphaDrSurf[0]*fUpwind[16]+0.3535533905932737*alphaDrSurf[4]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[9]; 
  Ghat[17] = 0.2828427124746191*alphaDrSurf[4]*fUpwind[18]+0.3162277660168379*alphaDrSurf[8]*fUpwind[17]+0.2258769757263128*alphaDrSurf[7]*fUpwind[17]+0.3535533905932737*alphaDrSurf[0]*fUpwind[17]+0.3162277660168379*alphaDrSurf[11]*fUpwind[14]+0.2258769757263128*alphaDrSurf[11]*fUpwind[13]+0.3535533905932737*alphaDrSurf[2]*fUpwind[13]+0.282842712474619*fUpwind[10]*alphaDrSurf[12]+0.3535533905932737*fUpwind[3]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[1]*fUpwind[10]+0.3535533905932737*fUpwind[6]*alphaDrSurf[7]+0.3162277660168379*alphaDrSurf[4]*fUpwind[5]; 
  Ghat[18] = 0.2258769757263128*alphaDrSurf[8]*fUpwind[18]+0.3162277660168379*alphaDrSurf[7]*fUpwind[18]+0.3535533905932737*alphaDrSurf[0]*fUpwind[18]+0.2828427124746191*alphaDrSurf[4]*fUpwind[17]+0.2258769757263128*alphaDrSurf[12]*fUpwind[14]+0.3535533905932737*alphaDrSurf[1]*fUpwind[14]+0.3162277660168379*alphaDrSurf[12]*fUpwind[13]+0.3535533905932737*fUpwind[3]*alphaDrSurf[12]+0.282842712474619*fUpwind[10]*alphaDrSurf[11]+0.3162277660168379*alphaDrSurf[2]*fUpwind[10]+0.3535533905932737*fUpwind[5]*alphaDrSurf[8]+0.3162277660168379*alphaDrSurf[4]*fUpwind[6]; 
  Ghat[19] = 0.3162277660168379*alphaDrSurf[8]*fUpwind[19]+0.3162277660168379*alphaDrSurf[7]*fUpwind[19]+0.3535533905932737*alphaDrSurf[0]*fUpwind[19]+0.3162277660168379*alphaDrSurf[12]*fUpwind[16]+0.3535533905932737*alphaDrSurf[1]*fUpwind[16]+0.3162277660168379*alphaDrSurf[11]*fUpwind[15]+0.3535533905932737*alphaDrSurf[2]*fUpwind[15]+0.3535533905932737*alphaDrSurf[4]*fUpwind[9]; 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 1.224744871391589*Ghat[0]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += -0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += -1.58113883008419*Ghat[0]*rdv2; 
  out[15] += -0.7071067811865475*Ghat[10]*rdv2; 
  out[16] += 1.224744871391589*Ghat[4]*rdv2; 
  out[17] += 1.224744871391589*Ghat[5]*rdv2; 
  out[18] += 1.224744871391589*Ghat[6]*rdv2; 
  out[19] += -0.7071067811865475*Ghat[11]*rdv2; 
  out[20] += -0.7071067811865475*Ghat[12]*rdv2; 
  out[21] += -0.7071067811865475*Ghat[13]*rdv2; 
  out[22] += -0.7071067811865475*Ghat[14]*rdv2; 
  out[23] += -0.7071067811865475*Ghat[15]*rdv2; 
  out[24] += -0.7071067811865475*Ghat[16]*rdv2; 
  out[25] += 1.224744871391589*Ghat[7]*rdv2; 
  out[26] += 1.224744871391589*Ghat[8]*rdv2; 
  out[27] += 1.224744871391589*Ghat[9]*rdv2; 
  out[28] += -1.58113883008419*Ghat[1]*rdv2; 
  out[29] += -1.58113883008419*Ghat[2]*rdv2; 
  out[30] += -1.58113883008419*Ghat[3]*rdv2; 
  out[31] += 1.224744871391589*Ghat[10]*rdv2; 
  out[32] += -0.7071067811865475*Ghat[17]*rdv2; 
  out[33] += -0.7071067811865475*Ghat[18]*rdv2; 
  out[34] += -0.7071067811865475*Ghat[19]*rdv2; 
  out[35] += 1.224744871391589*Ghat[11]*rdv2; 
  out[36] += 1.224744871391589*Ghat[12]*rdv2; 
  out[37] += 1.224744871391589*Ghat[13]*rdv2; 
  out[38] += 1.224744871391589*Ghat[14]*rdv2; 
  out[39] += 1.224744871391589*Ghat[15]*rdv2; 
  out[40] += 1.224744871391589*Ghat[16]*rdv2; 
  out[41] += -1.58113883008419*Ghat[4]*rdv2; 
  out[42] += -1.58113883008419*Ghat[5]*rdv2; 
  out[43] += -1.58113883008419*Ghat[6]*rdv2; 
  out[44] += 1.224744871391589*Ghat[17]*rdv2; 
  out[45] += 1.224744871391589*Ghat[18]*rdv2; 
  out[46] += 1.224744871391589*Ghat[19]*rdv2; 
  out[47] += -1.58113883008419*Ghat[10]*rdv2; 

  } 

  return 0.;

} 
