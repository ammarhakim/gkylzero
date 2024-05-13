#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_2x2v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_hyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_extem_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:     potentials, including external (scaled by appropriate factors).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 

  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 

  const double *phi = &field[0]; 

  const double *A0 = &field[4]; 
  const double *A1 = &field[8]; 

  double alpha[12] = {0.0}; 

  alpha[0] = 2.449489742783178*A0[2]*dx11*wv1-2.449489742783178*A1[1]*dx10*wv1-2.449489742783178*phi[2]*dx11; 
  alpha[1] = 2.449489742783178*A0[3]*dx11*wv1-2.449489742783178*phi[3]*dx11; 
  alpha[2] = -2.449489742783178*A1[3]*dx10*wv1; 
  alpha[3] = 0.7071067811865475*A0[2]*dv1*dx11-0.7071067811865475*A1[1]*dv1*dx10; 
  alpha[5] = 0.7071067811865475*A0[3]*dv1*dx11; 
  alpha[6] = -0.7071067811865475*A1[3]*dv1*dx10; 

  double fUpwindQuad_l[12] = {0.0};
  double fUpwindQuad_r[12] = {0.0};
  double fUpwind_l[12] = {0.0};;
  double fUpwind_r[12] = {0.0};
  double Ghat_l[12] = {0.0}; 
  double Ghat_r[12] = {0.0}; 

  if (0.4743416490252568*(alpha[6]+alpha[5])-0.4743416490252568*alpha[3]-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = hyb_2x2v_p1_surfx4_eval_quad_node_0_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3535533905932737*(alpha[2]+alpha[1]) > 0) { 
    fUpwindQuad_l[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = hyb_2x2v_p1_surfx4_eval_quad_node_1_l(fr); 
  } 
  if ((-0.4743416490252568*(alpha[6]+alpha[5]))+0.4743416490252568*alpha[3]-0.3535533905932737*(alpha[2]+alpha[1])+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = hyb_2x2v_p1_surfx4_eval_quad_node_2_l(fr); 
  } 
  if ((-0.4743416490252568*alpha[6])+0.4743416490252568*alpha[5]-0.4743416490252568*alpha[3]+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = hyb_2x2v_p1_surfx4_eval_quad_node_3_l(fr); 
  } 
  if (0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = hyb_2x2v_p1_surfx4_eval_quad_node_4_l(fr); 
  } 
  if (0.4743416490252568*alpha[6]-0.4743416490252568*alpha[5]+0.4743416490252568*alpha[3]+0.3535533905932737*alpha[2]-0.3535533905932737*alpha[1]+0.3535533905932737*alpha[0] > 0) { 
    fUpwindQuad_l[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = hyb_2x2v_p1_surfx4_eval_quad_node_5_l(fr); 
  } 
  if (0.4743416490252568*alpha[6]-0.4743416490252568*(alpha[5]+alpha[3])-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = hyb_2x2v_p1_surfx4_eval_quad_node_6_l(fr); 
  } 
  if (0.3535533905932737*(alpha[1]+alpha[0])-0.3535533905932737*alpha[2] > 0) { 
    fUpwindQuad_l[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = hyb_2x2v_p1_surfx4_eval_quad_node_7_l(fr); 
  } 
  if ((-0.4743416490252568*alpha[6])+0.4743416490252568*(alpha[5]+alpha[3])-0.3535533905932737*alpha[2]+0.3535533905932737*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = hyb_2x2v_p1_surfx4_eval_quad_node_8_l(fr); 
  } 
  if (0.3535533905932737*(alpha[2]+alpha[1]+alpha[0])-0.4743416490252568*(alpha[6]+alpha[5]+alpha[3]) > 0) { 
    fUpwindQuad_l[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = hyb_2x2v_p1_surfx4_eval_quad_node_9_l(fr); 
  } 
  if (0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = hyb_2x2v_p1_surfx4_eval_quad_node_10_l(fr); 
  } 
  if (0.4743416490252568*(alpha[6]+alpha[5]+alpha[3])+0.3535533905932737*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = hyb_2x2v_p1_surfx4_eval_quad_node_11_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alpha[6]*fUpwind_l[6]+0.3535533905932737*alpha[5]*fUpwind_l[5]+0.3535533905932737*alpha[3]*fUpwind_l[3]+0.3535533905932737*alpha[2]*fUpwind_l[2]+0.3535533905932737*alpha[1]*fUpwind_l[1]+0.3535533905932737*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alpha[6]*fUpwind_l[7]+0.3535533905932737*alpha[3]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[4]+0.3535533905932737*alpha[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.3535533905932737*alpha[5]*fUpwind_l[7]+0.3535533905932737*alpha[3]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_l[4]+0.3535533905932737*alpha[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.3162277660168379*alpha[6]*fUpwind_l[10]+0.3162277660168379*alpha[5]*fUpwind_l[9]+0.3162277660168379*alpha[3]*fUpwind_l[8]+0.3535533905932737*alpha[2]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[0]*alpha[3]; 
  Ghat_l[4] = 0.3535533905932737*alpha[3]*fUpwind_l[7]+0.3535533905932737*alpha[5]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_l[4]+0.3535533905932737*alpha[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alpha[2]; 
  Ghat_l[5] = 0.3162277660168379*alpha[6]*fUpwind_l[11]+0.3162277660168379*alpha[3]*fUpwind_l[9]+0.3162277660168379*alpha[5]*fUpwind_l[8]+0.3535533905932737*alpha[2]*fUpwind_l[7]+0.3535533905932737*fUpwind_l[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[1]*alpha[3]; 
  Ghat_l[6] = 0.3162277660168379*alpha[5]*fUpwind_l[11]+0.3162277660168379*alpha[3]*fUpwind_l[10]+0.3162277660168379*alpha[6]*fUpwind_l[8]+0.3535533905932737*alpha[1]*fUpwind_l[7]+0.3535533905932737*alpha[0]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[0]*alpha[6]+0.3535533905932737*fUpwind_l[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_l[3]+0.3535533905932737*fUpwind_l[2]*alpha[3]; 
  Ghat_l[7] = 0.3162277660168379*alpha[3]*fUpwind_l[11]+0.3162277660168379*alpha[5]*fUpwind_l[10]+0.3162277660168379*alpha[6]*fUpwind_l[9]+0.3535533905932737*alpha[0]*fUpwind_l[7]+0.3535533905932737*alpha[1]*fUpwind_l[6]+0.3535533905932737*fUpwind_l[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_l[4]; 
  Ghat_l[8] = 0.3535533905932737*alpha[2]*fUpwind_l[10]+0.3535533905932737*alpha[1]*fUpwind_l[9]+0.3535533905932737*alpha[0]*fUpwind_l[8]+0.3162277660168379*alpha[6]*fUpwind_l[6]+0.3162277660168379*alpha[5]*fUpwind_l[5]+0.3162277660168379*alpha[3]*fUpwind_l[3]; 
  Ghat_l[9] = 0.3535533905932737*alpha[2]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[9]+0.3535533905932737*alpha[1]*fUpwind_l[8]+0.3162277660168379*alpha[6]*fUpwind_l[7]+0.3162277660168379*alpha[3]*fUpwind_l[5]+0.3162277660168379*fUpwind_l[3]*alpha[5]; 
  Ghat_l[10] = 0.3535533905932737*alpha[1]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[10]+0.3535533905932737*alpha[2]*fUpwind_l[8]+0.3162277660168379*alpha[5]*fUpwind_l[7]+0.3162277660168379*alpha[3]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[3]*alpha[6]; 
  Ghat_l[11] = 0.3535533905932737*alpha[0]*fUpwind_l[11]+0.3535533905932737*alpha[1]*fUpwind_l[10]+0.3535533905932737*alpha[2]*fUpwind_l[9]+0.3162277660168379*alpha[3]*fUpwind_l[7]+0.3162277660168379*alpha[5]*fUpwind_l[6]+0.3162277660168379*fUpwind_l[5]*alpha[6]; 

  Ghat_r[0] = 0.3535533905932737*alpha[6]*fUpwind_r[6]+0.3535533905932737*alpha[5]*fUpwind_r[5]+0.3535533905932737*alpha[3]*fUpwind_r[3]+0.3535533905932737*alpha[2]*fUpwind_r[2]+0.3535533905932737*alpha[1]*fUpwind_r[1]+0.3535533905932737*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3535533905932737*alpha[6]*fUpwind_r[7]+0.3535533905932737*alpha[3]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[4]+0.3535533905932737*alpha[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.3535533905932737*alpha[5]*fUpwind_r[7]+0.3535533905932737*alpha[3]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[3]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_r[4]+0.3535533905932737*alpha[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.3162277660168379*alpha[6]*fUpwind_r[10]+0.3162277660168379*alpha[5]*fUpwind_r[9]+0.3162277660168379*alpha[3]*fUpwind_r[8]+0.3535533905932737*alpha[2]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[2]*alpha[6]+0.3535533905932737*alpha[1]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[1]*alpha[5]+0.3535533905932737*alpha[0]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[0]*alpha[3]; 
  Ghat_r[4] = 0.3535533905932737*alpha[3]*fUpwind_r[7]+0.3535533905932737*alpha[5]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[5]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_r[4]+0.3535533905932737*alpha[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alpha[2]; 
  Ghat_r[5] = 0.3162277660168379*alpha[6]*fUpwind_r[11]+0.3162277660168379*alpha[3]*fUpwind_r[9]+0.3162277660168379*alpha[5]*fUpwind_r[8]+0.3535533905932737*alpha[2]*fUpwind_r[7]+0.3535533905932737*fUpwind_r[4]*alpha[6]+0.3535533905932737*alpha[0]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[0]*alpha[5]+0.3535533905932737*alpha[1]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[1]*alpha[3]; 
  Ghat_r[6] = 0.3162277660168379*alpha[5]*fUpwind_r[11]+0.3162277660168379*alpha[3]*fUpwind_r[10]+0.3162277660168379*alpha[6]*fUpwind_r[8]+0.3535533905932737*alpha[1]*fUpwind_r[7]+0.3535533905932737*alpha[0]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[0]*alpha[6]+0.3535533905932737*fUpwind_r[4]*alpha[5]+0.3535533905932737*alpha[2]*fUpwind_r[3]+0.3535533905932737*fUpwind_r[2]*alpha[3]; 
  Ghat_r[7] = 0.3162277660168379*alpha[3]*fUpwind_r[11]+0.3162277660168379*alpha[5]*fUpwind_r[10]+0.3162277660168379*alpha[6]*fUpwind_r[9]+0.3535533905932737*alpha[0]*fUpwind_r[7]+0.3535533905932737*alpha[1]*fUpwind_r[6]+0.3535533905932737*fUpwind_r[1]*alpha[6]+0.3535533905932737*alpha[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[2]*alpha[5]+0.3535533905932737*alpha[3]*fUpwind_r[4]; 
  Ghat_r[8] = 0.3535533905932737*alpha[2]*fUpwind_r[10]+0.3535533905932737*alpha[1]*fUpwind_r[9]+0.3535533905932737*alpha[0]*fUpwind_r[8]+0.3162277660168379*alpha[6]*fUpwind_r[6]+0.3162277660168379*alpha[5]*fUpwind_r[5]+0.3162277660168379*alpha[3]*fUpwind_r[3]; 
  Ghat_r[9] = 0.3535533905932737*alpha[2]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[9]+0.3535533905932737*alpha[1]*fUpwind_r[8]+0.3162277660168379*alpha[6]*fUpwind_r[7]+0.3162277660168379*alpha[3]*fUpwind_r[5]+0.3162277660168379*fUpwind_r[3]*alpha[5]; 
  Ghat_r[10] = 0.3535533905932737*alpha[1]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[10]+0.3535533905932737*alpha[2]*fUpwind_r[8]+0.3162277660168379*alpha[5]*fUpwind_r[7]+0.3162277660168379*alpha[3]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[3]*alpha[6]; 
  Ghat_r[11] = 0.3535533905932737*alpha[0]*fUpwind_r[11]+0.3535533905932737*alpha[1]*fUpwind_r[10]+0.3535533905932737*alpha[2]*fUpwind_r[9]+0.3162277660168379*alpha[3]*fUpwind_r[7]+0.3162277660168379*alpha[5]*fUpwind_r[6]+0.3162277660168379*fUpwind_r[5]*alpha[6]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[4] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11; 
  out[8] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[9] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11; 
  out[12] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11; 
  out[16] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv11; 
  out[17] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv11; 
  out[18] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv11; 
  out[19] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv11; 
  out[20] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv11; 
  out[21] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv11; 
  out[22] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv11; 
  out[23] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv11; 
  out[24] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv11; 
  out[25] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv11; 
  out[26] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*dv11; 
  out[27] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*dv11; 
  out[28] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*dv11; 
  out[29] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*dv11; 
  out[30] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*dv11; 
  out[31] += (1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*dv11; 

  return 0.;

} 
