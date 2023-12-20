#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_2x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_nonuniformv_surfvx_2x2v_ser_p1(const double *w, const double *dxv, const double *vcoord, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // vcoord:    Discrete (DG) velocity coordinate.
  // field:     potential (scaled by appropriate factors).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double rdv2 = 2/dxv[2]; 
  const double *phi = &field[0]; 
  const double rdx2 = 2/dxv[0]; 
  const double rdy2 = 2/dxv[1]; 
  double alpha[12] = {0.0}; 

  alpha[0] = -2.449489742783178*phi[1]*rdx2; 
  alpha[2] = -2.449489742783178*phi[3]*rdx2; 

  double fUpwindQuad_l[12] = {0.0};
  double fUpwindQuad_r[12] = {0.0};
  double fUpwind_l[12] = {0.0};
  double fUpwind_r[12] = {0.0};
  double Ghat_l[12] = {0.0}; 
  double Ghat_r[12] = {0.0}; 

  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2] > 0) { 
    fUpwindQuad_l[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_r(fc); 
    fUpwindQuad_l[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_r(fc); 
    fUpwindQuad_l[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = hyb_2x2v_p1_surfx3_eval_quad_node_0_l(fr); 
    fUpwindQuad_l[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = hyb_2x2v_p1_surfx3_eval_quad_node_1_l(fr); 
    fUpwindQuad_l[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = hyb_2x2v_p1_surfx3_eval_quad_node_2_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2] > 0) { 
    fUpwindQuad_l[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_r(fc); 
    fUpwindQuad_l[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_r(fc); 
    fUpwindQuad_l[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = hyb_2x2v_p1_surfx3_eval_quad_node_3_l(fr); 
    fUpwindQuad_l[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = hyb_2x2v_p1_surfx3_eval_quad_node_4_l(fr); 
    fUpwindQuad_l[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = hyb_2x2v_p1_surfx3_eval_quad_node_5_l(fr); 
  } 
  if (0.3535533905932737*alpha[0]-0.3535533905932737*alpha[2] > 0) { 
    fUpwindQuad_l[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_r(fc); 
    fUpwindQuad_l[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_r(fc); 
    fUpwindQuad_l[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = hyb_2x2v_p1_surfx3_eval_quad_node_6_l(fr); 
    fUpwindQuad_l[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = hyb_2x2v_p1_surfx3_eval_quad_node_7_l(fr); 
    fUpwindQuad_l[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = hyb_2x2v_p1_surfx3_eval_quad_node_8_l(fr); 
  } 
  if (0.3535533905932737*(alpha[2]+alpha[0]) > 0) { 
    fUpwindQuad_l[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_r(fc); 
    fUpwindQuad_l[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_r(fc); 
    fUpwindQuad_l[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = hyb_2x2v_p1_surfx3_eval_quad_node_9_l(fr); 
    fUpwindQuad_l[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = hyb_2x2v_p1_surfx3_eval_quad_node_10_l(fr); 
    fUpwindQuad_l[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = hyb_2x2v_p1_surfx3_eval_quad_node_11_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_2x2v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.3535533905932737*alpha[2]*fUpwind_l[2]+0.3535533905932737*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.3535533905932737*alpha[2]*fUpwind_l[4]+0.3535533905932737*alpha[0]*fUpwind_l[1]; 
  Ghat_l[2] = 0.3535533905932737*alpha[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.3535533905932737*alpha[2]*fUpwind_l[6]+0.3535533905932737*alpha[0]*fUpwind_l[3]; 
  Ghat_l[4] = 0.3535533905932737*alpha[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alpha[2]; 
  Ghat_l[5] = 0.3535533905932737*alpha[2]*fUpwind_l[7]+0.3535533905932737*alpha[0]*fUpwind_l[5]; 
  Ghat_l[6] = 0.3535533905932737*alpha[0]*fUpwind_l[6]+0.3535533905932737*alpha[2]*fUpwind_l[3]; 
  Ghat_l[7] = 0.3535533905932737*alpha[0]*fUpwind_l[7]+0.3535533905932737*alpha[2]*fUpwind_l[5]; 
  Ghat_l[8] = 0.3535533905932737*alpha[2]*fUpwind_l[10]+0.3535533905932737*alpha[0]*fUpwind_l[8]; 
  Ghat_l[9] = 0.3535533905932737*alpha[2]*fUpwind_l[11]+0.3535533905932737*alpha[0]*fUpwind_l[9]; 
  Ghat_l[10] = 0.3535533905932737*alpha[0]*fUpwind_l[10]+0.3535533905932737*alpha[2]*fUpwind_l[8]; 
  Ghat_l[11] = 0.3535533905932737*alpha[0]*fUpwind_l[11]+0.3535533905932737*alpha[2]*fUpwind_l[9]; 

  Ghat_r[0] = 0.3535533905932737*alpha[2]*fUpwind_r[2]+0.3535533905932737*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.3535533905932737*alpha[2]*fUpwind_r[4]+0.3535533905932737*alpha[0]*fUpwind_r[1]; 
  Ghat_r[2] = 0.3535533905932737*alpha[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.3535533905932737*alpha[2]*fUpwind_r[6]+0.3535533905932737*alpha[0]*fUpwind_r[3]; 
  Ghat_r[4] = 0.3535533905932737*alpha[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alpha[2]; 
  Ghat_r[5] = 0.3535533905932737*alpha[2]*fUpwind_r[7]+0.3535533905932737*alpha[0]*fUpwind_r[5]; 
  Ghat_r[6] = 0.3535533905932737*alpha[0]*fUpwind_r[6]+0.3535533905932737*alpha[2]*fUpwind_r[3]; 
  Ghat_r[7] = 0.3535533905932737*alpha[0]*fUpwind_r[7]+0.3535533905932737*alpha[2]*fUpwind_r[5]; 
  Ghat_r[8] = 0.3535533905932737*alpha[2]*fUpwind_r[10]+0.3535533905932737*alpha[0]*fUpwind_r[8]; 
  Ghat_r[9] = 0.3535533905932737*alpha[2]*fUpwind_r[11]+0.3535533905932737*alpha[0]*fUpwind_r[9]; 
  Ghat_r[10] = 0.3535533905932737*alpha[0]*fUpwind_r[10]+0.3535533905932737*alpha[2]*fUpwind_r[8]; 
  Ghat_r[11] = 0.3535533905932737*alpha[0]*fUpwind_r[11]+0.3535533905932737*alpha[2]*fUpwind_r[9]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*rdv2; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*rdv2; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*rdv2; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*rdv2; 
  out[6] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdv2; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*rdv2; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*rdv2; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*rdv2; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*rdv2; 
  out[12] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*rdv2; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*rdv2; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*rdv2; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*rdv2; 
  out[16] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*rdv2; 
  out[17] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*rdv2; 
  out[18] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*rdv2; 
  out[19] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*rdv2; 
  out[20] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*rdv2; 
  out[21] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*rdv2; 
  out[22] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*rdv2; 
  out[23] += (1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*rdv2; 
  out[24] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*rdv2; 
  out[25] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*rdv2; 
  out[26] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*rdv2; 
  out[27] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*rdv2; 
  out[28] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*rdv2; 
  out[29] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*rdv2; 
  out[30] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*rdv2; 
  out[31] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*rdv2; 

  return 0.;

} 
