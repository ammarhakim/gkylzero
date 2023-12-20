#include <gkyl_vlasov_poisson_kernels.h> 
#include <gkyl_basis_hyb_2x3v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_hyb_2x3v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double vlasov_poisson_extem_nonuniformv_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *vcoord, const double *field, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // vcoord:    Discrete (DG) velocity coordinate.
  // field:     potentials, including external (scaled by appropriate factors).
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double rdv2 = 2/dxv[2]; 
  const double *phi = &field[0]; 
  const double rdx2 = 2/dxv[0]; 
  const double rdy2 = 2/dxv[1]; 
  const double *A0 = &field[4]; 
  const double *A1 = &field[8]; 
  const double *A2 = &field[12]; 
  double alpha[32] = {0.0}; 

  alpha[0] = (-1.224744871391589*A0[2]*vcoord[20]*rdy2)+1.224744871391589*A2[1]*vcoord[40]*rdx2+1.224744871391589*A1[1]*vcoord[20]*rdx2-3.464101615137754*phi[1]*rdx2; 
  alpha[1] = -1.224744871391589*A0[3]*vcoord[20]*rdy2; 
  alpha[2] = 1.224744871391589*A2[3]*vcoord[40]*rdx2+1.224744871391589*A1[3]*vcoord[20]*rdx2-3.464101615137754*phi[3]*rdx2; 
  alpha[3] = 1.224744871391589*A1[1]*vcoord[22]*rdx2-1.224744871391589*A0[2]*vcoord[22]*rdy2; 
  alpha[4] = 1.224744871391589*A2[1]*vcoord[43]*rdx2; 
  alpha[6] = -1.224744871391589*A0[3]*vcoord[22]*rdy2; 
  alpha[7] = 1.224744871391589*A1[3]*vcoord[22]*rdx2; 
  alpha[9] = 1.224744871391589*A2[3]*vcoord[43]*rdx2; 
  alpha[16] = 1.224744871391589*A1[1]*vcoord[28]*rdx2-1.224744871391589*A0[2]*vcoord[28]*rdy2; 
  alpha[17] = -1.224744871391589*A0[3]*vcoord[28]*rdy2; 
  alpha[18] = 1.224744871391589*A1[3]*vcoord[28]*rdx2; 
  alpha[24] = 1.224744871391589*A2[1]*vcoord[49]*rdx2; 
  alpha[26] = 1.224744871391589*A2[3]*vcoord[49]*rdx2; 

  double fUpwindQuad_l[36] = {0.0};
  double fUpwindQuad_r[36] = {0.0};
  double fUpwind_l[32] = {0.0};;
  double fUpwind_r[32] = {0.0};
  double Ghat_l[32] = {0.0}; 
  double Ghat_r[32] = {0.0}; 

  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]-0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]+0.3354101966249685*(alpha[9]+alpha[7]+alpha[6])-0.3354101966249685*(alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_r(fl); 
    fUpwindQuad_r[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_l(fc); 
    fUpwindQuad_r[0] = hyb_2x3v_p1_surfx3_eval_quad_node_0_l(fr); 
  } 
  if (0.2795084971874738*alpha[26]-0.2795084971874737*alpha[24]-0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]+0.3354101966249685*(alpha[7]+alpha[6])-0.3354101966249685*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_r(fl); 
    fUpwindQuad_r[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_l(fc); 
    fUpwindQuad_r[1] = hyb_2x3v_p1_surfx3_eval_quad_node_1_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]-0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]-0.3354101966249685*alpha[9]+0.3354101966249685*(alpha[7]+alpha[6]+alpha[4])-0.3354101966249685*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_r(fl); 
    fUpwindQuad_r[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_l(fc); 
    fUpwindQuad_r[2] = hyb_2x3v_p1_surfx3_eval_quad_node_2_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]+0.2795084971874738*(alpha[18]+alpha[17])-0.2795084971874737*alpha[16]+0.3354101966249685*alpha[9]-0.3354101966249685*alpha[4]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_r(fl); 
    fUpwindQuad_r[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_l(fc); 
    fUpwindQuad_r[3] = hyb_2x3v_p1_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (0.2795084971874738*alpha[26]-0.2795084971874737*alpha[24]+0.2795084971874738*(alpha[18]+alpha[17])-0.2795084971874737*alpha[16]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_r(fl); 
    fUpwindQuad_r[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_l(fc); 
    fUpwindQuad_r[4] = hyb_2x3v_p1_surfx3_eval_quad_node_4_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]+0.2795084971874738*(alpha[18]+alpha[17])-0.2795084971874737*alpha[16]-0.3354101966249685*alpha[9]+0.3354101966249685*alpha[4]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_r(fl); 
    fUpwindQuad_r[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_l(fc); 
    fUpwindQuad_r[5] = hyb_2x3v_p1_surfx3_eval_quad_node_5_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]-0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]+0.3354101966249685*alpha[9]-0.3354101966249685*(alpha[7]+alpha[6]+alpha[4])+0.3354101966249685*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_r(fl); 
    fUpwindQuad_r[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_l(fc); 
    fUpwindQuad_r[6] = hyb_2x3v_p1_surfx3_eval_quad_node_6_l(fr); 
  } 
  if (0.2795084971874738*alpha[26]-0.2795084971874737*alpha[24]-0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]-0.3354101966249685*(alpha[7]+alpha[6])+0.3354101966249685*alpha[3]-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_r(fl); 
    fUpwindQuad_r[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_l(fc); 
    fUpwindQuad_r[7] = hyb_2x3v_p1_surfx3_eval_quad_node_7_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]-0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]-0.3354101966249685*(alpha[9]+alpha[7]+alpha[6])+0.3354101966249685*(alpha[4]+alpha[3])-0.25*(alpha[2]+alpha[1])+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_r(fl); 
    fUpwindQuad_r[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpwindQuad_l[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_l(fc); 
    fUpwindQuad_r[8] = hyb_2x3v_p1_surfx3_eval_quad_node_8_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]+0.223606797749979*alpha[18]-0.223606797749979*alpha[17]+0.223606797749979*alpha[16]-0.3354101966249685*(alpha[9]+alpha[7])+0.3354101966249685*alpha[6]-0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_r(fl); 
    fUpwindQuad_r[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_r(fc); 
  } else { 
    fUpwindQuad_l[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_l(fc); 
    fUpwindQuad_r[9] = hyb_2x3v_p1_surfx3_eval_quad_node_9_l(fr); 
  } 
  if ((-0.2795084971874738*alpha[26])-0.2795084971874737*alpha[24]+0.223606797749979*alpha[18]-0.223606797749979*alpha[17]+0.223606797749979*alpha[16]-0.3354101966249685*alpha[7]+0.3354101966249685*alpha[6]-0.3354101966249685*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_r(fl); 
    fUpwindQuad_r[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_r(fc); 
  } else { 
    fUpwindQuad_l[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_l(fc); 
    fUpwindQuad_r[10] = hyb_2x3v_p1_surfx3_eval_quad_node_10_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]+0.223606797749979*alpha[18]-0.223606797749979*alpha[17]+0.223606797749979*alpha[16]+0.3354101966249685*alpha[9]-0.3354101966249685*alpha[7]+0.3354101966249685*(alpha[6]+alpha[4])-0.3354101966249685*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_r(fl); 
    fUpwindQuad_r[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_r(fc); 
  } else { 
    fUpwindQuad_l[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_l(fc); 
    fUpwindQuad_r[11] = hyb_2x3v_p1_surfx3_eval_quad_node_11_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]-0.2795084971874738*alpha[18]+0.2795084971874738*alpha[17]-0.2795084971874737*alpha[16]-0.3354101966249685*(alpha[9]+alpha[4])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_r(fl); 
    fUpwindQuad_r[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_r(fc); 
  } else { 
    fUpwindQuad_l[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_l(fc); 
    fUpwindQuad_r[12] = hyb_2x3v_p1_surfx3_eval_quad_node_12_l(fr); 
  } 
  if ((-0.2795084971874738*alpha[26])-0.2795084971874737*alpha[24]-0.2795084971874738*alpha[18]+0.2795084971874738*alpha[17]-0.2795084971874737*alpha[16]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_r(fl); 
    fUpwindQuad_r[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_r(fc); 
  } else { 
    fUpwindQuad_l[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_l(fc); 
    fUpwindQuad_r[13] = hyb_2x3v_p1_surfx3_eval_quad_node_13_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]-0.2795084971874738*alpha[18]+0.2795084971874738*alpha[17]-0.2795084971874737*alpha[16]+0.3354101966249685*(alpha[9]+alpha[4])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_r(fl); 
    fUpwindQuad_r[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_r(fc); 
  } else { 
    fUpwindQuad_l[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_l(fc); 
    fUpwindQuad_r[14] = hyb_2x3v_p1_surfx3_eval_quad_node_14_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]+0.223606797749979*alpha[18]-0.223606797749979*alpha[17]+0.223606797749979*alpha[16]-0.3354101966249685*alpha[9]+0.3354101966249685*alpha[7]-0.3354101966249685*(alpha[6]+alpha[4])+0.3354101966249685*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_r(fl); 
    fUpwindQuad_r[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_r(fc); 
  } else { 
    fUpwindQuad_l[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_l(fc); 
    fUpwindQuad_r[15] = hyb_2x3v_p1_surfx3_eval_quad_node_15_l(fr); 
  } 
  if ((-0.2795084971874738*alpha[26])-0.2795084971874737*alpha[24]+0.223606797749979*alpha[18]-0.223606797749979*alpha[17]+0.223606797749979*alpha[16]+0.3354101966249685*alpha[7]-0.3354101966249685*alpha[6]+0.3354101966249685*alpha[3]+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_r(fl); 
    fUpwindQuad_r[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_r(fc); 
  } else { 
    fUpwindQuad_l[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_l(fc); 
    fUpwindQuad_r[16] = hyb_2x3v_p1_surfx3_eval_quad_node_16_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]+0.223606797749979*alpha[18]-0.223606797749979*alpha[17]+0.223606797749979*alpha[16]+0.3354101966249685*(alpha[9]+alpha[7])-0.3354101966249685*alpha[6]+0.3354101966249685*(alpha[4]+alpha[3])+0.25*alpha[2]-0.25*alpha[1]+0.25*alpha[0] > 0) { 
    fUpwindQuad_l[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_r(fl); 
    fUpwindQuad_r[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_r(fc); 
  } else { 
    fUpwindQuad_l[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_l(fc); 
    fUpwindQuad_r[17] = hyb_2x3v_p1_surfx3_eval_quad_node_17_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]-0.223606797749979*alpha[18]+0.223606797749979*alpha[17]+0.223606797749979*alpha[16]+0.3354101966249685*(alpha[9]+alpha[7])-0.3354101966249685*(alpha[6]+alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_r(fl); 
    fUpwindQuad_r[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_r(fc); 
  } else { 
    fUpwindQuad_l[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_l(fc); 
    fUpwindQuad_r[18] = hyb_2x3v_p1_surfx3_eval_quad_node_18_l(fr); 
  } 
  if (0.2795084971874738*alpha[26]-0.2795084971874737*alpha[24]-0.223606797749979*alpha[18]+0.223606797749979*alpha[17]+0.223606797749979*alpha[16]+0.3354101966249685*alpha[7]-0.3354101966249685*(alpha[6]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_r(fl); 
    fUpwindQuad_r[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_r(fc); 
  } else { 
    fUpwindQuad_l[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_l(fc); 
    fUpwindQuad_r[19] = hyb_2x3v_p1_surfx3_eval_quad_node_19_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]-0.223606797749979*alpha[18]+0.223606797749979*alpha[17]+0.223606797749979*alpha[16]-0.3354101966249685*alpha[9]+0.3354101966249685*alpha[7]-0.3354101966249685*alpha[6]+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_r(fl); 
    fUpwindQuad_r[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_r(fc); 
  } else { 
    fUpwindQuad_l[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_l(fc); 
    fUpwindQuad_r[20] = hyb_2x3v_p1_surfx3_eval_quad_node_20_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]+0.2795084971874738*alpha[18]-0.2795084971874738*alpha[17]-0.2795084971874737*alpha[16]+0.3354101966249685*alpha[9]-0.3354101966249685*alpha[4]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_r(fl); 
    fUpwindQuad_r[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_r(fc); 
  } else { 
    fUpwindQuad_l[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_l(fc); 
    fUpwindQuad_r[21] = hyb_2x3v_p1_surfx3_eval_quad_node_21_l(fr); 
  } 
  if (0.2795084971874738*alpha[26]-0.2795084971874737*alpha[24]+0.2795084971874738*alpha[18]-0.2795084971874738*alpha[17]-0.2795084971874737*alpha[16]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_r(fl); 
    fUpwindQuad_r[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_r(fc); 
  } else { 
    fUpwindQuad_l[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_l(fc); 
    fUpwindQuad_r[22] = hyb_2x3v_p1_surfx3_eval_quad_node_22_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]+0.2795084971874738*alpha[18]-0.2795084971874738*alpha[17]-0.2795084971874737*alpha[16]-0.3354101966249685*alpha[9]+0.3354101966249685*alpha[4]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_r(fl); 
    fUpwindQuad_r[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_r(fc); 
  } else { 
    fUpwindQuad_l[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_l(fc); 
    fUpwindQuad_r[23] = hyb_2x3v_p1_surfx3_eval_quad_node_23_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]-0.223606797749979*alpha[18]+0.223606797749979*alpha[17]+0.223606797749979*alpha[16]+0.3354101966249685*alpha[9]-0.3354101966249685*alpha[7]+0.3354101966249685*alpha[6]-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_r(fl); 
    fUpwindQuad_r[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_r(fc); 
  } else { 
    fUpwindQuad_l[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_l(fc); 
    fUpwindQuad_r[24] = hyb_2x3v_p1_surfx3_eval_quad_node_24_l(fr); 
  } 
  if (0.2795084971874738*alpha[26]-0.2795084971874737*alpha[24]-0.223606797749979*alpha[18]+0.223606797749979*alpha[17]+0.223606797749979*alpha[16]-0.3354101966249685*alpha[7]+0.3354101966249685*(alpha[6]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_r(fl); 
    fUpwindQuad_r[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_r(fc); 
  } else { 
    fUpwindQuad_l[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_l(fc); 
    fUpwindQuad_r[25] = hyb_2x3v_p1_surfx3_eval_quad_node_25_l(fr); 
  } 
  if ((-0.223606797749979*alpha[26])+0.223606797749979*alpha[24]-0.223606797749979*alpha[18]+0.223606797749979*alpha[17]+0.223606797749979*alpha[16]-0.3354101966249685*(alpha[9]+alpha[7])+0.3354101966249685*(alpha[6]+alpha[4]+alpha[3])-0.25*alpha[2]+0.25*(alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_r(fl); 
    fUpwindQuad_r[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_r(fc); 
  } else { 
    fUpwindQuad_l[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_l(fc); 
    fUpwindQuad_r[26] = hyb_2x3v_p1_surfx3_eval_quad_node_26_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]+0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]-0.3354101966249685*(alpha[9]+alpha[7]+alpha[6]+alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_r(fl); 
    fUpwindQuad_r[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_r(fc); 
  } else { 
    fUpwindQuad_l[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_l(fc); 
    fUpwindQuad_r[27] = hyb_2x3v_p1_surfx3_eval_quad_node_27_l(fr); 
  } 
  if ((-0.2795084971874738*alpha[26])-0.2795084971874737*alpha[24]+0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]-0.3354101966249685*(alpha[7]+alpha[6]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_r(fl); 
    fUpwindQuad_r[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_r(fc); 
  } else { 
    fUpwindQuad_l[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_l(fc); 
    fUpwindQuad_r[28] = hyb_2x3v_p1_surfx3_eval_quad_node_28_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]+0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]+0.3354101966249685*alpha[9]-0.3354101966249685*(alpha[7]+alpha[6])+0.3354101966249685*alpha[4]-0.3354101966249685*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_r(fl); 
    fUpwindQuad_r[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_r(fc); 
  } else { 
    fUpwindQuad_l[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_l(fc); 
    fUpwindQuad_r[29] = hyb_2x3v_p1_surfx3_eval_quad_node_29_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]-0.2795084971874738*(alpha[18]+alpha[17])-0.2795084971874737*alpha[16]-0.3354101966249685*(alpha[9]+alpha[4])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_r(fl); 
    fUpwindQuad_r[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_r(fc); 
  } else { 
    fUpwindQuad_l[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_l(fc); 
    fUpwindQuad_r[30] = hyb_2x3v_p1_surfx3_eval_quad_node_30_l(fr); 
  } 
  if ((-0.2795084971874738*alpha[26])-0.2795084971874737*alpha[24]-0.2795084971874738*(alpha[18]+alpha[17])-0.2795084971874737*alpha[16]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_r(fl); 
    fUpwindQuad_r[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_r(fc); 
  } else { 
    fUpwindQuad_l[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_l(fc); 
    fUpwindQuad_r[31] = hyb_2x3v_p1_surfx3_eval_quad_node_31_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]-0.2795084971874738*(alpha[18]+alpha[17])-0.2795084971874737*alpha[16]+0.3354101966249685*(alpha[9]+alpha[4])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_r(fl); 
    fUpwindQuad_r[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_r(fc); 
  } else { 
    fUpwindQuad_l[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_l(fc); 
    fUpwindQuad_r[32] = hyb_2x3v_p1_surfx3_eval_quad_node_32_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]+0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]-0.3354101966249685*alpha[9]+0.3354101966249685*(alpha[7]+alpha[6])-0.3354101966249685*alpha[4]+0.3354101966249685*alpha[3]+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_r(fl); 
    fUpwindQuad_r[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_r(fc); 
  } else { 
    fUpwindQuad_l[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_l(fc); 
    fUpwindQuad_r[33] = hyb_2x3v_p1_surfx3_eval_quad_node_33_l(fr); 
  } 
  if ((-0.2795084971874738*alpha[26])-0.2795084971874737*alpha[24]+0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]+0.3354101966249685*(alpha[7]+alpha[6]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_r(fl); 
    fUpwindQuad_r[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_r(fc); 
  } else { 
    fUpwindQuad_l[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_l(fc); 
    fUpwindQuad_r[34] = hyb_2x3v_p1_surfx3_eval_quad_node_34_l(fr); 
  } 
  if (0.223606797749979*alpha[26]+0.223606797749979*alpha[24]+0.223606797749979*(alpha[18]+alpha[17])+0.223606797749979*alpha[16]+0.3354101966249685*(alpha[9]+alpha[7]+alpha[6]+alpha[4]+alpha[3])+0.25*(alpha[2]+alpha[1]+alpha[0]) > 0) { 
    fUpwindQuad_l[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_r(fl); 
    fUpwindQuad_r[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_r(fc); 
  } else { 
    fUpwindQuad_l[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_l(fc); 
    fUpwindQuad_r[35] = hyb_2x3v_p1_surfx3_eval_quad_node_35_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_l, fUpwind_l); 
  hyb_2x3v_p1_vdir_upwind_quad_to_modal(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] = 0.25*alpha[26]*fUpwind_l[26]+0.25*alpha[24]*fUpwind_l[24]+0.25*alpha[18]*fUpwind_l[18]+0.25*alpha[17]*fUpwind_l[17]+0.25*alpha[16]*fUpwind_l[16]+0.25*alpha[9]*fUpwind_l[9]+0.25*alpha[7]*fUpwind_l[7]+0.25*alpha[6]*fUpwind_l[6]+0.25*alpha[4]*fUpwind_l[4]+0.25*alpha[3]*fUpwind_l[3]+0.25*alpha[2]*fUpwind_l[2]+0.25*alpha[1]*fUpwind_l[1]+0.25*alpha[0]*fUpwind_l[0]; 
  Ghat_l[1] = 0.2500000000000001*alpha[26]*fUpwind_l[28]+0.2500000000000001*alpha[24]*fUpwind_l[25]+0.2500000000000001*alpha[18]*fUpwind_l[20]+0.2500000000000001*alpha[16]*fUpwind_l[17]+0.2500000000000001*fUpwind_l[16]*alpha[17]+0.25*alpha[9]*fUpwind_l[12]+0.25*alpha[7]*fUpwind_l[11]+0.25*alpha[4]*fUpwind_l[8]+0.25*alpha[3]*fUpwind_l[6]+0.25*fUpwind_l[3]*alpha[6]+0.25*alpha[2]*fUpwind_l[5]+0.25*alpha[0]*fUpwind_l[1]+0.25*fUpwind_l[0]*alpha[1]; 
  Ghat_l[2] = 0.2500000000000001*alpha[24]*fUpwind_l[26]+0.2500000000000001*fUpwind_l[24]*alpha[26]+0.2500000000000001*alpha[17]*fUpwind_l[20]+0.2500000000000001*alpha[16]*fUpwind_l[18]+0.2500000000000001*fUpwind_l[16]*alpha[18]+0.25*alpha[6]*fUpwind_l[11]+0.25*alpha[4]*fUpwind_l[9]+0.25*fUpwind_l[4]*alpha[9]+0.25*alpha[3]*fUpwind_l[7]+0.25*fUpwind_l[3]*alpha[7]+0.25*alpha[1]*fUpwind_l[5]+0.25*alpha[0]*fUpwind_l[2]+0.25*fUpwind_l[0]*alpha[2]; 
  Ghat_l[3] = 0.2500000000000001*alpha[26]*fUpwind_l[30]+0.2500000000000001*alpha[24]*fUpwind_l[27]+0.223606797749979*alpha[7]*fUpwind_l[18]+0.223606797749979*fUpwind_l[7]*alpha[18]+0.223606797749979*alpha[6]*fUpwind_l[17]+0.223606797749979*fUpwind_l[6]*alpha[17]+0.223606797749979*alpha[3]*fUpwind_l[16]+0.223606797749979*fUpwind_l[3]*alpha[16]+0.25*alpha[9]*fUpwind_l[14]+0.25*alpha[4]*fUpwind_l[10]+0.25*alpha[2]*fUpwind_l[7]+0.25*fUpwind_l[2]*alpha[7]+0.25*alpha[1]*fUpwind_l[6]+0.25*fUpwind_l[1]*alpha[6]+0.25*alpha[0]*fUpwind_l[3]+0.25*fUpwind_l[0]*alpha[3]; 
  Ghat_l[4] = 0.223606797749979*alpha[9]*fUpwind_l[26]+0.223606797749979*fUpwind_l[9]*alpha[26]+0.223606797749979*alpha[4]*fUpwind_l[24]+0.223606797749979*fUpwind_l[4]*alpha[24]+0.2500000000000001*alpha[18]*fUpwind_l[22]+0.2500000000000001*alpha[17]*fUpwind_l[21]+0.2500000000000001*alpha[16]*fUpwind_l[19]+0.25*alpha[7]*fUpwind_l[14]+0.25*alpha[6]*fUpwind_l[13]+0.25*alpha[3]*fUpwind_l[10]+0.25*alpha[2]*fUpwind_l[9]+0.25*fUpwind_l[2]*alpha[9]+0.25*alpha[1]*fUpwind_l[8]+0.25*alpha[0]*fUpwind_l[4]+0.25*fUpwind_l[0]*alpha[4]; 
  Ghat_l[5] = 0.25*alpha[24]*fUpwind_l[28]+0.25*fUpwind_l[25]*alpha[26]+0.25*alpha[16]*fUpwind_l[20]+0.25*alpha[17]*fUpwind_l[18]+0.25*fUpwind_l[17]*alpha[18]+0.25*alpha[4]*fUpwind_l[12]+0.25*alpha[3]*fUpwind_l[11]+0.25*fUpwind_l[8]*alpha[9]+0.25*alpha[6]*fUpwind_l[7]+0.25*fUpwind_l[6]*alpha[7]+0.25*alpha[0]*fUpwind_l[5]+0.25*alpha[1]*fUpwind_l[2]+0.25*fUpwind_l[1]*alpha[2]; 
  Ghat_l[6] = 0.25*alpha[26]*fUpwind_l[31]+0.25*alpha[24]*fUpwind_l[29]+0.223606797749979*alpha[7]*fUpwind_l[20]+0.223606797749979*fUpwind_l[11]*alpha[18]+0.223606797749979*alpha[3]*fUpwind_l[17]+0.223606797749979*fUpwind_l[3]*alpha[17]+0.223606797749979*alpha[6]*fUpwind_l[16]+0.223606797749979*fUpwind_l[6]*alpha[16]+0.25*alpha[9]*fUpwind_l[15]+0.25*alpha[4]*fUpwind_l[13]+0.25*alpha[2]*fUpwind_l[11]+0.25*fUpwind_l[5]*alpha[7]+0.25*alpha[0]*fUpwind_l[6]+0.25*fUpwind_l[0]*alpha[6]+0.25*alpha[1]*fUpwind_l[3]+0.25*fUpwind_l[1]*alpha[3]; 
  Ghat_l[7] = 0.25*alpha[24]*fUpwind_l[30]+0.25*alpha[26]*fUpwind_l[27]+0.223606797749979*alpha[6]*fUpwind_l[20]+0.223606797749979*alpha[3]*fUpwind_l[18]+0.223606797749979*fUpwind_l[3]*alpha[18]+0.223606797749979*fUpwind_l[11]*alpha[17]+0.223606797749979*alpha[7]*fUpwind_l[16]+0.223606797749979*fUpwind_l[7]*alpha[16]+0.25*alpha[4]*fUpwind_l[14]+0.25*alpha[1]*fUpwind_l[11]+0.25*alpha[9]*fUpwind_l[10]+0.25*alpha[0]*fUpwind_l[7]+0.25*fUpwind_l[0]*alpha[7]+0.25*fUpwind_l[5]*alpha[6]+0.25*alpha[2]*fUpwind_l[3]+0.25*fUpwind_l[2]*alpha[3]; 
  Ghat_l[8] = 0.223606797749979*alpha[9]*fUpwind_l[28]+0.223606797749979*fUpwind_l[12]*alpha[26]+0.223606797749979*alpha[4]*fUpwind_l[25]+0.223606797749979*fUpwind_l[8]*alpha[24]+0.25*alpha[18]*fUpwind_l[23]+0.25*alpha[16]*fUpwind_l[21]+0.25*alpha[17]*fUpwind_l[19]+0.25*alpha[7]*fUpwind_l[15]+0.25*alpha[3]*fUpwind_l[13]+0.25*alpha[2]*fUpwind_l[12]+0.25*alpha[6]*fUpwind_l[10]+0.25*fUpwind_l[5]*alpha[9]+0.25*alpha[0]*fUpwind_l[8]+0.25*alpha[1]*fUpwind_l[4]+0.25*fUpwind_l[1]*alpha[4]; 
  Ghat_l[9] = 0.223606797749979*alpha[4]*fUpwind_l[26]+0.223606797749979*fUpwind_l[4]*alpha[26]+0.223606797749979*alpha[9]*fUpwind_l[24]+0.223606797749979*fUpwind_l[9]*alpha[24]+0.25*alpha[17]*fUpwind_l[23]+0.25*alpha[16]*fUpwind_l[22]+0.25*alpha[18]*fUpwind_l[19]+0.25*alpha[6]*fUpwind_l[15]+0.25*alpha[3]*fUpwind_l[14]+0.25*alpha[1]*fUpwind_l[12]+0.25*alpha[7]*fUpwind_l[10]+0.25*alpha[0]*fUpwind_l[9]+0.25*fUpwind_l[0]*alpha[9]+0.25*alpha[2]*fUpwind_l[4]+0.25*fUpwind_l[2]*alpha[4]; 
  Ghat_l[10] = 0.223606797749979*alpha[9]*fUpwind_l[30]+0.223606797749979*alpha[4]*fUpwind_l[27]+0.223606797749979*fUpwind_l[14]*alpha[26]+0.223606797749979*fUpwind_l[10]*alpha[24]+0.223606797749979*alpha[7]*fUpwind_l[22]+0.223606797749979*alpha[6]*fUpwind_l[21]+0.223606797749979*alpha[3]*fUpwind_l[19]+0.223606797749979*fUpwind_l[14]*alpha[18]+0.223606797749979*fUpwind_l[13]*alpha[17]+0.223606797749979*fUpwind_l[10]*alpha[16]+0.25*alpha[2]*fUpwind_l[14]+0.25*alpha[1]*fUpwind_l[13]+0.25*alpha[0]*fUpwind_l[10]+0.25*alpha[7]*fUpwind_l[9]+0.25*fUpwind_l[7]*alpha[9]+0.25*alpha[6]*fUpwind_l[8]+0.25*alpha[3]*fUpwind_l[4]+0.25*fUpwind_l[3]*alpha[4]; 
  Ghat_l[11] = 0.2500000000000001*alpha[24]*fUpwind_l[31]+0.2500000000000001*alpha[26]*fUpwind_l[29]+0.223606797749979*alpha[3]*fUpwind_l[20]+0.223606797749979*alpha[6]*fUpwind_l[18]+0.223606797749979*fUpwind_l[6]*alpha[18]+0.223606797749979*alpha[7]*fUpwind_l[17]+0.223606797749979*fUpwind_l[7]*alpha[17]+0.223606797749979*fUpwind_l[11]*alpha[16]+0.25*alpha[4]*fUpwind_l[15]+0.25*alpha[9]*fUpwind_l[13]+0.25*alpha[0]*fUpwind_l[11]+0.25*alpha[1]*fUpwind_l[7]+0.25*fUpwind_l[1]*alpha[7]+0.25*alpha[2]*fUpwind_l[6]+0.25*fUpwind_l[2]*alpha[6]+0.25*alpha[3]*fUpwind_l[5]; 
  Ghat_l[12] = 0.223606797749979*alpha[4]*fUpwind_l[28]+0.223606797749979*fUpwind_l[8]*alpha[26]+0.223606797749979*alpha[9]*fUpwind_l[25]+0.223606797749979*fUpwind_l[12]*alpha[24]+0.2500000000000001*alpha[16]*fUpwind_l[23]+0.2500000000000001*alpha[17]*fUpwind_l[22]+0.2500000000000001*alpha[18]*fUpwind_l[21]+0.25*alpha[3]*fUpwind_l[15]+0.25*alpha[6]*fUpwind_l[14]+0.25*alpha[7]*fUpwind_l[13]+0.25*alpha[0]*fUpwind_l[12]+0.25*alpha[1]*fUpwind_l[9]+0.25*fUpwind_l[1]*alpha[9]+0.25*alpha[2]*fUpwind_l[8]+0.25*alpha[4]*fUpwind_l[5]; 
  Ghat_l[13] = 0.223606797749979*alpha[9]*fUpwind_l[31]+0.223606797749979*alpha[4]*fUpwind_l[29]+0.223606797749979*fUpwind_l[15]*alpha[26]+0.223606797749979*fUpwind_l[13]*alpha[24]+0.223606797749979*alpha[7]*fUpwind_l[23]+0.223606797749979*alpha[3]*fUpwind_l[21]+0.223606797749979*alpha[6]*fUpwind_l[19]+0.223606797749979*fUpwind_l[15]*alpha[18]+0.223606797749979*fUpwind_l[10]*alpha[17]+0.223606797749979*fUpwind_l[13]*alpha[16]+0.25*alpha[2]*fUpwind_l[15]+0.25*alpha[0]*fUpwind_l[13]+0.25*alpha[7]*fUpwind_l[12]+0.25*alpha[9]*fUpwind_l[11]+0.25*alpha[1]*fUpwind_l[10]+0.25*alpha[3]*fUpwind_l[8]+0.25*alpha[4]*fUpwind_l[6]+0.25*fUpwind_l[4]*alpha[6]; 
  Ghat_l[14] = 0.223606797749979*alpha[4]*fUpwind_l[30]+0.223606797749979*alpha[9]*fUpwind_l[27]+0.223606797749979*fUpwind_l[10]*alpha[26]+0.223606797749979*fUpwind_l[14]*alpha[24]+0.223606797749979*alpha[6]*fUpwind_l[23]+0.223606797749979*alpha[3]*fUpwind_l[22]+0.223606797749979*alpha[7]*fUpwind_l[19]+0.223606797749979*fUpwind_l[10]*alpha[18]+0.223606797749979*fUpwind_l[15]*alpha[17]+0.223606797749979*fUpwind_l[14]*alpha[16]+0.25*alpha[1]*fUpwind_l[15]+0.25*alpha[0]*fUpwind_l[14]+0.25*alpha[6]*fUpwind_l[12]+0.25*alpha[2]*fUpwind_l[10]+0.25*alpha[3]*fUpwind_l[9]+0.25*fUpwind_l[3]*alpha[9]+0.25*alpha[4]*fUpwind_l[7]+0.25*fUpwind_l[4]*alpha[7]; 
  Ghat_l[15] = 0.223606797749979*alpha[4]*fUpwind_l[31]+0.223606797749979*alpha[9]*fUpwind_l[29]+0.223606797749979*fUpwind_l[13]*alpha[26]+0.223606797749979*fUpwind_l[15]*alpha[24]+0.223606797749979*alpha[3]*fUpwind_l[23]+0.223606797749979*alpha[6]*fUpwind_l[22]+0.223606797749979*alpha[7]*fUpwind_l[21]+0.223606797749979*fUpwind_l[13]*alpha[18]+0.223606797749979*fUpwind_l[14]*alpha[17]+0.223606797749979*fUpwind_l[15]*alpha[16]+0.25*alpha[0]*fUpwind_l[15]+0.25*alpha[1]*fUpwind_l[14]+0.25*alpha[2]*fUpwind_l[13]+0.25*alpha[3]*fUpwind_l[12]+0.25*alpha[4]*fUpwind_l[11]+0.25*alpha[6]*fUpwind_l[9]+0.25*fUpwind_l[6]*alpha[9]+0.25*alpha[7]*fUpwind_l[8]; 
  Ghat_l[16] = 0.25*alpha[9]*fUpwind_l[22]+0.2500000000000001*alpha[4]*fUpwind_l[19]+0.159719141249985*alpha[18]*fUpwind_l[18]+0.2500000000000001*alpha[2]*fUpwind_l[18]+0.2500000000000001*fUpwind_l[2]*alpha[18]+0.159719141249985*alpha[17]*fUpwind_l[17]+0.2500000000000001*alpha[1]*fUpwind_l[17]+0.2500000000000001*fUpwind_l[1]*alpha[17]+0.159719141249985*alpha[16]*fUpwind_l[16]+0.25*alpha[0]*fUpwind_l[16]+0.25*fUpwind_l[0]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_l[7]+0.223606797749979*alpha[6]*fUpwind_l[6]+0.223606797749979*alpha[3]*fUpwind_l[3]; 
  Ghat_l[17] = 0.25*alpha[9]*fUpwind_l[23]+0.2500000000000001*alpha[4]*fUpwind_l[21]+0.159719141249985*alpha[18]*fUpwind_l[20]+0.2500000000000001*alpha[2]*fUpwind_l[20]+0.25*fUpwind_l[5]*alpha[18]+0.159719141249985*alpha[16]*fUpwind_l[17]+0.25*alpha[0]*fUpwind_l[17]+0.159719141249985*fUpwind_l[16]*alpha[17]+0.25*fUpwind_l[0]*alpha[17]+0.2500000000000001*alpha[1]*fUpwind_l[16]+0.2500000000000001*fUpwind_l[1]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_l[11]+0.223606797749979*alpha[3]*fUpwind_l[6]+0.223606797749979*fUpwind_l[3]*alpha[6]; 
  Ghat_l[18] = 0.2500000000000001*alpha[4]*fUpwind_l[22]+0.159719141249985*alpha[17]*fUpwind_l[20]+0.2500000000000001*alpha[1]*fUpwind_l[20]+0.25*alpha[9]*fUpwind_l[19]+0.159719141249985*alpha[16]*fUpwind_l[18]+0.25*alpha[0]*fUpwind_l[18]+0.159719141249985*fUpwind_l[16]*alpha[18]+0.25*fUpwind_l[0]*alpha[18]+0.25*fUpwind_l[5]*alpha[17]+0.2500000000000001*alpha[2]*fUpwind_l[16]+0.2500000000000001*fUpwind_l[2]*alpha[16]+0.223606797749979*alpha[6]*fUpwind_l[11]+0.223606797749979*alpha[3]*fUpwind_l[7]+0.223606797749979*fUpwind_l[3]*alpha[7]; 
  Ghat_l[19] = 0.223606797749979*fUpwind_l[22]*alpha[26]+0.223606797749979*fUpwind_l[19]*alpha[24]+0.159719141249985*alpha[18]*fUpwind_l[22]+0.2500000000000001*alpha[2]*fUpwind_l[22]+0.159719141249985*alpha[17]*fUpwind_l[21]+0.2500000000000001*alpha[1]*fUpwind_l[21]+0.159719141249985*alpha[16]*fUpwind_l[19]+0.25*alpha[0]*fUpwind_l[19]+0.25*alpha[9]*fUpwind_l[18]+0.25*fUpwind_l[9]*alpha[18]+0.25*fUpwind_l[8]*alpha[17]+0.2500000000000001*alpha[4]*fUpwind_l[16]+0.2500000000000001*fUpwind_l[4]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_l[14]+0.223606797749979*alpha[6]*fUpwind_l[13]+0.223606797749979*alpha[3]*fUpwind_l[10]; 
  Ghat_l[20] = 0.2500000000000001*alpha[4]*fUpwind_l[23]+0.25*alpha[9]*fUpwind_l[21]+0.159719141249985*alpha[16]*fUpwind_l[20]+0.25*alpha[0]*fUpwind_l[20]+0.159719141249985*alpha[17]*fUpwind_l[18]+0.2500000000000001*alpha[1]*fUpwind_l[18]+0.159719141249985*fUpwind_l[17]*alpha[18]+0.2500000000000001*fUpwind_l[1]*alpha[18]+0.2500000000000001*alpha[2]*fUpwind_l[17]+0.2500000000000001*fUpwind_l[2]*alpha[17]+0.25*fUpwind_l[5]*alpha[16]+0.223606797749979*alpha[3]*fUpwind_l[11]+0.223606797749979*alpha[6]*fUpwind_l[7]+0.223606797749979*fUpwind_l[6]*alpha[7]; 
  Ghat_l[21] = 0.223606797749979*fUpwind_l[23]*alpha[26]+0.223606797749979*fUpwind_l[21]*alpha[24]+0.159719141249985*alpha[18]*fUpwind_l[23]+0.2500000000000001*alpha[2]*fUpwind_l[23]+0.159719141249985*alpha[16]*fUpwind_l[21]+0.25*alpha[0]*fUpwind_l[21]+0.25*alpha[9]*fUpwind_l[20]+0.159719141249985*alpha[17]*fUpwind_l[19]+0.2500000000000001*alpha[1]*fUpwind_l[19]+0.2500000000000001*fUpwind_l[12]*alpha[18]+0.2500000000000001*alpha[4]*fUpwind_l[17]+0.2500000000000001*fUpwind_l[4]*alpha[17]+0.25*fUpwind_l[8]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_l[15]+0.223606797749979*alpha[3]*fUpwind_l[13]+0.223606797749979*alpha[6]*fUpwind_l[10]; 
  Ghat_l[22] = 0.223606797749979*fUpwind_l[19]*alpha[26]+0.223606797749979*fUpwind_l[22]*alpha[24]+0.159719141249985*alpha[17]*fUpwind_l[23]+0.2500000000000001*alpha[1]*fUpwind_l[23]+0.159719141249985*alpha[16]*fUpwind_l[22]+0.25*alpha[0]*fUpwind_l[22]+0.159719141249985*alpha[18]*fUpwind_l[19]+0.2500000000000001*alpha[2]*fUpwind_l[19]+0.2500000000000001*alpha[4]*fUpwind_l[18]+0.2500000000000001*fUpwind_l[4]*alpha[18]+0.2500000000000001*fUpwind_l[12]*alpha[17]+0.25*alpha[9]*fUpwind_l[16]+0.25*fUpwind_l[9]*alpha[16]+0.223606797749979*alpha[6]*fUpwind_l[15]+0.223606797749979*alpha[3]*fUpwind_l[14]+0.223606797749979*alpha[7]*fUpwind_l[10]; 
  Ghat_l[23] = 0.223606797749979*fUpwind_l[21]*alpha[26]+0.223606797749979*fUpwind_l[23]*alpha[24]+0.159719141249985*alpha[16]*fUpwind_l[23]+0.25*alpha[0]*fUpwind_l[23]+0.159719141249985*alpha[17]*fUpwind_l[22]+0.2500000000000001*alpha[1]*fUpwind_l[22]+0.159719141249985*alpha[18]*fUpwind_l[21]+0.2500000000000001*alpha[2]*fUpwind_l[21]+0.2500000000000001*alpha[4]*fUpwind_l[20]+0.25*fUpwind_l[8]*alpha[18]+0.25*alpha[9]*fUpwind_l[17]+0.25*fUpwind_l[9]*alpha[17]+0.2500000000000001*fUpwind_l[12]*alpha[16]+0.223606797749979*alpha[3]*fUpwind_l[15]+0.223606797749979*alpha[6]*fUpwind_l[14]+0.223606797749979*alpha[7]*fUpwind_l[13]; 
  Ghat_l[24] = 0.25*alpha[7]*fUpwind_l[30]+0.25*alpha[6]*fUpwind_l[29]+0.2500000000000001*alpha[3]*fUpwind_l[27]+0.159719141249985*alpha[26]*fUpwind_l[26]+0.2500000000000001*alpha[2]*fUpwind_l[26]+0.2500000000000001*fUpwind_l[2]*alpha[26]+0.2500000000000001*alpha[1]*fUpwind_l[25]+0.159719141249985*alpha[24]*fUpwind_l[24]+0.25*alpha[0]*fUpwind_l[24]+0.25*fUpwind_l[0]*alpha[24]+0.223606797749979*alpha[9]*fUpwind_l[9]+0.223606797749979*alpha[4]*fUpwind_l[4]; 
  Ghat_l[25] = 0.25*alpha[7]*fUpwind_l[31]+0.2500000000000001*alpha[3]*fUpwind_l[29]+0.159719141249985*alpha[26]*fUpwind_l[28]+0.2500000000000001*alpha[2]*fUpwind_l[28]+0.25*alpha[6]*fUpwind_l[27]+0.25*fUpwind_l[5]*alpha[26]+0.159719141249985*alpha[24]*fUpwind_l[25]+0.25*alpha[0]*fUpwind_l[25]+0.2500000000000001*alpha[1]*fUpwind_l[24]+0.2500000000000001*fUpwind_l[1]*alpha[24]+0.223606797749979*alpha[9]*fUpwind_l[12]+0.223606797749979*alpha[4]*fUpwind_l[8]; 
  Ghat_l[26] = 0.25*alpha[6]*fUpwind_l[31]+0.2500000000000001*alpha[3]*fUpwind_l[30]+0.2500000000000001*alpha[1]*fUpwind_l[28]+0.25*alpha[7]*fUpwind_l[27]+0.159719141249985*alpha[24]*fUpwind_l[26]+0.25*alpha[0]*fUpwind_l[26]+0.159719141249985*fUpwind_l[24]*alpha[26]+0.25*fUpwind_l[0]*alpha[26]+0.2500000000000001*alpha[2]*fUpwind_l[24]+0.2500000000000001*fUpwind_l[2]*alpha[24]+0.223606797749979*alpha[4]*fUpwind_l[9]+0.223606797749979*fUpwind_l[4]*alpha[9]; 
  Ghat_l[27] = 0.159719141249985*alpha[26]*fUpwind_l[30]+0.223606797749979*alpha[18]*fUpwind_l[30]+0.2500000000000001*alpha[2]*fUpwind_l[30]+0.223606797749979*alpha[17]*fUpwind_l[29]+0.2500000000000001*alpha[1]*fUpwind_l[29]+0.159719141249985*alpha[24]*fUpwind_l[27]+0.223606797749979*alpha[16]*fUpwind_l[27]+0.25*alpha[0]*fUpwind_l[27]+0.25*alpha[7]*fUpwind_l[26]+0.25*fUpwind_l[7]*alpha[26]+0.25*alpha[6]*fUpwind_l[25]+0.2500000000000001*alpha[3]*fUpwind_l[24]+0.2500000000000001*fUpwind_l[3]*alpha[24]+0.223606797749979*alpha[9]*fUpwind_l[14]+0.223606797749979*alpha[4]*fUpwind_l[10]; 
  Ghat_l[28] = 0.2500000000000001*alpha[3]*fUpwind_l[31]+0.25*alpha[6]*fUpwind_l[30]+0.25*alpha[7]*fUpwind_l[29]+0.159719141249985*alpha[24]*fUpwind_l[28]+0.25*alpha[0]*fUpwind_l[28]+0.2500000000000001*alpha[1]*fUpwind_l[26]+0.159719141249985*fUpwind_l[25]*alpha[26]+0.2500000000000001*fUpwind_l[1]*alpha[26]+0.2500000000000001*alpha[2]*fUpwind_l[25]+0.25*fUpwind_l[5]*alpha[24]+0.223606797749979*alpha[4]*fUpwind_l[12]+0.223606797749979*fUpwind_l[8]*alpha[9]; 
  Ghat_l[29] = 0.159719141249985*alpha[26]*fUpwind_l[31]+0.223606797749979*alpha[18]*fUpwind_l[31]+0.2500000000000001*alpha[2]*fUpwind_l[31]+0.159719141249985*alpha[24]*fUpwind_l[29]+0.223606797749979*alpha[16]*fUpwind_l[29]+0.25*alpha[0]*fUpwind_l[29]+0.25*alpha[7]*fUpwind_l[28]+0.223606797749979*alpha[17]*fUpwind_l[27]+0.2500000000000001*alpha[1]*fUpwind_l[27]+0.2500000000000001*fUpwind_l[11]*alpha[26]+0.2500000000000001*alpha[3]*fUpwind_l[25]+0.25*alpha[6]*fUpwind_l[24]+0.25*fUpwind_l[6]*alpha[24]+0.223606797749979*alpha[9]*fUpwind_l[15]+0.223606797749979*alpha[4]*fUpwind_l[13]; 
  Ghat_l[30] = 0.223606797749979*alpha[17]*fUpwind_l[31]+0.2500000000000001*alpha[1]*fUpwind_l[31]+0.159719141249985*alpha[24]*fUpwind_l[30]+0.223606797749979*alpha[16]*fUpwind_l[30]+0.25*alpha[0]*fUpwind_l[30]+0.25*alpha[6]*fUpwind_l[28]+0.159719141249985*alpha[26]*fUpwind_l[27]+0.223606797749979*alpha[18]*fUpwind_l[27]+0.2500000000000001*alpha[2]*fUpwind_l[27]+0.2500000000000001*alpha[3]*fUpwind_l[26]+0.2500000000000001*fUpwind_l[3]*alpha[26]+0.25*alpha[7]*fUpwind_l[24]+0.25*fUpwind_l[7]*alpha[24]+0.223606797749979*alpha[4]*fUpwind_l[14]+0.223606797749979*alpha[9]*fUpwind_l[10]; 
  Ghat_l[31] = 0.159719141249985*alpha[24]*fUpwind_l[31]+0.223606797749979*alpha[16]*fUpwind_l[31]+0.25*alpha[0]*fUpwind_l[31]+0.223606797749979*alpha[17]*fUpwind_l[30]+0.2500000000000001*alpha[1]*fUpwind_l[30]+0.159719141249985*alpha[26]*fUpwind_l[29]+0.223606797749979*alpha[18]*fUpwind_l[29]+0.2500000000000001*alpha[2]*fUpwind_l[29]+0.2500000000000001*alpha[3]*fUpwind_l[28]+0.25*alpha[6]*fUpwind_l[26]+0.25*fUpwind_l[6]*alpha[26]+0.25*alpha[7]*fUpwind_l[25]+0.2500000000000001*fUpwind_l[11]*alpha[24]+0.223606797749979*alpha[4]*fUpwind_l[15]+0.223606797749979*alpha[9]*fUpwind_l[13]; 

  Ghat_r[0] = 0.25*alpha[26]*fUpwind_r[26]+0.25*alpha[24]*fUpwind_r[24]+0.25*alpha[18]*fUpwind_r[18]+0.25*alpha[17]*fUpwind_r[17]+0.25*alpha[16]*fUpwind_r[16]+0.25*alpha[9]*fUpwind_r[9]+0.25*alpha[7]*fUpwind_r[7]+0.25*alpha[6]*fUpwind_r[6]+0.25*alpha[4]*fUpwind_r[4]+0.25*alpha[3]*fUpwind_r[3]+0.25*alpha[2]*fUpwind_r[2]+0.25*alpha[1]*fUpwind_r[1]+0.25*alpha[0]*fUpwind_r[0]; 
  Ghat_r[1] = 0.2500000000000001*alpha[26]*fUpwind_r[28]+0.2500000000000001*alpha[24]*fUpwind_r[25]+0.2500000000000001*alpha[18]*fUpwind_r[20]+0.2500000000000001*alpha[16]*fUpwind_r[17]+0.2500000000000001*fUpwind_r[16]*alpha[17]+0.25*alpha[9]*fUpwind_r[12]+0.25*alpha[7]*fUpwind_r[11]+0.25*alpha[4]*fUpwind_r[8]+0.25*alpha[3]*fUpwind_r[6]+0.25*fUpwind_r[3]*alpha[6]+0.25*alpha[2]*fUpwind_r[5]+0.25*alpha[0]*fUpwind_r[1]+0.25*fUpwind_r[0]*alpha[1]; 
  Ghat_r[2] = 0.2500000000000001*alpha[24]*fUpwind_r[26]+0.2500000000000001*fUpwind_r[24]*alpha[26]+0.2500000000000001*alpha[17]*fUpwind_r[20]+0.2500000000000001*alpha[16]*fUpwind_r[18]+0.2500000000000001*fUpwind_r[16]*alpha[18]+0.25*alpha[6]*fUpwind_r[11]+0.25*alpha[4]*fUpwind_r[9]+0.25*fUpwind_r[4]*alpha[9]+0.25*alpha[3]*fUpwind_r[7]+0.25*fUpwind_r[3]*alpha[7]+0.25*alpha[1]*fUpwind_r[5]+0.25*alpha[0]*fUpwind_r[2]+0.25*fUpwind_r[0]*alpha[2]; 
  Ghat_r[3] = 0.2500000000000001*alpha[26]*fUpwind_r[30]+0.2500000000000001*alpha[24]*fUpwind_r[27]+0.223606797749979*alpha[7]*fUpwind_r[18]+0.223606797749979*fUpwind_r[7]*alpha[18]+0.223606797749979*alpha[6]*fUpwind_r[17]+0.223606797749979*fUpwind_r[6]*alpha[17]+0.223606797749979*alpha[3]*fUpwind_r[16]+0.223606797749979*fUpwind_r[3]*alpha[16]+0.25*alpha[9]*fUpwind_r[14]+0.25*alpha[4]*fUpwind_r[10]+0.25*alpha[2]*fUpwind_r[7]+0.25*fUpwind_r[2]*alpha[7]+0.25*alpha[1]*fUpwind_r[6]+0.25*fUpwind_r[1]*alpha[6]+0.25*alpha[0]*fUpwind_r[3]+0.25*fUpwind_r[0]*alpha[3]; 
  Ghat_r[4] = 0.223606797749979*alpha[9]*fUpwind_r[26]+0.223606797749979*fUpwind_r[9]*alpha[26]+0.223606797749979*alpha[4]*fUpwind_r[24]+0.223606797749979*fUpwind_r[4]*alpha[24]+0.2500000000000001*alpha[18]*fUpwind_r[22]+0.2500000000000001*alpha[17]*fUpwind_r[21]+0.2500000000000001*alpha[16]*fUpwind_r[19]+0.25*alpha[7]*fUpwind_r[14]+0.25*alpha[6]*fUpwind_r[13]+0.25*alpha[3]*fUpwind_r[10]+0.25*alpha[2]*fUpwind_r[9]+0.25*fUpwind_r[2]*alpha[9]+0.25*alpha[1]*fUpwind_r[8]+0.25*alpha[0]*fUpwind_r[4]+0.25*fUpwind_r[0]*alpha[4]; 
  Ghat_r[5] = 0.25*alpha[24]*fUpwind_r[28]+0.25*fUpwind_r[25]*alpha[26]+0.25*alpha[16]*fUpwind_r[20]+0.25*alpha[17]*fUpwind_r[18]+0.25*fUpwind_r[17]*alpha[18]+0.25*alpha[4]*fUpwind_r[12]+0.25*alpha[3]*fUpwind_r[11]+0.25*fUpwind_r[8]*alpha[9]+0.25*alpha[6]*fUpwind_r[7]+0.25*fUpwind_r[6]*alpha[7]+0.25*alpha[0]*fUpwind_r[5]+0.25*alpha[1]*fUpwind_r[2]+0.25*fUpwind_r[1]*alpha[2]; 
  Ghat_r[6] = 0.25*alpha[26]*fUpwind_r[31]+0.25*alpha[24]*fUpwind_r[29]+0.223606797749979*alpha[7]*fUpwind_r[20]+0.223606797749979*fUpwind_r[11]*alpha[18]+0.223606797749979*alpha[3]*fUpwind_r[17]+0.223606797749979*fUpwind_r[3]*alpha[17]+0.223606797749979*alpha[6]*fUpwind_r[16]+0.223606797749979*fUpwind_r[6]*alpha[16]+0.25*alpha[9]*fUpwind_r[15]+0.25*alpha[4]*fUpwind_r[13]+0.25*alpha[2]*fUpwind_r[11]+0.25*fUpwind_r[5]*alpha[7]+0.25*alpha[0]*fUpwind_r[6]+0.25*fUpwind_r[0]*alpha[6]+0.25*alpha[1]*fUpwind_r[3]+0.25*fUpwind_r[1]*alpha[3]; 
  Ghat_r[7] = 0.25*alpha[24]*fUpwind_r[30]+0.25*alpha[26]*fUpwind_r[27]+0.223606797749979*alpha[6]*fUpwind_r[20]+0.223606797749979*alpha[3]*fUpwind_r[18]+0.223606797749979*fUpwind_r[3]*alpha[18]+0.223606797749979*fUpwind_r[11]*alpha[17]+0.223606797749979*alpha[7]*fUpwind_r[16]+0.223606797749979*fUpwind_r[7]*alpha[16]+0.25*alpha[4]*fUpwind_r[14]+0.25*alpha[1]*fUpwind_r[11]+0.25*alpha[9]*fUpwind_r[10]+0.25*alpha[0]*fUpwind_r[7]+0.25*fUpwind_r[0]*alpha[7]+0.25*fUpwind_r[5]*alpha[6]+0.25*alpha[2]*fUpwind_r[3]+0.25*fUpwind_r[2]*alpha[3]; 
  Ghat_r[8] = 0.223606797749979*alpha[9]*fUpwind_r[28]+0.223606797749979*fUpwind_r[12]*alpha[26]+0.223606797749979*alpha[4]*fUpwind_r[25]+0.223606797749979*fUpwind_r[8]*alpha[24]+0.25*alpha[18]*fUpwind_r[23]+0.25*alpha[16]*fUpwind_r[21]+0.25*alpha[17]*fUpwind_r[19]+0.25*alpha[7]*fUpwind_r[15]+0.25*alpha[3]*fUpwind_r[13]+0.25*alpha[2]*fUpwind_r[12]+0.25*alpha[6]*fUpwind_r[10]+0.25*fUpwind_r[5]*alpha[9]+0.25*alpha[0]*fUpwind_r[8]+0.25*alpha[1]*fUpwind_r[4]+0.25*fUpwind_r[1]*alpha[4]; 
  Ghat_r[9] = 0.223606797749979*alpha[4]*fUpwind_r[26]+0.223606797749979*fUpwind_r[4]*alpha[26]+0.223606797749979*alpha[9]*fUpwind_r[24]+0.223606797749979*fUpwind_r[9]*alpha[24]+0.25*alpha[17]*fUpwind_r[23]+0.25*alpha[16]*fUpwind_r[22]+0.25*alpha[18]*fUpwind_r[19]+0.25*alpha[6]*fUpwind_r[15]+0.25*alpha[3]*fUpwind_r[14]+0.25*alpha[1]*fUpwind_r[12]+0.25*alpha[7]*fUpwind_r[10]+0.25*alpha[0]*fUpwind_r[9]+0.25*fUpwind_r[0]*alpha[9]+0.25*alpha[2]*fUpwind_r[4]+0.25*fUpwind_r[2]*alpha[4]; 
  Ghat_r[10] = 0.223606797749979*alpha[9]*fUpwind_r[30]+0.223606797749979*alpha[4]*fUpwind_r[27]+0.223606797749979*fUpwind_r[14]*alpha[26]+0.223606797749979*fUpwind_r[10]*alpha[24]+0.223606797749979*alpha[7]*fUpwind_r[22]+0.223606797749979*alpha[6]*fUpwind_r[21]+0.223606797749979*alpha[3]*fUpwind_r[19]+0.223606797749979*fUpwind_r[14]*alpha[18]+0.223606797749979*fUpwind_r[13]*alpha[17]+0.223606797749979*fUpwind_r[10]*alpha[16]+0.25*alpha[2]*fUpwind_r[14]+0.25*alpha[1]*fUpwind_r[13]+0.25*alpha[0]*fUpwind_r[10]+0.25*alpha[7]*fUpwind_r[9]+0.25*fUpwind_r[7]*alpha[9]+0.25*alpha[6]*fUpwind_r[8]+0.25*alpha[3]*fUpwind_r[4]+0.25*fUpwind_r[3]*alpha[4]; 
  Ghat_r[11] = 0.2500000000000001*alpha[24]*fUpwind_r[31]+0.2500000000000001*alpha[26]*fUpwind_r[29]+0.223606797749979*alpha[3]*fUpwind_r[20]+0.223606797749979*alpha[6]*fUpwind_r[18]+0.223606797749979*fUpwind_r[6]*alpha[18]+0.223606797749979*alpha[7]*fUpwind_r[17]+0.223606797749979*fUpwind_r[7]*alpha[17]+0.223606797749979*fUpwind_r[11]*alpha[16]+0.25*alpha[4]*fUpwind_r[15]+0.25*alpha[9]*fUpwind_r[13]+0.25*alpha[0]*fUpwind_r[11]+0.25*alpha[1]*fUpwind_r[7]+0.25*fUpwind_r[1]*alpha[7]+0.25*alpha[2]*fUpwind_r[6]+0.25*fUpwind_r[2]*alpha[6]+0.25*alpha[3]*fUpwind_r[5]; 
  Ghat_r[12] = 0.223606797749979*alpha[4]*fUpwind_r[28]+0.223606797749979*fUpwind_r[8]*alpha[26]+0.223606797749979*alpha[9]*fUpwind_r[25]+0.223606797749979*fUpwind_r[12]*alpha[24]+0.2500000000000001*alpha[16]*fUpwind_r[23]+0.2500000000000001*alpha[17]*fUpwind_r[22]+0.2500000000000001*alpha[18]*fUpwind_r[21]+0.25*alpha[3]*fUpwind_r[15]+0.25*alpha[6]*fUpwind_r[14]+0.25*alpha[7]*fUpwind_r[13]+0.25*alpha[0]*fUpwind_r[12]+0.25*alpha[1]*fUpwind_r[9]+0.25*fUpwind_r[1]*alpha[9]+0.25*alpha[2]*fUpwind_r[8]+0.25*alpha[4]*fUpwind_r[5]; 
  Ghat_r[13] = 0.223606797749979*alpha[9]*fUpwind_r[31]+0.223606797749979*alpha[4]*fUpwind_r[29]+0.223606797749979*fUpwind_r[15]*alpha[26]+0.223606797749979*fUpwind_r[13]*alpha[24]+0.223606797749979*alpha[7]*fUpwind_r[23]+0.223606797749979*alpha[3]*fUpwind_r[21]+0.223606797749979*alpha[6]*fUpwind_r[19]+0.223606797749979*fUpwind_r[15]*alpha[18]+0.223606797749979*fUpwind_r[10]*alpha[17]+0.223606797749979*fUpwind_r[13]*alpha[16]+0.25*alpha[2]*fUpwind_r[15]+0.25*alpha[0]*fUpwind_r[13]+0.25*alpha[7]*fUpwind_r[12]+0.25*alpha[9]*fUpwind_r[11]+0.25*alpha[1]*fUpwind_r[10]+0.25*alpha[3]*fUpwind_r[8]+0.25*alpha[4]*fUpwind_r[6]+0.25*fUpwind_r[4]*alpha[6]; 
  Ghat_r[14] = 0.223606797749979*alpha[4]*fUpwind_r[30]+0.223606797749979*alpha[9]*fUpwind_r[27]+0.223606797749979*fUpwind_r[10]*alpha[26]+0.223606797749979*fUpwind_r[14]*alpha[24]+0.223606797749979*alpha[6]*fUpwind_r[23]+0.223606797749979*alpha[3]*fUpwind_r[22]+0.223606797749979*alpha[7]*fUpwind_r[19]+0.223606797749979*fUpwind_r[10]*alpha[18]+0.223606797749979*fUpwind_r[15]*alpha[17]+0.223606797749979*fUpwind_r[14]*alpha[16]+0.25*alpha[1]*fUpwind_r[15]+0.25*alpha[0]*fUpwind_r[14]+0.25*alpha[6]*fUpwind_r[12]+0.25*alpha[2]*fUpwind_r[10]+0.25*alpha[3]*fUpwind_r[9]+0.25*fUpwind_r[3]*alpha[9]+0.25*alpha[4]*fUpwind_r[7]+0.25*fUpwind_r[4]*alpha[7]; 
  Ghat_r[15] = 0.223606797749979*alpha[4]*fUpwind_r[31]+0.223606797749979*alpha[9]*fUpwind_r[29]+0.223606797749979*fUpwind_r[13]*alpha[26]+0.223606797749979*fUpwind_r[15]*alpha[24]+0.223606797749979*alpha[3]*fUpwind_r[23]+0.223606797749979*alpha[6]*fUpwind_r[22]+0.223606797749979*alpha[7]*fUpwind_r[21]+0.223606797749979*fUpwind_r[13]*alpha[18]+0.223606797749979*fUpwind_r[14]*alpha[17]+0.223606797749979*fUpwind_r[15]*alpha[16]+0.25*alpha[0]*fUpwind_r[15]+0.25*alpha[1]*fUpwind_r[14]+0.25*alpha[2]*fUpwind_r[13]+0.25*alpha[3]*fUpwind_r[12]+0.25*alpha[4]*fUpwind_r[11]+0.25*alpha[6]*fUpwind_r[9]+0.25*fUpwind_r[6]*alpha[9]+0.25*alpha[7]*fUpwind_r[8]; 
  Ghat_r[16] = 0.25*alpha[9]*fUpwind_r[22]+0.2500000000000001*alpha[4]*fUpwind_r[19]+0.159719141249985*alpha[18]*fUpwind_r[18]+0.2500000000000001*alpha[2]*fUpwind_r[18]+0.2500000000000001*fUpwind_r[2]*alpha[18]+0.159719141249985*alpha[17]*fUpwind_r[17]+0.2500000000000001*alpha[1]*fUpwind_r[17]+0.2500000000000001*fUpwind_r[1]*alpha[17]+0.159719141249985*alpha[16]*fUpwind_r[16]+0.25*alpha[0]*fUpwind_r[16]+0.25*fUpwind_r[0]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_r[7]+0.223606797749979*alpha[6]*fUpwind_r[6]+0.223606797749979*alpha[3]*fUpwind_r[3]; 
  Ghat_r[17] = 0.25*alpha[9]*fUpwind_r[23]+0.2500000000000001*alpha[4]*fUpwind_r[21]+0.159719141249985*alpha[18]*fUpwind_r[20]+0.2500000000000001*alpha[2]*fUpwind_r[20]+0.25*fUpwind_r[5]*alpha[18]+0.159719141249985*alpha[16]*fUpwind_r[17]+0.25*alpha[0]*fUpwind_r[17]+0.159719141249985*fUpwind_r[16]*alpha[17]+0.25*fUpwind_r[0]*alpha[17]+0.2500000000000001*alpha[1]*fUpwind_r[16]+0.2500000000000001*fUpwind_r[1]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_r[11]+0.223606797749979*alpha[3]*fUpwind_r[6]+0.223606797749979*fUpwind_r[3]*alpha[6]; 
  Ghat_r[18] = 0.2500000000000001*alpha[4]*fUpwind_r[22]+0.159719141249985*alpha[17]*fUpwind_r[20]+0.2500000000000001*alpha[1]*fUpwind_r[20]+0.25*alpha[9]*fUpwind_r[19]+0.159719141249985*alpha[16]*fUpwind_r[18]+0.25*alpha[0]*fUpwind_r[18]+0.159719141249985*fUpwind_r[16]*alpha[18]+0.25*fUpwind_r[0]*alpha[18]+0.25*fUpwind_r[5]*alpha[17]+0.2500000000000001*alpha[2]*fUpwind_r[16]+0.2500000000000001*fUpwind_r[2]*alpha[16]+0.223606797749979*alpha[6]*fUpwind_r[11]+0.223606797749979*alpha[3]*fUpwind_r[7]+0.223606797749979*fUpwind_r[3]*alpha[7]; 
  Ghat_r[19] = 0.223606797749979*fUpwind_r[22]*alpha[26]+0.223606797749979*fUpwind_r[19]*alpha[24]+0.159719141249985*alpha[18]*fUpwind_r[22]+0.2500000000000001*alpha[2]*fUpwind_r[22]+0.159719141249985*alpha[17]*fUpwind_r[21]+0.2500000000000001*alpha[1]*fUpwind_r[21]+0.159719141249985*alpha[16]*fUpwind_r[19]+0.25*alpha[0]*fUpwind_r[19]+0.25*alpha[9]*fUpwind_r[18]+0.25*fUpwind_r[9]*alpha[18]+0.25*fUpwind_r[8]*alpha[17]+0.2500000000000001*alpha[4]*fUpwind_r[16]+0.2500000000000001*fUpwind_r[4]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_r[14]+0.223606797749979*alpha[6]*fUpwind_r[13]+0.223606797749979*alpha[3]*fUpwind_r[10]; 
  Ghat_r[20] = 0.2500000000000001*alpha[4]*fUpwind_r[23]+0.25*alpha[9]*fUpwind_r[21]+0.159719141249985*alpha[16]*fUpwind_r[20]+0.25*alpha[0]*fUpwind_r[20]+0.159719141249985*alpha[17]*fUpwind_r[18]+0.2500000000000001*alpha[1]*fUpwind_r[18]+0.159719141249985*fUpwind_r[17]*alpha[18]+0.2500000000000001*fUpwind_r[1]*alpha[18]+0.2500000000000001*alpha[2]*fUpwind_r[17]+0.2500000000000001*fUpwind_r[2]*alpha[17]+0.25*fUpwind_r[5]*alpha[16]+0.223606797749979*alpha[3]*fUpwind_r[11]+0.223606797749979*alpha[6]*fUpwind_r[7]+0.223606797749979*fUpwind_r[6]*alpha[7]; 
  Ghat_r[21] = 0.223606797749979*fUpwind_r[23]*alpha[26]+0.223606797749979*fUpwind_r[21]*alpha[24]+0.159719141249985*alpha[18]*fUpwind_r[23]+0.2500000000000001*alpha[2]*fUpwind_r[23]+0.159719141249985*alpha[16]*fUpwind_r[21]+0.25*alpha[0]*fUpwind_r[21]+0.25*alpha[9]*fUpwind_r[20]+0.159719141249985*alpha[17]*fUpwind_r[19]+0.2500000000000001*alpha[1]*fUpwind_r[19]+0.2500000000000001*fUpwind_r[12]*alpha[18]+0.2500000000000001*alpha[4]*fUpwind_r[17]+0.2500000000000001*fUpwind_r[4]*alpha[17]+0.25*fUpwind_r[8]*alpha[16]+0.223606797749979*alpha[7]*fUpwind_r[15]+0.223606797749979*alpha[3]*fUpwind_r[13]+0.223606797749979*alpha[6]*fUpwind_r[10]; 
  Ghat_r[22] = 0.223606797749979*fUpwind_r[19]*alpha[26]+0.223606797749979*fUpwind_r[22]*alpha[24]+0.159719141249985*alpha[17]*fUpwind_r[23]+0.2500000000000001*alpha[1]*fUpwind_r[23]+0.159719141249985*alpha[16]*fUpwind_r[22]+0.25*alpha[0]*fUpwind_r[22]+0.159719141249985*alpha[18]*fUpwind_r[19]+0.2500000000000001*alpha[2]*fUpwind_r[19]+0.2500000000000001*alpha[4]*fUpwind_r[18]+0.2500000000000001*fUpwind_r[4]*alpha[18]+0.2500000000000001*fUpwind_r[12]*alpha[17]+0.25*alpha[9]*fUpwind_r[16]+0.25*fUpwind_r[9]*alpha[16]+0.223606797749979*alpha[6]*fUpwind_r[15]+0.223606797749979*alpha[3]*fUpwind_r[14]+0.223606797749979*alpha[7]*fUpwind_r[10]; 
  Ghat_r[23] = 0.223606797749979*fUpwind_r[21]*alpha[26]+0.223606797749979*fUpwind_r[23]*alpha[24]+0.159719141249985*alpha[16]*fUpwind_r[23]+0.25*alpha[0]*fUpwind_r[23]+0.159719141249985*alpha[17]*fUpwind_r[22]+0.2500000000000001*alpha[1]*fUpwind_r[22]+0.159719141249985*alpha[18]*fUpwind_r[21]+0.2500000000000001*alpha[2]*fUpwind_r[21]+0.2500000000000001*alpha[4]*fUpwind_r[20]+0.25*fUpwind_r[8]*alpha[18]+0.25*alpha[9]*fUpwind_r[17]+0.25*fUpwind_r[9]*alpha[17]+0.2500000000000001*fUpwind_r[12]*alpha[16]+0.223606797749979*alpha[3]*fUpwind_r[15]+0.223606797749979*alpha[6]*fUpwind_r[14]+0.223606797749979*alpha[7]*fUpwind_r[13]; 
  Ghat_r[24] = 0.25*alpha[7]*fUpwind_r[30]+0.25*alpha[6]*fUpwind_r[29]+0.2500000000000001*alpha[3]*fUpwind_r[27]+0.159719141249985*alpha[26]*fUpwind_r[26]+0.2500000000000001*alpha[2]*fUpwind_r[26]+0.2500000000000001*fUpwind_r[2]*alpha[26]+0.2500000000000001*alpha[1]*fUpwind_r[25]+0.159719141249985*alpha[24]*fUpwind_r[24]+0.25*alpha[0]*fUpwind_r[24]+0.25*fUpwind_r[0]*alpha[24]+0.223606797749979*alpha[9]*fUpwind_r[9]+0.223606797749979*alpha[4]*fUpwind_r[4]; 
  Ghat_r[25] = 0.25*alpha[7]*fUpwind_r[31]+0.2500000000000001*alpha[3]*fUpwind_r[29]+0.159719141249985*alpha[26]*fUpwind_r[28]+0.2500000000000001*alpha[2]*fUpwind_r[28]+0.25*alpha[6]*fUpwind_r[27]+0.25*fUpwind_r[5]*alpha[26]+0.159719141249985*alpha[24]*fUpwind_r[25]+0.25*alpha[0]*fUpwind_r[25]+0.2500000000000001*alpha[1]*fUpwind_r[24]+0.2500000000000001*fUpwind_r[1]*alpha[24]+0.223606797749979*alpha[9]*fUpwind_r[12]+0.223606797749979*alpha[4]*fUpwind_r[8]; 
  Ghat_r[26] = 0.25*alpha[6]*fUpwind_r[31]+0.2500000000000001*alpha[3]*fUpwind_r[30]+0.2500000000000001*alpha[1]*fUpwind_r[28]+0.25*alpha[7]*fUpwind_r[27]+0.159719141249985*alpha[24]*fUpwind_r[26]+0.25*alpha[0]*fUpwind_r[26]+0.159719141249985*fUpwind_r[24]*alpha[26]+0.25*fUpwind_r[0]*alpha[26]+0.2500000000000001*alpha[2]*fUpwind_r[24]+0.2500000000000001*fUpwind_r[2]*alpha[24]+0.223606797749979*alpha[4]*fUpwind_r[9]+0.223606797749979*fUpwind_r[4]*alpha[9]; 
  Ghat_r[27] = 0.159719141249985*alpha[26]*fUpwind_r[30]+0.223606797749979*alpha[18]*fUpwind_r[30]+0.2500000000000001*alpha[2]*fUpwind_r[30]+0.223606797749979*alpha[17]*fUpwind_r[29]+0.2500000000000001*alpha[1]*fUpwind_r[29]+0.159719141249985*alpha[24]*fUpwind_r[27]+0.223606797749979*alpha[16]*fUpwind_r[27]+0.25*alpha[0]*fUpwind_r[27]+0.25*alpha[7]*fUpwind_r[26]+0.25*fUpwind_r[7]*alpha[26]+0.25*alpha[6]*fUpwind_r[25]+0.2500000000000001*alpha[3]*fUpwind_r[24]+0.2500000000000001*fUpwind_r[3]*alpha[24]+0.223606797749979*alpha[9]*fUpwind_r[14]+0.223606797749979*alpha[4]*fUpwind_r[10]; 
  Ghat_r[28] = 0.2500000000000001*alpha[3]*fUpwind_r[31]+0.25*alpha[6]*fUpwind_r[30]+0.25*alpha[7]*fUpwind_r[29]+0.159719141249985*alpha[24]*fUpwind_r[28]+0.25*alpha[0]*fUpwind_r[28]+0.2500000000000001*alpha[1]*fUpwind_r[26]+0.159719141249985*fUpwind_r[25]*alpha[26]+0.2500000000000001*fUpwind_r[1]*alpha[26]+0.2500000000000001*alpha[2]*fUpwind_r[25]+0.25*fUpwind_r[5]*alpha[24]+0.223606797749979*alpha[4]*fUpwind_r[12]+0.223606797749979*fUpwind_r[8]*alpha[9]; 
  Ghat_r[29] = 0.159719141249985*alpha[26]*fUpwind_r[31]+0.223606797749979*alpha[18]*fUpwind_r[31]+0.2500000000000001*alpha[2]*fUpwind_r[31]+0.159719141249985*alpha[24]*fUpwind_r[29]+0.223606797749979*alpha[16]*fUpwind_r[29]+0.25*alpha[0]*fUpwind_r[29]+0.25*alpha[7]*fUpwind_r[28]+0.223606797749979*alpha[17]*fUpwind_r[27]+0.2500000000000001*alpha[1]*fUpwind_r[27]+0.2500000000000001*fUpwind_r[11]*alpha[26]+0.2500000000000001*alpha[3]*fUpwind_r[25]+0.25*alpha[6]*fUpwind_r[24]+0.25*fUpwind_r[6]*alpha[24]+0.223606797749979*alpha[9]*fUpwind_r[15]+0.223606797749979*alpha[4]*fUpwind_r[13]; 
  Ghat_r[30] = 0.223606797749979*alpha[17]*fUpwind_r[31]+0.2500000000000001*alpha[1]*fUpwind_r[31]+0.159719141249985*alpha[24]*fUpwind_r[30]+0.223606797749979*alpha[16]*fUpwind_r[30]+0.25*alpha[0]*fUpwind_r[30]+0.25*alpha[6]*fUpwind_r[28]+0.159719141249985*alpha[26]*fUpwind_r[27]+0.223606797749979*alpha[18]*fUpwind_r[27]+0.2500000000000001*alpha[2]*fUpwind_r[27]+0.2500000000000001*alpha[3]*fUpwind_r[26]+0.2500000000000001*fUpwind_r[3]*alpha[26]+0.25*alpha[7]*fUpwind_r[24]+0.25*fUpwind_r[7]*alpha[24]+0.223606797749979*alpha[4]*fUpwind_r[14]+0.223606797749979*alpha[9]*fUpwind_r[10]; 
  Ghat_r[31] = 0.159719141249985*alpha[24]*fUpwind_r[31]+0.223606797749979*alpha[16]*fUpwind_r[31]+0.25*alpha[0]*fUpwind_r[31]+0.223606797749979*alpha[17]*fUpwind_r[30]+0.2500000000000001*alpha[1]*fUpwind_r[30]+0.159719141249985*alpha[26]*fUpwind_r[29]+0.223606797749979*alpha[18]*fUpwind_r[29]+0.2500000000000001*alpha[2]*fUpwind_r[29]+0.2500000000000001*alpha[3]*fUpwind_r[28]+0.25*alpha[6]*fUpwind_r[26]+0.25*fUpwind_r[6]*alpha[26]+0.25*alpha[7]*fUpwind_r[25]+0.2500000000000001*fUpwind_r[11]*alpha[24]+0.223606797749979*alpha[4]*fUpwind_r[15]+0.223606797749979*alpha[9]*fUpwind_r[13]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*rdv2; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*rdv2; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*rdv2; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*rdv2; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*rdv2; 
  out[7] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[8] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdv2; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*rdv2; 
  out[10] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*rdv2; 
  out[11] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*rdv2; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*rdv2; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*rdv2; 
  out[14] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*rdv2; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*rdv2; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*rdv2; 
  out[17] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*rdv2; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*rdv2; 
  out[19] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*rdv2; 
  out[20] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*rdv2; 
  out[21] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*rdv2; 
  out[22] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*rdv2; 
  out[23] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*rdv2; 
  out[24] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*rdv2; 
  out[25] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*rdv2; 
  out[26] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*rdv2; 
  out[27] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*rdv2; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*rdv2; 
  out[29] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*rdv2; 
  out[30] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*rdv2; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*rdv2; 
  out[32] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*rdv2; 
  out[33] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*rdv2; 
  out[34] += (1.58113883008419*Ghat_l[2]-1.58113883008419*Ghat_r[2])*rdv2; 
  out[35] += (1.58113883008419*Ghat_l[3]-1.58113883008419*Ghat_r[3])*rdv2; 
  out[36] += (1.58113883008419*Ghat_l[4]-1.58113883008419*Ghat_r[4])*rdv2; 
  out[37] += (1.58113883008419*Ghat_l[5]-1.58113883008419*Ghat_r[5])*rdv2; 
  out[38] += (1.58113883008419*Ghat_l[6]-1.58113883008419*Ghat_r[6])*rdv2; 
  out[39] += (1.58113883008419*Ghat_l[7]-1.58113883008419*Ghat_r[7])*rdv2; 
  out[40] += (1.58113883008419*Ghat_l[8]-1.58113883008419*Ghat_r[8])*rdv2; 
  out[41] += (1.58113883008419*Ghat_l[9]-1.58113883008419*Ghat_r[9])*rdv2; 
  out[42] += (1.58113883008419*Ghat_l[10]-1.58113883008419*Ghat_r[10])*rdv2; 
  out[43] += (1.58113883008419*Ghat_l[11]-1.58113883008419*Ghat_r[11])*rdv2; 
  out[44] += (1.58113883008419*Ghat_l[12]-1.58113883008419*Ghat_r[12])*rdv2; 
  out[45] += (1.58113883008419*Ghat_l[13]-1.58113883008419*Ghat_r[13])*rdv2; 
  out[46] += (1.58113883008419*Ghat_l[14]-1.58113883008419*Ghat_r[14])*rdv2; 
  out[47] += (1.58113883008419*Ghat_l[15]-1.58113883008419*Ghat_r[15])*rdv2; 
  out[48] += (0.7071067811865475*Ghat_l[16]-0.7071067811865475*Ghat_r[16])*rdv2; 
  out[49] += (0.7071067811865475*Ghat_l[17]-0.7071067811865475*Ghat_r[17])*rdv2; 
  out[50] += (0.7071067811865475*Ghat_l[18]-0.7071067811865475*Ghat_r[18])*rdv2; 
  out[51] += -1.224744871391589*(Ghat_r[16]+Ghat_l[16])*rdv2; 
  out[52] += (0.7071067811865475*Ghat_l[19]-0.7071067811865475*Ghat_r[19])*rdv2; 
  out[53] += (0.7071067811865475*Ghat_l[20]-0.7071067811865475*Ghat_r[20])*rdv2; 
  out[54] += -1.224744871391589*(Ghat_r[17]+Ghat_l[17])*rdv2; 
  out[55] += -1.224744871391589*(Ghat_r[18]+Ghat_l[18])*rdv2; 
  out[56] += (0.7071067811865475*Ghat_l[21]-0.7071067811865475*Ghat_r[21])*rdv2; 
  out[57] += (0.7071067811865475*Ghat_l[22]-0.7071067811865475*Ghat_r[22])*rdv2; 
  out[58] += -1.224744871391589*(Ghat_r[19]+Ghat_l[19])*rdv2; 
  out[59] += -1.224744871391589*(Ghat_r[20]+Ghat_l[20])*rdv2; 
  out[60] += (0.7071067811865475*Ghat_l[23]-0.7071067811865475*Ghat_r[23])*rdv2; 
  out[61] += -1.224744871391589*(Ghat_r[21]+Ghat_l[21])*rdv2; 
  out[62] += -1.224744871391589*(Ghat_r[22]+Ghat_l[22])*rdv2; 
  out[63] += -1.224744871391589*(Ghat_r[23]+Ghat_l[23])*rdv2; 
  out[64] += (0.7071067811865475*Ghat_l[24]-0.7071067811865475*Ghat_r[24])*rdv2; 
  out[65] += (0.7071067811865475*Ghat_l[25]-0.7071067811865475*Ghat_r[25])*rdv2; 
  out[66] += (0.7071067811865475*Ghat_l[26]-0.7071067811865475*Ghat_r[26])*rdv2; 
  out[67] += -1.224744871391589*(Ghat_r[24]+Ghat_l[24])*rdv2; 
  out[68] += (0.7071067811865475*Ghat_l[27]-0.7071067811865475*Ghat_r[27])*rdv2; 
  out[69] += (0.7071067811865475*Ghat_l[28]-0.7071067811865475*Ghat_r[28])*rdv2; 
  out[70] += -1.224744871391589*(Ghat_r[25]+Ghat_l[25])*rdv2; 
  out[71] += -1.224744871391589*(Ghat_r[26]+Ghat_l[26])*rdv2; 
  out[72] += (0.7071067811865475*Ghat_l[29]-0.7071067811865475*Ghat_r[29])*rdv2; 
  out[73] += (0.7071067811865475*Ghat_l[30]-0.7071067811865475*Ghat_r[30])*rdv2; 
  out[74] += -1.224744871391589*(Ghat_r[27]+Ghat_l[27])*rdv2; 
  out[75] += -1.224744871391589*(Ghat_r[28]+Ghat_l[28])*rdv2; 
  out[76] += (0.7071067811865475*Ghat_l[31]-0.7071067811865475*Ghat_r[31])*rdv2; 
  out[77] += -1.224744871391589*(Ghat_r[29]+Ghat_l[29])*rdv2; 
  out[78] += -1.224744871391589*(Ghat_r[30]+Ghat_l[30])*rdv2; 
  out[79] += -1.224744871391589*(Ghat_r[31]+Ghat_l[31])*rdv2; 

  return 0.;

} 
