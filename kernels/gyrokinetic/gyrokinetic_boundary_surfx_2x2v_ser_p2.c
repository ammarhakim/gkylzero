#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_ser_4x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double *alpha_edge, const double *alpha_skin, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // alpha_edge: Surface expansion of phase space flux on the lower edges of the edge cell.
  // alpha_skin: Surface expansion of phase space flux on the lower edges of the skin cell.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];
  double wmu = w[3];
  double rdmu2 = 2.0/dxv[3];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wvparSq = w[2]*w[2];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[3]*w[3];
  double rdmu2Sq = rdmu2*rdmu2;

  double cflFreq = 0.0;

  if (edge == -1) { 

  const double *alphaR = &alpha_edge[0];
  double fUpOrdR[27] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = ser_4x_p2_surfx1_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = ser_4x_p2_surfx1_eval_quad_node_0_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = ser_4x_p2_surfx1_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = ser_4x_p2_surfx1_eval_quad_node_1_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = ser_4x_p2_surfx1_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = ser_4x_p2_surfx1_eval_quad_node_2_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.3952847075210471*alphaR[8])-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = ser_4x_p2_surfx1_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = ser_4x_p2_surfx1_eval_quad_node_3_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.3952847075210471*alphaR[8])-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = ser_4x_p2_surfx1_eval_quad_node_4_r(fskin); 
  } else { 
    fUpOrdR[4] = ser_4x_p2_surfx1_eval_quad_node_4_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.3952847075210471*alphaR[8])-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = ser_4x_p2_surfx1_eval_quad_node_5_r(fskin); 
  } else { 
    fUpOrdR[5] = ser_4x_p2_surfx1_eval_quad_node_5_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = ser_4x_p2_surfx1_eval_quad_node_6_r(fskin); 
  } else { 
    fUpOrdR[6] = ser_4x_p2_surfx1_eval_quad_node_6_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = ser_4x_p2_surfx1_eval_quad_node_7_r(fskin); 
  } else { 
    fUpOrdR[7] = ser_4x_p2_surfx1_eval_quad_node_7_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = ser_4x_p2_surfx1_eval_quad_node_8_r(fskin); 
  } else { 
    fUpOrdR[8] = ser_4x_p2_surfx1_eval_quad_node_8_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[9] = ser_4x_p2_surfx1_eval_quad_node_9_r(fskin); 
  } else { 
    fUpOrdR[9] = ser_4x_p2_surfx1_eval_quad_node_9_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[10] = ser_4x_p2_surfx1_eval_quad_node_10_r(fskin); 
  } else { 
    fUpOrdR[10] = ser_4x_p2_surfx1_eval_quad_node_10_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[11] = ser_4x_p2_surfx1_eval_quad_node_11_r(fskin); 
  } else { 
    fUpOrdR[11] = ser_4x_p2_surfx1_eval_quad_node_11_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0]-0.3952847075210471*alphaR[8];
  if (alphaR_n > 0.) {
    fUpOrdR[12] = ser_4x_p2_surfx1_eval_quad_node_12_r(fskin); 
  } else { 
    fUpOrdR[12] = ser_4x_p2_surfx1_eval_quad_node_12_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0]-0.3952847075210471*alphaR[8];
  if (alphaR_n > 0.) {
    fUpOrdR[13] = ser_4x_p2_surfx1_eval_quad_node_13_r(fskin); 
  } else { 
    fUpOrdR[13] = ser_4x_p2_surfx1_eval_quad_node_13_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0]-0.3952847075210471*alphaR[8];
  if (alphaR_n > 0.) {
    fUpOrdR[14] = ser_4x_p2_surfx1_eval_quad_node_14_r(fskin); 
  } else { 
    fUpOrdR[14] = ser_4x_p2_surfx1_eval_quad_node_14_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[15] = ser_4x_p2_surfx1_eval_quad_node_15_r(fskin); 
  } else { 
    fUpOrdR[15] = ser_4x_p2_surfx1_eval_quad_node_15_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[16] = ser_4x_p2_surfx1_eval_quad_node_16_r(fskin); 
  } else { 
    fUpOrdR[16] = ser_4x_p2_surfx1_eval_quad_node_16_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[17] = ser_4x_p2_surfx1_eval_quad_node_17_r(fskin); 
  } else { 
    fUpOrdR[17] = ser_4x_p2_surfx1_eval_quad_node_17_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[18] = ser_4x_p2_surfx1_eval_quad_node_18_r(fskin); 
  } else { 
    fUpOrdR[18] = ser_4x_p2_surfx1_eval_quad_node_18_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[19] = ser_4x_p2_surfx1_eval_quad_node_19_r(fskin); 
  } else { 
    fUpOrdR[19] = ser_4x_p2_surfx1_eval_quad_node_19_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[20] = ser_4x_p2_surfx1_eval_quad_node_20_r(fskin); 
  } else { 
    fUpOrdR[20] = ser_4x_p2_surfx1_eval_quad_node_20_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.3952847075210471*alphaR[8])+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[21] = ser_4x_p2_surfx1_eval_quad_node_21_r(fskin); 
  } else { 
    fUpOrdR[21] = ser_4x_p2_surfx1_eval_quad_node_21_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.3952847075210471*alphaR[8])+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[22] = ser_4x_p2_surfx1_eval_quad_node_22_r(fskin); 
  } else { 
    fUpOrdR[22] = ser_4x_p2_surfx1_eval_quad_node_22_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.3952847075210471*alphaR[8])+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[23] = ser_4x_p2_surfx1_eval_quad_node_23_r(fskin); 
  } else { 
    fUpOrdR[23] = ser_4x_p2_surfx1_eval_quad_node_23_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[24] = ser_4x_p2_surfx1_eval_quad_node_24_r(fskin); 
  } else { 
    fUpOrdR[24] = ser_4x_p2_surfx1_eval_quad_node_24_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[25] = ser_4x_p2_surfx1_eval_quad_node_25_r(fskin); 
  } else { 
    fUpOrdR[25] = ser_4x_p2_surfx1_eval_quad_node_25_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3162277660168378*alphaR[8]+0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[26] = ser_4x_p2_surfx1_eval_quad_node_26_r(fskin); 
  } else { 
    fUpOrdR[26] = ser_4x_p2_surfx1_eval_quad_node_26_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[20] = {0.};
  ser_4x_p2_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[20] = {0.}; 
  GhatR[0] = 0.3535533905932737*(alphaR[8]*fUpR[8]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.3535533905932737*alphaR[8]*fUpR[12]+0.3162277660168379*alphaR[1]*fUpR[7]+0.3535533905932737*(alphaR[2]*fUpR[4]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.3162277660168379*(alphaR[2]*fUpR[8]+fUpR[2]*alphaR[8])+0.3535533905932737*(alphaR[1]*fUpR[4]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.3535533905932737*(alphaR[8]*fUpR[14]+alphaR[2]*fUpR[6]+alphaR[1]*fUpR[5]+alphaR[0]*fUpR[3]); 
  GhatR[4] = 0.3162277660168379*(alphaR[2]*fUpR[12]+alphaR[1]*fUpR[11])+0.3162277660168379*fUpR[4]*alphaR[8]+0.3535533905932737*(alphaR[0]*fUpR[4]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[5] = 0.3535533905932737*alphaR[8]*fUpR[18]+0.3162277660168379*alphaR[1]*fUpR[13]+0.3535533905932737*(alphaR[2]*fUpR[10]+alphaR[0]*fUpR[5]+alphaR[1]*fUpR[3]); 
  GhatR[6] = 0.3162277660168379*alphaR[2]*fUpR[14]+0.3535533905932737*alphaR[1]*fUpR[10]+0.3162277660168379*fUpR[6]*alphaR[8]+0.3535533905932737*(alphaR[0]*fUpR[6]+alphaR[2]*fUpR[3]); 
  GhatR[7] = 0.3535533905932737*(alphaR[2]*fUpR[11]+alphaR[0]*fUpR[7])+0.3162277660168379*alphaR[1]*fUpR[1]; 
  GhatR[8] = 0.3535533905932737*alphaR[1]*fUpR[12]+0.2258769757263128*alphaR[8]*fUpR[8]+0.3535533905932737*(alphaR[0]*fUpR[8]+fUpR[0]*alphaR[8])+0.3162277660168379*alphaR[2]*fUpR[2]; 
  GhatR[9] = 0.3535533905932737*(alphaR[2]*fUpR[16]+alphaR[1]*fUpR[15]+alphaR[0]*fUpR[9]); 
  GhatR[10] = 0.3162277660168379*(alphaR[2]*fUpR[18]+alphaR[1]*fUpR[17]+alphaR[8]*fUpR[10])+0.3535533905932737*(alphaR[0]*fUpR[10]+alphaR[1]*fUpR[6]+alphaR[2]*fUpR[5]); 
  GhatR[11] = 0.3162277660168379*alphaR[8]*fUpR[11]+0.3535533905932737*(alphaR[0]*fUpR[11]+alphaR[2]*fUpR[7])+0.3162277660168379*alphaR[1]*fUpR[4]; 
  GhatR[12] = 0.2258769757263128*alphaR[8]*fUpR[12]+0.3535533905932737*(alphaR[0]*fUpR[12]+alphaR[1]*fUpR[8]+fUpR[1]*alphaR[8])+0.3162277660168379*alphaR[2]*fUpR[4]; 
  GhatR[13] = 0.3535533905932737*(alphaR[2]*fUpR[17]+alphaR[0]*fUpR[13])+0.3162277660168379*alphaR[1]*fUpR[5]; 
  GhatR[14] = 0.3535533905932737*alphaR[1]*fUpR[18]+0.2258769757263128*alphaR[8]*fUpR[14]+0.3535533905932737*(alphaR[0]*fUpR[14]+fUpR[3]*alphaR[8])+0.3162277660168379*alphaR[2]*fUpR[6]; 
  GhatR[15] = 0.3535533905932737*(alphaR[2]*fUpR[19]+alphaR[0]*fUpR[15]+alphaR[1]*fUpR[9]); 
  GhatR[16] = 0.3535533905932737*alphaR[1]*fUpR[19]+0.3162277660168379*alphaR[8]*fUpR[16]+0.3535533905932737*(alphaR[0]*fUpR[16]+alphaR[2]*fUpR[9]); 
  GhatR[17] = 0.3162277660168379*alphaR[8]*fUpR[17]+0.3535533905932737*(alphaR[0]*fUpR[17]+alphaR[2]*fUpR[13])+0.3162277660168379*alphaR[1]*fUpR[10]; 
  GhatR[18] = 0.2258769757263128*alphaR[8]*fUpR[18]+0.3535533905932737*(alphaR[0]*fUpR[18]+alphaR[1]*fUpR[14])+0.3162277660168379*alphaR[2]*fUpR[10]+0.3535533905932737*fUpR[5]*alphaR[8]; 
  GhatR[19] = 0.3162277660168379*alphaR[8]*fUpR[19]+0.3535533905932737*(alphaR[0]*fUpR[19]+alphaR[1]*fUpR[16]+alphaR[2]*fUpR[15]); 

  out[0] += -0.7071067811865475*GhatR[0]*rdx2; 
  out[1] += -1.224744871391589*GhatR[0]*rdx2; 
  out[2] += -0.7071067811865475*GhatR[1]*rdx2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdx2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdx2; 
  out[5] += -1.224744871391589*GhatR[1]*rdx2; 
  out[6] += -1.224744871391589*GhatR[2]*rdx2; 
  out[7] += -0.7071067811865475*GhatR[4]*rdx2; 
  out[8] += -1.224744871391589*GhatR[3]*rdx2; 
  out[9] += -0.7071067811865475*GhatR[5]*rdx2; 
  out[10] += -0.7071067811865475*GhatR[6]*rdx2; 
  out[11] += -1.58113883008419*GhatR[0]*rdx2; 
  out[12] += -0.7071067811865475*GhatR[7]*rdx2; 
  out[13] += -0.7071067811865475*GhatR[8]*rdx2; 
  out[14] += -0.7071067811865475*GhatR[9]*rdx2; 
  out[15] += -1.224744871391589*GhatR[4]*rdx2; 
  out[16] += -1.224744871391589*GhatR[5]*rdx2; 
  out[17] += -1.224744871391589*GhatR[6]*rdx2; 
  out[18] += -0.7071067811865475*GhatR[10]*rdx2; 
  out[19] += -1.58113883008419*GhatR[1]*rdx2; 
  out[20] += -1.224744871391589*GhatR[7]*rdx2; 
  out[21] += -1.58113883008419*GhatR[2]*rdx2; 
  out[22] += -0.7071067811865475*GhatR[11]*rdx2; 
  out[23] += -1.224744871391589*GhatR[8]*rdx2; 
  out[24] += -0.7071067811865475*GhatR[12]*rdx2; 
  out[25] += -1.58113883008419*GhatR[3]*rdx2; 
  out[26] += -0.7071067811865475*GhatR[13]*rdx2; 
  out[27] += -0.7071067811865475*GhatR[14]*rdx2; 
  out[28] += -1.224744871391589*GhatR[9]*rdx2; 
  out[29] += -0.7071067811865475*GhatR[15]*rdx2; 
  out[30] += -0.7071067811865475*GhatR[16]*rdx2; 
  out[31] += -1.224744871391589*GhatR[10]*rdx2; 
  out[32] += -1.58113883008419*GhatR[4]*rdx2; 
  out[33] += -1.224744871391589*GhatR[11]*rdx2; 
  out[34] += -1.224744871391589*GhatR[12]*rdx2; 
  out[35] += -1.58113883008419*GhatR[5]*rdx2; 
  out[36] += -1.224744871391589*GhatR[13]*rdx2; 
  out[37] += -1.58113883008419*GhatR[6]*rdx2; 
  out[38] += -0.7071067811865475*GhatR[17]*rdx2; 
  out[39] += -1.224744871391589*GhatR[14]*rdx2; 
  out[40] += -0.7071067811865475*GhatR[18]*rdx2; 
  out[41] += -1.224744871391589*GhatR[15]*rdx2; 
  out[42] += -1.224744871391589*GhatR[16]*rdx2; 
  out[43] += -0.7071067811865475*GhatR[19]*rdx2; 
  out[44] += -1.58113883008419*GhatR[10]*rdx2; 
  out[45] += -1.224744871391589*GhatR[17]*rdx2; 
  out[46] += -1.224744871391589*GhatR[18]*rdx2; 
  out[47] += -1.224744871391589*GhatR[19]*rdx2; 

  } else { 

  const double *alphaL = &alpha_skin[0];
  double fUpOrdL[27] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = ser_4x_p2_surfx1_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = ser_4x_p2_surfx1_eval_quad_node_0_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = ser_4x_p2_surfx1_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = ser_4x_p2_surfx1_eval_quad_node_1_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = ser_4x_p2_surfx1_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = ser_4x_p2_surfx1_eval_quad_node_2_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.3952847075210471*alphaL[8])-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = ser_4x_p2_surfx1_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = ser_4x_p2_surfx1_eval_quad_node_3_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.3952847075210471*alphaL[8])-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = ser_4x_p2_surfx1_eval_quad_node_4_r(fedge); 
  } else { 
    fUpOrdL[4] = ser_4x_p2_surfx1_eval_quad_node_4_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.3952847075210471*alphaL[8])-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = ser_4x_p2_surfx1_eval_quad_node_5_r(fedge); 
  } else { 
    fUpOrdL[5] = ser_4x_p2_surfx1_eval_quad_node_5_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = ser_4x_p2_surfx1_eval_quad_node_6_r(fedge); 
  } else { 
    fUpOrdL[6] = ser_4x_p2_surfx1_eval_quad_node_6_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = ser_4x_p2_surfx1_eval_quad_node_7_r(fedge); 
  } else { 
    fUpOrdL[7] = ser_4x_p2_surfx1_eval_quad_node_7_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = ser_4x_p2_surfx1_eval_quad_node_8_r(fedge); 
  } else { 
    fUpOrdL[8] = ser_4x_p2_surfx1_eval_quad_node_8_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[9] = ser_4x_p2_surfx1_eval_quad_node_9_r(fedge); 
  } else { 
    fUpOrdL[9] = ser_4x_p2_surfx1_eval_quad_node_9_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[10] = ser_4x_p2_surfx1_eval_quad_node_10_r(fedge); 
  } else { 
    fUpOrdL[10] = ser_4x_p2_surfx1_eval_quad_node_10_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[11] = ser_4x_p2_surfx1_eval_quad_node_11_r(fedge); 
  } else { 
    fUpOrdL[11] = ser_4x_p2_surfx1_eval_quad_node_11_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0]-0.3952847075210471*alphaL[8];
  if (alphaL_n > 0.) {
    fUpOrdL[12] = ser_4x_p2_surfx1_eval_quad_node_12_r(fedge); 
  } else { 
    fUpOrdL[12] = ser_4x_p2_surfx1_eval_quad_node_12_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0]-0.3952847075210471*alphaL[8];
  if (alphaL_n > 0.) {
    fUpOrdL[13] = ser_4x_p2_surfx1_eval_quad_node_13_r(fedge); 
  } else { 
    fUpOrdL[13] = ser_4x_p2_surfx1_eval_quad_node_13_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0]-0.3952847075210471*alphaL[8];
  if (alphaL_n > 0.) {
    fUpOrdL[14] = ser_4x_p2_surfx1_eval_quad_node_14_r(fedge); 
  } else { 
    fUpOrdL[14] = ser_4x_p2_surfx1_eval_quad_node_14_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[15] = ser_4x_p2_surfx1_eval_quad_node_15_r(fedge); 
  } else { 
    fUpOrdL[15] = ser_4x_p2_surfx1_eval_quad_node_15_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[16] = ser_4x_p2_surfx1_eval_quad_node_16_r(fedge); 
  } else { 
    fUpOrdL[16] = ser_4x_p2_surfx1_eval_quad_node_16_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[17] = ser_4x_p2_surfx1_eval_quad_node_17_r(fedge); 
  } else { 
    fUpOrdL[17] = ser_4x_p2_surfx1_eval_quad_node_17_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[18] = ser_4x_p2_surfx1_eval_quad_node_18_r(fedge); 
  } else { 
    fUpOrdL[18] = ser_4x_p2_surfx1_eval_quad_node_18_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[19] = ser_4x_p2_surfx1_eval_quad_node_19_r(fedge); 
  } else { 
    fUpOrdL[19] = ser_4x_p2_surfx1_eval_quad_node_19_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[20] = ser_4x_p2_surfx1_eval_quad_node_20_r(fedge); 
  } else { 
    fUpOrdL[20] = ser_4x_p2_surfx1_eval_quad_node_20_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.3952847075210471*alphaL[8])+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[21] = ser_4x_p2_surfx1_eval_quad_node_21_r(fedge); 
  } else { 
    fUpOrdL[21] = ser_4x_p2_surfx1_eval_quad_node_21_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.3952847075210471*alphaL[8])+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[22] = ser_4x_p2_surfx1_eval_quad_node_22_r(fedge); 
  } else { 
    fUpOrdL[22] = ser_4x_p2_surfx1_eval_quad_node_22_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.3952847075210471*alphaL[8])+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[23] = ser_4x_p2_surfx1_eval_quad_node_23_r(fedge); 
  } else { 
    fUpOrdL[23] = ser_4x_p2_surfx1_eval_quad_node_23_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[24] = ser_4x_p2_surfx1_eval_quad_node_24_r(fedge); 
  } else { 
    fUpOrdL[24] = ser_4x_p2_surfx1_eval_quad_node_24_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[25] = ser_4x_p2_surfx1_eval_quad_node_25_r(fedge); 
  } else { 
    fUpOrdL[25] = ser_4x_p2_surfx1_eval_quad_node_25_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3162277660168378*alphaL[8]+0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[26] = ser_4x_p2_surfx1_eval_quad_node_26_r(fedge); 
  } else { 
    fUpOrdL[26] = ser_4x_p2_surfx1_eval_quad_node_26_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[20] = {0.};
  ser_4x_p2_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[20] = {0.}; 
  GhatL[0] = 0.3535533905932737*(alphaL[8]*fUpL[8]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.3535533905932737*alphaL[8]*fUpL[12]+0.3162277660168379*alphaL[1]*fUpL[7]+0.3535533905932737*(alphaL[2]*fUpL[4]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.3162277660168379*(alphaL[2]*fUpL[8]+fUpL[2]*alphaL[8])+0.3535533905932737*(alphaL[1]*fUpL[4]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.3535533905932737*(alphaL[8]*fUpL[14]+alphaL[2]*fUpL[6]+alphaL[1]*fUpL[5]+alphaL[0]*fUpL[3]); 
  GhatL[4] = 0.3162277660168379*(alphaL[2]*fUpL[12]+alphaL[1]*fUpL[11])+0.3162277660168379*fUpL[4]*alphaL[8]+0.3535533905932737*(alphaL[0]*fUpL[4]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[5] = 0.3535533905932737*alphaL[8]*fUpL[18]+0.3162277660168379*alphaL[1]*fUpL[13]+0.3535533905932737*(alphaL[2]*fUpL[10]+alphaL[0]*fUpL[5]+alphaL[1]*fUpL[3]); 
  GhatL[6] = 0.3162277660168379*alphaL[2]*fUpL[14]+0.3535533905932737*alphaL[1]*fUpL[10]+0.3162277660168379*fUpL[6]*alphaL[8]+0.3535533905932737*(alphaL[0]*fUpL[6]+alphaL[2]*fUpL[3]); 
  GhatL[7] = 0.3535533905932737*(alphaL[2]*fUpL[11]+alphaL[0]*fUpL[7])+0.3162277660168379*alphaL[1]*fUpL[1]; 
  GhatL[8] = 0.3535533905932737*alphaL[1]*fUpL[12]+0.2258769757263128*alphaL[8]*fUpL[8]+0.3535533905932737*(alphaL[0]*fUpL[8]+fUpL[0]*alphaL[8])+0.3162277660168379*alphaL[2]*fUpL[2]; 
  GhatL[9] = 0.3535533905932737*(alphaL[2]*fUpL[16]+alphaL[1]*fUpL[15]+alphaL[0]*fUpL[9]); 
  GhatL[10] = 0.3162277660168379*(alphaL[2]*fUpL[18]+alphaL[1]*fUpL[17]+alphaL[8]*fUpL[10])+0.3535533905932737*(alphaL[0]*fUpL[10]+alphaL[1]*fUpL[6]+alphaL[2]*fUpL[5]); 
  GhatL[11] = 0.3162277660168379*alphaL[8]*fUpL[11]+0.3535533905932737*(alphaL[0]*fUpL[11]+alphaL[2]*fUpL[7])+0.3162277660168379*alphaL[1]*fUpL[4]; 
  GhatL[12] = 0.2258769757263128*alphaL[8]*fUpL[12]+0.3535533905932737*(alphaL[0]*fUpL[12]+alphaL[1]*fUpL[8]+fUpL[1]*alphaL[8])+0.3162277660168379*alphaL[2]*fUpL[4]; 
  GhatL[13] = 0.3535533905932737*(alphaL[2]*fUpL[17]+alphaL[0]*fUpL[13])+0.3162277660168379*alphaL[1]*fUpL[5]; 
  GhatL[14] = 0.3535533905932737*alphaL[1]*fUpL[18]+0.2258769757263128*alphaL[8]*fUpL[14]+0.3535533905932737*(alphaL[0]*fUpL[14]+fUpL[3]*alphaL[8])+0.3162277660168379*alphaL[2]*fUpL[6]; 
  GhatL[15] = 0.3535533905932737*(alphaL[2]*fUpL[19]+alphaL[0]*fUpL[15]+alphaL[1]*fUpL[9]); 
  GhatL[16] = 0.3535533905932737*alphaL[1]*fUpL[19]+0.3162277660168379*alphaL[8]*fUpL[16]+0.3535533905932737*(alphaL[0]*fUpL[16]+alphaL[2]*fUpL[9]); 
  GhatL[17] = 0.3162277660168379*alphaL[8]*fUpL[17]+0.3535533905932737*(alphaL[0]*fUpL[17]+alphaL[2]*fUpL[13])+0.3162277660168379*alphaL[1]*fUpL[10]; 
  GhatL[18] = 0.2258769757263128*alphaL[8]*fUpL[18]+0.3535533905932737*(alphaL[0]*fUpL[18]+alphaL[1]*fUpL[14])+0.3162277660168379*alphaL[2]*fUpL[10]+0.3535533905932737*fUpL[5]*alphaL[8]; 
  GhatL[19] = 0.3162277660168379*alphaL[8]*fUpL[19]+0.3535533905932737*(alphaL[0]*fUpL[19]+alphaL[1]*fUpL[16]+alphaL[2]*fUpL[15]); 

  out[0] += 0.7071067811865475*GhatL[0]*rdx2; 
  out[1] += -1.224744871391589*GhatL[0]*rdx2; 
  out[2] += 0.7071067811865475*GhatL[1]*rdx2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdx2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdx2; 
  out[5] += -1.224744871391589*GhatL[1]*rdx2; 
  out[6] += -1.224744871391589*GhatL[2]*rdx2; 
  out[7] += 0.7071067811865475*GhatL[4]*rdx2; 
  out[8] += -1.224744871391589*GhatL[3]*rdx2; 
  out[9] += 0.7071067811865475*GhatL[5]*rdx2; 
  out[10] += 0.7071067811865475*GhatL[6]*rdx2; 
  out[11] += 1.58113883008419*GhatL[0]*rdx2; 
  out[12] += 0.7071067811865475*GhatL[7]*rdx2; 
  out[13] += 0.7071067811865475*GhatL[8]*rdx2; 
  out[14] += 0.7071067811865475*GhatL[9]*rdx2; 
  out[15] += -1.224744871391589*GhatL[4]*rdx2; 
  out[16] += -1.224744871391589*GhatL[5]*rdx2; 
  out[17] += -1.224744871391589*GhatL[6]*rdx2; 
  out[18] += 0.7071067811865475*GhatL[10]*rdx2; 
  out[19] += 1.58113883008419*GhatL[1]*rdx2; 
  out[20] += -1.224744871391589*GhatL[7]*rdx2; 
  out[21] += 1.58113883008419*GhatL[2]*rdx2; 
  out[22] += 0.7071067811865475*GhatL[11]*rdx2; 
  out[23] += -1.224744871391589*GhatL[8]*rdx2; 
  out[24] += 0.7071067811865475*GhatL[12]*rdx2; 
  out[25] += 1.58113883008419*GhatL[3]*rdx2; 
  out[26] += 0.7071067811865475*GhatL[13]*rdx2; 
  out[27] += 0.7071067811865475*GhatL[14]*rdx2; 
  out[28] += -1.224744871391589*GhatL[9]*rdx2; 
  out[29] += 0.7071067811865475*GhatL[15]*rdx2; 
  out[30] += 0.7071067811865475*GhatL[16]*rdx2; 
  out[31] += -1.224744871391589*GhatL[10]*rdx2; 
  out[32] += 1.58113883008419*GhatL[4]*rdx2; 
  out[33] += -1.224744871391589*GhatL[11]*rdx2; 
  out[34] += -1.224744871391589*GhatL[12]*rdx2; 
  out[35] += 1.58113883008419*GhatL[5]*rdx2; 
  out[36] += -1.224744871391589*GhatL[13]*rdx2; 
  out[37] += 1.58113883008419*GhatL[6]*rdx2; 
  out[38] += 0.7071067811865475*GhatL[17]*rdx2; 
  out[39] += -1.224744871391589*GhatL[14]*rdx2; 
  out[40] += 0.7071067811865475*GhatL[18]*rdx2; 
  out[41] += -1.224744871391589*GhatL[15]*rdx2; 
  out[42] += -1.224744871391589*GhatL[16]*rdx2; 
  out[43] += 0.7071067811865475*GhatL[19]*rdx2; 
  out[44] += 1.58113883008419*GhatL[10]*rdx2; 
  out[45] += -1.224744871391589*GhatL[17]*rdx2; 
  out[46] += -1.224744871391589*GhatL[18]*rdx2; 
  out[47] += -1.224744871391589*GhatL[19]*rdx2; 

  } 

  return 2.5*rdx2*cflFreq; 

} 
