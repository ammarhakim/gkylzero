#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfy_3x2v_ser_p1(const double *w, const double *dxv, const double *alpha_edge, const double *alpha_skin, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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
  double wz = w[2];
  double rdz2 = 2.0/dxv[2];
  double wvpar = w[3];
  double rdvpar2 = 2.0/dxv[3];
  double wmu = w[4];
  double rdmu2 = 2.0/dxv[4];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wzSq = w[2]*w[2];
  double rdz2Sq = rdz2*rdz2;
  double wvparSq = w[3]*w[3];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[4]*w[4];
  double rdmu2Sq = rdmu2*rdmu2;

  double cflFreq = 0.0;

  if (edge == -1) { 

  const double *alphaR = &alpha_edge[24];
  double fUpOrdR[24] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.2236067977499786*alphaR[20]-0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]-0.3354101966249678*alphaR[11]+0.25*alphaR[9]+0.25*alphaR[8]+0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx2_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx2_eval_quad_node_0_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[20])+0.2795084971874732*alphaR[18]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.25*alphaR[12]+0.25*alphaR[9]+0.25*alphaR[8]+0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx2_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx2_eval_quad_node_1_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[20]-0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]+0.3354101966249678*alphaR[11]+0.25*alphaR[9]+0.25*alphaR[8]-0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx2_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx2_eval_quad_node_2_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[20]-0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]-0.3354101966249678*alphaR[11]-0.25*alphaR[9]-0.25*alphaR[8]+0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx2_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx2_eval_quad_node_3_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[20])+0.2795084971874732*alphaR[18]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]+0.25*alphaR[12]-0.25*alphaR[9]-0.25*alphaR[8]+0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx2_eval_quad_node_4_r(fskin); 
  } else { 
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx2_eval_quad_node_4_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[20]-0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]+0.3354101966249678*alphaR[11]-0.25*alphaR[9]-0.25*alphaR[8]-0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx2_eval_quad_node_5_r(fskin); 
  } else { 
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx2_eval_quad_node_5_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[20])+0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]+0.3354101966249678*alphaR[11]-0.25*alphaR[9]+0.25*alphaR[8]-0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx2_eval_quad_node_6_r(fskin); 
  } else { 
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx2_eval_quad_node_6_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[20]-0.2795084971874732*alphaR[18]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]+0.25*alphaR[12]-0.25*alphaR[9]+0.25*alphaR[8]-0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx2_eval_quad_node_7_r(fskin); 
  } else { 
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx2_eval_quad_node_7_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[20])+0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]-0.3354101966249678*alphaR[11]-0.25*alphaR[9]+0.25*alphaR[8]+0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx2_eval_quad_node_8_r(fskin); 
  } else { 
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx2_eval_quad_node_8_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[20])+0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]+0.3354101966249678*alphaR[11]+0.25*alphaR[9]-0.25*alphaR[8]-0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx2_eval_quad_node_9_r(fskin); 
  } else { 
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx2_eval_quad_node_9_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[20]-0.2795084971874732*alphaR[18]+0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.25*alphaR[12]+0.25*alphaR[9]-0.25*alphaR[8]-0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx2_eval_quad_node_10_r(fskin); 
  } else { 
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx2_eval_quad_node_10_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[20])+0.2236067977499786*alphaR[18]-0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]-0.3354101966249678*alphaR[11]+0.25*alphaR[9]-0.25*alphaR[8]+0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx2_eval_quad_node_11_r(fskin); 
  } else { 
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx2_eval_quad_node_11_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[20])-0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]+0.3354101966249678*alphaR[11]+0.25*alphaR[9]-0.25*alphaR[8]+0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx2_eval_quad_node_12_r(fskin); 
  } else { 
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx2_eval_quad_node_12_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[20]+0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]+0.25*alphaR[12]+0.25*alphaR[9]-0.25*alphaR[8]-0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx2_eval_quad_node_13_r(fskin); 
  } else { 
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx2_eval_quad_node_13_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[20])-0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]-0.3354101966249678*alphaR[11]+0.25*alphaR[9]-0.25*alphaR[8]-0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx2_eval_quad_node_14_r(fskin); 
  } else { 
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx2_eval_quad_node_14_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[20])-0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]+0.3354101966249678*alphaR[11]-0.25*alphaR[9]+0.25*alphaR[8]+0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx2_eval_quad_node_15_r(fskin); 
  } else { 
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx2_eval_quad_node_15_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[20]+0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.25*alphaR[12]-0.25*alphaR[9]+0.25*alphaR[8]-0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[16] = gkhyb_3x2v_p1_surfx2_eval_quad_node_16_r(fskin); 
  } else { 
    fUpOrdR[16] = gkhyb_3x2v_p1_surfx2_eval_quad_node_16_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[20])-0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]-0.3354101966249678*alphaR[11]-0.25*alphaR[9]+0.25*alphaR[8]-0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[17] = gkhyb_3x2v_p1_surfx2_eval_quad_node_17_r(fskin); 
  } else { 
    fUpOrdR[17] = gkhyb_3x2v_p1_surfx2_eval_quad_node_17_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[20]+0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]-0.3354101966249678*alphaR[11]-0.25*alphaR[9]-0.25*alphaR[8]-0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[18] = gkhyb_3x2v_p1_surfx2_eval_quad_node_18_r(fskin); 
  } else { 
    fUpOrdR[18] = gkhyb_3x2v_p1_surfx2_eval_quad_node_18_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[20])-0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]-0.25*alphaR[12]-0.25*alphaR[9]-0.25*alphaR[8]+0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[19] = gkhyb_3x2v_p1_surfx2_eval_quad_node_19_r(fskin); 
  } else { 
    fUpOrdR[19] = gkhyb_3x2v_p1_surfx2_eval_quad_node_19_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[20]+0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]-0.25*alphaR[12]+0.3354101966249678*alphaR[11]-0.25*alphaR[9]-0.25*alphaR[8]+0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[20] = gkhyb_3x2v_p1_surfx2_eval_quad_node_20_r(fskin); 
  } else { 
    fUpOrdR[20] = gkhyb_3x2v_p1_surfx2_eval_quad_node_20_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[20]+0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]-0.3354101966249678*alphaR[11]+0.25*alphaR[9]+0.25*alphaR[8]-0.3354101966249678*alphaR[7]-0.3354101966249678*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[21] = gkhyb_3x2v_p1_surfx2_eval_quad_node_21_r(fskin); 
  } else { 
    fUpOrdR[21] = gkhyb_3x2v_p1_surfx2_eval_quad_node_21_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[20])-0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[17]-0.2795084971874732*alphaR[16]+0.25*alphaR[12]+0.25*alphaR[9]+0.25*alphaR[8]+0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[22] = gkhyb_3x2v_p1_surfx2_eval_quad_node_22_r(fskin); 
  } else { 
    fUpOrdR[22] = gkhyb_3x2v_p1_surfx2_eval_quad_node_22_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[20]+0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[17]+0.2236067977499786*alphaR[16]+0.25*alphaR[12]+0.3354101966249678*alphaR[11]+0.25*alphaR[9]+0.25*alphaR[8]+0.3354101966249678*alphaR[7]+0.3354101966249678*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[23] = gkhyb_3x2v_p1_surfx2_eval_quad_node_23_r(fskin); 
  } else { 
    fUpOrdR[23] = gkhyb_3x2v_p1_surfx2_eval_quad_node_23_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[24] = {0.};
  gkhyb_3x2v_p1_xdir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[24] = {0.}; 
  GhatR[0] = 0.25*(alphaR[20]*fUpR[20]+alphaR[18]*fUpR[18]+alphaR[17]*fUpR[17]+alphaR[16]*fUpR[16]+alphaR[12]*fUpR[12]+alphaR[11]*fUpR[11]+alphaR[9]*fUpR[9]+alphaR[8]*fUpR[8]+alphaR[7]*fUpR[7]+alphaR[6]*fUpR[6]+alphaR[5]*fUpR[5]+alphaR[4]*fUpR[4]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.2500000000000001*(alphaR[18]*fUpR[20]+fUpR[18]*alphaR[20]+alphaR[16]*fUpR[17]+fUpR[16]*alphaR[17])+0.25*(alphaR[9]*fUpR[12]+fUpR[9]*alphaR[12]+alphaR[7]*fUpR[11]+fUpR[7]*alphaR[11]+alphaR[4]*fUpR[8]+fUpR[4]*alphaR[8]+alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6]+alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.2500000000000001*(alphaR[17]*fUpR[20]+fUpR[17]*alphaR[20]+alphaR[16]*fUpR[18]+fUpR[16]*alphaR[18])+0.25*(alphaR[8]*fUpR[12]+fUpR[8]*alphaR[12]+alphaR[6]*fUpR[11]+fUpR[6]*alphaR[11]+alphaR[4]*fUpR[9]+fUpR[4]*alphaR[9]+alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]+alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.223606797749979*(alphaR[11]*fUpR[20]+fUpR[11]*alphaR[20])+0.223606797749979*(alphaR[7]*fUpR[18]+fUpR[7]*alphaR[18]+alphaR[6]*fUpR[17]+fUpR[6]*alphaR[17])+0.223606797749979*(alphaR[3]*fUpR[16]+fUpR[3]*alphaR[16])+0.25*(alphaR[12]*fUpR[15]+alphaR[9]*fUpR[14]+alphaR[8]*fUpR[13]+alphaR[5]*fUpR[11]+fUpR[5]*alphaR[11]+alphaR[4]*fUpR[10]+alphaR[2]*fUpR[7]+fUpR[2]*alphaR[7]+alphaR[1]*fUpR[6]+fUpR[1]*alphaR[6]+alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]); 
  GhatR[4] = 0.2500000000000001*(alphaR[20]*fUpR[23]+alphaR[18]*fUpR[22]+alphaR[17]*fUpR[21]+alphaR[16]*fUpR[19])+0.25*(alphaR[11]*fUpR[15]+alphaR[7]*fUpR[14]+alphaR[6]*fUpR[13]+alphaR[5]*fUpR[12]+fUpR[5]*alphaR[12]+alphaR[3]*fUpR[10]+alphaR[2]*fUpR[9]+fUpR[2]*alphaR[9]+alphaR[1]*fUpR[8]+fUpR[1]*alphaR[8]+alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]); 
  GhatR[5] = 0.25*(alphaR[16]*fUpR[20]+fUpR[16]*alphaR[20]+alphaR[17]*fUpR[18]+fUpR[17]*alphaR[18]+alphaR[4]*fUpR[12]+fUpR[4]*alphaR[12]+alphaR[3]*fUpR[11]+fUpR[3]*alphaR[11]+alphaR[8]*fUpR[9]+fUpR[8]*alphaR[9]+alphaR[6]*fUpR[7]+fUpR[6]*alphaR[7]+alphaR[0]*fUpR[5]+fUpR[0]*alphaR[5]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[6] = 0.223606797749979*(alphaR[7]*fUpR[20]+fUpR[7]*alphaR[20])+0.223606797749979*(alphaR[11]*fUpR[18]+fUpR[11]*alphaR[18]+alphaR[3]*fUpR[17]+fUpR[3]*alphaR[17])+0.223606797749979*(alphaR[6]*fUpR[16]+fUpR[6]*alphaR[16])+0.25*(alphaR[9]*fUpR[15]+alphaR[12]*fUpR[14]+alphaR[4]*fUpR[13]+alphaR[2]*fUpR[11]+fUpR[2]*alphaR[11]+alphaR[8]*fUpR[10]+alphaR[5]*fUpR[7]+fUpR[5]*alphaR[7]+alphaR[0]*fUpR[6]+fUpR[0]*alphaR[6]+alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[7] = 0.223606797749979*(alphaR[6]*fUpR[20]+fUpR[6]*alphaR[20])+0.223606797749979*(alphaR[3]*fUpR[18]+fUpR[3]*alphaR[18]+alphaR[11]*fUpR[17]+fUpR[11]*alphaR[17])+0.223606797749979*(alphaR[7]*fUpR[16]+fUpR[7]*alphaR[16])+0.25*(alphaR[8]*fUpR[15]+alphaR[4]*fUpR[14]+alphaR[12]*fUpR[13]+alphaR[1]*fUpR[11]+fUpR[1]*alphaR[11]+alphaR[9]*fUpR[10]+alphaR[0]*fUpR[7]+fUpR[0]*alphaR[7]+alphaR[5]*fUpR[6]+fUpR[5]*alphaR[6]+alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 
  GhatR[8] = 0.25*(alphaR[18]*fUpR[23]+alphaR[20]*fUpR[22]+alphaR[16]*fUpR[21]+alphaR[17]*fUpR[19]+alphaR[7]*fUpR[15]+alphaR[11]*fUpR[14]+alphaR[3]*fUpR[13]+alphaR[2]*fUpR[12]+fUpR[2]*alphaR[12]+alphaR[6]*fUpR[10]+alphaR[5]*fUpR[9]+fUpR[5]*alphaR[9]+alphaR[0]*fUpR[8]+fUpR[0]*alphaR[8]+alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]); 
  GhatR[9] = 0.25*(alphaR[17]*fUpR[23]+alphaR[16]*fUpR[22]+alphaR[20]*fUpR[21]+alphaR[18]*fUpR[19]+alphaR[6]*fUpR[15]+alphaR[3]*fUpR[14]+alphaR[11]*fUpR[13]+alphaR[1]*fUpR[12]+fUpR[1]*alphaR[12]+alphaR[7]*fUpR[10]+alphaR[0]*fUpR[9]+fUpR[0]*alphaR[9]+alphaR[5]*fUpR[8]+fUpR[5]*alphaR[8]+alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]); 
  GhatR[10] = 0.223606797749979*alphaR[11]*fUpR[23]+0.223606797749979*(alphaR[7]*fUpR[22]+alphaR[6]*fUpR[21]+fUpR[15]*alphaR[20])+0.223606797749979*(alphaR[3]*fUpR[19]+fUpR[14]*alphaR[18]+fUpR[13]*alphaR[17])+0.223606797749979*fUpR[10]*alphaR[16]+0.25*(alphaR[5]*fUpR[15]+alphaR[2]*fUpR[14]+alphaR[1]*fUpR[13]+alphaR[11]*fUpR[12]+fUpR[11]*alphaR[12]+alphaR[0]*fUpR[10]+alphaR[7]*fUpR[9]+fUpR[7]*alphaR[9]+alphaR[6]*fUpR[8]+fUpR[6]*alphaR[8]+alphaR[3]*fUpR[4]+fUpR[3]*alphaR[4]); 
  GhatR[11] = 0.223606797749979*(alphaR[3]*fUpR[20]+fUpR[3]*alphaR[20])+0.223606797749979*(alphaR[6]*fUpR[18]+fUpR[6]*alphaR[18]+alphaR[7]*fUpR[17]+fUpR[7]*alphaR[17])+0.223606797749979*(alphaR[11]*fUpR[16]+fUpR[11]*alphaR[16])+0.25*(alphaR[4]*fUpR[15]+alphaR[8]*fUpR[14]+alphaR[9]*fUpR[13]+fUpR[10]*alphaR[12]+alphaR[0]*fUpR[11]+fUpR[0]*alphaR[11]+alphaR[1]*fUpR[7]+fUpR[1]*alphaR[7]+alphaR[2]*fUpR[6]+fUpR[2]*alphaR[6]+alphaR[3]*fUpR[5]+fUpR[3]*alphaR[5]); 
  GhatR[12] = 0.2500000000000001*(alphaR[16]*fUpR[23]+alphaR[17]*fUpR[22]+alphaR[18]*fUpR[21]+fUpR[19]*alphaR[20])+0.25*(alphaR[3]*fUpR[15]+alphaR[6]*fUpR[14]+alphaR[7]*fUpR[13]+alphaR[0]*fUpR[12]+fUpR[0]*alphaR[12]+fUpR[10]*alphaR[11]+alphaR[1]*fUpR[9]+fUpR[1]*alphaR[9]+alphaR[2]*fUpR[8]+fUpR[2]*alphaR[8]+alphaR[4]*fUpR[5]+fUpR[4]*alphaR[5]); 
  GhatR[13] = 0.223606797749979*alphaR[7]*fUpR[23]+0.223606797749979*(alphaR[11]*fUpR[22]+alphaR[3]*fUpR[21]+fUpR[14]*alphaR[20])+0.223606797749979*(alphaR[6]*fUpR[19]+fUpR[15]*alphaR[18]+fUpR[10]*alphaR[17])+0.223606797749979*fUpR[13]*alphaR[16]+0.25*(alphaR[2]*fUpR[15]+alphaR[5]*fUpR[14]+alphaR[0]*fUpR[13]+alphaR[7]*fUpR[12]+fUpR[7]*alphaR[12]+alphaR[9]*fUpR[11]+fUpR[9]*alphaR[11]+alphaR[1]*fUpR[10]+alphaR[3]*fUpR[8]+fUpR[3]*alphaR[8]+alphaR[4]*fUpR[6]+fUpR[4]*alphaR[6]); 
  GhatR[14] = 0.223606797749979*alphaR[6]*fUpR[23]+0.223606797749979*(alphaR[3]*fUpR[22]+alphaR[11]*fUpR[21]+fUpR[13]*alphaR[20])+0.223606797749979*(alphaR[7]*fUpR[19]+fUpR[10]*alphaR[18]+fUpR[15]*alphaR[17])+0.223606797749979*fUpR[14]*alphaR[16]+0.25*(alphaR[1]*fUpR[15]+alphaR[0]*fUpR[14]+alphaR[5]*fUpR[13]+alphaR[6]*fUpR[12]+fUpR[6]*alphaR[12]+alphaR[8]*fUpR[11]+fUpR[8]*alphaR[11]+alphaR[2]*fUpR[10]+alphaR[3]*fUpR[9]+fUpR[3]*alphaR[9]+alphaR[4]*fUpR[7]+fUpR[4]*alphaR[7]); 
  GhatR[15] = 0.223606797749979*alphaR[3]*fUpR[23]+0.223606797749979*(alphaR[6]*fUpR[22]+alphaR[7]*fUpR[21]+fUpR[10]*alphaR[20])+0.223606797749979*(alphaR[11]*fUpR[19]+fUpR[13]*alphaR[18]+fUpR[14]*alphaR[17])+0.223606797749979*fUpR[15]*alphaR[16]+0.25*(alphaR[0]*fUpR[15]+alphaR[1]*fUpR[14]+alphaR[2]*fUpR[13]+alphaR[3]*fUpR[12]+fUpR[3]*alphaR[12]+alphaR[4]*fUpR[11]+fUpR[4]*alphaR[11]+alphaR[5]*fUpR[10]+alphaR[6]*fUpR[9]+fUpR[6]*alphaR[9]+alphaR[7]*fUpR[8]+fUpR[7]*alphaR[8]); 
  GhatR[16] = 0.2500000000000001*alphaR[12]*fUpR[23]+0.25*(alphaR[9]*fUpR[22]+alphaR[8]*fUpR[21])+0.159719141249985*alphaR[20]*fUpR[20]+0.25*(alphaR[5]*fUpR[20]+fUpR[5]*alphaR[20])+0.2500000000000001*alphaR[4]*fUpR[19]+0.159719141249985*alphaR[18]*fUpR[18]+0.2500000000000001*(alphaR[2]*fUpR[18]+fUpR[2]*alphaR[18])+0.159719141249985*alphaR[17]*fUpR[17]+0.2500000000000001*(alphaR[1]*fUpR[17]+fUpR[1]*alphaR[17])+0.159719141249985*alphaR[16]*fUpR[16]+0.25*(alphaR[0]*fUpR[16]+fUpR[0]*alphaR[16])+0.223606797749979*(alphaR[11]*fUpR[11]+alphaR[7]*fUpR[7]+alphaR[6]*fUpR[6]+alphaR[3]*fUpR[3]); 
  GhatR[17] = 0.25*alphaR[9]*fUpR[23]+0.2500000000000001*(alphaR[12]*fUpR[22]+alphaR[4]*fUpR[21])+(0.159719141249985*alphaR[18]+0.2500000000000001*alphaR[2])*fUpR[20]+(0.159719141249985*fUpR[18]+0.2500000000000001*fUpR[2])*alphaR[20]+0.25*(alphaR[8]*fUpR[19]+alphaR[5]*fUpR[18]+fUpR[5]*alphaR[18])+(0.159719141249985*alphaR[16]+0.25*alphaR[0])*fUpR[17]+(0.159719141249985*fUpR[16]+0.25*fUpR[0])*alphaR[17]+0.2500000000000001*(alphaR[1]*fUpR[16]+fUpR[1]*alphaR[16])+0.223606797749979*(alphaR[7]*fUpR[11]+fUpR[7]*alphaR[11]+alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6]); 
  GhatR[18] = 0.25*alphaR[8]*fUpR[23]+0.2500000000000001*(alphaR[4]*fUpR[22]+alphaR[12]*fUpR[21])+(0.159719141249985*alphaR[17]+0.2500000000000001*alphaR[1])*fUpR[20]+(0.159719141249985*fUpR[17]+0.2500000000000001*fUpR[1])*alphaR[20]+0.25*alphaR[9]*fUpR[19]+(0.159719141249985*alphaR[16]+0.25*alphaR[0])*fUpR[18]+0.159719141249985*fUpR[16]*alphaR[18]+0.25*(fUpR[0]*alphaR[18]+alphaR[5]*fUpR[17]+fUpR[5]*alphaR[17])+0.2500000000000001*(alphaR[2]*fUpR[16]+fUpR[2]*alphaR[16])+0.223606797749979*(alphaR[6]*fUpR[11]+fUpR[6]*alphaR[11]+alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]); 
  GhatR[19] = (0.159719141249985*alphaR[20]+0.25*alphaR[5])*fUpR[23]+(0.159719141249985*alphaR[18]+0.2500000000000001*alphaR[2])*fUpR[22]+0.159719141249985*alphaR[17]*fUpR[21]+0.2500000000000001*(alphaR[1]*fUpR[21]+alphaR[12]*fUpR[20]+fUpR[12]*alphaR[20])+0.159719141249985*alphaR[16]*fUpR[19]+0.25*(alphaR[0]*fUpR[19]+alphaR[9]*fUpR[18]+fUpR[9]*alphaR[18]+alphaR[8]*fUpR[17]+fUpR[8]*alphaR[17])+0.2500000000000001*(alphaR[4]*fUpR[16]+fUpR[4]*alphaR[16])+0.223606797749979*(alphaR[11]*fUpR[15]+alphaR[7]*fUpR[14]+alphaR[6]*fUpR[13]+alphaR[3]*fUpR[10]); 
  GhatR[20] = 0.2500000000000001*alphaR[4]*fUpR[23]+0.25*(alphaR[8]*fUpR[22]+alphaR[9]*fUpR[21])+(0.159719141249985*alphaR[16]+0.25*alphaR[0])*fUpR[20]+(0.159719141249985*fUpR[16]+0.25*fUpR[0])*alphaR[20]+0.2500000000000001*alphaR[12]*fUpR[19]+(0.159719141249985*alphaR[17]+0.2500000000000001*alphaR[1])*fUpR[18]+0.159719141249985*fUpR[17]*alphaR[18]+0.2500000000000001*(fUpR[1]*alphaR[18]+alphaR[2]*fUpR[17]+fUpR[2]*alphaR[17])+0.25*(alphaR[5]*fUpR[16]+fUpR[5]*alphaR[16])+0.223606797749979*(alphaR[3]*fUpR[11]+fUpR[3]*alphaR[11]+alphaR[6]*fUpR[7]+fUpR[6]*alphaR[7]); 
  GhatR[21] = (0.159719141249985*alphaR[18]+0.2500000000000001*alphaR[2])*fUpR[23]+(0.159719141249985*alphaR[20]+0.25*alphaR[5])*fUpR[22]+0.159719141249985*alphaR[16]*fUpR[21]+0.25*(alphaR[0]*fUpR[21]+alphaR[9]*fUpR[20]+fUpR[9]*alphaR[20])+0.159719141249985*alphaR[17]*fUpR[19]+0.2500000000000001*(alphaR[1]*fUpR[19]+alphaR[12]*fUpR[18]+fUpR[12]*alphaR[18]+alphaR[4]*fUpR[17]+fUpR[4]*alphaR[17])+0.25*(alphaR[8]*fUpR[16]+fUpR[8]*alphaR[16])+0.223606797749979*(alphaR[7]*fUpR[15]+alphaR[11]*fUpR[14]+alphaR[3]*fUpR[13]+alphaR[6]*fUpR[10]); 
  GhatR[22] = (0.159719141249985*alphaR[17]+0.2500000000000001*alphaR[1])*fUpR[23]+(0.159719141249985*alphaR[16]+0.25*alphaR[0])*fUpR[22]+0.159719141249985*alphaR[20]*fUpR[21]+0.25*(alphaR[5]*fUpR[21]+alphaR[8]*fUpR[20]+fUpR[8]*alphaR[20])+0.159719141249985*alphaR[18]*fUpR[19]+0.2500000000000001*(alphaR[2]*fUpR[19]+alphaR[4]*fUpR[18]+fUpR[4]*alphaR[18]+alphaR[12]*fUpR[17]+fUpR[12]*alphaR[17])+0.25*(alphaR[9]*fUpR[16]+fUpR[9]*alphaR[16])+0.223606797749979*(alphaR[6]*fUpR[15]+alphaR[3]*fUpR[14]+alphaR[11]*fUpR[13]+alphaR[7]*fUpR[10]); 
  GhatR[23] = (0.159719141249985*alphaR[16]+0.25*alphaR[0])*fUpR[23]+(0.159719141249985*alphaR[17]+0.2500000000000001*alphaR[1])*fUpR[22]+0.159719141249985*alphaR[18]*fUpR[21]+0.2500000000000001*(alphaR[2]*fUpR[21]+alphaR[4]*fUpR[20])+(0.159719141249985*fUpR[19]+0.2500000000000001*fUpR[4])*alphaR[20]+0.25*(alphaR[5]*fUpR[19]+alphaR[8]*fUpR[18]+fUpR[8]*alphaR[18]+alphaR[9]*fUpR[17]+fUpR[9]*alphaR[17])+0.2500000000000001*(alphaR[12]*fUpR[16]+fUpR[12]*alphaR[16])+0.223606797749979*(alphaR[3]*fUpR[15]+alphaR[6]*fUpR[14]+alphaR[7]*fUpR[13]+fUpR[10]*alphaR[11]); 

  out[0] += -0.7071067811865475*GhatR[0]*rdy2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdy2; 
  out[2] += -1.224744871391589*GhatR[0]*rdy2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdy2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdy2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdy2; 
  out[6] += -1.224744871391589*GhatR[1]*rdy2; 
  out[7] += -0.7071067811865475*GhatR[5]*rdy2; 
  out[8] += -1.224744871391589*GhatR[2]*rdy2; 
  out[9] += -0.7071067811865475*GhatR[6]*rdy2; 
  out[10] += -1.224744871391589*GhatR[3]*rdy2; 
  out[11] += -0.7071067811865475*GhatR[7]*rdy2; 
  out[12] += -0.7071067811865475*GhatR[8]*rdy2; 
  out[13] += -1.224744871391589*GhatR[4]*rdy2; 
  out[14] += -0.7071067811865475*GhatR[9]*rdy2; 
  out[15] += -0.7071067811865475*GhatR[10]*rdy2; 
  out[16] += -1.224744871391589*GhatR[5]*rdy2; 
  out[17] += -1.224744871391589*GhatR[6]*rdy2; 
  out[18] += -0.7071067811865475*GhatR[11]*rdy2; 
  out[19] += -1.224744871391589*GhatR[7]*rdy2; 
  out[20] += -1.224744871391589*GhatR[8]*rdy2; 
  out[21] += -0.7071067811865475*GhatR[12]*rdy2; 
  out[22] += -1.224744871391589*GhatR[9]*rdy2; 
  out[23] += -0.7071067811865475*GhatR[13]*rdy2; 
  out[24] += -1.224744871391589*GhatR[10]*rdy2; 
  out[25] += -0.7071067811865475*GhatR[14]*rdy2; 
  out[26] += -1.224744871391589*GhatR[11]*rdy2; 
  out[27] += -1.224744871391589*GhatR[12]*rdy2; 
  out[28] += -1.224744871391589*GhatR[13]*rdy2; 
  out[29] += -0.7071067811865475*GhatR[15]*rdy2; 
  out[30] += -1.224744871391589*GhatR[14]*rdy2; 
  out[31] += -1.224744871391589*GhatR[15]*rdy2; 
  out[32] += -0.7071067811865475*GhatR[16]*rdy2; 
  out[33] += -0.7071067811865475*GhatR[17]*rdy2; 
  out[34] += -1.224744871391589*GhatR[16]*rdy2; 
  out[35] += -0.7071067811865475*GhatR[18]*rdy2; 
  out[36] += -0.7071067811865475*GhatR[19]*rdy2; 
  out[37] += -1.224744871391589*GhatR[17]*rdy2; 
  out[38] += -0.7071067811865475*GhatR[20]*rdy2; 
  out[39] += -1.224744871391589*GhatR[18]*rdy2; 
  out[40] += -0.7071067811865475*GhatR[21]*rdy2; 
  out[41] += -1.224744871391589*GhatR[19]*rdy2; 
  out[42] += -0.7071067811865475*GhatR[22]*rdy2; 
  out[43] += -1.224744871391589*GhatR[20]*rdy2; 
  out[44] += -1.224744871391589*GhatR[21]*rdy2; 
  out[45] += -0.7071067811865475*GhatR[23]*rdy2; 
  out[46] += -1.224744871391589*GhatR[22]*rdy2; 
  out[47] += -1.224744871391589*GhatR[23]*rdy2; 

  } else { 

  const double *alphaL = &alpha_skin[24];
  double fUpOrdL[24] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.2236067977499786*alphaL[20]-0.2236067977499786*alphaL[18]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.25*alphaL[12]-0.3354101966249678*alphaL[11]+0.25*alphaL[9]+0.25*alphaL[8]+0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx2_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx2_eval_quad_node_0_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[20])+0.2795084971874732*alphaL[18]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.25*alphaL[12]+0.25*alphaL[9]+0.25*alphaL[8]+0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx2_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx2_eval_quad_node_1_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[20]-0.2236067977499786*alphaL[18]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.25*alphaL[12]+0.3354101966249678*alphaL[11]+0.25*alphaL[9]+0.25*alphaL[8]-0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx2_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx2_eval_quad_node_2_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[20]-0.2236067977499786*alphaL[18]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.25*alphaL[12]-0.3354101966249678*alphaL[11]-0.25*alphaL[9]-0.25*alphaL[8]+0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx2_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx2_eval_quad_node_3_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[20])+0.2795084971874732*alphaL[18]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]+0.25*alphaL[12]-0.25*alphaL[9]-0.25*alphaL[8]+0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx2_eval_quad_node_4_r(fedge); 
  } else { 
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx2_eval_quad_node_4_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[20]-0.2236067977499786*alphaL[18]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.25*alphaL[12]+0.3354101966249678*alphaL[11]-0.25*alphaL[9]-0.25*alphaL[8]-0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx2_eval_quad_node_5_r(fedge); 
  } else { 
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx2_eval_quad_node_5_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[20])+0.2236067977499786*alphaL[18]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.25*alphaL[12]+0.3354101966249678*alphaL[11]-0.25*alphaL[9]+0.25*alphaL[8]-0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx2_eval_quad_node_6_r(fedge); 
  } else { 
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx2_eval_quad_node_6_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[20]-0.2795084971874732*alphaL[18]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]+0.25*alphaL[12]-0.25*alphaL[9]+0.25*alphaL[8]-0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx2_eval_quad_node_7_r(fedge); 
  } else { 
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx2_eval_quad_node_7_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[20])+0.2236067977499786*alphaL[18]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.25*alphaL[12]-0.3354101966249678*alphaL[11]-0.25*alphaL[9]+0.25*alphaL[8]+0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx2_eval_quad_node_8_r(fedge); 
  } else { 
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx2_eval_quad_node_8_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[20])+0.2236067977499786*alphaL[18]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.25*alphaL[12]+0.3354101966249678*alphaL[11]+0.25*alphaL[9]-0.25*alphaL[8]-0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx2_eval_quad_node_9_r(fedge); 
  } else { 
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx2_eval_quad_node_9_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[20]-0.2795084971874732*alphaL[18]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.25*alphaL[12]+0.25*alphaL[9]-0.25*alphaL[8]-0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx2_eval_quad_node_10_r(fedge); 
  } else { 
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx2_eval_quad_node_10_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[20])+0.2236067977499786*alphaL[18]-0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.25*alphaL[12]-0.3354101966249678*alphaL[11]+0.25*alphaL[9]-0.25*alphaL[8]+0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx2_eval_quad_node_11_r(fedge); 
  } else { 
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx2_eval_quad_node_11_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[20])-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.25*alphaL[12]+0.3354101966249678*alphaL[11]+0.25*alphaL[9]-0.25*alphaL[8]+0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx2_eval_quad_node_12_r(fedge); 
  } else { 
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx2_eval_quad_node_12_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[20]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]+0.25*alphaL[12]+0.25*alphaL[9]-0.25*alphaL[8]-0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx2_eval_quad_node_13_r(fedge); 
  } else { 
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx2_eval_quad_node_13_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[20])-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.25*alphaL[12]-0.3354101966249678*alphaL[11]+0.25*alphaL[9]-0.25*alphaL[8]-0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx2_eval_quad_node_14_r(fedge); 
  } else { 
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx2_eval_quad_node_14_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[20])-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.25*alphaL[12]+0.3354101966249678*alphaL[11]-0.25*alphaL[9]+0.25*alphaL[8]+0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx2_eval_quad_node_15_r(fedge); 
  } else { 
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx2_eval_quad_node_15_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[20]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.25*alphaL[12]-0.25*alphaL[9]+0.25*alphaL[8]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[16] = gkhyb_3x2v_p1_surfx2_eval_quad_node_16_r(fedge); 
  } else { 
    fUpOrdL[16] = gkhyb_3x2v_p1_surfx2_eval_quad_node_16_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[20])-0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.25*alphaL[12]-0.3354101966249678*alphaL[11]-0.25*alphaL[9]+0.25*alphaL[8]-0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[17] = gkhyb_3x2v_p1_surfx2_eval_quad_node_17_r(fedge); 
  } else { 
    fUpOrdL[17] = gkhyb_3x2v_p1_surfx2_eval_quad_node_17_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[20]+0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.25*alphaL[12]-0.3354101966249678*alphaL[11]-0.25*alphaL[9]-0.25*alphaL[8]-0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[18] = gkhyb_3x2v_p1_surfx2_eval_quad_node_18_r(fedge); 
  } else { 
    fUpOrdL[18] = gkhyb_3x2v_p1_surfx2_eval_quad_node_18_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[20])-0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.25*alphaL[12]-0.25*alphaL[9]-0.25*alphaL[8]+0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[19] = gkhyb_3x2v_p1_surfx2_eval_quad_node_19_r(fedge); 
  } else { 
    fUpOrdL[19] = gkhyb_3x2v_p1_surfx2_eval_quad_node_19_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[20]+0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]-0.25*alphaL[12]+0.3354101966249678*alphaL[11]-0.25*alphaL[9]-0.25*alphaL[8]+0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[20] = gkhyb_3x2v_p1_surfx2_eval_quad_node_20_r(fedge); 
  } else { 
    fUpOrdL[20] = gkhyb_3x2v_p1_surfx2_eval_quad_node_20_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[20]+0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.25*alphaL[12]-0.3354101966249678*alphaL[11]+0.25*alphaL[9]+0.25*alphaL[8]-0.3354101966249678*alphaL[7]-0.3354101966249678*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[21] = gkhyb_3x2v_p1_surfx2_eval_quad_node_21_r(fedge); 
  } else { 
    fUpOrdL[21] = gkhyb_3x2v_p1_surfx2_eval_quad_node_21_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[20])-0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]+0.25*alphaL[12]+0.25*alphaL[9]+0.25*alphaL[8]+0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[22] = gkhyb_3x2v_p1_surfx2_eval_quad_node_22_r(fedge); 
  } else { 
    fUpOrdL[22] = gkhyb_3x2v_p1_surfx2_eval_quad_node_22_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[20]+0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[17]+0.2236067977499786*alphaL[16]+0.25*alphaL[12]+0.3354101966249678*alphaL[11]+0.25*alphaL[9]+0.25*alphaL[8]+0.3354101966249678*alphaL[7]+0.3354101966249678*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[23] = gkhyb_3x2v_p1_surfx2_eval_quad_node_23_r(fedge); 
  } else { 
    fUpOrdL[23] = gkhyb_3x2v_p1_surfx2_eval_quad_node_23_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[24] = {0.};
  gkhyb_3x2v_p1_xdir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[24] = {0.}; 
  GhatL[0] = 0.25*(alphaL[20]*fUpL[20]+alphaL[18]*fUpL[18]+alphaL[17]*fUpL[17]+alphaL[16]*fUpL[16]+alphaL[12]*fUpL[12]+alphaL[11]*fUpL[11]+alphaL[9]*fUpL[9]+alphaL[8]*fUpL[8]+alphaL[7]*fUpL[7]+alphaL[6]*fUpL[6]+alphaL[5]*fUpL[5]+alphaL[4]*fUpL[4]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.2500000000000001*(alphaL[18]*fUpL[20]+fUpL[18]*alphaL[20]+alphaL[16]*fUpL[17]+fUpL[16]*alphaL[17])+0.25*(alphaL[9]*fUpL[12]+fUpL[9]*alphaL[12]+alphaL[7]*fUpL[11]+fUpL[7]*alphaL[11]+alphaL[4]*fUpL[8]+fUpL[4]*alphaL[8]+alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6]+alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.2500000000000001*(alphaL[17]*fUpL[20]+fUpL[17]*alphaL[20]+alphaL[16]*fUpL[18]+fUpL[16]*alphaL[18])+0.25*(alphaL[8]*fUpL[12]+fUpL[8]*alphaL[12]+alphaL[6]*fUpL[11]+fUpL[6]*alphaL[11]+alphaL[4]*fUpL[9]+fUpL[4]*alphaL[9]+alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]+alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.223606797749979*(alphaL[11]*fUpL[20]+fUpL[11]*alphaL[20])+0.223606797749979*(alphaL[7]*fUpL[18]+fUpL[7]*alphaL[18]+alphaL[6]*fUpL[17]+fUpL[6]*alphaL[17])+0.223606797749979*(alphaL[3]*fUpL[16]+fUpL[3]*alphaL[16])+0.25*(alphaL[12]*fUpL[15]+alphaL[9]*fUpL[14]+alphaL[8]*fUpL[13]+alphaL[5]*fUpL[11]+fUpL[5]*alphaL[11]+alphaL[4]*fUpL[10]+alphaL[2]*fUpL[7]+fUpL[2]*alphaL[7]+alphaL[1]*fUpL[6]+fUpL[1]*alphaL[6]+alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]); 
  GhatL[4] = 0.2500000000000001*(alphaL[20]*fUpL[23]+alphaL[18]*fUpL[22]+alphaL[17]*fUpL[21]+alphaL[16]*fUpL[19])+0.25*(alphaL[11]*fUpL[15]+alphaL[7]*fUpL[14]+alphaL[6]*fUpL[13]+alphaL[5]*fUpL[12]+fUpL[5]*alphaL[12]+alphaL[3]*fUpL[10]+alphaL[2]*fUpL[9]+fUpL[2]*alphaL[9]+alphaL[1]*fUpL[8]+fUpL[1]*alphaL[8]+alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]); 
  GhatL[5] = 0.25*(alphaL[16]*fUpL[20]+fUpL[16]*alphaL[20]+alphaL[17]*fUpL[18]+fUpL[17]*alphaL[18]+alphaL[4]*fUpL[12]+fUpL[4]*alphaL[12]+alphaL[3]*fUpL[11]+fUpL[3]*alphaL[11]+alphaL[8]*fUpL[9]+fUpL[8]*alphaL[9]+alphaL[6]*fUpL[7]+fUpL[6]*alphaL[7]+alphaL[0]*fUpL[5]+fUpL[0]*alphaL[5]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[6] = 0.223606797749979*(alphaL[7]*fUpL[20]+fUpL[7]*alphaL[20])+0.223606797749979*(alphaL[11]*fUpL[18]+fUpL[11]*alphaL[18]+alphaL[3]*fUpL[17]+fUpL[3]*alphaL[17])+0.223606797749979*(alphaL[6]*fUpL[16]+fUpL[6]*alphaL[16])+0.25*(alphaL[9]*fUpL[15]+alphaL[12]*fUpL[14]+alphaL[4]*fUpL[13]+alphaL[2]*fUpL[11]+fUpL[2]*alphaL[11]+alphaL[8]*fUpL[10]+alphaL[5]*fUpL[7]+fUpL[5]*alphaL[7]+alphaL[0]*fUpL[6]+fUpL[0]*alphaL[6]+alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[7] = 0.223606797749979*(alphaL[6]*fUpL[20]+fUpL[6]*alphaL[20])+0.223606797749979*(alphaL[3]*fUpL[18]+fUpL[3]*alphaL[18]+alphaL[11]*fUpL[17]+fUpL[11]*alphaL[17])+0.223606797749979*(alphaL[7]*fUpL[16]+fUpL[7]*alphaL[16])+0.25*(alphaL[8]*fUpL[15]+alphaL[4]*fUpL[14]+alphaL[12]*fUpL[13]+alphaL[1]*fUpL[11]+fUpL[1]*alphaL[11]+alphaL[9]*fUpL[10]+alphaL[0]*fUpL[7]+fUpL[0]*alphaL[7]+alphaL[5]*fUpL[6]+fUpL[5]*alphaL[6]+alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 
  GhatL[8] = 0.25*(alphaL[18]*fUpL[23]+alphaL[20]*fUpL[22]+alphaL[16]*fUpL[21]+alphaL[17]*fUpL[19]+alphaL[7]*fUpL[15]+alphaL[11]*fUpL[14]+alphaL[3]*fUpL[13]+alphaL[2]*fUpL[12]+fUpL[2]*alphaL[12]+alphaL[6]*fUpL[10]+alphaL[5]*fUpL[9]+fUpL[5]*alphaL[9]+alphaL[0]*fUpL[8]+fUpL[0]*alphaL[8]+alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]); 
  GhatL[9] = 0.25*(alphaL[17]*fUpL[23]+alphaL[16]*fUpL[22]+alphaL[20]*fUpL[21]+alphaL[18]*fUpL[19]+alphaL[6]*fUpL[15]+alphaL[3]*fUpL[14]+alphaL[11]*fUpL[13]+alphaL[1]*fUpL[12]+fUpL[1]*alphaL[12]+alphaL[7]*fUpL[10]+alphaL[0]*fUpL[9]+fUpL[0]*alphaL[9]+alphaL[5]*fUpL[8]+fUpL[5]*alphaL[8]+alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]); 
  GhatL[10] = 0.223606797749979*alphaL[11]*fUpL[23]+0.223606797749979*(alphaL[7]*fUpL[22]+alphaL[6]*fUpL[21]+fUpL[15]*alphaL[20])+0.223606797749979*(alphaL[3]*fUpL[19]+fUpL[14]*alphaL[18]+fUpL[13]*alphaL[17])+0.223606797749979*fUpL[10]*alphaL[16]+0.25*(alphaL[5]*fUpL[15]+alphaL[2]*fUpL[14]+alphaL[1]*fUpL[13]+alphaL[11]*fUpL[12]+fUpL[11]*alphaL[12]+alphaL[0]*fUpL[10]+alphaL[7]*fUpL[9]+fUpL[7]*alphaL[9]+alphaL[6]*fUpL[8]+fUpL[6]*alphaL[8]+alphaL[3]*fUpL[4]+fUpL[3]*alphaL[4]); 
  GhatL[11] = 0.223606797749979*(alphaL[3]*fUpL[20]+fUpL[3]*alphaL[20])+0.223606797749979*(alphaL[6]*fUpL[18]+fUpL[6]*alphaL[18]+alphaL[7]*fUpL[17]+fUpL[7]*alphaL[17])+0.223606797749979*(alphaL[11]*fUpL[16]+fUpL[11]*alphaL[16])+0.25*(alphaL[4]*fUpL[15]+alphaL[8]*fUpL[14]+alphaL[9]*fUpL[13]+fUpL[10]*alphaL[12]+alphaL[0]*fUpL[11]+fUpL[0]*alphaL[11]+alphaL[1]*fUpL[7]+fUpL[1]*alphaL[7]+alphaL[2]*fUpL[6]+fUpL[2]*alphaL[6]+alphaL[3]*fUpL[5]+fUpL[3]*alphaL[5]); 
  GhatL[12] = 0.2500000000000001*(alphaL[16]*fUpL[23]+alphaL[17]*fUpL[22]+alphaL[18]*fUpL[21]+fUpL[19]*alphaL[20])+0.25*(alphaL[3]*fUpL[15]+alphaL[6]*fUpL[14]+alphaL[7]*fUpL[13]+alphaL[0]*fUpL[12]+fUpL[0]*alphaL[12]+fUpL[10]*alphaL[11]+alphaL[1]*fUpL[9]+fUpL[1]*alphaL[9]+alphaL[2]*fUpL[8]+fUpL[2]*alphaL[8]+alphaL[4]*fUpL[5]+fUpL[4]*alphaL[5]); 
  GhatL[13] = 0.223606797749979*alphaL[7]*fUpL[23]+0.223606797749979*(alphaL[11]*fUpL[22]+alphaL[3]*fUpL[21]+fUpL[14]*alphaL[20])+0.223606797749979*(alphaL[6]*fUpL[19]+fUpL[15]*alphaL[18]+fUpL[10]*alphaL[17])+0.223606797749979*fUpL[13]*alphaL[16]+0.25*(alphaL[2]*fUpL[15]+alphaL[5]*fUpL[14]+alphaL[0]*fUpL[13]+alphaL[7]*fUpL[12]+fUpL[7]*alphaL[12]+alphaL[9]*fUpL[11]+fUpL[9]*alphaL[11]+alphaL[1]*fUpL[10]+alphaL[3]*fUpL[8]+fUpL[3]*alphaL[8]+alphaL[4]*fUpL[6]+fUpL[4]*alphaL[6]); 
  GhatL[14] = 0.223606797749979*alphaL[6]*fUpL[23]+0.223606797749979*(alphaL[3]*fUpL[22]+alphaL[11]*fUpL[21]+fUpL[13]*alphaL[20])+0.223606797749979*(alphaL[7]*fUpL[19]+fUpL[10]*alphaL[18]+fUpL[15]*alphaL[17])+0.223606797749979*fUpL[14]*alphaL[16]+0.25*(alphaL[1]*fUpL[15]+alphaL[0]*fUpL[14]+alphaL[5]*fUpL[13]+alphaL[6]*fUpL[12]+fUpL[6]*alphaL[12]+alphaL[8]*fUpL[11]+fUpL[8]*alphaL[11]+alphaL[2]*fUpL[10]+alphaL[3]*fUpL[9]+fUpL[3]*alphaL[9]+alphaL[4]*fUpL[7]+fUpL[4]*alphaL[7]); 
  GhatL[15] = 0.223606797749979*alphaL[3]*fUpL[23]+0.223606797749979*(alphaL[6]*fUpL[22]+alphaL[7]*fUpL[21]+fUpL[10]*alphaL[20])+0.223606797749979*(alphaL[11]*fUpL[19]+fUpL[13]*alphaL[18]+fUpL[14]*alphaL[17])+0.223606797749979*fUpL[15]*alphaL[16]+0.25*(alphaL[0]*fUpL[15]+alphaL[1]*fUpL[14]+alphaL[2]*fUpL[13]+alphaL[3]*fUpL[12]+fUpL[3]*alphaL[12]+alphaL[4]*fUpL[11]+fUpL[4]*alphaL[11]+alphaL[5]*fUpL[10]+alphaL[6]*fUpL[9]+fUpL[6]*alphaL[9]+alphaL[7]*fUpL[8]+fUpL[7]*alphaL[8]); 
  GhatL[16] = 0.2500000000000001*alphaL[12]*fUpL[23]+0.25*(alphaL[9]*fUpL[22]+alphaL[8]*fUpL[21])+0.159719141249985*alphaL[20]*fUpL[20]+0.25*(alphaL[5]*fUpL[20]+fUpL[5]*alphaL[20])+0.2500000000000001*alphaL[4]*fUpL[19]+0.159719141249985*alphaL[18]*fUpL[18]+0.2500000000000001*(alphaL[2]*fUpL[18]+fUpL[2]*alphaL[18])+0.159719141249985*alphaL[17]*fUpL[17]+0.2500000000000001*(alphaL[1]*fUpL[17]+fUpL[1]*alphaL[17])+0.159719141249985*alphaL[16]*fUpL[16]+0.25*(alphaL[0]*fUpL[16]+fUpL[0]*alphaL[16])+0.223606797749979*(alphaL[11]*fUpL[11]+alphaL[7]*fUpL[7]+alphaL[6]*fUpL[6]+alphaL[3]*fUpL[3]); 
  GhatL[17] = 0.25*alphaL[9]*fUpL[23]+0.2500000000000001*(alphaL[12]*fUpL[22]+alphaL[4]*fUpL[21])+(0.159719141249985*alphaL[18]+0.2500000000000001*alphaL[2])*fUpL[20]+(0.159719141249985*fUpL[18]+0.2500000000000001*fUpL[2])*alphaL[20]+0.25*(alphaL[8]*fUpL[19]+alphaL[5]*fUpL[18]+fUpL[5]*alphaL[18])+(0.159719141249985*alphaL[16]+0.25*alphaL[0])*fUpL[17]+(0.159719141249985*fUpL[16]+0.25*fUpL[0])*alphaL[17]+0.2500000000000001*(alphaL[1]*fUpL[16]+fUpL[1]*alphaL[16])+0.223606797749979*(alphaL[7]*fUpL[11]+fUpL[7]*alphaL[11]+alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6]); 
  GhatL[18] = 0.25*alphaL[8]*fUpL[23]+0.2500000000000001*(alphaL[4]*fUpL[22]+alphaL[12]*fUpL[21])+(0.159719141249985*alphaL[17]+0.2500000000000001*alphaL[1])*fUpL[20]+(0.159719141249985*fUpL[17]+0.2500000000000001*fUpL[1])*alphaL[20]+0.25*alphaL[9]*fUpL[19]+(0.159719141249985*alphaL[16]+0.25*alphaL[0])*fUpL[18]+0.159719141249985*fUpL[16]*alphaL[18]+0.25*(fUpL[0]*alphaL[18]+alphaL[5]*fUpL[17]+fUpL[5]*alphaL[17])+0.2500000000000001*(alphaL[2]*fUpL[16]+fUpL[2]*alphaL[16])+0.223606797749979*(alphaL[6]*fUpL[11]+fUpL[6]*alphaL[11]+alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]); 
  GhatL[19] = (0.159719141249985*alphaL[20]+0.25*alphaL[5])*fUpL[23]+(0.159719141249985*alphaL[18]+0.2500000000000001*alphaL[2])*fUpL[22]+0.159719141249985*alphaL[17]*fUpL[21]+0.2500000000000001*(alphaL[1]*fUpL[21]+alphaL[12]*fUpL[20]+fUpL[12]*alphaL[20])+0.159719141249985*alphaL[16]*fUpL[19]+0.25*(alphaL[0]*fUpL[19]+alphaL[9]*fUpL[18]+fUpL[9]*alphaL[18]+alphaL[8]*fUpL[17]+fUpL[8]*alphaL[17])+0.2500000000000001*(alphaL[4]*fUpL[16]+fUpL[4]*alphaL[16])+0.223606797749979*(alphaL[11]*fUpL[15]+alphaL[7]*fUpL[14]+alphaL[6]*fUpL[13]+alphaL[3]*fUpL[10]); 
  GhatL[20] = 0.2500000000000001*alphaL[4]*fUpL[23]+0.25*(alphaL[8]*fUpL[22]+alphaL[9]*fUpL[21])+(0.159719141249985*alphaL[16]+0.25*alphaL[0])*fUpL[20]+(0.159719141249985*fUpL[16]+0.25*fUpL[0])*alphaL[20]+0.2500000000000001*alphaL[12]*fUpL[19]+(0.159719141249985*alphaL[17]+0.2500000000000001*alphaL[1])*fUpL[18]+0.159719141249985*fUpL[17]*alphaL[18]+0.2500000000000001*(fUpL[1]*alphaL[18]+alphaL[2]*fUpL[17]+fUpL[2]*alphaL[17])+0.25*(alphaL[5]*fUpL[16]+fUpL[5]*alphaL[16])+0.223606797749979*(alphaL[3]*fUpL[11]+fUpL[3]*alphaL[11]+alphaL[6]*fUpL[7]+fUpL[6]*alphaL[7]); 
  GhatL[21] = (0.159719141249985*alphaL[18]+0.2500000000000001*alphaL[2])*fUpL[23]+(0.159719141249985*alphaL[20]+0.25*alphaL[5])*fUpL[22]+0.159719141249985*alphaL[16]*fUpL[21]+0.25*(alphaL[0]*fUpL[21]+alphaL[9]*fUpL[20]+fUpL[9]*alphaL[20])+0.159719141249985*alphaL[17]*fUpL[19]+0.2500000000000001*(alphaL[1]*fUpL[19]+alphaL[12]*fUpL[18]+fUpL[12]*alphaL[18]+alphaL[4]*fUpL[17]+fUpL[4]*alphaL[17])+0.25*(alphaL[8]*fUpL[16]+fUpL[8]*alphaL[16])+0.223606797749979*(alphaL[7]*fUpL[15]+alphaL[11]*fUpL[14]+alphaL[3]*fUpL[13]+alphaL[6]*fUpL[10]); 
  GhatL[22] = (0.159719141249985*alphaL[17]+0.2500000000000001*alphaL[1])*fUpL[23]+(0.159719141249985*alphaL[16]+0.25*alphaL[0])*fUpL[22]+0.159719141249985*alphaL[20]*fUpL[21]+0.25*(alphaL[5]*fUpL[21]+alphaL[8]*fUpL[20]+fUpL[8]*alphaL[20])+0.159719141249985*alphaL[18]*fUpL[19]+0.2500000000000001*(alphaL[2]*fUpL[19]+alphaL[4]*fUpL[18]+fUpL[4]*alphaL[18]+alphaL[12]*fUpL[17]+fUpL[12]*alphaL[17])+0.25*(alphaL[9]*fUpL[16]+fUpL[9]*alphaL[16])+0.223606797749979*(alphaL[6]*fUpL[15]+alphaL[3]*fUpL[14]+alphaL[11]*fUpL[13]+alphaL[7]*fUpL[10]); 
  GhatL[23] = (0.159719141249985*alphaL[16]+0.25*alphaL[0])*fUpL[23]+(0.159719141249985*alphaL[17]+0.2500000000000001*alphaL[1])*fUpL[22]+0.159719141249985*alphaL[18]*fUpL[21]+0.2500000000000001*(alphaL[2]*fUpL[21]+alphaL[4]*fUpL[20])+(0.159719141249985*fUpL[19]+0.2500000000000001*fUpL[4])*alphaL[20]+0.25*(alphaL[5]*fUpL[19]+alphaL[8]*fUpL[18]+fUpL[8]*alphaL[18]+alphaL[9]*fUpL[17]+fUpL[9]*alphaL[17])+0.2500000000000001*(alphaL[12]*fUpL[16]+fUpL[12]*alphaL[16])+0.223606797749979*(alphaL[3]*fUpL[15]+alphaL[6]*fUpL[14]+alphaL[7]*fUpL[13]+fUpL[10]*alphaL[11]); 

  out[0] += 0.7071067811865475*GhatL[0]*rdy2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdy2; 
  out[2] += -1.224744871391589*GhatL[0]*rdy2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdy2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdy2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdy2; 
  out[6] += -1.224744871391589*GhatL[1]*rdy2; 
  out[7] += 0.7071067811865475*GhatL[5]*rdy2; 
  out[8] += -1.224744871391589*GhatL[2]*rdy2; 
  out[9] += 0.7071067811865475*GhatL[6]*rdy2; 
  out[10] += -1.224744871391589*GhatL[3]*rdy2; 
  out[11] += 0.7071067811865475*GhatL[7]*rdy2; 
  out[12] += 0.7071067811865475*GhatL[8]*rdy2; 
  out[13] += -1.224744871391589*GhatL[4]*rdy2; 
  out[14] += 0.7071067811865475*GhatL[9]*rdy2; 
  out[15] += 0.7071067811865475*GhatL[10]*rdy2; 
  out[16] += -1.224744871391589*GhatL[5]*rdy2; 
  out[17] += -1.224744871391589*GhatL[6]*rdy2; 
  out[18] += 0.7071067811865475*GhatL[11]*rdy2; 
  out[19] += -1.224744871391589*GhatL[7]*rdy2; 
  out[20] += -1.224744871391589*GhatL[8]*rdy2; 
  out[21] += 0.7071067811865475*GhatL[12]*rdy2; 
  out[22] += -1.224744871391589*GhatL[9]*rdy2; 
  out[23] += 0.7071067811865475*GhatL[13]*rdy2; 
  out[24] += -1.224744871391589*GhatL[10]*rdy2; 
  out[25] += 0.7071067811865475*GhatL[14]*rdy2; 
  out[26] += -1.224744871391589*GhatL[11]*rdy2; 
  out[27] += -1.224744871391589*GhatL[12]*rdy2; 
  out[28] += -1.224744871391589*GhatL[13]*rdy2; 
  out[29] += 0.7071067811865475*GhatL[15]*rdy2; 
  out[30] += -1.224744871391589*GhatL[14]*rdy2; 
  out[31] += -1.224744871391589*GhatL[15]*rdy2; 
  out[32] += 0.7071067811865475*GhatL[16]*rdy2; 
  out[33] += 0.7071067811865475*GhatL[17]*rdy2; 
  out[34] += -1.224744871391589*GhatL[16]*rdy2; 
  out[35] += 0.7071067811865475*GhatL[18]*rdy2; 
  out[36] += 0.7071067811865475*GhatL[19]*rdy2; 
  out[37] += -1.224744871391589*GhatL[17]*rdy2; 
  out[38] += 0.7071067811865475*GhatL[20]*rdy2; 
  out[39] += -1.224744871391589*GhatL[18]*rdy2; 
  out[40] += 0.7071067811865475*GhatL[21]*rdy2; 
  out[41] += -1.224744871391589*GhatL[19]*rdy2; 
  out[42] += 0.7071067811865475*GhatL[22]*rdy2; 
  out[43] += -1.224744871391589*GhatL[20]*rdy2; 
  out[44] += -1.224744871391589*GhatL[21]*rdy2; 
  out[45] += 0.7071067811865475*GhatL[23]*rdy2; 
  out[46] += -1.224744871391589*GhatL[22]*rdy2; 
  out[47] += -1.224744871391589*GhatL[23]*rdy2; 

  } 

  return 1.5*rdy2*cflFreq; 

} 
