#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_ser_5x_p2_surfx4_eval_quad.h> 
#include <gkyl_basis_ser_5x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_3x2v_ser_p2(const double *w, const double *dxv, const double *alpha_edge, const double *alpha_skin, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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

  const double *alphaR = &alpha_edge[144];
  double fUpOrdR[81] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.4024922359499623*alphaR[39]+0.4024922359499623*alphaR[37]+0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]-0.3*alphaR[27]-0.3*alphaR[25]-0.3*alphaR[24]-0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]-0.3*alphaR[20]-0.3*alphaR[19]-0.603738353924943*alphaR[17]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]+0.45*alphaR[8]+0.45*alphaR[7]+0.45*alphaR[6]+0.45*alphaR[5]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = ser_5x_p2_surfx4_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = ser_5x_p2_surfx4_eval_quad_node_0_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]-0.3*alphaR[24]-0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]-0.3*alphaR[20]-0.3*alphaR[19]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[7]+0.45*alphaR[6]+0.45*alphaR[5]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = ser_5x_p2_surfx4_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = ser_5x_p2_surfx4_eval_quad_node_1_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])-0.4024922359499623*alphaR[37]+0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]+0.3*alphaR[27]+0.3*alphaR[25]-0.3*alphaR[24]-0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]-0.3*alphaR[20]-0.3*alphaR[19]+0.603738353924943*alphaR[17]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]-0.45*alphaR[8]+0.45*alphaR[7]+0.45*alphaR[6]+0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = ser_5x_p2_surfx4_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = ser_5x_p2_surfx4_eval_quad_node_2_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[39])-0.5031152949374518*alphaR[34]+0.375*alphaR[27]-0.3*alphaR[25]+0.375*alphaR[24]+0.375*alphaR[23]-0.3*alphaR[20]-0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]+0.45*alphaR[5]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = ser_5x_p2_surfx4_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = ser_5x_p2_surfx4_eval_quad_node_3_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[34])+0.375*alphaR[24]+0.375*alphaR[23]-0.3*alphaR[20]-0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[5]-0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = ser_5x_p2_surfx4_eval_quad_node_4_r(fskin); 
  } else { 
    fUpOrdR[4] = ser_5x_p2_surfx4_eval_quad_node_4_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[39]-0.5031152949374518*alphaR[34]-0.375*alphaR[27]+0.3*alphaR[25]+0.375*alphaR[24]+0.375*alphaR[23]-0.3*alphaR[20]-0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]+0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = ser_5x_p2_surfx4_eval_quad_node_5_r(fskin); 
  } else { 
    fUpOrdR[5] = ser_5x_p2_surfx4_eval_quad_node_5_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]+0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]-0.3*alphaR[27]-0.3*alphaR[25]-0.3*alphaR[24]-0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]-0.3*alphaR[20]-0.3*alphaR[19]+0.603738353924943*alphaR[17]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[7]-0.45*alphaR[6]+0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = ser_5x_p2_surfx4_eval_quad_node_6_r(fskin); 
  } else { 
    fUpOrdR[6] = ser_5x_p2_surfx4_eval_quad_node_6_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]-0.3*alphaR[24]-0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]-0.3*alphaR[20]-0.3*alphaR[19]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[7]-0.45*alphaR[6]+0.45*alphaR[5]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = ser_5x_p2_surfx4_eval_quad_node_7_r(fskin); 
  } else { 
    fUpOrdR[7] = ser_5x_p2_surfx4_eval_quad_node_7_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])+0.4024922359499623*alphaR[37]+0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]+0.3*alphaR[27]+0.3*alphaR[25]-0.3*alphaR[24]-0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]-0.3*alphaR[20]-0.3*alphaR[19]-0.603738353924943*alphaR[17]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]-0.45*alphaR[8]-0.45*alphaR[7]-0.45*alphaR[6]+0.45*alphaR[5]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = ser_5x_p2_surfx4_eval_quad_node_8_r(fskin); 
  } else { 
    fUpOrdR[8] = ser_5x_p2_surfx4_eval_quad_node_8_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]+0.4024922359499623*alphaR[37]-0.5031152949374518*alphaR[33]-0.3*alphaR[27]-0.3*alphaR[25]-0.3*alphaR[23]+0.375*alphaR[22]-0.3*alphaR[21]+0.375*alphaR[20]-0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]+0.45*alphaR[8]+0.45*alphaR[6]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[9] = ser_5x_p2_surfx4_eval_quad_node_9_r(fskin); 
  } else { 
    fUpOrdR[9] = ser_5x_p2_surfx4_eval_quad_node_9_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[33])-0.3*alphaR[23]+0.375*alphaR[22]-0.3*alphaR[21]+0.375*alphaR[20]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[6]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[10] = ser_5x_p2_surfx4_eval_quad_node_10_r(fskin); 
  } else { 
    fUpOrdR[10] = ser_5x_p2_surfx4_eval_quad_node_10_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])-0.4024922359499623*alphaR[37]-0.5031152949374518*alphaR[33]+0.3*alphaR[27]+0.3*alphaR[25]-0.3*alphaR[23]+0.375*alphaR[22]-0.3*alphaR[21]+0.375*alphaR[20]+0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]-0.45*alphaR[8]+0.45*alphaR[6]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[11] = ser_5x_p2_surfx4_eval_quad_node_11_r(fskin); 
  } else { 
    fUpOrdR[11] = ser_5x_p2_surfx4_eval_quad_node_11_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[39])+0.375*alphaR[27]-0.3*alphaR[25]+0.375*alphaR[23]+0.375*alphaR[20]-0.2795084971874732*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[12] = ser_5x_p2_surfx4_eval_quad_node_12_r(fskin); 
  } else { 
    fUpOrdR[12] = ser_5x_p2_surfx4_eval_quad_node_12_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.375*alphaR[23]+0.375*alphaR[20]-0.2795084971874732*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[13] = ser_5x_p2_surfx4_eval_quad_node_13_r(fskin); 
  } else { 
    fUpOrdR[13] = ser_5x_p2_surfx4_eval_quad_node_13_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[39]-0.375*alphaR[27]+0.3*alphaR[25]+0.375*alphaR[23]+0.375*alphaR[20]-0.2795084971874732*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[14] = ser_5x_p2_surfx4_eval_quad_node_14_r(fskin); 
  } else { 
    fUpOrdR[14] = ser_5x_p2_surfx4_eval_quad_node_14_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]+0.5031152949374518*alphaR[33]-0.3*alphaR[27]-0.3*alphaR[25]-0.3*alphaR[23]-0.375*alphaR[22]+0.3*alphaR[21]+0.375*alphaR[20]+0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[6]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[15] = ser_5x_p2_surfx4_eval_quad_node_15_r(fskin); 
  } else { 
    fUpOrdR[15] = ser_5x_p2_surfx4_eval_quad_node_15_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[33]-0.3*alphaR[23]-0.375*alphaR[22]+0.3*alphaR[21]+0.375*alphaR[20]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[6]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[16] = ser_5x_p2_surfx4_eval_quad_node_16_r(fskin); 
  } else { 
    fUpOrdR[16] = ser_5x_p2_surfx4_eval_quad_node_16_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])+0.4024922359499623*alphaR[37]+0.5031152949374518*alphaR[33]+0.3*alphaR[27]+0.3*alphaR[25]-0.3*alphaR[23]-0.375*alphaR[22]+0.3*alphaR[21]+0.375*alphaR[20]-0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]-0.45*alphaR[8]-0.45*alphaR[6]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[17] = ser_5x_p2_surfx4_eval_quad_node_17_r(fskin); 
  } else { 
    fUpOrdR[17] = ser_5x_p2_surfx4_eval_quad_node_17_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]+0.4024922359499623*alphaR[37]-0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]-0.3*alphaR[27]-0.3*alphaR[25]+0.3*alphaR[24]-0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]-0.3*alphaR[20]+0.3*alphaR[19]-0.603738353924943*alphaR[17]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[18] = ser_5x_p2_surfx4_eval_quad_node_18_r(fskin); 
  } else { 
    fUpOrdR[18] = ser_5x_p2_surfx4_eval_quad_node_18_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[34])+0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]+0.3*alphaR[24]-0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]-0.3*alphaR[20]+0.3*alphaR[19]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[19] = ser_5x_p2_surfx4_eval_quad_node_19_r(fskin); 
  } else { 
    fUpOrdR[19] = ser_5x_p2_surfx4_eval_quad_node_19_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])-0.4024922359499623*alphaR[37]-0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]+0.3*alphaR[27]+0.3*alphaR[25]+0.3*alphaR[24]-0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]-0.3*alphaR[20]+0.3*alphaR[19]+0.603738353924943*alphaR[17]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]-0.45*alphaR[8]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[20] = ser_5x_p2_surfx4_eval_quad_node_20_r(fskin); 
  } else { 
    fUpOrdR[20] = ser_5x_p2_surfx4_eval_quad_node_20_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[39])+0.5031152949374518*alphaR[34]+0.375*alphaR[27]-0.3*alphaR[25]-0.375*alphaR[24]+0.375*alphaR[23]-0.3*alphaR[20]+0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[21] = ser_5x_p2_surfx4_eval_quad_node_21_r(fskin); 
  } else { 
    fUpOrdR[21] = ser_5x_p2_surfx4_eval_quad_node_21_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[34]-0.375*alphaR[24]+0.375*alphaR[23]-0.3*alphaR[20]+0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[5]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[22] = ser_5x_p2_surfx4_eval_quad_node_22_r(fskin); 
  } else { 
    fUpOrdR[22] = ser_5x_p2_surfx4_eval_quad_node_22_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[39]+0.5031152949374518*alphaR[34]-0.375*alphaR[27]+0.3*alphaR[25]-0.375*alphaR[24]+0.375*alphaR[23]-0.3*alphaR[20]+0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[23] = ser_5x_p2_surfx4_eval_quad_node_23_r(fskin); 
  } else { 
    fUpOrdR[23] = ser_5x_p2_surfx4_eval_quad_node_23_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]-0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]-0.3*alphaR[27]-0.3*alphaR[25]+0.3*alphaR[24]-0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]-0.3*alphaR[20]+0.3*alphaR[19]+0.603738353924943*alphaR[17]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[8]+0.45*alphaR[7]-0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[24] = ser_5x_p2_surfx4_eval_quad_node_24_r(fskin); 
  } else { 
    fUpOrdR[24] = ser_5x_p2_surfx4_eval_quad_node_24_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[34])-0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]+0.3*alphaR[24]-0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]-0.3*alphaR[20]+0.3*alphaR[19]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[7]-0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[25] = ser_5x_p2_surfx4_eval_quad_node_25_r(fskin); 
  } else { 
    fUpOrdR[25] = ser_5x_p2_surfx4_eval_quad_node_25_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])+0.4024922359499623*alphaR[37]-0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]+0.3*alphaR[27]+0.3*alphaR[25]+0.3*alphaR[24]-0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]-0.3*alphaR[20]+0.3*alphaR[19]-0.603738353924943*alphaR[17]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]-0.45*alphaR[8]+0.45*alphaR[7]-0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]-0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[26] = ser_5x_p2_surfx4_eval_quad_node_26_r(fskin); 
  } else { 
    fUpOrdR[26] = ser_5x_p2_surfx4_eval_quad_node_26_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[37])-0.5031152949374518*alphaR[32]-0.3*alphaR[27]+0.375*alphaR[25]-0.3*alphaR[24]-0.3*alphaR[22]+0.375*alphaR[21]+0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[10]+0.45*alphaR[7]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[27] = ser_5x_p2_surfx4_eval_quad_node_27_r(fskin); 
  } else { 
    fUpOrdR[27] = ser_5x_p2_surfx4_eval_quad_node_27_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[32])-0.3*alphaR[24]-0.3*alphaR[22]+0.375*alphaR[21]+0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[7]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[28] = ser_5x_p2_surfx4_eval_quad_node_28_r(fskin); 
  } else { 
    fUpOrdR[28] = ser_5x_p2_surfx4_eval_quad_node_28_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[37]-0.5031152949374518*alphaR[32]+0.3*alphaR[27]-0.375*alphaR[25]-0.3*alphaR[24]-0.3*alphaR[22]+0.375*alphaR[21]+0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[7]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[29] = ser_5x_p2_surfx4_eval_quad_node_29_r(fskin); 
  } else { 
    fUpOrdR[29] = ser_5x_p2_surfx4_eval_quad_node_29_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.375*alphaR[27]+0.375*alphaR[25]+0.375*alphaR[24]+0.375*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[30] = ser_5x_p2_surfx4_eval_quad_node_30_r(fskin); 
  } else { 
    fUpOrdR[30] = ser_5x_p2_surfx4_eval_quad_node_30_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.375*alphaR[24]+0.375*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[31] = ser_5x_p2_surfx4_eval_quad_node_31_r(fskin); 
  } else { 
    fUpOrdR[31] = ser_5x_p2_surfx4_eval_quad_node_31_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.375*alphaR[27])-0.375*alphaR[25]+0.375*alphaR[24]+0.375*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[32] = ser_5x_p2_surfx4_eval_quad_node_32_r(fskin); 
  } else { 
    fUpOrdR[32] = ser_5x_p2_surfx4_eval_quad_node_32_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[37]+0.5031152949374518*alphaR[32]-0.3*alphaR[27]+0.375*alphaR[25]-0.3*alphaR[24]+0.3*alphaR[22]-0.375*alphaR[21]+0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[10]-0.45*alphaR[7]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[33] = ser_5x_p2_surfx4_eval_quad_node_33_r(fskin); 
  } else { 
    fUpOrdR[33] = ser_5x_p2_surfx4_eval_quad_node_33_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[32]-0.3*alphaR[24]+0.3*alphaR[22]-0.375*alphaR[21]+0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[7]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[34] = ser_5x_p2_surfx4_eval_quad_node_34_r(fskin); 
  } else { 
    fUpOrdR[34] = ser_5x_p2_surfx4_eval_quad_node_34_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[37])+0.5031152949374518*alphaR[32]+0.3*alphaR[27]-0.375*alphaR[25]-0.3*alphaR[24]+0.3*alphaR[22]-0.375*alphaR[21]+0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[10]-0.45*alphaR[7]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[35] = ser_5x_p2_surfx4_eval_quad_node_35_r(fskin); 
  } else { 
    fUpOrdR[35] = ser_5x_p2_surfx4_eval_quad_node_35_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[37])-0.3*alphaR[27]+0.375*alphaR[25]+0.375*alphaR[22]+0.375*alphaR[21]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[10]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[36] = ser_5x_p2_surfx4_eval_quad_node_36_r(fskin); 
  } else { 
    fUpOrdR[36] = ser_5x_p2_surfx4_eval_quad_node_36_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.375*alphaR[22]+0.375*alphaR[21]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]-0.2795084971874732*alphaR[11]-0.3354101966249678*alphaR[3]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[37] = ser_5x_p2_surfx4_eval_quad_node_37_r(fskin); 
  } else { 
    fUpOrdR[37] = ser_5x_p2_surfx4_eval_quad_node_37_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[37]+0.3*alphaR[27]-0.375*alphaR[25]+0.375*alphaR[22]+0.375*alphaR[21]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[10]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[38] = ser_5x_p2_surfx4_eval_quad_node_38_r(fskin); 
  } else { 
    fUpOrdR[38] = ser_5x_p2_surfx4_eval_quad_node_38_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.375*alphaR[27]+0.375*alphaR[25]-0.2795084971874732*alphaR[13]-0.2795084971874732*alphaR[12]-0.2795084971874732*alphaR[11]-0.3354101966249678*alphaR[4]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[39] = ser_5x_p2_surfx4_eval_quad_node_39_r(fskin); 
  } else { 
    fUpOrdR[39] = ser_5x_p2_surfx4_eval_quad_node_39_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[13])-0.2795084971874732*alphaR[12]-0.2795084971874732*alphaR[11]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[40] = ser_5x_p2_surfx4_eval_quad_node_40_r(fskin); 
  } else { 
    fUpOrdR[40] = ser_5x_p2_surfx4_eval_quad_node_40_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.375*alphaR[27])-0.375*alphaR[25]-0.2795084971874732*alphaR[13]-0.2795084971874732*alphaR[12]-0.2795084971874732*alphaR[11]+0.3354101966249678*alphaR[4]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[41] = ser_5x_p2_surfx4_eval_quad_node_41_r(fskin); 
  } else { 
    fUpOrdR[41] = ser_5x_p2_surfx4_eval_quad_node_41_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[37]-0.3*alphaR[27]+0.375*alphaR[25]-0.375*alphaR[22]-0.375*alphaR[21]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[10]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[42] = ser_5x_p2_surfx4_eval_quad_node_42_r(fskin); 
  } else { 
    fUpOrdR[42] = ser_5x_p2_surfx4_eval_quad_node_42_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.375*alphaR[22])-0.375*alphaR[21]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]-0.2795084971874732*alphaR[11]+0.3354101966249678*alphaR[3]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[43] = ser_5x_p2_surfx4_eval_quad_node_43_r(fskin); 
  } else { 
    fUpOrdR[43] = ser_5x_p2_surfx4_eval_quad_node_43_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[37])+0.3*alphaR[27]-0.375*alphaR[25]-0.375*alphaR[22]-0.375*alphaR[21]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[10]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[44] = ser_5x_p2_surfx4_eval_quad_node_44_r(fskin); 
  } else { 
    fUpOrdR[44] = ser_5x_p2_surfx4_eval_quad_node_44_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[37])+0.5031152949374518*alphaR[32]-0.3*alphaR[27]+0.375*alphaR[25]+0.3*alphaR[24]-0.3*alphaR[22]+0.375*alphaR[21]-0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[10]-0.45*alphaR[7]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[45] = ser_5x_p2_surfx4_eval_quad_node_45_r(fskin); 
  } else { 
    fUpOrdR[45] = ser_5x_p2_surfx4_eval_quad_node_45_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[32]+0.3*alphaR[24]-0.3*alphaR[22]+0.375*alphaR[21]-0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[7]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[46] = ser_5x_p2_surfx4_eval_quad_node_46_r(fskin); 
  } else { 
    fUpOrdR[46] = ser_5x_p2_surfx4_eval_quad_node_46_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[37]+0.5031152949374518*alphaR[32]+0.3*alphaR[27]-0.375*alphaR[25]+0.3*alphaR[24]-0.3*alphaR[22]+0.375*alphaR[21]-0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[10]-0.45*alphaR[7]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[47] = ser_5x_p2_surfx4_eval_quad_node_47_r(fskin); 
  } else { 
    fUpOrdR[47] = ser_5x_p2_surfx4_eval_quad_node_47_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.375*alphaR[27]+0.375*alphaR[25]-0.375*alphaR[24]-0.375*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[48] = ser_5x_p2_surfx4_eval_quad_node_48_r(fskin); 
  } else { 
    fUpOrdR[48] = ser_5x_p2_surfx4_eval_quad_node_48_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.375*alphaR[24])-0.375*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[49] = ser_5x_p2_surfx4_eval_quad_node_49_r(fskin); 
  } else { 
    fUpOrdR[49] = ser_5x_p2_surfx4_eval_quad_node_49_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.375*alphaR[27])-0.375*alphaR[25]-0.375*alphaR[24]-0.375*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[50] = ser_5x_p2_surfx4_eval_quad_node_50_r(fskin); 
  } else { 
    fUpOrdR[50] = ser_5x_p2_surfx4_eval_quad_node_50_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[37]-0.5031152949374518*alphaR[32]-0.3*alphaR[27]+0.375*alphaR[25]+0.3*alphaR[24]+0.3*alphaR[22]-0.375*alphaR[21]-0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[7]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[51] = ser_5x_p2_surfx4_eval_quad_node_51_r(fskin); 
  } else { 
    fUpOrdR[51] = ser_5x_p2_surfx4_eval_quad_node_51_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[32])+0.3*alphaR[24]+0.3*alphaR[22]-0.375*alphaR[21]-0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[7]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[52] = ser_5x_p2_surfx4_eval_quad_node_52_r(fskin); 
  } else { 
    fUpOrdR[52] = ser_5x_p2_surfx4_eval_quad_node_52_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[37])-0.5031152949374518*alphaR[32]+0.3*alphaR[27]-0.375*alphaR[25]+0.3*alphaR[24]+0.3*alphaR[22]-0.375*alphaR[21]-0.375*alphaR[19]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]-0.2795084971874732*alphaR[11]+0.45*alphaR[10]+0.45*alphaR[7]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[53] = ser_5x_p2_surfx4_eval_quad_node_53_r(fskin); 
  } else { 
    fUpOrdR[53] = ser_5x_p2_surfx4_eval_quad_node_53_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])+0.4024922359499623*alphaR[37]-0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]-0.3*alphaR[27]-0.3*alphaR[25]-0.3*alphaR[24]+0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]+0.3*alphaR[20]-0.3*alphaR[19]+0.603738353924943*alphaR[17]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]-0.45*alphaR[8]+0.45*alphaR[7]-0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[54] = ser_5x_p2_surfx4_eval_quad_node_54_r(fskin); 
  } else { 
    fUpOrdR[54] = ser_5x_p2_surfx4_eval_quad_node_54_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[34])-0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]-0.3*alphaR[24]+0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]+0.3*alphaR[20]-0.3*alphaR[19]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[7]-0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[55] = ser_5x_p2_surfx4_eval_quad_node_55_r(fskin); 
  } else { 
    fUpOrdR[55] = ser_5x_p2_surfx4_eval_quad_node_55_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]-0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]+0.3*alphaR[27]+0.3*alphaR[25]-0.3*alphaR[24]+0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]+0.3*alphaR[20]-0.3*alphaR[19]-0.603738353924943*alphaR[17]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[8]+0.45*alphaR[7]-0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[56] = ser_5x_p2_surfx4_eval_quad_node_56_r(fskin); 
  } else { 
    fUpOrdR[56] = ser_5x_p2_surfx4_eval_quad_node_56_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[39]+0.5031152949374518*alphaR[34]+0.375*alphaR[27]-0.3*alphaR[25]+0.375*alphaR[24]-0.375*alphaR[23]+0.3*alphaR[20]-0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[57] = ser_5x_p2_surfx4_eval_quad_node_57_r(fskin); 
  } else { 
    fUpOrdR[57] = ser_5x_p2_surfx4_eval_quad_node_57_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[34]+0.375*alphaR[24]-0.375*alphaR[23]+0.3*alphaR[20]-0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[5]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[58] = ser_5x_p2_surfx4_eval_quad_node_58_r(fskin); 
  } else { 
    fUpOrdR[58] = ser_5x_p2_surfx4_eval_quad_node_58_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[39])+0.5031152949374518*alphaR[34]-0.375*alphaR[27]+0.3*alphaR[25]+0.375*alphaR[24]-0.375*alphaR[23]+0.3*alphaR[20]-0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[59] = ser_5x_p2_surfx4_eval_quad_node_59_r(fskin); 
  } else { 
    fUpOrdR[59] = ser_5x_p2_surfx4_eval_quad_node_59_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])-0.4024922359499623*alphaR[37]-0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]-0.3*alphaR[27]-0.3*alphaR[25]-0.3*alphaR[24]+0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]+0.3*alphaR[20]-0.3*alphaR[19]-0.603738353924943*alphaR[17]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]-0.45*alphaR[8]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[60] = ser_5x_p2_surfx4_eval_quad_node_60_r(fskin); 
  } else { 
    fUpOrdR[60] = ser_5x_p2_surfx4_eval_quad_node_60_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[34])+0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]-0.3*alphaR[24]+0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]+0.3*alphaR[20]-0.3*alphaR[19]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[61] = ser_5x_p2_surfx4_eval_quad_node_61_r(fskin); 
  } else { 
    fUpOrdR[61] = ser_5x_p2_surfx4_eval_quad_node_61_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]+0.4024922359499623*alphaR[37]-0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]+0.3*alphaR[27]+0.3*alphaR[25]-0.3*alphaR[24]+0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]+0.3*alphaR[20]-0.3*alphaR[19]+0.603738353924943*alphaR[17]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[7]+0.45*alphaR[6]-0.45*alphaR[5]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]-0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[62] = ser_5x_p2_surfx4_eval_quad_node_62_r(fskin); 
  } else { 
    fUpOrdR[62] = ser_5x_p2_surfx4_eval_quad_node_62_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])+0.4024922359499623*alphaR[37]+0.5031152949374518*alphaR[33]-0.3*alphaR[27]-0.3*alphaR[25]+0.3*alphaR[23]+0.375*alphaR[22]-0.3*alphaR[21]-0.375*alphaR[20]+0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]-0.45*alphaR[8]-0.45*alphaR[6]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[63] = ser_5x_p2_surfx4_eval_quad_node_63_r(fskin); 
  } else { 
    fUpOrdR[63] = ser_5x_p2_surfx4_eval_quad_node_63_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[33]+0.3*alphaR[23]+0.375*alphaR[22]-0.3*alphaR[21]-0.375*alphaR[20]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[6]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[64] = ser_5x_p2_surfx4_eval_quad_node_64_r(fskin); 
  } else { 
    fUpOrdR[64] = ser_5x_p2_surfx4_eval_quad_node_64_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]+0.5031152949374518*alphaR[33]+0.3*alphaR[27]+0.3*alphaR[25]+0.3*alphaR[23]+0.375*alphaR[22]-0.3*alphaR[21]-0.375*alphaR[20]-0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[6]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[65] = ser_5x_p2_surfx4_eval_quad_node_65_r(fskin); 
  } else { 
    fUpOrdR[65] = ser_5x_p2_surfx4_eval_quad_node_65_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[39]+0.375*alphaR[27]-0.3*alphaR[25]-0.375*alphaR[23]-0.375*alphaR[20]-0.2795084971874732*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[66] = ser_5x_p2_surfx4_eval_quad_node_66_r(fskin); 
  } else { 
    fUpOrdR[66] = ser_5x_p2_surfx4_eval_quad_node_66_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.375*alphaR[23])-0.375*alphaR[20]-0.2795084971874732*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[67] = ser_5x_p2_surfx4_eval_quad_node_67_r(fskin); 
  } else { 
    fUpOrdR[67] = ser_5x_p2_surfx4_eval_quad_node_67_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[39])-0.375*alphaR[27]+0.3*alphaR[25]-0.375*alphaR[23]-0.375*alphaR[20]-0.2795084971874732*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[68] = ser_5x_p2_surfx4_eval_quad_node_68_r(fskin); 
  } else { 
    fUpOrdR[68] = ser_5x_p2_surfx4_eval_quad_node_68_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])-0.4024922359499623*alphaR[37]-0.5031152949374518*alphaR[33]-0.3*alphaR[27]-0.3*alphaR[25]+0.3*alphaR[23]-0.375*alphaR[22]+0.3*alphaR[21]-0.375*alphaR[20]-0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]-0.45*alphaR[8]+0.45*alphaR[6]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[69] = ser_5x_p2_surfx4_eval_quad_node_69_r(fskin); 
  } else { 
    fUpOrdR[69] = ser_5x_p2_surfx4_eval_quad_node_69_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[33])+0.3*alphaR[23]-0.375*alphaR[22]+0.3*alphaR[21]-0.375*alphaR[20]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[6]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[70] = ser_5x_p2_surfx4_eval_quad_node_70_r(fskin); 
  } else { 
    fUpOrdR[70] = ser_5x_p2_surfx4_eval_quad_node_70_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]+0.4024922359499623*alphaR[37]-0.5031152949374518*alphaR[33]+0.3*alphaR[27]+0.3*alphaR[25]+0.3*alphaR[23]-0.375*alphaR[22]+0.3*alphaR[21]-0.375*alphaR[20]+0.603738353924943*alphaR[17]+0.2236067977499786*alphaR[13]-0.2795084971874732*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]+0.45*alphaR[8]+0.45*alphaR[6]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[71] = ser_5x_p2_surfx4_eval_quad_node_71_r(fskin); 
  } else { 
    fUpOrdR[71] = ser_5x_p2_surfx4_eval_quad_node_71_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])+0.4024922359499623*alphaR[37]+0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]-0.3*alphaR[27]-0.3*alphaR[25]+0.3*alphaR[24]+0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]+0.3*alphaR[20]+0.3*alphaR[19]+0.603738353924943*alphaR[17]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]-0.45*alphaR[8]-0.45*alphaR[7]-0.45*alphaR[6]+0.45*alphaR[5]-0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[72] = ser_5x_p2_surfx4_eval_quad_node_72_r(fskin); 
  } else { 
    fUpOrdR[72] = ser_5x_p2_surfx4_eval_quad_node_72_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]+0.3*alphaR[24]+0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]+0.3*alphaR[20]+0.3*alphaR[19]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[7]-0.45*alphaR[6]+0.45*alphaR[5]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[73] = ser_5x_p2_surfx4_eval_quad_node_73_r(fskin); 
  } else { 
    fUpOrdR[73] = ser_5x_p2_surfx4_eval_quad_node_73_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]-0.4024922359499623*alphaR[37]+0.4024922359499623*alphaR[34]-0.4024922359499623*alphaR[33]-0.4024922359499623*alphaR[32]+0.3*alphaR[27]+0.3*alphaR[25]+0.3*alphaR[24]+0.3*alphaR[23]-0.3*alphaR[22]-0.3*alphaR[21]+0.3*alphaR[20]+0.3*alphaR[19]-0.603738353924943*alphaR[17]-0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]+0.45*alphaR[8]-0.45*alphaR[7]-0.45*alphaR[6]+0.45*alphaR[5]+0.3354101966249678*alphaR[4]-0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[74] = ser_5x_p2_surfx4_eval_quad_node_74_r(fskin); 
  } else { 
    fUpOrdR[74] = ser_5x_p2_surfx4_eval_quad_node_74_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5031152949374518*alphaR[39]-0.5031152949374518*alphaR[34]+0.375*alphaR[27]-0.3*alphaR[25]-0.375*alphaR[24]-0.375*alphaR[23]+0.3*alphaR[20]+0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[8]+0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[75] = ser_5x_p2_surfx4_eval_quad_node_75_r(fskin); 
  } else { 
    fUpOrdR[75] = ser_5x_p2_surfx4_eval_quad_node_75_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[34])-0.375*alphaR[24]-0.375*alphaR[23]+0.3*alphaR[20]+0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[5]+0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[76] = ser_5x_p2_surfx4_eval_quad_node_76_r(fskin); 
  } else { 
    fUpOrdR[76] = ser_5x_p2_surfx4_eval_quad_node_76_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5031152949374518*alphaR[39])-0.5031152949374518*alphaR[34]-0.375*alphaR[27]+0.3*alphaR[25]-0.375*alphaR[24]-0.375*alphaR[23]+0.3*alphaR[20]+0.3*alphaR[19]-0.2795084971874732*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[8]+0.45*alphaR[5]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[77] = ser_5x_p2_surfx4_eval_quad_node_77_r(fskin); 
  } else { 
    fUpOrdR[77] = ser_5x_p2_surfx4_eval_quad_node_77_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4024922359499623*alphaR[39])-0.4024922359499623*alphaR[37]+0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]-0.3*alphaR[27]-0.3*alphaR[25]+0.3*alphaR[24]+0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]+0.3*alphaR[20]+0.3*alphaR[19]-0.603738353924943*alphaR[17]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]-0.45*alphaR[10]-0.45*alphaR[8]+0.45*alphaR[7]+0.45*alphaR[6]+0.45*alphaR[5]-0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[78] = ser_5x_p2_surfx4_eval_quad_node_78_r(fskin); 
  } else { 
    fUpOrdR[78] = ser_5x_p2_surfx4_eval_quad_node_78_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]+0.3*alphaR[24]+0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]+0.3*alphaR[20]+0.3*alphaR[19]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[7]+0.45*alphaR[6]+0.45*alphaR[5]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[79] = ser_5x_p2_surfx4_eval_quad_node_79_r(fskin); 
  } else { 
    fUpOrdR[79] = ser_5x_p2_surfx4_eval_quad_node_79_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4024922359499623*alphaR[39]+0.4024922359499623*alphaR[37]+0.4024922359499623*alphaR[34]+0.4024922359499623*alphaR[33]+0.4024922359499623*alphaR[32]+0.3*alphaR[27]+0.3*alphaR[25]+0.3*alphaR[24]+0.3*alphaR[23]+0.3*alphaR[22]+0.3*alphaR[21]+0.3*alphaR[20]+0.3*alphaR[19]+0.603738353924943*alphaR[17]+0.603738353924943*alphaR[15]+0.2236067977499786*alphaR[13]+0.2236067977499786*alphaR[12]+0.2236067977499786*alphaR[11]+0.45*alphaR[10]+0.45*alphaR[8]+0.45*alphaR[7]+0.45*alphaR[6]+0.45*alphaR[5]+0.3354101966249678*alphaR[4]+0.3354101966249678*alphaR[3]+0.3354101966249678*alphaR[2]+0.3354101966249678*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[80] = ser_5x_p2_surfx4_eval_quad_node_80_r(fskin); 
  } else { 
    fUpOrdR[80] = ser_5x_p2_surfx4_eval_quad_node_80_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[48] = {0.};
  ser_5x_p2_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[48] = {0.}; 
  GhatR[0] = 0.25*(alphaR[39]*fUpR[39]+alphaR[37]*fUpR[37]+alphaR[34]*fUpR[34]+alphaR[33]*fUpR[33]+alphaR[32]*fUpR[32]+alphaR[27]*fUpR[27]+alphaR[25]*fUpR[25]+alphaR[24]*fUpR[24]+alphaR[23]*fUpR[23]+alphaR[22]*fUpR[22]+alphaR[21]*fUpR[21]+alphaR[20]*fUpR[20]+alphaR[19]*fUpR[19]+alphaR[17]*fUpR[17]+alphaR[15]*fUpR[15]+alphaR[13]*fUpR[13]+alphaR[12]*fUpR[12]+alphaR[11]*fUpR[11]+alphaR[10]*fUpR[10]+alphaR[8]*fUpR[8]+alphaR[7]*fUpR[7]+alphaR[6]*fUpR[6]+alphaR[5]*fUpR[5]+alphaR[4]*fUpR[4]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.2500000000000001*(alphaR[27]*fUpR[39]+fUpR[27]*alphaR[39])+0.223606797749979*(alphaR[17]*fUpR[37]+fUpR[17]*alphaR[37])+0.2500000000000001*(alphaR[24]*fUpR[34]+fUpR[24]*alphaR[34]+alphaR[22]*fUpR[33]+fUpR[22]*alphaR[33])+0.223606797749979*(alphaR[15]*fUpR[32]+fUpR[15]*alphaR[32])+0.223606797749979*(alphaR[8]*fUpR[25]+fUpR[8]*alphaR[25])+0.2500000000000001*(alphaR[13]*fUpR[23]+fUpR[13]*alphaR[23])+0.223606797749979*(alphaR[6]*fUpR[21]+fUpR[6]*alphaR[21])+0.2500000000000001*(alphaR[12]*fUpR[20]+fUpR[12]*alphaR[20])+0.223606797749979*(alphaR[5]*fUpR[19]+fUpR[5]*alphaR[19])+0.25*(alphaR[10]*fUpR[17]+fUpR[10]*alphaR[17]+alphaR[7]*fUpR[15]+fUpR[7]*alphaR[15])+0.223606797749979*(alphaR[1]*fUpR[11]+fUpR[1]*alphaR[11])+0.25*(alphaR[4]*fUpR[8]+fUpR[4]*alphaR[8]+alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6]+alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.2500000000000001*(alphaR[39]*fUpR[46]+alphaR[37]*fUpR[44]+alphaR[27]*fUpR[40]+alphaR[25]*fUpR[35]+alphaR[23]*fUpR[34]+fUpR[23]*alphaR[34])+0.223606797749979*(alphaR[15]*fUpR[33]+fUpR[15]*alphaR[33])+0.2500000000000001*(alphaR[21]*fUpR[32]+fUpR[21]*alphaR[32])+0.25*alphaR[17]*fUpR[31]+0.2500000000000001*(alphaR[13]*fUpR[24]+fUpR[13]*alphaR[24])+0.223606797749979*(alphaR[7]*fUpR[22]+fUpR[7]*alphaR[22]+alphaR[5]*fUpR[20]+fUpR[5]*alphaR[20])+0.2500000000000001*(alphaR[11]*fUpR[19]+fUpR[11]*alphaR[19])+0.25*(alphaR[10]*fUpR[18]+alphaR[8]*fUpR[16]+alphaR[6]*fUpR[15]+fUpR[6]*alphaR[15])+0.223606797749979*(alphaR[2]*fUpR[12]+fUpR[2]*alphaR[12])+0.25*(alphaR[4]*fUpR[9]+alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]+alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.223606797749979*(alphaR[17]*fUpR[39]+fUpR[17]*alphaR[39])+0.2500000000000001*(alphaR[25]*fUpR[37]+fUpR[25]*alphaR[37])+0.223606797749979*(alphaR[15]*fUpR[34]+fUpR[15]*alphaR[34])+0.2500000000000001*(alphaR[20]*fUpR[33]+fUpR[20]*alphaR[33]+alphaR[19]*fUpR[32]+fUpR[19]*alphaR[32])+0.223606797749979*(alphaR[10]*fUpR[27]+fUpR[10]*alphaR[27]+alphaR[7]*fUpR[24]+fUpR[7]*alphaR[24]+alphaR[6]*fUpR[23]+fUpR[6]*alphaR[23])+0.2500000000000001*(alphaR[12]*fUpR[22]+fUpR[12]*alphaR[22]+alphaR[11]*fUpR[21]+fUpR[11]*alphaR[21])+0.25*(alphaR[8]*fUpR[17]+fUpR[8]*alphaR[17]+alphaR[5]*fUpR[15]+fUpR[5]*alphaR[15])+0.223606797749979*(alphaR[3]*fUpR[13]+fUpR[3]*alphaR[13])+0.25*(alphaR[4]*fUpR[10]+fUpR[4]*alphaR[10]+alphaR[2]*fUpR[7]+fUpR[2]*alphaR[7]+alphaR[1]*fUpR[6]+fUpR[1]*alphaR[6]+alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]); 
  GhatR[4] = 0.2500000000000001*(alphaR[34]*fUpR[46]+alphaR[33]*fUpR[45]+alphaR[32]*fUpR[44])+0.223606797749979*alphaR[17]*fUpR[42]+0.2500000000000001*(alphaR[24]*fUpR[40]+alphaR[23]*fUpR[39]+fUpR[23]*alphaR[39]+alphaR[22]*fUpR[38]+alphaR[21]*fUpR[37]+fUpR[21]*alphaR[37]+alphaR[20]*fUpR[36]+alphaR[19]*fUpR[35])+0.25*alphaR[15]*fUpR[31]+0.223606797749979*(alphaR[10]*fUpR[30]+alphaR[8]*fUpR[28])+0.2500000000000001*(alphaR[13]*fUpR[27]+fUpR[13]*alphaR[27]+alphaR[12]*fUpR[26]+alphaR[11]*fUpR[25]+fUpR[11]*alphaR[25])+0.25*(alphaR[7]*fUpR[18]+alphaR[6]*fUpR[17]+fUpR[6]*alphaR[17]+alphaR[5]*fUpR[16])+0.223606797749979*alphaR[4]*fUpR[14]+0.25*(alphaR[3]*fUpR[10]+fUpR[3]*alphaR[10]+alphaR[2]*fUpR[9]+alphaR[1]*fUpR[8]+fUpR[1]*alphaR[8]+alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]); 
  GhatR[5] = 0.25*alphaR[27]*fUpR[46]+0.223606797749979*alphaR[17]*fUpR[44]+0.25*alphaR[39]*fUpR[40]+0.223606797749979*(fUpR[31]*alphaR[37]+alphaR[8]*fUpR[35])+0.25*(alphaR[13]*fUpR[34]+fUpR[13]*alphaR[34])+(0.2*alphaR[32]+0.223606797749979*alphaR[7])*fUpR[33]+0.2*fUpR[32]*alphaR[33]+0.223606797749979*(fUpR[7]*alphaR[33]+alphaR[6]*fUpR[32]+fUpR[6]*alphaR[32])+0.25*alphaR[10]*fUpR[31]+0.223606797749979*fUpR[16]*alphaR[25]+0.25*(alphaR[23]*fUpR[24]+fUpR[23]*alphaR[24])+0.223606797749979*(alphaR[15]*fUpR[22]+fUpR[15]*alphaR[22]+alphaR[15]*fUpR[21]+fUpR[15]*alphaR[21])+(0.2*alphaR[19]+0.223606797749979*alphaR[2])*fUpR[20]+0.2*fUpR[19]*alphaR[20]+0.223606797749979*(fUpR[2]*alphaR[20]+alphaR[1]*fUpR[19]+fUpR[1]*alphaR[19])+0.25*(alphaR[17]*fUpR[18]+alphaR[4]*fUpR[16]+alphaR[3]*fUpR[15]+fUpR[3]*alphaR[15])+0.223606797749979*(alphaR[5]*fUpR[12]+fUpR[5]*alphaR[12]+alphaR[5]*fUpR[11]+fUpR[5]*alphaR[11])+0.25*(alphaR[8]*fUpR[9]+alphaR[6]*fUpR[7]+fUpR[6]*alphaR[7]+alphaR[0]*fUpR[5]+fUpR[0]*alphaR[5]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[6] = (0.2*alphaR[37]+0.223606797749979*alphaR[10])*fUpR[39]+0.2*fUpR[37]*alphaR[39]+0.223606797749979*(fUpR[10]*alphaR[39]+alphaR[8]*fUpR[37]+fUpR[8]*alphaR[37])+(0.2*alphaR[32]+0.223606797749979*alphaR[7])*fUpR[34]+(0.2*fUpR[32]+0.223606797749979*fUpR[7])*alphaR[34]+0.25*(alphaR[12]*fUpR[33]+fUpR[12]*alphaR[33])+0.223606797749979*(alphaR[5]*fUpR[32]+fUpR[5]*alphaR[32])+0.223606797749979*(alphaR[17]*fUpR[27]+fUpR[17]*alphaR[27]+alphaR[17]*fUpR[25]+fUpR[17]*alphaR[25]+alphaR[15]*fUpR[24]+fUpR[15]*alphaR[24])+(0.2*alphaR[21]+0.223606797749979*alphaR[3])*fUpR[23]+(0.2*fUpR[21]+0.223606797749979*fUpR[3])*alphaR[23]+0.25*(alphaR[20]*fUpR[22]+fUpR[20]*alphaR[22])+0.223606797749979*(alphaR[1]*fUpR[21]+fUpR[1]*alphaR[21]+alphaR[15]*fUpR[19]+fUpR[15]*alphaR[19])+0.25*(alphaR[4]*fUpR[17]+fUpR[4]*alphaR[17]+alphaR[2]*fUpR[15]+fUpR[2]*alphaR[15])+0.223606797749979*(alphaR[6]*fUpR[13]+fUpR[6]*alphaR[13]+alphaR[6]*fUpR[11]+fUpR[6]*alphaR[11])+0.25*(alphaR[8]*fUpR[10]+fUpR[8]*alphaR[10]+alphaR[5]*fUpR[7]+fUpR[5]*alphaR[7]+alphaR[0]*fUpR[6]+fUpR[0]*alphaR[6]+alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[7] = 0.223606797749979*alphaR[17]*fUpR[46]+0.25*alphaR[25]*fUpR[44]+0.223606797749979*(alphaR[10]*fUpR[40]+fUpR[31]*alphaR[39])+0.25*fUpR[35]*alphaR[37]+(0.2*alphaR[33]+0.223606797749979*alphaR[6])*fUpR[34]+0.2*fUpR[33]*alphaR[34]+0.223606797749979*(fUpR[6]*alphaR[34]+alphaR[5]*fUpR[33]+fUpR[5]*alphaR[33])+0.25*(alphaR[11]*fUpR[32]+fUpR[11]*alphaR[32]+alphaR[8]*fUpR[31])+0.223606797749979*fUpR[18]*alphaR[27]+(0.2*alphaR[22]+0.223606797749979*alphaR[3])*fUpR[24]+0.2*fUpR[22]*alphaR[24]+0.223606797749979*(fUpR[3]*alphaR[24]+alphaR[15]*fUpR[23]+fUpR[15]*alphaR[23]+alphaR[2]*fUpR[22]+fUpR[2]*alphaR[22])+0.25*(alphaR[19]*fUpR[21]+fUpR[19]*alphaR[21])+0.223606797749979*(alphaR[15]*fUpR[20]+fUpR[15]*alphaR[20])+0.25*(alphaR[4]*fUpR[18]+fUpR[16]*alphaR[17]+alphaR[1]*fUpR[15]+fUpR[1]*alphaR[15])+0.223606797749979*(alphaR[7]*fUpR[13]+fUpR[7]*alphaR[13]+alphaR[7]*fUpR[12]+fUpR[7]*alphaR[12])+0.25*(fUpR[9]*alphaR[10]+alphaR[0]*fUpR[7]+fUpR[0]*alphaR[7]+alphaR[5]*fUpR[6]+fUpR[5]*alphaR[6]+alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 
  GhatR[8] = 0.25*(alphaR[24]*fUpR[46]+alphaR[22]*fUpR[45])+0.223606797749979*alphaR[15]*fUpR[44]+(0.2*alphaR[37]+0.223606797749979*alphaR[10])*fUpR[42]+0.25*(alphaR[34]*fUpR[40]+alphaR[13]*fUpR[39]+fUpR[13]*alphaR[39]+alphaR[33]*fUpR[38])+0.223606797749979*(alphaR[6]*fUpR[37]+fUpR[6]*alphaR[37])+0.25*alphaR[12]*fUpR[36]+0.223606797749979*alphaR[5]*fUpR[35]+fUpR[31]*(0.223606797749979*alphaR[32]+0.25*alphaR[7])+0.223606797749979*alphaR[17]*fUpR[30]+(0.2*alphaR[25]+0.223606797749979*alphaR[4])*fUpR[28]+0.25*(alphaR[23]*fUpR[27]+fUpR[23]*alphaR[27]+alphaR[20]*fUpR[26])+0.223606797749979*(alphaR[1]*fUpR[25]+fUpR[1]*alphaR[25]+alphaR[17]*fUpR[21]+fUpR[17]*alphaR[21]+fUpR[16]*alphaR[19])+0.25*(alphaR[15]*fUpR[18]+alphaR[3]*fUpR[17]+fUpR[3]*alphaR[17]+alphaR[2]*fUpR[16])+0.223606797749979*(alphaR[8]*(fUpR[14]+fUpR[11])+fUpR[8]*alphaR[11])+0.25*(alphaR[6]*fUpR[10]+fUpR[6]*alphaR[10]+alphaR[5]*fUpR[9]+alphaR[0]*fUpR[8]+fUpR[0]*alphaR[8]+alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]); 
  GhatR[9] = 0.223606797749979*alphaR[17]*fUpR[47]+0.25*alphaR[23]*fUpR[46]+0.223606797749979*alphaR[15]*fUpR[45]+0.25*alphaR[21]*fUpR[44]+0.223606797749979*(alphaR[10]*fUpR[43]+alphaR[8]*fUpR[41])+0.25*(alphaR[13]*fUpR[40]+alphaR[34]*fUpR[39]+fUpR[34]*alphaR[39])+0.223606797749979*alphaR[7]*fUpR[38]+0.25*(alphaR[32]*fUpR[37]+fUpR[32]*alphaR[37])+0.223606797749979*alphaR[5]*fUpR[36]+0.25*alphaR[11]*fUpR[35]+fUpR[31]*(0.223606797749979*alphaR[33]+0.25*alphaR[6])+0.223606797749979*alphaR[4]*fUpR[29]+0.25*(alphaR[24]*fUpR[27]+fUpR[24]*alphaR[27])+0.223606797749979*alphaR[2]*fUpR[26]+0.25*(alphaR[19]*fUpR[25]+fUpR[19]*alphaR[25])+0.223606797749979*(fUpR[18]*alphaR[22]+fUpR[16]*alphaR[20])+0.25*(alphaR[3]*fUpR[18]+alphaR[15]*fUpR[17]+fUpR[15]*alphaR[17]+alphaR[1]*fUpR[16])+0.223606797749979*fUpR[9]*alphaR[12]+0.25*(alphaR[7]*fUpR[10]+fUpR[7]*alphaR[10]+alphaR[0]*fUpR[9]+alphaR[5]*fUpR[8]+fUpR[5]*alphaR[8]+alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]); 
  GhatR[10] = 0.223606797749979*alphaR[15]*fUpR[46]+0.25*(alphaR[20]*fUpR[45]+alphaR[19]*fUpR[44])+0.2*alphaR[39]*fUpR[42]+0.223606797749979*(alphaR[8]*fUpR[42]+alphaR[7]*fUpR[40]+alphaR[6]*fUpR[39]+fUpR[6]*alphaR[39])+0.25*(alphaR[12]*fUpR[38]+alphaR[11]*fUpR[37]+fUpR[11]*alphaR[37]+alphaR[33]*fUpR[36]+alphaR[32]*fUpR[35])+fUpR[31]*(0.223606797749979*alphaR[34]+0.25*alphaR[5])+0.2*alphaR[27]*fUpR[30]+0.223606797749979*(alphaR[4]*fUpR[30]+alphaR[17]*fUpR[28]+alphaR[3]*fUpR[27]+fUpR[3]*alphaR[27])+0.25*(alphaR[22]*fUpR[26]+alphaR[21]*fUpR[25]+fUpR[21]*alphaR[25])+0.223606797749979*(fUpR[18]*alphaR[24]+alphaR[17]*fUpR[23]+fUpR[17]*alphaR[23])+0.25*(alphaR[2]*fUpR[18]+alphaR[1]*fUpR[17]+fUpR[1]*alphaR[17]+alphaR[15]*fUpR[16])+0.223606797749979*(alphaR[10]*(fUpR[14]+fUpR[13])+fUpR[10]*alphaR[13])+0.25*(alphaR[0]*fUpR[10]+fUpR[0]*alphaR[10]+alphaR[7]*fUpR[9]+alphaR[6]*fUpR[8]+fUpR[6]*alphaR[8]+alphaR[3]*fUpR[4]+fUpR[3]*alphaR[4]); 
  GhatR[11] = 0.223606797749979*alphaR[39]*fUpR[39]+0.159719141249985*alphaR[37]*fUpR[37]+0.25*(alphaR[10]*fUpR[37]+fUpR[10]*alphaR[37])+0.223606797749979*(alphaR[34]*fUpR[34]+alphaR[33]*fUpR[33])+0.159719141249985*alphaR[32]*fUpR[32]+0.25*(alphaR[7]*fUpR[32]+fUpR[7]*alphaR[32])+0.159719141249985*alphaR[25]*fUpR[25]+0.2500000000000001*(alphaR[4]*fUpR[25]+fUpR[4]*alphaR[25])+0.223606797749979*alphaR[23]*fUpR[23]+0.159719141249985*alphaR[21]*fUpR[21]+0.2500000000000001*(alphaR[3]*fUpR[21]+fUpR[3]*alphaR[21])+0.223606797749979*alphaR[20]*fUpR[20]+0.159719141249985*alphaR[19]*fUpR[19]+0.2500000000000001*(alphaR[2]*fUpR[19]+fUpR[2]*alphaR[19])+0.223606797749979*(alphaR[17]*fUpR[17]+alphaR[15]*fUpR[15])+0.159719141249985*alphaR[11]*fUpR[11]+0.25*(alphaR[0]*fUpR[11]+fUpR[0]*alphaR[11])+0.223606797749979*(alphaR[8]*fUpR[8]+alphaR[6]*fUpR[6]+alphaR[5]*fUpR[5]+alphaR[1]*fUpR[1]); 
  GhatR[12] = 0.2500000000000001*alphaR[17]*fUpR[45]+0.25*(alphaR[10]*fUpR[38]+alphaR[8]*fUpR[36])+0.223606797749979*alphaR[34]*fUpR[34]+0.159719141249985*alphaR[33]*fUpR[33]+0.25*(alphaR[6]*fUpR[33]+fUpR[6]*alphaR[33])+0.223606797749979*alphaR[32]*fUpR[32]+0.2500000000000001*alphaR[4]*fUpR[26]+0.223606797749979*alphaR[24]*fUpR[24]+0.159719141249985*alphaR[22]*fUpR[22]+0.2500000000000001*(alphaR[3]*fUpR[22]+fUpR[3]*alphaR[22])+0.159719141249985*alphaR[20]*fUpR[20]+0.2500000000000001*(alphaR[1]*fUpR[20]+fUpR[1]*alphaR[20])+0.223606797749979*(alphaR[19]*fUpR[19]+alphaR[15]*fUpR[15])+0.159719141249985*alphaR[12]*fUpR[12]+0.25*(alphaR[0]*fUpR[12]+fUpR[0]*alphaR[12])+0.223606797749979*(alphaR[7]*fUpR[7]+alphaR[5]*fUpR[5]+alphaR[2]*fUpR[2]); 
  GhatR[13] = 0.159719141249985*alphaR[39]*fUpR[39]+0.25*(alphaR[8]*fUpR[39]+fUpR[8]*alphaR[39])+0.223606797749979*alphaR[37]*fUpR[37]+0.159719141249985*alphaR[34]*fUpR[34]+0.25*(alphaR[5]*fUpR[34]+fUpR[5]*alphaR[34])+0.223606797749979*(alphaR[33]*fUpR[33]+alphaR[32]*fUpR[32])+0.159719141249985*alphaR[27]*fUpR[27]+0.2500000000000001*(alphaR[4]*fUpR[27]+fUpR[4]*alphaR[27])+0.159719141249985*alphaR[24]*fUpR[24]+0.2500000000000001*(alphaR[2]*fUpR[24]+fUpR[2]*alphaR[24])+0.159719141249985*alphaR[23]*fUpR[23]+0.2500000000000001*(alphaR[1]*fUpR[23]+fUpR[1]*alphaR[23])+0.223606797749979*(alphaR[22]*fUpR[22]+alphaR[21]*fUpR[21]+alphaR[17]*fUpR[17]+alphaR[15]*fUpR[15])+0.159719141249985*alphaR[13]*fUpR[13]+0.25*(alphaR[0]*fUpR[13]+fUpR[0]*alphaR[13])+0.223606797749979*(alphaR[10]*fUpR[10]+alphaR[7]*fUpR[7]+alphaR[6]*fUpR[6]+alphaR[3]*fUpR[3]); 
  GhatR[14] = 0.2500000000000001*alphaR[15]*fUpR[47]+0.25*(alphaR[7]*fUpR[43]+alphaR[6]*fUpR[42]+alphaR[5]*fUpR[41])+0.223606797749979*(alphaR[39]*fUpR[39]+alphaR[37]*fUpR[37])+0.2500000000000001*(alphaR[3]*fUpR[30]+alphaR[2]*fUpR[29]+alphaR[1]*fUpR[28])+0.223606797749979*(alphaR[27]*fUpR[27]+alphaR[25]*fUpR[25]+alphaR[17]*fUpR[17])+0.25*alphaR[0]*fUpR[14]+0.223606797749979*(alphaR[10]*fUpR[10]+alphaR[8]*fUpR[8]+alphaR[4]*fUpR[4]); 
  GhatR[15] = (0.2*alphaR[37]+0.223606797749979*alphaR[10])*fUpR[46]+(0.2*alphaR[39]+0.223606797749979*alphaR[8])*fUpR[44]+0.223606797749979*(alphaR[17]*fUpR[40]+fUpR[18]*alphaR[39]+fUpR[16]*alphaR[37]+alphaR[17]*fUpR[35])+(0.2*(alphaR[22]+alphaR[21])+0.223606797749979*alphaR[3])*fUpR[34]+(0.2*(fUpR[22]+fUpR[21])+0.223606797749979*fUpR[3])*alphaR[34]+(0.2*(alphaR[24]+alphaR[19])+0.223606797749979*alphaR[2])*fUpR[33]+(0.2*(fUpR[24]+fUpR[19])+0.223606797749979*fUpR[2])*alphaR[33]+(0.2*(alphaR[23]+alphaR[20])+0.223606797749979*alphaR[1])*fUpR[32]+(0.2*(fUpR[23]+fUpR[20])+0.223606797749979*fUpR[1])*alphaR[32]+(0.223606797749979*(alphaR[27]+alphaR[25])+0.25*alphaR[4])*fUpR[31]+0.223606797749979*(alphaR[6]*fUpR[24]+fUpR[6]*alphaR[24]+alphaR[7]*fUpR[23]+fUpR[7]*alphaR[23]+alphaR[5]*fUpR[22]+fUpR[5]*alphaR[22]+alphaR[5]*fUpR[21]+fUpR[5]*alphaR[21]+alphaR[7]*fUpR[20]+fUpR[7]*alphaR[20]+alphaR[6]*fUpR[19]+fUpR[6]*alphaR[19])+0.25*(alphaR[8]*fUpR[18]+fUpR[9]*alphaR[17]+alphaR[10]*fUpR[16])+(0.223606797749979*(alphaR[13]+alphaR[12]+alphaR[11])+0.25*alphaR[0])*fUpR[15]+0.223606797749979*(fUpR[13]+fUpR[12]+fUpR[11])*alphaR[15]+0.25*(fUpR[0]*alphaR[15]+alphaR[1]*fUpR[7]+fUpR[1]*alphaR[7]+alphaR[2]*fUpR[6]+fUpR[2]*alphaR[6]+alphaR[3]*fUpR[5]+fUpR[3]*alphaR[5]); 
  GhatR[16] = (0.2*alphaR[37]+0.223606797749979*alphaR[10])*fUpR[47]+0.2500000000000001*alphaR[13]*fUpR[46]+(0.2*alphaR[32]+0.223606797749979*alphaR[7])*fUpR[45]+(0.2*alphaR[33]+0.223606797749979*alphaR[6])*fUpR[44]+0.223606797749979*alphaR[17]*fUpR[43]+(0.2*alphaR[25]+0.223606797749979*alphaR[4])*fUpR[41]+0.2500000000000001*(alphaR[23]*fUpR[40]+alphaR[24]*fUpR[39]+fUpR[24]*alphaR[39])+0.223606797749979*(alphaR[15]*(fUpR[38]+fUpR[37])+fUpR[15]*alphaR[37])+(0.2*alphaR[19]+0.223606797749979*alphaR[2])*fUpR[36]+(0.2*alphaR[20]+0.223606797749979*alphaR[1])*fUpR[35]+0.2500000000000001*(alphaR[27]*fUpR[34]+fUpR[27]*alphaR[34])+0.223606797749979*(fUpR[18]*alphaR[33]+alphaR[17]*fUpR[32]+fUpR[17]*alphaR[32])+(0.223606797749979*(alphaR[22]+alphaR[21])+0.25*alphaR[3])*fUpR[31]+0.223606797749979*(alphaR[8]*fUpR[29]+alphaR[5]*(fUpR[26]+fUpR[25])+fUpR[5]*alphaR[25]+fUpR[9]*alphaR[20]+alphaR[8]*fUpR[19]+fUpR[8]*alphaR[19])+0.25*(alphaR[6]*fUpR[18]+alphaR[7]*fUpR[17]+fUpR[7]*alphaR[17])+0.223606797749979*(alphaR[12]+alphaR[11])*fUpR[16]+0.25*(alphaR[0]*fUpR[16]+alphaR[10]*fUpR[15]+fUpR[10]*alphaR[15]+alphaR[1]*fUpR[9]+alphaR[2]*fUpR[8]+fUpR[2]*alphaR[8]+alphaR[4]*fUpR[5]+fUpR[4]*alphaR[5]); 
  GhatR[17] = (0.2*alphaR[32]+0.223606797749979*alphaR[7])*fUpR[46]+0.2500000000000001*alphaR[12]*fUpR[45]+(0.2*alphaR[34]+0.223606797749979*alphaR[5])*fUpR[44]+0.2*(alphaR[27]+alphaR[25])*fUpR[42]+0.223606797749979*(alphaR[4]*fUpR[42]+alphaR[15]*fUpR[40])+(0.2*alphaR[21]+0.223606797749979*alphaR[3])*fUpR[39]+(0.2*(fUpR[30]+fUpR[21])+0.223606797749979*fUpR[3])*alphaR[39]+0.2500000000000001*alphaR[20]*fUpR[38]+(0.2*alphaR[23]+0.223606797749979*alphaR[1])*fUpR[37]+(0.2*(fUpR[28]+fUpR[23])+0.223606797749979*fUpR[1])*alphaR[37]+0.2500000000000001*alphaR[22]*fUpR[36]+0.223606797749979*(alphaR[15]*fUpR[35]+fUpR[18]*alphaR[34])+0.2500000000000001*fUpR[26]*alphaR[33]+0.223606797749979*fUpR[16]*alphaR[32]+(0.223606797749979*(alphaR[24]+alphaR[19])+0.25*alphaR[2])*fUpR[31]+0.223606797749979*(alphaR[8]*fUpR[30]+alphaR[10]*fUpR[28]+alphaR[6]*fUpR[27]+fUpR[6]*alphaR[27]+alphaR[6]*fUpR[25]+fUpR[6]*alphaR[25]+alphaR[10]*fUpR[23]+fUpR[10]*alphaR[23]+alphaR[8]*fUpR[21]+fUpR[8]*alphaR[21])+0.25*alphaR[5]*fUpR[18]+(0.223606797749979*(alphaR[13]+alphaR[11])+0.25*alphaR[0])*fUpR[17]+0.223606797749979*(fUpR[14]+fUpR[13]+fUpR[11])*alphaR[17]+0.25*(fUpR[0]*alphaR[17]+alphaR[7]*fUpR[16]+fUpR[9]*alphaR[15]+alphaR[1]*fUpR[10]+fUpR[1]*alphaR[10]+alphaR[3]*fUpR[8]+fUpR[3]*alphaR[8]+alphaR[4]*fUpR[6]+fUpR[4]*alphaR[6]); 
  GhatR[18] = (0.2*alphaR[39]+0.223606797749979*alphaR[8])*fUpR[47]+(0.2*alphaR[33]+0.223606797749979*alphaR[6])*fUpR[46]+(0.2*alphaR[34]+0.223606797749979*alphaR[5])*fUpR[45]+0.2500000000000001*alphaR[11]*fUpR[44]+0.2*alphaR[27]*fUpR[43]+0.223606797749979*(alphaR[4]*fUpR[43]+alphaR[17]*fUpR[41])+0.2*alphaR[22]*fUpR[40]+0.223606797749979*(alphaR[3]*fUpR[40]+alphaR[15]*fUpR[39]+fUpR[15]*alphaR[39])+(0.2*alphaR[24]+0.223606797749979*alphaR[2])*fUpR[38]+0.2500000000000001*(alphaR[19]*fUpR[37]+fUpR[19]*alphaR[37])+0.223606797749979*alphaR[15]*fUpR[36]+0.2500000000000001*alphaR[21]*fUpR[35]+0.223606797749979*(alphaR[17]*fUpR[34]+fUpR[17]*alphaR[34]+fUpR[16]*alphaR[33])+0.2500000000000001*(alphaR[25]*fUpR[32]+fUpR[25]*alphaR[32])+(0.223606797749979*(alphaR[23]+alphaR[20])+0.25*alphaR[1])*fUpR[31]+0.223606797749979*(alphaR[10]*fUpR[29]+alphaR[7]*fUpR[27]+fUpR[7]*alphaR[27]+alphaR[7]*fUpR[26]+alphaR[10]*fUpR[24]+fUpR[10]*alphaR[24]+fUpR[9]*alphaR[22])+0.223606797749979*(alphaR[13]+alphaR[12])*fUpR[18]+0.25*(alphaR[0]*fUpR[18]+alphaR[5]*fUpR[17]+fUpR[5]*alphaR[17]+alphaR[6]*fUpR[16]+alphaR[8]*fUpR[15]+fUpR[8]*alphaR[15]+alphaR[2]*fUpR[10]+fUpR[2]*alphaR[10]+alphaR[3]*fUpR[9]+alphaR[4]*fUpR[7]+fUpR[4]*alphaR[7]); 
  GhatR[19] = 0.223606797749979*alphaR[39]*fUpR[46]+(0.159719141249985*alphaR[37]+0.25*alphaR[10])*fUpR[44]+0.2500000000000001*fUpR[18]*alphaR[37]+(0.159719141249985*alphaR[25]+0.2500000000000001*alphaR[4])*fUpR[35]+0.223606797749979*(alphaR[23]*fUpR[34]+fUpR[23]*alphaR[34])+0.2*(alphaR[15]*fUpR[33]+fUpR[15]*alphaR[33])+(0.223606797749979*alphaR[22]+0.159719141249985*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[32]+(0.223606797749979*fUpR[22]+0.159719141249985*fUpR[21]+0.2500000000000001*fUpR[3])*alphaR[32]+0.223606797749979*alphaR[17]*fUpR[31]+0.25*(fUpR[9]*alphaR[25]+alphaR[7]*fUpR[21]+fUpR[7]*alphaR[21])+0.2*(alphaR[5]*fUpR[20]+fUpR[5]*alphaR[20])+(0.223606797749979*alphaR[12]+0.159719141249985*alphaR[11]+0.25*alphaR[0])*fUpR[19]+(0.223606797749979*fUpR[12]+0.159719141249985*fUpR[11]+0.25*fUpR[0])*alphaR[19]+0.223606797749979*(alphaR[8]*fUpR[16]+alphaR[6]*fUpR[15]+fUpR[6]*alphaR[15])+0.2500000000000001*(alphaR[2]*fUpR[11]+fUpR[2]*alphaR[11])+0.223606797749979*(alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]); 
  GhatR[20] = (0.223606797749979*alphaR[37]+0.25*alphaR[10])*fUpR[45]+0.2500000000000001*alphaR[17]*fUpR[38]+(0.223606797749979*alphaR[25]+0.2500000000000001*alphaR[4])*fUpR[36]+0.223606797749979*(alphaR[24]*fUpR[34]+fUpR[24]*alphaR[34])+(0.159719141249985*alphaR[22]+0.223606797749979*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[33]+(0.159719141249985*fUpR[22]+0.223606797749979*fUpR[21]+0.2500000000000001*fUpR[3])*alphaR[33]+0.2*(alphaR[15]*fUpR[32]+fUpR[15]*alphaR[32])+0.25*(alphaR[8]*fUpR[26]+alphaR[6]*fUpR[22]+fUpR[6]*alphaR[22])+(0.159719141249985*alphaR[12]+0.223606797749979*alphaR[11]+0.25*alphaR[0])*fUpR[20]+(0.159719141249985*fUpR[12]+0.223606797749979*fUpR[11]+0.25*fUpR[0])*alphaR[20]+0.2*(alphaR[5]*fUpR[19]+fUpR[5]*alphaR[19])+0.223606797749979*(alphaR[7]*fUpR[15]+fUpR[7]*alphaR[15])+0.2500000000000001*(alphaR[1]*fUpR[12]+fUpR[1]*alphaR[12])+0.223606797749979*(alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5]); 
  GhatR[21] = 0.2*(alphaR[17]*fUpR[39]+fUpR[17]*alphaR[39])+(0.223606797749979*alphaR[27]+0.159719141249985*alphaR[25]+0.2500000000000001*alphaR[4])*fUpR[37]+(0.223606797749979*fUpR[27]+0.159719141249985*fUpR[25]+0.2500000000000001*fUpR[4])*alphaR[37]+0.2*(alphaR[15]*fUpR[34]+fUpR[15]*alphaR[34])+0.223606797749979*(alphaR[20]*fUpR[33]+fUpR[20]*alphaR[33])+(0.223606797749979*alphaR[24]+0.159719141249985*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[32]+(0.223606797749979*fUpR[24]+0.159719141249985*fUpR[19]+0.2500000000000001*fUpR[2])*alphaR[32]+0.25*(alphaR[10]*fUpR[25]+fUpR[10]*alphaR[25])+0.2*(alphaR[6]*fUpR[23]+fUpR[6]*alphaR[23])+(0.223606797749979*alphaR[13]+0.159719141249985*alphaR[11]+0.25*alphaR[0])*fUpR[21]+(0.223606797749979*fUpR[13]+0.159719141249985*fUpR[11])*alphaR[21]+0.25*(fUpR[0]*alphaR[21]+alphaR[7]*fUpR[19]+fUpR[7]*alphaR[19])+0.223606797749979*(alphaR[8]*fUpR[17]+fUpR[8]*alphaR[17]+alphaR[5]*fUpR[15]+fUpR[5]*alphaR[15])+0.2500000000000001*(alphaR[3]*fUpR[11]+fUpR[3]*alphaR[11])+0.223606797749979*(alphaR[1]*fUpR[6]+fUpR[1]*alphaR[6]); 
  GhatR[22] = (0.223606797749979*alphaR[39]+0.25*alphaR[8])*fUpR[45]+0.223606797749979*alphaR[27]*fUpR[38]+0.2500000000000001*(alphaR[4]*fUpR[38]+alphaR[17]*fUpR[36])+0.2*(alphaR[15]*fUpR[34]+fUpR[15]*alphaR[34])+(0.223606797749979*alphaR[23]+0.159719141249985*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[33]+(0.223606797749979*fUpR[23]+0.159719141249985*fUpR[20]+0.2500000000000001*fUpR[1])*alphaR[33]+0.223606797749979*(alphaR[19]*fUpR[32]+fUpR[19]*alphaR[32])+0.25*alphaR[10]*fUpR[26]+0.2*(alphaR[7]*fUpR[24]+fUpR[7]*alphaR[24])+(0.223606797749979*alphaR[13]+0.159719141249985*alphaR[12]+0.25*alphaR[0])*fUpR[22]+(0.223606797749979*fUpR[13]+0.159719141249985*fUpR[12])*alphaR[22]+0.25*(fUpR[0]*alphaR[22]+alphaR[6]*fUpR[20]+fUpR[6]*alphaR[20])+0.223606797749979*(alphaR[5]*fUpR[15]+fUpR[5]*alphaR[15])+0.2500000000000001*(alphaR[3]*fUpR[12]+fUpR[3]*alphaR[12])+0.223606797749979*(alphaR[2]*fUpR[7]+fUpR[2]*alphaR[7]); 
  GhatR[23] = (0.159719141249985*alphaR[27]+0.223606797749979*alphaR[25]+0.2500000000000001*alphaR[4])*fUpR[39]+(0.159719141249985*fUpR[27]+0.223606797749979*fUpR[25]+0.2500000000000001*fUpR[4])*alphaR[39]+0.2*(alphaR[17]*fUpR[37]+fUpR[17]*alphaR[37])+(0.159719141249985*alphaR[24]+0.223606797749979*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[34]+(0.159719141249985*fUpR[24]+0.223606797749979*fUpR[19]+0.2500000000000001*fUpR[2])*alphaR[34]+0.223606797749979*(alphaR[22]*fUpR[33]+fUpR[22]*alphaR[33])+0.2*(alphaR[15]*fUpR[32]+fUpR[15]*alphaR[32])+0.25*(alphaR[8]*fUpR[27]+fUpR[8]*alphaR[27]+alphaR[5]*fUpR[24]+fUpR[5]*alphaR[24])+(0.159719141249985*alphaR[13]+0.223606797749979*alphaR[11]+0.25*alphaR[0])*fUpR[23]+(0.159719141249985*fUpR[13]+0.223606797749979*fUpR[11]+0.25*fUpR[0])*alphaR[23]+0.2*(alphaR[6]*fUpR[21]+fUpR[6]*alphaR[21])+0.223606797749979*(alphaR[10]*fUpR[17]+fUpR[10]*alphaR[17]+alphaR[7]*fUpR[15]+fUpR[7]*alphaR[15])+0.2500000000000001*(alphaR[1]*fUpR[13]+fUpR[1]*alphaR[13])+0.223606797749979*(alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6]); 
  GhatR[24] = (0.159719141249985*alphaR[39]+0.25*alphaR[8])*fUpR[46]+0.223606797749979*alphaR[37]*fUpR[44]+0.159719141249985*alphaR[27]*fUpR[40]+0.2500000000000001*(alphaR[4]*fUpR[40]+fUpR[16]*alphaR[39])+(0.159719141249985*alphaR[23]+0.223606797749979*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[34]+(0.159719141249985*fUpR[23]+0.223606797749979*fUpR[20]+0.2500000000000001*fUpR[1])*alphaR[34]+0.2*(alphaR[15]*fUpR[33]+fUpR[15]*alphaR[33])+0.223606797749979*(alphaR[21]*fUpR[32]+fUpR[21]*alphaR[32])+0.223606797749979*alphaR[17]*fUpR[31]+0.25*fUpR[9]*alphaR[27]+(0.159719141249985*alphaR[13]+0.223606797749979*alphaR[12]+0.25*alphaR[0])*fUpR[24]+(0.159719141249985*fUpR[13]+0.223606797749979*fUpR[12])*alphaR[24]+0.25*(fUpR[0]*alphaR[24]+alphaR[5]*fUpR[23]+fUpR[5]*alphaR[23])+0.2*(alphaR[7]*fUpR[22]+fUpR[7]*alphaR[22])+0.223606797749979*(alphaR[10]*fUpR[18]+alphaR[6]*fUpR[15]+fUpR[6]*alphaR[15])+0.2500000000000001*(alphaR[2]*fUpR[13]+fUpR[2]*alphaR[13])+0.223606797749979*(alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]); 
  GhatR[25] = 0.223606797749979*(alphaR[34]*fUpR[46]+alphaR[33]*fUpR[45])+(0.159719141249985*alphaR[32]+0.25*alphaR[7])*fUpR[44]+0.2*alphaR[17]*fUpR[42]+0.223606797749979*(alphaR[23]*fUpR[39]+fUpR[23]*alphaR[39])+(0.159719141249985*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[37]+(0.223606797749979*fUpR[30]+0.159719141249985*fUpR[21]+0.2500000000000001*fUpR[3])*alphaR[37]+0.223606797749979*alphaR[20]*fUpR[36]+0.159719141249985*alphaR[19]*fUpR[35]+0.2500000000000001*(alphaR[2]*fUpR[35]+fUpR[18]*alphaR[32])+0.223606797749979*alphaR[15]*fUpR[31]+0.2*alphaR[8]*fUpR[28]+(0.159719141249985*alphaR[11]+0.25*alphaR[0])*fUpR[25]+(0.223606797749979*fUpR[14]+0.159719141249985*fUpR[11])*alphaR[25]+0.25*(fUpR[0]*alphaR[25]+alphaR[10]*fUpR[21]+fUpR[10]*alphaR[21]+fUpR[9]*alphaR[19])+0.223606797749979*(alphaR[6]*fUpR[17]+fUpR[6]*alphaR[17]+alphaR[5]*fUpR[16])+0.2500000000000001*(alphaR[4]*fUpR[11]+fUpR[4]*alphaR[11])+0.223606797749979*(alphaR[1]*fUpR[8]+fUpR[1]*alphaR[8]); 
  GhatR[26] = 0.223606797749979*alphaR[34]*fUpR[46]+(0.159719141249985*alphaR[33]+0.25*alphaR[6])*fUpR[45]+0.223606797749979*(alphaR[32]*fUpR[44]+alphaR[24]*fUpR[40])+(0.159719141249985*alphaR[22]+0.2500000000000001*alphaR[3])*fUpR[38]+(0.159719141249985*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[36]+0.223606797749979*alphaR[19]*fUpR[35]+0.2500000000000001*(alphaR[17]*fUpR[33]+fUpR[17]*alphaR[33])+0.223606797749979*alphaR[15]*fUpR[31]+0.159719141249985*alphaR[12]*fUpR[26]+0.25*(alphaR[0]*fUpR[26]+alphaR[10]*fUpR[22]+fUpR[10]*alphaR[22]+alphaR[8]*fUpR[20]+fUpR[8]*alphaR[20])+0.223606797749979*(alphaR[7]*fUpR[18]+alphaR[5]*fUpR[16])+0.2500000000000001*(alphaR[4]*fUpR[12]+fUpR[4]*alphaR[12])+0.223606797749979*alphaR[2]*fUpR[9]; 
  GhatR[27] = (0.159719141249985*alphaR[34]+0.25*alphaR[5])*fUpR[46]+0.223606797749979*(alphaR[33]*fUpR[45]+alphaR[32]*fUpR[44])+0.2*alphaR[17]*fUpR[42]+(0.159719141249985*alphaR[24]+0.2500000000000001*alphaR[2])*fUpR[40]+(0.159719141249985*alphaR[23]+0.2500000000000001*alphaR[1])*fUpR[39]+(0.223606797749979*fUpR[28]+0.159719141249985*fUpR[23]+0.2500000000000001*fUpR[1])*alphaR[39]+0.223606797749979*(alphaR[22]*fUpR[38]+alphaR[21]*fUpR[37]+fUpR[21]*alphaR[37])+0.2500000000000001*fUpR[16]*alphaR[34]+0.223606797749979*alphaR[15]*fUpR[31]+0.2*alphaR[10]*fUpR[30]+(0.159719141249985*alphaR[13]+0.25*alphaR[0])*fUpR[27]+(0.223606797749979*fUpR[14]+0.159719141249985*fUpR[13])*alphaR[27]+0.25*(fUpR[0]*alphaR[27]+fUpR[9]*alphaR[24]+alphaR[8]*fUpR[23]+fUpR[8]*alphaR[23])+0.223606797749979*(alphaR[7]*fUpR[18]+alphaR[6]*fUpR[17]+fUpR[6]*alphaR[17])+0.2500000000000001*(alphaR[4]*fUpR[13]+fUpR[4]*alphaR[13])+0.223606797749979*(alphaR[3]*fUpR[10]+fUpR[3]*alphaR[10]); 
  GhatR[28] = (0.223606797749979*alphaR[32]+0.25*alphaR[7])*fUpR[47]+0.2500000000000001*alphaR[15]*fUpR[43]+(0.223606797749979*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[42]+(0.223606797749979*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[41]+0.223606797749979*(alphaR[27]*fUpR[39]+fUpR[27]*alphaR[39])+0.2*(alphaR[17]*fUpR[37]+fUpR[17]*alphaR[37])+0.25*(alphaR[6]*fUpR[30]+alphaR[5]*fUpR[29])+(0.223606797749979*alphaR[11]+0.25*alphaR[0])*fUpR[28]+0.2*(alphaR[8]*fUpR[25]+fUpR[8]*alphaR[25])+0.223606797749979*(alphaR[10]*fUpR[17]+fUpR[10]*alphaR[17])+0.2500000000000001*alphaR[1]*fUpR[14]+0.223606797749979*(alphaR[4]*fUpR[8]+fUpR[4]*alphaR[8]); 
  GhatR[29] = (0.223606797749979*alphaR[33]+0.25*alphaR[6])*fUpR[47]+0.223606797749979*(alphaR[39]*fUpR[46]+alphaR[37]*fUpR[44]+alphaR[22]*fUpR[43])+0.2500000000000001*(alphaR[3]*fUpR[43]+alphaR[15]*fUpR[42])+(0.223606797749979*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[41]+0.223606797749979*(alphaR[27]*fUpR[40]+alphaR[25]*fUpR[35])+0.223606797749979*alphaR[17]*fUpR[31]+0.25*alphaR[7]*fUpR[30]+0.223606797749979*alphaR[12]*fUpR[29]+0.25*(alphaR[0]*fUpR[29]+alphaR[5]*fUpR[28])+0.223606797749979*(alphaR[10]*fUpR[18]+alphaR[8]*fUpR[16])+0.2500000000000001*alphaR[2]*fUpR[14]+0.223606797749979*alphaR[4]*fUpR[9]; 
  GhatR[30] = (0.223606797749979*alphaR[34]+0.25*alphaR[5])*fUpR[47]+(0.223606797749979*alphaR[24]+0.2500000000000001*alphaR[2])*fUpR[43]+0.223606797749979*alphaR[23]*fUpR[42]+0.2500000000000001*(alphaR[1]*fUpR[42]+alphaR[15]*fUpR[41])+0.2*(alphaR[17]*fUpR[39]+fUpR[17]*alphaR[39])+0.223606797749979*(alphaR[25]*fUpR[37]+fUpR[25]*alphaR[37]+alphaR[13]*fUpR[30])+0.25*(alphaR[0]*fUpR[30]+alphaR[7]*fUpR[29]+alphaR[6]*fUpR[28])+0.2*(alphaR[10]*fUpR[27]+fUpR[10]*alphaR[27])+0.223606797749979*(alphaR[8]*fUpR[17]+fUpR[8]*alphaR[17])+0.2500000000000001*alphaR[3]*fUpR[14]+0.223606797749979*(alphaR[4]*fUpR[10]+fUpR[4]*alphaR[10]); 
  GhatR[31] = (0.2*(alphaR[27]+alphaR[25])+0.223606797749979*alphaR[4])*fUpR[47]+(0.2*(alphaR[22]+alphaR[21])+0.223606797749979*alphaR[3])*fUpR[46]+(0.2*(alphaR[24]+alphaR[19])+0.223606797749979*alphaR[2])*fUpR[45]+(0.2*(alphaR[23]+alphaR[20])+0.223606797749979*alphaR[1])*fUpR[44]+(0.2*alphaR[39]+0.223606797749979*alphaR[8])*fUpR[43]+(0.2*alphaR[37]+0.223606797749979*alphaR[10])*fUpR[41]+(0.2*alphaR[33]+0.223606797749979*alphaR[6])*fUpR[40]+(0.2*alphaR[32]+0.223606797749979*alphaR[7])*fUpR[39]+(0.2*fUpR[32]+0.223606797749979*fUpR[7])*alphaR[39]+(0.2*alphaR[34]+0.223606797749979*alphaR[5])*(fUpR[38]+fUpR[37])+(0.2*fUpR[34]+0.223606797749979*fUpR[5])*alphaR[37]+(0.2*alphaR[32]+0.223606797749979*alphaR[7])*fUpR[36]+0.2*alphaR[33]*fUpR[35]+0.223606797749979*(alphaR[6]*fUpR[35]+alphaR[10]*fUpR[34]+fUpR[10]*alphaR[34]+fUpR[9]*alphaR[33]+alphaR[8]*fUpR[32]+fUpR[8]*alphaR[32])+(0.223606797749979*(alphaR[13]+alphaR[12]+alphaR[11])+0.25*alphaR[0])*fUpR[31]+0.223606797749979*(alphaR[17]*fUpR[29]+alphaR[15]*fUpR[27]+fUpR[15]*alphaR[27]+alphaR[15]*(fUpR[26]+fUpR[25])+fUpR[15]*alphaR[25]+alphaR[17]*fUpR[24]+fUpR[17]*alphaR[24]+fUpR[18]*alphaR[23]+fUpR[16]*(alphaR[22]+alphaR[21])+fUpR[18]*alphaR[20]+alphaR[17]*fUpR[19]+fUpR[17]*alphaR[19])+0.25*(alphaR[1]*fUpR[18]+alphaR[2]*fUpR[17]+fUpR[2]*alphaR[17]+alphaR[3]*fUpR[16]+alphaR[4]*fUpR[15]+fUpR[4]*alphaR[15]+alphaR[5]*fUpR[10]+fUpR[5]*alphaR[10]+alphaR[6]*fUpR[9]+alphaR[7]*fUpR[8]+fUpR[7]*alphaR[8]); 
  GhatR[32] = 0.2*alphaR[17]*fUpR[46]+(0.223606797749979*alphaR[27]+0.159719141249985*alphaR[25]+0.2500000000000001*alphaR[4])*fUpR[44]+0.223606797749979*alphaR[37]*fUpR[40]+0.2*fUpR[31]*alphaR[39]+0.159719141249985*fUpR[35]*alphaR[37]+0.25*(fUpR[9]*alphaR[37]+alphaR[10]*fUpR[35])+(0.1788854381999831*alphaR[33]+0.2*alphaR[6])*fUpR[34]+0.1788854381999831*fUpR[33]*alphaR[34]+0.2*(fUpR[6]*alphaR[34]+alphaR[5]*fUpR[33]+fUpR[5]*alphaR[33])+(0.223606797749979*(alphaR[13]+alphaR[12])+0.159719141249985*alphaR[11]+0.25*alphaR[0])*fUpR[32]+(0.223606797749979*(fUpR[13]+fUpR[12])+0.159719141249985*fUpR[11]+0.25*fUpR[0])*alphaR[32]+0.223606797749979*alphaR[8]*fUpR[31]+0.2500000000000001*fUpR[18]*alphaR[25]+0.223606797749979*(alphaR[21]*fUpR[24]+fUpR[21]*alphaR[24])+0.2*(alphaR[15]*fUpR[23]+fUpR[15]*alphaR[23])+0.223606797749979*(alphaR[19]*fUpR[22]+fUpR[19]*alphaR[22])+(0.159719141249985*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[21]+(0.159719141249985*fUpR[19]+0.2500000000000001*fUpR[2])*alphaR[21]+0.2*(alphaR[15]*fUpR[20]+fUpR[15]*alphaR[20])+0.2500000000000001*(alphaR[3]*fUpR[19]+fUpR[3]*alphaR[19])+0.223606797749979*(fUpR[16]*alphaR[17]+alphaR[1]*fUpR[15]+fUpR[1]*alphaR[15])+0.25*(alphaR[7]*fUpR[11]+fUpR[7]*alphaR[11])+0.223606797749979*(alphaR[5]*fUpR[6]+fUpR[5]*alphaR[6]); 
  GhatR[33] = (0.223606797749979*(alphaR[27]+alphaR[25])+0.2500000000000001*alphaR[4])*fUpR[45]+fUpR[38]*(0.223606797749979*alphaR[39]+0.25*alphaR[8])+fUpR[36]*(0.223606797749979*alphaR[37]+0.25*alphaR[10])+(0.1788854381999831*alphaR[32]+0.2*alphaR[7])*fUpR[34]+(0.1788854381999831*fUpR[32]+0.2*fUpR[7])*alphaR[34]+(0.223606797749979*alphaR[13]+0.159719141249985*alphaR[12]+0.223606797749979*alphaR[11]+0.25*alphaR[0])*fUpR[33]+(0.223606797749979*fUpR[13]+0.159719141249985*fUpR[12]+0.223606797749979*fUpR[11]+0.25*fUpR[0])*alphaR[33]+0.2*(alphaR[5]*fUpR[32]+fUpR[5]*alphaR[32])+0.2500000000000001*alphaR[17]*fUpR[26]+0.2*(alphaR[15]*fUpR[24]+fUpR[15]*alphaR[24])+0.223606797749979*alphaR[22]*fUpR[23]+fUpR[22]*(0.223606797749979*alphaR[23]+0.159719141249985*alphaR[20]+0.2500000000000001*alphaR[1])+(0.159719141249985*fUpR[20]+0.2500000000000001*fUpR[1])*alphaR[22]+0.223606797749979*(alphaR[20]*fUpR[21]+fUpR[20]*alphaR[21])+0.2500000000000001*(alphaR[3]*fUpR[20]+fUpR[3]*alphaR[20])+0.2*(alphaR[15]*fUpR[19]+fUpR[15]*alphaR[19])+0.223606797749979*(alphaR[2]*fUpR[15]+fUpR[2]*alphaR[15])+0.25*(alphaR[6]*fUpR[12]+fUpR[6]*alphaR[12])+0.223606797749979*(alphaR[5]*fUpR[7]+fUpR[5]*alphaR[7]); 
  GhatR[34] = (0.159719141249985*alphaR[27]+0.223606797749979*alphaR[25]+0.2500000000000001*alphaR[4])*fUpR[46]+0.2*alphaR[17]*fUpR[44]+(0.159719141249985*alphaR[39]+0.25*alphaR[8])*fUpR[40]+(0.223606797749979*fUpR[35]+0.25*fUpR[9])*alphaR[39]+0.2*fUpR[31]*alphaR[37]+(0.159719141249985*alphaR[13]+0.223606797749979*(alphaR[12]+alphaR[11])+0.25*alphaR[0])*fUpR[34]+(0.159719141249985*fUpR[13]+0.223606797749979*(fUpR[12]+fUpR[11])+0.25*fUpR[0])*alphaR[34]+(0.1788854381999831*alphaR[32]+0.2*alphaR[7])*fUpR[33]+0.1788854381999831*fUpR[32]*alphaR[33]+0.2*(fUpR[7]*alphaR[33]+alphaR[6]*fUpR[32]+fUpR[6]*alphaR[32])+0.223606797749979*alphaR[10]*fUpR[31]+0.2500000000000001*fUpR[16]*alphaR[27]+(0.159719141249985*alphaR[23]+0.223606797749979*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[24]+(0.159719141249985*fUpR[23]+0.223606797749979*fUpR[20]+0.2500000000000001*fUpR[1])*alphaR[24]+(0.223606797749979*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[23]+(0.223606797749979*fUpR[19]+0.2500000000000001*fUpR[2])*alphaR[23]+0.2*(alphaR[15]*fUpR[22]+fUpR[15]*alphaR[22]+alphaR[15]*fUpR[21]+fUpR[15]*alphaR[21])+0.223606797749979*(alphaR[17]*fUpR[18]+alphaR[3]*fUpR[15]+fUpR[3]*alphaR[15])+0.25*(alphaR[5]*fUpR[13]+fUpR[5]*alphaR[13])+0.223606797749979*(alphaR[6]*fUpR[7]+fUpR[6]*alphaR[7]); 
  GhatR[35] = 0.2*alphaR[17]*fUpR[47]+0.223606797749979*alphaR[23]*fUpR[46]+0.2*alphaR[15]*fUpR[45]+(0.223606797749979*alphaR[22]+0.159719141249985*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[44]+0.223606797749979*alphaR[37]*fUpR[43]+0.2*alphaR[8]*fUpR[41]+0.223606797749979*(alphaR[34]*fUpR[39]+fUpR[34]*alphaR[39]+alphaR[32]*fUpR[38])+(0.159719141249985*alphaR[32]+0.25*alphaR[7])*fUpR[37]+(0.159719141249985*fUpR[32]+0.25*fUpR[7])*alphaR[37]+0.2*alphaR[5]*fUpR[36]+(0.223606797749979*alphaR[12]+0.159719141249985*alphaR[11]+0.25*alphaR[0])*fUpR[35]+0.2*fUpR[31]*alphaR[33]+0.25*(alphaR[10]*fUpR[32]+fUpR[10]*alphaR[32])+0.223606797749979*(alphaR[6]*fUpR[31]+alphaR[25]*fUpR[29]+alphaR[19]*fUpR[26])+(0.159719141249985*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[25]+0.159719141249985*fUpR[19]*alphaR[25]+0.2500000000000001*(fUpR[2]*alphaR[25]+fUpR[18]*alphaR[21])+0.2*fUpR[16]*alphaR[20]+0.2500000000000001*(alphaR[4]*fUpR[19]+fUpR[4]*alphaR[19])+0.223606797749979*(alphaR[15]*fUpR[17]+fUpR[15]*alphaR[17]+alphaR[1]*fUpR[16])+0.25*fUpR[9]*alphaR[11]+0.223606797749979*(alphaR[5]*fUpR[8]+fUpR[5]*alphaR[8]); 
  GhatR[36] = 0.223606797749979*alphaR[24]*fUpR[46]+(0.159719141249985*alphaR[22]+0.223606797749979*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[45]+0.2*alphaR[15]*fUpR[44]+0.223606797749979*alphaR[34]*fUpR[40]+(0.159719141249985*alphaR[33]+0.25*alphaR[6])*fUpR[38]+0.223606797749979*(alphaR[33]*fUpR[37]+fUpR[33]*alphaR[37])+(0.159719141249985*alphaR[12]+0.223606797749979*alphaR[11]+0.25*alphaR[0])*fUpR[36]+0.2*alphaR[5]*fUpR[35]+0.25*(alphaR[10]*fUpR[33]+fUpR[10]*alphaR[33])+fUpR[31]*(0.2*alphaR[32]+0.223606797749979*alphaR[7])+(0.159719141249985*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[26]+0.223606797749979*(alphaR[20]*fUpR[25]+fUpR[20]*alphaR[25])+0.2500000000000001*(alphaR[17]*fUpR[22]+fUpR[17]*alphaR[22]+alphaR[4]*fUpR[20]+fUpR[4]*alphaR[20])+0.2*fUpR[16]*alphaR[19]+0.223606797749979*(alphaR[15]*fUpR[18]+alphaR[2]*fUpR[16])+0.25*(alphaR[8]*fUpR[12]+fUpR[8]*alphaR[12])+0.223606797749979*alphaR[5]*fUpR[9]; 
  GhatR[37] = 0.2*alphaR[15]*fUpR[46]+0.223606797749979*alphaR[20]*fUpR[45]+(0.223606797749979*alphaR[24]+0.159719141249985*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[44]+(0.1788854381999831*alphaR[39]+0.2*alphaR[8])*fUpR[42]+0.223606797749979*alphaR[32]*fUpR[40]+0.2*(alphaR[6]*fUpR[39]+fUpR[6]*alphaR[39])+(0.223606797749979*alphaR[13]+0.159719141249985*alphaR[11]+0.25*alphaR[0])*fUpR[37]+(0.223606797749979*(fUpR[14]+fUpR[13])+0.159719141249985*fUpR[11]+0.25*fUpR[0])*alphaR[37]+0.223606797749979*alphaR[33]*fUpR[36]+(0.159719141249985*alphaR[32]+0.25*alphaR[7])*fUpR[35]+0.2*fUpR[31]*alphaR[34]+0.25*fUpR[9]*alphaR[32]+0.223606797749979*(alphaR[5]*fUpR[31]+alphaR[25]*fUpR[30])+0.2*alphaR[17]*fUpR[28]+0.223606797749979*(alphaR[21]*fUpR[27]+fUpR[21]*alphaR[27])+(0.159719141249985*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[25]+(0.159719141249985*fUpR[21]+0.2500000000000001*fUpR[3])*alphaR[25]+0.2*(alphaR[17]*fUpR[23]+fUpR[17]*alphaR[23])+0.2500000000000001*(alphaR[4]*fUpR[21]+fUpR[4]*alphaR[21]+fUpR[18]*alphaR[19])+0.223606797749979*(alphaR[1]*fUpR[17]+fUpR[1]*alphaR[17]+alphaR[15]*fUpR[16])+0.25*(alphaR[10]*fUpR[11]+fUpR[10]*alphaR[11])+0.223606797749979*(alphaR[6]*fUpR[8]+fUpR[6]*alphaR[8]); 
  GhatR[38] = 0.2*alphaR[15]*fUpR[46]+(0.223606797749979*alphaR[23]+0.159719141249985*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[45]+0.223606797749979*alphaR[19]*fUpR[44]+0.2*alphaR[7]*fUpR[40]+0.223606797749979*(alphaR[33]*fUpR[39]+fUpR[33]*alphaR[39])+(0.223606797749979*alphaR[13]+0.159719141249985*alphaR[12]+0.25*alphaR[0])*fUpR[38]+(0.159719141249985*alphaR[33]+0.25*alphaR[6])*fUpR[36]+0.223606797749979*alphaR[32]*fUpR[35]+0.2*fUpR[31]*alphaR[34]+0.25*(alphaR[8]*fUpR[33]+fUpR[8]*alphaR[33])+0.223606797749979*(alphaR[5]*fUpR[31]+alphaR[22]*fUpR[27]+fUpR[22]*alphaR[27])+(0.159719141249985*alphaR[22]+0.2500000000000001*alphaR[3])*fUpR[26]+0.2*fUpR[18]*alphaR[24]+0.2500000000000001*(alphaR[4]*fUpR[22]+fUpR[4]*alphaR[22]+alphaR[17]*fUpR[20]+fUpR[17]*alphaR[20])+0.223606797749979*(alphaR[2]*fUpR[18]+alphaR[15]*fUpR[16])+0.25*(alphaR[10]*fUpR[12]+fUpR[10]*alphaR[12])+0.223606797749979*alphaR[7]*fUpR[9]; 
  GhatR[39] = (0.159719141249985*alphaR[24]+0.223606797749979*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[46]+0.223606797749979*alphaR[22]*fUpR[45]+0.2*alphaR[15]*fUpR[44]+(0.1788854381999831*alphaR[37]+0.2*alphaR[10])*fUpR[42]+(0.159719141249985*alphaR[34]+0.25*alphaR[5])*fUpR[40]+(0.159719141249985*alphaR[13]+0.223606797749979*alphaR[11]+0.25*alphaR[0])*fUpR[39]+(0.223606797749979*fUpR[14]+0.159719141249985*fUpR[13]+0.223606797749979*fUpR[11]+0.25*fUpR[0])*alphaR[39]+0.223606797749979*alphaR[33]*fUpR[38]+0.2*(alphaR[6]*fUpR[37]+fUpR[6]*alphaR[37])+alphaR[34]*(0.223606797749979*fUpR[35]+0.25*fUpR[9])+fUpR[31]*(0.2*alphaR[32]+0.223606797749979*alphaR[7])+0.2*alphaR[17]*fUpR[30]+0.223606797749979*alphaR[27]*fUpR[28]+(0.159719141249985*alphaR[23]+0.2500000000000001*alphaR[1])*fUpR[27]+(0.159719141249985*fUpR[23]+0.2500000000000001*fUpR[1])*alphaR[27]+0.223606797749979*(alphaR[23]*fUpR[25]+fUpR[23]*alphaR[25])+0.2500000000000001*(fUpR[16]*alphaR[24]+alphaR[4]*fUpR[23]+fUpR[4]*alphaR[23])+0.2*(alphaR[17]*fUpR[21]+fUpR[17]*alphaR[21])+0.223606797749979*(alphaR[15]*fUpR[18]+alphaR[3]*fUpR[17]+fUpR[3]*alphaR[17])+0.25*(alphaR[8]*fUpR[13]+fUpR[8]*alphaR[13])+0.223606797749979*(alphaR[6]*fUpR[10]+fUpR[6]*alphaR[10]); 
  GhatR[40] = 0.2*alphaR[17]*fUpR[47]+(0.159719141249985*alphaR[23]+0.223606797749979*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[46]+0.2*alphaR[15]*fUpR[45]+0.223606797749979*alphaR[21]*fUpR[44]+0.2*alphaR[10]*fUpR[43]+0.223606797749979*alphaR[39]*fUpR[41]+(0.159719141249985*alphaR[13]+0.223606797749979*alphaR[12]+0.25*alphaR[0])*fUpR[40]+(0.159719141249985*alphaR[34]+0.25*alphaR[5])*fUpR[39]+(0.159719141249985*fUpR[34]+0.25*fUpR[5])*alphaR[39]+0.2*alphaR[7]*fUpR[38]+0.223606797749979*(alphaR[32]*fUpR[37]+fUpR[32]*alphaR[37]+alphaR[34]*fUpR[36])+0.25*(alphaR[8]*fUpR[34]+fUpR[8]*alphaR[34])+0.2*fUpR[31]*alphaR[33]+0.223606797749979*(alphaR[6]*fUpR[31]+alphaR[27]*fUpR[29])+(0.159719141249985*alphaR[24]+0.2500000000000001*alphaR[2])*fUpR[27]+(0.159719141249985*fUpR[24]+0.2500000000000001*fUpR[2])*alphaR[27]+0.223606797749979*alphaR[24]*fUpR[26]+0.2500000000000001*(alphaR[4]*fUpR[24]+fUpR[4]*alphaR[24]+fUpR[16]*alphaR[23])+0.2*fUpR[18]*alphaR[22]+0.223606797749979*(alphaR[3]*fUpR[18]+alphaR[15]*fUpR[17]+fUpR[15]*alphaR[17])+0.25*fUpR[9]*alphaR[13]+0.223606797749979*(alphaR[7]*fUpR[10]+fUpR[7]*alphaR[10]); 
  GhatR[41] = (0.223606797749979*(alphaR[22]+alphaR[21])+0.2500000000000001*alphaR[3])*fUpR[47]+0.223606797749979*alphaR[27]*fUpR[46]+0.2*alphaR[17]*fUpR[44]+(0.223606797749979*alphaR[33]+0.25*alphaR[6])*fUpR[43]+(0.223606797749979*alphaR[32]+0.25*alphaR[7])*fUpR[42]+(0.223606797749979*(alphaR[12]+alphaR[11])+0.25*alphaR[0])*fUpR[41]+0.223606797749979*alphaR[39]*fUpR[40]+0.2*(fUpR[31]*alphaR[37]+alphaR[8]*fUpR[35])+0.223606797749979*alphaR[10]*fUpR[31]+0.2500000000000001*alphaR[15]*fUpR[30]+(0.223606797749979*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[29]+(0.223606797749979*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[28]+0.2*fUpR[16]*alphaR[25]+0.223606797749979*(alphaR[17]*fUpR[18]+alphaR[4]*fUpR[16])+0.25*alphaR[5]*fUpR[14]+0.223606797749979*alphaR[8]*fUpR[9]; 
  GhatR[42] = (0.223606797749979*(alphaR[24]+alphaR[19])+0.2500000000000001*alphaR[2])*fUpR[47]+(0.223606797749979*alphaR[34]+0.25*alphaR[5])*fUpR[43]+(0.223606797749979*(alphaR[13]+alphaR[11])+0.25*alphaR[0])*fUpR[42]+(0.223606797749979*alphaR[32]+0.25*alphaR[7])*fUpR[41]+(0.1788854381999831*alphaR[37]+0.2*alphaR[10])*fUpR[39]+0.1788854381999831*fUpR[37]*alphaR[39]+0.2*(fUpR[10]*alphaR[39]+alphaR[8]*fUpR[37]+fUpR[8]*alphaR[37])+0.223606797749979*alphaR[23]*fUpR[30]+0.2500000000000001*(alphaR[1]*fUpR[30]+alphaR[15]*fUpR[29])+(0.223606797749979*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[28]+0.2*(alphaR[17]*fUpR[27]+fUpR[17]*alphaR[27]+alphaR[17]*fUpR[25]+fUpR[17]*alphaR[25])+0.223606797749979*(alphaR[4]*fUpR[17]+fUpR[4]*alphaR[17])+0.25*alphaR[6]*fUpR[14]+0.223606797749979*(alphaR[8]*fUpR[10]+fUpR[8]*alphaR[10]); 
  GhatR[43] = (0.223606797749979*(alphaR[23]+alphaR[20])+0.2500000000000001*alphaR[1])*fUpR[47]+0.2*alphaR[17]*fUpR[46]+0.223606797749979*alphaR[25]*fUpR[44]+(0.223606797749979*(alphaR[13]+alphaR[12])+0.25*alphaR[0])*fUpR[43]+(0.223606797749979*alphaR[34]+0.25*alphaR[5])*fUpR[42]+(0.223606797749979*alphaR[33]+0.25*alphaR[6])*fUpR[41]+0.2*(alphaR[10]*fUpR[40]+fUpR[31]*alphaR[39])+0.223606797749979*(fUpR[35]*alphaR[37]+alphaR[8]*fUpR[31])+(0.223606797749979*alphaR[24]+0.2500000000000001*alphaR[2])*fUpR[30]+0.223606797749979*alphaR[22]*fUpR[29]+0.2500000000000001*(alphaR[3]*fUpR[29]+alphaR[15]*fUpR[28])+0.2*fUpR[18]*alphaR[27]+0.223606797749979*(alphaR[4]*fUpR[18]+fUpR[16]*alphaR[17])+0.25*alphaR[7]*fUpR[14]+0.223606797749979*fUpR[9]*alphaR[10]; 
  GhatR[44] = (0.1788854381999831*alphaR[39]+0.2*alphaR[8])*fUpR[47]+(0.1788854381999831*alphaR[33]+0.2*alphaR[6])*fUpR[46]+(0.1788854381999831*alphaR[34]+0.2*alphaR[5])*fUpR[45]+(0.223606797749979*(alphaR[13]+alphaR[12])+0.159719141249985*alphaR[11]+0.25*alphaR[0])*fUpR[44]+0.223606797749979*alphaR[25]*fUpR[43]+0.2*alphaR[17]*fUpR[41]+0.223606797749979*alphaR[21]*fUpR[40]+0.2*(alphaR[15]*fUpR[39]+fUpR[15]*alphaR[39])+0.223606797749979*alphaR[19]*fUpR[38]+(0.223606797749979*alphaR[24]+0.159719141249985*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[37]+(0.223606797749979*(fUpR[29]+fUpR[24])+0.159719141249985*fUpR[19]+0.2500000000000001*fUpR[2])*alphaR[37]+0.2*alphaR[15]*fUpR[36]+(0.223606797749979*alphaR[22]+0.159719141249985*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[35]+0.2*(alphaR[17]*fUpR[34]+fUpR[17]*alphaR[34]+fUpR[16]*alphaR[33])+(0.223606797749979*alphaR[27]+0.159719141249985*alphaR[25]+0.2500000000000001*alphaR[4])*fUpR[32]+(0.223606797749979*(fUpR[27]+fUpR[26])+0.159719141249985*fUpR[25]+0.2500000000000001*fUpR[4])*alphaR[32]+(0.2*(alphaR[23]+alphaR[20])+0.223606797749979*alphaR[1])*fUpR[31]+0.25*(alphaR[7]*fUpR[25]+fUpR[7]*alphaR[25]+fUpR[9]*alphaR[21]+alphaR[10]*fUpR[19]+fUpR[10]*alphaR[19])+0.2500000000000001*alphaR[11]*fUpR[18]+0.223606797749979*(alphaR[5]*fUpR[17]+fUpR[5]*alphaR[17]+alphaR[6]*fUpR[16]+alphaR[8]*fUpR[15]+fUpR[8]*alphaR[15]); 
  GhatR[45] = (0.1788854381999831*alphaR[32]+0.2*alphaR[7])*fUpR[46]+(0.223606797749979*alphaR[13]+0.159719141249985*alphaR[12]+0.223606797749979*alphaR[11]+0.25*alphaR[0])*fUpR[45]+0.1788854381999831*alphaR[34]*fUpR[44]+0.2*(alphaR[5]*fUpR[44]+alphaR[15]*fUpR[40])+0.223606797749979*(alphaR[22]*fUpR[39]+fUpR[22]*alphaR[39])+(0.223606797749979*alphaR[23]+0.159719141249985*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[38]+0.223606797749979*(alphaR[20]*fUpR[37]+fUpR[20]*alphaR[37])+(0.159719141249985*alphaR[22]+0.223606797749979*alphaR[21]+0.2500000000000001*alphaR[3])*fUpR[36]+0.2*(alphaR[15]*fUpR[35]+fUpR[18]*alphaR[34])+(0.223606797749979*(alphaR[27]+alphaR[25])+0.2500000000000001*alphaR[4])*fUpR[33]+(0.223606797749979*fUpR[27]+0.159719141249985*fUpR[26]+0.223606797749979*fUpR[25]+0.2500000000000001*fUpR[4])*alphaR[33]+0.2*fUpR[16]*alphaR[32]+(0.2*(alphaR[24]+alphaR[19])+0.223606797749979*alphaR[2])*fUpR[31]+0.25*(alphaR[6]*fUpR[26]+alphaR[8]*fUpR[22]+fUpR[8]*alphaR[22]+alphaR[10]*fUpR[20]+fUpR[10]*alphaR[20])+0.223606797749979*alphaR[5]*fUpR[18]+0.2500000000000001*(alphaR[12]*fUpR[17]+fUpR[12]*alphaR[17])+0.223606797749979*(alphaR[7]*fUpR[16]+fUpR[9]*alphaR[15]); 
  GhatR[46] = (0.1788854381999831*alphaR[37]+0.2*alphaR[10])*fUpR[47]+(0.159719141249985*alphaR[13]+0.223606797749979*(alphaR[12]+alphaR[11])+0.25*alphaR[0])*fUpR[46]+(0.1788854381999831*alphaR[32]+0.2*alphaR[7])*fUpR[45]+0.1788854381999831*alphaR[33]*fUpR[44]+0.2*(alphaR[6]*fUpR[44]+alphaR[17]*fUpR[43])+0.223606797749979*alphaR[27]*fUpR[41]+(0.159719141249985*alphaR[23]+0.223606797749979*alphaR[20]+0.2500000000000001*alphaR[1])*fUpR[40]+(0.159719141249985*alphaR[24]+0.223606797749979*alphaR[19]+0.2500000000000001*alphaR[2])*fUpR[39]+(0.223606797749979*fUpR[29]+0.159719141249985*fUpR[24]+0.223606797749979*fUpR[19]+0.2500000000000001*fUpR[2])*alphaR[39]+0.2*(alphaR[15]*(fUpR[38]+fUpR[37])+fUpR[15]*alphaR[37])+0.223606797749979*(alphaR[24]*fUpR[36]+alphaR[23]*fUpR[35])+(0.159719141249985*alphaR[27]+0.223606797749979*alphaR[25]+0.2500000000000001*alphaR[4])*fUpR[34]+(0.159719141249985*fUpR[27]+0.223606797749979*(fUpR[26]+fUpR[25])+0.2500000000000001*fUpR[4])*alphaR[34]+0.2*(fUpR[18]*alphaR[33]+alphaR[17]*fUpR[32]+fUpR[17]*alphaR[32])+(0.2*(alphaR[22]+alphaR[21])+0.223606797749979*alphaR[3])*fUpR[31]+0.25*(alphaR[5]*fUpR[27]+fUpR[5]*alphaR[27]+alphaR[8]*fUpR[24]+fUpR[8]*alphaR[24]+fUpR[9]*alphaR[23])+0.223606797749979*(alphaR[6]*fUpR[18]+alphaR[7]*fUpR[17]+fUpR[7]*alphaR[17])+0.2500000000000001*alphaR[13]*fUpR[16]+0.223606797749979*(alphaR[10]*fUpR[15]+fUpR[10]*alphaR[15]); 
  GhatR[47] = (0.223606797749979*(alphaR[13]+alphaR[12]+alphaR[11])+0.25*alphaR[0])*fUpR[47]+(0.1788854381999831*alphaR[37]+0.2*alphaR[10])*fUpR[46]+(0.1788854381999831*alphaR[39]+0.2*alphaR[8])*fUpR[44]+(0.223606797749979*(alphaR[23]+alphaR[20])+0.2500000000000001*alphaR[1])*fUpR[43]+(0.223606797749979*(alphaR[24]+alphaR[19])+0.2500000000000001*alphaR[2])*fUpR[42]+(0.223606797749979*(alphaR[22]+alphaR[21])+0.2500000000000001*alphaR[3])*fUpR[41]+0.2*(alphaR[17]*fUpR[40]+fUpR[18]*alphaR[39]+fUpR[16]*alphaR[37]+alphaR[17]*fUpR[35])+0.223606797749979*(fUpR[30]*alphaR[34]+fUpR[29]*alphaR[33]+fUpR[28]*alphaR[32])+(0.2*(alphaR[27]+alphaR[25])+0.223606797749979*alphaR[4])*fUpR[31]+0.25*(alphaR[5]*fUpR[30]+alphaR[6]*fUpR[29]+alphaR[7]*fUpR[28])+0.223606797749979*(alphaR[8]*fUpR[18]+fUpR[9]*alphaR[17]+alphaR[10]*fUpR[16])+0.2500000000000001*fUpR[14]*alphaR[15]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[3] += -0.7071067811865475*GhatR[3]*rdvpar2; 
  out[4] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdvpar2; 
  out[6] += -0.7071067811865475*GhatR[5]*rdvpar2; 
  out[7] += -0.7071067811865475*GhatR[6]*rdvpar2; 
  out[8] += -0.7071067811865475*GhatR[7]*rdvpar2; 
  out[9] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[10] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[11] += -1.224744871391589*GhatR[3]*rdvpar2; 
  out[12] += -0.7071067811865475*GhatR[8]*rdvpar2; 
  out[13] += -0.7071067811865475*GhatR[9]*rdvpar2; 
  out[14] += -0.7071067811865475*GhatR[10]*rdvpar2; 
  out[15] += -1.224744871391589*GhatR[4]*rdvpar2; 
  out[16] += -0.7071067811865475*GhatR[11]*rdvpar2; 
  out[17] += -0.7071067811865475*GhatR[12]*rdvpar2; 
  out[18] += -0.7071067811865475*GhatR[13]*rdvpar2; 
  out[19] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[20] += -0.7071067811865475*GhatR[14]*rdvpar2; 
  out[21] += -0.7071067811865475*GhatR[15]*rdvpar2; 
  out[22] += -1.224744871391589*GhatR[5]*rdvpar2; 
  out[23] += -1.224744871391589*GhatR[6]*rdvpar2; 
  out[24] += -1.224744871391589*GhatR[7]*rdvpar2; 
  out[25] += -0.7071067811865475*GhatR[16]*rdvpar2; 
  out[26] += -0.7071067811865475*GhatR[17]*rdvpar2; 
  out[27] += -0.7071067811865475*GhatR[18]*rdvpar2; 
  out[28] += -1.224744871391589*GhatR[8]*rdvpar2; 
  out[29] += -1.224744871391589*GhatR[9]*rdvpar2; 
  out[30] += -1.224744871391589*GhatR[10]*rdvpar2; 
  out[31] += -0.7071067811865475*GhatR[19]*rdvpar2; 
  out[32] += -0.7071067811865475*GhatR[20]*rdvpar2; 
  out[33] += -0.7071067811865475*GhatR[21]*rdvpar2; 
  out[34] += -0.7071067811865475*GhatR[22]*rdvpar2; 
  out[35] += -0.7071067811865475*GhatR[23]*rdvpar2; 
  out[36] += -0.7071067811865475*GhatR[24]*rdvpar2; 
  out[37] += -1.224744871391589*GhatR[11]*rdvpar2; 
  out[38] += -1.224744871391589*GhatR[12]*rdvpar2; 
  out[39] += -1.224744871391589*GhatR[13]*rdvpar2; 
  out[40] += -1.58113883008419*GhatR[1]*rdvpar2; 
  out[41] += -1.58113883008419*GhatR[2]*rdvpar2; 
  out[42] += -1.58113883008419*GhatR[3]*rdvpar2; 
  out[43] += -0.7071067811865475*GhatR[25]*rdvpar2; 
  out[44] += -0.7071067811865475*GhatR[26]*rdvpar2; 
  out[45] += -0.7071067811865475*GhatR[27]*rdvpar2; 
  out[46] += -1.58113883008419*GhatR[4]*rdvpar2; 
  out[47] += -0.7071067811865475*GhatR[28]*rdvpar2; 
  out[48] += -0.7071067811865475*GhatR[29]*rdvpar2; 
  out[49] += -0.7071067811865475*GhatR[30]*rdvpar2; 
  out[50] += -1.224744871391589*GhatR[14]*rdvpar2; 
  out[51] += -1.224744871391589*GhatR[15]*rdvpar2; 
  out[52] += -0.7071067811865475*GhatR[31]*rdvpar2; 
  out[53] += -1.224744871391589*GhatR[16]*rdvpar2; 
  out[54] += -1.224744871391589*GhatR[17]*rdvpar2; 
  out[55] += -1.224744871391589*GhatR[18]*rdvpar2; 
  out[56] += -0.7071067811865475*GhatR[32]*rdvpar2; 
  out[57] += -0.7071067811865475*GhatR[33]*rdvpar2; 
  out[58] += -0.7071067811865475*GhatR[34]*rdvpar2; 
  out[59] += -1.224744871391589*GhatR[19]*rdvpar2; 
  out[60] += -1.224744871391589*GhatR[20]*rdvpar2; 
  out[61] += -1.224744871391589*GhatR[21]*rdvpar2; 
  out[62] += -1.224744871391589*GhatR[22]*rdvpar2; 
  out[63] += -1.224744871391589*GhatR[23]*rdvpar2; 
  out[64] += -1.224744871391589*GhatR[24]*rdvpar2; 
  out[65] += -1.58113883008419*GhatR[5]*rdvpar2; 
  out[66] += -1.58113883008419*GhatR[6]*rdvpar2; 
  out[67] += -1.58113883008419*GhatR[7]*rdvpar2; 
  out[68] += -0.7071067811865475*GhatR[35]*rdvpar2; 
  out[69] += -0.7071067811865475*GhatR[36]*rdvpar2; 
  out[70] += -0.7071067811865475*GhatR[37]*rdvpar2; 
  out[71] += -0.7071067811865475*GhatR[38]*rdvpar2; 
  out[72] += -0.7071067811865475*GhatR[39]*rdvpar2; 
  out[73] += -0.7071067811865475*GhatR[40]*rdvpar2; 
  out[74] += -1.224744871391589*GhatR[25]*rdvpar2; 
  out[75] += -1.224744871391589*GhatR[26]*rdvpar2; 
  out[76] += -1.224744871391589*GhatR[27]*rdvpar2; 
  out[77] += -1.58113883008419*GhatR[8]*rdvpar2; 
  out[78] += -1.58113883008419*GhatR[9]*rdvpar2; 
  out[79] += -1.58113883008419*GhatR[10]*rdvpar2; 
  out[80] += -0.7071067811865475*GhatR[41]*rdvpar2; 
  out[81] += -0.7071067811865475*GhatR[42]*rdvpar2; 
  out[82] += -0.7071067811865475*GhatR[43]*rdvpar2; 
  out[83] += -1.224744871391589*GhatR[28]*rdvpar2; 
  out[84] += -1.224744871391589*GhatR[29]*rdvpar2; 
  out[85] += -1.224744871391589*GhatR[30]*rdvpar2; 
  out[86] += -1.224744871391589*GhatR[31]*rdvpar2; 
  out[87] += -1.224744871391589*GhatR[32]*rdvpar2; 
  out[88] += -1.224744871391589*GhatR[33]*rdvpar2; 
  out[89] += -1.224744871391589*GhatR[34]*rdvpar2; 
  out[90] += -1.58113883008419*GhatR[15]*rdvpar2; 
  out[91] += -0.7071067811865475*GhatR[44]*rdvpar2; 
  out[92] += -0.7071067811865475*GhatR[45]*rdvpar2; 
  out[93] += -0.7071067811865475*GhatR[46]*rdvpar2; 
  out[94] += -1.224744871391589*GhatR[35]*rdvpar2; 
  out[95] += -1.224744871391589*GhatR[36]*rdvpar2; 
  out[96] += -1.224744871391589*GhatR[37]*rdvpar2; 
  out[97] += -1.224744871391589*GhatR[38]*rdvpar2; 
  out[98] += -1.224744871391589*GhatR[39]*rdvpar2; 
  out[99] += -1.224744871391589*GhatR[40]*rdvpar2; 
  out[100] += -1.58113883008419*GhatR[16]*rdvpar2; 
  out[101] += -1.58113883008419*GhatR[17]*rdvpar2; 
  out[102] += -1.58113883008419*GhatR[18]*rdvpar2; 
  out[103] += -0.7071067811865475*GhatR[47]*rdvpar2; 
  out[104] += -1.224744871391589*GhatR[41]*rdvpar2; 
  out[105] += -1.224744871391589*GhatR[42]*rdvpar2; 
  out[106] += -1.224744871391589*GhatR[43]*rdvpar2; 
  out[107] += -1.224744871391589*GhatR[44]*rdvpar2; 
  out[108] += -1.224744871391589*GhatR[45]*rdvpar2; 
  out[109] += -1.224744871391589*GhatR[46]*rdvpar2; 
  out[110] += -1.58113883008419*GhatR[31]*rdvpar2; 
  out[111] += -1.224744871391589*GhatR[47]*rdvpar2; 

  } else { 

  const double *alphaL = &alpha_skin[144];
  double fUpOrdL[81] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.4024922359499623*alphaL[39]+0.4024922359499623*alphaL[37]+0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]-0.3*alphaL[27]-0.3*alphaL[25]-0.3*alphaL[24]-0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]-0.3*alphaL[20]-0.3*alphaL[19]-0.603738353924943*alphaL[17]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]+0.45*alphaL[8]+0.45*alphaL[7]+0.45*alphaL[6]+0.45*alphaL[5]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = ser_5x_p2_surfx4_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = ser_5x_p2_surfx4_eval_quad_node_0_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]-0.3*alphaL[24]-0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]-0.3*alphaL[20]-0.3*alphaL[19]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[7]+0.45*alphaL[6]+0.45*alphaL[5]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = ser_5x_p2_surfx4_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = ser_5x_p2_surfx4_eval_quad_node_1_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])-0.4024922359499623*alphaL[37]+0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]+0.3*alphaL[27]+0.3*alphaL[25]-0.3*alphaL[24]-0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]-0.3*alphaL[20]-0.3*alphaL[19]+0.603738353924943*alphaL[17]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]-0.45*alphaL[8]+0.45*alphaL[7]+0.45*alphaL[6]+0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = ser_5x_p2_surfx4_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = ser_5x_p2_surfx4_eval_quad_node_2_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[39])-0.5031152949374518*alphaL[34]+0.375*alphaL[27]-0.3*alphaL[25]+0.375*alphaL[24]+0.375*alphaL[23]-0.3*alphaL[20]-0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[8]+0.45*alphaL[5]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = ser_5x_p2_surfx4_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = ser_5x_p2_surfx4_eval_quad_node_3_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[34])+0.375*alphaL[24]+0.375*alphaL[23]-0.3*alphaL[20]-0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[5]-0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = ser_5x_p2_surfx4_eval_quad_node_4_r(fedge); 
  } else { 
    fUpOrdL[4] = ser_5x_p2_surfx4_eval_quad_node_4_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[39]-0.5031152949374518*alphaL[34]-0.375*alphaL[27]+0.3*alphaL[25]+0.375*alphaL[24]+0.375*alphaL[23]-0.3*alphaL[20]-0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[8]+0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = ser_5x_p2_surfx4_eval_quad_node_5_r(fedge); 
  } else { 
    fUpOrdL[5] = ser_5x_p2_surfx4_eval_quad_node_5_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]+0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]-0.3*alphaL[27]-0.3*alphaL[25]-0.3*alphaL[24]-0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]-0.3*alphaL[20]-0.3*alphaL[19]+0.603738353924943*alphaL[17]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[7]-0.45*alphaL[6]+0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = ser_5x_p2_surfx4_eval_quad_node_6_r(fedge); 
  } else { 
    fUpOrdL[6] = ser_5x_p2_surfx4_eval_quad_node_6_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]-0.3*alphaL[24]-0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]-0.3*alphaL[20]-0.3*alphaL[19]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[7]-0.45*alphaL[6]+0.45*alphaL[5]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = ser_5x_p2_surfx4_eval_quad_node_7_r(fedge); 
  } else { 
    fUpOrdL[7] = ser_5x_p2_surfx4_eval_quad_node_7_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])+0.4024922359499623*alphaL[37]+0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]+0.3*alphaL[27]+0.3*alphaL[25]-0.3*alphaL[24]-0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]-0.3*alphaL[20]-0.3*alphaL[19]-0.603738353924943*alphaL[17]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]-0.45*alphaL[8]-0.45*alphaL[7]-0.45*alphaL[6]+0.45*alphaL[5]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = ser_5x_p2_surfx4_eval_quad_node_8_r(fedge); 
  } else { 
    fUpOrdL[8] = ser_5x_p2_surfx4_eval_quad_node_8_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]+0.4024922359499623*alphaL[37]-0.5031152949374518*alphaL[33]-0.3*alphaL[27]-0.3*alphaL[25]-0.3*alphaL[23]+0.375*alphaL[22]-0.3*alphaL[21]+0.375*alphaL[20]-0.603738353924943*alphaL[17]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]+0.45*alphaL[8]+0.45*alphaL[6]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[9] = ser_5x_p2_surfx4_eval_quad_node_9_r(fedge); 
  } else { 
    fUpOrdL[9] = ser_5x_p2_surfx4_eval_quad_node_9_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[33])-0.3*alphaL[23]+0.375*alphaL[22]-0.3*alphaL[21]+0.375*alphaL[20]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[6]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[10] = ser_5x_p2_surfx4_eval_quad_node_10_r(fedge); 
  } else { 
    fUpOrdL[10] = ser_5x_p2_surfx4_eval_quad_node_10_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])-0.4024922359499623*alphaL[37]-0.5031152949374518*alphaL[33]+0.3*alphaL[27]+0.3*alphaL[25]-0.3*alphaL[23]+0.375*alphaL[22]-0.3*alphaL[21]+0.375*alphaL[20]+0.603738353924943*alphaL[17]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]-0.45*alphaL[8]+0.45*alphaL[6]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[11] = ser_5x_p2_surfx4_eval_quad_node_11_r(fedge); 
  } else { 
    fUpOrdL[11] = ser_5x_p2_surfx4_eval_quad_node_11_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[39])+0.375*alphaL[27]-0.3*alphaL[25]+0.375*alphaL[23]+0.375*alphaL[20]-0.2795084971874732*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[8]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[12] = ser_5x_p2_surfx4_eval_quad_node_12_r(fedge); 
  } else { 
    fUpOrdL[12] = ser_5x_p2_surfx4_eval_quad_node_12_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.375*alphaL[23]+0.375*alphaL[20]-0.2795084971874732*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[13] = ser_5x_p2_surfx4_eval_quad_node_13_r(fedge); 
  } else { 
    fUpOrdL[13] = ser_5x_p2_surfx4_eval_quad_node_13_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[39]-0.375*alphaL[27]+0.3*alphaL[25]+0.375*alphaL[23]+0.375*alphaL[20]-0.2795084971874732*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[8]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[14] = ser_5x_p2_surfx4_eval_quad_node_14_r(fedge); 
  } else { 
    fUpOrdL[14] = ser_5x_p2_surfx4_eval_quad_node_14_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]+0.5031152949374518*alphaL[33]-0.3*alphaL[27]-0.3*alphaL[25]-0.3*alphaL[23]-0.375*alphaL[22]+0.3*alphaL[21]+0.375*alphaL[20]+0.603738353924943*alphaL[17]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[6]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[15] = ser_5x_p2_surfx4_eval_quad_node_15_r(fedge); 
  } else { 
    fUpOrdL[15] = ser_5x_p2_surfx4_eval_quad_node_15_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[33]-0.3*alphaL[23]-0.375*alphaL[22]+0.3*alphaL[21]+0.375*alphaL[20]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[6]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[16] = ser_5x_p2_surfx4_eval_quad_node_16_r(fedge); 
  } else { 
    fUpOrdL[16] = ser_5x_p2_surfx4_eval_quad_node_16_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])+0.4024922359499623*alphaL[37]+0.5031152949374518*alphaL[33]+0.3*alphaL[27]+0.3*alphaL[25]-0.3*alphaL[23]-0.375*alphaL[22]+0.3*alphaL[21]+0.375*alphaL[20]-0.603738353924943*alphaL[17]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]-0.45*alphaL[8]-0.45*alphaL[6]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[17] = ser_5x_p2_surfx4_eval_quad_node_17_r(fedge); 
  } else { 
    fUpOrdL[17] = ser_5x_p2_surfx4_eval_quad_node_17_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]+0.4024922359499623*alphaL[37]-0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]-0.3*alphaL[27]-0.3*alphaL[25]+0.3*alphaL[24]-0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]-0.3*alphaL[20]+0.3*alphaL[19]-0.603738353924943*alphaL[17]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[18] = ser_5x_p2_surfx4_eval_quad_node_18_r(fedge); 
  } else { 
    fUpOrdL[18] = ser_5x_p2_surfx4_eval_quad_node_18_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[34])+0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]+0.3*alphaL[24]-0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]-0.3*alphaL[20]+0.3*alphaL[19]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[19] = ser_5x_p2_surfx4_eval_quad_node_19_r(fedge); 
  } else { 
    fUpOrdL[19] = ser_5x_p2_surfx4_eval_quad_node_19_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])-0.4024922359499623*alphaL[37]-0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]+0.3*alphaL[27]+0.3*alphaL[25]+0.3*alphaL[24]-0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]-0.3*alphaL[20]+0.3*alphaL[19]+0.603738353924943*alphaL[17]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]-0.45*alphaL[8]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[20] = ser_5x_p2_surfx4_eval_quad_node_20_r(fedge); 
  } else { 
    fUpOrdL[20] = ser_5x_p2_surfx4_eval_quad_node_20_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[39])+0.5031152949374518*alphaL[34]+0.375*alphaL[27]-0.3*alphaL[25]-0.375*alphaL[24]+0.375*alphaL[23]-0.3*alphaL[20]+0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[8]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[21] = ser_5x_p2_surfx4_eval_quad_node_21_r(fedge); 
  } else { 
    fUpOrdL[21] = ser_5x_p2_surfx4_eval_quad_node_21_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[34]-0.375*alphaL[24]+0.375*alphaL[23]-0.3*alphaL[20]+0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[5]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[22] = ser_5x_p2_surfx4_eval_quad_node_22_r(fedge); 
  } else { 
    fUpOrdL[22] = ser_5x_p2_surfx4_eval_quad_node_22_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[39]+0.5031152949374518*alphaL[34]-0.375*alphaL[27]+0.3*alphaL[25]-0.375*alphaL[24]+0.375*alphaL[23]-0.3*alphaL[20]+0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[8]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[23] = ser_5x_p2_surfx4_eval_quad_node_23_r(fedge); 
  } else { 
    fUpOrdL[23] = ser_5x_p2_surfx4_eval_quad_node_23_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]-0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]-0.3*alphaL[27]-0.3*alphaL[25]+0.3*alphaL[24]-0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]-0.3*alphaL[20]+0.3*alphaL[19]+0.603738353924943*alphaL[17]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[8]+0.45*alphaL[7]-0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[24] = ser_5x_p2_surfx4_eval_quad_node_24_r(fedge); 
  } else { 
    fUpOrdL[24] = ser_5x_p2_surfx4_eval_quad_node_24_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[34])-0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]+0.3*alphaL[24]-0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]-0.3*alphaL[20]+0.3*alphaL[19]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[7]-0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[25] = ser_5x_p2_surfx4_eval_quad_node_25_r(fedge); 
  } else { 
    fUpOrdL[25] = ser_5x_p2_surfx4_eval_quad_node_25_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])+0.4024922359499623*alphaL[37]-0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]+0.3*alphaL[27]+0.3*alphaL[25]+0.3*alphaL[24]-0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]-0.3*alphaL[20]+0.3*alphaL[19]-0.603738353924943*alphaL[17]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]-0.45*alphaL[8]+0.45*alphaL[7]-0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]-0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[26] = ser_5x_p2_surfx4_eval_quad_node_26_r(fedge); 
  } else { 
    fUpOrdL[26] = ser_5x_p2_surfx4_eval_quad_node_26_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[37])-0.5031152949374518*alphaL[32]-0.3*alphaL[27]+0.375*alphaL[25]-0.3*alphaL[24]-0.3*alphaL[22]+0.375*alphaL[21]+0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]+0.45*alphaL[10]+0.45*alphaL[7]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[27] = ser_5x_p2_surfx4_eval_quad_node_27_r(fedge); 
  } else { 
    fUpOrdL[27] = ser_5x_p2_surfx4_eval_quad_node_27_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[32])-0.3*alphaL[24]-0.3*alphaL[22]+0.375*alphaL[21]+0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]+0.45*alphaL[7]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[28] = ser_5x_p2_surfx4_eval_quad_node_28_r(fedge); 
  } else { 
    fUpOrdL[28] = ser_5x_p2_surfx4_eval_quad_node_28_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[37]-0.5031152949374518*alphaL[32]+0.3*alphaL[27]-0.375*alphaL[25]-0.3*alphaL[24]-0.3*alphaL[22]+0.375*alphaL[21]+0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[7]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[29] = ser_5x_p2_surfx4_eval_quad_node_29_r(fedge); 
  } else { 
    fUpOrdL[29] = ser_5x_p2_surfx4_eval_quad_node_29_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.375*alphaL[27]+0.375*alphaL[25]+0.375*alphaL[24]+0.375*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[30] = ser_5x_p2_surfx4_eval_quad_node_30_r(fedge); 
  } else { 
    fUpOrdL[30] = ser_5x_p2_surfx4_eval_quad_node_30_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.375*alphaL[24]+0.375*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]-0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[31] = ser_5x_p2_surfx4_eval_quad_node_31_r(fedge); 
  } else { 
    fUpOrdL[31] = ser_5x_p2_surfx4_eval_quad_node_31_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.375*alphaL[27])-0.375*alphaL[25]+0.375*alphaL[24]+0.375*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[32] = ser_5x_p2_surfx4_eval_quad_node_32_r(fedge); 
  } else { 
    fUpOrdL[32] = ser_5x_p2_surfx4_eval_quad_node_32_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[37]+0.5031152949374518*alphaL[32]-0.3*alphaL[27]+0.375*alphaL[25]-0.3*alphaL[24]+0.3*alphaL[22]-0.375*alphaL[21]+0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]-0.45*alphaL[10]-0.45*alphaL[7]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[33] = ser_5x_p2_surfx4_eval_quad_node_33_r(fedge); 
  } else { 
    fUpOrdL[33] = ser_5x_p2_surfx4_eval_quad_node_33_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[32]-0.3*alphaL[24]+0.3*alphaL[22]-0.375*alphaL[21]+0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]-0.45*alphaL[7]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[34] = ser_5x_p2_surfx4_eval_quad_node_34_r(fedge); 
  } else { 
    fUpOrdL[34] = ser_5x_p2_surfx4_eval_quad_node_34_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[37])+0.5031152949374518*alphaL[32]+0.3*alphaL[27]-0.375*alphaL[25]-0.3*alphaL[24]+0.3*alphaL[22]-0.375*alphaL[21]+0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]+0.45*alphaL[10]-0.45*alphaL[7]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[35] = ser_5x_p2_surfx4_eval_quad_node_35_r(fedge); 
  } else { 
    fUpOrdL[35] = ser_5x_p2_surfx4_eval_quad_node_35_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[37])-0.3*alphaL[27]+0.375*alphaL[25]+0.375*alphaL[22]+0.375*alphaL[21]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]-0.2795084971874732*alphaL[11]+0.45*alphaL[10]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[36] = ser_5x_p2_surfx4_eval_quad_node_36_r(fedge); 
  } else { 
    fUpOrdL[36] = ser_5x_p2_surfx4_eval_quad_node_36_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.375*alphaL[22]+0.375*alphaL[21]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]-0.2795084971874732*alphaL[11]-0.3354101966249678*alphaL[3]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[37] = ser_5x_p2_surfx4_eval_quad_node_37_r(fedge); 
  } else { 
    fUpOrdL[37] = ser_5x_p2_surfx4_eval_quad_node_37_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[37]+0.3*alphaL[27]-0.375*alphaL[25]+0.375*alphaL[22]+0.375*alphaL[21]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]-0.2795084971874732*alphaL[11]-0.45*alphaL[10]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[38] = ser_5x_p2_surfx4_eval_quad_node_38_r(fedge); 
  } else { 
    fUpOrdL[38] = ser_5x_p2_surfx4_eval_quad_node_38_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.375*alphaL[27]+0.375*alphaL[25]-0.2795084971874732*alphaL[13]-0.2795084971874732*alphaL[12]-0.2795084971874732*alphaL[11]-0.3354101966249678*alphaL[4]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[39] = ser_5x_p2_surfx4_eval_quad_node_39_r(fedge); 
  } else { 
    fUpOrdL[39] = ser_5x_p2_surfx4_eval_quad_node_39_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[13])-0.2795084971874732*alphaL[12]-0.2795084971874732*alphaL[11]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[40] = ser_5x_p2_surfx4_eval_quad_node_40_r(fedge); 
  } else { 
    fUpOrdL[40] = ser_5x_p2_surfx4_eval_quad_node_40_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.375*alphaL[27])-0.375*alphaL[25]-0.2795084971874732*alphaL[13]-0.2795084971874732*alphaL[12]-0.2795084971874732*alphaL[11]+0.3354101966249678*alphaL[4]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[41] = ser_5x_p2_surfx4_eval_quad_node_41_r(fedge); 
  } else { 
    fUpOrdL[41] = ser_5x_p2_surfx4_eval_quad_node_41_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[37]-0.3*alphaL[27]+0.375*alphaL[25]-0.375*alphaL[22]-0.375*alphaL[21]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]-0.2795084971874732*alphaL[11]-0.45*alphaL[10]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[42] = ser_5x_p2_surfx4_eval_quad_node_42_r(fedge); 
  } else { 
    fUpOrdL[42] = ser_5x_p2_surfx4_eval_quad_node_42_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.375*alphaL[22])-0.375*alphaL[21]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]-0.2795084971874732*alphaL[11]+0.3354101966249678*alphaL[3]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[43] = ser_5x_p2_surfx4_eval_quad_node_43_r(fedge); 
  } else { 
    fUpOrdL[43] = ser_5x_p2_surfx4_eval_quad_node_43_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[37])+0.3*alphaL[27]-0.375*alphaL[25]-0.375*alphaL[22]-0.375*alphaL[21]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]-0.2795084971874732*alphaL[11]+0.45*alphaL[10]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[44] = ser_5x_p2_surfx4_eval_quad_node_44_r(fedge); 
  } else { 
    fUpOrdL[44] = ser_5x_p2_surfx4_eval_quad_node_44_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[37])+0.5031152949374518*alphaL[32]-0.3*alphaL[27]+0.375*alphaL[25]+0.3*alphaL[24]-0.3*alphaL[22]+0.375*alphaL[21]-0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]+0.45*alphaL[10]-0.45*alphaL[7]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[45] = ser_5x_p2_surfx4_eval_quad_node_45_r(fedge); 
  } else { 
    fUpOrdL[45] = ser_5x_p2_surfx4_eval_quad_node_45_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[32]+0.3*alphaL[24]-0.3*alphaL[22]+0.375*alphaL[21]-0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]-0.45*alphaL[7]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[46] = ser_5x_p2_surfx4_eval_quad_node_46_r(fedge); 
  } else { 
    fUpOrdL[46] = ser_5x_p2_surfx4_eval_quad_node_46_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[37]+0.5031152949374518*alphaL[32]+0.3*alphaL[27]-0.375*alphaL[25]+0.3*alphaL[24]-0.3*alphaL[22]+0.375*alphaL[21]-0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]-0.45*alphaL[10]-0.45*alphaL[7]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[47] = ser_5x_p2_surfx4_eval_quad_node_47_r(fedge); 
  } else { 
    fUpOrdL[47] = ser_5x_p2_surfx4_eval_quad_node_47_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.375*alphaL[27]+0.375*alphaL[25]-0.375*alphaL[24]-0.375*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[48] = ser_5x_p2_surfx4_eval_quad_node_48_r(fedge); 
  } else { 
    fUpOrdL[48] = ser_5x_p2_surfx4_eval_quad_node_48_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.375*alphaL[24])-0.375*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]+0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[49] = ser_5x_p2_surfx4_eval_quad_node_49_r(fedge); 
  } else { 
    fUpOrdL[49] = ser_5x_p2_surfx4_eval_quad_node_49_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.375*alphaL[27])-0.375*alphaL[25]-0.375*alphaL[24]-0.375*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[50] = ser_5x_p2_surfx4_eval_quad_node_50_r(fedge); 
  } else { 
    fUpOrdL[50] = ser_5x_p2_surfx4_eval_quad_node_50_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[37]-0.5031152949374518*alphaL[32]-0.3*alphaL[27]+0.375*alphaL[25]+0.3*alphaL[24]+0.3*alphaL[22]-0.375*alphaL[21]-0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[7]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[51] = ser_5x_p2_surfx4_eval_quad_node_51_r(fedge); 
  } else { 
    fUpOrdL[51] = ser_5x_p2_surfx4_eval_quad_node_51_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[32])+0.3*alphaL[24]+0.3*alphaL[22]-0.375*alphaL[21]-0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]+0.45*alphaL[7]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[52] = ser_5x_p2_surfx4_eval_quad_node_52_r(fedge); 
  } else { 
    fUpOrdL[52] = ser_5x_p2_surfx4_eval_quad_node_52_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[37])-0.5031152949374518*alphaL[32]+0.3*alphaL[27]-0.375*alphaL[25]+0.3*alphaL[24]+0.3*alphaL[22]-0.375*alphaL[21]-0.375*alphaL[19]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]-0.2795084971874732*alphaL[11]+0.45*alphaL[10]+0.45*alphaL[7]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[53] = ser_5x_p2_surfx4_eval_quad_node_53_r(fedge); 
  } else { 
    fUpOrdL[53] = ser_5x_p2_surfx4_eval_quad_node_53_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])+0.4024922359499623*alphaL[37]-0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]-0.3*alphaL[27]-0.3*alphaL[25]-0.3*alphaL[24]+0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]+0.3*alphaL[20]-0.3*alphaL[19]+0.603738353924943*alphaL[17]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]-0.45*alphaL[8]+0.45*alphaL[7]-0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[54] = ser_5x_p2_surfx4_eval_quad_node_54_r(fedge); 
  } else { 
    fUpOrdL[54] = ser_5x_p2_surfx4_eval_quad_node_54_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[34])-0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]-0.3*alphaL[24]+0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]+0.3*alphaL[20]-0.3*alphaL[19]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[7]-0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[55] = ser_5x_p2_surfx4_eval_quad_node_55_r(fedge); 
  } else { 
    fUpOrdL[55] = ser_5x_p2_surfx4_eval_quad_node_55_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]-0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]+0.3*alphaL[27]+0.3*alphaL[25]-0.3*alphaL[24]+0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]+0.3*alphaL[20]-0.3*alphaL[19]-0.603738353924943*alphaL[17]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[8]+0.45*alphaL[7]-0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[56] = ser_5x_p2_surfx4_eval_quad_node_56_r(fedge); 
  } else { 
    fUpOrdL[56] = ser_5x_p2_surfx4_eval_quad_node_56_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[39]+0.5031152949374518*alphaL[34]+0.375*alphaL[27]-0.3*alphaL[25]+0.375*alphaL[24]-0.375*alphaL[23]+0.3*alphaL[20]-0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[8]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[57] = ser_5x_p2_surfx4_eval_quad_node_57_r(fedge); 
  } else { 
    fUpOrdL[57] = ser_5x_p2_surfx4_eval_quad_node_57_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[34]+0.375*alphaL[24]-0.375*alphaL[23]+0.3*alphaL[20]-0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[5]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[58] = ser_5x_p2_surfx4_eval_quad_node_58_r(fedge); 
  } else { 
    fUpOrdL[58] = ser_5x_p2_surfx4_eval_quad_node_58_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[39])+0.5031152949374518*alphaL[34]-0.375*alphaL[27]+0.3*alphaL[25]+0.375*alphaL[24]-0.375*alphaL[23]+0.3*alphaL[20]-0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[8]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[59] = ser_5x_p2_surfx4_eval_quad_node_59_r(fedge); 
  } else { 
    fUpOrdL[59] = ser_5x_p2_surfx4_eval_quad_node_59_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])-0.4024922359499623*alphaL[37]-0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]-0.3*alphaL[27]-0.3*alphaL[25]-0.3*alphaL[24]+0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]+0.3*alphaL[20]-0.3*alphaL[19]-0.603738353924943*alphaL[17]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]-0.45*alphaL[8]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[60] = ser_5x_p2_surfx4_eval_quad_node_60_r(fedge); 
  } else { 
    fUpOrdL[60] = ser_5x_p2_surfx4_eval_quad_node_60_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[34])+0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]-0.3*alphaL[24]+0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]+0.3*alphaL[20]-0.3*alphaL[19]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[61] = ser_5x_p2_surfx4_eval_quad_node_61_r(fedge); 
  } else { 
    fUpOrdL[61] = ser_5x_p2_surfx4_eval_quad_node_61_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]+0.4024922359499623*alphaL[37]-0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]+0.3*alphaL[27]+0.3*alphaL[25]-0.3*alphaL[24]+0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]+0.3*alphaL[20]-0.3*alphaL[19]+0.603738353924943*alphaL[17]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[7]+0.45*alphaL[6]-0.45*alphaL[5]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]-0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[62] = ser_5x_p2_surfx4_eval_quad_node_62_r(fedge); 
  } else { 
    fUpOrdL[62] = ser_5x_p2_surfx4_eval_quad_node_62_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])+0.4024922359499623*alphaL[37]+0.5031152949374518*alphaL[33]-0.3*alphaL[27]-0.3*alphaL[25]+0.3*alphaL[23]+0.375*alphaL[22]-0.3*alphaL[21]-0.375*alphaL[20]+0.603738353924943*alphaL[17]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]-0.45*alphaL[8]-0.45*alphaL[6]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[63] = ser_5x_p2_surfx4_eval_quad_node_63_r(fedge); 
  } else { 
    fUpOrdL[63] = ser_5x_p2_surfx4_eval_quad_node_63_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[33]+0.3*alphaL[23]+0.375*alphaL[22]-0.3*alphaL[21]-0.375*alphaL[20]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[6]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[64] = ser_5x_p2_surfx4_eval_quad_node_64_r(fedge); 
  } else { 
    fUpOrdL[64] = ser_5x_p2_surfx4_eval_quad_node_64_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]+0.5031152949374518*alphaL[33]+0.3*alphaL[27]+0.3*alphaL[25]+0.3*alphaL[23]+0.375*alphaL[22]-0.3*alphaL[21]-0.375*alphaL[20]-0.603738353924943*alphaL[17]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[6]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[65] = ser_5x_p2_surfx4_eval_quad_node_65_r(fedge); 
  } else { 
    fUpOrdL[65] = ser_5x_p2_surfx4_eval_quad_node_65_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[39]+0.375*alphaL[27]-0.3*alphaL[25]-0.375*alphaL[23]-0.375*alphaL[20]-0.2795084971874732*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[8]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[66] = ser_5x_p2_surfx4_eval_quad_node_66_r(fedge); 
  } else { 
    fUpOrdL[66] = ser_5x_p2_surfx4_eval_quad_node_66_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.375*alphaL[23])-0.375*alphaL[20]-0.2795084971874732*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[67] = ser_5x_p2_surfx4_eval_quad_node_67_r(fedge); 
  } else { 
    fUpOrdL[67] = ser_5x_p2_surfx4_eval_quad_node_67_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[39])-0.375*alphaL[27]+0.3*alphaL[25]-0.375*alphaL[23]-0.375*alphaL[20]-0.2795084971874732*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[8]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[68] = ser_5x_p2_surfx4_eval_quad_node_68_r(fedge); 
  } else { 
    fUpOrdL[68] = ser_5x_p2_surfx4_eval_quad_node_68_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])-0.4024922359499623*alphaL[37]-0.5031152949374518*alphaL[33]-0.3*alphaL[27]-0.3*alphaL[25]+0.3*alphaL[23]-0.375*alphaL[22]+0.3*alphaL[21]-0.375*alphaL[20]-0.603738353924943*alphaL[17]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]-0.45*alphaL[8]+0.45*alphaL[6]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[69] = ser_5x_p2_surfx4_eval_quad_node_69_r(fedge); 
  } else { 
    fUpOrdL[69] = ser_5x_p2_surfx4_eval_quad_node_69_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[33])+0.3*alphaL[23]-0.375*alphaL[22]+0.3*alphaL[21]-0.375*alphaL[20]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[6]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[70] = ser_5x_p2_surfx4_eval_quad_node_70_r(fedge); 
  } else { 
    fUpOrdL[70] = ser_5x_p2_surfx4_eval_quad_node_70_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]+0.4024922359499623*alphaL[37]-0.5031152949374518*alphaL[33]+0.3*alphaL[27]+0.3*alphaL[25]+0.3*alphaL[23]-0.375*alphaL[22]+0.3*alphaL[21]-0.375*alphaL[20]+0.603738353924943*alphaL[17]+0.2236067977499786*alphaL[13]-0.2795084971874732*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]+0.45*alphaL[8]+0.45*alphaL[6]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[71] = ser_5x_p2_surfx4_eval_quad_node_71_r(fedge); 
  } else { 
    fUpOrdL[71] = ser_5x_p2_surfx4_eval_quad_node_71_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])+0.4024922359499623*alphaL[37]+0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]-0.3*alphaL[27]-0.3*alphaL[25]+0.3*alphaL[24]+0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]+0.3*alphaL[20]+0.3*alphaL[19]+0.603738353924943*alphaL[17]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]-0.45*alphaL[8]-0.45*alphaL[7]-0.45*alphaL[6]+0.45*alphaL[5]-0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[72] = ser_5x_p2_surfx4_eval_quad_node_72_r(fedge); 
  } else { 
    fUpOrdL[72] = ser_5x_p2_surfx4_eval_quad_node_72_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]+0.3*alphaL[24]+0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]+0.3*alphaL[20]+0.3*alphaL[19]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[7]-0.45*alphaL[6]+0.45*alphaL[5]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[73] = ser_5x_p2_surfx4_eval_quad_node_73_r(fedge); 
  } else { 
    fUpOrdL[73] = ser_5x_p2_surfx4_eval_quad_node_73_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]-0.4024922359499623*alphaL[37]+0.4024922359499623*alphaL[34]-0.4024922359499623*alphaL[33]-0.4024922359499623*alphaL[32]+0.3*alphaL[27]+0.3*alphaL[25]+0.3*alphaL[24]+0.3*alphaL[23]-0.3*alphaL[22]-0.3*alphaL[21]+0.3*alphaL[20]+0.3*alphaL[19]-0.603738353924943*alphaL[17]-0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]+0.45*alphaL[8]-0.45*alphaL[7]-0.45*alphaL[6]+0.45*alphaL[5]+0.3354101966249678*alphaL[4]-0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[74] = ser_5x_p2_surfx4_eval_quad_node_74_r(fedge); 
  } else { 
    fUpOrdL[74] = ser_5x_p2_surfx4_eval_quad_node_74_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5031152949374518*alphaL[39]-0.5031152949374518*alphaL[34]+0.375*alphaL[27]-0.3*alphaL[25]-0.375*alphaL[24]-0.375*alphaL[23]+0.3*alphaL[20]+0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[8]+0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[75] = ser_5x_p2_surfx4_eval_quad_node_75_r(fedge); 
  } else { 
    fUpOrdL[75] = ser_5x_p2_surfx4_eval_quad_node_75_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[34])-0.375*alphaL[24]-0.375*alphaL[23]+0.3*alphaL[20]+0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[5]+0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[76] = ser_5x_p2_surfx4_eval_quad_node_76_r(fedge); 
  } else { 
    fUpOrdL[76] = ser_5x_p2_surfx4_eval_quad_node_76_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5031152949374518*alphaL[39])-0.5031152949374518*alphaL[34]-0.375*alphaL[27]+0.3*alphaL[25]-0.375*alphaL[24]-0.375*alphaL[23]+0.3*alphaL[20]+0.3*alphaL[19]-0.2795084971874732*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[8]+0.45*alphaL[5]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[77] = ser_5x_p2_surfx4_eval_quad_node_77_r(fedge); 
  } else { 
    fUpOrdL[77] = ser_5x_p2_surfx4_eval_quad_node_77_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4024922359499623*alphaL[39])-0.4024922359499623*alphaL[37]+0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]-0.3*alphaL[27]-0.3*alphaL[25]+0.3*alphaL[24]+0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]+0.3*alphaL[20]+0.3*alphaL[19]-0.603738353924943*alphaL[17]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]-0.45*alphaL[10]-0.45*alphaL[8]+0.45*alphaL[7]+0.45*alphaL[6]+0.45*alphaL[5]-0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[78] = ser_5x_p2_surfx4_eval_quad_node_78_r(fedge); 
  } else { 
    fUpOrdL[78] = ser_5x_p2_surfx4_eval_quad_node_78_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]+0.3*alphaL[24]+0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]+0.3*alphaL[20]+0.3*alphaL[19]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[7]+0.45*alphaL[6]+0.45*alphaL[5]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[79] = ser_5x_p2_surfx4_eval_quad_node_79_r(fedge); 
  } else { 
    fUpOrdL[79] = ser_5x_p2_surfx4_eval_quad_node_79_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4024922359499623*alphaL[39]+0.4024922359499623*alphaL[37]+0.4024922359499623*alphaL[34]+0.4024922359499623*alphaL[33]+0.4024922359499623*alphaL[32]+0.3*alphaL[27]+0.3*alphaL[25]+0.3*alphaL[24]+0.3*alphaL[23]+0.3*alphaL[22]+0.3*alphaL[21]+0.3*alphaL[20]+0.3*alphaL[19]+0.603738353924943*alphaL[17]+0.603738353924943*alphaL[15]+0.2236067977499786*alphaL[13]+0.2236067977499786*alphaL[12]+0.2236067977499786*alphaL[11]+0.45*alphaL[10]+0.45*alphaL[8]+0.45*alphaL[7]+0.45*alphaL[6]+0.45*alphaL[5]+0.3354101966249678*alphaL[4]+0.3354101966249678*alphaL[3]+0.3354101966249678*alphaL[2]+0.3354101966249678*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[80] = ser_5x_p2_surfx4_eval_quad_node_80_r(fedge); 
  } else { 
    fUpOrdL[80] = ser_5x_p2_surfx4_eval_quad_node_80_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[48] = {0.};
  ser_5x_p2_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[48] = {0.}; 
  GhatL[0] = 0.25*(alphaL[39]*fUpL[39]+alphaL[37]*fUpL[37]+alphaL[34]*fUpL[34]+alphaL[33]*fUpL[33]+alphaL[32]*fUpL[32]+alphaL[27]*fUpL[27]+alphaL[25]*fUpL[25]+alphaL[24]*fUpL[24]+alphaL[23]*fUpL[23]+alphaL[22]*fUpL[22]+alphaL[21]*fUpL[21]+alphaL[20]*fUpL[20]+alphaL[19]*fUpL[19]+alphaL[17]*fUpL[17]+alphaL[15]*fUpL[15]+alphaL[13]*fUpL[13]+alphaL[12]*fUpL[12]+alphaL[11]*fUpL[11]+alphaL[10]*fUpL[10]+alphaL[8]*fUpL[8]+alphaL[7]*fUpL[7]+alphaL[6]*fUpL[6]+alphaL[5]*fUpL[5]+alphaL[4]*fUpL[4]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.2500000000000001*(alphaL[27]*fUpL[39]+fUpL[27]*alphaL[39])+0.223606797749979*(alphaL[17]*fUpL[37]+fUpL[17]*alphaL[37])+0.2500000000000001*(alphaL[24]*fUpL[34]+fUpL[24]*alphaL[34]+alphaL[22]*fUpL[33]+fUpL[22]*alphaL[33])+0.223606797749979*(alphaL[15]*fUpL[32]+fUpL[15]*alphaL[32])+0.223606797749979*(alphaL[8]*fUpL[25]+fUpL[8]*alphaL[25])+0.2500000000000001*(alphaL[13]*fUpL[23]+fUpL[13]*alphaL[23])+0.223606797749979*(alphaL[6]*fUpL[21]+fUpL[6]*alphaL[21])+0.2500000000000001*(alphaL[12]*fUpL[20]+fUpL[12]*alphaL[20])+0.223606797749979*(alphaL[5]*fUpL[19]+fUpL[5]*alphaL[19])+0.25*(alphaL[10]*fUpL[17]+fUpL[10]*alphaL[17]+alphaL[7]*fUpL[15]+fUpL[7]*alphaL[15])+0.223606797749979*(alphaL[1]*fUpL[11]+fUpL[1]*alphaL[11])+0.25*(alphaL[4]*fUpL[8]+fUpL[4]*alphaL[8]+alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6]+alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.2500000000000001*(alphaL[39]*fUpL[46]+alphaL[37]*fUpL[44]+alphaL[27]*fUpL[40]+alphaL[25]*fUpL[35]+alphaL[23]*fUpL[34]+fUpL[23]*alphaL[34])+0.223606797749979*(alphaL[15]*fUpL[33]+fUpL[15]*alphaL[33])+0.2500000000000001*(alphaL[21]*fUpL[32]+fUpL[21]*alphaL[32])+0.25*alphaL[17]*fUpL[31]+0.2500000000000001*(alphaL[13]*fUpL[24]+fUpL[13]*alphaL[24])+0.223606797749979*(alphaL[7]*fUpL[22]+fUpL[7]*alphaL[22]+alphaL[5]*fUpL[20]+fUpL[5]*alphaL[20])+0.2500000000000001*(alphaL[11]*fUpL[19]+fUpL[11]*alphaL[19])+0.25*(alphaL[10]*fUpL[18]+alphaL[8]*fUpL[16]+alphaL[6]*fUpL[15]+fUpL[6]*alphaL[15])+0.223606797749979*(alphaL[2]*fUpL[12]+fUpL[2]*alphaL[12])+0.25*(alphaL[4]*fUpL[9]+alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]+alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.223606797749979*(alphaL[17]*fUpL[39]+fUpL[17]*alphaL[39])+0.2500000000000001*(alphaL[25]*fUpL[37]+fUpL[25]*alphaL[37])+0.223606797749979*(alphaL[15]*fUpL[34]+fUpL[15]*alphaL[34])+0.2500000000000001*(alphaL[20]*fUpL[33]+fUpL[20]*alphaL[33]+alphaL[19]*fUpL[32]+fUpL[19]*alphaL[32])+0.223606797749979*(alphaL[10]*fUpL[27]+fUpL[10]*alphaL[27]+alphaL[7]*fUpL[24]+fUpL[7]*alphaL[24]+alphaL[6]*fUpL[23]+fUpL[6]*alphaL[23])+0.2500000000000001*(alphaL[12]*fUpL[22]+fUpL[12]*alphaL[22]+alphaL[11]*fUpL[21]+fUpL[11]*alphaL[21])+0.25*(alphaL[8]*fUpL[17]+fUpL[8]*alphaL[17]+alphaL[5]*fUpL[15]+fUpL[5]*alphaL[15])+0.223606797749979*(alphaL[3]*fUpL[13]+fUpL[3]*alphaL[13])+0.25*(alphaL[4]*fUpL[10]+fUpL[4]*alphaL[10]+alphaL[2]*fUpL[7]+fUpL[2]*alphaL[7]+alphaL[1]*fUpL[6]+fUpL[1]*alphaL[6]+alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]); 
  GhatL[4] = 0.2500000000000001*(alphaL[34]*fUpL[46]+alphaL[33]*fUpL[45]+alphaL[32]*fUpL[44])+0.223606797749979*alphaL[17]*fUpL[42]+0.2500000000000001*(alphaL[24]*fUpL[40]+alphaL[23]*fUpL[39]+fUpL[23]*alphaL[39]+alphaL[22]*fUpL[38]+alphaL[21]*fUpL[37]+fUpL[21]*alphaL[37]+alphaL[20]*fUpL[36]+alphaL[19]*fUpL[35])+0.25*alphaL[15]*fUpL[31]+0.223606797749979*(alphaL[10]*fUpL[30]+alphaL[8]*fUpL[28])+0.2500000000000001*(alphaL[13]*fUpL[27]+fUpL[13]*alphaL[27]+alphaL[12]*fUpL[26]+alphaL[11]*fUpL[25]+fUpL[11]*alphaL[25])+0.25*(alphaL[7]*fUpL[18]+alphaL[6]*fUpL[17]+fUpL[6]*alphaL[17]+alphaL[5]*fUpL[16])+0.223606797749979*alphaL[4]*fUpL[14]+0.25*(alphaL[3]*fUpL[10]+fUpL[3]*alphaL[10]+alphaL[2]*fUpL[9]+alphaL[1]*fUpL[8]+fUpL[1]*alphaL[8]+alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]); 
  GhatL[5] = 0.25*alphaL[27]*fUpL[46]+0.223606797749979*alphaL[17]*fUpL[44]+0.25*alphaL[39]*fUpL[40]+0.223606797749979*(fUpL[31]*alphaL[37]+alphaL[8]*fUpL[35])+0.25*(alphaL[13]*fUpL[34]+fUpL[13]*alphaL[34])+(0.2*alphaL[32]+0.223606797749979*alphaL[7])*fUpL[33]+0.2*fUpL[32]*alphaL[33]+0.223606797749979*(fUpL[7]*alphaL[33]+alphaL[6]*fUpL[32]+fUpL[6]*alphaL[32])+0.25*alphaL[10]*fUpL[31]+0.223606797749979*fUpL[16]*alphaL[25]+0.25*(alphaL[23]*fUpL[24]+fUpL[23]*alphaL[24])+0.223606797749979*(alphaL[15]*fUpL[22]+fUpL[15]*alphaL[22]+alphaL[15]*fUpL[21]+fUpL[15]*alphaL[21])+(0.2*alphaL[19]+0.223606797749979*alphaL[2])*fUpL[20]+0.2*fUpL[19]*alphaL[20]+0.223606797749979*(fUpL[2]*alphaL[20]+alphaL[1]*fUpL[19]+fUpL[1]*alphaL[19])+0.25*(alphaL[17]*fUpL[18]+alphaL[4]*fUpL[16]+alphaL[3]*fUpL[15]+fUpL[3]*alphaL[15])+0.223606797749979*(alphaL[5]*fUpL[12]+fUpL[5]*alphaL[12]+alphaL[5]*fUpL[11]+fUpL[5]*alphaL[11])+0.25*(alphaL[8]*fUpL[9]+alphaL[6]*fUpL[7]+fUpL[6]*alphaL[7]+alphaL[0]*fUpL[5]+fUpL[0]*alphaL[5]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[6] = (0.2*alphaL[37]+0.223606797749979*alphaL[10])*fUpL[39]+0.2*fUpL[37]*alphaL[39]+0.223606797749979*(fUpL[10]*alphaL[39]+alphaL[8]*fUpL[37]+fUpL[8]*alphaL[37])+(0.2*alphaL[32]+0.223606797749979*alphaL[7])*fUpL[34]+(0.2*fUpL[32]+0.223606797749979*fUpL[7])*alphaL[34]+0.25*(alphaL[12]*fUpL[33]+fUpL[12]*alphaL[33])+0.223606797749979*(alphaL[5]*fUpL[32]+fUpL[5]*alphaL[32])+0.223606797749979*(alphaL[17]*fUpL[27]+fUpL[17]*alphaL[27]+alphaL[17]*fUpL[25]+fUpL[17]*alphaL[25]+alphaL[15]*fUpL[24]+fUpL[15]*alphaL[24])+(0.2*alphaL[21]+0.223606797749979*alphaL[3])*fUpL[23]+(0.2*fUpL[21]+0.223606797749979*fUpL[3])*alphaL[23]+0.25*(alphaL[20]*fUpL[22]+fUpL[20]*alphaL[22])+0.223606797749979*(alphaL[1]*fUpL[21]+fUpL[1]*alphaL[21]+alphaL[15]*fUpL[19]+fUpL[15]*alphaL[19])+0.25*(alphaL[4]*fUpL[17]+fUpL[4]*alphaL[17]+alphaL[2]*fUpL[15]+fUpL[2]*alphaL[15])+0.223606797749979*(alphaL[6]*fUpL[13]+fUpL[6]*alphaL[13]+alphaL[6]*fUpL[11]+fUpL[6]*alphaL[11])+0.25*(alphaL[8]*fUpL[10]+fUpL[8]*alphaL[10]+alphaL[5]*fUpL[7]+fUpL[5]*alphaL[7]+alphaL[0]*fUpL[6]+fUpL[0]*alphaL[6]+alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[7] = 0.223606797749979*alphaL[17]*fUpL[46]+0.25*alphaL[25]*fUpL[44]+0.223606797749979*(alphaL[10]*fUpL[40]+fUpL[31]*alphaL[39])+0.25*fUpL[35]*alphaL[37]+(0.2*alphaL[33]+0.223606797749979*alphaL[6])*fUpL[34]+0.2*fUpL[33]*alphaL[34]+0.223606797749979*(fUpL[6]*alphaL[34]+alphaL[5]*fUpL[33]+fUpL[5]*alphaL[33])+0.25*(alphaL[11]*fUpL[32]+fUpL[11]*alphaL[32]+alphaL[8]*fUpL[31])+0.223606797749979*fUpL[18]*alphaL[27]+(0.2*alphaL[22]+0.223606797749979*alphaL[3])*fUpL[24]+0.2*fUpL[22]*alphaL[24]+0.223606797749979*(fUpL[3]*alphaL[24]+alphaL[15]*fUpL[23]+fUpL[15]*alphaL[23]+alphaL[2]*fUpL[22]+fUpL[2]*alphaL[22])+0.25*(alphaL[19]*fUpL[21]+fUpL[19]*alphaL[21])+0.223606797749979*(alphaL[15]*fUpL[20]+fUpL[15]*alphaL[20])+0.25*(alphaL[4]*fUpL[18]+fUpL[16]*alphaL[17]+alphaL[1]*fUpL[15]+fUpL[1]*alphaL[15])+0.223606797749979*(alphaL[7]*fUpL[13]+fUpL[7]*alphaL[13]+alphaL[7]*fUpL[12]+fUpL[7]*alphaL[12])+0.25*(fUpL[9]*alphaL[10]+alphaL[0]*fUpL[7]+fUpL[0]*alphaL[7]+alphaL[5]*fUpL[6]+fUpL[5]*alphaL[6]+alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 
  GhatL[8] = 0.25*(alphaL[24]*fUpL[46]+alphaL[22]*fUpL[45])+0.223606797749979*alphaL[15]*fUpL[44]+(0.2*alphaL[37]+0.223606797749979*alphaL[10])*fUpL[42]+0.25*(alphaL[34]*fUpL[40]+alphaL[13]*fUpL[39]+fUpL[13]*alphaL[39]+alphaL[33]*fUpL[38])+0.223606797749979*(alphaL[6]*fUpL[37]+fUpL[6]*alphaL[37])+0.25*alphaL[12]*fUpL[36]+0.223606797749979*alphaL[5]*fUpL[35]+fUpL[31]*(0.223606797749979*alphaL[32]+0.25*alphaL[7])+0.223606797749979*alphaL[17]*fUpL[30]+(0.2*alphaL[25]+0.223606797749979*alphaL[4])*fUpL[28]+0.25*(alphaL[23]*fUpL[27]+fUpL[23]*alphaL[27]+alphaL[20]*fUpL[26])+0.223606797749979*(alphaL[1]*fUpL[25]+fUpL[1]*alphaL[25]+alphaL[17]*fUpL[21]+fUpL[17]*alphaL[21]+fUpL[16]*alphaL[19])+0.25*(alphaL[15]*fUpL[18]+alphaL[3]*fUpL[17]+fUpL[3]*alphaL[17]+alphaL[2]*fUpL[16])+0.223606797749979*(alphaL[8]*(fUpL[14]+fUpL[11])+fUpL[8]*alphaL[11])+0.25*(alphaL[6]*fUpL[10]+fUpL[6]*alphaL[10]+alphaL[5]*fUpL[9]+alphaL[0]*fUpL[8]+fUpL[0]*alphaL[8]+alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]); 
  GhatL[9] = 0.223606797749979*alphaL[17]*fUpL[47]+0.25*alphaL[23]*fUpL[46]+0.223606797749979*alphaL[15]*fUpL[45]+0.25*alphaL[21]*fUpL[44]+0.223606797749979*(alphaL[10]*fUpL[43]+alphaL[8]*fUpL[41])+0.25*(alphaL[13]*fUpL[40]+alphaL[34]*fUpL[39]+fUpL[34]*alphaL[39])+0.223606797749979*alphaL[7]*fUpL[38]+0.25*(alphaL[32]*fUpL[37]+fUpL[32]*alphaL[37])+0.223606797749979*alphaL[5]*fUpL[36]+0.25*alphaL[11]*fUpL[35]+fUpL[31]*(0.223606797749979*alphaL[33]+0.25*alphaL[6])+0.223606797749979*alphaL[4]*fUpL[29]+0.25*(alphaL[24]*fUpL[27]+fUpL[24]*alphaL[27])+0.223606797749979*alphaL[2]*fUpL[26]+0.25*(alphaL[19]*fUpL[25]+fUpL[19]*alphaL[25])+0.223606797749979*(fUpL[18]*alphaL[22]+fUpL[16]*alphaL[20])+0.25*(alphaL[3]*fUpL[18]+alphaL[15]*fUpL[17]+fUpL[15]*alphaL[17]+alphaL[1]*fUpL[16])+0.223606797749979*fUpL[9]*alphaL[12]+0.25*(alphaL[7]*fUpL[10]+fUpL[7]*alphaL[10]+alphaL[0]*fUpL[9]+alphaL[5]*fUpL[8]+fUpL[5]*alphaL[8]+alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]); 
  GhatL[10] = 0.223606797749979*alphaL[15]*fUpL[46]+0.25*(alphaL[20]*fUpL[45]+alphaL[19]*fUpL[44])+0.2*alphaL[39]*fUpL[42]+0.223606797749979*(alphaL[8]*fUpL[42]+alphaL[7]*fUpL[40]+alphaL[6]*fUpL[39]+fUpL[6]*alphaL[39])+0.25*(alphaL[12]*fUpL[38]+alphaL[11]*fUpL[37]+fUpL[11]*alphaL[37]+alphaL[33]*fUpL[36]+alphaL[32]*fUpL[35])+fUpL[31]*(0.223606797749979*alphaL[34]+0.25*alphaL[5])+0.2*alphaL[27]*fUpL[30]+0.223606797749979*(alphaL[4]*fUpL[30]+alphaL[17]*fUpL[28]+alphaL[3]*fUpL[27]+fUpL[3]*alphaL[27])+0.25*(alphaL[22]*fUpL[26]+alphaL[21]*fUpL[25]+fUpL[21]*alphaL[25])+0.223606797749979*(fUpL[18]*alphaL[24]+alphaL[17]*fUpL[23]+fUpL[17]*alphaL[23])+0.25*(alphaL[2]*fUpL[18]+alphaL[1]*fUpL[17]+fUpL[1]*alphaL[17]+alphaL[15]*fUpL[16])+0.223606797749979*(alphaL[10]*(fUpL[14]+fUpL[13])+fUpL[10]*alphaL[13])+0.25*(alphaL[0]*fUpL[10]+fUpL[0]*alphaL[10]+alphaL[7]*fUpL[9]+alphaL[6]*fUpL[8]+fUpL[6]*alphaL[8]+alphaL[3]*fUpL[4]+fUpL[3]*alphaL[4]); 
  GhatL[11] = 0.223606797749979*alphaL[39]*fUpL[39]+0.159719141249985*alphaL[37]*fUpL[37]+0.25*(alphaL[10]*fUpL[37]+fUpL[10]*alphaL[37])+0.223606797749979*(alphaL[34]*fUpL[34]+alphaL[33]*fUpL[33])+0.159719141249985*alphaL[32]*fUpL[32]+0.25*(alphaL[7]*fUpL[32]+fUpL[7]*alphaL[32])+0.159719141249985*alphaL[25]*fUpL[25]+0.2500000000000001*(alphaL[4]*fUpL[25]+fUpL[4]*alphaL[25])+0.223606797749979*alphaL[23]*fUpL[23]+0.159719141249985*alphaL[21]*fUpL[21]+0.2500000000000001*(alphaL[3]*fUpL[21]+fUpL[3]*alphaL[21])+0.223606797749979*alphaL[20]*fUpL[20]+0.159719141249985*alphaL[19]*fUpL[19]+0.2500000000000001*(alphaL[2]*fUpL[19]+fUpL[2]*alphaL[19])+0.223606797749979*(alphaL[17]*fUpL[17]+alphaL[15]*fUpL[15])+0.159719141249985*alphaL[11]*fUpL[11]+0.25*(alphaL[0]*fUpL[11]+fUpL[0]*alphaL[11])+0.223606797749979*(alphaL[8]*fUpL[8]+alphaL[6]*fUpL[6]+alphaL[5]*fUpL[5]+alphaL[1]*fUpL[1]); 
  GhatL[12] = 0.2500000000000001*alphaL[17]*fUpL[45]+0.25*(alphaL[10]*fUpL[38]+alphaL[8]*fUpL[36])+0.223606797749979*alphaL[34]*fUpL[34]+0.159719141249985*alphaL[33]*fUpL[33]+0.25*(alphaL[6]*fUpL[33]+fUpL[6]*alphaL[33])+0.223606797749979*alphaL[32]*fUpL[32]+0.2500000000000001*alphaL[4]*fUpL[26]+0.223606797749979*alphaL[24]*fUpL[24]+0.159719141249985*alphaL[22]*fUpL[22]+0.2500000000000001*(alphaL[3]*fUpL[22]+fUpL[3]*alphaL[22])+0.159719141249985*alphaL[20]*fUpL[20]+0.2500000000000001*(alphaL[1]*fUpL[20]+fUpL[1]*alphaL[20])+0.223606797749979*(alphaL[19]*fUpL[19]+alphaL[15]*fUpL[15])+0.159719141249985*alphaL[12]*fUpL[12]+0.25*(alphaL[0]*fUpL[12]+fUpL[0]*alphaL[12])+0.223606797749979*(alphaL[7]*fUpL[7]+alphaL[5]*fUpL[5]+alphaL[2]*fUpL[2]); 
  GhatL[13] = 0.159719141249985*alphaL[39]*fUpL[39]+0.25*(alphaL[8]*fUpL[39]+fUpL[8]*alphaL[39])+0.223606797749979*alphaL[37]*fUpL[37]+0.159719141249985*alphaL[34]*fUpL[34]+0.25*(alphaL[5]*fUpL[34]+fUpL[5]*alphaL[34])+0.223606797749979*(alphaL[33]*fUpL[33]+alphaL[32]*fUpL[32])+0.159719141249985*alphaL[27]*fUpL[27]+0.2500000000000001*(alphaL[4]*fUpL[27]+fUpL[4]*alphaL[27])+0.159719141249985*alphaL[24]*fUpL[24]+0.2500000000000001*(alphaL[2]*fUpL[24]+fUpL[2]*alphaL[24])+0.159719141249985*alphaL[23]*fUpL[23]+0.2500000000000001*(alphaL[1]*fUpL[23]+fUpL[1]*alphaL[23])+0.223606797749979*(alphaL[22]*fUpL[22]+alphaL[21]*fUpL[21]+alphaL[17]*fUpL[17]+alphaL[15]*fUpL[15])+0.159719141249985*alphaL[13]*fUpL[13]+0.25*(alphaL[0]*fUpL[13]+fUpL[0]*alphaL[13])+0.223606797749979*(alphaL[10]*fUpL[10]+alphaL[7]*fUpL[7]+alphaL[6]*fUpL[6]+alphaL[3]*fUpL[3]); 
  GhatL[14] = 0.2500000000000001*alphaL[15]*fUpL[47]+0.25*(alphaL[7]*fUpL[43]+alphaL[6]*fUpL[42]+alphaL[5]*fUpL[41])+0.223606797749979*(alphaL[39]*fUpL[39]+alphaL[37]*fUpL[37])+0.2500000000000001*(alphaL[3]*fUpL[30]+alphaL[2]*fUpL[29]+alphaL[1]*fUpL[28])+0.223606797749979*(alphaL[27]*fUpL[27]+alphaL[25]*fUpL[25]+alphaL[17]*fUpL[17])+0.25*alphaL[0]*fUpL[14]+0.223606797749979*(alphaL[10]*fUpL[10]+alphaL[8]*fUpL[8]+alphaL[4]*fUpL[4]); 
  GhatL[15] = (0.2*alphaL[37]+0.223606797749979*alphaL[10])*fUpL[46]+(0.2*alphaL[39]+0.223606797749979*alphaL[8])*fUpL[44]+0.223606797749979*(alphaL[17]*fUpL[40]+fUpL[18]*alphaL[39]+fUpL[16]*alphaL[37]+alphaL[17]*fUpL[35])+(0.2*(alphaL[22]+alphaL[21])+0.223606797749979*alphaL[3])*fUpL[34]+(0.2*(fUpL[22]+fUpL[21])+0.223606797749979*fUpL[3])*alphaL[34]+(0.2*(alphaL[24]+alphaL[19])+0.223606797749979*alphaL[2])*fUpL[33]+(0.2*(fUpL[24]+fUpL[19])+0.223606797749979*fUpL[2])*alphaL[33]+(0.2*(alphaL[23]+alphaL[20])+0.223606797749979*alphaL[1])*fUpL[32]+(0.2*(fUpL[23]+fUpL[20])+0.223606797749979*fUpL[1])*alphaL[32]+(0.223606797749979*(alphaL[27]+alphaL[25])+0.25*alphaL[4])*fUpL[31]+0.223606797749979*(alphaL[6]*fUpL[24]+fUpL[6]*alphaL[24]+alphaL[7]*fUpL[23]+fUpL[7]*alphaL[23]+alphaL[5]*fUpL[22]+fUpL[5]*alphaL[22]+alphaL[5]*fUpL[21]+fUpL[5]*alphaL[21]+alphaL[7]*fUpL[20]+fUpL[7]*alphaL[20]+alphaL[6]*fUpL[19]+fUpL[6]*alphaL[19])+0.25*(alphaL[8]*fUpL[18]+fUpL[9]*alphaL[17]+alphaL[10]*fUpL[16])+(0.223606797749979*(alphaL[13]+alphaL[12]+alphaL[11])+0.25*alphaL[0])*fUpL[15]+0.223606797749979*(fUpL[13]+fUpL[12]+fUpL[11])*alphaL[15]+0.25*(fUpL[0]*alphaL[15]+alphaL[1]*fUpL[7]+fUpL[1]*alphaL[7]+alphaL[2]*fUpL[6]+fUpL[2]*alphaL[6]+alphaL[3]*fUpL[5]+fUpL[3]*alphaL[5]); 
  GhatL[16] = (0.2*alphaL[37]+0.223606797749979*alphaL[10])*fUpL[47]+0.2500000000000001*alphaL[13]*fUpL[46]+(0.2*alphaL[32]+0.223606797749979*alphaL[7])*fUpL[45]+(0.2*alphaL[33]+0.223606797749979*alphaL[6])*fUpL[44]+0.223606797749979*alphaL[17]*fUpL[43]+(0.2*alphaL[25]+0.223606797749979*alphaL[4])*fUpL[41]+0.2500000000000001*(alphaL[23]*fUpL[40]+alphaL[24]*fUpL[39]+fUpL[24]*alphaL[39])+0.223606797749979*(alphaL[15]*(fUpL[38]+fUpL[37])+fUpL[15]*alphaL[37])+(0.2*alphaL[19]+0.223606797749979*alphaL[2])*fUpL[36]+(0.2*alphaL[20]+0.223606797749979*alphaL[1])*fUpL[35]+0.2500000000000001*(alphaL[27]*fUpL[34]+fUpL[27]*alphaL[34])+0.223606797749979*(fUpL[18]*alphaL[33]+alphaL[17]*fUpL[32]+fUpL[17]*alphaL[32])+(0.223606797749979*(alphaL[22]+alphaL[21])+0.25*alphaL[3])*fUpL[31]+0.223606797749979*(alphaL[8]*fUpL[29]+alphaL[5]*(fUpL[26]+fUpL[25])+fUpL[5]*alphaL[25]+fUpL[9]*alphaL[20]+alphaL[8]*fUpL[19]+fUpL[8]*alphaL[19])+0.25*(alphaL[6]*fUpL[18]+alphaL[7]*fUpL[17]+fUpL[7]*alphaL[17])+0.223606797749979*(alphaL[12]+alphaL[11])*fUpL[16]+0.25*(alphaL[0]*fUpL[16]+alphaL[10]*fUpL[15]+fUpL[10]*alphaL[15]+alphaL[1]*fUpL[9]+alphaL[2]*fUpL[8]+fUpL[2]*alphaL[8]+alphaL[4]*fUpL[5]+fUpL[4]*alphaL[5]); 
  GhatL[17] = (0.2*alphaL[32]+0.223606797749979*alphaL[7])*fUpL[46]+0.2500000000000001*alphaL[12]*fUpL[45]+(0.2*alphaL[34]+0.223606797749979*alphaL[5])*fUpL[44]+0.2*(alphaL[27]+alphaL[25])*fUpL[42]+0.223606797749979*(alphaL[4]*fUpL[42]+alphaL[15]*fUpL[40])+(0.2*alphaL[21]+0.223606797749979*alphaL[3])*fUpL[39]+(0.2*(fUpL[30]+fUpL[21])+0.223606797749979*fUpL[3])*alphaL[39]+0.2500000000000001*alphaL[20]*fUpL[38]+(0.2*alphaL[23]+0.223606797749979*alphaL[1])*fUpL[37]+(0.2*(fUpL[28]+fUpL[23])+0.223606797749979*fUpL[1])*alphaL[37]+0.2500000000000001*alphaL[22]*fUpL[36]+0.223606797749979*(alphaL[15]*fUpL[35]+fUpL[18]*alphaL[34])+0.2500000000000001*fUpL[26]*alphaL[33]+0.223606797749979*fUpL[16]*alphaL[32]+(0.223606797749979*(alphaL[24]+alphaL[19])+0.25*alphaL[2])*fUpL[31]+0.223606797749979*(alphaL[8]*fUpL[30]+alphaL[10]*fUpL[28]+alphaL[6]*fUpL[27]+fUpL[6]*alphaL[27]+alphaL[6]*fUpL[25]+fUpL[6]*alphaL[25]+alphaL[10]*fUpL[23]+fUpL[10]*alphaL[23]+alphaL[8]*fUpL[21]+fUpL[8]*alphaL[21])+0.25*alphaL[5]*fUpL[18]+(0.223606797749979*(alphaL[13]+alphaL[11])+0.25*alphaL[0])*fUpL[17]+0.223606797749979*(fUpL[14]+fUpL[13]+fUpL[11])*alphaL[17]+0.25*(fUpL[0]*alphaL[17]+alphaL[7]*fUpL[16]+fUpL[9]*alphaL[15]+alphaL[1]*fUpL[10]+fUpL[1]*alphaL[10]+alphaL[3]*fUpL[8]+fUpL[3]*alphaL[8]+alphaL[4]*fUpL[6]+fUpL[4]*alphaL[6]); 
  GhatL[18] = (0.2*alphaL[39]+0.223606797749979*alphaL[8])*fUpL[47]+(0.2*alphaL[33]+0.223606797749979*alphaL[6])*fUpL[46]+(0.2*alphaL[34]+0.223606797749979*alphaL[5])*fUpL[45]+0.2500000000000001*alphaL[11]*fUpL[44]+0.2*alphaL[27]*fUpL[43]+0.223606797749979*(alphaL[4]*fUpL[43]+alphaL[17]*fUpL[41])+0.2*alphaL[22]*fUpL[40]+0.223606797749979*(alphaL[3]*fUpL[40]+alphaL[15]*fUpL[39]+fUpL[15]*alphaL[39])+(0.2*alphaL[24]+0.223606797749979*alphaL[2])*fUpL[38]+0.2500000000000001*(alphaL[19]*fUpL[37]+fUpL[19]*alphaL[37])+0.223606797749979*alphaL[15]*fUpL[36]+0.2500000000000001*alphaL[21]*fUpL[35]+0.223606797749979*(alphaL[17]*fUpL[34]+fUpL[17]*alphaL[34]+fUpL[16]*alphaL[33])+0.2500000000000001*(alphaL[25]*fUpL[32]+fUpL[25]*alphaL[32])+(0.223606797749979*(alphaL[23]+alphaL[20])+0.25*alphaL[1])*fUpL[31]+0.223606797749979*(alphaL[10]*fUpL[29]+alphaL[7]*fUpL[27]+fUpL[7]*alphaL[27]+alphaL[7]*fUpL[26]+alphaL[10]*fUpL[24]+fUpL[10]*alphaL[24]+fUpL[9]*alphaL[22])+0.223606797749979*(alphaL[13]+alphaL[12])*fUpL[18]+0.25*(alphaL[0]*fUpL[18]+alphaL[5]*fUpL[17]+fUpL[5]*alphaL[17]+alphaL[6]*fUpL[16]+alphaL[8]*fUpL[15]+fUpL[8]*alphaL[15]+alphaL[2]*fUpL[10]+fUpL[2]*alphaL[10]+alphaL[3]*fUpL[9]+alphaL[4]*fUpL[7]+fUpL[4]*alphaL[7]); 
  GhatL[19] = 0.223606797749979*alphaL[39]*fUpL[46]+(0.159719141249985*alphaL[37]+0.25*alphaL[10])*fUpL[44]+0.2500000000000001*fUpL[18]*alphaL[37]+(0.159719141249985*alphaL[25]+0.2500000000000001*alphaL[4])*fUpL[35]+0.223606797749979*(alphaL[23]*fUpL[34]+fUpL[23]*alphaL[34])+0.2*(alphaL[15]*fUpL[33]+fUpL[15]*alphaL[33])+(0.223606797749979*alphaL[22]+0.159719141249985*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[32]+(0.223606797749979*fUpL[22]+0.159719141249985*fUpL[21]+0.2500000000000001*fUpL[3])*alphaL[32]+0.223606797749979*alphaL[17]*fUpL[31]+0.25*(fUpL[9]*alphaL[25]+alphaL[7]*fUpL[21]+fUpL[7]*alphaL[21])+0.2*(alphaL[5]*fUpL[20]+fUpL[5]*alphaL[20])+(0.223606797749979*alphaL[12]+0.159719141249985*alphaL[11]+0.25*alphaL[0])*fUpL[19]+(0.223606797749979*fUpL[12]+0.159719141249985*fUpL[11]+0.25*fUpL[0])*alphaL[19]+0.223606797749979*(alphaL[8]*fUpL[16]+alphaL[6]*fUpL[15]+fUpL[6]*alphaL[15])+0.2500000000000001*(alphaL[2]*fUpL[11]+fUpL[2]*alphaL[11])+0.223606797749979*(alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]); 
  GhatL[20] = (0.223606797749979*alphaL[37]+0.25*alphaL[10])*fUpL[45]+0.2500000000000001*alphaL[17]*fUpL[38]+(0.223606797749979*alphaL[25]+0.2500000000000001*alphaL[4])*fUpL[36]+0.223606797749979*(alphaL[24]*fUpL[34]+fUpL[24]*alphaL[34])+(0.159719141249985*alphaL[22]+0.223606797749979*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[33]+(0.159719141249985*fUpL[22]+0.223606797749979*fUpL[21]+0.2500000000000001*fUpL[3])*alphaL[33]+0.2*(alphaL[15]*fUpL[32]+fUpL[15]*alphaL[32])+0.25*(alphaL[8]*fUpL[26]+alphaL[6]*fUpL[22]+fUpL[6]*alphaL[22])+(0.159719141249985*alphaL[12]+0.223606797749979*alphaL[11]+0.25*alphaL[0])*fUpL[20]+(0.159719141249985*fUpL[12]+0.223606797749979*fUpL[11]+0.25*fUpL[0])*alphaL[20]+0.2*(alphaL[5]*fUpL[19]+fUpL[5]*alphaL[19])+0.223606797749979*(alphaL[7]*fUpL[15]+fUpL[7]*alphaL[15])+0.2500000000000001*(alphaL[1]*fUpL[12]+fUpL[1]*alphaL[12])+0.223606797749979*(alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5]); 
  GhatL[21] = 0.2*(alphaL[17]*fUpL[39]+fUpL[17]*alphaL[39])+(0.223606797749979*alphaL[27]+0.159719141249985*alphaL[25]+0.2500000000000001*alphaL[4])*fUpL[37]+(0.223606797749979*fUpL[27]+0.159719141249985*fUpL[25]+0.2500000000000001*fUpL[4])*alphaL[37]+0.2*(alphaL[15]*fUpL[34]+fUpL[15]*alphaL[34])+0.223606797749979*(alphaL[20]*fUpL[33]+fUpL[20]*alphaL[33])+(0.223606797749979*alphaL[24]+0.159719141249985*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[32]+(0.223606797749979*fUpL[24]+0.159719141249985*fUpL[19]+0.2500000000000001*fUpL[2])*alphaL[32]+0.25*(alphaL[10]*fUpL[25]+fUpL[10]*alphaL[25])+0.2*(alphaL[6]*fUpL[23]+fUpL[6]*alphaL[23])+(0.223606797749979*alphaL[13]+0.159719141249985*alphaL[11]+0.25*alphaL[0])*fUpL[21]+(0.223606797749979*fUpL[13]+0.159719141249985*fUpL[11])*alphaL[21]+0.25*(fUpL[0]*alphaL[21]+alphaL[7]*fUpL[19]+fUpL[7]*alphaL[19])+0.223606797749979*(alphaL[8]*fUpL[17]+fUpL[8]*alphaL[17]+alphaL[5]*fUpL[15]+fUpL[5]*alphaL[15])+0.2500000000000001*(alphaL[3]*fUpL[11]+fUpL[3]*alphaL[11])+0.223606797749979*(alphaL[1]*fUpL[6]+fUpL[1]*alphaL[6]); 
  GhatL[22] = (0.223606797749979*alphaL[39]+0.25*alphaL[8])*fUpL[45]+0.223606797749979*alphaL[27]*fUpL[38]+0.2500000000000001*(alphaL[4]*fUpL[38]+alphaL[17]*fUpL[36])+0.2*(alphaL[15]*fUpL[34]+fUpL[15]*alphaL[34])+(0.223606797749979*alphaL[23]+0.159719141249985*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[33]+(0.223606797749979*fUpL[23]+0.159719141249985*fUpL[20]+0.2500000000000001*fUpL[1])*alphaL[33]+0.223606797749979*(alphaL[19]*fUpL[32]+fUpL[19]*alphaL[32])+0.25*alphaL[10]*fUpL[26]+0.2*(alphaL[7]*fUpL[24]+fUpL[7]*alphaL[24])+(0.223606797749979*alphaL[13]+0.159719141249985*alphaL[12]+0.25*alphaL[0])*fUpL[22]+(0.223606797749979*fUpL[13]+0.159719141249985*fUpL[12])*alphaL[22]+0.25*(fUpL[0]*alphaL[22]+alphaL[6]*fUpL[20]+fUpL[6]*alphaL[20])+0.223606797749979*(alphaL[5]*fUpL[15]+fUpL[5]*alphaL[15])+0.2500000000000001*(alphaL[3]*fUpL[12]+fUpL[3]*alphaL[12])+0.223606797749979*(alphaL[2]*fUpL[7]+fUpL[2]*alphaL[7]); 
  GhatL[23] = (0.159719141249985*alphaL[27]+0.223606797749979*alphaL[25]+0.2500000000000001*alphaL[4])*fUpL[39]+(0.159719141249985*fUpL[27]+0.223606797749979*fUpL[25]+0.2500000000000001*fUpL[4])*alphaL[39]+0.2*(alphaL[17]*fUpL[37]+fUpL[17]*alphaL[37])+(0.159719141249985*alphaL[24]+0.223606797749979*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[34]+(0.159719141249985*fUpL[24]+0.223606797749979*fUpL[19]+0.2500000000000001*fUpL[2])*alphaL[34]+0.223606797749979*(alphaL[22]*fUpL[33]+fUpL[22]*alphaL[33])+0.2*(alphaL[15]*fUpL[32]+fUpL[15]*alphaL[32])+0.25*(alphaL[8]*fUpL[27]+fUpL[8]*alphaL[27]+alphaL[5]*fUpL[24]+fUpL[5]*alphaL[24])+(0.159719141249985*alphaL[13]+0.223606797749979*alphaL[11]+0.25*alphaL[0])*fUpL[23]+(0.159719141249985*fUpL[13]+0.223606797749979*fUpL[11]+0.25*fUpL[0])*alphaL[23]+0.2*(alphaL[6]*fUpL[21]+fUpL[6]*alphaL[21])+0.223606797749979*(alphaL[10]*fUpL[17]+fUpL[10]*alphaL[17]+alphaL[7]*fUpL[15]+fUpL[7]*alphaL[15])+0.2500000000000001*(alphaL[1]*fUpL[13]+fUpL[1]*alphaL[13])+0.223606797749979*(alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6]); 
  GhatL[24] = (0.159719141249985*alphaL[39]+0.25*alphaL[8])*fUpL[46]+0.223606797749979*alphaL[37]*fUpL[44]+0.159719141249985*alphaL[27]*fUpL[40]+0.2500000000000001*(alphaL[4]*fUpL[40]+fUpL[16]*alphaL[39])+(0.159719141249985*alphaL[23]+0.223606797749979*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[34]+(0.159719141249985*fUpL[23]+0.223606797749979*fUpL[20]+0.2500000000000001*fUpL[1])*alphaL[34]+0.2*(alphaL[15]*fUpL[33]+fUpL[15]*alphaL[33])+0.223606797749979*(alphaL[21]*fUpL[32]+fUpL[21]*alphaL[32])+0.223606797749979*alphaL[17]*fUpL[31]+0.25*fUpL[9]*alphaL[27]+(0.159719141249985*alphaL[13]+0.223606797749979*alphaL[12]+0.25*alphaL[0])*fUpL[24]+(0.159719141249985*fUpL[13]+0.223606797749979*fUpL[12])*alphaL[24]+0.25*(fUpL[0]*alphaL[24]+alphaL[5]*fUpL[23]+fUpL[5]*alphaL[23])+0.2*(alphaL[7]*fUpL[22]+fUpL[7]*alphaL[22])+0.223606797749979*(alphaL[10]*fUpL[18]+alphaL[6]*fUpL[15]+fUpL[6]*alphaL[15])+0.2500000000000001*(alphaL[2]*fUpL[13]+fUpL[2]*alphaL[13])+0.223606797749979*(alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]); 
  GhatL[25] = 0.223606797749979*(alphaL[34]*fUpL[46]+alphaL[33]*fUpL[45])+(0.159719141249985*alphaL[32]+0.25*alphaL[7])*fUpL[44]+0.2*alphaL[17]*fUpL[42]+0.223606797749979*(alphaL[23]*fUpL[39]+fUpL[23]*alphaL[39])+(0.159719141249985*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[37]+(0.223606797749979*fUpL[30]+0.159719141249985*fUpL[21]+0.2500000000000001*fUpL[3])*alphaL[37]+0.223606797749979*alphaL[20]*fUpL[36]+0.159719141249985*alphaL[19]*fUpL[35]+0.2500000000000001*(alphaL[2]*fUpL[35]+fUpL[18]*alphaL[32])+0.223606797749979*alphaL[15]*fUpL[31]+0.2*alphaL[8]*fUpL[28]+(0.159719141249985*alphaL[11]+0.25*alphaL[0])*fUpL[25]+(0.223606797749979*fUpL[14]+0.159719141249985*fUpL[11])*alphaL[25]+0.25*(fUpL[0]*alphaL[25]+alphaL[10]*fUpL[21]+fUpL[10]*alphaL[21]+fUpL[9]*alphaL[19])+0.223606797749979*(alphaL[6]*fUpL[17]+fUpL[6]*alphaL[17]+alphaL[5]*fUpL[16])+0.2500000000000001*(alphaL[4]*fUpL[11]+fUpL[4]*alphaL[11])+0.223606797749979*(alphaL[1]*fUpL[8]+fUpL[1]*alphaL[8]); 
  GhatL[26] = 0.223606797749979*alphaL[34]*fUpL[46]+(0.159719141249985*alphaL[33]+0.25*alphaL[6])*fUpL[45]+0.223606797749979*(alphaL[32]*fUpL[44]+alphaL[24]*fUpL[40])+(0.159719141249985*alphaL[22]+0.2500000000000001*alphaL[3])*fUpL[38]+(0.159719141249985*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[36]+0.223606797749979*alphaL[19]*fUpL[35]+0.2500000000000001*(alphaL[17]*fUpL[33]+fUpL[17]*alphaL[33])+0.223606797749979*alphaL[15]*fUpL[31]+0.159719141249985*alphaL[12]*fUpL[26]+0.25*(alphaL[0]*fUpL[26]+alphaL[10]*fUpL[22]+fUpL[10]*alphaL[22]+alphaL[8]*fUpL[20]+fUpL[8]*alphaL[20])+0.223606797749979*(alphaL[7]*fUpL[18]+alphaL[5]*fUpL[16])+0.2500000000000001*(alphaL[4]*fUpL[12]+fUpL[4]*alphaL[12])+0.223606797749979*alphaL[2]*fUpL[9]; 
  GhatL[27] = (0.159719141249985*alphaL[34]+0.25*alphaL[5])*fUpL[46]+0.223606797749979*(alphaL[33]*fUpL[45]+alphaL[32]*fUpL[44])+0.2*alphaL[17]*fUpL[42]+(0.159719141249985*alphaL[24]+0.2500000000000001*alphaL[2])*fUpL[40]+(0.159719141249985*alphaL[23]+0.2500000000000001*alphaL[1])*fUpL[39]+(0.223606797749979*fUpL[28]+0.159719141249985*fUpL[23]+0.2500000000000001*fUpL[1])*alphaL[39]+0.223606797749979*(alphaL[22]*fUpL[38]+alphaL[21]*fUpL[37]+fUpL[21]*alphaL[37])+0.2500000000000001*fUpL[16]*alphaL[34]+0.223606797749979*alphaL[15]*fUpL[31]+0.2*alphaL[10]*fUpL[30]+(0.159719141249985*alphaL[13]+0.25*alphaL[0])*fUpL[27]+(0.223606797749979*fUpL[14]+0.159719141249985*fUpL[13])*alphaL[27]+0.25*(fUpL[0]*alphaL[27]+fUpL[9]*alphaL[24]+alphaL[8]*fUpL[23]+fUpL[8]*alphaL[23])+0.223606797749979*(alphaL[7]*fUpL[18]+alphaL[6]*fUpL[17]+fUpL[6]*alphaL[17])+0.2500000000000001*(alphaL[4]*fUpL[13]+fUpL[4]*alphaL[13])+0.223606797749979*(alphaL[3]*fUpL[10]+fUpL[3]*alphaL[10]); 
  GhatL[28] = (0.223606797749979*alphaL[32]+0.25*alphaL[7])*fUpL[47]+0.2500000000000001*alphaL[15]*fUpL[43]+(0.223606797749979*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[42]+(0.223606797749979*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[41]+0.223606797749979*(alphaL[27]*fUpL[39]+fUpL[27]*alphaL[39])+0.2*(alphaL[17]*fUpL[37]+fUpL[17]*alphaL[37])+0.25*(alphaL[6]*fUpL[30]+alphaL[5]*fUpL[29])+(0.223606797749979*alphaL[11]+0.25*alphaL[0])*fUpL[28]+0.2*(alphaL[8]*fUpL[25]+fUpL[8]*alphaL[25])+0.223606797749979*(alphaL[10]*fUpL[17]+fUpL[10]*alphaL[17])+0.2500000000000001*alphaL[1]*fUpL[14]+0.223606797749979*(alphaL[4]*fUpL[8]+fUpL[4]*alphaL[8]); 
  GhatL[29] = (0.223606797749979*alphaL[33]+0.25*alphaL[6])*fUpL[47]+0.223606797749979*(alphaL[39]*fUpL[46]+alphaL[37]*fUpL[44]+alphaL[22]*fUpL[43])+0.2500000000000001*(alphaL[3]*fUpL[43]+alphaL[15]*fUpL[42])+(0.223606797749979*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[41]+0.223606797749979*(alphaL[27]*fUpL[40]+alphaL[25]*fUpL[35])+0.223606797749979*alphaL[17]*fUpL[31]+0.25*alphaL[7]*fUpL[30]+0.223606797749979*alphaL[12]*fUpL[29]+0.25*(alphaL[0]*fUpL[29]+alphaL[5]*fUpL[28])+0.223606797749979*(alphaL[10]*fUpL[18]+alphaL[8]*fUpL[16])+0.2500000000000001*alphaL[2]*fUpL[14]+0.223606797749979*alphaL[4]*fUpL[9]; 
  GhatL[30] = (0.223606797749979*alphaL[34]+0.25*alphaL[5])*fUpL[47]+(0.223606797749979*alphaL[24]+0.2500000000000001*alphaL[2])*fUpL[43]+0.223606797749979*alphaL[23]*fUpL[42]+0.2500000000000001*(alphaL[1]*fUpL[42]+alphaL[15]*fUpL[41])+0.2*(alphaL[17]*fUpL[39]+fUpL[17]*alphaL[39])+0.223606797749979*(alphaL[25]*fUpL[37]+fUpL[25]*alphaL[37]+alphaL[13]*fUpL[30])+0.25*(alphaL[0]*fUpL[30]+alphaL[7]*fUpL[29]+alphaL[6]*fUpL[28])+0.2*(alphaL[10]*fUpL[27]+fUpL[10]*alphaL[27])+0.223606797749979*(alphaL[8]*fUpL[17]+fUpL[8]*alphaL[17])+0.2500000000000001*alphaL[3]*fUpL[14]+0.223606797749979*(alphaL[4]*fUpL[10]+fUpL[4]*alphaL[10]); 
  GhatL[31] = (0.2*(alphaL[27]+alphaL[25])+0.223606797749979*alphaL[4])*fUpL[47]+(0.2*(alphaL[22]+alphaL[21])+0.223606797749979*alphaL[3])*fUpL[46]+(0.2*(alphaL[24]+alphaL[19])+0.223606797749979*alphaL[2])*fUpL[45]+(0.2*(alphaL[23]+alphaL[20])+0.223606797749979*alphaL[1])*fUpL[44]+(0.2*alphaL[39]+0.223606797749979*alphaL[8])*fUpL[43]+(0.2*alphaL[37]+0.223606797749979*alphaL[10])*fUpL[41]+(0.2*alphaL[33]+0.223606797749979*alphaL[6])*fUpL[40]+(0.2*alphaL[32]+0.223606797749979*alphaL[7])*fUpL[39]+(0.2*fUpL[32]+0.223606797749979*fUpL[7])*alphaL[39]+(0.2*alphaL[34]+0.223606797749979*alphaL[5])*(fUpL[38]+fUpL[37])+(0.2*fUpL[34]+0.223606797749979*fUpL[5])*alphaL[37]+(0.2*alphaL[32]+0.223606797749979*alphaL[7])*fUpL[36]+0.2*alphaL[33]*fUpL[35]+0.223606797749979*(alphaL[6]*fUpL[35]+alphaL[10]*fUpL[34]+fUpL[10]*alphaL[34]+fUpL[9]*alphaL[33]+alphaL[8]*fUpL[32]+fUpL[8]*alphaL[32])+(0.223606797749979*(alphaL[13]+alphaL[12]+alphaL[11])+0.25*alphaL[0])*fUpL[31]+0.223606797749979*(alphaL[17]*fUpL[29]+alphaL[15]*fUpL[27]+fUpL[15]*alphaL[27]+alphaL[15]*(fUpL[26]+fUpL[25])+fUpL[15]*alphaL[25]+alphaL[17]*fUpL[24]+fUpL[17]*alphaL[24]+fUpL[18]*alphaL[23]+fUpL[16]*(alphaL[22]+alphaL[21])+fUpL[18]*alphaL[20]+alphaL[17]*fUpL[19]+fUpL[17]*alphaL[19])+0.25*(alphaL[1]*fUpL[18]+alphaL[2]*fUpL[17]+fUpL[2]*alphaL[17]+alphaL[3]*fUpL[16]+alphaL[4]*fUpL[15]+fUpL[4]*alphaL[15]+alphaL[5]*fUpL[10]+fUpL[5]*alphaL[10]+alphaL[6]*fUpL[9]+alphaL[7]*fUpL[8]+fUpL[7]*alphaL[8]); 
  GhatL[32] = 0.2*alphaL[17]*fUpL[46]+(0.223606797749979*alphaL[27]+0.159719141249985*alphaL[25]+0.2500000000000001*alphaL[4])*fUpL[44]+0.223606797749979*alphaL[37]*fUpL[40]+0.2*fUpL[31]*alphaL[39]+0.159719141249985*fUpL[35]*alphaL[37]+0.25*(fUpL[9]*alphaL[37]+alphaL[10]*fUpL[35])+(0.1788854381999831*alphaL[33]+0.2*alphaL[6])*fUpL[34]+0.1788854381999831*fUpL[33]*alphaL[34]+0.2*(fUpL[6]*alphaL[34]+alphaL[5]*fUpL[33]+fUpL[5]*alphaL[33])+(0.223606797749979*(alphaL[13]+alphaL[12])+0.159719141249985*alphaL[11]+0.25*alphaL[0])*fUpL[32]+(0.223606797749979*(fUpL[13]+fUpL[12])+0.159719141249985*fUpL[11]+0.25*fUpL[0])*alphaL[32]+0.223606797749979*alphaL[8]*fUpL[31]+0.2500000000000001*fUpL[18]*alphaL[25]+0.223606797749979*(alphaL[21]*fUpL[24]+fUpL[21]*alphaL[24])+0.2*(alphaL[15]*fUpL[23]+fUpL[15]*alphaL[23])+0.223606797749979*(alphaL[19]*fUpL[22]+fUpL[19]*alphaL[22])+(0.159719141249985*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[21]+(0.159719141249985*fUpL[19]+0.2500000000000001*fUpL[2])*alphaL[21]+0.2*(alphaL[15]*fUpL[20]+fUpL[15]*alphaL[20])+0.2500000000000001*(alphaL[3]*fUpL[19]+fUpL[3]*alphaL[19])+0.223606797749979*(fUpL[16]*alphaL[17]+alphaL[1]*fUpL[15]+fUpL[1]*alphaL[15])+0.25*(alphaL[7]*fUpL[11]+fUpL[7]*alphaL[11])+0.223606797749979*(alphaL[5]*fUpL[6]+fUpL[5]*alphaL[6]); 
  GhatL[33] = (0.223606797749979*(alphaL[27]+alphaL[25])+0.2500000000000001*alphaL[4])*fUpL[45]+fUpL[38]*(0.223606797749979*alphaL[39]+0.25*alphaL[8])+fUpL[36]*(0.223606797749979*alphaL[37]+0.25*alphaL[10])+(0.1788854381999831*alphaL[32]+0.2*alphaL[7])*fUpL[34]+(0.1788854381999831*fUpL[32]+0.2*fUpL[7])*alphaL[34]+(0.223606797749979*alphaL[13]+0.159719141249985*alphaL[12]+0.223606797749979*alphaL[11]+0.25*alphaL[0])*fUpL[33]+(0.223606797749979*fUpL[13]+0.159719141249985*fUpL[12]+0.223606797749979*fUpL[11]+0.25*fUpL[0])*alphaL[33]+0.2*(alphaL[5]*fUpL[32]+fUpL[5]*alphaL[32])+0.2500000000000001*alphaL[17]*fUpL[26]+0.2*(alphaL[15]*fUpL[24]+fUpL[15]*alphaL[24])+0.223606797749979*alphaL[22]*fUpL[23]+fUpL[22]*(0.223606797749979*alphaL[23]+0.159719141249985*alphaL[20]+0.2500000000000001*alphaL[1])+(0.159719141249985*fUpL[20]+0.2500000000000001*fUpL[1])*alphaL[22]+0.223606797749979*(alphaL[20]*fUpL[21]+fUpL[20]*alphaL[21])+0.2500000000000001*(alphaL[3]*fUpL[20]+fUpL[3]*alphaL[20])+0.2*(alphaL[15]*fUpL[19]+fUpL[15]*alphaL[19])+0.223606797749979*(alphaL[2]*fUpL[15]+fUpL[2]*alphaL[15])+0.25*(alphaL[6]*fUpL[12]+fUpL[6]*alphaL[12])+0.223606797749979*(alphaL[5]*fUpL[7]+fUpL[5]*alphaL[7]); 
  GhatL[34] = (0.159719141249985*alphaL[27]+0.223606797749979*alphaL[25]+0.2500000000000001*alphaL[4])*fUpL[46]+0.2*alphaL[17]*fUpL[44]+(0.159719141249985*alphaL[39]+0.25*alphaL[8])*fUpL[40]+(0.223606797749979*fUpL[35]+0.25*fUpL[9])*alphaL[39]+0.2*fUpL[31]*alphaL[37]+(0.159719141249985*alphaL[13]+0.223606797749979*(alphaL[12]+alphaL[11])+0.25*alphaL[0])*fUpL[34]+(0.159719141249985*fUpL[13]+0.223606797749979*(fUpL[12]+fUpL[11])+0.25*fUpL[0])*alphaL[34]+(0.1788854381999831*alphaL[32]+0.2*alphaL[7])*fUpL[33]+0.1788854381999831*fUpL[32]*alphaL[33]+0.2*(fUpL[7]*alphaL[33]+alphaL[6]*fUpL[32]+fUpL[6]*alphaL[32])+0.223606797749979*alphaL[10]*fUpL[31]+0.2500000000000001*fUpL[16]*alphaL[27]+(0.159719141249985*alphaL[23]+0.223606797749979*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[24]+(0.159719141249985*fUpL[23]+0.223606797749979*fUpL[20]+0.2500000000000001*fUpL[1])*alphaL[24]+(0.223606797749979*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[23]+(0.223606797749979*fUpL[19]+0.2500000000000001*fUpL[2])*alphaL[23]+0.2*(alphaL[15]*fUpL[22]+fUpL[15]*alphaL[22]+alphaL[15]*fUpL[21]+fUpL[15]*alphaL[21])+0.223606797749979*(alphaL[17]*fUpL[18]+alphaL[3]*fUpL[15]+fUpL[3]*alphaL[15])+0.25*(alphaL[5]*fUpL[13]+fUpL[5]*alphaL[13])+0.223606797749979*(alphaL[6]*fUpL[7]+fUpL[6]*alphaL[7]); 
  GhatL[35] = 0.2*alphaL[17]*fUpL[47]+0.223606797749979*alphaL[23]*fUpL[46]+0.2*alphaL[15]*fUpL[45]+(0.223606797749979*alphaL[22]+0.159719141249985*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[44]+0.223606797749979*alphaL[37]*fUpL[43]+0.2*alphaL[8]*fUpL[41]+0.223606797749979*(alphaL[34]*fUpL[39]+fUpL[34]*alphaL[39]+alphaL[32]*fUpL[38])+(0.159719141249985*alphaL[32]+0.25*alphaL[7])*fUpL[37]+(0.159719141249985*fUpL[32]+0.25*fUpL[7])*alphaL[37]+0.2*alphaL[5]*fUpL[36]+(0.223606797749979*alphaL[12]+0.159719141249985*alphaL[11]+0.25*alphaL[0])*fUpL[35]+0.2*fUpL[31]*alphaL[33]+0.25*(alphaL[10]*fUpL[32]+fUpL[10]*alphaL[32])+0.223606797749979*(alphaL[6]*fUpL[31]+alphaL[25]*fUpL[29]+alphaL[19]*fUpL[26])+(0.159719141249985*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[25]+0.159719141249985*fUpL[19]*alphaL[25]+0.2500000000000001*(fUpL[2]*alphaL[25]+fUpL[18]*alphaL[21])+0.2*fUpL[16]*alphaL[20]+0.2500000000000001*(alphaL[4]*fUpL[19]+fUpL[4]*alphaL[19])+0.223606797749979*(alphaL[15]*fUpL[17]+fUpL[15]*alphaL[17]+alphaL[1]*fUpL[16])+0.25*fUpL[9]*alphaL[11]+0.223606797749979*(alphaL[5]*fUpL[8]+fUpL[5]*alphaL[8]); 
  GhatL[36] = 0.223606797749979*alphaL[24]*fUpL[46]+(0.159719141249985*alphaL[22]+0.223606797749979*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[45]+0.2*alphaL[15]*fUpL[44]+0.223606797749979*alphaL[34]*fUpL[40]+(0.159719141249985*alphaL[33]+0.25*alphaL[6])*fUpL[38]+0.223606797749979*(alphaL[33]*fUpL[37]+fUpL[33]*alphaL[37])+(0.159719141249985*alphaL[12]+0.223606797749979*alphaL[11]+0.25*alphaL[0])*fUpL[36]+0.2*alphaL[5]*fUpL[35]+0.25*(alphaL[10]*fUpL[33]+fUpL[10]*alphaL[33])+fUpL[31]*(0.2*alphaL[32]+0.223606797749979*alphaL[7])+(0.159719141249985*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[26]+0.223606797749979*(alphaL[20]*fUpL[25]+fUpL[20]*alphaL[25])+0.2500000000000001*(alphaL[17]*fUpL[22]+fUpL[17]*alphaL[22]+alphaL[4]*fUpL[20]+fUpL[4]*alphaL[20])+0.2*fUpL[16]*alphaL[19]+0.223606797749979*(alphaL[15]*fUpL[18]+alphaL[2]*fUpL[16])+0.25*(alphaL[8]*fUpL[12]+fUpL[8]*alphaL[12])+0.223606797749979*alphaL[5]*fUpL[9]; 
  GhatL[37] = 0.2*alphaL[15]*fUpL[46]+0.223606797749979*alphaL[20]*fUpL[45]+(0.223606797749979*alphaL[24]+0.159719141249985*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[44]+(0.1788854381999831*alphaL[39]+0.2*alphaL[8])*fUpL[42]+0.223606797749979*alphaL[32]*fUpL[40]+0.2*(alphaL[6]*fUpL[39]+fUpL[6]*alphaL[39])+(0.223606797749979*alphaL[13]+0.159719141249985*alphaL[11]+0.25*alphaL[0])*fUpL[37]+(0.223606797749979*(fUpL[14]+fUpL[13])+0.159719141249985*fUpL[11]+0.25*fUpL[0])*alphaL[37]+0.223606797749979*alphaL[33]*fUpL[36]+(0.159719141249985*alphaL[32]+0.25*alphaL[7])*fUpL[35]+0.2*fUpL[31]*alphaL[34]+0.25*fUpL[9]*alphaL[32]+0.223606797749979*(alphaL[5]*fUpL[31]+alphaL[25]*fUpL[30])+0.2*alphaL[17]*fUpL[28]+0.223606797749979*(alphaL[21]*fUpL[27]+fUpL[21]*alphaL[27])+(0.159719141249985*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[25]+(0.159719141249985*fUpL[21]+0.2500000000000001*fUpL[3])*alphaL[25]+0.2*(alphaL[17]*fUpL[23]+fUpL[17]*alphaL[23])+0.2500000000000001*(alphaL[4]*fUpL[21]+fUpL[4]*alphaL[21]+fUpL[18]*alphaL[19])+0.223606797749979*(alphaL[1]*fUpL[17]+fUpL[1]*alphaL[17]+alphaL[15]*fUpL[16])+0.25*(alphaL[10]*fUpL[11]+fUpL[10]*alphaL[11])+0.223606797749979*(alphaL[6]*fUpL[8]+fUpL[6]*alphaL[8]); 
  GhatL[38] = 0.2*alphaL[15]*fUpL[46]+(0.223606797749979*alphaL[23]+0.159719141249985*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[45]+0.223606797749979*alphaL[19]*fUpL[44]+0.2*alphaL[7]*fUpL[40]+0.223606797749979*(alphaL[33]*fUpL[39]+fUpL[33]*alphaL[39])+(0.223606797749979*alphaL[13]+0.159719141249985*alphaL[12]+0.25*alphaL[0])*fUpL[38]+(0.159719141249985*alphaL[33]+0.25*alphaL[6])*fUpL[36]+0.223606797749979*alphaL[32]*fUpL[35]+0.2*fUpL[31]*alphaL[34]+0.25*(alphaL[8]*fUpL[33]+fUpL[8]*alphaL[33])+0.223606797749979*(alphaL[5]*fUpL[31]+alphaL[22]*fUpL[27]+fUpL[22]*alphaL[27])+(0.159719141249985*alphaL[22]+0.2500000000000001*alphaL[3])*fUpL[26]+0.2*fUpL[18]*alphaL[24]+0.2500000000000001*(alphaL[4]*fUpL[22]+fUpL[4]*alphaL[22]+alphaL[17]*fUpL[20]+fUpL[17]*alphaL[20])+0.223606797749979*(alphaL[2]*fUpL[18]+alphaL[15]*fUpL[16])+0.25*(alphaL[10]*fUpL[12]+fUpL[10]*alphaL[12])+0.223606797749979*alphaL[7]*fUpL[9]; 
  GhatL[39] = (0.159719141249985*alphaL[24]+0.223606797749979*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[46]+0.223606797749979*alphaL[22]*fUpL[45]+0.2*alphaL[15]*fUpL[44]+(0.1788854381999831*alphaL[37]+0.2*alphaL[10])*fUpL[42]+(0.159719141249985*alphaL[34]+0.25*alphaL[5])*fUpL[40]+(0.159719141249985*alphaL[13]+0.223606797749979*alphaL[11]+0.25*alphaL[0])*fUpL[39]+(0.223606797749979*fUpL[14]+0.159719141249985*fUpL[13]+0.223606797749979*fUpL[11]+0.25*fUpL[0])*alphaL[39]+0.223606797749979*alphaL[33]*fUpL[38]+0.2*(alphaL[6]*fUpL[37]+fUpL[6]*alphaL[37])+alphaL[34]*(0.223606797749979*fUpL[35]+0.25*fUpL[9])+fUpL[31]*(0.2*alphaL[32]+0.223606797749979*alphaL[7])+0.2*alphaL[17]*fUpL[30]+0.223606797749979*alphaL[27]*fUpL[28]+(0.159719141249985*alphaL[23]+0.2500000000000001*alphaL[1])*fUpL[27]+(0.159719141249985*fUpL[23]+0.2500000000000001*fUpL[1])*alphaL[27]+0.223606797749979*(alphaL[23]*fUpL[25]+fUpL[23]*alphaL[25])+0.2500000000000001*(fUpL[16]*alphaL[24]+alphaL[4]*fUpL[23]+fUpL[4]*alphaL[23])+0.2*(alphaL[17]*fUpL[21]+fUpL[17]*alphaL[21])+0.223606797749979*(alphaL[15]*fUpL[18]+alphaL[3]*fUpL[17]+fUpL[3]*alphaL[17])+0.25*(alphaL[8]*fUpL[13]+fUpL[8]*alphaL[13])+0.223606797749979*(alphaL[6]*fUpL[10]+fUpL[6]*alphaL[10]); 
  GhatL[40] = 0.2*alphaL[17]*fUpL[47]+(0.159719141249985*alphaL[23]+0.223606797749979*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[46]+0.2*alphaL[15]*fUpL[45]+0.223606797749979*alphaL[21]*fUpL[44]+0.2*alphaL[10]*fUpL[43]+0.223606797749979*alphaL[39]*fUpL[41]+(0.159719141249985*alphaL[13]+0.223606797749979*alphaL[12]+0.25*alphaL[0])*fUpL[40]+(0.159719141249985*alphaL[34]+0.25*alphaL[5])*fUpL[39]+(0.159719141249985*fUpL[34]+0.25*fUpL[5])*alphaL[39]+0.2*alphaL[7]*fUpL[38]+0.223606797749979*(alphaL[32]*fUpL[37]+fUpL[32]*alphaL[37]+alphaL[34]*fUpL[36])+0.25*(alphaL[8]*fUpL[34]+fUpL[8]*alphaL[34])+0.2*fUpL[31]*alphaL[33]+0.223606797749979*(alphaL[6]*fUpL[31]+alphaL[27]*fUpL[29])+(0.159719141249985*alphaL[24]+0.2500000000000001*alphaL[2])*fUpL[27]+(0.159719141249985*fUpL[24]+0.2500000000000001*fUpL[2])*alphaL[27]+0.223606797749979*alphaL[24]*fUpL[26]+0.2500000000000001*(alphaL[4]*fUpL[24]+fUpL[4]*alphaL[24]+fUpL[16]*alphaL[23])+0.2*fUpL[18]*alphaL[22]+0.223606797749979*(alphaL[3]*fUpL[18]+alphaL[15]*fUpL[17]+fUpL[15]*alphaL[17])+0.25*fUpL[9]*alphaL[13]+0.223606797749979*(alphaL[7]*fUpL[10]+fUpL[7]*alphaL[10]); 
  GhatL[41] = (0.223606797749979*(alphaL[22]+alphaL[21])+0.2500000000000001*alphaL[3])*fUpL[47]+0.223606797749979*alphaL[27]*fUpL[46]+0.2*alphaL[17]*fUpL[44]+(0.223606797749979*alphaL[33]+0.25*alphaL[6])*fUpL[43]+(0.223606797749979*alphaL[32]+0.25*alphaL[7])*fUpL[42]+(0.223606797749979*(alphaL[12]+alphaL[11])+0.25*alphaL[0])*fUpL[41]+0.223606797749979*alphaL[39]*fUpL[40]+0.2*(fUpL[31]*alphaL[37]+alphaL[8]*fUpL[35])+0.223606797749979*alphaL[10]*fUpL[31]+0.2500000000000001*alphaL[15]*fUpL[30]+(0.223606797749979*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[29]+(0.223606797749979*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[28]+0.2*fUpL[16]*alphaL[25]+0.223606797749979*(alphaL[17]*fUpL[18]+alphaL[4]*fUpL[16])+0.25*alphaL[5]*fUpL[14]+0.223606797749979*alphaL[8]*fUpL[9]; 
  GhatL[42] = (0.223606797749979*(alphaL[24]+alphaL[19])+0.2500000000000001*alphaL[2])*fUpL[47]+(0.223606797749979*alphaL[34]+0.25*alphaL[5])*fUpL[43]+(0.223606797749979*(alphaL[13]+alphaL[11])+0.25*alphaL[0])*fUpL[42]+(0.223606797749979*alphaL[32]+0.25*alphaL[7])*fUpL[41]+(0.1788854381999831*alphaL[37]+0.2*alphaL[10])*fUpL[39]+0.1788854381999831*fUpL[37]*alphaL[39]+0.2*(fUpL[10]*alphaL[39]+alphaL[8]*fUpL[37]+fUpL[8]*alphaL[37])+0.223606797749979*alphaL[23]*fUpL[30]+0.2500000000000001*(alphaL[1]*fUpL[30]+alphaL[15]*fUpL[29])+(0.223606797749979*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[28]+0.2*(alphaL[17]*fUpL[27]+fUpL[17]*alphaL[27]+alphaL[17]*fUpL[25]+fUpL[17]*alphaL[25])+0.223606797749979*(alphaL[4]*fUpL[17]+fUpL[4]*alphaL[17])+0.25*alphaL[6]*fUpL[14]+0.223606797749979*(alphaL[8]*fUpL[10]+fUpL[8]*alphaL[10]); 
  GhatL[43] = (0.223606797749979*(alphaL[23]+alphaL[20])+0.2500000000000001*alphaL[1])*fUpL[47]+0.2*alphaL[17]*fUpL[46]+0.223606797749979*alphaL[25]*fUpL[44]+(0.223606797749979*(alphaL[13]+alphaL[12])+0.25*alphaL[0])*fUpL[43]+(0.223606797749979*alphaL[34]+0.25*alphaL[5])*fUpL[42]+(0.223606797749979*alphaL[33]+0.25*alphaL[6])*fUpL[41]+0.2*(alphaL[10]*fUpL[40]+fUpL[31]*alphaL[39])+0.223606797749979*(fUpL[35]*alphaL[37]+alphaL[8]*fUpL[31])+(0.223606797749979*alphaL[24]+0.2500000000000001*alphaL[2])*fUpL[30]+0.223606797749979*alphaL[22]*fUpL[29]+0.2500000000000001*(alphaL[3]*fUpL[29]+alphaL[15]*fUpL[28])+0.2*fUpL[18]*alphaL[27]+0.223606797749979*(alphaL[4]*fUpL[18]+fUpL[16]*alphaL[17])+0.25*alphaL[7]*fUpL[14]+0.223606797749979*fUpL[9]*alphaL[10]; 
  GhatL[44] = (0.1788854381999831*alphaL[39]+0.2*alphaL[8])*fUpL[47]+(0.1788854381999831*alphaL[33]+0.2*alphaL[6])*fUpL[46]+(0.1788854381999831*alphaL[34]+0.2*alphaL[5])*fUpL[45]+(0.223606797749979*(alphaL[13]+alphaL[12])+0.159719141249985*alphaL[11]+0.25*alphaL[0])*fUpL[44]+0.223606797749979*alphaL[25]*fUpL[43]+0.2*alphaL[17]*fUpL[41]+0.223606797749979*alphaL[21]*fUpL[40]+0.2*(alphaL[15]*fUpL[39]+fUpL[15]*alphaL[39])+0.223606797749979*alphaL[19]*fUpL[38]+(0.223606797749979*alphaL[24]+0.159719141249985*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[37]+(0.223606797749979*(fUpL[29]+fUpL[24])+0.159719141249985*fUpL[19]+0.2500000000000001*fUpL[2])*alphaL[37]+0.2*alphaL[15]*fUpL[36]+(0.223606797749979*alphaL[22]+0.159719141249985*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[35]+0.2*(alphaL[17]*fUpL[34]+fUpL[17]*alphaL[34]+fUpL[16]*alphaL[33])+(0.223606797749979*alphaL[27]+0.159719141249985*alphaL[25]+0.2500000000000001*alphaL[4])*fUpL[32]+(0.223606797749979*(fUpL[27]+fUpL[26])+0.159719141249985*fUpL[25]+0.2500000000000001*fUpL[4])*alphaL[32]+(0.2*(alphaL[23]+alphaL[20])+0.223606797749979*alphaL[1])*fUpL[31]+0.25*(alphaL[7]*fUpL[25]+fUpL[7]*alphaL[25]+fUpL[9]*alphaL[21]+alphaL[10]*fUpL[19]+fUpL[10]*alphaL[19])+0.2500000000000001*alphaL[11]*fUpL[18]+0.223606797749979*(alphaL[5]*fUpL[17]+fUpL[5]*alphaL[17]+alphaL[6]*fUpL[16]+alphaL[8]*fUpL[15]+fUpL[8]*alphaL[15]); 
  GhatL[45] = (0.1788854381999831*alphaL[32]+0.2*alphaL[7])*fUpL[46]+(0.223606797749979*alphaL[13]+0.159719141249985*alphaL[12]+0.223606797749979*alphaL[11]+0.25*alphaL[0])*fUpL[45]+0.1788854381999831*alphaL[34]*fUpL[44]+0.2*(alphaL[5]*fUpL[44]+alphaL[15]*fUpL[40])+0.223606797749979*(alphaL[22]*fUpL[39]+fUpL[22]*alphaL[39])+(0.223606797749979*alphaL[23]+0.159719141249985*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[38]+0.223606797749979*(alphaL[20]*fUpL[37]+fUpL[20]*alphaL[37])+(0.159719141249985*alphaL[22]+0.223606797749979*alphaL[21]+0.2500000000000001*alphaL[3])*fUpL[36]+0.2*(alphaL[15]*fUpL[35]+fUpL[18]*alphaL[34])+(0.223606797749979*(alphaL[27]+alphaL[25])+0.2500000000000001*alphaL[4])*fUpL[33]+(0.223606797749979*fUpL[27]+0.159719141249985*fUpL[26]+0.223606797749979*fUpL[25]+0.2500000000000001*fUpL[4])*alphaL[33]+0.2*fUpL[16]*alphaL[32]+(0.2*(alphaL[24]+alphaL[19])+0.223606797749979*alphaL[2])*fUpL[31]+0.25*(alphaL[6]*fUpL[26]+alphaL[8]*fUpL[22]+fUpL[8]*alphaL[22]+alphaL[10]*fUpL[20]+fUpL[10]*alphaL[20])+0.223606797749979*alphaL[5]*fUpL[18]+0.2500000000000001*(alphaL[12]*fUpL[17]+fUpL[12]*alphaL[17])+0.223606797749979*(alphaL[7]*fUpL[16]+fUpL[9]*alphaL[15]); 
  GhatL[46] = (0.1788854381999831*alphaL[37]+0.2*alphaL[10])*fUpL[47]+(0.159719141249985*alphaL[13]+0.223606797749979*(alphaL[12]+alphaL[11])+0.25*alphaL[0])*fUpL[46]+(0.1788854381999831*alphaL[32]+0.2*alphaL[7])*fUpL[45]+0.1788854381999831*alphaL[33]*fUpL[44]+0.2*(alphaL[6]*fUpL[44]+alphaL[17]*fUpL[43])+0.223606797749979*alphaL[27]*fUpL[41]+(0.159719141249985*alphaL[23]+0.223606797749979*alphaL[20]+0.2500000000000001*alphaL[1])*fUpL[40]+(0.159719141249985*alphaL[24]+0.223606797749979*alphaL[19]+0.2500000000000001*alphaL[2])*fUpL[39]+(0.223606797749979*fUpL[29]+0.159719141249985*fUpL[24]+0.223606797749979*fUpL[19]+0.2500000000000001*fUpL[2])*alphaL[39]+0.2*(alphaL[15]*(fUpL[38]+fUpL[37])+fUpL[15]*alphaL[37])+0.223606797749979*(alphaL[24]*fUpL[36]+alphaL[23]*fUpL[35])+(0.159719141249985*alphaL[27]+0.223606797749979*alphaL[25]+0.2500000000000001*alphaL[4])*fUpL[34]+(0.159719141249985*fUpL[27]+0.223606797749979*(fUpL[26]+fUpL[25])+0.2500000000000001*fUpL[4])*alphaL[34]+0.2*(fUpL[18]*alphaL[33]+alphaL[17]*fUpL[32]+fUpL[17]*alphaL[32])+(0.2*(alphaL[22]+alphaL[21])+0.223606797749979*alphaL[3])*fUpL[31]+0.25*(alphaL[5]*fUpL[27]+fUpL[5]*alphaL[27]+alphaL[8]*fUpL[24]+fUpL[8]*alphaL[24]+fUpL[9]*alphaL[23])+0.223606797749979*(alphaL[6]*fUpL[18]+alphaL[7]*fUpL[17]+fUpL[7]*alphaL[17])+0.2500000000000001*alphaL[13]*fUpL[16]+0.223606797749979*(alphaL[10]*fUpL[15]+fUpL[10]*alphaL[15]); 
  GhatL[47] = (0.223606797749979*(alphaL[13]+alphaL[12]+alphaL[11])+0.25*alphaL[0])*fUpL[47]+(0.1788854381999831*alphaL[37]+0.2*alphaL[10])*fUpL[46]+(0.1788854381999831*alphaL[39]+0.2*alphaL[8])*fUpL[44]+(0.223606797749979*(alphaL[23]+alphaL[20])+0.2500000000000001*alphaL[1])*fUpL[43]+(0.223606797749979*(alphaL[24]+alphaL[19])+0.2500000000000001*alphaL[2])*fUpL[42]+(0.223606797749979*(alphaL[22]+alphaL[21])+0.2500000000000001*alphaL[3])*fUpL[41]+0.2*(alphaL[17]*fUpL[40]+fUpL[18]*alphaL[39]+fUpL[16]*alphaL[37]+alphaL[17]*fUpL[35])+0.223606797749979*(fUpL[30]*alphaL[34]+fUpL[29]*alphaL[33]+fUpL[28]*alphaL[32])+(0.2*(alphaL[27]+alphaL[25])+0.223606797749979*alphaL[4])*fUpL[31]+0.25*(alphaL[5]*fUpL[30]+alphaL[6]*fUpL[29]+alphaL[7]*fUpL[28])+0.223606797749979*(alphaL[8]*fUpL[18]+fUpL[9]*alphaL[17]+alphaL[10]*fUpL[16])+0.2500000000000001*fUpL[14]*alphaL[15]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[3] += 0.7071067811865475*GhatL[3]*rdvpar2; 
  out[4] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdvpar2; 
  out[6] += 0.7071067811865475*GhatL[5]*rdvpar2; 
  out[7] += 0.7071067811865475*GhatL[6]*rdvpar2; 
  out[8] += 0.7071067811865475*GhatL[7]*rdvpar2; 
  out[9] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[10] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[11] += -1.224744871391589*GhatL[3]*rdvpar2; 
  out[12] += 0.7071067811865475*GhatL[8]*rdvpar2; 
  out[13] += 0.7071067811865475*GhatL[9]*rdvpar2; 
  out[14] += 0.7071067811865475*GhatL[10]*rdvpar2; 
  out[15] += -1.224744871391589*GhatL[4]*rdvpar2; 
  out[16] += 0.7071067811865475*GhatL[11]*rdvpar2; 
  out[17] += 0.7071067811865475*GhatL[12]*rdvpar2; 
  out[18] += 0.7071067811865475*GhatL[13]*rdvpar2; 
  out[19] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[20] += 0.7071067811865475*GhatL[14]*rdvpar2; 
  out[21] += 0.7071067811865475*GhatL[15]*rdvpar2; 
  out[22] += -1.224744871391589*GhatL[5]*rdvpar2; 
  out[23] += -1.224744871391589*GhatL[6]*rdvpar2; 
  out[24] += -1.224744871391589*GhatL[7]*rdvpar2; 
  out[25] += 0.7071067811865475*GhatL[16]*rdvpar2; 
  out[26] += 0.7071067811865475*GhatL[17]*rdvpar2; 
  out[27] += 0.7071067811865475*GhatL[18]*rdvpar2; 
  out[28] += -1.224744871391589*GhatL[8]*rdvpar2; 
  out[29] += -1.224744871391589*GhatL[9]*rdvpar2; 
  out[30] += -1.224744871391589*GhatL[10]*rdvpar2; 
  out[31] += 0.7071067811865475*GhatL[19]*rdvpar2; 
  out[32] += 0.7071067811865475*GhatL[20]*rdvpar2; 
  out[33] += 0.7071067811865475*GhatL[21]*rdvpar2; 
  out[34] += 0.7071067811865475*GhatL[22]*rdvpar2; 
  out[35] += 0.7071067811865475*GhatL[23]*rdvpar2; 
  out[36] += 0.7071067811865475*GhatL[24]*rdvpar2; 
  out[37] += -1.224744871391589*GhatL[11]*rdvpar2; 
  out[38] += -1.224744871391589*GhatL[12]*rdvpar2; 
  out[39] += -1.224744871391589*GhatL[13]*rdvpar2; 
  out[40] += 1.58113883008419*GhatL[1]*rdvpar2; 
  out[41] += 1.58113883008419*GhatL[2]*rdvpar2; 
  out[42] += 1.58113883008419*GhatL[3]*rdvpar2; 
  out[43] += 0.7071067811865475*GhatL[25]*rdvpar2; 
  out[44] += 0.7071067811865475*GhatL[26]*rdvpar2; 
  out[45] += 0.7071067811865475*GhatL[27]*rdvpar2; 
  out[46] += 1.58113883008419*GhatL[4]*rdvpar2; 
  out[47] += 0.7071067811865475*GhatL[28]*rdvpar2; 
  out[48] += 0.7071067811865475*GhatL[29]*rdvpar2; 
  out[49] += 0.7071067811865475*GhatL[30]*rdvpar2; 
  out[50] += -1.224744871391589*GhatL[14]*rdvpar2; 
  out[51] += -1.224744871391589*GhatL[15]*rdvpar2; 
  out[52] += 0.7071067811865475*GhatL[31]*rdvpar2; 
  out[53] += -1.224744871391589*GhatL[16]*rdvpar2; 
  out[54] += -1.224744871391589*GhatL[17]*rdvpar2; 
  out[55] += -1.224744871391589*GhatL[18]*rdvpar2; 
  out[56] += 0.7071067811865475*GhatL[32]*rdvpar2; 
  out[57] += 0.7071067811865475*GhatL[33]*rdvpar2; 
  out[58] += 0.7071067811865475*GhatL[34]*rdvpar2; 
  out[59] += -1.224744871391589*GhatL[19]*rdvpar2; 
  out[60] += -1.224744871391589*GhatL[20]*rdvpar2; 
  out[61] += -1.224744871391589*GhatL[21]*rdvpar2; 
  out[62] += -1.224744871391589*GhatL[22]*rdvpar2; 
  out[63] += -1.224744871391589*GhatL[23]*rdvpar2; 
  out[64] += -1.224744871391589*GhatL[24]*rdvpar2; 
  out[65] += 1.58113883008419*GhatL[5]*rdvpar2; 
  out[66] += 1.58113883008419*GhatL[6]*rdvpar2; 
  out[67] += 1.58113883008419*GhatL[7]*rdvpar2; 
  out[68] += 0.7071067811865475*GhatL[35]*rdvpar2; 
  out[69] += 0.7071067811865475*GhatL[36]*rdvpar2; 
  out[70] += 0.7071067811865475*GhatL[37]*rdvpar2; 
  out[71] += 0.7071067811865475*GhatL[38]*rdvpar2; 
  out[72] += 0.7071067811865475*GhatL[39]*rdvpar2; 
  out[73] += 0.7071067811865475*GhatL[40]*rdvpar2; 
  out[74] += -1.224744871391589*GhatL[25]*rdvpar2; 
  out[75] += -1.224744871391589*GhatL[26]*rdvpar2; 
  out[76] += -1.224744871391589*GhatL[27]*rdvpar2; 
  out[77] += 1.58113883008419*GhatL[8]*rdvpar2; 
  out[78] += 1.58113883008419*GhatL[9]*rdvpar2; 
  out[79] += 1.58113883008419*GhatL[10]*rdvpar2; 
  out[80] += 0.7071067811865475*GhatL[41]*rdvpar2; 
  out[81] += 0.7071067811865475*GhatL[42]*rdvpar2; 
  out[82] += 0.7071067811865475*GhatL[43]*rdvpar2; 
  out[83] += -1.224744871391589*GhatL[28]*rdvpar2; 
  out[84] += -1.224744871391589*GhatL[29]*rdvpar2; 
  out[85] += -1.224744871391589*GhatL[30]*rdvpar2; 
  out[86] += -1.224744871391589*GhatL[31]*rdvpar2; 
  out[87] += -1.224744871391589*GhatL[32]*rdvpar2; 
  out[88] += -1.224744871391589*GhatL[33]*rdvpar2; 
  out[89] += -1.224744871391589*GhatL[34]*rdvpar2; 
  out[90] += 1.58113883008419*GhatL[15]*rdvpar2; 
  out[91] += 0.7071067811865475*GhatL[44]*rdvpar2; 
  out[92] += 0.7071067811865475*GhatL[45]*rdvpar2; 
  out[93] += 0.7071067811865475*GhatL[46]*rdvpar2; 
  out[94] += -1.224744871391589*GhatL[35]*rdvpar2; 
  out[95] += -1.224744871391589*GhatL[36]*rdvpar2; 
  out[96] += -1.224744871391589*GhatL[37]*rdvpar2; 
  out[97] += -1.224744871391589*GhatL[38]*rdvpar2; 
  out[98] += -1.224744871391589*GhatL[39]*rdvpar2; 
  out[99] += -1.224744871391589*GhatL[40]*rdvpar2; 
  out[100] += 1.58113883008419*GhatL[16]*rdvpar2; 
  out[101] += 1.58113883008419*GhatL[17]*rdvpar2; 
  out[102] += 1.58113883008419*GhatL[18]*rdvpar2; 
  out[103] += 0.7071067811865475*GhatL[47]*rdvpar2; 
  out[104] += -1.224744871391589*GhatL[41]*rdvpar2; 
  out[105] += -1.224744871391589*GhatL[42]*rdvpar2; 
  out[106] += -1.224744871391589*GhatL[43]*rdvpar2; 
  out[107] += -1.224744871391589*GhatL[44]*rdvpar2; 
  out[108] += -1.224744871391589*GhatL[45]*rdvpar2; 
  out[109] += -1.224744871391589*GhatL[46]*rdvpar2; 
  out[110] += 1.58113883008419*GhatL[31]*rdvpar2; 
  out[111] += -1.224744871391589*GhatL[47]*rdvpar2; 

  } 

  return 2.5*rdvpar2*cflFreq; 

} 
