#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfx_3x2v_ser_p1(const double *w, const double *dxv, const double *alpha_surf_l, const double *alpha_surf_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // alpha_surf_l: Surface expansion of phase space flux on the left.
  // alpha_surf_r: Surface expansion of phase space flux on the right.
  // fl,fc,fr: distribution function in left, center and right cells.
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

  const double *alphaL = &alpha_surf_l[0];
  const double *alphaR = &alpha_surf_r[0];
  double cflFreq = 0.0;
  double fUpOrdL[24] = {0.};
  double alphaL_n = 0.;

  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]+0.25*alphaL[9]+0.3354101966249678*alphaL[7]+0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx1_eval_quad_node_0_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]+0.25*alphaL[9]+0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx1_eval_quad_node_1_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]+0.25*alphaL[9]-0.3354101966249678*alphaL[7]+0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx1_eval_quad_node_2_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]-0.25*alphaL[9]+0.3354101966249678*alphaL[7]+0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx1_eval_quad_node_3_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.25*alphaL[9]+0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx1_eval_quad_node_4_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]-0.25*alphaL[9]-0.3354101966249678*alphaL[7]+0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx1_eval_quad_node_5_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.25*alphaL[9]-0.3354101966249678*alphaL[7]-0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx1_eval_quad_node_6_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[18])-0.2795084971874732*alphaL[16]-0.25*alphaL[9]-0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx1_eval_quad_node_7_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.25*alphaL[9]+0.3354101966249678*alphaL[7]-0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx1_eval_quad_node_8_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.25*alphaL[9]-0.3354101966249678*alphaL[7]-0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx1_eval_quad_node_9_r(fl); 
  } else { 
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx1_eval_quad_node_9_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[18])-0.2795084971874732*alphaL[16]+0.25*alphaL[9]-0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx1_eval_quad_node_10_r(fl); 
  } else { 
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx1_eval_quad_node_10_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.25*alphaL[9]+0.3354101966249678*alphaL[7]-0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx1_eval_quad_node_11_r(fl); 
  } else { 
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx1_eval_quad_node_11_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]+0.25*alphaL[9]+0.3354101966249678*alphaL[7]-0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx1_eval_quad_node_12_r(fl); 
  } else { 
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx1_eval_quad_node_12_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]+0.25*alphaL[9]-0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx1_eval_quad_node_13_r(fl); 
  } else { 
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx1_eval_quad_node_13_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]+0.25*alphaL[9]-0.3354101966249678*alphaL[7]-0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx1_eval_quad_node_14_r(fl); 
  } else { 
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx1_eval_quad_node_14_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]-0.25*alphaL[9]+0.3354101966249678*alphaL[7]-0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx1_eval_quad_node_15_r(fl); 
  } else { 
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx1_eval_quad_node_15_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.25*alphaL[9]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[16] = gkhyb_3x2v_p1_surfx1_eval_quad_node_16_r(fl); 
  } else { 
    fUpOrdL[16] = gkhyb_3x2v_p1_surfx1_eval_quad_node_16_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]-0.25*alphaL[9]-0.3354101966249678*alphaL[7]-0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[17] = gkhyb_3x2v_p1_surfx1_eval_quad_node_17_r(fl); 
  } else { 
    fUpOrdL[17] = gkhyb_3x2v_p1_surfx1_eval_quad_node_17_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.25*alphaL[9]-0.3354101966249678*alphaL[7]+0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[18] = gkhyb_3x2v_p1_surfx1_eval_quad_node_18_r(fl); 
  } else { 
    fUpOrdL[18] = gkhyb_3x2v_p1_surfx1_eval_quad_node_18_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[18])-0.2795084971874732*alphaL[16]-0.25*alphaL[9]+0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[19] = gkhyb_3x2v_p1_surfx1_eval_quad_node_19_r(fl); 
  } else { 
    fUpOrdL[19] = gkhyb_3x2v_p1_surfx1_eval_quad_node_19_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.25*alphaL[9]+0.3354101966249678*alphaL[7]+0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[20] = gkhyb_3x2v_p1_surfx1_eval_quad_node_20_r(fl); 
  } else { 
    fUpOrdL[20] = gkhyb_3x2v_p1_surfx1_eval_quad_node_20_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.25*alphaL[9]-0.3354101966249678*alphaL[7]+0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[21] = gkhyb_3x2v_p1_surfx1_eval_quad_node_21_r(fl); 
  } else { 
    fUpOrdL[21] = gkhyb_3x2v_p1_surfx1_eval_quad_node_21_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[18])-0.2795084971874732*alphaL[16]+0.25*alphaL[9]+0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[22] = gkhyb_3x2v_p1_surfx1_eval_quad_node_22_r(fl); 
  } else { 
    fUpOrdL[22] = gkhyb_3x2v_p1_surfx1_eval_quad_node_22_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.25*alphaL[9]+0.3354101966249678*alphaL[7]+0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[23] = gkhyb_3x2v_p1_surfx1_eval_quad_node_23_r(fl); 
  } else { 
    fUpOrdL[23] = gkhyb_3x2v_p1_surfx1_eval_quad_node_23_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[24] = {0.};
  gkhyb_3x2v_p1_xdir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[24] = {0.}; 
  GhatL[0] = 0.25*(alphaL[18]*fUpL[18]+alphaL[16]*fUpL[16]+alphaL[9]*fUpL[9]+alphaL[7]*fUpL[7]+alphaL[5]*fUpL[5]+alphaL[4]*fUpL[4]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.2500000000000001*(alphaL[18]*fUpL[20]+alphaL[16]*fUpL[17])+0.25*(alphaL[9]*fUpL[12]+alphaL[7]*fUpL[11]+alphaL[4]*fUpL[8]+alphaL[3]*fUpL[6]+alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.2500000000000001*(alphaL[16]*fUpL[18]+fUpL[16]*alphaL[18])+0.25*(alphaL[4]*fUpL[9]+fUpL[4]*alphaL[9]+alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]+alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.223606797749979*(alphaL[7]*fUpL[18]+fUpL[7]*alphaL[18])+0.223606797749979*(alphaL[3]*fUpL[16]+fUpL[3]*alphaL[16])+0.25*(alphaL[9]*fUpL[14]+alphaL[5]*fUpL[11]+alphaL[4]*fUpL[10]+alphaL[2]*fUpL[7]+fUpL[2]*alphaL[7]+alphaL[1]*fUpL[6]+alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]); 
  GhatL[4] = 0.2500000000000001*(alphaL[18]*fUpL[22]+alphaL[16]*fUpL[19])+0.25*(alphaL[7]*fUpL[14]+alphaL[5]*fUpL[12]+alphaL[3]*fUpL[10]+alphaL[2]*fUpL[9]+fUpL[2]*alphaL[9]+alphaL[1]*fUpL[8]+alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]); 
  GhatL[5] = 0.25*(alphaL[16]*fUpL[20]+fUpL[17]*alphaL[18]+alphaL[4]*fUpL[12]+alphaL[3]*fUpL[11]+fUpL[8]*alphaL[9]+fUpL[6]*alphaL[7]+alphaL[0]*fUpL[5]+fUpL[0]*alphaL[5]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[6] = 0.223606797749979*alphaL[7]*fUpL[20]+0.223606797749979*(fUpL[11]*alphaL[18]+alphaL[3]*fUpL[17])+0.223606797749979*fUpL[6]*alphaL[16]+0.25*(alphaL[9]*fUpL[15]+alphaL[4]*fUpL[13]+alphaL[2]*fUpL[11]+alphaL[5]*fUpL[7]+fUpL[5]*alphaL[7]+alphaL[0]*fUpL[6]+alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[7] = 0.223606797749979*(alphaL[3]*fUpL[18]+fUpL[3]*alphaL[18])+0.223606797749979*(alphaL[7]*fUpL[16]+fUpL[7]*alphaL[16])+0.25*(alphaL[4]*fUpL[14]+alphaL[1]*fUpL[11]+alphaL[9]*fUpL[10]+alphaL[0]*fUpL[7]+fUpL[0]*alphaL[7]+alphaL[5]*fUpL[6]+alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 
  GhatL[8] = 0.25*(alphaL[18]*fUpL[23]+alphaL[16]*fUpL[21]+alphaL[7]*fUpL[15]+alphaL[3]*fUpL[13]+alphaL[2]*fUpL[12]+alphaL[5]*fUpL[9]+fUpL[5]*alphaL[9]+alphaL[0]*fUpL[8]+alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]); 
  GhatL[9] = 0.25*(alphaL[16]*fUpL[22]+alphaL[18]*fUpL[19]+alphaL[3]*fUpL[14]+alphaL[1]*fUpL[12]+alphaL[7]*fUpL[10]+alphaL[0]*fUpL[9]+fUpL[0]*alphaL[9]+alphaL[5]*fUpL[8]+alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]); 
  GhatL[10] = 0.223606797749979*alphaL[7]*fUpL[22]+0.223606797749979*(alphaL[3]*fUpL[19]+fUpL[14]*alphaL[18])+0.223606797749979*fUpL[10]*alphaL[16]+0.25*(alphaL[5]*fUpL[15]+alphaL[2]*fUpL[14]+alphaL[1]*fUpL[13]+alphaL[0]*fUpL[10]+alphaL[7]*fUpL[9]+fUpL[7]*alphaL[9]+alphaL[3]*fUpL[4]+fUpL[3]*alphaL[4]); 
  GhatL[11] = 0.223606797749979*alphaL[3]*fUpL[20]+0.223606797749979*(fUpL[6]*alphaL[18]+alphaL[7]*fUpL[17])+0.223606797749979*fUpL[11]*alphaL[16]+0.25*(alphaL[4]*fUpL[15]+alphaL[9]*fUpL[13]+alphaL[0]*fUpL[11]+alphaL[1]*fUpL[7]+fUpL[1]*alphaL[7]+alphaL[2]*fUpL[6]+alphaL[3]*fUpL[5]+fUpL[3]*alphaL[5]); 
  GhatL[12] = 0.2500000000000001*(alphaL[16]*fUpL[23]+alphaL[18]*fUpL[21])+0.25*(alphaL[3]*fUpL[15]+alphaL[7]*fUpL[13]+alphaL[0]*fUpL[12]+alphaL[1]*fUpL[9]+fUpL[1]*alphaL[9]+alphaL[2]*fUpL[8]+alphaL[4]*fUpL[5]+fUpL[4]*alphaL[5]); 
  GhatL[13] = 0.223606797749979*alphaL[7]*fUpL[23]+0.223606797749979*alphaL[3]*fUpL[21]+0.223606797749979*fUpL[15]*alphaL[18]+0.223606797749979*fUpL[13]*alphaL[16]+0.25*(alphaL[2]*fUpL[15]+alphaL[5]*fUpL[14]+alphaL[0]*fUpL[13]+alphaL[7]*fUpL[12]+alphaL[9]*fUpL[11]+alphaL[1]*fUpL[10]+alphaL[3]*fUpL[8]+alphaL[4]*fUpL[6]); 
  GhatL[14] = 0.223606797749979*alphaL[3]*fUpL[22]+0.223606797749979*(alphaL[7]*fUpL[19]+fUpL[10]*alphaL[18])+0.223606797749979*fUpL[14]*alphaL[16]+0.25*(alphaL[1]*fUpL[15]+alphaL[0]*fUpL[14]+alphaL[5]*fUpL[13]+alphaL[2]*fUpL[10]+alphaL[3]*fUpL[9]+fUpL[3]*alphaL[9]+alphaL[4]*fUpL[7]+fUpL[4]*alphaL[7]); 
  GhatL[15] = 0.223606797749979*alphaL[3]*fUpL[23]+0.223606797749979*alphaL[7]*fUpL[21]+0.223606797749979*fUpL[13]*alphaL[18]+0.223606797749979*fUpL[15]*alphaL[16]+0.25*(alphaL[0]*fUpL[15]+alphaL[1]*fUpL[14]+alphaL[2]*fUpL[13]+alphaL[3]*fUpL[12]+alphaL[4]*fUpL[11]+alphaL[5]*fUpL[10]+fUpL[6]*alphaL[9]+alphaL[7]*fUpL[8]); 
  GhatL[16] = 0.25*(alphaL[9]*fUpL[22]+alphaL[5]*fUpL[20])+0.2500000000000001*alphaL[4]*fUpL[19]+0.159719141249985*alphaL[18]*fUpL[18]+0.2500000000000001*(alphaL[2]*fUpL[18]+fUpL[2]*alphaL[18]+alphaL[1]*fUpL[17])+0.159719141249985*alphaL[16]*fUpL[16]+0.25*(alphaL[0]*fUpL[16]+fUpL[0]*alphaL[16])+0.223606797749979*(alphaL[7]*fUpL[7]+alphaL[3]*fUpL[3]); 
  GhatL[17] = 0.25*alphaL[9]*fUpL[23]+0.2500000000000001*alphaL[4]*fUpL[21]+(0.159719141249985*alphaL[18]+0.2500000000000001*alphaL[2])*fUpL[20]+0.25*(alphaL[5]*fUpL[18]+fUpL[5]*alphaL[18])+(0.159719141249985*alphaL[16]+0.25*alphaL[0])*fUpL[17]+0.2500000000000001*(alphaL[1]*fUpL[16]+fUpL[1]*alphaL[16])+0.223606797749979*(alphaL[7]*fUpL[11]+alphaL[3]*fUpL[6]); 
  GhatL[18] = 0.2500000000000001*(alphaL[4]*fUpL[22]+alphaL[1]*fUpL[20])+0.25*alphaL[9]*fUpL[19]+(0.159719141249985*alphaL[16]+0.25*alphaL[0])*fUpL[18]+0.159719141249985*fUpL[16]*alphaL[18]+0.25*(fUpL[0]*alphaL[18]+alphaL[5]*fUpL[17])+0.2500000000000001*(alphaL[2]*fUpL[16]+fUpL[2]*alphaL[16])+0.223606797749979*(alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]); 
  GhatL[19] = 0.25*alphaL[5]*fUpL[23]+0.159719141249985*alphaL[18]*fUpL[22]+0.2500000000000001*(alphaL[2]*fUpL[22]+alphaL[1]*fUpL[21])+0.159719141249985*alphaL[16]*fUpL[19]+0.25*(alphaL[0]*fUpL[19]+alphaL[9]*fUpL[18]+fUpL[9]*alphaL[18])+0.2500000000000001*(alphaL[4]*fUpL[16]+fUpL[4]*alphaL[16])+0.223606797749979*(alphaL[7]*fUpL[14]+alphaL[3]*fUpL[10]); 
  GhatL[20] = 0.2500000000000001*alphaL[4]*fUpL[23]+0.25*alphaL[9]*fUpL[21]+(0.159719141249985*alphaL[16]+0.25*alphaL[0])*fUpL[20]+0.2500000000000001*alphaL[1]*fUpL[18]+0.159719141249985*fUpL[17]*alphaL[18]+0.2500000000000001*(fUpL[1]*alphaL[18]+alphaL[2]*fUpL[17])+0.25*(alphaL[5]*fUpL[16]+fUpL[5]*alphaL[16])+0.223606797749979*(alphaL[3]*fUpL[11]+fUpL[6]*alphaL[7]); 
  GhatL[21] = (0.159719141249985*alphaL[18]+0.2500000000000001*alphaL[2])*fUpL[23]+0.25*alphaL[5]*fUpL[22]+0.159719141249985*alphaL[16]*fUpL[21]+0.25*(alphaL[0]*fUpL[21]+alphaL[9]*fUpL[20])+0.2500000000000001*(alphaL[1]*fUpL[19]+fUpL[12]*alphaL[18]+alphaL[4]*fUpL[17])+0.25*fUpL[8]*alphaL[16]+0.223606797749979*(alphaL[7]*fUpL[15]+alphaL[3]*fUpL[13]); 
  GhatL[22] = 0.2500000000000001*alphaL[1]*fUpL[23]+0.159719141249985*alphaL[16]*fUpL[22]+0.25*(alphaL[0]*fUpL[22]+alphaL[5]*fUpL[21])+0.159719141249985*alphaL[18]*fUpL[19]+0.2500000000000001*(alphaL[2]*fUpL[19]+alphaL[4]*fUpL[18]+fUpL[4]*alphaL[18])+0.25*(alphaL[9]*fUpL[16]+fUpL[9]*alphaL[16])+0.223606797749979*(alphaL[3]*fUpL[14]+alphaL[7]*fUpL[10]); 
  GhatL[23] = (0.159719141249985*alphaL[16]+0.25*alphaL[0])*fUpL[23]+0.2500000000000001*alphaL[1]*fUpL[22]+0.159719141249985*alphaL[18]*fUpL[21]+0.2500000000000001*(alphaL[2]*fUpL[21]+alphaL[4]*fUpL[20])+0.25*(alphaL[5]*fUpL[19]+fUpL[8]*alphaL[18]+alphaL[9]*fUpL[17])+0.2500000000000001*fUpL[12]*alphaL[16]+0.223606797749979*(alphaL[3]*fUpL[15]+alphaL[7]*fUpL[13]); 

  double fUpOrdR[24] = {0.};
  double alphaR_n = 0.;

  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]+0.25*alphaR[9]+0.3354101966249678*alphaR[7]+0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx1_eval_quad_node_0_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[16]+0.25*alphaR[9]+0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx1_eval_quad_node_1_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]+0.25*alphaR[9]-0.3354101966249678*alphaR[7]+0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx1_eval_quad_node_2_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]-0.25*alphaR[9]+0.3354101966249678*alphaR[7]+0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx1_eval_quad_node_3_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[16]-0.25*alphaR[9]+0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx1_eval_quad_node_4_r(fc); 
  } else { 
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx1_eval_quad_node_4_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]-0.25*alphaR[9]-0.3354101966249678*alphaR[7]+0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx1_eval_quad_node_5_r(fc); 
  } else { 
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx1_eval_quad_node_5_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]-0.25*alphaR[9]-0.3354101966249678*alphaR[7]-0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx1_eval_quad_node_6_r(fc); 
  } else { 
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx1_eval_quad_node_6_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[18])-0.2795084971874732*alphaR[16]-0.25*alphaR[9]-0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx1_eval_quad_node_7_r(fc); 
  } else { 
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx1_eval_quad_node_7_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]-0.25*alphaR[9]+0.3354101966249678*alphaR[7]-0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx1_eval_quad_node_8_r(fc); 
  } else { 
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx1_eval_quad_node_8_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]+0.25*alphaR[9]-0.3354101966249678*alphaR[7]-0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx1_eval_quad_node_9_r(fc); 
  } else { 
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx1_eval_quad_node_9_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[18])-0.2795084971874732*alphaR[16]+0.25*alphaR[9]-0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx1_eval_quad_node_10_r(fc); 
  } else { 
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx1_eval_quad_node_10_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]+0.25*alphaR[9]+0.3354101966249678*alphaR[7]-0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx1_eval_quad_node_11_r(fc); 
  } else { 
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx1_eval_quad_node_11_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]+0.25*alphaR[9]+0.3354101966249678*alphaR[7]-0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx1_eval_quad_node_12_r(fc); 
  } else { 
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx1_eval_quad_node_12_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[16]+0.25*alphaR[9]-0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx1_eval_quad_node_13_r(fc); 
  } else { 
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx1_eval_quad_node_13_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]+0.25*alphaR[9]-0.3354101966249678*alphaR[7]-0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx1_eval_quad_node_14_r(fc); 
  } else { 
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx1_eval_quad_node_14_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]-0.25*alphaR[9]+0.3354101966249678*alphaR[7]-0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx1_eval_quad_node_15_r(fc); 
  } else { 
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx1_eval_quad_node_15_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[16]-0.25*alphaR[9]-0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[16] = gkhyb_3x2v_p1_surfx1_eval_quad_node_16_r(fc); 
  } else { 
    fUpOrdR[16] = gkhyb_3x2v_p1_surfx1_eval_quad_node_16_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]-0.25*alphaR[9]-0.3354101966249678*alphaR[7]-0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[17] = gkhyb_3x2v_p1_surfx1_eval_quad_node_17_r(fc); 
  } else { 
    fUpOrdR[17] = gkhyb_3x2v_p1_surfx1_eval_quad_node_17_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]-0.25*alphaR[9]-0.3354101966249678*alphaR[7]+0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[18] = gkhyb_3x2v_p1_surfx1_eval_quad_node_18_r(fc); 
  } else { 
    fUpOrdR[18] = gkhyb_3x2v_p1_surfx1_eval_quad_node_18_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[18])-0.2795084971874732*alphaR[16]-0.25*alphaR[9]+0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[19] = gkhyb_3x2v_p1_surfx1_eval_quad_node_19_r(fc); 
  } else { 
    fUpOrdR[19] = gkhyb_3x2v_p1_surfx1_eval_quad_node_19_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]-0.25*alphaR[9]+0.3354101966249678*alphaR[7]+0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[20] = gkhyb_3x2v_p1_surfx1_eval_quad_node_20_r(fc); 
  } else { 
    fUpOrdR[20] = gkhyb_3x2v_p1_surfx1_eval_quad_node_20_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]+0.25*alphaR[9]-0.3354101966249678*alphaR[7]+0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[21] = gkhyb_3x2v_p1_surfx1_eval_quad_node_21_r(fc); 
  } else { 
    fUpOrdR[21] = gkhyb_3x2v_p1_surfx1_eval_quad_node_21_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[18])-0.2795084971874732*alphaR[16]+0.25*alphaR[9]+0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[22] = gkhyb_3x2v_p1_surfx1_eval_quad_node_22_r(fc); 
  } else { 
    fUpOrdR[22] = gkhyb_3x2v_p1_surfx1_eval_quad_node_22_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]+0.25*alphaR[9]+0.3354101966249678*alphaR[7]+0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[23] = gkhyb_3x2v_p1_surfx1_eval_quad_node_23_r(fc); 
  } else { 
    fUpOrdR[23] = gkhyb_3x2v_p1_surfx1_eval_quad_node_23_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[24] = {0.};
  gkhyb_3x2v_p1_xdir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[24] = {0.}; 
  GhatR[0] = 0.25*(alphaR[18]*fUpR[18]+alphaR[16]*fUpR[16]+alphaR[9]*fUpR[9]+alphaR[7]*fUpR[7]+alphaR[5]*fUpR[5]+alphaR[4]*fUpR[4]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.2500000000000001*(alphaR[18]*fUpR[20]+alphaR[16]*fUpR[17])+0.25*(alphaR[9]*fUpR[12]+alphaR[7]*fUpR[11]+alphaR[4]*fUpR[8]+alphaR[3]*fUpR[6]+alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.2500000000000001*(alphaR[16]*fUpR[18]+fUpR[16]*alphaR[18])+0.25*(alphaR[4]*fUpR[9]+fUpR[4]*alphaR[9]+alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]+alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.223606797749979*(alphaR[7]*fUpR[18]+fUpR[7]*alphaR[18])+0.223606797749979*(alphaR[3]*fUpR[16]+fUpR[3]*alphaR[16])+0.25*(alphaR[9]*fUpR[14]+alphaR[5]*fUpR[11]+alphaR[4]*fUpR[10]+alphaR[2]*fUpR[7]+fUpR[2]*alphaR[7]+alphaR[1]*fUpR[6]+alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]); 
  GhatR[4] = 0.2500000000000001*(alphaR[18]*fUpR[22]+alphaR[16]*fUpR[19])+0.25*(alphaR[7]*fUpR[14]+alphaR[5]*fUpR[12]+alphaR[3]*fUpR[10]+alphaR[2]*fUpR[9]+fUpR[2]*alphaR[9]+alphaR[1]*fUpR[8]+alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]); 
  GhatR[5] = 0.25*(alphaR[16]*fUpR[20]+fUpR[17]*alphaR[18]+alphaR[4]*fUpR[12]+alphaR[3]*fUpR[11]+fUpR[8]*alphaR[9]+fUpR[6]*alphaR[7]+alphaR[0]*fUpR[5]+fUpR[0]*alphaR[5]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[6] = 0.223606797749979*alphaR[7]*fUpR[20]+0.223606797749979*(fUpR[11]*alphaR[18]+alphaR[3]*fUpR[17])+0.223606797749979*fUpR[6]*alphaR[16]+0.25*(alphaR[9]*fUpR[15]+alphaR[4]*fUpR[13]+alphaR[2]*fUpR[11]+alphaR[5]*fUpR[7]+fUpR[5]*alphaR[7]+alphaR[0]*fUpR[6]+alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[7] = 0.223606797749979*(alphaR[3]*fUpR[18]+fUpR[3]*alphaR[18])+0.223606797749979*(alphaR[7]*fUpR[16]+fUpR[7]*alphaR[16])+0.25*(alphaR[4]*fUpR[14]+alphaR[1]*fUpR[11]+alphaR[9]*fUpR[10]+alphaR[0]*fUpR[7]+fUpR[0]*alphaR[7]+alphaR[5]*fUpR[6]+alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 
  GhatR[8] = 0.25*(alphaR[18]*fUpR[23]+alphaR[16]*fUpR[21]+alphaR[7]*fUpR[15]+alphaR[3]*fUpR[13]+alphaR[2]*fUpR[12]+alphaR[5]*fUpR[9]+fUpR[5]*alphaR[9]+alphaR[0]*fUpR[8]+alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]); 
  GhatR[9] = 0.25*(alphaR[16]*fUpR[22]+alphaR[18]*fUpR[19]+alphaR[3]*fUpR[14]+alphaR[1]*fUpR[12]+alphaR[7]*fUpR[10]+alphaR[0]*fUpR[9]+fUpR[0]*alphaR[9]+alphaR[5]*fUpR[8]+alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]); 
  GhatR[10] = 0.223606797749979*alphaR[7]*fUpR[22]+0.223606797749979*(alphaR[3]*fUpR[19]+fUpR[14]*alphaR[18])+0.223606797749979*fUpR[10]*alphaR[16]+0.25*(alphaR[5]*fUpR[15]+alphaR[2]*fUpR[14]+alphaR[1]*fUpR[13]+alphaR[0]*fUpR[10]+alphaR[7]*fUpR[9]+fUpR[7]*alphaR[9]+alphaR[3]*fUpR[4]+fUpR[3]*alphaR[4]); 
  GhatR[11] = 0.223606797749979*alphaR[3]*fUpR[20]+0.223606797749979*(fUpR[6]*alphaR[18]+alphaR[7]*fUpR[17])+0.223606797749979*fUpR[11]*alphaR[16]+0.25*(alphaR[4]*fUpR[15]+alphaR[9]*fUpR[13]+alphaR[0]*fUpR[11]+alphaR[1]*fUpR[7]+fUpR[1]*alphaR[7]+alphaR[2]*fUpR[6]+alphaR[3]*fUpR[5]+fUpR[3]*alphaR[5]); 
  GhatR[12] = 0.2500000000000001*(alphaR[16]*fUpR[23]+alphaR[18]*fUpR[21])+0.25*(alphaR[3]*fUpR[15]+alphaR[7]*fUpR[13]+alphaR[0]*fUpR[12]+alphaR[1]*fUpR[9]+fUpR[1]*alphaR[9]+alphaR[2]*fUpR[8]+alphaR[4]*fUpR[5]+fUpR[4]*alphaR[5]); 
  GhatR[13] = 0.223606797749979*alphaR[7]*fUpR[23]+0.223606797749979*alphaR[3]*fUpR[21]+0.223606797749979*fUpR[15]*alphaR[18]+0.223606797749979*fUpR[13]*alphaR[16]+0.25*(alphaR[2]*fUpR[15]+alphaR[5]*fUpR[14]+alphaR[0]*fUpR[13]+alphaR[7]*fUpR[12]+alphaR[9]*fUpR[11]+alphaR[1]*fUpR[10]+alphaR[3]*fUpR[8]+alphaR[4]*fUpR[6]); 
  GhatR[14] = 0.223606797749979*alphaR[3]*fUpR[22]+0.223606797749979*(alphaR[7]*fUpR[19]+fUpR[10]*alphaR[18])+0.223606797749979*fUpR[14]*alphaR[16]+0.25*(alphaR[1]*fUpR[15]+alphaR[0]*fUpR[14]+alphaR[5]*fUpR[13]+alphaR[2]*fUpR[10]+alphaR[3]*fUpR[9]+fUpR[3]*alphaR[9]+alphaR[4]*fUpR[7]+fUpR[4]*alphaR[7]); 
  GhatR[15] = 0.223606797749979*alphaR[3]*fUpR[23]+0.223606797749979*alphaR[7]*fUpR[21]+0.223606797749979*fUpR[13]*alphaR[18]+0.223606797749979*fUpR[15]*alphaR[16]+0.25*(alphaR[0]*fUpR[15]+alphaR[1]*fUpR[14]+alphaR[2]*fUpR[13]+alphaR[3]*fUpR[12]+alphaR[4]*fUpR[11]+alphaR[5]*fUpR[10]+fUpR[6]*alphaR[9]+alphaR[7]*fUpR[8]); 
  GhatR[16] = 0.25*(alphaR[9]*fUpR[22]+alphaR[5]*fUpR[20])+0.2500000000000001*alphaR[4]*fUpR[19]+0.159719141249985*alphaR[18]*fUpR[18]+0.2500000000000001*(alphaR[2]*fUpR[18]+fUpR[2]*alphaR[18]+alphaR[1]*fUpR[17])+0.159719141249985*alphaR[16]*fUpR[16]+0.25*(alphaR[0]*fUpR[16]+fUpR[0]*alphaR[16])+0.223606797749979*(alphaR[7]*fUpR[7]+alphaR[3]*fUpR[3]); 
  GhatR[17] = 0.25*alphaR[9]*fUpR[23]+0.2500000000000001*alphaR[4]*fUpR[21]+(0.159719141249985*alphaR[18]+0.2500000000000001*alphaR[2])*fUpR[20]+0.25*(alphaR[5]*fUpR[18]+fUpR[5]*alphaR[18])+(0.159719141249985*alphaR[16]+0.25*alphaR[0])*fUpR[17]+0.2500000000000001*(alphaR[1]*fUpR[16]+fUpR[1]*alphaR[16])+0.223606797749979*(alphaR[7]*fUpR[11]+alphaR[3]*fUpR[6]); 
  GhatR[18] = 0.2500000000000001*(alphaR[4]*fUpR[22]+alphaR[1]*fUpR[20])+0.25*alphaR[9]*fUpR[19]+(0.159719141249985*alphaR[16]+0.25*alphaR[0])*fUpR[18]+0.159719141249985*fUpR[16]*alphaR[18]+0.25*(fUpR[0]*alphaR[18]+alphaR[5]*fUpR[17])+0.2500000000000001*(alphaR[2]*fUpR[16]+fUpR[2]*alphaR[16])+0.223606797749979*(alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]); 
  GhatR[19] = 0.25*alphaR[5]*fUpR[23]+0.159719141249985*alphaR[18]*fUpR[22]+0.2500000000000001*(alphaR[2]*fUpR[22]+alphaR[1]*fUpR[21])+0.159719141249985*alphaR[16]*fUpR[19]+0.25*(alphaR[0]*fUpR[19]+alphaR[9]*fUpR[18]+fUpR[9]*alphaR[18])+0.2500000000000001*(alphaR[4]*fUpR[16]+fUpR[4]*alphaR[16])+0.223606797749979*(alphaR[7]*fUpR[14]+alphaR[3]*fUpR[10]); 
  GhatR[20] = 0.2500000000000001*alphaR[4]*fUpR[23]+0.25*alphaR[9]*fUpR[21]+(0.159719141249985*alphaR[16]+0.25*alphaR[0])*fUpR[20]+0.2500000000000001*alphaR[1]*fUpR[18]+0.159719141249985*fUpR[17]*alphaR[18]+0.2500000000000001*(fUpR[1]*alphaR[18]+alphaR[2]*fUpR[17])+0.25*(alphaR[5]*fUpR[16]+fUpR[5]*alphaR[16])+0.223606797749979*(alphaR[3]*fUpR[11]+fUpR[6]*alphaR[7]); 
  GhatR[21] = (0.159719141249985*alphaR[18]+0.2500000000000001*alphaR[2])*fUpR[23]+0.25*alphaR[5]*fUpR[22]+0.159719141249985*alphaR[16]*fUpR[21]+0.25*(alphaR[0]*fUpR[21]+alphaR[9]*fUpR[20])+0.2500000000000001*(alphaR[1]*fUpR[19]+fUpR[12]*alphaR[18]+alphaR[4]*fUpR[17])+0.25*fUpR[8]*alphaR[16]+0.223606797749979*(alphaR[7]*fUpR[15]+alphaR[3]*fUpR[13]); 
  GhatR[22] = 0.2500000000000001*alphaR[1]*fUpR[23]+0.159719141249985*alphaR[16]*fUpR[22]+0.25*(alphaR[0]*fUpR[22]+alphaR[5]*fUpR[21])+0.159719141249985*alphaR[18]*fUpR[19]+0.2500000000000001*(alphaR[2]*fUpR[19]+alphaR[4]*fUpR[18]+fUpR[4]*alphaR[18])+0.25*(alphaR[9]*fUpR[16]+fUpR[9]*alphaR[16])+0.223606797749979*(alphaR[3]*fUpR[14]+alphaR[7]*fUpR[10]); 
  GhatR[23] = (0.159719141249985*alphaR[16]+0.25*alphaR[0])*fUpR[23]+0.2500000000000001*alphaR[1]*fUpR[22]+0.159719141249985*alphaR[18]*fUpR[21]+0.2500000000000001*(alphaR[2]*fUpR[21]+alphaR[4]*fUpR[20])+0.25*(alphaR[5]*fUpR[19]+fUpR[8]*alphaR[18]+alphaR[9]*fUpR[17])+0.2500000000000001*fUpR[12]*alphaR[16]+0.223606797749979*(alphaR[3]*fUpR[15]+alphaR[7]*fUpR[13]); 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdx2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdx2; 
  out[5] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdx2; 
  out[6] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 
  out[7] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdx2; 
  out[8] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdx2; 
  out[9] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdx2; 
  out[10] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdx2; 
  out[11] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdx2; 
  out[12] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdx2; 
  out[13] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdx2; 
  out[14] += (0.7071067811865475*GhatL[9]-0.7071067811865475*GhatR[9])*rdx2; 
  out[15] += (0.7071067811865475*GhatL[10]-0.7071067811865475*GhatR[10])*rdx2; 
  out[16] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdx2; 
  out[17] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdx2; 
  out[18] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdx2; 
  out[19] += (0.7071067811865475*GhatL[11]-0.7071067811865475*GhatR[11])*rdx2; 
  out[20] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdx2; 
  out[21] += ((-1.224744871391589*GhatR[9])-1.224744871391589*GhatL[9])*rdx2; 
  out[22] += (0.7071067811865475*GhatL[12]-0.7071067811865475*GhatR[12])*rdx2; 
  out[23] += ((-1.224744871391589*GhatR[10])-1.224744871391589*GhatL[10])*rdx2; 
  out[24] += (0.7071067811865475*GhatL[13]-0.7071067811865475*GhatR[13])*rdx2; 
  out[25] += (0.7071067811865475*GhatL[14]-0.7071067811865475*GhatR[14])*rdx2; 
  out[26] += ((-1.224744871391589*GhatR[11])-1.224744871391589*GhatL[11])*rdx2; 
  out[27] += ((-1.224744871391589*GhatR[12])-1.224744871391589*GhatL[12])*rdx2; 
  out[28] += ((-1.224744871391589*GhatR[13])-1.224744871391589*GhatL[13])*rdx2; 
  out[29] += ((-1.224744871391589*GhatR[14])-1.224744871391589*GhatL[14])*rdx2; 
  out[30] += (0.7071067811865475*GhatL[15]-0.7071067811865475*GhatR[15])*rdx2; 
  out[31] += ((-1.224744871391589*GhatR[15])-1.224744871391589*GhatL[15])*rdx2; 
  out[32] += (0.7071067811865475*GhatL[16]-0.7071067811865475*GhatR[16])*rdx2; 
  out[33] += ((-1.224744871391589*GhatR[16])-1.224744871391589*GhatL[16])*rdx2; 
  out[34] += (0.7071067811865475*GhatL[17]-0.7071067811865475*GhatR[17])*rdx2; 
  out[35] += (0.7071067811865475*GhatL[18]-0.7071067811865475*GhatR[18])*rdx2; 
  out[36] += (0.7071067811865475*GhatL[19]-0.7071067811865475*GhatR[19])*rdx2; 
  out[37] += ((-1.224744871391589*GhatR[17])-1.224744871391589*GhatL[17])*rdx2; 
  out[38] += ((-1.224744871391589*GhatR[18])-1.224744871391589*GhatL[18])*rdx2; 
  out[39] += (0.7071067811865475*GhatL[20]-0.7071067811865475*GhatR[20])*rdx2; 
  out[40] += ((-1.224744871391589*GhatR[19])-1.224744871391589*GhatL[19])*rdx2; 
  out[41] += (0.7071067811865475*GhatL[21]-0.7071067811865475*GhatR[21])*rdx2; 
  out[42] += (0.7071067811865475*GhatL[22]-0.7071067811865475*GhatR[22])*rdx2; 
  out[43] += ((-1.224744871391589*GhatR[20])-1.224744871391589*GhatL[20])*rdx2; 
  out[44] += ((-1.224744871391589*GhatR[21])-1.224744871391589*GhatL[21])*rdx2; 
  out[45] += ((-1.224744871391589*GhatR[22])-1.224744871391589*GhatL[22])*rdx2; 
  out[46] += (0.7071067811865475*GhatL[23]-0.7071067811865475*GhatR[23])*rdx2; 
  out[47] += ((-1.224744871391589*GhatR[23])-1.224744871391589*GhatL[23])*rdx2; 

  return 1.5*rdx2*cflFreq; 

} 
