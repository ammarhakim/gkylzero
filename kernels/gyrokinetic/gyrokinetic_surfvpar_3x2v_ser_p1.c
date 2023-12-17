#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double *alpha_surf_l, const double *alpha_surf_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  const double *alphaL = &alpha_surf_l[72];
  const double *alphaR = &alpha_surf_r[72];
  double cflFreq = 0.0;
  double fUpOrdL[16] = {0.};
  double alphaL_n = 0.;

  alphaL_n = (-0.25*alphaL[13])-0.25*alphaL[11]+0.25*alphaL[10]+0.25*alphaL[8]+0.25*alphaL[7]+0.25*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_r(fl); 
  } else { 
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]-0.25*alphaL[11]-0.25*alphaL[10]-0.25*alphaL[8]+0.25*alphaL[7]+0.25*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_r(fl); 
  } else { 
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]+0.25*alphaL[11]-0.25*alphaL[10]+0.25*alphaL[8]-0.25*alphaL[7]-0.25*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_r(fl); 
  } else { 
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])+0.25*alphaL[11]+0.25*alphaL[10]-0.25*alphaL[8]-0.25*alphaL[7]-0.25*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_r(fl); 
  } else { 
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])+0.25*alphaL[11]+0.25*alphaL[10]+0.25*alphaL[8]-0.25*alphaL[7]+0.25*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_r(fl); 
  } else { 
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]+0.25*alphaL[11]-0.25*alphaL[10]-0.25*alphaL[8]-0.25*alphaL[7]+0.25*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_r(fl); 
  } else { 
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]-0.25*alphaL[11]-0.25*alphaL[10]+0.25*alphaL[8]+0.25*alphaL[7]-0.25*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_r(fl); 
  } else { 
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])-0.25*alphaL[11]+0.25*alphaL[10]-0.25*alphaL[8]+0.25*alphaL[7]-0.25*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_r(fl); 
  } else { 
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]+0.25*alphaL[11]+0.25*alphaL[10]-0.25*alphaL[8]+0.25*alphaL[7]-0.25*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_r(fl); 
  } else { 
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])+0.25*alphaL[11]-0.25*alphaL[10]+0.25*alphaL[8]+0.25*alphaL[7]-0.25*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_r(fl); 
  } else { 
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])-0.25*alphaL[11]-0.25*alphaL[10]-0.25*alphaL[8]-0.25*alphaL[7]+0.25*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_r(fl); 
  } else { 
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]-0.25*alphaL[11]+0.25*alphaL[10]+0.25*alphaL[8]-0.25*alphaL[7]+0.25*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_r(fl); 
  } else { 
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]-0.25*alphaL[11]+0.25*alphaL[10]-0.25*alphaL[8]-0.25*alphaL[7]-0.25*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_r(fl); 
  } else { 
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])-0.25*alphaL[11]-0.25*alphaL[10]+0.25*alphaL[8]-0.25*alphaL[7]-0.25*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_r(fl); 
  } else { 
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])+0.25*alphaL[11]-0.25*alphaL[10]-0.25*alphaL[8]+0.25*alphaL[7]+0.25*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_r(fl); 
  } else { 
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]+0.25*alphaL[11]+0.25*alphaL[10]+0.25*alphaL[8]+0.25*alphaL[7]+0.25*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_r(fl); 
  } else { 
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[16] = {0.};
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[16] = {0.}; 
  GhatL[0] = 0.25*(alphaL[13]*fUpL[13]+alphaL[11]*fUpL[11]+alphaL[10]*fUpL[10]+alphaL[8]*fUpL[8]+alphaL[7]*fUpL[7]+alphaL[6]*fUpL[6]+alphaL[5]*fUpL[5]+alphaL[4]*fUpL[4]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.25*(alphaL[10]*fUpL[13]+fUpL[10]*alphaL[13]+alphaL[7]*fUpL[11]+fUpL[7]*alphaL[11]+alphaL[4]*fUpL[8]+fUpL[4]*alphaL[8]+alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6]+alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.25*(alphaL[13]*fUpL[15]+alphaL[10]*fUpL[14]+alphaL[8]*fUpL[12]+alphaL[6]*fUpL[11]+fUpL[6]*alphaL[11]+alphaL[4]*fUpL[9]+alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7]+alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.25*(alphaL[8]*fUpL[13]+fUpL[8]*alphaL[13]+alphaL[5]*fUpL[11]+fUpL[5]*alphaL[11]+alphaL[4]*fUpL[10]+fUpL[4]*alphaL[10]+alphaL[2]*fUpL[7]+fUpL[2]*alphaL[7]+alphaL[1]*fUpL[6]+fUpL[1]*alphaL[6]+alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]); 
  GhatL[4] = 0.25*(alphaL[11]*fUpL[15]+alphaL[7]*fUpL[14]+alphaL[6]*fUpL[13]+fUpL[6]*alphaL[13]+alphaL[5]*fUpL[12]+alphaL[3]*fUpL[10]+fUpL[3]*alphaL[10]+alphaL[2]*fUpL[9]+alphaL[1]*fUpL[8]+fUpL[1]*alphaL[8]+alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]); 
  GhatL[5] = 0.25*(alphaL[10]*fUpL[15]+alphaL[13]*fUpL[14]+alphaL[4]*fUpL[12]+alphaL[3]*fUpL[11]+fUpL[3]*alphaL[11]+alphaL[8]*fUpL[9]+alphaL[6]*fUpL[7]+fUpL[6]*alphaL[7]+alphaL[0]*fUpL[5]+fUpL[0]*alphaL[5]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[6] = 0.25*(alphaL[4]*fUpL[13]+fUpL[4]*alphaL[13]+alphaL[2]*fUpL[11]+fUpL[2]*alphaL[11]+alphaL[8]*fUpL[10]+fUpL[8]*alphaL[10]+alphaL[5]*fUpL[7]+fUpL[5]*alphaL[7]+alphaL[0]*fUpL[6]+fUpL[0]*alphaL[6]+alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[7] = 0.25*(alphaL[8]*fUpL[15]+alphaL[4]*fUpL[14]+fUpL[12]*alphaL[13]+alphaL[1]*fUpL[11]+fUpL[1]*alphaL[11]+fUpL[9]*alphaL[10]+alphaL[0]*fUpL[7]+fUpL[0]*alphaL[7]+alphaL[5]*fUpL[6]+fUpL[5]*alphaL[6]+alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 
  GhatL[8] = 0.25*(alphaL[7]*fUpL[15]+alphaL[11]*fUpL[14]+alphaL[3]*fUpL[13]+fUpL[3]*alphaL[13]+alphaL[2]*fUpL[12]+alphaL[6]*fUpL[10]+fUpL[6]*alphaL[10]+alphaL[5]*fUpL[9]+alphaL[0]*fUpL[8]+fUpL[0]*alphaL[8]+alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]); 
  GhatL[9] = 0.25*(alphaL[6]*fUpL[15]+alphaL[3]*fUpL[14]+alphaL[11]*fUpL[13]+fUpL[11]*alphaL[13]+alphaL[1]*fUpL[12]+alphaL[7]*fUpL[10]+fUpL[7]*alphaL[10]+alphaL[0]*fUpL[9]+alphaL[5]*fUpL[8]+fUpL[5]*alphaL[8]+alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]); 
  GhatL[10] = 0.25*(alphaL[5]*fUpL[15]+alphaL[2]*fUpL[14]+alphaL[1]*fUpL[13]+fUpL[1]*alphaL[13]+alphaL[11]*fUpL[12]+alphaL[0]*fUpL[10]+fUpL[0]*alphaL[10]+alphaL[7]*fUpL[9]+alphaL[6]*fUpL[8]+fUpL[6]*alphaL[8]+alphaL[3]*fUpL[4]+fUpL[3]*alphaL[4]); 
  GhatL[11] = 0.25*(alphaL[4]*fUpL[15]+alphaL[8]*fUpL[14]+fUpL[9]*alphaL[13]+alphaL[10]*fUpL[12]+alphaL[0]*fUpL[11]+fUpL[0]*alphaL[11]+alphaL[1]*fUpL[7]+fUpL[1]*alphaL[7]+alphaL[2]*fUpL[6]+fUpL[2]*alphaL[6]+alphaL[3]*fUpL[5]+fUpL[3]*alphaL[5]); 
  GhatL[12] = 0.25*(alphaL[3]*fUpL[15]+alphaL[6]*fUpL[14]+alphaL[7]*fUpL[13]+fUpL[7]*alphaL[13]+alphaL[0]*fUpL[12]+alphaL[10]*fUpL[11]+fUpL[10]*alphaL[11]+alphaL[1]*fUpL[9]+alphaL[2]*fUpL[8]+fUpL[2]*alphaL[8]+alphaL[4]*fUpL[5]+fUpL[4]*alphaL[5]); 
  GhatL[13] = 0.25*(alphaL[2]*fUpL[15]+alphaL[5]*fUpL[14]+alphaL[0]*fUpL[13]+fUpL[0]*alphaL[13]+alphaL[7]*fUpL[12]+fUpL[9]*alphaL[11]+alphaL[1]*fUpL[10]+fUpL[1]*alphaL[10]+alphaL[3]*fUpL[8]+fUpL[3]*alphaL[8]+alphaL[4]*fUpL[6]+fUpL[4]*alphaL[6]); 
  GhatL[14] = 0.25*(alphaL[1]*fUpL[15]+alphaL[0]*fUpL[14]+alphaL[5]*fUpL[13]+fUpL[5]*alphaL[13]+alphaL[6]*fUpL[12]+alphaL[8]*fUpL[11]+fUpL[8]*alphaL[11]+alphaL[2]*fUpL[10]+fUpL[2]*alphaL[10]+alphaL[3]*fUpL[9]+alphaL[4]*fUpL[7]+fUpL[4]*alphaL[7]); 
  GhatL[15] = 0.25*(alphaL[0]*fUpL[15]+alphaL[1]*fUpL[14]+alphaL[2]*fUpL[13]+fUpL[2]*alphaL[13]+alphaL[3]*fUpL[12]+alphaL[4]*fUpL[11]+fUpL[4]*alphaL[11]+alphaL[5]*fUpL[10]+fUpL[5]*alphaL[10]+alphaL[6]*fUpL[9]+alphaL[7]*fUpL[8]+fUpL[7]*alphaL[8]); 

  double fUpOrdR[16] = {0.};
  double alphaR_n = 0.;

  alphaR_n = (-0.25*alphaR[13])-0.25*alphaR[11]+0.25*alphaR[10]+0.25*alphaR[8]+0.25*alphaR[7]+0.25*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_r(fc); 
  } else { 
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]-0.25*alphaR[11]-0.25*alphaR[10]-0.25*alphaR[8]+0.25*alphaR[7]+0.25*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_r(fc); 
  } else { 
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]+0.25*alphaR[11]-0.25*alphaR[10]+0.25*alphaR[8]-0.25*alphaR[7]-0.25*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_r(fc); 
  } else { 
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])+0.25*alphaR[11]+0.25*alphaR[10]-0.25*alphaR[8]-0.25*alphaR[7]-0.25*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_r(fc); 
  } else { 
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])+0.25*alphaR[11]+0.25*alphaR[10]+0.25*alphaR[8]-0.25*alphaR[7]+0.25*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_r(fc); 
  } else { 
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]+0.25*alphaR[11]-0.25*alphaR[10]-0.25*alphaR[8]-0.25*alphaR[7]+0.25*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_r(fc); 
  } else { 
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]-0.25*alphaR[11]-0.25*alphaR[10]+0.25*alphaR[8]+0.25*alphaR[7]-0.25*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_r(fc); 
  } else { 
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])-0.25*alphaR[11]+0.25*alphaR[10]-0.25*alphaR[8]+0.25*alphaR[7]-0.25*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_r(fc); 
  } else { 
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]+0.25*alphaR[11]+0.25*alphaR[10]-0.25*alphaR[8]+0.25*alphaR[7]-0.25*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_r(fc); 
  } else { 
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])+0.25*alphaR[11]-0.25*alphaR[10]+0.25*alphaR[8]+0.25*alphaR[7]-0.25*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_r(fc); 
  } else { 
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])-0.25*alphaR[11]-0.25*alphaR[10]-0.25*alphaR[8]-0.25*alphaR[7]+0.25*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_r(fc); 
  } else { 
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]-0.25*alphaR[11]+0.25*alphaR[10]+0.25*alphaR[8]-0.25*alphaR[7]+0.25*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_r(fc); 
  } else { 
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]-0.25*alphaR[11]+0.25*alphaR[10]-0.25*alphaR[8]-0.25*alphaR[7]-0.25*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_r(fc); 
  } else { 
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])-0.25*alphaR[11]-0.25*alphaR[10]+0.25*alphaR[8]-0.25*alphaR[7]-0.25*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_r(fc); 
  } else { 
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])+0.25*alphaR[11]-0.25*alphaR[10]-0.25*alphaR[8]+0.25*alphaR[7]+0.25*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_r(fc); 
  } else { 
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]+0.25*alphaR[11]+0.25*alphaR[10]+0.25*alphaR[8]+0.25*alphaR[7]+0.25*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_r(fc); 
  } else { 
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[16] = {0.};
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[16] = {0.}; 
  GhatR[0] = 0.25*(alphaR[13]*fUpR[13]+alphaR[11]*fUpR[11]+alphaR[10]*fUpR[10]+alphaR[8]*fUpR[8]+alphaR[7]*fUpR[7]+alphaR[6]*fUpR[6]+alphaR[5]*fUpR[5]+alphaR[4]*fUpR[4]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.25*(alphaR[10]*fUpR[13]+fUpR[10]*alphaR[13]+alphaR[7]*fUpR[11]+fUpR[7]*alphaR[11]+alphaR[4]*fUpR[8]+fUpR[4]*alphaR[8]+alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6]+alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.25*(alphaR[13]*fUpR[15]+alphaR[10]*fUpR[14]+alphaR[8]*fUpR[12]+alphaR[6]*fUpR[11]+fUpR[6]*alphaR[11]+alphaR[4]*fUpR[9]+alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7]+alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.25*(alphaR[8]*fUpR[13]+fUpR[8]*alphaR[13]+alphaR[5]*fUpR[11]+fUpR[5]*alphaR[11]+alphaR[4]*fUpR[10]+fUpR[4]*alphaR[10]+alphaR[2]*fUpR[7]+fUpR[2]*alphaR[7]+alphaR[1]*fUpR[6]+fUpR[1]*alphaR[6]+alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]); 
  GhatR[4] = 0.25*(alphaR[11]*fUpR[15]+alphaR[7]*fUpR[14]+alphaR[6]*fUpR[13]+fUpR[6]*alphaR[13]+alphaR[5]*fUpR[12]+alphaR[3]*fUpR[10]+fUpR[3]*alphaR[10]+alphaR[2]*fUpR[9]+alphaR[1]*fUpR[8]+fUpR[1]*alphaR[8]+alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]); 
  GhatR[5] = 0.25*(alphaR[10]*fUpR[15]+alphaR[13]*fUpR[14]+alphaR[4]*fUpR[12]+alphaR[3]*fUpR[11]+fUpR[3]*alphaR[11]+alphaR[8]*fUpR[9]+alphaR[6]*fUpR[7]+fUpR[6]*alphaR[7]+alphaR[0]*fUpR[5]+fUpR[0]*alphaR[5]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[6] = 0.25*(alphaR[4]*fUpR[13]+fUpR[4]*alphaR[13]+alphaR[2]*fUpR[11]+fUpR[2]*alphaR[11]+alphaR[8]*fUpR[10]+fUpR[8]*alphaR[10]+alphaR[5]*fUpR[7]+fUpR[5]*alphaR[7]+alphaR[0]*fUpR[6]+fUpR[0]*alphaR[6]+alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[7] = 0.25*(alphaR[8]*fUpR[15]+alphaR[4]*fUpR[14]+fUpR[12]*alphaR[13]+alphaR[1]*fUpR[11]+fUpR[1]*alphaR[11]+fUpR[9]*alphaR[10]+alphaR[0]*fUpR[7]+fUpR[0]*alphaR[7]+alphaR[5]*fUpR[6]+fUpR[5]*alphaR[6]+alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 
  GhatR[8] = 0.25*(alphaR[7]*fUpR[15]+alphaR[11]*fUpR[14]+alphaR[3]*fUpR[13]+fUpR[3]*alphaR[13]+alphaR[2]*fUpR[12]+alphaR[6]*fUpR[10]+fUpR[6]*alphaR[10]+alphaR[5]*fUpR[9]+alphaR[0]*fUpR[8]+fUpR[0]*alphaR[8]+alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]); 
  GhatR[9] = 0.25*(alphaR[6]*fUpR[15]+alphaR[3]*fUpR[14]+alphaR[11]*fUpR[13]+fUpR[11]*alphaR[13]+alphaR[1]*fUpR[12]+alphaR[7]*fUpR[10]+fUpR[7]*alphaR[10]+alphaR[0]*fUpR[9]+alphaR[5]*fUpR[8]+fUpR[5]*alphaR[8]+alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]); 
  GhatR[10] = 0.25*(alphaR[5]*fUpR[15]+alphaR[2]*fUpR[14]+alphaR[1]*fUpR[13]+fUpR[1]*alphaR[13]+alphaR[11]*fUpR[12]+alphaR[0]*fUpR[10]+fUpR[0]*alphaR[10]+alphaR[7]*fUpR[9]+alphaR[6]*fUpR[8]+fUpR[6]*alphaR[8]+alphaR[3]*fUpR[4]+fUpR[3]*alphaR[4]); 
  GhatR[11] = 0.25*(alphaR[4]*fUpR[15]+alphaR[8]*fUpR[14]+fUpR[9]*alphaR[13]+alphaR[10]*fUpR[12]+alphaR[0]*fUpR[11]+fUpR[0]*alphaR[11]+alphaR[1]*fUpR[7]+fUpR[1]*alphaR[7]+alphaR[2]*fUpR[6]+fUpR[2]*alphaR[6]+alphaR[3]*fUpR[5]+fUpR[3]*alphaR[5]); 
  GhatR[12] = 0.25*(alphaR[3]*fUpR[15]+alphaR[6]*fUpR[14]+alphaR[7]*fUpR[13]+fUpR[7]*alphaR[13]+alphaR[0]*fUpR[12]+alphaR[10]*fUpR[11]+fUpR[10]*alphaR[11]+alphaR[1]*fUpR[9]+alphaR[2]*fUpR[8]+fUpR[2]*alphaR[8]+alphaR[4]*fUpR[5]+fUpR[4]*alphaR[5]); 
  GhatR[13] = 0.25*(alphaR[2]*fUpR[15]+alphaR[5]*fUpR[14]+alphaR[0]*fUpR[13]+fUpR[0]*alphaR[13]+alphaR[7]*fUpR[12]+fUpR[9]*alphaR[11]+alphaR[1]*fUpR[10]+fUpR[1]*alphaR[10]+alphaR[3]*fUpR[8]+fUpR[3]*alphaR[8]+alphaR[4]*fUpR[6]+fUpR[4]*alphaR[6]); 
  GhatR[14] = 0.25*(alphaR[1]*fUpR[15]+alphaR[0]*fUpR[14]+alphaR[5]*fUpR[13]+fUpR[5]*alphaR[13]+alphaR[6]*fUpR[12]+alphaR[8]*fUpR[11]+fUpR[8]*alphaR[11]+alphaR[2]*fUpR[10]+fUpR[2]*alphaR[10]+alphaR[3]*fUpR[9]+alphaR[4]*fUpR[7]+fUpR[4]*alphaR[7]); 
  GhatR[15] = 0.25*(alphaR[0]*fUpR[15]+alphaR[1]*fUpR[14]+alphaR[2]*fUpR[13]+fUpR[2]*alphaR[13]+alphaR[3]*fUpR[12]+alphaR[4]*fUpR[11]+fUpR[4]*alphaR[11]+alphaR[5]*fUpR[10]+fUpR[5]*alphaR[10]+alphaR[6]*fUpR[9]+alphaR[7]*fUpR[8]+fUpR[7]*alphaR[8]); 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdvpar2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdvpar2; 
  out[2] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdvpar2; 
  out[3] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdvpar2; 
  out[4] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdvpar2; 
  out[5] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdvpar2; 
  out[6] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdvpar2; 
  out[7] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdvpar2; 
  out[8] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdvpar2; 
  out[9] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdvpar2; 
  out[10] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdvpar2; 
  out[11] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdvpar2; 
  out[12] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdvpar2; 
  out[13] += (0.7071067811865475*GhatL[9]-0.7071067811865475*GhatR[9])*rdvpar2; 
  out[14] += (0.7071067811865475*GhatL[10]-0.7071067811865475*GhatR[10])*rdvpar2; 
  out[15] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdvpar2; 
  out[16] += (0.7071067811865475*GhatL[11]-0.7071067811865475*GhatR[11])*rdvpar2; 
  out[17] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdvpar2; 
  out[18] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdvpar2; 
  out[19] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdvpar2; 
  out[20] += (0.7071067811865475*GhatL[12]-0.7071067811865475*GhatR[12])*rdvpar2; 
  out[21] += (0.7071067811865475*GhatL[13]-0.7071067811865475*GhatR[13])*rdvpar2; 
  out[22] += (0.7071067811865475*GhatL[14]-0.7071067811865475*GhatR[14])*rdvpar2; 
  out[23] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdvpar2; 
  out[24] += ((-1.224744871391589*GhatR[9])-1.224744871391589*GhatL[9])*rdvpar2; 
  out[25] += ((-1.224744871391589*GhatR[10])-1.224744871391589*GhatL[10])*rdvpar2; 
  out[26] += ((-1.224744871391589*GhatR[11])-1.224744871391589*GhatL[11])*rdvpar2; 
  out[27] += (0.7071067811865475*GhatL[15]-0.7071067811865475*GhatR[15])*rdvpar2; 
  out[28] += ((-1.224744871391589*GhatR[12])-1.224744871391589*GhatL[12])*rdvpar2; 
  out[29] += ((-1.224744871391589*GhatR[13])-1.224744871391589*GhatL[13])*rdvpar2; 
  out[30] += ((-1.224744871391589*GhatR[14])-1.224744871391589*GhatL[14])*rdvpar2; 
  out[31] += ((-1.224744871391589*GhatR[15])-1.224744871391589*GhatL[15])*rdvpar2; 
  out[32] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdvpar2; 
  out[33] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdvpar2; 
  out[34] += (1.58113883008419*GhatL[2]-1.58113883008419*GhatR[2])*rdvpar2; 
  out[35] += (1.58113883008419*GhatL[3]-1.58113883008419*GhatR[3])*rdvpar2; 
  out[36] += (1.58113883008419*GhatL[4]-1.58113883008419*GhatR[4])*rdvpar2; 
  out[37] += (1.58113883008419*GhatL[5]-1.58113883008419*GhatR[5])*rdvpar2; 
  out[38] += (1.58113883008419*GhatL[6]-1.58113883008419*GhatR[6])*rdvpar2; 
  out[39] += (1.58113883008419*GhatL[7]-1.58113883008419*GhatR[7])*rdvpar2; 
  out[40] += (1.58113883008419*GhatL[8]-1.58113883008419*GhatR[8])*rdvpar2; 
  out[41] += (1.58113883008419*GhatL[9]-1.58113883008419*GhatR[9])*rdvpar2; 
  out[42] += (1.58113883008419*GhatL[10]-1.58113883008419*GhatR[10])*rdvpar2; 
  out[43] += (1.58113883008419*GhatL[11]-1.58113883008419*GhatR[11])*rdvpar2; 
  out[44] += (1.58113883008419*GhatL[12]-1.58113883008419*GhatR[12])*rdvpar2; 
  out[45] += (1.58113883008419*GhatL[13]-1.58113883008419*GhatR[13])*rdvpar2; 
  out[46] += (1.58113883008419*GhatL[14]-1.58113883008419*GhatR[14])*rdvpar2; 
  out[47] += (1.58113883008419*GhatL[15]-1.58113883008419*GhatR[15])*rdvpar2; 

  return 2.5*rdvpar2*cflFreq; 

} 
