#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_ser_4x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double *alpha_surf_l, const double *alpha_surf_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  const double *alphaL = &alpha_surf_l[20];
  const double *alphaR = &alpha_surf_r[20];
  double cflFreq = 0.0;
  double fUpOrdL[27] = {0.};
  double alphaL_n = 0.;

  alphaL_n = (-0.4242640687119286*alphaL[13])-0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[5]+0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = ser_4x_p2_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpOrdL[0] = ser_4x_p2_surfx2_eval_quad_node_0_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4242640687119282*alphaL[12])-0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = ser_4x_p2_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpOrdL[1] = ser_4x_p2_surfx2_eval_quad_node_1_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4242640687119286*alphaL[13]-0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[5]+0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = ser_4x_p2_surfx2_eval_quad_node_2_r(fl); 
  } else { 
    fUpOrdL[2] = ser_4x_p2_surfx2_eval_quad_node_2_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4242640687119286*alphaL[13])+0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[5]-0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = ser_4x_p2_surfx2_eval_quad_node_3_r(fl); 
  } else { 
    fUpOrdL[3] = ser_4x_p2_surfx2_eval_quad_node_3_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = ser_4x_p2_surfx2_eval_quad_node_4_r(fl); 
  } else { 
    fUpOrdL[4] = ser_4x_p2_surfx2_eval_quad_node_4_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4242640687119286*alphaL[13]+0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[5]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = ser_4x_p2_surfx2_eval_quad_node_5_r(fl); 
  } else { 
    fUpOrdL[5] = ser_4x_p2_surfx2_eval_quad_node_5_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4242640687119286*alphaL[13])-0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = ser_4x_p2_surfx2_eval_quad_node_6_r(fl); 
  } else { 
    fUpOrdL[6] = ser_4x_p2_surfx2_eval_quad_node_6_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4242640687119282*alphaL[12])+0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = ser_4x_p2_surfx2_eval_quad_node_7_r(fl); 
  } else { 
    fUpOrdL[7] = ser_4x_p2_surfx2_eval_quad_node_7_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4242640687119286*alphaL[13]-0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]-0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = ser_4x_p2_surfx2_eval_quad_node_8_r(fl); 
  } else { 
    fUpOrdL[8] = ser_4x_p2_surfx2_eval_quad_node_8_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5303300858899102*alphaL[13]+0.5303300858899102*alphaL[11]+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]-0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[9] = ser_4x_p2_surfx2_eval_quad_node_9_r(fl); 
  } else { 
    fUpOrdL[9] = ser_4x_p2_surfx2_eval_quad_node_9_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5303300858899102*alphaL[11]+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[10] = ser_4x_p2_surfx2_eval_quad_node_10_r(fl); 
  } else { 
    fUpOrdL[10] = ser_4x_p2_surfx2_eval_quad_node_10_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5303300858899102*alphaL[13])+0.5303300858899102*alphaL[11]+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[11] = ser_4x_p2_surfx2_eval_quad_node_11_r(fl); 
  } else { 
    fUpOrdL[11] = ser_4x_p2_surfx2_eval_quad_node_11_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5303300858899102*alphaL[13]-0.3952847075210471*alphaL[8]-0.3952847075210471*alphaL[7]-0.4743416490252568*alphaL[3]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[12] = ser_4x_p2_surfx2_eval_quad_node_12_r(fl); 
  } else { 
    fUpOrdL[12] = ser_4x_p2_surfx2_eval_quad_node_12_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.3952847075210471*alphaL[8])-0.3952847075210471*alphaL[7]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[13] = ser_4x_p2_surfx2_eval_quad_node_13_r(fl); 
  } else { 
    fUpOrdL[13] = ser_4x_p2_surfx2_eval_quad_node_13_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5303300858899102*alphaL[13])-0.3952847075210471*alphaL[8]-0.3952847075210471*alphaL[7]+0.4743416490252568*alphaL[3]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[14] = ser_4x_p2_surfx2_eval_quad_node_14_r(fl); 
  } else { 
    fUpOrdL[14] = ser_4x_p2_surfx2_eval_quad_node_14_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5303300858899102*alphaL[13]-0.5303300858899102*alphaL[11]+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[15] = ser_4x_p2_surfx2_eval_quad_node_15_r(fl); 
  } else { 
    fUpOrdL[15] = ser_4x_p2_surfx2_eval_quad_node_15_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5303300858899102*alphaL[11])+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[16] = ser_4x_p2_surfx2_eval_quad_node_16_r(fl); 
  } else { 
    fUpOrdL[16] = ser_4x_p2_surfx2_eval_quad_node_16_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5303300858899102*alphaL[13])-0.5303300858899102*alphaL[11]+0.3162277660168378*alphaL[8]-0.3952847075210471*alphaL[7]+0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[17] = ser_4x_p2_surfx2_eval_quad_node_17_r(fl); 
  } else { 
    fUpOrdL[17] = ser_4x_p2_surfx2_eval_quad_node_17_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4242640687119286*alphaL[13])+0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[18] = ser_4x_p2_surfx2_eval_quad_node_18_r(fl); 
  } else { 
    fUpOrdL[18] = ser_4x_p2_surfx2_eval_quad_node_18_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[19] = ser_4x_p2_surfx2_eval_quad_node_19_r(fl); 
  } else { 
    fUpOrdL[19] = ser_4x_p2_surfx2_eval_quad_node_19_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4242640687119286*alphaL[13]+0.4242640687119282*alphaL[12]-0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[5]-0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]-0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[20] = ser_4x_p2_surfx2_eval_quad_node_20_r(fl); 
  } else { 
    fUpOrdL[20] = ser_4x_p2_surfx2_eval_quad_node_20_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4242640687119286*alphaL[13])-0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[5]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[21] = ser_4x_p2_surfx2_eval_quad_node_21_r(fl); 
  } else { 
    fUpOrdL[21] = ser_4x_p2_surfx2_eval_quad_node_21_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5303300858899102*alphaL[12])-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[22] = ser_4x_p2_surfx2_eval_quad_node_22_r(fl); 
  } else { 
    fUpOrdL[22] = ser_4x_p2_surfx2_eval_quad_node_22_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4242640687119286*alphaL[13]-0.5303300858899102*alphaL[12]-0.3952847075210471*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[5]+0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[23] = ser_4x_p2_surfx2_eval_quad_node_23_r(fl); 
  } else { 
    fUpOrdL[23] = ser_4x_p2_surfx2_eval_quad_node_23_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.4242640687119286*alphaL[13])+0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]-0.6363961030678927*alphaL[5]+0.6363961030678927*alphaL[4]-0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[24] = ser_4x_p2_surfx2_eval_quad_node_24_r(fl); 
  } else { 
    fUpOrdL[24] = ser_4x_p2_surfx2_eval_quad_node_24_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[25] = ser_4x_p2_surfx2_eval_quad_node_25_r(fl); 
  } else { 
    fUpOrdL[25] = ser_4x_p2_surfx2_eval_quad_node_25_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4242640687119286*alphaL[13]+0.4242640687119282*alphaL[12]+0.4242640687119286*alphaL[11]+0.3162277660168378*alphaL[8]+0.3162277660168378*alphaL[7]+0.6363961030678927*alphaL[5]+0.6363961030678927*alphaL[4]+0.4743416490252568*alphaL[3]+0.4743416490252568*alphaL[2]+0.4743416490252568*alphaL[1]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[26] = ser_4x_p2_surfx2_eval_quad_node_26_r(fl); 
  } else { 
    fUpOrdL[26] = ser_4x_p2_surfx2_eval_quad_node_26_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[20] = {0.};
  ser_4x_p2_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[20] = {0.}; 
  GhatL[0] = 0.3535533905932737*(alphaL[13]*fUpL[13]+alphaL[12]*fUpL[12]+alphaL[11]*fUpL[11]+alphaL[8]*fUpL[8]+alphaL[7]*fUpL[7]+alphaL[5]*fUpL[5]+alphaL[4]*fUpL[4]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.3162277660168379*(alphaL[5]*fUpL[13]+fUpL[5]*alphaL[13])+0.3535533905932737*(alphaL[8]*fUpL[12]+fUpL[8]*alphaL[12])+0.3162277660168379*(alphaL[4]*fUpL[11]+fUpL[4]*alphaL[11])+0.3162277660168379*(alphaL[1]*fUpL[7]+fUpL[1]*alphaL[7])+0.3535533905932737*(alphaL[3]*fUpL[5]+fUpL[3]*alphaL[5]+alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.3535533905932737*alphaL[13]*fUpL[17]+0.3162277660168379*(alphaL[4]*fUpL[12]+fUpL[4]*alphaL[12])+0.3535533905932737*(alphaL[7]*fUpL[11]+fUpL[7]*alphaL[11]+alphaL[5]*fUpL[10])+0.3162277660168379*(alphaL[2]*fUpL[8]+fUpL[2]*alphaL[8])+0.3535533905932737*(alphaL[3]*fUpL[6]+alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.3535533905932737*(alphaL[12]*fUpL[18]+alphaL[11]*fUpL[17])+0.3162277660168379*alphaL[5]*fUpL[15]+0.3535533905932737*(alphaL[8]*fUpL[14]+alphaL[7]*fUpL[13]+fUpL[7]*alphaL[13]+alphaL[4]*fUpL[10])+0.3162277660168379*alphaL[3]*fUpL[9]+0.3535533905932737*(alphaL[2]*fUpL[6]+alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]+alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]); 
  GhatL[4] = 0.3162277660168379*alphaL[5]*fUpL[17]+0.3162277660168379*fUpL[10]*alphaL[13]+(0.2828427124746191*alphaL[11]+0.3162277660168379*alphaL[2])*fUpL[12]+0.2828427124746191*fUpL[11]*alphaL[12]+0.3162277660168379*(fUpL[2]*alphaL[12]+alphaL[1]*fUpL[11]+fUpL[1]*alphaL[11])+0.3535533905932737*alphaL[3]*fUpL[10]+0.3162277660168379*(alphaL[4]*fUpL[8]+fUpL[4]*alphaL[8]+alphaL[4]*fUpL[7]+fUpL[4]*alphaL[7])+0.3535533905932737*(alphaL[5]*fUpL[6]+alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[5] = 0.3535533905932737*alphaL[8]*fUpL[18]+0.3162277660168379*alphaL[4]*fUpL[17]+(0.2828427124746191*alphaL[13]+0.3162277660168379*alphaL[3])*fUpL[15]+0.3535533905932737*alphaL[12]*fUpL[14]+0.3162277660168379*(alphaL[1]*fUpL[13]+fUpL[1]*alphaL[13])+fUpL[10]*(0.3162277660168379*alphaL[11]+0.3535533905932737*alphaL[2])+0.3162277660168379*(alphaL[5]*(fUpL[9]+fUpL[7])+fUpL[5]*alphaL[7])+0.3535533905932737*(alphaL[4]*fUpL[6]+alphaL[0]*fUpL[5]+fUpL[0]*alphaL[5]+alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[6] = 0.3162277660168379*(alphaL[5]*fUpL[19]+alphaL[4]*fUpL[18])+0.3535533905932737*alphaL[7]*fUpL[17]+0.3162277660168379*(alphaL[3]*fUpL[16]+alphaL[2]*fUpL[14])+0.3535533905932737*(alphaL[11]*fUpL[13]+fUpL[11]*alphaL[13])+fUpL[10]*(0.3162277660168379*alphaL[12]+0.3535533905932737*alphaL[1])+0.3162277660168379*fUpL[6]*alphaL[8]+0.3535533905932737*(alphaL[0]*fUpL[6]+alphaL[4]*fUpL[5]+fUpL[4]*alphaL[5]+alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 
  GhatL[7] = 0.2258769757263128*alphaL[13]*fUpL[13]+0.3535533905932737*(alphaL[3]*fUpL[13]+fUpL[3]*alphaL[13])+0.3162277660168379*alphaL[12]*fUpL[12]+0.2258769757263128*alphaL[11]*fUpL[11]+0.3535533905932737*(alphaL[2]*fUpL[11]+fUpL[2]*alphaL[11])+0.2258769757263128*alphaL[7]*fUpL[7]+0.3535533905932737*(alphaL[0]*fUpL[7]+fUpL[0]*alphaL[7])+0.3162277660168379*(alphaL[5]*fUpL[5]+alphaL[4]*fUpL[4]+alphaL[1]*fUpL[1]); 
  GhatL[8] = 0.3535533905932737*(alphaL[5]*fUpL[18]+alphaL[3]*fUpL[14])+0.2258769757263128*alphaL[12]*fUpL[12]+0.3535533905932737*(alphaL[1]*fUpL[12]+fUpL[1]*alphaL[12])+0.3162277660168379*alphaL[11]*fUpL[11]+0.2258769757263128*alphaL[8]*fUpL[8]+0.3535533905932737*(alphaL[0]*fUpL[8]+fUpL[0]*alphaL[8])+0.3162277660168379*(alphaL[4]*fUpL[4]+alphaL[2]*fUpL[2]); 
  GhatL[9] = 0.3535533905932737*(alphaL[4]*fUpL[19]+alphaL[2]*fUpL[16]+alphaL[1]*fUpL[15])+0.3162277660168379*alphaL[13]*fUpL[13]+0.3535533905932737*alphaL[0]*fUpL[9]+0.3162277660168379*(alphaL[5]*fUpL[5]+alphaL[3]*fUpL[3]); 
  GhatL[10] = (0.282842712474619*alphaL[13]+0.3162277660168379*alphaL[3])*fUpL[19]+(0.282842712474619*alphaL[11]+0.3162277660168379*alphaL[2])*fUpL[18]+(0.282842712474619*alphaL[12]+0.3162277660168379*alphaL[1])*fUpL[17]+0.3162277660168379*(alphaL[5]*fUpL[16]+alphaL[4]*(fUpL[14]+fUpL[13])+fUpL[4]*alphaL[13]+fUpL[6]*alphaL[12]+alphaL[5]*fUpL[11]+fUpL[5]*alphaL[11])+0.3162277660168379*(alphaL[8]+alphaL[7])*fUpL[10]+0.3535533905932737*(alphaL[0]*fUpL[10]+alphaL[1]*fUpL[6]+alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5]+alphaL[3]*fUpL[4]+fUpL[3]*alphaL[4]); 
  GhatL[11] = 0.2258769757263128*alphaL[13]*fUpL[17]+0.3535533905932737*(alphaL[3]*fUpL[17]+fUpL[6]*alphaL[13])+0.2828427124746191*(alphaL[4]*fUpL[12]+fUpL[4]*alphaL[12])+(0.3162277660168379*alphaL[8]+0.2258769757263128*alphaL[7]+0.3535533905932737*alphaL[0])*fUpL[11]+(0.3162277660168379*fUpL[8]+0.2258769757263128*fUpL[7]+0.3535533905932737*fUpL[0])*alphaL[11]+0.3162277660168379*alphaL[5]*fUpL[10]+0.3535533905932737*(alphaL[2]*fUpL[7]+fUpL[2]*alphaL[7])+0.3162277660168379*(alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]); 
  GhatL[12] = 0.3162277660168379*alphaL[13]*fUpL[18]+0.3535533905932737*(alphaL[3]*fUpL[18]+alphaL[5]*fUpL[14])+(0.2258769757263128*alphaL[8]+0.3162277660168379*alphaL[7]+0.3535533905932737*alphaL[0])*fUpL[12]+(0.2258769757263128*fUpL[8]+0.3162277660168379*fUpL[7]+0.3535533905932737*fUpL[0])*alphaL[12]+0.2828427124746191*(alphaL[4]*fUpL[11]+fUpL[4]*alphaL[11])+0.3535533905932737*(alphaL[1]*fUpL[8]+fUpL[1]*alphaL[8])+0.3162277660168379*(alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]); 
  GhatL[13] = 0.3162277660168379*alphaL[12]*fUpL[18]+(0.2258769757263128*alphaL[11]+0.3535533905932737*alphaL[2])*fUpL[17]+0.2828427124746191*alphaL[5]*fUpL[15]+(0.2258769757263128*alphaL[7]+0.3535533905932737*alphaL[0])*fUpL[13]+(0.3162277660168379*fUpL[9]+0.2258769757263128*fUpL[7])*alphaL[13]+0.3535533905932737*(fUpL[0]*alphaL[13]+fUpL[6]*alphaL[11])+0.3162277660168379*alphaL[4]*fUpL[10]+0.3535533905932737*(alphaL[3]*fUpL[7]+fUpL[3]*alphaL[7])+0.3162277660168379*(alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]); 
  GhatL[14] = (0.2258769757263128*alphaL[12]+0.3535533905932737*alphaL[1])*fUpL[18]+0.3162277660168379*alphaL[11]*fUpL[17]+0.2258769757263128*alphaL[8]*fUpL[14]+0.3535533905932737*(alphaL[0]*fUpL[14]+alphaL[5]*fUpL[12]+fUpL[5]*alphaL[12])+0.3162277660168379*alphaL[4]*fUpL[10]+0.3535533905932737*(alphaL[3]*fUpL[8]+fUpL[3]*alphaL[8])+0.3162277660168379*alphaL[2]*fUpL[6]; 
  GhatL[15] = 0.3162277660168379*alphaL[11]*fUpL[19]+0.3535533905932737*(alphaL[2]*fUpL[19]+alphaL[4]*fUpL[16])+(0.3162277660168379*alphaL[7]+0.3535533905932737*alphaL[0])*fUpL[15]+0.2828427124746191*(alphaL[5]*fUpL[13]+fUpL[5]*alphaL[13])+0.3535533905932737*alphaL[1]*fUpL[9]+0.3162277660168379*(alphaL[3]*fUpL[5]+fUpL[3]*alphaL[5]); 
  GhatL[16] = (0.3162277660168379*alphaL[12]+0.3535533905932737*alphaL[1])*fUpL[19]+0.3162277660168379*(alphaL[13]*fUpL[17]+alphaL[8]*fUpL[16])+0.3535533905932737*(alphaL[0]*fUpL[16]+alphaL[4]*fUpL[15])+0.3162277660168379*alphaL[5]*fUpL[10]+0.3535533905932737*alphaL[2]*fUpL[9]+0.3162277660168379*alphaL[3]*fUpL[6]; 
  GhatL[17] = 0.2828427124746191*(alphaL[5]*fUpL[19]+alphaL[4]*fUpL[18])+(0.3162277660168379*alphaL[8]+0.2258769757263128*alphaL[7]+0.3535533905932737*alphaL[0])*fUpL[17]+0.3162277660168379*(alphaL[13]*fUpL[16]+alphaL[11]*fUpL[14])+(0.2258769757263128*alphaL[11]+0.3535533905932737*alphaL[2])*fUpL[13]+(0.2258769757263128*fUpL[11]+0.3535533905932737*fUpL[2])*alphaL[13]+0.282842712474619*fUpL[10]*alphaL[12]+0.3535533905932737*(alphaL[3]*fUpL[11]+fUpL[3]*alphaL[11])+0.3162277660168379*alphaL[1]*fUpL[10]+0.3535533905932737*fUpL[6]*alphaL[7]+0.3162277660168379*(alphaL[4]*fUpL[5]+fUpL[4]*alphaL[5]); 
  GhatL[18] = (0.2258769757263128*alphaL[8]+0.3162277660168379*alphaL[7]+0.3535533905932737*alphaL[0])*fUpL[18]+0.2828427124746191*alphaL[4]*fUpL[17]+(0.2258769757263128*alphaL[12]+0.3535533905932737*alphaL[1])*fUpL[14]+0.3162277660168379*(alphaL[12]*fUpL[13]+fUpL[12]*alphaL[13])+0.3535533905932737*(alphaL[3]*fUpL[12]+fUpL[3]*alphaL[12])+fUpL[10]*(0.282842712474619*alphaL[11]+0.3162277660168379*alphaL[2])+0.3535533905932737*(alphaL[5]*fUpL[8]+fUpL[5]*alphaL[8])+0.3162277660168379*alphaL[4]*fUpL[6]; 
  GhatL[19] = (0.3162277660168379*(alphaL[8]+alphaL[7])+0.3535533905932737*alphaL[0])*fUpL[19]+0.2828427124746191*alphaL[5]*fUpL[17]+(0.3162277660168379*alphaL[12]+0.3535533905932737*alphaL[1])*fUpL[16]+(0.3162277660168379*alphaL[11]+0.3535533905932737*alphaL[2])*fUpL[15]+fUpL[10]*(0.282842712474619*alphaL[13]+0.3162277660168379*alphaL[3])+0.3535533905932737*alphaL[4]*fUpL[9]+0.3162277660168379*alphaL[5]*fUpL[6]; 

  double fUpOrdR[27] = {0.};
  double alphaR_n = 0.;

  alphaR_n = (-0.4242640687119286*alphaR[13])-0.4242640687119282*alphaR[12]-0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[5]+0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = ser_4x_p2_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpOrdR[0] = ser_4x_p2_surfx2_eval_quad_node_0_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4242640687119282*alphaR[12])-0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = ser_4x_p2_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpOrdR[1] = ser_4x_p2_surfx2_eval_quad_node_1_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4242640687119286*alphaR[13]-0.4242640687119282*alphaR[12]-0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[5]+0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = ser_4x_p2_surfx2_eval_quad_node_2_r(fc); 
  } else { 
    fUpOrdR[2] = ser_4x_p2_surfx2_eval_quad_node_2_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4242640687119286*alphaR[13])+0.5303300858899102*alphaR[12]-0.3952847075210471*alphaR[8]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[5]-0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = ser_4x_p2_surfx2_eval_quad_node_3_r(fc); 
  } else { 
    fUpOrdR[3] = ser_4x_p2_surfx2_eval_quad_node_3_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5303300858899102*alphaR[12]-0.3952847075210471*alphaR[8]+0.3162277660168378*alphaR[7]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = ser_4x_p2_surfx2_eval_quad_node_4_r(fc); 
  } else { 
    fUpOrdR[4] = ser_4x_p2_surfx2_eval_quad_node_4_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4242640687119286*alphaR[13]+0.5303300858899102*alphaR[12]-0.3952847075210471*alphaR[8]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[5]+0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = ser_4x_p2_surfx2_eval_quad_node_5_r(fc); 
  } else { 
    fUpOrdR[5] = ser_4x_p2_surfx2_eval_quad_node_5_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4242640687119286*alphaR[13])-0.4242640687119282*alphaR[12]+0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[5]-0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = ser_4x_p2_surfx2_eval_quad_node_6_r(fc); 
  } else { 
    fUpOrdR[6] = ser_4x_p2_surfx2_eval_quad_node_6_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4242640687119282*alphaR[12])+0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = ser_4x_p2_surfx2_eval_quad_node_7_r(fc); 
  } else { 
    fUpOrdR[7] = ser_4x_p2_surfx2_eval_quad_node_7_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4242640687119286*alphaR[13]-0.4242640687119282*alphaR[12]+0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[5]-0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = ser_4x_p2_surfx2_eval_quad_node_8_r(fc); 
  } else { 
    fUpOrdR[8] = ser_4x_p2_surfx2_eval_quad_node_8_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5303300858899102*alphaR[13]+0.5303300858899102*alphaR[11]+0.3162277660168378*alphaR[8]-0.3952847075210471*alphaR[7]-0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[9] = ser_4x_p2_surfx2_eval_quad_node_9_r(fc); 
  } else { 
    fUpOrdR[9] = ser_4x_p2_surfx2_eval_quad_node_9_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5303300858899102*alphaR[11]+0.3162277660168378*alphaR[8]-0.3952847075210471*alphaR[7]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[10] = ser_4x_p2_surfx2_eval_quad_node_10_r(fc); 
  } else { 
    fUpOrdR[10] = ser_4x_p2_surfx2_eval_quad_node_10_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5303300858899102*alphaR[13])+0.5303300858899102*alphaR[11]+0.3162277660168378*alphaR[8]-0.3952847075210471*alphaR[7]+0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[11] = ser_4x_p2_surfx2_eval_quad_node_11_r(fc); 
  } else { 
    fUpOrdR[11] = ser_4x_p2_surfx2_eval_quad_node_11_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5303300858899102*alphaR[13]-0.3952847075210471*alphaR[8]-0.3952847075210471*alphaR[7]-0.4743416490252568*alphaR[3]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[12] = ser_4x_p2_surfx2_eval_quad_node_12_r(fc); 
  } else { 
    fUpOrdR[12] = ser_4x_p2_surfx2_eval_quad_node_12_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.3952847075210471*alphaR[8])-0.3952847075210471*alphaR[7]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[13] = ser_4x_p2_surfx2_eval_quad_node_13_r(fc); 
  } else { 
    fUpOrdR[13] = ser_4x_p2_surfx2_eval_quad_node_13_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5303300858899102*alphaR[13])-0.3952847075210471*alphaR[8]-0.3952847075210471*alphaR[7]+0.4743416490252568*alphaR[3]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[14] = ser_4x_p2_surfx2_eval_quad_node_14_r(fc); 
  } else { 
    fUpOrdR[14] = ser_4x_p2_surfx2_eval_quad_node_14_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5303300858899102*alphaR[13]-0.5303300858899102*alphaR[11]+0.3162277660168378*alphaR[8]-0.3952847075210471*alphaR[7]-0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[15] = ser_4x_p2_surfx2_eval_quad_node_15_r(fc); 
  } else { 
    fUpOrdR[15] = ser_4x_p2_surfx2_eval_quad_node_15_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5303300858899102*alphaR[11])+0.3162277660168378*alphaR[8]-0.3952847075210471*alphaR[7]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[16] = ser_4x_p2_surfx2_eval_quad_node_16_r(fc); 
  } else { 
    fUpOrdR[16] = ser_4x_p2_surfx2_eval_quad_node_16_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5303300858899102*alphaR[13])-0.5303300858899102*alphaR[11]+0.3162277660168378*alphaR[8]-0.3952847075210471*alphaR[7]+0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[17] = ser_4x_p2_surfx2_eval_quad_node_17_r(fc); 
  } else { 
    fUpOrdR[17] = ser_4x_p2_surfx2_eval_quad_node_17_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4242640687119286*alphaR[13])+0.4242640687119282*alphaR[12]-0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[5]-0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[18] = ser_4x_p2_surfx2_eval_quad_node_18_r(fc); 
  } else { 
    fUpOrdR[18] = ser_4x_p2_surfx2_eval_quad_node_18_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4242640687119282*alphaR[12]-0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[19] = ser_4x_p2_surfx2_eval_quad_node_19_r(fc); 
  } else { 
    fUpOrdR[19] = ser_4x_p2_surfx2_eval_quad_node_19_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4242640687119286*alphaR[13]+0.4242640687119282*alphaR[12]-0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[5]-0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[3]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[20] = ser_4x_p2_surfx2_eval_quad_node_20_r(fc); 
  } else { 
    fUpOrdR[20] = ser_4x_p2_surfx2_eval_quad_node_20_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4242640687119286*alphaR[13])-0.5303300858899102*alphaR[12]-0.3952847075210471*alphaR[8]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[5]-0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[21] = ser_4x_p2_surfx2_eval_quad_node_21_r(fc); 
  } else { 
    fUpOrdR[21] = ser_4x_p2_surfx2_eval_quad_node_21_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5303300858899102*alphaR[12])-0.3952847075210471*alphaR[8]+0.3162277660168378*alphaR[7]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[22] = ser_4x_p2_surfx2_eval_quad_node_22_r(fc); 
  } else { 
    fUpOrdR[22] = ser_4x_p2_surfx2_eval_quad_node_22_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4242640687119286*alphaR[13]-0.5303300858899102*alphaR[12]-0.3952847075210471*alphaR[8]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[5]+0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[23] = ser_4x_p2_surfx2_eval_quad_node_23_r(fc); 
  } else { 
    fUpOrdR[23] = ser_4x_p2_surfx2_eval_quad_node_23_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.4242640687119286*alphaR[13])+0.4242640687119282*alphaR[12]+0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[5]+0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[24] = ser_4x_p2_surfx2_eval_quad_node_24_r(fc); 
  } else { 
    fUpOrdR[24] = ser_4x_p2_surfx2_eval_quad_node_24_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4242640687119282*alphaR[12]+0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[25] = ser_4x_p2_surfx2_eval_quad_node_25_r(fc); 
  } else { 
    fUpOrdR[25] = ser_4x_p2_surfx2_eval_quad_node_25_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4242640687119286*alphaR[13]+0.4242640687119282*alphaR[12]+0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[8]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[5]+0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[3]+0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[26] = ser_4x_p2_surfx2_eval_quad_node_26_r(fc); 
  } else { 
    fUpOrdR[26] = ser_4x_p2_surfx2_eval_quad_node_26_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[20] = {0.};
  ser_4x_p2_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[20] = {0.}; 
  GhatR[0] = 0.3535533905932737*(alphaR[13]*fUpR[13]+alphaR[12]*fUpR[12]+alphaR[11]*fUpR[11]+alphaR[8]*fUpR[8]+alphaR[7]*fUpR[7]+alphaR[5]*fUpR[5]+alphaR[4]*fUpR[4]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.3162277660168379*(alphaR[5]*fUpR[13]+fUpR[5]*alphaR[13])+0.3535533905932737*(alphaR[8]*fUpR[12]+fUpR[8]*alphaR[12])+0.3162277660168379*(alphaR[4]*fUpR[11]+fUpR[4]*alphaR[11])+0.3162277660168379*(alphaR[1]*fUpR[7]+fUpR[1]*alphaR[7])+0.3535533905932737*(alphaR[3]*fUpR[5]+fUpR[3]*alphaR[5]+alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.3535533905932737*alphaR[13]*fUpR[17]+0.3162277660168379*(alphaR[4]*fUpR[12]+fUpR[4]*alphaR[12])+0.3535533905932737*(alphaR[7]*fUpR[11]+fUpR[7]*alphaR[11]+alphaR[5]*fUpR[10])+0.3162277660168379*(alphaR[2]*fUpR[8]+fUpR[2]*alphaR[8])+0.3535533905932737*(alphaR[3]*fUpR[6]+alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.3535533905932737*(alphaR[12]*fUpR[18]+alphaR[11]*fUpR[17])+0.3162277660168379*alphaR[5]*fUpR[15]+0.3535533905932737*(alphaR[8]*fUpR[14]+alphaR[7]*fUpR[13]+fUpR[7]*alphaR[13]+alphaR[4]*fUpR[10])+0.3162277660168379*alphaR[3]*fUpR[9]+0.3535533905932737*(alphaR[2]*fUpR[6]+alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]+alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]); 
  GhatR[4] = 0.3162277660168379*alphaR[5]*fUpR[17]+0.3162277660168379*fUpR[10]*alphaR[13]+(0.2828427124746191*alphaR[11]+0.3162277660168379*alphaR[2])*fUpR[12]+0.2828427124746191*fUpR[11]*alphaR[12]+0.3162277660168379*(fUpR[2]*alphaR[12]+alphaR[1]*fUpR[11]+fUpR[1]*alphaR[11])+0.3535533905932737*alphaR[3]*fUpR[10]+0.3162277660168379*(alphaR[4]*fUpR[8]+fUpR[4]*alphaR[8]+alphaR[4]*fUpR[7]+fUpR[4]*alphaR[7])+0.3535533905932737*(alphaR[5]*fUpR[6]+alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[5] = 0.3535533905932737*alphaR[8]*fUpR[18]+0.3162277660168379*alphaR[4]*fUpR[17]+(0.2828427124746191*alphaR[13]+0.3162277660168379*alphaR[3])*fUpR[15]+0.3535533905932737*alphaR[12]*fUpR[14]+0.3162277660168379*(alphaR[1]*fUpR[13]+fUpR[1]*alphaR[13])+fUpR[10]*(0.3162277660168379*alphaR[11]+0.3535533905932737*alphaR[2])+0.3162277660168379*(alphaR[5]*(fUpR[9]+fUpR[7])+fUpR[5]*alphaR[7])+0.3535533905932737*(alphaR[4]*fUpR[6]+alphaR[0]*fUpR[5]+fUpR[0]*alphaR[5]+alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[6] = 0.3162277660168379*(alphaR[5]*fUpR[19]+alphaR[4]*fUpR[18])+0.3535533905932737*alphaR[7]*fUpR[17]+0.3162277660168379*(alphaR[3]*fUpR[16]+alphaR[2]*fUpR[14])+0.3535533905932737*(alphaR[11]*fUpR[13]+fUpR[11]*alphaR[13])+fUpR[10]*(0.3162277660168379*alphaR[12]+0.3535533905932737*alphaR[1])+0.3162277660168379*fUpR[6]*alphaR[8]+0.3535533905932737*(alphaR[0]*fUpR[6]+alphaR[4]*fUpR[5]+fUpR[4]*alphaR[5]+alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 
  GhatR[7] = 0.2258769757263128*alphaR[13]*fUpR[13]+0.3535533905932737*(alphaR[3]*fUpR[13]+fUpR[3]*alphaR[13])+0.3162277660168379*alphaR[12]*fUpR[12]+0.2258769757263128*alphaR[11]*fUpR[11]+0.3535533905932737*(alphaR[2]*fUpR[11]+fUpR[2]*alphaR[11])+0.2258769757263128*alphaR[7]*fUpR[7]+0.3535533905932737*(alphaR[0]*fUpR[7]+fUpR[0]*alphaR[7])+0.3162277660168379*(alphaR[5]*fUpR[5]+alphaR[4]*fUpR[4]+alphaR[1]*fUpR[1]); 
  GhatR[8] = 0.3535533905932737*(alphaR[5]*fUpR[18]+alphaR[3]*fUpR[14])+0.2258769757263128*alphaR[12]*fUpR[12]+0.3535533905932737*(alphaR[1]*fUpR[12]+fUpR[1]*alphaR[12])+0.3162277660168379*alphaR[11]*fUpR[11]+0.2258769757263128*alphaR[8]*fUpR[8]+0.3535533905932737*(alphaR[0]*fUpR[8]+fUpR[0]*alphaR[8])+0.3162277660168379*(alphaR[4]*fUpR[4]+alphaR[2]*fUpR[2]); 
  GhatR[9] = 0.3535533905932737*(alphaR[4]*fUpR[19]+alphaR[2]*fUpR[16]+alphaR[1]*fUpR[15])+0.3162277660168379*alphaR[13]*fUpR[13]+0.3535533905932737*alphaR[0]*fUpR[9]+0.3162277660168379*(alphaR[5]*fUpR[5]+alphaR[3]*fUpR[3]); 
  GhatR[10] = (0.282842712474619*alphaR[13]+0.3162277660168379*alphaR[3])*fUpR[19]+(0.282842712474619*alphaR[11]+0.3162277660168379*alphaR[2])*fUpR[18]+(0.282842712474619*alphaR[12]+0.3162277660168379*alphaR[1])*fUpR[17]+0.3162277660168379*(alphaR[5]*fUpR[16]+alphaR[4]*(fUpR[14]+fUpR[13])+fUpR[4]*alphaR[13]+fUpR[6]*alphaR[12]+alphaR[5]*fUpR[11]+fUpR[5]*alphaR[11])+0.3162277660168379*(alphaR[8]+alphaR[7])*fUpR[10]+0.3535533905932737*(alphaR[0]*fUpR[10]+alphaR[1]*fUpR[6]+alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5]+alphaR[3]*fUpR[4]+fUpR[3]*alphaR[4]); 
  GhatR[11] = 0.2258769757263128*alphaR[13]*fUpR[17]+0.3535533905932737*(alphaR[3]*fUpR[17]+fUpR[6]*alphaR[13])+0.2828427124746191*(alphaR[4]*fUpR[12]+fUpR[4]*alphaR[12])+(0.3162277660168379*alphaR[8]+0.2258769757263128*alphaR[7]+0.3535533905932737*alphaR[0])*fUpR[11]+(0.3162277660168379*fUpR[8]+0.2258769757263128*fUpR[7]+0.3535533905932737*fUpR[0])*alphaR[11]+0.3162277660168379*alphaR[5]*fUpR[10]+0.3535533905932737*(alphaR[2]*fUpR[7]+fUpR[2]*alphaR[7])+0.3162277660168379*(alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]); 
  GhatR[12] = 0.3162277660168379*alphaR[13]*fUpR[18]+0.3535533905932737*(alphaR[3]*fUpR[18]+alphaR[5]*fUpR[14])+(0.2258769757263128*alphaR[8]+0.3162277660168379*alphaR[7]+0.3535533905932737*alphaR[0])*fUpR[12]+(0.2258769757263128*fUpR[8]+0.3162277660168379*fUpR[7]+0.3535533905932737*fUpR[0])*alphaR[12]+0.2828427124746191*(alphaR[4]*fUpR[11]+fUpR[4]*alphaR[11])+0.3535533905932737*(alphaR[1]*fUpR[8]+fUpR[1]*alphaR[8])+0.3162277660168379*(alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]); 
  GhatR[13] = 0.3162277660168379*alphaR[12]*fUpR[18]+(0.2258769757263128*alphaR[11]+0.3535533905932737*alphaR[2])*fUpR[17]+0.2828427124746191*alphaR[5]*fUpR[15]+(0.2258769757263128*alphaR[7]+0.3535533905932737*alphaR[0])*fUpR[13]+(0.3162277660168379*fUpR[9]+0.2258769757263128*fUpR[7])*alphaR[13]+0.3535533905932737*(fUpR[0]*alphaR[13]+fUpR[6]*alphaR[11])+0.3162277660168379*alphaR[4]*fUpR[10]+0.3535533905932737*(alphaR[3]*fUpR[7]+fUpR[3]*alphaR[7])+0.3162277660168379*(alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]); 
  GhatR[14] = (0.2258769757263128*alphaR[12]+0.3535533905932737*alphaR[1])*fUpR[18]+0.3162277660168379*alphaR[11]*fUpR[17]+0.2258769757263128*alphaR[8]*fUpR[14]+0.3535533905932737*(alphaR[0]*fUpR[14]+alphaR[5]*fUpR[12]+fUpR[5]*alphaR[12])+0.3162277660168379*alphaR[4]*fUpR[10]+0.3535533905932737*(alphaR[3]*fUpR[8]+fUpR[3]*alphaR[8])+0.3162277660168379*alphaR[2]*fUpR[6]; 
  GhatR[15] = 0.3162277660168379*alphaR[11]*fUpR[19]+0.3535533905932737*(alphaR[2]*fUpR[19]+alphaR[4]*fUpR[16])+(0.3162277660168379*alphaR[7]+0.3535533905932737*alphaR[0])*fUpR[15]+0.2828427124746191*(alphaR[5]*fUpR[13]+fUpR[5]*alphaR[13])+0.3535533905932737*alphaR[1]*fUpR[9]+0.3162277660168379*(alphaR[3]*fUpR[5]+fUpR[3]*alphaR[5]); 
  GhatR[16] = (0.3162277660168379*alphaR[12]+0.3535533905932737*alphaR[1])*fUpR[19]+0.3162277660168379*(alphaR[13]*fUpR[17]+alphaR[8]*fUpR[16])+0.3535533905932737*(alphaR[0]*fUpR[16]+alphaR[4]*fUpR[15])+0.3162277660168379*alphaR[5]*fUpR[10]+0.3535533905932737*alphaR[2]*fUpR[9]+0.3162277660168379*alphaR[3]*fUpR[6]; 
  GhatR[17] = 0.2828427124746191*(alphaR[5]*fUpR[19]+alphaR[4]*fUpR[18])+(0.3162277660168379*alphaR[8]+0.2258769757263128*alphaR[7]+0.3535533905932737*alphaR[0])*fUpR[17]+0.3162277660168379*(alphaR[13]*fUpR[16]+alphaR[11]*fUpR[14])+(0.2258769757263128*alphaR[11]+0.3535533905932737*alphaR[2])*fUpR[13]+(0.2258769757263128*fUpR[11]+0.3535533905932737*fUpR[2])*alphaR[13]+0.282842712474619*fUpR[10]*alphaR[12]+0.3535533905932737*(alphaR[3]*fUpR[11]+fUpR[3]*alphaR[11])+0.3162277660168379*alphaR[1]*fUpR[10]+0.3535533905932737*fUpR[6]*alphaR[7]+0.3162277660168379*(alphaR[4]*fUpR[5]+fUpR[4]*alphaR[5]); 
  GhatR[18] = (0.2258769757263128*alphaR[8]+0.3162277660168379*alphaR[7]+0.3535533905932737*alphaR[0])*fUpR[18]+0.2828427124746191*alphaR[4]*fUpR[17]+(0.2258769757263128*alphaR[12]+0.3535533905932737*alphaR[1])*fUpR[14]+0.3162277660168379*(alphaR[12]*fUpR[13]+fUpR[12]*alphaR[13])+0.3535533905932737*(alphaR[3]*fUpR[12]+fUpR[3]*alphaR[12])+fUpR[10]*(0.282842712474619*alphaR[11]+0.3162277660168379*alphaR[2])+0.3535533905932737*(alphaR[5]*fUpR[8]+fUpR[5]*alphaR[8])+0.3162277660168379*alphaR[4]*fUpR[6]; 
  GhatR[19] = (0.3162277660168379*(alphaR[8]+alphaR[7])+0.3535533905932737*alphaR[0])*fUpR[19]+0.2828427124746191*alphaR[5]*fUpR[17]+(0.3162277660168379*alphaR[12]+0.3535533905932737*alphaR[1])*fUpR[16]+(0.3162277660168379*alphaR[11]+0.3535533905932737*alphaR[2])*fUpR[15]+fUpR[10]*(0.282842712474619*alphaR[13]+0.3162277660168379*alphaR[3])+0.3535533905932737*alphaR[4]*fUpR[9]+0.3162277660168379*alphaR[5]*fUpR[6]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdy2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdy2; 
  out[2] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdy2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdy2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdy2; 
  out[5] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdy2; 
  out[6] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdy2; 
  out[7] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdy2; 
  out[8] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdy2; 
  out[9] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdy2; 
  out[10] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdy2; 
  out[11] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdy2; 
  out[12] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdy2; 
  out[13] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdy2; 
  out[14] += (0.7071067811865475*GhatL[9]-0.7071067811865475*GhatR[9])*rdy2; 
  out[15] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdy2; 
  out[16] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdy2; 
  out[17] += (0.7071067811865475*GhatL[10]-0.7071067811865475*GhatR[10])*rdy2; 
  out[18] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdy2; 
  out[19] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdy2; 
  out[20] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdy2; 
  out[21] += (0.7071067811865475*GhatL[11]-0.7071067811865475*GhatR[11])*rdy2; 
  out[22] += (1.58113883008419*GhatL[2]-1.58113883008419*GhatR[2])*rdy2; 
  out[23] += (0.7071067811865475*GhatL[12]-0.7071067811865475*GhatR[12])*rdy2; 
  out[24] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdy2; 
  out[25] += (0.7071067811865475*GhatL[13]-0.7071067811865475*GhatR[13])*rdy2; 
  out[26] += (1.58113883008419*GhatL[3]-1.58113883008419*GhatR[3])*rdy2; 
  out[27] += (0.7071067811865475*GhatL[14]-0.7071067811865475*GhatR[14])*rdy2; 
  out[28] += (0.7071067811865475*GhatL[15]-0.7071067811865475*GhatR[15])*rdy2; 
  out[29] += ((-1.224744871391589*GhatR[9])-1.224744871391589*GhatL[9])*rdy2; 
  out[30] += (0.7071067811865475*GhatL[16]-0.7071067811865475*GhatR[16])*rdy2; 
  out[31] += ((-1.224744871391589*GhatR[10])-1.224744871391589*GhatL[10])*rdy2; 
  out[32] += ((-1.224744871391589*GhatR[11])-1.224744871391589*GhatL[11])*rdy2; 
  out[33] += (1.58113883008419*GhatL[4]-1.58113883008419*GhatR[4])*rdy2; 
  out[34] += ((-1.224744871391589*GhatR[12])-1.224744871391589*GhatL[12])*rdy2; 
  out[35] += ((-1.224744871391589*GhatR[13])-1.224744871391589*GhatL[13])*rdy2; 
  out[36] += (1.58113883008419*GhatL[5]-1.58113883008419*GhatR[5])*rdy2; 
  out[37] += (0.7071067811865475*GhatL[17]-0.7071067811865475*GhatR[17])*rdy2; 
  out[38] += (1.58113883008419*GhatL[6]-1.58113883008419*GhatR[6])*rdy2; 
  out[39] += (0.7071067811865475*GhatL[18]-0.7071067811865475*GhatR[18])*rdy2; 
  out[40] += ((-1.224744871391589*GhatR[14])-1.224744871391589*GhatL[14])*rdy2; 
  out[41] += ((-1.224744871391589*GhatR[15])-1.224744871391589*GhatL[15])*rdy2; 
  out[42] += (0.7071067811865475*GhatL[19]-0.7071067811865475*GhatR[19])*rdy2; 
  out[43] += ((-1.224744871391589*GhatR[16])-1.224744871391589*GhatL[16])*rdy2; 
  out[44] += ((-1.224744871391589*GhatR[17])-1.224744871391589*GhatL[17])*rdy2; 
  out[45] += (1.58113883008419*GhatL[10]-1.58113883008419*GhatR[10])*rdy2; 
  out[46] += ((-1.224744871391589*GhatR[18])-1.224744871391589*GhatL[18])*rdy2; 
  out[47] += ((-1.224744871391589*GhatR[19])-1.224744871391589*GhatL[19])*rdy2; 

  return 2.5*rdy2*cflFreq; 

} 
