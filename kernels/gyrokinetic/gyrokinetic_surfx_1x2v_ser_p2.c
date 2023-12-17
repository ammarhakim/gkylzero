#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_ser_3x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double *alpha_surf_l, const double *alpha_surf_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // alpha_surf_l: Surface expansion of phase space flux on the left.
  // alpha_surf_r: Surface expansion of phase space flux on the right.
  // fl,fc,fr: distribution function in left, center and right cells.
  // out: output increment in center cell.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];
  double wmu = w[2];
  double rdmu2 = 2.0/dxv[2];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[2]*w[2];
  double rdmu2Sq = rdmu2*rdmu2;

  const double *alphaL = &alpha_surf_l[0];
  const double *alphaR = &alpha_surf_r[0];
  double cflFreq = 0.0;
  double fUpOrdL[9] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.5*alphaL[0]-0.6708203932499357*alphaL[1];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = ser_3x_p2_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpOrdL[0] = ser_3x_p2_surfx1_eval_quad_node_0_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0]-0.6708203932499357*alphaL[1];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = ser_3x_p2_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpOrdL[1] = ser_3x_p2_surfx1_eval_quad_node_1_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0]-0.6708203932499357*alphaL[1];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = ser_3x_p2_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpOrdL[2] = ser_3x_p2_surfx1_eval_quad_node_2_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = ser_3x_p2_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpOrdL[3] = ser_3x_p2_surfx1_eval_quad_node_3_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = ser_3x_p2_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpOrdL[4] = ser_3x_p2_surfx1_eval_quad_node_4_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = ser_3x_p2_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpOrdL[5] = ser_3x_p2_surfx1_eval_quad_node_5_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = ser_3x_p2_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpOrdL[6] = ser_3x_p2_surfx1_eval_quad_node_6_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = ser_3x_p2_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpOrdL[7] = ser_3x_p2_surfx1_eval_quad_node_7_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = ser_3x_p2_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpOrdL[8] = ser_3x_p2_surfx1_eval_quad_node_8_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[8] = {0.};
  ser_3x_p2_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[8] = {0.}; 
  GhatL[0] = 0.5*(alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.4472135954999579*alphaL[1]*fUpL[4]+0.5*(alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.5*(alphaL[1]*fUpL[3]+alphaL[0]*fUpL[2]); 
  GhatL[3] = 0.447213595499958*alphaL[1]*fUpL[6]+0.5*(alphaL[0]*fUpL[3]+alphaL[1]*fUpL[2]); 
  GhatL[4] = 0.5*alphaL[0]*fUpL[4]+0.4472135954999579*alphaL[1]*fUpL[1]; 
  GhatL[5] = 0.5000000000000001*alphaL[1]*fUpL[7]+0.5*alphaL[0]*fUpL[5]; 
  GhatL[6] = 0.5*alphaL[0]*fUpL[6]+0.447213595499958*alphaL[1]*fUpL[3]; 
  GhatL[7] = 0.5*alphaL[0]*fUpL[7]+0.5000000000000001*alphaL[1]*fUpL[5]; 

  double fUpOrdR[9] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.5*alphaR[0]-0.6708203932499357*alphaR[1];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = ser_3x_p2_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpOrdR[0] = ser_3x_p2_surfx1_eval_quad_node_0_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0]-0.6708203932499357*alphaR[1];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = ser_3x_p2_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpOrdR[1] = ser_3x_p2_surfx1_eval_quad_node_1_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0]-0.6708203932499357*alphaR[1];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = ser_3x_p2_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpOrdR[2] = ser_3x_p2_surfx1_eval_quad_node_2_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = ser_3x_p2_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpOrdR[3] = ser_3x_p2_surfx1_eval_quad_node_3_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = ser_3x_p2_surfx1_eval_quad_node_4_r(fc); 
  } else { 
    fUpOrdR[4] = ser_3x_p2_surfx1_eval_quad_node_4_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = ser_3x_p2_surfx1_eval_quad_node_5_r(fc); 
  } else { 
    fUpOrdR[5] = ser_3x_p2_surfx1_eval_quad_node_5_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = ser_3x_p2_surfx1_eval_quad_node_6_r(fc); 
  } else { 
    fUpOrdR[6] = ser_3x_p2_surfx1_eval_quad_node_6_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = ser_3x_p2_surfx1_eval_quad_node_7_r(fc); 
  } else { 
    fUpOrdR[7] = ser_3x_p2_surfx1_eval_quad_node_7_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = ser_3x_p2_surfx1_eval_quad_node_8_r(fc); 
  } else { 
    fUpOrdR[8] = ser_3x_p2_surfx1_eval_quad_node_8_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[8] = {0.};
  ser_3x_p2_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[8] = {0.}; 
  GhatR[0] = 0.5*(alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.4472135954999579*alphaR[1]*fUpR[4]+0.5*(alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.5*(alphaR[1]*fUpR[3]+alphaR[0]*fUpR[2]); 
  GhatR[3] = 0.447213595499958*alphaR[1]*fUpR[6]+0.5*(alphaR[0]*fUpR[3]+alphaR[1]*fUpR[2]); 
  GhatR[4] = 0.5*alphaR[0]*fUpR[4]+0.4472135954999579*alphaR[1]*fUpR[1]; 
  GhatR[5] = 0.5000000000000001*alphaR[1]*fUpR[7]+0.5*alphaR[0]*fUpR[5]; 
  GhatR[6] = 0.5*alphaR[0]*fUpR[6]+0.447213595499958*alphaR[1]*fUpR[3]; 
  GhatR[7] = 0.5*alphaR[0]*fUpR[7]+0.5000000000000001*alphaR[1]*fUpR[5]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdx2; 
  out[4] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 
  out[5] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdx2; 
  out[6] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdx2; 
  out[7] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdx2; 
  out[8] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdx2; 
  out[9] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdx2; 
  out[10] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdx2; 
  out[11] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdx2; 
  out[12] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdx2; 
  out[13] += (1.58113883008419*GhatL[2]-1.58113883008419*GhatR[2])*rdx2; 
  out[14] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdx2; 
  out[15] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdx2; 
  out[16] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdx2; 
  out[17] += (1.58113883008419*GhatL[3]-1.58113883008419*GhatR[3])*rdx2; 
  out[18] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdx2; 
  out[19] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdx2; 

  return 2.5*rdx2*cflFreq; 

} 
