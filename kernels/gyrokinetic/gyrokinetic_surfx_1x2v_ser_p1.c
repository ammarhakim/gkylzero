#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x2v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *alpha_surf_l, const double *alpha_surf_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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
  double fUpOrdL[6] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.5*alphaL[0]-0.6708203932499357*alphaL[1];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_1x2v_p1_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpOrdL[0] = gkhyb_1x2v_p1_surfx1_eval_quad_node_0_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_1x2v_p1_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpOrdL[1] = gkhyb_1x2v_p1_surfx1_eval_quad_node_1_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = gkhyb_1x2v_p1_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpOrdL[2] = gkhyb_1x2v_p1_surfx1_eval_quad_node_2_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0]-0.6708203932499357*alphaL[1];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = gkhyb_1x2v_p1_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpOrdL[3] = gkhyb_1x2v_p1_surfx1_eval_quad_node_3_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = gkhyb_1x2v_p1_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpOrdL[4] = gkhyb_1x2v_p1_surfx1_eval_quad_node_4_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = gkhyb_1x2v_p1_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpOrdL[5] = gkhyb_1x2v_p1_surfx1_eval_quad_node_5_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[6] = {0.};
  gkhyb_1x2v_p1_xdir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[6] = {0.}; 
  GhatL[0] = 0.5*(alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.4472135954999579*alphaL[1]*fUpL[4]+0.5*(alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.5*(alphaL[1]*fUpL[3]+alphaL[0]*fUpL[2]); 
  GhatL[3] = 0.447213595499958*alphaL[1]*fUpL[5]+0.5*(alphaL[0]*fUpL[3]+alphaL[1]*fUpL[2]); 
  GhatL[4] = 0.5*alphaL[0]*fUpL[4]+0.4472135954999579*alphaL[1]*fUpL[1]; 
  GhatL[5] = 0.5*alphaL[0]*fUpL[5]+0.447213595499958*alphaL[1]*fUpL[3]; 

  double fUpOrdR[6] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.5*alphaR[0]-0.6708203932499357*alphaR[1];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_1x2v_p1_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpOrdR[0] = gkhyb_1x2v_p1_surfx1_eval_quad_node_0_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_1x2v_p1_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpOrdR[1] = gkhyb_1x2v_p1_surfx1_eval_quad_node_1_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = gkhyb_1x2v_p1_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpOrdR[2] = gkhyb_1x2v_p1_surfx1_eval_quad_node_2_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0]-0.6708203932499357*alphaR[1];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = gkhyb_1x2v_p1_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpOrdR[3] = gkhyb_1x2v_p1_surfx1_eval_quad_node_3_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = gkhyb_1x2v_p1_surfx1_eval_quad_node_4_r(fc); 
  } else { 
    fUpOrdR[4] = gkhyb_1x2v_p1_surfx1_eval_quad_node_4_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = gkhyb_1x2v_p1_surfx1_eval_quad_node_5_r(fc); 
  } else { 
    fUpOrdR[5] = gkhyb_1x2v_p1_surfx1_eval_quad_node_5_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[6] = {0.};
  gkhyb_1x2v_p1_xdir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[6] = {0.}; 
  GhatR[0] = 0.5*(alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.4472135954999579*alphaR[1]*fUpR[4]+0.5*(alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.5*(alphaR[1]*fUpR[3]+alphaR[0]*fUpR[2]); 
  GhatR[3] = 0.447213595499958*alphaR[1]*fUpR[5]+0.5*(alphaR[0]*fUpR[3]+alphaR[1]*fUpR[2]); 
  GhatR[4] = 0.5*alphaR[0]*fUpR[4]+0.4472135954999579*alphaR[1]*fUpR[1]; 
  GhatR[5] = 0.5*alphaR[0]*fUpR[5]+0.447213595499958*alphaR[1]*fUpR[3]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdx2; 
  out[4] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 
  out[5] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdx2; 
  out[6] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdx2; 
  out[7] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdx2; 
  out[8] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdx2; 
  out[9] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdx2; 
  out[10] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdx2; 
  out[11] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdx2; 

  return 1.5*rdx2*cflFreq; 

} 
