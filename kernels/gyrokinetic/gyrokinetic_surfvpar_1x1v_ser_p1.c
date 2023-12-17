#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x1v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x1v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *alpha_surf_l, const double *alpha_surf_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;

  const double *alphaL = &alpha_surf_l[3];
  const double *alphaR = &alpha_surf_r[3];
  double cflFreq = 0.0;
  double fUpOrdL[2] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.7071067811865468*alphaL[0]-0.7071067811865468*alphaL[1];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_r(fl); 
  } else { 
    fUpOrdL[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.7071067811865468*alphaL[1]+0.7071067811865468*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_r(fl); 
  } else { 
    fUpOrdL[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[2] = {0.};
  gkhyb_1x1v_p1_vpardir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[2] = {0.}; 
  GhatL[0] = 0.7071067811865475*(alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.7071067811865475*(alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 

  double fUpOrdR[2] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.7071067811865468*alphaR[0]-0.7071067811865468*alphaR[1];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_r(fc); 
  } else { 
    fUpOrdR[0] = gkhyb_1x1v_p1_surfx2_eval_quad_node_0_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.7071067811865468*alphaR[1]+0.7071067811865468*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_r(fc); 
  } else { 
    fUpOrdR[1] = gkhyb_1x1v_p1_surfx2_eval_quad_node_1_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[2] = {0.};
  gkhyb_1x1v_p1_vpardir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[2] = {0.}; 
  GhatR[0] = 0.7071067811865475*(alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.7071067811865475*(alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdvpar2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdvpar2; 
  out[2] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdvpar2; 
  out[3] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdvpar2; 
  out[4] += (1.58113883008419*GhatL[0]-1.58113883008419*GhatR[0])*rdvpar2; 
  out[5] += (1.58113883008419*GhatL[1]-1.58113883008419*GhatR[1])*rdvpar2; 

  return 2.5*rdvpar2*cflFreq; 

} 
