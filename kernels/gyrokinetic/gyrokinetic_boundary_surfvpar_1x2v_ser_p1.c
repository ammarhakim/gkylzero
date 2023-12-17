#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x2v_p1_surfx2_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_ser_p1(const double *w, const double *dxv, const double *alpha_edge, const double *alpha_skin, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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

  double cflFreq = 0.0;

  if (edge == -1) { 

  const double *alphaR = &alpha_edge[6];
  double fUpOrdR[4] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.5*alphaR[3]-0.5*alphaR[2]-0.5*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5*alphaR[3])+0.5*alphaR[2]-0.5*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.5*alphaR[3])-0.5*alphaR[2]+0.5*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[3]+0.5*alphaR[2]+0.5*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[4] = {0.};
  gkhyb_1x2v_p1_vpardir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[4] = {0.}; 
  GhatR[0] = 0.5*(alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.5*(alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.5*(alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.5*(alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[4] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[5] += -0.7071067811865475*GhatR[3]*rdvpar2; 
  out[6] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[7] += -1.224744871391589*GhatR[3]*rdvpar2; 
  out[8] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[9] += -1.58113883008419*GhatR[1]*rdvpar2; 
  out[10] += -1.58113883008419*GhatR[2]*rdvpar2; 
  out[11] += -1.58113883008419*GhatR[3]*rdvpar2; 

  } else { 

  const double *alphaL = &alpha_skin[6];
  double fUpOrdL[4] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.5*alphaL[3]-0.5*alphaL[2]-0.5*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = gkhyb_1x2v_p1_surfx2_eval_quad_node_0_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5*alphaL[3])+0.5*alphaL[2]-0.5*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = gkhyb_1x2v_p1_surfx2_eval_quad_node_1_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.5*alphaL[3])-0.5*alphaL[2]+0.5*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = gkhyb_1x2v_p1_surfx2_eval_quad_node_2_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[3]+0.5*alphaL[2]+0.5*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = gkhyb_1x2v_p1_surfx2_eval_quad_node_3_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[4] = {0.};
  gkhyb_1x2v_p1_vpardir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[4] = {0.}; 
  GhatL[0] = 0.5*(alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.5*(alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.5*(alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.5*(alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[4] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[5] += 0.7071067811865475*GhatL[3]*rdvpar2; 
  out[6] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[7] += -1.224744871391589*GhatL[3]*rdvpar2; 
  out[8] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[9] += 1.58113883008419*GhatL[1]*rdvpar2; 
  out[10] += 1.58113883008419*GhatL[2]*rdvpar2; 
  out[11] += 1.58113883008419*GhatL[3]*rdvpar2; 

  } 

  return 2.5*rdvpar2*cflFreq; 

} 
