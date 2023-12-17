#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_ser_3x_p2_surfx2_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_1x2v_ser_p2(const double *w, const double *dxv, const double *alpha_edge, const double *alpha_skin, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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

  const double *alphaR = &alpha_edge[8];
  double fUpOrdR[9] = {0.};
  double alphaR_n = 0.;

  alphaR_n = (-0.6*alphaR[6])+0.4472135954999572*alphaR[4]+0.9*alphaR[3]-0.6708203932499357*alphaR[2]-0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]-0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.6*alphaR[6]+0.4472135954999572*alphaR[4]-0.9*alphaR[3]+0.6708203932499357*alphaR[2]-0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.75*alphaR[6]-0.5590169943749465*alphaR[4]-0.6708203932499357*alphaR[2]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0]-0.5590169943749465*alphaR[4];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(fskin); 
  } else { 
    fUpOrdR[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.75*alphaR[6])-0.5590169943749465*alphaR[4]+0.6708203932499357*alphaR[2]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(fskin); 
  } else { 
    fUpOrdR[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.6*alphaR[6])+0.4472135954999572*alphaR[4]-0.9*alphaR[3]-0.6708203932499357*alphaR[2]+0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(fskin); 
  } else { 
    fUpOrdR[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]+0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(fskin); 
  } else { 
    fUpOrdR[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.6*alphaR[6]+0.4472135954999572*alphaR[4]+0.9*alphaR[3]+0.6708203932499357*alphaR[2]+0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(fskin); 
  } else { 
    fUpOrdR[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[8] = {0.};
  ser_3x_p2_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[8] = {0.}; 
  GhatR[0] = 0.5*(alphaR[6]*fUpR[6]+alphaR[4]*fUpR[4]+alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.447213595499958*(alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6])+0.4472135954999579*(alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4])+0.5*(alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.447213595499958*alphaR[3]*fUpR[7]+0.5000000000000001*(alphaR[4]*fUpR[6]+fUpR[4]*alphaR[6])+0.4472135954999579*alphaR[2]*fUpR[5]+0.5*(alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.4*alphaR[6]*fUpR[7]+0.447213595499958*(alphaR[2]*fUpR[7]+alphaR[1]*fUpR[6]+fUpR[1]*alphaR[6])+0.4472135954999579*(alphaR[3]*(fUpR[5]+fUpR[4])+fUpR[3]*alphaR[4])+0.5*(alphaR[0]*fUpR[3]+fUpR[0]*alphaR[3]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[4] = 0.31943828249997*alphaR[6]*fUpR[6]+0.5000000000000001*(alphaR[2]*fUpR[6]+fUpR[2]*alphaR[6])+0.31943828249997*alphaR[4]*fUpR[4]+0.5*(alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4])+0.4472135954999579*(alphaR[3]*fUpR[3]+alphaR[1]*fUpR[1]); 
  GhatR[5] = 0.5000000000000001*alphaR[1]*fUpR[7]+0.4472135954999579*alphaR[6]*fUpR[6]+0.5*alphaR[0]*fUpR[5]+0.4472135954999579*(alphaR[3]*fUpR[3]+alphaR[2]*fUpR[2]); 
  GhatR[6] = 0.4*alphaR[3]*fUpR[7]+(0.31943828249997*alphaR[4]+0.5*alphaR[0])*fUpR[6]+(0.4472135954999579*fUpR[5]+0.31943828249997*fUpR[4]+0.5*fUpR[0])*alphaR[6]+0.5000000000000001*(alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4])+0.447213595499958*(alphaR[1]*fUpR[3]+fUpR[1]*alphaR[3]); 
  GhatR[7] = (0.4472135954999579*alphaR[4]+0.5*alphaR[0])*fUpR[7]+0.4*(alphaR[3]*fUpR[6]+fUpR[3]*alphaR[6])+0.5000000000000001*alphaR[1]*fUpR[5]+0.447213595499958*(alphaR[2]*fUpR[3]+fUpR[2]*alphaR[3]); 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[4] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[5] += -0.7071067811865475*GhatR[3]*rdvpar2; 
  out[6] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[7] += -0.7071067811865475*GhatR[4]*rdvpar2; 
  out[8] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[9] += -0.7071067811865475*GhatR[5]*rdvpar2; 
  out[10] += -1.224744871391589*GhatR[3]*rdvpar2; 
  out[11] += -1.224744871391589*GhatR[4]*rdvpar2; 
  out[12] += -1.58113883008419*GhatR[1]*rdvpar2; 
  out[13] += -0.7071067811865475*GhatR[6]*rdvpar2; 
  out[14] += -1.58113883008419*GhatR[2]*rdvpar2; 
  out[15] += -0.7071067811865475*GhatR[7]*rdvpar2; 
  out[16] += -1.224744871391589*GhatR[5]*rdvpar2; 
  out[17] += -1.224744871391589*GhatR[6]*rdvpar2; 
  out[18] += -1.58113883008419*GhatR[3]*rdvpar2; 
  out[19] += -1.224744871391589*GhatR[7]*rdvpar2; 

  } else { 

  const double *alphaL = &alpha_skin[8];
  double fUpOrdL[9] = {0.};
  double alphaL_n = 0.;

  alphaL_n = (-0.6*alphaL[6])+0.4472135954999572*alphaL[4]+0.9*alphaL[3]-0.6708203932499357*alphaL[2]-0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = ser_3x_p2_surfx2_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = ser_3x_p2_surfx2_eval_quad_node_0_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]-0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = ser_3x_p2_surfx2_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = ser_3x_p2_surfx2_eval_quad_node_1_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.6*alphaL[6]+0.4472135954999572*alphaL[4]-0.9*alphaL[3]+0.6708203932499357*alphaL[2]-0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = ser_3x_p2_surfx2_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = ser_3x_p2_surfx2_eval_quad_node_2_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.75*alphaL[6]-0.5590169943749465*alphaL[4]-0.6708203932499357*alphaL[2]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = ser_3x_p2_surfx2_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = ser_3x_p2_surfx2_eval_quad_node_3_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0]-0.5590169943749465*alphaL[4];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = ser_3x_p2_surfx2_eval_quad_node_4_r(fedge); 
  } else { 
    fUpOrdL[4] = ser_3x_p2_surfx2_eval_quad_node_4_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.75*alphaL[6])-0.5590169943749465*alphaL[4]+0.6708203932499357*alphaL[2]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = ser_3x_p2_surfx2_eval_quad_node_5_r(fedge); 
  } else { 
    fUpOrdL[5] = ser_3x_p2_surfx2_eval_quad_node_5_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.6*alphaL[6])+0.4472135954999572*alphaL[4]-0.9*alphaL[3]-0.6708203932499357*alphaL[2]+0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = ser_3x_p2_surfx2_eval_quad_node_6_r(fedge); 
  } else { 
    fUpOrdL[6] = ser_3x_p2_surfx2_eval_quad_node_6_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]+0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = ser_3x_p2_surfx2_eval_quad_node_7_r(fedge); 
  } else { 
    fUpOrdL[7] = ser_3x_p2_surfx2_eval_quad_node_7_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.6*alphaL[6]+0.4472135954999572*alphaL[4]+0.9*alphaL[3]+0.6708203932499357*alphaL[2]+0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = ser_3x_p2_surfx2_eval_quad_node_8_r(fedge); 
  } else { 
    fUpOrdL[8] = ser_3x_p2_surfx2_eval_quad_node_8_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[8] = {0.};
  ser_3x_p2_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[8] = {0.}; 
  GhatL[0] = 0.5*(alphaL[6]*fUpL[6]+alphaL[4]*fUpL[4]+alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.447213595499958*(alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6])+0.4472135954999579*(alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4])+0.5*(alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.447213595499958*alphaL[3]*fUpL[7]+0.5000000000000001*(alphaL[4]*fUpL[6]+fUpL[4]*alphaL[6])+0.4472135954999579*alphaL[2]*fUpL[5]+0.5*(alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.4*alphaL[6]*fUpL[7]+0.447213595499958*(alphaL[2]*fUpL[7]+alphaL[1]*fUpL[6]+fUpL[1]*alphaL[6])+0.4472135954999579*(alphaL[3]*(fUpL[5]+fUpL[4])+fUpL[3]*alphaL[4])+0.5*(alphaL[0]*fUpL[3]+fUpL[0]*alphaL[3]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[4] = 0.31943828249997*alphaL[6]*fUpL[6]+0.5000000000000001*(alphaL[2]*fUpL[6]+fUpL[2]*alphaL[6])+0.31943828249997*alphaL[4]*fUpL[4]+0.5*(alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4])+0.4472135954999579*(alphaL[3]*fUpL[3]+alphaL[1]*fUpL[1]); 
  GhatL[5] = 0.5000000000000001*alphaL[1]*fUpL[7]+0.4472135954999579*alphaL[6]*fUpL[6]+0.5*alphaL[0]*fUpL[5]+0.4472135954999579*(alphaL[3]*fUpL[3]+alphaL[2]*fUpL[2]); 
  GhatL[6] = 0.4*alphaL[3]*fUpL[7]+(0.31943828249997*alphaL[4]+0.5*alphaL[0])*fUpL[6]+(0.4472135954999579*fUpL[5]+0.31943828249997*fUpL[4]+0.5*fUpL[0])*alphaL[6]+0.5000000000000001*(alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4])+0.447213595499958*(alphaL[1]*fUpL[3]+fUpL[1]*alphaL[3]); 
  GhatL[7] = (0.4472135954999579*alphaL[4]+0.5*alphaL[0])*fUpL[7]+0.4*(alphaL[3]*fUpL[6]+fUpL[3]*alphaL[6])+0.5000000000000001*alphaL[1]*fUpL[5]+0.447213595499958*(alphaL[2]*fUpL[3]+fUpL[2]*alphaL[3]); 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[4] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[5] += 0.7071067811865475*GhatL[3]*rdvpar2; 
  out[6] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[7] += 0.7071067811865475*GhatL[4]*rdvpar2; 
  out[8] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[9] += 0.7071067811865475*GhatL[5]*rdvpar2; 
  out[10] += -1.224744871391589*GhatL[3]*rdvpar2; 
  out[11] += -1.224744871391589*GhatL[4]*rdvpar2; 
  out[12] += 1.58113883008419*GhatL[1]*rdvpar2; 
  out[13] += 0.7071067811865475*GhatL[6]*rdvpar2; 
  out[14] += 1.58113883008419*GhatL[2]*rdvpar2; 
  out[15] += 0.7071067811865475*GhatL[7]*rdvpar2; 
  out[16] += -1.224744871391589*GhatL[5]*rdvpar2; 
  out[17] += -1.224744871391589*GhatL[6]*rdvpar2; 
  out[18] += 1.58113883008419*GhatL[3]*rdvpar2; 
  out[19] += -1.224744871391589*GhatL[7]*rdvpar2; 

  } 

  return 2.5*rdvpar2*cflFreq; 

} 
