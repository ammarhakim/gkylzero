#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_ser_3x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_3x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential .
  // apar: parallel component of magnetic vector potential.
  // apardot: time derivative of Apar.
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

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[3];
  const double *b_z = &b_i[6];

  double hamil[20] = {0.}; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+2.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[7] = 2.0*(bmag[2]*wmu+phi[2]*q_); 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 
  hamil[13] = (1.154700538379251*bmag[2])/rdmu2; 

  double BstarZdBmag[20] = {0.}; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobtot_inv[2]+cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_))/q_; 
  BstarZdBmag[1] = (0.2828427124746191*(1.732050807568877*(b_y[2]*(10.0*jacobtot_inv[2]+11.18033988749895*jacobtot_inv[0])+5.0*b_y[1]*jacobtot_inv[1])*m_*rdx2*wvpar+(4.47213595499958*(cmag[1]*jacobtot_inv[2]+jacobtot_inv[1]*cmag[2])+5.0*(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*(2.23606797749979*jacobtot_inv[1]*b_y[2]+jacobtot_inv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*(b_y[2]*(2.0*jacobtot_inv[2]+2.23606797749979*jacobtot_inv[0])+b_y[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2*wvpar+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobtot_inv[2]+7.0*(5.0*jacobtot_inv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobtot_inv[1]))*q_))/q_; 
  BstarZdBmag[11] = (1.414213562373095*(b_y[1]*jacobtot_inv[2]+2.0*jacobtot_inv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0;

  if (edge == -1) { 

  double alphaR[8] = {0.}; 
  alphaR[0] = (0.25*(hamil[8]*(8.660254037844387*BstarZdBmag[11]+6.708203932499369*BstarZdBmag[4]+3.872983346207417*BstarZdBmag[2])+hamil[2]*(3.872983346207417*BstarZdBmag[7]+3.0*BstarZdBmag[1]+1.732050807568877*BstarZdBmag[0]))*rdvpar2)/m_; 
  alphaR[1] = (0.25*(3.872983346207417*hamil[2]*BstarZdBmag[11]+(8.660254037844386*BstarZdBmag[7]+6.708203932499369*BstarZdBmag[1]+3.872983346207417*BstarZdBmag[0])*hamil[8]+hamil[2]*(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2]))*rdvpar2)/m_; 
  alphaR[4] = (0.5*hamil[8]*(3.872983346207417*BstarZdBmag[11]+3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double fUpOrdR[9] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.4472135954999572*alphaR[4]-0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = ser_3x_p2_surfx1_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = ser_3x_p2_surfx1_eval_quad_node_0_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]-0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = ser_3x_p2_surfx1_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = ser_3x_p2_surfx1_eval_quad_node_1_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]-0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = ser_3x_p2_surfx1_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = ser_3x_p2_surfx1_eval_quad_node_2_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0]-0.5590169943749465*alphaR[4];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = ser_3x_p2_surfx1_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = ser_3x_p2_surfx1_eval_quad_node_3_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0]-0.5590169943749465*alphaR[4];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = ser_3x_p2_surfx1_eval_quad_node_4_r(fskin); 
  } else { 
    fUpOrdR[4] = ser_3x_p2_surfx1_eval_quad_node_4_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0]-0.5590169943749465*alphaR[4];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = ser_3x_p2_surfx1_eval_quad_node_5_r(fskin); 
  } else { 
    fUpOrdR[5] = ser_3x_p2_surfx1_eval_quad_node_5_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]+0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = ser_3x_p2_surfx1_eval_quad_node_6_r(fskin); 
  } else { 
    fUpOrdR[6] = ser_3x_p2_surfx1_eval_quad_node_6_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]+0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = ser_3x_p2_surfx1_eval_quad_node_7_r(fskin); 
  } else { 
    fUpOrdR[7] = ser_3x_p2_surfx1_eval_quad_node_7_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]+0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = ser_3x_p2_surfx1_eval_quad_node_8_r(fskin); 
  } else { 
    fUpOrdR[8] = ser_3x_p2_surfx1_eval_quad_node_8_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[8] = {0.};
  ser_3x_p2_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[8] = {0.}; 
  GhatR[0] = 0.5*alphaR[4]*fUpR[4]+0.5*alphaR[1]*fUpR[1]+0.5*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.4472135954999579*alphaR[1]*fUpR[4]+0.4472135954999579*fUpR[1]*alphaR[4]+0.5*alphaR[0]*fUpR[1]+0.5*fUpR[0]*alphaR[1]; 
  GhatR[2] = 0.5000000000000001*alphaR[4]*fUpR[6]+0.5*alphaR[1]*fUpR[3]+0.5*alphaR[0]*fUpR[2]; 
  GhatR[3] = 0.447213595499958*alphaR[1]*fUpR[6]+0.4472135954999579*fUpR[3]*alphaR[4]+0.5*alphaR[0]*fUpR[3]+0.5*alphaR[1]*fUpR[2]; 
  GhatR[4] = 0.31943828249997*alphaR[4]*fUpR[4]+0.5*alphaR[0]*fUpR[4]+0.5*fUpR[0]*alphaR[4]+0.4472135954999579*alphaR[1]*fUpR[1]; 
  GhatR[5] = 0.5000000000000001*alphaR[1]*fUpR[7]+0.5*alphaR[0]*fUpR[5]; 
  GhatR[6] = 0.31943828249997*alphaR[4]*fUpR[6]+0.5*alphaR[0]*fUpR[6]+0.5000000000000001*fUpR[2]*alphaR[4]+0.447213595499958*alphaR[1]*fUpR[3]; 
  GhatR[7] = 0.4472135954999579*alphaR[4]*fUpR[7]+0.5*alphaR[0]*fUpR[7]+0.5000000000000001*alphaR[1]*fUpR[5]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdx2; 
  out[1] += -1.224744871391589*GhatR[0]*rdx2; 
  out[2] += -0.7071067811865475*GhatR[1]*rdx2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdx2; 
  out[4] += -1.224744871391589*GhatR[1]*rdx2; 
  out[5] += -1.224744871391589*GhatR[2]*rdx2; 
  out[6] += -0.7071067811865475*GhatR[3]*rdx2; 
  out[7] += -1.58113883008419*GhatR[0]*rdx2; 
  out[8] += -0.7071067811865475*GhatR[4]*rdx2; 
  out[9] += -0.7071067811865475*GhatR[5]*rdx2; 
  out[10] += -1.224744871391589*GhatR[3]*rdx2; 
  out[11] += -1.58113883008419*GhatR[1]*rdx2; 
  out[12] += -1.224744871391589*GhatR[4]*rdx2; 
  out[13] += -1.58113883008419*GhatR[2]*rdx2; 
  out[14] += -0.7071067811865475*GhatR[6]*rdx2; 
  out[15] += -1.224744871391589*GhatR[5]*rdx2; 
  out[16] += -0.7071067811865475*GhatR[7]*rdx2; 
  out[17] += -1.58113883008419*GhatR[3]*rdx2; 
  out[18] += -1.224744871391589*GhatR[6]*rdx2; 
  out[19] += -1.224744871391589*GhatR[7]*rdx2; 

  } else { 

  double alphaL[8] = {0.}; 
  alphaL[0] = (0.25*(hamil[8]*(8.660254037844387*BstarZdBmag[11]-6.708203932499369*BstarZdBmag[4]+3.872983346207417*BstarZdBmag[2])+hamil[2]*(3.872983346207417*BstarZdBmag[7]-3.0*BstarZdBmag[1]+1.732050807568877*BstarZdBmag[0]))*rdvpar2)/m_; 
  alphaL[1] = (0.25*(3.872983346207417*hamil[2]*BstarZdBmag[11]+(8.660254037844386*BstarZdBmag[7]-6.708203932499369*BstarZdBmag[1]+3.872983346207417*BstarZdBmag[0])*hamil[8]+hamil[2]*(1.732050807568877*BstarZdBmag[2]-3.0*BstarZdBmag[4]))*rdvpar2)/m_; 
  alphaL[4] = (0.5*hamil[8]*(3.872983346207417*BstarZdBmag[11]-3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2])*rdvpar2)/m_; 

  double fUpOrdL[9] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.4472135954999572*alphaL[4]-0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = ser_3x_p2_surfx1_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = ser_3x_p2_surfx1_eval_quad_node_0_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]-0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = ser_3x_p2_surfx1_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = ser_3x_p2_surfx1_eval_quad_node_1_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]-0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = ser_3x_p2_surfx1_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = ser_3x_p2_surfx1_eval_quad_node_2_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0]-0.5590169943749465*alphaL[4];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = ser_3x_p2_surfx1_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = ser_3x_p2_surfx1_eval_quad_node_3_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0]-0.5590169943749465*alphaL[4];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = ser_3x_p2_surfx1_eval_quad_node_4_r(fedge); 
  } else { 
    fUpOrdL[4] = ser_3x_p2_surfx1_eval_quad_node_4_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0]-0.5590169943749465*alphaL[4];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = ser_3x_p2_surfx1_eval_quad_node_5_r(fedge); 
  } else { 
    fUpOrdL[5] = ser_3x_p2_surfx1_eval_quad_node_5_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]+0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = ser_3x_p2_surfx1_eval_quad_node_6_r(fedge); 
  } else { 
    fUpOrdL[6] = ser_3x_p2_surfx1_eval_quad_node_6_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]+0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = ser_3x_p2_surfx1_eval_quad_node_7_r(fedge); 
  } else { 
    fUpOrdL[7] = ser_3x_p2_surfx1_eval_quad_node_7_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]+0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = ser_3x_p2_surfx1_eval_quad_node_8_r(fedge); 
  } else { 
    fUpOrdL[8] = ser_3x_p2_surfx1_eval_quad_node_8_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[8] = {0.};
  ser_3x_p2_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[8] = {0.}; 
  GhatL[0] = 0.5*alphaL[4]*fUpL[4]+0.5*alphaL[1]*fUpL[1]+0.5*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.4472135954999579*alphaL[1]*fUpL[4]+0.4472135954999579*fUpL[1]*alphaL[4]+0.5*alphaL[0]*fUpL[1]+0.5*fUpL[0]*alphaL[1]; 
  GhatL[2] = 0.5000000000000001*alphaL[4]*fUpL[6]+0.5*alphaL[1]*fUpL[3]+0.5*alphaL[0]*fUpL[2]; 
  GhatL[3] = 0.447213595499958*alphaL[1]*fUpL[6]+0.4472135954999579*fUpL[3]*alphaL[4]+0.5*alphaL[0]*fUpL[3]+0.5*alphaL[1]*fUpL[2]; 
  GhatL[4] = 0.31943828249997*alphaL[4]*fUpL[4]+0.5*alphaL[0]*fUpL[4]+0.5*fUpL[0]*alphaL[4]+0.4472135954999579*alphaL[1]*fUpL[1]; 
  GhatL[5] = 0.5000000000000001*alphaL[1]*fUpL[7]+0.5*alphaL[0]*fUpL[5]; 
  GhatL[6] = 0.31943828249997*alphaL[4]*fUpL[6]+0.5*alphaL[0]*fUpL[6]+0.5000000000000001*fUpL[2]*alphaL[4]+0.447213595499958*alphaL[1]*fUpL[3]; 
  GhatL[7] = 0.4472135954999579*alphaL[4]*fUpL[7]+0.5*alphaL[0]*fUpL[7]+0.5000000000000001*alphaL[1]*fUpL[5]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdx2; 
  out[1] += -1.224744871391589*GhatL[0]*rdx2; 
  out[2] += 0.7071067811865475*GhatL[1]*rdx2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdx2; 
  out[4] += -1.224744871391589*GhatL[1]*rdx2; 
  out[5] += -1.224744871391589*GhatL[2]*rdx2; 
  out[6] += 0.7071067811865475*GhatL[3]*rdx2; 
  out[7] += 1.58113883008419*GhatL[0]*rdx2; 
  out[8] += 0.7071067811865475*GhatL[4]*rdx2; 
  out[9] += 0.7071067811865475*GhatL[5]*rdx2; 
  out[10] += -1.224744871391589*GhatL[3]*rdx2; 
  out[11] += 1.58113883008419*GhatL[1]*rdx2; 
  out[12] += -1.224744871391589*GhatL[4]*rdx2; 
  out[13] += 1.58113883008419*GhatL[2]*rdx2; 
  out[14] += 0.7071067811865475*GhatL[6]*rdx2; 
  out[15] += -1.224744871391589*GhatL[5]*rdx2; 
  out[16] += 0.7071067811865475*GhatL[7]*rdx2; 
  out[17] += 1.58113883008419*GhatL[3]*rdx2; 
  out[18] += -1.224744871391589*GhatL[6]*rdx2; 
  out[19] += -1.224744871391589*GhatL[7]*rdx2; 

  } 

  return 5.0*rdx2*cflFreq; 

} 
