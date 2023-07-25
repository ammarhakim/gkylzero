#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_1x2v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_gkhyb_1x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  double hamil[12] = {0.}; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+2.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 

  double BstarZdBmag[12] = {0.}; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*jacobtot_inv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_))/q_; 
  BstarZdBmag[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobtot_inv[1]*m_*rdx2*wvpar+(cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*jacobtot_inv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*b_y[1]*jacobtot_inv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0;

  if (edge == -1) { 

  double alphaR[6] = {0.}; 
  alphaR[0] = (0.25*((6.708203932499369*BstarZdBmag[4]+3.872983346207417*BstarZdBmag[2])*hamil[8]+(3.0*BstarZdBmag[1]+1.732050807568877*BstarZdBmag[0])*hamil[2])*rdvpar2)/m_; 
  alphaR[1] = (0.25*((6.708203932499369*BstarZdBmag[1]+3.872983346207417*BstarZdBmag[0])*hamil[8]+hamil[2]*(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2]))*rdvpar2)/m_; 
  alphaR[4] = (0.5*(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[2])*hamil[8]*rdvpar2)/m_; 

  double fUpOrdR[6] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.4472135954999572*alphaR[4]-0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_1x2v_p1_surfx1_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = gkhyb_1x2v_p1_surfx1_eval_quad_node_0_l(fedge); 
  } 
  cflFreq += -0.1875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0]-0.5590169943749465*alphaR[4];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_1x2v_p1_surfx1_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = gkhyb_1x2v_p1_surfx1_eval_quad_node_1_l(fedge); 
  } 
  cflFreq += -0.1875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]+0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = gkhyb_1x2v_p1_surfx1_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = gkhyb_1x2v_p1_surfx1_eval_quad_node_2_l(fedge); 
  } 
  cflFreq += -0.1875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]-0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = gkhyb_1x2v_p1_surfx1_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = gkhyb_1x2v_p1_surfx1_eval_quad_node_3_l(fedge); 
  } 
  cflFreq += -0.1875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.5*alphaR[0]-0.5590169943749465*alphaR[4];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = gkhyb_1x2v_p1_surfx1_eval_quad_node_4_r(fskin); 
  } else { 
    fUpOrdR[4] = gkhyb_1x2v_p1_surfx1_eval_quad_node_4_l(fedge); 
  } 
  cflFreq += -0.1875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.4472135954999572*alphaR[4]+0.6708203932499357*alphaR[1]+0.5*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = gkhyb_1x2v_p1_surfx1_eval_quad_node_5_r(fskin); 
  } else { 
    fUpOrdR[5] = gkhyb_1x2v_p1_surfx1_eval_quad_node_5_l(fedge); 
  } 
  cflFreq += -0.1875*rdx2*(alphaR_n-fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[6] = {0.};
  gkhyb_1x2v_p1_xdir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[6] = {0.}; 
  GhatR[0] = 0.5*alphaR[4]*fUpR[4]+0.5*alphaR[1]*fUpR[1]+0.5*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.4472135954999579*alphaR[1]*fUpR[4]+0.4472135954999579*fUpR[1]*alphaR[4]+0.5*alphaR[0]*fUpR[1]+0.5*fUpR[0]*alphaR[1]; 
  GhatR[2] = 0.5000000000000001*alphaR[4]*fUpR[5]+0.5*alphaR[1]*fUpR[3]+0.5*alphaR[0]*fUpR[2]; 
  GhatR[3] = 0.447213595499958*alphaR[1]*fUpR[5]+0.4472135954999579*fUpR[3]*alphaR[4]+0.5*alphaR[0]*fUpR[3]+0.5*alphaR[1]*fUpR[2]; 
  GhatR[4] = 0.31943828249997*alphaR[4]*fUpR[4]+0.5*alphaR[0]*fUpR[4]+0.5*fUpR[0]*alphaR[4]+0.4472135954999579*alphaR[1]*fUpR[1]; 
  GhatR[5] = 0.31943828249997*alphaR[4]*fUpR[5]+0.5*alphaR[0]*fUpR[5]+0.5000000000000001*fUpR[2]*alphaR[4]+0.447213595499958*alphaR[1]*fUpR[3]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdx2; 
  out[1] += -1.224744871391589*GhatR[0]*rdx2; 
  out[2] += -0.7071067811865475*GhatR[1]*rdx2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdx2; 
  out[4] += -1.224744871391589*GhatR[1]*rdx2; 
  out[5] += -1.224744871391589*GhatR[2]*rdx2; 
  out[6] += -0.7071067811865475*GhatR[3]*rdx2; 
  out[7] += -1.224744871391589*GhatR[3]*rdx2; 
  out[8] += -0.7071067811865475*GhatR[4]*rdx2; 
  out[9] += -1.224744871391589*GhatR[4]*rdx2; 
  out[10] += -0.7071067811865475*GhatR[5]*rdx2; 
  out[11] += -1.224744871391589*GhatR[5]*rdx2; 

  } else { 

  double alphaL[6] = {0.}; 
  alphaL[0] = -(0.25*((6.708203932499369*BstarZdBmag[4]-3.872983346207417*BstarZdBmag[2])*hamil[8]+(3.0*BstarZdBmag[1]-1.732050807568877*BstarZdBmag[0])*hamil[2])*rdvpar2)/m_; 
  alphaL[1] = -(0.25*((6.708203932499369*BstarZdBmag[1]-3.872983346207417*BstarZdBmag[0])*hamil[8]+hamil[2]*(3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[2]))*rdvpar2)/m_; 
  alphaL[4] = -(0.5*(3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[2])*hamil[8]*rdvpar2)/m_; 

  double fUpOrdL[6] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.4472135954999572*alphaL[4]-0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_1x2v_p1_surfx1_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = gkhyb_1x2v_p1_surfx1_eval_quad_node_0_l(fskin); 
  } 
  cflFreq += -0.1875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0]-0.5590169943749465*alphaL[4];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_1x2v_p1_surfx1_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = gkhyb_1x2v_p1_surfx1_eval_quad_node_1_l(fskin); 
  } 
  cflFreq += -0.1875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]+0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = gkhyb_1x2v_p1_surfx1_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = gkhyb_1x2v_p1_surfx1_eval_quad_node_2_l(fskin); 
  } 
  cflFreq += -0.1875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]-0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = gkhyb_1x2v_p1_surfx1_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = gkhyb_1x2v_p1_surfx1_eval_quad_node_3_l(fskin); 
  } 
  cflFreq += -0.1875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.5*alphaL[0]-0.5590169943749465*alphaL[4];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = gkhyb_1x2v_p1_surfx1_eval_quad_node_4_r(fedge); 
  } else { 
    fUpOrdL[4] = gkhyb_1x2v_p1_surfx1_eval_quad_node_4_l(fskin); 
  } 
  cflFreq += -0.1875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.4472135954999572*alphaL[4]+0.6708203932499357*alphaL[1]+0.5*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = gkhyb_1x2v_p1_surfx1_eval_quad_node_5_r(fedge); 
  } else { 
    fUpOrdL[5] = gkhyb_1x2v_p1_surfx1_eval_quad_node_5_l(fskin); 
  } 
  cflFreq += -0.1875*rdx2*(alphaL_n-fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[6] = {0.};
  gkhyb_1x2v_p1_xdir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[6] = {0.}; 
  GhatL[0] = 0.5*alphaL[4]*fUpL[4]+0.5*alphaL[1]*fUpL[1]+0.5*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.4472135954999579*alphaL[1]*fUpL[4]+0.4472135954999579*fUpL[1]*alphaL[4]+0.5*alphaL[0]*fUpL[1]+0.5*fUpL[0]*alphaL[1]; 
  GhatL[2] = 0.5000000000000001*alphaL[4]*fUpL[5]+0.5*alphaL[1]*fUpL[3]+0.5*alphaL[0]*fUpL[2]; 
  GhatL[3] = 0.447213595499958*alphaL[1]*fUpL[5]+0.4472135954999579*fUpL[3]*alphaL[4]+0.5*alphaL[0]*fUpL[3]+0.5*alphaL[1]*fUpL[2]; 
  GhatL[4] = 0.31943828249997*alphaL[4]*fUpL[4]+0.5*alphaL[0]*fUpL[4]+0.5*fUpL[0]*alphaL[4]+0.4472135954999579*alphaL[1]*fUpL[1]; 
  GhatL[5] = 0.31943828249997*alphaL[4]*fUpL[5]+0.5*alphaL[0]*fUpL[5]+0.5000000000000001*fUpL[2]*alphaL[4]+0.447213595499958*alphaL[1]*fUpL[3]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdx2; 
  out[1] += -1.224744871391589*GhatL[0]*rdx2; 
  out[2] += 0.7071067811865475*GhatL[1]*rdx2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdx2; 
  out[4] += -1.224744871391589*GhatL[1]*rdx2; 
  out[5] += -1.224744871391589*GhatL[2]*rdx2; 
  out[6] += 0.7071067811865475*GhatL[3]*rdx2; 
  out[7] += -1.224744871391589*GhatL[3]*rdx2; 
  out[8] += 0.7071067811865475*GhatL[4]*rdx2; 
  out[9] += -1.224744871391589*GhatL[4]*rdx2; 
  out[10] += 0.7071067811865475*GhatL[5]*rdx2; 
  out[11] += -1.224744871391589*GhatL[5]*rdx2; 

  } 

  return cflFreq; 

} 
