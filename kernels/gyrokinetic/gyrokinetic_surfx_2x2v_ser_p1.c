#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_2x2v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_gkhyb_2x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential .
  // apar: parallel component of magnetic vector potential.
  // apardot: time derivative of Apar.
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

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[4];
  const double *b_z = &b_i[8];

  double hamil[24] = {0.}; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[16] = (0.5962847939999438*m_)/rdvpar2Sq; 

  double BstarXdBmag[24] = {0.}; 
  BstarXdBmag[0] = 0.8660254037844386*((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[3]+(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[2])*rdy2; 
  BstarXdBmag[1] = 0.1732050807568877*((9.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0])*apar[3]+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[2])*rdy2; 

  double alphaL[12] = {0.}; 
  alphaL[0] = (0.125*((((12.72792206135786*b_z[1]-7.348469228349534*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_z[0]-7.348469228349534*b_z[1]))*hamil[5]+((4.242640687119286*b_z[0]-7.348469228349534*b_z[1])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_z[1]-2.449489742783178*b_z[0]))*hamil[2])*m_*rdy2+(2.449489742783178*BstarXdBmag[0]-4.242640687119286*BstarXdBmag[1])*hamil[3]*q_*rdvpar2))/(m_*q_); 
  alphaL[2] = -(0.125*(9.48683298050514*BstarXdBmag[1]-5.477225575051662*BstarXdBmag[0])*hamil[16]*rdvpar2)/m_; 

  double alphaR[12] = {0.}; 
  alphaR[0] = -(0.125*((((12.72792206135786*b_z[1]+7.348469228349534*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(7.348469228349534*b_z[1]+4.242640687119286*b_z[0]))*hamil[5]+((7.348469228349534*b_z[1]+4.242640687119286*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_z[1]+2.449489742783178*b_z[0]))*hamil[2])*m_*rdy2+((-4.242640687119286*BstarXdBmag[1])-2.449489742783178*BstarXdBmag[0])*hamil[3]*q_*rdvpar2))/(m_*q_); 
  alphaR[2] = (0.125*(9.48683298050514*BstarXdBmag[1]+5.477225575051662*BstarXdBmag[0])*hamil[16]*rdvpar2)/m_; 

  double cflFreq = 0.0;
  double fUpOrdL[12] = {0.};
  double alphaL_n = 0.;

  alphaL_n = 0.3535533905932734*alphaL[0]-0.4743416490252568*alphaL[2];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_2x2v_p1_surfx1_eval_quad_node_0_r(fl); 
  } else { 
    fUpOrdL[0] = gkhyb_2x2v_p1_surfx1_eval_quad_node_0_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_2x2v_p1_surfx1_eval_quad_node_1_r(fl); 
  } else { 
    fUpOrdL[1] = gkhyb_2x2v_p1_surfx1_eval_quad_node_1_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = gkhyb_2x2v_p1_surfx1_eval_quad_node_2_r(fl); 
  } else { 
    fUpOrdL[2] = gkhyb_2x2v_p1_surfx1_eval_quad_node_2_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0]-0.4743416490252568*alphaL[2];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = gkhyb_2x2v_p1_surfx1_eval_quad_node_3_r(fl); 
  } else { 
    fUpOrdL[3] = gkhyb_2x2v_p1_surfx1_eval_quad_node_3_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = gkhyb_2x2v_p1_surfx1_eval_quad_node_4_r(fl); 
  } else { 
    fUpOrdL[4] = gkhyb_2x2v_p1_surfx1_eval_quad_node_4_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = gkhyb_2x2v_p1_surfx1_eval_quad_node_5_r(fl); 
  } else { 
    fUpOrdL[5] = gkhyb_2x2v_p1_surfx1_eval_quad_node_5_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0]-0.4743416490252568*alphaL[2];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = gkhyb_2x2v_p1_surfx1_eval_quad_node_6_r(fl); 
  } else { 
    fUpOrdL[6] = gkhyb_2x2v_p1_surfx1_eval_quad_node_6_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = gkhyb_2x2v_p1_surfx1_eval_quad_node_7_r(fl); 
  } else { 
    fUpOrdL[7] = gkhyb_2x2v_p1_surfx1_eval_quad_node_7_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = gkhyb_2x2v_p1_surfx1_eval_quad_node_8_r(fl); 
  } else { 
    fUpOrdL[8] = gkhyb_2x2v_p1_surfx1_eval_quad_node_8_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0]-0.4743416490252568*alphaL[2];
  if (alphaL_n > 0.) {
    fUpOrdL[9] = gkhyb_2x2v_p1_surfx1_eval_quad_node_9_r(fl); 
  } else { 
    fUpOrdL[9] = gkhyb_2x2v_p1_surfx1_eval_quad_node_9_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[10] = gkhyb_2x2v_p1_surfx1_eval_quad_node_10_r(fl); 
  } else { 
    fUpOrdL[10] = gkhyb_2x2v_p1_surfx1_eval_quad_node_10_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.4743416490252568*alphaL[2]+0.3535533905932734*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[11] = gkhyb_2x2v_p1_surfx1_eval_quad_node_11_r(fl); 
  } else { 
    fUpOrdL[11] = gkhyb_2x2v_p1_surfx1_eval_quad_node_11_l(fc); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[12] = {0.};
  gkhyb_2x2v_p1_xdir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[12] = {0.}; 
  GhatL[0] = 0.3535533905932737*alphaL[2]*fUpL[2]+0.3535533905932737*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.3535533905932737*alphaL[2]*fUpL[4]+0.3535533905932737*alphaL[0]*fUpL[1]; 
  GhatL[2] = 0.3162277660168379*alphaL[2]*fUpL[8]+0.3535533905932737*alphaL[0]*fUpL[2]+0.3535533905932737*fUpL[0]*alphaL[2]; 
  GhatL[3] = 0.3535533905932737*alphaL[2]*fUpL[6]+0.3535533905932737*alphaL[0]*fUpL[3]; 
  GhatL[4] = 0.3162277660168379*alphaL[2]*fUpL[9]+0.3535533905932737*alphaL[0]*fUpL[4]+0.3535533905932737*fUpL[1]*alphaL[2]; 
  GhatL[5] = 0.3535533905932737*alphaL[2]*fUpL[7]+0.3535533905932737*alphaL[0]*fUpL[5]; 
  GhatL[6] = 0.3162277660168379*alphaL[2]*fUpL[10]+0.3535533905932737*alphaL[0]*fUpL[6]+0.3535533905932737*alphaL[2]*fUpL[3]; 
  GhatL[7] = 0.3162277660168379*alphaL[2]*fUpL[11]+0.3535533905932737*alphaL[0]*fUpL[7]+0.3535533905932737*alphaL[2]*fUpL[5]; 
  GhatL[8] = 0.3535533905932737*alphaL[0]*fUpL[8]+0.3162277660168379*alphaL[2]*fUpL[2]; 
  GhatL[9] = 0.3535533905932737*alphaL[0]*fUpL[9]+0.3162277660168379*alphaL[2]*fUpL[4]; 
  GhatL[10] = 0.3535533905932737*alphaL[0]*fUpL[10]+0.3162277660168379*alphaL[2]*fUpL[6]; 
  GhatL[11] = 0.3535533905932737*alphaL[0]*fUpL[11]+0.3162277660168379*alphaL[2]*fUpL[7]; 

  double fUpOrdR[12] = {0.};
  double alphaR_n = 0.;

  alphaR_n = 0.3535533905932734*alphaR[0]-0.4743416490252568*alphaR[2];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_2x2v_p1_surfx1_eval_quad_node_0_r(fc); 
  } else { 
    fUpOrdR[0] = gkhyb_2x2v_p1_surfx1_eval_quad_node_0_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_2x2v_p1_surfx1_eval_quad_node_1_r(fc); 
  } else { 
    fUpOrdR[1] = gkhyb_2x2v_p1_surfx1_eval_quad_node_1_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = gkhyb_2x2v_p1_surfx1_eval_quad_node_2_r(fc); 
  } else { 
    fUpOrdR[2] = gkhyb_2x2v_p1_surfx1_eval_quad_node_2_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0]-0.4743416490252568*alphaR[2];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = gkhyb_2x2v_p1_surfx1_eval_quad_node_3_r(fc); 
  } else { 
    fUpOrdR[3] = gkhyb_2x2v_p1_surfx1_eval_quad_node_3_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = gkhyb_2x2v_p1_surfx1_eval_quad_node_4_r(fc); 
  } else { 
    fUpOrdR[4] = gkhyb_2x2v_p1_surfx1_eval_quad_node_4_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = gkhyb_2x2v_p1_surfx1_eval_quad_node_5_r(fc); 
  } else { 
    fUpOrdR[5] = gkhyb_2x2v_p1_surfx1_eval_quad_node_5_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0]-0.4743416490252568*alphaR[2];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = gkhyb_2x2v_p1_surfx1_eval_quad_node_6_r(fc); 
  } else { 
    fUpOrdR[6] = gkhyb_2x2v_p1_surfx1_eval_quad_node_6_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = gkhyb_2x2v_p1_surfx1_eval_quad_node_7_r(fc); 
  } else { 
    fUpOrdR[7] = gkhyb_2x2v_p1_surfx1_eval_quad_node_7_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = gkhyb_2x2v_p1_surfx1_eval_quad_node_8_r(fc); 
  } else { 
    fUpOrdR[8] = gkhyb_2x2v_p1_surfx1_eval_quad_node_8_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0]-0.4743416490252568*alphaR[2];
  if (alphaR_n > 0.) {
    fUpOrdR[9] = gkhyb_2x2v_p1_surfx1_eval_quad_node_9_r(fc); 
  } else { 
    fUpOrdR[9] = gkhyb_2x2v_p1_surfx1_eval_quad_node_9_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[10] = gkhyb_2x2v_p1_surfx1_eval_quad_node_10_r(fc); 
  } else { 
    fUpOrdR[10] = gkhyb_2x2v_p1_surfx1_eval_quad_node_10_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[11] = gkhyb_2x2v_p1_surfx1_eval_quad_node_11_r(fc); 
  } else { 
    fUpOrdR[11] = gkhyb_2x2v_p1_surfx1_eval_quad_node_11_l(fr); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[12] = {0.};
  gkhyb_2x2v_p1_xdir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[12] = {0.}; 
  GhatR[0] = 0.3535533905932737*alphaR[2]*fUpR[2]+0.3535533905932737*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.3535533905932737*alphaR[2]*fUpR[4]+0.3535533905932737*alphaR[0]*fUpR[1]; 
  GhatR[2] = 0.3162277660168379*alphaR[2]*fUpR[8]+0.3535533905932737*alphaR[0]*fUpR[2]+0.3535533905932737*fUpR[0]*alphaR[2]; 
  GhatR[3] = 0.3535533905932737*alphaR[2]*fUpR[6]+0.3535533905932737*alphaR[0]*fUpR[3]; 
  GhatR[4] = 0.3162277660168379*alphaR[2]*fUpR[9]+0.3535533905932737*alphaR[0]*fUpR[4]+0.3535533905932737*fUpR[1]*alphaR[2]; 
  GhatR[5] = 0.3535533905932737*alphaR[2]*fUpR[7]+0.3535533905932737*alphaR[0]*fUpR[5]; 
  GhatR[6] = 0.3162277660168379*alphaR[2]*fUpR[10]+0.3535533905932737*alphaR[0]*fUpR[6]+0.3535533905932737*alphaR[2]*fUpR[3]; 
  GhatR[7] = 0.3162277660168379*alphaR[2]*fUpR[11]+0.3535533905932737*alphaR[0]*fUpR[7]+0.3535533905932737*alphaR[2]*fUpR[5]; 
  GhatR[8] = 0.3535533905932737*alphaR[0]*fUpR[8]+0.3162277660168379*alphaR[2]*fUpR[2]; 
  GhatR[9] = 0.3535533905932737*alphaR[0]*fUpR[9]+0.3162277660168379*alphaR[2]*fUpR[4]; 
  GhatR[10] = 0.3535533905932737*alphaR[0]*fUpR[10]+0.3162277660168379*alphaR[2]*fUpR[6]; 
  GhatR[11] = 0.3535533905932737*alphaR[0]*fUpR[11]+0.3162277660168379*alphaR[2]*fUpR[7]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdx2; 
  out[1] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdx2; 
  out[2] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdx2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdx2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdx2; 
  out[5] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdx2; 
  out[6] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdx2; 
  out[7] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdx2; 
  out[8] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdx2; 
  out[9] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdx2; 
  out[10] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdx2; 
  out[11] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdx2; 
  out[12] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdx2; 
  out[13] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdx2; 
  out[14] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdx2; 
  out[15] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdx2; 
  out[16] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdx2; 
  out[17] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdx2; 
  out[18] += (0.7071067811865475*GhatL[9]-0.7071067811865475*GhatR[9])*rdx2; 
  out[19] += (0.7071067811865475*GhatL[10]-0.7071067811865475*GhatR[10])*rdx2; 
  out[20] += ((-1.224744871391589*GhatR[9])-1.224744871391589*GhatL[9])*rdx2; 
  out[21] += ((-1.224744871391589*GhatR[10])-1.224744871391589*GhatL[10])*rdx2; 
  out[22] += (0.7071067811865475*GhatL[11]-0.7071067811865475*GhatR[11])*rdx2; 
  out[23] += ((-1.224744871391589*GhatR[11])-1.224744871391589*GhatL[11])*rdx2; 

  return 3.0*rdx2*cflFreq; 

} 
