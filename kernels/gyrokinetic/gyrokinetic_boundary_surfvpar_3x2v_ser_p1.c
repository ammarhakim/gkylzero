#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_surfx4_eval_quad.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  double hamil[48] = {0.}; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(4.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+4.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*(bmag[3]*wmu+phi[3]*q_); 
  hamil[4] = (3.265986323710906*m_*wvpar)/rdvpar2; 
  hamil[5] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*(bmag[5]*wmu+phi[5]*q_); 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[12] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[14] = (1.154700538379252*bmag[3])/rdmu2; 
  hamil[16] = 2.0*phi[7]*q_; 
  hamil[21] = (1.154700538379252*bmag[5])/rdmu2; 
  hamil[32] = (0.8432740427115681*m_)/rdvpar2Sq; 

  double BstarXdBmag[48] = {0.}; 
  BstarXdBmag[0] = -(0.3061862178478971*(4.0*(jacobtot_inv[1]*b_y[5]+jacobtot_inv[0]*b_y[3])*m_*rdz2*wvpar-1.414213562373095*((b_z[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_z[5]+b_z[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_z[3])*apar[7]+(b_z[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_z[5]+b_z[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_z[3])*apar[6]+(apar[2]*b_z[5]+b_z[3]*apar[4])*jacobtot_inv[5]+apar[4]*(jacobtot_inv[3]*b_z[5]+b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])+apar[2]*(b_z[3]*jacobtot_inv[3]+b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0]))*q_*rdy2))/q_; 
  BstarXdBmag[1] = -(0.06123724356957942*(20.0*(jacobtot_inv[0]*b_y[5]+jacobtot_inv[1]*b_y[3])*m_*rdz2*wvpar-1.414213562373095*((9.0*(b_z[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_z[5])+5.0*(b_z[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_z[3]))*apar[7]+5.0*(b_z[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_z[5]+b_z[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_z[3])*apar[6]+(9.0*apar[4]*b_z[5]+5.0*apar[2]*b_z[3])*jacobtot_inv[5]+5.0*apar[2]*jacobtot_inv[3]*b_z[5]+(5.0*b_z[3]*jacobtot_inv[3]+9.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0])*apar[4]+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[2])*q_*rdy2))/q_; 
  BstarXdBmag[3] = -(0.06123724356957942*(20.0*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])*m_*rdz2*wvpar-1.414213562373095*((9.0*(b_z[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_z[5])+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[7]+(9.0*(b_z[5]*jacobtot_inv[5]+b_z[3]*jacobtot_inv[3])+5.0*(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0]))*apar[6]+5.0*((b_z[0]*apar[4]+b_z[1]*apar[2])*jacobtot_inv[5]+(jacobtot_inv[0]*apar[4]+jacobtot_inv[1]*apar[2])*b_z[5]+(b_z[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_z[3])*apar[4]+apar[2]*(b_z[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_z[3])))*q_*rdy2))/q_; 
  BstarXdBmag[4] = -(0.7071067811865475*(jacobtot_inv[1]*b_y[5]+jacobtot_inv[0]*b_y[3])*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[7] = -(0.01224744871391589*(100.0*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])*m_*rdz2*wvpar-1.414213562373095*((81.0*b_z[5]*jacobtot_inv[5]+5.0*(9.0*(b_z[3]*jacobtot_inv[3]+b_z[1]*jacobtot_inv[1])+5.0*b_z[0]*jacobtot_inv[0]))*apar[7]+5.0*((9.0*(b_z[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_z[5])+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[6]+(9.0*b_z[1]*apar[4]+5.0*b_z[0]*apar[2])*jacobtot_inv[5]+(9.0*jacobtot_inv[1]*apar[4]+5.0*jacobtot_inv[0]*apar[2])*b_z[5]+5.0*((b_z[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_z[3])*apar[4]+apar[2]*(b_z[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_z[3]))))*q_*rdy2))/q_; 
  BstarXdBmag[9] = -(0.7071067811865475*(jacobtot_inv[0]*b_y[5]+jacobtot_inv[1]*b_y[3])*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[11] = -(0.7071067811865475*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[18] = -(0.7071067811865475*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])*m_*rdz2)/(q_*rdvpar2); 

  double BstarYdBmag[48] = {0.}; 
  BstarYdBmag[0] = (0.3061862178478971*(4.0*m_*((jacobtot_inv[1]*b_x[5]+jacobtot_inv[0]*b_x[3])*rdz2-1.0*(jacobtot_inv[3]*b_z[5]+jacobtot_inv[0]*b_z[1])*rdx2)*wvpar-1.414213562373095*(2.0*(apar[1]*b_z[5]+b_z[1]*apar[5])*jacobtot_inv[5]+(2.0*jacobtot_inv[1]*apar[5]+apar[0]*jacobtot_inv[3]+jacobtot_inv[0]*apar[3])*b_z[5]+(b_z[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_z[3])*apar[5]+(apar[1]*b_z[3]+b_z[1]*apar[3])*jacobtot_inv[3]+2.0*apar[1]*b_z[1]*jacobtot_inv[1]+jacobtot_inv[0]*(apar[0]*b_z[1]+b_z[0]*apar[1]))*q_*rdx2))/q_; 
  BstarYdBmag[1] = (0.3061862178478971*(4.0*m_*((jacobtot_inv[0]*b_x[5]+jacobtot_inv[1]*b_x[3])*rdz2-1.0*(b_z[5]*jacobtot_inv[5]+b_z[1]*jacobtot_inv[1])*rdx2)*wvpar-1.414213562373095*((apar[0]*b_z[5]+b_z[0]*apar[5]+apar[1]*b_z[3]+b_z[1]*apar[3])*jacobtot_inv[5]+(2.0*(jacobtot_inv[0]*apar[5]+apar[1]*jacobtot_inv[3])+jacobtot_inv[1]*apar[3])*b_z[5]+(2.0*b_z[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_z[3])*apar[5]+(apar[0]*b_z[1]+b_z[0]*apar[1])*jacobtot_inv[1]+2.0*jacobtot_inv[0]*apar[1]*b_z[1])*q_*rdx2))/q_; 
  BstarYdBmag[2] = -0.4330127018922193*((2.0*(b_z[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_z[5])+b_z[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_z[3])*apar[7]+(jacobtot_inv[0]*b_z[5]+b_z[1]*jacobtot_inv[3])*apar[6]+b_z[5]*(2.0*apar[4]*jacobtot_inv[5]+apar[2]*jacobtot_inv[3])+(b_z[3]*jacobtot_inv[3]+2.0*b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[4]+jacobtot_inv[0]*b_z[1]*apar[2])*rdx2; 
  BstarYdBmag[3] = (0.06123724356957942*(20.0*m_*((b_x[5]*jacobtot_inv[5]+b_x[3]*jacobtot_inv[3])*rdz2-1.0*(jacobtot_inv[0]*b_z[5]+b_z[1]*jacobtot_inv[3])*rdx2)*wvpar-1.414213562373095*(2.0*(9.0*apar[5]*b_z[5]+5.0*apar[1]*b_z[1])*jacobtot_inv[5]+(9.0*apar[3]*jacobtot_inv[3]+5.0*(2.0*apar[1]*jacobtot_inv[1]+apar[0]*jacobtot_inv[0]))*b_z[5]+(9.0*b_z[3]*jacobtot_inv[3]+5.0*(2.0*b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0]))*apar[5]+5.0*((apar[0]*b_z[1]+b_z[0]*apar[1])*jacobtot_inv[3]+jacobtot_inv[0]*(apar[1]*b_z[3]+b_z[1]*apar[3])))*q_*rdx2))/q_; 
  BstarYdBmag[4] = (0.7071067811865475*m_*((jacobtot_inv[1]*b_x[5]+jacobtot_inv[0]*b_x[3])*rdz2-1.0*(jacobtot_inv[3]*b_z[5]+jacobtot_inv[0]*b_z[1])*rdx2))/(q_*rdvpar2); 
  BstarYdBmag[6] = -0.4330127018922193*((b_z[0]*jacobtot_inv[5]+2.0*(jacobtot_inv[0]*b_z[5]+b_z[1]*jacobtot_inv[3])+jacobtot_inv[1]*b_z[3])*apar[7]+(b_z[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_z[5])*apar[6]+(apar[2]*b_z[5]+b_z[3]*apar[4])*jacobtot_inv[5]+apar[4]*(2.0*jacobtot_inv[3]*b_z[5]+b_z[0]*jacobtot_inv[1]+2.0*jacobtot_inv[0]*b_z[1])+b_z[1]*jacobtot_inv[1]*apar[2])*rdx2; 
  BstarYdBmag[7] = (0.06123724356957942*(20.0*m_*((b_x[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_x[5])*rdz2-1.0*(b_z[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_z[5])*rdx2)*wvpar-1.414213562373095*((9.0*(apar[3]*b_z[5]+b_z[3]*apar[5])+5.0*(apar[0]*b_z[1]+b_z[0]*apar[1]))*jacobtot_inv[5]+(18.0*jacobtot_inv[3]*apar[5]+5.0*(apar[0]*jacobtot_inv[1]+2.0*jacobtot_inv[0]*apar[1]))*b_z[5]+5.0*((b_z[0]*jacobtot_inv[1]+2.0*jacobtot_inv[0]*b_z[1])*apar[5]+2.0*apar[1]*b_z[1]*jacobtot_inv[3]+jacobtot_inv[1]*(apar[1]*b_z[3]+b_z[1]*apar[3])))*q_*rdx2))/q_; 
  BstarYdBmag[8] = -0.08660254037844387*((9.0*(2.0*b_z[5]*jacobtot_inv[5]+b_z[3]*jacobtot_inv[3])+5.0*(2.0*b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0]))*apar[7]+(9.0*jacobtot_inv[3]*b_z[5]+5.0*jacobtot_inv[0]*b_z[1])*apar[6]+5.0*(2.0*b_z[1]*apar[4]*jacobtot_inv[5]+(2.0*jacobtot_inv[1]*apar[4]+jacobtot_inv[0]*apar[2])*b_z[5]+(b_z[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_z[3])*apar[4]+b_z[1]*apar[2]*jacobtot_inv[3]))*rdx2; 
  BstarYdBmag[9] = (0.7071067811865475*m_*((jacobtot_inv[0]*b_x[5]+jacobtot_inv[1]*b_x[3])*rdz2-1.0*(b_z[5]*jacobtot_inv[5]+b_z[1]*jacobtot_inv[1])*rdx2))/(q_*rdvpar2); 
  BstarYdBmag[11] = (0.7071067811865475*m_*((b_x[5]*jacobtot_inv[5]+b_x[3]*jacobtot_inv[3])*rdz2-1.0*(jacobtot_inv[0]*b_z[5]+b_z[1]*jacobtot_inv[3])*rdx2))/(q_*rdvpar2); 
  BstarYdBmag[16] = -0.08660254037844387*((9.0*(b_z[3]*jacobtot_inv[5]+2.0*jacobtot_inv[3]*b_z[5])+5.0*(b_z[0]*jacobtot_inv[1]+2.0*jacobtot_inv[0]*b_z[1]))*apar[7]+(9.0*b_z[5]*jacobtot_inv[5]+5.0*b_z[1]*jacobtot_inv[1])*apar[6]+5.0*((b_z[0]*apar[4]+b_z[1]*apar[2])*jacobtot_inv[5]+(2.0*jacobtot_inv[0]*apar[4]+jacobtot_inv[1]*apar[2])*b_z[5]+(2.0*b_z[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_z[3])*apar[4]))*rdx2; 
  BstarYdBmag[18] = (0.7071067811865475*m_*((b_x[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_x[5])*rdz2-1.0*(b_z[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_z[5])*rdx2))/(q_*rdvpar2); 

  double BstarZdBmag[48] = {0.}; 
  BstarZdBmag[0] = (0.7071067811865475*(1.732050807568877*(jacobtot_inv[3]*b_y[5]+jacobtot_inv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[5]*jacobtot_inv[5]+cmag[3]*jacobtot_inv[3]+cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_))/q_; 
  BstarZdBmag[1] = (0.7071067811865475*(1.732050807568877*(b_y[5]*jacobtot_inv[5]+b_y[1]*jacobtot_inv[1])*m_*rdx2*wvpar+(cmag[3]*jacobtot_inv[5]+jacobtot_inv[3]*cmag[5]+cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[3] = (0.7071067811865475*(1.732050807568877*(jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3])*m_*rdx2*wvpar+(cmag[1]*jacobtot_inv[5]+jacobtot_inv[1]*cmag[5]+cmag[0]*jacobtot_inv[3]+jacobtot_inv[0]*cmag[3])*q_))/q_; 
  BstarZdBmag[4] = (0.7071067811865475*(jacobtot_inv[3]*b_y[5]+jacobtot_inv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[7] = (0.7071067811865475*(1.732050807568877*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5])*m_*rdx2*wvpar+(cmag[0]*jacobtot_inv[5]+jacobtot_inv[0]*cmag[5]+cmag[1]*jacobtot_inv[3]+jacobtot_inv[1]*cmag[3])*q_))/q_; 
  BstarZdBmag[9] = (0.7071067811865475*(b_y[5]*jacobtot_inv[5]+b_y[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[11] = (0.7071067811865475*(jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[18] = (0.7071067811865475*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0;

  if (edge == -1) { 

  double alphaR[16] = {0.}; 
  alphaR[0] = -(0.125*((hamil[7]*(3.0*BstarZdBmag[9]+1.732050807568877*BstarZdBmag[1])+hamil[3]*(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[0]))*rdz2+(hamil[16]*(3.0*BstarYdBmag[18]+1.732050807568877*BstarYdBmag[7])+3.0*(hamil[8]*BstarYdBmag[11]+hamil[6]*BstarYdBmag[9])+1.732050807568877*(BstarYdBmag[3]*hamil[8]+BstarYdBmag[1]*hamil[6])+hamil[2]*(3.0*BstarYdBmag[4]+1.732050807568877*BstarYdBmag[0]))*rdy2+(hamil[7]*(3.0*BstarXdBmag[11]+1.732050807568877*BstarXdBmag[3])+hamil[1]*(3.0*BstarXdBmag[4]+1.732050807568877*BstarXdBmag[0]))*rdx2+11.31370849898477*apardot[0]*q_))/m_; 
  alphaR[1] = -(0.125*((3.0*hamil[3]*BstarZdBmag[9]+(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[0])*hamil[7]+1.732050807568877*BstarZdBmag[1]*hamil[3])*rdz2+(3.0*hamil[8]*BstarYdBmag[18]+(3.0*BstarYdBmag[11]+1.732050807568877*BstarYdBmag[3])*hamil[16]+3.0*hamil[2]*BstarYdBmag[9]+1.732050807568877*BstarYdBmag[7]*hamil[8]+(3.0*BstarYdBmag[4]+1.732050807568877*BstarYdBmag[0])*hamil[6]+1.732050807568877*BstarYdBmag[1]*hamil[2])*rdy2+(3.0*(hamil[7]*BstarXdBmag[18]+hamil[1]*BstarXdBmag[9])+1.732050807568877*(BstarXdBmag[7]*hamil[7]+BstarXdBmag[1]*hamil[1]))*rdx2+11.31370849898477*apardot[1]*q_))/m_; 
  alphaR[2] = -(0.125*(((3.0*BstarZdBmag[9]+1.732050807568877*BstarZdBmag[1])*hamil[16]+(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[0])*hamil[8])*rdz2+1.732050807568877*(BstarYdBmag[16]*hamil[16]+BstarYdBmag[8]*hamil[8]+BstarYdBmag[6]*hamil[6]+BstarYdBmag[2]*hamil[2])*rdy2+((3.0*BstarXdBmag[11]+1.732050807568877*BstarXdBmag[3])*hamil[16]+(3.0*BstarXdBmag[4]+1.732050807568877*BstarXdBmag[0])*hamil[6])*rdx2+11.31370849898477*apardot[2]*q_))/m_; 
  alphaR[3] = -(0.125*((3.0*(hamil[7]*BstarZdBmag[18]+hamil[3]*BstarZdBmag[11])+1.732050807568877*(BstarZdBmag[7]*hamil[7]+BstarZdBmag[3]*hamil[3]))*rdz2+(3.0*hamil[6]*BstarYdBmag[18]+(3.0*BstarYdBmag[9]+1.732050807568877*BstarYdBmag[1])*hamil[16]+3.0*hamil[2]*BstarYdBmag[11]+(3.0*BstarYdBmag[4]+1.732050807568877*BstarYdBmag[0])*hamil[8]+1.732050807568877*(hamil[6]*BstarYdBmag[7]+hamil[2]*BstarYdBmag[3]))*rdy2+(3.0*hamil[1]*BstarXdBmag[11]+(3.0*BstarXdBmag[4]+1.732050807568877*BstarXdBmag[0])*hamil[7]+1.732050807568877*hamil[1]*BstarXdBmag[3])*rdx2+11.31370849898477*apardot[3]*q_))/m_; 
  alphaR[4] = -(0.125*(((3.0*BstarZdBmag[9]+1.732050807568877*BstarZdBmag[1])*hamil[21]+(3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[0])*hamil[14])*rdz2+((3.0*BstarXdBmag[11]+1.732050807568877*BstarXdBmag[3])*hamil[21]+(3.0*BstarXdBmag[4]+1.732050807568877*BstarXdBmag[0])*hamil[12])*rdx2))/m_; 
  alphaR[5] = -(0.125*(((3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[0])*hamil[16]+hamil[8]*(3.0*BstarZdBmag[9]+1.732050807568877*BstarZdBmag[1]))*rdz2+1.732050807568877*(BstarYdBmag[8]*hamil[16]+hamil[8]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[6]+hamil[2]*BstarYdBmag[6])*rdy2+(hamil[16]*(3.0*BstarXdBmag[18]+1.732050807568877*BstarXdBmag[7])+hamil[6]*(3.0*BstarXdBmag[9]+1.732050807568877*BstarXdBmag[1]))*rdx2+11.31370849898477*apardot[4]*q_))/m_; 
  alphaR[6] = -(0.125*((3.0*(hamil[3]*BstarZdBmag[18]+hamil[7]*BstarZdBmag[11])+1.732050807568877*(BstarZdBmag[3]*hamil[7]+hamil[3]*BstarZdBmag[7]))*rdz2+(3.0*hamil[2]*BstarYdBmag[18]+(3.0*BstarYdBmag[4]+1.732050807568877*BstarYdBmag[0])*hamil[16]+3.0*(hamil[6]*BstarYdBmag[11]+hamil[8]*BstarYdBmag[9])+1.732050807568877*(BstarYdBmag[1]*hamil[8]+hamil[2]*BstarYdBmag[7]+BstarYdBmag[3]*hamil[6]))*rdy2+(3.0*(hamil[1]*BstarXdBmag[18]+hamil[7]*BstarXdBmag[9])+1.732050807568877*(BstarXdBmag[1]*hamil[7]+hamil[1]*BstarXdBmag[7]))*rdx2+11.31370849898477*apardot[5]*q_))/m_; 
  alphaR[7] = -(0.125*((hamil[16]*(3.0*BstarZdBmag[18]+1.732050807568877*BstarZdBmag[7])+hamil[8]*(3.0*BstarZdBmag[11]+1.732050807568877*BstarZdBmag[3]))*rdz2+1.732050807568877*(BstarYdBmag[6]*hamil[16]+hamil[6]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[8]+hamil[2]*BstarYdBmag[8])*rdy2+((3.0*BstarXdBmag[4]+1.732050807568877*BstarXdBmag[0])*hamil[16]+hamil[6]*(3.0*BstarXdBmag[11]+1.732050807568877*BstarXdBmag[3]))*rdx2+11.31370849898477*apardot[6]*q_))/m_; 
  alphaR[8] = -(0.125*(((3.0*BstarZdBmag[4]+1.732050807568877*BstarZdBmag[0])*hamil[21]+(3.0*BstarZdBmag[9]+1.732050807568877*BstarZdBmag[1])*hamil[14])*rdz2+((3.0*BstarXdBmag[18]+1.732050807568877*BstarXdBmag[7])*hamil[21]+(3.0*BstarXdBmag[9]+1.732050807568877*BstarXdBmag[1])*hamil[12])*rdx2))/m_; 
  alphaR[10] = -(0.125*(((3.0*BstarZdBmag[18]+1.732050807568877*BstarZdBmag[7])*hamil[21]+(3.0*BstarZdBmag[11]+1.732050807568877*BstarZdBmag[3])*hamil[14])*rdz2+((3.0*BstarXdBmag[4]+1.732050807568877*BstarXdBmag[0])*hamil[21]+(3.0*BstarXdBmag[11]+1.732050807568877*BstarXdBmag[3])*hamil[12])*rdx2))/m_; 
  alphaR[11] = -(0.125*((3.0*hamil[8]*BstarZdBmag[18]+(3.0*BstarZdBmag[11]+1.732050807568877*BstarZdBmag[3])*hamil[16]+1.732050807568877*BstarZdBmag[7]*hamil[8])*rdz2+1.732050807568877*(BstarYdBmag[2]*hamil[16]+hamil[2]*BstarYdBmag[16]+BstarYdBmag[6]*hamil[8]+hamil[6]*BstarYdBmag[8])*rdy2+(3.0*hamil[6]*BstarXdBmag[18]+(3.0*BstarXdBmag[9]+1.732050807568877*BstarXdBmag[1])*hamil[16]+1.732050807568877*hamil[6]*BstarXdBmag[7])*rdx2+11.31370849898477*apardot[7]*q_))/m_; 
  alphaR[13] = -(0.125*(((3.0*BstarZdBmag[11]+1.732050807568877*BstarZdBmag[3])*hamil[21]+hamil[14]*(3.0*BstarZdBmag[18]+1.732050807568877*BstarZdBmag[7]))*rdz2+((3.0*BstarXdBmag[9]+1.732050807568877*BstarXdBmag[1])*hamil[21]+hamil[12]*(3.0*BstarXdBmag[18]+1.732050807568877*BstarXdBmag[7]))*rdx2))/m_; 

  double fUpOrdR[16] = {0.};
  double alphaR_n = 0.;

  alphaR_n = (-0.25*alphaR[13])-0.25*alphaR[11]+0.25*alphaR[10]+0.25*alphaR[8]+0.25*alphaR[7]+0.25*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]-0.25*alphaR[11]-0.25*alphaR[10]-0.25*alphaR[8]+0.25*alphaR[7]+0.25*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]+0.25*alphaR[11]-0.25*alphaR[10]+0.25*alphaR[8]-0.25*alphaR[7]-0.25*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])+0.25*alphaR[11]+0.25*alphaR[10]-0.25*alphaR[8]-0.25*alphaR[7]-0.25*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])+0.25*alphaR[11]+0.25*alphaR[10]+0.25*alphaR[8]-0.25*alphaR[7]+0.25*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_r(fskin); 
  } else { 
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]+0.25*alphaR[11]-0.25*alphaR[10]-0.25*alphaR[8]-0.25*alphaR[7]+0.25*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_r(fskin); 
  } else { 
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]-0.25*alphaR[11]-0.25*alphaR[10]+0.25*alphaR[8]+0.25*alphaR[7]-0.25*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_r(fskin); 
  } else { 
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])-0.25*alphaR[11]+0.25*alphaR[10]-0.25*alphaR[8]+0.25*alphaR[7]-0.25*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_r(fskin); 
  } else { 
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]+0.25*alphaR[11]+0.25*alphaR[10]-0.25*alphaR[8]+0.25*alphaR[7]-0.25*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_r(fskin); 
  } else { 
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])+0.25*alphaR[11]-0.25*alphaR[10]+0.25*alphaR[8]+0.25*alphaR[7]-0.25*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_r(fskin); 
  } else { 
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])-0.25*alphaR[11]-0.25*alphaR[10]-0.25*alphaR[8]-0.25*alphaR[7]+0.25*alphaR[6]-0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_r(fskin); 
  } else { 
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]-0.25*alphaR[11]+0.25*alphaR[10]+0.25*alphaR[8]-0.25*alphaR[7]+0.25*alphaR[6]-0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_r(fskin); 
  } else { 
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]-0.25*alphaR[11]+0.25*alphaR[10]-0.25*alphaR[8]-0.25*alphaR[7]-0.25*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_r(fskin); 
  } else { 
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])-0.25*alphaR[11]-0.25*alphaR[10]+0.25*alphaR[8]-0.25*alphaR[7]-0.25*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_r(fskin); 
  } else { 
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = (-0.25*alphaR[13])+0.25*alphaR[11]-0.25*alphaR[10]-0.25*alphaR[8]+0.25*alphaR[7]+0.25*alphaR[6]+0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_r(fskin); 
  } else { 
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 
  alphaR_n = 0.25*alphaR[13]+0.25*alphaR[11]+0.25*alphaR[10]+0.25*alphaR[8]+0.25*alphaR[7]+0.25*alphaR[6]+0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_r(fskin); 
  } else { 
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_l(fedge); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[16] = {0.};
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[16] = {0.}; 
  GhatR[0] = 0.25*alphaR[13]*fUpR[13]+0.25*alphaR[11]*fUpR[11]+0.25*alphaR[10]*fUpR[10]+0.25*alphaR[8]*fUpR[8]+0.25*alphaR[7]*fUpR[7]+0.25*alphaR[6]*fUpR[6]+0.25*alphaR[5]*fUpR[5]+0.25*alphaR[4]*fUpR[4]+0.25*alphaR[3]*fUpR[3]+0.25*alphaR[2]*fUpR[2]+0.25*alphaR[1]*fUpR[1]+0.25*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.25*alphaR[10]*fUpR[13]+0.25*fUpR[10]*alphaR[13]+0.25*alphaR[7]*fUpR[11]+0.25*fUpR[7]*alphaR[11]+0.25*alphaR[4]*fUpR[8]+0.25*fUpR[4]*alphaR[8]+0.25*alphaR[3]*fUpR[6]+0.25*fUpR[3]*alphaR[6]+0.25*alphaR[2]*fUpR[5]+0.25*fUpR[2]*alphaR[5]+0.25*alphaR[0]*fUpR[1]+0.25*fUpR[0]*alphaR[1]; 
  GhatR[2] = 0.25*alphaR[13]*fUpR[15]+0.25*alphaR[10]*fUpR[14]+0.25*alphaR[8]*fUpR[12]+0.25*alphaR[6]*fUpR[11]+0.25*fUpR[6]*alphaR[11]+0.25*alphaR[4]*fUpR[9]+0.25*alphaR[3]*fUpR[7]+0.25*fUpR[3]*alphaR[7]+0.25*alphaR[1]*fUpR[5]+0.25*fUpR[1]*alphaR[5]+0.25*alphaR[0]*fUpR[2]+0.25*fUpR[0]*alphaR[2]; 
  GhatR[3] = 0.25*alphaR[8]*fUpR[13]+0.25*fUpR[8]*alphaR[13]+0.25*alphaR[5]*fUpR[11]+0.25*fUpR[5]*alphaR[11]+0.25*alphaR[4]*fUpR[10]+0.25*fUpR[4]*alphaR[10]+0.25*alphaR[2]*fUpR[7]+0.25*fUpR[2]*alphaR[7]+0.25*alphaR[1]*fUpR[6]+0.25*fUpR[1]*alphaR[6]+0.25*alphaR[0]*fUpR[3]+0.25*fUpR[0]*alphaR[3]; 
  GhatR[4] = 0.25*alphaR[11]*fUpR[15]+0.25*alphaR[7]*fUpR[14]+0.25*alphaR[6]*fUpR[13]+0.25*fUpR[6]*alphaR[13]+0.25*alphaR[5]*fUpR[12]+0.25*alphaR[3]*fUpR[10]+0.25*fUpR[3]*alphaR[10]+0.25*alphaR[2]*fUpR[9]+0.25*alphaR[1]*fUpR[8]+0.25*fUpR[1]*alphaR[8]+0.25*alphaR[0]*fUpR[4]+0.25*fUpR[0]*alphaR[4]; 
  GhatR[5] = 0.25*alphaR[10]*fUpR[15]+0.25*alphaR[13]*fUpR[14]+0.25*alphaR[4]*fUpR[12]+0.25*alphaR[3]*fUpR[11]+0.25*fUpR[3]*alphaR[11]+0.25*alphaR[8]*fUpR[9]+0.25*alphaR[6]*fUpR[7]+0.25*fUpR[6]*alphaR[7]+0.25*alphaR[0]*fUpR[5]+0.25*fUpR[0]*alphaR[5]+0.25*alphaR[1]*fUpR[2]+0.25*fUpR[1]*alphaR[2]; 
  GhatR[6] = 0.25*alphaR[4]*fUpR[13]+0.25*fUpR[4]*alphaR[13]+0.25*alphaR[2]*fUpR[11]+0.25*fUpR[2]*alphaR[11]+0.25*alphaR[8]*fUpR[10]+0.25*fUpR[8]*alphaR[10]+0.25*alphaR[5]*fUpR[7]+0.25*fUpR[5]*alphaR[7]+0.25*alphaR[0]*fUpR[6]+0.25*fUpR[0]*alphaR[6]+0.25*alphaR[1]*fUpR[3]+0.25*fUpR[1]*alphaR[3]; 
  GhatR[7] = 0.25*alphaR[8]*fUpR[15]+0.25*alphaR[4]*fUpR[14]+0.25*fUpR[12]*alphaR[13]+0.25*alphaR[1]*fUpR[11]+0.25*fUpR[1]*alphaR[11]+0.25*fUpR[9]*alphaR[10]+0.25*alphaR[0]*fUpR[7]+0.25*fUpR[0]*alphaR[7]+0.25*alphaR[5]*fUpR[6]+0.25*fUpR[5]*alphaR[6]+0.25*alphaR[2]*fUpR[3]+0.25*fUpR[2]*alphaR[3]; 
  GhatR[8] = 0.25*alphaR[7]*fUpR[15]+0.25*alphaR[11]*fUpR[14]+0.25*alphaR[3]*fUpR[13]+0.25*fUpR[3]*alphaR[13]+0.25*alphaR[2]*fUpR[12]+0.25*alphaR[6]*fUpR[10]+0.25*fUpR[6]*alphaR[10]+0.25*alphaR[5]*fUpR[9]+0.25*alphaR[0]*fUpR[8]+0.25*fUpR[0]*alphaR[8]+0.25*alphaR[1]*fUpR[4]+0.25*fUpR[1]*alphaR[4]; 
  GhatR[9] = 0.25*alphaR[6]*fUpR[15]+0.25*alphaR[3]*fUpR[14]+0.25*alphaR[11]*fUpR[13]+0.25*fUpR[11]*alphaR[13]+0.25*alphaR[1]*fUpR[12]+0.25*alphaR[7]*fUpR[10]+0.25*fUpR[7]*alphaR[10]+0.25*alphaR[0]*fUpR[9]+0.25*alphaR[5]*fUpR[8]+0.25*fUpR[5]*alphaR[8]+0.25*alphaR[2]*fUpR[4]+0.25*fUpR[2]*alphaR[4]; 
  GhatR[10] = 0.25*alphaR[5]*fUpR[15]+0.25*alphaR[2]*fUpR[14]+0.25*alphaR[1]*fUpR[13]+0.25*fUpR[1]*alphaR[13]+0.25*alphaR[11]*fUpR[12]+0.25*alphaR[0]*fUpR[10]+0.25*fUpR[0]*alphaR[10]+0.25*alphaR[7]*fUpR[9]+0.25*alphaR[6]*fUpR[8]+0.25*fUpR[6]*alphaR[8]+0.25*alphaR[3]*fUpR[4]+0.25*fUpR[3]*alphaR[4]; 
  GhatR[11] = 0.25*alphaR[4]*fUpR[15]+0.25*alphaR[8]*fUpR[14]+0.25*fUpR[9]*alphaR[13]+0.25*alphaR[10]*fUpR[12]+0.25*alphaR[0]*fUpR[11]+0.25*fUpR[0]*alphaR[11]+0.25*alphaR[1]*fUpR[7]+0.25*fUpR[1]*alphaR[7]+0.25*alphaR[2]*fUpR[6]+0.25*fUpR[2]*alphaR[6]+0.25*alphaR[3]*fUpR[5]+0.25*fUpR[3]*alphaR[5]; 
  GhatR[12] = 0.25*alphaR[3]*fUpR[15]+0.25*alphaR[6]*fUpR[14]+0.25*alphaR[7]*fUpR[13]+0.25*fUpR[7]*alphaR[13]+0.25*alphaR[0]*fUpR[12]+0.25*alphaR[10]*fUpR[11]+0.25*fUpR[10]*alphaR[11]+0.25*alphaR[1]*fUpR[9]+0.25*alphaR[2]*fUpR[8]+0.25*fUpR[2]*alphaR[8]+0.25*alphaR[4]*fUpR[5]+0.25*fUpR[4]*alphaR[5]; 
  GhatR[13] = 0.25*alphaR[2]*fUpR[15]+0.25*alphaR[5]*fUpR[14]+0.25*alphaR[0]*fUpR[13]+0.25*fUpR[0]*alphaR[13]+0.25*alphaR[7]*fUpR[12]+0.25*fUpR[9]*alphaR[11]+0.25*alphaR[1]*fUpR[10]+0.25*fUpR[1]*alphaR[10]+0.25*alphaR[3]*fUpR[8]+0.25*fUpR[3]*alphaR[8]+0.25*alphaR[4]*fUpR[6]+0.25*fUpR[4]*alphaR[6]; 
  GhatR[14] = 0.25*alphaR[1]*fUpR[15]+0.25*alphaR[0]*fUpR[14]+0.25*alphaR[5]*fUpR[13]+0.25*fUpR[5]*alphaR[13]+0.25*alphaR[6]*fUpR[12]+0.25*alphaR[8]*fUpR[11]+0.25*fUpR[8]*alphaR[11]+0.25*alphaR[2]*fUpR[10]+0.25*fUpR[2]*alphaR[10]+0.25*alphaR[3]*fUpR[9]+0.25*alphaR[4]*fUpR[7]+0.25*fUpR[4]*alphaR[7]; 
  GhatR[15] = 0.25*alphaR[0]*fUpR[15]+0.25*alphaR[1]*fUpR[14]+0.25*alphaR[2]*fUpR[13]+0.25*fUpR[2]*alphaR[13]+0.25*alphaR[3]*fUpR[12]+0.25*alphaR[4]*fUpR[11]+0.25*fUpR[4]*alphaR[11]+0.25*alphaR[5]*fUpR[10]+0.25*fUpR[5]*alphaR[10]+0.25*alphaR[6]*fUpR[9]+0.25*alphaR[7]*fUpR[8]+0.25*fUpR[7]*alphaR[8]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdvpar2; 
  out[1] += -0.7071067811865475*GhatR[1]*rdvpar2; 
  out[2] += -0.7071067811865475*GhatR[2]*rdvpar2; 
  out[3] += -0.7071067811865475*GhatR[3]*rdvpar2; 
  out[4] += -1.224744871391589*GhatR[0]*rdvpar2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdvpar2; 
  out[6] += -0.7071067811865475*GhatR[5]*rdvpar2; 
  out[7] += -0.7071067811865475*GhatR[6]*rdvpar2; 
  out[8] += -0.7071067811865475*GhatR[7]*rdvpar2; 
  out[9] += -1.224744871391589*GhatR[1]*rdvpar2; 
  out[10] += -1.224744871391589*GhatR[2]*rdvpar2; 
  out[11] += -1.224744871391589*GhatR[3]*rdvpar2; 
  out[12] += -0.7071067811865475*GhatR[8]*rdvpar2; 
  out[13] += -0.7071067811865475*GhatR[9]*rdvpar2; 
  out[14] += -0.7071067811865475*GhatR[10]*rdvpar2; 
  out[15] += -1.224744871391589*GhatR[4]*rdvpar2; 
  out[16] += -0.7071067811865475*GhatR[11]*rdvpar2; 
  out[17] += -1.224744871391589*GhatR[5]*rdvpar2; 
  out[18] += -1.224744871391589*GhatR[6]*rdvpar2; 
  out[19] += -1.224744871391589*GhatR[7]*rdvpar2; 
  out[20] += -0.7071067811865475*GhatR[12]*rdvpar2; 
  out[21] += -0.7071067811865475*GhatR[13]*rdvpar2; 
  out[22] += -0.7071067811865475*GhatR[14]*rdvpar2; 
  out[23] += -1.224744871391589*GhatR[8]*rdvpar2; 
  out[24] += -1.224744871391589*GhatR[9]*rdvpar2; 
  out[25] += -1.224744871391589*GhatR[10]*rdvpar2; 
  out[26] += -1.224744871391589*GhatR[11]*rdvpar2; 
  out[27] += -0.7071067811865475*GhatR[15]*rdvpar2; 
  out[28] += -1.224744871391589*GhatR[12]*rdvpar2; 
  out[29] += -1.224744871391589*GhatR[13]*rdvpar2; 
  out[30] += -1.224744871391589*GhatR[14]*rdvpar2; 
  out[31] += -1.224744871391589*GhatR[15]*rdvpar2; 
  out[32] += -1.58113883008419*GhatR[0]*rdvpar2; 
  out[33] += -1.58113883008419*GhatR[1]*rdvpar2; 
  out[34] += -1.58113883008419*GhatR[2]*rdvpar2; 
  out[35] += -1.58113883008419*GhatR[3]*rdvpar2; 
  out[36] += -1.58113883008419*GhatR[4]*rdvpar2; 
  out[37] += -1.58113883008419*GhatR[5]*rdvpar2; 
  out[38] += -1.58113883008419*GhatR[6]*rdvpar2; 
  out[39] += -1.58113883008419*GhatR[7]*rdvpar2; 
  out[40] += -1.58113883008419*GhatR[8]*rdvpar2; 
  out[41] += -1.58113883008419*GhatR[9]*rdvpar2; 
  out[42] += -1.58113883008419*GhatR[10]*rdvpar2; 
  out[43] += -1.58113883008419*GhatR[11]*rdvpar2; 
  out[44] += -1.58113883008419*GhatR[12]*rdvpar2; 
  out[45] += -1.58113883008419*GhatR[13]*rdvpar2; 
  out[46] += -1.58113883008419*GhatR[14]*rdvpar2; 
  out[47] += -1.58113883008419*GhatR[15]*rdvpar2; 

  } else { 

  double alphaL[16] = {0.}; 
  alphaL[0] = (0.125*((hamil[7]*(3.0*BstarZdBmag[9]-1.732050807568877*BstarZdBmag[1])+hamil[3]*(3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[0]))*rdz2+(hamil[16]*(3.0*BstarYdBmag[18]-1.732050807568877*BstarYdBmag[7])+3.0*(hamil[8]*BstarYdBmag[11]+hamil[6]*BstarYdBmag[9])-1.732050807568877*(BstarYdBmag[3]*hamil[8]+BstarYdBmag[1]*hamil[6])+hamil[2]*(3.0*BstarYdBmag[4]-1.732050807568877*BstarYdBmag[0]))*rdy2+(hamil[7]*(3.0*BstarXdBmag[11]-1.732050807568877*BstarXdBmag[3])+hamil[1]*(3.0*BstarXdBmag[4]-1.732050807568877*BstarXdBmag[0]))*rdx2-11.31370849898477*apardot[0]*q_))/m_; 
  alphaL[1] = (0.125*((3.0*hamil[3]*BstarZdBmag[9]+(3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[0])*hamil[7]-1.732050807568877*BstarZdBmag[1]*hamil[3])*rdz2+(3.0*hamil[8]*BstarYdBmag[18]+(3.0*BstarYdBmag[11]-1.732050807568877*BstarYdBmag[3])*hamil[16]+3.0*hamil[2]*BstarYdBmag[9]-1.732050807568877*BstarYdBmag[7]*hamil[8]+(3.0*BstarYdBmag[4]-1.732050807568877*BstarYdBmag[0])*hamil[6]-1.732050807568877*BstarYdBmag[1]*hamil[2])*rdy2+(3.0*(hamil[7]*BstarXdBmag[18]+hamil[1]*BstarXdBmag[9])-1.732050807568877*(BstarXdBmag[7]*hamil[7]+BstarXdBmag[1]*hamil[1]))*rdx2-11.31370849898477*apardot[1]*q_))/m_; 
  alphaL[2] = (0.125*(((3.0*BstarZdBmag[9]-1.732050807568877*BstarZdBmag[1])*hamil[16]+(3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[0])*hamil[8])*rdz2-1.732050807568877*(BstarYdBmag[16]*hamil[16]+BstarYdBmag[8]*hamil[8]+BstarYdBmag[6]*hamil[6]+BstarYdBmag[2]*hamil[2])*rdy2+((3.0*BstarXdBmag[11]-1.732050807568877*BstarXdBmag[3])*hamil[16]+(3.0*BstarXdBmag[4]-1.732050807568877*BstarXdBmag[0])*hamil[6])*rdx2-11.31370849898477*apardot[2]*q_))/m_; 
  alphaL[3] = (0.125*((3.0*(hamil[7]*BstarZdBmag[18]+hamil[3]*BstarZdBmag[11])-1.732050807568877*(BstarZdBmag[7]*hamil[7]+BstarZdBmag[3]*hamil[3]))*rdz2+(3.0*hamil[6]*BstarYdBmag[18]+(3.0*BstarYdBmag[9]-1.732050807568877*BstarYdBmag[1])*hamil[16]+3.0*hamil[2]*BstarYdBmag[11]+(3.0*BstarYdBmag[4]-1.732050807568877*BstarYdBmag[0])*hamil[8]-1.732050807568877*(hamil[6]*BstarYdBmag[7]+hamil[2]*BstarYdBmag[3]))*rdy2+(3.0*hamil[1]*BstarXdBmag[11]+(3.0*BstarXdBmag[4]-1.732050807568877*BstarXdBmag[0])*hamil[7]-1.732050807568877*hamil[1]*BstarXdBmag[3])*rdx2-11.31370849898477*apardot[3]*q_))/m_; 
  alphaL[4] = (0.125*(((3.0*BstarZdBmag[9]-1.732050807568877*BstarZdBmag[1])*hamil[21]+(3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[0])*hamil[14])*rdz2+((3.0*BstarXdBmag[11]-1.732050807568877*BstarXdBmag[3])*hamil[21]+(3.0*BstarXdBmag[4]-1.732050807568877*BstarXdBmag[0])*hamil[12])*rdx2))/m_; 
  alphaL[5] = (0.125*(((3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[0])*hamil[16]+hamil[8]*(3.0*BstarZdBmag[9]-1.732050807568877*BstarZdBmag[1]))*rdz2-1.732050807568877*(BstarYdBmag[8]*hamil[16]+hamil[8]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[6]+hamil[2]*BstarYdBmag[6])*rdy2+(hamil[16]*(3.0*BstarXdBmag[18]-1.732050807568877*BstarXdBmag[7])+hamil[6]*(3.0*BstarXdBmag[9]-1.732050807568877*BstarXdBmag[1]))*rdx2-11.31370849898477*apardot[4]*q_))/m_; 
  alphaL[6] = (0.125*((3.0*(hamil[3]*BstarZdBmag[18]+hamil[7]*BstarZdBmag[11])-1.732050807568877*(BstarZdBmag[3]*hamil[7]+hamil[3]*BstarZdBmag[7]))*rdz2+(3.0*hamil[2]*BstarYdBmag[18]+(3.0*BstarYdBmag[4]-1.732050807568877*BstarYdBmag[0])*hamil[16]+3.0*(hamil[6]*BstarYdBmag[11]+hamil[8]*BstarYdBmag[9])-1.732050807568877*(BstarYdBmag[1]*hamil[8]+hamil[2]*BstarYdBmag[7]+BstarYdBmag[3]*hamil[6]))*rdy2+(3.0*(hamil[1]*BstarXdBmag[18]+hamil[7]*BstarXdBmag[9])-1.732050807568877*(BstarXdBmag[1]*hamil[7]+hamil[1]*BstarXdBmag[7]))*rdx2-11.31370849898477*apardot[5]*q_))/m_; 
  alphaL[7] = (0.125*((hamil[16]*(3.0*BstarZdBmag[18]-1.732050807568877*BstarZdBmag[7])+hamil[8]*(3.0*BstarZdBmag[11]-1.732050807568877*BstarZdBmag[3]))*rdz2-1.732050807568877*(BstarYdBmag[6]*hamil[16]+hamil[6]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[8]+hamil[2]*BstarYdBmag[8])*rdy2+((3.0*BstarXdBmag[4]-1.732050807568877*BstarXdBmag[0])*hamil[16]+hamil[6]*(3.0*BstarXdBmag[11]-1.732050807568877*BstarXdBmag[3]))*rdx2-11.31370849898477*apardot[6]*q_))/m_; 
  alphaL[8] = (0.125*(((3.0*BstarZdBmag[4]-1.732050807568877*BstarZdBmag[0])*hamil[21]+(3.0*BstarZdBmag[9]-1.732050807568877*BstarZdBmag[1])*hamil[14])*rdz2+((3.0*BstarXdBmag[18]-1.732050807568877*BstarXdBmag[7])*hamil[21]+(3.0*BstarXdBmag[9]-1.732050807568877*BstarXdBmag[1])*hamil[12])*rdx2))/m_; 
  alphaL[10] = (0.125*(((3.0*BstarZdBmag[18]-1.732050807568877*BstarZdBmag[7])*hamil[21]+(3.0*BstarZdBmag[11]-1.732050807568877*BstarZdBmag[3])*hamil[14])*rdz2+((3.0*BstarXdBmag[4]-1.732050807568877*BstarXdBmag[0])*hamil[21]+(3.0*BstarXdBmag[11]-1.732050807568877*BstarXdBmag[3])*hamil[12])*rdx2))/m_; 
  alphaL[11] = (0.125*((3.0*hamil[8]*BstarZdBmag[18]+(3.0*BstarZdBmag[11]-1.732050807568877*BstarZdBmag[3])*hamil[16]-1.732050807568877*BstarZdBmag[7]*hamil[8])*rdz2-1.732050807568877*(BstarYdBmag[2]*hamil[16]+hamil[2]*BstarYdBmag[16]+BstarYdBmag[6]*hamil[8]+hamil[6]*BstarYdBmag[8])*rdy2+(3.0*hamil[6]*BstarXdBmag[18]+(3.0*BstarXdBmag[9]-1.732050807568877*BstarXdBmag[1])*hamil[16]-1.732050807568877*hamil[6]*BstarXdBmag[7])*rdx2-11.31370849898477*apardot[7]*q_))/m_; 
  alphaL[13] = (0.125*(((3.0*BstarZdBmag[11]-1.732050807568877*BstarZdBmag[3])*hamil[21]+hamil[14]*(3.0*BstarZdBmag[18]-1.732050807568877*BstarZdBmag[7]))*rdz2+((3.0*BstarXdBmag[9]-1.732050807568877*BstarXdBmag[1])*hamil[21]+hamil[12]*(3.0*BstarXdBmag[18]-1.732050807568877*BstarXdBmag[7]))*rdx2))/m_; 

  double fUpOrdL[16] = {0.};
  double alphaL_n = 0.;

  alphaL_n = (-0.25*alphaL[13])-0.25*alphaL[11]+0.25*alphaL[10]+0.25*alphaL[8]+0.25*alphaL[7]+0.25*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx4_eval_quad_node_0_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]-0.25*alphaL[11]-0.25*alphaL[10]-0.25*alphaL[8]+0.25*alphaL[7]+0.25*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx4_eval_quad_node_1_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]+0.25*alphaL[11]-0.25*alphaL[10]+0.25*alphaL[8]-0.25*alphaL[7]-0.25*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx4_eval_quad_node_2_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])+0.25*alphaL[11]+0.25*alphaL[10]-0.25*alphaL[8]-0.25*alphaL[7]-0.25*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx4_eval_quad_node_3_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])+0.25*alphaL[11]+0.25*alphaL[10]+0.25*alphaL[8]-0.25*alphaL[7]+0.25*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_r(fedge); 
  } else { 
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx4_eval_quad_node_4_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]+0.25*alphaL[11]-0.25*alphaL[10]-0.25*alphaL[8]-0.25*alphaL[7]+0.25*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_r(fedge); 
  } else { 
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx4_eval_quad_node_5_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]-0.25*alphaL[11]-0.25*alphaL[10]+0.25*alphaL[8]+0.25*alphaL[7]-0.25*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_r(fedge); 
  } else { 
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx4_eval_quad_node_6_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])-0.25*alphaL[11]+0.25*alphaL[10]-0.25*alphaL[8]+0.25*alphaL[7]-0.25*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_r(fedge); 
  } else { 
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx4_eval_quad_node_7_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]+0.25*alphaL[11]+0.25*alphaL[10]-0.25*alphaL[8]+0.25*alphaL[7]-0.25*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_r(fedge); 
  } else { 
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx4_eval_quad_node_8_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])+0.25*alphaL[11]-0.25*alphaL[10]+0.25*alphaL[8]+0.25*alphaL[7]-0.25*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_r(fedge); 
  } else { 
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx4_eval_quad_node_9_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])-0.25*alphaL[11]-0.25*alphaL[10]-0.25*alphaL[8]-0.25*alphaL[7]+0.25*alphaL[6]-0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_r(fedge); 
  } else { 
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx4_eval_quad_node_10_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]-0.25*alphaL[11]+0.25*alphaL[10]+0.25*alphaL[8]-0.25*alphaL[7]+0.25*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_r(fedge); 
  } else { 
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx4_eval_quad_node_11_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]-0.25*alphaL[11]+0.25*alphaL[10]-0.25*alphaL[8]-0.25*alphaL[7]-0.25*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_r(fedge); 
  } else { 
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx4_eval_quad_node_12_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])-0.25*alphaL[11]-0.25*alphaL[10]+0.25*alphaL[8]-0.25*alphaL[7]-0.25*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_r(fedge); 
  } else { 
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx4_eval_quad_node_13_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = (-0.25*alphaL[13])+0.25*alphaL[11]-0.25*alphaL[10]-0.25*alphaL[8]+0.25*alphaL[7]+0.25*alphaL[6]+0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_r(fedge); 
  } else { 
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx4_eval_quad_node_14_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 
  alphaL_n = 0.25*alphaL[13]+0.25*alphaL[11]+0.25*alphaL[10]+0.25*alphaL[8]+0.25*alphaL[7]+0.25*alphaL[6]+0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_r(fedge); 
  } else { 
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx4_eval_quad_node_15_l(fskin); 
  } 
  cflFreq = fmax(cflFreq, fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[16] = {0.};
  gkhyb_3x2v_p1_vpardir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[16] = {0.}; 
  GhatL[0] = 0.25*alphaL[13]*fUpL[13]+0.25*alphaL[11]*fUpL[11]+0.25*alphaL[10]*fUpL[10]+0.25*alphaL[8]*fUpL[8]+0.25*alphaL[7]*fUpL[7]+0.25*alphaL[6]*fUpL[6]+0.25*alphaL[5]*fUpL[5]+0.25*alphaL[4]*fUpL[4]+0.25*alphaL[3]*fUpL[3]+0.25*alphaL[2]*fUpL[2]+0.25*alphaL[1]*fUpL[1]+0.25*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.25*alphaL[10]*fUpL[13]+0.25*fUpL[10]*alphaL[13]+0.25*alphaL[7]*fUpL[11]+0.25*fUpL[7]*alphaL[11]+0.25*alphaL[4]*fUpL[8]+0.25*fUpL[4]*alphaL[8]+0.25*alphaL[3]*fUpL[6]+0.25*fUpL[3]*alphaL[6]+0.25*alphaL[2]*fUpL[5]+0.25*fUpL[2]*alphaL[5]+0.25*alphaL[0]*fUpL[1]+0.25*fUpL[0]*alphaL[1]; 
  GhatL[2] = 0.25*alphaL[13]*fUpL[15]+0.25*alphaL[10]*fUpL[14]+0.25*alphaL[8]*fUpL[12]+0.25*alphaL[6]*fUpL[11]+0.25*fUpL[6]*alphaL[11]+0.25*alphaL[4]*fUpL[9]+0.25*alphaL[3]*fUpL[7]+0.25*fUpL[3]*alphaL[7]+0.25*alphaL[1]*fUpL[5]+0.25*fUpL[1]*alphaL[5]+0.25*alphaL[0]*fUpL[2]+0.25*fUpL[0]*alphaL[2]; 
  GhatL[3] = 0.25*alphaL[8]*fUpL[13]+0.25*fUpL[8]*alphaL[13]+0.25*alphaL[5]*fUpL[11]+0.25*fUpL[5]*alphaL[11]+0.25*alphaL[4]*fUpL[10]+0.25*fUpL[4]*alphaL[10]+0.25*alphaL[2]*fUpL[7]+0.25*fUpL[2]*alphaL[7]+0.25*alphaL[1]*fUpL[6]+0.25*fUpL[1]*alphaL[6]+0.25*alphaL[0]*fUpL[3]+0.25*fUpL[0]*alphaL[3]; 
  GhatL[4] = 0.25*alphaL[11]*fUpL[15]+0.25*alphaL[7]*fUpL[14]+0.25*alphaL[6]*fUpL[13]+0.25*fUpL[6]*alphaL[13]+0.25*alphaL[5]*fUpL[12]+0.25*alphaL[3]*fUpL[10]+0.25*fUpL[3]*alphaL[10]+0.25*alphaL[2]*fUpL[9]+0.25*alphaL[1]*fUpL[8]+0.25*fUpL[1]*alphaL[8]+0.25*alphaL[0]*fUpL[4]+0.25*fUpL[0]*alphaL[4]; 
  GhatL[5] = 0.25*alphaL[10]*fUpL[15]+0.25*alphaL[13]*fUpL[14]+0.25*alphaL[4]*fUpL[12]+0.25*alphaL[3]*fUpL[11]+0.25*fUpL[3]*alphaL[11]+0.25*alphaL[8]*fUpL[9]+0.25*alphaL[6]*fUpL[7]+0.25*fUpL[6]*alphaL[7]+0.25*alphaL[0]*fUpL[5]+0.25*fUpL[0]*alphaL[5]+0.25*alphaL[1]*fUpL[2]+0.25*fUpL[1]*alphaL[2]; 
  GhatL[6] = 0.25*alphaL[4]*fUpL[13]+0.25*fUpL[4]*alphaL[13]+0.25*alphaL[2]*fUpL[11]+0.25*fUpL[2]*alphaL[11]+0.25*alphaL[8]*fUpL[10]+0.25*fUpL[8]*alphaL[10]+0.25*alphaL[5]*fUpL[7]+0.25*fUpL[5]*alphaL[7]+0.25*alphaL[0]*fUpL[6]+0.25*fUpL[0]*alphaL[6]+0.25*alphaL[1]*fUpL[3]+0.25*fUpL[1]*alphaL[3]; 
  GhatL[7] = 0.25*alphaL[8]*fUpL[15]+0.25*alphaL[4]*fUpL[14]+0.25*fUpL[12]*alphaL[13]+0.25*alphaL[1]*fUpL[11]+0.25*fUpL[1]*alphaL[11]+0.25*fUpL[9]*alphaL[10]+0.25*alphaL[0]*fUpL[7]+0.25*fUpL[0]*alphaL[7]+0.25*alphaL[5]*fUpL[6]+0.25*fUpL[5]*alphaL[6]+0.25*alphaL[2]*fUpL[3]+0.25*fUpL[2]*alphaL[3]; 
  GhatL[8] = 0.25*alphaL[7]*fUpL[15]+0.25*alphaL[11]*fUpL[14]+0.25*alphaL[3]*fUpL[13]+0.25*fUpL[3]*alphaL[13]+0.25*alphaL[2]*fUpL[12]+0.25*alphaL[6]*fUpL[10]+0.25*fUpL[6]*alphaL[10]+0.25*alphaL[5]*fUpL[9]+0.25*alphaL[0]*fUpL[8]+0.25*fUpL[0]*alphaL[8]+0.25*alphaL[1]*fUpL[4]+0.25*fUpL[1]*alphaL[4]; 
  GhatL[9] = 0.25*alphaL[6]*fUpL[15]+0.25*alphaL[3]*fUpL[14]+0.25*alphaL[11]*fUpL[13]+0.25*fUpL[11]*alphaL[13]+0.25*alphaL[1]*fUpL[12]+0.25*alphaL[7]*fUpL[10]+0.25*fUpL[7]*alphaL[10]+0.25*alphaL[0]*fUpL[9]+0.25*alphaL[5]*fUpL[8]+0.25*fUpL[5]*alphaL[8]+0.25*alphaL[2]*fUpL[4]+0.25*fUpL[2]*alphaL[4]; 
  GhatL[10] = 0.25*alphaL[5]*fUpL[15]+0.25*alphaL[2]*fUpL[14]+0.25*alphaL[1]*fUpL[13]+0.25*fUpL[1]*alphaL[13]+0.25*alphaL[11]*fUpL[12]+0.25*alphaL[0]*fUpL[10]+0.25*fUpL[0]*alphaL[10]+0.25*alphaL[7]*fUpL[9]+0.25*alphaL[6]*fUpL[8]+0.25*fUpL[6]*alphaL[8]+0.25*alphaL[3]*fUpL[4]+0.25*fUpL[3]*alphaL[4]; 
  GhatL[11] = 0.25*alphaL[4]*fUpL[15]+0.25*alphaL[8]*fUpL[14]+0.25*fUpL[9]*alphaL[13]+0.25*alphaL[10]*fUpL[12]+0.25*alphaL[0]*fUpL[11]+0.25*fUpL[0]*alphaL[11]+0.25*alphaL[1]*fUpL[7]+0.25*fUpL[1]*alphaL[7]+0.25*alphaL[2]*fUpL[6]+0.25*fUpL[2]*alphaL[6]+0.25*alphaL[3]*fUpL[5]+0.25*fUpL[3]*alphaL[5]; 
  GhatL[12] = 0.25*alphaL[3]*fUpL[15]+0.25*alphaL[6]*fUpL[14]+0.25*alphaL[7]*fUpL[13]+0.25*fUpL[7]*alphaL[13]+0.25*alphaL[0]*fUpL[12]+0.25*alphaL[10]*fUpL[11]+0.25*fUpL[10]*alphaL[11]+0.25*alphaL[1]*fUpL[9]+0.25*alphaL[2]*fUpL[8]+0.25*fUpL[2]*alphaL[8]+0.25*alphaL[4]*fUpL[5]+0.25*fUpL[4]*alphaL[5]; 
  GhatL[13] = 0.25*alphaL[2]*fUpL[15]+0.25*alphaL[5]*fUpL[14]+0.25*alphaL[0]*fUpL[13]+0.25*fUpL[0]*alphaL[13]+0.25*alphaL[7]*fUpL[12]+0.25*fUpL[9]*alphaL[11]+0.25*alphaL[1]*fUpL[10]+0.25*fUpL[1]*alphaL[10]+0.25*alphaL[3]*fUpL[8]+0.25*fUpL[3]*alphaL[8]+0.25*alphaL[4]*fUpL[6]+0.25*fUpL[4]*alphaL[6]; 
  GhatL[14] = 0.25*alphaL[1]*fUpL[15]+0.25*alphaL[0]*fUpL[14]+0.25*alphaL[5]*fUpL[13]+0.25*fUpL[5]*alphaL[13]+0.25*alphaL[6]*fUpL[12]+0.25*alphaL[8]*fUpL[11]+0.25*fUpL[8]*alphaL[11]+0.25*alphaL[2]*fUpL[10]+0.25*fUpL[2]*alphaL[10]+0.25*alphaL[3]*fUpL[9]+0.25*alphaL[4]*fUpL[7]+0.25*fUpL[4]*alphaL[7]; 
  GhatL[15] = 0.25*alphaL[0]*fUpL[15]+0.25*alphaL[1]*fUpL[14]+0.25*alphaL[2]*fUpL[13]+0.25*fUpL[2]*alphaL[13]+0.25*alphaL[3]*fUpL[12]+0.25*alphaL[4]*fUpL[11]+0.25*fUpL[4]*alphaL[11]+0.25*alphaL[5]*fUpL[10]+0.25*fUpL[5]*alphaL[10]+0.25*alphaL[6]*fUpL[9]+0.25*alphaL[7]*fUpL[8]+0.25*fUpL[7]*alphaL[8]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdvpar2; 
  out[1] += 0.7071067811865475*GhatL[1]*rdvpar2; 
  out[2] += 0.7071067811865475*GhatL[2]*rdvpar2; 
  out[3] += 0.7071067811865475*GhatL[3]*rdvpar2; 
  out[4] += -1.224744871391589*GhatL[0]*rdvpar2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdvpar2; 
  out[6] += 0.7071067811865475*GhatL[5]*rdvpar2; 
  out[7] += 0.7071067811865475*GhatL[6]*rdvpar2; 
  out[8] += 0.7071067811865475*GhatL[7]*rdvpar2; 
  out[9] += -1.224744871391589*GhatL[1]*rdvpar2; 
  out[10] += -1.224744871391589*GhatL[2]*rdvpar2; 
  out[11] += -1.224744871391589*GhatL[3]*rdvpar2; 
  out[12] += 0.7071067811865475*GhatL[8]*rdvpar2; 
  out[13] += 0.7071067811865475*GhatL[9]*rdvpar2; 
  out[14] += 0.7071067811865475*GhatL[10]*rdvpar2; 
  out[15] += -1.224744871391589*GhatL[4]*rdvpar2; 
  out[16] += 0.7071067811865475*GhatL[11]*rdvpar2; 
  out[17] += -1.224744871391589*GhatL[5]*rdvpar2; 
  out[18] += -1.224744871391589*GhatL[6]*rdvpar2; 
  out[19] += -1.224744871391589*GhatL[7]*rdvpar2; 
  out[20] += 0.7071067811865475*GhatL[12]*rdvpar2; 
  out[21] += 0.7071067811865475*GhatL[13]*rdvpar2; 
  out[22] += 0.7071067811865475*GhatL[14]*rdvpar2; 
  out[23] += -1.224744871391589*GhatL[8]*rdvpar2; 
  out[24] += -1.224744871391589*GhatL[9]*rdvpar2; 
  out[25] += -1.224744871391589*GhatL[10]*rdvpar2; 
  out[26] += -1.224744871391589*GhatL[11]*rdvpar2; 
  out[27] += 0.7071067811865475*GhatL[15]*rdvpar2; 
  out[28] += -1.224744871391589*GhatL[12]*rdvpar2; 
  out[29] += -1.224744871391589*GhatL[13]*rdvpar2; 
  out[30] += -1.224744871391589*GhatL[14]*rdvpar2; 
  out[31] += -1.224744871391589*GhatL[15]*rdvpar2; 
  out[32] += 1.58113883008419*GhatL[0]*rdvpar2; 
  out[33] += 1.58113883008419*GhatL[1]*rdvpar2; 
  out[34] += 1.58113883008419*GhatL[2]*rdvpar2; 
  out[35] += 1.58113883008419*GhatL[3]*rdvpar2; 
  out[36] += 1.58113883008419*GhatL[4]*rdvpar2; 
  out[37] += 1.58113883008419*GhatL[5]*rdvpar2; 
  out[38] += 1.58113883008419*GhatL[6]*rdvpar2; 
  out[39] += 1.58113883008419*GhatL[7]*rdvpar2; 
  out[40] += 1.58113883008419*GhatL[8]*rdvpar2; 
  out[41] += 1.58113883008419*GhatL[9]*rdvpar2; 
  out[42] += 1.58113883008419*GhatL[10]*rdvpar2; 
  out[43] += 1.58113883008419*GhatL[11]*rdvpar2; 
  out[44] += 1.58113883008419*GhatL[12]*rdvpar2; 
  out[45] += 1.58113883008419*GhatL[13]*rdvpar2; 
  out[46] += 1.58113883008419*GhatL[14]*rdvpar2; 
  out[47] += 1.58113883008419*GhatL[15]*rdvpar2; 

  } 

  return 5.0*rdvpar2*cflFreq; 

} 
