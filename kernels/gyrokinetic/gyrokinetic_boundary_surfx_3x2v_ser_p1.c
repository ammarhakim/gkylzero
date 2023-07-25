#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_surfx1_eval_quad.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double gyrokinetic_boundary_surfx_3x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
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

  double cflFreq = 0.0;

  if (edge == -1) { 

  double alphaR[24] = {0.}; 
  alphaR[0] = (0.0625*(m_*((((12.72792206135786*b_y[5]+7.348469228349534*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(7.348469228349534*b_y[5]+4.242640687119286*b_y[3])+(12.72792206135786*b_y[1]+7.348469228349534*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(7.348469228349534*b_y[1]+4.242640687119286*b_y[0]))*hamil[7]+hamil[3]*((7.348469228349534*b_y[5]+4.242640687119286*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(4.242640687119286*b_y[5]+2.449489742783178*b_y[3])+(7.348469228349534*b_y[1]+4.242640687119286*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_y[1]+2.449489742783178*b_y[0])))*rdz2+((((-12.72792206135786*b_z[1])-7.348469228349534*b_z[0])*jacobtot_inv[5]+((-12.72792206135786*jacobtot_inv[1])-7.348469228349534*jacobtot_inv[0])*b_z[5]+((-7.348469228349534*b_z[1])-4.242640687119286*b_z[0])*jacobtot_inv[3]+((-7.348469228349534*jacobtot_inv[1])-4.242640687119286*jacobtot_inv[0])*b_z[3])*hamil[16]+(((-7.348469228349534*b_z[1])-4.242640687119286*b_z[0])*jacobtot_inv[5]+((-7.348469228349534*jacobtot_inv[1])-4.242640687119286*jacobtot_inv[0])*b_z[5]+((-4.242640687119286*b_z[1])-2.449489742783178*b_z[0])*jacobtot_inv[3]+((-4.242640687119286*jacobtot_inv[1])-2.449489742783178*jacobtot_inv[0])*b_z[3])*hamil[8]+(((-12.72792206135786*b_z[5])-7.348469228349534*b_z[3])*jacobtot_inv[5]+jacobtot_inv[3]*((-7.348469228349534*b_z[5])-4.242640687119286*b_z[3])+((-12.72792206135786*b_z[1])-7.348469228349534*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*((-7.348469228349534*b_z[1])-4.242640687119286*b_z[0]))*hamil[6]+hamil[2]*(((-7.348469228349534*b_z[5])-4.242640687119286*b_z[3])*jacobtot_inv[5]+jacobtot_inv[3]*((-4.242640687119286*b_z[5])-2.449489742783178*b_z[3])+((-7.348469228349534*b_z[1])-4.242640687119286*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*((-4.242640687119286*b_z[1])-2.449489742783178*b_z[0])))*rdy2)+((13.41640786499874*BstarXdBmag[9]+7.745966692414834*BstarXdBmag[4])*hamil[32]+(6.0*BstarXdBmag[1]+3.464101615137754*BstarXdBmag[0])*hamil[4])*q_*rdvpar2))/(m_*q_); 
  alphaR[1] = (0.0625*(((12.72792206135786*b_y[5]+7.348469228349534*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(7.348469228349534*b_y[5]+4.242640687119286*b_y[3])+(12.72792206135786*b_y[1]+7.348469228349534*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(7.348469228349534*b_y[1]+4.242640687119286*b_y[0]))*hamil[16]+((7.348469228349534*b_y[5]+4.242640687119286*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(4.242640687119286*b_y[5]+2.449489742783178*b_y[3])+(7.348469228349534*b_y[1]+4.242640687119286*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_y[1]+2.449489742783178*b_y[0]))*hamil[8])*rdz2)/q_; 
  alphaR[2] = (0.0125*(m_*((((63.63961030678928*b_y[1]+36.74234614174767*b_y[0])*jacobtot_inv[5]+(63.63961030678928*jacobtot_inv[1]+36.74234614174767*jacobtot_inv[0])*b_y[5]+(36.74234614174767*b_y[1]+21.21320343559643*b_y[0])*jacobtot_inv[3]+(36.74234614174767*jacobtot_inv[1]+21.21320343559643*jacobtot_inv[0])*b_y[3])*hamil[7]+hamil[3]*((36.74234614174767*b_y[1]+21.21320343559643*b_y[0])*jacobtot_inv[5]+(36.74234614174767*jacobtot_inv[1]+21.21320343559643*jacobtot_inv[0])*b_y[5]+(21.21320343559643*b_y[1]+12.24744871391589*b_y[0])*jacobtot_inv[3]+(21.21320343559643*jacobtot_inv[1]+12.24744871391589*jacobtot_inv[0])*b_y[3]))*rdz2+((((-114.5512985522207*b_z[5])-66.13622305514579*b_z[3])*jacobtot_inv[5]+jacobtot_inv[3]*((-66.13622305514579*b_z[5])-38.18376618407357*b_z[3])+((-63.63961030678928*b_z[1])-36.74234614174767*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*((-36.74234614174767*b_z[1])-21.21320343559643*b_z[0]))*hamil[16]+(((-66.13622305514579*b_z[5])-38.18376618407357*b_z[3])*jacobtot_inv[5]+jacobtot_inv[3]*((-38.18376618407357*b_z[5])-22.0454076850486*b_z[3])+((-36.74234614174767*b_z[1])-21.21320343559643*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*((-21.21320343559643*b_z[1])-12.24744871391589*b_z[0]))*hamil[8]+(((-63.63961030678928*b_z[1])-36.74234614174767*b_z[0])*jacobtot_inv[5]+((-63.63961030678928*jacobtot_inv[1])-36.74234614174767*jacobtot_inv[0])*b_z[5]+((-36.74234614174767*b_z[1])-21.21320343559643*b_z[0])*jacobtot_inv[3]+((-36.74234614174767*jacobtot_inv[1])-21.21320343559643*jacobtot_inv[0])*b_z[3])*hamil[6]+hamil[2]*(((-36.74234614174767*b_z[1])-21.21320343559643*b_z[0])*jacobtot_inv[5]+((-36.74234614174767*jacobtot_inv[1])-21.21320343559643*jacobtot_inv[0])*b_z[5]+((-21.21320343559643*b_z[1])-12.24744871391589*b_z[0])*jacobtot_inv[3]+((-21.21320343559643*jacobtot_inv[1])-12.24744871391589*jacobtot_inv[0])*b_z[3]))*rdy2)+((67.0820393249937*BstarXdBmag[18]+38.72983346207418*BstarXdBmag[11])*hamil[32]+hamil[4]*(30.0*BstarXdBmag[7]+17.32050807568877*BstarXdBmag[3]))*q_*rdvpar2))/(m_*q_); 
  alphaR[3] = (0.125*((6.708203932499369*BstarXdBmag[1]+3.872983346207417*BstarXdBmag[0])*hamil[32]+hamil[4]*(3.0*BstarXdBmag[9]+1.732050807568877*BstarXdBmag[4]))*rdvpar2)/m_; 
  alphaR[4] = (0.0625*(((12.72792206135786*b_y[5]+7.348469228349534*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(7.348469228349534*b_y[5]+4.242640687119286*b_y[3])+(12.72792206135786*b_y[1]+7.348469228349534*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(7.348469228349534*b_y[1]+4.242640687119286*b_y[0]))*hamil[21]+((7.348469228349534*b_y[5]+4.242640687119286*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(4.242640687119286*b_y[5]+2.449489742783178*b_y[3])+(7.348469228349534*b_y[1]+4.242640687119286*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_y[1]+2.449489742783178*b_y[0]))*hamil[14])*rdz2)/q_; 
  alphaR[5] = (0.0625*(((12.72792206135786*b_y[1]+7.348469228349534*b_y[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[1]+7.348469228349534*jacobtot_inv[0])*b_y[5]+(7.348469228349534*b_y[1]+4.242640687119286*b_y[0])*jacobtot_inv[3]+(7.348469228349534*jacobtot_inv[1]+4.242640687119286*jacobtot_inv[0])*b_y[3])*hamil[16]+((7.348469228349534*b_y[1]+4.242640687119286*b_y[0])*jacobtot_inv[5]+(7.348469228349534*jacobtot_inv[1]+4.242640687119286*jacobtot_inv[0])*b_y[5]+(4.242640687119286*b_y[1]+2.449489742783178*b_y[0])*jacobtot_inv[3]+(4.242640687119286*jacobtot_inv[1]+2.449489742783178*jacobtot_inv[0])*b_y[3])*hamil[8])*rdz2)/q_; 
  alphaR[7] = (0.125*((6.708203932499369*BstarXdBmag[7]+3.872983346207417*BstarXdBmag[3])*hamil[32]+hamil[4]*(3.0*BstarXdBmag[18]+1.732050807568877*BstarXdBmag[11]))*rdvpar2)/m_; 
  alphaR[9] = (0.0625*(((12.72792206135786*b_y[1]+7.348469228349534*b_y[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[1]+7.348469228349534*jacobtot_inv[0])*b_y[5]+(7.348469228349534*b_y[1]+4.242640687119286*b_y[0])*jacobtot_inv[3]+(7.348469228349534*jacobtot_inv[1]+4.242640687119286*jacobtot_inv[0])*b_y[3])*hamil[21]+((7.348469228349534*b_y[1]+4.242640687119286*b_y[0])*jacobtot_inv[5]+(7.348469228349534*jacobtot_inv[1]+4.242640687119286*jacobtot_inv[0])*b_y[5]+(4.242640687119286*b_y[1]+2.449489742783178*b_y[0])*jacobtot_inv[3]+(4.242640687119286*jacobtot_inv[1]+2.449489742783178*jacobtot_inv[0])*b_y[3])*hamil[14])*rdz2)/q_; 
  alphaR[16] = (0.25*(3.0*BstarXdBmag[9]+1.732050807568877*BstarXdBmag[4])*hamil[32]*rdvpar2)/m_; 
  alphaR[18] = (0.05*(15.0*BstarXdBmag[18]+8.660254037844387*BstarXdBmag[11])*hamil[32]*rdvpar2)/m_; 

  double fUpOrdR[24] = {0.};
  double alphaR_n = 0.;

  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]+0.25*alphaR[9]+0.3354101966249678*alphaR[7]+0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx1_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx1_eval_quad_node_0_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[16]+0.25*alphaR[9]+0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx1_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx1_eval_quad_node_1_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]+0.25*alphaR[9]-0.3354101966249678*alphaR[7]+0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx1_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx1_eval_quad_node_2_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]-0.25*alphaR[9]+0.3354101966249678*alphaR[7]+0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx1_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx1_eval_quad_node_3_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[16]-0.25*alphaR[9]+0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx1_eval_quad_node_4_r(fskin); 
  } else { 
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx1_eval_quad_node_4_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]-0.25*alphaR[9]-0.3354101966249678*alphaR[7]+0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx1_eval_quad_node_5_r(fskin); 
  } else { 
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx1_eval_quad_node_5_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]-0.25*alphaR[9]-0.3354101966249678*alphaR[7]-0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx1_eval_quad_node_6_r(fskin); 
  } else { 
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx1_eval_quad_node_6_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[18])-0.2795084971874732*alphaR[16]-0.25*alphaR[9]-0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx1_eval_quad_node_7_r(fskin); 
  } else { 
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx1_eval_quad_node_7_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]-0.25*alphaR[9]+0.3354101966249678*alphaR[7]-0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx1_eval_quad_node_8_r(fskin); 
  } else { 
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx1_eval_quad_node_8_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]+0.25*alphaR[9]-0.3354101966249678*alphaR[7]-0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx1_eval_quad_node_9_r(fskin); 
  } else { 
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx1_eval_quad_node_9_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[18])-0.2795084971874732*alphaR[16]+0.25*alphaR[9]-0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx1_eval_quad_node_10_r(fskin); 
  } else { 
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx1_eval_quad_node_10_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]+0.25*alphaR[9]+0.3354101966249678*alphaR[7]-0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]-0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx1_eval_quad_node_11_r(fskin); 
  } else { 
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx1_eval_quad_node_11_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]+0.25*alphaR[9]+0.3354101966249678*alphaR[7]-0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx1_eval_quad_node_12_r(fskin); 
  } else { 
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx1_eval_quad_node_12_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[16]+0.25*alphaR[9]-0.25*alphaR[5]-0.25*alphaR[4]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx1_eval_quad_node_13_r(fskin); 
  } else { 
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx1_eval_quad_node_13_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]+0.25*alphaR[9]-0.3354101966249678*alphaR[7]-0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx1_eval_quad_node_14_r(fskin); 
  } else { 
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx1_eval_quad_node_14_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]-0.25*alphaR[9]+0.3354101966249678*alphaR[7]-0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx1_eval_quad_node_15_r(fskin); 
  } else { 
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx1_eval_quad_node_15_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2795084971874732*alphaR[18]-0.2795084971874732*alphaR[16]-0.25*alphaR[9]-0.25*alphaR[5]+0.25*alphaR[4]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[16] = gkhyb_3x2v_p1_surfx1_eval_quad_node_16_r(fskin); 
  } else { 
    fUpOrdR[16] = gkhyb_3x2v_p1_surfx1_eval_quad_node_16_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2236067977499786*alphaR[18])+0.2236067977499786*alphaR[16]-0.25*alphaR[9]-0.3354101966249678*alphaR[7]-0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]-0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[17] = gkhyb_3x2v_p1_surfx1_eval_quad_node_17_r(fskin); 
  } else { 
    fUpOrdR[17] = gkhyb_3x2v_p1_surfx1_eval_quad_node_17_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]-0.25*alphaR[9]-0.3354101966249678*alphaR[7]+0.25*alphaR[5]-0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[18] = gkhyb_3x2v_p1_surfx1_eval_quad_node_18_r(fskin); 
  } else { 
    fUpOrdR[18] = gkhyb_3x2v_p1_surfx1_eval_quad_node_18_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[18])-0.2795084971874732*alphaR[16]-0.25*alphaR[9]+0.25*alphaR[5]-0.25*alphaR[4]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[19] = gkhyb_3x2v_p1_surfx1_eval_quad_node_19_r(fskin); 
  } else { 
    fUpOrdR[19] = gkhyb_3x2v_p1_surfx1_eval_quad_node_19_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]-0.25*alphaR[9]+0.3354101966249678*alphaR[7]+0.25*alphaR[5]-0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[20] = gkhyb_3x2v_p1_surfx1_eval_quad_node_20_r(fskin); 
  } else { 
    fUpOrdR[20] = gkhyb_3x2v_p1_surfx1_eval_quad_node_20_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]+0.25*alphaR[9]-0.3354101966249678*alphaR[7]+0.25*alphaR[5]+0.25*alphaR[4]-0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[21] = gkhyb_3x2v_p1_surfx1_eval_quad_node_21_r(fskin); 
  } else { 
    fUpOrdR[21] = gkhyb_3x2v_p1_surfx1_eval_quad_node_21_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = (-0.2795084971874732*alphaR[18])-0.2795084971874732*alphaR[16]+0.25*alphaR[9]+0.25*alphaR[5]+0.25*alphaR[4]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[22] = gkhyb_3x2v_p1_surfx1_eval_quad_node_22_r(fskin); 
  } else { 
    fUpOrdR[22] = gkhyb_3x2v_p1_surfx1_eval_quad_node_22_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 
  alphaR_n = 0.2236067977499786*alphaR[18]+0.2236067977499786*alphaR[16]+0.25*alphaR[9]+0.3354101966249678*alphaR[7]+0.25*alphaR[5]+0.25*alphaR[4]+0.3354101966249678*alphaR[3]+0.25*alphaR[2]+0.25*alphaR[1]+0.25*alphaR[0];
  if (alphaR_n > 0.) {
    fUpOrdR[23] = gkhyb_3x2v_p1_surfx1_eval_quad_node_23_r(fskin); 
  } else { 
    fUpOrdR[23] = gkhyb_3x2v_p1_surfx1_eval_quad_node_23_l(fedge); 
  } 
  cflFreq += -0.046875*rdx2*(alphaR_n-fabs(alphaR_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[24] = {0.};
  gkhyb_3x2v_p1_xdir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[24] = {0.}; 
  GhatR[0] = 0.25*alphaR[18]*fUpR[18]+0.25*alphaR[16]*fUpR[16]+0.25*alphaR[9]*fUpR[9]+0.25*alphaR[7]*fUpR[7]+0.25*alphaR[5]*fUpR[5]+0.25*alphaR[4]*fUpR[4]+0.25*alphaR[3]*fUpR[3]+0.25*alphaR[2]*fUpR[2]+0.25*alphaR[1]*fUpR[1]+0.25*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.2500000000000001*alphaR[18]*fUpR[20]+0.2500000000000001*alphaR[16]*fUpR[17]+0.25*alphaR[9]*fUpR[12]+0.25*alphaR[7]*fUpR[11]+0.25*alphaR[4]*fUpR[8]+0.25*alphaR[3]*fUpR[6]+0.25*alphaR[2]*fUpR[5]+0.25*fUpR[2]*alphaR[5]+0.25*alphaR[0]*fUpR[1]+0.25*fUpR[0]*alphaR[1]; 
  GhatR[2] = 0.2500000000000001*alphaR[16]*fUpR[18]+0.2500000000000001*fUpR[16]*alphaR[18]+0.25*alphaR[4]*fUpR[9]+0.25*fUpR[4]*alphaR[9]+0.25*alphaR[3]*fUpR[7]+0.25*fUpR[3]*alphaR[7]+0.25*alphaR[1]*fUpR[5]+0.25*fUpR[1]*alphaR[5]+0.25*alphaR[0]*fUpR[2]+0.25*fUpR[0]*alphaR[2]; 
  GhatR[3] = 0.223606797749979*alphaR[7]*fUpR[18]+0.223606797749979*fUpR[7]*alphaR[18]+0.223606797749979*alphaR[3]*fUpR[16]+0.223606797749979*fUpR[3]*alphaR[16]+0.25*alphaR[9]*fUpR[14]+0.25*alphaR[5]*fUpR[11]+0.25*alphaR[4]*fUpR[10]+0.25*alphaR[2]*fUpR[7]+0.25*fUpR[2]*alphaR[7]+0.25*alphaR[1]*fUpR[6]+0.25*alphaR[0]*fUpR[3]+0.25*fUpR[0]*alphaR[3]; 
  GhatR[4] = 0.2500000000000001*alphaR[18]*fUpR[22]+0.2500000000000001*alphaR[16]*fUpR[19]+0.25*alphaR[7]*fUpR[14]+0.25*alphaR[5]*fUpR[12]+0.25*alphaR[3]*fUpR[10]+0.25*alphaR[2]*fUpR[9]+0.25*fUpR[2]*alphaR[9]+0.25*alphaR[1]*fUpR[8]+0.25*alphaR[0]*fUpR[4]+0.25*fUpR[0]*alphaR[4]; 
  GhatR[5] = 0.25*alphaR[16]*fUpR[20]+0.25*fUpR[17]*alphaR[18]+0.25*alphaR[4]*fUpR[12]+0.25*alphaR[3]*fUpR[11]+0.25*fUpR[8]*alphaR[9]+0.25*fUpR[6]*alphaR[7]+0.25*alphaR[0]*fUpR[5]+0.25*fUpR[0]*alphaR[5]+0.25*alphaR[1]*fUpR[2]+0.25*fUpR[1]*alphaR[2]; 
  GhatR[6] = 0.223606797749979*alphaR[7]*fUpR[20]+0.223606797749979*fUpR[11]*alphaR[18]+0.223606797749979*alphaR[3]*fUpR[17]+0.223606797749979*fUpR[6]*alphaR[16]+0.25*alphaR[9]*fUpR[15]+0.25*alphaR[4]*fUpR[13]+0.25*alphaR[2]*fUpR[11]+0.25*alphaR[5]*fUpR[7]+0.25*fUpR[5]*alphaR[7]+0.25*alphaR[0]*fUpR[6]+0.25*alphaR[1]*fUpR[3]+0.25*fUpR[1]*alphaR[3]; 
  GhatR[7] = 0.223606797749979*alphaR[3]*fUpR[18]+0.223606797749979*fUpR[3]*alphaR[18]+0.223606797749979*alphaR[7]*fUpR[16]+0.223606797749979*fUpR[7]*alphaR[16]+0.25*alphaR[4]*fUpR[14]+0.25*alphaR[1]*fUpR[11]+0.25*alphaR[9]*fUpR[10]+0.25*alphaR[0]*fUpR[7]+0.25*fUpR[0]*alphaR[7]+0.25*alphaR[5]*fUpR[6]+0.25*alphaR[2]*fUpR[3]+0.25*fUpR[2]*alphaR[3]; 
  GhatR[8] = 0.25*alphaR[18]*fUpR[23]+0.25*alphaR[16]*fUpR[21]+0.25*alphaR[7]*fUpR[15]+0.25*alphaR[3]*fUpR[13]+0.25*alphaR[2]*fUpR[12]+0.25*alphaR[5]*fUpR[9]+0.25*fUpR[5]*alphaR[9]+0.25*alphaR[0]*fUpR[8]+0.25*alphaR[1]*fUpR[4]+0.25*fUpR[1]*alphaR[4]; 
  GhatR[9] = 0.25*alphaR[16]*fUpR[22]+0.25*alphaR[18]*fUpR[19]+0.25*alphaR[3]*fUpR[14]+0.25*alphaR[1]*fUpR[12]+0.25*alphaR[7]*fUpR[10]+0.25*alphaR[0]*fUpR[9]+0.25*fUpR[0]*alphaR[9]+0.25*alphaR[5]*fUpR[8]+0.25*alphaR[2]*fUpR[4]+0.25*fUpR[2]*alphaR[4]; 
  GhatR[10] = 0.223606797749979*alphaR[7]*fUpR[22]+0.223606797749979*alphaR[3]*fUpR[19]+0.223606797749979*fUpR[14]*alphaR[18]+0.223606797749979*fUpR[10]*alphaR[16]+0.25*alphaR[5]*fUpR[15]+0.25*alphaR[2]*fUpR[14]+0.25*alphaR[1]*fUpR[13]+0.25*alphaR[0]*fUpR[10]+0.25*alphaR[7]*fUpR[9]+0.25*fUpR[7]*alphaR[9]+0.25*alphaR[3]*fUpR[4]+0.25*fUpR[3]*alphaR[4]; 
  GhatR[11] = 0.223606797749979*alphaR[3]*fUpR[20]+0.223606797749979*fUpR[6]*alphaR[18]+0.223606797749979*alphaR[7]*fUpR[17]+0.223606797749979*fUpR[11]*alphaR[16]+0.25*alphaR[4]*fUpR[15]+0.25*alphaR[9]*fUpR[13]+0.25*alphaR[0]*fUpR[11]+0.25*alphaR[1]*fUpR[7]+0.25*fUpR[1]*alphaR[7]+0.25*alphaR[2]*fUpR[6]+0.25*alphaR[3]*fUpR[5]+0.25*fUpR[3]*alphaR[5]; 
  GhatR[12] = 0.2500000000000001*alphaR[16]*fUpR[23]+0.2500000000000001*alphaR[18]*fUpR[21]+0.25*alphaR[3]*fUpR[15]+0.25*alphaR[7]*fUpR[13]+0.25*alphaR[0]*fUpR[12]+0.25*alphaR[1]*fUpR[9]+0.25*fUpR[1]*alphaR[9]+0.25*alphaR[2]*fUpR[8]+0.25*alphaR[4]*fUpR[5]+0.25*fUpR[4]*alphaR[5]; 
  GhatR[13] = 0.223606797749979*alphaR[7]*fUpR[23]+0.223606797749979*alphaR[3]*fUpR[21]+0.223606797749979*fUpR[15]*alphaR[18]+0.223606797749979*fUpR[13]*alphaR[16]+0.25*alphaR[2]*fUpR[15]+0.25*alphaR[5]*fUpR[14]+0.25*alphaR[0]*fUpR[13]+0.25*alphaR[7]*fUpR[12]+0.25*alphaR[9]*fUpR[11]+0.25*alphaR[1]*fUpR[10]+0.25*alphaR[3]*fUpR[8]+0.25*alphaR[4]*fUpR[6]; 
  GhatR[14] = 0.223606797749979*alphaR[3]*fUpR[22]+0.223606797749979*alphaR[7]*fUpR[19]+0.223606797749979*fUpR[10]*alphaR[18]+0.223606797749979*fUpR[14]*alphaR[16]+0.25*alphaR[1]*fUpR[15]+0.25*alphaR[0]*fUpR[14]+0.25*alphaR[5]*fUpR[13]+0.25*alphaR[2]*fUpR[10]+0.25*alphaR[3]*fUpR[9]+0.25*fUpR[3]*alphaR[9]+0.25*alphaR[4]*fUpR[7]+0.25*fUpR[4]*alphaR[7]; 
  GhatR[15] = 0.223606797749979*alphaR[3]*fUpR[23]+0.223606797749979*alphaR[7]*fUpR[21]+0.223606797749979*fUpR[13]*alphaR[18]+0.223606797749979*fUpR[15]*alphaR[16]+0.25*alphaR[0]*fUpR[15]+0.25*alphaR[1]*fUpR[14]+0.25*alphaR[2]*fUpR[13]+0.25*alphaR[3]*fUpR[12]+0.25*alphaR[4]*fUpR[11]+0.25*alphaR[5]*fUpR[10]+0.25*fUpR[6]*alphaR[9]+0.25*alphaR[7]*fUpR[8]; 
  GhatR[16] = 0.25*alphaR[9]*fUpR[22]+0.25*alphaR[5]*fUpR[20]+0.2500000000000001*alphaR[4]*fUpR[19]+0.159719141249985*alphaR[18]*fUpR[18]+0.2500000000000001*alphaR[2]*fUpR[18]+0.2500000000000001*fUpR[2]*alphaR[18]+0.2500000000000001*alphaR[1]*fUpR[17]+0.159719141249985*alphaR[16]*fUpR[16]+0.25*alphaR[0]*fUpR[16]+0.25*fUpR[0]*alphaR[16]+0.223606797749979*alphaR[7]*fUpR[7]+0.223606797749979*alphaR[3]*fUpR[3]; 
  GhatR[17] = 0.25*alphaR[9]*fUpR[23]+0.2500000000000001*alphaR[4]*fUpR[21]+0.159719141249985*alphaR[18]*fUpR[20]+0.2500000000000001*alphaR[2]*fUpR[20]+0.25*alphaR[5]*fUpR[18]+0.25*fUpR[5]*alphaR[18]+0.159719141249985*alphaR[16]*fUpR[17]+0.25*alphaR[0]*fUpR[17]+0.2500000000000001*alphaR[1]*fUpR[16]+0.2500000000000001*fUpR[1]*alphaR[16]+0.223606797749979*alphaR[7]*fUpR[11]+0.223606797749979*alphaR[3]*fUpR[6]; 
  GhatR[18] = 0.2500000000000001*alphaR[4]*fUpR[22]+0.2500000000000001*alphaR[1]*fUpR[20]+0.25*alphaR[9]*fUpR[19]+0.159719141249985*alphaR[16]*fUpR[18]+0.25*alphaR[0]*fUpR[18]+0.159719141249985*fUpR[16]*alphaR[18]+0.25*fUpR[0]*alphaR[18]+0.25*alphaR[5]*fUpR[17]+0.2500000000000001*alphaR[2]*fUpR[16]+0.2500000000000001*fUpR[2]*alphaR[16]+0.223606797749979*alphaR[3]*fUpR[7]+0.223606797749979*fUpR[3]*alphaR[7]; 
  GhatR[19] = 0.25*alphaR[5]*fUpR[23]+0.159719141249985*alphaR[18]*fUpR[22]+0.2500000000000001*alphaR[2]*fUpR[22]+0.2500000000000001*alphaR[1]*fUpR[21]+0.159719141249985*alphaR[16]*fUpR[19]+0.25*alphaR[0]*fUpR[19]+0.25*alphaR[9]*fUpR[18]+0.25*fUpR[9]*alphaR[18]+0.2500000000000001*alphaR[4]*fUpR[16]+0.2500000000000001*fUpR[4]*alphaR[16]+0.223606797749979*alphaR[7]*fUpR[14]+0.223606797749979*alphaR[3]*fUpR[10]; 
  GhatR[20] = 0.2500000000000001*alphaR[4]*fUpR[23]+0.25*alphaR[9]*fUpR[21]+0.159719141249985*alphaR[16]*fUpR[20]+0.25*alphaR[0]*fUpR[20]+0.2500000000000001*alphaR[1]*fUpR[18]+0.159719141249985*fUpR[17]*alphaR[18]+0.2500000000000001*fUpR[1]*alphaR[18]+0.2500000000000001*alphaR[2]*fUpR[17]+0.25*alphaR[5]*fUpR[16]+0.25*fUpR[5]*alphaR[16]+0.223606797749979*alphaR[3]*fUpR[11]+0.223606797749979*fUpR[6]*alphaR[7]; 
  GhatR[21] = 0.159719141249985*alphaR[18]*fUpR[23]+0.2500000000000001*alphaR[2]*fUpR[23]+0.25*alphaR[5]*fUpR[22]+0.159719141249985*alphaR[16]*fUpR[21]+0.25*alphaR[0]*fUpR[21]+0.25*alphaR[9]*fUpR[20]+0.2500000000000001*alphaR[1]*fUpR[19]+0.2500000000000001*fUpR[12]*alphaR[18]+0.2500000000000001*alphaR[4]*fUpR[17]+0.25*fUpR[8]*alphaR[16]+0.223606797749979*alphaR[7]*fUpR[15]+0.223606797749979*alphaR[3]*fUpR[13]; 
  GhatR[22] = 0.2500000000000001*alphaR[1]*fUpR[23]+0.159719141249985*alphaR[16]*fUpR[22]+0.25*alphaR[0]*fUpR[22]+0.25*alphaR[5]*fUpR[21]+0.159719141249985*alphaR[18]*fUpR[19]+0.2500000000000001*alphaR[2]*fUpR[19]+0.2500000000000001*alphaR[4]*fUpR[18]+0.2500000000000001*fUpR[4]*alphaR[18]+0.25*alphaR[9]*fUpR[16]+0.25*fUpR[9]*alphaR[16]+0.223606797749979*alphaR[3]*fUpR[14]+0.223606797749979*alphaR[7]*fUpR[10]; 
  GhatR[23] = 0.159719141249985*alphaR[16]*fUpR[23]+0.25*alphaR[0]*fUpR[23]+0.2500000000000001*alphaR[1]*fUpR[22]+0.159719141249985*alphaR[18]*fUpR[21]+0.2500000000000001*alphaR[2]*fUpR[21]+0.2500000000000001*alphaR[4]*fUpR[20]+0.25*alphaR[5]*fUpR[19]+0.25*fUpR[8]*alphaR[18]+0.25*alphaR[9]*fUpR[17]+0.2500000000000001*fUpR[12]*alphaR[16]+0.223606797749979*alphaR[3]*fUpR[15]+0.223606797749979*alphaR[7]*fUpR[13]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdx2; 
  out[1] += -1.224744871391589*GhatR[0]*rdx2; 
  out[2] += -0.7071067811865475*GhatR[1]*rdx2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdx2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdx2; 
  out[5] += -0.7071067811865475*GhatR[4]*rdx2; 
  out[6] += -1.224744871391589*GhatR[1]*rdx2; 
  out[7] += -1.224744871391589*GhatR[2]*rdx2; 
  out[8] += -0.7071067811865475*GhatR[5]*rdx2; 
  out[9] += -1.224744871391589*GhatR[3]*rdx2; 
  out[10] += -0.7071067811865475*GhatR[6]*rdx2; 
  out[11] += -0.7071067811865475*GhatR[7]*rdx2; 
  out[12] += -1.224744871391589*GhatR[4]*rdx2; 
  out[13] += -0.7071067811865475*GhatR[8]*rdx2; 
  out[14] += -0.7071067811865475*GhatR[9]*rdx2; 
  out[15] += -0.7071067811865475*GhatR[10]*rdx2; 
  out[16] += -1.224744871391589*GhatR[5]*rdx2; 
  out[17] += -1.224744871391589*GhatR[6]*rdx2; 
  out[18] += -1.224744871391589*GhatR[7]*rdx2; 
  out[19] += -0.7071067811865475*GhatR[11]*rdx2; 
  out[20] += -1.224744871391589*GhatR[8]*rdx2; 
  out[21] += -1.224744871391589*GhatR[9]*rdx2; 
  out[22] += -0.7071067811865475*GhatR[12]*rdx2; 
  out[23] += -1.224744871391589*GhatR[10]*rdx2; 
  out[24] += -0.7071067811865475*GhatR[13]*rdx2; 
  out[25] += -0.7071067811865475*GhatR[14]*rdx2; 
  out[26] += -1.224744871391589*GhatR[11]*rdx2; 
  out[27] += -1.224744871391589*GhatR[12]*rdx2; 
  out[28] += -1.224744871391589*GhatR[13]*rdx2; 
  out[29] += -1.224744871391589*GhatR[14]*rdx2; 
  out[30] += -0.7071067811865475*GhatR[15]*rdx2; 
  out[31] += -1.224744871391589*GhatR[15]*rdx2; 
  out[32] += -0.7071067811865475*GhatR[16]*rdx2; 
  out[33] += -1.224744871391589*GhatR[16]*rdx2; 
  out[34] += -0.7071067811865475*GhatR[17]*rdx2; 
  out[35] += -0.7071067811865475*GhatR[18]*rdx2; 
  out[36] += -0.7071067811865475*GhatR[19]*rdx2; 
  out[37] += -1.224744871391589*GhatR[17]*rdx2; 
  out[38] += -1.224744871391589*GhatR[18]*rdx2; 
  out[39] += -0.7071067811865475*GhatR[20]*rdx2; 
  out[40] += -1.224744871391589*GhatR[19]*rdx2; 
  out[41] += -0.7071067811865475*GhatR[21]*rdx2; 
  out[42] += -0.7071067811865475*GhatR[22]*rdx2; 
  out[43] += -1.224744871391589*GhatR[20]*rdx2; 
  out[44] += -1.224744871391589*GhatR[21]*rdx2; 
  out[45] += -1.224744871391589*GhatR[22]*rdx2; 
  out[46] += -0.7071067811865475*GhatR[23]*rdx2; 
  out[47] += -1.224744871391589*GhatR[23]*rdx2; 

  } else { 

  double alphaL[24] = {0.}; 
  alphaL[0] = -(0.0625*(m_*((((12.72792206135786*b_y[5]-7.348469228349534*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(4.242640687119286*b_y[3]-7.348469228349534*b_y[5])+(12.72792206135786*b_y[1]-7.348469228349534*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_y[0]-7.348469228349534*b_y[1]))*hamil[7]+hamil[3]*((4.242640687119286*b_y[3]-7.348469228349534*b_y[5])*jacobtot_inv[5]+jacobtot_inv[3]*(4.242640687119286*b_y[5]-2.449489742783178*b_y[3])+(4.242640687119286*b_y[0]-7.348469228349534*b_y[1])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_y[1]-2.449489742783178*b_y[0])))*rdz2+(((7.348469228349534*b_z[0]-12.72792206135786*b_z[1])*jacobtot_inv[5]+(7.348469228349534*jacobtot_inv[0]-12.72792206135786*jacobtot_inv[1])*b_z[5]+(7.348469228349534*b_z[1]-4.242640687119286*b_z[0])*jacobtot_inv[3]+(7.348469228349534*jacobtot_inv[1]-4.242640687119286*jacobtot_inv[0])*b_z[3])*hamil[16]+((7.348469228349534*b_z[1]-4.242640687119286*b_z[0])*jacobtot_inv[5]+(7.348469228349534*jacobtot_inv[1]-4.242640687119286*jacobtot_inv[0])*b_z[5]+(2.449489742783178*b_z[0]-4.242640687119286*b_z[1])*jacobtot_inv[3]+(2.449489742783178*jacobtot_inv[0]-4.242640687119286*jacobtot_inv[1])*b_z[3])*hamil[8]+((7.348469228349534*b_z[3]-12.72792206135786*b_z[5])*jacobtot_inv[5]+jacobtot_inv[3]*(7.348469228349534*b_z[5]-4.242640687119286*b_z[3])+(7.348469228349534*b_z[0]-12.72792206135786*b_z[1])*jacobtot_inv[1]+jacobtot_inv[0]*(7.348469228349534*b_z[1]-4.242640687119286*b_z[0]))*hamil[6]+hamil[2]*((7.348469228349534*b_z[5]-4.242640687119286*b_z[3])*jacobtot_inv[5]+jacobtot_inv[3]*(2.449489742783178*b_z[3]-4.242640687119286*b_z[5])+(7.348469228349534*b_z[1]-4.242640687119286*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(2.449489742783178*b_z[0]-4.242640687119286*b_z[1])))*rdy2)+((13.41640786499874*BstarXdBmag[9]-7.745966692414834*BstarXdBmag[4])*hamil[32]+(6.0*BstarXdBmag[1]-3.464101615137754*BstarXdBmag[0])*hamil[4])*q_*rdvpar2))/(m_*q_); 
  alphaL[1] = -(0.0625*(((12.72792206135786*b_y[5]-7.348469228349534*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(4.242640687119286*b_y[3]-7.348469228349534*b_y[5])+(12.72792206135786*b_y[1]-7.348469228349534*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_y[0]-7.348469228349534*b_y[1]))*hamil[16]+((4.242640687119286*b_y[3]-7.348469228349534*b_y[5])*jacobtot_inv[5]+jacobtot_inv[3]*(4.242640687119286*b_y[5]-2.449489742783178*b_y[3])+(4.242640687119286*b_y[0]-7.348469228349534*b_y[1])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_y[1]-2.449489742783178*b_y[0]))*hamil[8])*rdz2)/q_; 
  alphaL[2] = -(0.0125*(m_*((((63.63961030678928*b_y[1]-36.74234614174767*b_y[0])*jacobtot_inv[5]+(63.63961030678928*jacobtot_inv[1]-36.74234614174767*jacobtot_inv[0])*b_y[5]+(21.21320343559643*b_y[0]-36.74234614174767*b_y[1])*jacobtot_inv[3]+(21.21320343559643*jacobtot_inv[0]-36.74234614174767*jacobtot_inv[1])*b_y[3])*hamil[7]+hamil[3]*((21.21320343559643*b_y[0]-36.74234614174767*b_y[1])*jacobtot_inv[5]+(21.21320343559643*jacobtot_inv[0]-36.74234614174767*jacobtot_inv[1])*b_y[5]+(21.21320343559643*b_y[1]-12.24744871391589*b_y[0])*jacobtot_inv[3]+(21.21320343559643*jacobtot_inv[1]-12.24744871391589*jacobtot_inv[0])*b_y[3]))*rdz2+(((66.13622305514579*b_z[3]-114.5512985522207*b_z[5])*jacobtot_inv[5]+jacobtot_inv[3]*(66.13622305514579*b_z[5]-38.18376618407357*b_z[3])+(36.74234614174767*b_z[0]-63.63961030678928*b_z[1])*jacobtot_inv[1]+jacobtot_inv[0]*(36.74234614174767*b_z[1]-21.21320343559643*b_z[0]))*hamil[16]+((66.13622305514579*b_z[5]-38.18376618407357*b_z[3])*jacobtot_inv[5]+jacobtot_inv[3]*(22.0454076850486*b_z[3]-38.18376618407357*b_z[5])+(36.74234614174767*b_z[1]-21.21320343559643*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(12.24744871391589*b_z[0]-21.21320343559643*b_z[1]))*hamil[8]+((36.74234614174767*b_z[0]-63.63961030678928*b_z[1])*jacobtot_inv[5]+(36.74234614174767*jacobtot_inv[0]-63.63961030678928*jacobtot_inv[1])*b_z[5]+(36.74234614174767*b_z[1]-21.21320343559643*b_z[0])*jacobtot_inv[3]+(36.74234614174767*jacobtot_inv[1]-21.21320343559643*jacobtot_inv[0])*b_z[3])*hamil[6]+hamil[2]*((36.74234614174767*b_z[1]-21.21320343559643*b_z[0])*jacobtot_inv[5]+(36.74234614174767*jacobtot_inv[1]-21.21320343559643*jacobtot_inv[0])*b_z[5]+(12.24744871391589*b_z[0]-21.21320343559643*b_z[1])*jacobtot_inv[3]+(12.24744871391589*jacobtot_inv[0]-21.21320343559643*jacobtot_inv[1])*b_z[3]))*rdy2)+((67.0820393249937*BstarXdBmag[18]-38.72983346207418*BstarXdBmag[11])*hamil[32]+hamil[4]*(30.0*BstarXdBmag[7]-17.32050807568877*BstarXdBmag[3]))*q_*rdvpar2))/(m_*q_); 
  alphaL[3] = -(0.125*((6.708203932499369*BstarXdBmag[1]-3.872983346207417*BstarXdBmag[0])*hamil[32]+hamil[4]*(3.0*BstarXdBmag[9]-1.732050807568877*BstarXdBmag[4]))*rdvpar2)/m_; 
  alphaL[4] = -(0.0625*(((12.72792206135786*b_y[5]-7.348469228349534*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(4.242640687119286*b_y[3]-7.348469228349534*b_y[5])+(12.72792206135786*b_y[1]-7.348469228349534*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_y[0]-7.348469228349534*b_y[1]))*hamil[21]+((4.242640687119286*b_y[3]-7.348469228349534*b_y[5])*jacobtot_inv[5]+jacobtot_inv[3]*(4.242640687119286*b_y[5]-2.449489742783178*b_y[3])+(4.242640687119286*b_y[0]-7.348469228349534*b_y[1])*jacobtot_inv[1]+jacobtot_inv[0]*(4.242640687119286*b_y[1]-2.449489742783178*b_y[0]))*hamil[14])*rdz2)/q_; 
  alphaL[5] = -(0.0625*(((12.72792206135786*b_y[1]-7.348469228349534*b_y[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[1]-7.348469228349534*jacobtot_inv[0])*b_y[5]+(4.242640687119286*b_y[0]-7.348469228349534*b_y[1])*jacobtot_inv[3]+(4.242640687119286*jacobtot_inv[0]-7.348469228349534*jacobtot_inv[1])*b_y[3])*hamil[16]+((4.242640687119286*b_y[0]-7.348469228349534*b_y[1])*jacobtot_inv[5]+(4.242640687119286*jacobtot_inv[0]-7.348469228349534*jacobtot_inv[1])*b_y[5]+(4.242640687119286*b_y[1]-2.449489742783178*b_y[0])*jacobtot_inv[3]+(4.242640687119286*jacobtot_inv[1]-2.449489742783178*jacobtot_inv[0])*b_y[3])*hamil[8])*rdz2)/q_; 
  alphaL[7] = -(0.125*((6.708203932499369*BstarXdBmag[7]-3.872983346207417*BstarXdBmag[3])*hamil[32]+hamil[4]*(3.0*BstarXdBmag[18]-1.732050807568877*BstarXdBmag[11]))*rdvpar2)/m_; 
  alphaL[9] = -(0.0625*(((12.72792206135786*b_y[1]-7.348469228349534*b_y[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[1]-7.348469228349534*jacobtot_inv[0])*b_y[5]+(4.242640687119286*b_y[0]-7.348469228349534*b_y[1])*jacobtot_inv[3]+(4.242640687119286*jacobtot_inv[0]-7.348469228349534*jacobtot_inv[1])*b_y[3])*hamil[21]+((4.242640687119286*b_y[0]-7.348469228349534*b_y[1])*jacobtot_inv[5]+(4.242640687119286*jacobtot_inv[0]-7.348469228349534*jacobtot_inv[1])*b_y[5]+(4.242640687119286*b_y[1]-2.449489742783178*b_y[0])*jacobtot_inv[3]+(4.242640687119286*jacobtot_inv[1]-2.449489742783178*jacobtot_inv[0])*b_y[3])*hamil[14])*rdz2)/q_; 
  alphaL[16] = -(0.25*(3.0*BstarXdBmag[9]-1.732050807568877*BstarXdBmag[4])*hamil[32]*rdvpar2)/m_; 
  alphaL[18] = -(0.05*(15.0*BstarXdBmag[18]-8.660254037844387*BstarXdBmag[11])*hamil[32]*rdvpar2)/m_; 

  double fUpOrdL[24] = {0.};
  double alphaL_n = 0.;

  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]+0.25*alphaL[9]+0.3354101966249678*alphaL[7]+0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx1_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx1_eval_quad_node_0_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]+0.25*alphaL[9]+0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx1_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx1_eval_quad_node_1_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]+0.25*alphaL[9]-0.3354101966249678*alphaL[7]+0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx1_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx1_eval_quad_node_2_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]-0.25*alphaL[9]+0.3354101966249678*alphaL[7]+0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx1_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx1_eval_quad_node_3_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.25*alphaL[9]+0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx1_eval_quad_node_4_r(fedge); 
  } else { 
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx1_eval_quad_node_4_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]-0.25*alphaL[9]-0.3354101966249678*alphaL[7]+0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx1_eval_quad_node_5_r(fedge); 
  } else { 
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx1_eval_quad_node_5_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.25*alphaL[9]-0.3354101966249678*alphaL[7]-0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx1_eval_quad_node_6_r(fedge); 
  } else { 
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx1_eval_quad_node_6_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[18])-0.2795084971874732*alphaL[16]-0.25*alphaL[9]-0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx1_eval_quad_node_7_r(fedge); 
  } else { 
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx1_eval_quad_node_7_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.25*alphaL[9]+0.3354101966249678*alphaL[7]-0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx1_eval_quad_node_8_r(fedge); 
  } else { 
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx1_eval_quad_node_8_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.25*alphaL[9]-0.3354101966249678*alphaL[7]-0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx1_eval_quad_node_9_r(fedge); 
  } else { 
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx1_eval_quad_node_9_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[18])-0.2795084971874732*alphaL[16]+0.25*alphaL[9]-0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx1_eval_quad_node_10_r(fedge); 
  } else { 
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx1_eval_quad_node_10_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.25*alphaL[9]+0.3354101966249678*alphaL[7]-0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx1_eval_quad_node_11_r(fedge); 
  } else { 
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx1_eval_quad_node_11_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]+0.25*alphaL[9]+0.3354101966249678*alphaL[7]-0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx1_eval_quad_node_12_r(fedge); 
  } else { 
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx1_eval_quad_node_12_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]+0.25*alphaL[9]-0.25*alphaL[5]-0.25*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx1_eval_quad_node_13_r(fedge); 
  } else { 
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx1_eval_quad_node_13_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]+0.25*alphaL[9]-0.3354101966249678*alphaL[7]-0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx1_eval_quad_node_14_r(fedge); 
  } else { 
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx1_eval_quad_node_14_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]-0.25*alphaL[9]+0.3354101966249678*alphaL[7]-0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx1_eval_quad_node_15_r(fedge); 
  } else { 
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx1_eval_quad_node_15_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.25*alphaL[9]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[16] = gkhyb_3x2v_p1_surfx1_eval_quad_node_16_r(fedge); 
  } else { 
    fUpOrdL[16] = gkhyb_3x2v_p1_surfx1_eval_quad_node_16_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2236067977499786*alphaL[18])+0.2236067977499786*alphaL[16]-0.25*alphaL[9]-0.3354101966249678*alphaL[7]-0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]-0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[17] = gkhyb_3x2v_p1_surfx1_eval_quad_node_17_r(fedge); 
  } else { 
    fUpOrdL[17] = gkhyb_3x2v_p1_surfx1_eval_quad_node_17_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.25*alphaL[9]-0.3354101966249678*alphaL[7]+0.25*alphaL[5]-0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[18] = gkhyb_3x2v_p1_surfx1_eval_quad_node_18_r(fedge); 
  } else { 
    fUpOrdL[18] = gkhyb_3x2v_p1_surfx1_eval_quad_node_18_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[18])-0.2795084971874732*alphaL[16]-0.25*alphaL[9]+0.25*alphaL[5]-0.25*alphaL[4]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[19] = gkhyb_3x2v_p1_surfx1_eval_quad_node_19_r(fedge); 
  } else { 
    fUpOrdL[19] = gkhyb_3x2v_p1_surfx1_eval_quad_node_19_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]-0.25*alphaL[9]+0.3354101966249678*alphaL[7]+0.25*alphaL[5]-0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[20] = gkhyb_3x2v_p1_surfx1_eval_quad_node_20_r(fedge); 
  } else { 
    fUpOrdL[20] = gkhyb_3x2v_p1_surfx1_eval_quad_node_20_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.25*alphaL[9]-0.3354101966249678*alphaL[7]+0.25*alphaL[5]+0.25*alphaL[4]-0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[21] = gkhyb_3x2v_p1_surfx1_eval_quad_node_21_r(fedge); 
  } else { 
    fUpOrdL[21] = gkhyb_3x2v_p1_surfx1_eval_quad_node_21_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = (-0.2795084971874732*alphaL[18])-0.2795084971874732*alphaL[16]+0.25*alphaL[9]+0.25*alphaL[5]+0.25*alphaL[4]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[22] = gkhyb_3x2v_p1_surfx1_eval_quad_node_22_r(fedge); 
  } else { 
    fUpOrdL[22] = gkhyb_3x2v_p1_surfx1_eval_quad_node_22_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 
  alphaL_n = 0.2236067977499786*alphaL[18]+0.2236067977499786*alphaL[16]+0.25*alphaL[9]+0.3354101966249678*alphaL[7]+0.25*alphaL[5]+0.25*alphaL[4]+0.3354101966249678*alphaL[3]+0.25*alphaL[2]+0.25*alphaL[1]+0.25*alphaL[0];
  if (alphaL_n > 0.) {
    fUpOrdL[23] = gkhyb_3x2v_p1_surfx1_eval_quad_node_23_r(fedge); 
  } else { 
    fUpOrdL[23] = gkhyb_3x2v_p1_surfx1_eval_quad_node_23_l(fskin); 
  } 
  cflFreq += -0.046875*rdx2*(alphaL_n-fabs(alphaL_n)); 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[24] = {0.};
  gkhyb_3x2v_p1_xdir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[24] = {0.}; 
  GhatL[0] = 0.25*alphaL[18]*fUpL[18]+0.25*alphaL[16]*fUpL[16]+0.25*alphaL[9]*fUpL[9]+0.25*alphaL[7]*fUpL[7]+0.25*alphaL[5]*fUpL[5]+0.25*alphaL[4]*fUpL[4]+0.25*alphaL[3]*fUpL[3]+0.25*alphaL[2]*fUpL[2]+0.25*alphaL[1]*fUpL[1]+0.25*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.2500000000000001*alphaL[18]*fUpL[20]+0.2500000000000001*alphaL[16]*fUpL[17]+0.25*alphaL[9]*fUpL[12]+0.25*alphaL[7]*fUpL[11]+0.25*alphaL[4]*fUpL[8]+0.25*alphaL[3]*fUpL[6]+0.25*alphaL[2]*fUpL[5]+0.25*fUpL[2]*alphaL[5]+0.25*alphaL[0]*fUpL[1]+0.25*fUpL[0]*alphaL[1]; 
  GhatL[2] = 0.2500000000000001*alphaL[16]*fUpL[18]+0.2500000000000001*fUpL[16]*alphaL[18]+0.25*alphaL[4]*fUpL[9]+0.25*fUpL[4]*alphaL[9]+0.25*alphaL[3]*fUpL[7]+0.25*fUpL[3]*alphaL[7]+0.25*alphaL[1]*fUpL[5]+0.25*fUpL[1]*alphaL[5]+0.25*alphaL[0]*fUpL[2]+0.25*fUpL[0]*alphaL[2]; 
  GhatL[3] = 0.223606797749979*alphaL[7]*fUpL[18]+0.223606797749979*fUpL[7]*alphaL[18]+0.223606797749979*alphaL[3]*fUpL[16]+0.223606797749979*fUpL[3]*alphaL[16]+0.25*alphaL[9]*fUpL[14]+0.25*alphaL[5]*fUpL[11]+0.25*alphaL[4]*fUpL[10]+0.25*alphaL[2]*fUpL[7]+0.25*fUpL[2]*alphaL[7]+0.25*alphaL[1]*fUpL[6]+0.25*alphaL[0]*fUpL[3]+0.25*fUpL[0]*alphaL[3]; 
  GhatL[4] = 0.2500000000000001*alphaL[18]*fUpL[22]+0.2500000000000001*alphaL[16]*fUpL[19]+0.25*alphaL[7]*fUpL[14]+0.25*alphaL[5]*fUpL[12]+0.25*alphaL[3]*fUpL[10]+0.25*alphaL[2]*fUpL[9]+0.25*fUpL[2]*alphaL[9]+0.25*alphaL[1]*fUpL[8]+0.25*alphaL[0]*fUpL[4]+0.25*fUpL[0]*alphaL[4]; 
  GhatL[5] = 0.25*alphaL[16]*fUpL[20]+0.25*fUpL[17]*alphaL[18]+0.25*alphaL[4]*fUpL[12]+0.25*alphaL[3]*fUpL[11]+0.25*fUpL[8]*alphaL[9]+0.25*fUpL[6]*alphaL[7]+0.25*alphaL[0]*fUpL[5]+0.25*fUpL[0]*alphaL[5]+0.25*alphaL[1]*fUpL[2]+0.25*fUpL[1]*alphaL[2]; 
  GhatL[6] = 0.223606797749979*alphaL[7]*fUpL[20]+0.223606797749979*fUpL[11]*alphaL[18]+0.223606797749979*alphaL[3]*fUpL[17]+0.223606797749979*fUpL[6]*alphaL[16]+0.25*alphaL[9]*fUpL[15]+0.25*alphaL[4]*fUpL[13]+0.25*alphaL[2]*fUpL[11]+0.25*alphaL[5]*fUpL[7]+0.25*fUpL[5]*alphaL[7]+0.25*alphaL[0]*fUpL[6]+0.25*alphaL[1]*fUpL[3]+0.25*fUpL[1]*alphaL[3]; 
  GhatL[7] = 0.223606797749979*alphaL[3]*fUpL[18]+0.223606797749979*fUpL[3]*alphaL[18]+0.223606797749979*alphaL[7]*fUpL[16]+0.223606797749979*fUpL[7]*alphaL[16]+0.25*alphaL[4]*fUpL[14]+0.25*alphaL[1]*fUpL[11]+0.25*alphaL[9]*fUpL[10]+0.25*alphaL[0]*fUpL[7]+0.25*fUpL[0]*alphaL[7]+0.25*alphaL[5]*fUpL[6]+0.25*alphaL[2]*fUpL[3]+0.25*fUpL[2]*alphaL[3]; 
  GhatL[8] = 0.25*alphaL[18]*fUpL[23]+0.25*alphaL[16]*fUpL[21]+0.25*alphaL[7]*fUpL[15]+0.25*alphaL[3]*fUpL[13]+0.25*alphaL[2]*fUpL[12]+0.25*alphaL[5]*fUpL[9]+0.25*fUpL[5]*alphaL[9]+0.25*alphaL[0]*fUpL[8]+0.25*alphaL[1]*fUpL[4]+0.25*fUpL[1]*alphaL[4]; 
  GhatL[9] = 0.25*alphaL[16]*fUpL[22]+0.25*alphaL[18]*fUpL[19]+0.25*alphaL[3]*fUpL[14]+0.25*alphaL[1]*fUpL[12]+0.25*alphaL[7]*fUpL[10]+0.25*alphaL[0]*fUpL[9]+0.25*fUpL[0]*alphaL[9]+0.25*alphaL[5]*fUpL[8]+0.25*alphaL[2]*fUpL[4]+0.25*fUpL[2]*alphaL[4]; 
  GhatL[10] = 0.223606797749979*alphaL[7]*fUpL[22]+0.223606797749979*alphaL[3]*fUpL[19]+0.223606797749979*fUpL[14]*alphaL[18]+0.223606797749979*fUpL[10]*alphaL[16]+0.25*alphaL[5]*fUpL[15]+0.25*alphaL[2]*fUpL[14]+0.25*alphaL[1]*fUpL[13]+0.25*alphaL[0]*fUpL[10]+0.25*alphaL[7]*fUpL[9]+0.25*fUpL[7]*alphaL[9]+0.25*alphaL[3]*fUpL[4]+0.25*fUpL[3]*alphaL[4]; 
  GhatL[11] = 0.223606797749979*alphaL[3]*fUpL[20]+0.223606797749979*fUpL[6]*alphaL[18]+0.223606797749979*alphaL[7]*fUpL[17]+0.223606797749979*fUpL[11]*alphaL[16]+0.25*alphaL[4]*fUpL[15]+0.25*alphaL[9]*fUpL[13]+0.25*alphaL[0]*fUpL[11]+0.25*alphaL[1]*fUpL[7]+0.25*fUpL[1]*alphaL[7]+0.25*alphaL[2]*fUpL[6]+0.25*alphaL[3]*fUpL[5]+0.25*fUpL[3]*alphaL[5]; 
  GhatL[12] = 0.2500000000000001*alphaL[16]*fUpL[23]+0.2500000000000001*alphaL[18]*fUpL[21]+0.25*alphaL[3]*fUpL[15]+0.25*alphaL[7]*fUpL[13]+0.25*alphaL[0]*fUpL[12]+0.25*alphaL[1]*fUpL[9]+0.25*fUpL[1]*alphaL[9]+0.25*alphaL[2]*fUpL[8]+0.25*alphaL[4]*fUpL[5]+0.25*fUpL[4]*alphaL[5]; 
  GhatL[13] = 0.223606797749979*alphaL[7]*fUpL[23]+0.223606797749979*alphaL[3]*fUpL[21]+0.223606797749979*fUpL[15]*alphaL[18]+0.223606797749979*fUpL[13]*alphaL[16]+0.25*alphaL[2]*fUpL[15]+0.25*alphaL[5]*fUpL[14]+0.25*alphaL[0]*fUpL[13]+0.25*alphaL[7]*fUpL[12]+0.25*alphaL[9]*fUpL[11]+0.25*alphaL[1]*fUpL[10]+0.25*alphaL[3]*fUpL[8]+0.25*alphaL[4]*fUpL[6]; 
  GhatL[14] = 0.223606797749979*alphaL[3]*fUpL[22]+0.223606797749979*alphaL[7]*fUpL[19]+0.223606797749979*fUpL[10]*alphaL[18]+0.223606797749979*fUpL[14]*alphaL[16]+0.25*alphaL[1]*fUpL[15]+0.25*alphaL[0]*fUpL[14]+0.25*alphaL[5]*fUpL[13]+0.25*alphaL[2]*fUpL[10]+0.25*alphaL[3]*fUpL[9]+0.25*fUpL[3]*alphaL[9]+0.25*alphaL[4]*fUpL[7]+0.25*fUpL[4]*alphaL[7]; 
  GhatL[15] = 0.223606797749979*alphaL[3]*fUpL[23]+0.223606797749979*alphaL[7]*fUpL[21]+0.223606797749979*fUpL[13]*alphaL[18]+0.223606797749979*fUpL[15]*alphaL[16]+0.25*alphaL[0]*fUpL[15]+0.25*alphaL[1]*fUpL[14]+0.25*alphaL[2]*fUpL[13]+0.25*alphaL[3]*fUpL[12]+0.25*alphaL[4]*fUpL[11]+0.25*alphaL[5]*fUpL[10]+0.25*fUpL[6]*alphaL[9]+0.25*alphaL[7]*fUpL[8]; 
  GhatL[16] = 0.25*alphaL[9]*fUpL[22]+0.25*alphaL[5]*fUpL[20]+0.2500000000000001*alphaL[4]*fUpL[19]+0.159719141249985*alphaL[18]*fUpL[18]+0.2500000000000001*alphaL[2]*fUpL[18]+0.2500000000000001*fUpL[2]*alphaL[18]+0.2500000000000001*alphaL[1]*fUpL[17]+0.159719141249985*alphaL[16]*fUpL[16]+0.25*alphaL[0]*fUpL[16]+0.25*fUpL[0]*alphaL[16]+0.223606797749979*alphaL[7]*fUpL[7]+0.223606797749979*alphaL[3]*fUpL[3]; 
  GhatL[17] = 0.25*alphaL[9]*fUpL[23]+0.2500000000000001*alphaL[4]*fUpL[21]+0.159719141249985*alphaL[18]*fUpL[20]+0.2500000000000001*alphaL[2]*fUpL[20]+0.25*alphaL[5]*fUpL[18]+0.25*fUpL[5]*alphaL[18]+0.159719141249985*alphaL[16]*fUpL[17]+0.25*alphaL[0]*fUpL[17]+0.2500000000000001*alphaL[1]*fUpL[16]+0.2500000000000001*fUpL[1]*alphaL[16]+0.223606797749979*alphaL[7]*fUpL[11]+0.223606797749979*alphaL[3]*fUpL[6]; 
  GhatL[18] = 0.2500000000000001*alphaL[4]*fUpL[22]+0.2500000000000001*alphaL[1]*fUpL[20]+0.25*alphaL[9]*fUpL[19]+0.159719141249985*alphaL[16]*fUpL[18]+0.25*alphaL[0]*fUpL[18]+0.159719141249985*fUpL[16]*alphaL[18]+0.25*fUpL[0]*alphaL[18]+0.25*alphaL[5]*fUpL[17]+0.2500000000000001*alphaL[2]*fUpL[16]+0.2500000000000001*fUpL[2]*alphaL[16]+0.223606797749979*alphaL[3]*fUpL[7]+0.223606797749979*fUpL[3]*alphaL[7]; 
  GhatL[19] = 0.25*alphaL[5]*fUpL[23]+0.159719141249985*alphaL[18]*fUpL[22]+0.2500000000000001*alphaL[2]*fUpL[22]+0.2500000000000001*alphaL[1]*fUpL[21]+0.159719141249985*alphaL[16]*fUpL[19]+0.25*alphaL[0]*fUpL[19]+0.25*alphaL[9]*fUpL[18]+0.25*fUpL[9]*alphaL[18]+0.2500000000000001*alphaL[4]*fUpL[16]+0.2500000000000001*fUpL[4]*alphaL[16]+0.223606797749979*alphaL[7]*fUpL[14]+0.223606797749979*alphaL[3]*fUpL[10]; 
  GhatL[20] = 0.2500000000000001*alphaL[4]*fUpL[23]+0.25*alphaL[9]*fUpL[21]+0.159719141249985*alphaL[16]*fUpL[20]+0.25*alphaL[0]*fUpL[20]+0.2500000000000001*alphaL[1]*fUpL[18]+0.159719141249985*fUpL[17]*alphaL[18]+0.2500000000000001*fUpL[1]*alphaL[18]+0.2500000000000001*alphaL[2]*fUpL[17]+0.25*alphaL[5]*fUpL[16]+0.25*fUpL[5]*alphaL[16]+0.223606797749979*alphaL[3]*fUpL[11]+0.223606797749979*fUpL[6]*alphaL[7]; 
  GhatL[21] = 0.159719141249985*alphaL[18]*fUpL[23]+0.2500000000000001*alphaL[2]*fUpL[23]+0.25*alphaL[5]*fUpL[22]+0.159719141249985*alphaL[16]*fUpL[21]+0.25*alphaL[0]*fUpL[21]+0.25*alphaL[9]*fUpL[20]+0.2500000000000001*alphaL[1]*fUpL[19]+0.2500000000000001*fUpL[12]*alphaL[18]+0.2500000000000001*alphaL[4]*fUpL[17]+0.25*fUpL[8]*alphaL[16]+0.223606797749979*alphaL[7]*fUpL[15]+0.223606797749979*alphaL[3]*fUpL[13]; 
  GhatL[22] = 0.2500000000000001*alphaL[1]*fUpL[23]+0.159719141249985*alphaL[16]*fUpL[22]+0.25*alphaL[0]*fUpL[22]+0.25*alphaL[5]*fUpL[21]+0.159719141249985*alphaL[18]*fUpL[19]+0.2500000000000001*alphaL[2]*fUpL[19]+0.2500000000000001*alphaL[4]*fUpL[18]+0.2500000000000001*fUpL[4]*alphaL[18]+0.25*alphaL[9]*fUpL[16]+0.25*fUpL[9]*alphaL[16]+0.223606797749979*alphaL[3]*fUpL[14]+0.223606797749979*alphaL[7]*fUpL[10]; 
  GhatL[23] = 0.159719141249985*alphaL[16]*fUpL[23]+0.25*alphaL[0]*fUpL[23]+0.2500000000000001*alphaL[1]*fUpL[22]+0.159719141249985*alphaL[18]*fUpL[21]+0.2500000000000001*alphaL[2]*fUpL[21]+0.2500000000000001*alphaL[4]*fUpL[20]+0.25*alphaL[5]*fUpL[19]+0.25*fUpL[8]*alphaL[18]+0.25*alphaL[9]*fUpL[17]+0.2500000000000001*fUpL[12]*alphaL[16]+0.223606797749979*alphaL[3]*fUpL[15]+0.223606797749979*alphaL[7]*fUpL[13]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdx2; 
  out[1] += -1.224744871391589*GhatL[0]*rdx2; 
  out[2] += 0.7071067811865475*GhatL[1]*rdx2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdx2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdx2; 
  out[5] += 0.7071067811865475*GhatL[4]*rdx2; 
  out[6] += -1.224744871391589*GhatL[1]*rdx2; 
  out[7] += -1.224744871391589*GhatL[2]*rdx2; 
  out[8] += 0.7071067811865475*GhatL[5]*rdx2; 
  out[9] += -1.224744871391589*GhatL[3]*rdx2; 
  out[10] += 0.7071067811865475*GhatL[6]*rdx2; 
  out[11] += 0.7071067811865475*GhatL[7]*rdx2; 
  out[12] += -1.224744871391589*GhatL[4]*rdx2; 
  out[13] += 0.7071067811865475*GhatL[8]*rdx2; 
  out[14] += 0.7071067811865475*GhatL[9]*rdx2; 
  out[15] += 0.7071067811865475*GhatL[10]*rdx2; 
  out[16] += -1.224744871391589*GhatL[5]*rdx2; 
  out[17] += -1.224744871391589*GhatL[6]*rdx2; 
  out[18] += -1.224744871391589*GhatL[7]*rdx2; 
  out[19] += 0.7071067811865475*GhatL[11]*rdx2; 
  out[20] += -1.224744871391589*GhatL[8]*rdx2; 
  out[21] += -1.224744871391589*GhatL[9]*rdx2; 
  out[22] += 0.7071067811865475*GhatL[12]*rdx2; 
  out[23] += -1.224744871391589*GhatL[10]*rdx2; 
  out[24] += 0.7071067811865475*GhatL[13]*rdx2; 
  out[25] += 0.7071067811865475*GhatL[14]*rdx2; 
  out[26] += -1.224744871391589*GhatL[11]*rdx2; 
  out[27] += -1.224744871391589*GhatL[12]*rdx2; 
  out[28] += -1.224744871391589*GhatL[13]*rdx2; 
  out[29] += -1.224744871391589*GhatL[14]*rdx2; 
  out[30] += 0.7071067811865475*GhatL[15]*rdx2; 
  out[31] += -1.224744871391589*GhatL[15]*rdx2; 
  out[32] += 0.7071067811865475*GhatL[16]*rdx2; 
  out[33] += -1.224744871391589*GhatL[16]*rdx2; 
  out[34] += 0.7071067811865475*GhatL[17]*rdx2; 
  out[35] += 0.7071067811865475*GhatL[18]*rdx2; 
  out[36] += 0.7071067811865475*GhatL[19]*rdx2; 
  out[37] += -1.224744871391589*GhatL[17]*rdx2; 
  out[38] += -1.224744871391589*GhatL[18]*rdx2; 
  out[39] += 0.7071067811865475*GhatL[20]*rdx2; 
  out[40] += -1.224744871391589*GhatL[19]*rdx2; 
  out[41] += 0.7071067811865475*GhatL[21]*rdx2; 
  out[42] += 0.7071067811865475*GhatL[22]*rdx2; 
  out[43] += -1.224744871391589*GhatL[20]*rdx2; 
  out[44] += -1.224744871391589*GhatL[21]*rdx2; 
  out[45] += -1.224744871391589*GhatL[22]*rdx2; 
  out[46] += 0.7071067811865475*GhatL[23]*rdx2; 
  out[47] += -1.224744871391589*GhatL[23]*rdx2; 

  } 

  return cflFreq; 

} 
