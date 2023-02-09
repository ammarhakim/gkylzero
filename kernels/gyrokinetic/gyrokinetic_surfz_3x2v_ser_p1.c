#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_gkhyb_3x2v_p1_surfx3_eval_quad.h> 
#include <gkyl_basis_gkhyb_3x2v_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH void gyrokinetic_surfz_3x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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

  double BstarZdBmag[48] = {0.}; 
  BstarZdBmag[0] = (0.7071067811865475*(1.732050807568877*(jacobtot_inv[3]*b_y[5]+jacobtot_inv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[5]*jacobtot_inv[5]+cmag[3]*jacobtot_inv[3]+cmag[1]*jacobtot_inv[1]+cmag[0]*jacobtot_inv[0])*q_))/q_; 
  BstarZdBmag[1] = (0.7071067811865475*(1.732050807568877*(b_y[5]*jacobtot_inv[5]+b_y[1]*jacobtot_inv[1])*m_*rdx2*wvpar+(cmag[3]*jacobtot_inv[5]+jacobtot_inv[3]*cmag[5]+cmag[0]*jacobtot_inv[1]+jacobtot_inv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[3] = (0.7071067811865475*(1.732050807568877*(jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3])*m_*rdx2*wvpar+(cmag[1]*jacobtot_inv[5]+jacobtot_inv[1]*cmag[5]+cmag[0]*jacobtot_inv[3]+jacobtot_inv[0]*cmag[3])*q_))/q_; 
  BstarZdBmag[4] = (0.7071067811865475*(jacobtot_inv[3]*b_y[5]+jacobtot_inv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[7] = (0.7071067811865475*(1.732050807568877*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5])*m_*rdx2*wvpar+(cmag[0]*jacobtot_inv[5]+jacobtot_inv[0]*cmag[5]+cmag[1]*jacobtot_inv[3]+jacobtot_inv[1]*cmag[3])*q_))/q_; 
  BstarZdBmag[9] = (0.7071067811865475*(b_y[5]*jacobtot_inv[5]+b_y[1]*jacobtot_inv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[11] = (0.7071067811865475*(jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[18] = (0.7071067811865475*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5])*m_*rdx2)/(q_*rdvpar2); 

  double alphaL[24] = {0.}; 
  alphaL[0] = -(0.0625*(m_*((((12.72792206135786*b_x[3]-7.348469228349534*b_x[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[3]-7.348469228349534*jacobtot_inv[0])*b_x[5]-7.348469228349534*(b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+4.242640687119286*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1]))*hamil[16]+((12.72792206135786*b_x[5]-7.348469228349534*b_x[1])*jacobtot_inv[5]-7.348469228349534*jacobtot_inv[1]*b_x[5]+(12.72792206135786*b_x[3]-7.348469228349534*b_x[0])*jacobtot_inv[3]-7.348469228349534*jacobtot_inv[0]*b_x[3]+4.242640687119286*(b_x[1]*jacobtot_inv[1]+b_x[0]*jacobtot_inv[0]))*hamil[8]+((4.242640687119286*b_x[0]-7.348469228349534*b_x[3])*jacobtot_inv[5]+(4.242640687119286*jacobtot_inv[0]-7.348469228349534*jacobtot_inv[3])*b_x[5]+4.242640687119286*(b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])-2.449489742783178*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1]))*hamil[6]+hamil[2]*((4.242640687119286*b_x[1]-7.348469228349534*b_x[5])*jacobtot_inv[5]+4.242640687119286*jacobtot_inv[1]*b_x[5]+(4.242640687119286*b_x[0]-7.348469228349534*b_x[3])*jacobtot_inv[3]+4.242640687119286*jacobtot_inv[0]*b_x[3]-2.449489742783178*(b_x[1]*jacobtot_inv[1]+b_x[0]*jacobtot_inv[0])))*rdy2+(((7.348469228349534*b_y[1]-12.72792206135786*b_y[5])*jacobtot_inv[5]+7.348469228349534*jacobtot_inv[1]*b_y[5]+(7.348469228349534*b_y[0]-12.72792206135786*b_y[3])*jacobtot_inv[3]+7.348469228349534*jacobtot_inv[0]*b_y[3]-4.242640687119286*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[7]+hamil[1]*((7.348469228349534*b_y[5]-4.242640687119286*b_y[1])*jacobtot_inv[5]-4.242640687119286*jacobtot_inv[1]*b_y[5]+(7.348469228349534*b_y[3]-4.242640687119286*b_y[0])*jacobtot_inv[3]-4.242640687119286*jacobtot_inv[0]*b_y[3]+2.449489742783178*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0])))*rdx2)+((13.41640786499874*BstarZdBmag[11]-7.745966692414834*BstarZdBmag[4])*hamil[32]+(6.0*BstarZdBmag[3]-3.464101615137754*BstarZdBmag[0])*hamil[4])*q_*rdvpar2))/(m_*q_); 
  alphaL[1] = -(0.0125*(m_*((((114.5512985522207*b_x[5]-66.13622305514579*b_x[1])*jacobtot_inv[5]-66.13622305514579*jacobtot_inv[1]*b_x[5]+(63.63961030678928*b_x[3]-36.74234614174767*b_x[0])*jacobtot_inv[3]-36.74234614174767*jacobtot_inv[0]*b_x[3]+38.18376618407357*b_x[1]*jacobtot_inv[1]+21.21320343559643*b_x[0]*jacobtot_inv[0])*hamil[16]+((63.63961030678928*b_x[3]-36.74234614174767*b_x[0])*jacobtot_inv[5]+(63.63961030678928*jacobtot_inv[3]-36.74234614174767*jacobtot_inv[0])*b_x[5]-36.74234614174767*(b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+21.21320343559643*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1]))*hamil[8]+((38.18376618407357*b_x[1]-66.13622305514579*b_x[5])*jacobtot_inv[5]+38.18376618407357*jacobtot_inv[1]*b_x[5]+(21.21320343559643*b_x[0]-36.74234614174767*b_x[3])*jacobtot_inv[3]+21.21320343559643*jacobtot_inv[0]*b_x[3]-22.0454076850486*b_x[1]*jacobtot_inv[1]-12.24744871391589*b_x[0]*jacobtot_inv[0])*hamil[6]+hamil[2]*((21.21320343559643*b_x[0]-36.74234614174767*b_x[3])*jacobtot_inv[5]+(21.21320343559643*jacobtot_inv[0]-36.74234614174767*jacobtot_inv[3])*b_x[5]+21.21320343559643*(b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])-12.24744871391589*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1])))*rdy2+(((36.74234614174767*b_y[0]-63.63961030678928*b_y[3])*jacobtot_inv[5]+(36.74234614174767*jacobtot_inv[0]-63.63961030678928*jacobtot_inv[3])*b_y[5]+36.74234614174767*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])-21.21320343559643*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[7]+hamil[1]*((36.74234614174767*b_y[3]-21.21320343559643*b_y[0])*jacobtot_inv[5]+(36.74234614174767*jacobtot_inv[3]-21.21320343559643*jacobtot_inv[0])*b_y[5]-21.21320343559643*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])+12.24744871391589*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1])))*rdx2)+((67.0820393249937*BstarZdBmag[18]-38.72983346207418*BstarZdBmag[9])*hamil[32]+hamil[4]*(30.0*BstarZdBmag[7]-17.32050807568877*BstarZdBmag[1]))*q_*rdvpar2))/(m_*q_); 
  alphaL[2] = (0.0625*(((12.72792206135786*b_y[5]-7.348469228349534*b_y[1])*jacobtot_inv[5]-7.348469228349534*jacobtot_inv[1]*b_y[5]+(12.72792206135786*b_y[3]-7.348469228349534*b_y[0])*jacobtot_inv[3]-7.348469228349534*jacobtot_inv[0]*b_y[3]+4.242640687119286*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[16]+((4.242640687119286*b_y[1]-7.348469228349534*b_y[5])*jacobtot_inv[5]+4.242640687119286*jacobtot_inv[1]*b_y[5]+(4.242640687119286*b_y[0]-7.348469228349534*b_y[3])*jacobtot_inv[3]+4.242640687119286*jacobtot_inv[0]*b_y[3]-2.449489742783178*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[6])*rdx2)/q_; 
  alphaL[3] = -(0.125*((6.708203932499369*BstarZdBmag[3]-3.872983346207417*BstarZdBmag[0])*hamil[32]+hamil[4]*(3.0*BstarZdBmag[11]-1.732050807568877*BstarZdBmag[4]))*rdvpar2)/m_; 
  alphaL[4] = (0.0625*(((12.72792206135786*b_y[5]-7.348469228349534*b_y[1])*jacobtot_inv[5]-7.348469228349534*jacobtot_inv[1]*b_y[5]+(12.72792206135786*b_y[3]-7.348469228349534*b_y[0])*jacobtot_inv[3]-7.348469228349534*jacobtot_inv[0]*b_y[3]+4.242640687119286*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[21]+((4.242640687119286*b_y[1]-7.348469228349534*b_y[5])*jacobtot_inv[5]+4.242640687119286*jacobtot_inv[1]*b_y[5]+(4.242640687119286*b_y[0]-7.348469228349534*b_y[3])*jacobtot_inv[3]+4.242640687119286*jacobtot_inv[0]*b_y[3]-2.449489742783178*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[12])*rdx2)/q_; 
  alphaL[5] = (0.0625*(((12.72792206135786*b_y[3]-7.348469228349534*b_y[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[3]-7.348469228349534*jacobtot_inv[0])*b_y[5]-7.348469228349534*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])+4.242640687119286*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[16]+((4.242640687119286*b_y[0]-7.348469228349534*b_y[3])*jacobtot_inv[5]+(4.242640687119286*jacobtot_inv[0]-7.348469228349534*jacobtot_inv[3])*b_y[5]+4.242640687119286*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])-2.449489742783178*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[6])*rdx2)/q_; 
  alphaL[6] = -(0.125*((6.708203932499369*BstarZdBmag[7]-3.872983346207417*BstarZdBmag[1])*hamil[32]+hamil[4]*(3.0*BstarZdBmag[18]-1.732050807568877*BstarZdBmag[9]))*rdvpar2)/m_; 
  alphaL[8] = (0.0625*(((12.72792206135786*b_y[3]-7.348469228349534*b_y[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[3]-7.348469228349534*jacobtot_inv[0])*b_y[5]-7.348469228349534*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])+4.242640687119286*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[21]+((4.242640687119286*b_y[0]-7.348469228349534*b_y[3])*jacobtot_inv[5]+(4.242640687119286*jacobtot_inv[0]-7.348469228349534*jacobtot_inv[3])*b_y[5]+4.242640687119286*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])-2.449489742783178*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[12])*rdx2)/q_; 
  alphaL[16] = -(0.25*(3.0*BstarZdBmag[11]-1.732050807568877*BstarZdBmag[4])*hamil[32]*rdvpar2)/m_; 
  alphaL[17] = -(0.05*(15.0*BstarZdBmag[18]-8.660254037844387*BstarZdBmag[9])*hamil[32]*rdvpar2)/m_; 

  double alphaR[24] = {0.}; 
  alphaR[0] = (0.0625*(m_*((((12.72792206135786*b_x[3]+7.348469228349534*b_x[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[3]+7.348469228349534*jacobtot_inv[0])*b_x[5]+7.348469228349534*(b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+4.242640687119286*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1]))*hamil[16]+((12.72792206135786*b_x[5]+7.348469228349534*b_x[1])*jacobtot_inv[5]+7.348469228349534*jacobtot_inv[1]*b_x[5]+(12.72792206135786*b_x[3]+7.348469228349534*b_x[0])*jacobtot_inv[3]+7.348469228349534*jacobtot_inv[0]*b_x[3]+4.242640687119286*(b_x[1]*jacobtot_inv[1]+b_x[0]*jacobtot_inv[0]))*hamil[8]+((7.348469228349534*b_x[3]+4.242640687119286*b_x[0])*jacobtot_inv[5]+(7.348469228349534*jacobtot_inv[3]+4.242640687119286*jacobtot_inv[0])*b_x[5]+4.242640687119286*(b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+2.449489742783178*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1]))*hamil[6]+hamil[2]*((7.348469228349534*b_x[5]+4.242640687119286*b_x[1])*jacobtot_inv[5]+4.242640687119286*jacobtot_inv[1]*b_x[5]+(7.348469228349534*b_x[3]+4.242640687119286*b_x[0])*jacobtot_inv[3]+4.242640687119286*jacobtot_inv[0]*b_x[3]+2.449489742783178*(b_x[1]*jacobtot_inv[1]+b_x[0]*jacobtot_inv[0])))*rdy2+((((-12.72792206135786*b_y[5])-7.348469228349534*b_y[1])*jacobtot_inv[5]-7.348469228349534*jacobtot_inv[1]*b_y[5]+((-12.72792206135786*b_y[3])-7.348469228349534*b_y[0])*jacobtot_inv[3]-7.348469228349534*jacobtot_inv[0]*b_y[3]-4.242640687119286*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[7]+hamil[1]*(((-7.348469228349534*b_y[5])-4.242640687119286*b_y[1])*jacobtot_inv[5]-4.242640687119286*jacobtot_inv[1]*b_y[5]+((-7.348469228349534*b_y[3])-4.242640687119286*b_y[0])*jacobtot_inv[3]-4.242640687119286*jacobtot_inv[0]*b_y[3]-2.449489742783178*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0])))*rdx2)+((13.41640786499874*BstarZdBmag[11]+7.745966692414834*BstarZdBmag[4])*hamil[32]+(6.0*BstarZdBmag[3]+3.464101615137754*BstarZdBmag[0])*hamil[4])*q_*rdvpar2))/(m_*q_); 
  alphaR[1] = (0.0125*(m_*((((114.5512985522207*b_x[5]+66.13622305514579*b_x[1])*jacobtot_inv[5]+66.13622305514579*jacobtot_inv[1]*b_x[5]+(63.63961030678928*b_x[3]+36.74234614174767*b_x[0])*jacobtot_inv[3]+36.74234614174767*jacobtot_inv[0]*b_x[3]+38.18376618407357*b_x[1]*jacobtot_inv[1]+21.21320343559643*b_x[0]*jacobtot_inv[0])*hamil[16]+((63.63961030678928*b_x[3]+36.74234614174767*b_x[0])*jacobtot_inv[5]+(63.63961030678928*jacobtot_inv[3]+36.74234614174767*jacobtot_inv[0])*b_x[5]+36.74234614174767*(b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+21.21320343559643*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1]))*hamil[8]+((66.13622305514579*b_x[5]+38.18376618407357*b_x[1])*jacobtot_inv[5]+38.18376618407357*jacobtot_inv[1]*b_x[5]+(36.74234614174767*b_x[3]+21.21320343559643*b_x[0])*jacobtot_inv[3]+21.21320343559643*jacobtot_inv[0]*b_x[3]+22.0454076850486*b_x[1]*jacobtot_inv[1]+12.24744871391589*b_x[0]*jacobtot_inv[0])*hamil[6]+hamil[2]*((36.74234614174767*b_x[3]+21.21320343559643*b_x[0])*jacobtot_inv[5]+(36.74234614174767*jacobtot_inv[3]+21.21320343559643*jacobtot_inv[0])*b_x[5]+21.21320343559643*(b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+12.24744871391589*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1])))*rdy2+((((-63.63961030678928*b_y[3])-36.74234614174767*b_y[0])*jacobtot_inv[5]+((-63.63961030678928*jacobtot_inv[3])-36.74234614174767*jacobtot_inv[0])*b_y[5]-36.74234614174767*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])-21.21320343559643*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[7]+hamil[1]*(((-36.74234614174767*b_y[3])-21.21320343559643*b_y[0])*jacobtot_inv[5]+((-36.74234614174767*jacobtot_inv[3])-21.21320343559643*jacobtot_inv[0])*b_y[5]-21.21320343559643*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])-12.24744871391589*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1])))*rdx2)+((67.0820393249937*BstarZdBmag[18]+38.72983346207418*BstarZdBmag[9])*hamil[32]+hamil[4]*(30.0*BstarZdBmag[7]+17.32050807568877*BstarZdBmag[1]))*q_*rdvpar2))/(m_*q_); 
  alphaR[2] = -(0.0625*(((12.72792206135786*b_y[5]+7.348469228349534*b_y[1])*jacobtot_inv[5]+7.348469228349534*jacobtot_inv[1]*b_y[5]+(12.72792206135786*b_y[3]+7.348469228349534*b_y[0])*jacobtot_inv[3]+7.348469228349534*jacobtot_inv[0]*b_y[3]+4.242640687119286*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[16]+((7.348469228349534*b_y[5]+4.242640687119286*b_y[1])*jacobtot_inv[5]+4.242640687119286*jacobtot_inv[1]*b_y[5]+(7.348469228349534*b_y[3]+4.242640687119286*b_y[0])*jacobtot_inv[3]+4.242640687119286*jacobtot_inv[0]*b_y[3]+2.449489742783178*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[6])*rdx2)/q_; 
  alphaR[3] = (0.125*((6.708203932499369*BstarZdBmag[3]+3.872983346207417*BstarZdBmag[0])*hamil[32]+hamil[4]*(3.0*BstarZdBmag[11]+1.732050807568877*BstarZdBmag[4]))*rdvpar2)/m_; 
  alphaR[4] = -(0.0625*(((12.72792206135786*b_y[5]+7.348469228349534*b_y[1])*jacobtot_inv[5]+7.348469228349534*jacobtot_inv[1]*b_y[5]+(12.72792206135786*b_y[3]+7.348469228349534*b_y[0])*jacobtot_inv[3]+7.348469228349534*jacobtot_inv[0]*b_y[3]+4.242640687119286*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[21]+((7.348469228349534*b_y[5]+4.242640687119286*b_y[1])*jacobtot_inv[5]+4.242640687119286*jacobtot_inv[1]*b_y[5]+(7.348469228349534*b_y[3]+4.242640687119286*b_y[0])*jacobtot_inv[3]+4.242640687119286*jacobtot_inv[0]*b_y[3]+2.449489742783178*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[12])*rdx2)/q_; 
  alphaR[5] = -(0.0625*(((12.72792206135786*b_y[3]+7.348469228349534*b_y[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[3]+7.348469228349534*jacobtot_inv[0])*b_y[5]+7.348469228349534*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])+4.242640687119286*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[16]+((7.348469228349534*b_y[3]+4.242640687119286*b_y[0])*jacobtot_inv[5]+(7.348469228349534*jacobtot_inv[3]+4.242640687119286*jacobtot_inv[0])*b_y[5]+4.242640687119286*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])+2.449489742783178*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[6])*rdx2)/q_; 
  alphaR[6] = (0.125*((6.708203932499369*BstarZdBmag[7]+3.872983346207417*BstarZdBmag[1])*hamil[32]+hamil[4]*(3.0*BstarZdBmag[18]+1.732050807568877*BstarZdBmag[9]))*rdvpar2)/m_; 
  alphaR[8] = -(0.0625*(((12.72792206135786*b_y[3]+7.348469228349534*b_y[0])*jacobtot_inv[5]+(12.72792206135786*jacobtot_inv[3]+7.348469228349534*jacobtot_inv[0])*b_y[5]+7.348469228349534*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])+4.242640687119286*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[21]+((7.348469228349534*b_y[3]+4.242640687119286*b_y[0])*jacobtot_inv[5]+(7.348469228349534*jacobtot_inv[3]+4.242640687119286*jacobtot_inv[0])*b_y[5]+4.242640687119286*(b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3])+2.449489742783178*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[12])*rdx2)/q_; 
  alphaR[16] = (0.25*(3.0*BstarZdBmag[11]+1.732050807568877*BstarZdBmag[4])*hamil[32]*rdvpar2)/m_; 
  alphaR[17] = (0.05*(15.0*BstarZdBmag[18]+8.660254037844387*BstarZdBmag[9])*hamil[32]*rdvpar2)/m_; 

  double fUpOrdL[24] = {0.};
  if ((-1.0*alphaL[17])+alphaL[16]+1.118033988749897*alphaL[8]+1.5*alphaL[6]+1.118033988749897*alphaL[5]-1.118033988749897*alphaL[4]-1.5*alphaL[3]-1.118033988749897*alphaL[2]-1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx3_eval_quad_node_0_r(fl); 
  } else { 
    fUpOrdL[0] = gkhyb_3x2v_p1_surfx3_eval_quad_node_0_l(fc); 
  } 
  if (alphaL[17]-1.0*alphaL[16]+0.8944271909999175*alphaL[8]+0.8944271909999175*alphaL[5]-0.8944271909999175*alphaL[4]-0.8944271909999175*alphaL[2]-0.8944271909999175*alphaL[1]+0.8944271909999175*alphaL[0] > 0.) {
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx3_eval_quad_node_1_r(fl); 
  } else { 
    fUpOrdL[1] = gkhyb_3x2v_p1_surfx3_eval_quad_node_1_l(fc); 
  } 
  if ((-1.0*alphaL[17])+alphaL[16]+1.118033988749897*alphaL[8]-1.5*alphaL[6]+1.118033988749897*alphaL[5]-1.118033988749897*alphaL[4]+1.5*alphaL[3]-1.118033988749897*alphaL[2]-1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx3_eval_quad_node_2_r(fl); 
  } else { 
    fUpOrdL[2] = gkhyb_3x2v_p1_surfx3_eval_quad_node_2_l(fc); 
  } 
  if ((-1.0*alphaL[17])+alphaL[16]-1.118033988749897*alphaL[8]+1.5*alphaL[6]+1.118033988749897*alphaL[5]+1.118033988749897*alphaL[4]-1.5*alphaL[3]-1.118033988749897*alphaL[2]-1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx3_eval_quad_node_3_r(fl); 
  } else { 
    fUpOrdL[3] = gkhyb_3x2v_p1_surfx3_eval_quad_node_3_l(fc); 
  } 
  if (alphaL[17]-1.0*alphaL[16]-0.8944271909999175*alphaL[8]+0.8944271909999175*alphaL[5]+0.8944271909999175*alphaL[4]-0.8944271909999175*alphaL[2]-0.8944271909999175*alphaL[1]+0.8944271909999175*alphaL[0] > 0.) {
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx3_eval_quad_node_4_r(fl); 
  } else { 
    fUpOrdL[4] = gkhyb_3x2v_p1_surfx3_eval_quad_node_4_l(fc); 
  } 
  if ((-1.0*alphaL[17])+alphaL[16]-1.118033988749897*alphaL[8]-1.5*alphaL[6]+1.118033988749897*alphaL[5]+1.118033988749897*alphaL[4]+1.5*alphaL[3]-1.118033988749897*alphaL[2]-1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx3_eval_quad_node_5_r(fl); 
  } else { 
    fUpOrdL[5] = gkhyb_3x2v_p1_surfx3_eval_quad_node_5_l(fc); 
  } 
  if ((-1.0*alphaL[17])+alphaL[16]+1.118033988749897*alphaL[8]+1.5*alphaL[6]-1.118033988749897*alphaL[5]-1.118033988749897*alphaL[4]-1.5*alphaL[3]+1.118033988749897*alphaL[2]-1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx3_eval_quad_node_6_r(fl); 
  } else { 
    fUpOrdL[6] = gkhyb_3x2v_p1_surfx3_eval_quad_node_6_l(fc); 
  } 
  if (alphaL[17]-1.0*alphaL[16]+0.8944271909999175*alphaL[8]-0.8944271909999175*alphaL[5]-0.8944271909999175*alphaL[4]+0.8944271909999175*alphaL[2]-0.8944271909999175*alphaL[1]+0.8944271909999175*alphaL[0] > 0.) {
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx3_eval_quad_node_7_r(fl); 
  } else { 
    fUpOrdL[7] = gkhyb_3x2v_p1_surfx3_eval_quad_node_7_l(fc); 
  } 
  if ((-1.0*alphaL[17])+alphaL[16]+1.118033988749897*alphaL[8]-1.5*alphaL[6]-1.118033988749897*alphaL[5]-1.118033988749897*alphaL[4]+1.5*alphaL[3]+1.118033988749897*alphaL[2]-1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx3_eval_quad_node_8_r(fl); 
  } else { 
    fUpOrdL[8] = gkhyb_3x2v_p1_surfx3_eval_quad_node_8_l(fc); 
  } 
  if ((-1.0*alphaL[17])+alphaL[16]-1.118033988749897*alphaL[8]+1.5*alphaL[6]-1.118033988749897*alphaL[5]+1.118033988749897*alphaL[4]-1.5*alphaL[3]+1.118033988749897*alphaL[2]-1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx3_eval_quad_node_9_r(fl); 
  } else { 
    fUpOrdL[9] = gkhyb_3x2v_p1_surfx3_eval_quad_node_9_l(fc); 
  } 
  if (alphaL[17]-1.0*alphaL[16]-0.8944271909999175*alphaL[8]-0.8944271909999175*alphaL[5]+0.8944271909999175*alphaL[4]+0.8944271909999175*alphaL[2]-0.8944271909999175*alphaL[1]+0.8944271909999175*alphaL[0] > 0.) {
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx3_eval_quad_node_10_r(fl); 
  } else { 
    fUpOrdL[10] = gkhyb_3x2v_p1_surfx3_eval_quad_node_10_l(fc); 
  } 
  if ((-1.0*alphaL[17])+alphaL[16]-1.118033988749897*alphaL[8]-1.5*alphaL[6]-1.118033988749897*alphaL[5]+1.118033988749897*alphaL[4]+1.5*alphaL[3]+1.118033988749897*alphaL[2]-1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx3_eval_quad_node_11_r(fl); 
  } else { 
    fUpOrdL[11] = gkhyb_3x2v_p1_surfx3_eval_quad_node_11_l(fc); 
  } 
  if (alphaL[17]+alphaL[16]-1.118033988749897*alphaL[8]-1.5*alphaL[6]-1.118033988749897*alphaL[5]-1.118033988749897*alphaL[4]-1.5*alphaL[3]-1.118033988749897*alphaL[2]+1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx3_eval_quad_node_12_r(fl); 
  } else { 
    fUpOrdL[12] = gkhyb_3x2v_p1_surfx3_eval_quad_node_12_l(fc); 
  } 
  if ((-1.0*alphaL[17])-1.0*alphaL[16]-0.8944271909999175*alphaL[8]-0.8944271909999175*alphaL[5]-0.8944271909999175*alphaL[4]-0.8944271909999175*alphaL[2]+0.8944271909999175*alphaL[1]+0.8944271909999175*alphaL[0] > 0.) {
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx3_eval_quad_node_13_r(fl); 
  } else { 
    fUpOrdL[13] = gkhyb_3x2v_p1_surfx3_eval_quad_node_13_l(fc); 
  } 
  if (alphaL[17]+alphaL[16]-1.118033988749897*alphaL[8]+1.5*alphaL[6]-1.118033988749897*alphaL[5]-1.118033988749897*alphaL[4]+1.5*alphaL[3]-1.118033988749897*alphaL[2]+1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx3_eval_quad_node_14_r(fl); 
  } else { 
    fUpOrdL[14] = gkhyb_3x2v_p1_surfx3_eval_quad_node_14_l(fc); 
  } 
  if (alphaL[17]+alphaL[16]+1.118033988749897*alphaL[8]-1.5*alphaL[6]-1.118033988749897*alphaL[5]+1.118033988749897*alphaL[4]-1.5*alphaL[3]-1.118033988749897*alphaL[2]+1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx3_eval_quad_node_15_r(fl); 
  } else { 
    fUpOrdL[15] = gkhyb_3x2v_p1_surfx3_eval_quad_node_15_l(fc); 
  } 
  if ((-1.0*alphaL[17])-1.0*alphaL[16]+0.8944271909999175*alphaL[8]-0.8944271909999175*alphaL[5]+0.8944271909999175*alphaL[4]-0.8944271909999175*alphaL[2]+0.8944271909999175*alphaL[1]+0.8944271909999175*alphaL[0] > 0.) {
    fUpOrdL[16] = gkhyb_3x2v_p1_surfx3_eval_quad_node_16_r(fl); 
  } else { 
    fUpOrdL[16] = gkhyb_3x2v_p1_surfx3_eval_quad_node_16_l(fc); 
  } 
  if (alphaL[17]+alphaL[16]+1.118033988749897*alphaL[8]+1.5*alphaL[6]-1.118033988749897*alphaL[5]+1.118033988749897*alphaL[4]+1.5*alphaL[3]-1.118033988749897*alphaL[2]+1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[17] = gkhyb_3x2v_p1_surfx3_eval_quad_node_17_r(fl); 
  } else { 
    fUpOrdL[17] = gkhyb_3x2v_p1_surfx3_eval_quad_node_17_l(fc); 
  } 
  if (alphaL[17]+alphaL[16]-1.118033988749897*alphaL[8]-1.5*alphaL[6]+1.118033988749897*alphaL[5]-1.118033988749897*alphaL[4]-1.5*alphaL[3]+1.118033988749897*alphaL[2]+1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[18] = gkhyb_3x2v_p1_surfx3_eval_quad_node_18_r(fl); 
  } else { 
    fUpOrdL[18] = gkhyb_3x2v_p1_surfx3_eval_quad_node_18_l(fc); 
  } 
  if ((-1.0*alphaL[17])-1.0*alphaL[16]-0.8944271909999175*alphaL[8]+0.8944271909999175*alphaL[5]-0.8944271909999175*alphaL[4]+0.8944271909999175*alphaL[2]+0.8944271909999175*alphaL[1]+0.8944271909999175*alphaL[0] > 0.) {
    fUpOrdL[19] = gkhyb_3x2v_p1_surfx3_eval_quad_node_19_r(fl); 
  } else { 
    fUpOrdL[19] = gkhyb_3x2v_p1_surfx3_eval_quad_node_19_l(fc); 
  } 
  if (alphaL[17]+alphaL[16]-1.118033988749897*alphaL[8]+1.5*alphaL[6]+1.118033988749897*alphaL[5]-1.118033988749897*alphaL[4]+1.5*alphaL[3]+1.118033988749897*alphaL[2]+1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[20] = gkhyb_3x2v_p1_surfx3_eval_quad_node_20_r(fl); 
  } else { 
    fUpOrdL[20] = gkhyb_3x2v_p1_surfx3_eval_quad_node_20_l(fc); 
  } 
  if (alphaL[17]+alphaL[16]+1.118033988749897*alphaL[8]-1.5*alphaL[6]+1.118033988749897*alphaL[5]+1.118033988749897*alphaL[4]-1.5*alphaL[3]+1.118033988749897*alphaL[2]+1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[21] = gkhyb_3x2v_p1_surfx3_eval_quad_node_21_r(fl); 
  } else { 
    fUpOrdL[21] = gkhyb_3x2v_p1_surfx3_eval_quad_node_21_l(fc); 
  } 
  if ((-1.0*alphaL[17])-1.0*alphaL[16]+0.8944271909999175*alphaL[8]+0.8944271909999175*alphaL[5]+0.8944271909999175*alphaL[4]+0.8944271909999175*alphaL[2]+0.8944271909999175*alphaL[1]+0.8944271909999175*alphaL[0] > 0.) {
    fUpOrdL[22] = gkhyb_3x2v_p1_surfx3_eval_quad_node_22_r(fl); 
  } else { 
    fUpOrdL[22] = gkhyb_3x2v_p1_surfx3_eval_quad_node_22_l(fc); 
  } 
  if (alphaL[17]+alphaL[16]+1.118033988749897*alphaL[8]+1.5*alphaL[6]+1.118033988749897*alphaL[5]+1.118033988749897*alphaL[4]+1.5*alphaL[3]+1.118033988749897*alphaL[2]+1.118033988749897*alphaL[1]+1.118033988749897*alphaL[0] > 0.) {
    fUpOrdL[23] = gkhyb_3x2v_p1_surfx3_eval_quad_node_23_r(fl); 
  } else { 
    fUpOrdL[23] = gkhyb_3x2v_p1_surfx3_eval_quad_node_23_l(fc); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[24] = {0.};
  gkhyb_3x2v_p1_xdir_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[48] = {0.}; 
  GhatL[0] = 0.25*alphaL[17]*fUpL[17]+0.25*alphaL[16]*fUpL[16]+0.25*alphaL[8]*fUpL[8]+0.25*alphaL[6]*fUpL[6]+0.25*alphaL[5]*fUpL[5]+0.25*alphaL[4]*fUpL[4]+0.25*alphaL[3]*fUpL[3]+0.25*alphaL[2]*fUpL[2]+0.25*alphaL[1]*fUpL[1]+0.25*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.2500000000000001*alphaL[16]*fUpL[17]+0.2500000000000001*fUpL[16]*alphaL[17]+0.25*alphaL[4]*fUpL[8]+0.25*fUpL[4]*alphaL[8]+0.25*alphaL[3]*fUpL[6]+0.25*fUpL[3]*alphaL[6]+0.25*alphaL[2]*fUpL[5]+0.25*fUpL[2]*alphaL[5]+0.25*alphaL[0]*fUpL[1]+0.25*fUpL[0]*alphaL[1]; 
  GhatL[2] = 0.2500000000000001*alphaL[17]*fUpL[20]+0.2500000000000001*alphaL[16]*fUpL[18]+0.25*alphaL[8]*fUpL[12]+0.25*alphaL[6]*fUpL[11]+0.25*alphaL[4]*fUpL[9]+0.25*alphaL[3]*fUpL[7]+0.25*alphaL[1]*fUpL[5]+0.25*fUpL[1]*alphaL[5]+0.25*alphaL[0]*fUpL[2]+0.25*fUpL[0]*alphaL[2]; 
  GhatL[3] = 0.223606797749979*alphaL[6]*fUpL[17]+0.223606797749979*fUpL[6]*alphaL[17]+0.223606797749979*alphaL[3]*fUpL[16]+0.223606797749979*fUpL[3]*alphaL[16]+0.25*alphaL[8]*fUpL[13]+0.25*alphaL[5]*fUpL[11]+0.25*alphaL[4]*fUpL[10]+0.25*alphaL[2]*fUpL[7]+0.25*alphaL[1]*fUpL[6]+0.25*fUpL[1]*alphaL[6]+0.25*alphaL[0]*fUpL[3]+0.25*fUpL[0]*alphaL[3]; 
  GhatL[4] = 0.2500000000000001*alphaL[17]*fUpL[21]+0.2500000000000001*alphaL[16]*fUpL[19]+0.25*alphaL[6]*fUpL[13]+0.25*alphaL[5]*fUpL[12]+0.25*alphaL[3]*fUpL[10]+0.25*alphaL[2]*fUpL[9]+0.25*alphaL[1]*fUpL[8]+0.25*fUpL[1]*alphaL[8]+0.25*alphaL[0]*fUpL[4]+0.25*fUpL[0]*alphaL[4]; 
  GhatL[5] = 0.25*alphaL[16]*fUpL[20]+0.25*alphaL[17]*fUpL[18]+0.25*alphaL[4]*fUpL[12]+0.25*alphaL[3]*fUpL[11]+0.25*alphaL[8]*fUpL[9]+0.25*alphaL[6]*fUpL[7]+0.25*alphaL[0]*fUpL[5]+0.25*fUpL[0]*alphaL[5]+0.25*alphaL[1]*fUpL[2]+0.25*fUpL[1]*alphaL[2]; 
  GhatL[6] = 0.223606797749979*alphaL[3]*fUpL[17]+0.223606797749979*fUpL[3]*alphaL[17]+0.223606797749979*alphaL[6]*fUpL[16]+0.223606797749979*fUpL[6]*alphaL[16]+0.25*alphaL[4]*fUpL[13]+0.25*alphaL[2]*fUpL[11]+0.25*alphaL[8]*fUpL[10]+0.25*alphaL[5]*fUpL[7]+0.25*alphaL[0]*fUpL[6]+0.25*fUpL[0]*alphaL[6]+0.25*alphaL[1]*fUpL[3]+0.25*fUpL[1]*alphaL[3]; 
  GhatL[7] = 0.223606797749979*alphaL[6]*fUpL[20]+0.223606797749979*alphaL[3]*fUpL[18]+0.223606797749979*fUpL[11]*alphaL[17]+0.223606797749979*fUpL[7]*alphaL[16]+0.25*alphaL[8]*fUpL[15]+0.25*alphaL[4]*fUpL[14]+0.25*alphaL[1]*fUpL[11]+0.25*alphaL[0]*fUpL[7]+0.25*alphaL[5]*fUpL[6]+0.25*fUpL[5]*alphaL[6]+0.25*alphaL[2]*fUpL[3]+0.25*fUpL[2]*alphaL[3]; 
  GhatL[8] = 0.25*alphaL[16]*fUpL[21]+0.25*alphaL[17]*fUpL[19]+0.25*alphaL[3]*fUpL[13]+0.25*alphaL[2]*fUpL[12]+0.25*alphaL[6]*fUpL[10]+0.25*alphaL[5]*fUpL[9]+0.25*alphaL[0]*fUpL[8]+0.25*fUpL[0]*alphaL[8]+0.25*alphaL[1]*fUpL[4]+0.25*fUpL[1]*alphaL[4]; 
  GhatL[9] = 0.25*alphaL[17]*fUpL[23]+0.25*alphaL[16]*fUpL[22]+0.25*alphaL[6]*fUpL[15]+0.25*alphaL[3]*fUpL[14]+0.25*alphaL[1]*fUpL[12]+0.25*alphaL[0]*fUpL[9]+0.25*alphaL[5]*fUpL[8]+0.25*fUpL[5]*alphaL[8]+0.25*alphaL[2]*fUpL[4]+0.25*fUpL[2]*alphaL[4]; 
  GhatL[10] = 0.223606797749979*alphaL[6]*fUpL[21]+0.223606797749979*alphaL[3]*fUpL[19]+0.223606797749979*fUpL[13]*alphaL[17]+0.223606797749979*fUpL[10]*alphaL[16]+0.25*alphaL[5]*fUpL[15]+0.25*alphaL[2]*fUpL[14]+0.25*alphaL[1]*fUpL[13]+0.25*alphaL[0]*fUpL[10]+0.25*alphaL[6]*fUpL[8]+0.25*fUpL[6]*alphaL[8]+0.25*alphaL[3]*fUpL[4]+0.25*fUpL[3]*alphaL[4]; 
  GhatL[11] = 0.223606797749979*alphaL[3]*fUpL[20]+0.223606797749979*alphaL[6]*fUpL[18]+0.223606797749979*fUpL[7]*alphaL[17]+0.223606797749979*fUpL[11]*alphaL[16]+0.25*alphaL[4]*fUpL[15]+0.25*alphaL[8]*fUpL[14]+0.25*alphaL[0]*fUpL[11]+0.25*alphaL[1]*fUpL[7]+0.25*alphaL[2]*fUpL[6]+0.25*fUpL[2]*alphaL[6]+0.25*alphaL[3]*fUpL[5]+0.25*fUpL[3]*alphaL[5]; 
  GhatL[12] = 0.2500000000000001*alphaL[16]*fUpL[23]+0.2500000000000001*alphaL[17]*fUpL[22]+0.25*alphaL[3]*fUpL[15]+0.25*alphaL[6]*fUpL[14]+0.25*alphaL[0]*fUpL[12]+0.25*alphaL[1]*fUpL[9]+0.25*alphaL[2]*fUpL[8]+0.25*fUpL[2]*alphaL[8]+0.25*alphaL[4]*fUpL[5]+0.25*fUpL[4]*alphaL[5]; 
  GhatL[13] = 0.223606797749979*alphaL[3]*fUpL[21]+0.223606797749979*alphaL[6]*fUpL[19]+0.223606797749979*fUpL[10]*alphaL[17]+0.223606797749979*fUpL[13]*alphaL[16]+0.25*alphaL[2]*fUpL[15]+0.25*alphaL[5]*fUpL[14]+0.25*alphaL[0]*fUpL[13]+0.25*alphaL[1]*fUpL[10]+0.25*alphaL[3]*fUpL[8]+0.25*fUpL[3]*alphaL[8]+0.25*alphaL[4]*fUpL[6]+0.25*fUpL[4]*alphaL[6]; 
  GhatL[14] = 0.223606797749979*alphaL[6]*fUpL[23]+0.223606797749979*alphaL[3]*fUpL[22]+0.223606797749979*fUpL[15]*alphaL[17]+0.223606797749979*fUpL[14]*alphaL[16]+0.25*alphaL[1]*fUpL[15]+0.25*alphaL[0]*fUpL[14]+0.25*alphaL[5]*fUpL[13]+0.25*alphaL[6]*fUpL[12]+0.25*alphaL[8]*fUpL[11]+0.25*alphaL[2]*fUpL[10]+0.25*alphaL[3]*fUpL[9]+0.25*alphaL[4]*fUpL[7]; 
  GhatL[15] = 0.223606797749979*alphaL[3]*fUpL[23]+0.223606797749979*alphaL[6]*fUpL[22]+0.223606797749979*fUpL[14]*alphaL[17]+0.223606797749979*fUpL[15]*alphaL[16]+0.25*alphaL[0]*fUpL[15]+0.25*alphaL[1]*fUpL[14]+0.25*alphaL[2]*fUpL[13]+0.25*alphaL[3]*fUpL[12]+0.25*alphaL[4]*fUpL[11]+0.25*alphaL[5]*fUpL[10]+0.25*alphaL[6]*fUpL[9]+0.25*fUpL[7]*alphaL[8]; 
  GhatL[16] = 0.25*alphaL[8]*fUpL[21]+0.25*alphaL[5]*fUpL[20]+0.2500000000000001*alphaL[4]*fUpL[19]+0.2500000000000001*alphaL[2]*fUpL[18]+0.159719141249985*alphaL[17]*fUpL[17]+0.2500000000000001*alphaL[1]*fUpL[17]+0.2500000000000001*fUpL[1]*alphaL[17]+0.159719141249985*alphaL[16]*fUpL[16]+0.25*alphaL[0]*fUpL[16]+0.25*fUpL[0]*alphaL[16]+0.223606797749979*alphaL[6]*fUpL[6]+0.223606797749979*alphaL[3]*fUpL[3]; 
  GhatL[17] = 0.2500000000000001*alphaL[4]*fUpL[21]+0.2500000000000001*alphaL[2]*fUpL[20]+0.25*alphaL[8]*fUpL[19]+0.25*alphaL[5]*fUpL[18]+0.159719141249985*alphaL[16]*fUpL[17]+0.25*alphaL[0]*fUpL[17]+0.159719141249985*fUpL[16]*alphaL[17]+0.25*fUpL[0]*alphaL[17]+0.2500000000000001*alphaL[1]*fUpL[16]+0.2500000000000001*fUpL[1]*alphaL[16]+0.223606797749979*alphaL[3]*fUpL[6]+0.223606797749979*fUpL[3]*alphaL[6]; 
  GhatL[18] = 0.25*alphaL[8]*fUpL[23]+0.2500000000000001*alphaL[4]*fUpL[22]+0.159719141249985*alphaL[17]*fUpL[20]+0.2500000000000001*alphaL[1]*fUpL[20]+0.159719141249985*alphaL[16]*fUpL[18]+0.25*alphaL[0]*fUpL[18]+0.25*alphaL[5]*fUpL[17]+0.25*fUpL[5]*alphaL[17]+0.2500000000000001*alphaL[2]*fUpL[16]+0.2500000000000001*fUpL[2]*alphaL[16]+0.223606797749979*alphaL[6]*fUpL[11]+0.223606797749979*alphaL[3]*fUpL[7]; 
  GhatL[19] = 0.25*alphaL[5]*fUpL[23]+0.2500000000000001*alphaL[2]*fUpL[22]+0.159719141249985*alphaL[17]*fUpL[21]+0.2500000000000001*alphaL[1]*fUpL[21]+0.159719141249985*alphaL[16]*fUpL[19]+0.25*alphaL[0]*fUpL[19]+0.25*alphaL[8]*fUpL[17]+0.25*fUpL[8]*alphaL[17]+0.2500000000000001*alphaL[4]*fUpL[16]+0.2500000000000001*fUpL[4]*alphaL[16]+0.223606797749979*alphaL[6]*fUpL[13]+0.223606797749979*alphaL[3]*fUpL[10]; 
  GhatL[20] = 0.2500000000000001*alphaL[4]*fUpL[23]+0.25*alphaL[8]*fUpL[22]+0.159719141249985*alphaL[16]*fUpL[20]+0.25*alphaL[0]*fUpL[20]+0.159719141249985*alphaL[17]*fUpL[18]+0.2500000000000001*alphaL[1]*fUpL[18]+0.2500000000000001*alphaL[2]*fUpL[17]+0.2500000000000001*fUpL[2]*alphaL[17]+0.25*alphaL[5]*fUpL[16]+0.25*fUpL[5]*alphaL[16]+0.223606797749979*alphaL[3]*fUpL[11]+0.223606797749979*alphaL[6]*fUpL[7]; 
  GhatL[21] = 0.2500000000000001*alphaL[2]*fUpL[23]+0.25*alphaL[5]*fUpL[22]+0.159719141249985*alphaL[16]*fUpL[21]+0.25*alphaL[0]*fUpL[21]+0.159719141249985*alphaL[17]*fUpL[19]+0.2500000000000001*alphaL[1]*fUpL[19]+0.2500000000000001*alphaL[4]*fUpL[17]+0.2500000000000001*fUpL[4]*alphaL[17]+0.25*alphaL[8]*fUpL[16]+0.25*fUpL[8]*alphaL[16]+0.223606797749979*alphaL[3]*fUpL[13]+0.223606797749979*alphaL[6]*fUpL[10]; 
  GhatL[22] = 0.159719141249985*alphaL[17]*fUpL[23]+0.2500000000000001*alphaL[1]*fUpL[23]+0.159719141249985*alphaL[16]*fUpL[22]+0.25*alphaL[0]*fUpL[22]+0.25*alphaL[5]*fUpL[21]+0.25*alphaL[8]*fUpL[20]+0.2500000000000001*alphaL[2]*fUpL[19]+0.2500000000000001*alphaL[4]*fUpL[18]+0.2500000000000001*fUpL[12]*alphaL[17]+0.25*fUpL[9]*alphaL[16]+0.223606797749979*alphaL[6]*fUpL[15]+0.223606797749979*alphaL[3]*fUpL[14]; 
  GhatL[23] = 0.159719141249985*alphaL[16]*fUpL[23]+0.25*alphaL[0]*fUpL[23]+0.159719141249985*alphaL[17]*fUpL[22]+0.2500000000000001*alphaL[1]*fUpL[22]+0.2500000000000001*alphaL[2]*fUpL[21]+0.2500000000000001*alphaL[4]*fUpL[20]+0.25*alphaL[5]*fUpL[19]+0.25*alphaL[8]*fUpL[18]+0.25*fUpL[9]*alphaL[17]+0.2500000000000001*fUpL[12]*alphaL[16]+0.223606797749979*alphaL[3]*fUpL[15]+0.223606797749979*alphaL[6]*fUpL[14]; 

  double fUpOrdR[24] = {0.};
  if ((-1.0*alphaR[17])+alphaR[16]+1.118033988749897*alphaR[8]+1.5*alphaR[6]+1.118033988749897*alphaR[5]-1.118033988749897*alphaR[4]-1.5*alphaR[3]-1.118033988749897*alphaR[2]-1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx3_eval_quad_node_0_r(fc); 
  } else { 
    fUpOrdR[0] = gkhyb_3x2v_p1_surfx3_eval_quad_node_0_l(fr); 
  } 
  if (alphaR[17]-1.0*alphaR[16]+0.8944271909999175*alphaR[8]+0.8944271909999175*alphaR[5]-0.8944271909999175*alphaR[4]-0.8944271909999175*alphaR[2]-0.8944271909999175*alphaR[1]+0.8944271909999175*alphaR[0] > 0.) {
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx3_eval_quad_node_1_r(fc); 
  } else { 
    fUpOrdR[1] = gkhyb_3x2v_p1_surfx3_eval_quad_node_1_l(fr); 
  } 
  if ((-1.0*alphaR[17])+alphaR[16]+1.118033988749897*alphaR[8]-1.5*alphaR[6]+1.118033988749897*alphaR[5]-1.118033988749897*alphaR[4]+1.5*alphaR[3]-1.118033988749897*alphaR[2]-1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx3_eval_quad_node_2_r(fc); 
  } else { 
    fUpOrdR[2] = gkhyb_3x2v_p1_surfx3_eval_quad_node_2_l(fr); 
  } 
  if ((-1.0*alphaR[17])+alphaR[16]-1.118033988749897*alphaR[8]+1.5*alphaR[6]+1.118033988749897*alphaR[5]+1.118033988749897*alphaR[4]-1.5*alphaR[3]-1.118033988749897*alphaR[2]-1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx3_eval_quad_node_3_r(fc); 
  } else { 
    fUpOrdR[3] = gkhyb_3x2v_p1_surfx3_eval_quad_node_3_l(fr); 
  } 
  if (alphaR[17]-1.0*alphaR[16]-0.8944271909999175*alphaR[8]+0.8944271909999175*alphaR[5]+0.8944271909999175*alphaR[4]-0.8944271909999175*alphaR[2]-0.8944271909999175*alphaR[1]+0.8944271909999175*alphaR[0] > 0.) {
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx3_eval_quad_node_4_r(fc); 
  } else { 
    fUpOrdR[4] = gkhyb_3x2v_p1_surfx3_eval_quad_node_4_l(fr); 
  } 
  if ((-1.0*alphaR[17])+alphaR[16]-1.118033988749897*alphaR[8]-1.5*alphaR[6]+1.118033988749897*alphaR[5]+1.118033988749897*alphaR[4]+1.5*alphaR[3]-1.118033988749897*alphaR[2]-1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx3_eval_quad_node_5_r(fc); 
  } else { 
    fUpOrdR[5] = gkhyb_3x2v_p1_surfx3_eval_quad_node_5_l(fr); 
  } 
  if ((-1.0*alphaR[17])+alphaR[16]+1.118033988749897*alphaR[8]+1.5*alphaR[6]-1.118033988749897*alphaR[5]-1.118033988749897*alphaR[4]-1.5*alphaR[3]+1.118033988749897*alphaR[2]-1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx3_eval_quad_node_6_r(fc); 
  } else { 
    fUpOrdR[6] = gkhyb_3x2v_p1_surfx3_eval_quad_node_6_l(fr); 
  } 
  if (alphaR[17]-1.0*alphaR[16]+0.8944271909999175*alphaR[8]-0.8944271909999175*alphaR[5]-0.8944271909999175*alphaR[4]+0.8944271909999175*alphaR[2]-0.8944271909999175*alphaR[1]+0.8944271909999175*alphaR[0] > 0.) {
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx3_eval_quad_node_7_r(fc); 
  } else { 
    fUpOrdR[7] = gkhyb_3x2v_p1_surfx3_eval_quad_node_7_l(fr); 
  } 
  if ((-1.0*alphaR[17])+alphaR[16]+1.118033988749897*alphaR[8]-1.5*alphaR[6]-1.118033988749897*alphaR[5]-1.118033988749897*alphaR[4]+1.5*alphaR[3]+1.118033988749897*alphaR[2]-1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx3_eval_quad_node_8_r(fc); 
  } else { 
    fUpOrdR[8] = gkhyb_3x2v_p1_surfx3_eval_quad_node_8_l(fr); 
  } 
  if ((-1.0*alphaR[17])+alphaR[16]-1.118033988749897*alphaR[8]+1.5*alphaR[6]-1.118033988749897*alphaR[5]+1.118033988749897*alphaR[4]-1.5*alphaR[3]+1.118033988749897*alphaR[2]-1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx3_eval_quad_node_9_r(fc); 
  } else { 
    fUpOrdR[9] = gkhyb_3x2v_p1_surfx3_eval_quad_node_9_l(fr); 
  } 
  if (alphaR[17]-1.0*alphaR[16]-0.8944271909999175*alphaR[8]-0.8944271909999175*alphaR[5]+0.8944271909999175*alphaR[4]+0.8944271909999175*alphaR[2]-0.8944271909999175*alphaR[1]+0.8944271909999175*alphaR[0] > 0.) {
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx3_eval_quad_node_10_r(fc); 
  } else { 
    fUpOrdR[10] = gkhyb_3x2v_p1_surfx3_eval_quad_node_10_l(fr); 
  } 
  if ((-1.0*alphaR[17])+alphaR[16]-1.118033988749897*alphaR[8]-1.5*alphaR[6]-1.118033988749897*alphaR[5]+1.118033988749897*alphaR[4]+1.5*alphaR[3]+1.118033988749897*alphaR[2]-1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx3_eval_quad_node_11_r(fc); 
  } else { 
    fUpOrdR[11] = gkhyb_3x2v_p1_surfx3_eval_quad_node_11_l(fr); 
  } 
  if (alphaR[17]+alphaR[16]-1.118033988749897*alphaR[8]-1.5*alphaR[6]-1.118033988749897*alphaR[5]-1.118033988749897*alphaR[4]-1.5*alphaR[3]-1.118033988749897*alphaR[2]+1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx3_eval_quad_node_12_r(fc); 
  } else { 
    fUpOrdR[12] = gkhyb_3x2v_p1_surfx3_eval_quad_node_12_l(fr); 
  } 
  if ((-1.0*alphaR[17])-1.0*alphaR[16]-0.8944271909999175*alphaR[8]-0.8944271909999175*alphaR[5]-0.8944271909999175*alphaR[4]-0.8944271909999175*alphaR[2]+0.8944271909999175*alphaR[1]+0.8944271909999175*alphaR[0] > 0.) {
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx3_eval_quad_node_13_r(fc); 
  } else { 
    fUpOrdR[13] = gkhyb_3x2v_p1_surfx3_eval_quad_node_13_l(fr); 
  } 
  if (alphaR[17]+alphaR[16]-1.118033988749897*alphaR[8]+1.5*alphaR[6]-1.118033988749897*alphaR[5]-1.118033988749897*alphaR[4]+1.5*alphaR[3]-1.118033988749897*alphaR[2]+1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx3_eval_quad_node_14_r(fc); 
  } else { 
    fUpOrdR[14] = gkhyb_3x2v_p1_surfx3_eval_quad_node_14_l(fr); 
  } 
  if (alphaR[17]+alphaR[16]+1.118033988749897*alphaR[8]-1.5*alphaR[6]-1.118033988749897*alphaR[5]+1.118033988749897*alphaR[4]-1.5*alphaR[3]-1.118033988749897*alphaR[2]+1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx3_eval_quad_node_15_r(fc); 
  } else { 
    fUpOrdR[15] = gkhyb_3x2v_p1_surfx3_eval_quad_node_15_l(fr); 
  } 
  if ((-1.0*alphaR[17])-1.0*alphaR[16]+0.8944271909999175*alphaR[8]-0.8944271909999175*alphaR[5]+0.8944271909999175*alphaR[4]-0.8944271909999175*alphaR[2]+0.8944271909999175*alphaR[1]+0.8944271909999175*alphaR[0] > 0.) {
    fUpOrdR[16] = gkhyb_3x2v_p1_surfx3_eval_quad_node_16_r(fc); 
  } else { 
    fUpOrdR[16] = gkhyb_3x2v_p1_surfx3_eval_quad_node_16_l(fr); 
  } 
  if (alphaR[17]+alphaR[16]+1.118033988749897*alphaR[8]+1.5*alphaR[6]-1.118033988749897*alphaR[5]+1.118033988749897*alphaR[4]+1.5*alphaR[3]-1.118033988749897*alphaR[2]+1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[17] = gkhyb_3x2v_p1_surfx3_eval_quad_node_17_r(fc); 
  } else { 
    fUpOrdR[17] = gkhyb_3x2v_p1_surfx3_eval_quad_node_17_l(fr); 
  } 
  if (alphaR[17]+alphaR[16]-1.118033988749897*alphaR[8]-1.5*alphaR[6]+1.118033988749897*alphaR[5]-1.118033988749897*alphaR[4]-1.5*alphaR[3]+1.118033988749897*alphaR[2]+1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[18] = gkhyb_3x2v_p1_surfx3_eval_quad_node_18_r(fc); 
  } else { 
    fUpOrdR[18] = gkhyb_3x2v_p1_surfx3_eval_quad_node_18_l(fr); 
  } 
  if ((-1.0*alphaR[17])-1.0*alphaR[16]-0.8944271909999175*alphaR[8]+0.8944271909999175*alphaR[5]-0.8944271909999175*alphaR[4]+0.8944271909999175*alphaR[2]+0.8944271909999175*alphaR[1]+0.8944271909999175*alphaR[0] > 0.) {
    fUpOrdR[19] = gkhyb_3x2v_p1_surfx3_eval_quad_node_19_r(fc); 
  } else { 
    fUpOrdR[19] = gkhyb_3x2v_p1_surfx3_eval_quad_node_19_l(fr); 
  } 
  if (alphaR[17]+alphaR[16]-1.118033988749897*alphaR[8]+1.5*alphaR[6]+1.118033988749897*alphaR[5]-1.118033988749897*alphaR[4]+1.5*alphaR[3]+1.118033988749897*alphaR[2]+1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[20] = gkhyb_3x2v_p1_surfx3_eval_quad_node_20_r(fc); 
  } else { 
    fUpOrdR[20] = gkhyb_3x2v_p1_surfx3_eval_quad_node_20_l(fr); 
  } 
  if (alphaR[17]+alphaR[16]+1.118033988749897*alphaR[8]-1.5*alphaR[6]+1.118033988749897*alphaR[5]+1.118033988749897*alphaR[4]-1.5*alphaR[3]+1.118033988749897*alphaR[2]+1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[21] = gkhyb_3x2v_p1_surfx3_eval_quad_node_21_r(fc); 
  } else { 
    fUpOrdR[21] = gkhyb_3x2v_p1_surfx3_eval_quad_node_21_l(fr); 
  } 
  if ((-1.0*alphaR[17])-1.0*alphaR[16]+0.8944271909999175*alphaR[8]+0.8944271909999175*alphaR[5]+0.8944271909999175*alphaR[4]+0.8944271909999175*alphaR[2]+0.8944271909999175*alphaR[1]+0.8944271909999175*alphaR[0] > 0.) {
    fUpOrdR[22] = gkhyb_3x2v_p1_surfx3_eval_quad_node_22_r(fc); 
  } else { 
    fUpOrdR[22] = gkhyb_3x2v_p1_surfx3_eval_quad_node_22_l(fr); 
  } 
  if (alphaR[17]+alphaR[16]+1.118033988749897*alphaR[8]+1.5*alphaR[6]+1.118033988749897*alphaR[5]+1.118033988749897*alphaR[4]+1.5*alphaR[3]+1.118033988749897*alphaR[2]+1.118033988749897*alphaR[1]+1.118033988749897*alphaR[0] > 0.) {
    fUpOrdR[23] = gkhyb_3x2v_p1_surfx3_eval_quad_node_23_r(fc); 
  } else { 
    fUpOrdR[23] = gkhyb_3x2v_p1_surfx3_eval_quad_node_23_l(fr); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[24] = {0.};
  gkhyb_3x2v_p1_xdir_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[48] = {0.}; 
  GhatR[0] = 0.25*alphaR[17]*fUpR[17]+0.25*alphaR[16]*fUpR[16]+0.25*alphaR[8]*fUpR[8]+0.25*alphaR[6]*fUpR[6]+0.25*alphaR[5]*fUpR[5]+0.25*alphaR[4]*fUpR[4]+0.25*alphaR[3]*fUpR[3]+0.25*alphaR[2]*fUpR[2]+0.25*alphaR[1]*fUpR[1]+0.25*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.2500000000000001*alphaR[16]*fUpR[17]+0.2500000000000001*fUpR[16]*alphaR[17]+0.25*alphaR[4]*fUpR[8]+0.25*fUpR[4]*alphaR[8]+0.25*alphaR[3]*fUpR[6]+0.25*fUpR[3]*alphaR[6]+0.25*alphaR[2]*fUpR[5]+0.25*fUpR[2]*alphaR[5]+0.25*alphaR[0]*fUpR[1]+0.25*fUpR[0]*alphaR[1]; 
  GhatR[2] = 0.2500000000000001*alphaR[17]*fUpR[20]+0.2500000000000001*alphaR[16]*fUpR[18]+0.25*alphaR[8]*fUpR[12]+0.25*alphaR[6]*fUpR[11]+0.25*alphaR[4]*fUpR[9]+0.25*alphaR[3]*fUpR[7]+0.25*alphaR[1]*fUpR[5]+0.25*fUpR[1]*alphaR[5]+0.25*alphaR[0]*fUpR[2]+0.25*fUpR[0]*alphaR[2]; 
  GhatR[3] = 0.223606797749979*alphaR[6]*fUpR[17]+0.223606797749979*fUpR[6]*alphaR[17]+0.223606797749979*alphaR[3]*fUpR[16]+0.223606797749979*fUpR[3]*alphaR[16]+0.25*alphaR[8]*fUpR[13]+0.25*alphaR[5]*fUpR[11]+0.25*alphaR[4]*fUpR[10]+0.25*alphaR[2]*fUpR[7]+0.25*alphaR[1]*fUpR[6]+0.25*fUpR[1]*alphaR[6]+0.25*alphaR[0]*fUpR[3]+0.25*fUpR[0]*alphaR[3]; 
  GhatR[4] = 0.2500000000000001*alphaR[17]*fUpR[21]+0.2500000000000001*alphaR[16]*fUpR[19]+0.25*alphaR[6]*fUpR[13]+0.25*alphaR[5]*fUpR[12]+0.25*alphaR[3]*fUpR[10]+0.25*alphaR[2]*fUpR[9]+0.25*alphaR[1]*fUpR[8]+0.25*fUpR[1]*alphaR[8]+0.25*alphaR[0]*fUpR[4]+0.25*fUpR[0]*alphaR[4]; 
  GhatR[5] = 0.25*alphaR[16]*fUpR[20]+0.25*alphaR[17]*fUpR[18]+0.25*alphaR[4]*fUpR[12]+0.25*alphaR[3]*fUpR[11]+0.25*alphaR[8]*fUpR[9]+0.25*alphaR[6]*fUpR[7]+0.25*alphaR[0]*fUpR[5]+0.25*fUpR[0]*alphaR[5]+0.25*alphaR[1]*fUpR[2]+0.25*fUpR[1]*alphaR[2]; 
  GhatR[6] = 0.223606797749979*alphaR[3]*fUpR[17]+0.223606797749979*fUpR[3]*alphaR[17]+0.223606797749979*alphaR[6]*fUpR[16]+0.223606797749979*fUpR[6]*alphaR[16]+0.25*alphaR[4]*fUpR[13]+0.25*alphaR[2]*fUpR[11]+0.25*alphaR[8]*fUpR[10]+0.25*alphaR[5]*fUpR[7]+0.25*alphaR[0]*fUpR[6]+0.25*fUpR[0]*alphaR[6]+0.25*alphaR[1]*fUpR[3]+0.25*fUpR[1]*alphaR[3]; 
  GhatR[7] = 0.223606797749979*alphaR[6]*fUpR[20]+0.223606797749979*alphaR[3]*fUpR[18]+0.223606797749979*fUpR[11]*alphaR[17]+0.223606797749979*fUpR[7]*alphaR[16]+0.25*alphaR[8]*fUpR[15]+0.25*alphaR[4]*fUpR[14]+0.25*alphaR[1]*fUpR[11]+0.25*alphaR[0]*fUpR[7]+0.25*alphaR[5]*fUpR[6]+0.25*fUpR[5]*alphaR[6]+0.25*alphaR[2]*fUpR[3]+0.25*fUpR[2]*alphaR[3]; 
  GhatR[8] = 0.25*alphaR[16]*fUpR[21]+0.25*alphaR[17]*fUpR[19]+0.25*alphaR[3]*fUpR[13]+0.25*alphaR[2]*fUpR[12]+0.25*alphaR[6]*fUpR[10]+0.25*alphaR[5]*fUpR[9]+0.25*alphaR[0]*fUpR[8]+0.25*fUpR[0]*alphaR[8]+0.25*alphaR[1]*fUpR[4]+0.25*fUpR[1]*alphaR[4]; 
  GhatR[9] = 0.25*alphaR[17]*fUpR[23]+0.25*alphaR[16]*fUpR[22]+0.25*alphaR[6]*fUpR[15]+0.25*alphaR[3]*fUpR[14]+0.25*alphaR[1]*fUpR[12]+0.25*alphaR[0]*fUpR[9]+0.25*alphaR[5]*fUpR[8]+0.25*fUpR[5]*alphaR[8]+0.25*alphaR[2]*fUpR[4]+0.25*fUpR[2]*alphaR[4]; 
  GhatR[10] = 0.223606797749979*alphaR[6]*fUpR[21]+0.223606797749979*alphaR[3]*fUpR[19]+0.223606797749979*fUpR[13]*alphaR[17]+0.223606797749979*fUpR[10]*alphaR[16]+0.25*alphaR[5]*fUpR[15]+0.25*alphaR[2]*fUpR[14]+0.25*alphaR[1]*fUpR[13]+0.25*alphaR[0]*fUpR[10]+0.25*alphaR[6]*fUpR[8]+0.25*fUpR[6]*alphaR[8]+0.25*alphaR[3]*fUpR[4]+0.25*fUpR[3]*alphaR[4]; 
  GhatR[11] = 0.223606797749979*alphaR[3]*fUpR[20]+0.223606797749979*alphaR[6]*fUpR[18]+0.223606797749979*fUpR[7]*alphaR[17]+0.223606797749979*fUpR[11]*alphaR[16]+0.25*alphaR[4]*fUpR[15]+0.25*alphaR[8]*fUpR[14]+0.25*alphaR[0]*fUpR[11]+0.25*alphaR[1]*fUpR[7]+0.25*alphaR[2]*fUpR[6]+0.25*fUpR[2]*alphaR[6]+0.25*alphaR[3]*fUpR[5]+0.25*fUpR[3]*alphaR[5]; 
  GhatR[12] = 0.2500000000000001*alphaR[16]*fUpR[23]+0.2500000000000001*alphaR[17]*fUpR[22]+0.25*alphaR[3]*fUpR[15]+0.25*alphaR[6]*fUpR[14]+0.25*alphaR[0]*fUpR[12]+0.25*alphaR[1]*fUpR[9]+0.25*alphaR[2]*fUpR[8]+0.25*fUpR[2]*alphaR[8]+0.25*alphaR[4]*fUpR[5]+0.25*fUpR[4]*alphaR[5]; 
  GhatR[13] = 0.223606797749979*alphaR[3]*fUpR[21]+0.223606797749979*alphaR[6]*fUpR[19]+0.223606797749979*fUpR[10]*alphaR[17]+0.223606797749979*fUpR[13]*alphaR[16]+0.25*alphaR[2]*fUpR[15]+0.25*alphaR[5]*fUpR[14]+0.25*alphaR[0]*fUpR[13]+0.25*alphaR[1]*fUpR[10]+0.25*alphaR[3]*fUpR[8]+0.25*fUpR[3]*alphaR[8]+0.25*alphaR[4]*fUpR[6]+0.25*fUpR[4]*alphaR[6]; 
  GhatR[14] = 0.223606797749979*alphaR[6]*fUpR[23]+0.223606797749979*alphaR[3]*fUpR[22]+0.223606797749979*fUpR[15]*alphaR[17]+0.223606797749979*fUpR[14]*alphaR[16]+0.25*alphaR[1]*fUpR[15]+0.25*alphaR[0]*fUpR[14]+0.25*alphaR[5]*fUpR[13]+0.25*alphaR[6]*fUpR[12]+0.25*alphaR[8]*fUpR[11]+0.25*alphaR[2]*fUpR[10]+0.25*alphaR[3]*fUpR[9]+0.25*alphaR[4]*fUpR[7]; 
  GhatR[15] = 0.223606797749979*alphaR[3]*fUpR[23]+0.223606797749979*alphaR[6]*fUpR[22]+0.223606797749979*fUpR[14]*alphaR[17]+0.223606797749979*fUpR[15]*alphaR[16]+0.25*alphaR[0]*fUpR[15]+0.25*alphaR[1]*fUpR[14]+0.25*alphaR[2]*fUpR[13]+0.25*alphaR[3]*fUpR[12]+0.25*alphaR[4]*fUpR[11]+0.25*alphaR[5]*fUpR[10]+0.25*alphaR[6]*fUpR[9]+0.25*fUpR[7]*alphaR[8]; 
  GhatR[16] = 0.25*alphaR[8]*fUpR[21]+0.25*alphaR[5]*fUpR[20]+0.2500000000000001*alphaR[4]*fUpR[19]+0.2500000000000001*alphaR[2]*fUpR[18]+0.159719141249985*alphaR[17]*fUpR[17]+0.2500000000000001*alphaR[1]*fUpR[17]+0.2500000000000001*fUpR[1]*alphaR[17]+0.159719141249985*alphaR[16]*fUpR[16]+0.25*alphaR[0]*fUpR[16]+0.25*fUpR[0]*alphaR[16]+0.223606797749979*alphaR[6]*fUpR[6]+0.223606797749979*alphaR[3]*fUpR[3]; 
  GhatR[17] = 0.2500000000000001*alphaR[4]*fUpR[21]+0.2500000000000001*alphaR[2]*fUpR[20]+0.25*alphaR[8]*fUpR[19]+0.25*alphaR[5]*fUpR[18]+0.159719141249985*alphaR[16]*fUpR[17]+0.25*alphaR[0]*fUpR[17]+0.159719141249985*fUpR[16]*alphaR[17]+0.25*fUpR[0]*alphaR[17]+0.2500000000000001*alphaR[1]*fUpR[16]+0.2500000000000001*fUpR[1]*alphaR[16]+0.223606797749979*alphaR[3]*fUpR[6]+0.223606797749979*fUpR[3]*alphaR[6]; 
  GhatR[18] = 0.25*alphaR[8]*fUpR[23]+0.2500000000000001*alphaR[4]*fUpR[22]+0.159719141249985*alphaR[17]*fUpR[20]+0.2500000000000001*alphaR[1]*fUpR[20]+0.159719141249985*alphaR[16]*fUpR[18]+0.25*alphaR[0]*fUpR[18]+0.25*alphaR[5]*fUpR[17]+0.25*fUpR[5]*alphaR[17]+0.2500000000000001*alphaR[2]*fUpR[16]+0.2500000000000001*fUpR[2]*alphaR[16]+0.223606797749979*alphaR[6]*fUpR[11]+0.223606797749979*alphaR[3]*fUpR[7]; 
  GhatR[19] = 0.25*alphaR[5]*fUpR[23]+0.2500000000000001*alphaR[2]*fUpR[22]+0.159719141249985*alphaR[17]*fUpR[21]+0.2500000000000001*alphaR[1]*fUpR[21]+0.159719141249985*alphaR[16]*fUpR[19]+0.25*alphaR[0]*fUpR[19]+0.25*alphaR[8]*fUpR[17]+0.25*fUpR[8]*alphaR[17]+0.2500000000000001*alphaR[4]*fUpR[16]+0.2500000000000001*fUpR[4]*alphaR[16]+0.223606797749979*alphaR[6]*fUpR[13]+0.223606797749979*alphaR[3]*fUpR[10]; 
  GhatR[20] = 0.2500000000000001*alphaR[4]*fUpR[23]+0.25*alphaR[8]*fUpR[22]+0.159719141249985*alphaR[16]*fUpR[20]+0.25*alphaR[0]*fUpR[20]+0.159719141249985*alphaR[17]*fUpR[18]+0.2500000000000001*alphaR[1]*fUpR[18]+0.2500000000000001*alphaR[2]*fUpR[17]+0.2500000000000001*fUpR[2]*alphaR[17]+0.25*alphaR[5]*fUpR[16]+0.25*fUpR[5]*alphaR[16]+0.223606797749979*alphaR[3]*fUpR[11]+0.223606797749979*alphaR[6]*fUpR[7]; 
  GhatR[21] = 0.2500000000000001*alphaR[2]*fUpR[23]+0.25*alphaR[5]*fUpR[22]+0.159719141249985*alphaR[16]*fUpR[21]+0.25*alphaR[0]*fUpR[21]+0.159719141249985*alphaR[17]*fUpR[19]+0.2500000000000001*alphaR[1]*fUpR[19]+0.2500000000000001*alphaR[4]*fUpR[17]+0.2500000000000001*fUpR[4]*alphaR[17]+0.25*alphaR[8]*fUpR[16]+0.25*fUpR[8]*alphaR[16]+0.223606797749979*alphaR[3]*fUpR[13]+0.223606797749979*alphaR[6]*fUpR[10]; 
  GhatR[22] = 0.159719141249985*alphaR[17]*fUpR[23]+0.2500000000000001*alphaR[1]*fUpR[23]+0.159719141249985*alphaR[16]*fUpR[22]+0.25*alphaR[0]*fUpR[22]+0.25*alphaR[5]*fUpR[21]+0.25*alphaR[8]*fUpR[20]+0.2500000000000001*alphaR[2]*fUpR[19]+0.2500000000000001*alphaR[4]*fUpR[18]+0.2500000000000001*fUpR[12]*alphaR[17]+0.25*fUpR[9]*alphaR[16]+0.223606797749979*alphaR[6]*fUpR[15]+0.223606797749979*alphaR[3]*fUpR[14]; 
  GhatR[23] = 0.159719141249985*alphaR[16]*fUpR[23]+0.25*alphaR[0]*fUpR[23]+0.159719141249985*alphaR[17]*fUpR[22]+0.2500000000000001*alphaR[1]*fUpR[22]+0.2500000000000001*alphaR[2]*fUpR[21]+0.2500000000000001*alphaR[4]*fUpR[20]+0.25*alphaR[5]*fUpR[19]+0.25*alphaR[8]*fUpR[18]+0.25*fUpR[9]*alphaR[17]+0.2500000000000001*fUpR[12]*alphaR[16]+0.223606797749979*alphaR[3]*fUpR[15]+0.223606797749979*alphaR[6]*fUpR[14]; 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdz2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdz2; 
  out[2] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdz2; 
  out[3] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdz2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdz2; 
  out[5] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdz2; 
  out[6] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdz2; 
  out[7] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdz2; 
  out[8] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdz2; 
  out[9] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdz2; 
  out[10] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdz2; 
  out[11] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdz2; 
  out[12] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdz2; 
  out[13] += (0.7071067811865475*GhatL[9]-0.7071067811865475*GhatR[9])*rdz2; 
  out[14] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdz2; 
  out[15] += (0.7071067811865475*GhatL[10]-0.7071067811865475*GhatR[10])*rdz2; 
  out[16] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdz2; 
  out[17] += (0.7071067811865475*GhatL[11]-0.7071067811865475*GhatR[11])*rdz2; 
  out[18] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdz2; 
  out[19] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdz2; 
  out[20] += (0.7071067811865475*GhatL[12]-0.7071067811865475*GhatR[12])*rdz2; 
  out[21] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdz2; 
  out[22] += ((-1.224744871391589*GhatR[9])-1.224744871391589*GhatL[9])*rdz2; 
  out[23] += (0.7071067811865475*GhatL[13]-0.7071067811865475*GhatR[13])*rdz2; 
  out[24] += (0.7071067811865475*GhatL[14]-0.7071067811865475*GhatR[14])*rdz2; 
  out[25] += ((-1.224744871391589*GhatR[10])-1.224744871391589*GhatL[10])*rdz2; 
  out[26] += ((-1.224744871391589*GhatR[11])-1.224744871391589*GhatL[11])*rdz2; 
  out[27] += ((-1.224744871391589*GhatR[12])-1.224744871391589*GhatL[12])*rdz2; 
  out[28] += (0.7071067811865475*GhatL[15]-0.7071067811865475*GhatR[15])*rdz2; 
  out[29] += ((-1.224744871391589*GhatR[13])-1.224744871391589*GhatL[13])*rdz2; 
  out[30] += ((-1.224744871391589*GhatR[14])-1.224744871391589*GhatL[14])*rdz2; 
  out[31] += ((-1.224744871391589*GhatR[15])-1.224744871391589*GhatL[15])*rdz2; 
  out[32] += (0.7071067811865475*GhatL[16]-0.7071067811865475*GhatR[16])*rdz2; 
  out[33] += (0.7071067811865475*GhatL[17]-0.7071067811865475*GhatR[17])*rdz2; 
  out[34] += (0.7071067811865475*GhatL[18]-0.7071067811865475*GhatR[18])*rdz2; 
  out[35] += ((-1.224744871391589*GhatR[16])-1.224744871391589*GhatL[16])*rdz2; 
  out[36] += (0.7071067811865475*GhatL[19]-0.7071067811865475*GhatR[19])*rdz2; 
  out[37] += (0.7071067811865475*GhatL[20]-0.7071067811865475*GhatR[20])*rdz2; 
  out[38] += ((-1.224744871391589*GhatR[17])-1.224744871391589*GhatL[17])*rdz2; 
  out[39] += ((-1.224744871391589*GhatR[18])-1.224744871391589*GhatL[18])*rdz2; 
  out[40] += (0.7071067811865475*GhatL[21]-0.7071067811865475*GhatR[21])*rdz2; 
  out[41] += (0.7071067811865475*GhatL[22]-0.7071067811865475*GhatR[22])*rdz2; 
  out[42] += ((-1.224744871391589*GhatR[19])-1.224744871391589*GhatL[19])*rdz2; 
  out[43] += ((-1.224744871391589*GhatR[20])-1.224744871391589*GhatL[20])*rdz2; 
  out[44] += (0.7071067811865475*GhatL[23]-0.7071067811865475*GhatR[23])*rdz2; 
  out[45] += ((-1.224744871391589*GhatR[21])-1.224744871391589*GhatL[21])*rdz2; 
  out[46] += ((-1.224744871391589*GhatR[22])-1.224744871391589*GhatL[22])*rdz2; 
  out[47] += ((-1.224744871391589*GhatR[23])-1.224744871391589*GhatL[23])*rdz2; 

} 
