#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_no_by_vol_2x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out) 
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
  // fin: Distribution function.
  // out: output increment.

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
  hamil[0] = 2.0*m_*wvparSq+2.0*bmag[0]*wmu+(0.6666666666666666*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*bmag[1]*wmu+2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[16] = (0.5962847939999438*m_)/rdvpar2Sq; 

  double alphax[24] = {0.}; 
  alphax[0] = (((-0.4330127018922193*b_z[0]*jacobtot_inv[1]*hamil[5])-0.4330127018922193*jacobtot_inv[0]*b_z[1]*hamil[5]-0.4330127018922193*b_z[1]*jacobtot_inv[1]*hamil[2]-0.4330127018922193*b_z[0]*jacobtot_inv[0]*hamil[2])*rdx2*rdy2)/q_; 
  alphax[1] = (((-0.7794228634059945*b_z[1]*jacobtot_inv[1]*hamil[5])-0.4330127018922193*b_z[0]*jacobtot_inv[0]*hamil[5]-0.4330127018922193*b_z[0]*jacobtot_inv[1]*hamil[2]-0.4330127018922193*jacobtot_inv[0]*b_z[1]*hamil[2])*rdx2*rdy2)/q_; 

  double alphay[24] = {0.}; 
  alphay[0] = (((-0.9682458365518543*jacobtot_inv[0]*b_z[1]*hamil[16])+0.4330127018922193*b_z[1]*hamil[1]*jacobtot_inv[1]+0.4330127018922193*b_z[0]*jacobtot_inv[0]*hamil[1])*rdx2*rdy2)/q_-(0.75*jacobtot_inv[0]*b_z[1]*hamil[3]*rdvpar2*rdx2*rdy2*wvpar)/q_; 
  alphay[1] = (((-0.9682458365518543*b_z[1]*jacobtot_inv[1]*hamil[16])+0.4330127018922193*b_z[0]*hamil[1]*jacobtot_inv[1]+0.4330127018922193*jacobtot_inv[0]*b_z[1]*hamil[1])*rdx2*rdy2)/q_-(0.75*b_z[1]*jacobtot_inv[1]*hamil[3]*rdvpar2*rdx2*rdy2*wvpar)/q_; 
  alphay[2] = ((0.4330127018922193*b_z[1]*jacobtot_inv[1]*hamil[5]+0.4330127018922193*b_z[0]*jacobtot_inv[0]*hamil[5])*rdx2*rdy2)/q_; 
  alphay[3] = (-(1.677050983124842*jacobtot_inv[0]*b_z[1]*hamil[16]*rdvpar2*rdx2*rdy2*wvpar)/q_)-(0.4330127018922193*jacobtot_inv[0]*b_z[1]*hamil[3]*rdx2*rdy2)/q_; 
  alphay[4] = ((0.4330127018922193*b_z[1]*jacobtot_inv[1]*hamil[8]+0.4330127018922193*b_z[0]*jacobtot_inv[0]*hamil[8])*rdx2*rdy2)/q_; 
  alphay[5] = ((0.4330127018922193*b_z[0]*jacobtot_inv[1]*hamil[5]+0.4330127018922193*jacobtot_inv[0]*b_z[1]*hamil[5])*rdx2*rdy2)/q_; 
  alphay[6] = (-(1.677050983124842*b_z[1]*jacobtot_inv[1]*hamil[16]*rdvpar2*rdx2*rdy2*wvpar)/q_)-(0.4330127018922193*b_z[1]*jacobtot_inv[1]*hamil[3]*rdx2*rdy2)/q_; 
  alphay[8] = ((0.4330127018922193*b_z[0]*jacobtot_inv[1]*hamil[8]+0.4330127018922193*jacobtot_inv[0]*b_z[1]*hamil[8])*rdx2*rdy2)/q_; 
  alphay[16] = -(0.8660254037844386*jacobtot_inv[0]*b_z[1]*hamil[16]*rdx2*rdy2)/q_; 
  alphay[17] = -(0.8660254037844387*b_z[1]*jacobtot_inv[1]*hamil[16]*rdx2*rdy2)/q_; 

  double alphavpar[24] = {0.}; 
  alphavpar[0] = ((0.75*b_z[1]*jacobtot_inv[1]*hamil[5]+0.75*jacobtot_inv[0]*b_z[1]*hamil[2])*rdvpar2*rdx2*rdy2*wvpar)/q_; 
  alphavpar[1] = ((0.75*jacobtot_inv[0]*b_z[1]*hamil[5]+0.75*b_z[1]*jacobtot_inv[1]*hamil[2])*rdvpar2*rdx2*rdy2*wvpar)/q_; 
  alphavpar[3] = ((0.4330127018922193*b_z[1]*jacobtot_inv[1]*hamil[5]+0.4330127018922193*jacobtot_inv[0]*b_z[1]*hamil[2])*rdx2*rdy2)/q_; 
  alphavpar[6] = ((0.4330127018922193*jacobtot_inv[0]*b_z[1]*hamil[5]+0.4330127018922193*b_z[1]*jacobtot_inv[1]*hamil[2])*rdx2*rdy2)/q_; 

  out[1] += 0.4330127018922193*(alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.4330127018922193*(alphay[17]*fin[17]+alphay[16]*fin[16]+alphay[8]*fin[8]+alphay[6]*fin[6]+alphay[5]*fin[5]+alphay[4]*fin[4]+alphay[3]*fin[3]+alphay[2]*fin[2]+alphay[1]*fin[1]+alphay[0]*fin[0]); 
  out[3] += 0.4330127018922193*(alphavpar[6]*fin[6]+alphavpar[3]*fin[3]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[5] += 0.4330127018922194*(alphay[16]*fin[17]+fin[16]*alphay[17])+0.4330127018922193*(alphay[4]*fin[8]+fin[4]*alphay[8]+alphay[3]*fin[6]+fin[3]*alphay[6]+(alphay[2]+alphax[1])*fin[5]+fin[2]*(alphay[5]+alphax[0])+alphay[0]*fin[1]+fin[0]*alphay[1]); 
  out[6] += 0.4330127018922193*((alphavpar[3]+alphax[1])*fin[6]+fin[3]*(alphavpar[6]+alphax[0])+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[7] += 0.3872983346207417*(alphay[6]*fin[17]+fin[6]*alphay[17])+0.3872983346207416*(alphay[3]*fin[16]+fin[3]*alphay[16])+0.4330127018922193*(alphay[8]*fin[13]+(alphavpar[6]+alphay[5])*fin[11]+alphay[4]*fin[10]+(alphavpar[3]+alphay[2])*fin[7]+alphay[1]*fin[6]+fin[1]*alphay[6]+alphavpar[1]*fin[5]+alphay[0]*fin[3]+fin[0]*alphay[3]+alphavpar[0]*fin[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*fin[8]+alphax[0]*fin[4]); 
  out[9] += 0.4330127018922194*(alphay[17]*fin[21]+alphay[16]*fin[19])+0.4330127018922193*(alphay[6]*fin[13]+alphay[5]*fin[12]+alphay[3]*fin[10]+alphay[2]*fin[9]+alphay[1]*fin[8]+fin[1]*alphay[8]+alphay[0]*fin[4]+fin[0]*alphay[4]); 
  out[10] += 0.4330127018922193*(alphavpar[6]*fin[13]+alphavpar[3]*fin[10]+alphavpar[1]*fin[8]+alphavpar[0]*fin[4]); 
  out[11] += 0.3872983346207417*(alphay[3]*fin[17]+fin[3]*alphay[17])+0.3872983346207416*(alphay[6]*fin[16]+fin[6]*alphay[16])+0.4330127018922193*(alphay[4]*fin[13]+(alphavpar[3]+alphay[2]+alphax[1])*fin[11]+alphay[8]*fin[10]+(alphavpar[6]+alphay[5]+alphax[0])*fin[7]+alphay[0]*fin[6]+fin[0]*alphay[6]+alphavpar[0]*fin[5]+alphay[1]*fin[3]+fin[1]*alphay[3]+alphavpar[1]*fin[2]); 
  out[12] += 0.4330127018922193*(alphay[16]*fin[21]+alphay[17]*fin[19]+alphay[3]*fin[13]+(alphay[2]+alphax[1])*fin[12]+alphay[6]*fin[10]+(alphay[5]+alphax[0])*fin[9]+alphay[0]*fin[8]+fin[0]*alphay[8]+alphay[1]*fin[4]+fin[1]*alphay[4]); 
  out[13] += 0.4330127018922193*((alphavpar[3]+alphax[1])*fin[13]+(alphavpar[6]+alphax[0])*fin[10]+alphavpar[0]*fin[8]+alphavpar[1]*fin[4]); 
  out[14] += 0.3872983346207416*alphay[6]*fin[21]+0.3872983346207417*(alphay[3]*fin[19]+fin[13]*alphay[17])+0.3872983346207416*fin[10]*alphay[16]+0.4330127018922193*((alphavpar[6]+alphay[5])*fin[15]+(alphavpar[3]+alphay[2])*fin[14]+alphay[1]*fin[13]+alphavpar[1]*fin[12]+alphay[0]*fin[10]+alphavpar[0]*fin[9]+alphay[6]*fin[8]+fin[6]*alphay[8]+alphay[3]*fin[4]+fin[3]*alphay[4]); 
  out[15] += 0.3872983346207416*alphay[3]*fin[21]+0.3872983346207417*(alphay[6]*fin[19]+fin[10]*alphay[17])+0.3872983346207416*fin[13]*alphay[16]+0.4330127018922193*((alphavpar[3]+alphay[2]+alphax[1])*fin[15]+(alphavpar[6]+alphay[5]+alphax[0])*fin[14]+alphay[0]*fin[13]+alphavpar[0]*fin[12]+alphay[1]*fin[10]+alphavpar[1]*fin[9]+alphay[3]*fin[8]+fin[3]*alphay[8]+alphay[4]*fin[6]+fin[4]*alphay[6]); 
  out[16] += 0.8660254037844387*alphavpar[6]*fin[17]+0.8660254037844386*alphavpar[3]*fin[16]+0.9682458365518543*(alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3]); 
  out[17] += (0.8660254037844386*alphavpar[3]+0.4330127018922193*alphax[1])*fin[17]+(0.8660254037844387*alphavpar[6]+0.4330127018922194*alphax[0])*fin[16]+0.9682458365518543*(alphavpar[0]*fin[6]+fin[0]*alphavpar[6]+alphavpar[1]*fin[3]+fin[1]*alphavpar[3]); 
  out[18] += 0.4330127018922194*alphay[8]*fin[21]+(0.8660254037844387*alphavpar[6]+0.4330127018922194*alphay[5])*fin[20]+0.4330127018922193*alphay[4]*fin[19]+(0.8660254037844386*alphavpar[3]+0.4330127018922193*alphay[2])*fin[18]+0.276641667586244*alphay[17]*fin[17]+0.4330127018922193*(alphay[1]*fin[17]+fin[1]*alphay[17])+0.276641667586244*alphay[16]*fin[16]+0.4330127018922194*(alphay[0]*fin[16]+fin[0]*alphay[16])+0.9682458365518543*(alphavpar[1]*fin[11]+alphavpar[0]*fin[7])+0.3872983346207417*alphay[6]*fin[6]+0.9682458365518543*fin[5]*alphavpar[6]+0.3872983346207417*alphay[3]*fin[3]+0.9682458365518543*fin[2]*alphavpar[3]; 
  out[19] += 0.8660254037844387*alphavpar[6]*fin[21]+0.8660254037844386*alphavpar[3]*fin[19]+0.9682458365518543*(alphavpar[1]*fin[13]+alphavpar[0]*fin[10]+alphavpar[6]*fin[8]+alphavpar[3]*fin[4]); 
  out[20] += 0.4330127018922193*alphay[4]*fin[21]+(0.8660254037844386*alphavpar[3]+0.4330127018922193*(alphay[2]+alphax[1]))*fin[20]+0.4330127018922194*alphay[8]*fin[19]+(0.8660254037844387*alphavpar[6]+0.4330127018922194*(alphay[5]+alphax[0]))*fin[18]+(0.276641667586244*alphay[16]+0.4330127018922194*alphay[0])*fin[17]+(0.276641667586244*fin[16]+0.4330127018922194*fin[0])*alphay[17]+0.4330127018922193*(alphay[1]*fin[16]+fin[1]*alphay[16])+0.9682458365518543*(alphavpar[0]*fin[11]+alphavpar[1]*fin[7])+0.3872983346207416*(alphay[3]*fin[6]+fin[3]*alphay[6])+0.9682458365518543*(fin[2]*alphavpar[6]+alphavpar[3]*fin[5]); 
  out[21] += (0.8660254037844386*alphavpar[3]+0.4330127018922193*alphax[1])*fin[21]+(0.8660254037844387*alphavpar[6]+0.4330127018922194*alphax[0])*fin[19]+0.9682458365518543*(alphavpar[0]*fin[13]+alphavpar[1]*fin[10]+alphavpar[3]*fin[8]+fin[4]*alphavpar[6]); 
  out[22] += (0.8660254037844387*alphavpar[6]+0.4330127018922194*alphay[5])*fin[23]+(0.8660254037844386*alphavpar[3]+0.4330127018922193*alphay[2])*fin[22]+(0.276641667586244*alphay[17]+0.4330127018922193*alphay[1])*fin[21]+0.276641667586244*alphay[16]*fin[19]+0.4330127018922194*(alphay[0]*fin[19]+alphay[8]*fin[17]+fin[8]*alphay[17])+0.4330127018922193*(alphay[4]*fin[16]+fin[4]*alphay[16])+0.9682458365518543*(alphavpar[1]*fin[15]+alphavpar[0]*fin[14])+0.3872983346207416*alphay[6]*fin[13]+0.9682458365518543*alphavpar[6]*fin[12]+0.3872983346207416*alphay[3]*fin[10]+0.9682458365518543*alphavpar[3]*fin[9]; 
  out[23] += (0.8660254037844386*alphavpar[3]+0.4330127018922193*(alphay[2]+alphax[1]))*fin[23]+(0.8660254037844387*alphavpar[6]+0.4330127018922194*(alphay[5]+alphax[0]))*fin[22]+(0.276641667586244*alphay[16]+0.4330127018922194*alphay[0])*fin[21]+0.276641667586244*alphay[17]*fin[19]+0.4330127018922193*(alphay[1]*fin[19]+alphay[4]*fin[17]+fin[4]*alphay[17])+0.4330127018922194*(alphay[8]*fin[16]+fin[8]*alphay[16])+0.9682458365518543*(alphavpar[0]*fin[15]+alphavpar[1]*fin[14])+0.3872983346207417*alphay[3]*fin[13]+0.9682458365518543*alphavpar[3]*fin[12]+0.3872983346207417*alphay[6]*fin[10]+0.9682458365518543*alphavpar[6]*fin[9]; 

  return 0.; 
} 