#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH double gyrokinetic_vol_2x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double* apardot, const double *fin, double* GKYL_RESTRICT out) 
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

  double BstarYdBmag[24] = {0.}; 
  BstarYdBmag[0] = -(0.8660254037844386*rdx2*(2.0*jacobtot_inv[0]*b_z[1]*m_*wvpar+(2.0*apar[1]*b_z[1]*jacobtot_inv[1]+jacobtot_inv[0]*(apar[0]*b_z[1]+b_z[0]*apar[1]))*q_))/q_; 
  BstarYdBmag[1] = -(0.8660254037844386*rdx2*(2.0*b_z[1]*jacobtot_inv[1]*m_*wvpar+((apar[0]*b_z[1]+b_z[0]*apar[1])*jacobtot_inv[1]+2.0*jacobtot_inv[0]*apar[1]*b_z[1])*q_))/q_; 
  BstarYdBmag[2] = -0.8660254037844386*((2.0*b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[3]+jacobtot_inv[0]*b_z[1]*apar[2])*rdx2; 
  BstarYdBmag[3] = -(1.0*jacobtot_inv[0]*b_z[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[5] = -0.8660254037844386*((b_z[0]*jacobtot_inv[1]+2.0*jacobtot_inv[0]*b_z[1])*apar[3]+b_z[1]*jacobtot_inv[1]*apar[2])*rdx2; 
  BstarYdBmag[6] = -(1.0*b_z[1]*jacobtot_inv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 

  double alphax[24] = {0.}; 
  alphax[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[0]*hamil[3]*rdvpar2)/m_-(0.25*((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*hamil[5]+(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*hamil[2])*rdy2)/q_); 
  alphax[1] = 1.732050807568877*rdx2*((((-0.25*(b_z[0]*jacobtot_inv[0]*hamil[5]+(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*hamil[2]))-0.45*b_z[1]*jacobtot_inv[1]*hamil[5])*rdy2)/q_+(0.25*BstarXdBmag[1]*hamil[3]*rdvpar2)/m_); 
  alphax[3] = (0.9682458365518543*BstarXdBmag[0]*hamil[16]*rdvpar2*rdx2)/m_; 
  alphax[6] = (0.9682458365518543*BstarXdBmag[1]*hamil[16]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.5809475019311124*alphax[6]-0.3354101966249685*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.5809475019311124*alphax[6])+0.3354101966249685*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.5809475019311124*alphax[6]-0.3354101966249685*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.5809475019311124*alphax[6])+0.3354101966249685*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.5809475019311124*alphax[6]-0.3354101966249685*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.5809475019311124*alphax[6])+0.3354101966249685*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.5809475019311124*alphax[6]-0.3354101966249685*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.5809475019311124*alphax[6])+0.3354101966249685*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.5809475019311124*alphax[6])-0.3354101966249685*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.5809475019311124*alphax[6]+0.3354101966249685*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.5809475019311124*alphax[6])-0.3354101966249685*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.5809475019311124*alphax[6]+0.3354101966249685*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.5809475019311124*alphax[6])-0.3354101966249685*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.5809475019311124*alphax[6]+0.3354101966249685*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.5809475019311124*alphax[6])-0.3354101966249685*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.5809475019311124*alphax[6]+0.3354101966249685*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 

  double alphay[24] = {0.}; 
  alphay[0] = 0.4330127018922193*((hamil[1]*(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*rdx2)/q_+((2.23606797749979*BstarYdBmag[3]*hamil[16]+BstarYdBmag[0]*hamil[3])*rdvpar2)/m_)*rdy2; 
  alphay[1] = 0.4330127018922193*((hamil[1]*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*rdx2)/q_+((2.23606797749979*BstarYdBmag[6]*hamil[16]+BstarYdBmag[1]*hamil[3])*rdvpar2)/m_)*rdy2; 
  alphay[2] = 0.4330127018922193*(((b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*hamil[5]*rdx2)/q_+(BstarYdBmag[2]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[3] = (0.4330127018922193*(2.23606797749979*BstarYdBmag[0]*hamil[16]+BstarYdBmag[3]*hamil[3])*rdvpar2*rdy2)/m_; 
  alphay[4] = (0.4330127018922193*(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*hamil[8]*rdx2*rdy2)/q_; 
  alphay[5] = 0.4330127018922193*(((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*hamil[5]*rdx2)/q_+(hamil[3]*BstarYdBmag[5]*rdvpar2)/m_)*rdy2; 
  alphay[6] = (0.4330127018922193*(2.23606797749979*BstarYdBmag[1]*hamil[16]+hamil[3]*BstarYdBmag[6])*rdvpar2*rdy2)/m_; 
  alphay[7] = (0.9682458365518543*BstarYdBmag[2]*hamil[16]*rdvpar2*rdy2)/m_; 
  alphay[8] = (0.4330127018922193*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*hamil[8]*rdx2*rdy2)/q_; 
  alphay[11] = (0.9682458365518543*BstarYdBmag[5]*hamil[16]*rdvpar2*rdy2)/m_; 
  alphay[16] = (0.8660254037844386*BstarYdBmag[3]*hamil[16]*rdvpar2*rdy2)/m_; 
  alphay[17] = (0.8660254037844387*BstarYdBmag[6]*hamil[16]*rdvpar2*rdy2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*((-0.223606797749979*alphay[17])+0.223606797749979*alphay[16]-0.5809475019311124*alphay[11]+0.25*alphay[8]+0.5809475019311124*alphay[7]+0.3354101966249685*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.2795084971874738*alphay[17]-0.2795084971874737*alphay[16]+0.25*alphay[8]+0.4330127018922193*alphay[5]-0.25*alphay[4]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.223606797749979*alphay[17])+0.223606797749979*alphay[16]+0.5809475019311124*alphay[11]+0.25*alphay[8]-0.5809475019311124*alphay[7]-0.3354101966249685*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.223606797749979*alphay[17])+0.223606797749979*alphay[16]-0.5809475019311124*alphay[11]-0.25*alphay[8]+0.5809475019311124*alphay[7]+0.3354101966249685*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.2795084971874738*alphay[17]-0.2795084971874737*alphay[16]-0.25*alphay[8]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.223606797749979*alphay[17])+0.223606797749979*alphay[16]+0.5809475019311124*alphay[11]-0.25*alphay[8]-0.5809475019311124*alphay[7]-0.3354101966249685*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.223606797749979*alphay[17]+0.223606797749979*alphay[16]+0.5809475019311124*alphay[11]-0.25*alphay[8]+0.5809475019311124*alphay[7]-0.3354101966249685*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.2795084971874738*alphay[17])-0.2795084971874737*alphay[16]-0.25*alphay[8]-0.4330127018922193*alphay[5]-0.25*alphay[4]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.223606797749979*alphay[17]+0.223606797749979*alphay[16]-0.5809475019311124*alphay[11]-0.25*alphay[8]-0.5809475019311124*alphay[7]+0.3354101966249685*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.223606797749979*alphay[17]+0.223606797749979*alphay[16]+0.5809475019311124*alphay[11]+0.25*alphay[8]+0.5809475019311124*alphay[7]-0.3354101966249685*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.2795084971874738*alphay[17])-0.2795084971874737*alphay[16]+0.25*alphay[8]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.223606797749979*alphay[17]+0.223606797749979*alphay[16]-0.5809475019311124*alphay[11]+0.25*alphay[8]-0.5809475019311124*alphay[7]+0.3354101966249685*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -1.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.223606797749979*alphay[17])+0.223606797749979*alphay[16]+0.5809475019311124*alphay[11]+0.25*alphay[8]-0.5809475019311124*alphay[7]+0.3354101966249685*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.2795084971874738*alphay[17]-0.2795084971874737*alphay[16]+0.25*alphay[8]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.223606797749979*alphay[17])+0.223606797749979*alphay[16]-0.5809475019311124*alphay[11]+0.25*alphay[8]+0.5809475019311124*alphay[7]-0.3354101966249685*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.223606797749979*alphay[17])+0.223606797749979*alphay[16]+0.5809475019311124*alphay[11]-0.25*alphay[8]-0.5809475019311124*alphay[7]+0.3354101966249685*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.2795084971874738*alphay[17]-0.2795084971874737*alphay[16]-0.25*alphay[8]-0.4330127018922193*alphay[5]+0.25*alphay[4]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.223606797749979*alphay[17])+0.223606797749979*alphay[16]-0.5809475019311124*alphay[11]-0.25*alphay[8]+0.5809475019311124*alphay[7]-0.3354101966249685*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.223606797749979*alphay[17]+0.223606797749979*alphay[16]-0.5809475019311124*alphay[11]-0.25*alphay[8]-0.5809475019311124*alphay[7]-0.3354101966249685*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.2795084971874738*alphay[17])-0.2795084971874737*alphay[16]-0.25*alphay[8]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.223606797749979*alphay[17]+0.223606797749979*alphay[16]+0.5809475019311124*alphay[11]-0.25*alphay[8]+0.5809475019311124*alphay[7]+0.3354101966249685*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.223606797749979*alphay[17]+0.223606797749979*alphay[16]-0.5809475019311124*alphay[11]+0.25*alphay[8]-0.5809475019311124*alphay[7]-0.3354101966249685*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.2795084971874738*alphay[17])-0.2795084971874737*alphay[16]+0.25*alphay[8]+0.4330127018922193*alphay[5]+0.25*alphay[4]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.223606797749979*alphay[17]+0.223606797749979*alphay[16]+0.5809475019311124*alphay[11]+0.25*alphay[8]+0.5809475019311124*alphay[7]+0.3354101966249685*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 1.5*(alphaR+fabs(alphaR)); 

  double alphavpar[24] = {0.}; 
  alphavpar[0] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[1]*hamil[5]+BstarYdBmag[0]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[1]*rdx2))/m_; 
  alphavpar[1] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[0]*hamil[5]+BstarYdBmag[1]*hamil[2])*rdy2+BstarXdBmag[1]*hamil[1]*rdx2))/m_; 
  alphavpar[2] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[5]*hamil[5]+BstarYdBmag[2]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[5]*rdx2))/m_; 
  alphavpar[3] = -(0.4330127018922193*(hamil[5]*BstarYdBmag[6]+hamil[2]*BstarYdBmag[3])*rdvpar2*rdy2)/m_; 
  alphavpar[4] = -(0.4330127018922193*BstarXdBmag[0]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[2]*hamil[5]+hamil[2]*BstarYdBmag[5])*rdy2+BstarXdBmag[1]*hamil[5]*rdx2))/m_; 
  alphavpar[6] = -(0.4330127018922193*(hamil[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil[5])*rdvpar2*rdy2)/m_; 
  alphavpar[8] = -(0.4330127018922193*BstarXdBmag[1]*hamil[8]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  out[1] += 0.4330127018922193*(alphax[6]*fin[6]+alphax[3]*fin[3]+alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.4330127018922193*(alphay[17]*fin[17]+alphay[16]*fin[16]+alphay[11]*fin[11]+alphay[8]*fin[8]+alphay[7]*fin[7]+alphay[6]*fin[6]+alphay[5]*fin[5]+alphay[4]*fin[4]+alphay[3]*fin[3]+alphay[2]*fin[2]+alphay[1]*fin[1]+alphay[0]*fin[0]); 
  out[3] += 0.4330127018922193*(alphavpar[8]*fin[8]+alphavpar[6]*fin[6]+alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[5] += 0.05*(8.660254037844387*(alphay[16]*fin[17]+fin[16]*alphay[17])+8.660254037844386*((alphay[7]+alphax[6])*fin[11]+fin[7]*alphay[11]+alphay[4]*fin[8]+fin[4]*alphay[8]+alphax[3]*fin[7]+alphay[3]*fin[6]+fin[3]*alphay[6]+(alphay[2]+alphax[1])*fin[5]+fin[2]*(alphay[5]+alphax[0])+alphay[0]*fin[1]+fin[0]*alphay[1])); 
  out[6] += 0.05*(7.745966692414834*(alphax[6]*fin[17]+alphax[3]*fin[16])+8.660254037844386*(alphavpar[4]*fin[8]+fin[4]*alphavpar[8]+(alphavpar[3]+alphax[1])*fin[6]+fin[1]*alphax[6]+fin[3]*alphavpar[6]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphax[0]*fin[3]+fin[0]*alphax[3]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1])); 
  out[7] += 0.05*(7.745966692414834*(alphay[11]*fin[20]+alphay[7]*fin[18]+alphay[6]*fin[17]+fin[6]*alphay[17]+alphay[3]*fin[16]+fin[3]*alphay[16])+8.660254037844386*(alphay[8]*fin[13]+alphavpar[8]*fin[12]+(alphavpar[6]+alphay[5])*fin[11]+fin[5]*alphay[11]+alphay[4]*fin[10]+alphavpar[4]*fin[9]+(alphavpar[3]+alphay[2])*fin[7]+fin[2]*alphay[7]+alphay[1]*fin[6]+fin[1]*alphay[6]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphay[0]*fin[3]+fin[0]*alphay[3]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2])); 
  out[8] += 0.4330127018922193*(alphax[6]*fin[13]+alphax[3]*fin[10]+alphax[1]*fin[8]+alphax[0]*fin[4]); 
  out[9] += 0.05*(8.660254037844387*(alphay[17]*fin[21]+alphay[16]*fin[19])+8.660254037844386*(alphay[11]*fin[15]+alphay[7]*fin[14]+alphay[6]*fin[13]+alphay[5]*fin[12]+alphay[3]*fin[10]+alphay[2]*fin[9]+alphay[1]*fin[8]+fin[1]*alphay[8]+alphay[0]*fin[4]+fin[0]*alphay[4])); 
  out[10] += 0.4330127018922193*(alphavpar[6]*fin[13]+alphavpar[5]*fin[12]+alphavpar[3]*fin[10]+alphavpar[2]*fin[9]+alphavpar[1]*fin[8]+fin[1]*alphavpar[8]+alphavpar[0]*fin[4]+fin[0]*alphavpar[4]); 
  out[11] += 0.05*(7.745966692414834*((alphay[7]+alphax[6])*fin[20]+(alphay[11]+alphax[3])*fin[18]+alphay[3]*fin[17]+fin[3]*alphay[17]+alphay[6]*fin[16]+fin[6]*alphay[16])+8.660254037844386*(alphay[4]*fin[13]+alphavpar[4]*fin[12]+(alphavpar[3]+alphay[2]+alphax[1])*fin[11]+fin[2]*alphay[11]+alphay[8]*fin[10]+alphavpar[8]*fin[9]+(alphavpar[6]+alphay[5]+alphax[0])*fin[7]+fin[5]*alphay[7]+alphay[0]*fin[6]+fin[0]*alphay[6]+fin[5]*(alphax[6]+alphavpar[0])+fin[0]*alphavpar[5]+alphay[1]*fin[3]+fin[1]*alphay[3]+fin[2]*(alphax[3]+alphavpar[1])+fin[1]*alphavpar[2])); 
  out[12] += 0.4330127018922193*(alphay[16]*fin[21]+alphay[17]*fin[19]+(alphay[7]+alphax[6])*fin[15]+(alphay[11]+alphax[3])*fin[14]+alphay[3]*fin[13]+(alphay[2]+alphax[1])*fin[12]+alphay[6]*fin[10]+(alphay[5]+alphax[0])*fin[9]+alphay[0]*fin[8]+fin[0]*alphay[8]+alphay[1]*fin[4]+fin[1]*alphay[4]); 
  out[13] += 0.05*(7.745966692414834*(alphax[6]*fin[21]+alphax[3]*fin[19])+8.660254037844386*((alphavpar[3]+alphax[1])*fin[13]+alphavpar[2]*fin[12]+(alphavpar[6]+alphax[0])*fin[10]+alphavpar[5]*fin[9]+(alphax[6]+alphavpar[0])*fin[8]+fin[0]*alphavpar[8]+(alphax[3]+alphavpar[1])*fin[4]+fin[1]*alphavpar[4])); 
  out[14] += 0.05*(7.745966692414834*(alphay[11]*fin[23]+alphay[7]*fin[22]+alphay[6]*fin[21]+alphay[3]*fin[19]+fin[13]*alphay[17]+fin[10]*alphay[16])+8.660254037844386*((alphavpar[6]+alphay[5])*fin[15]+(alphavpar[3]+alphay[2])*fin[14]+alphay[1]*fin[13]+(alphay[11]+alphavpar[1])*fin[12]+alphay[0]*fin[10]+(alphay[7]+alphavpar[0])*fin[9]+(alphay[6]+alphavpar[5])*fin[8]+fin[6]*alphay[8]+fin[5]*alphavpar[8]+(alphay[3]+alphavpar[2])*fin[4]+fin[3]*alphay[4]+fin[2]*alphavpar[4])); 
  out[15] += 0.05*(7.745966692414834*((alphay[7]+alphax[6])*fin[23]+(alphay[11]+alphax[3])*fin[22]+alphay[3]*fin[21]+alphay[6]*fin[19]+fin[10]*alphay[17]+fin[13]*alphay[16])+8.660254037844386*((alphavpar[3]+alphay[2]+alphax[1])*fin[15]+(alphavpar[6]+alphay[5]+alphax[0])*fin[14]+alphay[0]*fin[13]+(alphay[7]+alphax[6]+alphavpar[0])*fin[12]+fin[9]*alphay[11]+alphay[1]*fin[10]+(alphax[3]+alphavpar[1])*fin[9]+(alphay[3]+alphavpar[2])*fin[8]+fin[3]*alphay[8]+fin[2]*alphavpar[8]+alphay[4]*fin[6]+fin[4]*alphay[6]+alphavpar[4]*fin[5]+fin[4]*alphavpar[5])); 
  out[16] += 0.05*(17.32050807568877*alphavpar[6]*fin[17]+17.32050807568877*alphavpar[3]*fin[16]+19.36491673103709*(alphavpar[8]*fin[13]+alphavpar[5]*fin[11]+alphavpar[4]*fin[10]+alphavpar[2]*fin[7]+alphavpar[1]*fin[6]+fin[1]*alphavpar[6]+alphavpar[0]*fin[3]+fin[0]*alphavpar[3])); 
  out[17] += 0.05*((17.32050807568877*alphavpar[3]+8.660254037844386*alphax[1])*fin[17]+(17.32050807568877*alphavpar[6]+8.660254037844387*alphax[0])*fin[16]+19.36491673103708*(alphavpar[4]*fin[13]+alphavpar[2]*fin[11]+alphavpar[8]*fin[10]+alphavpar[5]*fin[7])+7.745966692414834*alphax[6]*fin[6]+19.36491673103708*(alphavpar[0]*fin[6]+fin[0]*alphavpar[6])+7.745966692414834*alphax[3]*fin[3]+19.36491673103708*(alphavpar[1]*fin[3]+fin[1]*alphavpar[3])); 
  out[18] += 0.007142857142857143*(60.62177826491071*alphay[8]*fin[21]+(121.2435565298214*alphavpar[6]+60.62177826491071*alphay[5])*fin[20]+60.6217782649107*alphay[4]*fin[19]+(121.2435565298214*alphavpar[3]+60.6217782649107*alphay[2])*fin[18]+38.72983346207417*alphay[17]*fin[17]+60.6217782649107*(alphay[1]*fin[17]+fin[1]*alphay[17])+38.72983346207417*alphay[16]*fin[16]+60.62177826491071*(alphay[0]*fin[16]+fin[0]*alphay[16])+135.5544171172596*(alphavpar[8]*fin[15]+alphavpar[4]*fin[14])+(54.22176684690384*alphay[11]+135.5544171172596*alphavpar[1])*fin[11]+(54.22176684690384*alphay[7]+135.5544171172596*alphavpar[0])*fin[7]+54.22176684690384*alphay[6]*fin[6]+135.5544171172596*(alphavpar[5]*fin[6]+fin[5]*alphavpar[6])+54.22176684690384*alphay[3]*fin[3]+135.5544171172596*(alphavpar[2]*fin[3]+fin[2]*alphavpar[3])); 
  out[19] += 0.05*(17.32050807568877*alphavpar[6]*fin[21]+17.32050807568877*alphavpar[3]*fin[19]+19.36491673103708*(alphavpar[5]*fin[15]+alphavpar[2]*fin[14]+alphavpar[1]*fin[13]+alphavpar[0]*fin[10]+alphavpar[6]*fin[8]+fin[6]*alphavpar[8]+alphavpar[3]*fin[4]+fin[3]*alphavpar[4])); 
  out[20] += 0.007142857142857143*(60.6217782649107*alphay[4]*fin[21]+(121.2435565298214*alphavpar[3]+60.6217782649107*(alphay[2]+alphax[1]))*fin[20]+60.62177826491071*alphay[8]*fin[19]+(121.2435565298214*alphavpar[6]+60.62177826491071*(alphay[5]+alphax[0]))*fin[18]+(38.72983346207417*alphay[16]+60.62177826491071*alphay[0])*fin[17]+(38.72983346207417*fin[16]+60.62177826491071*fin[0])*alphay[17]+60.6217782649107*(alphay[1]*fin[16]+fin[1]*alphay[16])+135.5544171172596*(alphavpar[4]*fin[15]+alphavpar[8]*fin[14])+(54.22176684690384*(alphay[7]+alphax[6])+135.5544171172596*alphavpar[0])*fin[11]+fin[7]*(54.22176684690384*(alphay[11]+alphax[3])+135.5544171172596*alphavpar[1])+(54.22176684690384*alphay[3]+135.5544171172596*alphavpar[2])*fin[6]+54.22176684690384*fin[3]*alphay[6]+135.5544171172596*(fin[2]*alphavpar[6]+alphavpar[3]*fin[5]+fin[3]*alphavpar[5])); 
  out[21] += 0.05*((17.32050807568877*alphavpar[3]+8.660254037844386*alphax[1])*fin[21]+(17.32050807568877*alphavpar[6]+8.660254037844387*alphax[0])*fin[19]+19.36491673103709*(alphavpar[2]*fin[15]+alphavpar[5]*fin[14])+(7.745966692414834*alphax[6]+19.36491673103709*alphavpar[0])*fin[13]+7.745966692414834*alphax[3]*fin[10]+19.36491673103709*(alphavpar[1]*fin[10]+alphavpar[3]*fin[8]+fin[3]*alphavpar[8]+alphavpar[4]*fin[6]+fin[4]*alphavpar[6])); 
  out[22] += 0.007142857142857143*((121.2435565298214*alphavpar[6]+60.62177826491071*alphay[5])*fin[23]+(121.2435565298214*alphavpar[3]+60.6217782649107*alphay[2])*fin[22]+(38.72983346207417*alphay[17]+60.6217782649107*alphay[1])*fin[21]+38.72983346207417*alphay[16]*fin[19]+60.62177826491071*(alphay[0]*fin[19]+alphay[8]*fin[17]+fin[8]*alphay[17])+60.6217782649107*(alphay[4]*fin[16]+fin[4]*alphay[16])+(54.22176684690384*alphay[11]+135.5544171172596*alphavpar[1])*fin[15]+(54.22176684690384*alphay[7]+135.5544171172596*alphavpar[0])*fin[14]+54.22176684690384*alphay[6]*fin[13]+135.5544171172596*(alphavpar[5]*fin[13]+alphavpar[6]*fin[12]+alphavpar[8]*fin[11])+54.22176684690384*alphay[3]*fin[10]+135.5544171172596*(alphavpar[2]*fin[10]+alphavpar[3]*fin[9]+alphavpar[4]*fin[7])); 
  out[23] += 0.007142857142857143*((121.2435565298214*alphavpar[3]+60.6217782649107*(alphay[2]+alphax[1]))*fin[23]+(121.2435565298214*alphavpar[6]+60.62177826491071*(alphay[5]+alphax[0]))*fin[22]+(38.72983346207417*alphay[16]+60.62177826491071*alphay[0])*fin[21]+38.72983346207417*alphay[17]*fin[19]+60.6217782649107*(alphay[1]*fin[19]+alphay[4]*fin[17]+fin[4]*alphay[17])+60.62177826491071*(alphay[8]*fin[16]+fin[8]*alphay[16])+(54.22176684690384*(alphay[7]+alphax[6])+135.5544171172596*alphavpar[0])*fin[15]+(54.22176684690384*(alphay[11]+alphax[3])+135.5544171172596*alphavpar[1])*fin[14]+54.22176684690384*alphay[3]*fin[13]+135.5544171172596*(alphavpar[2]*fin[13]+alphavpar[3]*fin[12]+alphavpar[4]*fin[11])+54.22176684690384*alphay[6]*fin[10]+135.5544171172596*(alphavpar[5]*fin[10]+alphavpar[6]*fin[9]+fin[7]*alphavpar[8])); 

  return cflFreq; 
} 
GKYL_CU_DH double gyrokinetic_step2_vol_2x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *apardot, const double *fin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // apardot: time derivative of Apar.
  // fIn: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[2]; 
  out[3] += -(0.8660254037844386*(apardot[3]*fin[5]+apardot[2]*fin[2]+apardot[1]*fin[1]+apardot[0]*fin[0])*q_*rdvpar2)/m_; 
  out[6] += -(0.8660254037844386*(apardot[2]*fin[5]+fin[2]*apardot[3]+apardot[0]*fin[1]+fin[0]*apardot[1])*q_*rdvpar2)/m_; 
  out[7] += -(0.8660254037844386*(apardot[1]*fin[5]+fin[1]*apardot[3]+apardot[0]*fin[2]+fin[0]*apardot[2])*q_*rdvpar2)/m_; 
  out[10] += -(0.8660254037844386*(apardot[3]*fin[12]+apardot[2]*fin[9]+apardot[1]*fin[8]+apardot[0]*fin[4])*q_*rdvpar2)/m_; 
  out[11] += -(0.8660254037844386*(apardot[0]*fin[5]+fin[0]*apardot[3]+apardot[1]*fin[2]+fin[1]*apardot[2])*q_*rdvpar2)/m_; 
  out[13] += -(0.8660254037844386*(apardot[2]*fin[12]+apardot[3]*fin[9]+apardot[0]*fin[8]+apardot[1]*fin[4])*q_*rdvpar2)/m_; 
  out[14] += -(0.8660254037844386*(apardot[1]*fin[12]+apardot[0]*fin[9]+apardot[3]*fin[8]+apardot[2]*fin[4])*q_*rdvpar2)/m_; 
  out[15] += -(0.8660254037844386*(apardot[0]*fin[12]+apardot[1]*fin[9]+apardot[2]*fin[8]+apardot[3]*fin[4])*q_*rdvpar2)/m_; 
  out[16] += -(1.936491673103709*(apardot[3]*fin[11]+apardot[2]*fin[7]+apardot[1]*fin[6]+apardot[0]*fin[3])*q_*rdvpar2)/m_; 
  out[17] += -(1.936491673103709*(apardot[2]*fin[11]+apardot[3]*fin[7]+apardot[0]*fin[6]+apardot[1]*fin[3])*q_*rdvpar2)/m_; 
  out[18] += -(1.936491673103709*(apardot[1]*fin[11]+apardot[0]*fin[7]+apardot[3]*fin[6]+apardot[2]*fin[3])*q_*rdvpar2)/m_; 
  out[19] += -(1.936491673103709*(apardot[3]*fin[15]+apardot[2]*fin[14]+apardot[1]*fin[13]+apardot[0]*fin[10])*q_*rdvpar2)/m_; 
  out[20] += -(1.936491673103709*(apardot[0]*fin[11]+apardot[1]*fin[7]+apardot[2]*fin[6]+apardot[3]*fin[3])*q_*rdvpar2)/m_; 
  out[21] += -(1.936491673103709*(apardot[2]*fin[15]+apardot[3]*fin[14]+apardot[0]*fin[13]+apardot[1]*fin[10])*q_*rdvpar2)/m_; 
  out[22] += -(1.936491673103709*(apardot[1]*fin[15]+apardot[0]*fin[14]+apardot[3]*fin[13]+apardot[2]*fin[10])*q_*rdvpar2)/m_; 
  out[23] += -(1.936491673103709*(apardot[0]*fin[15]+apardot[1]*fin[14]+apardot[2]*fin[13]+apardot[3]*fin[10])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.0625*(0.5*apardot[3]-0.5*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*(0.5*apardot[3]-0.5*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*((-0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*((-0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*(0.5*(apardot[1]+apardot[0])-0.5*(apardot[3]+apardot[2]))*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*(0.5*(apardot[1]+apardot[0])-0.5*(apardot[3]+apardot[2]))*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.03125*(apardot[3]+apardot[2]+apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.03125*(apardot[3]+apardot[2]+apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -2.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.0625*(0.5*apardot[3]-0.5*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*(0.5*apardot[3]-0.5*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*((-0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*((-0.5*apardot[3])+0.5*apardot[2]-0.5*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*(0.5*(apardot[1]+apardot[0])-0.5*(apardot[3]+apardot[2]))*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*(0.5*(apardot[1]+apardot[0])-0.5*(apardot[3]+apardot[2]))*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.03125*(apardot[3]+apardot[2]+apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.03125*(apardot[3]+apardot[2]+apardot[1]+apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 2.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
