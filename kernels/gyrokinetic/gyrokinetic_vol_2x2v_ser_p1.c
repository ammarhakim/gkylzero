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

  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 

  double BstarXdBmag[16]; 
  BstarXdBmag[0] = 0.8660254037844386*((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[3]+(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[2])*rdy2; 
  BstarXdBmag[1] = 0.1732050807568877*((9.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0])*apar[3]+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[2])*rdy2; 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = -(0.8660254037844386*rdx2*(2.0*jacobtot_inv[0]*b_z[1]*m_*wvpar+(2.0*apar[1]*b_z[1]*jacobtot_inv[1]+jacobtot_inv[0]*(apar[0]*b_z[1]+b_z[0]*apar[1]))*q_))/q_; 
  BstarYdBmag[1] = -(0.8660254037844386*rdx2*(2.0*b_z[1]*jacobtot_inv[1]*m_*wvpar+((apar[0]*b_z[1]+b_z[0]*apar[1])*jacobtot_inv[1]+2.0*jacobtot_inv[0]*apar[1]*b_z[1])*q_))/q_; 
  BstarYdBmag[2] = -0.8660254037844386*((2.0*b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[3]+jacobtot_inv[0]*b_z[1]*apar[2])*rdx2; 
  BstarYdBmag[3] = -(1.0*jacobtot_inv[0]*b_z[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[5] = -0.8660254037844386*((b_z[0]*jacobtot_inv[1]+2.0*jacobtot_inv[0]*b_z[1])*apar[3]+b_z[1]*jacobtot_inv[1]*apar[2])*rdx2; 
  BstarYdBmag[6] = -(1.0*b_z[1]*jacobtot_inv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[0]*hamil[3]*rdvpar2)/m_-(0.25*((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*hamil[5]+(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*hamil[2])*rdy2)/q_); 
  alphax[1] = 1.732050807568877*rdx2*((((-0.25*(b_z[0]*jacobtot_inv[0]*hamil[5]+(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*hamil[2]))-0.45*b_z[1]*jacobtot_inv[1]*hamil[5])*rdy2)/q_+(0.25*BstarXdBmag[1]*hamil[3]*rdvpar2)/m_); 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  double alphay[16]; 
  alphay[0] = 0.4330127018922193*((hamil[1]*(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*rdx2)/q_+(BstarYdBmag[0]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[1] = 0.4330127018922193*((hamil[1]*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*rdx2)/q_+(BstarYdBmag[1]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[2] = 0.4330127018922193*(((b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*hamil[5]*rdx2)/q_+(BstarYdBmag[2]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[3] = (0.4330127018922193*BstarYdBmag[3]*hamil[3]*rdvpar2*rdy2)/m_; 
  alphay[4] = (0.4330127018922193*(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*hamil[8]*rdx2*rdy2)/q_; 
  alphay[5] = 0.4330127018922193*(((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*hamil[5]*rdx2)/q_+(hamil[3]*BstarYdBmag[5]*rdvpar2)/m_)*rdy2; 
  alphay[6] = (0.4330127018922193*hamil[3]*BstarYdBmag[6]*rdvpar2*rdy2)/m_; 
  alphay[8] = (0.4330127018922193*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*hamil[8]*rdx2*rdy2)/q_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  double alphavpar[16]; 
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
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  out[1] += 0.4330127018922193*(alphax[1]*fin[1]+alphax[0]*fin[0]); 
  out[2] += 0.4330127018922193*(alphay[8]*fin[8]+alphay[6]*fin[6]+alphay[5]*fin[5]+alphay[4]*fin[4]+alphay[3]*fin[3]+alphay[2]*fin[2]+alphay[1]*fin[1]+alphay[0]*fin[0]); 
  out[3] += 0.4330127018922193*(alphavpar[8]*fin[8]+alphavpar[6]*fin[6]+alphavpar[5]*fin[5]+alphavpar[4]*fin[4]+alphavpar[3]*fin[3]+alphavpar[2]*fin[2]+alphavpar[1]*fin[1]+alphavpar[0]*fin[0]); 
  out[5] += 0.4330127018922193*(alphay[4]*fin[8]+fin[4]*alphay[8]+alphay[3]*fin[6]+fin[3]*alphay[6]+(alphay[2]+alphax[1])*fin[5]+fin[2]*(alphay[5]+alphax[0])+alphay[0]*fin[1]+fin[0]*alphay[1]); 
  out[6] += 0.4330127018922193*(alphavpar[4]*fin[8]+fin[4]*alphavpar[8]+(alphavpar[3]+alphax[1])*fin[6]+fin[3]*alphavpar[6]+alphavpar[2]*fin[5]+fin[2]*alphavpar[5]+alphax[0]*fin[3]+alphavpar[0]*fin[1]+fin[0]*alphavpar[1]); 
  out[7] += 0.4330127018922193*(alphay[8]*fin[13]+alphavpar[8]*fin[12]+(alphavpar[6]+alphay[5])*fin[11]+alphay[4]*fin[10]+alphavpar[4]*fin[9]+(alphavpar[3]+alphay[2])*fin[7]+alphay[1]*fin[6]+fin[1]*alphay[6]+alphavpar[1]*fin[5]+fin[1]*alphavpar[5]+alphay[0]*fin[3]+fin[0]*alphay[3]+alphavpar[0]*fin[2]+fin[0]*alphavpar[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*fin[8]+alphax[0]*fin[4]); 
  out[9] += 0.4330127018922193*(alphay[6]*fin[13]+alphay[5]*fin[12]+alphay[3]*fin[10]+alphay[2]*fin[9]+alphay[1]*fin[8]+fin[1]*alphay[8]+alphay[0]*fin[4]+fin[0]*alphay[4]); 
  out[10] += 0.4330127018922193*(alphavpar[6]*fin[13]+alphavpar[5]*fin[12]+alphavpar[3]*fin[10]+alphavpar[2]*fin[9]+alphavpar[1]*fin[8]+fin[1]*alphavpar[8]+alphavpar[0]*fin[4]+fin[0]*alphavpar[4]); 
  out[11] += 0.4330127018922193*(alphay[4]*fin[13]+alphavpar[4]*fin[12]+(alphavpar[3]+alphay[2]+alphax[1])*fin[11]+alphay[8]*fin[10]+alphavpar[8]*fin[9]+(alphavpar[6]+alphay[5]+alphax[0])*fin[7]+alphay[0]*fin[6]+fin[0]*alphay[6]+alphavpar[0]*fin[5]+fin[0]*alphavpar[5]+alphay[1]*fin[3]+fin[1]*alphay[3]+alphavpar[1]*fin[2]+fin[1]*alphavpar[2]); 
  out[12] += 0.4330127018922193*(alphay[3]*fin[13]+(alphay[2]+alphax[1])*fin[12]+alphay[6]*fin[10]+(alphay[5]+alphax[0])*fin[9]+alphay[0]*fin[8]+fin[0]*alphay[8]+alphay[1]*fin[4]+fin[1]*alphay[4]); 
  out[13] += 0.4330127018922193*((alphavpar[3]+alphax[1])*fin[13]+alphavpar[2]*fin[12]+(alphavpar[6]+alphax[0])*fin[10]+alphavpar[5]*fin[9]+alphavpar[0]*fin[8]+fin[0]*alphavpar[8]+alphavpar[1]*fin[4]+fin[1]*alphavpar[4]); 
  out[14] += 0.4330127018922193*((alphavpar[6]+alphay[5])*fin[15]+(alphavpar[3]+alphay[2])*fin[14]+alphay[1]*fin[13]+alphavpar[1]*fin[12]+alphay[0]*fin[10]+alphavpar[0]*fin[9]+(alphay[6]+alphavpar[5])*fin[8]+fin[6]*alphay[8]+fin[5]*alphavpar[8]+(alphay[3]+alphavpar[2])*fin[4]+fin[3]*alphay[4]+fin[2]*alphavpar[4]); 
  out[15] += 0.4330127018922193*((alphavpar[3]+alphay[2]+alphax[1])*fin[15]+(alphavpar[6]+alphay[5]+alphax[0])*fin[14]+alphay[0]*fin[13]+alphavpar[0]*fin[12]+alphay[1]*fin[10]+alphavpar[1]*fin[9]+(alphay[3]+alphavpar[2])*fin[8]+fin[3]*alphay[8]+fin[2]*alphavpar[8]+alphay[4]*fin[6]+fin[4]*alphay[6]+alphavpar[4]*fin[5]+fin[4]*alphavpar[5]); 

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
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.0625*(0.5*apardot[3]-0.4999999999999999*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*(0.5*apardot[3]-0.4999999999999999*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*((-0.5*apardot[3])+0.4999999999999999*apardot[2]-0.4999999999999999*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*((-0.5*apardot[3])+0.4999999999999999*apardot[2]-0.4999999999999999*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*((-0.5*apardot[3])-0.4999999999999999*apardot[2]+0.4999999999999999*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*((-0.5*apardot[3])-0.4999999999999999*apardot[2]+0.4999999999999999*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*(0.5*apardot[3]+0.4999999999999999*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  alphaL = -(0.0625*(0.5*apardot[3]+0.4999999999999999*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-fabs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.0625*(0.5*apardot[3]-0.4999999999999999*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*(0.5*apardot[3]-0.4999999999999999*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*((-0.5*apardot[3])+0.4999999999999999*apardot[2]-0.4999999999999999*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*((-0.5*apardot[3])+0.4999999999999999*apardot[2]-0.4999999999999999*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*((-0.5*apardot[3])-0.4999999999999999*apardot[2]+0.4999999999999999*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*((-0.5*apardot[3])-0.4999999999999999*apardot[2]+0.4999999999999999*apardot[1]+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*(0.5*apardot[3]+0.4999999999999999*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 
  alphaR = -(0.0625*(0.5*apardot[3]+0.4999999999999999*(apardot[2]+apardot[1])+0.5*apardot[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+fabs(alphaR)); 

  return cflFreq; 
} 
